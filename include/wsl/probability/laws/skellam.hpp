#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Loi de Skellam : K = N₁ − N₂, N₁ ~ Poisson(μ₁) et N₂ ~ Poisson(μ₂) indépendantes, sur ℤ.
// P(K = k) = e^{−(μ₁ + μ₂)} (μ₁/μ₂)^{k/2} I_|k|(2√(μ₁μ₂)).
//
//   wsl::probability::laws::skellam K({.mean1 = 3.0, .mean2 = 1.5});
//   K.pmf(1); K.cdf(-2); K.quantile(0.9); K.mean(); K(g);
//
// Masse en logarithme double-double : −(μ₁ + μ₂) + (k/2) log(μ₁/μ₂) + log I_|k|(z) (detail de
// special_functions) ; z = 2√(μ₁μ₂) n'est pas un double en général : log I_n(z_hi + z_lo) =
// log I_n(z_hi) + z_lo (I_{n+1}/I_n + n/z) au premier ordre (le second est sous 2^-100).
// Répartition : la queue du côté de k sommée depuis k (termes décroissants au-delà du
// mode), en double-double ; l'autre par complément. Quantile exact (detail/discrete.hpp).
// Tirage : N₁ − N₂, deux tirages de Poisson (même exactitude que poisson).

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/discrete.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/laws/poisson.hpp>
#include <wsl/random/detail/portable_math.hpp>

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_integer I = std::int64_t>
    class skellam : public law_interface
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = poisson<I>::sampling;

        struct parameters
        {
            double mean1, mean2;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.mean1 = mu1_, .mean2 = mu2_}; }

        constexpr skellam(parameters p) : mu1_(p.mean1), mu2_(p.mean2)
        {
            namespace sd = detail::sd;
            if (!(mu1_ > 0) || !detail::is_finite(mu1_) || !(mu2_ > 0) || !detail::is_finite(mu2_))
                throw std::domain_error("wsl::probability::laws::skellam : il faut mean1 > 0 et mean2 > 0, finis");
            auto const [s, s_error] = sd::two_sum(mu1_, mu2_);
            minus_s_ = {-s, -s_error};
            // r = ½ log(μ₁/μ₂), z = 2√(μ₁μ₂) en double-double
            auto const l = detail::sum(detail::log_precise(mu1_), detail::negated(detail::log_precise(mu2_)));
            half_log_ratio_ = {sd::fmul(l.hi, 0.5), sd::fmul(l.lo, 0.5)};
            auto const z = detail::unscaled(detail::sqrt(detail::normalized(mu1_, 0.0) * detail::normalized(mu2_, 0.0)));
            z_hi_ = sd::fmul(z.hi, 2.0);
            z_lo_ = sd::fmul(z.lo, 2.0);
            auto const [d, d_error] = sd::two_sum(mu1_, -mu2_);
            center_ = static_cast<I>(detail::floor_small(d + d_error));
        }

        constexpr double mean1() const noexcept { return mu1_; }
        constexpr double mean2() const noexcept { return mu2_; }
        constexpr interval<I> support() const noexcept { return {std::numeric_limits<I>::min(), std::numeric_limits<I>::max()}; }

        constexpr double mean() const noexcept { return mu1_ - mu2_; }
        constexpr double variance() const noexcept { return mu1_ + mu2_; }
        // (μ₁ − μ₂)/(μ₁ + μ₂)^{3/2}
        constexpr double skewness() const noexcept
        {
            namespace sd = detail::sd;
            auto const [d, d_error] = sd::two_sum(mu1_, -mu2_);
            if (d == 0)
                return 0;
            bool const negative = d < 0;
            detail::scaled const s = detail::normalized(-minus_s_.hi, -minus_s_.lo);
            double const r = detail::rounded(detail::normalized(negative ? -d : d, negative ? -d_error : d_error) / (s * detail::sqrt(s)));
            return negative ? -r : r;
        }
        constexpr double excess_kurtosis() const noexcept
        {
            return detail::rounded(detail::normalized(1.0, 0.0) / detail::normalized(-minus_s_.hi, -minus_s_.lo));
        }
        constexpr I median() const noexcept { return quantile(0.5); }

        constexpr double pmf(I k) const noexcept
        {
            auto const l = log_mass(k);
            return l.hi > -1500 ? detail::rounded(detail::exp_scaled(l.hi, l.lo)) : 0.0;
        }
        constexpr double log_pmf(I k) const noexcept
        {
            auto const l = log_mass(k);
            return l.hi + l.lo;
        }
        constexpr double cdf(I k) const noexcept
        {
            auto const [v, is_f] = distribution(k);
            return is_f ? detail::rounded(v) : detail::sd::one_minus(v);
        }
        constexpr double ccdf(I k) const noexcept
        {
            auto const [v, is_f] = distribution(k);
            return is_f ? detail::sd::one_minus(v) : detail::rounded(v);
        }

        constexpr I quantile(double u) const
        {
            auto const& X = *this;
            if (!(u >= 0 && u <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            constexpr I low = std::numeric_limits<I>::min(), high = std::numeric_limits<I>::max();
            if (u == 0)
                return low;
            if (u == 1)
                return high;
            double const sigma = wsl::random::detail::pm::sqrt(-X.minus_s_.hi);
            double const guess = (X.mu1_ - X.mu2_) + detail::sd::fmul(detail::normal_quantile_guess(u), sigma);
            I const span = static_cast<I>(detail::sd::fmul(60.0, sigma)) + 60;   // queues : décroissance au moins géométrique
            I const lo = X.center_ - span * 40, hi = X.center_ + span * 40;
            auto const ge = [&](I k) {
                auto const [v, is_f] = X.distribution(k);
                return detail::at_least(v, is_f, u);
            };
            return detail::smallest_true(static_cast<I>(guess < double(lo) ? double(lo) : guess > double(hi) ? double(hi) : guess), lo, hi, ge);
        }
        constexpr I cquantile(double q) const
        {
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::probability::cquantile : q doit être dans [0, 1]");
            return quantile(1 - q);
        }

        template <std::uniform_random_bit_generator G>
        constexpr I sample(G& g) const
        {
            auto const& X = *this;
            I const a = detail::native_sample(poisson<I>({.mean = X.mu1_}), g);
            I const b = detail::native_sample(poisson<I>({.mean = X.mu2_}), g);
            return a - b;
        }
        template <std::uniform_random_bit_generator G>
        constexpr I sample(G& g, approximate_t) const
        {
            auto const& X = *this;
            I const a = wsl::probability::sample(poisson<I>({.mean = X.mu1_}), g, wsl::approximate);
            I const b = wsl::probability::sample(poisson<I>({.mean = X.mu2_}), g, wsl::approximate);
            return a - b;
        }

        friend constexpr bool operator==(skellam const& a, skellam const& b) noexcept { return a.mu1_ == b.mu1_ && a.mu2_ == b.mu2_; }

    private:
        // log I_n(z), n >= 0 entier ; mémoire des trois derniers (les queues avancent d'un pas).
        struct bessel_cache
        {
            double n[3] = {-1, -1, -1};
            detail::sd::double_double v[3]{};
            int next = 0;
        };
        constexpr detail::sd::double_double log_bessel_at_hi(double n, bessel_cache& c) const noexcept
        {
            for (int i = 0; i < 3; ++i)
                if (c.n[i] == n)
                    return c.v[i];
            auto const v = wsl::random::detail::sd::log_cyl_bessel_i(n, z_hi_);
            c.n[c.next] = n;
            c.v[c.next] = v;
            c.next = (c.next + 1) % 3;
            return v;
        }
        // log I_n(z_hi + z_lo) : correction du premier ordre par la dérivée logarithmique.
        constexpr detail::sd::double_double log_bessel(double n, bessel_cache& c) const noexcept
        {
            namespace sd = detail::sd;
            auto const l = log_bessel_at_hi(n, c);
            if (z_lo_ == 0)
                return l;
            auto const up = log_bessel_at_hi(n + 1, c);
            double const ratio = wsl::random::detail::pm::exp((up.hi - l.hi) + (up.lo - l.lo));   // I_{n+1}/I_n
            double const derivative = ratio + n / z_hi_;
            return detail::sum(l, {sd::fmul(z_lo_, derivative), 0.0});
        }
        constexpr detail::sd::double_double log_mass(I k, bessel_cache& c) const noexcept
        {
            double const n = k < 0 ? -double(k) : double(k);
            auto l = log_bessel(n, c);
            l = detail::sum(l, minus_s_);
            return detail::sum(l, detail::mul({double(k), 0.0}, half_log_ratio_));
        }
        constexpr detail::sd::double_double log_mass(I k) const noexcept
        {
            bessel_cache c;
            return log_mass(k, c);
        }

        // F(k) (is_f) si k < c, sinon 1 − F(k) : la queue sommée depuis k, termes décroissants.
        struct part
        {
            detail::scaled value;
            bool is_f;
        };
        constexpr part distribution(I k) const noexcept
        {
            bool const lower = k < center_;
            if (lower ? k == std::numeric_limits<I>::min() : k == std::numeric_limits<I>::max())
                return {{0.0, 0.0, 0}, lower};
            I j = lower ? k : k + 1;
            bessel_cache c;
            auto const base = log_mass(j, c);
            if (!(base.hi > -1500))                                         // queue sous 2^-2000
                return {{0.0, 0.0, 0}, lower};
            detail::sd::double_double total{0.0, 0.0};
            for (std::int64_t step = 0; step < (std::int64_t{1} << 40); ++step)
            {
                auto const e = detail::sum(log_mass(j, c), {-base.hi, -base.lo});   // relatif au premier terme
                if (!(e.hi > -1500))
                    break;
                auto const w = detail::unscaled(detail::exp_scaled(e.hi, e.lo));
                total = detail::sum(total, w);
                bool const away = lower ? j <= center_ : j > center_;
                if (step > 0 && away && w.hi < detail::sd::fmul(total.hi, 0x1p-112))
                    break;
                if (lower ? j == std::numeric_limits<I>::min() : j == std::numeric_limits<I>::max())
                    break;
                j = lower ? j - 1 : j + 1;
            }
            return {detail::exp_scaled(base.hi, base.lo) * detail::normalized(total.hi, total.lo), lower};
        }

        double mu1_, mu2_;
        detail::sd::double_double minus_s_{}, half_log_ratio_{};
        double z_hi_ = 0, z_lo_ = 0;
        I center_ = 0;
    };

}
