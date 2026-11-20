#pragma once

// Loi gaussienne discrète sur ℤ : P(X = k) = e^{−(k − μ)²/(2σ²)}/Z, Z = Σ_j e^{−(j − μ)²/(2σ²)},
// μ réel, σ > 0.
//
//   wsl::random::discrete_gaussian X({.location = 0.3, .scale = 2.5});
//   X.pmf(1); X.cdf(-2); X.quantile(0.9); X.mean(); X(g);
//
// Z : σ >= 1, sommation de Poisson, σ√(2π) (1 + 2 Σ_{m>=1} e^{−2π²σ²m²} cos 2πmμ) (deux
// termes suffisent : e^{−8π²} < 2^-113) ; σ < 1, somme directe. Moyenne et variance de même
// (μ − (4πσ³√(2π)/Z) Σ m e^{−2π²σ²m²} sin 2πmμ ; Σ (j − μ)² e^{…} par sa transformée).
// Répartition : la queue du côté de k, sommée depuis k (termes décroissants), en
// double-double, O(σ) termes ; l'autre par complément. Quantile exact (detail/discrete.hpp).
//
// Tirage (approché, bit-exact) : rejet depuis une loi de Laplace discrète centrée en
// c = ⌊μ⌉, d'échelle t = ⌊σ⌋ + 1 (Canonne, Kamath et Steinke, 2020), acceptation
// e^{−(k − μ)²/(2σ²) + |k − c|/t − σ²/(2t²) − |μ − c|/t} ; taux d'acceptation >= 0,4.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/trigonometry.hpp>
#include <wsl/random/detail/discrete.hpp>
#include <wsl/random/geometric.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    namespace detail
    {
    }

    template <law_integer I = std::int64_t>
    class discrete_gaussian : public law_interface
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            double location = 0, scale = 1;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.location = location(), .scale = scale()}; }

        constexpr discrete_gaussian() : discrete_gaussian(parameters{}) {}

        constexpr discrete_gaussian(parameters q) : mu_(q.location), sigma_(q.scale)
        {
            if (!detail::is_finite(mu_) || !(detail::magnitude(mu_) < 0x1p52) || !(sigma_ > 0) || !(sigma_ < 0x1p40))
                throw std::domain_error("wsl::random::discrete_gaussian : il faut |location| < 2^52 et 0 < scale < 2^40");
            namespace sd = detail::sd;
            auto const [s2, s2_error] = sd::two_product(sigma_, sigma_);
            two_s2_ = {sd::fmul(s2, 2.0), sd::fmul(s2_error, 2.0)};
            center_ = static_cast<I>(mu_ < 0 ? -detail::floor_small(0.5 - mu_) : detail::floor_small(mu_ + 0.5));
            moments m = sums();
            log_z_ = detail::log(m.z);
            mean_ = m.mean;
            variance_ = m.variance;
        }

        constexpr double location() const noexcept { return mu_; }
        constexpr double scale() const noexcept { return sigma_; }
        constexpr interval<I> support() const noexcept { return {std::numeric_limits<I>::min(), std::numeric_limits<I>::max()}; }
        constexpr double mean() const noexcept { return mean_; }
        constexpr double variance() const noexcept { return variance_; }
        constexpr I mode() const noexcept { return center_; }

        constexpr double pmf(I k) const noexcept
        {
            auto const& X = *this;
            auto const l = X.log_mass(k);
            if (!(l.hi > -1500))
                return 0;
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }
        constexpr double log_pmf(I k) const noexcept
        {
            auto const& X = *this;
            auto const l = X.log_mass(k);
            return l.hi + l.lo;
        }
        constexpr double cdf(I k) const noexcept
        {
            auto const& X = *this;
            auto const [v, is_f] = X.distribution(k);
            return is_f ? detail::rounded(v) : detail::sd::one_minus(v);
        }
        constexpr double ccdf(I k) const noexcept
        {
            auto const& X = *this;
            auto const [v, is_f] = X.distribution(k);
            return is_f ? detail::sd::one_minus(v) : detail::rounded(v);
        }

        // Plus petit k tel que P(X <= k) >= u.
        constexpr I quantile(double u) const
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            if (!(u >= 0 && u <= 1))
                throw std::domain_error("wsl::random::quantile : p doit être dans [0, 1]");
            constexpr I low = std::numeric_limits<I>::min(), high = std::numeric_limits<I>::max();
            if (u == 0)
                return low;
            if (u == 1)
                return high;
            double const guess = X.mu_ + sd::fmul(detail::normal_quantile_guess(u), X.sigma_);
            I const span = static_cast<I>(sd::fmul(40.0, X.sigma_)) + 40;   // hors de μ ± 40σ + 40 : masse < 2^-1100
            I const lo = X.center_ - span, hi = X.center_ + span;
            auto const ge = [&](I k) {
                auto const [v, is_f] = X.distribution(k);
                return detail::at_least(v, is_f, u);
            };
            return detail::smallest_true(static_cast<I>(guess < double(lo) ? double(lo) : guess > double(hi) ? double(hi) : guess), lo, hi, ge);
        }

        template <std::uniform_random_bit_generator G>
        constexpr I sample(G& g) const
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            double const t = detail::floor_small(X.sigma_) + 1;
            double const delta = X.mu_ - double(X.center_);                 // |δ| <= ½, exact
            double const shift = sd::fmul(X.sigma_, X.sigma_) / sd::fmul(2.0, sd::fmul(t, t)) + detail::magnitude(delta) / t;
            geometric<I> const magnitude({.p = -wsl::special::expm1(-1.0 / t)});
            for (;;)
            {
                bool const negative = bernoulli_trial(0.5, g);
                I const d = wsl::random::sample(magnitude, g);
                if (negative && d == 0)
                    continue;
                double const dd = negative ? -double(d) : double(d);
                double const e = dd - delta;
                double const exponent = (double(d) / t - shift) - sd::fmul(e, e) / sd::fmul(2.0, sd::fmul(X.sigma_, X.sigma_));
                double const accept = wsl::special::exp(exponent);
                if (bernoulli_trial(accept < 1 ? accept : 1.0, g))
                    return X.center_ + (negative ? -d : d);
            }
        }

        friend constexpr bool operator==(discrete_gaussian const& a, discrete_gaussian const& b) noexcept { return a.mu_ == b.mu_ && a.sigma_ == b.sigma_; }

    private:
        // −(k − μ)²/(2σ²), en double-double
        constexpr detail::sd::double_double exponent(I k) const noexcept
        {
            auto const [d, d_error] = detail::sd::two_sum(double(k), -mu_);
            auto const sq = detail::mul({d, d_error}, {d, d_error});
            auto const q = detail::sd::quotient_dd(sq, two_s2_);
            return {-q.hi, -q.lo};
        }
        constexpr detail::sd::double_double log_mass(I k) const noexcept { return detail::sum(exponent(k), {-log_z_.hi, -log_z_.lo}); }

        struct moments
        {
            detail::scaled z;
            double mean, variance;
        };
        constexpr moments sums() const noexcept
        {
            namespace sd = detail::sd;
            if (sigma_ >= 1)                                                // Poisson : m = 1, 2
            {
                double const pi = detail::pi_hi;
                double c = 0, sn = 0, v = 0;                                // Σ 2 e^{…} cos, Σ m e^{…} sin, Σ 2 (1 − 4π²σ²m²) e^{…} cos
                for (int m = 1; m <= 2; ++m)
                {
                    double const md = m;
                    double const w = sd::fmul(sd::fmul(2.0, sd::fmul(pi, pi)), sd::fmul(sd::fmul(sigma_, sigma_), sd::fmul(md, md)));
                    if (w > 800)
                        break;
                    double const e = wsl::special::exp(-w);
                    double const f = mu_ - detail::floor_small(mu_);       // exact
                    double const angle = sd::fmul(sd::fmul(2.0, pi), sd::fmul(md, f));
                    double const cs = wsl::special::cos(angle), si = wsl::special::sin(angle);
                    c += sd::fmul(2.0, sd::fmul(e, cs));
                    sn += sd::fmul(md, sd::fmul(e, si));
                    v += sd::fmul(2.0, sd::fmul(1.0 - sd::fmul(2.0, w), sd::fmul(e, cs)));
                }
                detail::scaled const base = detail::normalized(sigma_, 0.0) * detail::normalized(detail::sqrt2_hi, detail::sqrt2_lo)
                                          * detail::sqrt(detail::normalized(detail::pi_hi, detail::pi_lo));       // σ√(2π)
                auto const [one_c, one_c_error] = sd::fast_two_sum(1.0, c);
                detail::scaled const z = base * detail::normalized(one_c, one_c_error);
                // E[X] − μ = −4πσ² Σ m e^{…} sin/(1 + c) ; E[(X − μ)²] = σ² (1 + v)/(1 + c)
                double const shift = -sd::fmul(sd::fmul(4.0, pi), sd::fmul(sd::fmul(sigma_, sigma_), sn)) / (1 + c);
                double const second = detail::rounded(detail::normalized(sigma_, 0.0) * detail::normalized(sigma_, 0.0)
                                                      * detail::normalized(1 + v, 0.0) / detail::normalized(one_c, one_c_error));
                return {z, mu_ + shift, second - sd::fmul(shift, shift)};
            }
            // somme directe autour de c
            detail::sd::double_double z{0.0, 0.0}, first{0.0, 0.0}, second{0.0, 0.0};
            for (I j = center_ - 60; j <= center_ + 60; ++j)
            {
                auto const e = exponent(j);
                if (!(e.hi > -1500))
                    continue;
                auto const w = detail::unscaled(detail::exp_scaled(e.hi, e.lo));
                auto const [d, d_error] = sd::two_sum(double(j), -mu_);
                z = detail::sum(z, w);
                first = detail::sum(first, detail::mul(w, {d, d_error}));
                second = detail::sum(second, detail::mul(w, detail::mul({d, d_error}, {d, d_error})));
            }
            auto const m1 = sd::quotient_dd(first, z);
            auto const m2 = sd::quotient_dd(second, z);
            auto const var = detail::sum(m2, detail::negated(detail::mul(m1, m1)));
            auto const mean = detail::sum({mu_, 0.0}, m1);
            return {detail::normalized(z.hi, z.lo), mean.hi + mean.lo, var.hi + var.lo};
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
            I const span = static_cast<I>(detail::sd::fmul(40.0, sigma_)) + 40;   // hors de c ± span : masse < 2^-1100
            if (lower ? k < center_ - span : k > center_ + span)
                return {{0.0, 0.0, 0}, lower};
            I j = lower ? k : k + 1;
            auto const base = exponent(j);
            detail::sd::double_double total{0.0, 0.0};
            for (std::int64_t step = 0; step < (std::int64_t{1} << 40); ++step)
            {
                auto const e = detail::sum(exponent(j), {-base.hi, -base.lo});   // relatif au premier terme
                if (!(e.hi > -1500))
                    break;
                auto const w = detail::unscaled(detail::exp_scaled(e.hi, e.lo));
                total = detail::sum(total, w);
                if (step > 0 && w.hi < detail::sd::fmul(total.hi, 0x1p-112))
                    break;
                j = lower ? j - 1 : j + 1;
            }
            auto const l = detail::sum(base, {-log_z_.hi, -log_z_.lo});
            if (!(l.hi > -1500))
                return {{0.0, 0.0, 0}, lower};
            return {detail::exp_scaled(l.hi, l.lo) * detail::normalized(total.hi, total.lo), lower};
        }

        double mu_, sigma_;
        detail::sd::double_double two_s2_{};
        I center_ = 0;
        detail::sd::double_double log_z_{};
        double mean_ = 0, variance_ = 0;
    };
}
