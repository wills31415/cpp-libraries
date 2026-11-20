#pragma once

// Lois de Gumbel (valeurs extrêmes de type I, maximum) et de Fréchet (type II).
//
//   wsl::random::gumbel G({.location = 0, .scale = 1});               // P(X <= x) = e^{−e^{−(x−μ)/β}}
//   wsl::random::frechet F({.location = 0, .scale = 2, .shape = 3});   // P(X <= x) = e^{−((x−m)/s)^{−α}}, x > m
//   G.cdf(1.0); F.quantile(0.99); G.mean(); F(g);
//
// Noyau commun : y = (x − μ)/β (Gumbel) ou α log((x − m)/s) (Fréchet, log du rapport par
// log1p près de 1), en double-double ; t = e^{−y} ; répartition e^{−t}, survie −expm1(−t)
// (t minuscule : t), densité e^{−y − t} (fois α/(x − m) pour Fréchet), un arrondi.
// Quantiles : y = −log(−log p) (cquantile : −log(−log1p(−q))), puis μ + β y ou
// m + s e^{y/α}. Moments de Fréchet par Γ(1 − 1/α), variance sans compensation
// (detail/gamma_moments.hpp).
//
// Tirage (approché, bit-exact) : μ − β log E, m + s E^{−1/α}, E exponentielle exacte.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/gamma_moments.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/exponential.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    namespace detail
    {
        // P(Y <= y) = e^{−e^{−y}} (ou la survie), y en double-double.
        constexpr double gumbel_distribution(sd::double_double y, bool upper) noexcept
        {
            if (!(y.hi == y.hi))
                return y.hi;
            if (y.hi < -7.32)                                               // t > 1500 : e^{−t} < 2^-2000
                return upper ? 1.0 : 0.0;
            if (y.hi > 1500)
                return upper ? 0.0 : 1.0;
            scaled const ts = exp_scaled(-y.hi, -y.lo);                     // t = e^{−y}
            if (ts.e < -900)                                                // 1 − e^{−t} = t (1 − t/2 …)
                return upper ? rounded(ts) : 1.0;
            auto const t = unscaled(ts);
            if (upper)
                return t.hi > 40 ? 1.0 : -expm1_dd(-t.hi, -t.lo);
            return rounded(exp_scaled(-t.hi, -t.lo));
        }

        // −y − e^{−y}, en double-double (−∞ si e^{−y} déborde).
        constexpr sd::double_double gumbel_log_kernel(sd::double_double y) noexcept
        {
            if (y.hi < -709)
                return {-std::numeric_limits<double>::infinity(), 0.0};
            sd::double_double t{0.0, 0.0};
            if (y.hi < 1500)
                t = unscaled(exp_scaled(-y.hi, -y.lo));
            return sum({-y.hi, -y.lo}, {-t.hi, -t.lo});
        }

        // −log(−log p) (p dans (0, 1)) ; upper : −log(−log(1 − q)).
        constexpr sd::double_double gumbel_quantile_kernel(double p, bool upper) noexcept
        {
            auto const l = upper ? log1p_dd(-p) : log_precise(p);           // < 0
            auto const ll = log(normalized(-l.hi, -l.lo));
            return {-ll.hi, -ll.lo};
        }

    }

    template <law_real T = double>
    class gumbel : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T location = 0, scale = 1;
        };

        constexpr gumbel() : gumbel(parameters{}) {}

        constexpr gumbel(parameters p) : mu_(p.location), beta_(p.scale)
        {
            if (!detail::is_finite(mu_) || !(beta_ > 0) || !detail::is_finite(beta_))
                throw std::domain_error("wsl::random::gumbel : il faut location fini et scale > 0, fini");
        }

        constexpr T location() const noexcept { return mu_; }
        constexpr T scale() const noexcept { return beta_; }
        constexpr interval<T> support() const noexcept { return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()}; }
        constexpr T mean() const noexcept { return static_cast<T>(detail::affine(mu_, beta_, detail::sd::euler_gamma_hi, detail::sd::euler_gamma_lo)); }   // μ + βγ
        constexpr T mode() const noexcept { return mu_; }
        // μ − β log log 2
        constexpr T median() const noexcept
        {
            auto const l = detail::log(detail::normalized(detail::sd::ln2_hi, detail::sd::ln2_lo));
            return static_cast<T>(detail::affine(mu_, beta_, -l.hi, -l.lo));
        }
        // β² π²/6
        constexpr T variance() const noexcept
        {
            detail::scaled const b = detail::normalized(double(beta_), 0.0);
            return static_cast<T>(detail::rounded(b * b * detail::normalized(detail::pi_squared_over_3_hi, detail::pi_squared_over_3_lo) / 2.0));
        }
        constexpr T skewness() const noexcept { return static_cast<T>(detail::gumbel_skewness_hi); }
        constexpr T excess_kurtosis() const noexcept { return static_cast<T>(2.4); }
        // log β + γ + 1
        constexpr T entropy() const noexcept
        {
            auto h = detail::sum(detail::log_precise(double(beta_)), {detail::sd::euler_gamma_hi, detail::sd::euler_gamma_lo});
            h = detail::sum(h, {1.0, 0.0});
            return static_cast<T>(h.hi + h.lo);
        }

        constexpr T pdf(T x) const noexcept
        {
            auto const& G = *this;
            if (!(x == x))
                return x;
            auto const l = detail::gumbel_log_kernel(detail::standardize(x, G.mu_, G.beta_));
            if (!(l.hi > -1500))
                return 0;
            return static_cast<T>(detail::rounded(detail::exp_scaled(l.hi, l.lo) / double(G.beta_)));
        }
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& G = *this;
            if (!(x == x))
                return x;
            auto const l = detail::gumbel_log_kernel(detail::standardize(x, G.mu_, G.beta_));
            auto const lb = detail::log_precise(double(G.beta_));
            auto const r = detail::sum(l, {-lb.hi, -lb.lo});
            return static_cast<T>(r.hi + r.lo);
        }
        constexpr T cdf(T x) const noexcept { auto const& G = *this; return x == x ? static_cast<T>(detail::gumbel_distribution(detail::standardize(x, G.mu_, G.beta_), false)) : x; }
        constexpr T ccdf(T x) const noexcept { auto const& G = *this; return x == x ? static_cast<T>(detail::gumbel_distribution(detail::standardize(x, G.mu_, G.beta_), true)) : x; }

        constexpr T quantile(double p) const { auto const& G = *this; return G.inverse(p, false); }
        constexpr T cquantile(double q) const { auto const& G = *this; return G.inverse(q, true); }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            double const e = wsl::random::sample(exponential<double>({.rate = 1.0}), g);
            auto const l = detail::log_precise(e);
            return static_cast<T>(detail::affine(X.mu_, X.beta_, -l.hi, -l.lo));
        }

        friend constexpr bool operator==(gumbel const&, gumbel const&) = default;

    private:
        constexpr T inverse(double t, bool upper) const
        {
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::random::quantile : probabilité hors de [0, 1]");
            if (t == 0)
                return upper ? std::numeric_limits<T>::infinity() : -std::numeric_limits<T>::infinity();
            if (t == 1)
                return upper ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity();
            auto const y = detail::gumbel_quantile_kernel(t, upper);
            return static_cast<T>(detail::affine(mu_, beta_, y.hi, y.lo));
        }

        T mu_, beta_;
    };

    template <law_real T = double>
    class frechet : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T location = 0, scale = 1, shape;
        };

        constexpr frechet(parameters p) : alpha_(p.shape), s_(p.scale), m_(p.location)
        {
            if (!(alpha_ > 0) || !detail::is_finite(alpha_) || !(s_ > 0) || !detail::is_finite(s_) || !detail::is_finite(m_))
                throw std::domain_error("wsl::random::frechet : il faut shape > 0 et scale > 0, finis, location fini");
        }

        constexpr T shape() const noexcept { return alpha_; }
        constexpr T scale() const noexcept { return s_; }
        constexpr T location() const noexcept { return m_; }
        constexpr interval<T> support() const noexcept { return {m_, std::numeric_limits<T>::infinity()}; }

        // m + s Γ(1 − 1/α), α > 1
        constexpr T mean() const noexcept
        {
            if (!(double(alpha_) > 1))
                return std::numeric_limits<T>::infinity();
            auto const lg = detail::lgamma_one_plus(-1.0 / double(alpha_), -inv_lo());
            auto const l = detail::sum(lg, detail::log_precise(double(s_)));
            return static_cast<T>(shifted(l));
        }
        // s² [Γ(1 − 2/α) − Γ(1 − 1/α)²], α > 2
        constexpr T variance() const noexcept
        {
            if (!(double(alpha_) > 2))
                return std::numeric_limits<T>::infinity();
            return static_cast<T>(detail::gamma_variance(-1.0 / double(alpha_), -inv_lo(), detail::log_precise(double(s_))));
        }
        // m + s (log 2)^{−1/α}
        constexpr T median() const noexcept
        {
            auto const l2 = detail::log(detail::normalized(detail::sd::ln2_hi, detail::sd::ln2_lo));
            auto const l = detail::sum(detail::sd::quotient_dd({-l2.hi, -l2.lo}, {double(alpha_), 0.0}), detail::log_precise(double(s_)));
            return static_cast<T>(shifted(l));
        }
        // m + s (α/(1 + α))^{1/α}
        constexpr T mode() const noexcept
        {
            auto const r = detail::sd::quotient_dd({1.0, 0.0}, {double(alpha_), 0.0});
            auto const lp = detail::log1p_positive(r.hi, r.lo);              // log(1 + 1/α)
            auto const l = detail::sum(detail::sd::quotient_dd({-lp.hi, -lp.lo}, {double(alpha_), 0.0}), detail::log_precise(double(s_)));
            return static_cast<T>(shifted(l));
        }
        // 1 + γ/α + γ + log(s/α)
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            auto h = detail::sum({1.0, 0.0}, {fp::euler_gamma_hi, fp::euler_gamma_lo});
            h = detail::sum(h, fp::quotient_dd({fp::euler_gamma_hi, fp::euler_gamma_lo}, {double(alpha_), 0.0}));
            h = detail::sum(h, detail::log_precise(double(s_)));
            auto const la = detail::log_precise(double(alpha_));
            h = detail::sum(h, {-la.hi, -la.lo});
            return static_cast<T>(h.hi + h.lo);
        }

        // (α/(x − m)) e^{−y − t}, y = α log((x − m)/s)
        constexpr T pdf(T x) const noexcept
        {
            auto const& F = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > double(F.m_)) || v == std::numeric_limits<double>::infinity())
                return 0;
            auto const l = F.log_density_dd(v);
            return static_cast<T>(detail::rounded_exp(l));
        }
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& F = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > double(F.m_)) || v == std::numeric_limits<double>::infinity())
                return -std::numeric_limits<T>::infinity();
            auto const l = F.log_density_dd(v);
            return static_cast<T>(l.hi + l.lo);
        }
        constexpr T cdf(T x) const noexcept { auto const& F = *this; return F.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& F = *this; return F.distribution(x, true); }

        // m + s e^{y/α}, y = −log(−log p)
        constexpr T quantile(double p) const { auto const& F = *this; return F.inverse(p, false); }
        constexpr T cquantile(double q) const { auto const& F = *this; return F.inverse(q, true); }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& F = *this;
            double const e = wsl::random::sample(exponential<double>({.rate = 1.0}), g);
            double const w = wsl::special::exp(-wsl::special::log(e) / double(F.alpha_));
            return static_cast<T>(double(F.m_) + detail::sd::fmul(double(F.s_), w));
        }

        friend constexpr bool operator==(frechet const&, frechet const&) = default;

    private:
        // 1/α − fl(1/α)
        constexpr double inv_lo() const noexcept
        {
            double const a = alpha_, r = 1.0 / a;
            auto const [m, m_error] = detail::sd::two_product(r, a);
            return ((1.0 - m) - m_error) / a;
        }

        // m + e^{l}, un arrondi.
        constexpr double shifted(detail::sd::double_double l) const noexcept
        {
            if (l.hi > 1400)
                return std::numeric_limits<double>::infinity();
            if (l.hi < -1500)
                return m_;
            return detail::affine_scaled(double(m_), false, detail::exp_scaled(l.hi, l.lo));
        }

        // y = α log((x − m)/s), le log du rapport par log1p près de 1.
        constexpr detail::sd::double_double reduced(double x) const noexcept
        {
            namespace fp = detail::sd;
            auto const [d_hi, d_lo] = fp::two_sum(x, -double(m_));
            detail::scaled const r = detail::normalized(d_hi, d_lo) / double(s_);
            fp::double_double l{};
            if (r.e == 0 || (r.e == -1 && r.hi > 1.0))
            {
                auto const [r_hi, r_lo] = detail::unscaled(r);
                auto const [g_hi, g_lo] = fp::fast_two_sum(r_hi - 1.0, r_lo);
                if (g_hi >= 0)
                    l = detail::log1p_positive(g_hi, g_lo);
                else
                {
                    auto const [h, lo] = detail::log1p_dd(g_hi);
                    l = {h, lo + g_lo / (1.0 + g_hi)};
                }
            }
            else
                l = detail::log(r);
            return detail::mul({double(alpha_), 0.0}, l);
        }

        // log α − log(x − m) − y − e^{−y}
        constexpr detail::sd::double_double log_density_dd(double x) const noexcept
        {
            namespace fp = detail::sd;
            auto const y = reduced(x);
            auto l = detail::gumbel_log_kernel(y);
            if (!(l.hi > -std::numeric_limits<double>::infinity()))
                return l;
            auto const [d_hi, d_lo] = fp::two_sum(x, -double(m_));
            auto const ld = detail::log(detail::normalized(d_hi, d_lo));
            l = detail::sum(l, detail::log_precise(double(alpha_)));
            return detail::sum(l, {-ld.hi, -ld.lo});
        }

        constexpr T distribution(T x, bool upper) const noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > double(m_)))
                return upper ? T(1) : T(0);
            if (v == std::numeric_limits<double>::infinity())
                return upper ? T(0) : T(1);
            return static_cast<T>(detail::gumbel_distribution(reduced(v), upper));
        }

        constexpr T inverse(double t, bool upper) const
        {
            namespace fp = detail::sd;
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::random::quantile : probabilité hors de [0, 1]");
            if (t == 0)
                return upper ? std::numeric_limits<T>::infinity() : m_;
            if (t == 1)
                return upper ? m_ : std::numeric_limits<T>::infinity();
            auto const y = detail::gumbel_quantile_kernel(t, upper);
            auto const l = detail::sum(fp::quotient_dd(y, {double(alpha_), 0.0}), detail::log_precise(double(s_)));
            return static_cast<T>(shifted(l));
        }

        T alpha_, s_, m_;
    };
}
