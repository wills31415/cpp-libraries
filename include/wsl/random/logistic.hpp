#pragma once

// Loi logistique de position μ et d'échelle s > 0 : P(X <= x) = 1/(1 + e^{−z}),
// z = (x − μ)/s.
//
//   wsl::logistic L({.location = 0, .scale = 1});
//   cdf(L, 2.0); quantile(L, 0.9); sample(L, g);
//
// z en double-double ; queue e^{−|z|}/(1 + e^{−|z|}), densité e^{−|z|}/(s (1 + e^{−|z|})²) ;
// quantile μ + s logit(p), logit(p) = log1p((2p − 1)/(1 − p)) près de ½ (sans
// compensation), log p − log1p(−p) ailleurs ; un seul arrondi.
//
// Tirage (approché, bit-exact) : μ + s logit(U).

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl
{
    namespace detail
    {
        // logit(p) = log(p/(1 − p)), p dans (0, 1), en double-double.
        constexpr sd::double_double logit(double p) noexcept
        {
            if (p >= 0.25 && p <= 0.75)                                     // log1p((2p − 1)/(1 − p))
            {
                auto const [q_hi, q_lo] = sd::two_sum(1.0, -p);             // exact
                sd::double_double const r = sd::quotient_dd({sd::fmul(p, 2.0) - 1.0, 0.0}, {q_hi, q_lo});   // 2p − 1 exact
                if (r.hi >= 0)
                    return log1p_positive(r.hi, r.lo);
                auto const [h, l] = log1p_dd(r.hi);
                return {h, l + r.lo / (1.0 + r.hi)};
            }
            auto const lp = log_precise(p);
            auto const lq = p < 0.5 ? log1p_dd(-p) : log(normalized(1.0 - p, 0.0));   // 1 − p exact pour p > ½
            return sum(lp, {-lq.hi, -lq.lo});
        }
    }

    template <law_real T = double>
    class logistic
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T location, scale;
        };

        constexpr logistic(parameters p) : mu_(p.location), s_(p.scale)
        {
            if (!detail::is_finite(mu_) || !(s_ > 0) || !detail::is_finite(s_))
                throw std::domain_error("wsl::logistic : il faut location fini et scale > 0, fini");
        }

        constexpr T location() const noexcept { return mu_; }
        constexpr T scale() const noexcept { return s_; }
        constexpr interval<T> support() const noexcept
        {
            return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()};
        }
        constexpr T mean() const noexcept { return mu_; }
        constexpr T median() const noexcept { return mu_; }
        constexpr T mode() const noexcept { return mu_; }
        // s² π²/3
        constexpr T variance() const noexcept
        {
            detail::scaled const s = detail::normalized(double(s_), 0.0);
            return static_cast<T>(detail::rounded(s * s * detail::normalized(detail::pi_squared_over_3_hi, detail::pi_squared_over_3_lo)));
        }
        constexpr T skewness() const noexcept { return 0; }
        constexpr T excess_kurtosis() const noexcept { return static_cast<T>(1.2); }
        // log s + 2
        constexpr T entropy() const noexcept
        {
            auto const e = detail::sum(detail::log_precise(double(s_)), {2.0, 0.0});
            return static_cast<T>(e.hi + e.lo);
        }

        // e^{−|z|} / (s (1 + e^{−|z|})²)
        friend constexpr T pdf(logistic const& X, T x) noexcept
        {
            if (!(x == x))
                return x;
            auto const [e, one_plus] = X.parts(x);
            if (e.hi == 0)
                return 0;
            return static_cast<T>(detail::rounded(e / (one_plus * one_plus) / double(X.s_)));
        }
        friend constexpr T log_pdf(logistic const& X, T x) noexcept
        {
            if (!(x == x))
                return x;
            auto const [z_hi, z_lo] = detail::standardize(x, X.mu_, X.s_);
            double const a = detail::magnitude(z_hi), a_lo = z_hi < 0 ? -z_lo : z_lo;
            if (!(a < 1500))                                                 // −|z| − log s
            {
                auto const e = detail::sum({a, a_lo}, detail::log_precise(double(X.s_)));
                return static_cast<T>(-(e.hi + e.lo));
            }
            auto const lp = detail::log(X.parts(x).one_plus);              // −|z| − 2 log(1 + e^{−|z|}) − log s
            auto l = detail::sum({-a, -a_lo}, {detail::sd::fmul(lp.hi, -2.0), detail::sd::fmul(lp.lo, -2.0)});
            auto const ls = detail::log_precise(double(X.s_));
            l = detail::sum(l, {-ls.hi, -ls.lo});
            return static_cast<T>(l.hi + l.lo);
        }

        friend constexpr T cdf(logistic const& X, T x) noexcept { return X.distribution(x, false); }
        friend constexpr T ccdf(logistic const& X, T x) noexcept { return X.distribution(x, true); }

        friend constexpr T quantile(logistic const& X, double p)
        {
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return -std::numeric_limits<T>::infinity();
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            if (p == 0.5)
                return X.mu_;
            auto const [l_hi, l_lo] = detail::logit(p);
            return static_cast<T>(detail::affine(X.mu_, X.s_, l_hi, l_lo));
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(logistic const& X, G& g)
        {
            double const u = detail::uniform_positive(g);
            if (u == 0.5)
                return X.mu_;
            auto const [l_hi, l_lo] = detail::logit(u);
            return static_cast<T>(detail::affine(X.mu_, X.s_, l_hi, l_lo));
        }

        friend constexpr bool operator==(logistic const&, logistic const&) = default;

    private:
        struct exp_parts
        {
            detail::scaled e, one_plus;                                     // e^{−|z|}, 1 + e^{−|z|}
        };
        constexpr exp_parts parts(T x) const noexcept
        {
            auto const [z_hi, z_lo] = detail::standardize(x, mu_, s_);
            double const a = detail::magnitude(z_hi), a_lo = z_hi < 0 ? -z_lo : z_lo;
            if (!(a < 1500))
                return {{0.0, 0.0, 0}, detail::normalized(1.0, 0.0)};
            detail::scaled const e = detail::exp_scaled(-a, -a_lo);
            if (e.e < -200)
                return {e, detail::normalized(1.0, 0.0)};
            auto const [e_hi, e_lo] = detail::unscaled(e);
            auto const [o_hi, o_lo] = detail::sum({1.0, 0.0}, {e_hi, e_lo});
            return {e, detail::normalized(o_hi, o_lo)};
        }

        // queue e^{−|z|}/(1 + e^{−|z|})
        constexpr T distribution(T x, bool upper) const noexcept
        {
            namespace fp = detail::sd;
            if (!(x == x))
                return x;
            auto const [z_hi, z_lo] = detail::standardize(x, mu_, s_);
            if (z_hi == 0)
                return T(0.5);
            bool const tail_side = (z_hi < 0) != upper;
            auto const [e, one_plus] = parts(x);
            if (e.hi == 0)
                return tail_side ? T(0) : T(1);
            detail::scaled const tail = e / one_plus;
            return static_cast<T>(tail_side ? detail::rounded(tail) : fp::one_minus(tail));
        }

        T mu_, s_;
    };
}
