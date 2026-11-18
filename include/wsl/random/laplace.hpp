#pragma once

// Loi de Laplace (double exponentielle) de position μ et d'échelle b > 0 :
// densité e^{−|x − μ|/b}/(2b).
//
//   wsl::laplace L({.location = 0, .scale = 1});
//   cdf(L, 1.5); quantile(L, 0.9); sample(L, g);
//
// z = (x − μ)/b en double-double ; répartition e^{z}/2 (z < 0), 1 − e^{−z}/2 sinon ;
// quantile μ + b log(2p) (p < ½), μ − b log(2(1 − p)) ; un seul arrondi.
//
// Tirage exact : μ ± b E, E exponentielle tirée exactement, signe uniforme, puis
// l'arrondi correct (comme pour la normale).

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
    template <law_real T = double>
    class laplace
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            T location, scale;
        };

        constexpr laplace(parameters p) : mu_(p.location), b_(p.scale)
        {
            if (!detail::is_finite(mu_) || !(b_ > 0) || !detail::is_finite(b_))
                throw std::domain_error("wsl::laplace : il faut location fini et scale > 0, fini");
        }

        constexpr T location() const noexcept { return mu_; }
        constexpr T scale() const noexcept { return b_; }
        constexpr interval<T> support() const noexcept
        {
            return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()};
        }
        constexpr T mean() const noexcept { return mu_; }
        constexpr T median() const noexcept { return mu_; }
        constexpr T mode() const noexcept { return mu_; }
        constexpr T variance() const noexcept { return static_cast<T>(detail::sd::fmul(2.0, detail::sd::fmul(double(b_), double(b_)))); }
        constexpr T skewness() const noexcept { return 0; }
        constexpr T excess_kurtosis() const noexcept { return 3; }
        // 1 + log 2 + log b
        constexpr T entropy() const noexcept
        {
            auto e = detail::sum(detail::log_precise(double(b_)), {detail::sd::ln2_hi, detail::sd::ln2_lo});
            e = detail::sum(e, {1.0, 0.0});
            return static_cast<T>(e.hi + e.lo);
        }

        // e^{−|z|}/(2b)
        friend constexpr T pdf(laplace const& X, T x) noexcept
        {
            if (!(x == x))
                return x;
            auto const [z_hi, z_lo] = detail::standardize(x, X.mu_, X.b_);
            double const a = detail::magnitude(z_hi), a_lo = z_hi < 0 ? -z_lo : z_lo;
            if (!(a < 1500))
                return 0;
            detail::scaled const e = detail::exp_scaled(-a, -a_lo);
            return static_cast<T>(detail::rounded(e / double(X.b_) / 2.0));
        }
        // −|z| − log(2b)
        friend constexpr T log_pdf(laplace const& X, T x) noexcept
        {
            if (!(x == x))
                return x;
            auto const [z_hi, z_lo] = detail::standardize(x, X.mu_, X.b_);
            if (!detail::is_finite(z_hi) && detail::is_finite(x))           // |z| > 2^600 : −|z| (le reste, négligeable)
                return static_cast<T>(-detail::rounded(detail::standardize_scaled(x, X.mu_, X.b_).z));
            double const a = detail::magnitude(z_hi), a_lo = z_hi < 0 ? -z_lo : z_lo;
            auto e = detail::sum(detail::log_precise(double(X.b_)), {detail::sd::ln2_hi, detail::sd::ln2_lo});
            e = detail::sum(e, {a, a_lo});
            return static_cast<T>(-(e.hi + e.lo));
        }

        friend constexpr T cdf(laplace const& X, T x) noexcept { return X.distribution(x, false); }
        friend constexpr T ccdf(laplace const& X, T x) noexcept { return X.distribution(x, true); }

        // μ + b log(2p) (p < ½), μ − b log(2(1 − p))
        friend constexpr T quantile(laplace const& X, double p)
        {
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return -std::numeric_limits<T>::infinity();
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            if (p == 0.5)
                return X.mu_;
            bool const lower = p < 0.5;
            double const q = detail::sd::fmul(lower ? p : 1.0 - p, 2.0);   // exact, dans (0, 1)
            auto const [l_hi, l_lo] = detail::log_precise(q);               // < 0
            detail::scaled const t = detail::normalized(-l_hi, -l_lo) * detail::normalized(double(X.b_), 0.0);
            return static_cast<T>(detail::affine_scaled(X.mu_, lower, t));
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        friend constexpr T cquantile(laplace const& X, double q)
        {
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return quantile(X, 1.0 - q);                                // 1 − q exact
            if (q == 0)
                return std::numeric_limits<T>::infinity();
            auto const [l_hi, l_lo] = detail::log_precise(detail::sd::fmul(q, 2.0));   // μ − b log(2q)
            detail::scaled const t = detail::normalized(-l_hi, -l_lo) * detail::normalized(double(X.b_), 0.0);
            return static_cast<T>(detail::affine_scaled(X.mu_, false, t));
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(laplace const& X, G& g)
        {
            detail::lazy_real v = detail::exact_exponential(g);
            bool const negative = (random_bits(g) >> 63) != 0;
            return detail::round_lazy<T>(X.mu_, X.b_, 1.0, negative, v, g);
        }

        friend constexpr bool operator==(laplace const&, laplace const&) = default;

    private:
        // queue P(X > μ + b|z|) = e^{−|z|}/2
        constexpr T distribution(T x, bool upper) const noexcept
        {
            namespace fp = detail::sd;
            if (!(x == x))
                return x;
            auto const [z_hi, z_lo] = detail::standardize(x, mu_, b_);
            if (z_hi == 0)
                return T(0.5);
            bool const tail_side = (z_hi < 0) != upper;
            double const a = detail::magnitude(z_hi), a_lo = z_hi < 0 ? -z_lo : z_lo;
            if (!(a < 1500))
                return tail_side ? T(0) : T(1);
            detail::scaled tail = detail::exp_scaled(-a, -a_lo);
            tail.e -= 1;
            return static_cast<T>(tail_side ? detail::rounded(tail) : fp::one_minus(tail));
        }

        T mu_, b_;
    };
}
