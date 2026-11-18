#pragma once

// Loi exponentielle de paramètre λ > 0 (rate) : densité λ e^{−λx} sur [0, +∞).
//
//   wsl::exponential X({.rate = 2});
//   X.mean(); cdf(X, 1.0); quantile(X, 0.5); sample(X, g); sample(X, g, wsl::fast);
//
// Tirage exact : E exponentielle de paramètre 1 tirée exactement (von Neumann,
// réel paresseux), puis l'arrondi correct de E/λ — la valeur tirée suit donc la
// loi exacte arrondie au flottant le plus proche. Tirage rapide (wsl::fast) :
// inversion, −log(U)/λ, deux arrondis.
//
// Densité, répartition, survie, quantile : λx en double-double (produit exact),
// exp, expm1 et log1p à argument double-double, un seul arrondi final ; erreur
// mesurée sous 0,51 ulp (tests/test_random_precision.cpp).

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl
{
    template <law_real T = double>
    class exponential
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            T rate;
        };

        constexpr exponential(parameters p) : rate_(p.rate)
        {
            if (!(rate_ > 0) || !detail::is_finite(rate_))
                throw std::domain_error("wsl::exponential : il faut rate > 0, fini");
        }

        constexpr T rate() const noexcept { return rate_; }
        constexpr interval<T> support() const noexcept { return {T(0), std::numeric_limits<T>::infinity()}; }

        constexpr T mean() const noexcept { return static_cast<T>(1.0 / double(rate_)); }
        constexpr T variance() const noexcept
        {
            return static_cast<T>(detail::rounded(detail::normalized(1.0, 0.0) / double(rate_) / double(rate_)));
        }
        constexpr T median() const noexcept
        {
            return static_cast<T>(detail::rounded(detail::normalized(detail::sd::ln2_hi, detail::sd::ln2_lo) / double(rate_)));
        }
        constexpr T mode() const noexcept { return 0; }
        constexpr T skewness() const noexcept { return 2; }
        constexpr T excess_kurtosis() const noexcept { return 6; }
        // 1 − log λ.
        constexpr T entropy() const noexcept
        {
            auto const [hi, lo] = detail::log_precise(double(rate_));
            auto const [s, s_error] = detail::sd::two_sum(1.0, -hi);
            return static_cast<T>(s + (s_error - lo));
        }

        friend constexpr T pdf(exponential const& X, T x) noexcept
        {
            double const l = X.rate_, v = x;
            if (!(v >= 0))
                return v == v ? T(0) : x;
            double const y = l * v;                                         // λx arrondi : pour choisir le régime
            if (y > 1460)                                                   // λ e^{−λx} < 2^-1075, +∞ compris
                return 0;
            if (y < 0x1p-900)                                               // λ (1 − λx) s'arrondit à λ
                return X.rate_;
            auto const [p, p_error] = detail::balanced_product(l, v);
            return static_cast<T>(detail::rounded(detail::exp_scaled(-p, -p_error) * detail::normalized(l, 0.0)));
        }

        // log λ − λx.
        friend constexpr T log_pdf(exponential const& X, T x) noexcept
        {
            double const l = X.rate_, v = x;
            if (!(v >= 0))
                return v == v ? -std::numeric_limits<T>::infinity() : x;
            double const y = l * v;
            if (y > 0x1p100)                                                // log λ sous le demi-ulp
                return static_cast<T>(-y);
            auto const [hi, lo] = detail::log_precise(l);
            if (y < 0x1p-900)
                return static_cast<T>(hi + lo);
            auto const [p, p_error] = detail::balanced_product(l, v);
            auto const [s, s_error] = detail::sd::two_sum(hi, -p);
            return static_cast<T>(s + ((s_error + lo) - p_error));
        }

        // 1 − e^{−λx} = −expm1(−λx).
        friend constexpr T cdf(exponential const& X, T x) noexcept
        {
            double const l = X.rate_, v = x;
            if (!(v > 0))
                return v == v ? T(0) : x;
            double const y = l * v;
            if (y > 40)                                                     // e^{−40} < 2^-57
                return 1;
            if (y < 0x1p-900)                                               // λx − (λx)²/2 s'arrondit à λx
                return static_cast<T>(y);
            auto const [p, p_error] = detail::balanced_product(l, v);
            return static_cast<T>(-detail::expm1_dd(-p, -p_error));
        }

        // e^{−λx}.
        friend constexpr T ccdf(exponential const& X, T x) noexcept
        {
            double const l = X.rate_, v = x;
            if (!(v > 0))
                return v == v ? T(1) : x;
            double const y = l * v;
            if (y > 746)                                                    // < 2^-1076
                return 0;
            if (y < 0x1p-900)
                return 1;
            auto const [p, p_error] = detail::balanced_product(l, v);
            return static_cast<T>(detail::rounded(detail::exp_scaled(-p, -p_error)));
        }

        // −log(1 − p)/λ.
        friend constexpr T quantile(exponential const& X, double p)
        {
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return 0;
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            auto const [hi, lo] = detail::log1p_dd(-p);
            return static_cast<T>(detail::rounded(detail::normalized(-hi, -lo) / double(X.rate_)));
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        friend constexpr T cquantile(exponential const& X, double q)
        {
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return quantile(X, 1.0 - q);                                // 1 − q exact
            if (q == 0)
                return std::numeric_limits<T>::infinity();
            auto const [hi, lo] = detail::log_precise(q);                   // x = −log q / λ
            return static_cast<T>(detail::rounded(detail::normalized(-hi, -lo) / double(X.rate_)));
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(exponential const& X, G& g)
        {
            detail::lazy_real v = detail::exact_exponential(g);
            return detail::round_lazy<T>(0.0, 1.0, double(X.rate_), false, v, g);
        }

        // Inversion : −log(U)/λ, U uniforme dans (0, 1) (approché, bit-exact).
        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(exponential const& X, G& g, fast_t)
        {
            return static_cast<T>(-wsl::special::log(detail::uniform_positive(g)) / double(X.rate_));
        }

        friend constexpr bool operator==(exponential const&, exponential const&) = default;

    private:
        T rate_;
    };
}
