#pragma once

// Loi exponentielle de paramètre λ > 0 (rate) : densité λ e^{−λx} sur [0, +∞).
//
//   wsl::random::exponential X({.rate = 2});
//   X.mean(); X.cdf(1.0); X.quantile(0.5); X(g); X(g, wsl::random::fast);
//
// Tirage exact : E exponentielle de paramètre 1 tirée exactement (von Neumann,
// réel paresseux), puis l'arrondi correct de E/λ — la valeur tirée suit donc la
// loi exacte arrondie au flottant le plus proche. Tirage rapide (wsl::random::fast) :
// inversion, −log(U)/λ, deux arrondis.
//
// Densité, répartition, survie, quantile : λx en double-double (produit exact),
// exp, expm1 et log1p à argument double-double, un seul arrondi final ; erreur
// mesurée sous 0,51 ulp (tests/test_random_precision.cpp).

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/rational.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    template <law_real T = double>
    class exponential : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::exact;

        // Intensité .rate ou échelle .scale (l'une des deux au plus ; aucune : λ = 1), double ou
        // wr::ratio(a, b) : exactes toutes deux (λ rationnel, voir detail/rational.hpp).
        struct parameters
        {
            std::optional<ratio> rate = std::nullopt, scale = std::nullopt;
        };
        using param_type = parameters;
        constexpr param_type param() const
        {
            auto const [r, s] = detail::parameters_of(lambda_);
            return {.rate = r, .scale = s};
        }

        constexpr exponential() : exponential(parameters{}) {}
        constexpr exponential(parameters p)
            : lambda_(detail::rate_of(p.rate, p.scale, "wsl::random::exponential : il faut .rate ou .scale (pas les deux), > 0 et finis")),
              rate_(lambda_.as_double())
        {
        }

        // λ et 1/λ arrondis (exacts quand ils sont des doubles).
        constexpr T rate() const noexcept { return static_cast<T>(rate_ != 0 ? rate_ : detail::rounded(lambda_.value())); }
        constexpr T scale() const noexcept { return static_cast<T>(detail::rounded(per_lambda(detail::normalized(1.0, 0.0)))); }
        constexpr detail::exact_rate const& exact_rate() const noexcept { return lambda_; }
        constexpr interval<T> support() const noexcept { return {T(0), std::numeric_limits<T>::infinity()}; }

        constexpr T mean() const noexcept
        {
            return static_cast<T>(rate_ != 0 ? 1.0 / rate_ : detail::rounded(per_lambda(detail::normalized(1.0, 0.0))));
        }
        constexpr T variance() const noexcept
        {
            return static_cast<T>(detail::rounded(per_lambda(per_lambda(detail::normalized(1.0, 0.0)))));
        }
        constexpr T median() const noexcept
        {
            return static_cast<T>(detail::rounded(per_lambda(detail::normalized(detail::sd::ln2_hi, detail::sd::ln2_lo))));
        }
        constexpr T mode() const noexcept { return 0; }
        constexpr T skewness() const noexcept { return 2; }
        constexpr T excess_kurtosis() const noexcept { return 6; }
        // 1 − log λ.
        constexpr T entropy() const noexcept
        {
            auto const [hi, lo] = log_lambda();
            auto const [s, s_error] = detail::sd::two_sum(1.0, -hi);
            return static_cast<T>(s + (s_error - lo));
        }

        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v >= 0))
                return v == v ? T(0) : x;
            double const y = X.approximate_product(v);                      // λx arrondi : pour choisir le régime
            if (y > 1460)                                                   // λ e^{−λx} < 2^-1075, +∞ compris
                return 0;
            if (y < 0x1p-900)                                               // λ (1 − λx) s'arrondit à λ
                return X.rate();
            auto const [p, p_error] = X.product(v);
            return static_cast<T>(detail::rounded(detail::exp_scaled(-p, -p_error) * X.lambda_scaled()));
        }

        // log λ − λx.
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v >= 0))
                return v == v ? -std::numeric_limits<T>::infinity() : x;
            double const y = X.approximate_product(v);
            if (y > 0x1p100)                                                // log λ sous le demi-ulp
                return static_cast<T>(-y);
            auto const [hi, lo] = X.log_lambda();
            if (y < 0x1p-900)
                return static_cast<T>(hi + lo);
            auto const [p, p_error] = X.product(v);
            auto const [s, s_error] = detail::sd::two_sum(hi, -p);
            return static_cast<T>(s + ((s_error + lo) - p_error));
        }

        // 1 − e^{−λx} = −expm1(−λx).
        constexpr T cdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v > 0))
                return v == v ? T(0) : x;
            double const y = X.approximate_product(v);
            if (y > 40)                                                     // e^{−40} < 2^-57
                return 1;
            if (y < 0x1p-900)                                               // λx − (λx)²/2 s'arrondit à λx
                return static_cast<T>(y);
            auto const [p, p_error] = X.product(v);
            return static_cast<T>(-detail::expm1_dd(-p, -p_error));
        }

        // e^{−λx}.
        constexpr T ccdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v > 0))
                return v == v ? T(1) : x;
            double const y = X.approximate_product(v);
            if (y > 746)                                                    // < 2^-1076
                return 0;
            if (y < 0x1p-900)
                return 1;
            auto const [p, p_error] = X.product(v);
            return static_cast<T>(detail::rounded(detail::exp_scaled(-p, -p_error)));
        }

        // −log(1 − p)/λ.
        constexpr T quantile(double p) const
        {
            auto const& X = *this;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::random::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return 0;
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            auto const [hi, lo] = detail::log1p_dd(-p);
            return static_cast<T>(detail::rounded(X.per_lambda(detail::normalized(-hi, -lo))));
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        constexpr T cquantile(double q) const
        {
            auto const& X = *this;
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::random::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return wsl::random::quantile(X, 1.0 - q);                                // 1 − q exact
            if (q == 0)
                return std::numeric_limits<T>::infinity();
            auto const [hi, lo] = detail::log_precise(q);                   // x = −log q / λ
            return static_cast<T>(detail::rounded(X.per_lambda(detail::normalized(-hi, -lo))));
        }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            detail::lazy_real v = detail::exact_exponential(g);
            if (X.rate_ != 0)
                return detail::round_lazy<T>(0.0, 1.0, X.rate_, false, v, g);
            auto const f = X.lambda_.sampling_factors();                     // V/λ = V σ / d, exact
            return detail::round_lazy<T>(0.0, f.sigma, f.divisor, false, v, g);
        }

        // Inversion : −log(U)/λ, U uniforme dans (0, 1) (approché, bit-exact).
        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g, fast_t) const
        {
            auto const& X = *this;
            double const e = -wsl::special::log(detail::uniform_positive(g));
            return static_cast<T>(X.rate_ != 0 ? e / X.rate_ : detail::rounded(X.per_lambda(detail::normalized(e, 0.0))));
        }

        friend constexpr bool operator==(exponential const& a, exponential const& b) noexcept { return a.lambda_ == b.lambda_; }

    private:
        // λx en double-double (x > 0 fini), et arrondi (x > 0, +∞ compris) ; voie rapide si λ est un double.
        constexpr detail::sd::double_double product(double x) const noexcept
        {
            if (rate_ != 0)
            {
                auto const [p, p_error] = detail::balanced_product(rate_, x);
                return {p, p_error};
            }
            return detail::unscaled(lambda_.times(detail::normalized(x, 0.0)));
        }
        constexpr double approximate_product(double x) const noexcept
        {
            if (rate_ != 0)
                return rate_ * x;
            return detail::is_finite(x) ? detail::rounded(lambda_.times(detail::normalized(x, 0.0))) : x;
        }
        constexpr detail::scaled lambda_scaled() const noexcept { return rate_ != 0 ? detail::normalized(rate_, 0.0) : lambda_.value(); }
        constexpr detail::scaled per_lambda(detail::scaled const& v) const noexcept { return rate_ != 0 ? v / rate_ : lambda_.per(v); }
        constexpr detail::sd::double_double log_lambda() const noexcept { return rate_ != 0 ? detail::log_precise(rate_) : lambda_.log(); }

        detail::exact_rate lambda_;
        double rate_ = 0;                                                   // λ s'il est un double (sinon 0)
    };
}
