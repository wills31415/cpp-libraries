#pragma once

// Premières lois : uniforme continue, uniforme discrète, Bernoulli.
//
//   wsl::uniform X({.low = 2, .high = 6});              // [2, 6), double par défaut
//   wsl::uniform<float> Y({.low = 0, .high = 1});
//   wsl::discrete_uniform D({.low = 1, .high = 6});     // un dé : {1, …, 6}
//   wsl::bernoulli B({.p = 0.3});
//
//   X.mean(); X.variance(); cdf(X, 3.0); quantile(X, 0.25); sample(X, g);
//
// Les calculs se font en double (puis un arrondi vers float pour les lois en
// float), en double-double là où c'est nécessaire pour rester sous 1 ulp.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <concepts>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace wsl
{
    template <class T>
    concept law_real = std::same_as<T, float> || std::same_as<T, double>;

    template <class I>
    concept law_integer = std::integral<I> && !std::same_as<I, bool> && sizeof(I) <= 8;

    // -----------------------------------------------------------------------
    // Uniforme continue sur [low, high).
    //
    // Tirage : low + U (high − low), U uniforme exact sur [0, 1), arrondi une
    // fois (high exclu) : approximate (l'arrondi du produit n'est pas celui d'un
    // réel uniforme exact sur [low, high)).
    // -----------------------------------------------------------------------
    template <law_real T = double>
    class uniform
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct bounds
        {
            T low, high;
        };

        constexpr uniform(bounds b) : low_(b.low), high_(b.high)
        {
            if (!(low_ < high_) || !detail::is_finite(low_) || !detail::is_finite(high_))
                throw std::domain_error("wsl::uniform : il faut low < high, finis");
        }

        constexpr T low() const noexcept { return low_; }
        constexpr T high() const noexcept { return high_; }
        constexpr interval<T> support() const noexcept { return {low_, high_}; }

        constexpr T mean() const noexcept
        {
            double const a = low_, b = high_;
            return static_cast<T>(detail::is_finite(a + b) ? (a + b) * 0.5 : detail::sd::fmul(a, 0.5) + detail::sd::fmul(b, 0.5));
        }
        constexpr T median() const noexcept { return mean(); }

        // (b − a)² / 12, b − a en double-double, mis à l'échelle (débordement voulu
        // seulement si le résultat déborde).
        constexpr T variance() const noexcept
        {
            namespace sd = detail::sd;
            auto const [d, d_error, doubled] = width();
            int const e = sd::exponent_of(d);
            double const dn = sd::scale_by_power_of_two(d, -e), dl = sd::scale_by_power_of_two(d_error, -e);
            auto const [s, s_error] = sd::two_product(dn, dn);
            double const s_lo = s_error + sd::fmul(2.0, sd::fmul(dn, dl));
            return static_cast<T>(sd::scale_by_power_of_two(detail::divide(s, s_lo, 12.0, 0.0), 2 * (e + doubled)));
        }
        constexpr T skewness() const noexcept { return 0; }
        constexpr T excess_kurtosis() const noexcept { return static_cast<T>(-1.2); }
        constexpr T entropy() const noexcept
        {
            auto const [d, d_error, doubled] = width();
            return static_cast<T>(detail::log_dd(d, d_error) + (doubled ? detail::sd::ln2_hi : 0.0));
        }

        friend constexpr T pdf(uniform const& X, T x) noexcept
        {
            if (!(x >= X.low_ && x < X.high_))
                return 0;
            auto const [d, d_error, doubled] = X.width();
            return static_cast<T>(detail::sd::scale_by_power_of_two(detail::divide(1.0, 0.0, d, d_error), -doubled));
        }
        friend constexpr T log_pdf(uniform const& X, T x) noexcept
        {
            if (!(x >= X.low_ && x < X.high_))
                return -std::numeric_limits<T>::infinity();
            return -X.entropy();
        }
        friend constexpr T cdf(uniform const& X, T x) noexcept
        {
            if (!(x > X.low_))
                return x == x ? T(0) : x;
            if (!(x < X.high_))
                return 1;
            return static_cast<T>(detail::ratio(x, X.low_, X.high_));
        }
        // P(X > x) = (high − x)/(high − low), sans annulation.
        friend constexpr T ccdf(uniform const& X, T x) noexcept
        {
            if (!(x > X.low_))
                return x == x ? T(1) : x;
            if (!(x < X.high_))
                return 0;
            return static_cast<T>(detail::ratio(-static_cast<double>(x), -static_cast<double>(X.high_), -static_cast<double>(X.low_)));
        }
        friend constexpr T quantile(uniform const& X, double p)
        {
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::quantile : p doit être dans [0, 1]");
            if (p == 1)
                return X.high_;
            return X.clamp(detail::interpolate(X.low_, X.high_, p));
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(uniform const& X, G& g)
        {
            return X.clamp(detail::interpolate(X.low_, X.high_, uniform_real<double>(g)));
        }

        friend constexpr bool operator==(uniform const&, uniform const&) = default;

    private:
        // high − low = (hi + lo) · 2^doubled exactement ; sur les moitiés (doubled =
        // 1) quand la différence déborderait.
        struct difference
        {
            double hi, lo;
            int doubled;
        };
        constexpr difference width() const noexcept
        {
            double const a = low_, b = high_;
            if (detail::is_finite(b - a))
            {
                auto const [d, e] = detail::sd::two_sum(b, -a);
                return {d, e, 0};
            }
            auto const [d, e] = detail::sd::two_sum(detail::sd::fmul(b, 0.5), -detail::sd::fmul(a, 0.5));
            return {d, e, 1};
        }
        // Valeur arrondie dans [low, high) (l'arrondi peut atteindre high).
        constexpr T clamp(double r) const noexcept
        {
            auto const v = static_cast<T>(r);
            return v < high_ ? (v < low_ ? low_ : v) : detail::next_below(high_);
        }

        T low_, high_;
    };

    // -----------------------------------------------------------------------
    // Uniforme discrète sur {low, …, high} (bornes comprises). Tirage exact.
    // -----------------------------------------------------------------------
    template <law_integer I = std::int64_t>
    class discrete_uniform
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::exact;

        struct bounds
        {
            I low, high;
        };

        constexpr discrete_uniform(bounds b) : low_(b.low), high_(b.high)
        {
            if (high_ < low_)
                throw std::domain_error("wsl::discrete_uniform : il faut low <= high");
        }

        constexpr I low() const noexcept { return low_; }
        constexpr I high() const noexcept { return high_; }
        constexpr interval<I> support() const noexcept { return {low_, high_}; }

        constexpr double mean() const noexcept
        {
            double const a = static_cast<double>(low_), b = static_cast<double>(high_);
            return detail::is_finite(a + b) ? (a + b) * 0.5 : detail::sd::fmul(a, 0.5) + detail::sd::fmul(b, 0.5);
        }
        // Médiane inférieure (pour un nombre pair de valeurs, toute valeur entre
        // les deux du milieu en est une).
        constexpr I median() const noexcept { return static_cast<I>(low_ + static_cast<I>(span() / 2)); }
        // (n² − 1)/12, n = high − low + 1, en double-double.
        constexpr double variance() const noexcept
        {
            namespace sd = detail::sd;
            auto const [n, n_lo] = count();
            if (n > 0x1p500)
                return sd::fmul(n, n) / 12.0;                                // impossible ici : n <= 2^64
            auto const [s, s_error] = sd::two_product(n, n);
            auto const [t, t_error] = sd::two_sum(s, -1.0);
            return detail::divide(t, (t_error + s_error) + sd::fmul(2.0, sd::fmul(n, n_lo)), 12.0, 0.0);
        }
        constexpr double skewness() const noexcept { return 0; }
        constexpr double entropy() const noexcept
        {
            auto const [n, n_lo] = count();
            return detail::log_dd(n, n_lo);
        }

        friend constexpr double pmf(discrete_uniform const& X, I k) noexcept
        {
            if (k < X.low_ || k > X.high_)
                return 0;
            auto const [n, n_lo] = X.count();
            return detail::divide(1.0, 0.0, n, n_lo);
        }
        friend constexpr double log_pmf(discrete_uniform const& X, I k) noexcept
        {
            if (k < X.low_ || k > X.high_)
                return -std::numeric_limits<double>::infinity();
            return -X.entropy();
        }
        // P(X <= k) = (k − low + 1)/n.
        friend constexpr double cdf(discrete_uniform const& X, I k) noexcept
        {
            if (k < X.low_)
                return 0;
            if (k >= X.high_)
                return 1;
            auto const [m, m_lo] = split(static_cast<std::uint64_t>(X.offset(k)) + 1);
            auto const [n, n_lo] = X.count();
            return detail::divide(m, m_lo, n, n_lo);
        }
        // P(X > k) = (high − k)/n.
        friend constexpr double ccdf(discrete_uniform const& X, I k) noexcept
        {
            if (k < X.low_)
                return 1;
            if (k >= X.high_)
                return 0;
            auto const [m, m_lo] = split(X.span() - static_cast<std::uint64_t>(X.offset(k)));
            auto const [n, n_lo] = X.count();
            return detail::divide(m, m_lo, n, n_lo);
        }
        // Quantile exact : plus petit k tel que P(X <= k) >= p, pour la vraie loi,
        // soit low + ⌈p n⌉ − 1 (p = m 2^e, produit entier m n sur 128 bits). La cdf
        // calculée étant arrondie, quantile(cdf(k)) vaut k, ou k + 1 quand cdf(k)
        // a été arrondie vers le haut. (« Le plus petit k dont la cdf calculée
        // atteint p » ne serait pas tenable : au-delà de 2^53 valeurs, des milliers
        // de k voisins ont la même cdf arrondie.)
        friend constexpr I quantile(discrete_uniform const& X, double p)
        {
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return X.low_;
            if (p == 1)
                return X.high_;
            std::uint64_t const b = std::bit_cast<std::uint64_t>(p);
            auto const biased = static_cast<int>(b >> 52);
            std::uint64_t const m = biased == 0 ? (b & 0x000fffffffffffff) : ((b & 0x000fffffffffffff) | 0x0010000000000000);
            int const e = biased == 0 ? -1074 : biased - 1075;               // p = m 2^e, e < 0
            detail::uint128 const n = detail::uint128{X.span()} + 1;
            detail::uint128 const product = detail::uint128{m} * n;          // < 2^117
            int const s = -e;
            detail::uint128 const whole = s >= 128 ? 0 : product >> s;
            bool const exact = s >= 128 ? product == 0 : (product & ((detail::uint128{1} << s) - 1)) == 0;
            detail::uint128 const ceiling = whole + (exact ? 0 : 1);         // ⌈p n⌉ >= 1
            using U = std::make_unsigned_t<I>;
            return static_cast<I>(static_cast<U>(static_cast<U>(X.low_) + static_cast<U>(static_cast<std::uint64_t>(ceiling - 1))));
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr I sample(discrete_uniform const& X, G& g)
        {
            return uniform_int(g, X.low_, X.high_);
        }

        friend constexpr bool operator==(discrete_uniform const&, discrete_uniform const&) = default;

    private:
        using U = std::make_unsigned_t<I>;
        constexpr std::uint64_t span() const noexcept                       // high − low, sans débordement
        {
            return static_cast<std::uint64_t>(static_cast<U>(static_cast<U>(high_) - static_cast<U>(low_)));
        }
        constexpr std::uint64_t offset(I k) const noexcept
        {
            return static_cast<std::uint64_t>(static_cast<U>(static_cast<U>(k) - static_cast<U>(low_)));
        }
        struct exact_double
        {
            double hi, lo;
        };
        // Entier de 64 bits en double-double exact : exact en double sous 2^53 ;
        // au-delà, les 53 bits de tête (11 bits nuls en queue) et le reste.
        static constexpr exact_double split(std::uint64_t v) noexcept
        {
            if (v < (std::uint64_t{1} << 53))
                return {static_cast<double>(v), 0.0};
            return {static_cast<double>(v & ~std::uint64_t{0x7ff}), static_cast<double>(v & 0x7ff)};
        }
        // n = high − low + 1, jusqu'à 2^64.
        constexpr exact_double count() const noexcept
        {
            if (span() == ~std::uint64_t{0})
                return {0x1p64, 0.0};
            return split(span() + 1);
        }

        I low_, high_;
    };

    // -----------------------------------------------------------------------
    // Bernoulli : vrai avec probabilité p. Tirage exact.
    // -----------------------------------------------------------------------
    class bernoulli
    {
    public:
        using value_type = bool;
        static constexpr exactness sampling = exactness::exact;

        struct probability
        {
            double p;
        };

        constexpr bernoulli(probability parameter) : p_(parameter.p)
        {
            if (!(p_ >= 0 && p_ <= 1))
                throw std::domain_error("wsl::bernoulli : p doit être dans [0, 1]");
        }

        constexpr double p() const noexcept { return p_; }
        constexpr interval<bool> support() const noexcept { return {p_ == 1, p_ != 0}; }

        constexpr double mean() const noexcept { return p_; }
        // p (1 − p) = p − p², en double-double.
        constexpr double variance() const noexcept
        {
            namespace sd = detail::sd;
            if (p_ < 0x1p-480)                                               // p² sous le demi-ulp de p, et
                return p_;                                                   // produit exact loin des sous-normaux
            auto const [s, s_error] = sd::two_product(p_, p_);
            auto const [d, d_error] = sd::two_sum(p_, -s);
            return d + (d_error - s_error);
        }
        // −p log p − (1 − p) log(1 − p), log(1 − p) par log1p.
        constexpr double entropy() const noexcept
        {
            if (p_ == 0 || p_ == 1)
                return 0;
            namespace sp = wsl::special;
            return -(detail::sd::fmul(p_, sp::log(p_)) + detail::sd::fmul(1.0 - p_, sp::log1p(-p_)));
        }

        friend constexpr double pmf(bernoulli const& X, bool k) noexcept { return k ? X.p_ : 1.0 - X.p_; }
        friend constexpr double cdf(bernoulli const& X, bool k) noexcept { return k ? 1.0 : 1.0 - X.p_; }
        friend constexpr double ccdf(bernoulli const& X, bool k) noexcept { return k ? 0.0 : X.p_; }

        template <std::uniform_random_bit_generator G>
        friend constexpr bool sample(bernoulli const& X, G& g)
        {
            return bernoulli_trial(g, X.p_);
        }

        friend constexpr bool operator==(bernoulli const&, bernoulli const&) = default;

    private:
        double p_;
    };
}
