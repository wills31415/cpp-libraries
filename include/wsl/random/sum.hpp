#pragma once

// Algèbre des lois : somme de deux variables indépendantes (voir docs/random.md §4.6).
//
//   auto S = wsl::random::normal({.mean = 1, .sd = 3}) + wsl::random::normal({.mean = 2, .sd = 4});   // wsl::random::normal : N(3, 5)
//   auto P = wsl::random::poisson({.mean = 2.0}) + wsl::random::binomial({.trials = 10, .p = 0.3});  // convolution
//   pmf(P, 4); cdf(P, 6); quantile(P, 0.9); P.mean(); sample(P, g);
//
// X + Y est la loi de X + Y, X et Y indépendantes : X + X est donc la somme de deux
// copies indépendantes (de loi celle de 2X seulement pour une loi stable).
//
// Identités, par la fonction sum_closure(X, Y) trouvée par ADL :
//   - rendant une loi : toujours vraie, X + Y est de ce type (normale + normale,
//     Cauchy + Cauchy, Poisson + Poisson, χ² + χ²), paramètres arrondis une fois ;
//   - rendant std::optional<R> : vraie sous condition (gamma + gamma de même
//     intensité, exponentielles, binomiales de même p) ; X + Y est independent_sum,
//     dont les fonctions sont celles de R quand la condition tient.
// Sinon independent_sum compose : moments par les cumulants (κ₁, κ₂, κ₃, κ₄
// s'ajoutent) ; lois discrètes : masse, répartition et survie par convolution exacte
// (sommes finies de termes positifs, en double-double, O(longueur du support utile)),
// quantile exact de cette répartition. Une somme continue sans identité n'a ni
// densité ni répartition (pas de convolution numérique) : ses fonctions rendent NaN,
// et quantile lève std::domain_error, quand la condition d'une identité échoue.
//
// Tirage : x + y, x et y tirés dans cet ordre (exact pour des lois discrètes exactes ;
// arrondi une fois pour des réels).

#include <wsl/random/affine.hpp>
#include <wsl/random/binomial.hpp>
#include <wsl/random/cauchy.hpp>
#include <wsl/random/chi_squared.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/discrete.hpp>
#include <wsl/random/exponential.hpp>
#include <wsl/random/gamma.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/random/poisson.hpp>
#include <wsl/special.hpp>

#include <algorithm>
#include <concepts>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace wsl::random
{
    // -----------------------------------------------------------------------
    // Identités.
    // -----------------------------------------------------------------------
    template <law_real T>
    constexpr normal<T> sum_closure(normal<T> const& X, normal<T> const& Y)
    {
        return normal<T>({.mean = static_cast<T>(double(X.mean()) + double(Y.mean())), .sd = static_cast<T>(wsl::special::hypot(double(X.sd()), double(Y.sd())))});
    }
    template <law_real T>
    constexpr cauchy<T> sum_closure(cauchy<T> const& X, cauchy<T> const& Y)
    {
        return cauchy<T>({.location = static_cast<T>(double(X.location()) + double(Y.location())), .scale = static_cast<T>(double(X.scale()) + double(Y.scale()))});
    }
    template <law_integer I>
    constexpr poisson<I> sum_closure(poisson<I> const& X, poisson<I> const& Y)
    {
        return poisson<I>({.mean = X.mean() + Y.mean()});
    }
    template <law_real T>
    constexpr chi_squared<T> sum_closure(chi_squared<T> const& X, chi_squared<T> const& Y)
    {
        return chi_squared<T>({.degrees_of_freedom = static_cast<T>(double(X.degrees_of_freedom()) + double(Y.degrees_of_freedom()))});
    }
    // E(λ) = gamma(1, λ).
    template <law_real T>
    constexpr std::optional<gamma<T>> sum_closure(exponential<T> const& X, exponential<T> const& Y) noexcept
    {
        if (!(X.rate() == Y.rate()))
            return std::nullopt;
        return gamma<T>({.shape = T(2), .rate = X.rate()});
    }
    template <law_real T>
    constexpr std::optional<gamma<T>> sum_closure(gamma<T> const& X, exponential<T> const& Y) noexcept
    {
        return sum_closure(X, gamma<T>({.shape = T(1), .rate = Y.rate()}));
    }
    template <law_real T>
    constexpr std::optional<gamma<T>> sum_closure(exponential<T> const& X, gamma<T> const& Y) noexcept
    {
        return sum_closure(Y, gamma<T>({.shape = T(1), .rate = X.rate()}));
    }
    template <law_integer I>
    constexpr std::optional<binomial<I>> sum_closure(binomial<I> const& X, binomial<I> const& Y) noexcept
    {
        if (!(X.p() == Y.p()) || X.trials() > std::numeric_limits<I>::max() - Y.trials())
            return std::nullopt;
        return binomial<I>({.trials = static_cast<I>(X.trials() + Y.trials()), .p = X.p()});
    }

    namespace detail
    {
        template <class R>
        struct optional_value
        {
            using type = void;
        };
        template <class R>
        struct optional_value<std::optional<R>>
        {
            using type = R;
        };

        // Somme saturée d'entiers.
        template <std::integral I>
        constexpr I saturated_sum(I a, I b) noexcept
        {
            if (b > 0 && a > std::numeric_limits<I>::max() - b)
                return std::numeric_limits<I>::max();
            if (b < 0 && a < std::numeric_limits<I>::min() - b)
                return std::numeric_limits<I>::min();
            return static_cast<I>(a + b);
        }
        template <std::integral I>
        constexpr I saturated_difference(I a, I b) noexcept
        {
            if (b < 0 && a > std::numeric_limits<I>::max() + b)
                return std::numeric_limits<I>::max();
            if (b > 0 && a < std::numeric_limits<I>::min() + b)
                return std::numeric_limits<I>::min();
            return static_cast<I>(a - b);
        }
    }

    template <class L1, class L2>
    concept has_sum_closure = requires(L1 const& x, L2 const& y) { sum_closure(x, y); };

    template <class L1, class L2>
    using sum_closure_t = decltype(sum_closure(std::declval<L1 const&>(), std::declval<L2 const&>()));

    // Identité toujours vraie : X + Y est de ce type.
    template <class L1, class L2>
    concept sum_closed = has_sum_closure<L1, L2> && probability_law<sum_closure_t<L1, L2>>;

    // -----------------------------------------------------------------------
    // Somme indépendante générique.
    // -----------------------------------------------------------------------
    template <probability_law L1, probability_law L2>
        requires std::same_as<typename L1::value_type, typename L2::value_type>
    class independent_sum
    {
        template <class A, class B>
        struct family_of
        {
            using type = void;
        };
        template <class A, class B>
            requires has_sum_closure<A, B>
        struct family_of<A, B>
        {
            using type = typename detail::optional_value<sum_closure_t<A, B>>::type;
        };

    public:
        using value_type = typename L1::value_type;
        using family_type = typename family_of<L1, L2>::type;              // void sans identité conditionnelle
        static constexpr bool continuous = std::is_floating_point_v<value_type>;
        static constexpr bool has_family = !std::is_void_v<family_type>;
        using moment_type = std::conditional_t<continuous, value_type, double>;
        static constexpr exactness sampling =
            !continuous && L1::sampling == exactness::exact && L2::sampling == exactness::exact ? exactness::exact : exactness::approximate;

        constexpr independent_sum(L1 const& x, L2 const& y) : x_(x), y_(y)
        {
            if constexpr (has_family)
                family_ = sum_closure(x_, y_);
        }

        constexpr L1 const& first() const noexcept { return x_; }
        constexpr L2 const& second() const noexcept { return y_; }
        // La loi de la famille quand l'identité conditionnelle tient.
        constexpr auto const& family() const noexcept requires has_family { return family_; }

        constexpr interval<value_type> support() const noexcept
        {
            auto const a = x_.support(), b = y_.support();
            return {add(a.lower, b.lower), add(a.upper, b.upper)};
        }

        // Cumulants : κ₁ et κ₂ s'ajoutent…
        constexpr moment_type mean() const noexcept requires has_mean<L1> && has_mean<L2>
        {
            if constexpr (has_family && has_mean<family_type>)
                if (family_)
                    return static_cast<moment_type>(family_->mean());
            return static_cast<moment_type>(double(x_.mean()) + double(y_.mean()));
        }
        constexpr moment_type variance() const noexcept requires has_variance<L1> && has_variance<L2>
        {
            if constexpr (has_family && has_variance<family_type>)
                if (family_)
                    return static_cast<moment_type>(family_->variance());
            return static_cast<moment_type>(double(x_.variance()) + double(y_.variance()));
        }
        // … ainsi que κ₃ = γ₁ σ³ et κ₄ = γ₂ σ⁴.
        constexpr moment_type skewness() const noexcept requires has_variance<L1> && has_variance<L2> && has_skewness<L1> && has_skewness<L2>
        {
            namespace fp = detail::sd;
            if constexpr (has_family && has_skewness<family_type>)
                if (family_)
                    return static_cast<moment_type>(family_->skewness());
            double const v1 = x_.variance(), v2 = y_.variance(), v = v1 + v2;
            double const k3 = fp::fmul(double(x_.skewness()), fp::fmul(v1, wsl::special::sqrt(v1))) + fp::fmul(double(y_.skewness()), fp::fmul(v2, wsl::special::sqrt(v2)));
            return static_cast<moment_type>(k3 / (v * wsl::special::sqrt(v)));
        }
        constexpr moment_type excess_kurtosis() const noexcept requires has_variance<L1> && has_variance<L2> && has_excess_kurtosis<L1> && has_excess_kurtosis<L2>
        {
            namespace fp = detail::sd;
            if constexpr (has_family && has_excess_kurtosis<family_type>)
                if (family_)
                    return static_cast<moment_type>(family_->excess_kurtosis());
            double const v1 = x_.variance(), v2 = y_.variance(), v = v1 + v2;
            double const k4 = fp::fmul(double(x_.excess_kurtosis()), fp::fmul(v1, v1)) + fp::fmul(double(y_.excess_kurtosis()), fp::fmul(v2, v2));
            return static_cast<moment_type>(k4 / fp::fmul(v, v));
        }

        // ---- Lois continues : seulement par une identité conditionnelle ----
        friend constexpr value_type pdf(independent_sum const& S, value_type x) noexcept requires continuous && has_family && has_pdf<family_type>
        {
            return S.family_ ? pdf(*S.family_, x) : std::numeric_limits<value_type>::quiet_NaN();
        }
        friend constexpr value_type log_pdf(independent_sum const& S, value_type x) noexcept
            requires continuous && has_family && requires(family_type const& F, value_type v) { log_pdf(F, v); }
        {
            return S.family_ ? log_pdf(*S.family_, x) : std::numeric_limits<value_type>::quiet_NaN();
        }
        friend constexpr value_type cdf(independent_sum const& S, value_type x) noexcept requires continuous && has_family && has_cdf<family_type>
        {
            return S.family_ ? cdf(*S.family_, x) : std::numeric_limits<value_type>::quiet_NaN();
        }
        friend constexpr value_type ccdf(independent_sum const& S, value_type x) noexcept requires continuous && has_family && has_ccdf<family_type>
        {
            return S.family_ ? ccdf(*S.family_, x) : std::numeric_limits<value_type>::quiet_NaN();
        }
        friend constexpr value_type quantile(independent_sum const& S, double p) requires continuous && has_family && has_quantile<family_type>
        {
            if (!S.family_)
                throw std::domain_error("wsl::random::quantile : somme indépendante sans forme connue");
            return quantile(*S.family_, p);
        }
        friend constexpr value_type cquantile(independent_sum const& S, double q) requires continuous && has_family && has_cquantile<family_type>
        {
            if (!S.family_)
                throw std::domain_error("wsl::random::cquantile : somme indépendante sans forme connue");
            return cquantile(*S.family_, q);
        }

        // ---- Lois discrètes : convolution ----
        static constexpr bool convolvable = !continuous && has_pmf<L1> && has_pmf<L2> && has_cdf<L1> && has_cdf<L2> && has_ccdf<L1> && has_ccdf<L2>;

        // Σⱼ P(X = j) P(Y = k − j)
        friend constexpr double pmf(independent_sum const& S, value_type k) noexcept requires convolvable
        {
            if constexpr (has_family)
                if (S.family_)
                    return pmf(*S.family_, k);
            auto const [lo, hi] = S.range(k);
            detail::sd::double_double s{0.0, 0.0};
            for (value_type j = lo; j <= hi; ++j)
            {
                s = detail::sum(s, detail::product(double(pmf(S.x_, j)), double(pmf(S.y_, static_cast<value_type>(k - j)))));
                if (j == hi)
                    break;
            }
            return s.hi + s.lo;
        }
        // P(X <= k − sup Y − 1) + Σⱼ P(X = j) P(Y <= k − j)
        friend constexpr double cdf(independent_sum const& S, value_type k) noexcept requires convolvable
        {
            if constexpr (has_family)
                if (S.family_)
                    return cdf(*S.family_, k);
            auto const s = S.cdf_dd(k);
            return clamp01(s.hi + s.lo);
        }
        // Σⱼ P(X = j) P(Y > k − j) + P(X > k − inf Y)
        friend constexpr double ccdf(independent_sum const& S, value_type k) noexcept requires convolvable
        {
            if constexpr (has_family)
                if (S.family_)
                    return ccdf(*S.family_, k);
            auto const s = S.ccdf_dd(k);
            return clamp01(s.hi + s.lo);
        }
        // Plus petit k tel que P(X + Y <= k) >= p.
        friend constexpr value_type quantile(independent_sum const& S, double p) requires convolvable
        {
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::random::quantile : p doit être dans [0, 1]");
            if constexpr (has_family && has_quantile<family_type>)
                if (S.family_)
                    return quantile(*S.family_, p);
            auto const [low, high] = S.support();
            if (p == 0)
                return low;
            if (p == 1)
                return high;
            value_type guess = low;
            if constexpr (has_quantile<L1> && has_quantile<L2>)
                guess = detail::saturated_sum(quantile(S.x_, p), quantile(S.y_, p));
            double const q = 1.0 - p;                                       // utilisé pour p > ½ : exact
            auto const ge = [&](value_type k) { return p <= 0.5 ? detail::at_least(S.cdf_dd(k), p) : detail::at_most(S.ccdf_dd(k), q); };
            return detail::smallest_true(guess, low, high, ge);
        }

        template <std::uniform_random_bit_generator G>
            requires samplable<L1, G> && samplable<L2, G>
        friend constexpr value_type sample(independent_sum const& S, G& g)
        {
            value_type const a = sample(S.x_, g);
            value_type const b = sample(S.y_, g);
            return add(a, b);
        }
        template <std::uniform_random_bit_generator G>
            requires samplable<L1, G> && samplable<L2, G>
        friend constexpr value_type sample(independent_sum const& S, G& g, fast_t)
        {
            value_type const a = sample(S.x_, g, fast);
            value_type const b = sample(S.y_, g, fast);
            return add(a, b);
        }

        friend constexpr bool operator==(independent_sum const& a, independent_sum const& b) noexcept { return a.x_ == b.x_ && a.y_ == b.y_; }

    private:
        static constexpr double clamp01(double v) noexcept { return v < 0 ? 0.0 : v > 1 ? 1.0 : v; }

        static constexpr value_type add(value_type a, value_type b) noexcept
        {
            if constexpr (continuous)
                return a + b;
            else
                return detail::saturated_sum(a, b);
        }

        // j tels que P(X = j) P(Y = k − j) peut être non nul : [max(inf X, k − sup Y), min(sup X, k − inf Y)].
        constexpr interval<value_type> range(value_type k) const noexcept
        {
            auto const a = x_.support(), b = y_.support();
            value_type const lo = std::max(a.lower, detail::saturated_difference(k, b.upper));
            value_type const hi = std::min(a.upper, detail::saturated_difference(k, b.lower));
            return {lo, hi};
        }

        // Somme de termes positifs, en double-double.
        template <class F>
        constexpr detail::sd::double_double convolve(value_type k, F f) const
        {
            auto const [lo, hi] = range(k);
            detail::sd::double_double s{0.0, 0.0};
            for (value_type j = lo; j <= hi; ++j)
            {
                s = detail::sum(s, detail::product(double(pmf(x_, j)), f(static_cast<value_type>(k - j))));
                if (j == hi)
                    break;
            }
            return s;
        }
        constexpr detail::sd::double_double cdf_dd(value_type k) const
        {
            auto s = convolve(k, [&](value_type m) { return double(cdf(y_, m)); });
            auto const b = y_.support().upper;
            if (b < std::numeric_limits<value_type>::max())                // j < k − sup Y : P(Y <= k − j) = 1
            {
                value_type const edge = detail::saturated_difference(detail::saturated_difference(k, b), value_type(1));
                if (edge >= x_.support().lower)
                    s = detail::sum(s, {double(cdf(x_, edge)), 0.0});
            }
            return s;
        }
        constexpr detail::sd::double_double ccdf_dd(value_type k) const
        {
            auto s = convolve(k, [&](value_type m) { return double(ccdf(y_, m)); });
            value_type const edge = detail::saturated_difference(k, y_.support().lower);   // j > k − inf Y : P(Y > k − j) = 1
            if (edge < x_.support().upper)
                s = detail::sum(s, {double(ccdf(x_, edge)), 0.0});
            return s;
        }

        L1 x_;
        L2 y_;
        struct no_family
        {
            friend constexpr bool operator==(no_family, no_family) noexcept = default;
        };
        [[no_unique_address]] std::conditional_t<has_family, std::optional<std::conditional_t<has_family, family_type, int>>, no_family> family_{};
    };

    // X + Y, X et Y indépendantes : la loi de la famille si l'identité est toujours vraie.
    template <probability_law L1, probability_law L2>
        requires std::same_as<typename L1::value_type, typename L2::value_type>
    constexpr auto operator+(L1 const& X, L2 const& Y)
    {
        if constexpr (sum_closed<L1, L2>)
            return sum_closure(X, Y);
        else
            return independent_sum<L1, L2>(X, Y);
    }
    // X − Y = X + (−Y), lois réelles.
    template <real_law L1, real_law L2>
        requires std::same_as<typename L1::value_type, typename L2::value_type>
    constexpr auto operator-(L1 const& X, L2 const& Y)
    {
        return X + (-Y);
    }
}
