#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Algèbre des lois : image d'une loi par une fonction, Y = f(X) (voir design/random.md §4.6).
//
//   namespace sp = wsl::special;
//   auto Q = wsl::probability::pushforward(wsl::probability::normal({.mean = 0, .stddev = 1}), [](double x) { return x * x; });   // tirage seul
//   auto L = wsl::probability::pushforward(wsl::probability::normal({.mean = 0, .stddev = 1}),
//                             wsl::probability::increasing([](double x) { return sp::exp(x); },          // f
//                                             [](double y) { return sp::log(y); },          // f⁻¹
//                                             [](double y) { return 1 / y; }));             // (f⁻¹)'
//   L.cdf(2.0); L.quantile(0.9); L.pdf(1.0);       // la log-normale, par composition
//
// Fonction quelconque : seul le tirage f(x) est hérité ; support : tout le type.
// Bijection monotone déclarée (increasing ou decreasing, f, f⁻¹ et facultativement la
// dérivée de f⁻¹) : répartition F_X(f⁻¹(y)) (survie si f décroît), quantiles
// f(Q_X(p)) (f(cQ_X(p)) si f décroît), médiane, support f(support) ; densité
// f_X(f⁻¹(y)) |(f⁻¹)'(y)| et log-densité si la dérivée est donnée. Composition : la
// précision et la reproductibilité au bit près dépendent de celles de f (wsl::special
// les garantit, <cmath> non). f est appliquée aux bornes du support (±∞ comprises).

#include <wsl/probability/affine.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/special.hpp>

#include <concepts>
#include <functional>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>

namespace wsl::inline WSL_DEFAULTS::probability
{
    struct no_derivative
    {
    };

    // Bijection monotone : f, sa réciproque et, facultativement, la dérivée de la réciproque.
    template <bool Increasing, class F, class G, class D = no_derivative>
    struct monotone_map
    {
        static constexpr bool increasing = Increasing;
        static constexpr bool has_derivative = !std::is_same_v<D, no_derivative>;
        F forward;
        G inverse;
        D inverse_derivative;
    };

    template <class F, class G>
    constexpr monotone_map<true, F, G> increasing(F f, G inverse) { return {f, inverse, {}}; }
    template <class F, class G, class D>
    constexpr monotone_map<true, F, G, D> increasing(F f, G inverse, D inverse_derivative) { return {f, inverse, inverse_derivative}; }
    template <class F, class G>
    constexpr monotone_map<false, F, G> decreasing(F f, G inverse) { return {f, inverse, {}}; }
    template <class F, class G, class D>
    constexpr monotone_map<false, F, G, D> decreasing(F f, G inverse, D inverse_derivative) { return {f, inverse, inverse_derivative}; }

    namespace detail
    {
        template <class T>
        struct is_monotone_map : std::false_type
        {
        };
        template <bool I, class F, class G, class D>
        struct is_monotone_map<monotone_map<I, F, G, D>> : std::true_type
        {
        };

        template <class T>
        constexpr interval<T> whole_type() noexcept
        {
            if constexpr (std::is_floating_point_v<T>)
                return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()};
            else
                return {std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max()};
        }
    }

    // Image par une fonction quelconque : tirage seul.
    template <probability_law L, class F>
        requires std::invocable<F const&, typename L::value_type> && std::is_arithmetic_v<std::invoke_result_t<F const&, typename L::value_type>>
    class image : public law_interface
    {
    public:
        using value_type = std::invoke_result_t<F const&, typename L::value_type>;
        static constexpr exactness sampling = exactness::approximate;

        constexpr image(L const& law, F f) : law_(law), f_(std::move(f)) {}

        constexpr L const& law() const noexcept { return law_; }
        constexpr interval<value_type> support() const noexcept { return detail::whole_type<value_type>(); }

        template <std::uniform_random_bit_generator G>
            requires samplable<L, G>
        constexpr value_type sample(G& g) const
        {
            auto const& Y = *this;
            return std::invoke(Y.f_, wsl::probability::sample(Y.law_, g));
        }
        template <std::uniform_random_bit_generator G>
            requires samplable<L, G>
        constexpr value_type sample(G& g, fast_t) const
        {
            auto const& Y = *this;
            return std::invoke(Y.f_, wsl::probability::sample(Y.law_, g, fast));
        }

    private:
        L law_;
        F f_;
    };

    // Image par une bijection monotone déclarée.
    template <real_law L, class M>
        requires detail::is_monotone_map<M>::value
    class monotone_image : public law_interface
    {
        using T = typename L::value_type;

    public:
        using value_type = std::invoke_result_t<decltype(M::forward) const&, T>;
        static_assert(std::is_floating_point_v<value_type>, "wsl::probability::pushforward : f doit être à valeurs réelles");
        static constexpr exactness sampling = exactness::approximate;
        static constexpr bool increasing = M::increasing;

        constexpr monotone_image(L const& law, M m) : law_(law), m_(std::move(m)) {}

        constexpr L const& law() const noexcept { return law_; }
        constexpr interval<value_type> support() const noexcept
        {
            auto const [lo, up] = law_.support();
            value_type const a = forward(lo), b = forward(up);
            return increasing ? interval<value_type>{a, b} : interval<value_type>{b, a};
        }
        constexpr value_type median() const requires has_median<L> { return forward(law_.median()); }

        constexpr value_type cdf(value_type y) const requires has_cdf<L> && has_ccdf<L>
        {
            auto const& Y = *this;
            T const x = Y.inverse(y);
            return static_cast<value_type>(increasing ? wsl::probability::cdf(Y.law_, x) : wsl::probability::ccdf(Y.law_, x));
        }
        constexpr value_type ccdf(value_type y) const requires has_cdf<L> && has_ccdf<L>
        {
            auto const& Y = *this;
            T const x = Y.inverse(y);
            return static_cast<value_type>(increasing ? wsl::probability::ccdf(Y.law_, x) : wsl::probability::cdf(Y.law_, x));
        }
        constexpr value_type quantile(double p) const requires has_quantile<L> && has_cquantile<L>
        {
            auto const& Y = *this;
            return Y.forward(increasing ? wsl::probability::quantile(Y.law_, p) : wsl::probability::cquantile(Y.law_, p));
        }
        constexpr value_type cquantile(double q) const requires has_quantile<L> && has_cquantile<L>
        {
            auto const& Y = *this;
            return Y.forward(increasing ? wsl::probability::cquantile(Y.law_, q) : wsl::probability::quantile(Y.law_, q));
        }
        // f_X(f⁻¹(y)) |(f⁻¹)'(y)|, 0 hors de l'image du support.
        constexpr value_type pdf(value_type y) const requires M::has_derivative && has_pdf<L>
        {
            auto const& Y = *this;
            auto const [lo, up] = Y.support();
            if (!(y >= lo && y <= up))
                return y == y ? value_type(0) : y;
            double const d = std::invoke(Y.m_.inverse_derivative, y);
            return static_cast<value_type>(detail::sd::fmul(double(wsl::probability::pdf(Y.law_, Y.inverse(y))), d < 0 ? -d : d));
        }
        constexpr value_type log_pdf(value_type y) const
            requires M::has_derivative && requires(L const& X, T x) { wsl::probability::log_pdf(X, x); }
        {
            auto const& Y = *this;
            auto const [lo, up] = Y.support();
            if (!(y >= lo && y <= up))
                return y == y ? -std::numeric_limits<value_type>::infinity() : y;
            double const d = std::invoke(Y.m_.inverse_derivative, y);
            return static_cast<value_type>(double(wsl::probability::log_pdf(Y.law_, Y.inverse(y))) + wsl::special::log(d < 0 ? -d : d));
        }

        template <std::uniform_random_bit_generator G>
            requires samplable<L, G>
        constexpr value_type sample(G& g) const
        {
            auto const& Y = *this;
            return Y.forward(wsl::probability::sample(Y.law_, g));
        }
        template <std::uniform_random_bit_generator G>
            requires samplable<L, G>
        constexpr value_type sample(G& g, fast_t) const
        {
            auto const& Y = *this;
            return Y.forward(wsl::probability::sample(Y.law_, g, fast));
        }

    private:
        constexpr value_type forward(T x) const { return static_cast<value_type>(std::invoke(m_.forward, x)); }
        constexpr T inverse(value_type y) const { return static_cast<T>(std::invoke(m_.inverse, y)); }

        L law_;
        M m_;
    };

    // Loi de f(X).
    template <probability_law L, class F>
    constexpr auto pushforward(L const& X, F f)
    {
        if constexpr (detail::is_monotone_map<F>::value)
            return monotone_image<L, F>(X, std::move(f));
        else
            return image<L, F>(X, std::move(f));
    }
}
