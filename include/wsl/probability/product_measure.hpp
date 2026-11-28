#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Algèbre des lois : loi jointe de variables indépendantes (loi produit, voir
// design/random.md §4.6).
//
//   wsl::probability::product_measure J(wsl::probability::normal({.mean = 0, .stddev = 1}), wsl::probability::exponential({.rate = 2}));
//   auto [x, y] = J(g);                       // tirés dans cet ordre
//   J.pdf({0.5, 1.0}); J.cdf({0.0, 1.0}); J.mean(); J.entropy();
//
// value_type : std::tuple des valeurs des composantes. Densité (lois toutes continues)
// ou masse (toutes discrètes) : produit, arrondi une fois ; log-densité : somme ;
// répartition jointe P(X₁ <= x₁, …) et survie jointe P(X₁ > x₁, …) : produits ;
// entropie : somme ; moyenne, variance : tuples. Support : le pavé des supports.
// Tirage exact si celui de chaque composante l'est.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/law.hpp>

#include <cstddef>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>

namespace wsl::inline WSL_DEFAULTS::probability
{
    template <probability_law... L>
        requires(sizeof...(L) >= 1)
    class product_measure : public law_interface
    {
        static constexpr bool all_continuous = (std::is_floating_point_v<typename L::value_type> && ...);
        static constexpr bool all_discrete = (std::is_integral_v<typename L::value_type> && ...);

    public:
        using value_type = std::tuple<typename L::value_type...>;
        static constexpr exactness sampling = ((L::sampling == exactness::exact) && ...) ? exactness::exact : exactness::approximate;

        constexpr explicit product_measure(L const&... laws) : laws_(laws...) {}

        static constexpr std::size_t size() noexcept { return sizeof...(L); }
        template <std::size_t I>
        constexpr auto const& component() const noexcept { return std::get<I>(laws_); }

        constexpr interval<value_type> support() const noexcept
        {
            return std::apply([](auto const&... c) { return interval<value_type>{value_type{c.support().lower...}, value_type{c.support().upper...}}; }, laws_);
        }

        constexpr auto mean() const noexcept requires(has_mean<L> && ...)
        {
            return std::apply([](auto const&... c) { return std::tuple{c.mean()...}; }, laws_);
        }
        constexpr auto variance() const noexcept requires(has_variance<L> && ...)
        {
            return std::apply([](auto const&... c) { return std::tuple{c.variance()...}; }, laws_);
        }
        // Σ H(Xᵢ)
        // Type commun des entropies des composantes (T des lois continues, double des discrètes).
        constexpr auto entropy() const noexcept requires(has_entropy<L> && ...)
        {
            using R = std::common_type_t<decltype(std::declval<L const&>().entropy())...>;
            detail::sd::double_double s{0.0, 0.0};
            std::apply([&](auto const&... c) { ((s = detail::sum(s, {double(c.entropy()), 0.0})), ...); }, laws_);
            return static_cast<R>(s.hi + s.lo);
        }

        constexpr double pdf(value_type const& x) const noexcept requires all_continuous && (has_pdf<L> && ...)
        {
            auto const& J = *this;
            return J.product(x, [](auto const& c, auto v) { return double(wsl::probability::pdf(c, v)); });
        }
        constexpr double pmf(value_type const& x) const noexcept requires all_discrete && (has_pmf<L> && ...)
        {
            auto const& J = *this;
            return J.product(x, [](auto const& c, auto v) { return double(wsl::probability::pmf(c, v)); });
        }
        constexpr double log_pdf(value_type const& x) const noexcept
            requires all_continuous && (requires(L const& c, typename L::value_type v) { wsl::probability::log_pdf(c, v); } && ...)
        {
            auto const& J = *this;
            return J.total(x, [](auto const& c, auto v) { return double(wsl::probability::log_pdf(c, v)); });
        }
        constexpr double log_pmf(value_type const& x) const noexcept
            requires all_discrete && (requires(L const& c, typename L::value_type v) { wsl::probability::log_pmf(c, v); } && ...)
        {
            auto const& J = *this;
            return J.total(x, [](auto const& c, auto v) { return double(wsl::probability::log_pmf(c, v)); });
        }
        // P(X₁ <= x₁, …, Xₙ <= xₙ)
        constexpr double cdf(value_type const& x) const noexcept requires(has_cdf<L> && ...)
        {
            auto const& J = *this;
            return J.product(x, [](auto const& c, auto v) { return double(wsl::probability::cdf(c, v)); });
        }
        // P(X₁ > x₁, …, Xₙ > xₙ)
        constexpr double ccdf(value_type const& x) const noexcept requires(has_ccdf<L> && ...)
        {
            auto const& J = *this;
            return J.product(x, [](auto const& c, auto v) { return double(wsl::probability::ccdf(c, v)); });
        }

        template <std::uniform_random_bit_generator G>
            requires(samplable<L, G> && ...)
        constexpr value_type sample(G& g) const
        {
            auto const& J = *this;
            return std::apply([&g](auto const&... c) { return value_type{detail::native_sample(c, g)...}; }, J.laws_);   // accolades : dans l'ordre
        }
        template <std::uniform_random_bit_generator G>
            requires(samplable<L, G> && ...)
        constexpr value_type sample(G& g, approximate_t) const
        {
            auto const& J = *this;
            return std::apply([&g](auto const&... c) { return value_type{wsl::probability::sample(c, g, wsl::approximate)...}; }, J.laws_);
        }

        friend constexpr bool operator==(product_measure const& a, product_measure const& b) noexcept { return a.laws_ == b.laws_; }

    private:
        // Π f(Xᵢ, xᵢ), à l'échelle, arrondi une fois (0 si un facteur est nul).
        template <class F>
        constexpr double product(value_type const& x, F f) const noexcept
        {
            double result = 1;
            [&]<std::size_t... I>(std::index_sequence<I...>) {
                detail::scaled p = detail::normalized(1.0, 0.0);
                bool zero = false, special = false;
                auto const step = [&](double v) {
                    if (v == 0)
                        zero = true;
                    else if (!(v > 0) || !detail::is_finite(v))
                    {
                        special = true;
                        result = v;
                    }
                    else
                        p = p * detail::normalized(v, 0.0);
                };
                (step(f(std::get<I>(laws_), std::get<I>(x))), ...);
                if (!special)
                    result = zero ? 0.0 : detail::rounded(p);
            }(std::index_sequence_for<L...>{});
            return result;
        }
        // Σ f(Xᵢ, xᵢ), en double-double.
        template <class F>
        constexpr double total(value_type const& x, F f) const noexcept
        {
            detail::sd::double_double s{0.0, 0.0};
            bool infinite = false;
            double special = 0;
            [&]<std::size_t... I>(std::index_sequence<I...>) {
                auto const step = [&](double v) {
                    if (!detail::is_finite(v))
                    {
                        infinite = true;
                        special = special == 0 || v != v ? v : special + v;
                    }
                    else
                        s = detail::sum(s, {v, 0.0});
                };
                (step(f(std::get<I>(laws_), std::get<I>(x))), ...);
            }(std::index_sequence_for<L...>{});
            return infinite ? special : s.hi + s.lo;
        }

        std::tuple<L...> laws_;
    };

    template <class... L>
    product_measure(L...) -> product_measure<L...>;
}
