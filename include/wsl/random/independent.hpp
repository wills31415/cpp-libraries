#pragma once

// Algèbre des lois : loi jointe de variables indépendantes (loi produit, voir
// docs/random.md §4.6).
//
//   wsl::random::independent J(wsl::random::normal({.mean = 0, .sd = 1}), wsl::random::exponential({.rate = 2}));
//   auto [x, y] = J(g);                       // tirés dans cet ordre
//   J.pdf({0.5, 1.0}); J.cdf({0.0, 1.0}); J.mean(); J.entropy();
//
// value_type : std::tuple des valeurs des composantes. Densité (lois toutes continues)
// ou masse (toutes discrètes) : produit, arrondi une fois ; log-densité : somme ;
// répartition jointe P(X₁ <= x₁, …) et survie jointe P(X₁ > x₁, …) : produits ;
// entropie : somme ; moyenne, variance : tuples. Support : le pavé des supports.
// Tirage exact si celui de chaque composante l'est.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/law.hpp>

#include <cstddef>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>

namespace wsl::random
{
    template <probability_law... L>
        requires(sizeof...(L) >= 1)
    class independent : public law_interface
    {
        static constexpr bool all_continuous = (std::is_floating_point_v<typename L::value_type> && ...);
        static constexpr bool all_discrete = (std::is_integral_v<typename L::value_type> && ...);

    public:
        using value_type = std::tuple<typename L::value_type...>;
        static constexpr exactness sampling = ((L::sampling == exactness::exact) && ...) ? exactness::exact : exactness::approximate;

        constexpr explicit independent(L const&... laws) : laws_(laws...) {}

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
        constexpr double entropy() const noexcept requires(has_entropy<L> && ...)
        {
            detail::sd::double_double s{0.0, 0.0};
            std::apply([&](auto const&... c) { ((s = detail::sum(s, {double(c.entropy()), 0.0})), ...); }, laws_);
            return s.hi + s.lo;
        }

        constexpr double pdf(value_type const& x) const noexcept requires all_continuous && (has_pdf<L> && ...)
        {
            auto const& J = *this;
            return J.product(x, [](auto const& c, auto v) { return double(wsl::random::pdf(c, v)); });
        }
        constexpr double pmf(value_type const& x) const noexcept requires all_discrete && (has_pmf<L> && ...)
        {
            auto const& J = *this;
            return J.product(x, [](auto const& c, auto v) { return double(wsl::random::pmf(c, v)); });
        }
        constexpr double log_pdf(value_type const& x) const noexcept
            requires all_continuous && (requires(L const& c, typename L::value_type v) { wsl::random::log_pdf(c, v); } && ...)
        {
            auto const& J = *this;
            return J.total(x, [](auto const& c, auto v) { return double(wsl::random::log_pdf(c, v)); });
        }
        constexpr double log_pmf(value_type const& x) const noexcept
            requires all_discrete && (requires(L const& c, typename L::value_type v) { wsl::random::log_pmf(c, v); } && ...)
        {
            auto const& J = *this;
            return J.total(x, [](auto const& c, auto v) { return double(wsl::random::log_pmf(c, v)); });
        }
        // P(X₁ <= x₁, …, Xₙ <= xₙ)
        constexpr double cdf(value_type const& x) const noexcept requires(has_cdf<L> && ...)
        {
            auto const& J = *this;
            return J.product(x, [](auto const& c, auto v) { return double(wsl::random::cdf(c, v)); });
        }
        // P(X₁ > x₁, …, Xₙ > xₙ)
        constexpr double ccdf(value_type const& x) const noexcept requires(has_ccdf<L> && ...)
        {
            auto const& J = *this;
            return J.product(x, [](auto const& c, auto v) { return double(wsl::random::ccdf(c, v)); });
        }

        template <std::uniform_random_bit_generator G>
            requires(samplable<L, G> && ...)
        constexpr value_type sample(G& g) const
        {
            auto const& J = *this;
            return std::apply([&g](auto const&... c) { return value_type{wsl::random::sample(c, g)...}; }, J.laws_);   // accolades : dans l'ordre
        }
        template <std::uniform_random_bit_generator G>
            requires(samplable<L, G> && ...)
        constexpr value_type sample(G& g, fast_t) const
        {
            auto const& J = *this;
            return std::apply([&g](auto const&... c) { return value_type{wsl::random::sample(c, g, fast)...}; }, J.laws_);
        }

        friend constexpr bool operator==(independent const& a, independent const& b) noexcept { return a.laws_ == b.laws_; }

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
    independent(L...) -> independent<L...>;
}
