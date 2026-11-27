#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Processus aléatoires : trajectoires paresseuses (voir design/random.md §14).
//
//   for (auto x : wsl::probability::trajectory(P, g) | std::views::take(100)) …
//
// Un processus P fournit un « marcheur » (P.walker()), objet qui tire l'état initial
// (first(g)) puis chaque état suivant (next(g)). trajectory(P, g) en fait une vue
// std::ranges (input_range, infinie) qui tient le moteur par référence : il doit lui
// survivre, et elle l'avance.
//
// Contrat de reproductibilité : la valeur d'indice k ne dépend que de l'état initial du
// moteur et de k ; une valeur n'est tirée qu'au moment où on la lit (déréférencement),
// jamais par avance : views::take(k) consomme exactement ce que demandent les k valeurs,
// et tronquer une trajectoire ne change pas les valeurs gardées. Les restrictions finies
// (P.steps(n)…) sont des lois ordinaires dont le tirage suit le même marcheur : leurs
// valeurs sont les n + 1 premières de la trajectoire.

#include <wsl/probability/law.hpp>

#include <concepts>
#include <cstddef>
#include <iterator>
#include <random>
#include <ranges>
#include <type_traits>
#include <utility>

namespace wsl::inline WSL_DEFAULTS::probability::stochastic_processes
{
    template <class W, std::uniform_random_bit_generator G>
    class process_view : public std::ranges::view_interface<process_view<W, G>>
    {
    public:
        using value_type = std::remove_cvref_t<decltype(std::declval<W&>().first(std::declval<G&>()))>;

        constexpr process_view(W walker, G& g) : walker_(std::move(walker)), g_(&g) {}

        class iterator
        {
        public:
            using value_type = typename process_view::value_type;
            using difference_type = std::ptrdiff_t;
            using iterator_concept = std::input_iterator_tag;

            constexpr iterator() = default;
            constexpr explicit iterator(process_view* v) : v_(v) {}

            constexpr value_type const& operator*() const { return v_->current(); }
            constexpr iterator& operator++()
            {
                v_->advance();
                return *this;
            }
            constexpr void operator++(int) { ++*this; }

            friend constexpr bool operator==(iterator const&, std::unreachable_sentinel_t) noexcept { return false; }

        private:
            process_view* v_ = nullptr;
        };

        constexpr iterator begin() { return iterator(this); }
        constexpr std::unreachable_sentinel_t end() const noexcept { return {}; }

    private:
        // L'état d'indice courant, tiré au premier accès.
        constexpr value_type const& current()
        {
            if (!started_)
            {
                state_ = walker_.first(*g_);
                started_ = true;
                pending_ = 0;
            }
            while (pending_ > 0)
            {
                state_ = walker_.next(*g_);
                --pending_;
            }
            return state_;
        }
        // Passer à l'indice suivant : rien n'est tiré avant la lecture.
        constexpr void advance()
        {
            if (!started_)
            {
                state_ = walker_.first(*g_);                                // l'état suivant dépend de lui
                started_ = true;
            }
            else if (pending_ > 0)
            {
                state_ = walker_.next(*g_);
                --pending_;
            }
            ++pending_;
        }

        W walker_;
        G* g_;
        value_type state_{};
        bool started_ = false;
        std::size_t pending_ = 0;
    };

    // Trajectoire paresseuse de P, tirée sur g.
    template <class P, std::uniform_random_bit_generator G>
    constexpr auto trajectory(P const& process, G& g)
    {
        return process_view<decltype(process.walker()), G>(process.walker(), g);
    }

}

// Outils internes : dans wsl::probability::detail (un seul « detail » visible depuis probability).
namespace wsl::inline WSL_DEFAULTS::probability::detail
{
        // Marcheur des tirages indépendants : chaque état est un tirage de la loi.
        template <class L>
        struct sample_walker
        {
            L law;
            template <class G>
            constexpr auto first(G& g) const
            {
                return wsl::probability::sample(law, g);
            }
            template <class G>
            constexpr auto next(G& g) const
            {
                return wsl::probability::sample(law, g);
            }
        };
}

namespace wsl::inline WSL_DEFAULTS::probability::stochastic_processes
{

    // Tirages indépendants de X, en vue paresseuse (même contrat que trajectory : une valeur
    // n'est tirée qu'à la lecture ; les n premières valeurs sont celles de n appels
    // sample(X, g) successifs).
    //
    //   for (double x : wsl::probability::samples(X, g) | std::views::take(1000)) …
    template <class L, std::uniform_random_bit_generator G>
        requires samplable<L, G>
    constexpr auto samples(L const& law, G& g)
    {
        return process_view<detail::sample_walker<L>, G>(detail::sample_walker<L>{law}, g);
    }

    template <class P>
    concept stochastic_process = requires(P const& p, std::mt19937_64& g) {
        typename P::state_type;
        { trajectory(p, g) } -> std::ranges::input_range;
    };
}
