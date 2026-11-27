#pragma once

// Modes de calcul, les mêmes dans toutes les bibliothèques (design/interfaces.md, décisions
// 24 à 27). Deux axes indépendants, deux noms chacun ; le défaut est toujours la
// permission, les étiquettes servent à exiger davantage.
//
//   reproductibilité : deterministic (défaut : mêmes bits sur une configuration)
//                      portable      (mêmes bits sur toute configuration IEEE 754)
//   algorithme       : approximate   (défaut : erreur bornée et documentée)
//                      exact         (le résultat de l'algorithme idéal)
//
// Une étiquette s'emploie en argument ou en enveloppe :
//
//   auto z = law(g, wsl::approximate);
//   y[i] = wsl::portable(A[i, j] * x[j] + b[i]);
//
// portable et exact sont des exigences : une fonction qui ne peut pas les honorer les refuse
// à la compilation. deterministic et approximate sont des permissions, toujours acceptées :
// la fonction rend ce qu'elle a de mieux. wsl::provides<F, Mode> dit si F a une vraie
// variante pour une permission. Une fonction de l'utilisateur dans une expression portable
// doit être signée wsl::assume_portable(f) : la bibliothèque ne peut pas inspecter son corps.
//
// Défaut du projet : WSL_PORTABLE et WSL_EXACT (options CMake du même nom, propagées à tout
// consommateur ; une option globale renforce une garantie, jamais ne l'affaiblit). Le
// défaut vit dans un espace de noms en ligne qui en porte le nom : deux parties d'un
// programme compilées avec des défauts différents ont des entités distinctes, et un
// mélange devient une erreur à l'édition de liens au lieu d'une violation silencieuse de
// la règle de la définition unique.

#include <concepts>
#include <type_traits>
#include <utility>

#if !defined(WSL_PORTABLE)
#    define WSL_PORTABLE 0
#endif
#if !defined(WSL_EXACT)
#    define WSL_EXACT 0
#endif

// Espace de noms en ligne des entités qui dépendent du défaut : namespace wsl::inline
// WSL_DEFAULTS { … }.
#if WSL_PORTABLE && WSL_EXACT
#    define WSL_DEFAULTS portable_exact
#elif WSL_PORTABLE
#    define WSL_DEFAULTS portable_approximate
#elif WSL_EXACT
#    define WSL_DEFAULTS deterministic_exact
#else
#    define WSL_DEFAULTS deterministic_approximate
#endif

namespace wsl::inline WSL_DEFAULTS
{
    // Une étiquette est appelable : tag(e) met l'expression e dans ce mode, par la
    // fonction apply_mode(tag, e) que la bibliothèque de e fournit (trouvée par ADL).
    template <class Tag>
    struct mode_tag
    {
        template <class E>
        constexpr auto operator()(E&& e) const -> decltype(apply_mode(std::declval<Tag const&>(), std::forward<E>(e)))
        {
            return apply_mode(static_cast<Tag const&>(*this), std::forward<E>(e));
        }
    };

    struct deterministic_t : mode_tag<deterministic_t>
    {
    };
    struct portable_t : mode_tag<portable_t>
    {
    };
    struct approximate_t : mode_tag<approximate_t>
    {
    };
    struct exact_t : mode_tag<exact_t>
    {
    };

    inline constexpr deterministic_t deterministic{};
    inline constexpr portable_t portable{};
    inline constexpr approximate_t approximate{};
    inline constexpr exact_t exact{};

    template <class T>
    concept reproducibility_mode = std::same_as<T, deterministic_t> || std::same_as<T, portable_t>;
    template <class T>
    concept algorithm_mode = std::same_as<T, approximate_t> || std::same_as<T, exact_t>;
    template <class T>
    concept mode = reproducibility_mode<T> || algorithm_mode<T>;

    // Exigence ou permission (décision 25).
    template <mode M>
    inline constexpr bool is_requirement = std::same_as<M, portable_t> || std::same_as<M, exact_t>;

    // Défauts du projet (décision 27).
    using default_reproducibility = std::conditional_t<WSL_PORTABLE != 0, portable_t, deterministic_t>;
    using default_algorithm = std::conditional_t<WSL_EXACT != 0, exact_t, approximate_t>;

    // Mode résultant des étiquettes passées à une fonction : au plus une par axe, le défaut
    // du projet pour un axe sans étiquette (f(x), f(x, wsl::portable), f(x, wsl::exact,
    // wsl::portable)…).
    namespace detail::modes
    {
        template <class Default, template <class> class Is, class... M>
        struct select
        {
            using type = Default;
        };
        template <class Default, template <class> class Is, class First, class... Rest>
        struct select<Default, Is, First, Rest...>
        {
            using type = std::conditional_t<Is<First>::value, First, typename select<Default, Is, Rest...>::type>;
        };
        template <class T>
        struct is_reproducibility : std::bool_constant<reproducibility_mode<T>>
        {
        };
        template <class T>
        struct is_algorithm : std::bool_constant<algorithm_mode<T>>
        {
        };
    }

    // Un jeu d'étiquettes valide : au plus une par axe.
    template <class... M>
    concept mode_set = (mode<M> && ...) && (reproducibility_mode<M> + ... + 0) <= 1
                    && (algorithm_mode<M> + ... + 0) <= 1;

    template <mode... M>
        requires((reproducibility_mode<M> + ... + 0) <= 1 && (algorithm_mode<M> + ... + 0) <= 1)
    using reproducibility_of =
        typename detail::modes::select<default_reproducibility, detail::modes::is_reproducibility, M...>::type;

    template <mode... M>
        requires((reproducibility_mode<M> + ... + 0) <= 1 && (algorithm_mode<M> + ... + 0) <= 1)
    using algorithm_of = typename detail::modes::select<default_algorithm, detail::modes::is_algorithm, M...>::type;

    // Vrai si exact est demandé explicitement parmi les étiquettes M (une fonction qui ne
    // sait pas être exacte refuse alors l'appel, décision 25).
    template <mode... M>
    inline constexpr bool exact_requested = (std::same_as<M, exact_t> || ...);

    // wsl::provides<F, M> : F a une vraie variante pour le mode M (une permission sans
    // variante propre est acceptée, mais rend le résultat le plus fort). Chaque
    // bibliothèque spécialise ce trait pour ses fonctions.
    template <class F, mode M>
    inline constexpr bool provides = false;

    // Signature de l'utilisateur : f est portable (aucun produit fusionnable, aucune
    // fonction de la libm) ; la garantie repose sur lui (comme std::assume_aligned).
    template <class F>
    struct assume_portable_fn
    {
        F function;

        template <class... A>
        constexpr decltype(auto) operator()(A&&... a) const
            noexcept(noexcept(function(std::forward<A>(a)...)))
        {
            return function(std::forward<A>(a)...);
        }
    };

    template <class F>
    constexpr assume_portable_fn<std::decay_t<F>> assume_portable(F&& f)
    {
        return {std::forward<F>(f)};
    }

    template <class T>
    inline constexpr bool is_assumed_portable = false;
    template <class F>
    inline constexpr bool is_assumed_portable<assume_portable_fn<F>> = true;
}
