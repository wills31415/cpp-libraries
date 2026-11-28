#pragma once

// Interfaces publiques des fonctions élémentaires dans le mode des étiquettes (design/
// interfaces.md, décisions 24 à 27, 38) : déterministe par défaut, wsl::portable sur demande ;
// wsl::approximate accepté. wsl::exact n'est accepté que par les fonctions correctement
// arrondies (EXACT vrai) : les autres le refusent à la compilation (décision 25).
//
//   WSL_EF_UNARY(nom, conversion_float, EXACT)   nom(double, M...), nom(float, M...), nom(entier, M...)
//   WSL_EF_BINARY(nom, EXACT)                    nom(double, double, M...), nom(float, float, M...),
//                                                arguments mêlés en double
//
// L'implémentation est detail::nom<R>(…), R le mode de reproductibilité.

#include <wsl/modes.hpp>

#include <concepts>
#include <type_traits>

// Refus d'exact : demandé par une étiquette, ou par le défaut du projet (WSL_EXACT, mode strict
// provisoire de la décision 41).
#define WSL_EF_REFUSE_EXACT(name, EXACT)                                                                       \
    static_assert(EXACT || !std::is_same_v<algorithm_of<M...>, exact_t>,                                       \
                  "wsl::elementary_functions::" #name " : pas de version exacte (erreur < 0,51 ulp, pas "      \
                  "l'arrondi correct ; design/interfaces.md, décisions 25 et 41)")

#define WSL_EF_UNARY(name, to_single, EXACT)                                                                   \
    template <class... M>                                                                                      \
        requires mode_set<M...>                                                                                \
    constexpr double name(double x, M...) noexcept                                                             \
    {                                                                                                          \
        WSL_EF_REFUSE_EXACT(name, EXACT);                                                                      \
        return detail::name<reproducibility_of<M...>>(x);                                                      \
    }                                                                                                          \
    template <class... M>                                                                                      \
        requires mode_set<M...>                                                                                \
    constexpr float name(float x, M... m) noexcept                                                             \
    {                                                                                                          \
        return to_single(name(static_cast<double>(x), m...));                                                  \
    }                                                                                                          \
    template <std::integral I, class... M>                                                                     \
        requires mode_set<M...>                                                                                \
    constexpr double name(I x, M... m) noexcept                                                                \
    {                                                                                                          \
        return name(static_cast<double>(x), m...);                                                             \
    }

#define WSL_EF_BINARY(name, to_single, EXACT)                                                                  \
    template <class... M>                                                                                      \
        requires mode_set<M...>                                                                                \
    constexpr double name(double x, double y, M...) noexcept                                                   \
    {                                                                                                          \
        WSL_EF_REFUSE_EXACT(name, EXACT);                                                                      \
        return detail::name<reproducibility_of<M...>>(x, y);                                                   \
    }                                                                                                          \
    template <class... M>                                                                                      \
        requires mode_set<M...>                                                                                \
    constexpr float name(float x, float y, M... m) noexcept                                                    \
    {                                                                                                          \
        return to_single(name(static_cast<double>(x), static_cast<double>(y), m...));                          \
    }                                                                                                          \
    template <class A, class B, class... M>                                                                    \
        requires detail::mixed_arithmetic<A, B> && mode_set<M...>                                              \
    constexpr double name(A x, B y, M... m) noexcept                                                           \
    {                                                                                                          \
        return name(static_cast<double>(x), static_cast<double>(y), m...);                                     \
    }
