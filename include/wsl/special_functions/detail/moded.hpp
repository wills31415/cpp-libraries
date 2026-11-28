#pragma once

// Interfaces publiques des fonctions spéciales dans le mode des étiquettes (design/
// interfaces.md, décisions 24 à 27, 38) : déterministe par défaut, wsl::portable sur demande ;
// wsl::approximate accepté ; wsl::exact refusé (aucune n'est correctement arrondie).
//
//   WSL_SF_FUNCTION(nom, conversion_float, arité)   nom(double…, M...), nom(float…, M...),
//                                                    arguments entiers ou mêlés en double
//
// L'implémentation est detail::nom<R>(…), R le mode de reproductibilité.

#include <wsl/modes.hpp>

#include <concepts>
#include <type_traits>

namespace wsl::inline WSL_DEFAULTS::special_functions::detail
{
    // Arguments acceptés par la surcharge en double : arithmétiques, ni tous float (leur
    // surcharge), ni long double.
    template <class... A>
    concept double_arguments = (std::is_arithmetic_v<A> && ...) && (!std::is_same_v<A, long double> && ...)
                            && !(std::is_same_v<A, float> && ...) && !(std::is_same_v<A, double> && ...);
}

// Refus d'exact : demandé par une étiquette, ou par le défaut du projet (WSL_EXACT, mode strict
// provisoire de la décision 41).
#define WSL_SF_REFUSE_EXACT(name)                                                                              \
    static_assert(!std::is_same_v<algorithm_of<M...>, exact_t>, "wsl::special_functions::" #name               \
                  " : pas de version exacte (design/interfaces.md, décisions 25 et 41)")

#define WSL_SF_FUNCTION_1(name, to_single)                                                                     \
    template <class... M>                                                                                      \
        requires mode_set<M...>                                                                                \
    constexpr double name(double x, M...) noexcept                                                             \
    {                                                                                                          \
        WSL_SF_REFUSE_EXACT(name);                                                                             \
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

#define WSL_SF_FUNCTION_2(name, to_single)                                                                     \
    template <class... M>                                                                                      \
        requires mode_set<M...>                                                                                \
    constexpr double name(double a, double x, M...) noexcept                                                   \
    {                                                                                                          \
        WSL_SF_REFUSE_EXACT(name);                                                                             \
        return detail::name<reproducibility_of<M...>>(a, x);                                                   \
    }                                                                                                          \
    template <class... M>                                                                                      \
        requires mode_set<M...>                                                                                \
    constexpr float name(float a, float x, M... m) noexcept                                                    \
    {                                                                                                          \
        return to_single(name(static_cast<double>(a), static_cast<double>(x), m...));                          \
    }                                                                                                          \
    template <class A, class B, class... M>                                                                    \
        requires detail::double_arguments<A, B> && mode_set<M...>                                              \
    constexpr double name(A a, B x, M... m) noexcept                                                           \
    {                                                                                                          \
        return name(static_cast<double>(a), static_cast<double>(x), m...);                                     \
    }

#define WSL_SF_FUNCTION_3(name, to_single)                                                                     \
    template <class... M>                                                                                      \
        requires mode_set<M...>                                                                                \
    constexpr double name(double a, double b, double x, M...) noexcept                                         \
    {                                                                                                          \
        WSL_SF_REFUSE_EXACT(name);                                                                             \
        return detail::name<reproducibility_of<M...>>(a, b, x);                                                \
    }                                                                                                          \
    template <class... M>                                                                                      \
        requires mode_set<M...>                                                                                \
    constexpr float name(float a, float b, float x, M... m) noexcept                                           \
    {                                                                                                          \
        return to_single(name(static_cast<double>(a), static_cast<double>(b), static_cast<double>(x), m...));  \
    }                                                                                                          \
    template <class A, class B, class C, class... M>                                                           \
        requires detail::double_arguments<A, B, C> && mode_set<M...>                                           \
    constexpr double name(A a, B b, C x, M... m) noexcept                                                      \
    {                                                                                                          \
        return name(static_cast<double>(a), static_cast<double>(b), static_cast<double>(x), m...);             \
    }
