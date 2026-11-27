#pragma once

// Marqueur commun des expressions paresseuses (design/mdarray.md, V8, L2).

#include <concepts>
#include <cstddef>
#include <type_traits>

namespace wsl::detail::mdarray
{
    struct expression_tag
    {
    };

    template <class T>
    concept expression = std::derived_from<std::remove_cvref_t<T>, expression_tag>;

    // Paquet de W valeurs (voies) : W positions calculées côte à côte, chacune par les
    // mêmes opérations que seule (indexed.hpp, expression.hpp).
    template <class T, std::size_t W>
    struct lane_pack
    {
        T v[W];
    };

    // Expressions indicées (monde des axes, design/mdarray.md V3, V4).
    struct indexed_tag
    {
    };

    template <class T>
    concept indexed = std::derived_from<std::remove_cvref_t<T>, indexed_tag>;

    // +, −, × élément par élément. En mode rapide, leur fusion en fma est permise
    // (I9) : marquée explicitement, puisqu'elles viennent de fonctions distinctes que
    // le compilateur ne fusionnerait pas d'elles-mêmes. En mode reproductible, le produit
    // passe par une barrière (mode.hpp) : aucune fusion possible.
    struct add_op
    {
        template <class A, class B>
        constexpr auto operator()(A const& a, B const& b) const -> decltype(a + b)
        {
#if defined(__clang__)
#    pragma clang fp contract(fast)
#endif
            return a + b;
        }
    };
    struct subtract_op
    {
        template <class A, class B>
        constexpr auto operator()(A const& a, B const& b) const -> decltype(a - b)
        {
#if defined(__clang__)
#    pragma clang fp contract(fast)
#endif
            return a - b;
        }
    };
    struct multiply_op
    {
        template <class A, class B>
        constexpr auto operator()(A const& a, B const& b) const -> decltype(a * b)
        {
#if defined(__clang__)
#    pragma clang fp contract(fast)
#endif
            return a * b;
        }
    };

    // Modes de calcul (mode.hpp) : réécriture d'un arbre, mode par défaut d'une
    // instruction non marquée. Inclure <wsl/mdarray.hpp>.
    template <class Mode, class T>
    constexpr auto with_mode(T&& t);
    template <class Mode, class Op>
    constexpr auto with_mode_op(Op op);
    template <class T>
    constexpr decltype(auto) default_mode(T&& t);

    // m[i, j] avec au moins un axe : défini dans indexed.hpp.
    template <class Source, class... Args>
    constexpr auto make_indexed(Source&& source, Args const&... args);
}
