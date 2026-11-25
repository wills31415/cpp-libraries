#pragma once

// Marqueur commun des expressions paresseuses (docs/mdarray.md, V8, L2).

#include <concepts>
#include <type_traits>

namespace wsl::detail::mdarray
{
    struct expression_tag
    {
    };

    template <class T>
    concept expression = std::derived_from<std::remove_cvref_t<T>, expression_tag>;

    // Expressions indicées (monde des axes, docs/mdarray.md V3, V4).
    struct indexed_tag
    {
    };

    template <class T>
    concept indexed = std::derived_from<std::remove_cvref_t<T>, indexed_tag>;

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
