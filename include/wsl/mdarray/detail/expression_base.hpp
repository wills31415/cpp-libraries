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
}
