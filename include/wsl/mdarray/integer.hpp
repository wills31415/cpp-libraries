#pragma once

// Arithmétique entière générale (docs/mdarray.md, D15) : division euclidienne,
// reste dans {0, …, k − 1} (noms de Rust). `/` et `%` de C++ tronquent vers zéro :
// (−1) / 4 == 0 et (−1) % 4 == −1, alors que div_euclid(−1, 4) == −1 et
// rem_euclid(−1, 4) == 3 (bord périodique, tore).

#include <wsl/mdarray/detail/precondition.hpp>

#include <concepts>
#include <type_traits>

namespace wsl
{
    namespace detail::mdarray
    {
        // Entiers acceptés comme indices : ni bool ni caractères (comme std::extents).
        template <class T>
        concept index_integer = std::integral<T>
            && !std::same_as<std::remove_cv_t<T>, bool>
            && !std::same_as<std::remove_cv_t<T>, char>
            && !std::same_as<std::remove_cv_t<T>, wchar_t>
            && !std::same_as<std::remove_cv_t<T>, char8_t>
            && !std::same_as<std::remove_cv_t<T>, char16_t>
            && !std::same_as<std::remove_cv_t<T>, char32_t>;
    }

    // Reste euclidien de a par k > 0 : dans {0, …, k − 1}.
    template <detail::mdarray::index_integer T>
    constexpr T rem_euclid(T a, T k) noexcept
    {
        detail::mdarray::precondition(k > 0, "rem_euclid : diviseur non positif");
        T const r = static_cast<T>(a % k);
        if constexpr (std::is_signed_v<T>)
            return r < 0 ? static_cast<T>(r + k) : r;
        else
            return r;
    }

    // Quotient euclidien de a par k > 0 : a == k * div_euclid(a, k) + rem_euclid(a, k).
    template <detail::mdarray::index_integer T>
    constexpr T div_euclid(T a, T k) noexcept
    {
        detail::mdarray::precondition(k > 0, "div_euclid : diviseur non positif");
        T const q = static_cast<T>(a / k);
        if constexpr (std::is_signed_v<T>)
            return a % k < 0 ? static_cast<T>(q - 1) : q;
        else
            return q;
    }
}
