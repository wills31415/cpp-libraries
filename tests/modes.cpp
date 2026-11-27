// Modes de calcul (design/interfaces.md, décisions 24 à 27) : étiquettes, défauts, signature
// de l'utilisateur.
#include <wsl/modes.hpp>

#include <cassert>
#include <cstdio>
#include <type_traits>

namespace
{
    struct some_function
    {
    };
}

// Une bibliothèque déclare une vraie variante par spécialisation.
template <>
inline constexpr bool wsl::provides<some_function, wsl::approximate_t> = true;

int main()
{
    // Les deux axes, leurs concepts, exigences et permissions.
    static_assert(wsl::reproducibility_mode<wsl::deterministic_t> && wsl::reproducibility_mode<wsl::portable_t>);
    static_assert(wsl::algorithm_mode<wsl::approximate_t> && wsl::algorithm_mode<wsl::exact_t>);
    static_assert(!wsl::reproducibility_mode<wsl::exact_t> && !wsl::algorithm_mode<wsl::portable_t>);
    static_assert(!wsl::mode<int>);
    static_assert(wsl::is_requirement<wsl::portable_t> && wsl::is_requirement<wsl::exact_t>);
    static_assert(!wsl::is_requirement<wsl::deterministic_t> && !wsl::is_requirement<wsl::approximate_t>);

    // Défauts : la permission, sauf option globale (WSL_PORTABLE, WSL_EXACT).
#if WSL_PORTABLE
    static_assert(std::is_same_v<wsl::default_reproducibility, wsl::portable_t>);
#else
    static_assert(std::is_same_v<wsl::default_reproducibility, wsl::deterministic_t>);
    // L'espace de noms en ligne porte le défaut.
    static_assert(std::is_same_v<wsl::deterministic_approximate::portable_t, wsl::portable_t>);
#endif
#if WSL_EXACT
    static_assert(std::is_same_v<wsl::default_algorithm, wsl::exact_t>);
#else
    static_assert(std::is_same_v<wsl::default_algorithm, wsl::approximate_t>);
#endif

    // wsl::provides : faux sauf déclaration.
    static_assert(wsl::provides<some_function, wsl::approximate_t>);
    static_assert(!wsl::provides<some_function, wsl::deterministic_t>);

    // Signature de l'utilisateur : la fonction reste appelable, et se reconnaît.
    auto f = wsl::assume_portable([](double x, double y) { return x + y; });
    static_assert(wsl::is_assumed_portable<decltype(f)>);
    static_assert(!wsl::is_assumed_portable<some_function>);
    assert(f(1.0, 2.0) == 3.0);
    constexpr auto g = wsl::assume_portable([](int x) { return 2 * x; });
    static_assert(g(21) == 42);

    std::puts("modes : OK");
}
