#pragma once

#include <wsl/modes.hpp>

// Préconditions de wsl::mdarray (design/mdarray.md, P3) : aucune vérification en
// Release, arrêt avec message quand NDEBUG n'est pas défini, choix forcé par
// WSL_ASSERTIONS (0 ou 1) ; toujours vérifiées en évaluation constante, où une
// violation est une erreur de compilation.
//
// WSL_ASSERTIONS doit valoir la même chose dans tout le programme (comme NDEBUG
// pour assert).

#include <cstdio>
#include <cstdlib>

#if !defined(WSL_ASSERTIONS)
#    if defined(NDEBUG)
#        define WSL_ASSERTIONS 0
#    else
#        define WSL_ASSERTIONS 1
#    endif
#endif

namespace wsl::inline WSL_DEFAULTS::detail::mdarray
{
    // Non constexpr exprès : l'appeler pendant une évaluation constante la rend
    // invalide, d'où l'erreur de compilation.
    [[noreturn]] inline void precondition_failed(char const* message) noexcept
    {
        std::fprintf(stderr, "wsl::mdarray : précondition violée : %s\n", message);
        std::abort();
    }

    constexpr void precondition(bool condition, char const* message) noexcept
    {
        if consteval
        {
            if (!condition)
                precondition_failed(message);
        }
        else
        {
            if constexpr (WSL_ASSERTIONS)
            {
                if (!condition)
                    precondition_failed(message);
            }
        }
    }
}
