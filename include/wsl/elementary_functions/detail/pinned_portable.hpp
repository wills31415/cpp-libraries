#pragma once

// Transition (refonte, passe 4b) : les fonctions élémentaires, figées en mode portable, pour
// les appelants pas encore convertis aux modes (special_functions, <wsl/special>, le mode
// portable de mdarray). À retirer à la fin de la passe 4b.

#include <wsl/elementary_functions/algebraic.hpp>
#include <wsl/elementary_functions/exponential.hpp>
#include <wsl/elementary_functions/trigonometric.hpp>
#include <wsl/modes.hpp>

namespace wsl::inline WSL_DEFAULTS::elementary_functions::detail::pinned_portable
{
#define WSL_EF_PINNED(name)                                                                                    \
    template <class... A>                                                                                      \
    constexpr auto name(A... a) noexcept -> decltype(elementary_functions::name(a..., wsl::portable))         \
    {                                                                                                          \
        return elementary_functions::name(a..., wsl::portable);                                                \
    }
    WSL_EF_PINNED(exp)
    WSL_EF_PINNED(exp2)
    WSL_EF_PINNED(expm1)
    WSL_EF_PINNED(log)
    WSL_EF_PINNED(log2)
    WSL_EF_PINNED(log10)
    WSL_EF_PINNED(log1p)
    WSL_EF_PINNED(pow)
    WSL_EF_PINNED(sqrt)
    WSL_EF_PINNED(cbrt)
    WSL_EF_PINNED(hypot)
    WSL_EF_PINNED(sin)
    WSL_EF_PINNED(cos)
    WSL_EF_PINNED(tan)
    WSL_EF_PINNED(atan)
    WSL_EF_PINNED(atan2)
#undef WSL_EF_PINNED

}
