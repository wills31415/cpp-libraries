#pragma once

// Transition (refonte, passe 4b) : les fonctions spéciales, figées en mode portable, pour les
// appelants pas encore convertis aux modes (<wsl/special>, donc random et probability ; le
// mode portable de mdarray). À retirer avec la passe 5b.

#include <wsl/elementary_functions/detail/pinned_portable.hpp>
#include <wsl/modes.hpp>
#include <wsl/special_functions.hpp>

namespace wsl::inline WSL_DEFAULTS::special_functions::detail::pinned_portable
{
    using namespace elementary_functions::detail::pinned_portable;

#define WSL_SF_PINNED(name)                                                                                    \
    template <class... A>                                                                                      \
    constexpr auto name(A... a) noexcept -> decltype(special_functions::name(a..., wsl::portable))            \
    {                                                                                                          \
        return special_functions::name(a..., wsl::portable);                                                   \
    }
    WSL_SF_PINNED(erf)
    WSL_SF_PINNED(erfc)
    WSL_SF_PINNED(erf_inv)
    WSL_SF_PINNED(erfc_inv)
    WSL_SF_PINNED(lgamma)
    WSL_SF_PINNED(tgamma)
    WSL_SF_PINNED(digamma)
    WSL_SF_PINNED(beta)
    WSL_SF_PINNED(lbeta)
    WSL_SF_PINNED(gamma_p)
    WSL_SF_PINNED(gamma_q)
    WSL_SF_PINNED(gamma_p_inv)
    WSL_SF_PINNED(gamma_q_inv)
    WSL_SF_PINNED(ibeta)
    WSL_SF_PINNED(ibetac)
    WSL_SF_PINNED(ibeta_inv)
    WSL_SF_PINNED(ibetac_inv)
    WSL_SF_PINNED(cyl_bessel_i)
#undef WSL_SF_PINNED

}
