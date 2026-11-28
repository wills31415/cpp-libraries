#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Décomposition de Loader : celle de wsl::special (detail/loader.hpp).

#include <wsl/random/detail/portable_math.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/special_functions/detail/loader.hpp>

namespace wsl::inline WSL_DEFAULTS::probability::detail
{
    namespace sd = random::detail::sd;
    namespace xd = random::detail::xd;

    using xd::deviance;
    using xd::digamma_approximate;
    using xd::log_beta;
    using xd::log_gamma;
    using xd::stirling_error;
    using xd::stirling_error_derivative;
}
