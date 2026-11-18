#pragma once

// Décomposition de Loader : celle de wsl::special (detail/loader.hpp).

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/special/detail/loader.hpp>

namespace wsl::detail
{
    namespace sd = wsl::special::detail;
    namespace xd = wsl::special::detail::extended;

    using xd::deviance;
    using xd::digamma_approximate;
    using xd::log_beta;
    using xd::log_gamma;
    using xd::stirling_error;
    using xd::stirling_error_derivative;
}
