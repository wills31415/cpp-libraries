#pragma once

// Décomposition de Loader : celle de wsl::special (detail/loader.hpp).

#include <wsl/special.hpp>   // transition (étape 4) : noms wsl::special
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/special_functions/detail/loader.hpp>

namespace wsl::random::detail
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
