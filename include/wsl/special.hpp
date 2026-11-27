#pragma once

// Transition (refonte, étape 4) : wsl::special est découpée en
// wsl::floating_point_arithmetic, wsl::elementary_functions et wsl::special_functions
// (design/interfaces.md, décision 37). Cet en-tête réexporte leurs noms sous wsl::special
// pour random et mdarray, jusqu'à leur propre refonte (étapes 5 et 6) ; il disparaîtra
// ensuite. Tout nouveau code emploie les nouveaux domaines.

#include <wsl/elementary_functions.hpp>
#include <wsl/special_functions/detail/pinned_portable.hpp>
#include <wsl/floating_point_arithmetic.hpp>
#include <wsl/modes.hpp>
#include <wsl/special_functions.hpp>

namespace wsl::inline WSL_DEFAULTS::special
{
    using namespace special_functions::detail::pinned_portable;   // figées en mode portable

    namespace detail
    {
        using namespace floating_point_arithmetic::detail;
        using namespace elementary_functions::detail;
        using namespace special_functions::detail;
    }
}
