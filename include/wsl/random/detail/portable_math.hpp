#pragma once

// Le calcul flottant des tirages et des lois, toujours portable (design/interfaces.md, décision
// 42 : les lois n'ont pour l'instant que leur version portable) :
//   sd  primitives de floating_point_arithmetic et noyaux de elementary_functions et
//       special_functions (sd::two_sum, sd::fmul, sd::scaled…) ;
//   xd  le socle double-double de special_functions (detail::extended) ;
//   pm  les fonctions publiques, figées en mode portable et approché (pm::exp, pm::lgamma…).

#include <wsl/elementary_functions.hpp>
#include <wsl/floating_point_arithmetic.hpp>
#include <wsl/modes.hpp>
#include <wsl/special_functions.hpp>
#include <wsl/special_functions/detail/extended.hpp>
#include <wsl/special_functions/detail/loader.hpp>
#include <wsl/special_functions/detail/pinned_portable.hpp>

namespace wsl::inline WSL_DEFAULTS::random::detail
{
    namespace sd
    {
        using namespace floating_point_arithmetic::detail;
        using namespace elementary_functions::detail;
        using namespace special_functions::detail;
    }
    namespace xd = special_functions::detail::extended;
    namespace pm = special_functions::detail::pinned_portable;
}
