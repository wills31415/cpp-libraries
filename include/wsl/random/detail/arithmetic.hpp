#pragma once

// Arithmétique étendue des lois : celle de wsl::special (detail/extended.hpp).

#include <wsl/special/detail/extended.hpp>

namespace wsl::detail
{
    namespace sd = wsl::special::detail;
    namespace xd = wsl::special::detail::extended;

    using xd::balanced_product;
    using xd::ceiling;
    using xd::divide;
    using xd::exact_split;
    using xd::exp_scaled;
    using xd::expm1_dd;
    using xd::expm1_double_double;
    using xd::floor_nonnegative;
    using xd::half_log_2pi_hi;
    using xd::half_log_2pi_lo;
    using xd::half_sqrt_pi;
    using xd::interpolate;
    using xd::inv_sqrt2_hi;
    using xd::inv_sqrt2_lo;
    using xd::inv_pi_hi;
    using xd::inv_pi_lo;
    using xd::inv_sqrt_2pi_hi;
    using xd::log_4pi_hi;
    using xd::log_4pi_lo;
    using xd::inv_sqrt_2pi_lo;
    using xd::inv_sqrt_pi;
    using xd::is_finite;
    using xd::log;
    using xd::log1p_dd;
    using xd::log1p_positive;
    using xd::log_dd;
    using xd::log_precise;
    using xd::magnitude;
    using xd::max_double;
    using xd::mul;
    using xd::next_below;
    using xd::normal_entropy_hi;
    using xd::normal_entropy_lo;
    using xd::normalized;
    using xd::one_twelfth_hi;
    using xd::one_twelfth_lo;
    using xd::product;
    using xd::ratio;
    using xd::reciprocal;
    using xd::rounded;
    using xd::scaled;
    using xd::sqrt;
    using xd::sqrt2_hi;
    using xd::sqrt2_lo;
    using xd::sum;
    using xd::times;
    using xd::unscaled;
}
