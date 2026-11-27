#pragma once

// Arithmétique étendue des lois : celle de wsl::special (detail/extended.hpp), et ses
// compléments communs (opposé, somme, différence et complément à l'échelle, e^l arrondi…).

#include <wsl/special.hpp>   // transition (étape 4) : noms wsl::special
#include <wsl/special_functions/detail/extended.hpp>

#include <cstdint>
#include <limits>

namespace wsl::random::detail
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
    using xd::pi_squared_over_3_hi;
    using xd::sqrt_half_pi_hi;
    using xd::sqrt_half_pi_lo;
    using xd::four_minus_pi_half_hi;
    using xd::four_minus_pi_half_lo;
    using xd::sqrt_two_ln2_hi;
    using xd::sqrt_two_ln2_lo;
    using xd::gumbel_skewness_hi;
    using xd::gumbel_skewness_lo;
    using xd::rayleigh_skewness_hi;
    using xd::rayleigh_skewness_lo;
    using xd::rayleigh_kurtosis_hi;
    using xd::rayleigh_kurtosis_lo;
    using xd::pi_squared_over_3_lo;
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

    // -----------------------------------------------------------------------
    // Compléments communs aux lois (double-double, nombres à l'échelle)
    // -----------------------------------------------------------------------

    // −v.
    constexpr sd::double_double negated(sd::double_double v) noexcept { return {-v.hi, -v.lo}; }

    // v d, d double.
    constexpr sd::double_double times_double(sd::double_double v, double d) noexcept { return mul(v, {d, 0.0}); }

    // log v, v > 0 en double-double.
    constexpr sd::double_double log_positive(sd::double_double v) noexcept { return log(normalized(v.hi, v.lo)); }

    // e^l arrondi une fois (0 sous e^-1500, +∞ au-delà de e^1400, NaN conservé).
    constexpr double rounded_exp(sd::double_double l) noexcept
    {
        if (!(l.hi > -1500))
            return l.hi == l.hi ? 0.0 : l.hi;
        if (l.hi > 1400)
            return std::numeric_limits<double>::infinity();
        return rounded(exp_scaled(l.hi, l.lo));
    }

    // a + b, a, b >= 0 à l'échelle (zéro : hi = 0).
    constexpr scaled add_scaled(scaled const& a, scaled const& b) noexcept
    {
        if (a.hi == 0)
            return b;
        if (b.hi == 0)
            return a;
        scaled const& big = a.e >= b.e ? a : b;
        scaled const& small = a.e >= b.e ? b : a;
        if (big.e - small.e > 120)
            return big;
        int const shift = small.e - big.e;
        auto const t = sum({big.hi, big.lo}, {sd::scale_by_power_of_two(small.hi, shift), sd::scale_by_power_of_two(small.lo, shift)});
        return normalized(t.hi, t.lo, big.e);
    }

    // a − b, a > b >= 0 à l'échelle.
    constexpr scaled scaled_difference(scaled const& a, scaled const& b) noexcept
    {
        if (b.hi == 0 || a.e - b.e > 120)
            return a;
        int const shift = b.e - a.e;
        auto const s = sum({a.hi, a.lo}, {-sd::scale_by_power_of_two(b.hi, shift), -sd::scale_by_power_of_two(b.lo, shift)});
        if (!(s.hi > 0))
            return {0.0, 0.0, 0};
        return normalized(s.hi, s.lo, a.e);
    }

    // 1 − v, 0 <= v <= 1 à l'échelle (v <= 0,7 en pratique : résultat normal).
    constexpr scaled scaled_complement(scaled const& v) noexcept
    {
        if (v.hi == 0)
            return {1.0, 0.0, 0};
        auto const [w_hi, w_lo] = unscaled(v);
        auto const [s, s_error] = sd::two_sum(1.0, -w_hi);
        auto const [hi, lo] = sd::fast_two_sum(s, s_error - w_lo);
        if (!(hi > 0))
            return {0.0, 0.0, 0};
        return normalized(hi, lo);
    }

    // ⌊x⌋, |x| < 2^52 (exact)
    constexpr double floor_small(double x) noexcept
    {
        double const t = static_cast<double>(static_cast<std::int64_t>(x));
        return t > x ? t - 1 : t;
    }
}
