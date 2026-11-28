#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Trigonométrie en double-double des lois : π, π/2, sin u (u dans [0, π]), atan v.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/portable_math.hpp>

namespace wsl::inline WSL_DEFAULTS::probability::detail
{
    inline constexpr double pi_hi = 0x1.921fb54442d18p+1, pi_lo = 0x1.1a62633145c07p-53;

    inline constexpr sd::double_double half_pi_dd{0x1.921fb54442d18p+0, 0x1.1a62633145c07p-54};

    // sin(h), |h| <= π/2 + ε, h en double-double (à 2^-70 près).
    constexpr sd::double_double sine_dd(double h_hi, double h_lo) noexcept
    {
        namespace fp = sd;
        double const a = magnitude(h_hi);
        if (a == 0)
            return {h_hi, h_lo};
        if (a <= 0x1.921fb54442d18p-1)                                  // |h| <= π/4
            return fp::sin_of(fp::sincos_prepare(h_hi, h_lo));
        bool const negative = h_hi < 0;                                 // sin|h| = cos(π/2 − |h|)
        auto const r = sum({fp::fmul(pi_hi, 0.5), fp::fmul(pi_lo, 0.5)}, {-a, negative ? h_lo : -h_lo});
        auto const c = fp::cos_of(fp::sincos_prepare(r.hi, r.lo));
        return negative ? fp::double_double{-c.hi, -c.lo} : c;
    }

    // sin u, 0 <= u <= π (+ε), en double-double.
    constexpr sd::double_double sin_any(sd::double_double u) noexcept
    {
        if (u.hi <= half_pi_dd.hi)
            return sine_dd(u.hi, u.lo);
        auto const r = sum({pi_hi, pi_lo}, negated(u));
        return sine_dd(r.hi, r.lo);
    }

    // atan v, en double-double.
    constexpr sd::double_double atan_dd(sd::double_double v) noexcept
    {
        bool const negative = v.hi < 0;
        if (negative)
            v = negated(v);
        sd::double_double r{};
        if (v.hi < 0x1p-27)
            r = {v.hi, v.lo - sd::fmul(sd::fmul(v.hi, v.hi), v.hi) / 3};
        else if (v.hi <= 1)
            r = sd::atan_kernel(v.hi, v.lo);
        else
        {
            auto const w = sd::quotient_dd({1.0, 0.0}, v);
            sd::double_double const a = w.hi < 0x1p-27 ? sd::double_double{w.hi, w.lo - sd::fmul(sd::fmul(w.hi, w.hi), w.hi) / 3} : sd::atan_kernel(w.hi, w.lo);
            r = sum(half_pi_dd, negated(a));
        }
        return negative ? negated(r) : r;
    }
}
