#pragma once

// Quantiles des lois discrètes : plus petit k tel que F(k) >= p, F connue par sa
// plus petite part (F ou 1 − F, en double-double à exposant séparé, à 2^-64 près),
// comparée à p exactement — le quantile est celui de la vraie loi, sauf si F(k)
// tombe à 2^-60 près (relativement) de p.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/special.hpp>

#include <limits>

namespace wsl::detail
{
    // signe de v − d, v >= 0 à l'échelle, d >= 0 en double-double.
    constexpr int compare(scaled const& v, sd::double_double d) noexcept
    {
        if (d.hi == 0)
            return v.hi == 0 ? 0 : 1;
        if (v.hi == 0)
            return -1;
        scaled const r = v / normalized(d.hi, d.lo);                        // v/d = (hi + lo) 2^e
        if (r.e != 0)
            return r.e < 0 ? -1 : 1;
        if (r.hi != 1.0)
            return r.hi < 1.0 ? -1 : 1;
        return r.lo < 0 ? -1 : r.lo > 0 ? 1 : 0;
    }

    // F >= p, value = F si is_f, sinon value = 1 − F.
    constexpr bool at_least(scaled const& value, bool is_f, double p) noexcept
    {
        if (is_f)
            return compare(value, {p, 0.0}) >= 0;
        auto const [q_hi, q_lo] = sd::two_sum(1.0, -p);                     // 1 − p, exact
        return compare(value, {q_hi, q_lo}) <= 0;                           // 1 − F <= 1 − p
    }

    // Plus petit k de [low, high] tel que ge(k), ge croissante et ge(high) vrai :
    // recherche exponentielle depuis guess, puis dichotomie.
    template <class I, class Ge>
    constexpr I smallest_true(I guess, I low, I high, Ge ge)
    {
        guess = guess < low ? low : guess > high ? high : guess;
        I below = low, above = high;                                         // ge(above) ; ge faux sous below
        bool below_known_false = false;
        if (ge(guess))
        {
            above = guess;
            I step = 1;
            while (above > low)
            {
                I const probe = above - low > step ? static_cast<I>(above - step) : low;
                if (!ge(probe))
                {
                    below = probe;
                    below_known_false = true;
                    break;
                }
                above = probe;
                if (step > std::numeric_limits<I>::max() / 2)
                    break;
                step = static_cast<I>(step * 2);
            }
            if (!below_known_false)
                return above;                                               // above == low
        }
        else
        {
            below = guess;
            below_known_false = true;
            I step = 1;
            for (;;)
            {
                if (high - below <= step)
                {
                    above = high;
                    break;
                }
                I const probe = static_cast<I>(below + step);
                if (ge(probe))
                {
                    above = probe;
                    break;
                }
                below = probe;
                if (step > std::numeric_limits<I>::max() / 2)
                {
                    above = high;
                    break;
                }
                step = static_cast<I>(step * 2);
            }
        }
        while (above - below > 1)                                            // ge(below) faux, ge(above) vrai
        {
            I const middle = static_cast<I>(below + (above - below) / 2);
            if (ge(middle))
                above = middle;
            else
                below = middle;
        }
        return above;
    }

    // Φ⁻¹(p), p dans (0, 1), approché (point de départ des recherches).
    constexpr double normal_quantile_guess(double p) noexcept
    {
        return -sd::fmul(sqrt2_hi, wsl::special::erfc_inv(2.0 * p));
    }
}
