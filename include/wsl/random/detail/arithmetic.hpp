#pragma once

// Arithmétique soignée commune aux lois : calculs en double-double arrondis une
// seule fois, mis à l'échelle pour rester dans le domaine où les produits exacts
// le sont (facteurs sous 2^995, produits loin des sous-normaux : voir
// wsl/special/detail/fp.hpp). Bit-exact comme le reste.

#include <wsl/special.hpp>

#include <bit>
#include <cstdint>
#include <limits>

namespace wsl::detail
{
    namespace sd = wsl::special::detail;

    inline constexpr double max_double = std::numeric_limits<double>::max();

    constexpr bool is_finite(double x) noexcept { return x - x == 0; }

    // Le flottant immédiatement inférieur à x (x fini, > plus petit double).
    constexpr double next_below(double x) noexcept
    {
        if (x == 0)
            return -0x1p-1074;
        std::uint64_t const b = std::bit_cast<std::uint64_t>(x);
        return std::bit_cast<double>(x > 0 ? b - 1 : b + 1);
    }

    constexpr float next_below(float x) noexcept
    {
        if (x == 0)
            return -0x1p-149f;
        std::uint32_t const b = std::bit_cast<std::uint32_t>(x);
        return std::bit_cast<float>(x > 0 ? b - 1 : b + 1);
    }

    // (n_hi + n_lo) / (d_hi + d_lo), d > 0, 0 <= n, arrondi une fois (à 2^-100 près
    // avant l'arrondi). Mise à l'échelle : d dans [1, 2), n remonté s'il est
    // minuscule, pour que le produit exact q·d soit valable.
    constexpr double divide(double n_hi, double n_lo, double d_hi, double d_lo) noexcept
    {
        if (n_hi == 0)
            return 0.0;
        int const e = sd::exponent_of(d_hi);
        double const dh = sd::scale_by_power_of_two(d_hi, -e);
        double const dl = sd::scale_by_power_of_two(d_lo, -e);
        int shift = -e;
        if (sd::exponent_of(n_hi) - e < -800)                               // quotient minuscule : remonté
            shift += 400;
        double const nh = sd::scale_by_power_of_two(n_hi, shift);
        double const nl = sd::scale_by_power_of_two(n_lo, shift);
        double const q = nh / dh;
        auto const [p, p_error] = sd::two_product(q, dh);
        double const r = q + (((nh - p) - p_error) + (nl - sd::fmul(q, dl))) / dh;
        return sd::scale_by_power_of_two(r, -(shift + e));
    }

    // (x − a) / (b − a), a <= x <= b, a < b finis : arrondi une fois. Les
    // différences, qui peuvent déborder, sont prises sur les moitiés (exactes).
    constexpr double ratio(double x, double a, double b) noexcept
    {
        if (!is_finite(b - a))
        {
            x = sd::fmul(x, 0.5);
            a = sd::fmul(a, 0.5);
            b = sd::fmul(b, 0.5);
        }
        auto const [n, n_error] = sd::two_sum(x, -a);
        auto const [d, d_error] = sd::two_sum(b, -a);
        return divide(n, n_error, d, d_error);
    }

    // a + p (b − a), p dans [0, 1], a < b finis : arrondi une fois.
    constexpr double interpolate(double a, double b, double p) noexcept
    {
        double factor = 1.0;
        if (!is_finite(b - a))
        {
            a = sd::fmul(a, 0.5);
            b = sd::fmul(b, 0.5);
            factor = 2.0;
        }
        auto const [d, d_error] = sd::two_sum(b, -a);                        // b − a exactement
        int const e = sd::exponent_of(d);
        double const dn = sd::scale_by_power_of_two(d, -e);                  // dans [1, 2)
        double const dl = sd::scale_by_power_of_two(d_error, -e);
        int shift = e;
        double pp = p;
        if (p != 0 && p < 0x1p-800)                                          // produit loin des sous-normaux
        {
            pp = sd::fmul(p, 0x1p400);
            shift -= 400;
        }
        auto const [m, m_error] = sd::two_product(pp, dn);
        double const tail = m_error + sd::fmul(pp, dl);
        auto const [s, s_error] = sd::two_sum(a, sd::scale_by_power_of_two(m, shift));
        return sd::fmul(s + (s_error + sd::scale_by_power_of_two(tail, shift)), factor);
    }

    // log(hi + lo), hi > 0 : log(hi) + lo/hi, à moins d'un ulp.
    constexpr double log_dd(double hi, double lo) noexcept
    {
        return wsl::special::log(hi) + lo / hi;
    }
}
