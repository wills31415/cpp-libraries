#pragma once

#include <wsl/modes.hpp>

// Fonctions algébriques bit-exactes : même résultat au bit près sur toute
// plateforme (voir detail/fp.hpp pour les conditions), constexpr.
//
//   sqrt     arrondi correct (celui d'IEEE 754), évaluable à la compilation
//   cbrt     erreur < 0,51 ulp
//   hypot    erreur < 0,51 ulp, sans débordement intermédiaire
//
// Cas limites ceux de <cmath> (Annexe F de C) : sqrt(±0) = ±0, sqrt(x < 0) = NaN,
// sqrt(+∞) = +∞ ; cbrt(±0) = ±0, cbrt(±∞) = ±∞ ; hypot(±∞, y) = +∞ même si y est
// NaN, hypot(x, ±0) = |x| ; NaN donne NaN.

#include <wsl/floating_point_arithmetic/detail/primitives.hpp>

#include <concepts>
#include <cstdint>
#include <type_traits>

// Couche du dessous : les primitives flottantes.
namespace wsl::inline WSL_DEFAULTS::elementary_functions::detail
{
    using namespace floating_point_arithmetic::detail;
}

namespace wsl::inline WSL_DEFAULTS::elementary_functions
{
    namespace detail
    {
        // Racine carrée correctement arrondie, en arithmétique entière : ne sert
        // qu'à l'évaluation constante, où l'instruction matérielle n'est pas
        // disponible. Même résultat qu'elle, puisque l'arrondi correct est unique.
        constexpr double sqrt_by_integers(double x) noexcept     // x > 0 fini
        {
            std::uint64_t const b = bits(x);
            int exponent = static_cast<int>(b >> 52);
            std::uint64_t mantissa = b & 0x000fffffffffffff;
            if (exponent == 0)                                          // sous-normal : normaliser
            {
                int const shift = std::countl_zero(mantissa) - 11;
                mantissa <<= shift;
                exponent = 1 - shift;
            }
            mantissa |= 0x0010000000000000;
            int e = exponent - 1075;                                    // x = m · 2^e, m sur 53 bits
            uint128 m = mantissa;
            if (e & 1)                                                  // exposant pair
            {
                m <<= 1;
                --e;
            }
            // n = m · 2^54 dans [2^106, 2^108) : sa racine entière R a 54 bits, un
            // de plus que le résultat, pour l'arrondi.
            uint128 const n = m << 54;
            uint128 root = 0, remainder = 0;
            for (int i = 106; i >= 0; i -= 2)                           // racine chiffre par chiffre
            {
                remainder = (remainder << 2) | ((n >> i) & 3);
                uint128 const trial = (root << 2) | 1;
                root <<= 1;
                if (remainder >= trial)
                {
                    remainder -= trial;
                    root |= 1;
                }
            }
            // Arrondi au plus proche sur le bit de plus. Égalité impossible : une
            // racine de double n'est jamais à mi-chemin de deux doubles.
            std::uint64_t result = static_cast<std::uint64_t>(root >> 1) + static_cast<std::uint64_t>(root & 1);
            int result_exponent = (e - 54) / 2 + 1;                     // √(n 2^(e−54)) = R 2^((e−54)/2)
            if (result >> 53)
            {
                result >>= 1;
                ++result_exponent;
            }
            return fmul(static_cast<double>(result), power_of_two(result_exponent));
        }
    }

    constexpr double sqrt(double x) noexcept
    {
        if (!std::is_constant_evaluated())
            return __builtin_sqrt(x);                                   // IEEE 754 : arrondi correct
        if (!(x == x) || x == 0 || x == detail::infinity)
            return x;
        if (x < 0)
            return detail::quiet_nan;
        return detail::sqrt_by_integers(x);
    }

    // -----------------------------------------------------------------------
    // cbrt
    //
    // x = 2^(3q + s) z avec z dans [1, 8) : cbrt(x) = 2^q cbrt(z). Approximation
    // initiale par les bits (5 bits justes), deux itérations de Halley
    // (y ← y (y³ + 2z)/(2y³ + z), convergence cubique : 2^-45), puis un pas de
    // Newton à résidu exact (z − y³ en double-double) : 2^-89, un seul arrondi.
    // -----------------------------------------------------------------------
    constexpr double cbrt(double x) noexcept
    {
        using namespace detail;
        if (!(x == x) || x == 0 || x == infinity || x == -infinity)
            return x;
        bool const negative = x < 0;
        double ax = negative ? -x : x;
        int adjust = 0;
        if (ax < 0x1p-1022)                                             // sous-normal : × 2^54, exact
        {
            ax = fmul(ax, 0x1p54);
            adjust = -18;
        }

        std::uint64_t const b = bits(ax);
        int const exponent = static_cast<int>(b >> 52) - 1023;
        int const q = (exponent >= 0 ? exponent : exponent - 2) / 3;    // division par défaut
        int const s = exponent - 3 * q;                                 // dans {0, 1, 2}
        double const z = from_bits((b & 0x000fffffffffffff) | (static_cast<std::uint64_t>(1023 + s) << 52));

        double y = from_bits(bits(z) / 3 + 0x2a9f789300000000);
        for (int i = 0; i < 2; ++i)                                     // Halley
        {
            double const y3 = fmul(y * y, y);
            y = y * (y3 + fmul(2.0, z)) / (fmul(2.0, y3) + z);
        }

        // Pas de Newton exact : y ← y + (z − y³)/(3y²).
        auto const [a, a_error] = two_product(y, y);                    // y² exactement
        auto const [c, c_error] = two_product(a, y);
        double const residual = ((z - c) - c_error) - fmul(a_error, y); // z − y³, à 2^-106 près
        y = y + residual / fmul(3.0, a);

        double const result = fmul(y, power_of_two(q + adjust));        // exact
        return negative ? -result : result;
    }

    // -----------------------------------------------------------------------
    // hypot
    //
    // √(a² + b²), a = max(|x|, |y|), b = min : mise à l'échelle par une puissance
    // de 2 (exacte) quand a est très grand ou très petit, a² + b² en
    // double-double, puis √ corrigée par son résidu exact : un seul arrondi.
    // -----------------------------------------------------------------------
    constexpr double hypot(double x, double y) noexcept
    {
        using namespace detail;
        double ax = x < 0 ? -x : x;
        double ay = y < 0 ? -y : y;
        if (ax == infinity || ay == infinity)
            return infinity;                                            // même si l'autre est NaN
        if (!(ax == ax) || !(ay == ay))
            return x + y;
        double a = ax > ay ? ax : ay;
        double b = ax > ay ? ay : ax;
        if (b == 0)
            return a;
        if (b < fmul(a, 0x1p-54))                                       // b²/a² < 2^-108 : √ arrondie à a
            return a;

        // a ramené dans [1, 2) par une puissance de 2 (exact ; b >= a 2^-54 reste
        // normal) : ni débordement, ni carré sous-normal dans les produits exacts.
        int const exponent = exponent_of(a);
        a = scale_by_power_of_two(a, -exponent);
        b = scale_by_power_of_two(b, -exponent);

        auto const [a2, a2_error] = two_product(a, a);
        auto const [b2, b2_error] = two_product(b, b);
        auto const [s, s_error] = two_sum(a2, b2);
        double const s_lo = (s_error + a2_error) + b2_error;

        double const h = sqrt(s);
        auto const [hh, hh_error] = two_product(h, h);
        double const residual = ((s - hh) - hh_error) + s_lo;           // s − h², à 2^-106 près
        double const result = h + residual / fmul(2.0, h);
        return scale_by_power_of_two(result, exponent);                 // un seul arrondi s'il est sous-normal
    }

    // -----------------------------------------------------------------------
    // float : calcul en double, puis un seul arrondi vers float.
    // -----------------------------------------------------------------------
    constexpr float sqrt(float x) noexcept
    {
        // √ en double puis arrondi en float : arrondi correct (double précision
        // plus que double de celle du float, pas de double arrondi possible).
        return static_cast<float>(sqrt(static_cast<double>(x)));
    }
    constexpr float cbrt(float x) noexcept { return static_cast<float>(cbrt(static_cast<double>(x))); }
    constexpr float hypot(float x, float y) noexcept
    {
        return static_cast<float>(hypot(static_cast<double>(x), static_cast<double>(y)));
    }

    template <std::integral I> constexpr double sqrt(I x) noexcept { return sqrt(static_cast<double>(x)); }
    template <std::integral I> constexpr double cbrt(I x) noexcept { return cbrt(static_cast<double>(x)); }

    long double sqrt(long double) = delete;
    long double cbrt(long double) = delete;
    long double hypot(long double, long double) = delete;

    // Arguments mêlés, comme <cmath>.
    template <class A, class B> requires detail::mixed_arithmetic<A, B>
    constexpr double hypot(A x, B y) noexcept { return hypot(static_cast<double>(x), static_cast<double>(y)); }
}
