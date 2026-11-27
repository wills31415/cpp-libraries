#pragma once

#include <wsl/modes.hpp>

// Fonctions trigonométriques bit-exactes : même résultat au bit près sur toute
// plateforme (voir detail/fp.hpp pour les conditions), constexpr.
//
//   sin, cos, tan, atan, atan2
//
// En double, erreur < 0,51 ulp sur tout le domaine, arguments immenses compris ;
// en float, calcul en double puis un seul arrondi. Cas limites ceux de <cmath>
// (Annexe F de C) : sin(±0) = ±0, tan(±0) = ±0, cos(±0) = 1, atan(±0) = ±0,
// atan(±∞) = ±π/2 ; sin, cos, tan de ±∞ : NaN ; atan2 selon les signes et les
// infinis (voir atan2) ; NaN donne NaN.
//
// Règle d'écriture : tout produit dont le résultat entre dans une addition ou une
// soustraction passe par fmul (voir detail/fp.hpp).

#include <wsl/floating_point_arithmetic/detail/primitives.hpp>
#include <wsl/elementary_functions/detail/trigonometric_tables.hpp>

#include <bit>
#include <concepts>
#include <cstdint>

// Couche du dessous : les primitives flottantes.
namespace wsl::inline WSL_DEFAULTS::elementary_functions::detail
{
    using namespace floating_point_arithmetic::detail;
}

namespace wsl::inline WSL_DEFAULTS::elementary_functions
{
    namespace detail
    {
        // -------------------------------------------------------------------
        // Réduction : x = (4n + quadrant) π/2 + r, |r| <= π/4 (à un epsilon
        // près), r en double-double, précis à 2^-70 près relativement même
        // quand r est minuscule (il existe des doubles dont le reste vaut à peine
        // 2^-61 ; le pire connu est 6381956970095103 · 2^797).
        // -------------------------------------------------------------------
        struct trig_reduction
        {
            int quadrant;
            double r_hi, r_lo;
        };

        // |x| < 2^20 : Cody–Waite, π/2 en cinq morceaux de 33 bits (2^-165 près) :
        // chaque k pio2_n est exact (|k| < 2^20, 20 + 33 bits), et les trois
        // premiers sont retranchés exactement (two_sum) ; les deux derniers, de
        // l'ordre de 2^-99 au plus, ne pèsent que sur la queue.
        constexpr trig_reduction reduce_trig_medium(double x) noexcept
        {
            constexpr double shift = 0x1.8p52;
            double const t = fmul(x, inv_pio2) + shift;
            auto const k = static_cast<std::int64_t>(bits(t) - bits(shift));   // round(2x/π)
            double const kd = t - shift;

            double const a = x - fmul(kd, pio2_1);                  // exact (Sterbenz)
            auto const [t2, e2] = two_sum(a, -fmul(kd, pio2_2));
            auto const [t3, e3] = two_sum(t2, -fmul(kd, pio2_3));
            double const lo = (e2 + e3) - (fmul(kd, pio2_4) + fmul(kd, pio2_5));
            auto const [r_hi, r_lo] = two_sum(t3, lo);
            return {static_cast<int>(k & 3), r_hi, r_lo};
        }

        // |x| >= 2^20 : Payne–Hanek. x = m 2^e (m entier de 53 bits) ; x · 2/π
        // modulo 4 ne dépend que des bits de 2/π de rang >= e − 1 (les précédents
        // donnent des multiples de 4). On multiplie m par une fenêtre de 256 bits
        // de 2/π, en entiers 128 bits : produit exact, partie entière (le
        // quadrant) et fraction (2^-254 près) exactes, donc déterministes.
        constexpr std::uint64_t two_over_pi_chunk(int position) noexcept
        {
            int const word = position >> 6, offset = position & 63;
            if (offset == 0)
                return two_over_pi_bits[word];
            return (two_over_pi_bits[word] << offset) | (two_over_pi_bits[word + 1] >> (64 - offset));
        }

        constexpr trig_reduction reduce_trig_large(double ax) noexcept   // ax >= 2^20, fini
        {
            std::uint64_t const b = bits(ax);
            int const e = static_cast<int>(b >> 52) - 1075;               // ax = m 2^e
            std::uint64_t const m = (b & 0x000fffffffffffff) | 0x0010000000000000;

            // Fenêtre : bits de rang e − 1 à e + 254 de 2/π (position e + 62 dans
            // la suite précédée d'un mot nul). x · 2/π ≡ m · W · 2^-254 (mod 4).
            int const position = e + 62;
            std::uint64_t const w0 = two_over_pi_chunk(position);
            std::uint64_t const w1 = two_over_pi_chunk(position + 64);
            std::uint64_t const w2 = two_over_pi_chunk(position + 128);
            std::uint64_t const w3 = two_over_pi_chunk(position + 192);

            uint128 product = uint128{m} * w3;
            std::uint64_t const l0 = static_cast<std::uint64_t>(product);
            product = uint128{m} * w2 + (product >> 64);
            std::uint64_t l1 = static_cast<std::uint64_t>(product);
            product = uint128{m} * w1 + (product >> 64);
            std::uint64_t l2 = static_cast<std::uint64_t>(product);
            product = uint128{m} * w0 + (product >> 64);
            std::uint64_t l3 = static_cast<std::uint64_t>(product);   // bits 192 … 255 (le reste : multiples de 4)

            int quadrant = static_cast<int>(l3 >> 62);
            l3 &= (std::uint64_t{1} << 62) - 1;                         // fraction : 254 bits
            std::uint64_t l0f = l0;
            bool negative = false;
            if ((l3 >> 61) != 0)                                        // fraction >= 1/2 : quadrant suivant, reste négatif
            {
                ++quadrant;
                negative = true;
                // 2^254 − fraction, sur quatre mots.
                l0f = ~l0f + 1;
                std::uint64_t carry = l0f == 0 ? 1 : 0;
                l1 = ~l1 + carry;
                carry = (carry != 0 && l1 == 0) ? 1 : 0;
                l2 = ~l2 + carry;
                carry = (carry != 0 && l2 == 0) ? 1 : 0;
                l3 = (~l3 + carry) & ((std::uint64_t{1} << 62) - 1);
            }

            // Fraction en double-double : les 106 bits qui suivent le premier bit à 1.
            std::uint64_t const limbs[4] = {l3, l2, l1, l0f};
            int first = 0;
            while (first < 3 && limbs[first] == 0)
                ++first;
            std::uint64_t const a0 = limbs[first];
            std::uint64_t const a1 = first + 1 < 4 ? limbs[first + 1] : 0;
            std::uint64_t const a2 = first + 2 < 4 ? limbs[first + 2] : 0;
            int const c = std::countl_zero(a0);
            std::uint64_t const top = c == 0 ? a0 : (a0 << c) | (a1 >> (64 - c));
            std::uint64_t const next = c == 0 ? a1 : (a1 << c) | (a2 >> (64 - c));
            int const weight = 64 * (3 - first) - c - 254;              // poids du bit de tête : 2^(weight + 63)
            double const f_hi = fmul(static_cast<double>(top >> 11), power_of_two(weight + 11));   // exacts
            double const f_lo = fmul(static_cast<double>(((top & 0x7ff) << 42) | (next >> 22)), power_of_two(weight - 42));

            // r = fraction · π/2.
            auto const [p, p_error] = two_product(f_hi, pio2_hi);
            double const lo = p_error + (fmul(f_hi, pio2_lo) + fmul(f_lo, pio2_hi));
            auto const [r_hi, r_lo] = fast_two_sum(p, lo);
            if (negative)
                return {quadrant & 3, -r_hi, -r_lo};
            return {quadrant & 3, r_hi, r_lo};
        }

        // ax >= 0 fini.
        constexpr trig_reduction reduce_trig(double ax) noexcept
        {
            if (ax <= 0x1.921fb54442d18p-1)                             // <= π/4 : rien à réduire
                return {0, ax, 0.0};
            if (ax < 0x1p20)
                return reduce_trig_medium(ax);
            return reduce_trig_large(ax);
        }

        // -------------------------------------------------------------------
        // sin(r) et cos(r), |r| <= π/4 + ε, en double-double.
        //
        // r = c + d, c = i/64 (exact), |d| <= 1/128, et
        //     sin(r) = S + C d + S (cos d − 1) + C (sin d − d)
        //     cos(r) = C − S d + C (cos d − 1) − S (sin d − d)
        // avec S = sin c, C = cos c tabulés en double-double. Le terme d'ordre 1
        // est formé exactement (produit de Dekker, two_sum) ; sin d − d et
        // cos d − 1 en Taylor (restes < 2^-74 et 2^-91). sin_of et cos_of ne
        // calculent que ce qu'on leur demande : sin et cos n'en appellent qu'un.
        // -------------------------------------------------------------------
        struct sincos_point
        {
            bool negative;
            double s_hi, s_lo, c_hi, c_lo;      // sin c, cos c
            double d, d_lo, sm, cm;             // d, sin d − d, cos d − 1
        };

        constexpr sincos_point sincos_prepare(double r_hi, double r_lo) noexcept
        {
            bool const negative = r_hi < 0;
            double const a_hi = negative ? -r_hi : r_hi;
            double const a_lo = negative ? -r_lo : r_lo;

            constexpr double shift = 0x1.8p52;
            double const t = fmul(a_hi, 64.0) + shift;                 // produit exact
            auto const i = static_cast<int>(bits(t) - bits(shift));     // dans [0, 51]
            double const d = a_hi - fmul(t - shift, 0x1p-6);           // exact
            double const d_lo = a_lo;

            double const d2 = d * d;
            // sin d − d = d³ (−1/6 + d²/120 − d⁴/5040)
            double const sm = fmul(d * d2, -0x1.5555555555555p-3 + fmul(d2, 0x1.1111111111111p-7
                                        + fmul(d2, -0x1.a01a01a01a01ap-13)));
            // cos(d + d_lo) − 1 = d² (−1/2 + d²/24 − d⁴/720 + d⁶/40320) − d d_lo
            double const cm = fmul(d2, -0x1p-1 + fmul(d2, 0x1.5555555555555p-5 + fmul(d2, -0x1.6c16c16c16c17p-10
                                        + fmul(d2, 0x1.a01a01a01a01ap-16)))) - fmul(d, d_lo);

            auto const& entry = sincos_table[i];
            return {negative, entry[0], entry[1], entry[2], entry[3], d, d_lo, sm, cm};
        }

        // sin(r) = S + C d + S (cos d − 1) + C (sin d − d)
        constexpr double_double sin_of(sincos_point const& u) noexcept
        {
            auto const [p, p_error] = two_product(u.c_hi, u.d);
            auto const [s0, s_error] = two_sum(u.s_hi, p);
            double const tail = ((s_error + p_error) + (u.s_lo + (fmul(u.c_hi, u.d_lo) + fmul(u.c_lo, u.d))))
                              + (fmul(u.s_hi, u.cm) + fmul(u.c_hi, u.sm));
            auto const [hi, lo] = fast_two_sum(s0, tail);
            return u.negative ? double_double{-hi, -lo} : double_double{hi, lo};
        }

        // cos(r) = C − S d + C (cos d − 1) − S (sin d − d)
        constexpr double_double cos_of(sincos_point const& u) noexcept
        {
            auto const [q, q_error] = two_product(u.s_hi, u.d);
            auto const [c0, c_error] = two_sum(u.c_hi, -q);
            double const tail = ((c_error - q_error) + (u.c_lo - (fmul(u.s_hi, u.d_lo) + fmul(u.s_lo, u.d))))
                              + (fmul(u.c_hi, u.cm) - fmul(u.s_hi, u.sm));
            auto const [hi, lo] = fast_two_sum(c0, tail);
            return {hi, lo};
        }

        // n / d en double-double, n et d en double-double.
        constexpr double divide(double_double n, double_double d) noexcept
        {
            double const q = n.hi / d.hi;
            auto const [p, p_error] = two_product(q, d.hi);
            double const q_lo = (((n.hi - p) - p_error) + (n.lo - fmul(q, d.lo))) / d.hi;
            return q + q_lo;
        }

        // -------------------------------------------------------------------
        // atan(v), v = v_hi + v_lo dans [0, 1], en double-double.
        //
        // c = i/64 au plus près de v, t = (v − c)/(1 + v c) en double-double,
        // |t| <= 1/128, et atan(v) = atan(c) + atan(t), atan(c) tabulé en
        // double-double, atan(t) − t en Taylor (reste < 2^-73 relativement à t).
        // -------------------------------------------------------------------
        constexpr double_double atan_kernel(double v_hi, double v_lo) noexcept
        {
            constexpr double shift = 0x1.8p52;
            double const t = fmul(v_hi, 64.0) + shift;
            auto const i = static_cast<int>(bits(t) - bits(shift));     // dans [0, 64]
            double t_hi = v_hi, t_lo = v_lo;
            if (i != 0)
            {
                double const c = fmul(t - shift, 0x1p-6);
                double const n = v_hi - c;                              // exact (Sterbenz)
                auto const [m, m_error] = two_product(v_hi, c);
                auto const [den, den_error] = two_sum(1.0, m);
                double const den_lo = (den_error + m_error) + fmul(v_lo, c);
                t_hi = n / den;
                auto const [p, p_error] = two_product(t_hi, den);
                t_lo = (((n - p) - p_error) + (v_lo - fmul(t_hi, den_lo))) / den;
            }
            double const t2 = t_hi * t_hi;
            // atan t − t = t³ (−1/3 + t²/5 − t⁴/7 + t⁶/9)
            double const at = fmul(t_hi * t2, -0x1.5555555555555p-2 + fmul(t2, 0x1.999999999999ap-3
                                        + fmul(t2, -0x1.2492492492492p-3 + fmul(t2, 0x1.c71c71c71c71cp-4))));
            auto const& entry = atan_table[i];
            auto const [s, s_error] = two_sum(entry[0], t_hi);
            auto const [hi, lo] = fast_two_sum(s, (s_error + entry[1]) + (t_lo + at));
            return {hi, lo};
        }

        // (c_hi + c_lo) − a, en double-double.
        constexpr double_double subtract_from(double c_hi, double c_lo, double_double a) noexcept
        {
            auto const [u, u_error] = two_sum(c_hi, -a.hi);
            auto const [hi, lo] = fast_two_sum(u, (u_error + c_lo) - a.lo);
            return {hi, lo};
        }

        constexpr bool is_infinite(double x) noexcept { return x == infinity || x == -infinity; }
    }

    // -----------------------------------------------------------------------
    // sin, cos, tan
    // -----------------------------------------------------------------------
    constexpr double sin(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return detail::opaque(x + x);
        if (is_infinite(x))
            return quiet_nan;
        double const ax = x < 0 ? -x : x;
        if (ax < 0x1p-27)                                               // sin x = x (1 − x²/6), x²/6 < 2^-55
            return x;
        auto const [quadrant, r_hi, r_lo] = reduce_trig(ax);
        auto const point = sincos_prepare(r_hi, r_lo);
        auto const [hi, lo] = (quadrant & 1) == 0 ? sin_of(point) : cos_of(point);   // seul celui qui sert
        double const y = quadrant < 2 ? hi + lo : -(hi + lo);
        return x < 0 ? -y : y;
    }

    constexpr double cos(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return detail::opaque(x + x);
        if (is_infinite(x))
            return quiet_nan;
        double const ax = x < 0 ? -x : x;
        if (ax < 0x1p-27)                                               // cos x = 1 − x²/2, x²/2 < 2^-55
            return 1.0;
        auto const [quadrant, r_hi, r_lo] = reduce_trig(ax);
        auto const point = sincos_prepare(r_hi, r_lo);
        auto const [hi, lo] = (quadrant & 1) == 0 ? cos_of(point) : sin_of(point);
        return (quadrant == 0 || quadrant == 3) ? hi + lo : -(hi + lo);
    }

    // tan = sin/cos, ou −cos/sin dans les quadrants impairs, la division en
    // double-double.
    constexpr double tan(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return detail::opaque(x + x);
        if (is_infinite(x))
            return quiet_nan;
        double const ax = x < 0 ? -x : x;
        if (ax < 0x1p-27)                                               // tan x = x (1 + x²/3), x²/3 < 2^-55
            return x;
        auto const [quadrant, r_hi, r_lo] = reduce_trig(ax);
        auto const point = sincos_prepare(r_hi, r_lo);
        double_double const s = sin_of(point), c = cos_of(point);
        double const y = (quadrant & 1) == 0 ? divide(s, c) : -divide(c, s);
        return x < 0 ? -y : y;
    }

    // -----------------------------------------------------------------------
    // atan, atan2
    // -----------------------------------------------------------------------
    constexpr double atan(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return detail::opaque(x + x);
        double const ax = x < 0 ? -x : x;
        if (ax < 0x1p-27)                                               // atan x = x (1 − x²/3)
            return x;
        double y = 0;
        if (ax > 0x1p60)                                                // π/2 − 1/x, 1/x sous le demi-ulp
            y = pio2_hi;
        else if (ax <= 1)
        {
            auto const [hi, lo] = atan_kernel(ax, 0.0);
            y = hi + lo;
        }
        else                                                            // π/2 − atan(1/x), 1/x en double-double
        {
            double const v = 1.0 / ax;
            auto const [m, m_error] = two_product(v, ax);
            double const v_lo = ((1.0 - m) - m_error) / ax;
            auto const [hi, lo] = subtract_from(pio2_hi, pio2_lo, atan_kernel(v, v_lo));
            y = hi + lo;
        }
        return x < 0 ? -y : y;
    }

    // Angle de (x, y), dans [−π, π]. Cas à part de l'Annexe F de C :
    //   atan2(±0, −0) = ±π, atan2(±0, +0) = ±0, atan2(±0, x < 0) = ±π,
    //   atan2(±0, x > 0) = ±0, atan2(y ≠ 0, ±0) = ±π/2,
    //   atan2(±y, −∞) = ±π, atan2(±y, +∞) = ±0 (y fini > 0),
    //   atan2(±∞, x fini) = ±π/2, atan2(±∞, −∞) = ±3π/4, atan2(±∞, +∞) = ±π/4.
    // Sinon, rapport q <= 1 des valeurs absolues en double-double, atan(q),
    // puis π/2 − ·, π − · selon l'octant, toujours en double-double.
    constexpr double atan2(double y, double x) noexcept
    {
        using namespace detail;
        if (!(x == x) || !(y == y))
            return x + y;
        bool const y_negative = (bits(y) >> 63) != 0;
        bool const x_negative = (bits(x) >> 63) != 0;
        auto const with_sign = [y_negative](double v) { return y_negative ? -v : v; };

        if (y == 0)
            return x_negative ? with_sign(pi_hi) : y;
        if (x == 0)
            return with_sign(pio2_hi);
        if (is_infinite(x))
        {
            if (is_infinite(y))
                return with_sign(x_negative ? three_pio4 : pio4);
            return with_sign(x_negative ? pi_hi : 0.0);
        }
        if (is_infinite(y))
            return with_sign(pio2_hi);

        // Mise à l'échelle par une puissance de 2 (seul le rapport compte) : ni
        // débordement (le découpage de Veltkamp de d multiplie par 2^27 + 1 :
        // d doit rester sous 2^995), ni produits sous-normaux dans le produit de
        // Dekker q·d ≈ n, qui y deviendrait inexact. Si la plus petite valeur
        // reste minuscule quand la plus grande interdit de remonter, leur
        // rapport est sous 2^-1500 : nul ou presque en double, le résultat est
        // exact sans correction.
        double ay = y_negative ? -y : y;
        double ax = x_negative ? -x : x;
        double const largest = ay > ax ? ay : ax;
        double const smallest = ay > ax ? ax : ay;
        if (largest > 0x1p900)
        {
            ay = fmul(ay, 0x1p-200);                                    // exacts, mais entrent dans n − p
            ax = fmul(ax, 0x1p-200);
        }
        else if (smallest < 0x1p-800 && largest < 0x1p700)
        {
            ay = fmul(ay, 0x1p200);
            ax = fmul(ax, 0x1p200);
        }

        bool const swap = ay > ax;
        double const n = swap ? ax : ay;
        double const d = swap ? ay : ax;
        double const q = n / d;                                         // dans [0, 1]
        double q_lo = 0;
        if (n >= 0x1p-800)                                              // produits partiels normaux
        {
            auto const [p, p_error] = two_product(q, d);
            q_lo = ((n - p) - p_error) / d;
        }
        double_double a = atan_kernel(q, q_lo);
        if (swap)
            a = subtract_from(pio2_hi, pio2_lo, a);
        if (x_negative)
            a = subtract_from(pi_hi, pi_lo, a);
        return with_sign(a.hi + a.lo);
    }

    // -----------------------------------------------------------------------
    // float : calcul en double, puis un seul arrondi vers float.
    // -----------------------------------------------------------------------
    constexpr float sin(float x) noexcept { return static_cast<float>(sin(static_cast<double>(x))); }
    constexpr float cos(float x) noexcept { return static_cast<float>(cos(static_cast<double>(x))); }
    constexpr float tan(float x) noexcept { return static_cast<float>(tan(static_cast<double>(x))); }
    constexpr float atan(float x) noexcept { return static_cast<float>(atan(static_cast<double>(x))); }
    constexpr float atan2(float y, float x) noexcept
    {
        return static_cast<float>(atan2(static_cast<double>(y), static_cast<double>(x)));
    }

    template <std::integral I> constexpr double sin(I x) noexcept { return sin(static_cast<double>(x)); }
    template <std::integral I> constexpr double cos(I x) noexcept { return cos(static_cast<double>(x)); }
    template <std::integral I> constexpr double tan(I x) noexcept { return tan(static_cast<double>(x)); }
    template <std::integral I> constexpr double atan(I x) noexcept { return atan(static_cast<double>(x)); }

    long double sin(long double) = delete;
    long double cos(long double) = delete;
    long double tan(long double) = delete;
    long double atan(long double) = delete;
    long double atan2(long double, long double) = delete;

    // Arguments mêlés, comme <cmath>.
    template <class A, class B> requires detail::mixed_arithmetic<A, B>
    constexpr double atan2(A y, B x) noexcept { return atan2(static_cast<double>(y), static_cast<double>(x)); }
}
