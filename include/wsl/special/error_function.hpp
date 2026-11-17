#pragma once

// Fonction d'erreur bit-exacte : même résultat au bit près sur toute plateforme
// (voir detail/fp.hpp pour les conditions), constexpr.
//
//   erf, erfc, erf_inv, erfc_inv
//
// Cas limites ceux de <cmath> (Annexe F de C) : erf(±0) = ±0, erf(±∞) = ±1,
// erfc(+∞) = +0, erfc(−∞) = 2 ; NaN donne NaN. Inverses : erf_inv(±0) = ±0,
// erf_inv(±1) = ±∞, erf_inv hors de [−1, 1] : NaN ; erfc_inv(0) = +∞,
// erfc_inv(1) = +0, erfc_inv(2) = −∞, erfc_inv hors de [0, 2] : NaN.
//
// Règle d'écriture : tout produit dont le résultat entre dans une addition ou une
// soustraction passe par fmul (voir detail/fp.hpp).

#include <wsl/special/algebraic.hpp>
#include <wsl/special/detail/erf_tables.hpp>
#include <wsl/special/detail/fp.hpp>
#include <wsl/special/elementary.hpp>

#include <concepts>
#include <cstddef>
#include <cstdint>

namespace wsl::special
{
    namespace detail
    {
        // erf(x) pour 2^-900 <= |x| < 1/4, en double-double :
        //     erf(x) = (2/√π) x (1 + x² Q(x²)),   Q(u) = −1/3 + u (1/10 − …),
        // à 2^-63 près relativement : le terme principal (2/√π) x et la correction
        // principale p · x² · (−1/3) sont formés exactement (produits de Dekker,
        // x² et −1/3 en double-double), le reste ne pesant que 2^-11. (Avec x²
        // arrondi, l'erreur relative était de 2^-59 : assez pour erf, pas pour
        // erf_inv, qui la reporte telle quelle sur son résultat.)
        constexpr double_double erf_small(double x) noexcept
        {
            auto const [p, p_error] = two_product(x, two_over_sqrt_pi_hi);
            if (x > -0x1p-300 && x < 0x1p-300)                              // x² Q < 2^-600 : sans effet, et x², x³
            {                                                               // tomberaient dans les sous-normaux
                auto const [hi, lo] = fast_two_sum(p, p_error + fmul(x, two_over_sqrt_pi_lo));
                return {hi, lo};
            }
            auto const [u, u_lo] = two_product(x, x);
            constexpr int terms = sizeof(erf_small_q) / sizeof(erf_small_q[0]);
            double rest = erf_small_q[terms - 1];                           // q_2 + q_3 u + …
            for (int n = terms - 2; n >= 1; --n)
                rest = erf_small_q[n] + fmul(u, rest);
            double const q1 = erf_small_q[0];                               // −1/3 arrondi
            auto const [m, m_error] = two_product(u, q1);                   // u q_1, exact
            double const uq_lo = (m_error + (fmul(u_lo, q1) + fmul(u, erf_small_q1_lo))) + fmul(u * u, rest);
            auto const [c, c_error] = two_product(p, m);                    // p · x² Q, terme de tête exact
            auto const [s, s_error] = two_sum(p, c);
            double const tail = (s_error + (p_error + fmul(x, two_over_sqrt_pi_lo)))
                              + (c_error + fmul(p, uq_lo));
            auto const [hi, lo] = fast_two_sum(s, tail);
            return {hi, lo};
        }

        // F(x) = e^{x²} erfc(x), lisse et lentement variable, pour x dans
        // [1/4, 28), en double-double : par morceaux (largeur 1/8 sous 4, 1/2 au
        // delà), Taylor de degré 15 autour du milieu c (reste < 2^-66, vérifié par
        // le générateur). t = x − c est exact (Sterbenz) ; le terme a_1 t est
        // formé exactement, le reste ne pesant que quelques pour cent.
        constexpr double_double erfc_factor(double x) noexcept
        {
            int i = 0;
            double c = 0;
            if (x < 4.0)
            {
                i = static_cast<int>(fmul(x - 0.25, 8.0));
                c = fmul(static_cast<double>(2 * i + 1), 0x1p-4) + 0.25;   // exact
            }
            else
            {
                int const j = static_cast<int>(fmul(x - 4.0, 2.0));
                i = erfc_fine_pieces + j;
                c = fmul(static_cast<double>(2 * j + 1), 0x1p-2) + 4.0;    // exact
            }
            double const t = x - c;                                         // exact
            auto const& a = erfc_table[i];
            constexpr int last = 4 + erfc_degree - 2;                       // a_15
            double h = a[last];
            for (int k = last - 1; k >= 4; --k)                             // a_2 + a_3 t + … + a_15 t^13
                h = a[k] + fmul(t, h);
            auto const [q, q_error] = two_product(a[2], t);
            auto const [s, s_error] = two_sum(a[0], q);
            double const tail = ((s_error + a[1]) + (q_error + fmul(a[3], t))) + fmul(t * t, h);
            auto const [hi, lo] = fast_two_sum(s, tail);
            return {hi, lo};
        }

        // erfc(x) = 2^e (hi + lo), hi dans [1, 2), pour x dans [1/4, 28) :
        // e^{-x²} F(x), x² exact (produit de Dekker), e^{-x²} par le noyau de exp à
        // argument double-double, le produit en double-double.
        struct scaled_double_double
        {
            double hi, lo;
            int e;
        };

        constexpr scaled_double_double erfc_scaled(double x) noexcept
        {
            auto const [s, s_lo] = two_product(x, x);
            auto const [k, r_hi, r_lo] = reduce_exp(-s);
            auto const [t_hi, tail, e] = exp_parts_dd(k, r_hi, r_lo - s_lo);
            auto const [f_hi, f_lo] = erfc_factor(x);
            auto const [m, m_error] = two_product(t_hi, f_hi);
            double const m_tail = m_error + (fmul(t_hi, f_lo) + fmul(tail, f_hi));
            int const shift = 1023 - static_cast<int>(bits(m) >> 52);       // m · 2^shift dans [1, 2)
            double const scale = power_of_two(shift);
            return {fmul(m, scale), fmul(m_tail, scale), e - shift};
        }

        // erfc(x) en double-double, pour x dans [1/4, 6] (résultat normal).
        constexpr double_double erfc_dd(double x) noexcept
        {
            auto const [hi, lo, e] = erfc_scaled(x);
            double const scale = power_of_two(e);
            return {fmul(hi, scale), fmul(lo, scale)};
        }
    }

    // erf : série pour |x| < 1/4, 1 − erfc(|x|) au-delà (l'annulation ne coûte
    // qu'un facteur 2,6 au plus, sur une erreur de 2^-59), 1 dès que erfc(|x|) est
    // sous le demi-ulp de 1.
    constexpr double erf(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return detail::opaque(x + x);
        if (x == 0)
            return x;                                                       // garde le signe de −0
        double const ax = x < 0 ? -x : x;
        if (ax >= 6.0)                                                      // erfc(6) < 2^-54
            return x < 0 ? -1.0 : 1.0;
        if (ax < 0.25)
        {
            if (ax < 0x1p-900)                                              // produits partiels normaux
            {
                auto const [hi, lo] = erf_small(fmul(x, 0x1p200));          // exact
                return fmul(hi + lo, 0x1p-200);
            }
            auto const [hi, lo] = erf_small(x);
            return hi + lo;
        }
        auto const [c_hi, c_lo] = erfc_dd(ax);
        auto const [u, u_error] = two_sum(1.0, -c_hi);
        double const y = u + (u_error - c_lo);
        return x < 0 ? -y : y;
    }

    // erfc : 1 − erf(x) pour |x| < 1/4, e^{-x²} F(x) pour x >= 1/4 (un seul
    // arrondi, résultats sous-normaux compris), 2 − erfc(−x) pour x <= −1/4.
    constexpr double erfc(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return detail::opaque(x + x);
        if (x >= 28.0)                                                      // erfc(28) < 2^-1100 ; +∞ compris
            return 0.0;
        double const ax = x < 0 ? -x : x;
        if (ax < 0.25)
        {
            if (ax < 0x1p-60)                                               // |erf(x)| < 2^-59
                return 1.0;
            auto const [e_hi, e_lo] = erf_small(x);
            auto const [u, u_error] = two_sum(1.0, -e_hi);
            return u + (u_error - e_lo);
        }
        if (x > 0)
        {
            auto const [hi, lo, e] = erfc_scaled(x);
            return scale_exp(hi, lo, e);
        }
        if (ax >= 6.0)                                                      // 2 − erfc(6) s'arrondit à 2 ; −∞ compris
            return 2.0;
        auto const [c_hi, c_lo] = erfc_dd(ax);
        auto const [u, u_error] = two_sum(2.0, -c_hi);
        return u + (u_error - c_lo);
    }

    namespace detail
    {
        // Polynôme aux coefficients rangés du degré le plus haut au plus bas.
        template <std::size_t N>
        constexpr double horner(double const (&c)[N], double u) noexcept
        {
            double p = c[0];
            for (std::size_t i = 1; i < N; ++i)
                p = c[i] + fmul(u, p);
            return p;
        }

        // Approximation initiale de erf⁻¹, à 2^-35 près relativement, pour |y| < 1,
        // z = 1 − |y| : erf⁻¹(y) = |y| P, P en w = −log(z (2 − z)) au centre
        // (w <= 5), en s = √w dans les queues (six morceaux). z est fourni à part :
        // près de 1, calculer 1 − |y| ne garderait rien.
        constexpr double erfinv_initial(double ay, double z) noexcept
        {
            double const w = -log(z * (2.0 - z));
            if (w <= 5.0)
                return fmul(ay, horner(erfinv_central, fmul(w, 0.4) - 1.0));
            double const s = sqrt(w);
            int i = 0;
            while (i < 5 && s > erfinv_tail_bounds[i + 1])
                ++i;
            double const a = erfinv_tail_bounds[i], b = erfinv_tail_bounds[i + 1];
            double u = (s - fmul(a + b, 0.5)) / fmul(b - a, 0.5);
            u = u < -1.0 ? -1.0 : (u > 1.0 ? 1.0 : u);
            return fmul(ay, horner(erfinv_tail[i], u));
        }

        // Un pas de Newton sur erf, pour |y| <= 1/2 : le résidu erf(x) − |y| en
        // double-double (erf_small ou 1 − erfc), d'où x − résidu / erf'(x) ; parti
        // de 2^-35, le pas laisse 2^-70.
        constexpr double erfinv_refine_central(double ay, double x) noexcept
        {
            double e_hi = 0, e_lo = 0;
            if (x < 0.25)
            {
                auto const [hi, lo] = erf_small(x);
                e_hi = hi;
                e_lo = lo;
            }
            else
            {
                auto const [c_hi, c_lo] = erfc_dd(x);
                auto const [u, u_error] = two_sum(1.0, -c_hi);
                e_hi = u;
                e_lo = u_error - c_lo;
            }
            double const residual = (e_hi - ay) + e_lo;                       // e_hi − |y| exact
            double const slope = fmul(two_over_sqrt_pi_hi, exp(-(x * x)));
            return x - residual / slope;
        }

        // Un pas de Newton sur log erfc, pour z dans (0, 1/2] (sous-normaux
        // compris) : (log erfc)'(x) = −(2/√π)/F(x), d'où
        //     x ← x + log(erfc(x)/z) · F(x) · √π/2,
        // le rapport erfc(x)/z formé en double-double par mises à l'échelle
        // (erfc(x) = 2^e (h + l), z = 2^ez m), sans quitter le domaine normal.
        constexpr double erfcinv_refine(double z, double x) noexcept
        {
            x = x < 0.25 ? 0.25 : (x > 27.9 ? 27.9 : x);
            auto const [h, l, e] = erfc_scaled(x);
            double zs = z;
            int ez = 0;
            if (zs < 0x1p-1022)
            {
                zs = fmul(zs, 0x1p64);                                          // exact
                ez = -64;
            }
            std::uint64_t const bz = bits(zs);
            ez += static_cast<int>(bz >> 52) - 1023;
            double const m = from_bits((bz & 0x000fffffffffffff) | 0x3ff0000000000000);
            double const scale = power_of_two(e - ez);
            double const q = ((fmul(h, scale) - m) + fmul(l, scale)) / m;    // erfc(x)/z − 1
            double const f = erfc_factor(x).hi;
            return x + fmul(fmul(log1p(q), f), half_sqrt_pi_hi);
        }
    }

    // erf⁻¹ : |y| <= 1/2, approximation et un pas de Newton sur erf ; au-delà,
    // erfc⁻¹(1 − |y|), 1 − |y| exact (Sterbenz).
    constexpr double erf_inv(double y) noexcept
    {
        using namespace detail;
        if (!(y == y))
            return detail::opaque(y + y);
        double const ay = y < 0 ? -y : y;
        if (ay > 1.0)
            return quiet_nan;
        if (ay == 1.0)
            return y < 0 ? -infinity : infinity;
        if (y == 0)
            return y;                                                           // garde le signe de −0
        double x = 0;
        if (ay < 0x1p-60)                                                       // (√π/2) y (1 + π y²/12), y² négligeable
        {
            double const ys = ay < 0x1p-900 ? fmul(ay, 0x1p200) : ay;           // produits partiels normaux
            auto const [p, p_error] = two_product(ys, half_sqrt_pi_hi);
            double const r = p + (p_error + fmul(ys, half_sqrt_pi_lo));
            x = ay < 0x1p-900 ? fmul(r, 0x1p-200) : r;
        }
        else if (ay <= 0.5)
            x = erfinv_refine_central(ay, erfinv_initial(ay, 1.0 - ay));
        else
            x = erfcinv_refine(1.0 - ay, erfinv_initial(ay, 1.0 - ay));
        return y < 0 ? -x : x;
    }

    // erfc⁻¹ : z <= 1/2 par Newton sur log erfc ; 1/2 < z < 3/2, erf⁻¹(1 − z) ;
    // z >= 3/2, −erfc⁻¹(2 − z) (exacts : Sterbenz).
    constexpr double erfc_inv(double z) noexcept
    {
        using namespace detail;
        if (!(z == z))
            return detail::opaque(z + z);
        if (z < 0 || z > 2)
            return quiet_nan;
        if (z == 0)
            return infinity;
        if (z == 2)
            return -infinity;
        if (z <= 0.5)
            return erfcinv_refine(z, erfinv_initial(1.0 - z, z));
        if (z < 1.5)
            return erf_inv(1.0 - z);
        double const w = 2.0 - z;
        return -erfcinv_refine(w, erfinv_initial(1.0 - w, w));
    }

    // float : calcul en double, puis un seul arrondi vers float.
    constexpr float erf(float x) noexcept { return static_cast<float>(erf(static_cast<double>(x))); }
    constexpr float erfc(float x) noexcept { return static_cast<float>(erfc(static_cast<double>(x))); }

    template <std::integral I> constexpr double erf(I x) noexcept { return erf(static_cast<double>(x)); }
    template <std::integral I> constexpr double erfc(I x) noexcept { return erfc(static_cast<double>(x)); }
    constexpr float erf_inv(float y) noexcept { return static_cast<float>(erf_inv(static_cast<double>(y))); }
    constexpr float erfc_inv(float z) noexcept { return static_cast<float>(erfc_inv(static_cast<double>(z))); }

    long double erf(long double) = delete;
    long double erfc(long double) = delete;
    long double erf_inv(long double) = delete;
    long double erfc_inv(long double) = delete;
}
