#pragma once

// Fonction d'erreur bit-exacte : même résultat au bit près sur toute plateforme
// (voir detail/fp.hpp pour les conditions), constexpr.
//
//   erf, erfc
//
// Cas limites ceux de <cmath> (Annexe F de C) : erf(±0) = ±0, erf(±∞) = ±1,
// erfc(+∞) = +0, erfc(−∞) = 2 ; NaN donne NaN.
//
// Règle d'écriture : tout produit dont le résultat entre dans une addition ou une
// soustraction passe par fmul (voir detail/fp.hpp).

#include <wsl/special/detail/erf_tables.hpp>
#include <wsl/special/detail/fp.hpp>
#include <wsl/special/elementary.hpp>

#include <concepts>
#include <cstdint>

namespace wsl::special
{
    namespace detail
    {
        // erf(x) pour 2^-900 <= |x| < 1/4, en double-double :
        //     erf(x) = (2/√π) x (1 + x² Q(x²)),
        // le terme principal (2/√π) x formé exactement (produit de Dekker, ses
        // produits partiels restant normaux), x² Q(x²) ≈ −x²/3 ne pesant que 2 %.
        constexpr double_double erf_small(double x) noexcept
        {
            auto const [p, p_error] = two_product(x, two_over_sqrt_pi_hi);
            double const u = x * x;
            constexpr int terms = sizeof(erf_small_q) / sizeof(erf_small_q[0]);
            double q = erf_small_q[terms - 1];
            for (int n = terms - 2; n >= 0; --n)
                q = erf_small_q[n] + fmul(u, q);
            double const tail = (p_error + fmul(x, two_over_sqrt_pi_lo)) + fmul(p, fmul(u, q));
            auto const [hi, lo] = fast_two_sum(p, tail);
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
            return x + x;
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
            return x + x;
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

    // float : calcul en double, puis un seul arrondi vers float.
    constexpr float erf(float x) noexcept { return static_cast<float>(erf(static_cast<double>(x))); }
    constexpr float erfc(float x) noexcept { return static_cast<float>(erfc(static_cast<double>(x))); }

    template <std::integral I> constexpr double erf(I x) noexcept { return erf(static_cast<double>(x)); }
    template <std::integral I> constexpr double erfc(I x) noexcept { return erfc(static_cast<double>(x)); }

    long double erf(long double) = delete;
    long double erfc(long double) = delete;
}
