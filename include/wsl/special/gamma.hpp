#pragma once

// Fonction gamma bit-exacte : même résultat au bit près sur toute plateforme
// (voir detail/fp.hpp pour les conditions), constexpr.
//
//   lgamma(x) = log |Γ(x)|      tgamma(x) = Γ(x)
//
// Cas limites ceux de <cmath> (Annexe F de C) : lgamma(1) = lgamma(2) = +0,
// lgamma de ±0, des entiers négatifs et de ±∞ : +∞ ; tgamma(±0) = ±∞,
// tgamma(entier négatif) = NaN, tgamma(−∞) = NaN, tgamma(+∞) = +∞ ; NaN donne NaN.
//
// Précision : voir docs/special.md. Pour x < 0, près des zéros de lgamma
// (x ≈ −2,457, −2,748, …), la formule de réflexion retranche deux termes presque
// égaux : l'erreur y est petite en valeur absolue, pas en ulp.
//
// Règle d'écriture : tout produit dont le résultat entre dans une addition ou une
// soustraction passe par fmul (voir detail/fp.hpp).

#include <wsl/special/detail/fp.hpp>
#include <wsl/special/detail/gamma_tables.hpp>
#include <wsl/special/elementary.hpp>
#include <wsl/special/trigonometric.hpp>

#include <concepts>
#include <cstdint>

namespace wsl::special
{
    namespace detail
    {
        // log x en double-double, x > 0 fini (sous-normal compris).
        constexpr double_double log_dd(double x) noexcept
        {
            std::uint64_t ix = bits(x);
            if (x < 0x1p-1022)
                ix = bits(fmul(x, 0x1p52)) - (std::uint64_t{52} << 52);
            return log_accurate(ix);
        }

        // lgamma(y_hi + y_lo) pour y dans [0,4375 ; 10,0625) : table, Taylor de
        // degré 26 autour de c = k/8 (reste < 2^-66, vérifié par le générateur).
        // t = y_hi − c est exact (Sterbenz) ; le terme a_1 t est formé exactement.
        constexpr double_double lgamma_piece(double y_hi, double y_lo) noexcept
        {
            constexpr double shift = 0x1.8p52;
            double const u = fmul(y_hi, 8.0) + shift;                       // produit exact
            auto k = static_cast<int>(bits(u) - bits(shift));              // round(8 y)
            // Près des zéros 1 et 2, le morceau centré sur le zéro couvre aussi ses
            // voisins (|t| <= 3/16) : a_0 = 0, pas d'annulation entre a_0 et a_1 t.
            if (k == 7 || k == 9)
                k = 8;
            else if (k == 15 || k == 17)
                k = 16;
            double const t = y_hi - fmul(static_cast<double>(k), 0.125);    // exact
            // a_0 + a_1 t + a_2 t² exacts (a_k en double-double, produits de Dekker) :
            // près des zéros, a_2 t² pèse jusqu'à 10 % du résultat. Le reste, t³ (…),
            // moins de 0,3 %.
            auto const& a = lgamma_table[k - lgamma_first];
            constexpr int last = 6 + lgamma_degree - 3;                     // a_26
            double h = a[last];
            for (int n = last - 1; n >= 6; --n)                             // a_3 + a_4 t + … + a_26 t^23
                h = a[n] + fmul(t, h);
            auto const [q1, q1_error] = two_product(a[2], t);
            auto const [tt, tt_error] = two_product(t, t);
            auto const [q2, q2_error] = two_product(a[4], tt);
            auto const [s1, e1] = two_sum(a[0], q1);
            auto const [s2, e2] = two_sum(s1, q2);
            double const tail = ((e1 + e2) + (a[1] + (q1_error + q2_error)))
                              + ((fmul(a[3], t) + (fmul(a[4], tt_error) + fmul(a[5], tt))) + fmul(a[2], y_lo))
                              + fmul(t * tt, h);
            auto const [hi, lo] = fast_two_sum(s2, tail);
            return {hi, lo};
        }

        // Stirling en double-double, y dans [10,0625 ; 2^900) :
        //     lgamma(y) = y log y − y − (log y)/2 + log(2π)/2 + Σ B_2k / (2k (2k − 1) y^(2k−1)),
        // y log y exact (produit de Dekker sur log y en double-double), les termes
        // de tête sommés exactement.
        constexpr double_double lgamma_stirling(double y) noexcept
        {
            auto const [l_hi, l_lo] = log_dd(y);
            auto const [p, p_error] = two_product(y, l_hi);
            auto const [s1, e1] = two_sum(p, -y);
            auto const [s2, e2] = two_sum(s1, fmul(-0.5, l_hi));            // −(log y)/2 : exact
            auto const [s3, e3] = two_sum(s2, half_log_two_pi_hi);
            double const z = 1.0 / y;
            double const z2 = z * z;
            constexpr int terms = sizeof(stirling) / sizeof(stirling[0]);
            double series = stirling[terms - 1];
            for (int n = terms - 2; n >= 0; --n)
                series = stirling[n] + fmul(z2, series);
            double const tail = ((e1 + e2) + (e3 + p_error))
                              + ((fmul(y, l_lo) - fmul(0.5, l_lo)) + half_log_two_pi_lo) + fmul(z, series);
            auto const [hi, lo] = fast_two_sum(s3, tail);
            return {hi, lo};
        }

        // lgamma(y) en double-double, y > 0, 2^-54 <= y <= seuil de débordement.
        constexpr double_double lgamma_positive(double y) noexcept
        {
            if (y < 0.5625)                                                  // lgamma(1 + y) − log y
            {
                auto const [u, u_lo] = two_sum(1.0, y);
                auto const [g_hi, g_lo] = lgamma_piece(u, u_lo);
                auto const [l_hi, l_lo] = log_dd(y);
                auto const [s, s_error] = two_sum(g_hi, -l_hi);
                auto const [hi, lo] = fast_two_sum(s, (s_error + g_lo) - l_lo);
                return {hi, lo};
            }
            if (y < 10.0625)
                return lgamma_piece(y, 0.0);
            if (y < 0x1p900)
                return lgamma_stirling(y);
            // y immense : lgamma(y) = y (log y − 1) à moins d'un ulp près, les autres
            // termes (quelques centaines) étant bien en dessous. Mise à l'échelle par
            // 2^-100 : le produit de Dekker n'admet pas de facteur au-delà de 2^995.
            auto const [l_hi, l_lo] = log_dd(y);
            auto const [m, m_error] = two_sum(l_hi, -1.0);
            double const ys = fmul(y, 0x1p-100);                             // exact
            auto const [p, p_error] = two_product(ys, m);
            double const tail = p_error + fmul(ys, m_error + l_lo);
            return {fmul(p, 0x1p100), fmul(tail, 0x1p100)};
        }

        // |sin(π r)| en double-double, |r| <= 1/2, r ≠ 0 : π r en double-double,
        // puis sin(π |r|) si |r| <= 1/4, cos(π (1/2 − |r|)) sinon (exact : Sterbenz).
        constexpr double_double sin_pi_abs(double r) noexcept
        {
            double const a = r < 0 ? -r : r;
            bool const near_half = a > 0.25;
            double const w = near_half ? 0.5 - a : a;
            auto const [p, p_error] = two_product(w, pi_hi);
            double const p_lo = p_error + fmul(w, pi_lo);
            auto const point = sincos_prepare(p, p_lo);
            return near_half ? cos_of(point) : sin_of(point);
        }

        // Réflexion, x < 0 non entier (−2^52 < x < −2^-54) :
        //     log |Γ(x)| = log π − log |sin(π x)| − lgamma(1 − x),
        // tout en double-double. 1 − x = y_hi + y_lo (two_sum) ; lgamma(y_hi + y_lo)
        // ≈ lgamma(y_hi) + ψ(y_hi) y_lo, ψ(y) ≈ log y − 1/(2y) suffisant pour cette
        // correction de l'ordre de 2^-53.
        struct reflection
        {
            double_double log_abs;
            bool negative;                                                   // signe de Γ(x)
        };

        // x entier ? (|x| < 2^52 ; au-delà, tout double est entier.)
        constexpr bool is_integer(double x) noexcept
        {
            return static_cast<double>(static_cast<std::int64_t>(x)) == x;
        }

        constexpr reflection lgamma_reflection(double x) noexcept
        {
            // x = n + r, n entier, |r| <= 1/2 : troncature (exacte sous 2^52), puis
            // repli du reste. (L'astuce du décalage n'arrondirait à l'entier que
            // sous 2^51 : au-delà, un double peut être demi-entier.)
            double n = static_cast<double>(static_cast<std::int64_t>(x));
            double r = x - n;                                                // exact, dans (−1, 0)
            if (r < -0.5)
            {
                r += 1.0;                                                    // exact (Sterbenz)
                n -= 1.0;
            }
            auto const [s_hi, s_lo] = sin_pi_abs(r);
            auto const [ls_hi, ls_lo0] = log_dd(s_hi);
            double const ls_lo = ls_lo0 + s_lo / s_hi;

            auto const [y_hi, y_lo] = two_sum(1.0, -x);
            auto const [g_hi, g_lo0] = lgamma_positive(y_hi);
            double const digamma = log(y_hi) - 0.5 / y_hi;
            double const g_lo = g_lo0 + fmul(digamma, y_lo);

            auto const [a, a_error] = two_sum(log_pi_hi, -ls_hi);
            auto const [b, b_error] = two_sum(a, -g_hi);
            auto const [hi, lo] = fast_two_sum(b, ((a_error + b_error) + log_pi_lo) - (ls_lo + g_lo));

            // Γ(x) < 0 exactement quand floor(x) est impair (r ≠ 0).
            auto const floor_x = static_cast<std::int64_t>(r > 0 ? n : n - 1);
            return {{hi, lo}, (floor_x & 1) != 0};
        }

        constexpr double exp_dd(double_double a) noexcept
        {
            auto const [k, r_hi, r_lo] = reduce_exp(a.hi);
            return exp_finish_dd(k, r_hi, r_lo + a.lo);
        }
    }

    constexpr double lgamma(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return x + x;
        if (x == infinity || x == -infinity || x == 0)
            return infinity;
        double const ax = x < 0 ? -x : x;
        if (ax < 0x1p-54)                                                    // −log|x| − γx, γx sous le demi-ulp
            return -log(ax);
        if (x > 0)
        {
            if (x > 0x1.754d9278b51a7p+1014)                                 // lgamma(x) > DBL_MAX
                return infinity;
            if (x == 1.0 || x == 2.0)
                return 0.0;
            auto const [hi, lo] = lgamma_positive(x);
            return hi + lo;
        }
        if (x <= -0x1p52 || is_integer(x))                                   // entier négatif : pôle
            return infinity;
        auto const [value, negative] = lgamma_reflection(x);
        return value.hi + value.lo;
    }

    constexpr double tgamma(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return x + x;
        if (x == infinity)
            return x;
        if (x == -infinity)
            return quiet_nan;
        if (x == 0)
            return 1.0 / x;                                                  // ±∞ selon le signe du zéro
        double const ax = x < 0 ? -x : x;
        if (ax < 0x1p-54)                                                    // Γ(x) = 1/x − γ + O(x)
        {
            if (ax < 0x1p-60)                                                // γ < 0,003 ulp de 1/x
                return 1.0 / x;
            double const q = 1.0 / x;                                        // 1/x en double-double
            auto const [p, p_error] = two_product(q, x);
            double const q_lo = ((1.0 - p) - p_error) / x;
            double const euler = -lgamma_table[8 - lgamma_first][2];         // γ = −ψ(1)
            return q + (q_lo - euler);
        }
        if (x > 0)
        {
            if (x > 0x1.573fae561f647p+7)                                    // Γ(x) > DBL_MAX
                return infinity;
            return exp_dd(lgamma_positive(x));
        }
        if (x <= -0x1p52 || is_integer(x))                                   // entier négatif
            return quiet_nan;
        auto const [value, negative] = lgamma_reflection(x);
        double const magnitude = value.hi < -800.0 ? 0.0 : exp_dd(value);  // |Γ(x)| < 2^-1150 : 0
        return negative ? -magnitude : magnitude;
    }

    // float : calcul en double, puis un seul arrondi vers float.
    constexpr float lgamma(float x) noexcept { return detail::to_float(lgamma(static_cast<double>(x))); }
    constexpr float tgamma(float x) noexcept
    {
        double const y = tgamma(static_cast<double>(x));
        if (y <= -0x1.ffffffp+127)                                           // < −FLT_MAX − ulp/2 : −∞
            return static_cast<float>(-detail::infinity);
        return detail::to_float(y);
    }

    template <std::integral I> constexpr double lgamma(I x) noexcept { return lgamma(static_cast<double>(x)); }
    template <std::integral I> constexpr double tgamma(I x) noexcept { return tgamma(static_cast<double>(x)); }

    long double lgamma(long double) = delete;
    long double tgamma(long double) = delete;
}
