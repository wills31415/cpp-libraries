#pragma once

#include <wsl/modes.hpp>
#include <wsl/special_functions/detail/moded.hpp>

// Fonction gamma bit-exacte : même résultat au bit près sur toute plateforme
// (voir detail/fp.hpp pour les conditions), constexpr.
//
//   lgamma(x) = log |Γ(x)|      tgamma(x) = Γ(x)      digamma(x) = ψ(x) = Γ'(x)/Γ(x)
//
// Cas limites ceux de <cmath> (Annexe F de C) : lgamma(1) = lgamma(2) = +0,
// lgamma de ±0, des entiers négatifs et de ±∞ : +∞ ; tgamma(±0) = ±∞,
// tgamma(entier négatif) = NaN, tgamma(−∞) = NaN, tgamma(+∞) = +∞ ; NaN donne NaN.
// digamma(+0) = −∞, digamma(−0) = +∞, digamma(entier négatif) = NaN (pôles),
// digamma(+∞) = +∞, digamma(−∞) = NaN.
//
// Précision : voir design/special.md. Pour x < 0, près des zéros de lgamma
// (x ≈ −2,457, −2,748, …), la formule de réflexion retranche deux termes presque
// égaux : l'erreur y est petite en valeur absolue, pas en ulp.
//
// Règle d'écriture (design/interfaces.md, décision 38) : tout produit dont le résultat entre
// dans une addition ou une soustraction passe par mul<R> ou, s'il entre directement dans une
// somme exacte ou un découpage, par unfused_mul ; tout appel porte son mode (<R>).

#include <wsl/floating_point_arithmetic/detail/primitives.hpp>
#include <wsl/special_functions/detail/digamma_tables.hpp>
#include <wsl/special_functions/detail/gamma_tables.hpp>
#include <wsl/elementary_functions/exponential.hpp>
#include <wsl/elementary_functions/trigonometric.hpp>

#include <concepts>
#include <cstdint>

#include <wsl/elementary_functions.hpp>

// Couches du dessous : primitives flottantes et fonctions élémentaires.
namespace wsl::inline WSL_DEFAULTS::special_functions::detail
{
    using namespace floating_point_arithmetic::detail;
    using namespace elementary_functions::detail;
}

namespace wsl::inline WSL_DEFAULTS::special_functions
{
    namespace detail
    {
        // log x en double-double, x > 0 fini (sous-normal compris).
        template <reproducibility_mode R = portable_t>
        constexpr double_double log_dd(double x) noexcept
        {
            std::uint64_t ix = bits(x);
            if (x < 0x1p-1022)
                ix = bits(mul<R>(x, 0x1p52)) - (std::uint64_t{52} << 52);
            return log_accurate<R>(ix);
        }

        // lgamma(y_hi + y_lo) pour y dans [0,4375 ; 10,0625) : table, Taylor de
        // degré 26 autour de c = k/8 (reste < 2^-66, vérifié par le générateur).
        // t = y_hi − c est exact (Sterbenz) ; le terme a_1 t est formé exactement.
        template <reproducibility_mode R = portable_t>
        constexpr double_double lgamma_piece(double y_hi, double y_lo) noexcept
        {
            constexpr double shift = 0x1.8p52;
            double const u = mul<R>(y_hi, 8.0) + shift;                       // produit exact
            auto k = static_cast<int>(bits(u) - bits(shift));              // round(8 y)
            // Près des zéros 1 et 2, le morceau centré sur le zéro couvre aussi ses
            // voisins (|t| <= 3/16) : a_0 = 0, pas d'annulation entre a_0 et a_1 t.
            if (k == 7 || k == 9)
                k = 8;
            else if (k == 15 || k == 17)
                k = 16;
            double const t = y_hi - mul<R>(static_cast<double>(k), 0.125);    // exact
            // a_0 + a_1 t + a_2 t² exacts (a_k en double-double, produits de Dekker) :
            // près des zéros, a_2 t² pèse jusqu'à 10 % du résultat. Le reste, t³ (…),
            // moins de 0,3 %.
            auto const& a = lgamma_table[k - lgamma_first];
            constexpr int last = 6 + lgamma_degree - 3;                     // a_26
            double h = a[last];
            for (int n = last - 1; n >= 6; --n)                             // a_3 + a_4 t + … + a_26 t^23
                h = a[n] + mul<R>(t, h);
            auto const [q1, q1_error] = two_product(a[2], t);
            auto const [tt, tt_error] = two_product(t, t);
            auto const [q2, q2_error] = two_product(a[4], tt);
            auto const [s1, e1] = two_sum(a[0], q1);
            auto const [s2, e2] = two_sum(s1, q2);
            double const tail = ((e1 + e2) + (a[1] + (q1_error + q2_error)))
                              + ((mul<R>(a[3], t) + (mul<R>(a[4], tt_error) + mul<R>(a[5], tt))) + mul<R>(a[2], y_lo))
                              + mul<R>(t * tt, h);
            auto const [hi, lo] = fast_two_sum(s2, tail);
            return {hi, lo};
        }

        // Stirling en double-double, y dans [10,0625 ; 2^900) :
        //     lgamma(y) = y log y − y − (log y)/2 + log(2π)/2 + Σ B_2k / (2k (2k − 1) y^(2k−1)),
        // y log y exact (produit de Dekker sur log y en double-double), les termes
        // de tête sommés exactement.
        template <reproducibility_mode R = portable_t>
        constexpr double_double lgamma_stirling(double y) noexcept
        {
            auto const [l_hi, l_lo] = log_dd<R>(y);
            auto const [p, p_error] = two_product(y, l_hi);
            auto const [s1, e1] = two_sum(p, -y);
            auto const [s2, e2] = two_sum(s1, unfused_mul(-0.5, l_hi));            // −(log y)/2 : exact
            auto const [s3, e3] = two_sum(s2, half_log_two_pi_hi);
            double const z = 1.0 / y;
            double const z2 = z * z;
            constexpr int terms = sizeof(stirling) / sizeof(stirling[0]);
            double series = stirling[terms - 1];
            for (int n = terms - 2; n >= 0; --n)
                series = stirling[n] + mul<R>(z2, series);
            double const tail = ((e1 + e2) + (e3 + p_error))
                              + ((mul<R>(y, l_lo) - mul<R>(0.5, l_lo)) + half_log_two_pi_lo) + mul<R>(z, series);
            auto const [hi, lo] = fast_two_sum(s3, tail);
            return {hi, lo};
        }

        // lgamma(y) en double-double, y > 0, 2^-54 <= y <= seuil de débordement.
        template <reproducibility_mode R = portable_t>
        constexpr double_double lgamma_positive(double y) noexcept
        {
            if (y < 0.5625)                                                  // lgamma(1 + y) − log y
            {
                auto const [u, u_lo] = two_sum(1.0, y);
                auto const [g_hi, g_lo] = lgamma_piece<R>(u, u_lo);
                auto const [l_hi, l_lo] = log_dd<R>(y);
                auto const [s, s_error] = two_sum(g_hi, -l_hi);
                auto const [hi, lo] = fast_two_sum(s, (s_error + g_lo) - l_lo);
                return {hi, lo};
            }
            if (y < 10.0625)
                return lgamma_piece<R>(y, 0.0);
            if (y < 0x1p900)
                return lgamma_stirling<R>(y);
            // y immense : lgamma(y) = y (log y − 1) à moins d'un ulp près, les autres
            // termes (quelques centaines) étant bien en dessous. Mise à l'échelle par
            // 2^-100 : le produit de Dekker n'admet pas de facteur au-delà de 2^995.
            auto const [l_hi, l_lo] = log_dd<R>(y);
            auto const [m, m_error] = two_sum(l_hi, -1.0);
            double const ys = mul<R>(y, 0x1p-100);                             // exact
            auto const [p, p_error] = two_product(ys, m);
            double const tail = p_error + mul<R>(ys, m_error + l_lo);
            return {mul<R>(p, 0x1p100), mul<R>(tail, 0x1p100)};
        }

        // |sin(π r)| en double-double, |r| <= 1/2, r ≠ 0 : π r en double-double,
        // puis sin(π |r|) si |r| <= 1/4, cos(π (1/2 − |r|)) sinon (exact : Sterbenz).
        template <reproducibility_mode R = portable_t>
        constexpr double_double sin_pi_abs(double r) noexcept
        {
            double const a = r < 0 ? -r : r;
            bool const near_half = a > 0.25;
            double const w = near_half ? 0.5 - a : a;
            auto const [p, p_error] = two_product(w, pi_hi);
            double const p_lo = p_error + mul<R>(w, pi_lo);
            auto const point = sincos_prepare<R>(p, p_lo);
            return near_half ? cos_of<R>(point) : sin_of<R>(point);
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
        template <reproducibility_mode R = portable_t>
        constexpr bool is_integer(double x) noexcept
        {
            return static_cast<double>(static_cast<std::int64_t>(x)) == x;
        }

        template <reproducibility_mode R = portable_t>
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
            auto const [s_hi, s_lo] = sin_pi_abs<R>(r);
            auto const [ls_hi, ls_lo0] = log_dd<R>(s_hi);
            double const ls_lo = ls_lo0 + s_lo / s_hi;

            auto const [y_hi, y_lo] = two_sum(1.0, -x);
            auto const [g_hi, g_lo0] = lgamma_positive<R>(y_hi);
            double const digamma = detail::log<R>(y_hi) - 0.5 / y_hi;
            double const g_lo = g_lo0 + mul<R>(digamma, y_lo);

            auto const [a, a_error] = two_sum(log_pi_hi, -ls_hi);
            auto const [b, b_error] = two_sum(a, -g_hi);
            auto const [hi, lo] = fast_two_sum(b, ((a_error + b_error) + log_pi_lo) - (ls_lo + g_lo));

            // Γ(x) < 0 exactement quand floor(x) est impair (r ≠ 0).
            auto const floor_x = static_cast<std::int64_t>(r > 0 ? n : n - 1);
            return {{hi, lo}, (floor_x & 1) != 0};
        }

        // ---------------------------------------------------------------------
        // digamma
        // ---------------------------------------------------------------------

        // 1/(a_hi + a_lo) en double-double, |1/a| < 2^990.
        template <reproducibility_mode R = portable_t>
        constexpr double_double reciprocal_dd(double a_hi, double a_lo) noexcept
        {
            double const q = 1.0 / a_hi;
            auto const [p, p_error] = two_product(q, a_hi);
            auto const [hi, lo] = fast_two_sum(q, (((1.0 - p) - p_error) - mul<R>(q, a_lo)) / a_hi);
            return {hi, lo};
        }

        template <reproducibility_mode R = portable_t>
        constexpr double_double add_dd(double_double a, double_double b) noexcept
        {
            auto const [s, s_error] = two_sum(a.hi, b.hi);
            auto const [hi, lo] = fast_two_sum(s, s_error + (a.lo + b.lo));
            return {hi, lo};
        }

        // ψ(y), y >= 10, en double-double : log y − 1/(2y) − Σ B_2k / (2k y^(2k)), la
        // série (moins de 1/1200 du tout) en double.
        template <reproducibility_mode R = portable_t>
        constexpr double_double digamma_asymptotic_dd(double y) noexcept
        {
            auto const [l_hi, l_lo] = log_dd<R>(y);
            double const h = 0.5 / y;
            double h_lo = 0;                                                // 1/(2y) = h + h_lo
            if (y < 0x1p900)                                                // au-delà, négligeable (et Dekker
            {                                                               // n'admet pas de facteur > 2^995)
                auto const [p, p_error] = two_product(h, y);
                h_lo = ((0.5 - p) - p_error) / y;
            }
            double const z = 1.0 / mul<R>(y, y);
            constexpr int terms = sizeof(digamma_asymptotic) / sizeof(digamma_asymptotic[0]);
            double series = digamma_asymptotic[terms - 1];
            for (int k = terms - 2; k >= 0; --k)
                series = digamma_asymptotic[k] + mul<R>(z, series);
            series = mul<R>(z, series);
            auto const [s, s_error] = two_sum(l_hi, -h);
            auto const [hi, lo] = fast_two_sum(s, s_error + ((l_lo - h_lo) - series));
            return {hi, lo};
        }

        // ψ(x), 2^-990 <= x < 2^1023, en double-double, précis en valeur absolue
        // (à 2^-63 près relativement à ψ(x + n)) : x >= 10 asymptotique ; sinon
        //     ψ(x) = ψ(s) − Σ_{i<n} 1/(x + i),   s = x + n dans [10, 11),
        // s en double-double (ψ(s_hi + s_lo) ≈ ψ(s_hi) + ψ'(s_hi) s_lo), la somme exacte à
        // 2^-100 près. Près du zéro x₀ la compensation coûte la précision relative :
        // digamma passe alors par le développement autour de x₀.
        template <reproducibility_mode R = portable_t>
        constexpr double_double digamma_positive_dd(double x) noexcept
        {
            if (x >= digamma_asymptotic_from)
                return digamma_asymptotic_dd<R>(x);
            int const n = 10 - static_cast<int>(x);
            auto const [s_hi, s_lo] = two_sum(x, static_cast<double>(n));
            double_double result = digamma_asymptotic_dd<R>(s_hi);
            double const w = 1.0 / s_hi;                                    // ψ'(s) ≈ 1/s + 1/(2s²) + 1/(6s³)
            result = add_dd<R>(result, {mul<R>(s_lo, w + mul<R>(mul<R>(w, w), 0.5 + mul<R>(w, 0x1.5555555555555p-3))), 0.0});
            for (int i = n - 1; i >= 0; --i)                                // les plus petits termes d'abord
            {
                auto const [f_hi, f_lo] = two_sum(x, static_cast<double>(i));
                auto const [r_hi, r_lo] = reciprocal_dd<R>(f_hi, f_lo);
                result = add_dd<R>(result, {-r_hi, -r_lo});
            }
            return result;
        }

        // ψ(x₀ + t), |t| <= 1/8 : (c_1 + q t) t, q = c_2 + t (c_3 + …), c_1 et c_2 en
        // double-double, t = x − x₀ à 2^-160 près.
        template <reproducibility_mode R = portable_t>
        constexpr double digamma_near_root(double x) noexcept
        {
            auto const [t_hi, t_lo0] = two_sum(x - digamma_root_hi, -digamma_root_lo);   // x − x₀_hi : exact
            double const t_lo = t_lo0 - digamma_root_lo2;
            constexpr auto const& c = digamma_root_taylor;
            constexpr int last = digamma_root_degree + 1;                   // indice de c_24
            double h = c[last];
            for (int k = last - 1; k >= 4; --k)                             // c_3 + c_4 t + …
                h = c[k] + mul<R>(t_hi, h);
            auto const [q_hi, q_lo0] = two_sum(c[2], unfused_mul(t_hi, h));
            double const q_lo = q_lo0 + c[3];
            auto const [m, m_error] = two_product(q_hi, t_hi);
            auto const [u_hi, u_lo0] = two_sum(c[0], m);
            double const u_lo = u_lo0 + ((c[1] + m_error) + (mul<R>(q_lo, t_hi) + mul<R>(q_hi, t_lo)));
            auto const [v, v_error] = two_product(u_hi, t_hi);
            return v + (v_error + (mul<R>(u_lo, t_hi) + mul<R>(u_hi, t_lo)));
        }

        // ψ(x), x < 0 non entier, −2^52 < x <= −2^-54, en double-double :
        //     ψ(x) = ψ(y) + 1/y + π cot(π y),   y = −x
        // (réflexion ψ(1 − x) − π cot(π x), avec ψ(1 + y) = ψ(y) + 1/y). y = n + r,
        // |r| <= 1/2 exact ; cot(π r) = cos(π r)/sin(π r), les deux en double-double.
        template <reproducibility_mode R = portable_t>
        constexpr double_double digamma_reflection_dd(double x) noexcept
        {
            double const y = -x;
            double n = static_cast<double>(static_cast<std::int64_t>(y));
            double r = y - n;                                                // exact, dans [0, 1)
            if (r > 0.5)
                r -= 1.0;                                                    // exact
            auto const [s_hi, s_lo] = sin_pi_abs<R>(r);                        // |sin(π r)|
            auto const [c_hi, c_lo] = sin_pi_abs<R>(0.5 - (r < 0 ? -r : r));   // cos(π r) >= 0
            auto const [m, m_error] = two_product(pi_hi, c_hi);
            double const m_lo = m_error + (mul<R>(pi_hi, c_lo) + mul<R>(pi_lo, c_hi));
            double const q = m / s_hi;
            auto const [p, p_error] = two_product(q, s_hi);
            double const q_lo = (((m - p) - p_error) + (m_lo - mul<R>(q, s_lo))) / s_hi;
            double_double const cot = r < 0 ? double_double{-q, -q_lo} : double_double{q, q_lo};
            double_double result = digamma_positive_dd<R>(y);
            result = add_dd<R>(result, reciprocal_dd<R>(y, 0.0));
            return add_dd<R>(result, cot);
        }

        template <reproducibility_mode R = portable_t>
        constexpr double exp_dd(double_double a) noexcept
        {
            auto const [k, r_hi, r_lo] = reduce_exp<R>(a.hi);
            return exp_finish_dd<R>(k, r_hi, r_lo + a.lo);
        }
    }

    namespace detail
    {
    template <reproducibility_mode R>
    constexpr double lgamma(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return detail::opaque(x + x);
        if (x == infinity || x == -infinity || x == 0)
            return infinity;
        double const ax = x < 0 ? -x : x;
        if (ax < 0x1p-54)                                                    // −log|x| − γx, γx sous le demi-ulp
            return -detail::log<R>(ax);
        if (x > 0)
        {
            if (x > 0x1.754d9278b51a7p+1014)                                 // lgamma(x) > DBL_MAX
                return infinity;
            if (x == 1.0 || x == 2.0)
                return 0.0;
            auto const [hi, lo] = lgamma_positive<R>(x);
            return hi + lo;
        }
        if (x <= -0x1p52 || is_integer<R>(x))                                   // entier négatif : pôle
            return infinity;
        auto const [value, negative] = lgamma_reflection<R>(x);
        return value.hi + value.lo;
    }
    }

    // ψ(x) = Γ'(x)/Γ(x). Erreur mesurée sous 0,51 ulp pour x > 0 ; pour x < 0, près des
    // zéros de ψ (un par intervalle (−n − 1, −n)), la réflexion retranche des termes
    // presque égaux : l'erreur y est petite en valeur absolue, pas en ulp.
    namespace detail
    {
    template <reproducibility_mode R>
    constexpr double digamma(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return opaque(x + x);
        if (x == infinity)
            return x;
        if (x == 0)
            return (bits(x) >> 63) != 0 ? infinity : -infinity;
        if (x > 0)
        {
            if (x < 0x1p-990)                                                // −1/x − γ, γ x sous 2^-990
                return -1.0 / x;
            double const t = x - digamma_root_hi;
            if (t > -digamma_root_half_width && t < digamma_root_half_width)
                return digamma_near_root<R>(x);
            auto const [hi, lo] = digamma_positive_dd<R>(x);
            return hi + lo;
        }
        if (x <= -0x1p52 || x == -infinity || is_integer<R>(x))                // pôles (−∞ : NaN)
            return quiet_nan;
        if (x > -0x1p-54)                                                    // −1/x − γ + O(x)
        {
            if (x > -0x1p-990)
                return -1.0 / x;
            auto const [r_hi, r_lo] = reciprocal_dd<R>(x, 0.0);
            auto const [s, s_error] = two_sum(-r_hi, -euler_gamma_hi);
            return s + (s_error - (r_lo + euler_gamma_lo));
        }
        auto const [hi, lo] = digamma_reflection_dd<R>(x);
        return hi + lo;
    }
    }

    namespace detail
    {
    template <reproducibility_mode R>
    constexpr double tgamma(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return detail::opaque(x + x);
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
            return exp_dd<R>(lgamma_positive<R>(x));
        }
        if (x <= -0x1p52 || is_integer<R>(x))                                   // entier négatif
            return quiet_nan;
        auto const [value, negative] = lgamma_reflection<R>(x);
        double const magnitude = value.hi < -800.0 ? 0.0 : exp_dd<R>(value);  // |Γ(x)| < 2^-1150 : 0
        return negative ? -magnitude : magnitude;
    }
    }

    // float : calcul en double, puis un seul arrondi vers float.


    namespace detail
    {
        constexpr float tgamma_to_float(double y) noexcept
        {
            if (y <= -0x1.ffffffp+127)                                       // < −FLT_MAX − ulp/2 : −∞
                return static_cast<float>(-infinity);
            return to_float(y);
        }
    }

    // float : calcul en double, puis un seul arrondi vers float (tgamma : −∞ sous −FLT_MAX).
    WSL_SF_FUNCTION_1(lgamma, detail::to_float)
    WSL_SF_FUNCTION_1(digamma, detail::to_float)
    WSL_SF_FUNCTION_1(tgamma, detail::tgamma_to_float)

    long double lgamma(long double) = delete;
    long double tgamma(long double) = delete;

    long double digamma(long double) = delete;
}
