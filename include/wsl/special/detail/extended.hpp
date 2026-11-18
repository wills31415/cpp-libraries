#pragma once

// Arithmétique étendue (double-double à exposant séparé, exp, expm1, log, log1p à
// argument double-double…), commune aux fonctions incomplètes et aux lois de
// wsl::random. Arithmétique soignée : calculs en double-double arrondis une
// seule fois, mis à l'échelle pour rester dans le domaine où les produits exacts
// le sont (facteurs sous 2^995, produits loin des sous-normaux : voir
// wsl/special/detail/fp.hpp). Bit-exact comme le reste.

#include <wsl/special/algebraic.hpp>
#include <wsl/special/detail/fp.hpp>
#include <wsl/special/elementary.hpp>
#include <wsl/special/error_function.hpp>
#include <wsl/special/gamma.hpp>

#include <bit>
#include <cstdint>
#include <limits>

namespace wsl::special::detail::extended
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

    // ---------------------------------------------------------------------
    // Constantes en double-double (valeurs de mpmath, vérifiées par
    // tools/random/law_references.py).
    // ---------------------------------------------------------------------
    inline constexpr double inv_sqrt2_hi = 0x1.6a09e667f3bcdp-1, inv_sqrt2_lo = -0x1.bdd3413b26456p-55;
    inline constexpr double sqrt2_hi = 0x1.6a09e667f3bcdp+0, sqrt2_lo = -0x1.bdd3413b26456p-54;
    inline constexpr double inv_sqrt_2pi_hi = 0x1.9884533d43651p-2, inv_sqrt_2pi_lo = -0x1.cbc0d30ebfd15p-56;
    inline constexpr double half_log_2pi_hi = 0x1.d67f1c864beb5p-1, half_log_2pi_lo = -0x1.65b5a1b7ff5dfp-55;
    inline constexpr double normal_entropy_hi = 0x1.6b3f8e4325f5ap+0, normal_entropy_lo = 0x1.4d252f2400510p-54;  // (1 + log 2π)/2
    inline constexpr double half_sqrt_pi = 0x1.c5bf891b4ef6bp-1;
    inline constexpr double inv_sqrt_pi = 0x1.20dd750429b6dp-1;
    inline constexpr double inv_pi_hi = 0x1.45f306dc9c883p-2, inv_pi_lo = -0x1.6b01ec5417056p-56;
    inline constexpr double log_4pi_hi = 0x1.43f89a3f0edd6p+1, log_4pi_lo = 0x1.053cd734e6a31p-54;
    inline constexpr double one_twelfth_hi = 0x1.5555555555555p-4, one_twelfth_lo = 0x1.5555555555555p-58;

    constexpr double magnitude(double x) noexcept { return x < 0 ? -x : x; }

    // a·b exactement (produit et erreur), |a·b| dans [2^-900, 2^900] : les
    // facteurs sont rééquilibrés (exactement) pour que le produit de Dekker le soit.
    constexpr sd::sum_error balanced_product(double a, double b) noexcept
    {
        int const shift = (sd::exponent_of(magnitude(a)) - sd::exponent_of(magnitude(b))) / 2;
        return sd::two_product(sd::scale_by_power_of_two(a, -shift), sd::scale_by_power_of_two(b, shift));
    }

    // a·b en double-double, a, b > 0 finis : exact si |ab| est dans [2^-900,
    // 2^900] ; arrondi sinon (produit minuscule ou immense, sans erreur rendue).
    constexpr sd::double_double product(double a, double b) noexcept
    {
        double const approx = magnitude(sd::fmul(a, b));
        if (approx >= 0x1p-900 && approx <= 0x1p900)
        {
            auto const [p, p_error] = balanced_product(a, b);
            return {p, p_error};
        }
        return {sd::fmul(a, b), 0.0};
    }

    // Entier de 64 bits en somme exacte de deux doubles (le second nul sous 2^53).
    constexpr sd::double_double exact_split(std::uint64_t v) noexcept
    {
        if (v < (std::uint64_t{1} << 53))
            return {static_cast<double>(v), 0.0};
        return {static_cast<double>(v & ~std::uint64_t{0x7ff}), static_cast<double>(v & 0x7ff)};
    }

    // Produit et somme en double-double (produit exact dans [2^-900, 2^900]).
    constexpr sd::double_double mul(sd::double_double a, sd::double_double b) noexcept
    {
        bool const negative = (a.hi < 0) != (b.hi < 0);
        auto const [p, p_error] = product(magnitude(a.hi), magnitude(b.hi));
        double const ph = negative ? -p : p, pe = negative ? -p_error : p_error;
        auto const [hi, lo] = sd::fast_two_sum(ph, pe + (sd::fmul(a.hi, b.lo) + sd::fmul(a.lo, b.hi)));
        return {hi, lo};
    }

    constexpr sd::double_double sum(sd::double_double a, sd::double_double b) noexcept
    {
        auto const [s, s_error] = sd::two_sum(a.hi, b.hi);
        auto const [hi, lo] = sd::fast_two_sum(s, s_error + (a.lo + b.lo));
        return {hi, lo};
    }

    // 1/n en double-double, n entier > 0 exact.
    constexpr sd::double_double reciprocal(double n) noexcept
    {
        double const q = 1.0 / n;
        auto const [p, p_error] = sd::two_product(q, n);
        return {q, ((1.0 - p) - p_error) / n};
    }

    // (k_hi + k_lo)(l_hi + l_lo), en double-double : k entier exact, l quelconque.
    constexpr sd::double_double times(sd::double_double k, sd::double_double l) noexcept
    {
        auto const [p, p_error] = product(k.hi, magnitude(l.hi));
        double const signed_p = l.hi < 0 ? -p : p, signed_error = l.hi < 0 ? -p_error : p_error;
        auto const [hi, lo] = sd::fast_two_sum(signed_p, signed_error + (sd::fmul(k.hi, l.lo) + sd::fmul(k.lo, l.hi)));
        return {hi, lo};
    }

    // ⌊x⌋, 0 <= x < 2^63.
    constexpr double floor_nonnegative(double x) noexcept
    {
        return static_cast<double>(static_cast<std::uint64_t>(x));
    }

    // ⌈hi + lo⌉, 0 <= hi + lo < 2^63, en entier.
    constexpr std::uint64_t ceiling(double hi, double lo) noexcept
    {
        auto const n = static_cast<std::uint64_t>(hi);                      // ⌊hi⌋, exact en double
        if (static_cast<double>(n) < hi)                                    // hi non entier : lo ne fait
            return n + 1;                                                   // franchir aucun entier
        // hi entier : ⌈hi + lo⌉ = hi + ⌈lo⌉ (lo peut dépasser 1 quand hi >= 2^54).
        auto const m = static_cast<std::int64_t>(lo);                       // troncature vers 0
        std::int64_t const c = (lo > 0 && static_cast<double>(m) < lo) ? m + 1 : m;
        return static_cast<std::uint64_t>(static_cast<std::int64_t>(n) + c);
    }

    // log(x) en double-double, x > 0 fini (sous-normal compris), à 2^-72 près
    // relativement.
    constexpr sd::double_double log_precise(double x) noexcept
    {
        return sd::log_dd(x);
    }

    // ---------------------------------------------------------------------
    // Double-double à exposant séparé : (hi + lo) · 2^e, hi dans [1, 2). Ni
    // débordement ni sous-normal dans les calculs intermédiaires ; un seul arrondi
    // à la fin (rounded), sous-normaux compris.
    // ---------------------------------------------------------------------
    struct scaled
    {
        double hi, lo;
        int e;
    };

    // hi > 0.
    constexpr scaled normalized(double hi, double lo, int e = 0) noexcept
    {
        int const eh = sd::exponent_of(hi);
        return {sd::scale_by_power_of_two(hi, -eh), sd::scale_by_power_of_two(lo, -eh), e + eh};
    }

    constexpr scaled operator*(scaled const& a, scaled const& b) noexcept
    {
        auto const [p, p_error] = sd::two_product(a.hi, b.hi);             // dans [1, 4) : exact
        auto const [hi, lo] = sd::fast_two_sum(p, p_error + (sd::fmul(a.hi, b.lo) + sd::fmul(a.lo, b.hi)));
        return normalized(hi, lo, a.e + b.e);
    }

    // a / d, d > 0 fini.
    constexpr scaled operator/(scaled const& a, double d) noexcept
    {
        int const ed = sd::exponent_of(d);
        double const m = sd::scale_by_power_of_two(d, -ed);                 // dans [1, 2)
        double const q = a.hi / m;
        auto const [p, p_error] = sd::two_product(q, m);
        auto const [hi, lo] = sd::fast_two_sum(q, (((a.hi - p) - p_error) + a.lo) / m);
        return normalized(hi, lo, a.e - ed);
    }

    // a / b.
    constexpr scaled operator/(scaled const& a, scaled const& b) noexcept
    {
        double const q = a.hi / b.hi;                                       // dans (1/2, 2)
        auto const [p, p_error] = sd::two_product(q, b.hi);
        auto const [hi, lo] = sd::fast_two_sum(q, ((((a.hi - p) - p_error) + a.lo) - sd::fmul(q, b.lo)) / b.hi);
        return normalized(hi, lo, a.e - b.e);
    }

    // √a.
    constexpr scaled sqrt(scaled a) noexcept
    {
        if (a.e % 2 != 0)                                                   // exposant pair
        {
            a.hi = sd::fmul(a.hi, 2.0);
            a.lo = sd::fmul(a.lo, 2.0);
            --a.e;
        }
        double const s = wsl::special::sqrt(a.hi);
        auto const [p, p_error] = sd::two_product(s, s);
        auto const [hi, lo] = sd::fast_two_sum(s, (((a.hi - p) - p_error) + a.lo) / (s + s));
        return normalized(hi, lo, a.e / 2);
    }

    // log a, en double-double.
    constexpr sd::double_double log(scaled const& a) noexcept
    {
        auto const [l_hi, l_lo] = log_precise(a.hi);
        double const ed = static_cast<double>(a.e);
        auto const [m, m_error] = sd::two_product(ed, sd::ln2_hi);
        auto const [s, s_error] = sd::two_sum(m, l_hi);
        auto const [hi, lo] = sd::fast_two_sum(s, (s_error + m_error) + ((l_lo + a.lo / a.hi) + sd::fmul(ed, sd::ln2_lo)));
        return {hi, lo};
    }

    // Au plus proche : un arrondi, sous-normaux compris.
    constexpr double rounded(scaled const& a) noexcept
    {
        if (a.e > 1023)
            return std::numeric_limits<double>::infinity();
        if (a.e < -1076)                                                    // < 2^-1075 : vers 0
            return 0.0;
        return sd::scale_exp(a.hi, a.lo, a.e);
    }

    // En double-double ordinaire (exact tant que le résultat est normal).
    constexpr sd::double_double unscaled(scaled const& a) noexcept
    {
        return {sd::scale_by_power_of_two(a.hi, a.e), sd::scale_by_power_of_two(a.lo, a.e)};
    }

    // ---------------------------------------------------------------------
    // Fonctions à argument double-double.
    // ---------------------------------------------------------------------

    // e^x, x = x_hi + x_lo dans [−1500, 1400], à 2^-70 près relativement, à
    // l'échelle. Au-dessous de −1400 (hors de la réduction de Cody–Waite),
    // e^x = 2^-1024 e^{x + 1024 ln 2}. Réduction de wsl::special (x = k ln2/128 + r),
    // puis e^r = 1 + r_hi + q, q = t + r t + r²/2 + … + r⁷/5040 (t = r_lo + x_lo,
    // reste < 2^-80) et T (1 + r_hi + q) avec T r_hi exact : le noyau de
    // wsl::special::exp, à 2^-60 près, ne suffit pas pour rester sous 0,501 ulp
    // après les produits qui suivent.
    constexpr scaled exp_scaled(double x_hi, double x_lo) noexcept
    {
        int extra = 0;
        if (x_hi < -1400)
        {
            auto const [s, s_error] = sd::two_sum(x_hi, sd::fmul(sd::ln2_hi, 1024.0));
            auto const [hi, lo] = sd::fast_two_sum(s, (s_error + x_lo) + sd::fmul(sd::ln2_lo, 1024.0));
            x_hi = hi;
            x_lo = lo;
            extra = -1024;
        }
        auto const [k, r, r_lo] = sd::reduce_exp(x_hi);
        double const t = r_lo + x_lo;
        double const r2 = sd::fmul(r, r);
        double const q = (t + sd::fmul(r, t))
                       + sd::fmul(r2, 0.5 + sd::fmul(r, 0x1.5555555555555p-3 + sd::fmul(r, 0x1.5555555555555p-5
                                  + sd::fmul(r, 0x1.1111111111111p-7 + sd::fmul(r, 0x1.6c16c16c16c17p-10 + sd::fmul(r, 0x1.a01a01a01a01ap-13))))));
        int const j = sd::exp_index(k);
        double const t_hi = sd::exp_table[j][0], t_lo = sd::exp_table[j][1];
        auto const [m, m_error] = sd::two_product(t_hi, r);
        auto const [s, s_error] = sd::fast_two_sum(t_hi, m);               // |T| > |T r|
        double const tail = s_error + (m_error + ((t_lo + sd::fmul(t_hi, q)) + sd::fmul(t_lo, r)));
        return normalized(s, tail, sd::exp_exponent(k) + extra);
    }

    // e^x − 1, x = x_hi + x_lo dans [−40, 1], en double-double : la méthode de
    // wsl::special::expm1, l'argument réduit gardé en double-double.
    constexpr sd::double_double expm1_double_double(double x_hi, double x_lo) noexcept
    {
        if (magnitude(x_hi) < 0x1p-54)                                      // e^x − 1 = x + x²/2 + …
        {
            if (magnitude(x_hi) < 0x1p-500)                                 // x² sous 2^-1000 : sans effet
            {
                auto const [h, l] = sd::fast_two_sum(x_hi, x_lo);
                return {h, l};
            }
            auto const [h, l] = sd::fast_two_sum(x_hi, x_lo + sd::fmul(sd::fmul(x_hi, x_hi), 0.5));
            return {h, l};
        }
        auto const [k, a_hi, a_lo] = sd::reduce_exp(x_hi);
        auto const [r_hi, r_lo] = sd::two_sum(a_hi, a_lo + x_lo);
        int const j = sd::exp_index(k);
        int const e = sd::exp_exponent(k);
        double const t_hi = sd::exp_table[j][0];
        double const t_lo = sd::exp_table[j][1];
        double const r = r_hi;
        double const rr = r * r;
        double const a = 0x1p-1 + sd::fmul(r, 0x1.5555555555555p-3);
        double const b = 0x1.5555555555555p-5 + sd::fmul(r, 0x1.1111111111111p-7);
        double const c = 0x1.6c16c16c16c17p-10 + sd::fmul(r, 0x1.a01a01a01a01ap-13);
        double const p_lo = r_lo + sd::fmul(rr, (a + sd::fmul(rr, b)) + sd::fmul(rr * rr, c));
        double const scale = sd::power_of_two(e);
        auto const [s, s_error] = sd::two_sum(sd::fmul(scale, t_hi), -1.0);
        auto const [m, m_error] = sd::two_product(t_hi, r_hi);
        auto const [s2, s2_error] = sd::two_sum(s, sd::fmul(scale, m));
        double const rest = sd::fmul(scale, (m_error + sd::fmul(t_hi, p_lo)) + sd::fmul(t_lo, 1.0 + (r_hi + p_lo)));
        auto const [h, l] = sd::fast_two_sum(s2, (s2_error + s_error) + rest);
        return {h, l};
    }

    // e^x − 1 arrondi une fois.
    constexpr double expm1_dd(double x_hi, double x_lo) noexcept
    {
        auto const [hi, lo] = expm1_double_double(x_hi, x_lo);
        return hi + lo;
    }

    // log(1 + r) en double-double, r = r_hi + r_lo >= 0, à 2^-69 près
    // relativement.
    constexpr sd::double_double log1p_positive(double r_hi, double r_lo) noexcept
    {
        if (r_hi == 0)
            return {0.0, 0.0};
        if (r_hi < 0x1p-10)
        {
            auto const [h, l] = sd::fast_two_sum(r_hi, sd::log1p_minus_identity(r_hi) + r_lo / (1.0 + r_hi));
            return {h, l};
        }
        auto const [u, u_error] = sd::two_sum(1.0, r_hi);
        auto const [hi, lo] = log_precise(u);
        auto const [h, l] = sd::fast_two_sum(hi, lo + (u_error + r_lo) / u);
        return {h, l};
    }

    // log(1 + x) en double-double, x dans (−1, 0], à 2^-69 près relativement.
    constexpr sd::double_double log1p_dd(double x) noexcept
    {
        if (x == 0)
            return {0.0, 0.0};
        if (x > -0x1p-10)
        {
            auto const [h, l] = sd::fast_two_sum(x, sd::log1p_minus_identity(x));
            return {h, l};
        }
        auto const [u, u_error] = sd::two_sum(1.0, x);
        auto const [hi, lo] = sd::log_accurate(sd::bits(u));
        auto const [h, l] = sd::fast_two_sum(hi, lo + u_error / u);
        return {h, l};
    }
}
