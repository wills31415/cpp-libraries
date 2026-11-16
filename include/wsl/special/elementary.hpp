#pragma once

// Fonctions élémentaires bit-exactes : même résultat au bit près sur toute
// plateforme (voir detail/fp.hpp pour les conditions), constexpr.
//
//   exp, exp2, expm1     log, log2, log10, log1p
//
// En double, erreur < 0,51 ulp ; en float, calcul en double puis un seul arrondi
// (erreur < 0,5 + 2^-28 ulp). Les erreurs maximales mesurées contre mpmath sont
// vérifiées par tests/test_special_elementary.cpp. Cas limites ceux de <cmath>
// (Annexe F de C) : exp(±0) = 1, exp(+∞) = +∞, exp(−∞) = +0, expm1(±0) = ±0,
// expm1(−∞) = −1 ; log(1) = +0, log(±0) = −∞, log(+∞) = +∞, log(x < 0) = NaN,
// log1p(±0) = ±0, log1p(−1) = −∞ ; NaN donne NaN.
//
// Règle d'écriture : tout produit dont le résultat entre dans une addition ou une
// soustraction passe par fmul (voir detail/fp.hpp).

#include <wsl/special/detail/fp.hpp>
#include <wsl/special/detail/tables.hpp>

#include <concepts>
#include <cstdint>

namespace wsl::special
{
    namespace detail
    {
        // -------------------------------------------------------------------
        // Noyaux de la famille exponentielle.
        //
        // x = k ln2/N + r, |r| <= ln2/(2N), N = 128, k = N e + j :
        //     exp(x) = 2^e · 2^(j/N) · exp(r).
        // 2^(j/N) vient d'une table en double-double (T_hi + T_lo). Réduction de
        // Cody–Waite : k ln2/N est exact (ln2/N sur 35 bits, |k| < 2^18), et
        // x − k ln2_hi aussi (lemme de Sterbenz). r est rendu en double-double
        // (r_hi + r_lo), dont seul expm1 a besoin.
        // -------------------------------------------------------------------
        struct exp_reduction
        {
            std::int64_t k;
            double r_hi, r_lo;
        };

        constexpr exp_reduction reduce_exp(double x) noexcept
        {
            double const z = fmul(x, exp_inv_ln2_n);
            auto const k = static_cast<std::int64_t>(z < 0 ? z - 0.5 : z + 0.5);
            double const kd = static_cast<double>(k);
            double const a = x - fmul(kd, exp_ln2_n_hi);            // exact
            double const b = fmul(kd, exp_ln2_n_lo);
            double const r_hi = a - b;
            return {k, r_hi, (a - r_hi) - b};                       // Fast2Sum : |a| >= |b|
        }

        constexpr int exp_index(std::int64_t k) noexcept
        {
            return static_cast<int>(k & ((std::int64_t{1} << exp_table_bits) - 1));
        }

        constexpr int exp_exponent(std::int64_t k) noexcept
        {
            return static_cast<int>(k >> exp_table_bits);           // division entière par défaut
        }

        // 2^e · (t_hi + tail), avec t_hi + tail dans [0,99 ; 2) : un seul arrondi,
        // résultats sous-normaux compris.
        constexpr double scale_exp(double t_hi, double tail, int e) noexcept
        {
            if (e >= -1022)
            {
                double const y = t_hi + tail;
                if (e > 1023)                                       // résultat juste sous DBL_MAX
                    return y * power_of_two(1023) * 2.0;
                return y * power_of_two(e);
            }

            // Résultat sous-normal. Arrondir y à 53 bits puis le ramener dans la
            // zone sous-normale, de précision moindre, arrondirait deux fois
            // (jusqu'à 0,75 ulp). On travaille donc sur y = résultat · 2^1022, dans
            // (0, 1) : 1 + y arrondit au bit voulu (2^-52, soit 2^-1074 une fois
            // remis à l'échelle), et l'erreur de chaque étape est reportée, si bien
            // que le résultat n'est arrondi qu'une fois. La remise à l'échelle
            // finale est exacte.
            double const scale = power_of_two(e + 1022);             // e + 1022 dans [-58, -1]
            double const a = fmul(scale, t_hi);                     // exacts : puissance de 2
            double const b = fmul(scale, tail);
            double const y = a + b;
            double const y_error = (a - y) + b;                     // exacte, |a| >= |b|
            double const high = 1.0 + y;
            double const low = ((1.0 - high) + y) + y_error;
            return ((high + low) - 1.0) * 0x1p-1022;
        }

        // 2^e 2^(j/N) exp(r) pour k = N e + j. exp(r) − 1 : Taylor de degré 5
        // (reste < 2^-60 sur l'intervalle). T_hi + (T_lo + T_hi p) n'est arrondi
        // qu'une fois, d'où une erreur à peine supérieure à 1/2 ulp.
        constexpr double exp_finish(std::int64_t k, double r) noexcept
        {
            // exp(r) − 1 = r + r²/2 + r³/6 + r⁴/24 + r⁵/120
            double const r2 = r * r;
            double const p = r + fmul(r2, 0x1p-1 + fmul(r, 0x1.5555555555555p-3
                                       + fmul(r, 0x1.5555555555555p-5 + fmul(r, 0x1.1111111111111p-7))));
            int const j = exp_index(k);
            double const t_hi = exp_table[j][0];
            return scale_exp(t_hi, exp_table[j][1] + fmul(t_hi, p), exp_exponent(k));
        }

        // -------------------------------------------------------------------
        // Noyau de la famille logarithme.
        //
        // x = 2^e z avec z dans [√2/2, √2), puis z = c (1 + r) avec c = j/512 :
        //     log(x) = e ln2 + log(c) + log1p(r),  |r| < 2^-9,5.
        //
        // - Réduction sans branche, en arithmétique entière sur les bits de x : z
        //   garde la mantisse de x, e est l'exposant décalé de sorte que z tombe
        //   dans [√2/2, √2). Près de 1, e = 0 : pas d'annulation entre e ln2 et
        //   log(z).
        // - j = round(512 z) par l'astuce du décalage (2^52 · 1,5) : sans
        //   conversion.
        // - z − c est exact (Sterbenz). r est obtenu en double-double sans
        //   division : r_hi = (z − c)·(1/c), 1/c tabulé, puis le résidu
        //   z − c − r_hi·c, exact par un produit de Dekker (c n'a que 10 bits),
        //   donne r_lo.
        // - log(c) est tabulé en deux parties, hi multiple de 2^-42 comme e ln2_hi :
        //   w = e ln2_hi + log(c)_hi est exact. Puis |w| >= |r| (ou w = 0 quand
        //   c = 1), d'où une somme exacte w + r_hi en trois opérations (Fast2Sum).
        // - log1p(r) − r : Taylor de degré 7 (reste < 2^-69 relativement à r),
        //   évalué par le schéma d'Estrin, dont les termes se calculent en
        //   parallèle.
        // Le logarithme est rendu en double-double (hi + lo) : log l'arrondit en
        // une addition, log2, log10 et log1p le prolongent.
        // -------------------------------------------------------------------

        // log1p(r) − r pour |r| < 2^-9,5 : Taylor de degré 7, Estrin.
        //     r² [(−1/2 + r/3) + r² (−1/4 + r/5) + r⁴ (−1/6 + r/7)]
        constexpr double log1p_minus_identity(double r) noexcept
        {
            double const rr = r * r;
            double const p01 = -0x1p-1 + fmul(r, 0x1.5555555555555p-2);
            double const p23 = -0x1p-2 + fmul(r, 0x1.999999999999ap-3);
            double const p45 = -0x1.5555555555555p-3 + fmul(r, 0x1.2492492492492p-3);
            return fmul(rr, (p01 + fmul(rr, p23)) + fmul(rr * rr, p45));
        }

        // Cas à part : renvoie vrai, avec le résultat, si x n'est pas un réel > 0
        // fini ; sinon met dans ix les bits de x, remis à l'échelle s'il est
        // sous-normal. Un seul test pour tous les cas : x <= 0 (le bit de signe
        // rend les bits très grands), sous-normal, +∞ ou NaN.
        constexpr bool log_special_case(double x, std::uint64_t& ix, double& result) noexcept
        {
            ix = bits(x);
            if (ix - 0x0010000000000000 >= 0x7fe0000000000000) [[unlikely]]
            {
                if (!(x == x))
                    result = x + x;                                 // NaN
                else if (x == 0)
                    result = -infinity;
                else if (x < 0)
                    result = quiet_nan;
                else if (x == infinity)
                    result = x;
                else
                {
                    ix = bits(x * 0x1p52) - (std::uint64_t{52} << 52);   // sous-normal : x · 2^52
                    return false;
                }
                return true;
            }
            return false;
        }

        constexpr double_double log_kernel(std::uint64_t ix) noexcept
        {
            constexpr std::uint64_t half_sqrt2 = 0x3fe6a09e667f3bcd;   // bits de √2/2
            std::uint64_t const offset = ix - half_sqrt2;
            auto const e = static_cast<int>(static_cast<std::int64_t>(offset) >> 52);
            double const z = from_bits(ix - (offset & (std::uint64_t{0xfff} << 52)));   // dans [√2/2, √2)

            constexpr double shift = 0x1.8p52;                      // ulp 1 : t − shift = round(512 z)
            constexpr double scale = static_cast<double>(1 << log_table_bits);
            double const t = fmul(z, scale) + shift;                // produit exact
            auto const j = static_cast<int>(bits(t) - bits(shift)); // dans [362, 724]
            double const c = fmul(t - shift, 1.0 / scale);          // exact, mais entre dans z − c
            double const d = z - c;                                 // exact

            auto const& entry = log_table[j - log_table_first];
            double const inverse_c = entry[2];
            double const r_hi = fmul(d, inverse_c);
            auto const [r1, r2] = split(r_hi);
            double const r_lo = fmul((d - fmul(r1, c)) - fmul(r2, c), inverse_c);

            double const q = log1p_minus_identity(r_hi);

            double const ed = static_cast<double>(e);
            double const w = fmul(ed, log_ln2_hi) + entry[0];      // exact
            double const s = w + r_hi;
            double const s_error = (w - s) + r_hi;                  // exact : |w| >= |r_hi| ou w = 0
            double const tail = (fmul(ed, log_ln2_lo) + entry[1]) + r_lo;
            return {s, (s_error + tail) + q};
        }

        // log(x) · (c_hi + c_lo), le produit en double-double, arrondi une fois.
        constexpr double log_times(double x, double c_hi, double c_lo) noexcept
        {
            std::uint64_t ix = 0;
            double special = 0;
            if (log_special_case(x, ix, special))
                return special;                                     // ±∞, NaN : inchangés (c > 0)
            auto const [hi, lo] = log_kernel(ix);
            auto const [p, p_error] = two_product(hi, c_hi);
            return p + ((p_error + fmul(hi, c_lo)) + fmul(lo, c_hi));
        }
    }

    // -----------------------------------------------------------------------
    // exp, exp2, expm1
    // -----------------------------------------------------------------------
    constexpr double exp(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return x + x;                                           // NaN
        if (x > 0x1.62e42fefa39efp+9)                               // > ln(DBL_MAX) : +∞
            return infinity;
        if (x < -746.0)                                             // exp(x) < 2^-1076 : +0
            return 0.0;
        auto const [k, r, r_lo] = reduce_exp(x);
        return exp_finish(k, r);
    }

    // 2^x : x = k/N + f, réduction exacte (N = 128 est une puissance de 2), puis
    // 2^f = exp(f ln2), |f ln2| <= ln2/(2N) : le noyau de exp. 2^n est exact pour
    // tout entier n.
    constexpr double exp2(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return x + x;
        if (x >= 1024.0)
            return infinity;
        if (x < -1080.0)                                            // 2^x < 2^-1076 : +0
            return 0.0;
        constexpr double shift = 0x1.8p52;
        double const t = fmul(x, 128.0) + shift;                    // produit exact
        auto const k = static_cast<std::int64_t>(bits(t) - bits(shift));   // round(128 x)
        double const f = x - fmul(t - shift, 0x1p-7);               // exact (Sterbenz)
        return exp_finish(k, fmul(f, ln2_hi));
    }

    // e^x − 1, précis aussi quand le résultat est petit. Même réduction que exp,
    // r gardé en double-double, exp(r) − 1 en Taylor de degré 7 (reste < 2^-84) :
    //     2^e T (1 + p) − 1 = (2^e T_hi − 1) + 2^e T_hi r_hi + 2^e [T_hi p_lo + T_lo (1 + p)],
    // les deux premiers termes ajoutés exactement (two_sum, produit de Dekker), si
    // bien que les annulations entre eux ne coûtent rien.
    constexpr double expm1(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return x + x;
        if (x > 0x1.62e42fefa39efp+9)
            return infinity;
        if (x < -38.0)                                              // e^x < 2^-54 : −1 + e^x s'arrondit à −1
            return -1.0;
        if (x == 0)
            return x;                                               // garde le signe de −0

        auto const [k, r_hi, r_lo] = reduce_exp(x);
        int const j = exp_index(k);
        int const e = exp_exponent(k);                              // dans [-55, 1024]
        double const t_hi = exp_table[j][0];
        double const t_lo = exp_table[j][1];

        // exp(r) − 1 = r + r² [(1/2 + r/6) + r² (1/24 + r/120) + r⁴ (1/720 + r/5040)] = r_hi + p_lo
        double const r = r_hi;
        double const rr = r * r;
        double const a = 0x1p-1 + fmul(r, 0x1.5555555555555p-3);
        double const b = 0x1.5555555555555p-5 + fmul(r, 0x1.1111111111111p-7);
        double const c = 0x1.6c16c16c16c17p-10 + fmul(r, 0x1.a01a01a01a01ap-13);
        double const p_lo = r_lo + fmul(rr, (a + fmul(rr, b)) + fmul(rr * rr, c));

        if (e > 60)                                                 // e^x > 2^60 : le −1 est négligeable
            return scale_exp(t_hi, t_lo + fmul(t_hi, r_hi + p_lo), e);

        double const scale = power_of_two(e);
        auto const [s, s_error] = two_sum(fmul(scale, t_hi), -1.0);
        auto const [m, m_error] = two_product(t_hi, r_hi);
        auto const [s2, s2_error] = two_sum(s, fmul(scale, m));
        double const rest = fmul(scale, (m_error + fmul(t_hi, p_lo)) + fmul(t_lo, 1.0 + (r_hi + p_lo)));
        return s2 + ((s2_error + s_error) + rest);
    }

    // -----------------------------------------------------------------------
    // log, log2, log10, log1p
    // -----------------------------------------------------------------------
    constexpr double log(double x) noexcept
    {
        std::uint64_t ix = 0;
        double special = 0;
        if (detail::log_special_case(x, ix, special))
            return special;
        auto const [hi, lo] = detail::log_kernel(ix);
        return hi + lo;
    }

    // log(x) · (1/ln 2) en double-double : log2(2^n) = n exactement.
    constexpr double log2(double x) noexcept
    {
        return detail::log_times(x, detail::inv_ln2_hi, detail::inv_ln2_lo);
    }

    constexpr double log10(double x) noexcept
    {
        return detail::log_times(x, detail::inv_ln10_hi, detail::inv_ln10_lo);
    }

    // log(1 + x), précis aussi pour x petit.
    // - |x| < 2^-10 : x + (log1p(x) − x) directement, avec le polynôme de log ; x
    //   reste exact, un seul arrondi. (Passer par 1 + x ne garderait que quelques
    //   bits de x, et la correction aurait presque la taille du résultat.)
    // - sinon : 1 + x = u + u_error exactement (two_sum), puis
    //   log(1 + x) = log(u) + log1p(u_error/u) ≈ log(u) + u_error/u, le terme
    //   négligé étant inférieur à 2^-106 relativement.
    constexpr double log1p(double x) noexcept
    {
        using namespace detail;
        if (!(x == x))
            return x + x;
        if (x == 0)
            return x;                                               // garde le signe de −0
        if (x < -1.0)
            return quiet_nan;
        if (x == -1.0)
            return -infinity;
        if (x == infinity)
            return x;
        if (x > -0x1p-10 && x < 0x1p-10)
            return x + log1p_minus_identity(x);
        auto const [u, u_error] = two_sum(1.0, x);                  // u >= 2^-53 : normal
        auto const [hi, lo] = log_kernel(bits(u));
        return hi + (lo + u_error / u);
    }

    // -----------------------------------------------------------------------
    // float : calcul en double, puis un seul arrondi vers float. L'erreur du
    // double (< 0,51 ulp double, soit < 2^-29 ulp float) ne change l'arrondi que
    // si la valeur exacte est à moins de 2^-29 ulp d'un milieu de deux floats :
    // le résultat est presque toujours l'arrondi correct, et toujours à moins de
    // 0,5 + 2^-28 ulp.
    // -----------------------------------------------------------------------
    namespace detail
    {
        constexpr float to_float(double y) noexcept
        {
            if (y >= 0x1.ffffffp+127)                               // ≥ FLT_MAX + ulp/2 : +∞
                return static_cast<float>(infinity);
            return static_cast<float>(y);
        }
    }

    constexpr float exp(float x) noexcept { return detail::to_float(exp(static_cast<double>(x))); }
    constexpr float exp2(float x) noexcept { return detail::to_float(exp2(static_cast<double>(x))); }
    constexpr float expm1(float x) noexcept { return detail::to_float(expm1(static_cast<double>(x))); }
    constexpr float log(float x) noexcept { return static_cast<float>(log(static_cast<double>(x))); }
    constexpr float log2(float x) noexcept { return static_cast<float>(log2(static_cast<double>(x))); }
    constexpr float log10(float x) noexcept { return static_cast<float>(log10(static_cast<double>(x))); }
    constexpr float log1p(float x) noexcept { return static_cast<float>(log1p(static_cast<double>(x))); }

    // Entiers : comme <cmath>, traités en double.
    template <std::integral I> constexpr double exp(I x) noexcept { return exp(static_cast<double>(x)); }
    template <std::integral I> constexpr double exp2(I x) noexcept { return exp2(static_cast<double>(x)); }
    template <std::integral I> constexpr double expm1(I x) noexcept { return expm1(static_cast<double>(x)); }
    template <std::integral I> constexpr double log(I x) noexcept { return log(static_cast<double>(x)); }
    template <std::integral I> constexpr double log2(I x) noexcept { return log2(static_cast<double>(x)); }
    template <std::integral I> constexpr double log10(I x) noexcept { return log10(static_cast<double>(x)); }
    template <std::integral I> constexpr double log1p(I x) noexcept { return log1p(static_cast<double>(x)); }

    // long double n'a pas la même taille partout (64, 80 ou 128 bits) : aucun
    // résultat bit-exact n'est possible.
    long double exp(long double) = delete;
    long double exp2(long double) = delete;
    long double expm1(long double) = delete;
    long double log(long double) = delete;
    long double log2(long double) = delete;
    long double log10(long double) = delete;
    long double log1p(long double) = delete;
}
