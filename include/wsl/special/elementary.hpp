#pragma once

// Fonctions élémentaires bit-exactes : même résultat au bit près sur toute
// plateforme (voir detail/fp.hpp pour les conditions), constexpr.
//
//   exp(x)   erreur < 0,51 ulp en double ; float calculé en double puis arrondi
//   log(x)   erreur < 0,51 ulp en double ; float calculé en double puis arrondi
//
// Les erreurs maximales mesurées contre mpmath sont vérifiées par
// tests/test_special_elementary.cpp. Cas limites ceux de <cmath> (Annexe F de C) :
// exp(±0) = 1, exp(+∞) = +∞, exp(−∞) = +0 ; log(1) = +0, log(±0) = −∞,
// log(+∞) = +∞, log(x < 0) = NaN ; NaN donne NaN.
//
// Règle d'écriture : tout produit dont le résultat entre dans une addition ou une
// soustraction passe par fmul (voir detail/fp.hpp).

#include <wsl/special/detail/fp.hpp>
#include <wsl/special/detail/tables.hpp>

#include <concepts>
#include <cstdint>

namespace wsl::special
{
    // -----------------------------------------------------------------------
    // exp
    //
    // x = k ln2/N + r, |r| <= ln2/(2N), N = 128, k = N e + j :
    //     exp(x) = 2^e · 2^(j/N) · exp(r).
    // 2^(j/N) vient d'une table en double-double (T_hi + T_lo), exp(r) − 1 d'un
    // polynôme de Taylor de degré 5 (reste < 2^-60 sur l'intervalle). Le résultat
    // T_hi + (T_lo + T_hi p) n'est arrondi qu'une fois à la fin, d'où une erreur
    // à peine supérieure à 1/2 ulp. Réduction de Cody–Waite : k ln2/N est exact
    // (ln2/N sur 35 bits, |k| < 2^18), et x − k ln2_hi aussi (lemme de Sterbenz).
    // -----------------------------------------------------------------------
    constexpr double exp(double x) noexcept
    {
        using namespace detail;

        if (!(x == x))
            return x + x;                                   // NaN
        if (x > 0x1.62e42fefa39efp+9)                       // > ln(DBL_MAX) : +∞
            return infinity;
        if (x < -746.0)                                     // exp(x) < 2^-1076 : +0
            return 0.0;

        double const z = fmul(x, exp_inv_ln2_n);
        auto const k = static_cast<std::int64_t>(z < 0 ? z - 0.5 : z + 0.5);
        double const kd = static_cast<double>(k);
        double const r = (x - fmul(kd, exp_ln2_n_hi)) - fmul(kd, exp_ln2_n_lo);

        constexpr std::int64_t n = std::int64_t{1} << exp_table_bits;
        auto const j = static_cast<int>(k & (n - 1));
        auto const e = static_cast<int>(k >> exp_table_bits);      // division entière par défaut

        // exp(r) − 1 = r + r²/2 + r³/6 + r⁴/24 + r⁵/120
        double const r2 = r * r;
        double const p = r + fmul(r2, 0x1p-1 + fmul(r, 0x1.5555555555555p-3
                                   + fmul(r, 0x1.5555555555555p-5 + fmul(r, 0x1.1111111111111p-7))));

        double const t_hi = exp_table[j][0];
        double const tail = exp_table[j][1] + fmul(t_hi, p);

        if (e >= -1022)
        {
            double const y = t_hi + tail;                   // 2^(j/N) exp(r), dans [0,99 ; 2)
            if (e > 1023)                                   // résultat juste sous DBL_MAX
                return y * power_of_two(1023) * 2.0;
            return y * power_of_two(e);
        }

        // Résultat sous-normal. Arrondir y à 53 bits puis le ramener dans la zone
        // sous-normale, de précision moindre, arrondirait deux fois (jusqu'à 0,75
        // ulp). On travaille donc sur y = résultat · 2^1022, dans (0, 1) : 1 + y
        // arrondit au bit voulu (2^-52, soit 2^-1074 une fois remis à l'échelle),
        // et l'erreur de chaque étape est reportée, si bien que le résultat n'est
        // arrondi qu'une fois. La remise à l'échelle finale est exacte.
        double const scale = power_of_two(e + 1022);         // e + 1022 dans [-55, -1]
        double const a = fmul(scale, t_hi);                 // exacts : puissance de 2
        double const b = fmul(scale, tail);
        double const y = a + b;
        double const y_error = (a - y) + b;                 // exacte, |a| >= |b|
        double const high = 1.0 + y;
        double const low = ((1.0 - high) + y) + y_error;
        return ((high + low) - 1.0) * 0x1p-1022;
    }

    // -----------------------------------------------------------------------
    // log
    //
    // x = 2^e z avec z dans [√2/2, √2), puis z = c (1 + r) avec c = j/512 :
    //     log(x) = e ln2 + log(c) + log1p(r),  |r| < 2^-9,5.
    //
    // - Réduction sans branche, en arithmétique entière sur les bits de x : z garde
    //   la mantisse de x, e est l'exposant décalé de sorte que z tombe dans
    //   [√2/2, √2). Près de 1, e = 0 : pas d'annulation entre e ln2 et log(z).
    // - j = round(512 z) par l'astuce du décalage (2^52 · 1,5) : sans conversion.
    // - z − c est exact (Sterbenz). r est obtenu en double-double sans division :
    //   r_hi = (z − c)·(1/c), 1/c tabulé, puis le résidu z − c − r_hi·c, exact par
    //   un produit de Dekker (c n'a que 10 bits), donne r_lo.
    // - log(c) est tabulé en deux parties, hi multiple de 2^-42 comme e ln2_hi :
    //   w = e ln2_hi + log(c)_hi est exact. Puis |w| >= |r| (ou w = 0 quand c = 1),
    //   d'où une somme exacte w + r_hi en trois opérations (Fast2Sum).
    // - log1p(r) − r : Taylor de degré 7 (reste < 2^-69 relativement à r), évalué
    //   par le schéma d'Estrin, dont les termes se calculent en parallèle.
    // Le résultat n'est arrondi qu'une fois, à la dernière addition.
    // -----------------------------------------------------------------------
    constexpr double log(double x) noexcept
    {
        using namespace detail;

        std::uint64_t ix = bits(x);
        // Un seul test pour tous les cas à part : x <= 0 (le bit de signe rend ix
        // très grand), sous-normal, +∞ ou NaN.
        if (ix - 0x0010000000000000 >= 0x7fe0000000000000) [[unlikely]]
        {
            if (!(x == x))
                return x + x;                               // NaN
            if (x == 0)
                return -infinity;
            if (x < 0)
                return quiet_nan;
            if (x == infinity)
                return x;
            ix = bits(x * 0x1p52) - (std::uint64_t{52} << 52);   // sous-normal : x · 2^52, exposant corrigé
        }

        constexpr std::uint64_t half_sqrt2 = 0x3fe6a09e667f3bcd;       // bits de √2/2
        std::uint64_t const offset = ix - half_sqrt2;
        auto const e = static_cast<int>(static_cast<std::int64_t>(offset) >> 52);
        double const z = from_bits(ix - (offset & (std::uint64_t{0xfff} << 52)));   // dans [√2/2, √2)

        constexpr double shift = 0x1.8p52;                  // ulp 1 : t − shift = round(512 z)
        constexpr double scale = static_cast<double>(1 << log_table_bits);
        double const t = fmul(z, scale) + shift;            // produit exact
        auto const j = static_cast<int>(bits(t) - bits(shift));   // dans [362, 724]
        double const c = (t - shift) * (1.0 / scale);
        double const d = z - c;                             // exact

        auto const& entry = log_table[j - log_table_first];
        double const inverse_c = entry[2];
        double const r_hi = fmul(d, inverse_c);
        auto const [r1, r2] = split(r_hi);
        double const r_lo = fmul((d - fmul(r1, c)) - fmul(r2, c), inverse_c);

        // log1p(r) − r = r² [(−1/2 + r/3) + r² (−1/4 + r/5) + r⁴ (−1/6 + r/7)]
        double const r = r_hi;
        double const rr = r * r;
        double const p01 = -0x1p-1 + fmul(r, 0x1.5555555555555p-2);
        double const p23 = -0x1p-2 + fmul(r, 0x1.999999999999ap-3);
        double const p45 = -0x1.5555555555555p-3 + fmul(r, 0x1.2492492492492p-3);
        double const q = fmul(rr, (p01 + fmul(rr, p23)) + fmul(rr * rr, p45));

        double const ed = static_cast<double>(e);
        double const w = fmul(ed, log_ln2_hi) + entry[0];  // exact
        double const s = w + r_hi;
        double const s_error = (w - s) + r_hi;              // exact : |w| >= |r_hi| ou w = 0
        double const tail = (fmul(ed, log_ln2_lo) + entry[1]) + r_lo;
        return s + ((s_error + tail) + q);
    }

    // -----------------------------------------------------------------------
    // float : calcul en double, puis un seul arrondi vers float. L'erreur du
    // double (< 0,51 ulp double, soit < 2^-29 ulp float) ne change l'arrondi que
    // si la valeur exacte est à moins de 2^-29 ulp d'un milieu de deux floats :
    // le résultat est presque toujours l'arrondi correct, et toujours à moins de
    // 0,5 + 2^-29 ulp.
    // -----------------------------------------------------------------------
    constexpr float exp(float x) noexcept
    {
        double const y = exp(static_cast<double>(x));
        if (y >= 0x1.ffffffp+127)                           // ≥ FLT_MAX + ulp/2 : +∞
            return static_cast<float>(detail::infinity);
        return static_cast<float>(y);
    }

    constexpr float log(float x) noexcept
    {
        return static_cast<float>(log(static_cast<double>(x)));
    }

    // Entiers : comme <cmath>, traités en double.
    template <std::integral I>
    constexpr double exp(I x) noexcept { return exp(static_cast<double>(x)); }

    template <std::integral I>
    constexpr double log(I x) noexcept { return log(static_cast<double>(x)); }

    // long double n'a pas la même taille partout (64, 80 ou 128 bits) : aucun
    // résultat bit-exact n'est possible.
    long double exp(long double) = delete;
    long double log(long double) = delete;
}
