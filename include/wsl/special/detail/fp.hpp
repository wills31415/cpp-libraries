#pragma once

// Socle de wsl::special : conditions de la garantie bit-exacte, et opérations
// flottantes élémentaires dont le résultat est imposé par IEEE 754.
//
// Toute fonction de wsl::special est écrite avec les seules opérations + − × ÷,
// dont la norme IEEE 754 impose l'arrondi correct : leur résultat est donc le même
// sur toute plateforme. Trois choses pourraient le changer, et sont exclues ici :
//
//   - les optimisations « rapides » (-ffast-math) : refusées à la compilation ;
//   - le calcul en précision étendue (x87) : refusé à la compilation ;
//   - la contraction de a*b + c en une fma, qu'un compilateur a le droit de faire
//     (Clang par défaut sur ARM64, GCC dès que la cible a une fma, et les deux
//     sous -ffp-contract=fast, même avec #pragma) : tout produit dont le résultat
//     entre dans une addition passe par fmul, qui la rend impossible.
//
// Le mode d'arrondi est supposé être celui par défaut (au plus proche, pair).

#include <bit>
#include <cfloat>
#include <cstdint>
#include <limits>
#include <type_traits>

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
#error "wsl::special : -ffast-math (ou -ffinite-math-only) rend la garantie bit-exacte impossible"
#endif

#if !defined(FLT_EVAL_METHOD) || FLT_EVAL_METHOD != 0
#error "wsl::special : le calcul flottant doit se faire dans le type des opérandes (FLT_EVAL_METHOD == 0)"
#endif

#if !defined(__GNUC__)
#error "wsl::special : seuls GCC et Clang sont pris en charge (barrière de contraction en assembleur)"
#endif

static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<float>::is_iec559,
              "wsl::special : float et double doivent être les formats binary32 et binary64 d'IEEE 754");

namespace wsl::special::detail
{
    // Rend la valeur opaque pour l'optimiseur : rien ne peut être fusionné à
    // travers. Sans effet en évaluation constante, qui ne fusionne jamais.
    template <class T>
    constexpr T opaque(T value) noexcept
    {
        if (!std::is_constant_evaluated())
        {
#if defined(__x86_64__)
            asm("" : "+x"(value));
#elif defined(__aarch64__)
            asm("" : "+w"(value));
#else
            asm("" : "+m"(value));
#endif
        }
        return value;
    }

    // Produit destiné à une addition : arrondi une fois, jamais fusionné.
    template <class T>
    constexpr T fmul(T a, T b) noexcept
    {
        return opaque(a * b);
    }

    // Somme et erreur exactes : a + b == sum + error exactement (Knuth).
    struct sum_error
    {
        double sum, error;
    };

    constexpr sum_error two_sum(double a, double b) noexcept
    {
        double const sum = a + b;
        double const b_virtual = sum - a;
        double const a_virtual = sum - b_virtual;
        return {sum, (a - a_virtual) + (b - b_virtual)};
    }

    // Découpage de Veltkamp : x == high + low exactement, high sur 26 bits
    // significatifs et low sur 27, si bien que les produits de deux moitiés sont
    // exacts.
    struct halves
    {
        double high, low;
    };

    constexpr halves split(double x) noexcept
    {
        double const t = fmul(x, 0x1.0000002p27);   // 2^27 + 1
        double const high = t - (t - x);
        return {high, x - high};
    }

    // Produit et erreur exacts : a · b == product + error exactement.
    //
    // Avec une fma matérielle (ARM64 ; x86-64 compilé avec -mfma, -march=haswell
    // ou x86-64-v3…), l'erreur est fma(a, b, −product) : 2 opérations. Sinon, et en
    // évaluation constante, le produit de Dekker : 17 opérations. L'erreur d'un
    // produit étant unique, les deux rendent les mêmes bits — tant que Dekker est
    // exact : facteurs sous 2^995 (le découpage multiplie par 2^27 + 1) et
    // exposants de a et b de somme >= −970 (pas de produit partiel sous-normal).
    // Chaque appel doit le garantir ; les empreintes de test_special_bitexact,
    // relevées sur des cibles avec et sans fma, le vérifient.
    //
    // WSL_SPECIAL_NO_FMA force Dekker partout.
#if (defined(__FMA__) || defined(__FP_FAST_FMA) || defined(__aarch64__)) && !defined(WSL_SPECIAL_NO_FMA)
    inline constexpr bool hardware_fma = true;
#else
    inline constexpr bool hardware_fma = false;
#endif

    constexpr sum_error two_product(double a, double b) noexcept
    {
        double const product = fmul(a, b);
        if constexpr (hardware_fma)
        {
            if (!std::is_constant_evaluated())
                return {product, __builtin_fma(a, b, -product)};
        }
        auto const [a_high, a_low] = split(a);
        auto const [b_high, b_low] = split(b);
        double const error = (((fmul(a_high, b_high) - product) + fmul(a_high, b_low)) + fmul(a_low, b_high))
                           + fmul(a_low, b_low);
        return {product, error};
    }

    // Entier non signé de 128 bits (extension de GCC et Clang, d'où __extension__
    // sous -Wpedantic) : réduction de Payne–Hanek, en arithmétique entière exacte.
    __extension__ typedef unsigned __int128 uint128;

    // Somme exacte quand |a| >= |b| (ou a = 0) : trois opérations (Dekker).
    constexpr sum_error fast_two_sum(double a, double b) noexcept
    {
        double const sum = a + b;
        return {sum, (a - sum) + b};
    }

    // Valeur en double-double : hi + lo, |lo| bien plus petit que |hi|.
    struct double_double
    {
        double hi, lo;
    };

    constexpr std::uint64_t bits(double x) noexcept { return std::bit_cast<std::uint64_t>(x); }
    constexpr double from_bits(std::uint64_t b) noexcept { return std::bit_cast<double>(b); }

    // 2^e, pour e dans [-1022, 1023] (exposant d'un double normal).
    constexpr double power_of_two(int e) noexcept
    {
        return from_bits(static_cast<std::uint64_t>(e + 1023) << 52);
    }

    // Exposant e de x > 0 fini : 2^e <= x < 2^(e+1), sous-normaux compris.
    constexpr int exponent_of(double x) noexcept
    {
        std::uint64_t const b = bits(x);
        if ((b >> 52) == 0)
            return (63 - std::countl_zero(b)) - 1074;
        return static_cast<int>(b >> 52) - 1023;
    }

    // x · 2^n, n quelconque, en plusieurs étapes si n sort de [−1022, 1023]. Les
    // étapes de montée sont exactes jusqu'au débordement éventuel (alors ±∞, voulu) ;
    // en descente, chaque étape intermédiaire reste normale tant que le résultat
    // n'est pas minuscule, si bien qu'il n'y a qu'un arrondi, celui de la dernière,
    // si le résultat est sous-normal.
    constexpr double scale_by_power_of_two(double x, int n) noexcept
    {
        while (n > 1023)
        {
            x = fmul(x, power_of_two(1023));
            n -= 1023;
        }
        while (n < -1022)
        {
            x = fmul(x, power_of_two(-1022));
            n += 1022;
        }
        return fmul(x, power_of_two(n));
    }

    inline constexpr double infinity = std::numeric_limits<double>::infinity();
    inline constexpr double quiet_nan = std::numeric_limits<double>::quiet_NaN();
}
