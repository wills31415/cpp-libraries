// Preuve de la garantie bit-exacte de wsl::special.
//
// Pour chaque fonction, empreinte (FNV-1a 64 bits) des résultats sur 2^20 entrées
// fixées, couvrant tout le domaine. Les empreintes attendues ont été relevées une
// fois : toute plateforme, tout compilateur et toute option de compilation doivent
// les reproduire. Ce test tourne sur les environnements de la CI (x86-64 et ARM64,
// GCC et Clang, avec et sans fma) : c'est lui qui prouve la garantie, au lieu de
// la supposer.
//
// Une empreinte qui change sur une seule plateforme : la garantie est cassée, à
// corriger. Une empreinte qui change partout, après une modification voulue d'un
// algorithme : c'est une rupture, à noter au CHANGELOG, puis mettre à jour ici.

#include <wsl/special.hpp>

#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <utility>

namespace sp = wsl::special;

namespace
{
    // Entrées déterministes (SplitMix64), indépendantes de la bibliothèque standard.
    struct splitmix64
    {
        std::uint64_t state;
        std::uint64_t operator()()
        {
            std::uint64_t z = (state += 0x9e3779b97f4a7c15);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
            z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
            return z ^ (z >> 31);
        }
    };

    struct fnv1a
    {
        std::uint64_t value = 0xcbf29ce484222325;
        void add(std::uint64_t bits)
        {
            for (int i = 0; i < 8; ++i)
            {
                value ^= (bits >> (8 * i)) & 0xff;
                value *= 0x100000001b3;
            }
        }
    };

    constexpr int count = 1 << 20;

    // Réel uniforme de [low, high]. Les entrées doivent elles aussi être les mêmes
    // partout : le produit passe par la barrière de wsl::special, sans quoi
    // -ffp-contract=fast en ferait une fma et changerait les entrées, pas les
    // fonctions testées.
    double uniform(splitmix64& g, double low, double high)
    {
        double const u = double(g() >> 11) * 0x1p-53;
        return low + wsl::special::detail::fmul(high - low, u);
    }

    template <class F>
    std::uint64_t fingerprint(F f, int n = count)
    {
        fnv1a h;
        splitmix64 g{2026};
        for (int i = 0; i < n; ++i)
            h.add(f(g));
        return h.value;
    }

    bool check(char const* name, std::uint64_t actual, std::uint64_t expected)
    {
        bool const ok = actual == expected;
        std::printf("%-16s %016llx %s\n", name, static_cast<unsigned long long>(actual),
                    ok ? "ok" : "DIFFÉRENTE de l'empreinte attendue");
        std::fflush(stdout);
        return ok;
    }
}

// Chaque tirage est une instruction à part : deux appels à g() dans une même
// expression seraient évalués dans un ordre que la norme ne fixe pas.
int main()
{
    // exp : moitié uniforme sur tout le domaine, moitié petites valeurs de tous
    // exposants (bits aléatoires sous 2), de signe aléatoire.
    auto const exp_double = fingerprint([](splitmix64& g) {
        std::uint64_t const choice = g();
        double x = 0;
        if (choice & 1)
            x = uniform(g, -746, 710);
        else
        {
            x = std::bit_cast<double>(g() & 0x3fefffffffffffff);
            if (choice & 2)
                x = -x;
        }
        return std::bit_cast<std::uint64_t>(sp::exp(x));
    });
    // log : moitié bits aléatoires (tous les doubles finis positifs, sous-normaux
    // compris), moitié près de 1.
    auto const log_double = fingerprint([](splitmix64& g) {
        std::uint64_t const choice = g();
        double const x = (choice & 1) ? std::bit_cast<double>(g() & 0x7fefffffffffffff) : uniform(g, 0.5, 2.0);
        return std::bit_cast<std::uint64_t>(sp::log(x));
    });
    auto const exp_float = fingerprint([](splitmix64& g) {
        float const x = static_cast<float>(uniform(g, -104, 89));
        return std::uint64_t{std::bit_cast<std::uint32_t>(sp::exp(x))};
    });
    auto const log_float = fingerprint([](splitmix64& g) {
        float const x = std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0x7f7fffff));
        return std::uint64_t{std::bit_cast<std::uint32_t>(sp::log(x))};
    });

    // Nouvelles fonctions : mêmes principes. Petites valeurs signées de tous
    // exposants : bits aléatoires sous 1 (sous 2 pour exp2).
    auto const small = [](splitmix64& g, std::uint64_t mask) {
        std::uint64_t const sign = g();
        double const x = std::bit_cast<double>(g() & mask);
        return (sign & 1) ? -x : x;
    };
    auto const exp2_input = [&](splitmix64& g) {
        std::uint64_t const choice = g();
        return (choice & 1) ? uniform(g, -1080, 1024) : small(g, 0x3fefffffffffffff);
    };
    auto const expm1_input = [&](splitmix64& g) {
        std::uint64_t const choice = g();
        return (choice & 1) ? uniform(g, -40, 710) : small(g, 0x3fdfffffffffffff);
    };
    auto const positive_input = [&](splitmix64& g) {
        std::uint64_t const choice = g();
        return (choice & 1) ? std::bit_cast<double>(g() & 0x7fefffffffffffff) : uniform(g, 0.5, 2.0);
    };
    auto const log1p_input = [&](splitmix64& g) {
        std::uint64_t const choice = g();
        return (choice & 1) ? uniform(g, -1, 1e6) : small(g, 0x3fdfffffffffffff);
    };
    auto const of_double = [](double y) { return std::bit_cast<std::uint64_t>(y); };
    auto const of_float = [](float y) { return std::uint64_t{std::bit_cast<std::uint32_t>(y)}; };

    auto const exp2_double = fingerprint([&](splitmix64& g) { return of_double(sp::exp2(exp2_input(g))); });
    auto const expm1_double = fingerprint([&](splitmix64& g) { return of_double(sp::expm1(expm1_input(g))); });
    auto const log2_double = fingerprint([&](splitmix64& g) { return of_double(sp::log2(positive_input(g))); });
    auto const log10_double = fingerprint([&](splitmix64& g) { return of_double(sp::log10(positive_input(g))); });
    auto const log1p_double = fingerprint([&](splitmix64& g) { return of_double(sp::log1p(log1p_input(g))); });
    auto const exp2_float = fingerprint([&](splitmix64& g) { return of_float(sp::exp2(static_cast<float>(uniform(g, -155, 130)))); });
    auto const expm1_float = fingerprint([&](splitmix64& g) { return of_float(sp::expm1(static_cast<float>(expm1_input(g) / 8))); });
    auto const log2_float = fingerprint([&](splitmix64& g) { return of_float(sp::log2(std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0x7f7fffff)))); });
    auto const log10_float = fingerprint([&](splitmix64& g) { return of_float(sp::log10(std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0x7f7fffff)))); });
    auto const log1p_float = fingerprint([&](splitmix64& g) { return of_float(sp::log1p(static_cast<float>(log1p_input(g)))); });

    // pow : x quelconque (bits aléatoires) ou proche de 1, y uniforme ; et bases
    // négatives à exposant entier.
    auto const pow_double = fingerprint([&](splitmix64& g) {
        std::uint64_t const choice = g();
        double x = 0, y = 0;
        if (choice & 1)
        {
            x = std::bit_cast<double>(g() & 0x7fefffffffffffff);
            y = uniform(g, -3, 3);
        }
        else if (choice & 2)
        {
            x = uniform(g, 0.999, 1.001);
            y = uniform(g, -5e5, 5e5);
        }
        else
        {
            x = uniform(g, -10, 10);
            y = static_cast<double>(static_cast<std::int64_t>(g() % 601) - 300);
        }
        return of_double(sp::pow(x, y));
    });
    auto const pow_float = fingerprint([&](splitmix64& g) {
        float const x = std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0x7f7fffff));
        float const y = static_cast<float>(uniform(g, -3, 3));
        return of_float(sp::pow(x, y));
    });

    // Trigonométrie : petits, moyens et immenses arguments (bits aléatoires).
    auto const trig_input = [&](splitmix64& g) {
        std::uint64_t const choice = g() % 3;
        if (choice == 0)
            return uniform(g, -10, 10);
        if (choice == 1)
            return uniform(g, -0x1p20, 0x1p20);
        return std::bit_cast<double>(g() & 0xffefffffffffffff);           // tous les doubles finis
    };
    auto const sin_double = fingerprint([&](splitmix64& g) { return of_double(sp::sin(trig_input(g))); });
    auto const cos_double = fingerprint([&](splitmix64& g) { return of_double(sp::cos(trig_input(g))); });
    auto const tan_double = fingerprint([&](splitmix64& g) { return of_double(sp::tan(trig_input(g))); });
    auto const atan_double = fingerprint([&](splitmix64& g) { return of_double(sp::atan(std::bit_cast<double>(g() & 0xffefffffffffffff))); });
    auto const atan2_double = fingerprint([&](splitmix64& g) {
        double const y = std::bit_cast<double>(g() & 0xffefffffffffffff);
        double const x = std::bit_cast<double>(g() & 0xffefffffffffffff);
        return of_double(sp::atan2(y, x));
    });
    auto const sin_float = fingerprint([&](splitmix64& g) { return of_float(sp::sin(std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0xff7fffff)))); });
    auto const cos_float = fingerprint([&](splitmix64& g) { return of_float(sp::cos(std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0xff7fffff)))); });
    auto const tan_float = fingerprint([&](splitmix64& g) { return of_float(sp::tan(std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0xff7fffff)))); });
    auto const atan_float = fingerprint([&](splitmix64& g) { return of_float(sp::atan(std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0xff7fffff)))); });
    auto const atan2_float = fingerprint([&](splitmix64& g) {
        float const y = std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0xff7fffff));
        float const x = std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0xff7fffff));
        return of_float(sp::atan2(y, x));
    });

    // Fonctions algébriques : tous les doubles finis (bits aléatoires).
    auto const any_double = [](splitmix64& g) { return std::bit_cast<double>(g() & 0xffefffffffffffff); };
    auto const cbrt_double = fingerprint([&](splitmix64& g) { return of_double(sp::cbrt(any_double(g))); });
    auto const hypot_double = fingerprint([&](splitmix64& g) {
        double const x = any_double(g);
        double const y = any_double(g);
        return of_double(sp::hypot(x, y));
    });
    auto const cbrt_float = fingerprint([&](splitmix64& g) { return of_float(sp::cbrt(std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0xff7fffff)))); });
    auto const hypot_float = fingerprint([&](splitmix64& g) {
        float const x = std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0xff7fffff));
        float const y = std::bit_cast<float>(static_cast<std::uint32_t>(g() & 0xff7fffff));
        return of_float(sp::hypot(x, y));
    });

    // Fonction d'erreur : toute la zone utile (|x| < 30), et de petites valeurs.
    auto const erf_input = [&](splitmix64& g) {
        std::uint64_t const choice = g();
        return (choice & 1) ? uniform(g, -8, 30) : small(g, 0x3fdfffffffffffff);
    };
    auto const erf_double = fingerprint([&](splitmix64& g) { return of_double(sp::erf(erf_input(g))); });
    auto const erfc_double = fingerprint([&](splitmix64& g) { return of_double(sp::erfc(erf_input(g))); });
    auto const erf_float = fingerprint([&](splitmix64& g) { return of_float(sp::erf(static_cast<float>(erf_input(g)))); });
    auto const erfc_float = fingerprint([&](splitmix64& g) { return of_float(sp::erfc(static_cast<float>(erf_input(g)))); });

    // Gamma : positifs et négatifs (non entiers en pratique), tous les ordres de grandeur.
    auto const gamma_input = [&](splitmix64& g) {
        std::uint64_t const choice = g() % 3;
        if (choice == 0)
            return uniform(g, -40, 40);
        if (choice == 1)
            return uniform(g, -190, 180);
        return std::bit_cast<double>(g() & 0xffefffffffffffff);
    };
    auto const lgamma_double = fingerprint([&](splitmix64& g) { return of_double(sp::lgamma(gamma_input(g))); });
    auto const tgamma_double = fingerprint([&](splitmix64& g) { return of_double(sp::tgamma(gamma_input(g))); });
    auto const lgamma_float = fingerprint([&](splitmix64& g) { return of_float(sp::lgamma(static_cast<float>(uniform(g, -40, 1e6)))); });
    auto const tgamma_float = fingerprint([&](splitmix64& g) { return of_float(sp::tgamma(static_cast<float>(uniform(g, -40, 36)))); });
    auto const digamma_double = fingerprint([&](splitmix64& g) { return of_double(sp::digamma(gamma_input(g))); });
    auto const digamma_float = fingerprint([&](splitmix64& g) { return of_float(sp::digamma(static_cast<float>(uniform(g, -40, 1e6)))); });
    auto const igamma_input = [&](splitmix64& g) {
        double const a = std::bit_cast<double>((std::uint64_t{1023 - 30 + g() % 64} << 52) | (g() >> 12));   // 2^-30 … 2^34
        double const x = g() % 2 == 0 ? a + wsl::special::detail::fmul(uniform(g, -8, 8), sp::sqrt(a))
                                      : a * std::bit_cast<double>((std::uint64_t{1023 - 6 + g() % 12} << 52) | (g() >> 12));
        return std::pair{a, x < 0 ? -x : x};
    };
    auto const gamma_p_double = fingerprint([&](splitmix64& g) { auto const [a, x] = igamma_input(g); return of_double(sp::gamma_p(a, x)); }, 1 << 15);
    auto const gamma_q_double = fingerprint([&](splitmix64& g) { auto const [a, x] = igamma_input(g); return of_double(sp::gamma_q(a, x)); }, 1 << 15);
    auto const gamma_q_float = fingerprint([&](splitmix64& g) { auto const [a, x] = igamma_input(g); return of_float(sp::gamma_q(static_cast<float>(a), static_cast<float>(x))); }, 1 << 15);
    auto const gamma_inverse_input = [&](splitmix64& g) {
        double const a = std::bit_cast<double>((std::uint64_t{1023 - 20 + g() % 44} << 52) | (g() >> 12));   // 2^-20 … 2^24
        double const t = g() % 2 == 0 ? uniform(g, 0, 1) : std::bit_cast<double>((std::uint64_t{1023 - 1 - g() % 900} << 52) | (g() >> 12));
        return std::pair{a, t};
    };
    auto const gamma_p_inv_double = fingerprint([&](splitmix64& g) { auto const [a, t] = gamma_inverse_input(g); return of_double(sp::gamma_p_inv(a, t)); }, 1 << 12);
    auto const gamma_q_inv_double = fingerprint([&](splitmix64& g) { auto const [a, t] = gamma_inverse_input(g); return of_double(sp::gamma_q_inv(a, t)); }, 1 << 12);

    // Inverses de la fonction d'erreur : tout (−1, 1) et tout (0, 2), extrémités comprises.
    auto const erf_inv_double = fingerprint([&](splitmix64& g) {
        std::uint64_t const choice = g();
        double const y = (choice & 1) ? uniform(g, -1, 1) : 1.0 - std::bit_cast<double>(g() & 0x3fefffffffffffff);
        return of_double(sp::erf_inv((choice & 2) ? -y : y));
    });
    auto const erfc_inv_double = fingerprint([&](splitmix64& g) {
        std::uint64_t const choice = g();
        double const z = (choice & 1) ? uniform(g, 0, 2) : std::bit_cast<double>(g() & 0x3fefffffffffffff);
        return of_double(sp::erfc_inv(z));
    });
    auto const erf_inv_float = fingerprint([&](splitmix64& g) { return of_float(sp::erf_inv(static_cast<float>(uniform(g, -1, 1)))); });
    auto const erfc_inv_float = fingerprint([&](splitmix64& g) { return of_float(sp::erfc_inv(static_cast<float>(uniform(g, 0, 2)))); });

    bool ok = true;
    ok &= check("exp double", exp_double, 0xdc7c71a0d207c234);
    ok &= check("log double", log_double, 0x75529ab2341c6559);
    ok &= check("exp float", exp_float, 0xe70dae8bac4c54dc);
    ok &= check("log float", log_float, 0xfd2df93b55a28f5e);
    ok &= check("exp2 double", exp2_double, 0x57b305f423cd4aa9);
    ok &= check("expm1 double", expm1_double, 0x33d154266c84cd52);
    ok &= check("log2 double", log2_double, 0x7bf1c6c8661a8ae2);
    ok &= check("log10 double", log10_double, 0xc5088a1926d99501);
    ok &= check("log1p double", log1p_double, 0x6275560ffbe3d9a9);
    ok &= check("exp2 float", exp2_float, 0xd61d09d825ea0e17);
    ok &= check("expm1 float", expm1_float, 0x4454c409717fbdf3);
    ok &= check("log2 float", log2_float, 0x386f09ee6688ef03);
    ok &= check("log10 float", log10_float, 0x2feaf388eb377657);
    ok &= check("log1p float", log1p_float, 0x79935e48bb26c27c);
    ok &= check("pow double", pow_double, 0x979387e8dfa097bb);
    ok &= check("pow float", pow_float, 0x554a23153033bb57);
    ok &= check("sin double", sin_double, 0x698975b623670a55);
    ok &= check("cos double", cos_double, 0x73c98d0ba6c856d1);
    ok &= check("tan double", tan_double, 0xa6ff7c3298679e4c);
    ok &= check("atan double", atan_double, 0xc0eaff1c7939caa3);
    ok &= check("atan2 double", atan2_double, 0xc11e0b5fb09d6dbd);
    ok &= check("sin float", sin_float, 0xf3dd2fa428f21680);
    ok &= check("cos float", cos_float, 0x216da6f1a9bb6224);
    ok &= check("tan float", tan_float, 0xf847c6dbd98a670d);
    ok &= check("atan float", atan_float, 0x46044c8b18db7619);
    ok &= check("atan2 float", atan2_float, 0x8138cf20810c3d95);
    ok &= check("cbrt double", cbrt_double, 0x130f4c4a6e6bc014);
    ok &= check("hypot double", hypot_double, 0x6571c9f1c46dc8d9);
    ok &= check("cbrt float", cbrt_float, 0x2f8352bd44cdbd49);
    ok &= check("hypot float", hypot_float, 0x72174db6500aa335);
    ok &= check("erf double", erf_double, 0xac731c85c5f4811d);
    ok &= check("erfc double", erfc_double, 0xcca7cee2fa853ae9);
    ok &= check("erf float", erf_float, 0x8d664252d9a17955);
    ok &= check("erfc float", erfc_float, 0x8cef7e01c22f4892);
    ok &= check("lgamma double", lgamma_double, 0x4397c43ee6b5bcb9);
    ok &= check("tgamma double", tgamma_double, 0x8bb4f651e6380bcb);
    ok &= check("lgamma float", lgamma_float, 0xea350830d6044a96);
    ok &= check("tgamma float", tgamma_float, 0x99520856f145ad49);
    ok &= check("digamma double", digamma_double, 0x3a9ec0df558f01bc);
    ok &= check("digamma float", digamma_float, 0x260ca54d39c41921);
    ok &= check("gamma_p double", gamma_p_double, 0xe7b7bbc1b3a35053);
    ok &= check("gamma_q double", gamma_q_double, 0x9782b6a6da56fe8d);
    ok &= check("gamma_q float", gamma_q_float, 0x84058102c180cb27);
    ok &= check("gamma_p_inv dbl", gamma_p_inv_double, 0x5cf887df114e68d5);
    ok &= check("gamma_q_inv dbl", gamma_q_inv_double, 0xd762c5cf914f3b15);
    ok &= check("erf_inv double", erf_inv_double, 0xf3178cd1f5ba8e43);
    ok &= check("erfc_inv double", erfc_inv_double, 0x4f276a282a6d14f2);
    ok &= check("erf_inv float", erf_inv_float, 0x11aceac4e44c6603);
    ok &= check("erfc_inv float", erfc_inv_float, 0xc04b1a3cf1bb8efb);
    assert(ok);
    return 0;
}
