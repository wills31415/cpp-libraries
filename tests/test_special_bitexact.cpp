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
#include <cstdint>
#include <cstdio>

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
    std::uint64_t fingerprint(F f)
    {
        fnv1a h;
        splitmix64 g{2026};
        for (int i = 0; i < count; ++i)
            h.add(f(g));
        return h.value;
    }

    bool check(char const* name, std::uint64_t actual, std::uint64_t expected)
    {
        bool const ok = actual == expected;
        std::printf("%-12s %016llx %s\n", name, static_cast<unsigned long long>(actual),
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
    assert(ok);
    return 0;
}
