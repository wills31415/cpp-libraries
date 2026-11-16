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

    bool ok = true;
    ok &= check("exp double", exp_double, 0xdc7c71a0d207c234);
    ok &= check("log double", log_double, 0x75529ab2341c6559);
    ok &= check("exp float", exp_float, 0xe70dae8bac4c54dc);
    ok &= check("log float", log_float, 0xfd2df93b55a28f5e);
    assert(ok);
    return 0;
}
