#pragma once
// Tirages entiers exactement uniformes, par rejet.
//
// Ces fonctions ne dépendent que des sorties 64 bits du générateur : si elles
// sont uniformes et indépendantes, les tirages sont exactement uniformes, quelle
// que soit l'implémentation de la bibliothèque standard (on n'utilise ni
// std::uniform_int_distribution ni std::shuffle, dont l'algorithme n'est pas
// imposé par la norme).

#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <utility>

namespace wsl {

// Entier signé sur 128 bits (extension GCC/Clang, d'où __extension__ sous -Wpedantic).
__extension__ typedef __int128 int128;

// Entier uniforme dans [0, bound), bound > 0. On rejette les x < 2^64 mod bound :
// il reste 2^64 - (2^64 mod bound) valeurs, un multiple de bound.
inline std::uint64_t uniform_below(std::uint64_t bound, std::mt19937_64& rng)
{
    if (bound == 0)
        throw std::invalid_argument("uniform_below : borne nulle");
    const std::uint64_t threshold = (0 - bound) % bound;   // 2^64 mod bound
    while (true) {
        const std::uint64_t x = rng();
        if (x >= threshold)
            return x % bound;
    }
}

// Même chose sur 128 bits (deux sorties du générateur par essai).
inline int128 uniform_below_128(int128 bound, std::mt19937_64& rng)
{
    __extension__ typedef unsigned __int128 U;
    if (bound <= 0)
        throw std::invalid_argument("uniform_below_128 : borne non positive");
    const U b = static_cast<U>(bound);
    if (b <= U{UINT64_MAX})
        return static_cast<int128>(uniform_below(static_cast<std::uint64_t>(b), rng));
    const U threshold = (U{0} - b) % b;   // 2^128 mod b
    while (true) {
        const U x = (U{rng()} << 64) | U{rng()};
        if (x >= threshold)
            return static_cast<int128>(x % b);
    }
}

// Mélange de Fisher–Yates : chacune des k! permutations avec probabilité 1/k!.
template <class Container>
void exact_shuffle(Container& c, std::mt19937_64& rng)
{
    for (std::size_t i = c.size(); i > 1; --i)
        std::swap(c[i - 1], c[static_cast<std::size_t>(uniform_below(i, rng))]);
}

// Générateur dont tout l'état (19 968 bits) est initialisé par std::random_device.
// Une graine de 64 bits n'atteint que 2^64 états de départ, moins que les
// ≈ 2^72,5 grilles 9×9 : le premier tirage d'une graine fixée ne peut pas couvrir
// toutes les grilles. Ce générateur n'a pas cette limite (mais n'est pas reproductible).
inline std::mt19937_64 seeded_engine()
{
    // 624 mots de 32 bits = la taille de l'état du Mersenne Twister 64 bits.
    std::random_device device;
    std::array<std::uint32_t, 624> words{};
    for (auto& w : words)
        w = device();
    std::seed_seq seq(words.begin(), words.end());
    return std::mt19937_64(seq);
}

}  // namespace wsl
