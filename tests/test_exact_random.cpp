#include <wsl/exact_random.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>

using namespace wsl;

int main()
{
    std::mt19937_64 rng(2028);

    // Bornes et uniformité grossière : 7 classes, 10 000 tirages attendus chacune.
    {
        std::array<std::size_t, 7> counts{};
        for (int k = 0; k < 70000; ++k)
            ++counts[uniform_below(7, rng)];
        for (auto c : counts)
            assert(c > 9500 && c < 10500);
    }

    // Borne 1 : toujours 0. Borne maximale : jamais de rejet infini.
    for (int k = 0; k < 100; ++k)
    {
        assert(uniform_below(1, rng) == 0);
        (void)uniform_below(UINT64_MAX, rng);
    }

    // Borne nulle refusée.
    {
        bool thrown = false;
        try { (void)uniform_below(0, rng); }
        catch (std::invalid_argument const &) { thrown = true; }
        assert(thrown);
    }

    // 128 bits : sous la borne, et au-delà de 2^64 pour une grande borne.
    {
        int128 const bound = int128{1} << 100;
        bool above_64 = false;
        for (int k = 0; k < 1000; ++k)
        {
            int128 const x = uniform_below_128(bound, rng);
            assert(x >= 0 && x < bound);
            above_64 = above_64 || x > int128{UINT64_MAX};
        }
        assert(above_64);
        // Petite borne : délègue à la version 64 bits.
        for (int k = 0; k < 100; ++k)
            assert(uniform_below_128(3, rng) < 3);

        bool thrown = false;
        try { (void)uniform_below_128(-5, rng); }
        catch (std::invalid_argument const &) { thrown = true; }
        assert(thrown);
    }

    // Fisher–Yates : une permutation, et les 3! = 6 permutations également fréquentes.
    {
        std::map<std::array<int, 3>, std::size_t> seen;
        for (int k = 0; k < 60000; ++k)
        {
            std::array<int, 3> values{0, 1, 2};
            exact_shuffle(values, rng);
            ++seen[values];
        }
        assert(seen.size() == 6);
        for (auto const &[permutation, count] : seen)
        {
            assert(std::ranges::is_permutation(permutation, std::array{0, 1, 2}));
            assert(count > 9500 && count < 10500);
        }

        std::array<int, 0> none{};
        exact_shuffle(none, rng);   // rien à faire, sans erreur
    }

    // Générateur entièrement initialisé : deux générateurs diffèrent.
    {
        auto a = seeded_engine();
        auto b = seeded_engine();
        assert(a() != b() || a() != b());
    }

    return 0;
}
