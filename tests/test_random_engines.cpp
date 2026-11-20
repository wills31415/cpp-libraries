// Moteurs de wsl::random : sorties comparées à des sources indépendantes de ce
// code (tests/random_engine_vectors.hpp : NumPy, matrices sur GF(2)), propriétés
// d'interface, sauts et avances.

#include <wsl/random.hpp>

#include "random_engine_vectors.hpp"

#include <cassert>
#include <cstdint>
#include <random>
#include <span>

namespace v = engine_vectors;

namespace
{
    template <class G>
    bool same_outputs(G g, std::span<std::uint64_t const> expected)
    {
        for (std::uint64_t e : expected)
            if (g() != e)
                return false;
        return true;
    }

    constexpr wsl::random::detail::uint128 u128(std::uint64_t const (&w)[2])
    {
        return (wsl::random::detail::uint128{w[0]} << 64) | w[1];
    }
}

// Interface standard, et utilisables en évaluation constante.
static_assert(std::uniform_random_bit_generator<wsl::random::splitmix64>);
static_assert(std::uniform_random_bit_generator<wsl::random::xoshiro256pp>);
static_assert(std::uniform_random_bit_generator<wsl::random::pcg64>);
static_assert(std::uniform_random_bit_generator<wsl::random::philox4x64>);

constexpr std::uint64_t first_output(auto g) { return g(); }
static_assert(first_output(wsl::random::splitmix64(0)) == 0xe220a8397b1dcdaf);
static_assert(first_output(wsl::random::xoshiro256pp({1, 2, 3, 4})) == 41943041);   // rotl(1 + 4, 23) + 1
static_assert(first_output(wsl::random::philox4x64(7)) == first_output(wsl::random::philox4x64(7)));

int main()
{
    // splitmix64
    assert(same_outputs(wsl::random::splitmix64(0), v::splitmix64_seed_0));
    assert(same_outputs(wsl::random::splitmix64(1234567), v::splitmix64_seed_1234567));
    {
        wsl::random::splitmix64 a(5), b(5);
        a.discard(1000);
        for (int i = 0; i < 1000; ++i)
            b();
        assert(a == b && a() == b());
    }

    // xoshiro256++ : état brut, initialisation par splitmix64, sauts (vérifiés par GF(2)).
    assert(same_outputs(wsl::random::xoshiro256pp({1, 2, 3, 4}), v::xoshiro_state_1234));
    {
        wsl::random::xoshiro256pp const seeded(42);
        auto const& s = seeded.state();
        for (int i = 0; i < 4; ++i)
            assert(s[static_cast<std::size_t>(i)] == v::xoshiro_seed_42_state[i]);
        assert(same_outputs(seeded, v::xoshiro_seed_42));
    }
    {
        wsl::random::xoshiro256pp g({1, 2, 3, 4}), h({1, 2, 3, 4});
        g.jump();
        h.long_jump();
        for (int i = 0; i < 4; ++i)
        {
            assert(g.state()[static_cast<std::size_t>(i)] == v::xoshiro_1234_jump_state[i]);
            assert(h.state()[static_cast<std::size_t>(i)] == v::xoshiro_1234_long_jump_state[i]);
        }
    }

    // pcg64 : état brut de NumPy, avance de 10^20, initialisation (42, 54).
    wsl::random::pcg64::raw_state const raw{u128(v::pcg64_state), u128(v::pcg64_inc)};
    assert(same_outputs(wsl::random::pcg64(raw), v::pcg64_outputs));
    {
        wsl::random::pcg64 g(raw);
        wsl::random::detail::uint128 delta = 1;
        for (int i = 0; i < 20; ++i)
            delta *= 10;
        g.advance(delta);
        assert(same_outputs(g, v::pcg64_after_advance_1e20));
        wsl::random::pcg64 a(raw), b(raw);
        a.discard(777);
        for (int i = 0; i < 777; ++i)
            b();
        assert(a == b);
    }
    assert(same_outputs(wsl::random::pcg64(42, 54), v::pcg64_seed_42_stream_54));

    // philox4x64 : clé et compteur, propagation de retenue, accès direct.
    wsl::random::philox4x64::key_type const key{v::philox_key[0], v::philox_key[1]};
    assert(same_outputs(wsl::random::philox4x64(key, {v::philox_counter[0], v::philox_counter[1], v::philox_counter[2], v::philox_counter[3]}),
                        v::philox_outputs));
    assert(same_outputs(wsl::random::philox4x64(key, {v::philox_carry_counter[0], v::philox_carry_counter[1],
                                              v::philox_carry_counter[2], v::philox_carry_counter[3]}),
                        v::philox_carry_outputs));
    {
        for (unsigned long long n : {0ull, 1ull, 3ull, 4ull, 5ull, 17ull, 1000ull})
        {
            wsl::random::philox4x64 a(9), b(9);
            b();                                                            // discard depuis un bloc entamé
            a();
            a.discard(n);
            for (unsigned long long i = 0; i < n; ++i)
                b();
            assert(a == b && a() == b());
        }
        auto const block = wsl::random::philox4x64::block(key, {5, 0, 0, 0});
        assert(block[0] == v::philox_outputs[0] && block[3] == v::philox_outputs[3]);
    }

    // Avec <random> : un moteur wsl alimente une distribution standard.
    wsl::random::xoshiro256pp g(1);
    std::uniform_int_distribution<int> die(1, 6);
    int const roll = die(g);
    assert(roll >= 1 && roll <= 6);
    return 0;
}
