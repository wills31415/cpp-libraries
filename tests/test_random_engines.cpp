// Moteurs de wsl::random : sorties comparées à des sources indépendantes de ce
// code (tests/random_engine_vectors.hpp : NumPy, matrices sur GF(2)), propriétés
// d'interface, sauts et avances.

#include <wsl/random.hpp>

#include "random_engine_vectors.hpp"

#include <cassert>
#include <cstdint>
#include <random>
#include <span>
#include <iomanip>
#include <sstream>
#include <string>

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

// RandomNumberEngine : E(), E(s), E(q), seed(), seed(s), seed(q), ==, <<, >>.
template <class E>
void check_engine(std::uint64_t expected_from_seq)
{
    static_assert(std::uniform_random_bit_generator<E> && std::is_nothrow_default_constructible_v<E>);
    E const d = {};                                                         // E() non explicite
    assert(d == E(0));
    std::seed_seq q1{1, 2, 3}, q2{1, 2, 3};
    E a(q1);                                                                // valeur portable : seed_seq est spécifiée
    E b;
    b.seed(q2);
    assert(a == b && a() == expected_from_seq && b() == expected_from_seq);
    b.seed();
    assert(b == E());
    b.seed(7);
    assert(b == E(7));
    // aller-retour par un flux, en tout point d'un bloc (philox)
    for (int k = 0; k < 6; ++k)
    {
        std::stringstream ss;
        ss << std::hex << std::setfill('#') << a;                           // format du flux rétabli
        assert((ss.flags() & std::ios::basefield) == std::ios::hex && ss.fill() == '#');
        E c;
        ss >> c;
        assert(!ss.fail() && c == a);
        for (int i = 0; i < 5; ++i)
            assert(c() == a());
        a.discard(static_cast<unsigned long long>(k));
    }
    std::wstringstream ws;
    ws << a;
    E w;
    ws >> w;
    assert(w == a);
    // lecture invalide : moteur intact, failbit
    std::stringstream bad("abc");
    E e = a;
    bad >> e;
    assert(bad.fail() && e == a);
}

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

    // ---- RandomNumberEngine complets -------------------------------------------------
    {
        check_engine<wsl::random::splitmix64>(0x7f48bd39a96b09ab);
        check_engine<wsl::random::xoshiro256pp>(0x7adf7674778b678d);
        check_engine<wsl::random::pcg64>(0x93e5d4564e9e31fb);
        check_engine<wsl::random::philox4x64>(0x02accfa08a6e1d7e);
        // états invalides refusés : xoshiro nul, incrément pcg pair, position philox > 4
        std::stringstream z("0 0 0 0"), even("1 2 3 4"), position("1 2 3 4 5 6 7");
        wsl::random::xoshiro256pp x(5);
        wsl::random::pcg64 c(5);
        wsl::random::philox4x64 f(5);
        z >> x;
        even >> c;
        position >> f;
        assert(z.fail() && x == wsl::random::xoshiro256pp(5) && even.fail() && c == wsl::random::pcg64(5) && position.fail() && f == wsl::random::philox4x64(5));
        // philox : discard en O(1), retenue sur les quatre mots du compteur
        wsl::random::philox4x64 s1(9), s2(9);
        s1();
        s2();
        s1.discard(4003);
        for (int i = 0; i < 4003; ++i)
            s2();
        assert(s1 == s2 && s1() == s2());
        wsl::random::philox4x64 h({1, 2}, {~0ull, ~0ull, 5, 0}), k({1, 2}, {3, 0, 6, 0});
        h.discard(16);                                                      // 2^128 − 1 + 4 blocs
        assert(h == k && h() == k());
    }
    return 0;
}
