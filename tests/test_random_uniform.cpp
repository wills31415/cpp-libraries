// Primitives uniformes de wsl::random : cas limites par des générateurs factices
// (mots choisis), lois vérifiées statistiquement (graines fixes : aucun test
// aléatoire), et empreintes des suites rendues, qui doivent être les mêmes sur
// toute plateforme.

#include <wsl/random.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <vector>

namespace
{
    // Générateur factice : rend les mots donnés, puis le dernier indéfiniment.
    struct scripted
    {
        using result_type = std::uint64_t;
        static constexpr result_type min() { return 0; }
        static constexpr result_type max() { return ~result_type{0}; }
        std::vector<std::uint64_t> words;
        std::size_t next = 0;
        result_type operator()() { return words[std::min(next++, words.size() - 1)]; }
    };

    // χ² des effectifs observés contre des effectifs attendus égaux.
    double chi_square(std::vector<double> const& counts)
    {
        double total = 0;
        for (double c : counts)
            total += c;
        double const expected = total / double(counts.size());
        double chi = 0;
        for (double c : counts)
            chi += (c - expected) * (c - expected) / expected;
        return chi;
    }

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

    bool check(char const* name, std::uint64_t actual, std::uint64_t expected)
    {
        bool const ok = actual == expected;
        std::printf("%-22s %016llx %s\n", name, static_cast<unsigned long long>(actual),
                    ok ? "ok" : "DIFFÉRENTE de l'empreinte attendue");
        std::fflush(stdout);
        return ok;
    }

    template <class F>
    std::uint64_t fingerprint(F f, int count = 1 << 16)
    {
        fnv1a h;
        for (int i = 0; i < count; ++i)
            h.add(f());
        return h.value;
    }
}

// Évaluation constante : un mélange calculé à la compilation.
constexpr std::array<int, 6> shuffled_at_compile_time()
{
    std::array<int, 6> a{0, 1, 2, 3, 4, 5};
    wsl::random::xoshiro256pp g(3);
    wsl::random::shuffle(a, g);
    return a;
}
static_assert(std::ranges::is_permutation(shuffled_at_compile_time(), std::array{0, 1, 2, 3, 4, 5}));
static_assert([] { wsl::random::pcg64 g(1); return wsl::random::uniform_below(10, g) < 10; }());

int main()
{
    // ---- random_bits : générateurs de 32 bits et de plage non binaire --------
    {
        std::mt19937 a(5), b(5);
        std::uint64_t const high = b(), low = b();
        assert(wsl::random::random_bits(a) == ((high << 32) | low));                // premier appel en poids fort
        std::minstd_rand m(1);                                              // plage [1, 2^31 − 2] : rejet
        double ones = 0;
        for (int i = 0; i < 20000; ++i)
            ones += std::popcount(wsl::random::random_bits(m));
        assert(std::abs(ones / (20000.0 * 64) - 0.5) < 0.005);
    }

    wsl::random::xoshiro256pp g(2026);

    // ---- uniform_below ----------------------------------------------------------
    {
        for (int i = 0; i < 100; ++i)
            assert(wsl::random::uniform_below(1, g) == 0);
        std::vector<double> counts(7);
        for (int i = 0; i < 70000; ++i)
            ++counts[wsl::random::uniform_below(7, g)];
        assert(chi_square(counts) < 30);                                    // 6 degrés de liberté
        // n = 2^63 + 1 : un simple modulo donnerait les petites valeurs deux fois
        // plus souvent ; Lemire, exactement la moitié des tirages sous n/2.
        std::uint64_t const n = (std::uint64_t{1} << 63) + 1;
        int below_half = 0;
        for (int i = 0; i < 100000; ++i)
            below_half += wsl::random::uniform_below(n, g) < n / 2;
        assert(std::abs(below_half / 100000.0 - 0.5) < 0.006);
        bool thrown = false;
        try { (void)wsl::random::uniform_below(0, g); } catch (std::domain_error const&) { thrown = true; }
        assert(thrown);
    }

    // ---- uniform_int ------------------------------------------------------------
    {
        std::array<int, 256> seen{};
        for (int i = 0; i < 100000; ++i)
            ++seen[static_cast<std::size_t>(wsl::random::uniform_int<std::int8_t>(-128, 127, g) + 128)];
        assert(std::ranges::all_of(seen, [](int c) { return c > 250; }));
        for (int i = 0; i < 1000; ++i)
        {
            auto const x = wsl::random::uniform_int(-5LL, -3LL, g);
            assert(x >= -5 && x <= -3);
            (void)wsl::random::uniform_int(std::numeric_limits<std::int64_t>::min(), std::numeric_limits<std::int64_t>::max(), g);
        }
        assert(wsl::random::uniform_int(42u, 42u, g) == 42u);
        bool thrown = false;
        try { (void)wsl::random::uniform_int(3, 2, g); } catch (std::domain_error const&) { thrown = true; }
        assert(thrown);
    }

    // ---- uniform_real : cas limites par mots choisis ----------------------------
    {
        scripted s{{0x8000000000000000}};                                   // U = 0,1000… : 1/2 exactement
        assert(wsl::random::uniform_real(s) == 0.5);
        scripted half_ones{{0xc000000000000000}};                           // 0,11 : 3/4
        assert(wsl::random::uniform_real(half_ones) == 0.75);
        scripted most{{~std::uint64_t{0}}};                                 // 0,111… : le plus grand double < 1
        assert(wsl::random::uniform_real(most) == std::nextafter(1.0, 0.0));
        // 16 mots nuls puis un 1 en tête : U ∈ [2^-1025, 2^-1024), sous-normal.
        std::vector<std::uint64_t> words(16, 0);
        words.push_back(0x8000000000000000);
        scripted tiny{words};
        assert(wsl::random::uniform_real(tiny) == std::ldexp(1.0, -1025));
        scripted zero{{0}};                                                 // U < 2^-1074 : 0
        assert(wsl::random::uniform_real(zero) == 0.0 && zero.next == 17);
        scripted f{{0x8000000000000000}};
        assert(wsl::random::uniform_real<float>(f) == 0.5f);
        scripted fz{{0}};
        assert(wsl::random::uniform_real<float>(fz) == 0.0f && fz.next == 3);       // 2^-192 < 2^-149
    }

    // ---- uniform_real : loi --------------------------------------------------------
    {
        constexpr int count = 400000;
        std::array<int, 12> below{};                                        // U < 2^-k
        double sum = 0, largest = 0;
        for (int i = 0; i < count; ++i)
        {
            double const u = wsl::random::uniform_real(g);
            assert(u >= 0 && u < 1);
            sum += u;
            largest = std::max(largest, u);
            for (int k = 0; k < 12; ++k)
                below[static_cast<std::size_t>(k)] += u < std::ldexp(1.0, -k - 1);
        }
        assert(std::abs(sum / count - 0.5) < 0.003);
        for (int k = 0; k < 12; ++k)
        {
            double const expected = count * std::ldexp(1.0, -k - 1);
            assert(std::abs(below[static_cast<std::size_t>(k)] - expected) < 5 * std::sqrt(expected) + 2);
        }
        for (int i = 0; i < 100000; ++i)
        {
            float const u = wsl::random::uniform_real<float>(g);
            double const v = wsl::random::uniform_real_fast(g);
            float const w = wsl::random::uniform_real_fast<float>(g);
            assert(u >= 0 && u < 1 && v >= 0 && v < 1 && w >= 0 && w < 1);
            assert(v * 0x1p53 == std::floor(v * 0x1p53));                   // multiples de 2^-53
        }
    }

    // ---- bernoulli ----------------------------------------------------------------
    {
        int hits = 0;
        for (int i = 0; i < 200000; ++i)
            hits += wsl::random::bernoulli_trial(0.3, g);
        assert(std::abs(hits / 200000.0 - 0.3) < 0.004);
        assert(!wsl::random::bernoulli_trial(0.0, g) && wsl::random::bernoulli_trial(1.0, g));
        scripted zeros{{0}};
        assert(wsl::random::bernoulli_trial(0x1p-1074, zeros));                           // U = 0 < p
        scripted ones{{~std::uint64_t{0}}};
        assert(!wsl::random::bernoulli_trial(std::nextafter(1.0, 0.0), ones));            // U = 1 − 2^-64·… >= p
        scripted equal{{0xc000000000000000, 0}};                            // U = 0,11000… = p
        assert(!wsl::random::bernoulli_trial(0.75, equal) && equal.next == 1);            // U >= p, décidé sans autre mot
        scripted below_p{{0xbfffffffffffffff}};
        assert(wsl::random::bernoulli_trial(0.75, below_p));
        int third = 0;
        for (int i = 0; i < 300000; ++i)
            third += wsl::random::bernoulli_trial(1, 3, g);
        assert(std::abs(third / 300000.0 - 1.0 / 3) < 0.004);
        bool thrown = false;
        try { (void)wsl::random::bernoulli_trial(1.5, g); } catch (std::domain_error const&) { thrown = true; }
        assert(thrown);
    }

    // ---- shuffle -----------------------------------------------------------------
    {
        std::map<std::array<int, 4>, double> permutations;
        for (int i = 0; i < 240000; ++i)
        {
            std::array<int, 4> a{0, 1, 2, 3};
            wsl::random::shuffle(a, g);
            ++permutations[a];
        }
        assert(permutations.size() == 24);
        std::vector<double> counts;
        for (auto const& [p, c] : permutations)
            counts.push_back(c);
        assert(chi_square(counts) < 50);                                    // 23 degrés de liberté
        std::vector<int> empty;
        wsl::random::shuffle(empty, g);
        std::vector<int> single{7};
        wsl::random::shuffle(single, g);
        assert(single[0] == 7);
    }

    // ---- empreintes : mêmes suites sur toute plateforme --------------------------
    bool ok = true;
    {
        wsl::random::xoshiro256pp x(7);
        ok &= check("uniform_below", fingerprint([&] { return wsl::random::uniform_below(1000003, x); }), 0xc31c297194c07ec7);
        ok &= check("uniform_int", fingerprint([&] { return static_cast<std::uint64_t>(wsl::random::uniform_int(-50, 50, x)); }), 0xe68e1093afb2279b);
        ok &= check("uniform_real double", fingerprint([&] { return std::bit_cast<std::uint64_t>(wsl::random::uniform_real(x)); }), 0x3dac0cd7502899fb);
        ok &= check("uniform_real float", fingerprint([&] { return std::uint64_t{std::bit_cast<std::uint32_t>(wsl::random::uniform_real<float>(x))}; }), 0x4713094695c22526);
        ok &= check("uniform_real_fast", fingerprint([&] { return std::bit_cast<std::uint64_t>(wsl::random::uniform_real_fast(x)); }), 0xf767bee0632d4554);
        ok &= check("bernoulli", fingerprint([&] { return std::uint64_t{wsl::random::bernoulli_trial(0.3, x)}; }), 0xb356cd399bd84ee4);
        ok &= check("shuffle", fingerprint([&] {
            std::array<std::uint64_t, 52> deck{};
            for (std::size_t i = 0; i < deck.size(); ++i)
                deck[i] = i;
            wsl::random::shuffle(deck, x);
            return deck[0] | deck[1] << 8 | deck[51] << 16;
        }, 4096), 0x0713aac5d0827c70);
        std::mt19937 mt(7);
        std::minstd_rand ms(7);
        ok &= check("random_bits mt19937", fingerprint([&] { return wsl::random::random_bits(mt); }), 0x89135d00cef863bb);
        ok &= check("random_bits minstd", fingerprint([&] { return wsl::random::random_bits(ms); }), 0xad7600aa69d20a1c);
        wsl::random::philox4x64 ph(7);
        wsl::random::pcg64 pc(7);
        ok &= check("uniform_real philox", fingerprint([&] { return std::bit_cast<std::uint64_t>(wsl::random::uniform_real(ph)); }), 0x7b280f1d68239a86);
        ok &= check("uniform_below pcg64", fingerprint([&] { return wsl::random::uniform_below(6, pc); }), 0xbe14811296869dc7);
    }
    assert(ok);
    return 0;
}
