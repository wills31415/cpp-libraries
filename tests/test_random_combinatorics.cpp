// Objets combinatoires de wsl::random : permutations (uniformes, Ewens, Mallows),
// sous-ensembles, combinaisons, mots, partitions d'entiers. Masses sommant à 1,
// χ² contre la loi exacte sur tous les objets d'une petite taille (graines fixes),
// choix exacts à poids flottants, empreintes.

#include <wsl/random.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <vector>

namespace
{
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

    template <class F>
    std::uint64_t fingerprint(F f, int count = 1 << 12)
    {
        fnv1a h;
        for (int i = 0; i < count; ++i)
            f(h);
        return h.value;
    }

    bool check(char const* name, std::uint64_t actual, std::uint64_t expected)
    {
        bool const ok = actual == expected;
        std::printf("%-22s %016llx %s\n", name, static_cast<unsigned long long>(actual), ok ? "ok" : "DIFFÉRENTE de l'empreinte attendue");
        std::fflush(stdout);
        return ok;
    }

    template <class F>
    bool throws_domain_error(F f)
    {
        try { f(); } catch (std::domain_error const&) { return true; }
        return false;
    }

    bool close(double a, double b, double relative) { return std::abs(a - b) <= relative * std::abs(b); }

    // χ² de draws tirages de law contre sa masse, sur tous les objets de objects ;
    // écart réduit (χ² − ddl)/√(2 ddl), < 4 attendu.
    template <class L, class G, class V>
    double chi_square(L const& law, G& g, std::vector<V> const& objects, int draws)
    {
        std::map<V, int> counts;
        for (int i = 0; i < draws; ++i)
            ++counts[sample(law, g)];
        double chi = 0, total = 0;
        for (auto const& o : objects)
        {
            double const e = draws * pmf(law, o);
            total += pmf(law, o);
            auto const it = counts.find(o);
            double const c = it == counts.end() ? 0 : it->second;
            chi += (c - e) * (c - e) / e;
        }
        assert(close(total, 1, 1e-13));
        std::size_t seen = 0;
        for (auto const& [o, c] : counts)
            seen += static_cast<std::size_t>(c);
        assert(seen == static_cast<std::size_t>(draws));                     // aucun objet hors de la liste
        double const dof = double(objects.size()) - 1;
        return (chi - dof) / std::sqrt(2 * dof);
    }

    std::vector<wsl::permutation> all_permutations(std::size_t n)
    {
        std::vector<wsl::permutation> r;
        wsl::permutation p(n);
        do
            r.push_back(p);
        while (p.next());
        return r;
    }

    std::vector<std::vector<std::size_t>> all_subsets(std::size_t n, std::size_t k = SIZE_MAX)
    {
        std::vector<std::vector<std::size_t>> r;
        for (std::uint64_t m = 0; m < (std::uint64_t{1} << n); ++m)
        {
            std::vector<std::size_t> s;
            for (std::size_t i = 0; i < n; ++i)
                if ((m >> i) & 1)
                    s.push_back(i);
            if (k == SIZE_MAX || s.size() == k)
                r.push_back(s);
        }
        return r;
    }

    void partitions(std::size_t n, std::size_t largest, std::vector<std::size_t>& current, std::vector<std::vector<std::size_t>>& out)
    {
        if (n == 0)
        {
            out.push_back(current);
            return;
        }
        for (std::size_t d = std::min(n, largest); d >= 1; --d)
        {
            current.push_back(d);
            partitions(n - d, d, current, out);
            current.pop_back();
        }
    }
}

static_assert(wsl::probability_law<wsl::ewens> && !wsl::has_support<wsl::ewens> && wsl::has_pmf<wsl::mallows>);
static_assert(wsl::uniform_permutation::sampling == wsl::exactness::exact && wsl::mallows::sampling == wsl::exactness::exact);
static_assert(wsl::has_support<wsl::normal<>> && wsl::has_pmf<wsl::uniform_integer_partition>);

int main()
{
    wsl::xoshiro256pp g(20261118);
    constexpr int draws = 60000;

    // ---- Choix exacts à poids flottants ---------------------------------------------
    {
        int hits = 0;                                                        // 128 bits : 0,1/(0,1 + 3)
        for (int i = 0; i < draws; ++i)
            hits += wsl::detail::choose(g, 0.1, 3.0);
        assert(std::abs(hits / double(draws) - 0.1 / 3.1) < 0.003);
        hits = 0;                                                            // grands entiers : 2^-600 contre 2^600
        for (int i = 0; i < 2000; ++i)
            hits += wsl::detail::choose(g, 0x1p-600, 0x1p600);
        assert(hits == 0);
        hits = 0;                                                            // a/b, b non dyadique
        for (int i = 0; i < draws; ++i)
            hits += wsl::detail::fraction(g, 1.0, 2.5);
        assert(std::abs(hits / double(draws) - 0.4) < 0.006);
        assert(wsl::detail::choose(g, 1.0, 0.0) && !wsl::detail::choose(g, 0.0, 1.0) && wsl::detail::fraction(g, 2.0, 2.0));
    }

    // ---- Permutations ----------------------------------------------------------------
    {
        auto const s4 = all_permutations(4);
        assert(chi_square(wsl::uniform_permutation({.size = 4}), g, s4, draws) < 4);
        for (double theta : {0.3, 2.5, 40.0})
            assert(chi_square(wsl::ewens({.size = 4, .theta = theta}), g, s4, draws) < 4);
        for (double q : {0.01, 0.4, 0.999, 1.0, 2.5, 300.0})                // 0,999 : uniforme acceptée ; sinon géométrique
            assert(chi_square(wsl::mallows({.size = 4, .q = q}), g, s4, draws) < 4);
        auto const s6 = all_permutations(6);
        wsl::ewens const E({.size = 6, .theta = 2.5});
        wsl::mallows const M({.size = 6, .q = 0.4});
        double te = 0, tm = 0, ce = 0, im = 0;
        for (auto const& p : s6)
        {
            te += pmf(E, p);
            tm += pmf(M, p);
            ce += pmf(E, p) * double(p.cycle_count());
            im += pmf(M, p) * double(p.inversions());
        }
        assert(close(te, 1, 1e-14) && close(tm, 1, 1e-14) && close(ce, E.mean_cycles(), 1e-14) && close(im, M.mean_inversions(), 1e-13));
        assert(close(wsl::mallows({.size = 6, .q = 2.5}).mean_inversions(), 15 - M.mean_inversions(), 1e-13));   // q ↔ 1/q : inv ↔ 15 − inv
        // Ewens θ = 1 : uniforme.
        assert(close(pmf(wsl::ewens({.size = 6, .theta = 1.0}), s6[17]), 1.0 / 720, 1e-15) && close(pmf(wsl::uniform_permutation({.size = 6}), s6[3]), 1.0 / 720, 1e-15));
        // Grandes tailles : validité et statistiques simples.
        wsl::ewens const big({.size = 2000, .theta = 3.0});
        double cycles = 0;
        for (int i = 0; i < 200; ++i)
            cycles += double(sample(big, g).cycle_count());
        assert(std::abs(cycles / 200 - big.mean_cycles()) < 0.6);
        wsl::mallows const near({.size = 300, .q = 1 - 1e-4});
        double inv = 0;
        for (int i = 0; i < 200; ++i)
            inv += double(sample(near, g).inversions());
        assert(std::abs(inv / 200 - near.mean_inversions()) / near.mean_inversions() < 0.02);
        assert(pmf(E, wsl::permutation(5)) == 0 && log_pmf(M, wsl::permutation(7)) == -INFINITY);
        assert(throws_domain_error([] { wsl::ewens X({.size = 3, .theta = 0.0}); }) && throws_domain_error([] { wsl::mallows X({.size = 3, .q = -1.0}); }));
    }

    // ---- Sous-ensembles, combinaisons, mots -------------------------------------------
    {
        assert(chi_square(wsl::uniform_subset({.size = 4}), g, all_subsets(4), draws) < 4);
        assert(chi_square(wsl::bernoulli_subset({.size = 4, .p = 0.3}), g, all_subsets(4), draws) < 4);
        assert(chi_square(wsl::uniform_combination({.size = 6, .count = 3}), g, all_subsets(6, 3), draws) < 4);
        std::vector<std::vector<std::size_t>> words;
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b)
                for (std::size_t c = 0; c < 3; ++c)
                    words.push_back({a, b, c});
        assert(chi_square(wsl::uniform_word({.length = 3, .alphabet = 3}), g, words, draws) < 4);
        wsl::uniform_combination const C({.size = 1000000, .count = 50});
        auto const s = sample(C, g);
        assert(s.size() == 50 && pmf(C, s) > 0 && close(C.entropy(), std::lgamma(1000001.0) - std::lgamma(51.0) - std::lgamma(999951.0), 1e-12));
        assert(pmf(C, {1, 0}) == 0 && pmf(wsl::uniform_subset({.size = 3}), {0, 3}) == 0 && pmf(wsl::uniform_word({.length = 2, .alphabet = 2}), {0, 2}) == 0);
        assert(throws_domain_error([] { wsl::uniform_combination X({.size = 3, .count = 4}); }));
    }

    // ---- Partitions d'entiers ------------------------------------------------------------
    {
        std::vector<std::vector<std::size_t>> p8;
        std::vector<std::size_t> current;
        partitions(8, 8, current, p8);
        assert(p8.size() == 22);
        assert(chi_square(wsl::uniform_integer_partition({.n = 8}), g, p8, draws) < 4);
        assert(wsl::uniform_integer_partition({.n = 100}).count() == 190569292);
        assert(close(wsl::uniform_integer_partition({.n = 1000}).count(), 24061467864032622473692149727991.0, 1e-15));
        wsl::uniform_integer_partition const P({.n = 500});
        auto const q = sample(P, g);
        std::size_t total = 0;
        for (std::size_t v : q)
            total += v;
        assert(total == 500 && std::is_sorted(q.rbegin(), q.rend()) && pmf(P, q) > 0 && pmf(P, {3, 4}) == 0);
    }

    // ---- Empreintes ----------------------------------------------------------------
    bool ok = true;
    {
        wsl::xoshiro256pp x(42);
        auto const add_perm = [](fnv1a& h, wsl::permutation const& p) { for (std::size_t v : p) h.add(v); };
        auto const add_vec = [](fnv1a& h, std::vector<std::size_t> const& v) { h.add(v.size()); for (std::size_t e : v) h.add(e); };
        ok &= check("uniform_permutation", fingerprint([&](fnv1a& h) { add_perm(h, sample(wsl::uniform_permutation({.size = 20}), x)); }), 0xc24e2f0c1702f325);
        ok &= check("ewens", fingerprint([&](fnv1a& h) { add_perm(h, sample(wsl::ewens({.size = 20, .theta = 1.7}), x)); }), 0x16e0ba54eb024465);
        ok &= check("mallows", fingerprint([&](fnv1a& h) { add_perm(h, sample(wsl::mallows({.size = 20, .q = 0.8}), x)); }), 0x1dbf2a90e4629485);
        ok &= check("mallows (q > 1)", fingerprint([&](fnv1a& h) { add_perm(h, sample(wsl::mallows({.size = 20, .q = 1.05}), x)); }), 0xd7a97235d3517c65);
        ok &= check("subsets", fingerprint([&](fnv1a& h) {
            add_vec(h, sample(wsl::uniform_subset({.size = 70}), x));
            add_vec(h, sample(wsl::bernoulli_subset({.size = 30, .p = 0.2}), x));
            add_vec(h, sample(wsl::uniform_combination({.size = 100, .count = 7}), x));
            add_vec(h, sample(wsl::uniform_word({.length = 9, .alphabet = 26}), x));
        }), 0x0ed6bf2a62e76d9b);
        wsl::uniform_integer_partition const P({.n = 60});
        ok &= check("integer_partition", fingerprint([&](fnv1a& h) { add_vec(h, sample(P, x)); }, 1 << 10), 0x8d0ee620a5ae13cc);
    }
    assert(ok);
    std::printf("test_random_combinatorics : ok\n");
}
