// Objets combinatoires de wsl::random : permutations (uniformes, Ewens, Mallows),
// sous-ensembles, combinaisons, mots, partitions d'entiers. Masses sommant à 1,
// χ² contre la loi exacte sur tous les objets d'une petite taille (graines fixes),
// choix exacts à poids flottants, empreintes.

#include <wsl/probability.hpp>
#include <wsl/random.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <vector>

// Les noms de l'ancien wsl::random, réunis (refonte, étape 5a) : moteurs, primitives, lois,
// processus, mesures de référence.
namespace wr
{
    using namespace wsl::random;
    using namespace wsl::random::engines;
    using namespace wsl::probability;
    using namespace wsl::probability::laws;
    using namespace wsl::probability::stochastic_processes;
    using namespace wsl::probability::stochastic_processes::processes;
    using namespace wsl::measure_theory;
    namespace detail = wsl::probability::detail;

    // Le tirage natif, l'ancien défaut (exact pour une loi exacte, son seul tirage sinon) :
    // celui que mesurent les empreintes de ces tests.
    template <class X, class G>
    constexpr auto native(X const& x, G& g)
    {
        if constexpr (requires { wsl::probability::detail::native_sample(x, g); })
            return wsl::probability::detail::native_sample(x, g);
        else if constexpr (requires { x.law(); })
            return wsl::probability::detail::native_sample(x.law(), g);
        else
            return x(g, wsl::exact);
    }
    template <class X, class G>
    constexpr auto native(X const& x, std::size_t n, G& g)
    {
        std::vector<decltype(native(x, g))> out;
        out.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            out.push_back(native(x, g));
        return out;
    }
}


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
    template <class L, class G, class V, class P>
    double chi_square(L const& law, G& g, std::vector<V> const& objects, int draws, P probability);

    template <class L, class G, class V>
    double chi_square(L const& law, G& g, std::vector<V> const& objects, int draws)
    {
        return chi_square(law, g, objects, draws, [&](V const& o) { return wr::pmf(law, o); });
    }

    template <class L, class G, class V, class P>
    double chi_square(L const& law, G& g, std::vector<V> const& objects, int draws, P probability)
    {
        std::map<V, int> counts;
        for (int i = 0; i < draws; ++i)
            ++counts[wr::native(law, g)];
        double chi = 0, total = 0;
        for (auto const& o : objects)
        {
            double const e = draws * probability(o);
            total += probability(o);
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

    // Graphes sur n sommets, arêtes parmi les C(n, 2) paires (ordre lexicographique).
    std::vector<wr::edge_list> all_graphs(std::size_t n, std::size_t edges = SIZE_MAX)
    {
        wr::edge_list pairs;
        for (std::size_t u = 0; u < n; ++u)
            for (std::size_t v = u + 1; v < n; ++v)
                pairs.emplace_back(u, v);
        std::vector<wr::edge_list> r;
        for (std::uint64_t m = 0; m < (std::uint64_t{1} << pairs.size()); ++m)
        {
            wr::edge_list e;
            for (std::size_t i = 0; i < pairs.size(); ++i)
                if ((m >> i) & 1)
                    e.push_back(pairs[i]);
            if (edges == SIZE_MAX || e.size() == edges)
                r.push_back(e);
        }
        return r;
    }

    std::vector<std::size_t> degrees(wr::edge_list const& e, std::size_t n)
    {
        std::vector<std::size_t> d(n, 0);
        for (auto const& [u, v] : e)
        {
            ++d[u];
            ++d[v];
        }
        return d;
    }

    // Arbres binaires à n nœuds, numérotés en préordre.
    std::vector<wr::binary_tree> all_binary_trees(std::size_t n)
    {
        constexpr std::size_t none = wr::binary_tree::none;
        if (n == 0)
            return {wr::binary_tree{}};
        std::vector<wr::binary_tree> r;
        for (std::size_t l = 0; l < n; ++l)
            for (auto const& left : all_binary_trees(l))
                for (auto const& right : all_binary_trees(n - 1 - l))
                {
                    wr::binary_tree t;
                    t.children.push_back({l > 0 ? 1 : none, n - 1 - l > 0 ? 1 + l : none});
                    for (auto c : left.children)
                        t.children.push_back({c[0] == none ? none : c[0] + 1, c[1] == none ? none : c[1] + 1});
                    for (auto c : right.children)
                        t.children.push_back({c[0] == none ? none : c[0] + 1 + l, c[1] == none ? none : c[1] + 1 + l});
                    r.push_back(t);
                }
        return r;
    }

    // Mots à croissance restreinte de longueur n (partitions d'ensemble).
    void restricted_growth(std::size_t n, std::size_t blocks, std::vector<std::size_t>& a, std::vector<std::vector<std::size_t>>& out)
    {
        if (a.size() == n)
        {
            out.push_back(a);
            return;
        }
        for (std::size_t v = 0; v <= blocks; ++v)
        {
            a.push_back(v);
            restricted_growth(n, v == blocks ? blocks + 1 : blocks, a, out);
            a.pop_back();
        }
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

static_assert(wr::probability_law<wr::ewens> && !wr::has_support<wr::ewens> && wr::has_pmf<wr::mallows>);
static_assert(wr::uniform_permutation::sampling == wr::exactness::exact && wr::mallows::sampling == wr::exactness::exact);
static_assert(wr::has_support<wr::normal<>> && wr::has_pmf<wr::uniform_integer_partition>);

int main()
{
    wr::xoshiro256pp g(20261118);
    constexpr int draws = 60000;

    // ---- Choix exacts à poids flottants ---------------------------------------------
    {
        int hits = 0;                                                        // 128 bits : 0,1/(0,1 + 3)
        for (int i = 0; i < draws; ++i)
            hits += wsl::probability::detail::choose(g, 0.1, 3.0);
        assert(std::abs(hits / double(draws) - 0.1 / 3.1) < 0.003);
        hits = 0;                                                            // grands entiers : 2^-600 contre 2^600
        for (int i = 0; i < 2000; ++i)
            hits += wsl::probability::detail::choose(g, 0x1p-600, 0x1p600);
        assert(hits == 0);
        hits = 0;                                                            // a/b, b non dyadique
        for (int i = 0; i < draws; ++i)
            hits += wsl::probability::detail::fraction(g, 1.0, 2.5);
        assert(std::abs(hits / double(draws) - 0.4) < 0.006);
        assert(wsl::probability::detail::choose(g, 1.0, 0.0) && !wsl::probability::detail::choose(g, 0.0, 1.0) && wsl::probability::detail::fraction(g, 2.0, 2.0));
    }

    // ---- Permutations ----------------------------------------------------------------
    {
        auto const s4 = all_permutations(4);
        assert(chi_square(wr::uniform_permutation({.size = 4}), g, s4, draws) < 4);
        for (double theta : {0.3, 2.5, 40.0})
            assert(chi_square(wr::ewens({.size = 4, .theta = theta}), g, s4, draws) < 4);
        for (double q : {0.01, 0.4, 0.999, 1.0, 2.5, 300.0})                // 0,999 : uniforme acceptée ; sinon géométrique
            assert(chi_square(wr::mallows({.size = 4, .q = q}), g, s4, draws) < 4);
        auto const s6 = all_permutations(6);
        wr::ewens const E({.size = 6, .theta = 2.5});
        wr::mallows const M({.size = 6, .q = 0.4});
        double te = 0, tm = 0, ce = 0, im = 0;
        for (auto const& p : s6)
        {
            te += wr::pmf(E, p);
            tm += wr::pmf(M, p);
            ce += wr::pmf(E, p) * double(p.cycle_count());
            im += wr::pmf(M, p) * double(p.inversions());
        }
        assert(close(te, 1, 1e-14) && close(tm, 1, 1e-14) && close(ce, E.mean_cycles(), 1e-14) && close(im, M.mean_inversions(), 1e-13));
        assert(close(wr::mallows({.size = 6, .q = 2.5}).mean_inversions(), 15 - M.mean_inversions(), 1e-13));   // q ↔ 1/q : inv ↔ 15 − inv
        // Ewens θ = 1 : uniforme.
        assert(close(wr::pmf(wr::ewens({.size = 6, .theta = 1.0}), s6[17]), 1.0 / 720, 1e-15) && close(wr::pmf(wr::uniform_permutation({.size = 6}), s6[3]), 1.0 / 720, 1e-15));
        // Grandes tailles : validité et statistiques simples.
        wr::ewens const big({.size = 2000, .theta = 3.0});
        double cycles = 0;
        for (int i = 0; i < 200; ++i)
            cycles += double(wr::native(big, g).cycle_count());
        assert(std::abs(cycles / 200 - big.mean_cycles()) < 0.6);
        wr::mallows const near({.size = 300, .q = 1 - 1e-4});
        double inv = 0;
        for (int i = 0; i < 200; ++i)
            inv += double(wr::native(near, g).inversions());
        assert(std::abs(inv / 200 - near.mean_inversions()) / near.mean_inversions() < 0.02);
        assert(wr::pmf(E, wsl::permutation(5)) == 0 && wr::log_pmf(M, wsl::permutation(7)) == -INFINITY);
        assert(throws_domain_error([] { wr::ewens X({.size = 3, .theta = 0.0}); }) && throws_domain_error([] { wr::mallows X({.size = 3, .q = -1.0}); }));
    }

    // ---- Sous-ensembles, combinaisons, mots -------------------------------------------
    {
        assert(chi_square(wr::uniform_subset({.size = 4}), g, all_subsets(4), draws) < 4);
        assert(chi_square(wr::bernoulli_subset({.size = 4, .p = 0.3}), g, all_subsets(4), draws) < 4);
        assert(chi_square(wr::uniform_combination({.size = 6, .count = 3}), g, all_subsets(6, 3), draws) < 4);
        std::vector<std::vector<std::size_t>> words;
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b)
                for (std::size_t c = 0; c < 3; ++c)
                    words.push_back({a, b, c});
        assert(chi_square(wr::uniform_word({.length = 3, .alphabet = 3}), g, words, draws) < 4);
        wr::uniform_combination const C({.size = 1000000, .count = 50});
        auto const s = wr::native(C, g);
        assert(s.size() == 50 && wr::pmf(C, s) > 0 && close(C.entropy(), std::lgamma(1000001.0) - std::lgamma(51.0) - std::lgamma(999951.0), 1e-12));
        assert(wr::pmf(C, {1, 0}) == 0 && wr::pmf(wr::uniform_subset({.size = 3}), {0, 3}) == 0 && wr::pmf(wr::uniform_word({.length = 2, .alphabet = 2}), {0, 2}) == 0);
        assert(throws_domain_error([] { wr::uniform_combination X({.size = 3, .count = 4}); }));
    }

    // ---- Partitions d'entiers ------------------------------------------------------------
    {
        std::vector<std::vector<std::size_t>> p8;
        std::vector<std::size_t> current;
        partitions(8, 8, current, p8);
        assert(p8.size() == 22);
        assert(chi_square(wr::uniform_integer_partition({.total = 8}), g, p8, draws) < 4);
        assert(wr::uniform_integer_partition({.total = 100}).count() == 190569292);
        assert(close(wr::uniform_integer_partition({.total = 1000}).count(), 24061467864032622473692149727991.0, 1e-15));
        wr::uniform_integer_partition const P({.total = 500});
        auto const q = wr::native(P, g);
        std::size_t total = 0;
        for (std::size_t v : q)
            total += v;
        assert(total == 500 && std::is_sorted(q.rbegin(), q.rend()) && wr::pmf(P, q) > 0 && wr::pmf(P, {3, 4}) == 0);
    }

    // ---- Graphes ------------------------------------------------------------------------
    {
        assert(chi_square(wr::erdos_renyi({.vertices = 4, .p = 0.35}), g, all_graphs(4), draws) < 4);
        assert(chi_square(wr::uniform_graph({.vertices = 5, .edges = 3}), g, all_graphs(5, 3), draws) < 4);
        std::vector<wr::edge_list> regular;                                 // 2-réguliers sur 6 sommets : 70
        for (auto const& e : all_graphs(6, 6))
            if (std::ranges::all_of(degrees(e, 6), [](std::size_t d) { return d == 2; }))
                regular.push_back(e);
        assert(regular.size() == 70);
        assert(chi_square(wr::random_regular_graph({.vertices = 6, .degree = 2}), g, regular, draws, [](wr::edge_list const&) { return 1.0 / 70; }) < 4);
        // Barabási–Albert, n = 4, m = 1 : {01, x2, y3}, P(x) = ½, P(y = x) = ½, sinon ¼.
        std::vector<wr::edge_list> ba;
        std::vector<double> pr;
        for (std::size_t x : {0, 1})
            for (std::size_t y : {0, 1, 2})
            {
                wr::edge_list e{{0, 1}, {x, 2}, {y, 3}};
                std::ranges::sort(e);
                ba.push_back(e);
                pr.push_back(0.5 * (y == x ? 0.5 : 0.25));
            }
        assert((chi_square(wr::barabasi_albert({.vertices = 4, .edges = 1}), g, ba, draws, [&](wr::edge_list const& e) { return pr[static_cast<std::size_t>(std::ranges::find(ba, e) - ba.begin())]; }) < 4));
        wr::barabasi_albert const B({.vertices = 2000, .edges = 3});
        auto const be = wr::native(B, g);
        assert(be.size() == B.edge_count() && std::ranges::adjacent_find(be) == be.end() && wsl::probability::detail::is_edge_list(be, 2000));
        auto const re = wr::native(wr::random_regular_graph({.vertices = 500, .degree = 3}), g);
        assert(std::ranges::all_of(degrees(re, 500), [](std::size_t d) { return d == 3; }) && wsl::probability::detail::is_edge_list(re, 500));
        auto const er = wr::native(wr::erdos_renyi({.vertices = 300, .p = 0.1}), g);
        assert(std::abs(double(er.size()) - 4485) < 6 * std::sqrt(4485 * 0.9) && wr::pmf(wr::erdos_renyi({.vertices = 3, .p = 0.5}), {{0, 1}, {0, 2}}) == 0.125);
        assert(throws_domain_error([] { wr::random_regular_graph X({.vertices = 5, .degree = 3}); }) && throws_domain_error([] { wr::uniform_graph X({.vertices = 4, .edges = 7}); }));
    }

    // ---- Arbres ------------------------------------------------------------------------
    {
        std::vector<wr::edge_list> trees;
        for (auto const& e : all_graphs(4, 3))
            if (wsl::probability::detail::is_tree(e, 4))
                trees.push_back(e);
        assert(trees.size() == 16);                                          // Cayley : 4²
        assert(chi_square(wr::uniform_tree({.vertices = 4}), g, trees, draws) < 4);
        auto const binary = all_binary_trees(4);
        assert(binary.size() == 14);                                         // Catalan C₄
        assert(chi_square(wr::uniform_binary_tree({.vertices = 4}), g, binary, draws) < 4);
        wr::uniform_tree const T({.vertices = 1000});
        auto const t = T(g);
        assert(wsl::probability::detail::is_tree(t, 1000) && wr::pmf(T, t) == 0 && wr::log_pmf(T, t) == -T.entropy());   // n^{−(n−2)} sous-déborde
        wr::uniform_binary_tree const A({.vertices = 1000});
        auto const a = wr::native(A, g);
        assert(a.size() == 1000 && wr::log_pmf(A, a) == -A.entropy() && close(A.entropy(), std::lgamma(2001.0) - 2 * std::lgamma(1001.0) - std::log(1001.0), 1e-13));
        assert(wr::pmf(wr::uniform_tree({.vertices = 2}), {{0, 1}}) == 1 && wr::pmf(wr::uniform_tree({.vertices = 4}), {{0, 1}, {1, 2}, {0, 2}}) == 0);
    }

    // ---- Partitions d'ensemble, compositions --------------------------------------------
    {
        std::vector<std::vector<std::size_t>> rgs, current_rgs;
        std::vector<std::size_t> a;
        restricted_growth(5, 0, a, rgs);
        assert(rgs.size() == 52);                                            // Bell B₅
        assert(chi_square(wr::uniform_set_partition({.size = 5}), g, rgs, draws) < 4);
        assert(wr::uniform_set_partition({.size = 20}).count() == 51724158235372.0);
        std::vector<std::vector<std::size_t>> comps, comps3;
        for (std::uint64_t m = 0; m < 32; ++m)                               // coupures parmi 1 … 5
        {
            std::vector<std::size_t> c;
            std::size_t last = 0;
            for (std::size_t i = 1; i < 6; ++i)
                if ((m >> (i - 1)) & 1)
                {
                    c.push_back(i - last);
                    last = i;
                }
            c.push_back(6 - last);
            comps.push_back(c);
            if (c.size() == 3)
                comps3.push_back(c);
        }
        assert(chi_square(wr::uniform_composition({.total = 6}), g, comps, draws) < 4);
        assert(comps3.size() == 10 && chi_square(wr::uniform_composition({.total = 6, .parts = 3}), g, comps3, draws) < 4);
        auto const sp = wr::native(wr::uniform_set_partition({.size = 150}), g);
        assert(wr::pmf(wr::uniform_set_partition({.size = 150}), sp) > 0 && wr::pmf(wr::uniform_set_partition({.size = 3}), {0, 2, 1}) == 0);
    }

    // ---- Empreintes ----------------------------------------------------------------
    bool ok = true;
    {
        wr::xoshiro256pp x(42);
        auto const add_perm = [](fnv1a& h, wsl::permutation const& p) { for (std::size_t v : p) h.add(v); };
        auto const add_vec = [](fnv1a& h, std::vector<std::size_t> const& v) { h.add(v.size()); for (std::size_t e : v) h.add(e); };
        ok &= check("uniform_permutation", fingerprint([&](fnv1a& h) { add_perm(h, wr::native(wr::uniform_permutation({.size = 20}), x)); }), 0xc24e2f0c1702f325);
        ok &= check("ewens", fingerprint([&](fnv1a& h) { add_perm(h, wr::native(wr::ewens({.size = 20, .theta = 1.7}), x)); }), 0x16e0ba54eb024465);
        ok &= check("mallows", fingerprint([&](fnv1a& h) { add_perm(h, wr::native(wr::mallows({.size = 20, .q = 0.8}), x)); }), 0x1dbf2a90e4629485);
        ok &= check("mallows (q > 1)", fingerprint([&](fnv1a& h) { add_perm(h, wr::native(wr::mallows({.size = 20, .q = 1.05}), x)); }), 0xd7a97235d3517c65);
        ok &= check("subsets", fingerprint([&](fnv1a& h) {
            add_vec(h, wr::native(wr::uniform_subset({.size = 70}), x));
            add_vec(h, wr::native(wr::bernoulli_subset({.size = 30, .p = 0.2}), x));
            add_vec(h, wr::native(wr::uniform_combination({.size = 100, .count = 7}), x));
            add_vec(h, wr::native(wr::uniform_word({.length = 9, .alphabet = 26}), x));
        }), 0x0ed6bf2a62e76d9b);
        auto const add_edges = [](fnv1a& h, wr::edge_list const& e) { h.add(e.size()); for (auto [u, v] : e) h.add((std::uint64_t(u) << 32) | v); };
        ok &= check("graphs", fingerprint([&](fnv1a& h) {
            add_edges(h, wr::native(wr::erdos_renyi({.vertices = 12, .p = 0.3}), x));
            add_edges(h, wr::native(wr::uniform_graph({.vertices = 12, .edges = 15}), x));
            add_edges(h, wr::native(wr::barabasi_albert({.vertices = 20, .edges = 2}), x));
            add_edges(h, wr::native(wr::random_regular_graph({.vertices = 12, .degree = 3}), x));
        }, 1 << 10), 0x553dd3266e70997e);
        ok &= check("trees", fingerprint([&](fnv1a& h) {
            add_edges(h, wr::native(wr::uniform_tree({.vertices = 30}), x));
            for (auto c : wr::native(wr::uniform_binary_tree({.vertices = 20}), x).children)
                h.add((std::uint64_t(c[0]) << 32) ^ c[1]);
        }, 1 << 10), 0xd9ba05da9d874029);
        wr::uniform_set_partition const SP({.size = 25});
        ok &= check("set_partition", fingerprint([&](fnv1a& h) { add_vec(h, wr::native(SP, x)); add_vec(h, wr::native(wr::uniform_composition({.total = 15, .parts = 4}), x)); }, 1 << 10), 0xb24e8c428ebea7a8);
        wr::uniform_integer_partition const P({.total = 60});
        ok &= check("integer_partition", fingerprint([&](fnv1a& h) { add_vec(h, wr::native(P, x)); }, 1 << 10), 0x85ac6c433bff617f);
    }
    assert(ok);
    std::printf("test_random_combinatorics : ok\n");
}
