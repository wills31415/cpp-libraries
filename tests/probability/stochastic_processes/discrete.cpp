// Processus de wsl::random en temps discret : loi catégorielle exacte, vues paresseuses
// (contrat de reproductibilité), marches, chaînes de Markov, urnes de Pólya, restaurant
// chinois ; χ² contre les lois exactes (graines fixes), empreintes.

#include <wsl/probability.hpp>
#include <wsl/random.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <ranges>
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

    // χ² de draws valeurs tirées par draw() contre probability, sur objects ; écart
    // réduit (χ² − ddl)/√(2 ddl), < 4 attendu. Les objets de probabilité nulle doivent
    // ne jamais sortir.
    template <class V, class D, class P>
    double chi_square(std::vector<V> const& objects, int draws, D draw, P probability)
    {
        std::map<V, int> counts;
        for (int i = 0; i < draws; ++i)
            ++counts[draw()];
        double chi = 0, total = 0;
        int classes = 0;
        for (auto const& o : objects)
        {
            double const p = probability(o);
            total += p;
            auto const it = counts.find(o);
            double const c = it == counts.end() ? 0 : it->second;
            if (p == 0)
            {
                assert(c == 0);
                continue;
            }
            double const e = draws * p;
            chi += (c - e) * (c - e) / e;
            ++classes;
        }
        assert(close(total, 1, 1e-12));
        int seen = 0;
        for (auto const& [o, c] : counts)
            seen += c;
        assert(seen == draws);                                               // rien hors de la liste
        double const dof = classes - 1.0;
        return (chi - dof) / std::sqrt(2 * dof);
    }

    // Suites de longueur n sur {0, …, k − 1}.
    std::vector<std::vector<std::size_t>> sequences(std::size_t n, std::size_t k)
    {
        std::vector<std::vector<std::size_t>> r{{}};
        for (std::size_t i = 0; i < n; ++i)
        {
            std::vector<std::vector<std::size_t>> next;
            for (auto const& s : r)
                for (std::size_t c = 0; c < k; ++c)
                {
                    auto t = s;
                    t.push_back(c);
                    next.push_back(t);
                }
            r = next;
        }
        return r;
    }

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
}

static_assert(wr::stochastic_process<wr::markov_chain> && wr::stochastic_process<wr::chinese_restaurant> && wr::stochastic_process<wr::polya_urn>);
static_assert(wr::stochastic_process<wr::random_walk<wr::rademacher>>);
static_assert(wr::probability_law<wr::categorical<>> && wr::categorical<>::sampling == wr::exactness::exact && wr::has_cdf<wr::categorical<>>);
static_assert(std::ranges::input_range<decltype(wr::trajectory(std::declval<wr::markov_chain const&>(), std::declval<wr::xoshiro256pp&>()))>);

int main()
{
    wr::xoshiro256pp g(20261118);
    constexpr int draws = 60000;

    // ---- Loi catégorielle ------------------------------------------------------------
    {
        wr::categorical const C({.p = {0.1, 1e-3, 2.5, 0.0, 7.0}});
        std::vector<std::size_t> const k{0, 1, 2, 3, 4};
        assert(chi_square(k, draws, [&] { return wr::sample(C, g); }, [&](std::size_t i) { return wr::pmf(C, i); }) < 4);
        double const total = 0.1 + 1e-3 + 2.5 + 7.0;
        assert(close(wr::pmf(C, 2), 2.5 / total, 1e-15) && wr::pmf(C, 3) == 0 && wr::pmf(C, 7) == 0 && wr::cdf(C, 4) == 1 && close(wr::cdf(C, 1) + wr::ccdf(C, 1), 1, 1e-15));
        assert(C.mode() == 4 && close(C.mean(), (1e-3 + 5.0 + 28.0) / total, 1e-15));
        wr::categorical const Q({.p = {1.0, 1.0, 2.0}});                    // quantile exact aux frontières
        assert(wr::quantile(Q, 0.25) == 0 && wr::quantile(Q, std::nextafter(0.25, 1.0)) == 1 && wr::quantile(Q, 0.5) == 1 && wr::quantile(Q, 0.5000000001) == 2 && wr::quantile(Q, 0.0) == 0);
        assert(wr::quantile(wr::categorical({.p = {0.0, 0.0, 3.0}}), 0.0) == 2);
        wr::categorical const B({.p = {0x1p600, 0x1p-600, 0x1p599}});      // grands entiers (pas d'alias)
        int above = 0;
        for (int i = 0; i < 30000; ++i)
        {
            auto const x = wr::sample(B, g);
            assert(x != 1);
            above += x == 2;
        }
        assert(std::abs(above / 30000.0 - 1.0 / 3) < 0.012);
        assert(close(C.entropy(), [&] { double h = 0; for (double w : C.weights()) if (w > 0) h -= w / total * std::log(w / total); return h; }(), 1e-14));
        assert(throws_domain_error([] { wr::categorical X({.p = {0.0, 0.0}}); }) && throws_domain_error([] { wr::categorical X({.p = {1.0, -1.0}}); }));
    }

    // ---- Vues paresseuses : contrat de reproductibilité ----------------------------------
    {
        wr::markov_chain const M({.transition = {{0.9, 0.1}, {0.5, 0.5}}, .start = 0});
        for (std::size_t k : {1, 2, 7, 30})
        {
            wr::xoshiro256pp a(5), b(5);
            auto const path = wr::sample(M.steps(k - 1), a);                    // k états
            std::vector<std::size_t> lazy;
            for (std::size_t x : wr::trajectory(M, b) | std::views::take(k))
                lazy.push_back(x);
            assert(path == lazy && a() == b());                             // mêmes valeurs, même consommation
        }
        wr::xoshiro256pp a(9), b(9);
        {
            auto view = wr::trajectory(M, a);                              // rien n'est tiré sans lecture
            (void)view;
        }
        assert(a() == b());
        wr::xoshiro256pp c(11), d(11);                                     // take_while : la valeur qui arrête est tirée
        wr::chinese_restaurant const R({.theta = 2.0});
        std::vector<std::size_t> head;
        for (std::size_t t : wr::trajectory(R, c) | std::views::take_while([](std::size_t t) { return t < 3; }))
            head.push_back(t);
        auto const ref = wr::sample(R.at(head.size() + 1), d);
        assert(std::equal(head.begin(), head.end(), ref.begin()) && ref.back() >= 3 && c() == d());
    }

    // ---- Marche aléatoire ----------------------------------------------------------------
    {
        wr::random_walk const W(wr::rademacher({.p = 0.5}));
        std::vector<std::int64_t> ends;
        for (std::int64_t v = -10; v <= 10; v += 2)
            ends.push_back(v);
        wr::binomial const B({.trials = 10, .p = 0.5});
        assert(chi_square(ends, draws, [&] { return wr::sample(W.steps(10), g).back(); }, [&](std::int64_t v) { return wr::pmf(B, (v + 10) / 2); }) < 4);
        assert(W.mean(10) == 0 && W.variance(10) == 10);
        wr::random_walk const V(wr::rademacher({.p = 0.3}), {.start = 5});
        auto const p = wr::sample(V.steps(6), g);
        assert(p[0] == 5 && close(wr::pmf(V.steps(6), p), std::pow(0.3, double((p.back() - 5 + 6) / 2)) * std::pow(0.7, double(6 - (p.back() - 5 + 6) / 2)), 1e-14));
        assert(wr::pmf(V.steps(2), {5, 7, 8}) == 0 && wr::pmf(V.steps(2), {4, 5, 6}) == 0);
        wr::random_walk const P(wr::poisson({.mean = 2.0}));              // pas quelconques entiers
        assert(close(P.mean(100), 200, 1e-15) && wr::sample(P.steps(100), g).back() > 100);
    }

    // ---- Chaînes de Markov -------------------------------------------------------------
    {
        wr::markov_chain const M({.transition = {{0.5, 0.3, 0.2}, {0.1, 0.8, 0.1}, {0.25, 0.25, 0.5}}, .start = 2});
        auto const at5 = M.at(5);
        std::vector<std::size_t> const states{0, 1, 2};
        assert(chi_square(states, draws, [&] { return wr::sample(M.steps(5), g).back(); }, [&](std::size_t i) { return wr::pmf(at5, i); }) < 4);
        auto const pi = M.stationary(), far = M.at(200);
        for (std::size_t i = 0; i < 3; ++i)
            assert(std::abs(wr::pmf(pi, i) - wr::pmf(far, i)) < 1e-14);
        double s = 0;                                                       // π P = π
        for (std::size_t j = 0; j < 3; ++j)
        {
            double v = 0;
            for (std::size_t i = 0; i < 3; ++i)
                v += wr::pmf(pi, i) * M.probability(i, j);
            s += std::abs(v - wr::pmf(pi, j));
        }
        assert(s < 1e-14);
        auto const paths = sequences(4, 3);                                 // masses des trajectoires : somme 1
        double t = 0;
        for (auto const& q : paths)
            t += wr::pmf(M.steps(3), q);
        assert(close(t, 1, 1e-14) && wr::pmf(M.steps(3), paths[0]) == 0);       // x₀ = 0 impossible
        assert(chi_square(paths, draws, [&] { return wr::sample(M.steps(3), g); }, [&](std::vector<std::size_t> const& q) { return wr::pmf(M.steps(3), q); }) < 4);
        wr::markov_chain const L({.transition = {{1, 0}, {0.4, 0.6}}, .initial = {0.25, 0.75}});
        assert(L.is_absorbing(0) && !L.is_absorbing(1) && close(wr::pmf(L.at(3), 1), 0.75 * 0.216, 1e-15));
        assert(throws_domain_error([] { (void)wr::markov_chain({.transition = {{1, 0}, {0, 1}}, .start = 0}).stationary(); }));   // non unique
        assert(throws_domain_error([] { wr::markov_chain X({.transition = {{0.5, 0.4}, {0.5, 0.5}}, .start = 0}); }));
        assert(throws_domain_error([] { wr::markov_chain X({.transition = {{1, 0}, {0, 1}}, .start = 2}); }));
    }

    // ---- Urnes de Pólya ----------------------------------------------------------------------
    {
        wr::polya_urn const U({.balls = {1, 2}, .added = 1});
        auto const D = U.at(6);                                             // couleurs après 6 tirages
        std::vector<std::vector<std::int64_t>> counts;
        for (std::int64_t a = 0; a <= 6; ++a)
            counts.push_back({a, 6 - a});
        assert(chi_square(counts, draws, [&] {
            auto const c = wr::sample(U.steps(5), g);
            std::int64_t const zeros = std::ranges::count(c, std::size_t{0});
            return std::vector<std::int64_t>{zeros, 6 - zeros};
        }, [&](std::vector<std::int64_t> const& k) { return wr::pmf(D, k); }) < 4);
        auto const seq = sequences(4, 2);
        double t = 0;
        for (auto const& q : seq)
            t += wr::pmf(U.steps(3), q);
        assert(close(t, 1, 1e-14) && close(wr::pmf(U.steps(3), {0, 1, 1, 0}), wr::pmf(U.steps(3), {1, 1, 0, 0}), 1e-15));   // échangeable
        wr::dirichlet_multinomial const M({.trials = 4, .alpha = {0.5, 1.5, 2.0}});
        std::vector<std::vector<std::int64_t>> ks;
        for (std::int64_t a = 0; a <= 4; ++a)
            for (std::int64_t b = 0; a + b <= 4; ++b)
                ks.push_back({a, b, 4 - a - b});
        assert(chi_square(ks, draws, [&] { return wr::sample(M, g); }, [&](std::vector<std::int64_t> const& k) { return wr::pmf(M, k); }) < 4);
        assert(close(M.mean()[2], 2.0, 1e-15) && close(M.variance()[0], 4 * 0.125 * 0.875 * 8 / 5, 1e-14));
        wr::polya_urn const R({.balls = {0.3, 0.7, 1.1}, .added = 0.25});  // poids non entiers : exacts
        assert(chi_square(sequences(3, 3), draws, [&] { return wr::sample(R.steps(2), g); }, [&](std::vector<std::size_t> const& q) { return wr::pmf(R.steps(2), q); }) < 4);
        assert(throws_domain_error([] { (void)wr::polya_urn({.balls = {1, 1}, .added = 0}).at(3); }));
    }

    // ---- Restaurant chinois, partitions d'Ewens–Pitman -----------------------------------------
    {
        std::vector<std::vector<std::size_t>> rgs;
        std::vector<std::size_t> a;
        restricted_growth(5, 0, a, rgs);
        for (auto [theta, d] : {std::pair{1.5, 0.0}, std::pair{1.0, 0.5}, std::pair{-0.25, 0.5}, std::pair{40.0, 0.1}})
        {
            wr::chinese_restaurant const R({.theta = theta, .discount = d});
            auto const P = R.at(5);
            assert(chi_square(rgs, draws, [&] {
                std::vector<std::size_t> s;
                for (std::size_t t : wr::trajectory(R, g) | std::views::take(5))
                    s.push_back(t);
                return s;
            }, [&](std::vector<std::size_t> const& q) { return wr::pmf(P, q); }) < 4);
        }
        wr::ewens_partition const E({.size = 200, .theta = 2.0});
        double k = 0;
        for (int i = 0; i < 500; ++i)
        {
            auto const s = wr::sample(E, g);
            k += double(*std::ranges::max_element(s) + 1);
        }
        assert(std::abs(k / 500 - E.mean_blocks()) < 0.3 && close(E.mean_blocks(), wr::ewens({.size = 200, .theta = 2.0}).mean_cycles(), 1e-15));
        wr::ewens_partition const PY({.size = 200, .theta = 1.0, .discount = 0.5});
        double kp = 0;
        for (int i = 0; i < 500; ++i)
        {
            auto const s = wr::sample(PY, g);
            kp += double(*std::ranges::max_element(s) + 1);
        }
        assert(std::abs(kp / 500 - PY.mean_blocks()) < 0.8);
        assert(wr::pmf(E, std::vector<std::size_t>(200, 0)) > 0 && wr::pmf(wr::ewens_partition({.size = 3, .theta = 1.0}), {0, 2, 1}) == 0);
        assert(throws_domain_error([] { wr::chinese_restaurant X({.theta = -0.6, .discount = 0.5}); }));
    }

    // ---- Empreintes ----------------------------------------------------------------
    bool ok = true;
    {
        wr::xoshiro256pp x(42);
        wr::categorical const C({.p = {0.1, 0.25, 0.3, 0.05, 0.3}});
        wr::markov_chain const M({.transition = {{0.5, 0.3, 0.2}, {0.1, 0.8, 0.1}, {0.25, 0.25, 0.5}}, .start = 0});
        auto const take = [&](fnv1a& h, auto const& process, std::size_t n) { for (auto v : wr::trajectory(process, x) | std::views::take(n)) h.add(static_cast<std::uint64_t>(v)); };
        ok &= check("categorical", fingerprint([&](fnv1a& h) { h.add(wr::sample(C, x)); }), 0x75204fe3d6fc7447);
        ok &= check("random_walk", fingerprint([&](fnv1a& h) { take(h, wr::random_walk(wr::rademacher({.p = 0.6})), 20); }, 1 << 9), 0xe141a7ef8f3eb38a);
        ok &= check("markov_chain", fingerprint([&](fnv1a& h) { take(h, M, 20); }, 1 << 9), 0xbd7b75b6f3dea6a4);
        ok &= check("polya_urn", fingerprint([&](fnv1a& h) { take(h, wr::polya_urn({.balls = {1, 2, 0.5}, .added = 1.5}), 20); }, 1 << 9), 0x0c3d0a6d77d909e4);
        ok &= check("chinese_restaurant", fingerprint([&](fnv1a& h) {
            take(h, wr::chinese_restaurant({.theta = 1.5}), 20);
            take(h, wr::chinese_restaurant({.theta = 0.5, .discount = 0.3}), 20);
        }, 1 << 9), 0xb288345378897f45);
    }
    assert(ok);
    std::printf("test_random_processes : ok\n");
}
