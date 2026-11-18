// Processus de wsl::random en temps continu (option B : exacts point par point) : sommes
// paresseuses exactes (contrôle par raffinement), brownien, pont, Ornstein–Uhlenbeck,
// brownien géométrique, Poisson (homogène, inhomogène, composé), chaîne de Markov en temps
// continu, marches à pas réels ; Kolmogorov–Smirnov et χ² contre les lois exactes,
// contrat de préfixe, empreintes.

#include <wsl/random.hpp>

#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace
{
    namespace d = wsl::detail;

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
    std::uint64_t fingerprint(F f, int count = 1 << 9)
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

    std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

    template <class F>
    bool throws_domain_error(F f)
    {
        try { f(); } catch (std::domain_error const&) { return true; }
        return false;
    }

    // √n · écart maximal de la répartition empirique de v à F, < 1,95 au seuil de 0,1 %.
    template <class F>
    double kolmogorov_smirnov(std::vector<double> v, F F_)
    {
        std::ranges::sort(v);
        double worst = 0, n = double(v.size());
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            double const f = F_(v[i]);
            worst = std::max({worst, std::abs(f - double(i) / n), std::abs(f - double(i + 1) / n)});
        }
        return worst * std::sqrt(n);
    }

    template <class L>
    double ks_against(std::vector<double> const& v, L const& law) { return kolmogorov_smirnov(v, [&](double x) { return double(cdf(law, x)); }); }
}

int main()
{
    wsl::xoshiro256pp g(20261118);
    constexpr int n = 20000;

    // ---- Sommes paresseuses : arrondi correct, contrôlé par raffinement ---------------------
    {
        int bad = 0;
        for (int trial = 0; trial < 1500; ++trial)
        {
            d::lazy_sum s;
            int const m = 1 + int(wsl::uniform_below(g, 30));
            for (int i = 0; i < m; ++i)
                switch (wsl::uniform_below(g, 4))
                {
                case 0: s.add(d::lazy_sample(wsl::normal({.mean = 0.3, .sd = 1e-3 * (i + 1)}), g), g); break;
                case 1: s.add(d::lazy_sample(wsl::exponential({.rate = 2.5}), g), g); break;
                case 2: s.add(d::lazy_sample(wsl::laplace({.location = -1.0, .scale = 0.7}), g), g); break;
                default: s.add(d::lazy_sample(wsl::uniform({.low = -2.0, .high = 3.0}), g), g); break;
                }
            d::dyadic_value const offset = d::to_dyadic(trial % 3 == 0 ? 0.0 : -0.3 * m);   // compensations
            d::dyadic_value const alpha = d::to_dyadic(trial % 2 ? 0.5 : 0.0), beta = d::to_dyadic(trial % 5 == 0 ? -1.0 / 3 : 1.0);
            std::size_t const k = std::size_t(m) / 2;
            double const r = s.rounded<double>(g, offset, alpha, k, beta);
            float const rf = s.rounded<float>(g, offset, alpha, k, beta);
            bad += !s.confirms(g, r, 2, offset, alpha, k, beta) || !s.confirms(g, rf, 0, offset, alpha, k, beta);
        }
        assert(bad == 0);
        int tight = 0;                                                       // valeur visée à un milieu de deux doubles
        for (int trial = 0; trial < 300; ++trial)
        {
            d::lazy_sum s;
            for (int i = 0; i < 5; ++i)
                s.add(d::lazy_sample(wsl::normal({.mean = 0.0, .sd = 1.0}), g), g);
            double const v = s.rounded(g);
            d::dyadic_value const half = d::to_dyadic((std::nextafter(v, INFINITY) - v) / 2);
            double const r = s.rounded<double>(g, half, {}, 0, d::to_dyadic(1.0));
            tight += !s.confirms(g, r, 2, half, {}, 0, d::to_dyadic(1.0));
        }
        assert(tight == 0);
    }

    // ---- Mouvement brownien --------------------------------------------------------------
    {
        wsl::brownian_motion const B({.drift = 0.5, .volatility = 2.0, .start = 1.0});
        // Contrat de préfixe : trajectoire à pas 1/4 = grille {0, 1/4, …}.
        for (std::size_t k : {1, 5, 40})
        {
            wsl::xoshiro256pp a(7), b(7);
            std::vector<double> times;
            for (std::size_t i = 0; i < k; ++i)
                times.push_back(0.25 * double(i));
            auto const grid = sample(B.on_grid(times), a);
            std::vector<double> lazy;
            for (auto [t, v] : wsl::trajectory(B, b, {.step = 0.25}) | std::views::take(k))
                lazy.push_back(v);
            assert(grid == lazy && a() == b());
        }
        std::vector<double> at1, increments;
        for (int i = 0; i < n; ++i)
        {
            auto const x = sample(B.on_grid({1.0, 3.0}), g);
            at1.push_back(x[0]);
            increments.push_back(x[1] - x[0]);
        }
        assert(ks_against(at1, B.at(1.0)) < 1.95);
        assert(ks_against(increments, wsl::normal({.mean = 1.0, .sd = 2.0 * std::sqrt(2.0)})) < 1.95);   // accroissement indépendant
        assert(sample(B.on_grid({0.0}), g)[0] == 1.0);
        std::vector<double> fast;
        for (int i = 0; i < n; ++i)
        {
            double last = 0;
            for (auto [t, v] : wsl::trajectory(B, g, {.step = 0.5}, wsl::fast) | std::views::take(3))
                last = v;
            fast.push_back(last);
        }
        assert(ks_against(fast, B.at(1.0)) < 1.95);
        assert(throws_domain_error([] { wsl::brownian_motion X({.volatility = 0.0}); }) && throws_domain_error([&] { (void)B.on_grid({1.0, 0.5}); }));
    }

    // ---- Pont brownien, Ornstein–Uhlenbeck, brownien géométrique -------------------------------
    {
        wsl::brownian_bridge const D({.from = {0, 1}, .to = {2, -1}, .volatility = 1.5});
        std::vector<double> mid;
        for (int i = 0; i < n; ++i)
        {
            auto const x = sample(D.on_grid({0.0, 0.5, 2.0}), g);
            assert(x[0] == 1 && x[2] == -1);                                 // extrémités exactes
            mid.push_back(x[1]);
        }
        assert(ks_against(mid, D.at(0.5)) < 1.95);
        wsl::ornstein_uhlenbeck const O({.theta = 1.5, .mean = 1.0, .volatility = 0.8, .start = -2.0});
        std::vector<double> o1, o5;
        for (int i = 0; i < n; ++i)
        {
            auto const x = sample(O.on_grid({0.4, 5.0}), g);
            o1.push_back(x[0]);
            o5.push_back(x[1]);
        }
        assert(ks_against(o1, O.at(0.4)) < 1.95 && ks_against(o5, O.at(5.0)) < 1.95 && ks_against(o5, O.stationary()) < 1.95);
        std::vector<double> lazy;
        wsl::xoshiro256pp a(3), b(3);
        for (auto [t, v] : wsl::trajectory(O, a, {.step = 0.5}) | std::views::take(4))
            lazy.push_back(v);
        assert(lazy == sample(O.on_grid({0.0, 0.5, 1.0, 1.5}), b));
        assert(throws_domain_error([&] { (void)sample(O.on_grid({300.0}), g); }));   // θ t > 350
        wsl::geometric_brownian_motion const S({.drift = 0.05, .volatility = 0.2, .start = 100});
        std::vector<double> s1;
        for (int i = 0; i < n; ++i)
            s1.push_back(sample(S.on_grid({1.0}), g)[0]);
        assert(ks_against(s1, S.at(1.0)) < 1.95);
    }

    // ---- Processus de Poisson -----------------------------------------------------------------
    {
        wsl::poisson_process const N({.rate = 2.0});
        std::vector<double> gaps;
        std::vector<int> counts(40, 0);
        for (int i = 0; i < n; ++i)
        {
            auto const t = sample(N.on(3.0), g);
            assert(std::ranges::is_sorted(t) && (t.empty() || t.back() <= 3.0));
            ++counts[std::min<std::size_t>(t.size(), 39)];
            if (t.size() >= 2)
                gaps.push_back(t[1] - t[0]);
        }
        double chi = 0;
        int classes = 0;
        for (int k = 0; k < 39; ++k)                                         // χ² de N(3) contre poisson(6)
        {
            double const e = n * pmf(N.at(3.0), k);
            if (e < 20)
                continue;
            chi += (counts[k] - e) * (counts[k] - e) / e;
            ++classes;
        }
        assert((chi - (classes - 1)) / std::sqrt(2.0 * (classes - 1)) < 4);
        std::vector<double> first;
        for (int i = 0; i < n; ++i)
            for (double t : wsl::trajectory(N, g) | std::views::take(1))
                first.push_back(t);
        assert(ks_against(first, wsl::exponential({.rate = 2.0})) < 1.95);
        wsl::inhomogeneous_poisson_process const M([](double t) { return 1 + std::sin(t); }, 2.0);
        double total = 0;
        for (int i = 0; i < 4000; ++i)
            total += double(sample(M.on(5.0), g).size());
        double const expected = 5 + (1 - std::cos(5.0));                    // ∫₀⁵ (1 + sin t) dt
        assert(std::abs(total / 4000 - expected) < 0.15);
        wsl::compound_poisson const C(2.0, wsl::exponential({.rate = 0.5}));
        double s = 0;
        for (int i = 0; i < 4000; ++i)
        {
            double value = 0;
            for (auto [t, v] : wsl::trajectory(C, g) | std::views::take_while([](auto p) { return p.first <= 3.0; }))
                value = v;
            s += value;
        }
        assert(std::abs(s / 4000 - C.mean(3.0)) < 0.35 && C.mean(3.0) == 12 && C.variance(3.0) == 48);
        wsl::compound_poisson const A(1.0, wsl::gamma({.shape = 2.0, .rate = 1.0}));   // sauts non paresseux
        static_assert(decltype(A)::sampling == wsl::exactness::approximate && wsl::compound_poisson<wsl::normal<>>::sampling == wsl::exactness::exact);
    }

    // ---- Chaîne de Markov en temps continu ------------------------------------------------------
    {
        wsl::continuous_markov_chain const Q({.generator = {{-1, 1, 0}, {2, -3, 1}, {0, 0.5, -0.5}}, .start = 0});
        std::vector<int> counts(3, 0);
        for (int i = 0; i < n; ++i)
        {
            std::size_t state = 0;
            for (auto [t, x] : wsl::trajectory(Q, g) | std::views::take_while([](auto p) { return p.first <= 1.5; }))
                state = x;
            ++counts[state];
        }
        auto const law = Q.at(1.5);
        double chi = 0;
        for (std::size_t k = 0; k < 3; ++k)
        {
            double const e = n * pmf(law, k);
            chi += (counts[k] - e) * (counts[k] - e) / e;
        }
        assert(chi < 13.8);                                                  // 2 ddl, 0,1 %
        wsl::continuous_markov_chain const T({.generator = {{-1, 1}, {0, 0}}, .start = 0});   // absorbant
        std::vector<std::pair<double, std::size_t>> path;
        for (auto p : wsl::trajectory(T, g) | std::views::take(4))
            path.push_back(p);
        assert(path[1].second == 1 && path[2].first == INFINITY && path[3].second == 1);
        assert(throws_domain_error([] { wsl::continuous_markov_chain X({.generator = {{-1, 2}, {1, -1}}, .start = 0}); }));
    }

    // ---- Marches à pas réels ------------------------------------------------------------------
    {
        wsl::random_walk const W(wsl::normal({.mean = 0.0, .sd = 1.0}));
        std::vector<double> x10, e5;
        for (int i = 0; i < n; ++i)
            x10.push_back(sample(W.steps(10), g).back());
        assert(ks_against(x10, wsl::normal({.mean = 0.0, .sd = std::sqrt(10.0)})) < 1.95);
        wsl::random_walk const E(wsl::exponential({.rate = 2.0}));
        for (int i = 0; i < n; ++i)
            e5.push_back(sample(E.steps(5), g).back());
        assert(ks_against(e5, wsl::gamma({.shape = 5.0, .rate = 2.0})) < 1.95);   // somme exacte d'exponentielles
        wsl::random_walk const G(wsl::gamma({.shape = 2.0, .rate = 1.0}), {.start = 1.0});   // pas non paresseux
        double s = 0;
        for (int i = 0; i < 4000; ++i)
            s += sample(G.steps(4), g).back();
        assert(std::abs(s / 4000 - 9) < 0.15);
    }

    // ---- Empreintes ----------------------------------------------------------------
    bool ok = true;
    {
        wsl::xoshiro256pp x(42);
        auto const take = [&](fnv1a& h, auto&& view) { for (auto v : view) { if constexpr (requires { v.first; }) { h.add(bits(double(v.first))); h.add(bits(double(v.second))); } else h.add(bits(double(v))); } };
        wsl::brownian_motion const B({.drift = 0.1, .volatility = 1.3});
        wsl::ornstein_uhlenbeck const O({.theta = 0.7, .mean = 1.0, .volatility = 0.5, .start = 0.0});
        ok &= check("brownian", fingerprint([&](fnv1a& h) { take(h, wsl::trajectory(B, x, {.step = 0.1}) | std::views::take(20)); }), 0x74f897663e90c36e);
        ok &= check("bridge", fingerprint([&](fnv1a& h) { for (double v : sample(wsl::brownian_bridge({.from = {0, 0}, .to = {1, 1}}).on_grid({0.1, 0.4, 0.9}), x)) h.add(bits(v)); }), 0x637ed1f3e7bb1e29);
        ok &= check("ornstein_uhlenbeck", fingerprint([&](fnv1a& h) { take(h, wsl::trajectory(O, x, {.step = 0.25}) | std::views::take(20)); }), 0xc1712e764addec70);
        ok &= check("poisson", fingerprint([&](fnv1a& h) { take(h, wsl::trajectory(wsl::poisson_process({.rate = 3.0}), x) | std::views::take(20)); }), 0x05e3b438010ccb3a);
        ok &= check("compound_poisson", fingerprint([&](fnv1a& h) { take(h, wsl::trajectory(wsl::compound_poisson(1.5, wsl::laplace({.location = 0.0, .scale = 1.0})), x) | std::views::take(10)); }), 0x20e17e9986f507f7);
        wsl::continuous_markov_chain const Q({.generator = {{-1, 1, 0}, {2, -3, 1}, {0, 0.5, -0.5}}, .start = 0});
        ok &= check("continuous_markov", fingerprint([&](fnv1a& h) { take(h, wsl::trajectory(Q, x) | std::views::take(10)); }), 0xa46709f22568617a);
        ok &= check("real_walk", fingerprint([&](fnv1a& h) { take(h, wsl::trajectory(wsl::random_walk(wsl::uniform({.low = -1.0, .high = 2.0})), x) | std::views::take(20)); }), 0x96f2cfeff63723cb);
    }
    assert(ok);
    std::printf("test_random_continuous : ok\n");
}
