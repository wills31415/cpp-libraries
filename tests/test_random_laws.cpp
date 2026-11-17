// Lois de wsl::random : concept et attributs, paramètres, valeurs exactes,
// cas extrêmes, lois vérifiées statistiquement (graines fixes), adaptateur vers
// <random>, empreintes des tirages.

#include <wsl/random.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
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
    std::uint64_t fingerprint(F f, int count = 1 << 16)
    {
        fnv1a h;
        for (int i = 0; i < count; ++i)
            h.add(f());
        return h.value;
    }

    bool check(char const* name, std::uint64_t actual, std::uint64_t expected)
    {
        bool const ok = actual == expected;
        std::printf("%-22s %016llx %s\n", name, static_cast<unsigned long long>(actual),
                    ok ? "ok" : "DIFFÉRENTE de l'empreinte attendue");
        std::fflush(stdout);
        return ok;
    }

    template <class F>
    bool throws_domain_error(F f)
    {
        try { f(); } catch (std::domain_error const&) { return true; }
        return false;
    }

    // Kolmogorov–Smirnov : écart maximal entre la répartition empirique et cdf.
    template <class L>
    double kolmogorov_smirnov(L const& law, std::vector<double> sample)
    {
        std::ranges::sort(sample);
        double worst = 0;
        double const n = double(sample.size());
        for (std::size_t i = 0; i < sample.size(); ++i)
        {
            double const f = cdf(law, sample[i]);
            worst = std::max({worst, std::abs(f - double(i) / n), std::abs(f - double(i + 1) / n)});
        }
        return worst * std::sqrt(n);                                         // < 1,95 au seuil de 0,1 %
    }
}

// ---- Concept, attributs présents ou absents, déduction du type ------------------
static_assert(wsl::probability_law<wsl::uniform<>> && wsl::probability_law<wsl::uniform<float>>);
static_assert(wsl::probability_law<wsl::discrete_uniform<>> && wsl::probability_law<wsl::bernoulli>);
static_assert(wsl::has_mean<wsl::uniform<>> && wsl::has_entropy<wsl::uniform<>> && wsl::has_quantile<wsl::uniform<>>);
static_assert(wsl::has_pdf<wsl::uniform<>> && !wsl::has_pmf<wsl::uniform<>>);
static_assert(wsl::has_pmf<wsl::discrete_uniform<>> && !wsl::has_pdf<wsl::discrete_uniform<>>);
static_assert(!wsl::has_quantile<wsl::bernoulli> && !wsl::has_mode<wsl::uniform<>>);
static_assert(wsl::samplable<wsl::uniform<>, wsl::xoshiro256pp> && wsl::samplable<wsl::bernoulli, std::mt19937>);
static_assert(wsl::uniform<>::sampling == wsl::exactness::approximate && wsl::discrete_uniform<>::sampling == wsl::exactness::exact);

static_assert([] { wsl::uniform X({.low = 2, .high = 6}); return std::is_same_v<decltype(X), wsl::uniform<double>>; }());
static_assert(wsl::uniform({.low = 2.0, .high = 6.0}).mean() == 4.0);
static_assert(cdf(wsl::uniform({.low = 2.0, .high = 6.0}), 3.0) == 0.25);
static_assert(quantile(wsl::discrete_uniform({.low = 1, .high = 6}), 0.5) == 3);
static_assert([] { wsl::pcg64 g(1); auto const x = sample(wsl::discrete_uniform({.low = 1, .high = 6}), g); return x >= 1 && x <= 6; }());

int main()
{
    wsl::xoshiro256pp g(2026);

    // ---- uniform --------------------------------------------------------------
    {
        wsl::uniform const X({.low = 2, .high = 6});
        assert(X.mean() == 4 && X.median() == 4 && X.variance() == 16.0 / 12 && X.skewness() == 0);
        assert(X.entropy() == wsl::special::log(4.0));
        assert(pdf(X, 3.0) == 0.25 && pdf(X, 6.0) == 0 && pdf(X, 1.0) == 0 && log_pdf(X, 7.0) == -INFINITY);
        assert(cdf(X, 3.0) == 0.25 && cdf(X, 1.0) == 0 && cdf(X, 9.0) == 1 && ccdf(X, 3.0) == 0.75);
        assert(quantile(X, 0.25) == 3 && quantile(X, 0.0) == 2 && quantile(X, 1.0) == 6);
        assert(cdf(X, 0.1 + 2) + ccdf(X, 0.1 + 2) == 1);
        assert(throws_domain_error([&] { (void)quantile(X, 1.5); }));
        assert(throws_domain_error([] { wsl::uniform Y({.low = 1, .high = 1}); }));
        assert(throws_domain_error([] { wsl::uniform Y({.low = 0, .high = INFINITY}); }));

        // Intervalle immense : b − a déborde, rien d'infini ni de NaN là où le résultat est fini.
        double const big = std::numeric_limits<double>::max();
        wsl::uniform const H({.low = -big, .high = big});
        assert(H.mean() == 0 && std::isinf(H.variance()) && std::abs(H.entropy() - 710.4758600739439) < 1e-12);
        assert(pdf(H, 0.0) > 0 && std::isfinite(pdf(H, 0.0)) && cdf(H, 0.0) == 0.5 && quantile(H, 0.5) == 0);
        assert(quantile(H, 0.75) == big / 2 && ccdf(H, big / 2) == 0.25);

        // Valeurs tirées dans [low, high), loi uniforme (Kolmogorov–Smirnov).
        std::vector<double> values;
        for (int i = 0; i < 20000; ++i)
        {
            double const x = sample(X, g);
            assert(x >= 2 && x < 6);
            values.push_back(x);
        }
        assert(kolmogorov_smirnov(X, values) < 1.95);
        wsl::uniform const tiny({.low = 1.0, .high = std::nextafter(1.0, 2.0)});   // un seul flottant possible
        for (int i = 0; i < 100; ++i)
            assert(sample(tiny, g) == 1.0);
        wsl::uniform<float> const F({.low = -1, .high = 1});
        for (int i = 0; i < 10000; ++i)
        {
            float const x = sample(F, g);
            assert(x >= -1 && x < 1);
        }
        assert(F.mean() == 0 && cdf(F, 0.5f) == 0.75f);
    }

    // ---- discrete_uniform -----------------------------------------------------
    {
        wsl::discrete_uniform const D({.low = 1, .high = 6});
        assert(D.mean() == 3.5 && D.median() == 3 && D.variance() == 35.0 / 12);
        assert(pmf(D, 3) == 1.0 / 6 && pmf(D, 7) == 0 && cdf(D, 3) == 0.5 && ccdf(D, 3) == 0.5 && cdf(D, 6) == 1);
        for (int k = 1; k <= 6; ++k)
        {
            // Quantile exact ; la cdf arrondie le décale au plus d'un cran.
            long long const q = quantile(D, cdf(D, k));
            assert(q == k || q == k + 1);
        }
        assert(quantile(D, 1e-300) == 1 && quantile(D, 1.0) == 6 && quantile(D, 0.5) == 3 && quantile(D, 0.25) == 2);
        assert(quantile(D, 0.5 + 0x1p-52) == 4);                            // juste au-dessus de 3/6
        std::array<int, 7> counts{};
        for (int i = 0; i < 60000; ++i)
            ++counts[static_cast<std::size_t>(sample(D, g))];
        for (int k = 1; k <= 6; ++k)
            assert(std::abs(counts[static_cast<std::size_t>(k)] - 10000) < 450);

        // Pleine plage de 64 bits : n = 2^64.
        wsl::discrete_uniform const W({.low = std::numeric_limits<std::int64_t>::min(), .high = std::numeric_limits<std::int64_t>::max()});
        assert(pmf(W, 0) == 0x1p-64 && cdf(W, -1) == 0.5 && W.variance() == 0x1p128 / 12 && W.entropy() == 64 * std::log(2.0));
        assert(quantile(W, 0.5) == -1 && quantile(W, 0x1p-64) == std::numeric_limits<std::int64_t>::min());
        wsl::discrete_uniform<std::uint8_t> const B({.low = 0, .high = 255});
        assert(quantile(B, 0.5) == 127 && pmf(B, 255) == 1.0 / 256);
        assert(throws_domain_error([] { wsl::discrete_uniform E({.low = 2, .high = 1}); }));
    }

    // ---- bernoulli ------------------------------------------------------------
    {
        wsl::bernoulli const B({.p = 0.3});
        assert(B.mean() == 0.3 && std::abs(B.variance() - 0.21) < 1e-16 && pmf(B, true) == 0.3 && cdf(B, false) == 0.7);
        assert(std::abs(B.entropy() - 0.6108643020548935) < 1e-15);
        int hits = 0;
        for (int i = 0; i < 100000; ++i)
            hits += sample(B, g);
        assert(std::abs(hits / 100000.0 - 0.3) < 0.005);
        assert((wsl::bernoulli({.p = 0}).support() == wsl::interval<bool>{false, false}));
        assert(wsl::bernoulli({.p = 1e-300}).variance() == 1e-300 && wsl::bernoulli({.p = 0}).entropy() == 0);
        assert(throws_domain_error([] { wsl::bernoulli E({.p = 1.5}); }));
    }

    // ---- Tirage rapide demandé : à défaut de méthode propre, le tirage par défaut ----
    {
        wsl::xoshiro256pp a(5), b(5);
        wsl::uniform const X({.low = 0, .high = 1});
        for (int i = 0; i < 100; ++i)
            assert(sample(X, a, wsl::fast) == sample(X, b));
    }

    // ---- Adaptateur vers <random> -------------------------------------------------
    {
        auto d = wsl::as_distribution(wsl::discrete_uniform({.low = 1, .high = 6}));
        std::vector<long long> rolls(1000);
        std::ranges::generate(rolls, [&] { return d(g); });
        assert(std::ranges::all_of(rolls, [](long long x) { return x >= 1 && x <= 6; }));
        assert(d.min() == 1 && d.max() == 6 && d == d);
        decltype(d)::param_type const other{wsl::discrete_uniform({.low = 10, .high = 12})};
        assert(d(g, other) >= 10);
        auto u = wsl::as_distribution(wsl::uniform({.low = 0.0, .high = 1.0}));
        assert(u.min() == 0.0 && u.max() == 1.0 && u(g) < 1.0);
    }

    // ---- Empreintes ---------------------------------------------------------------
    bool ok = true;
    {
        wsl::xoshiro256pp x(11);
        wsl::uniform const U({.low = -3.5, .high = 1e10});
        wsl::uniform<float> const F({.low = 0, .high = 7});
        wsl::discrete_uniform const D({.low = -1000, .high = 1000000});
        wsl::bernoulli const B({.p = 0.125});
        ok &= check("uniform", fingerprint([&] { return std::bit_cast<std::uint64_t>(sample(U, x)); }), 0xe9028c88f6ddde7f);
        ok &= check("uniform float", fingerprint([&] { return std::uint64_t{std::bit_cast<std::uint32_t>(sample(F, x))}; }), 0x53a47c22e845da0b);
        ok &= check("discrete_uniform", fingerprint([&] { return static_cast<std::uint64_t>(sample(D, x)); }), 0x3f058c88c4d229e8);
        ok &= check("bernoulli", fingerprint([&] { return std::uint64_t{sample(B, x)}; }), 0x88e84634bf3325c5);
        ok &= check("uniform quantile", fingerprint([&] {
            double const p = wsl::uniform_real(x);
            return std::bit_cast<std::uint64_t>(quantile(U, p)) ^ std::bit_cast<std::uint64_t>(cdf(U, quantile(U, p)));
        }), 0x40887993489ea177);
    }
    assert(ok);
    return 0;
}
