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

static_assert(wsl::probability_law<wsl::exponential<>> && wsl::probability_law<wsl::normal<float>>);
static_assert(wsl::has_pdf<wsl::exponential<>> && wsl::has_mode<wsl::exponential<>> && wsl::has_quantile<wsl::normal<>>);
static_assert(wsl::exponential<>::sampling == wsl::exactness::exact && wsl::normal<>::sampling == wsl::exactness::exact);
static_assert([] { wsl::normal X({.mean = 0, .sd = 1}); return std::is_same_v<decltype(X), wsl::normal<double>>; }());
static_assert(cdf(wsl::exponential({.rate = 1.0}), 0.0) == 0 && cdf(wsl::normal({.mean = 2.0, .sd = 3.0}), 2.0) == 0.5);
static_assert(quantile(wsl::normal({.mean = 2.0, .sd = 3.0}), 0.5) == 2.0);
// Tirages exacts en évaluation constante (réels paresseux, grands entiers, voie lente comprise).
static_assert([] { wsl::pcg64 g(3); return sample(wsl::exponential({.rate = 2.0}), g) > 0; }());
static_assert([] { wsl::pcg64 g(4); double const x = sample(wsl::normal({.mean = 0.0, .sd = 0x1p-900}), g); return x == x; }());

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

    // ---- exponential ----------------------------------------------------------
    {
        wsl::exponential const E({.rate = 2});
        assert(E.mean() == 0.5 && E.variance() == 0.25 && E.mode() == 0 && E.skewness() == 2 && E.excess_kurtosis() == 6);
        assert(E.median() == quantile(E, 0.5) && std::abs(E.median() - std::log(2.0) / 2) < 1e-16);
        assert(std::abs(E.entropy() - (1 - std::log(2.0))) < 1e-16);
        assert(pdf(E, 0.0) == 2 && pdf(E, -1.0) == 0 && pdf(E, INFINITY) == 0 && log_pdf(E, -1.0) == -INFINITY);
        assert(cdf(E, 0.0) == 0 && ccdf(E, 0.0) == 1 && cdf(E, INFINITY) == 1 && ccdf(E, INFINITY) == 0);
        assert(quantile(E, 0.0) == 0 && quantile(E, 1.0) == INFINITY && std::isnan(cdf(E, NAN)));
        assert(cdf(E, 1e-300) == 2e-300 && ccdf(E, 360.0) > 0 && ccdf(E, 360.0) < 1e-300);
        assert(throws_domain_error([] { wsl::exponential X({.rate = 0}); }));
        assert(throws_domain_error([] { wsl::exponential X({.rate = INFINITY}); }));
        assert(throws_domain_error([&] { (void)quantile(E, -0.1); }));
        // Paramètres extrêmes : rien de NaN, densité finie.
        wsl::exponential const big({.rate = 0x1p1000}), small({.rate = 0x1p-1000});
        assert(pdf(big, 0x1p-1000) > 0 && std::isfinite(pdf(big, 0x1p-1000)) && cdf(small, 0x1p1000) > 0.63);

        std::vector<double> exact, fast;
        for (int i = 0; i < 20000; ++i)
        {
            exact.push_back(sample(E, g));
            fast.push_back(sample(E, g, wsl::fast));
        }
        assert(std::ranges::all_of(exact, [](double x) { return x >= 0; }));
        assert(kolmogorov_smirnov(E, exact) < 1.95 && kolmogorov_smirnov(E, fast) < 1.95);
        wsl::exponential<float> const F({.rate = 0.5f});
        std::vector<double> floats;
        for (int i = 0; i < 20000; ++i)
            floats.push_back(sample(F, g));
        assert(kolmogorov_smirnov(wsl::exponential({.rate = 0.5}), floats) < 1.95);
    }

    // ---- normal ---------------------------------------------------------------
    {
        wsl::normal const N({.mean = 1, .sd = 2});
        assert(N.mean() == 1 && N.median() == 1 && N.mode() == 1 && N.variance() == 4 && N.skewness() == 0 && N.excess_kurtosis() == 0);
        assert(std::abs(N.entropy() - (std::log(2.0) + 0.5 * (1 + std::log(2 * 3.141592653589793)))) < 1e-15);
        assert(std::abs(pdf(N, 1.0) - 0.19947114020071635) < 1e-17 && cdf(N, 1.0) == 0.5 && quantile(N, 0.5) == 1);
        for (double d : {0.375, 1.0, 3.0, 10.0, 40.0})                 // 1 ± d exacts
            assert(cdf(N, 1 - d) == ccdf(N, 1 + d));                        // symétrie, exacte
        assert(std::abs(cdf(N, 3.0) - 0.8413447460685429) < 1e-16 && std::abs(quantile(N, 0.8413447460685429) - 3) < 1e-14);
        assert(quantile(N, 0.0) == -INFINITY && quantile(N, 1.0) == INFINITY && cdf(N, -INFINITY) == 0 && pdf(N, INFINITY) == 0);
        assert(std::isnan(pdf(N, NAN)) && std::isnan(cdf(N, NAN)) && log_pdf(N, INFINITY) == -INFINITY);
        assert(cdf(N, -74.0) > 0 && cdf(N, -74.0) < 1e-307 && ccdf(N, 1e6) == 0);
        assert(throws_domain_error([] { wsl::normal X({.mean = 0, .sd = 0}); }));
        assert(throws_domain_error([] { wsl::normal X({.mean = INFINITY, .sd = 1}); }));
        // Paramètres extrêmes.
        double const big = std::numeric_limits<double>::max();
        wsl::normal const H({.mean = 0, .sd = big}), T({.mean = 0, .sd = 1e-300});
        assert(std::abs(cdf(H, big) - 0.8413447460685429) < 1e-16 && quantile(H, 0.975) == INFINITY && std::isfinite(quantile(H, 0.75)));
        assert(std::isfinite(pdf(T, 0.0)) && pdf(T, 0.0) > 1e299 && cdf(T, 1e-300) == cdf(N, 3.0));
        assert(quantile(wsl::normal({.mean = big, .sd = big}), 0.01) < big);   // μ + σ z sur les moitiés

        std::vector<double> exact, fast;
        int tail = 0;
        for (int i = 0; i < 40000; ++i)
        {
            exact.push_back(sample(N, g));
            fast.push_back(sample(N, g, wsl::fast));
            tail += std::abs(fast.back() - 1) > 2 * 3.6541528853610088;     // queue de la ziggourat
        }
        assert(kolmogorov_smirnov(N, exact) < 1.95 && kolmogorov_smirnov(N, fast) < 1.95);
        assert(tail > 0 && tail < 40);                                      // P(|Z| > R) = 2,6·10⁻⁴
        wsl::normal<float> const F({.mean = -3, .sd = 0.25f});
        std::vector<double> floats;
        for (int i = 0; i < 20000; ++i)
            floats.push_back(sample(F, g));
        assert(kolmogorov_smirnov(wsl::normal({.mean = -3.0, .sd = 0.25}), floats) < 1.95);
        // Tirage exact par la voie lente (σ hors de la voie rapide) : même loi.
        std::vector<double> slow;
        for (int i = 0; i < 300; ++i)
            slow.push_back(sample(T, g));
        assert(kolmogorov_smirnov(T, slow) < 1.95);
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
        wsl::exponential const E({.rate = 0.75});
        wsl::exponential<float> const EF({.rate = 3});
        wsl::normal const N({.mean = 10, .sd = 0.5});
        wsl::normal<float> const NF({.mean = 0, .sd = 1});
        auto bits = [](auto v) { if constexpr (sizeof(v) == 8) return std::bit_cast<std::uint64_t>(v); else return std::uint64_t{std::bit_cast<std::uint32_t>(v)}; };
        ok &= check("exponential", fingerprint([&] { return bits(sample(E, x)); }), 0x93410ac575f17399);
        ok &= check("exponential float", fingerprint([&] { return bits(sample(EF, x)); }), 0x707f86aa03e3abe5);
        ok &= check("exponential fast", fingerprint([&] { return bits(sample(E, x, wsl::fast)); }), 0xf0883bbac7b578e0);
        ok &= check("normal", fingerprint([&] { return bits(sample(N, x)); }), 0xb505dc877374dd0a);
        ok &= check("normal float", fingerprint([&] { return bits(sample(NF, x)); }), 0x1c5f8b545c8869d7);
        ok &= check("normal fast", fingerprint([&] { return bits(sample(N, x, wsl::fast)); }), 0xd599f5eebf8f926d);
        ok &= check("normal functions", fingerprint([&] {
            double const p = wsl::uniform_real(x), y = quantile(N, p);
            return bits(y) ^ bits(cdf(N, y)) ^ (bits(ccdf(N, y)) << 1) ^ (bits(pdf(N, y)) << 2) ^ (bits(log_pdf(N, y)) << 3);
        }, 1 << 13), 0x2c59eb64ebd4bc77);
        ok &= check("exponential functions", fingerprint([&] {
            double const p = wsl::uniform_real(x), y = quantile(E, p);
            return bits(y) ^ bits(cdf(E, y)) ^ (bits(ccdf(E, y)) << 1) ^ (bits(pdf(E, y)) << 2) ^ (bits(log_pdf(E, y)) << 3);
        }, 1 << 13), 0x5c920e5f73d7dffd);
    }
    assert(ok);
    return 0;
}
