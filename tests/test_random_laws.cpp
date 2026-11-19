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
#include <numeric>
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

    // Kolmogorov–Smirnov contre une répartition F donnée (lois sans cdf pour l'instant).
    template <class L, class G, class F>
    double kolmogorov_smirnov_against(L const& law, G& g, F F_, int n = 20000)
    {
        std::vector<double> v;
        for (int i = 0; i < n; ++i)
            v.push_back(double(sample(law, g)));
        std::ranges::sort(v);
        double worst = 0;
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            double const f = F_(v[i]);
            worst = std::max({worst, std::abs(f - double(i) / n), std::abs(f - double(i + 1) / n)});
        }
        return worst * std::sqrt(double(n));
    }

    // χ² d'ajustement d'un tirage discret à pmf, classes d'effectif attendu >= 20
    // (queues regroupées) : écart réduit (χ² − ddl)/√(2 ddl), < 4 attendu.
    template <class L, class G>
    double chi_square(L const& law, G& g, int draws, long long low, long long high)
    {
        std::vector<long long> counts(static_cast<std::size_t>(high - low + 1));
        long long below = 0, above = 0;
        for (int i = 0; i < draws; ++i)
        {
            long long const k = sample(law, g);
            if (k < low)
                ++below;
            else if (k > high)
                ++above;
            else
                ++counts[static_cast<std::size_t>(k - low)];
        }
        double chi = 0, expected_below = 0, expected_above = 0;
        int classes = 0;
        double pooled_expected = 0, pooled_count = 0;
        for (long long k = low; k <= high; ++k)
        {
            pooled_expected += draws * pmf(law, k);
            pooled_count += double(counts[static_cast<std::size_t>(k - low)]);
            if (pooled_expected >= 20)
            {
                chi += (pooled_count - pooled_expected) * (pooled_count - pooled_expected) / pooled_expected;
                ++classes;
                pooled_expected = pooled_count = 0;
            }
        }
        for (long long k = low - 1; k >= 0 && k >= low - 100000; --k)
            expected_below += draws * pmf(law, k);
        double const total = draws - [&] { double e = 0; for (long long k = low; k <= high; ++k) e += draws * pmf(law, k); return e; }();
        expected_above = total - expected_below;
        for (auto [c, e] : {std::pair{double(below + 0) + 0.0, expected_below}, {double(above) + pooled_count, expected_above + pooled_expected}})
            if (e >= 20)
            {
                chi += (c - e) * (c - e) / e;
                ++classes;
            }
        int const freedom = classes - 1;
        return (chi - freedom) / std::sqrt(2.0 * freedom);
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

static_assert(wsl::probability_law<wsl::geometric<>> && wsl::probability_law<wsl::poisson<int>> && wsl::probability_law<wsl::binomial<>>);
static_assert(wsl::has_pmf<wsl::geometric<>> && wsl::has_cdf<wsl::geometric<>> && wsl::has_quantile<wsl::geometric<>>);
static_assert(wsl::has_pmf<wsl::poisson<>> && wsl::has_cdf<wsl::poisson<>> && wsl::has_ccdf<wsl::binomial<>> && wsl::has_quantile<wsl::binomial<>>);
static_assert(wsl::geometric<>::sampling == wsl::exactness::approximate && wsl::poisson<>::sampling == wsl::exactness::approximate);
static_assert([] { wsl::binomial X({.trials = 20, .p = 0.3}); return std::is_same_v<decltype(X), wsl::binomial<std::int64_t>>; }());
static_assert(pmf(wsl::geometric({.p = 0.5}), 1) == 0.25 && quantile(wsl::geometric({.p = 0.5}), 0.75) == 1);
static_assert([] { wsl::pcg64 g(5); return sample(wsl::poisson({.mean = 30.0}), g) > 0 && sample(wsl::binomial({.trials = 100, .p = 0.5}), g) > 0; }());

static_assert(wsl::probability_law<wsl::gamma<>> && wsl::probability_law<wsl::beta<float>>);
static_assert(wsl::probability_law<wsl::chi_squared<>> && wsl::probability_law<wsl::student_t<>> && wsl::probability_law<wsl::lognormal<float>> && wsl::probability_law<wsl::cauchy<>>);
static_assert(!wsl::has_mean<wsl::cauchy<>> && !wsl::has_variance<wsl::cauchy<>> && wsl::has_median<wsl::cauchy<>> && wsl::has_quantile<wsl::student_t<>>);
static_assert([] { wsl::student_t X({.degrees_of_freedom = 3}); wsl::cauchy C({.location = 0, .scale = 1}); return std::is_same_v<decltype(X), wsl::student_t<double>> && std::is_same_v<decltype(C), wsl::cauchy<double>>; }());
static_assert(quantile(wsl::cauchy({.location = 2.0, .scale = 3.0}), 0.5) == 2 && cdf(wsl::student_t({.degrees_of_freedom = 4.0}), 0.0) == 0.5);
static_assert(wsl::probability_law<wsl::fisher_f<>> && wsl::probability_law<wsl::weibull<float>> && wsl::probability_law<wsl::laplace<>> && wsl::probability_law<wsl::logistic<>> && wsl::probability_law<wsl::pareto<>>);
static_assert(wsl::laplace<>::sampling == wsl::exactness::exact && wsl::has_entropy<wsl::fisher_f<>> && wsl::has_quantile<wsl::pareto<>>);
static_assert(cdf(wsl::laplace({.location = 1.0, .scale = 2.0}), 1.0) == 0.5 && quantile(wsl::logistic({.location = 3.0, .scale = 1.0}), 0.5) == 3);
static_assert(wsl::has_pdf<wsl::gamma<>> && wsl::has_cdf<wsl::gamma<>> && wsl::has_quantile<wsl::gamma<>> && wsl::has_mode<wsl::beta<>> && wsl::has_entropy<wsl::beta<>> && wsl::has_entropy<wsl::gamma<>>);
static_assert([] { wsl::gamma X({.shape = 2, .rate = 3}); wsl::gamma Y({.shape = 2, .scale = 3}); return std::is_same_v<decltype(X), wsl::gamma<double>> && std::is_same_v<decltype(Y), wsl::gamma<double>>; }());
static_assert([] { wsl::beta X({.alpha = 2, .beta = 5}); return std::is_same_v<decltype(X), wsl::beta<double>>; }());
static_assert(wsl::gamma({.shape = 3.0, .rate = 2.0}).mean() == 1.5 && wsl::gamma({.shape = 3.0, .scale = 2.0}).mean() == 6);
static_assert([] { wsl::pcg64 g(6); return sample(wsl::gamma({.shape = 0.5, .rate = 1.0}), g) > 0 && sample(wsl::beta({.alpha = 2.0, .beta = 0.5}), g) < 1; }());

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

    // ---- geometric ------------------------------------------------------------
    {
        wsl::geometric const G({.p = 0.5});
        assert(G.mean() == 1 && G.variance() == 2 && G.mode() == 0 && G.median() == 0 && G.entropy() == 2 * std::log(2.0));
        assert(pmf(G, 0) == 0.5 && pmf(G, -1) == 0 && cdf(G, 1) == 0.75 && ccdf(G, 1) == 0.25 && cdf(G, -1) == 0);
        assert(quantile(G, 0.75) == 1 && quantile(G, std::nextafter(0.75, 1.0)) == 2 && quantile(G, 0.0) == 0);
        assert(quantile(G, 1.0) == std::numeric_limits<std::int64_t>::max() && log_pmf(G, 3) == -4 * std::log(2.0));
        wsl::geometric const one({.p = 1});
        assert(pmf(one, 0) == 1 && pmf(one, 1) == 0 && cdf(one, 0) == 1 && one.mean() == 0 && sample(one, g) == 0);
        wsl::geometric const tiny({.p = 0x1p-60});
        assert(tiny.mean() == 0x1p60 - 1 && cdf(tiny, std::numeric_limits<std::int64_t>::max()) > 0.99);
        assert(throws_domain_error([] { wsl::geometric X({.p = 0}); }) && throws_domain_error([] { wsl::geometric X({.p = 1.5}); }));
        assert(std::abs(chi_square(wsl::geometric({.p = 0.3}), g, 200000, 0, 30)) < 4);
        assert(std::abs(chi_square(wsl::geometric({.p = 0.001}), g, 200000, 0, 5000)) < 4);
    }

    // ---- poisson --------------------------------------------------------------
    {
        wsl::poisson const P({.mean = 3.5});
        assert(P.mean() == 3.5 && P.variance() == 3.5 && P.mode() == 3 && pmf(P, -1) == 0 && log_pmf(P, 0) == -3.5);
        assert(pmf(P, 0) == wsl::special::exp(-3.5));
        double total = 0;
        for (int k = 0; k < 60; ++k)
            total += pmf(P, k);
        assert(std::abs(total - 1) < 1e-15);
        assert(pmf(wsl::poisson({.mean = 1e300}), 5) == 0 && pmf(wsl::poisson({.mean = 1e-300}), 1) == 1e-300);
        assert(throws_domain_error([] { wsl::poisson X({.mean = 0}); }) && throws_domain_error([] { wsl::poisson X({.mean = INFINITY}); }));
        for (double mu : {0.7, 3.5, 9.9, 10.0, 50.0, 1e4, 1e9})              // inversion, puis PTRS
        {
            double const sd = std::sqrt(mu);
            auto const low = static_cast<long long>(std::max(0.0, mu - 6 * sd)), high = static_cast<long long>(mu + 6 * sd + 6);
            assert(std::abs(chi_square(wsl::poisson({.mean = mu}), g, 200000, low, high)) < 4);
        }
    }

    // ---- binomial -------------------------------------------------------------
    {
        wsl::binomial const B({.trials = 20, .p = 0.3});
        assert(B.trials() == 20 && B.mean() == 6 && std::abs(B.variance() - 4.2) < 1e-15 && B.mode() == 6);
        assert(pmf(B, 21) == 0 && pmf(B, -1) == 0 && pmf(B, 20) == wsl::special::pow(0.3, 20.0));
        double total = 0;
        for (int k = 0; k <= 20; ++k)
            total += pmf(B, k);
        assert(std::abs(total - 1) < 1e-15);
        assert(pmf(wsl::binomial({.trials = 5, .p = 0.0}), 0) == 1 && pmf(wsl::binomial({.trials = 5, .p = 1.0}), 5) == 1);
        assert(sample(wsl::binomial({.trials = 5, .p = 1.0}), g) == 5 && sample(wsl::binomial({.trials = 0, .p = 0.5}), g) == 0);
        assert(throws_domain_error([] { wsl::binomial X({.trials = -1, .p = 0.5}); }) && throws_domain_error([] { wsl::binomial X({.trials = 3, .p = 2}); }));
        for (auto [n, p] : {std::pair{20LL, 0.3}, {1000LL, 0.004}, {1000LL, 0.3}, {1000LL, 0.9}, {100000000LL, 0.5}, {50LL, 0.97}})
        {
            double const mean = double(n) * p, sd = std::sqrt(mean * (1 - p));
            auto const low = static_cast<long long>(std::max(0.0, mean - 6 * sd)), high = std::min(n, static_cast<long long>(mean + 6 * sd + 6));
            assert(std::abs(chi_square(wsl::binomial({.trials = n, .p = p}), g, 200000, low, high)) < 4);
        }
    }

    // ---- gamma ----------------------------------------------------------------
    {
        wsl::gamma const X({.shape = 2.5, .rate = 2}), Y({.shape = 2.5, .scale = 0.5});
        assert(X.mean() == 1.25 && X.variance() == 0.625 && X.mode() == 0.75 && X.shape() == 2.5 && X.rate() == 2 && Y.scale() == 0.5);
        assert(Y.mean() == 1.25 && Y.variance() == 0.625 && pdf(X, 1.0) == pdf(Y, 1.0) && log_pdf(X, 0.3) == log_pdf(Y, 0.3));
        assert(pdf(X, -1.0) == 0 && pdf(X, 0.0) == 0 && pdf(X, INFINITY) == 0 && std::isnan(pdf(X, NAN)) && log_pdf(X, -1.0) == -INFINITY);
        assert(pdf(wsl::gamma({.shape = 1.0, .rate = 3.0}), 0.0) == 3 && pdf(wsl::gamma({.shape = 0.5, .rate = 3.0}), 0.0) == INFINITY);
        assert(pdf(wsl::gamma({.shape = 1.0, .rate = 3.0}), 0.5) == pdf(wsl::exponential({.rate = 3.0}), 0.5));
        assert(throws_domain_error([] { wsl::gamma Z({.shape = 0, .rate = 1}); }) && throws_domain_error([] { wsl::gamma Z({.shape = 1, .scale = -1}); }));
        assert(std::isfinite(pdf(wsl::gamma({.shape = 1e-300, .rate = 1.0}), 1e-300)) && pdf(wsl::gamma({.shape = 1e300, .rate = 1.0}), 1e300) > 0);
        // Répartitions en forme close : forme entière (Erlang), χ²₁ = gamma(1/2, taux 1/2).
        auto const erlang3 = [](double x) { double const t = 2 * x; return 1 - std::exp(-t) * (1 + t + t * t / 2); };
        assert(kolmogorov_smirnov_against(wsl::gamma({.shape = 3.0, .rate = 2.0}), g, erlang3) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::gamma({.shape = 3.0, .scale = 0.5}), g, erlang3) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::gamma({.shape = 0.5, .rate = 0.5}), g, [](double x) { return std::erf(std::sqrt(x / 2)); }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::gamma<float>({.shape = 1, .rate = 4}), g, [](double x) { return 1 - std::exp(-4 * x); }) < 1.95);
        double s = 0;
        wsl::gamma const big({.shape = 400.0, .rate = 0.5});
        for (int i = 0; i < 20000; ++i)
            s += sample(big, g);
        assert(std::abs(s / 20000 - 800) < 4 * std::sqrt(1600.0 / 20000));
    }

    // ---- beta -----------------------------------------------------------------
    {
        wsl::beta const B({.alpha = 2, .beta = 5});
        assert(B.a() == 2 && B.b() == 5 && std::abs(B.mean() - 2.0 / 7) < 1e-16 && B.mode() == 0.2 && std::abs(B.variance() - 10.0 / 392) < 1e-17);
        assert(pdf(B, 0.0) == 0 && pdf(B, 1.0) == 0 && pdf(B, -0.5) == 0 && std::abs(pdf(B, 0.3) - 2.1609) < 1e-15);
        assert(pdf(wsl::beta({.alpha = 1.0, .beta = 3.0}), 0.0) == 3 && pdf(wsl::beta({.alpha = 0.5, .beta = 3.0}), 0.0) == INFINITY);
        assert(pdf(wsl::beta({.alpha = 1.0, .beta = 1.0}), 0.7) == 1 && std::isnan(wsl::beta({.alpha = 0.5, .beta = 0.5}).mode()));
        assert(throws_domain_error([] { wsl::beta Z({.alpha = 0, .beta = 1}); }));
        double const pi = 3.141592653589793;
        assert(kolmogorov_smirnov_against(wsl::beta({.alpha = 0.5, .beta = 0.5}), g, [&](double x) { return 2 / pi * std::asin(std::sqrt(x)); }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::beta({.alpha = 3.0, .beta = 1.0}), g, [](double x) { return x * x * x; }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::beta({.alpha = 1.0, .beta = 0.2}), g, [](double x) { return 1 - std::pow(1 - x, 0.2); }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::beta({.alpha = 0.05, .beta = 1.0}), g, [](double x) { return std::pow(x, 0.05); }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::beta<float>({.alpha = 2, .beta = 1}), g, [](double x) { return x * x; }) < 1.95);
    }

    // ---- Répartitions et quantiles (phase 4) : cohérences ---------------------------
    {
        wsl::poisson const P({.mean = 3.5});
        wsl::binomial const B({.trials = 20, .p = 0.3});
        double sp = 0, sb = 0;
        for (int k = 0; k <= 20; ++k)
        {
            sp += pmf(P, k);
            sb += pmf(B, k);
            assert(std::abs(cdf(P, k) - sp) < 1e-15 && std::abs(cdf(B, k) - sb) < 1e-15);
            assert(std::abs(cdf(P, k) + ccdf(P, k) - 1) < 1e-15 && std::abs(cdf(B, k) + ccdf(B, k) - 1) < 1e-15);
        }
        assert(cdf(P, -1) == 0 && ccdf(P, -1) == 1 && cdf(B, 20) == 1 && ccdf(B, 25) == 0);
        for (double u : {1e-300, 1e-9, 0.1, 0.5, 0.75, 0.999999})
        {
            long long const k = quantile(P, u), j = quantile(B, u);          // quantile ∘ cdf : plus petit k
            assert(cdf(P, k) >= u && (k == 0 || cdf(P, k - 1) < u) && cdf(B, j) >= u && (j == 0 || cdf(B, j - 1) < u));
        }
        assert(quantile(P, 0.0) == 0 && quantile(B, 1.0) == 20 && quantile(wsl::binomial({.trials = 5, .p = 1.0}), 0.3) == 5);
        wsl::gamma const G({.shape = 2.5, .rate = 3});
        wsl::beta const Be({.alpha = 2, .beta = 5});
        for (double u : {1e-200, 1e-6, 0.3, 0.5, 0.9, 1 - 1e-10})
        {
            double const x = quantile(G, u), y = quantile(Be, u);
            assert(std::abs(cdf(G, x) / u - 1) < 1e-14 && std::abs(cdf(Be, y) / u - 1) < 1e-14);
        }
        assert(cdf(G, 0.0) == 0 && ccdf(G, INFINITY) == 0 && quantile(G, 1.0) == INFINITY && cdf(Be, 1.0) == 1 && quantile(Be, 0.0) == 0);
        assert(std::abs(wsl::gamma({.shape = 1.0, .rate = 2.0}).entropy() - wsl::exponential({.rate = 2.0}).entropy()) < 1e-16);
        assert(std::abs(wsl::beta({.alpha = 1.0, .beta = 1.0}).entropy()) < 1e-20);
        assert(std::abs(cdf(wsl::gamma({.shape = 1.0, .rate = 2.0}), 0.7) - cdf(wsl::exponential({.rate = 2.0}), 0.7)) < 1e-16);
    }

    // ---- Lois dérivées -----------------------------------------------------------
    {
        wsl::chi_squared const C({.degrees_of_freedom = 4});
        assert(C.mean() == 4 && C.variance() == 8 && C.mode() == 2 && std::abs(cdf(C, 2.0) - (1 - 2 * std::exp(-1.0))) < 1e-16);
        assert(kolmogorov_smirnov(C, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(sample(C, g)); return v; }()) < 1.95);
        wsl::student_t const T1({.degrees_of_freedom = 1}), T5({.degrees_of_freedom = 5});
        assert(std::isnan(T1.mean()) && T5.mean() == 0 && std::abs(T5.variance() - 5.0 / 3) < 1e-16 && std::isinf(wsl::student_t({.degrees_of_freedom = 1.5}).variance()));
        // t à 1 degré de liberté = Cauchy standard.
        wsl::cauchy const K({.location = 0, .scale = 1});
        for (double t : {-1e6, -3.0, -0.25, 0.5, 2.0, 1e12})
        {
            assert(std::abs(cdf(T1, t) - cdf(K, t)) <= 2e-16 * cdf(K, t) && std::abs(pdf(T1, t) / pdf(K, t) - 1) < 3e-16);
            assert(std::abs(ccdf(T5, t) - cdf(T5, -t)) <= 2e-16 * ccdf(T5, t));                   // symétrie
        }
        assert(quantile(T5, 0.25) == -quantile(T5, 0.75));                  // 1 − 0,25 exact
        for (double u : {1e-300, 1e-8, 0.2, 0.5, 0.975})
        {
            double const t = quantile(T5, u);
            assert(u == 0.5 ? t == 0 : std::abs(cdf(T5, t) / u - 1) < 1e-14);
        }
        assert(kolmogorov_smirnov(T5, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(sample(T5, g)); return v; }()) < 1.95);
        assert(kolmogorov_smirnov(wsl::student_t({.degrees_of_freedom = 0.7}), [&] { std::vector<double> v; wsl::student_t const T({.degrees_of_freedom = 0.7}); for (int i = 0; i < 20000; ++i) v.push_back(sample(T, g)); return v; }()) < 1.95);
        wsl::lognormal const L({.mu = 0.3, .sigma = 0.5});
        assert(std::abs(L.median() - std::exp(0.3)) < 3e-16 && std::abs(L.mean() - std::exp(0.3 + 0.125)) < 1e-15 && std::abs(cdf(L, L.median()) - 0.5) < 2e-16);
        assert(cdf(L, 0.0) == 0 && pdf(L, -1.0) == 0 && quantile(L, 0.0) == 0 && quantile(L, 1.0) == INFINITY);
        assert(kolmogorov_smirnov(L, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(sample(L, g)); return v; }()) < 1.95);
        assert(kolmogorov_smirnov(L, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(sample(L, g, wsl::fast)); return v; }()) < 1.95);
        wsl::cauchy const C2({.location = 1, .scale = 2});
        assert(cdf(C2, 3.0) == 0.75 && quantile(C2, 0.75) == 3 && C2.median() == 1 && std::abs(C2.entropy() - std::log(8 * 3.141592653589793)) < 1e-15);
        assert(kolmogorov_smirnov(C2, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(sample(C2, g)); return v; }()) < 1.95);
        assert(throws_domain_error([] { wsl::student_t X({.degrees_of_freedom = 0}); }) && throws_domain_error([] { wsl::cauchy X({.location = 0, .scale = -1}); }));
        assert(throws_domain_error([] { wsl::chi_squared X({.degrees_of_freedom = -2}); }) && throws_domain_error([] { wsl::lognormal X({.mu = 0, .sigma = 0}); }));
    }

    // ---- Lois dérivées, second lot --------------------------------------------------
    {
        auto draws = [&](auto const& law) { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(double(sample(law, g))); return v; };
        wsl::fisher_f const F({.d1 = 5, .d2 = 10});
        assert(F.mean() == 1.25 && std::isnan(wsl::fisher_f({.d1 = 3, .d2 = 2}).mean()) && std::isinf(wsl::fisher_f({.d1 = 3, .d2 = 3}).variance()));
        assert(pdf(wsl::fisher_f({.d1 = 2, .d2 = 5}), 0.0) == 1 && pdf(F, 0.0) == 0 && pdf(wsl::fisher_f({.d1 = 1, .d2 = 5}), 0.0) == INFINITY);
        // F(1, ν) = t(ν)² : P(F <= t²) = P(|T| <= t).
        wsl::student_t const T7({.degrees_of_freedom = 7});
        wsl::fisher_f const F17({.d1 = 1, .d2 = 7});
        for (double t : {0.5, 1.0, 2.5, 9.0})
            assert(std::abs(cdf(F17, t * t) - (1 - 2 * ccdf(T7, t))) < 1e-15);
        assert(kolmogorov_smirnov(F, draws(F)) < 1.95 && kolmogorov_smirnov(wsl::fisher_f({.d1 = 0.6, .d2 = 1.5}), draws(wsl::fisher_f({.d1 = 0.6, .d2 = 1.5}))) < 1.95);
        wsl::weibull const W1({.shape = 1, .scale = 2});                    // = exponentielle de taux 1/2
        wsl::exponential const E2({.rate = 0.5});
        for (double x : {0.1, 1.0, 5.0, 30.0})
            assert(std::abs(cdf(W1, x) - cdf(E2, x)) <= 2e-16 * cdf(E2, x) && std::abs(pdf(W1, x) / pdf(E2, x) - 1) < 3e-16);
        wsl::weibull const W({.shape = 1.5, .scale = 2});
        assert(std::abs(W.median() - quantile(W, 0.5)) == 0 && kolmogorov_smirnov(W, draws(W)) < 1.95);
        assert(wsl::weibull({.shape = 1e4, .scale = 1}).variance() > 0);     // compensation traitée
        wsl::laplace const L({.location = 1, .scale = 2});
        assert(L.variance() == 8 && cdf(L, 1.0) == 0.5 && quantile(L, 0.5) == 1 && kolmogorov_smirnov(L, draws(L)) < 1.95);
        wsl::logistic const Lg({.location = 1, .scale = 2});
        assert(cdf(Lg, 1.0) == 0.5 && std::abs(quantile(Lg, 0.75) - (1 + 2 * std::log(3.0))) < 1e-15 && kolmogorov_smirnov(Lg, draws(Lg)) < 1.95);
        wsl::pareto const P({.scale = 1, .shape = 2.5});
        assert(cdf(P, 1.0) == 0 && pdf(P, 0.5) == 0 && quantile(P, 0.0) == 1 && std::isinf(wsl::pareto({.scale = 1, .shape = 1}).mean()));
        assert(kolmogorov_smirnov(P, draws(P)) < 1.95);
        assert(throws_domain_error([] { wsl::fisher_f X({.d1 = 0, .d2 = 1}); }) && throws_domain_error([] { wsl::weibull X({.shape = 1, .scale = 0}); }));
        assert(throws_domain_error([] { wsl::laplace X({.location = 0, .scale = -1}); }) && throws_domain_error([] { wsl::pareto X({.scale = 0, .shape = 1}); }));
    }

    // ---- Compléments : Gumbel, Fréchet, Rayleigh, triangulaire, arcsinus, Erlang ----
    {
        auto draws = [&](auto const& law) { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(double(sample(law, g))); return v; };
        wsl::gumbel const Gu({.location = 1, .scale = 2});
        wsl::frechet const Fr({.shape = 3, .scale = 2, .location = 0.5});
        wsl::rayleigh const Ra({.scale = 2});
        wsl::triangular const Tr({.low = 0, .mode = 1, .high = 3});
        wsl::arcsine const Ar({.low = -1, .high = 2});
        wsl::erlang const Er({.shape = 4, .rate = 1.5});
        assert(kolmogorov_smirnov(Gu, draws(Gu)) < 1.95 && kolmogorov_smirnov(Fr, draws(Fr)) < 1.95 && kolmogorov_smirnov(Ra, draws(Ra)) < 1.95);
        assert(kolmogorov_smirnov(Tr, draws(Tr)) < 1.95 && kolmogorov_smirnov(Ar, draws(Ar)) < 1.95 && kolmogorov_smirnov(Er, draws(Er)) < 1.95);
        static_assert(wsl::erlang<>::sampling == wsl::exactness::exact && wsl::has_cquantile<wsl::frechet<>> && wsl::has_entropy<wsl::arcsine<float>>);
        wsl::arcsine const A01({.low = 0, .high = 1});                     // = bêta(½, ½)
        wsl::beta const B({.alpha = 0.5, .beta = 0.5});
        for (double x : {1e-9, 0.1, 0.5, 0.73, 0.999})
            assert(std::abs(cdf(A01, x) - cdf(B, x)) <= 2e-16 * cdf(B, x) && std::abs(pdf(A01, x) / pdf(B, x) - 1) < 1e-15);
        wsl::erlang const E1({.shape = 1, .rate = 2});                      // = exponentielle
        for (double x : {0.01, 1.0, 7.0})
            assert(cdf(E1, x) == cdf(wsl::gamma({.shape = 1.0, .rate = 2.0}), x) && std::abs(cdf(E1, x) / cdf(wsl::exponential({.rate = 2.0}), x) - 1) < 3e-16);
        std::vector<double> y;                                              // α log((X − m)/s) : Gumbel(0, 1)
        for (double v : draws(Fr))
            y.push_back(3 * std::log((v - 0.5) / 2));
        assert(kolmogorov_smirnov(wsl::gumbel({.location = 0.0, .scale = 1.0}), y) < 1.95);
        assert(Tr.median() == quantile(Tr, 0.5) && Ar.median() == 0.5 && quantile(Ra, 0.5) == Ra.median() && cquantile(Gu, 0.5) == quantile(Gu, 0.5));
        assert(pdf(Fr, 0.5) == 0 && cdf(Fr, 0.4) == 0 && pdf(Tr, 3.5) == 0 && pdf(Ar, -1.0) == INFINITY && quantile(Tr, 1.0) == 3);
        assert(throws_domain_error([] { wsl::triangular X({.low = 0, .mode = 4, .high = 3}); }) && throws_domain_error([] { wsl::erlang X({.shape = 0, .rate = 1.0}); }));
        assert(throws_domain_error([] { wsl::frechet X({.shape = 0, .scale = 1}); }) && throws_domain_error([] { wsl::arcsine X({.low = 1, .high = 1}); }));
    }

    // ---- Inverse-gaussienne ----
    {
        for (auto [mu, lambda] : {std::pair{1.0, 3.0}, std::pair{2.0, 0.05}, std::pair{1e-2, 1e3}})
        {
            wsl::inverse_gaussian const W({.mean = mu, .shape = lambda});
            std::vector<double> v;
            for (int i = 0; i < 5000; ++i)
                v.push_back(sample(W, g));
            assert(kolmogorov_smirnov(W, v) < 1.95);
            double const m = std::accumulate(v.begin(), v.end(), 0.0) / double(v.size());
            assert(std::abs(m - mu) < 8 * std::sqrt(W.variance() / double(v.size())));
            for (double p : {1e-200, 1e-8, 0.25, 0.5, 0.75, 1 - 1e-9})   // quantiles et répartition se répondent
            {
                double const x = quantile(W, p), y = cquantile(W, p);
                assert(std::abs(cdf(W, x) / p - 1) < 1e-14 || std::abs(ccdf(W, x) / (1 - p) - 1) < 1e-14);
                assert(std::abs(ccdf(W, y) / p - 1) < 1e-14 || std::abs(cdf(W, y) / (1 - p) - 1) < 1e-14);
            }
            assert(W.median() == quantile(W, 0.5) && cdf(W, 0.0) == 0 && ccdf(W, INFINITY) == 0 && pdf(W, -1.0) == 0);
            double area = 0;                                                // Simpson sur [μ, 1,5 μ]
            for (int i = 0; i < 2000; ++i)
            {
                double const x0 = mu * (1 + i / 4000.0), x1 = mu * (1 + (i + 1) / 4000.0);
                area += (x1 - x0) / 6 * (pdf(W, x0) + 4 * pdf(W, (x0 + x1) / 2) + pdf(W, x1));
            }
            assert(std::abs((cdf(W, 1.5 * mu) - cdf(W, mu)) / area - 1) < 1e-9);
        }
        wsl::inverse_gaussian const W({.mean = 1.0, .shape = 0.5});        // queue droite : ccdf(x) ~ e^{−λx/(2μ²)}
        assert(ccdf(W, 1e3) > 0 && ccdf(W, 1e3) < 1e-100 && ccdf(W, 3e3) == 0);
        assert(throws_domain_error([] { wsl::inverse_gaussian X({.mean = 0.0, .shape = 1.0}); }) && throws_domain_error([] { wsl::inverse_gaussian X({.mean = 1.0, .shape = INFINITY}); }));
    }

    // ---- Répartition de von Mises ----
    {
        for (double kappa : {0.0, 0.3, 4.0, 200.0, 1e6})
        {
            wsl::von_mises const V({.mu = 1.0, .kappa = kappa});
            std::vector<double> v;
            for (int i = 0; i < 2000; ++i)
                v.push_back(sample(V, g));
            assert(kolmogorov_smirnov(V, v) < 1.95);
            assert(cdf(V, 1.0) == 0.5 && ccdf(V, 1.0) == 0.5 && quantile(V, 0.5) == 1.0);
            assert(cdf(V, 1.0 - 3.2) == 0 && cdf(V, 1.0 + 3.2) == 1 && quantile(V, 0.0) == 1.0 - 3.141592653589793);
            for (double p : {1e-300, 1e-9, 0.1, 0.4, 0.75, 1 - 1e-12})
            {
                double const x = quantile(V, p);                            // F(x⁻) <= p <= F(x⁺), à un ulp près
                double const below = std::nextafter(std::nextafter(x, -INFINITY), -INFINITY), above = std::nextafter(std::nextafter(x, INFINITY), INFINITY);
                if (p <= 0.5)
                    assert(cdf(V, below) <= p && cdf(V, above) >= p);
                else
                    assert(ccdf(V, below) >= 1 - p && ccdf(V, above) <= 1 - p);
                assert(p < 0.5 || cquantile(V, 1 - p) == x);
            }
            for (double x : {-1.5, 0.2, 0.99, 1.7, 3.9})                    // F(x) + S(x) = 1, et symétrie autour de μ
                assert(std::abs(cdf(V, x) + ccdf(V, x) - 1) < 1e-15 && cdf(V, 1.0 - (x - 1.0)) == ccdf(V, x));
        }
        wsl::von_mises const U({.mu = 0.0, .kappa = 0.0});                 // κ = 0 : uniforme sur [−π, π]
        for (double x : {-3.0, -1.0, 0.5, 2.0})
            assert(std::abs(cdf(U, x) - (x + 3.141592653589793) / (2 * 3.141592653589793)) < 1e-16);
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
        }, 1 << 13), 0x16ac969b0974669b);
        ok &= check("exponential functions", fingerprint([&] {
            double const p = wsl::uniform_real(x), y = quantile(E, p);
            return bits(y) ^ bits(cdf(E, y)) ^ (bits(ccdf(E, y)) << 1) ^ (bits(pdf(E, y)) << 2) ^ (bits(log_pdf(E, y)) << 3);
        }, 1 << 13), 0x7fd190e1632d6e27);
        wsl::geometric const Gf({.p = 0.05});
        wsl::poisson const Ps({.mean = 4.25}), Pl({.mean = 812.5});
        wsl::binomial const Bs({.trials = 30, .p = 0.2}), Bl({.trials = 100000, .p = 0.7});
        ok &= check("geometric", fingerprint([&] { return static_cast<std::uint64_t>(sample(Gf, x)); }), 0xf0cebc4b03f0655d);
        ok &= check("poisson inversion", fingerprint([&] { return static_cast<std::uint64_t>(sample(Ps, x)); }), 0xab8d234f51d7332e);
        ok &= check("poisson PTRS", fingerprint([&] { return static_cast<std::uint64_t>(sample(Pl, x)); }), 0x4bfc82d01b403440);
        ok &= check("binomial inversion", fingerprint([&] { return static_cast<std::uint64_t>(sample(Bs, x)); }), 0x758aa2af1860fcd7);
        ok &= check("binomial BTRS", fingerprint([&] { return static_cast<std::uint64_t>(sample(Bl, x)); }), 0xc68cde3109166ca4);
        ok &= check("discrete functions", fingerprint([&] {
            auto const k = static_cast<long long>(wsl::uniform_below(x, 2000));
            return bits(pmf(Gf, k)) ^ (bits(cdf(Gf, k)) << 1) ^ (bits(pmf(Pl, k)) << 2) ^ (bits(log_pmf(Pl, k)) << 3)
                 ^ (bits(pmf(Bl, k + 69000)) << 4) ^ (bits(log_pmf(Bs, k % 31)) << 5) ^ static_cast<std::uint64_t>(quantile(Gf, wsl::uniform_real(x)));
        }, 1 << 13), 0x994f5c2864bf1649);
        wsl::gamma const Ga({.shape = 3.75, .rate = 1.5}), Ga_s({.shape = 0.4, .scale = 2.0});
        wsl::gamma<float> const Ga_f({.shape = 9, .rate = 1});
        wsl::beta const Be_a({.alpha = 2.5, .beta = 7.0}), Be_s({.alpha = 0.3, .beta = 0.8});
        ok &= check("gamma", fingerprint([&] { return bits(sample(Ga, x)); }), 0xf1ec39e9eacde03c);
        ok &= check("gamma small shape", fingerprint([&] { return bits(sample(Ga_s, x)); }), 0x965ae0d5f9d27c99);
        ok &= check("gamma float", fingerprint([&] { return bits(sample(Ga_f, x)); }), 0xe510068f6c674400);
        ok &= check("beta", fingerprint([&] { return bits(sample(Be_a, x)); }), 0x78abaeb973605116);
        ok &= check("beta small shapes", fingerprint([&] { return bits(sample(Be_s, x)); }), 0xa2504ea65588f633);
        ok &= check("gamma beta functions", fingerprint([&] {
            double const u = wsl::uniform_real(x), y = 20 * u;
            return bits(pdf(Ga, y)) ^ (bits(log_pdf(Ga_s, y)) << 1) ^ (bits(pdf(Be_a, u)) << 2) ^ (bits(log_pdf(Be_s, u)) << 3);
        }, 1 << 13), 0x9d4f36e7cb669b74);
        ok &= check("phase 4 functions", fingerprint([&] {
            double const u = wsl::uniform_real(x);
            auto const k = static_cast<long long>(wsl::uniform_below(x, 2000));
            return bits(cdf(Ga, 10 * u)) ^ (bits(quantile(Ga_s, u)) << 1) ^ (bits(cdf(Be_a, u)) << 2) ^ (bits(quantile(Be_s, u)) << 3)
                 ^ (bits(cdf(Pl, k)) << 4) ^ static_cast<std::uint64_t>(quantile(Bl, u)) ^ (bits(Ga.entropy() + Be_s.entropy()) << 5);
        }, 1 << 10), 0x1a38363ea9e3931d);
        wsl::chi_squared const Cq({.degrees_of_freedom = 7});
        wsl::student_t const St({.degrees_of_freedom = 3.5});
        wsl::lognormal const Ln({.mu = -1, .sigma = 0.8});
        wsl::cauchy const Cy({.location = 5, .scale = 0.5});
        ok &= check("chi_squared", fingerprint([&] { return bits(sample(Cq, x)); }), 0x8b01e521b2059e12);
        ok &= check("student_t", fingerprint([&] { return bits(sample(St, x)); }), 0xb833b44187ecf66a);
        ok &= check("lognormal", fingerprint([&] { return bits(sample(Ln, x)); }), 0xc6350ba1601077cb);
        ok &= check("cauchy", fingerprint([&] { return bits(sample(Cy, x)); }), 0xf82fd9546762bb11);
        ok &= check("derived functions", fingerprint([&] {
            namespace fp = wsl::special::detail;                            // entrées protégées de la contraction
            double const u = wsl::uniform_real(x), t = fp::fmul(20.0, u) - 10;
            return bits(cdf(St, t)) ^ (bits(quantile(St, u)) << 1) ^ (bits(pdf(Ln, fp::fmul(3.0, u))) << 2) ^ (bits(quantile(Ln, u)) << 3)
                 ^ (bits(cdf(Cy, t)) << 4) ^ (bits(quantile(Cy, u)) << 5) ^ (bits(cdf(Cq, fp::fmul(10.0, u))) << 6);
        }, 1 << 10), 0x6e71d9a455e07fde);
        wsl::fisher_f const Ff({.d1 = 4, .d2 = 9});
        wsl::weibull const Wb({.shape = 2.2, .scale = 3});
        wsl::laplace const Lp({.location = -2, .scale = 0.75});
        wsl::logistic const Lo({.location = 4, .scale = 1.5});
        wsl::pareto const Pa({.scale = 2, .shape = 1.7});
        ok &= check("fisher_f", fingerprint([&] { return bits(sample(Ff, x)); }), 0x2d8f53b32c9947ee);
        ok &= check("weibull", fingerprint([&] { return bits(sample(Wb, x)); }), 0x6b695ec07d8c4fb3);
        ok &= check("laplace", fingerprint([&] { return bits(sample(Lp, x)); }), 0x514ae68e514769b0);
        ok &= check("logistic", fingerprint([&] { return bits(sample(Lo, x)); }), 0x0dc308ffc659ed1d);
        ok &= check("pareto", fingerprint([&] { return bits(sample(Pa, x)); }), 0x6a50d8d1248f5578);
        ok &= check("second functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::uniform_real(x), t = fp::fmul(20.0, u) - 10, y = fp::fmul(8.0, u);
            return bits(cdf(Ff, y)) ^ (bits(quantile(Ff, u)) << 1) ^ (bits(pdf(Wb, y)) << 2) ^ (bits(quantile(Wb, u)) << 3)
                 ^ (bits(cdf(Lp, t)) << 4) ^ (bits(quantile(Lo, u)) << 5) ^ (bits(ccdf(Pa, y + 2)) << 6);
        }, 1 << 10), 0x572d7f7b8ca7f626);
        wsl::gumbel const Gu({.location = -1, .scale = 0.5});
        wsl::frechet const Fr({.shape = 2.5, .scale = 1.5});
        wsl::rayleigh const Ra({.scale = 0.7});
        wsl::triangular const Tr({.low = -1, .mode = 0, .high = 4});
        wsl::arcsine const Ar({.low = 0, .high = 1});
        wsl::erlang const Er({.shape = 3, .rate = 0.5});
        ok &= check("complements", fingerprint([&] {
            return bits(sample(Gu, x)) ^ (bits(sample(Fr, x)) << 1) ^ (bits(sample(Ra, x)) << 2) ^ (bits(sample(Tr, x)) << 3) ^ (bits(sample(Ar, x)) << 4) ^ (bits(sample(Er, x)) << 5);
        }, 1 << 12), 0x8048e3cc57673536);
        ok &= check("complement functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::uniform_real(x), t = fp::fmul(8.0, u);
            return bits(cdf(Gu, t)) ^ (bits(quantile(Fr, u)) << 1) ^ (bits(pdf(Ra, t)) << 2) ^ (bits(cquantile(Tr, u)) << 3) ^ (bits(cdf(Ar, u)) << 4) ^ (bits(pdf(Fr, t + 0.1)) << 5);
        }, 1 << 10), 0x31ef6ddc23202445);
        wsl::inverse_gaussian const Ig({.mean = 1.5, .shape = 0.8});
        ok &= check("inverse_gaussian", fingerprint([&] { return bits(sample(Ig, x)); }), 0xb65cd37eac96c4a3);
        ok &= check("inverse_gaussian functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::uniform_real(x), t = fp::fmul(40.0, u);
            return bits(cdf(Ig, t)) ^ (bits(ccdf(Ig, t)) << 1) ^ (bits(pdf(Ig, t)) << 2) ^ (bits(quantile(Ig, u)) << 3) ^ (bits(cquantile(Ig, u)) << 4);
        }, 1 << 8), 0x2c052dfac3fa0be4);
        wsl::von_mises const Vm({.mu = 0.25, .kappa = 7.5});
        ok &= check("von_mises functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::uniform_real(x), t = fp::fmul(6.0, u) - 2.75;
            return bits(cdf(Vm, t)) ^ (bits(ccdf(Vm, t)) << 1) ^ (bits(quantile(Vm, u)) << 2) ^ (bits(cquantile(Vm, u)) << 3);
        }, 1 << 8), 0xedcd34abd3d9dd7f);
    }
    assert(ok);
    return 0;
}
