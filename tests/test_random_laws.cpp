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
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <ranges>
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
            double const f = wsl::random::cdf(law, sample[i]);
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
            v.push_back(double(wsl::random::sample(law, g)));
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
            long long const k = wsl::random::sample(law, g);
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
            pooled_expected += draws * wsl::random::pmf(law, k);
            pooled_count += double(counts[static_cast<std::size_t>(k - low)]);
            if (pooled_expected >= 20)
            {
                chi += (pooled_count - pooled_expected) * (pooled_count - pooled_expected) / pooled_expected;
                ++classes;
                pooled_expected = pooled_count = 0;
            }
        }
        for (long long k = low - 1; k >= 0 && k >= low - 100000; --k)
            expected_below += draws * wsl::random::pmf(law, k);
        double const total = draws - [&] { double e = 0; for (long long k = low; k <= high; ++k) e += draws * wsl::random::pmf(law, k); return e; }();
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
static_assert(wsl::random::probability_law<wsl::random::uniform<>> && wsl::random::probability_law<wsl::random::uniform<float>>);
static_assert(wsl::random::probability_law<wsl::random::discrete_uniform<>> && wsl::random::probability_law<wsl::random::bernoulli>);
static_assert(wsl::random::has_mean<wsl::random::uniform<>> && wsl::random::has_entropy<wsl::random::uniform<>> && wsl::random::has_quantile<wsl::random::uniform<>>);
static_assert(wsl::random::has_pdf<wsl::random::uniform<>> && !wsl::random::has_pmf<wsl::random::uniform<>>);
static_assert(wsl::random::has_pmf<wsl::random::discrete_uniform<>> && !wsl::random::has_pdf<wsl::random::discrete_uniform<>>);
static_assert(!wsl::random::has_quantile<wsl::random::bernoulli> && !wsl::random::has_mode<wsl::random::uniform<>>);
static_assert(wsl::random::samplable<wsl::random::uniform<>, wsl::random::xoshiro256pp> && wsl::random::samplable<wsl::random::bernoulli, std::mt19937>);
static_assert(wsl::random::uniform<>::sampling == wsl::random::exactness::approximate && wsl::random::discrete_uniform<>::sampling == wsl::random::exactness::exact);

static_assert([] { wsl::random::uniform X({.low = 2, .high = 6}); return std::is_same_v<decltype(X), wsl::random::uniform<double>>; }());
static_assert(wsl::random::uniform({.low = 2.0, .high = 6.0}).mean() == 4.0);
static_assert(wsl::random::cdf(wsl::random::uniform({.low = 2.0, .high = 6.0}), 3.0) == 0.25);
static_assert(wsl::random::quantile(wsl::random::discrete_uniform({.low = 1, .high = 6}), 0.5) == 3);
static_assert([] { wsl::random::pcg64 g(1); auto const x = wsl::random::sample(wsl::random::discrete_uniform({.low = 1, .high = 6}), g); return x >= 1 && x <= 6; }());

static_assert(wsl::random::probability_law<wsl::random::exponential<>> && wsl::random::probability_law<wsl::random::normal<float>>);
static_assert(wsl::random::has_pdf<wsl::random::exponential<>> && wsl::random::has_mode<wsl::random::exponential<>> && wsl::random::has_quantile<wsl::random::normal<>>);
static_assert(wsl::random::exponential<>::sampling == wsl::random::exactness::exact && wsl::random::normal<>::sampling == wsl::random::exactness::exact);
static_assert([] { wsl::random::normal X({.mean = 0, .sd = 1}); return std::is_same_v<decltype(X), wsl::random::normal<double>>; }());
static_assert(wsl::random::cdf(wsl::random::exponential({.rate = 1.0}), 0.0) == 0 && wsl::random::cdf(wsl::random::normal({.mean = 2.0, .sd = 3.0}), 2.0) == 0.5);
static_assert(wsl::random::quantile(wsl::random::normal({.mean = 2.0, .sd = 3.0}), 0.5) == 2.0);
// Tirages exacts en évaluation constante (réels paresseux, grands entiers, voie lente comprise).
static_assert([] { wsl::random::pcg64 g(3); return wsl::random::sample(wsl::random::exponential({.rate = 2.0}), g) > 0; }());
static_assert([] { wsl::random::pcg64 g(4); double const x = wsl::random::sample(wsl::random::normal({.mean = 0.0, .sd = 0x1p-900}), g); return x == x; }());

static_assert(wsl::random::probability_law<wsl::random::geometric<>> && wsl::random::probability_law<wsl::random::poisson<int>> && wsl::random::probability_law<wsl::random::binomial<>>);
static_assert(wsl::random::has_pmf<wsl::random::geometric<>> && wsl::random::has_cdf<wsl::random::geometric<>> && wsl::random::has_quantile<wsl::random::geometric<>>);
static_assert(wsl::random::has_pmf<wsl::random::poisson<>> && wsl::random::has_cdf<wsl::random::poisson<>> && wsl::random::has_ccdf<wsl::random::binomial<>> && wsl::random::has_quantile<wsl::random::binomial<>>);
static_assert(wsl::random::geometric<>::sampling == wsl::random::exactness::approximate && wsl::random::poisson<>::sampling == wsl::random::exactness::approximate);
static_assert([] { wsl::random::binomial X({.trials = 20, .p = 0.3}); return std::is_same_v<decltype(X), wsl::random::binomial<std::int64_t>>; }());
static_assert(wsl::random::pmf(wsl::random::geometric({.p = 0.5}), 1) == 0.25 && wsl::random::quantile(wsl::random::geometric({.p = 0.5}), 0.75) == 1);
static_assert([] { wsl::random::pcg64 g(5); return wsl::random::sample(wsl::random::poisson({.mean = 30.0}), g) > 0 && wsl::random::sample(wsl::random::binomial({.trials = 100, .p = 0.5}), g) > 0; }());

static_assert(wsl::random::probability_law<wsl::random::gamma<>> && wsl::random::probability_law<wsl::random::beta<float>>);
static_assert(wsl::random::probability_law<wsl::random::chi_squared<>> && wsl::random::probability_law<wsl::random::student_t<>> && wsl::random::probability_law<wsl::random::lognormal<float>> && wsl::random::probability_law<wsl::random::cauchy<>>);
static_assert(!wsl::random::has_mean<wsl::random::cauchy<>> && !wsl::random::has_variance<wsl::random::cauchy<>> && wsl::random::has_median<wsl::random::cauchy<>> && wsl::random::has_quantile<wsl::random::student_t<>>);
static_assert([] { wsl::random::student_t X({.degrees_of_freedom = 3}); wsl::random::cauchy C({.location = 0, .scale = 1}); return std::is_same_v<decltype(X), wsl::random::student_t<double>> && std::is_same_v<decltype(C), wsl::random::cauchy<double>>; }());
static_assert(wsl::random::quantile(wsl::random::cauchy({.location = 2.0, .scale = 3.0}), 0.5) == 2 && wsl::random::cdf(wsl::random::student_t({.degrees_of_freedom = 4.0}), 0.0) == 0.5);
static_assert(wsl::random::probability_law<wsl::random::fisher_f<>> && wsl::random::probability_law<wsl::random::weibull<float>> && wsl::random::probability_law<wsl::random::laplace<>> && wsl::random::probability_law<wsl::random::logistic<>> && wsl::random::probability_law<wsl::random::pareto<>>);
static_assert(wsl::random::laplace<>::sampling == wsl::random::exactness::exact && wsl::random::has_entropy<wsl::random::fisher_f<>> && wsl::random::has_quantile<wsl::random::pareto<>>);
static_assert(wsl::random::cdf(wsl::random::laplace({.location = 1.0, .scale = 2.0}), 1.0) == 0.5 && wsl::random::quantile(wsl::random::logistic({.location = 3.0, .scale = 1.0}), 0.5) == 3);
static_assert(wsl::random::has_pdf<wsl::random::gamma<>> && wsl::random::has_cdf<wsl::random::gamma<>> && wsl::random::has_quantile<wsl::random::gamma<>> && wsl::random::has_mode<wsl::random::beta<>> && wsl::random::has_entropy<wsl::random::beta<>> && wsl::random::has_entropy<wsl::random::gamma<>>);
static_assert([] { wsl::random::gamma X({.shape = 2, .rate = 3}); wsl::random::gamma Y({.shape = 2, .scale = 3}); return std::is_same_v<decltype(X), wsl::random::gamma<double>> && std::is_same_v<decltype(Y), wsl::random::gamma<double>>; }());
static_assert([] { wsl::random::beta X({.alpha = 2, .beta = 5}); return std::is_same_v<decltype(X), wsl::random::beta<double>>; }());
static_assert(wsl::random::gamma({.shape = 3.0, .rate = 2.0}).mean() == 1.5 && wsl::random::gamma({.shape = 3.0, .scale = 2.0}).mean() == 6);
static_assert([] { wsl::random::pcg64 g(6); return wsl::random::sample(wsl::random::gamma({.shape = 0.5, .rate = 1.0}), g) > 0 && wsl::random::sample(wsl::random::beta({.alpha = 2.0, .beta = 0.5}), g) < 1; }());

int main()
{
    wsl::random::xoshiro256pp g(2026);

    // ---- uniform --------------------------------------------------------------
    {
        wsl::random::uniform const X({.low = 2, .high = 6});
        assert(X.mean() == 4 && X.median() == 4 && X.variance() == 16.0 / 12 && X.skewness() == 0);
        assert(X.entropy() == wsl::special::log(4.0));
        assert(wsl::random::pdf(X, 3.0) == 0.25 && wsl::random::pdf(X, 6.0) == 0 && wsl::random::pdf(X, 1.0) == 0 && wsl::random::log_pdf(X, 7.0) == -INFINITY);
        assert(wsl::random::cdf(X, 3.0) == 0.25 && wsl::random::cdf(X, 1.0) == 0 && wsl::random::cdf(X, 9.0) == 1 && wsl::random::ccdf(X, 3.0) == 0.75);
        assert(wsl::random::quantile(X, 0.25) == 3 && wsl::random::quantile(X, 0.0) == 2 && wsl::random::quantile(X, 1.0) == 6);
        assert(wsl::random::cdf(X, 0.1 + 2) + wsl::random::ccdf(X, 0.1 + 2) == 1);
        assert(throws_domain_error([&] { (void)wsl::random::quantile(X, 1.5); }));
        assert(throws_domain_error([] { wsl::random::uniform Y({.low = 1, .high = 1}); }));
        assert(throws_domain_error([] { wsl::random::uniform Y({.low = 0, .high = INFINITY}); }));

        // Intervalle immense : b − a déborde, rien d'infini ni de NaN là où le résultat est fini.
        double const big = std::numeric_limits<double>::max();
        wsl::random::uniform const H({.low = -big, .high = big});
        assert(H.mean() == 0 && std::isinf(H.variance()) && std::abs(H.entropy() - 710.4758600739439) < 1e-12);
        assert(wsl::random::pdf(H, 0.0) > 0 && std::isfinite(wsl::random::pdf(H, 0.0)) && wsl::random::cdf(H, 0.0) == 0.5 && wsl::random::quantile(H, 0.5) == 0);
        assert(wsl::random::quantile(H, 0.75) == big / 2 && wsl::random::ccdf(H, big / 2) == 0.25);

        // Valeurs tirées dans [low, high), loi uniforme (Kolmogorov–Smirnov).
        std::vector<double> values;
        for (int i = 0; i < 20000; ++i)
        {
            double const x = wsl::random::sample(X, g);
            assert(x >= 2 && x < 6);
            values.push_back(x);
        }
        assert(kolmogorov_smirnov(X, values) < 1.95);
        wsl::random::uniform const tiny({.low = 1.0, .high = std::nextafter(1.0, 2.0)});   // un seul flottant possible
        for (int i = 0; i < 100; ++i)
            assert(wsl::random::sample(tiny, g) == 1.0);
        wsl::random::uniform<float> const F({.low = -1, .high = 1});
        for (int i = 0; i < 10000; ++i)
        {
            float const x = wsl::random::sample(F, g);
            assert(x >= -1 && x < 1);
        }
        assert(F.mean() == 0 && wsl::random::cdf(F, 0.5f) == 0.75f);
    }

    // ---- discrete_uniform -----------------------------------------------------
    {
        wsl::random::discrete_uniform const D({.low = 1, .high = 6});
        assert(D.mean() == 3.5 && D.median() == 3 && D.variance() == 35.0 / 12);
        assert(wsl::random::pmf(D, 3) == 1.0 / 6 && wsl::random::pmf(D, 7) == 0 && wsl::random::cdf(D, 3) == 0.5 && wsl::random::ccdf(D, 3) == 0.5 && wsl::random::cdf(D, 6) == 1);
        for (int k = 1; k <= 6; ++k)
        {
            // Quantile exact ; la cdf arrondie le décale au plus d'un cran.
            long long const q = wsl::random::quantile(D, wsl::random::cdf(D, k));
            assert(q == k || q == k + 1);
        }
        assert(wsl::random::quantile(D, 1e-300) == 1 && wsl::random::quantile(D, 1.0) == 6 && wsl::random::quantile(D, 0.5) == 3 && wsl::random::quantile(D, 0.25) == 2);
        assert(wsl::random::quantile(D, 0.5 + 0x1p-52) == 4);                            // juste au-dessus de 3/6
        std::array<int, 7> counts{};
        for (int i = 0; i < 60000; ++i)
            ++counts[static_cast<std::size_t>(wsl::random::sample(D, g))];
        for (int k = 1; k <= 6; ++k)
            assert(std::abs(counts[static_cast<std::size_t>(k)] - 10000) < 450);

        // Pleine plage de 64 bits : n = 2^64.
        wsl::random::discrete_uniform const W({.low = std::numeric_limits<std::int64_t>::min(), .high = std::numeric_limits<std::int64_t>::max()});
        assert(wsl::random::pmf(W, 0) == 0x1p-64 && wsl::random::cdf(W, -1) == 0.5 && W.variance() == 0x1p128 / 12 && W.entropy() == 64 * std::log(2.0));
        assert(wsl::random::quantile(W, 0.5) == -1 && wsl::random::quantile(W, 0x1p-64) == std::numeric_limits<std::int64_t>::min());
        wsl::random::discrete_uniform<std::uint8_t> const B({.low = 0, .high = 255});
        assert(wsl::random::quantile(B, 0.5) == 127 && wsl::random::pmf(B, 255) == 1.0 / 256);
        assert(throws_domain_error([] { wsl::random::discrete_uniform E({.low = 2, .high = 1}); }));
    }

    // ---- bernoulli ------------------------------------------------------------
    {
        wsl::random::bernoulli const B({.p = 0.3});
        assert(B.mean() == 0.3 && std::abs(B.variance() - 0.21) < 1e-16 && wsl::random::pmf(B, true) == 0.3 && wsl::random::cdf(B, false) == 0.7);
        assert(std::abs(B.entropy() - 0.6108643020548935) < 1e-15);
        int hits = 0;
        for (int i = 0; i < 100000; ++i)
            hits += wsl::random::sample(B, g);
        assert(std::abs(hits / 100000.0 - 0.3) < 0.005);
        assert((wsl::random::bernoulli({.p = 0}).support() == wsl::random::interval<bool>{false, false}));
        assert(wsl::random::bernoulli({.p = 1e-300}).variance() == 1e-300 && wsl::random::bernoulli({.p = 0}).entropy() == 0);
        assert(throws_domain_error([] { wsl::random::bernoulli E({.p = 1.5}); }));
    }

    // ---- exponential ----------------------------------------------------------
    {
        wsl::random::exponential const E({.rate = 2});
        assert(E.mean() == 0.5 && E.variance() == 0.25 && E.mode() == 0 && E.skewness() == 2 && E.excess_kurtosis() == 6);
        assert(E.median() == wsl::random::quantile(E, 0.5) && std::abs(E.median() - std::log(2.0) / 2) < 1e-16);
        assert(std::abs(E.entropy() - (1 - std::log(2.0))) < 1e-16);
        assert(wsl::random::pdf(E, 0.0) == 2 && wsl::random::pdf(E, -1.0) == 0 && wsl::random::pdf(E, INFINITY) == 0 && wsl::random::log_pdf(E, -1.0) == -INFINITY);
        assert(wsl::random::cdf(E, 0.0) == 0 && wsl::random::ccdf(E, 0.0) == 1 && wsl::random::cdf(E, INFINITY) == 1 && wsl::random::ccdf(E, INFINITY) == 0);
        assert(wsl::random::quantile(E, 0.0) == 0 && wsl::random::quantile(E, 1.0) == INFINITY && std::isnan(wsl::random::cdf(E, NAN)));
        assert(wsl::random::cdf(E, 1e-300) == 2e-300 && wsl::random::ccdf(E, 360.0) > 0 && wsl::random::ccdf(E, 360.0) < 1e-300);
        assert(throws_domain_error([] { wsl::random::exponential X({.rate = 0}); }));
        assert(throws_domain_error([] { wsl::random::exponential X({.rate = INFINITY}); }));
        assert(throws_domain_error([&] { (void)wsl::random::quantile(E, -0.1); }));
        // Paramètres extrêmes : rien de NaN, densité finie.
        wsl::random::exponential const big({.rate = 0x1p1000}), small({.rate = 0x1p-1000});
        assert(wsl::random::pdf(big, 0x1p-1000) > 0 && std::isfinite(wsl::random::pdf(big, 0x1p-1000)) && wsl::random::cdf(small, 0x1p1000) > 0.63);

        std::vector<double> exact, fast;
        for (int i = 0; i < 20000; ++i)
        {
            exact.push_back(wsl::random::sample(E, g));
            fast.push_back(wsl::random::sample(E, g, wsl::random::fast));
        }
        assert(std::ranges::all_of(exact, [](double x) { return x >= 0; }));
        assert(kolmogorov_smirnov(E, exact) < 1.95 && kolmogorov_smirnov(E, fast) < 1.95);
        wsl::random::exponential<float> const F({.rate = 0.5f});
        std::vector<double> floats;
        for (int i = 0; i < 20000; ++i)
            floats.push_back(wsl::random::sample(F, g));
        assert(kolmogorov_smirnov(wsl::random::exponential({.rate = 0.5}), floats) < 1.95);
    }

    // ---- normal ---------------------------------------------------------------
    {
        wsl::random::normal const N({.mean = 1, .sd = 2});
        assert(N.mean() == 1 && N.median() == 1 && N.mode() == 1 && N.variance() == 4 && N.skewness() == 0 && N.excess_kurtosis() == 0);
        assert(std::abs(N.entropy() - (std::log(2.0) + 0.5 * (1 + std::log(2 * 3.141592653589793)))) < 1e-15);
        assert(std::abs(wsl::random::pdf(N, 1.0) - 0.19947114020071635) < 1e-17 && wsl::random::cdf(N, 1.0) == 0.5 && wsl::random::quantile(N, 0.5) == 1);
        for (double d : {0.375, 1.0, 3.0, 10.0, 40.0})                 // 1 ± d exacts
            assert(wsl::random::cdf(N, 1 - d) == wsl::random::ccdf(N, 1 + d));                        // symétrie, exacte
        assert(std::abs(wsl::random::cdf(N, 3.0) - 0.8413447460685429) < 1e-16 && std::abs(wsl::random::quantile(N, 0.8413447460685429) - 3) < 1e-14);
        assert(wsl::random::quantile(N, 0.0) == -INFINITY && wsl::random::quantile(N, 1.0) == INFINITY && wsl::random::cdf(N, -INFINITY) == 0 && wsl::random::pdf(N, INFINITY) == 0);
        assert(std::isnan(wsl::random::pdf(N, NAN)) && std::isnan(wsl::random::cdf(N, NAN)) && wsl::random::log_pdf(N, INFINITY) == -INFINITY);
        assert(wsl::random::cdf(N, -74.0) > 0 && wsl::random::cdf(N, -74.0) < 1e-307 && wsl::random::ccdf(N, 1e6) == 0);
        assert(throws_domain_error([] { wsl::random::normal X({.mean = 0, .sd = 0}); }));
        assert(throws_domain_error([] { wsl::random::normal X({.mean = INFINITY, .sd = 1}); }));
        // Paramètres extrêmes.
        double const big = std::numeric_limits<double>::max();
        wsl::random::normal const H({.mean = 0, .sd = big}), T({.mean = 0, .sd = 1e-300});
        assert(std::abs(wsl::random::cdf(H, big) - 0.8413447460685429) < 1e-16 && wsl::random::quantile(H, 0.975) == INFINITY && std::isfinite(wsl::random::quantile(H, 0.75)));
        assert(std::isfinite(wsl::random::pdf(T, 0.0)) && wsl::random::pdf(T, 0.0) > 1e299 && wsl::random::cdf(T, 1e-300) == wsl::random::cdf(N, 3.0));
        assert(wsl::random::quantile(wsl::random::normal({.mean = big, .sd = big}), 0.01) < big);   // μ + σ z sur les moitiés

        std::vector<double> exact, fast;
        int tail = 0;
        for (int i = 0; i < 40000; ++i)
        {
            exact.push_back(wsl::random::sample(N, g));
            fast.push_back(wsl::random::sample(N, g, wsl::random::fast));
            tail += std::abs(fast.back() - 1) > 2 * 3.6541528853610088;     // queue de la ziggourat
        }
        assert(kolmogorov_smirnov(N, exact) < 1.95 && kolmogorov_smirnov(N, fast) < 1.95);
        assert(tail > 0 && tail < 40);                                      // P(|Z| > R) = 2,6·10⁻⁴
        wsl::random::normal<float> const F({.mean = -3, .sd = 0.25f});
        std::vector<double> floats;
        for (int i = 0; i < 20000; ++i)
            floats.push_back(wsl::random::sample(F, g));
        assert(kolmogorov_smirnov(wsl::random::normal({.mean = -3.0, .sd = 0.25}), floats) < 1.95);
        // Tirage exact par la voie lente (σ hors de la voie rapide) : même loi.
        std::vector<double> slow;
        for (int i = 0; i < 300; ++i)
            slow.push_back(wsl::random::sample(T, g));
        assert(kolmogorov_smirnov(T, slow) < 1.95);
    }

    // ---- geometric ------------------------------------------------------------
    {
        wsl::random::geometric const G({.p = 0.5});
        assert(G.mean() == 1 && G.variance() == 2 && G.mode() == 0 && G.median() == 0 && G.entropy() == 2 * std::log(2.0));
        assert(wsl::random::pmf(G, 0) == 0.5 && wsl::random::pmf(G, -1) == 0 && wsl::random::cdf(G, 1) == 0.75 && wsl::random::ccdf(G, 1) == 0.25 && wsl::random::cdf(G, -1) == 0);
        assert(wsl::random::quantile(G, 0.75) == 1 && wsl::random::quantile(G, std::nextafter(0.75, 1.0)) == 2 && wsl::random::quantile(G, 0.0) == 0);
        assert(wsl::random::quantile(G, 1.0) == std::numeric_limits<std::int64_t>::max() && wsl::random::log_pmf(G, 3) == -4 * std::log(2.0));
        wsl::random::geometric const one({.p = 1});
        assert(wsl::random::pmf(one, 0) == 1 && wsl::random::pmf(one, 1) == 0 && wsl::random::cdf(one, 0) == 1 && one.mean() == 0 && wsl::random::sample(one, g) == 0);
        wsl::random::geometric const tiny({.p = 0x1p-60});
        assert(tiny.mean() == 0x1p60 - 1 && wsl::random::cdf(tiny, std::numeric_limits<std::int64_t>::max()) > 0.99);
        assert(throws_domain_error([] { wsl::random::geometric X({.p = 0}); }) && throws_domain_error([] { wsl::random::geometric X({.p = 1.5}); }));
        assert(std::abs(chi_square(wsl::random::geometric({.p = 0.3}), g, 200000, 0, 30)) < 4);
        assert(std::abs(chi_square(wsl::random::geometric({.p = 0.001}), g, 200000, 0, 5000)) < 4);
    }

    // ---- poisson --------------------------------------------------------------
    {
        wsl::random::poisson const P({.mean = 3.5});
        assert(P.mean() == 3.5 && P.variance() == 3.5 && P.mode() == 3 && wsl::random::pmf(P, -1) == 0 && wsl::random::log_pmf(P, 0) == -3.5);
        assert(wsl::random::pmf(P, 0) == wsl::special::exp(-3.5));
        double total = 0;
        for (int k = 0; k < 60; ++k)
            total += wsl::random::pmf(P, k);
        assert(std::abs(total - 1) < 1e-15);
        assert(wsl::random::pmf(wsl::random::poisson({.mean = 1e300}), 5) == 0 && wsl::random::pmf(wsl::random::poisson({.mean = 1e-300}), 1) == 1e-300);
        assert(throws_domain_error([] { wsl::random::poisson X({.mean = 0}); }) && throws_domain_error([] { wsl::random::poisson X({.mean = INFINITY}); }));
        for (double mu : {0.7, 3.5, 9.9, 10.0, 50.0, 1e4, 1e9})              // inversion, puis PTRS
        {
            double const sd = std::sqrt(mu);
            auto const low = static_cast<long long>(std::max(0.0, mu - 6 * sd)), high = static_cast<long long>(mu + 6 * sd + 6);
            assert(std::abs(chi_square(wsl::random::poisson({.mean = mu}), g, 200000, low, high)) < 4);
        }
    }

    // ---- binomial -------------------------------------------------------------
    {
        wsl::random::binomial const B({.trials = 20, .p = 0.3});
        assert(B.trials() == 20 && B.mean() == 6 && std::abs(B.variance() - 4.2) < 1e-15 && B.mode() == 6);
        assert(wsl::random::pmf(B, 21) == 0 && wsl::random::pmf(B, -1) == 0 && wsl::random::pmf(B, 20) == wsl::special::pow(0.3, 20.0));
        double total = 0;
        for (int k = 0; k <= 20; ++k)
            total += wsl::random::pmf(B, k);
        assert(std::abs(total - 1) < 1e-15);
        assert(wsl::random::pmf(wsl::random::binomial({.trials = 5, .p = 0.0}), 0) == 1 && wsl::random::pmf(wsl::random::binomial({.trials = 5, .p = 1.0}), 5) == 1);
        assert(wsl::random::sample(wsl::random::binomial({.trials = 5, .p = 1.0}), g) == 5 && wsl::random::sample(wsl::random::binomial({.trials = 0, .p = 0.5}), g) == 0);
        assert(throws_domain_error([] { wsl::random::binomial X({.trials = -1, .p = 0.5}); }) && throws_domain_error([] { wsl::random::binomial X({.trials = 3, .p = 2}); }));
        for (auto [n, p] : {std::pair{20LL, 0.3}, {1000LL, 0.004}, {1000LL, 0.3}, {1000LL, 0.9}, {100000000LL, 0.5}, {50LL, 0.97}})
        {
            double const mean = double(n) * p, sd = std::sqrt(mean * (1 - p));
            auto const low = static_cast<long long>(std::max(0.0, mean - 6 * sd)), high = std::min(n, static_cast<long long>(mean + 6 * sd + 6));
            assert(std::abs(chi_square(wsl::random::binomial({.trials = n, .p = p}), g, 200000, low, high)) < 4);
        }
    }

    // ---- gamma ----------------------------------------------------------------
    {
        wsl::random::gamma const X({.shape = 2.5, .rate = 2}), Y({.shape = 2.5, .scale = 0.5});
        assert(X.mean() == 1.25 && X.variance() == 0.625 && X.mode() == 0.75 && X.shape() == 2.5 && X.rate() == 2 && Y.scale() == 0.5);
        assert(Y.mean() == 1.25 && Y.variance() == 0.625 && wsl::random::pdf(X, 1.0) == wsl::random::pdf(Y, 1.0) && wsl::random::log_pdf(X, 0.3) == wsl::random::log_pdf(Y, 0.3));
        assert(wsl::random::pdf(X, -1.0) == 0 && wsl::random::pdf(X, 0.0) == 0 && wsl::random::pdf(X, INFINITY) == 0 && std::isnan(wsl::random::pdf(X, NAN)) && wsl::random::log_pdf(X, -1.0) == -INFINITY);
        assert(wsl::random::pdf(wsl::random::gamma({.shape = 1.0, .rate = 3.0}), 0.0) == 3 && wsl::random::pdf(wsl::random::gamma({.shape = 0.5, .rate = 3.0}), 0.0) == INFINITY);
        assert(wsl::random::pdf(wsl::random::gamma({.shape = 1.0, .rate = 3.0}), 0.5) == wsl::random::pdf(wsl::random::exponential({.rate = 3.0}), 0.5));
        assert(throws_domain_error([] { wsl::random::gamma Z({.shape = 0, .rate = 1}); }) && throws_domain_error([] { wsl::random::gamma Z({.shape = 1, .scale = -1}); }));
        assert(std::isfinite(wsl::random::pdf(wsl::random::gamma({.shape = 1e-300, .rate = 1.0}), 1e-300)) && wsl::random::pdf(wsl::random::gamma({.shape = 1e300, .rate = 1.0}), 1e300) > 0);
        // Répartitions en forme close : forme entière (Erlang), χ²₁ = gamma(1/2, taux 1/2).
        auto const erlang3 = [](double x) { double const t = 2 * x; return 1 - std::exp(-t) * (1 + t + t * t / 2); };
        assert(kolmogorov_smirnov_against(wsl::random::gamma({.shape = 3.0, .rate = 2.0}), g, erlang3) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::random::gamma({.shape = 3.0, .scale = 0.5}), g, erlang3) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::random::gamma({.shape = 0.5, .rate = 0.5}), g, [](double x) { return std::erf(std::sqrt(x / 2)); }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::random::gamma<float>({.shape = 1, .rate = 4}), g, [](double x) { return 1 - std::exp(-4 * x); }) < 1.95);
        double s = 0;
        wsl::random::gamma const big({.shape = 400.0, .rate = 0.5});
        for (int i = 0; i < 20000; ++i)
            s += wsl::random::sample(big, g);
        assert(std::abs(s / 20000 - 800) < 4 * std::sqrt(1600.0 / 20000));
    }

    // ---- beta -----------------------------------------------------------------
    {
        wsl::random::beta const B({.alpha = 2, .beta = 5});
        assert(B.a() == 2 && B.b() == 5 && std::abs(B.mean() - 2.0 / 7) < 1e-16 && B.mode() == 0.2 && std::abs(B.variance() - 10.0 / 392) < 1e-17);
        assert(wsl::random::pdf(B, 0.0) == 0 && wsl::random::pdf(B, 1.0) == 0 && wsl::random::pdf(B, -0.5) == 0 && std::abs(wsl::random::pdf(B, 0.3) - 2.1609) < 1e-15);
        assert(wsl::random::pdf(wsl::random::beta({.alpha = 1.0, .beta = 3.0}), 0.0) == 3 && wsl::random::pdf(wsl::random::beta({.alpha = 0.5, .beta = 3.0}), 0.0) == INFINITY);
        assert(wsl::random::pdf(wsl::random::beta({.alpha = 1.0, .beta = 1.0}), 0.7) == 1 && std::isnan(wsl::random::beta({.alpha = 0.5, .beta = 0.5}).mode()));
        assert(throws_domain_error([] { wsl::random::beta Z({.alpha = 0, .beta = 1}); }));
        double const pi = 3.141592653589793;
        assert(kolmogorov_smirnov_against(wsl::random::beta({.alpha = 0.5, .beta = 0.5}), g, [&](double x) { return 2 / pi * std::asin(std::sqrt(x)); }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::random::beta({.alpha = 3.0, .beta = 1.0}), g, [](double x) { return x * x * x; }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::random::beta({.alpha = 1.0, .beta = 0.2}), g, [](double x) { return 1 - std::pow(1 - x, 0.2); }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::random::beta({.alpha = 0.05, .beta = 1.0}), g, [](double x) { return std::pow(x, 0.05); }) < 1.95);
        assert(kolmogorov_smirnov_against(wsl::random::beta<float>({.alpha = 2, .beta = 1}), g, [](double x) { return x * x; }) < 1.95);
    }

    // ---- Répartitions et quantiles (phase 4) : cohérences ---------------------------
    {
        wsl::random::poisson const P({.mean = 3.5});
        wsl::random::binomial const B({.trials = 20, .p = 0.3});
        double sp = 0, sb = 0;
        for (int k = 0; k <= 20; ++k)
        {
            sp += wsl::random::pmf(P, k);
            sb += wsl::random::pmf(B, k);
            assert(std::abs(wsl::random::cdf(P, k) - sp) < 1e-15 && std::abs(wsl::random::cdf(B, k) - sb) < 1e-15);
            assert(std::abs(wsl::random::cdf(P, k) + wsl::random::ccdf(P, k) - 1) < 1e-15 && std::abs(wsl::random::cdf(B, k) + wsl::random::ccdf(B, k) - 1) < 1e-15);
        }
        assert(wsl::random::cdf(P, -1) == 0 && wsl::random::ccdf(P, -1) == 1 && wsl::random::cdf(B, 20) == 1 && wsl::random::ccdf(B, 25) == 0);
        for (double u : {1e-300, 1e-9, 0.1, 0.5, 0.75, 0.999999})
        {
            long long const k = wsl::random::quantile(P, u), j = wsl::random::quantile(B, u);          // quantile ∘ cdf : plus petit k
            assert(wsl::random::cdf(P, k) >= u && (k == 0 || wsl::random::cdf(P, k - 1) < u) && wsl::random::cdf(B, j) >= u && (j == 0 || wsl::random::cdf(B, j - 1) < u));
        }
        assert(wsl::random::quantile(P, 0.0) == 0 && wsl::random::quantile(B, 1.0) == 20 && wsl::random::quantile(wsl::random::binomial({.trials = 5, .p = 1.0}), 0.3) == 5);
        wsl::random::gamma const G({.shape = 2.5, .rate = 3});
        wsl::random::beta const Be({.alpha = 2, .beta = 5});
        for (double u : {1e-200, 1e-6, 0.3, 0.5, 0.9, 1 - 1e-10})
        {
            double const x = wsl::random::quantile(G, u), y = wsl::random::quantile(Be, u);
            assert(std::abs(wsl::random::cdf(G, x) / u - 1) < 1e-14 && std::abs(wsl::random::cdf(Be, y) / u - 1) < 1e-14);
        }
        assert(wsl::random::cdf(G, 0.0) == 0 && wsl::random::ccdf(G, INFINITY) == 0 && wsl::random::quantile(G, 1.0) == INFINITY && wsl::random::cdf(Be, 1.0) == 1 && wsl::random::quantile(Be, 0.0) == 0);
        assert(std::abs(wsl::random::gamma({.shape = 1.0, .rate = 2.0}).entropy() - wsl::random::exponential({.rate = 2.0}).entropy()) < 1e-16);
        assert(std::abs(wsl::random::beta({.alpha = 1.0, .beta = 1.0}).entropy()) < 1e-20);
        assert(std::abs(wsl::random::cdf(wsl::random::gamma({.shape = 1.0, .rate = 2.0}), 0.7) - wsl::random::cdf(wsl::random::exponential({.rate = 2.0}), 0.7)) < 1e-16);
    }

    // ---- Lois dérivées -----------------------------------------------------------
    {
        wsl::random::chi_squared const C({.degrees_of_freedom = 4});
        assert(C.mean() == 4 && C.variance() == 8 && C.mode() == 2 && std::abs(wsl::random::cdf(C, 2.0) - (1 - 2 * std::exp(-1.0))) < 1e-16);
        assert(kolmogorov_smirnov(C, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wsl::random::sample(C, g)); return v; }()) < 1.95);
        wsl::random::student_t const T1({.degrees_of_freedom = 1}), T5({.degrees_of_freedom = 5});
        assert(std::isnan(T1.mean()) && T5.mean() == 0 && std::abs(T5.variance() - 5.0 / 3) < 1e-16 && std::isinf(wsl::random::student_t({.degrees_of_freedom = 1.5}).variance()));
        // t à 1 degré de liberté = Cauchy standard.
        wsl::random::cauchy const K({.location = 0, .scale = 1});
        for (double t : {-1e6, -3.0, -0.25, 0.5, 2.0, 1e12})
        {
            assert(std::abs(wsl::random::cdf(T1, t) - wsl::random::cdf(K, t)) <= 2e-16 * wsl::random::cdf(K, t) && std::abs(wsl::random::pdf(T1, t) / wsl::random::pdf(K, t) - 1) < 3e-16);
            assert(std::abs(wsl::random::ccdf(T5, t) - wsl::random::cdf(T5, -t)) <= 2e-16 * wsl::random::ccdf(T5, t));                   // symétrie
        }
        assert(wsl::random::quantile(T5, 0.25) == -wsl::random::quantile(T5, 0.75));                  // 1 − 0,25 exact
        for (double u : {1e-300, 1e-8, 0.2, 0.5, 0.975})
        {
            double const t = wsl::random::quantile(T5, u);
            assert(u == 0.5 ? t == 0 : std::abs(wsl::random::cdf(T5, t) / u - 1) < 1e-14);
        }
        assert(kolmogorov_smirnov(T5, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wsl::random::sample(T5, g)); return v; }()) < 1.95);
        assert(kolmogorov_smirnov(wsl::random::student_t({.degrees_of_freedom = 0.7}), [&] { std::vector<double> v; wsl::random::student_t const T({.degrees_of_freedom = 0.7}); for (int i = 0; i < 20000; ++i) v.push_back(wsl::random::sample(T, g)); return v; }()) < 1.95);
        wsl::random::lognormal const L({.mu = 0.3, .sigma = 0.5});
        assert(std::abs(L.median() - std::exp(0.3)) < 3e-16 && std::abs(L.mean() - std::exp(0.3 + 0.125)) < 1e-15 && std::abs(wsl::random::cdf(L, L.median()) - 0.5) < 2e-16);
        assert(wsl::random::cdf(L, 0.0) == 0 && wsl::random::pdf(L, -1.0) == 0 && wsl::random::quantile(L, 0.0) == 0 && wsl::random::quantile(L, 1.0) == INFINITY);
        assert(kolmogorov_smirnov(L, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wsl::random::sample(L, g)); return v; }()) < 1.95);
        assert(kolmogorov_smirnov(L, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wsl::random::sample(L, g, wsl::random::fast)); return v; }()) < 1.95);
        wsl::random::cauchy const C2({.location = 1, .scale = 2});
        assert(wsl::random::cdf(C2, 3.0) == 0.75 && wsl::random::quantile(C2, 0.75) == 3 && C2.median() == 1 && std::abs(C2.entropy() - std::log(8 * 3.141592653589793)) < 1e-15);
        assert(kolmogorov_smirnov(C2, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wsl::random::sample(C2, g)); return v; }()) < 1.95);
        assert(throws_domain_error([] { wsl::random::student_t X({.degrees_of_freedom = 0}); }) && throws_domain_error([] { wsl::random::cauchy X({.location = 0, .scale = -1}); }));
        assert(throws_domain_error([] { wsl::random::chi_squared X({.degrees_of_freedom = -2}); }) && throws_domain_error([] { wsl::random::lognormal X({.mu = 0, .sigma = 0}); }));
    }

    // ---- Lois dérivées, second lot --------------------------------------------------
    {
        auto draws = [&](auto const& law) { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(double(wsl::random::sample(law, g))); return v; };
        wsl::random::fisher_f const F({.d1 = 5, .d2 = 10});
        assert(F.mean() == 1.25 && std::isnan(wsl::random::fisher_f({.d1 = 3, .d2 = 2}).mean()) && std::isinf(wsl::random::fisher_f({.d1 = 3, .d2 = 3}).variance()));
        assert(wsl::random::pdf(wsl::random::fisher_f({.d1 = 2, .d2 = 5}), 0.0) == 1 && wsl::random::pdf(F, 0.0) == 0 && wsl::random::pdf(wsl::random::fisher_f({.d1 = 1, .d2 = 5}), 0.0) == INFINITY);
        // F(1, ν) = t(ν)² : P(F <= t²) = P(|T| <= t).
        wsl::random::student_t const T7({.degrees_of_freedom = 7});
        wsl::random::fisher_f const F17({.d1 = 1, .d2 = 7});
        for (double t : {0.5, 1.0, 2.5, 9.0})
            assert(std::abs(wsl::random::cdf(F17, t * t) - (1 - 2 * wsl::random::ccdf(T7, t))) < 1e-15);
        assert(kolmogorov_smirnov(F, draws(F)) < 1.95 && kolmogorov_smirnov(wsl::random::fisher_f({.d1 = 0.6, .d2 = 1.5}), draws(wsl::random::fisher_f({.d1 = 0.6, .d2 = 1.5}))) < 1.95);
        wsl::random::weibull const W1({.shape = 1, .scale = 2});                    // = exponentielle de taux 1/2
        wsl::random::exponential const E2({.rate = 0.5});
        for (double x : {0.1, 1.0, 5.0, 30.0})
            assert(std::abs(wsl::random::cdf(W1, x) - wsl::random::cdf(E2, x)) <= 2e-16 * wsl::random::cdf(E2, x) && std::abs(wsl::random::pdf(W1, x) / wsl::random::pdf(E2, x) - 1) < 3e-16);
        wsl::random::weibull const W({.shape = 1.5, .scale = 2});
        assert(std::abs(W.median() - wsl::random::quantile(W, 0.5)) == 0 && kolmogorov_smirnov(W, draws(W)) < 1.95);
        assert(wsl::random::weibull({.shape = 1e4, .scale = 1}).variance() > 0);     // compensation traitée
        wsl::random::laplace const L({.location = 1, .scale = 2});
        assert(L.variance() == 8 && wsl::random::cdf(L, 1.0) == 0.5 && wsl::random::quantile(L, 0.5) == 1 && kolmogorov_smirnov(L, draws(L)) < 1.95);
        wsl::random::logistic const Lg({.location = 1, .scale = 2});
        assert(wsl::random::cdf(Lg, 1.0) == 0.5 && std::abs(wsl::random::quantile(Lg, 0.75) - (1 + 2 * std::log(3.0))) < 1e-15 && kolmogorov_smirnov(Lg, draws(Lg)) < 1.95);
        wsl::random::pareto const P({.scale = 1, .shape = 2.5});
        assert(wsl::random::cdf(P, 1.0) == 0 && wsl::random::pdf(P, 0.5) == 0 && wsl::random::quantile(P, 0.0) == 1 && std::isinf(wsl::random::pareto({.scale = 1, .shape = 1}).mean()));
        assert(kolmogorov_smirnov(P, draws(P)) < 1.95);
        assert(throws_domain_error([] { wsl::random::fisher_f X({.d1 = 0, .d2 = 1}); }) && throws_domain_error([] { wsl::random::weibull X({.shape = 1, .scale = 0}); }));
        assert(throws_domain_error([] { wsl::random::laplace X({.location = 0, .scale = -1}); }) && throws_domain_error([] { wsl::random::pareto X({.scale = 0, .shape = 1}); }));
    }

    // ---- Compléments : Gumbel, Fréchet, Rayleigh, triangulaire, arcsinus, Erlang ----
    {
        auto draws = [&](auto const& law) { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(double(wsl::random::sample(law, g))); return v; };
        wsl::random::gumbel const Gu({.location = 1, .scale = 2});
        wsl::random::frechet const Fr({.shape = 3, .scale = 2, .location = 0.5});
        wsl::random::rayleigh const Ra({.scale = 2});
        wsl::random::triangular const Tr({.low = 0, .mode = 1, .high = 3});
        wsl::random::arcsine const Ar({.low = -1, .high = 2});
        wsl::random::erlang const Er({.shape = 4, .rate = 1.5});
        assert(kolmogorov_smirnov(Gu, draws(Gu)) < 1.95 && kolmogorov_smirnov(Fr, draws(Fr)) < 1.95 && kolmogorov_smirnov(Ra, draws(Ra)) < 1.95);
        assert(kolmogorov_smirnov(Tr, draws(Tr)) < 1.95 && kolmogorov_smirnov(Ar, draws(Ar)) < 1.95 && kolmogorov_smirnov(Er, draws(Er)) < 1.95);
        static_assert(wsl::random::erlang<>::sampling == wsl::random::exactness::exact && wsl::random::has_cquantile<wsl::random::frechet<>> && wsl::random::has_entropy<wsl::random::arcsine<float>>);
        wsl::random::arcsine const A01({.low = 0, .high = 1});                     // = bêta(½, ½)
        wsl::random::beta const B({.alpha = 0.5, .beta = 0.5});
        for (double x : {1e-9, 0.1, 0.5, 0.73, 0.999})
            assert(std::abs(wsl::random::cdf(A01, x) - wsl::random::cdf(B, x)) <= 2e-16 * wsl::random::cdf(B, x) && std::abs(wsl::random::pdf(A01, x) / wsl::random::pdf(B, x) - 1) < 1e-15);
        wsl::random::erlang const E1({.shape = 1, .rate = 2});                      // = exponentielle
        for (double x : {0.01, 1.0, 7.0})
            assert(wsl::random::cdf(E1, x) == wsl::random::cdf(wsl::random::gamma({.shape = 1.0, .rate = 2.0}), x) && std::abs(wsl::random::cdf(E1, x) / wsl::random::cdf(wsl::random::exponential({.rate = 2.0}), x) - 1) < 3e-16);
        std::vector<double> y;                                              // α log((X − m)/s) : Gumbel(0, 1)
        for (double v : draws(Fr))
            y.push_back(3 * std::log((v - 0.5) / 2));
        assert(kolmogorov_smirnov(wsl::random::gumbel({.location = 0.0, .scale = 1.0}), y) < 1.95);
        assert(Tr.median() == wsl::random::quantile(Tr, 0.5) && Ar.median() == 0.5 && wsl::random::quantile(Ra, 0.5) == Ra.median() && wsl::random::cquantile(Gu, 0.5) == wsl::random::quantile(Gu, 0.5));
        assert(wsl::random::pdf(Fr, 0.5) == 0 && wsl::random::cdf(Fr, 0.4) == 0 && wsl::random::pdf(Tr, 3.5) == 0 && wsl::random::pdf(Ar, -1.0) == INFINITY && wsl::random::quantile(Tr, 1.0) == 3);
        assert(throws_domain_error([] { wsl::random::triangular X({.low = 0, .mode = 4, .high = 3}); }) && throws_domain_error([] { wsl::random::erlang X({.shape = 0, .rate = 1.0}); }));
        assert(throws_domain_error([] { wsl::random::frechet X({.shape = 0, .scale = 1}); }) && throws_domain_error([] { wsl::random::arcsine X({.low = 1, .high = 1}); }));
    }

    // ---- Inverse-gaussienne ----
    {
        for (auto [mu, lambda] : {std::pair{1.0, 3.0}, std::pair{2.0, 0.05}, std::pair{1e-2, 1e3}})
        {
            wsl::random::inverse_gaussian const W({.mean = mu, .shape = lambda});
            std::vector<double> v;
            for (int i = 0; i < 5000; ++i)
                v.push_back(wsl::random::sample(W, g));
            assert(kolmogorov_smirnov(W, v) < 1.95);
            double const m = std::accumulate(v.begin(), v.end(), 0.0) / double(v.size());
            assert(std::abs(m - mu) < 8 * std::sqrt(W.variance() / double(v.size())));
            for (double p : {1e-200, 1e-8, 0.25, 0.5, 0.75, 1 - 1e-9})   // quantiles et répartition se répondent
            {
                double const x = wsl::random::quantile(W, p), y = wsl::random::cquantile(W, p);
                assert(std::abs(wsl::random::cdf(W, x) / p - 1) < 1e-14 || std::abs(wsl::random::ccdf(W, x) / (1 - p) - 1) < 1e-14);
                assert(std::abs(wsl::random::ccdf(W, y) / p - 1) < 1e-14 || std::abs(wsl::random::cdf(W, y) / (1 - p) - 1) < 1e-14);
            }
            assert(W.median() == wsl::random::quantile(W, 0.5) && wsl::random::cdf(W, 0.0) == 0 && wsl::random::ccdf(W, INFINITY) == 0 && wsl::random::pdf(W, -1.0) == 0);
            double area = 0;                                                // Simpson sur [μ, 1,5 μ]
            for (int i = 0; i < 2000; ++i)
            {
                double const x0 = mu * (1 + i / 4000.0), x1 = mu * (1 + (i + 1) / 4000.0);
                area += (x1 - x0) / 6 * (wsl::random::pdf(W, x0) + 4 * wsl::random::pdf(W, (x0 + x1) / 2) + wsl::random::pdf(W, x1));
            }
            assert(std::abs((wsl::random::cdf(W, 1.5 * mu) - wsl::random::cdf(W, mu)) / area - 1) < 1e-9);
        }
        wsl::random::inverse_gaussian const W({.mean = 1.0, .shape = 0.5});        // queue droite : wsl::random::ccdf(x) ~ e^{−λx/(2μ²)}
        assert(wsl::random::ccdf(W, 1e3) > 0 && wsl::random::ccdf(W, 1e3) < 1e-100 && wsl::random::ccdf(W, 3e3) == 0);
        assert(throws_domain_error([] { wsl::random::inverse_gaussian X({.mean = 0.0, .shape = 1.0}); }) && throws_domain_error([] { wsl::random::inverse_gaussian X({.mean = 1.0, .shape = INFINITY}); }));
    }

    // ---- Répartition de von Mises ----
    {
        for (double kappa : {0.0, 0.3, 4.0, 200.0, 1e6})
        {
            wsl::random::von_mises const V({.mu = 1.0, .kappa = kappa});
            std::vector<double> v;
            for (int i = 0; i < 2000; ++i)
                v.push_back(wsl::random::sample(V, g));
            assert(kolmogorov_smirnov(V, v) < 1.95);
            assert(wsl::random::cdf(V, 1.0) == 0.5 && wsl::random::ccdf(V, 1.0) == 0.5 && wsl::random::quantile(V, 0.5) == 1.0);
            assert(wsl::random::cdf(V, 1.0 - 3.2) == 0 && wsl::random::cdf(V, 1.0 + 3.2) == 1 && wsl::random::quantile(V, 0.0) == 1.0 - 3.141592653589793);
            for (double p : {1e-300, 1e-9, 0.1, 0.4, 0.75, 1 - 1e-12})
            {
                double const x = wsl::random::quantile(V, p);                            // F(x⁻) <= p <= F(x⁺), à un ulp près
                double const below = std::nextafter(std::nextafter(x, -INFINITY), -INFINITY), above = std::nextafter(std::nextafter(x, INFINITY), INFINITY);
                if (p <= 0.5)
                    assert(wsl::random::cdf(V, below) <= p && wsl::random::cdf(V, above) >= p);
                else
                    assert(wsl::random::ccdf(V, below) >= 1 - p && wsl::random::ccdf(V, above) <= 1 - p);
                assert(p < 0.5 || wsl::random::cquantile(V, 1 - p) == x);
            }
            for (double x : {-1.5, 0.2, 0.99, 1.7, 3.9})                    // F(x) + S(x) = 1, et symétrie autour de μ
                assert(std::abs(wsl::random::cdf(V, x) + wsl::random::ccdf(V, x) - 1) < 1e-15 && wsl::random::cdf(V, 1.0 - (x - 1.0)) == wsl::random::ccdf(V, x));
        }
        wsl::random::von_mises const U({.mu = 0.0, .kappa = 0.0});                 // κ = 0 : uniforme sur [−π, π]
        for (double x : {-3.0, -1.0, 0.5, 2.0})
            assert(std::abs(wsl::random::cdf(U, x) - (x + 3.141592653589793) / (2 * 3.141592653589793)) < 1e-16);
    }

    // ---- Binomiale négative, hypergéométrique ----
    {
        wsl::random::negative_binomial const B({.successes = 3.5, .p = 0.4});
        assert(std::abs(chi_square(B, g, 50000, 0, 30)) < 4);
        assert(std::abs(chi_square(wsl::random::negative_binomial({.successes = 0.3, .p = 0.05}), g, 50000, 0, 60)) < 4);
        wsl::random::negative_binomial const G1({.successes = 1.0, .p = 0.3});     // r = 1 : géométrique
        for (long long k : {0LL, 3LL, 40LL})
            assert(std::abs(wsl::random::pmf(G1, k) / wsl::random::pmf(wsl::random::geometric({.p = 0.3}), k) - 1) < 4e-16 && std::abs(wsl::random::cdf(G1, k) / wsl::random::cdf(wsl::random::geometric({.p = 0.3}), k) - 1) < 4e-16);
        for (double u : {1e-12, 0.3, 0.999999})
        {
            long long const k = wsl::random::quantile(B, u);
            assert(wsl::random::cdf(B, k) >= u && (k == 0 || wsl::random::cdf(B, k - 1) < u));
        }
        wsl::random::hypergeometric const H({.population = 50, .successes = 20, .draws = 10});
        assert(std::abs(chi_square(H, g, 50000, 0, 10)) < 4);
        wsl::random::hypergeometric const H2({.population = 100000, .successes = 30000, .draws = 2000});
        assert(std::abs(chi_square(H2, g, 30000, 520, 680)) < 4);
        double total = 0;
        for (long long k = 0; k <= 10; ++k)
            total += wsl::random::pmf(H, k);
        assert(std::abs(total - 1) < 1e-15 && H.mode() == 4 && H.mean() == 4 && wsl::random::cdf(H, 10) == 1 && wsl::random::ccdf(H, -1) == 1);
        for (double u : {1e-9, 0.5, 0.97})
        {
            long long const k = wsl::random::quantile(H2, u);
            assert(wsl::random::cdf(H2, k) >= u && wsl::random::cdf(H2, k - 1) < u);
        }
        wsl::random::hypergeometric const H3({.population = 30, .successes = 29, .draws = 15});   // support {14, 15}
        assert(H3.low() == 14 && wsl::random::pmf(H3, 13) == 0 && std::abs(wsl::random::pmf(H3, 14) - 0.5) < 1e-16);
        assert(throws_domain_error([] { wsl::random::hypergeometric X({.population = 10, .successes = 11, .draws = 2}); }) && throws_domain_error([] { wsl::random::negative_binomial X({.successes = 2, .p = 0}); }));
    }

    // ---- Bêta-binomiale, logarithmique, gaussienne discrète ----
    {
        wsl::random::beta_binomial const B({.trials = 20, .alpha = 2, .beta = 3});
        assert(std::abs(chi_square(B, g, 50000, 0, 20)) < 4);
        assert(std::abs(chi_square(wsl::random::beta_binomial({.trials = 30, .alpha = 0.3, .beta = 0.4}), g, 50000, 0, 30)) < 4);
        wsl::random::beta_binomial const U({.trials = 9, .alpha = 1, .beta = 1});   // α = β = 1 : uniforme sur {0 … n}
        for (long long k = 0; k <= 9; ++k)
            assert(std::abs(wsl::random::pmf(U, k) - 0.1) < 2e-17);
        wsl::random::logarithmic const L({.p = 0.9});
        assert(std::abs(chi_square(L, g, 50000, 1, 40)) < 4 && std::abs(chi_square(wsl::random::logarithmic({.p = 0.05}), g, 50000, 1, 4)) < 4);
        for (double u : {1e-6, 0.5, 0.9, 0.999999})
        {
            long long const k = wsl::random::quantile(L, u);
            assert(wsl::random::cdf(L, k) >= u && (k == 1 || wsl::random::cdf(L, k - 1) < u));
        }
        wsl::random::discrete_gaussian const D({.mu = 0.3, .sigma = 2.5});
        assert(std::abs(chi_square(D, g, 50000, -9, 10)) < 4);
        assert(std::abs(chi_square(wsl::random::discrete_gaussian({.mu = -0.5, .sigma = 0.3}), g, 50000, -2, 1)) < 4);
        assert(std::abs(chi_square(wsl::random::discrete_gaussian({.mu = 7.25, .sigma = 40.0}), g, 50000, -110, 125)) < 4);
        wsl::random::discrete_gaussian const E({.mu = 2.0, .sigma = 3.0});        // μ entier : symétrie, moyenne μ
        assert(E.mean() == 2.0 && wsl::random::cdf(E, 1) == wsl::random::ccdf(E, 2) && wsl::random::pmf(E, 0) == wsl::random::pmf(E, 4));
        for (double u : {1e-12, 0.3, 0.5, 0.97})
        {
            long long const k = wsl::random::quantile(D, u);
            assert(wsl::random::cdf(D, k) >= u && wsl::random::cdf(D, k - 1) < u);
        }
        assert(throws_domain_error([] { wsl::random::logarithmic X({.p = 1.0}); }) && throws_domain_error([] { wsl::random::discrete_gaussian X({.mu = 0, .sigma = 0}); }));
    }

    // ---- Normale asymétrique ----
    {
        for (double a : {0.0, 1.0, -3.0, 20.0})
        {
            wsl::random::skew_normal const X({.location = 0.5, .scale = 2.0, .shape = a});
            std::vector<double> v;
            for (int i = 0; i < 3000; ++i)
                v.push_back(wsl::random::sample(X, g));
            assert(kolmogorov_smirnov(X, v) < 1.95);
            for (double p : {1e-200, 1e-6, 0.3, 0.5, 0.9})
            {
                double const x = wsl::random::quantile(X, p);
                assert(wsl::random::cdf(X, std::nextafter(x, -INFINITY)) <= p * (1 + 1e-15) && wsl::random::cdf(X, std::nextafter(x, INFINITY)) >= p * (1 - 1e-15));
            }
        }
        wsl::random::skew_normal const Z({.shape = 0.0}), One({.shape = 1.0});    // α = 0 : normale ; α = 1 : F = Φ²
        wsl::random::normal const N({.mean = 0.0, .sd = 1.0});
        for (double x : {-20.0, -3.0, 0.4, 5.0})
        {
            assert(std::abs(wsl::random::cdf(Z, x) / wsl::random::cdf(N, x) - 1) < 3e-16 && std::abs(wsl::random::pdf(Z, x) / wsl::random::pdf(N, x) - 1) < 3e-16);
            assert(std::abs(wsl::random::cdf(One, x) / (wsl::random::cdf(N, x) * wsl::random::cdf(N, x)) - 1) < 1e-15);
        }
        assert(throws_domain_error([] { wsl::random::skew_normal X({.scale = 0.0}); }));
    }

    // ---- Lois stables ----
    {
        auto draws = [&](auto const& law, int n) { std::vector<double> v; for (int i = 0; i < n; ++i) v.push_back(wsl::random::sample(law, g)); return v; };
        wsl::random::stable const Gs({.alpha = 2.0, .beta = 0.4, .scale = 1.5});   // normale de variance 2γ²
        wsl::random::stable const Ca({.alpha = 1.0, .beta = 0.0, .scale = 0.5, .location = 1.0});
        wsl::random::stable const Lv({.alpha = 0.5, .beta = 1.0});                  // Lévy : densité √(1/2π) y^{−3/2} e^{−1/(2y)}
        assert(kolmogorov_smirnov(Gs, draws(Gs, 5000)) < 1.95 && kolmogorov_smirnov(Ca, draws(Ca, 5000)) < 1.95 && kolmogorov_smirnov(Lv, draws(Lv, 3000)) < 1.95);
        for (auto [a, b] : {std::pair{1.5, 0.5}, std::pair{0.7, -0.4}, std::pair{1.0, 0.8}})
        {
            wsl::random::stable const X({.alpha = a, .beta = b, .scale = 2.0, .location = -1.0});
            assert(kolmogorov_smirnov(X, draws(X, 200)) < 1.95);                // répartition coûteuse (~1 ms)
            for (double p : {1e-6, 0.8})
            {
                double const x = wsl::random::quantile(X, p);
                assert(wsl::random::cdf(X, std::nextafter(x, -INFINITY)) <= p * (1 + 1e-14) && wsl::random::cdf(X, std::nextafter(x, INFINITY)) >= p * (1 - 1e-14));
            }
        }
        wsl::random::normal const N({.mean = 0.0, .sd = 1.4142135623730951 * 1.5});
        for (double x : {-3.0, 0.5, 4.0})
            assert(std::abs(wsl::random::cdf(Gs, x) / wsl::random::cdf(N, x) - 1) < 1e-15);
        for (double y : {0.1, 1.0, 30.0})
            assert(std::abs(wsl::random::pdf(Lv, y) / (0.3989422804014327 * std::pow(y, -1.5) * std::exp(-0.5 / y)) - 1) < 1e-14);
        assert(wsl::random::cdf(Lv, -1.0) == 0 && wsl::random::pdf(Lv, 0.0) == 0 && Lv.support().lower == 0);
        assert(throws_domain_error([] { wsl::random::stable X({.alpha = 2.5}); }) && throws_domain_error([] { wsl::random::stable X({.alpha = 1.0, .beta = 1.5}); }));
    }

    // ---- Objets-fonctions, densité, vue des tirages ----
    {
        namespace wr = wsl::random;
        wr::normal const N({.mean = 1.0, .sd = 2.0});
        wr::poisson const P({.mean = 3.5});
        static_assert(std::is_same_v<wr::reference_measure_t<wr::normal<>>, wr::lebesgue> && std::is_same_v<wr::reference_measure_t<wr::poisson<>>, wr::counting>);
        static_assert(wr::has_density<wr::normal<>> && wr::has_density<wr::poisson<>> && !wr::has_mean<wr::cauchy<>>);
        assert(wr::pdf(N, 0.5) == wsl::random::pdf(N, 0.5) && wr::density(N, 0.5) == wsl::random::pdf(N, 0.5) && wr::density(P, 2) == wsl::random::pmf(P, 2));
        assert(wr::mean(N) == N.mean() && wr::quantile(N, 0.9) == wsl::random::quantile(N, 0.9) && wr::log_density(P, 4) == wsl::random::log_pmf(P, 4));
        std::vector<double> xs{-1.0, 0.0, 2.5};
        auto const densities = xs | std::views::transform(std::bind_front(wr::pdf, N));   // passable en argument
        assert(std::ranges::equal(densities, std::vector<double>{wsl::random::pdf(N, -1.0), wsl::random::pdf(N, 0.0), wsl::random::pdf(N, 2.5)}));
        wr::xoshiro256pp a(7), b(7);
        std::vector<double> lazy;
        for (double v : wr::samples(N, a) | std::views::take(100))
            lazy.push_back(v);
        for (double v : lazy)
            assert(v == wr::sample(N, b));
        assert(a == b);                                                     // même consommation
        // membres et appel X(g), X(n, g) : mêmes valeurs et même consommation que sample
        assert(N.pdf(0.5) == wr::pdf(N, 0.5) && N.density(0.5) == N.pdf(0.5) && P.density(2) == P.pmf(2) && P.cdf(3) == wr::cdf(P, 3));
        wr::xoshiro256pp c(11), d(11), e(11);
        std::vector<double> const batch = N(5, c);
        for (double v : batch)
            assert(v == N(d));
        for (double v : batch)
            assert(v == wr::sample(N, e));
        assert(c == d && d == e && N(c, wr::fast) == wr::sample(N, d, wr::fast));
    }

    // ---- Tirage rapide demandé : à défaut de méthode propre, le tirage par défaut ----
    {
        wsl::random::xoshiro256pp a(5), b(5);
        wsl::random::uniform const X({.low = 0, .high = 1});
        for (int i = 0; i < 100; ++i)
            assert(wsl::random::sample(X, a, wsl::random::fast) == wsl::random::sample(X, b));
    }

    // ---- Adaptateur vers <random> -------------------------------------------------
    {
        auto d = wsl::random::as_distribution(wsl::random::discrete_uniform({.low = 1, .high = 6}));
        std::vector<long long> rolls(1000);
        std::ranges::generate(rolls, [&] { return d(g); });
        assert(std::ranges::all_of(rolls, [](long long x) { return x >= 1 && x <= 6; }));
        assert(d.min() == 1 && d.max() == 6 && d == d);
        decltype(d)::param_type const other{wsl::random::discrete_uniform({.low = 10, .high = 12})};
        assert(d(g, other) >= 10);
        auto u = wsl::random::as_distribution(wsl::random::uniform({.low = 0.0, .high = 1.0}));
        assert(u.min() == 0.0 && u.max() == 1.0 && u(g) < 1.0);
    }

    // ---- Empreintes ---------------------------------------------------------------
    bool ok = true;
    {
        wsl::random::xoshiro256pp x(11);
        wsl::random::uniform const U({.low = -3.5, .high = 1e10});
        wsl::random::uniform<float> const F({.low = 0, .high = 7});
        wsl::random::discrete_uniform const D({.low = -1000, .high = 1000000});
        wsl::random::bernoulli const B({.p = 0.125});
        ok &= check("uniform", fingerprint([&] { return std::bit_cast<std::uint64_t>(wsl::random::sample(U, x)); }), 0xe9028c88f6ddde7f);
        ok &= check("uniform float", fingerprint([&] { return std::uint64_t{std::bit_cast<std::uint32_t>(wsl::random::sample(F, x))}; }), 0x53a47c22e845da0b);
        ok &= check("discrete_uniform", fingerprint([&] { return static_cast<std::uint64_t>(wsl::random::sample(D, x)); }), 0x3f058c88c4d229e8);
        ok &= check("bernoulli", fingerprint([&] { return std::uint64_t{wsl::random::sample(B, x)}; }), 0x88e84634bf3325c5);
        ok &= check("uniform quantile", fingerprint([&] {
            double const p = wsl::random::uniform_real(x);
            return std::bit_cast<std::uint64_t>(wsl::random::quantile(U, p)) ^ std::bit_cast<std::uint64_t>(wsl::random::cdf(U, wsl::random::quantile(U, p)));
        }), 0x40887993489ea177);
        wsl::random::exponential const E({.rate = 0.75});
        wsl::random::exponential<float> const EF({.rate = 3});
        wsl::random::normal const N({.mean = 10, .sd = 0.5});
        wsl::random::normal<float> const NF({.mean = 0, .sd = 1});
        auto bits = [](auto v) { if constexpr (sizeof(v) == 8) return std::bit_cast<std::uint64_t>(v); else return std::uint64_t{std::bit_cast<std::uint32_t>(v)}; };
        ok &= check("exponential", fingerprint([&] { return bits(wsl::random::sample(E, x)); }), 0x93410ac575f17399);
        ok &= check("exponential float", fingerprint([&] { return bits(wsl::random::sample(EF, x)); }), 0x707f86aa03e3abe5);
        ok &= check("exponential fast", fingerprint([&] { return bits(wsl::random::sample(E, x, wsl::random::fast)); }), 0xf0883bbac7b578e0);
        ok &= check("normal", fingerprint([&] { return bits(wsl::random::sample(N, x)); }), 0xb505dc877374dd0a);
        ok &= check("normal float", fingerprint([&] { return bits(wsl::random::sample(NF, x)); }), 0x1c5f8b545c8869d7);
        ok &= check("normal fast", fingerprint([&] { return bits(wsl::random::sample(N, x, wsl::random::fast)); }), 0xd599f5eebf8f926d);
        ok &= check("normal functions", fingerprint([&] {
            double const p = wsl::random::uniform_real(x), y = wsl::random::quantile(N, p);
            return bits(y) ^ bits(wsl::random::cdf(N, y)) ^ (bits(wsl::random::ccdf(N, y)) << 1) ^ (bits(wsl::random::pdf(N, y)) << 2) ^ (bits(wsl::random::log_pdf(N, y)) << 3);
        }, 1 << 13), 0x16ac969b0974669b);
        ok &= check("exponential functions", fingerprint([&] {
            double const p = wsl::random::uniform_real(x), y = wsl::random::quantile(E, p);
            return bits(y) ^ bits(wsl::random::cdf(E, y)) ^ (bits(wsl::random::ccdf(E, y)) << 1) ^ (bits(wsl::random::pdf(E, y)) << 2) ^ (bits(wsl::random::log_pdf(E, y)) << 3);
        }, 1 << 13), 0x7fd190e1632d6e27);
        wsl::random::geometric const Gf({.p = 0.05});
        wsl::random::poisson const Ps({.mean = 4.25}), Pl({.mean = 812.5});
        wsl::random::binomial const Bs({.trials = 30, .p = 0.2}), Bl({.trials = 100000, .p = 0.7});
        ok &= check("geometric", fingerprint([&] { return static_cast<std::uint64_t>(wsl::random::sample(Gf, x)); }), 0xf0cebc4b03f0655d);
        ok &= check("poisson inversion", fingerprint([&] { return static_cast<std::uint64_t>(wsl::random::sample(Ps, x)); }), 0xab8d234f51d7332e);
        ok &= check("poisson PTRS", fingerprint([&] { return static_cast<std::uint64_t>(wsl::random::sample(Pl, x)); }), 0x4bfc82d01b403440);
        ok &= check("binomial inversion", fingerprint([&] { return static_cast<std::uint64_t>(wsl::random::sample(Bs, x)); }), 0x758aa2af1860fcd7);
        ok &= check("binomial BTRS", fingerprint([&] { return static_cast<std::uint64_t>(wsl::random::sample(Bl, x)); }), 0xc68cde3109166ca4);
        ok &= check("discrete functions", fingerprint([&] {
            auto const k = static_cast<long long>(wsl::random::uniform_below(2000, x));
            return bits(wsl::random::pmf(Gf, k)) ^ (bits(wsl::random::cdf(Gf, k)) << 1) ^ (bits(wsl::random::pmf(Pl, k)) << 2) ^ (bits(wsl::random::log_pmf(Pl, k)) << 3)
                 ^ (bits(wsl::random::pmf(Bl, k + 69000)) << 4) ^ (bits(wsl::random::log_pmf(Bs, k % 31)) << 5) ^ static_cast<std::uint64_t>(wsl::random::quantile(Gf, wsl::random::uniform_real(x)));
        }, 1 << 13), 0x994f5c2864bf1649);
        wsl::random::gamma const Ga({.shape = 3.75, .rate = 1.5}), Ga_s({.shape = 0.4, .scale = 2.0});
        wsl::random::gamma<float> const Ga_f({.shape = 9, .rate = 1});
        wsl::random::beta const Be_a({.alpha = 2.5, .beta = 7.0}), Be_s({.alpha = 0.3, .beta = 0.8});
        ok &= check("gamma", fingerprint([&] { return bits(wsl::random::sample(Ga, x)); }), 0xf1ec39e9eacde03c);
        ok &= check("gamma small shape", fingerprint([&] { return bits(wsl::random::sample(Ga_s, x)); }), 0x965ae0d5f9d27c99);
        ok &= check("gamma float", fingerprint([&] { return bits(wsl::random::sample(Ga_f, x)); }), 0xe510068f6c674400);
        ok &= check("beta", fingerprint([&] { return bits(wsl::random::sample(Be_a, x)); }), 0x78abaeb973605116);
        ok &= check("beta small shapes", fingerprint([&] { return bits(wsl::random::sample(Be_s, x)); }), 0xa2504ea65588f633);
        ok &= check("gamma beta functions", fingerprint([&] {
            double const u = wsl::random::uniform_real(x), y = 20 * u;
            return bits(wsl::random::pdf(Ga, y)) ^ (bits(wsl::random::log_pdf(Ga_s, y)) << 1) ^ (bits(wsl::random::pdf(Be_a, u)) << 2) ^ (bits(wsl::random::log_pdf(Be_s, u)) << 3);
        }, 1 << 13), 0x9d4f36e7cb669b74);
        ok &= check("phase 4 functions", fingerprint([&] {
            double const u = wsl::random::uniform_real(x);
            auto const k = static_cast<long long>(wsl::random::uniform_below(2000, x));
            return bits(wsl::random::cdf(Ga, 10 * u)) ^ (bits(wsl::random::quantile(Ga_s, u)) << 1) ^ (bits(wsl::random::cdf(Be_a, u)) << 2) ^ (bits(wsl::random::quantile(Be_s, u)) << 3)
                 ^ (bits(wsl::random::cdf(Pl, k)) << 4) ^ static_cast<std::uint64_t>(wsl::random::quantile(Bl, u)) ^ (bits(Ga.entropy() + Be_s.entropy()) << 5);
        }, 1 << 10), 0x1a38363ea9e3931d);
        wsl::random::chi_squared const Cq({.degrees_of_freedom = 7});
        wsl::random::student_t const St({.degrees_of_freedom = 3.5});
        wsl::random::lognormal const Ln({.mu = -1, .sigma = 0.8});
        wsl::random::cauchy const Cy({.location = 5, .scale = 0.5});
        ok &= check("chi_squared", fingerprint([&] { return bits(wsl::random::sample(Cq, x)); }), 0x8b01e521b2059e12);
        ok &= check("student_t", fingerprint([&] { return bits(wsl::random::sample(St, x)); }), 0xb833b44187ecf66a);
        ok &= check("lognormal", fingerprint([&] { return bits(wsl::random::sample(Ln, x)); }), 0xc6350ba1601077cb);
        ok &= check("cauchy", fingerprint([&] { return bits(wsl::random::sample(Cy, x)); }), 0xf82fd9546762bb11);
        ok &= check("derived functions", fingerprint([&] {
            namespace fp = wsl::special::detail;                            // entrées protégées de la contraction
            double const u = wsl::random::uniform_real(x), t = fp::fmul(20.0, u) - 10;
            return bits(wsl::random::cdf(St, t)) ^ (bits(wsl::random::quantile(St, u)) << 1) ^ (bits(wsl::random::pdf(Ln, fp::fmul(3.0, u))) << 2) ^ (bits(wsl::random::quantile(Ln, u)) << 3)
                 ^ (bits(wsl::random::cdf(Cy, t)) << 4) ^ (bits(wsl::random::quantile(Cy, u)) << 5) ^ (bits(wsl::random::cdf(Cq, fp::fmul(10.0, u))) << 6);
        }, 1 << 10), 0x6e71d9a455e07fde);
        wsl::random::fisher_f const Ff({.d1 = 4, .d2 = 9});
        wsl::random::weibull const Wb({.shape = 2.2, .scale = 3});
        wsl::random::laplace const Lp({.location = -2, .scale = 0.75});
        wsl::random::logistic const Lo({.location = 4, .scale = 1.5});
        wsl::random::pareto const Pa({.scale = 2, .shape = 1.7});
        ok &= check("fisher_f", fingerprint([&] { return bits(wsl::random::sample(Ff, x)); }), 0x2d8f53b32c9947ee);
        ok &= check("weibull", fingerprint([&] { return bits(wsl::random::sample(Wb, x)); }), 0x6b695ec07d8c4fb3);
        ok &= check("laplace", fingerprint([&] { return bits(wsl::random::sample(Lp, x)); }), 0x514ae68e514769b0);
        ok &= check("logistic", fingerprint([&] { return bits(wsl::random::sample(Lo, x)); }), 0x0dc308ffc659ed1d);
        ok &= check("pareto", fingerprint([&] { return bits(wsl::random::sample(Pa, x)); }), 0x6a50d8d1248f5578);
        ok &= check("second functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::random::uniform_real(x), t = fp::fmul(20.0, u) - 10, y = fp::fmul(8.0, u);
            return bits(wsl::random::cdf(Ff, y)) ^ (bits(wsl::random::quantile(Ff, u)) << 1) ^ (bits(wsl::random::pdf(Wb, y)) << 2) ^ (bits(wsl::random::quantile(Wb, u)) << 3)
                 ^ (bits(wsl::random::cdf(Lp, t)) << 4) ^ (bits(wsl::random::quantile(Lo, u)) << 5) ^ (bits(wsl::random::ccdf(Pa, y + 2)) << 6);
        }, 1 << 10), 0x572d7f7b8ca7f626);
        wsl::random::gumbel const Gu({.location = -1, .scale = 0.5});
        wsl::random::frechet const Fr({.shape = 2.5, .scale = 1.5});
        wsl::random::rayleigh const Ra({.scale = 0.7});
        wsl::random::triangular const Tr({.low = -1, .mode = 0, .high = 4});
        wsl::random::arcsine const Ar({.low = 0, .high = 1});
        wsl::random::erlang const Er({.shape = 3, .rate = 0.5});
        ok &= check("complements", fingerprint([&] {
            return bits(wsl::random::sample(Gu, x)) ^ (bits(wsl::random::sample(Fr, x)) << 1) ^ (bits(wsl::random::sample(Ra, x)) << 2) ^ (bits(wsl::random::sample(Tr, x)) << 3) ^ (bits(wsl::random::sample(Ar, x)) << 4) ^ (bits(wsl::random::sample(Er, x)) << 5);
        }, 1 << 12), 0x8048e3cc57673536);
        ok &= check("complement functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::random::uniform_real(x), t = fp::fmul(8.0, u);
            return bits(wsl::random::cdf(Gu, t)) ^ (bits(wsl::random::quantile(Fr, u)) << 1) ^ (bits(wsl::random::pdf(Ra, t)) << 2) ^ (bits(wsl::random::cquantile(Tr, u)) << 3) ^ (bits(wsl::random::cdf(Ar, u)) << 4) ^ (bits(wsl::random::pdf(Fr, t + 0.1)) << 5);
        }, 1 << 10), 0x31ef6ddc23202445);
        wsl::random::inverse_gaussian const Ig({.mean = 1.5, .shape = 0.8});
        ok &= check("inverse_gaussian", fingerprint([&] { return bits(wsl::random::sample(Ig, x)); }), 0xb65cd37eac96c4a3);
        ok &= check("inverse_gaussian functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::random::uniform_real(x), t = fp::fmul(40.0, u);
            return bits(wsl::random::cdf(Ig, t)) ^ (bits(wsl::random::ccdf(Ig, t)) << 1) ^ (bits(wsl::random::pdf(Ig, t)) << 2) ^ (bits(wsl::random::quantile(Ig, u)) << 3) ^ (bits(wsl::random::cquantile(Ig, u)) << 4);
        }, 1 << 8), 0x2c052dfac3fa0be4);
        wsl::random::von_mises const Vm({.mu = 0.25, .kappa = 7.5});
        ok &= check("von_mises functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::random::uniform_real(x), t = fp::fmul(6.0, u) - 2.75;
            return bits(wsl::random::cdf(Vm, t)) ^ (bits(wsl::random::ccdf(Vm, t)) << 1) ^ (bits(wsl::random::quantile(Vm, u)) << 2) ^ (bits(wsl::random::cquantile(Vm, u)) << 3);
        }, 1 << 8), 0xedcd34abd3d9dd7f);
        wsl::random::negative_binomial const Nb({.successes = 2.5, .p = 0.3});
        wsl::random::hypergeometric const Hg({.population = 500, .successes = 120, .draws = 60});
        ok &= check("discrete lot 3", fingerprint([&] { return bits(double(wsl::random::sample(Nb, x))) ^ (bits(double(wsl::random::sample(Hg, x))) << 1); }, 1 << 12), 0x1abf82d977d052a7);
        ok &= check("discrete lot 3 functions", fingerprint([&] {
            double const u = wsl::random::uniform_real(x);
            long long const k = static_cast<long long>(wsl::random::uniform_below(30, x));
            return bits(wsl::random::pmf(Nb, k)) ^ (bits(wsl::random::cdf(Nb, k)) << 1) ^ (bits(double(wsl::random::quantile(Nb, u))) << 2) ^ (bits(wsl::random::pmf(Hg, k)) << 3) ^ (bits(wsl::random::ccdf(Hg, k)) << 4) ^ (bits(double(wsl::random::quantile(Hg, u))) << 5);
        }, 1 << 9), 0x2fa168fb952bbce9);
        wsl::random::beta_binomial const Bb({.trials = 40, .alpha = 0.7, .beta = 2.5});
        wsl::random::logarithmic const Lg({.p = 0.95});
        wsl::random::discrete_gaussian const Dg({.mu = -1.3, .sigma = 4.2});
        ok &= check("discrete lot 4", fingerprint([&] { return bits(double(wsl::random::sample(Bb, x))) ^ (bits(double(wsl::random::sample(Lg, x))) << 1) ^ (bits(double(wsl::random::sample(Dg, x))) << 2); }, 1 << 12), 0x5715c99366cd16cf);
        ok &= check("discrete lot 4 functions", fingerprint([&] {
            double const u = wsl::random::uniform_real(x);
            long long const k = static_cast<long long>(wsl::random::uniform_below(40, x));
            return bits(wsl::random::pmf(Bb, k)) ^ (bits(wsl::random::cdf(Bb, k)) << 1) ^ (bits(double(wsl::random::quantile(Bb, u))) << 2) ^ (bits(wsl::random::pmf(Lg, k)) << 3) ^ (bits(wsl::random::ccdf(Lg, k)) << 4)
                 ^ (bits(double(wsl::random::quantile(Lg, u))) << 5) ^ (bits(wsl::random::cdf(Dg, k - 20)) << 6) ^ (bits(double(wsl::random::quantile(Dg, u))) << 7);
        }, 1 << 9), 0x26b8cef1aa562a21);
        wsl::random::skew_normal const Sn({.location = -0.5, .scale = 1.5, .shape = 3.0});
        ok &= check("skew_normal", fingerprint([&] { return bits(wsl::random::sample(Sn, x)); }, 1 << 12), 0x646d26c92778b5d5);
        ok &= check("skew_normal functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::random::uniform_real(x), t = fp::fmul(10.0, u) - 4;
            return bits(wsl::random::pdf(Sn, t)) ^ (bits(wsl::random::cdf(Sn, t)) << 1) ^ (bits(wsl::random::ccdf(Sn, t)) << 2) ^ (bits(wsl::random::quantile(Sn, u)) << 3);
        }, 1 << 7), 0x564a5825a3a680fd);
        wsl::random::stable const S1({.alpha = 1.3, .beta = -0.6, .scale = 0.8, .location = 0.2});
        wsl::random::stable const S2({.alpha = 1.0, .beta = 0.5, .scale = 3.0});
        ok &= check("stable", fingerprint([&] { return bits(wsl::random::sample(S1, x)) ^ (bits(wsl::random::sample(S2, x)) << 1); }, 1 << 12), 0x3e15bab78fbb0172);
        ok &= check("stable functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::random::uniform_real(x), t = fp::fmul(10.0, u) - 5;
            return bits(wsl::random::pdf(S1, t)) ^ (bits(wsl::random::cdf(S1, t)) << 1) ^ (bits(wsl::random::ccdf(S2, t)) << 2);
        }, 1 << 5), 0xba852547a6b18b91);
    }
    assert(ok);
    return 0;
}
