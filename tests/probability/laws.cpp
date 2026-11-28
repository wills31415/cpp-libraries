// Lois de wsl::random : concept et attributs, paramètres, valeurs exactes,
// cas extrêmes, lois vérifiées statistiquement (graines fixes), adaptateur vers
// <random>, empreintes des tirages.

#include <wsl/probability.hpp>
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
            double const f = wr::cdf(law, sample[i]);
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
            v.push_back(double(wr::native(law, g)));
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
            long long const k = wr::native(law, g);
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
            pooled_expected += draws * wr::pmf(law, k);
            pooled_count += double(counts[static_cast<std::size_t>(k - low)]);
            if (pooled_expected >= 20)
            {
                chi += (pooled_count - pooled_expected) * (pooled_count - pooled_expected) / pooled_expected;
                ++classes;
                pooled_expected = pooled_count = 0;
            }
        }
        for (long long k = low - 1; k >= low - 100000; --k)              // support négatif compris
        {
            double const p = wr::pmf(law, k);
            expected_below += draws * p;
            if (p < 1e-300 && k < low - 50)
                break;
        }
        double const total = draws - [&] { double e = 0; for (long long k = low; k <= high; ++k) e += draws * wr::pmf(law, k); return e; }();
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
static_assert(wr::probability_law<wr::uniform<>> && wr::probability_law<wr::uniform<float>>);
static_assert(wr::probability_law<wr::discrete_uniform<>> && wr::probability_law<wr::bernoulli>);
static_assert(wr::has_mean<wr::uniform<>> && wr::has_entropy<wr::uniform<>> && wr::has_quantile<wr::uniform<>>);
static_assert(wr::has_pdf<wr::uniform<>> && !wr::has_pmf<wr::uniform<>>);
static_assert(wr::has_pmf<wr::discrete_uniform<>> && !wr::has_pdf<wr::discrete_uniform<>>);
static_assert(!wr::has_quantile<wr::bernoulli> && !wr::has_mode<wr::uniform<>>);
static_assert(wr::samplable<wr::uniform<>, wr::xoshiro256pp> && wr::samplable<wr::bernoulli, std::mt19937>);
static_assert(wr::uniform<>::sampling == wr::exactness::approximate && wr::discrete_uniform<>::sampling == wr::exactness::exact);

static_assert([] { wr::uniform X({.lower = 2, .upper = 6}); return std::is_same_v<decltype(X), wr::uniform<double>>; }());
static_assert(wr::uniform({.lower = 2.0, .upper = 6.0}).mean() == 4.0);
static_assert(wr::cdf(wr::uniform({.lower = 2.0, .upper = 6.0}), 3.0) == 0.25);
static_assert(wr::quantile(wr::discrete_uniform({.lower = 1, .upper = 6}), 0.5) == 3);
static_assert([] { wr::pcg64 g(1); auto const x = wr::native(wr::discrete_uniform({.lower = 1, .upper = 6}), g); return x >= 1 && x <= 6; }());

static_assert(wr::probability_law<wr::exponential<>> && wr::probability_law<wr::normal<float>>);
static_assert(wr::has_pdf<wr::exponential<>> && wr::has_mode<wr::exponential<>> && wr::has_quantile<wr::normal<>>);
static_assert(wr::exponential<>::sampling == wr::exactness::exact && wr::normal<>::sampling == wr::exactness::exact);
static_assert([] { wr::normal X({.mean = 0, .stddev = 1}); return std::is_same_v<decltype(X), wr::normal<double>>; }());
static_assert(wr::cdf(wr::exponential({.rate = 1.0}), 0.0) == 0 && wr::cdf(wr::normal({.mean = 2.0, .stddev = 3.0}), 2.0) == 0.5);
static_assert(wr::quantile(wr::normal({.mean = 2.0, .stddev = 3.0}), 0.5) == 2.0);
// Tirages exacts en évaluation constante (réels paresseux, grands entiers, voie lente comprise).
static_assert([] { wr::pcg64 g(3); return wr::native(wr::exponential({.rate = 2.0}), g) > 0; }());
static_assert([] { wr::pcg64 g(4); double const x = wr::native(wr::normal({.mean = 0.0, .stddev = 0x1p-900}), g); return x == x; }());

static_assert(wr::probability_law<wr::geometric<>> && wr::probability_law<wr::poisson<int>> && wr::probability_law<wr::binomial<>>);
static_assert(wr::has_pmf<wr::geometric<>> && wr::has_cdf<wr::geometric<>> && wr::has_quantile<wr::geometric<>>);
static_assert(wr::has_pmf<wr::poisson<>> && wr::has_cdf<wr::poisson<>> && wr::has_ccdf<wr::binomial<>> && wr::has_quantile<wr::binomial<>>);
static_assert(wr::geometric<>::sampling == wr::exactness::approximate && wr::poisson<>::sampling == wr::exactness::approximate);
static_assert([] { wr::binomial X({.trials = 20, .p = 0.3}); return std::is_same_v<decltype(X), wr::binomial<std::int64_t>>; }());
static_assert(wr::pmf(wr::geometric({.p = 0.5}), 1) == 0.25 && wr::quantile(wr::geometric({.p = 0.5}), 0.75) == 1);
static_assert([] { wr::pcg64 g(5); return wr::native(wr::poisson({.mean = 30.0}), g) > 0 && wr::native(wr::binomial({.trials = 100, .p = 0.5}), g) > 0; }());

static_assert(wr::probability_law<wr::gamma<>> && wr::probability_law<wr::beta<float>>);
static_assert(wr::probability_law<wr::chi_squared<>> && wr::probability_law<wr::student_t<>> && wr::probability_law<wr::lognormal<float>> && wr::probability_law<wr::cauchy<>>);
static_assert(!wr::has_mean<wr::cauchy<>> && !wr::has_variance<wr::cauchy<>> && wr::has_median<wr::cauchy<>> && wr::has_quantile<wr::student_t<>>);
static_assert([] { wr::student_t X({.degrees_of_freedom = 3}); wr::cauchy C({.location = 0, .scale = 1}); return std::is_same_v<decltype(X), wr::student_t<double>> && std::is_same_v<decltype(C), wr::cauchy<double>>; }());
static_assert(wr::quantile(wr::cauchy({.location = 2.0, .scale = 3.0}), 0.5) == 2 && wr::cdf(wr::student_t({.degrees_of_freedom = 4.0}), 0.0) == 0.5);
static_assert(wr::probability_law<wr::fisher_f<>> && wr::probability_law<wr::weibull<float>> && wr::probability_law<wr::laplace<>> && wr::probability_law<wr::logistic<>> && wr::probability_law<wr::pareto<>>);
static_assert(wr::laplace<>::sampling == wr::exactness::exact && wr::has_entropy<wr::fisher_f<>> && wr::has_quantile<wr::pareto<>>);
static_assert(wr::cdf(wr::laplace({.location = 1.0, .scale = 2.0}), 1.0) == 0.5 && wr::quantile(wr::logistic({.location = 3.0, .scale = 1.0}), 0.5) == 3);
static_assert(wr::has_pdf<wr::gamma<>> && wr::has_cdf<wr::gamma<>> && wr::has_quantile<wr::gamma<>> && wr::has_mode<wr::beta<>> && wr::has_entropy<wr::beta<>> && wr::has_entropy<wr::gamma<>>);
static_assert([] { wr::gamma X({.rate = 3, .shape = 2}); wr::gamma Y({.scale = 3, .shape = 2}); return std::is_same_v<decltype(X), wr::gamma<double>> && std::is_same_v<decltype(Y), wr::gamma<double>>; }());
static_assert([] { wr::beta X({.alpha = 2, .beta = 5}); return std::is_same_v<decltype(X), wr::beta<double>>; }());
static_assert(wr::gamma({.rate = 2.0, .shape = 3.0}).mean() == 1.5 && wr::gamma({.scale = 2.0, .shape = 3.0}).mean() == 6);
static_assert([] { wr::pcg64 g(6); return wr::native(wr::gamma({.rate = 1.0, .shape = 0.5}), g) > 0 && wr::native(wr::beta({.alpha = 2.0, .beta = 0.5}), g) < 1; }());

int main()
{
    wr::xoshiro256pp g(2026);

    // ---- uniform --------------------------------------------------------------
    {
        wr::uniform const X({.lower = 2, .upper = 6});
        assert(X.mean() == 4 && X.median() == 4 && X.variance() == 16.0 / 12 && X.skewness() == 0);
        assert(X.entropy() == wsl::special_functions::detail::pinned_portable::log(4.0));
        assert(wr::pdf(X, 3.0) == 0.25 && wr::pdf(X, 6.0) == 0 && wr::pdf(X, 1.0) == 0 && wr::log_pdf(X, 7.0) == -INFINITY);
        assert(wr::cdf(X, 3.0) == 0.25 && wr::cdf(X, 1.0) == 0 && wr::cdf(X, 9.0) == 1 && wr::ccdf(X, 3.0) == 0.75);
        assert(wr::quantile(X, 0.25) == 3 && wr::quantile(X, 0.0) == 2 && wr::quantile(X, 1.0) == 6);
        assert(wr::cdf(X, 0.1 + 2) + wr::ccdf(X, 0.1 + 2) == 1);
        assert(throws_domain_error([&] { (void)wr::quantile(X, 1.5); }));
        assert(throws_domain_error([] { wr::uniform Y({.lower = 1, .upper = 1}); }));
        assert(throws_domain_error([] { wr::uniform Y({.lower = 0, .upper = INFINITY}); }));

        // Intervalle immense : b − a déborde, rien d'infini ni de NaN là où le résultat est fini.
        double const big = std::numeric_limits<double>::max();
        wr::uniform const H({.lower = -big, .upper = big});
        assert(H.mean() == 0 && std::isinf(H.variance()) && std::abs(H.entropy() - 710.4758600739439) < 1e-12);
        assert(wr::pdf(H, 0.0) > 0 && std::isfinite(wr::pdf(H, 0.0)) && wr::cdf(H, 0.0) == 0.5 && wr::quantile(H, 0.5) == 0);
        assert(wr::quantile(H, 0.75) == big / 2 && wr::ccdf(H, big / 2) == 0.25);

        // Valeurs tirées dans [low, high), loi uniforme (Kolmogorov–Smirnov).
        std::vector<double> values;
        for (int i = 0; i < 20000; ++i)
        {
            double const x = wr::native(X, g);
            assert(x >= 2 && x < 6);
            values.push_back(x);
        }
        assert(kolmogorov_smirnov(X, values) < 1.95);
        wr::uniform const tiny({.lower = 1.0, .upper = std::nextafter(1.0, 2.0)});   // un seul flottant possible
        for (int i = 0; i < 100; ++i)
            assert(wr::native(tiny, g) == 1.0);
        wr::uniform<float> const F({.lower = -1, .upper = 1});
        for (int i = 0; i < 10000; ++i)
        {
            float const x = F(g);
            assert(x >= -1 && x < 1);
        }
        assert(F.mean() == 0 && wr::cdf(F, 0.5f) == 0.75f);
    }

    // ---- discrete_uniform -----------------------------------------------------
    {
        wr::discrete_uniform const D({.lower = 1, .upper = 6});
        assert(D.mean() == 3.5 && D.median() == 3 && D.variance() == 35.0 / 12);
        assert(wr::pmf(D, 3) == 1.0 / 6 && wr::pmf(D, 7) == 0 && wr::cdf(D, 3) == 0.5 && wr::ccdf(D, 3) == 0.5 && wr::cdf(D, 6) == 1);
        for (int k = 1; k <= 6; ++k)
        {
            // Quantile exact ; la cdf arrondie le décale au plus d'un cran.
            long long const q = wr::quantile(D, wr::cdf(D, k));
            assert(q == k || q == k + 1);
        }
        assert(wr::quantile(D, 1e-300) == 1 && wr::quantile(D, 1.0) == 6 && wr::quantile(D, 0.5) == 3 && wr::quantile(D, 0.25) == 2);
        assert(wr::quantile(D, 0.5 + 0x1p-52) == 4);                            // juste au-dessus de 3/6
        std::array<int, 7> counts{};
        for (int i = 0; i < 60000; ++i)
            ++counts[static_cast<std::size_t>(D(g))];
        for (int k = 1; k <= 6; ++k)
            assert(std::abs(counts[static_cast<std::size_t>(k)] - 10000) < 450);

        // Pleine plage de 64 bits : n = 2^64.
        wr::discrete_uniform const W({.lower = std::numeric_limits<std::int64_t>::min(), .upper = std::numeric_limits<std::int64_t>::max()});
        assert(wr::pmf(W, 0) == 0x1p-64 && wr::cdf(W, -1) == 0.5 && W.variance() == 0x1p128 / 12 && W.entropy() == 64 * std::log(2.0));
        assert(wr::quantile(W, 0.5) == -1 && wr::quantile(W, 0x1p-64) == std::numeric_limits<std::int64_t>::min());
        wr::discrete_uniform<std::uint8_t> const B({.lower = 0, .upper = 255});
        assert(wr::quantile(B, 0.5) == 127 && wr::pmf(B, 255) == 1.0 / 256);
        assert(throws_domain_error([] { wr::discrete_uniform E({.lower = 2, .upper = 1}); }));
    }

    // ---- bernoulli ------------------------------------------------------------
    {
        wr::bernoulli const B({.p = 0.3});
        assert(B.mean() == 0.3 && std::abs(B.variance() - 0.21) < 1e-16 && wr::pmf(B, true) == 0.3 && wr::cdf(B, false) == 0.7);
        assert(std::abs(B.entropy() - 0.6108643020548935) < 1e-15);
        int hits = 0;
        for (int i = 0; i < 100000; ++i)
            hits += wr::native(B, g);
        assert(std::abs(hits / 100000.0 - 0.3) < 0.005);
        assert((wr::bernoulli({.p = 0}).support() == wr::interval<bool>{false, false}));
        assert(wr::bernoulli({.p = 1e-300}).variance() == 1e-300 && wr::bernoulli({.p = 0}).entropy() == 0);
        assert(throws_domain_error([] { wr::bernoulli E({.p = 1.5}); }));
    }

    // ---- exponential ----------------------------------------------------------
    {
        wr::exponential const E({.rate = 2});
        assert(E.mean() == 0.5 && E.variance() == 0.25 && E.mode() == 0 && E.skewness() == 2 && E.excess_kurtosis() == 6);
        assert(E.median() == wr::quantile(E, 0.5) && std::abs(E.median() - std::log(2.0) / 2) < 1e-16);
        assert(std::abs(E.entropy() - (1 - std::log(2.0))) < 1e-16);
        assert(wr::pdf(E, 0.0) == 2 && wr::pdf(E, -1.0) == 0 && wr::pdf(E, INFINITY) == 0 && wr::log_pdf(E, -1.0) == -INFINITY);
        assert(wr::cdf(E, 0.0) == 0 && wr::ccdf(E, 0.0) == 1 && wr::cdf(E, INFINITY) == 1 && wr::ccdf(E, INFINITY) == 0);
        assert(wr::quantile(E, 0.0) == 0 && wr::quantile(E, 1.0) == INFINITY && std::isnan(wr::cdf(E, NAN)));
        assert(wr::cdf(E, 1e-300) == 2e-300 && wr::ccdf(E, 360.0) > 0 && wr::ccdf(E, 360.0) < 1e-300);
        assert(throws_domain_error([] { wr::exponential X({.rate = 0}); }));
        assert(throws_domain_error([] { wr::exponential X({.rate = INFINITY}); }));
        assert(throws_domain_error([&] { (void)wr::quantile(E, -0.1); }));
        // Paramètres extrêmes : rien de NaN, densité finie.
        wr::exponential const big({.rate = 0x1p1000}), small({.rate = 0x1p-1000});
        assert(wr::pdf(big, 0x1p-1000) > 0 && std::isfinite(wr::pdf(big, 0x1p-1000)) && wr::cdf(small, 0x1p1000) > 0.63);

        std::vector<double> exact, fast;
        for (int i = 0; i < 20000; ++i)
        {
            exact.push_back(wr::native(E, g));
            fast.push_back(wr::sample(E, g, wsl::approximate));
        }
        assert(std::ranges::all_of(exact, [](double x) { return x >= 0; }));
        assert(kolmogorov_smirnov(E, exact) < 1.95 && kolmogorov_smirnov(E, fast) < 1.95);
        wr::exponential<float> const F({.rate = 0.5f});
        std::vector<double> floats;
        for (int i = 0; i < 20000; ++i)
            floats.push_back(F(g));
        assert(kolmogorov_smirnov(wr::exponential({.rate = 0.5}), floats) < 1.95);
    }

    // ---- normal ---------------------------------------------------------------
    {
        wr::normal const N({.mean = 1, .stddev = 2});
        assert(N.mean() == 1 && N.median() == 1 && N.mode() == 1 && N.variance() == 4 && N.skewness() == 0 && N.excess_kurtosis() == 0);
        assert(std::abs(N.entropy() - (std::log(2.0) + 0.5 * (1 + std::log(2 * 3.141592653589793)))) < 1e-15);
        assert(std::abs(wr::pdf(N, 1.0) - 0.19947114020071635) < 1e-17 && wr::cdf(N, 1.0) == 0.5 && wr::quantile(N, 0.5) == 1);
        for (double d : {0.375, 1.0, 3.0, 10.0, 40.0})                 // 1 ± d exacts
            assert(wr::cdf(N, 1 - d) == wr::ccdf(N, 1 + d));                        // symétrie, exacte
        assert(std::abs(wr::cdf(N, 3.0) - 0.8413447460685429) < 1e-16 && std::abs(wr::quantile(N, 0.8413447460685429) - 3) < 1e-14);
        assert(wr::quantile(N, 0.0) == -INFINITY && wr::quantile(N, 1.0) == INFINITY && wr::cdf(N, -INFINITY) == 0 && wr::pdf(N, INFINITY) == 0);
        assert(std::isnan(wr::pdf(N, NAN)) && std::isnan(wr::cdf(N, NAN)) && wr::log_pdf(N, INFINITY) == -INFINITY);
        assert(wr::cdf(N, -74.0) > 0 && wr::cdf(N, -74.0) < 1e-307 && wr::ccdf(N, 1e6) == 0);
        assert(throws_domain_error([] { wr::normal X({.mean = 0, .stddev = 0}); }));
        assert(throws_domain_error([] { wr::normal X({.mean = INFINITY, .stddev = 1}); }));
        // Paramètres extrêmes.
        double const big = std::numeric_limits<double>::max();
        wr::normal const H({.mean = 0, .stddev = big}), T({.mean = 0, .stddev = 1e-300});
        assert(std::abs(wr::cdf(H, big) - 0.8413447460685429) < 1e-16 && wr::quantile(H, 0.975) == INFINITY && std::isfinite(wr::quantile(H, 0.75)));
        assert(std::isfinite(wr::pdf(T, 0.0)) && wr::pdf(T, 0.0) > 1e299 && wr::cdf(T, 1e-300) == wr::cdf(N, 3.0));
        assert(wr::quantile(wr::normal({.mean = big, .stddev = big}), 0.01) < big);   // μ + σ z sur les moitiés

        std::vector<double> exact, fast;
        int tail = 0;
        for (int i = 0; i < 40000; ++i)
        {
            exact.push_back(wr::native(N, g));
            fast.push_back(wr::sample(N, g, wsl::approximate));
            tail += std::abs(fast.back() - 1) > 2 * 3.6541528853610088;     // queue de la ziggourat
        }
        assert(kolmogorov_smirnov(N, exact) < 1.95 && kolmogorov_smirnov(N, fast) < 1.95);
        assert(tail > 0 && tail < 40);                                      // P(|Z| > R) = 2,6·10⁻⁴
        wr::normal<float> const F({.mean = -3, .stddev = 0.25f});
        std::vector<double> floats;
        for (int i = 0; i < 20000; ++i)
            floats.push_back(F(g));
        assert(kolmogorov_smirnov(wr::normal({.mean = -3.0, .stddev = 0.25}), floats) < 1.95);
        // Tirage exact par la voie lente (σ hors de la voie rapide) : même loi.
        std::vector<double> slow;
        for (int i = 0; i < 300; ++i)
            slow.push_back(T(g));
        assert(kolmogorov_smirnov(T, slow) < 1.95);
    }

    // ---- geometric ------------------------------------------------------------
    {
        wr::geometric const G({.p = 0.5});
        assert(G.mean() == 1 && G.variance() == 2 && G.mode() == 0 && G.median() == 0 && G.entropy() == 2 * std::log(2.0));
        assert(wr::pmf(G, 0) == 0.5 && wr::pmf(G, -1) == 0 && wr::cdf(G, 1) == 0.75 && wr::ccdf(G, 1) == 0.25 && wr::cdf(G, -1) == 0);
        assert(wr::quantile(G, 0.75) == 1 && wr::quantile(G, std::nextafter(0.75, 1.0)) == 2 && wr::quantile(G, 0.0) == 0);
        assert(wr::quantile(G, 1.0) == std::numeric_limits<std::int64_t>::max() && wr::log_pmf(G, 3) == -4 * std::log(2.0));
        wr::geometric const one({.p = 1});
        assert(wr::pmf(one, 0) == 1 && wr::pmf(one, 1) == 0 && wr::cdf(one, 0) == 1 && one.mean() == 0 && wr::native(one, g) == 0);
        wr::geometric const tiny({.p = 0x1p-60});
        assert(tiny.mean() == 0x1p60 - 1 && wr::cdf(tiny, std::numeric_limits<std::int64_t>::max()) > 0.99);
        assert(throws_domain_error([] { wr::geometric X({.p = 0}); }) && throws_domain_error([] { wr::geometric X({.p = 1.5}); }));
        assert(std::abs(chi_square(wr::geometric({.p = 0.3}), g, 200000, 0, 30)) < 4);
        assert(std::abs(chi_square(wr::geometric({.p = 0.001}), g, 200000, 0, 5000)) < 4);
    }

    // ---- poisson --------------------------------------------------------------
    {
        wr::poisson const P({.mean = 3.5});
        assert(P.mean() == 3.5 && P.variance() == 3.5 && P.mode() == 3 && wr::pmf(P, -1) == 0 && wr::log_pmf(P, 0) == -3.5);
        assert(wr::pmf(P, 0) == wsl::special_functions::detail::pinned_portable::exp(-3.5));
        double total = 0;
        for (int k = 0; k < 60; ++k)
            total += wr::pmf(P, k);
        assert(std::abs(total - 1) < 1e-15);
        assert(wr::pmf(wr::poisson({.mean = 1e300}), 5) == 0 && wr::pmf(wr::poisson({.mean = 1e-300}), 1) == 1e-300);
        assert(throws_domain_error([] { wr::poisson X({.mean = 0}); }) && throws_domain_error([] { wr::poisson X({.mean = INFINITY}); }));
        for (double mu : {0.7, 3.5, 9.9, 10.0, 50.0, 1e4, 1e9})              // inversion, puis PTRS
        {
            double const sd = std::sqrt(mu);
            auto const low = static_cast<long long>(std::max(0.0, mu - 6 * sd)), high = static_cast<long long>(mu + 6 * sd + 6);
            assert(std::abs(chi_square(wr::poisson({.mean = mu}), g, 200000, low, high)) < 4);
        }
    }

    // ---- binomial -------------------------------------------------------------
    {
        wr::binomial const B({.trials = 20, .p = 0.3});
        assert(B.trials() == 20 && B.mean() == 6 && std::abs(B.variance() - 4.2) < 1e-15 && B.mode() == 6);
        assert(wr::pmf(B, 21) == 0 && wr::pmf(B, -1) == 0 && wr::pmf(B, 20) == wsl::special_functions::detail::pinned_portable::pow(0.3, 20.0));
        double total = 0;
        for (int k = 0; k <= 20; ++k)
            total += wr::pmf(B, k);
        assert(std::abs(total - 1) < 1e-15);
        assert(wr::pmf(wr::binomial({.trials = 5, .p = 0.0}), 0) == 1 && wr::pmf(wr::binomial({.trials = 5, .p = 1.0}), 5) == 1);
        assert(wr::native(wr::binomial({.trials = 5, .p = 1.0}), g) == 5 && wr::native(wr::binomial({.trials = 0, .p = 0.5}), g) == 0);
        assert(throws_domain_error([] { wr::binomial X({.trials = -1, .p = 0.5}); }) && throws_domain_error([] { wr::binomial X({.trials = 3, .p = 2}); }));
        for (auto [n, p] : {std::pair{20LL, 0.3}, {1000LL, 0.004}, {1000LL, 0.3}, {1000LL, 0.9}, {100000000LL, 0.5}, {50LL, 0.97}})
        {
            double const mean = double(n) * p, sd = std::sqrt(mean * (1 - p));
            auto const low = static_cast<long long>(std::max(0.0, mean - 6 * sd)), high = std::min(n, static_cast<long long>(mean + 6 * sd + 6));
            assert(std::abs(chi_square(wr::binomial({.trials = n, .p = p}), g, 200000, low, high)) < 4);
        }
    }

    // ---- gamma ----------------------------------------------------------------
    {
        wr::gamma const X({.rate = 2, .shape = 2.5}), Y({.scale = 0.5, .shape = 2.5});
        assert(X.mean() == 1.25 && X.variance() == 0.625 && X.mode() == 0.75 && X.shape() == 2.5 && X.rate() == 2 && Y.scale() == 0.5);
        assert(Y.mean() == 1.25 && Y.variance() == 0.625 && wr::pdf(X, 1.0) == wr::pdf(Y, 1.0) && wr::log_pdf(X, 0.3) == wr::log_pdf(Y, 0.3));
        assert(wr::pdf(X, -1.0) == 0 && wr::pdf(X, 0.0) == 0 && wr::pdf(X, INFINITY) == 0 && std::isnan(wr::pdf(X, NAN)) && wr::log_pdf(X, -1.0) == -INFINITY);
        assert(wr::pdf(wr::gamma({.rate = 3.0, .shape = 1.0}), 0.0) == 3 && wr::pdf(wr::gamma({.rate = 3.0, .shape = 0.5}), 0.0) == INFINITY);
        assert(wr::pdf(wr::gamma({.rate = 3.0, .shape = 1.0}), 0.5) == wr::pdf(wr::exponential({.rate = 3.0}), 0.5));
        assert(throws_domain_error([] { wr::gamma Z({.rate = 1, .shape = 0}); }) && throws_domain_error([] { wr::gamma Z({.scale = -1, .shape = 1}); }));
        assert(std::isfinite(wr::pdf(wr::gamma({.rate = 1.0, .shape = 1e-300}), 1e-300)) && wr::pdf(wr::gamma({.rate = 1.0, .shape = 1e300}), 1e300) > 0);
        // Répartitions en forme close : forme entière (Erlang), χ²₁ = gamma(1/2, taux 1/2).
        auto const erlang3 = [](double x) { double const t = 2 * x; return 1 - std::exp(-t) * (1 + t + t * t / 2); };
        assert(kolmogorov_smirnov_against(wr::gamma({.rate = 2.0, .shape = 3.0}), g, erlang3) < 1.95);
        assert(kolmogorov_smirnov_against(wr::gamma({.scale = 0.5, .shape = 3.0}), g, erlang3) < 1.95);
        assert(kolmogorov_smirnov_against(wr::gamma({.rate = 0.5, .shape = 0.5}), g, [](double x) { return std::erf(std::sqrt(x / 2)); }) < 1.95);
        assert(kolmogorov_smirnov_against(wr::gamma<float>({.rate = 4, .shape = 1}), g, [](double x) { return 1 - std::exp(-4 * x); }) < 1.95);
        double s = 0;
        wr::gamma const big({.rate = 0.5, .shape = 400.0});
        for (int i = 0; i < 20000; ++i)
            s += wr::native(big, g);
        assert(std::abs(s / 20000 - 800) < 4 * std::sqrt(1600.0 / 20000));
    }

    // ---- beta -----------------------------------------------------------------
    {
        wr::beta const B({.alpha = 2, .beta = 5});
        assert(B.a() == 2 && B.b() == 5 && std::abs(B.mean() - 2.0 / 7) < 1e-16 && B.mode() == 0.2 && std::abs(B.variance() - 10.0 / 392) < 1e-17);
        assert(wr::pdf(B, 0.0) == 0 && wr::pdf(B, 1.0) == 0 && wr::pdf(B, -0.5) == 0 && std::abs(wr::pdf(B, 0.3) - 2.1609) < 1e-15);
        assert(wr::pdf(wr::beta({.alpha = 1.0, .beta = 3.0}), 0.0) == 3 && wr::pdf(wr::beta({.alpha = 0.5, .beta = 3.0}), 0.0) == INFINITY);
        assert(wr::pdf(wr::beta({.alpha = 1.0, .beta = 1.0}), 0.7) == 1 && std::isnan(wr::beta({.alpha = 0.5, .beta = 0.5}).mode()));
        assert(throws_domain_error([] { wr::beta Z({.alpha = 0, .beta = 1}); }));
        double const pi = 3.141592653589793;
        assert(kolmogorov_smirnov_against(wr::beta({.alpha = 0.5, .beta = 0.5}), g, [&](double x) { return 2 / pi * std::asin(std::sqrt(x)); }) < 1.95);
        assert(kolmogorov_smirnov_against(wr::beta({.alpha = 3.0, .beta = 1.0}), g, [](double x) { return x * x * x; }) < 1.95);
        assert(kolmogorov_smirnov_against(wr::beta({.alpha = 1.0, .beta = 0.2}), g, [](double x) { return 1 - std::pow(1 - x, 0.2); }) < 1.95);
        assert(kolmogorov_smirnov_against(wr::beta({.alpha = 0.05, .beta = 1.0}), g, [](double x) { return std::pow(x, 0.05); }) < 1.95);
        assert(kolmogorov_smirnov_against(wr::beta<float>({.alpha = 2, .beta = 1}), g, [](double x) { return x * x; }) < 1.95);
    }

    // ---- Répartitions et quantiles (phase 4) : cohérences ---------------------------
    {
        wr::poisson const P({.mean = 3.5});
        wr::binomial const B({.trials = 20, .p = 0.3});
        double sp = 0, sb = 0;
        for (int k = 0; k <= 20; ++k)
        {
            sp += wr::pmf(P, k);
            sb += wr::pmf(B, k);
            assert(std::abs(wr::cdf(P, k) - sp) < 1e-15 && std::abs(wr::cdf(B, k) - sb) < 1e-15);
            assert(std::abs(wr::cdf(P, k) + wr::ccdf(P, k) - 1) < 1e-15 && std::abs(wr::cdf(B, k) + wr::ccdf(B, k) - 1) < 1e-15);
        }
        assert(wr::cdf(P, -1) == 0 && wr::ccdf(P, -1) == 1 && wr::cdf(B, 20) == 1 && wr::ccdf(B, 25) == 0);
        for (double u : {1e-300, 1e-9, 0.1, 0.5, 0.75, 0.999999})
        {
            long long const k = wr::quantile(P, u), j = wr::quantile(B, u);          // quantile ∘ cdf : plus petit k
            assert(wr::cdf(P, k) >= u && (k == 0 || wr::cdf(P, k - 1) < u) && wr::cdf(B, j) >= u && (j == 0 || wr::cdf(B, j - 1) < u));
        }
        assert(wr::quantile(P, 0.0) == 0 && wr::quantile(B, 1.0) == 20 && wr::quantile(wr::binomial({.trials = 5, .p = 1.0}), 0.3) == 5);
        wr::gamma const G({.rate = 3, .shape = 2.5});
        wr::beta const Be({.alpha = 2, .beta = 5});
        for (double u : {1e-200, 1e-6, 0.3, 0.5, 0.9, 1 - 1e-10})
        {
            double const x = wr::quantile(G, u), y = wr::quantile(Be, u);
            assert(std::abs(wr::cdf(G, x) / u - 1) < 1e-14 && std::abs(wr::cdf(Be, y) / u - 1) < 1e-14);
        }
        assert(wr::cdf(G, 0.0) == 0 && wr::ccdf(G, INFINITY) == 0 && wr::quantile(G, 1.0) == INFINITY && wr::cdf(Be, 1.0) == 1 && wr::quantile(Be, 0.0) == 0);
        assert(std::abs(wr::gamma({.rate = 2.0, .shape = 1.0}).entropy() - wr::exponential({.rate = 2.0}).entropy()) < 1e-16);
        assert(std::abs(wr::beta({.alpha = 1.0, .beta = 1.0}).entropy()) < 1e-20);
        assert(std::abs(wr::cdf(wr::gamma({.rate = 2.0, .shape = 1.0}), 0.7) - wr::cdf(wr::exponential({.rate = 2.0}), 0.7)) < 1e-16);
    }

    // ---- Lois dérivées -----------------------------------------------------------
    {
        wr::chi_squared const C({.degrees_of_freedom = 4});
        assert(C.mean() == 4 && C.variance() == 8 && C.mode() == 2 && std::abs(wr::cdf(C, 2.0) - (1 - 2 * std::exp(-1.0))) < 1e-16);
        assert(kolmogorov_smirnov(C, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wr::native(C, g)); return v; }()) < 1.95);
        wr::student_t const T1({.degrees_of_freedom = 1}), T5({.degrees_of_freedom = 5});
        assert(std::isnan(T1.mean()) && T5.mean() == 0 && std::abs(T5.variance() - 5.0 / 3) < 1e-16 && std::isinf(wr::student_t({.degrees_of_freedom = 1.5}).variance()));
        // t à 1 degré de liberté = Cauchy standard.
        wr::cauchy const K({.location = 0, .scale = 1});
        for (double t : {-1e6, -3.0, -0.25, 0.5, 2.0, 1e12})
        {
            assert(std::abs(wr::cdf(T1, t) - wr::cdf(K, t)) <= 2e-16 * wr::cdf(K, t) && std::abs(wr::pdf(T1, t) / wr::pdf(K, t) - 1) < 3e-16);
            assert(std::abs(wr::ccdf(T5, t) - wr::cdf(T5, -t)) <= 2e-16 * wr::ccdf(T5, t));                   // symétrie
        }
        assert(wr::quantile(T5, 0.25) == -wr::quantile(T5, 0.75));                  // 1 − 0,25 exact
        for (double u : {1e-300, 1e-8, 0.2, 0.5, 0.975})
        {
            double const t = wr::quantile(T5, u);
            assert(u == 0.5 ? t == 0 : std::abs(wr::cdf(T5, t) / u - 1) < 1e-14);
        }
        assert(kolmogorov_smirnov(T5, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wr::native(T5, g)); return v; }()) < 1.95);
        assert(kolmogorov_smirnov(wr::student_t({.degrees_of_freedom = 0.7}), [&] { std::vector<double> v; wr::student_t const T({.degrees_of_freedom = 0.7}); for (int i = 0; i < 20000; ++i) v.push_back(T(g)); return v; }()) < 1.95);
        wr::lognormal const L({.log_mean = 0.3, .log_stddev = 0.5});
        assert(std::abs(L.median() - std::exp(0.3)) < 3e-16 && std::abs(L.mean() - std::exp(0.3 + 0.125)) < 1e-15 && std::abs(wr::cdf(L, L.median()) - 0.5) < 2e-16);
        assert(wr::cdf(L, 0.0) == 0 && wr::pdf(L, -1.0) == 0 && wr::quantile(L, 0.0) == 0 && wr::quantile(L, 1.0) == INFINITY);
        assert(kolmogorov_smirnov(L, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wr::native(L, g)); return v; }()) < 1.95);
        assert(kolmogorov_smirnov(L, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wr::sample(L, g, wsl::approximate)); return v; }()) < 1.95);
        wr::cauchy const C2({.location = 1, .scale = 2});
        assert(wr::cdf(C2, 3.0) == 0.75 && wr::quantile(C2, 0.75) == 3 && C2.median() == 1 && std::abs(C2.entropy() - std::log(8 * 3.141592653589793)) < 1e-15);
        assert(kolmogorov_smirnov(C2, [&] { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(wr::native(C2, g)); return v; }()) < 1.95);
        assert(throws_domain_error([] { wr::student_t X({.degrees_of_freedom = 0}); }) && throws_domain_error([] { wr::cauchy X({.location = 0, .scale = -1}); }));
        assert(throws_domain_error([] { wr::chi_squared X({.degrees_of_freedom = -2}); }) && throws_domain_error([] { wr::lognormal X({.log_mean = 0, .log_stddev = 0}); }));
    }

    // ---- Lois dérivées, second lot --------------------------------------------------
    {
        auto draws = [&](auto const& law) { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(double(wr::native(law, g))); return v; };
        wr::fisher_f const F({.d1 = 5, .d2 = 10});
        assert(F.mean() == 1.25 && wr::fisher_f({.d1 = 3, .d2 = 2}).mean() == INFINITY && std::isinf(wr::fisher_f({.d1 = 3, .d2 = 3}).variance()));
        assert(wr::pdf(wr::fisher_f({.d1 = 2, .d2 = 5}), 0.0) == 1 && wr::pdf(F, 0.0) == 0 && wr::pdf(wr::fisher_f({.d1 = 1, .d2 = 5}), 0.0) == INFINITY);
        // F(1, ν) = t(ν)² : P(F <= t²) = P(|T| <= t).
        wr::student_t const T7({.degrees_of_freedom = 7});
        wr::fisher_f const F17({.d1 = 1, .d2 = 7});
        for (double t : {0.5, 1.0, 2.5, 9.0})
            assert(std::abs(wr::cdf(F17, t * t) - (1 - 2 * wr::ccdf(T7, t))) < 1e-15);
        assert(kolmogorov_smirnov(F, draws(F)) < 1.95 && kolmogorov_smirnov(wr::fisher_f({.d1 = 0.6, .d2 = 1.5}), draws(wr::fisher_f({.d1 = 0.6, .d2 = 1.5}))) < 1.95);
        wr::weibull const W1({.scale = 2, .shape = 1});                    // = exponentielle de taux 1/2
        wr::exponential const E2({.rate = 0.5});
        for (double x : {0.1, 1.0, 5.0, 30.0})
            assert(std::abs(wr::cdf(W1, x) - wr::cdf(E2, x)) <= 2e-16 * wr::cdf(E2, x) && std::abs(wr::pdf(W1, x) / wr::pdf(E2, x) - 1) < 3e-16);
        wr::weibull const W({.scale = 2, .shape = 1.5});
        assert(std::abs(W.median() - wr::quantile(W, 0.5)) == 0 && kolmogorov_smirnov(W, draws(W)) < 1.95);
        assert(wr::weibull({.scale = 1, .shape = 1e4}).variance() > 0);     // compensation traitée
        wr::laplace const L({.location = 1, .scale = 2});
        assert(L.variance() == 8 && wr::cdf(L, 1.0) == 0.5 && wr::quantile(L, 0.5) == 1 && kolmogorov_smirnov(L, draws(L)) < 1.95);
        wr::logistic const Lg({.location = 1, .scale = 2});
        assert(wr::cdf(Lg, 1.0) == 0.5 && std::abs(wr::quantile(Lg, 0.75) - (1 + 2 * std::log(3.0))) < 1e-15 && kolmogorov_smirnov(Lg, draws(Lg)) < 1.95);
        wr::pareto const P({.scale = 1, .shape = 2.5});
        assert(wr::cdf(P, 1.0) == 0 && wr::pdf(P, 0.5) == 0 && wr::quantile(P, 0.0) == 1 && std::isinf(wr::pareto({.scale = 1, .shape = 1}).mean()));
        assert(kolmogorov_smirnov(P, draws(P)) < 1.95);
        assert(throws_domain_error([] { wr::fisher_f X({.d1 = 0, .d2 = 1}); }) && throws_domain_error([] { wr::weibull X({.scale = 0, .shape = 1}); }));
        assert(throws_domain_error([] { wr::laplace X({.location = 0, .scale = -1}); }) && throws_domain_error([] { wr::pareto X({.scale = 0, .shape = 1}); }));
    }

    // ---- Compléments : Gumbel, Fréchet, Rayleigh, triangulaire, arcsinus, Erlang ----
    {
        auto draws = [&](auto const& law) { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(double(wr::native(law, g))); return v; };
        wr::gumbel const Gu({.location = 1, .scale = 2});
        wr::frechet const Fr({.location = 0.5, .scale = 2, .shape = 3});
        wr::rayleigh const Ra({.scale = 2});
        wr::triangular const Tr({.lower = 0, .mode = 1, .upper = 3});
        wr::arcsine const Ar({.lower = -1, .upper = 2});
        wr::erlang const Er({.rate = 1.5, .shape = 4});
        assert(kolmogorov_smirnov(Gu, draws(Gu)) < 1.95 && kolmogorov_smirnov(Fr, draws(Fr)) < 1.95 && kolmogorov_smirnov(Ra, draws(Ra)) < 1.95);
        assert(kolmogorov_smirnov(Tr, draws(Tr)) < 1.95 && kolmogorov_smirnov(Ar, draws(Ar)) < 1.95 && kolmogorov_smirnov(Er, draws(Er)) < 1.95);
        static_assert(wr::erlang<>::sampling == wr::exactness::exact && wr::has_cquantile<wr::frechet<>> && wr::has_entropy<wr::arcsine<float>>);
        wr::arcsine const A01({.lower = 0, .upper = 1});                     // = bêta(½, ½)
        wr::beta const B({.alpha = 0.5, .beta = 0.5});
        for (double x : {1e-9, 0.1, 0.5, 0.73, 0.999})
            assert(std::abs(wr::cdf(A01, x) - wr::cdf(B, x)) <= 2e-16 * wr::cdf(B, x) && std::abs(wr::pdf(A01, x) / wr::pdf(B, x) - 1) < 1e-15);
        wr::erlang const E1({.rate = 2, .shape = 1});                      // = exponentielle
        for (double x : {0.01, 1.0, 7.0})
            assert(wr::cdf(E1, x) == wr::cdf(wr::gamma({.rate = 2.0, .shape = 1.0}), x) && std::abs(wr::cdf(E1, x) / wr::cdf(wr::exponential({.rate = 2.0}), x) - 1) < 3e-16);
        std::vector<double> y;                                              // α log((X − m)/s) : Gumbel(0, 1)
        for (double v : draws(Fr))
            y.push_back(3 * std::log((v - 0.5) / 2));
        assert(kolmogorov_smirnov(wr::gumbel({.location = 0.0, .scale = 1.0}), y) < 1.95);
        assert(Tr.median() == wr::quantile(Tr, 0.5) && Ar.median() == 0.5 && wr::quantile(Ra, 0.5) == Ra.median() && wr::cquantile(Gu, 0.5) == wr::quantile(Gu, 0.5));
        assert(wr::pdf(Fr, 0.5) == 0 && wr::cdf(Fr, 0.4) == 0 && wr::pdf(Tr, 3.5) == 0 && wr::pdf(Ar, -1.0) == INFINITY && wr::quantile(Tr, 1.0) == 3);
        assert(throws_domain_error([] { wr::triangular X({.lower = 0, .mode = 4, .upper = 3}); }) && throws_domain_error([] { wr::erlang X({.rate = 1.0, .shape = 0}); }));
        assert(throws_domain_error([] { wr::frechet X({.scale = 1, .shape = 0}); }) && throws_domain_error([] { wr::arcsine X({.lower = 1, .upper = 1}); }));
    }

    // ---- Inverse-gaussienne ----
    {
        for (auto [mu, lambda] : {std::pair{1.0, 3.0}, std::pair{2.0, 0.05}, std::pair{1e-2, 1e3}})
        {
            wr::inverse_gaussian const W({.mean = mu, .shape = lambda});
            std::vector<double> v;
            for (int i = 0; i < 5000; ++i)
                v.push_back(wr::native(W, g));
            assert(kolmogorov_smirnov(W, v) < 1.95);
            double const m = std::accumulate(v.begin(), v.end(), 0.0) / double(v.size());
            assert(std::abs(m - mu) < 8 * std::sqrt(W.variance() / double(v.size())));
            for (double p : {1e-200, 1e-8, 0.25, 0.5, 0.75, 1 - 1e-9})   // quantiles et répartition se répondent
            {
                double const x = wr::quantile(W, p), y = wr::cquantile(W, p);
                assert(std::abs(wr::cdf(W, x) / p - 1) < 1e-14 || std::abs(wr::ccdf(W, x) / (1 - p) - 1) < 1e-14);
                assert(std::abs(wr::ccdf(W, y) / p - 1) < 1e-14 || std::abs(wr::cdf(W, y) / (1 - p) - 1) < 1e-14);
            }
            assert(W.median() == wr::quantile(W, 0.5) && wr::cdf(W, 0.0) == 0 && wr::ccdf(W, INFINITY) == 0 && wr::pdf(W, -1.0) == 0);
            double area = 0;                                                // Simpson sur [μ, 1,5 μ]
            for (int i = 0; i < 2000; ++i)
            {
                double const x0 = mu * (1 + i / 4000.0), x1 = mu * (1 + (i + 1) / 4000.0);
                area += (x1 - x0) / 6 * (wr::pdf(W, x0) + 4 * wr::pdf(W, (x0 + x1) / 2) + wr::pdf(W, x1));
            }
            assert(std::abs((wr::cdf(W, 1.5 * mu) - wr::cdf(W, mu)) / area - 1) < 1e-9);
        }
        wr::inverse_gaussian const W({.mean = 1.0, .shape = 0.5});        // queue droite : wr::ccdf(x) ~ e^{−λx/(2μ²)}
        assert(wr::ccdf(W, 1e3) > 0 && wr::ccdf(W, 1e3) < 1e-100 && wr::ccdf(W, 3e3) == 0);
        assert(throws_domain_error([] { wr::inverse_gaussian X({.mean = 0.0, .shape = 1.0}); }) && throws_domain_error([] { wr::inverse_gaussian X({.mean = 1.0, .shape = INFINITY}); }));
    }

    // ---- Répartition de von Mises ----
    {
        for (double kappa : {0.0, 0.3, 4.0, 200.0, 1e6})
        {
            wr::von_mises const V({.location = 1.0, .concentration = kappa});
            std::vector<double> v;
            for (int i = 0; i < 2000; ++i)
                v.push_back(V(g));
            assert(kolmogorov_smirnov(V, v) < 1.95);
            assert(wr::cdf(V, 1.0) == 0.5 && wr::ccdf(V, 1.0) == 0.5 && wr::quantile(V, 0.5) == 1.0);
            assert(wr::cdf(V, 1.0 - 3.2) == 0 && wr::cdf(V, 1.0 + 3.2) == 1 && wr::quantile(V, 0.0) == 1.0 - 3.141592653589793);
            for (double p : {1e-300, 1e-9, 0.1, 0.4, 0.75, 1 - 1e-12})
            {
                double const x = wr::quantile(V, p);                            // F(x⁻) <= p <= F(x⁺), à un ulp près
                double const below = std::nextafter(std::nextafter(x, -INFINITY), -INFINITY), above = std::nextafter(std::nextafter(x, INFINITY), INFINITY);
                if (p <= 0.5)
                    assert(wr::cdf(V, below) <= p && wr::cdf(V, above) >= p);
                else
                    assert(wr::ccdf(V, below) >= 1 - p && wr::ccdf(V, above) <= 1 - p);
                assert(p < 0.5 || wr::cquantile(V, 1 - p) == x);
            }
            for (double x : {-1.5, 0.2, 0.99, 1.7, 3.9})                    // F(x) + S(x) = 1, et symétrie autour de μ
                assert(std::abs(wr::cdf(V, x) + wr::ccdf(V, x) - 1) < 1e-15 && wr::cdf(V, 1.0 - (x - 1.0)) == wr::ccdf(V, x));
        }
        wr::von_mises const U({.location = 0.0, .concentration = 0.0});                 // κ = 0 : uniforme sur [−π, π]
        for (double x : {-3.0, -1.0, 0.5, 2.0})
            assert(std::abs(wr::cdf(U, x) - (x + 3.141592653589793) / (2 * 3.141592653589793)) < 1e-16);
    }

    // ---- Binomiale négative, hypergéométrique ----
    {
        wr::negative_binomial const B({.successes = 3.5, .p = 0.4});
        assert(std::abs(chi_square(B, g, 50000, 0, 30)) < 4);
        assert(std::abs(chi_square(wr::negative_binomial({.successes = 0.3, .p = 0.05}), g, 50000, 0, 60)) < 4);
        wr::negative_binomial const G1({.successes = 1.0, .p = 0.3});     // r = 1 : géométrique
        for (long long k : {0LL, 3LL, 40LL})
            assert(std::abs(wr::pmf(G1, k) / wr::pmf(wr::geometric({.p = 0.3}), k) - 1) < 4e-16 && std::abs(wr::cdf(G1, k) / wr::cdf(wr::geometric({.p = 0.3}), k) - 1) < 4e-16);
        for (double u : {1e-12, 0.3, 0.999999})
        {
            long long const k = wr::quantile(B, u);
            assert(wr::cdf(B, k) >= u && (k == 0 || wr::cdf(B, k - 1) < u));
        }
        wr::hypergeometric const H({.population = 50, .successes = 20, .draws = 10});
        assert(std::abs(chi_square(H, g, 50000, 0, 10)) < 4);
        wr::hypergeometric const H2({.population = 100000, .successes = 30000, .draws = 2000});
        assert(std::abs(chi_square(H2, g, 30000, 520, 680)) < 4);
        double total = 0;
        for (long long k = 0; k <= 10; ++k)
            total += wr::pmf(H, k);
        assert(std::abs(total - 1) < 1e-15 && H.mode() == 4 && H.mean() == 4 && wr::cdf(H, 10) == 1 && wr::ccdf(H, -1) == 1);
        for (double u : {1e-9, 0.5, 0.97})
        {
            long long const k = wr::quantile(H2, u);
            assert(wr::cdf(H2, k) >= u && wr::cdf(H2, k - 1) < u);
        }
        wr::hypergeometric const H3({.population = 30, .successes = 29, .draws = 15});   // support {14, 15}
        assert(H3.lower() == 14 && wr::pmf(H3, 13) == 0 && std::abs(wr::pmf(H3, 14) - 0.5) < 1e-16);
        assert(throws_domain_error([] { wr::hypergeometric X({.population = 10, .successes = 11, .draws = 2}); }) && throws_domain_error([] { wr::negative_binomial X({.successes = 2, .p = 0}); }));
    }

    // ---- Bêta-binomiale, logarithmique, gaussienne discrète ----
    {
        wr::beta_binomial const B({.trials = 20, .alpha = 2, .beta = 3});
        assert(std::abs(chi_square(B, g, 50000, 0, 20)) < 4);
        assert(std::abs(chi_square(wr::beta_binomial({.trials = 30, .alpha = 0.3, .beta = 0.4}), g, 50000, 0, 30)) < 4);
        wr::beta_binomial const U({.trials = 9, .alpha = 1, .beta = 1});   // α = β = 1 : uniforme sur {0 … n}
        for (long long k = 0; k <= 9; ++k)
            assert(std::abs(wr::pmf(U, k) - 0.1) < 2e-17);
        wr::logarithmic const L({.p = 0.9});
        assert(std::abs(chi_square(L, g, 50000, 1, 40)) < 4 && std::abs(chi_square(wr::logarithmic({.p = 0.05}), g, 50000, 1, 4)) < 4);
        for (double u : {1e-6, 0.5, 0.9, 0.999999})
        {
            long long const k = wr::quantile(L, u);
            assert(wr::cdf(L, k) >= u && (k == 1 || wr::cdf(L, k - 1) < u));
        }
        wr::discrete_gaussian const D({.location = 0.3, .scale = 2.5});
        assert(std::abs(chi_square(D, g, 50000, -9, 10)) < 4);
        assert(std::abs(chi_square(wr::discrete_gaussian({.location = -0.5, .scale = 0.3}), g, 50000, -2, 1)) < 4);
        assert(std::abs(chi_square(wr::discrete_gaussian({.location = 7.25, .scale = 40.0}), g, 50000, -110, 125)) < 4);
        wr::discrete_gaussian const E({.location = 2.0, .scale = 3.0});        // μ entier : symétrie, moyenne μ
        assert(E.mean() == 2.0 && wr::cdf(E, 1) == wr::ccdf(E, 2) && wr::pmf(E, 0) == wr::pmf(E, 4));
        for (double u : {1e-12, 0.3, 0.5, 0.97})
        {
            long long const k = wr::quantile(D, u);
            assert(wr::cdf(D, k) >= u && wr::cdf(D, k - 1) < u);
        }
        assert(throws_domain_error([] { wr::logarithmic X({.p = 1.0}); }) && throws_domain_error([] { wr::discrete_gaussian X({.location = 0, .scale = 0}); }));
    }

    // ---- Normale asymétrique ----
    {
        for (double a : {0.0, 1.0, -3.0, 20.0})
        {
            wr::skew_normal const X({.location = 0.5, .scale = 2.0, .shape = a});
            std::vector<double> v;
            for (int i = 0; i < 3000; ++i)
                v.push_back(wr::native(X, g));
            assert(kolmogorov_smirnov(X, v) < 1.95);
            for (double p : {1e-200, 1e-6, 0.3, 0.5, 0.9})
            {
                double const x = wr::quantile(X, p);
                assert(wr::cdf(X, std::nextafter(x, -INFINITY)) <= p * (1 + 1e-15) && wr::cdf(X, std::nextafter(x, INFINITY)) >= p * (1 - 1e-15));
            }
        }
        wr::skew_normal const Z({.shape = 0.0}), One({.shape = 1.0});    // α = 0 : normale ; α = 1 : F = Φ²
        wr::normal const N({.mean = 0.0, .stddev = 1.0});
        for (double x : {-20.0, -3.0, 0.4, 5.0})
        {
            assert(std::abs(wr::cdf(Z, x) / wr::cdf(N, x) - 1) < 3e-16 && std::abs(wr::pdf(Z, x) / wr::pdf(N, x) - 1) < 3e-16);
            assert(std::abs(wr::cdf(One, x) / (wr::cdf(N, x) * wr::cdf(N, x)) - 1) < 1e-15);
        }
        assert(throws_domain_error([] { wr::skew_normal X({.scale = 0.0}); }));
    }

    // ---- Lois stables ----
    {
        auto draws = [&](auto const& law, int n) { std::vector<double> v; for (int i = 0; i < n; ++i) v.push_back(wr::native(law, g)); return v; };
        wr::stable const Gs({.scale = 1.5, .alpha = 2.0, .beta = 0.4});   // normale de variance 2γ²
        wr::stable const Ca({.location = 1.0, .scale = 0.5, .alpha = 1.0, .beta = 0.0});
        wr::stable const Lv({.alpha = 0.5, .beta = 1.0});                  // Lévy : densité √(1/2π) y^{−3/2} e^{−1/(2y)}
        assert(kolmogorov_smirnov(Gs, draws(Gs, 5000)) < 1.95 && kolmogorov_smirnov(Ca, draws(Ca, 5000)) < 1.95 && kolmogorov_smirnov(Lv, draws(Lv, 3000)) < 1.95);
        for (auto [a, b] : {std::pair{1.5, 0.5}, std::pair{0.7, -0.4}, std::pair{1.0, 0.8}})
        {
            wr::stable const X({.location = -1.0, .scale = 2.0, .alpha = a, .beta = b});
            assert(kolmogorov_smirnov(X, draws(X, 200)) < 1.95);                // répartition coûteuse (~1 ms)
            for (double p : {1e-6, 0.8})
            {
                double const x = wr::quantile(X, p);
                assert(wr::cdf(X, std::nextafter(x, -INFINITY)) <= p * (1 + 1e-14) && wr::cdf(X, std::nextafter(x, INFINITY)) >= p * (1 - 1e-14));
            }
        }
        wr::normal const N({.mean = 0.0, .stddev = 1.4142135623730951 * 1.5});
        for (double x : {-3.0, 0.5, 4.0})
            assert(std::abs(wr::cdf(Gs, x) / wr::cdf(N, x) - 1) < 1e-15);
        for (double y : {0.1, 1.0, 30.0})
            assert(std::abs(wr::pdf(Lv, y) / (0.3989422804014327 * std::pow(y, -1.5) * std::exp(-0.5 / y)) - 1) < 1e-14);
        assert(wr::cdf(Lv, -1.0) == 0 && wr::pdf(Lv, 0.0) == 0 && Lv.support().lower == 0);
        assert(throws_domain_error([] { wr::stable X({.alpha = 2.5}); }) && throws_domain_error([] { wr::stable X({.alpha = 1.0, .beta = 1.5}); }));
    }

    // ---- Objets-fonctions, densité, vue des tirages ----
    {
                wr::normal const N({.mean = 1.0, .stddev = 2.0});
        wr::poisson const P({.mean = 3.5});
        static_assert(std::is_same_v<wr::reference_measure_t<wr::normal<>>, wr::lebesgue> && std::is_same_v<wr::reference_measure_t<wr::poisson<>>, wr::counting>);
        static_assert(wr::has_density<wr::normal<>> && wr::has_density<wr::poisson<>> && !wr::has_mean<wr::cauchy<>>);
        assert(wr::pdf(N, 0.5) == wr::pdf(N, 0.5) && wr::density(N, 0.5) == wr::pdf(N, 0.5) && wr::density(P, 2) == wr::pmf(P, 2));
        assert(wr::mean(N) == N.mean() && wr::quantile(N, 0.9) == wr::quantile(N, 0.9) && wr::log_density(P, 4) == wr::log_pmf(P, 4));
        std::vector<double> xs{-1.0, 0.0, 2.5};
        auto const densities = xs | std::views::transform(std::bind_front(wr::pdf, N));   // passable en argument
        assert(std::ranges::equal(densities, std::vector<double>{wr::pdf(N, -1.0), wr::pdf(N, 0.0), wr::pdf(N, 2.5)}));
        wr::xoshiro256pp a(7), b(7);
        std::vector<double> lazy;
        for (double v : wr::samples(N, a) | std::views::take(100))
            lazy.push_back(v);
        for (double v : lazy)
            assert(v == wr::native(N, b));
        assert(a == b);                                                     // même consommation
        // membres et appel X(g), X(n, g) : mêmes valeurs et même consommation que sample
        assert(N.pdf(0.5) == wr::pdf(N, 0.5) && N.density(0.5) == N.pdf(0.5) && P.density(2) == P.pmf(2) && P.cdf(3) == wr::cdf(P, 3));
        wr::xoshiro256pp c(11), d(11), e(11);
        std::vector<double> const batch = N(5, c);
        for (double v : batch)
            assert(v == N(d));
        for (double v : batch)
            assert(v == wr::sample(N, e));
        assert(c == d && d == e && N(c, wsl::approximate) == wr::sample(N, d, wsl::approximate));
        // Le défaut est approché (décision 27) : X(g) == X(g, wsl::approximate) ; wsl::exact
        // donne le tirage natif d'une loi exacte ; l'axe de la reproductibilité est sans effet.
        wr::xoshiro256pp g1(5), g2(5), g3(5), g4(5);
        assert(N(g1) == N(g2, wsl::approximate) && N(g3, wsl::exact) == wr::native(N, g4) && g1 == g2 && g3 == g4);
        assert(N(g1, wsl::portable) == N(g2, wsl::deterministic, wsl::approximate));
    }

    // ---- Tirage rapide demandé : à défaut de méthode propre, le tirage par défaut ----
    {
        wr::xoshiro256pp a(5), b(5);
        wr::uniform const X({.lower = 0, .upper = 1});
        for (int i = 0; i < 100; ++i)
            assert(wr::sample(X, a, wsl::approximate) == wr::native(X, b));
    }

    // ---- Adaptateur vers <random> -------------------------------------------------
    {
        auto d = wr::as_distribution(wr::discrete_uniform({.lower = 1, .upper = 6}));
        std::vector<long long> rolls(1000);
        std::ranges::generate(rolls, [&] { return d(g); });
        assert(std::ranges::all_of(rolls, [](long long x) { return x >= 1 && x <= 6; }));
        assert(d.min() == 1 && d.max() == 6 && d == d);
        decltype(d)::param_type const other{wr::discrete_uniform({.lower = 10, .upper = 12})};
        assert(d(g, other) >= 10);
        auto u = wr::as_distribution(wr::uniform({.lower = 0.0, .upper = 1.0}));
        assert(u.min() == 0.0 && u.max() == 1.0 && u(g) < 1.0);
    }

    // ---- Empreintes ---------------------------------------------------------------
    bool ok = true;
    {
        wr::xoshiro256pp x(11);
        wr::uniform const U({.lower = -3.5, .upper = 1e10});
        wr::uniform<float> const F({.lower = 0, .upper = 7});
        wr::discrete_uniform const D({.lower = -1000, .upper = 1000000});
        wr::bernoulli const B({.p = 0.125});
        ok &= check("uniform", fingerprint([&] { return std::bit_cast<std::uint64_t>(U(x)); }), 0xe9028c88f6ddde7f);
        ok &= check("uniform float", fingerprint([&] { return std::uint64_t{std::bit_cast<std::uint32_t>(F(x))}; }), 0x53a47c22e845da0b);
        ok &= check("discrete_uniform", fingerprint([&] { return static_cast<std::uint64_t>(D(x)); }), 0x3f058c88c4d229e8);
        ok &= check("bernoulli", fingerprint([&] { return std::uint64_t{wr::native(B, x)}; }), 0x88e84634bf3325c5);
        ok &= check("uniform quantile", fingerprint([&] {
            double const p = wr::uniform_real(x, wsl::exact);
            return std::bit_cast<std::uint64_t>(wr::quantile(U, p)) ^ std::bit_cast<std::uint64_t>(wr::cdf(U, wr::quantile(U, p)));
        }), 0x40887993489ea177);
        wr::exponential const E({.rate = 0.75});
        wr::exponential<float> const EF({.rate = 3});
        wr::normal const N({.mean = 10, .stddev = 0.5});
        wr::normal<float> const NF({.mean = 0, .stddev = 1});
        auto bits = [](auto v) { if constexpr (sizeof(v) == 8) return std::bit_cast<std::uint64_t>(v); else return std::uint64_t{std::bit_cast<std::uint32_t>(v)}; };
        ok &= check("exponential", fingerprint([&] { return bits(wr::native(E, x)); }), 0x93410ac575f17399);
        ok &= check("exponential float", fingerprint([&] { return bits(wr::native(EF, x)); }), 0x707f86aa03e3abe5);
        ok &= check("exponential fast", fingerprint([&] { return bits(wr::sample(E, x, wsl::approximate)); }), 0xf0883bbac7b578e0);
        ok &= check("normal", fingerprint([&] { return bits(wr::native(N, x)); }), 0xb505dc877374dd0a);
        ok &= check("normal float", fingerprint([&] { return bits(wr::native(NF, x)); }), 0x1c5f8b545c8869d7);
        ok &= check("normal fast", fingerprint([&] { return bits(wr::sample(N, x, wsl::approximate)); }), 0xd599f5eebf8f926d);
        ok &= check("normal functions", fingerprint([&] {
            double const p = wr::uniform_real(x, wsl::exact), y = wr::quantile(N, p);
            return bits(y) ^ bits(wr::cdf(N, y)) ^ (bits(wr::ccdf(N, y)) << 1) ^ (bits(wr::pdf(N, y)) << 2) ^ (bits(wr::log_pdf(N, y)) << 3);
        }, 1 << 13), 0x16ac969b0974669b);
        ok &= check("exponential functions", fingerprint([&] {
            double const p = wr::uniform_real(x, wsl::exact), y = wr::quantile(E, p);
            return bits(y) ^ bits(wr::cdf(E, y)) ^ (bits(wr::ccdf(E, y)) << 1) ^ (bits(wr::pdf(E, y)) << 2) ^ (bits(wr::log_pdf(E, y)) << 3);
        }, 1 << 13), 0x7fd190e1632d6e27);
        wr::geometric const Gf({.p = 0.05});
        wr::poisson const Ps({.mean = 4.25}), Pl({.mean = 812.5});
        wr::binomial const Bs({.trials = 30, .p = 0.2}), Bl({.trials = 100000, .p = 0.7});
        ok &= check("geometric", fingerprint([&] { return static_cast<std::uint64_t>(wr::native(Gf, x)); }), 0xf0cebc4b03f0655d);
        ok &= check("poisson inversion", fingerprint([&] { return static_cast<std::uint64_t>(wr::native(Ps, x)); }), 0xab8d234f51d7332e);
        ok &= check("poisson PTRS", fingerprint([&] { return static_cast<std::uint64_t>(wr::native(Pl, x)); }), 0x4bfc82d01b403440);
        ok &= check("binomial inversion", fingerprint([&] { return static_cast<std::uint64_t>(wr::native(Bs, x)); }), 0x758aa2af1860fcd7);
        ok &= check("binomial BTRS", fingerprint([&] { return static_cast<std::uint64_t>(wr::native(Bl, x)); }), 0xc68cde3109166ca4);
        ok &= check("discrete functions", fingerprint([&] {
            auto const k = static_cast<long long>(wr::uniform_below(2000, x));
            return bits(wr::pmf(Gf, k)) ^ (bits(wr::cdf(Gf, k)) << 1) ^ (bits(wr::pmf(Pl, k)) << 2) ^ (bits(wr::log_pmf(Pl, k)) << 3)
                 ^ (bits(wr::pmf(Bl, k + 69000)) << 4) ^ (bits(wr::log_pmf(Bs, k % 31)) << 5) ^ static_cast<std::uint64_t>(wr::quantile(Gf, wr::uniform_real(x, wsl::exact)));
        }, 1 << 13), 0x994f5c2864bf1649);
        wr::gamma const Ga({.rate = 1.5, .shape = 3.75}), Ga_s({.scale = 2.0, .shape = 0.4});
        wr::gamma<float> const Ga_f({.rate = 1, .shape = 9});
        wr::beta const Be_a({.alpha = 2.5, .beta = 7.0}), Be_s({.alpha = 0.3, .beta = 0.8});
        ok &= check("gamma", fingerprint([&] { return bits(wr::native(Ga, x)); }), 0xf1ec39e9eacde03c);
        ok &= check("gamma small shape", fingerprint([&] { return bits(wr::native(Ga_s, x)); }), 0x965ae0d5f9d27c99);
        ok &= check("gamma float", fingerprint([&] { return bits(wr::native(Ga_f, x)); }), 0xe510068f6c674400);
        ok &= check("beta", fingerprint([&] { return bits(wr::native(Be_a, x)); }), 0x78abaeb973605116);
        ok &= check("beta small shapes", fingerprint([&] { return bits(wr::native(Be_s, x)); }), 0xa2504ea65588f633);
        ok &= check("gamma beta functions", fingerprint([&] {
            double const u = wr::uniform_real(x, wsl::exact), y = 20 * u;
            return bits(wr::pdf(Ga, y)) ^ (bits(wr::log_pdf(Ga_s, y)) << 1) ^ (bits(wr::pdf(Be_a, u)) << 2) ^ (bits(wr::log_pdf(Be_s, u)) << 3);
        }, 1 << 13), 0x9d4f36e7cb669b74);
        ok &= check("phase 4 functions", fingerprint([&] {
            double const u = wr::uniform_real(x, wsl::exact);
            auto const k = static_cast<long long>(wr::uniform_below(2000, x));
            return bits(wr::cdf(Ga, 10 * u)) ^ (bits(wr::quantile(Ga_s, u)) << 1) ^ (bits(wr::cdf(Be_a, u)) << 2) ^ (bits(wr::quantile(Be_s, u)) << 3)
                 ^ (bits(wr::cdf(Pl, k)) << 4) ^ static_cast<std::uint64_t>(wr::quantile(Bl, u)) ^ (bits(Ga.entropy() + Be_s.entropy()) << 5);
        }, 1 << 10), 0x1a38363ea9e3931d);
        wr::chi_squared const Cq({.degrees_of_freedom = 7});
        wr::student_t const St({.degrees_of_freedom = 3.5});
        wr::lognormal const Ln({.log_mean = -1, .log_stddev = 0.8});
        wr::cauchy const Cy({.location = 5, .scale = 0.5});
        ok &= check("chi_squared", fingerprint([&] { return bits(wr::native(Cq, x)); }), 0x8b01e521b2059e12);
        ok &= check("student_t", fingerprint([&] { return bits(wr::native(St, x)); }), 0xb833b44187ecf66a);
        ok &= check("lognormal", fingerprint([&] { return bits(wr::native(Ln, x)); }), 0xc6350ba1601077cb);
        ok &= check("cauchy", fingerprint([&] { return bits(wr::native(Cy, x)); }), 0xf82fd9546762bb11);
        ok &= check("derived functions", fingerprint([&] {
            namespace fp = wsl::random::detail::sd;                            // entrées protégées de la contraction
            double const u = wr::uniform_real(x, wsl::exact), t = fp::fmul(20.0, u) - 10;
            return bits(wr::cdf(St, t)) ^ (bits(wr::quantile(St, u)) << 1) ^ (bits(wr::pdf(Ln, fp::fmul(3.0, u))) << 2) ^ (bits(wr::quantile(Ln, u)) << 3)
                 ^ (bits(wr::cdf(Cy, t)) << 4) ^ (bits(wr::quantile(Cy, u)) << 5) ^ (bits(wr::cdf(Cq, fp::fmul(10.0, u))) << 6);
        }, 1 << 10), 0x6e71d9a455e07fde);
        wr::fisher_f const Ff({.d1 = 4, .d2 = 9});
        wr::weibull const Wb({.scale = 3, .shape = 2.2});
        wr::laplace const Lp({.location = -2, .scale = 0.75});
        wr::logistic const Lo({.location = 4, .scale = 1.5});
        wr::pareto const Pa({.scale = 2, .shape = 1.7});
        ok &= check("fisher_f", fingerprint([&] { return bits(wr::native(Ff, x)); }), 0x2d8f53b32c9947ee);
        ok &= check("weibull", fingerprint([&] { return bits(wr::native(Wb, x)); }), 0x6b695ec07d8c4fb3);
        ok &= check("laplace", fingerprint([&] { return bits(wr::native(Lp, x)); }), 0x514ae68e514769b0);
        ok &= check("logistic", fingerprint([&] { return bits(wr::native(Lo, x)); }), 0x0dc308ffc659ed1d);
        ok &= check("pareto", fingerprint([&] { return bits(wr::native(Pa, x)); }), 0x6a50d8d1248f5578);
        ok &= check("second functions", fingerprint([&] {
            namespace fp = wsl::random::detail::sd;
            double const u = wr::uniform_real(x, wsl::exact), t = fp::fmul(20.0, u) - 10, y = fp::fmul(8.0, u);
            return bits(wr::cdf(Ff, y)) ^ (bits(wr::quantile(Ff, u)) << 1) ^ (bits(wr::pdf(Wb, y)) << 2) ^ (bits(wr::quantile(Wb, u)) << 3)
                 ^ (bits(wr::cdf(Lp, t)) << 4) ^ (bits(wr::quantile(Lo, u)) << 5) ^ (bits(wr::ccdf(Pa, y + 2)) << 6);
        }, 1 << 10), 0x572d7f7b8ca7f626);
        wr::gumbel const Gu({.location = -1, .scale = 0.5});
        wr::frechet const Fr({.scale = 1.5, .shape = 2.5});
        wr::rayleigh const Ra({.scale = 0.7});
        wr::triangular const Tr({.lower = -1, .mode = 0, .upper = 4});
        wr::arcsine const Ar({.lower = 0, .upper = 1});
        wr::erlang const Er({.rate = 0.5, .shape = 3});
        ok &= check("complements", fingerprint([&] {
            return bits(wr::native(Gu, x)) ^ (bits(wr::native(Fr, x)) << 1) ^ (bits(wr::native(Ra, x)) << 2) ^ (bits(wr::native(Tr, x)) << 3) ^ (bits(wr::native(Ar, x)) << 4) ^ (bits(wr::native(Er, x)) << 5);
        }, 1 << 12), 0x8048e3cc57673536);
        ok &= check("complement functions", fingerprint([&] {
            namespace fp = wsl::random::detail::sd;
            double const u = wr::uniform_real(x, wsl::exact), t = fp::fmul(8.0, u);
            return bits(wr::cdf(Gu, t)) ^ (bits(wr::quantile(Fr, u)) << 1) ^ (bits(wr::pdf(Ra, t)) << 2) ^ (bits(wr::cquantile(Tr, u)) << 3) ^ (bits(wr::cdf(Ar, u)) << 4) ^ (bits(wr::pdf(Fr, t + 0.1)) << 5);
        }, 1 << 10), 0x31ef6ddc23202445);
        wr::inverse_gaussian const Ig({.mean = 1.5, .shape = 0.8});
        ok &= check("inverse_gaussian", fingerprint([&] { return bits(wr::native(Ig, x)); }), 0xb65cd37eac96c4a3);
        ok &= check("inverse_gaussian functions", fingerprint([&] {
            namespace fp = wsl::random::detail::sd;
            double const u = wr::uniform_real(x, wsl::exact), t = fp::fmul(40.0, u);
            return bits(wr::cdf(Ig, t)) ^ (bits(wr::ccdf(Ig, t)) << 1) ^ (bits(wr::pdf(Ig, t)) << 2) ^ (bits(wr::quantile(Ig, u)) << 3) ^ (bits(wr::cquantile(Ig, u)) << 4);
        }, 1 << 8), 0x2c052dfac3fa0be4);
        wr::von_mises const Vm({.location = 0.25, .concentration = 7.5});
        ok &= check("von_mises functions", fingerprint([&] {
            namespace fp = wsl::random::detail::sd;
            double const u = wr::uniform_real(x, wsl::exact), t = fp::fmul(6.0, u) - 2.75;
            return bits(wr::cdf(Vm, t)) ^ (bits(wr::ccdf(Vm, t)) << 1) ^ (bits(wr::quantile(Vm, u)) << 2) ^ (bits(wr::cquantile(Vm, u)) << 3);
        }, 1 << 8), 0xedcd34abd3d9dd7f);
        wr::negative_binomial const Nb({.successes = 2.5, .p = 0.3});
        wr::hypergeometric const Hg({.population = 500, .successes = 120, .draws = 60});
        ok &= check("discrete lot 3", fingerprint([&] { return bits(double(wr::native(Nb, x))) ^ (bits(double(wr::native(Hg, x))) << 1); }, 1 << 12), 0x1abf82d977d052a7);
        ok &= check("discrete lot 3 functions", fingerprint([&] {
            double const u = wr::uniform_real(x, wsl::exact);
            long long const k = static_cast<long long>(wr::uniform_below(30, x));
            return bits(wr::pmf(Nb, k)) ^ (bits(wr::cdf(Nb, k)) << 1) ^ (bits(double(wr::quantile(Nb, u))) << 2) ^ (bits(wr::pmf(Hg, k)) << 3) ^ (bits(wr::ccdf(Hg, k)) << 4) ^ (bits(double(wr::quantile(Hg, u))) << 5);
        }, 1 << 9), 0x2fa168fb952bbce9);
        wr::beta_binomial const Bb({.trials = 40, .alpha = 0.7, .beta = 2.5});
        wr::logarithmic const Lg({.p = 0.95});
        wr::discrete_gaussian const Dg({.location = -1.3, .scale = 4.2});
        ok &= check("discrete lot 4", fingerprint([&] { return bits(double(wr::native(Bb, x))) ^ (bits(double(wr::native(Lg, x))) << 1) ^ (bits(double(wr::native(Dg, x))) << 2); }, 1 << 12), 0x5715c99366cd16cf);
        ok &= check("discrete lot 4 functions", fingerprint([&] {
            double const u = wr::uniform_real(x, wsl::exact);
            long long const k = static_cast<long long>(wr::uniform_below(40, x));
            return bits(wr::pmf(Bb, k)) ^ (bits(wr::cdf(Bb, k)) << 1) ^ (bits(double(wr::quantile(Bb, u))) << 2) ^ (bits(wr::pmf(Lg, k)) << 3) ^ (bits(wr::ccdf(Lg, k)) << 4)
                 ^ (bits(double(wr::quantile(Lg, u))) << 5) ^ (bits(wr::cdf(Dg, k - 20)) << 6) ^ (bits(double(wr::quantile(Dg, u))) << 7);
        }, 1 << 9), 0x26b8cef1aa562a21);
        wr::skew_normal const Sn({.location = -0.5, .scale = 1.5, .shape = 3.0});
        ok &= check("skew_normal", fingerprint([&] { return bits(wr::native(Sn, x)); }, 1 << 12), 0x646d26c92778b5d5);
        ok &= check("skew_normal functions", fingerprint([&] {
            namespace fp = wsl::random::detail::sd;
            double const u = wr::uniform_real(x, wsl::exact), t = fp::fmul(10.0, u) - 4;
            return bits(wr::pdf(Sn, t)) ^ (bits(wr::cdf(Sn, t)) << 1) ^ (bits(wr::ccdf(Sn, t)) << 2) ^ (bits(wr::quantile(Sn, u)) << 3);
        }, 1 << 7), 0x564a5825a3a680fd);
        wr::stable const S1({.location = 0.2, .scale = 0.8, .alpha = 1.3, .beta = -0.6});
        wr::stable const S2({.scale = 3.0, .alpha = 1.0, .beta = 0.5});
        ok &= check("stable", fingerprint([&] { return bits(wr::native(S1, x)) ^ (bits(wr::native(S2, x)) << 1); }, 1 << 12), 0x3e15bab78fbb0172);
        ok &= check("stable functions", fingerprint([&] {
            namespace fp = wsl::random::detail::sd;
            double const u = wr::uniform_real(x, wsl::exact), t = fp::fmul(10.0, u) - 5;
            return bits(wr::pdf(S1, t)) ^ (bits(wr::cdf(S1, t)) << 1) ^ (bits(wr::ccdf(S2, t)) << 2);
        }, 1 << 5), 0xba852547a6b18b91);
        wr::normal const NV({.mean = 0.5, .variance = 2.0});                 // √2 : tirage exact par les carrés
        wr::normal<float> const NVF({.mean = -1.0f, .variance = 3.0f});
        ok &= check("normal variance", fingerprint([&] { return bits(wr::native(NV, x)) ^ (std::uint64_t{std::bit_cast<std::uint32_t>(wr::native(NVF, x))} << 1); }, 1 << 12), 0x2df7ca19b86ab85d);
        ok &= check("normal variance fns", fingerprint([&] {
            namespace fp = wsl::random::detail::sd;
            double const u = wr::uniform_real(x, wsl::exact), t = fp::fmul(12.0, u) - 6;
            return bits(NV.pdf(t)) ^ (bits(NV.cdf(t)) << 1) ^ (bits(NV.ccdf(t)) << 2) ^ (bits(NV.quantile(u)) << 3);
        }, 1 << 12), 0xe9abcd992c0dde7d);
        wr::skellam const SK({.mean1 = 3.5, .mean2 = 1.25});
        ok &= check("skellam", fingerprint([&] { return static_cast<std::uint64_t>(wr::native(SK, x)); }, 1 << 12), 0x3603842a3cbc0d13);
        ok &= check("skellam functions", fingerprint([&] {
            long long const k = static_cast<long long>(wr::uniform_below(30, x)) - 12;
            return bits(SK.pmf(k)) ^ (bits(SK.cdf(k)) << 1) ^ (bits(SK.ccdf(k)) << 2);
        }, 1 << 10), 0x01f121fa1149564f);
        wr::zipf const ZP({.size = 1000, .exponent = 1.3});
        ok &= check("zipf", fingerprint([&] { return static_cast<std::uint64_t>(wr::native(ZP, x)); }, 1 << 12), 0xcf20e17a6221804c);
        ok &= check("zipf functions", fingerprint([&] {
            long long const k = 1 + static_cast<long long>(wr::uniform_below(1000, x));
            return bits(ZP.pmf(k)) ^ (bits(ZP.cdf(k)) << 1) ^ (bits(ZP.ccdf(k)) << 2);
        }, 1 << 10), 0x458210ef1666a63d);
        wr::rice const RC({.scale = 1.5, .noncentrality = 2.0});
        ok &= check("rice", fingerprint([&] { return bits(wr::native(RC, x)); }, 1 << 12), 0xba4941ebb4520dcc);
        ok &= check("rice functions", fingerprint([&] {
            double const u = wr::uniform_real(x, wsl::exact), t = wsl::floating_point_arithmetic::detail::fmul(8.0, u);   // entrées protégées de la contraction
            return bits(RC.pdf(t)) ^ (bits(RC.cdf(t)) << 1) ^ (bits(RC.ccdf(t)) << 2) ^ (bits(RC.quantile(u)) << 3);
        }, 1 << 8), 0x36a6e864e0bde893);
    }
    assert(ok);

    // ---- param() : L(X.param()) == X pour toutes les lois à paramètres -------------------
    {
        auto const round_trip = [](auto const& X) { return std::remove_cvref_t<decltype(X)>(X.param()) == X; };
        assert(round_trip(wr::uniform({.lower = -1.0, .upper = 2.0})) && round_trip(wr::discrete_uniform({.lower = 1, .upper = 6})));
        assert(round_trip(wr::bernoulli({.p = 0.3})) && round_trip(wr::binomial({.trials = 10, .p = 0.3})) && round_trip(wr::geometric({.p = 0.2})));
        assert(round_trip(wr::poisson({.mean = 3.0})) && round_trip(wr::negative_binomial({.successes = 2.5, .p = 0.4})) && round_trip(wr::logarithmic({.p = 0.6})));
        assert(round_trip(wr::hypergeometric({.population = 50, .successes = 20, .draws = 10})) && round_trip(wr::beta_binomial({.trials = 8, .alpha = 2.0, .beta = 3.0})));
        assert(round_trip(wr::categorical({.p = {0.2, 0.5, 0.3}})) && round_trip(wr::discrete_gaussian({.location = 1.5, .scale = 2.0})) && round_trip(wr::rademacher({.p = 0.4})));
        assert(round_trip(wr::normal({.mean = 1.0, .stddev = 2.0})) && round_trip(wr::lognormal({.log_mean = 0.5, .log_stddev = 0.3})) && round_trip(wr::exponential({.rate = 2.0})));
        assert(round_trip(wr::gamma({.scale = 3.0, .shape = 2.0})) && round_trip(wr::gamma({.rate = wr::ratio(3, 5), .shape = 2.0})) && round_trip(wr::exponential({.scale = 3.0})));
        assert(round_trip(wr::erlang({.rate = 1.5, .shape = 3})) && round_trip(wr::beta({.alpha = 2.0, .beta = 3.0})) && round_trip(wr::chi_squared({.degrees_of_freedom = 3.0})));
        assert(round_trip(wr::student_t({.degrees_of_freedom = 4.0})) && round_trip(wr::fisher_f({.d1 = 3.0, .d2 = 7.0})) && round_trip(wr::cauchy({.location = 1.0, .scale = 2.0})));
        assert(round_trip(wr::laplace({.location = 1.0, .scale = 2.0})) && round_trip(wr::logistic({.location = 1.0, .scale = 2.0})) && round_trip(wr::gumbel({.location = 1.0, .scale = 2.0})));
        assert(round_trip(wr::frechet({.location = 1.0, .scale = 2.0, .shape = 3.0})) && round_trip(wr::weibull({.scale = 2.0, .shape = 1.5})) && round_trip(wr::pareto({.scale = 2.0, .shape = 2.5})));
        assert(round_trip(wr::rayleigh({.scale = 2.0})) && round_trip(wr::triangular({.lower = 0.0, .mode = 1.0, .upper = 3.0})) && round_trip(wr::arcsine({.lower = -1.0, .upper = 1.0})));
        assert(round_trip(wr::inverse_gaussian({.mean = 2.0, .shape = 3.0})) && round_trip(wr::skew_normal({.location = 1.0, .scale = 2.0, .shape = -3.0})));
        assert(round_trip(wr::stable({.location = 1.0, .scale = 2.0, .alpha = 1.5, .beta = 0.5})) && round_trip(wr::von_mises({.location = 0.5, .concentration = 4.0})));
        assert(round_trip(wr::von_mises_fisher<double, 3>({.mean_direction = {0.0, 0.0, 1.0}, .concentration = 10.0})) && round_trip(wr::dirichlet<double, 3>({.alpha = {1.0, 2.0, 3.0}})));
        assert(round_trip(wr::multinomial<std::int64_t, 3>({.trials = 10, .p = {0.2, 0.3, 0.5}})) && round_trip(wr::multivariate_normal<double, 2>({.mean = {0.0, 1.0}, .covariance = {{1.0, 0.3}, {0.3, 2.0}}})));
        assert(round_trip(wr::multivariate_student_t<double, 2>({.location = {0.0, 1.0}, .scale = {{1.0, 0.3}, {0.3, 2.0}}, .degrees_of_freedom = 4.0})));
        assert(round_trip(wr::wishart<double, 2>({.scale = {{1.0, 0.2}, {0.2, 2.0}}, .degrees_of_freedom = 5.0})) && round_trip(wr::inverse_wishart<double, 2>({.scale = {{1.0, 0.2}, {0.2, 2.0}}, .degrees_of_freedom = 6.0})));
        assert(round_trip(wr::uniform_sphere({.dimension = 3})) && round_trip(wr::uniform_ball({.dimension = 3})) && round_trip(wr::ewens({.size = 10, .theta = 2.0})) && round_trip(wr::mallows({.size = 10, .q = 0.5})));
        assert(round_trip(wr::uniform_subset({.size = 10})) && round_trip(wr::bernoulli_subset({.size = 10, .p = 0.3})) && round_trip(wr::uniform_combination({.size = 10, .count = 3})));
        assert(round_trip(wr::uniform_word({.length = 5, .alphabet = 4})) && round_trip(wr::uniform_set_partition({.size = 6})) && round_trip(wr::uniform_composition({.total = 10, .parts = 3})));
        assert(round_trip(wr::uniform_integer_partition({.total = 10})) && round_trip(wr::dirichlet_multinomial({.trials = 10, .alpha = {1.0, 2.0}})) && round_trip(wr::ewens_partition({.size = 10, .theta = 1.5})));
        assert(round_trip(wr::erdos_renyi({.vertices = 10, .p = 0.3})) && round_trip(wr::uniform_graph({.vertices = 10, .edges = 12})) && round_trip(wr::barabasi_albert({.vertices = 10, .edges = 2})));
        assert(round_trip(wr::random_regular_graph({.vertices = 10, .degree = 3})) && round_trip(wr::uniform_tree({.vertices = 10})) && round_trip(wr::uniform_binary_tree({.vertices = 10})));
        static_assert(std::is_same_v<wr::normal<>::param_type, wr::normal<>::parameters> && wr::normal<>().stddev() == 1 && wr::cauchy<>().scale() == 1);
    }

    // ---- Constantes indéfinies (décision 3b) : +∞ si l'intégrale diverge vers +∞, NaN si elle
    // n'a pas de sens (∞ − ∞, 0/0) ; types : T pour les lois continues, double pour les discrètes.
    {
        constexpr double inf = INFINITY;
        auto const nan = [](auto v) { return std::isnan(v); };
        // lois dégénérées (variance nulle) : asymétrie et kurtosis 0/0
        assert(nan(wr::binomial({.trials = 5, .p = 0.0}).skewness()) && nan(wr::binomial({.trials = 0, .p = 0.3}).excess_kurtosis()));
        assert(nan(wr::geometric({.p = 1.0}).skewness()) && nan(wr::geometric({.p = 1.0}).excess_kurtosis()));
        assert(nan(wr::negative_binomial({.successes = 2.0, .p = 1.0}).skewness()) && nan(wr::discrete_uniform({.lower = 3, .upper = 3}).skewness()));
        // queues lourdes
        assert(wr::fisher_f({.d1 = 3.0, .d2 = 2.0}).mean() == inf && nan(wr::fisher_f({.d1 = 3.0, .d2 = 2.0}).variance()));
        assert(wr::fisher_f({.d1 = 3.0, .d2 = 5.0}).skewness() == inf && wr::fisher_f({.d1 = 3.0, .d2 = 7.0}).excess_kurtosis() == inf && nan(wr::fisher_f({.d1 = 3.0, .d2 = 4.0}).skewness()));
        assert(wr::pareto({.shape = 0.5}).mean() == inf && nan(wr::pareto({.shape = 0.5}).variance()) && wr::pareto({.shape = 1.5}).variance() == inf);
        assert(wr::pareto({.shape = 2.5}).skewness() == inf && wr::pareto({.shape = 3.5}).excess_kurtosis() == inf && nan(wr::pareto({.shape = 1.5}).skewness()));
        assert(wr::frechet({.shape = 0.5}).mean() == inf && nan(wr::frechet({.shape = 0.5}).variance()) && wr::frechet({.shape = 1.5}).variance() == inf);
        assert(nan(wr::student_t({.degrees_of_freedom = 1.0}).mean()) && wr::student_t({.degrees_of_freedom = 1.5}).variance() == inf && wr::student_t({.degrees_of_freedom = 3.5}).excess_kurtosis() == inf);
        assert(wr::stable({.alpha = 0.5, .beta = 1.0}).mean() == inf && wr::stable({.alpha = 0.5, .beta = -1.0}).mean() == -inf && nan(wr::stable({.alpha = 0.5}).mean()));
        assert(nan(wr::stable({.alpha = 0.5}).variance()) && wr::stable({.alpha = 1.5}).variance() == inf && wr::stable({.alpha = 1.5}).mean() == 0);
        auto const M = wr::inverse_wishart<double, 2>({.scale = {{1.0, 0.2}, {0.2, 2.0}}, .degrees_of_freedom = 2.5}).mean();
        assert((M[0, 0] == inf && M[1, 1] == inf && nan(M[0, 1])));
        // types des constantes
        static_assert(std::is_same_v<decltype(wr::von_mises<float>({.concentration = 1.0f}).entropy()), float>);
        static_assert(std::is_same_v<decltype(wr::multivariate_normal<float, 2>({.mean = {0.0f, 0.0f}, .covariance = {{1.0f, 0.0f}, {0.0f, 1.0f}}}).entropy()), float>);
        static_assert(std::is_same_v<decltype(wr::product_measure(wr::normal<float>(), wr::exponential<float>()).entropy()), float>);
        static_assert(std::is_same_v<decltype(wr::product_measure(wr::normal<float>(), wr::geometric({.p = 0.3})).entropy()), double>);
        static_assert(std::is_same_v<decltype(wr::geometric({.p = 0.3}).mean()), double> && std::is_same_v<decltype(wr::geometric({.p = 0.3}).median()), std::int64_t>);
        static_assert(noexcept(std::declval<wr::stable<> const&>().median()) && noexcept(std::declval<wr::geometric<> const&>().median()));
    }

    // ---- Normale par sa variance (décision 5, N4) -----------------------------------------
    {
        namespace fp = wsl::random::detail::sd;
        auto const close = [](double a, double b, double tol) { return std::abs(a - b) <= tol * std::max(std::abs(a), std::abs(b)); };
        // forme canonique : la variance exacte
        assert(wr::normal({.mean = 1.0, .variance = 4.0}) == wr::normal({.mean = 1.0, .stddev = 2.0}) && wr::normal({.variance = 4.0}).param().stddev);
        assert(!(wr::normal({.variance = 2.0}) == wr::normal({.stddev = wsl::special_functions::detail::pinned_portable::sqrt(2.0)})));      // √2 arrondi : autre loi
        wr::normal const V({.mean = 0.5, .variance = 2.0});
        assert(V.variance() == 2 && V.stddev() == wsl::special_functions::detail::pinned_portable::sqrt(2.0) && !V.param().stddev && *V.param().variance == 2 && wr::normal(V.param()) == V);
        assert(throws_domain_error([] { (void)wr::normal({.stddev = 1.0, .variance = 1.0}); }) && throws_domain_error([] { (void)wr::normal({.variance = -1.0}); }));
        // algèbre : les variances s'additionnent exactement
        assert(wr::convolution(wr::normal(), wr::normal()) == wr::normal({.variance = 2.0}));
        assert(wr::affine(wr::normal(), {.scale = 3}) == wr::normal({.stddev = 3.0}) && wr::affine(wr::normal({.variance = 2.0}), {.scale = 0.5}) == wr::normal({.variance = 0.5}));
        // fonctions : à 1 ulp de celles de σ = √2 arrondi (l'écart de σ est sous 2^-53)
        wr::normal const R({.mean = 0.5, .stddev = wsl::special_functions::detail::pinned_portable::sqrt(2.0)});
        for (double t = -8; t <= 8; t += 0.37)
        {
            assert(close(V.cdf(t), R.cdf(t), 4e-16 * (1 + std::abs(t))) && close(V.pdf(t), R.pdf(t), 4e-16 * (1 + t * t)));
            double const p = V.cdf(t);
            if (p > 1e-300 && p < 0.9)
                assert(close(V.quantile(p), t, 1e-14 * (1 + std::abs(t))));
        }
        // arrondi exact par les carrés, contre μ + σ t en double-double (σ = √v à 2^-104) :
        // désaccord possible seulement à moins de 2^-100 d'un milieu (probabilité ~2^-50 par point)
        wr::xoshiro256pp h(77);
        auto const square = wsl::probability::detail::exact_square_of(wsl::probability::detail::normalized(2.0, 0.0));
        wsl::probability::detail::scaled const sigma = wsl::probability::detail::sqrt(wsl::probability::detail::normalized(2.0, 0.0));
        for (int i = 0; i < 20000; ++i)
        {
            std::uint64_t const w = (h() >> 11) >> (h() % 50);
            int const b = static_cast<int>(h() % 80);
            bool const negative = (h() & 1) != 0;
            double const mu = i % 3 == 0 ? 0.0 : fp::fmul(wr::uniform_real(h, wsl::exact) - 0.5, 1e3);
            double const exact = wsl::probability::detail::round_root_affine<double>(mu, square, wsl::special_functions::detail::pinned_portable::sqrt(2.0), negative, wsl::probability::detail::big(w), b);
            if (w == 0)
            {
                assert(exact == mu);
                continue;
            }
            auto const t = sigma * wsl::probability::detail::normalized(static_cast<double>(w), 0.0);
            wsl::probability::detail::scaled const st{t.hi, t.lo, t.e - b};
            double const reference = wsl::probability::detail::affine_scaled(mu, negative, st);
            assert(exact == reference);
        }
        assert(wsl::probability::detail::lazy_sample(V, h).root && !wsl::probability::detail::lazy_sample(R, h).root);   // somme exacte : terme √v V
        assert(kolmogorov_smirnov(V, [&] { std::vector<double> d; for (int i = 0; i < 20000; ++i) d.push_back(V(h)); return d; }()) < 1.95);
    }

    // ---- Skellam ----------------------------------------------------------------------------------
    {
        wr::skellam const S({.mean1 = 3.0, .mean2 = 1.5});
        assert(S.mean() == 1.5 && S.variance() == 4.5 && wr::skellam(S.param()) == S);
        assert(std::abs(chi_square(S, g, 50000, -6, 10)) < 4);
        assert(std::abs(chi_square(wr::skellam({.mean1 = 0.4, .mean2 = 6.0}), g, 50000, -14, 3)) < 4);
        assert(std::abs(chi_square(wr::skellam({.mean1 = 200.0, .mean2 = 180.0}), g, 50000, -40, 80)) < 4);
        // symétrie : P_{μ₁,μ₂}(k) = P_{μ₂,μ₁}(−k), au bit près
        wr::skellam const R({.mean1 = 1.5, .mean2 = 3.0});
        for (long long k = -10; k <= 10; ++k)
            assert(S.pmf(k) == R.pmf(-k) && std::abs(S.cdf(k) - R.ccdf(-k - 1)) <= 2.3e-16 * S.cdf(k));
        for (double u : {1e-12, 0.3, 0.5, 0.97})
        {
            long long const k = S.quantile(u);
            assert(S.cdf(k) >= u && S.cdf(k - 1) < u);
        }
        assert(throws_domain_error([] { wr::skellam X({.mean1 = 0.0, .mean2 = 1.0}); }));
        assert(wr::geometric({.p = 0.999}).quantile(0x1p-1073) == 0 && wr::geometric({.p = 0.5}).quantile(0x1p-1074) == 0);   // u sous-normal
    }

    // ---- Zipf ---------------------------------------------------------------------------------------
    {
        wr::zipf const Z({.size = 50, .exponent = 1.1});
        assert(wr::zipf(Z.param()) == Z && Z.mode() == 1 && Z.cdf(50) == 1 && Z.ccdf(0) == 1);
        assert(std::abs(chi_square(Z, g, 50000, 1, 50)) < 4);
        assert(std::abs(chi_square(wr::zipf({.size = 10, .exponent = 0.0}), g, 50000, 1, 10)) < 4);       // uniforme
        assert(std::abs(chi_square(wr::zipf({.size = 1L << 40, .exponent = 2.5}), g, 50000, 1, 30)) < 4);
        assert(std::abs(chi_square(wr::zipf({.size = 200, .exponent = 1.0}), g, 50000, 1, 200)) < 4);
        for (double u : {1e-12, 0.3, 0.5, 0.97})
        {
            long long const k = Z.quantile(u);
            assert(Z.cdf(k) >= u && (k == 1 || Z.cdf(k - 1) < u));
        }
        assert(throws_domain_error([] { wr::zipf X({.size = 0, .exponent = 1.0}); }) && throws_domain_error([] { wr::zipf X({.size = 10, .exponent = -1.0}); }));
    }

    // ---- Rice ---------------------------------------------------------------------------------------
    {
        auto draws = [&](auto const& law) { std::vector<double> v; for (int i = 0; i < 20000; ++i) v.push_back(double(wr::native(law, g))); return v; };
        for (auto [s, n] : {std::pair{1.0, 0.0}, {1.5, 2.0}, {0.5, 0.2}, {2.0, 40.0}})
        {
            wr::rice const R({.scale = s, .noncentrality = n});
            assert(wr::rice(R.param()) == R && kolmogorov_smirnov(R, draws(R)) < 1.95);
            for (double u : {1e-300, 1e-9, 0.3, 0.5, 0.97})
            {
                double const y = R.quantile(u), z = R.cquantile(u);
                // encadrement (le conditionnement de F dans les queues interdit un aller-retour serré)
                assert(R.cdf(std::nextafter(y, 0.0)) <= u && u <= R.cdf(std::nextafter(y, 1e300)));
                assert(R.ccdf(std::nextafter(z, 1e300)) <= u && u <= R.ccdf(std::nextafter(z, 0.0)));
            }
        }
        wr::rice const R0({.scale = 2.0, .noncentrality = 0.0});
        wr::rayleigh const Ra({.scale = 2.0});
        for (double t : {1e-200, 0.3, 2.0, 9.0, 70.0})
            assert(R0.pdf(t) == Ra.pdf(t) && R0.cdf(t) == Ra.cdf(t) && R0.ccdf(t) == Ra.ccdf(t));
        assert(R0.mean() == Ra.mean() && R0.variance() == Ra.variance());
        assert(throws_domain_error([] { wr::rice X({.scale = 0.0, .noncentrality = 1.0}); }) && throws_domain_error([] { wr::rice X({.scale = 1.0, .noncentrality = -1.0}); }));
    }
    return 0;
}
