// Lois multivariées de wsl::random : Dirichlet, multinomiale, normale multivariée,
// uniformes sur la sphère et la boule ; valeurs (mpmath), moments empiriques,
// marginales (Kolmogorov–Smirnov), empreintes des tirages.

#include <wsl/random.hpp>

#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <stdexcept>
#include <type_traits>
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

    std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

    template <class F>
    bool throws_domain_error(F f)
    {
        try { f(); } catch (std::domain_error const&) { return true; }
        return false;
    }

    // |a − b| <= n ulp(b)
    bool ulps(double a, double b, double n = 2) { return std::abs(a - b) <= n * (std::nextafter(std::abs(b), INFINITY) - std::abs(b)); }
    bool close(double a, double b, double relative) { return std::abs(a - b) <= relative * std::abs(b); }

    // √n · écart maximal de la répartition empirique de v à F.
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
}

static_assert(std::is_same_v<wsl::dirichlet<double, 3>::value_type, std::array<double, 3>>);
static_assert(std::is_same_v<wsl::dirichlet<float>::value_type, std::vector<float>>);
static_assert(std::is_same_v<wsl::multinomial<std::int64_t, 4>::value_type, std::array<std::int64_t, 4>>);
static_assert(wsl::probability_law<wsl::multivariate_normal<double, 2>> && wsl::probability_law<wsl::uniform_sphere<>>);
static_assert(wsl::has_pdf<wsl::dirichlet<double, 3>> && wsl::has_pmf<wsl::multinomial<>> && wsl::has_entropy<wsl::uniform_ball<double, 4>>);
static_assert(!wsl::has_cdf<wsl::multivariate_normal<>>);

int main()
{
    wsl::xoshiro256pp g(20261118);
    constexpr int n = 20000;

    // ---- Dirichlet -----------------------------------------------------------------
    {
        wsl::dirichlet<double, 3> const D({.alpha = {1.0, 2.0, 3.0}});
        assert(D.dimension() == 3 && ulps(D.mean()[0], 1.0 / 6, 1) && D.mean()[2] == 0.5);
        assert(ulps(D.variance()[0], 5.0 / 252, 1) && ulps(D.covariance()(0, 1), -2.0 / 252, 1) && D.covariance()(1, 0) == D.covariance()(0, 1));
        assert(ulps(pdf(D, {0.2, 0.3, 0.5}), 4.5) && ulps(log_pdf(D, {0.2, 0.3, 0.5}), std::log(4.5), 4));
        assert(ulps(D.entropy(), -1.2443445622221007));
        assert(pdf(D, {0.2, 0.3, 0.6}) == 0 && pdf(D, {-0.1, 0.6, 0.5}) == 0 && log_pdf(D, {0.2, 0.3, 0.6}) == -INFINITY);
        assert(D.marginal(1) == wsl::beta({.alpha = 2.0, .beta = 4.0}) && std::isnan(D.mode()[0]));
        assert((wsl::dirichlet<double, 3>({.alpha = {2.0, 3.0, 5.0}}).mode()[2] == 4.0 / 7));
        wsl::dirichlet<double, wsl::dynamic> const E({.alpha = std::vector<double>(4, 0.5)});
        assert(E.dimension() == 4 && ulps(E.entropy(), -2.4831289505409809));
        // Moments empiriques et marginale (bêta).
        std::vector<double> first;
        double m[3] = {}, c01 = 0;
        for (int i = 0; i < n; ++i)
        {
            auto const x = sample(D, g);
            assert(std::abs(x[0] + x[1] + x[2] - 1) < 1e-15 && x[0] >= 0 && x[1] >= 0);
            for (int j = 0; j < 3; ++j)
                m[j] += x[j];
            c01 += (x[0] - 1.0 / 6) * (x[1] - 1.0 / 3);
            first.push_back(x[0]);
        }
        assert(std::abs(m[0] / n - 1.0 / 6) < 0.004 && std::abs(m[2] / n - 0.5) < 0.005 && std::abs(c01 / n + 2.0 / 252) < 0.001);
        auto const B = D.marginal(0);
        assert(kolmogorov_smirnov(first, [&](double v) { return cdf(B, v); }) < 1.95);
        assert(throws_domain_error([] { wsl::dirichlet<double, 2> X({.alpha = {1.0, 0.0}}); }));
        assert(throws_domain_error([] { wsl::dirichlet<double> X({.alpha = {1.0}}); }));
    }

    // ---- Multinomiale --------------------------------------------------------------
    {
        wsl::multinomial<std::int64_t, 3> const M({.trials = 10, .p = {0.2, 0.3, 0.5}});
        assert(ulps(pmf(M, {2, 3, 5}), 0.08505) && pmf(M, {2, 3, 4}) == 0 && pmf(M, {-1, 6, 5}) == 0 && log_pmf(M, {11, 0, -1}) == -INFINITY);
        assert(M.mean()[2] == 5 && ulps(M.covariance()(0, 1), -0.6, 1) && ulps(M.variance()[1], 2.1, 1) && M.marginal(1).p() == 0.3);
        double total = 0;                                                   // Σ des masses = 1
        for (std::int64_t a = 0; a <= 10; ++a)
            for (std::int64_t b = 0; a + b <= 10; ++b)
                total += pmf(M, {a, b, 10 - a - b});
        assert(close(total, 1, 1e-14));
        double mean1 = 0;
        std::vector<double> c;
        for (int i = 0; i < n; ++i)
        {
            auto const k = sample(M, g);
            assert(k[0] + k[1] + k[2] == 10 && k[0] >= 0 && k[1] >= 0 && k[2] >= 0);
            mean1 += double(k[1]);
        }
        assert(std::abs(mean1 / n - 3) < 0.04);
        wsl::multinomial<std::int64_t> const L({.trials = 1000000000, .p = {0.25, 0.25, 0.5}});
        auto const k = sample(L, g);
        assert(k[0] + k[1] + k[2] == 1000000000 && std::abs(double(k[2]) - 5e8) < 2e5);
        assert(pmf(L, {250000000, 250000000, 500000000}) > 0);
        assert(throws_domain_error([] { wsl::multinomial<std::int64_t, 2> X({.trials = 3, .p = {0.2, 0.3}}); }));
        assert(throws_domain_error([] { wsl::multinomial<std::int64_t, 2> X({.trials = -1, .p = {0.5, 0.5}}); }));
    }

    // ---- Normale multivariée -------------------------------------------------------
    {
        wsl::multivariate_normal<double, 2> const X({.mean = {1, 2}, .covariance = {{1, 0.5}, {0.5, 2}}});
        assert(ulps(pdf(X, {1.5, 1.0}), 0.067941140344700179) && ulps(X.entropy(), 3.1176849603770568));
        assert(ulps(log_pdf(X, {1.5, 1.0}), std::log(0.067941140344700179), 4) && X.marginal(1).sd() == std::sqrt(2.0) && X.marginal(0).mean() == 1);
        double s[2] = {}, c = 0;
        std::vector<double> second;
        for (int i = 0; i < n; ++i)
        {
            auto const y = sample(X, g);
            s[0] += y[0];
            s[1] += y[1];
            c += (y[0] - 1) * (y[1] - 2);
            second.push_back(y[1]);
        }
        assert(std::abs(s[0] / n - 1) < 0.03 && std::abs(s[1] / n - 2) < 0.04 && std::abs(c / n - 0.5) < 0.04);
        auto const N1 = X.marginal(1);
        assert(kolmogorov_smirnov(second, [&](double v) { return cdf(N1, v); }) < 1.95);
        // Dimension dynamique : diagonale, densité produit des marginales.
        wsl::multivariate_normal<double> const D({.mean = {0, 1, -1}, .covariance = {{4, 0, 0}, {0, 1, 0}, {0, 0, 9}}});
        double const p = pdf(wsl::normal({.mean = 0.0, .sd = 2.0}), 0.5) * pdf(wsl::normal({.mean = 1.0, .sd = 1.0}), 0.0) * pdf(wsl::normal({.mean = -1.0, .sd = 3.0}), 2.0);
        assert(close(pdf(D, {0.5, 0.0, 2.0}), p, 1e-15) && D.mahalanobis_squared({2.0, 1.0, -1.0}) == 1);
        assert(throws_domain_error([] { wsl::multivariate_normal<double, 2> Y({.mean = {0, 0}, .covariance = {{1, 2}, {2, 1}}}); }));   // non définie positive
        assert(throws_domain_error([] { wsl::multivariate_normal<double, 2> Y({.mean = {0, 0}, .covariance = {{1, 0.1}, {0.2, 1}}}); })); // non symétrique
        assert(throws_domain_error([] { wsl::multivariate_normal<double> Y({.mean = {0, 0}, .covariance = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}); }));
        assert(throws_domain_error([] { wsl::square_matrix<double, 2> M{{1, 0}, {0}}; }));
    }

    // ---- Sphère et boule -----------------------------------------------------------
    {
        wsl::uniform_sphere<double, 3> const S;
        assert(ulps(S.surface_area(), 4 * std::numbers::pi) && ulps(pdf(S, {0.0, 0.0, 1.0}), 1 / (4 * std::numbers::pi)) && pdf(S, {0.0, 0.0, 0.9}) == 0);
        std::vector<double> z;                                              // Archimède : x₃ uniforme sur [−1, 1]
        for (int i = 0; i < n; ++i)
        {
            auto const x = sample(S, g);
            assert(std::abs(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1) < 4e-16);
            z.push_back(x[2]);
        }
        assert(kolmogorov_smirnov(z, [](double v) { return (v + 1) / 2; }) < 1.95);
        wsl::uniform_ball<double> const B({.dimension = 10});
        assert(ulps(B.volume(), 2.5501640398773454) && pdf(B, std::vector<double>(10, 0.0)) == 1 / B.volume() && pdf(B, std::vector<double>(10, 0.5)) == 0);
        std::vector<double> r;                                              // |x|ⁿ uniforme sur [0, 1]
        for (int i = 0; i < n; ++i)
        {
            auto const x = sample(B, g);
            double q = 0;
            for (double v : x)
                q += v * v;
            r.push_back(std::pow(q, 5.0));
        }
        assert(kolmogorov_smirnov(r, [](double v) { return v; }) < 1.95);
        wsl::uniform_sphere<float> const F({.dimension = 2});
        assert(F.variance()[0] == 0.5f && sample(F, g).size() == 2);
        assert(throws_domain_error([] { wsl::uniform_sphere<double, 3> X({.dimension = 2}); }));
    }

    // ---- Empreintes ----------------------------------------------------------------
    bool ok = true;
    {
        wsl::xoshiro256pp x(42);
        wsl::dirichlet<double, 3> const D({.alpha = {0.5, 2.0, 7.0}});
        wsl::multinomial<std::int64_t> const M({.trials = 50, .p = {0.1, 0.2, 0.3, 0.4}});
        wsl::multivariate_normal<double, 3> const X({.mean = {1, -1, 0}, .covariance = {{2, 0.3, 0.1}, {0.3, 1, -0.2}, {0.1, -0.2, 0.5}}});
        wsl::uniform_sphere<double, 4> const S;
        wsl::uniform_ball<double, 3> const B;
        ok &= check("dirichlet", fingerprint([&](fnv1a& h) { for (double v : sample(D, x)) h.add(bits(v)); }), 0xd64c83a93edccedb);
        ok &= check("multinomial", fingerprint([&](fnv1a& h) { for (auto v : sample(M, x)) h.add(static_cast<std::uint64_t>(v)); }), 0x9160e76ae8ecf5e9);
        ok &= check("multivariate_normal", fingerprint([&](fnv1a& h) { for (double v : sample(X, x)) h.add(bits(v)); }), 0x3b220932212ad366);
        ok &= check("mvn (fast)", fingerprint([&](fnv1a& h) { for (double v : sample(X, x, wsl::fast)) h.add(bits(v)); }), 0x885b367a151c3b43);
        ok &= check("uniform_sphere", fingerprint([&](fnv1a& h) { for (double v : sample(S, x)) h.add(bits(v)); }), 0xd472b57b4039a242);
        ok &= check("uniform_ball", fingerprint([&](fnv1a& h) { for (double v : sample(B, x)) h.add(bits(v)); }), 0x9a77e07637a1c221);
        ok &= check("functions", fingerprint([&](fnv1a& h) {
            namespace fp = wsl::special::detail;
            double const u = wsl::uniform_real(x), v = fp::fmul(0.5, u);
            h.add(bits(pdf(D, {v, 0.25, 0.75 - v})) ^ (bits(log_pdf(X, {u, v, -u})) << 1) ^ (bits(pmf(M, {5, 10, 15, 20})) << 2));
        }, 1 << 9), 0x22fd3630ae225505);
    }
    assert(ok);
    std::printf("test_random_multivariate : ok\n");
}
