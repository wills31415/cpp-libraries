// Lois multivariées de wsl::random : Dirichlet, multinomiale, normale multivariée,
// uniformes sur la sphère et la boule ; valeurs (mpmath), moments empiriques,
// marginales (Kolmogorov–Smirnov), empreintes des tirages.

#include <wsl/random.hpp>

#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <complex>
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

static_assert(std::is_same_v<wsl::random::dirichlet<double, 3>::value_type, std::array<double, 3>>);
static_assert(std::is_same_v<wsl::random::dirichlet<float>::value_type, std::vector<float>>);
static_assert(std::is_same_v<wsl::random::multinomial<std::int64_t, 4>::value_type, std::array<std::int64_t, 4>>);
static_assert(wsl::random::probability_law<wsl::random::multivariate_normal<double, 2>> && wsl::random::probability_law<wsl::random::uniform_sphere<>>);
static_assert(wsl::random::has_pdf<wsl::random::dirichlet<double, 3>> && wsl::random::has_pmf<wsl::random::multinomial<>> && wsl::random::has_entropy<wsl::random::uniform_ball<double, 4>>);
static_assert(!wsl::random::has_cdf<wsl::random::multivariate_normal<>>);

int main()
{
    wsl::random::xoshiro256pp g(20261118);
    constexpr int n = 20000;

    // ---- Dirichlet -----------------------------------------------------------------
    {
        wsl::random::dirichlet<double, 3> const D({.alpha = {1.0, 2.0, 3.0}});
        assert(D.dimension() == 3 && ulps(D.mean()[0], 1.0 / 6, 1) && D.mean()[2] == 0.5);
        assert(ulps(D.variance()[0], 5.0 / 252, 1) && ulps(D.covariance()(0, 1), -2.0 / 252, 1) && D.covariance()(1, 0) == D.covariance()(0, 1));
        assert(ulps(pdf(D, {0.2, 0.3, 0.5}), 4.5) && ulps(log_pdf(D, {0.2, 0.3, 0.5}), std::log(4.5), 4));
        assert(ulps(D.entropy(), -1.2443445622221007));
        assert(pdf(D, {0.2, 0.3, 0.6}) == 0 && pdf(D, {-0.1, 0.6, 0.5}) == 0 && log_pdf(D, {0.2, 0.3, 0.6}) == -INFINITY);
        assert(D.marginal(1) == wsl::random::beta({.alpha = 2.0, .beta = 4.0}) && std::isnan(D.mode()[0]));
        assert((wsl::random::dirichlet<double, 3>({.alpha = {2.0, 3.0, 5.0}}).mode()[2] == 4.0 / 7));
        wsl::random::dirichlet<double, wsl::random::dynamic> const E({.alpha = std::vector<double>(4, 0.5)});
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
        assert(throws_domain_error([] { wsl::random::dirichlet<double, 2> X({.alpha = {1.0, 0.0}}); }));
        assert(throws_domain_error([] { wsl::random::dirichlet<double> X({.alpha = {1.0}}); }));
    }

    // ---- Multinomiale --------------------------------------------------------------
    {
        wsl::random::multinomial<std::int64_t, 3> const M({.trials = 10, .p = {0.2, 0.3, 0.5}});
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
        wsl::random::multinomial<std::int64_t> const L({.trials = 1000000000, .p = {0.25, 0.25, 0.5}});
        auto const k = sample(L, g);
        assert(k[0] + k[1] + k[2] == 1000000000 && std::abs(double(k[2]) - 5e8) < 2e5);
        assert(pmf(L, {250000000, 250000000, 500000000}) > 0);
        assert(throws_domain_error([] { wsl::random::multinomial<std::int64_t, 2> X({.trials = 3, .p = {0.2, 0.3}}); }));
        assert(throws_domain_error([] { wsl::random::multinomial<std::int64_t, 2> X({.trials = -1, .p = {0.5, 0.5}}); }));
    }

    // ---- Normale multivariée -------------------------------------------------------
    {
        wsl::random::multivariate_normal<double, 2> const X({.mean = {1, 2}, .covariance = {{1, 0.5}, {0.5, 2}}});
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
        wsl::random::multivariate_normal<double> const D({.mean = {0, 1, -1}, .covariance = {{4, 0, 0}, {0, 1, 0}, {0, 0, 9}}});
        double const p = pdf(wsl::random::normal({.mean = 0.0, .sd = 2.0}), 0.5) * pdf(wsl::random::normal({.mean = 1.0, .sd = 1.0}), 0.0) * pdf(wsl::random::normal({.mean = -1.0, .sd = 3.0}), 2.0);
        assert(close(pdf(D, {0.5, 0.0, 2.0}), p, 1e-15) && D.mahalanobis_squared({2.0, 1.0, -1.0}) == 1);
        assert(throws_domain_error([] { wsl::random::multivariate_normal<double, 2> Y({.mean = {0, 0}, .covariance = {{1, 2}, {2, 1}}}); }));   // non définie positive
        assert(throws_domain_error([] { wsl::random::multivariate_normal<double, 2> Y({.mean = {0, 0}, .covariance = {{1, 0.1}, {0.2, 1}}}); })); // non symétrique
        assert(throws_domain_error([] { wsl::random::multivariate_normal<double> Y({.mean = {0, 0}, .covariance = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}); }));
        assert(throws_domain_error([] { wsl::random::square_matrix<double, 2> M{{1, 0}, {0}}; }));
    }

    // ---- Sphère et boule -----------------------------------------------------------
    {
        wsl::random::uniform_sphere<double, 3> const S;
        assert(ulps(S.surface_area(), 4 * std::numbers::pi) && ulps(pdf(S, {0.0, 0.0, 1.0}), 1 / (4 * std::numbers::pi)) && pdf(S, {0.0, 0.0, 0.9}) == 0);
        std::vector<double> z;                                              // Archimède : x₃ uniforme sur [−1, 1]
        for (int i = 0; i < n; ++i)
        {
            auto const x = sample(S, g);
            assert(std::abs(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1) < 4e-16);
            z.push_back(x[2]);
        }
        assert(kolmogorov_smirnov(z, [](double v) { return (v + 1) / 2; }) < 1.95);
        wsl::random::uniform_ball<double> const B({.dimension = 10});
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
        wsl::random::uniform_sphere<float> const F({.dimension = 2});
        assert(F.variance()[0] == 0.5f && sample(F, g).size() == 2);
        assert(throws_domain_error([] { wsl::random::uniform_sphere<double, 3> X({.dimension = 2}); }));
    }

    // ---- Student multivariée --------------------------------------------------------
    {
        wsl::random::multivariate_student_t<double, 2> const X({.degrees_of_freedom = 4, .location = {0, 1}, .scale = {{1, 0.3}, {0.3, 2}}});
        assert(ulps(pdf(X, {0.5, 0.5}), 0.082453654979996737) && ulps(X.entropy(), 3.6614286874386147));
        assert(X.covariance()(0, 0) == 2 && ulps(X.covariance()(0, 1), 0.6, 1) && X.mean()[1] == 1);
        assert(std::isinf(wsl::random::multivariate_student_t<double, 2>({.degrees_of_freedom = 1.5, .location = {0, 0}, .scale = {{1, 0}, {0, 1}}}).covariance()(0, 0)));
        assert(std::isnan(wsl::random::multivariate_student_t<double, 2>({.degrees_of_freedom = 1.0, .location = {0, 0}, .scale = {{1, 0}, {0, 1}}}).mean()[0]));
        wsl::random::multivariate_student_t<double, 1> const One({.degrees_of_freedom = 3, .location = {0}, .scale = {{1}}});   // = Student
        for (double t : {-4.0, 0.0, 0.7, 12.0})
            assert(close(pdf(One, {t}), pdf(wsl::random::student_t({.degrees_of_freedom = 3.0}), t), 2e-15));
        std::vector<double> first;
        for (int i = 0; i < n; ++i)
            first.push_back(sample(X, g)[1]);
        auto const M1 = X.marginal(1);
        assert(kolmogorov_smirnov(first, [&](double v) { return cdf(M1, v); }) < 1.95);
    }

    // ---- Wishart, inverse-Wishart --------------------------------------------------
    {
        wsl::random::square_matrix<double, 2> const A{{4, 1}, {1, 9}};
        wsl::random::wishart<double, 2> const W({.degrees_of_freedom = 5, .scale = {{1, 0.2}, {0.2, 2}}});
        assert(ulps(log_pdf(W, A), -6.6844907153650172, 2) && ulps(W.entropy(), 7.8199650891973836, 2));
        assert(W.mean()(1, 1) == 10 && ulps(W.variance()(0, 1), 10.2, 1) && W.mode()(0, 0) == 2 && log_pdf(W, wsl::random::square_matrix<double, 2>{{1, 2}, {2, 1}}) == -INFINITY);
        wsl::random::wishart<double, 1> const W1({.degrees_of_freedom = 3, .scale = {{2}}});   // = gamma(3/2, échelle 4)
        for (double x : {0.5, 3.0, 20.0})
            assert(close(log_pdf(W1, wsl::random::square_matrix<double, 1>{{x}}), log_pdf(wsl::random::gamma({.shape = 1.5, .scale = 4.0}), x), 1e-14));
        double m01 = 0;
        std::vector<double> d;
        for (int i = 0; i < n; ++i)
        {
            auto const w = sample(W, g);
            assert(w(0, 1) == w(1, 0) && w(0, 0) > 0 && w(0, 0) * w(1, 1) > w(0, 1) * w(0, 1));
            m01 += w(0, 1);
            d.push_back(w(1, 1));
        }
        assert(std::abs(m01 / n - 1) < 0.06);
        auto const G = W.marginal(1);
        assert(kolmogorov_smirnov(d, [&](double v) { return cdf(G, v); }) < 1.95);
        wsl::random::inverse_wishart<double, 2> const I({.degrees_of_freedom = 6, .scale = {{1, 0.2}, {0.2, 2}}});
        assert(ulps(log_pdf(I, A), -19.926453791435178, 2) && std::abs(I.entropy() - -0.83656422870532905) < 1e-15);
        assert(ulps(I.mean()(1, 1), 2.0 / 3, 1) && ulps(I.mode()(0, 0), 1.0 / 9, 1));
        double i00 = 0;
        for (int i = 0; i < n; ++i)
        {
            auto const x = sample(I, g);
            assert(x(0, 1) == x(1, 0) && x(0, 0) > 0);
            i00 += x(0, 0);
        }
        assert(std::abs(i00 / n - 1.0 / 3) < 0.02);
        assert(throws_domain_error([] { wsl::random::wishart<double, 2> Y({.degrees_of_freedom = 0.9, .scale = {{1, 0}, {0, 1}}}); }));
        assert(throws_domain_error([] { wsl::random::inverse_wishart<double, 2> Y({.degrees_of_freedom = 3, .scale = {{1, 2}, {2, 1}}}); }));
    }

    // ---- Matrices aléatoires -------------------------------------------------------
    {
        wsl::random::haar_orthogonal<double, 3> const O;
        std::vector<double> z;
        int positive = 0;
        for (int i = 0; i < n; ++i)
        {
            auto const q = sample(O, g);
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                {
                    double s = 0;
                    for (int k = 0; k < 3; ++k)
                        s += q(k, a) * q(k, b);
                    assert(std::abs(s - (a == b)) < 1e-15);
                }
            double const det = q(0, 0) * (q(1, 1) * q(2, 2) - q(1, 2) * q(2, 1)) - q(0, 1) * (q(1, 0) * q(2, 2) - q(1, 2) * q(2, 0)) + q(0, 2) * (q(1, 0) * q(2, 1) - q(1, 1) * q(2, 0));
            positive += det > 0;
            z.push_back(q(2, 0));                                           // colonne uniforme sur S² : Archimède
        }
        assert(kolmogorov_smirnov(z, [](double v) { return (v + 1) / 2; }) < 1.95 && std::abs(positive / double(n) - 0.5) < 0.015);
        wsl::random::haar_unitary<double> const U({.dimension = 4});
        std::vector<double> u2;
        for (int i = 0; i < n; ++i)
        {
            auto const u = sample(U, g);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                {
                    std::complex<double> s = 0;
                    for (int k = 0; k < 4; ++k)
                        s += std::conj(u(k, a)) * u(k, b);
                    assert(std::abs(s - double(a == b)) < 2e-15);
                }
            u2.push_back(std::norm(u(0, 0)));                               // |u₀₀|² ~ bêta(1, n − 1)
        }
        assert(kolmogorov_smirnov(u2, [](double t) { return 1 - std::pow(1 - t, 3); }) < 1.95);
        wsl::random::goe<double, 3> const H;
        assert(ulps(pdf(H, wsl::random::square_matrix<double, 3>{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}), 0.011402639350693397, 1));
        assert(log_pdf(H, wsl::random::square_matrix<double, 3>{{0, 1, 0}, {0, 0, 0}, {0, 0, 0}}) == -INFINITY);
        double v00 = 0, v01 = 0;
        for (int i = 0; i < n; ++i)
        {
            auto const h = sample(H, g);
            assert(h(0, 2) == h(2, 0));
            v00 += h(0, 0) * h(0, 0);
            v01 += h(0, 1) * h(0, 1);
        }
        assert(std::abs(v00 / n - 1) < 0.04 && std::abs(v01 / n - 0.5) < 0.02);
        wsl::random::gue<double, 2> const K;
        assert(ulps(pdf(K, wsl::random::square_matrix<std::complex<double>, 2>{{0, 0}, {0, 0}}), 0.050660591821168886, 1));
        assert(log_pdf(K, wsl::random::square_matrix<std::complex<double>, 2>{{{0, 1}, 0}, {0, 0}}) == -INFINITY);   // diagonale non réelle
        double vr = 0;
        for (int i = 0; i < n; ++i)
        {
            auto const k = sample(K, g);
            assert(k(0, 1) == std::conj(k(1, 0)) && k(0, 0).imag() == 0);
            vr += k(0, 1).real() * k(0, 1).real();
        }
        assert(std::abs(vr / n - 0.5) < 0.02);
    }

    // ---- von Mises, von Mises–Fisher -----------------------------------------------
    {
        for (double kappa : {0.3, 4.0, 50.0})
        {
            wsl::random::von_mises const V({.mu = 1.0, .kappa = kappa});
            // Répartition par intégration (trapèzes, 2^16 pas) de la densité.
            int const steps = 1 << 16;
            double const lo = V.support().lower, h = (V.support().upper - lo) / steps;
            std::vector<double> F(steps + 1, 0.0);
            for (int i = 1; i <= steps; ++i)
                F[i] = F[i - 1] + h * (pdf(V, lo + (i - 1) * h) + pdf(V, lo + i * h)) / 2;
            assert(std::abs(F[steps] - 1) < 1e-9);
            std::vector<double> v;
            double c = 0;
            for (int i = 0; i < n; ++i)
            {
                double const x = sample(V, g);
                assert(x >= V.support().lower - 1e-15 && x <= V.support().upper + 1e-15);
                v.push_back(x);
                c += std::cos(x - 1);
            }
            assert(std::abs(c / n - V.mean_resultant_length()) < 0.02);
            assert(kolmogorov_smirnov(v, [&](double x) {
                double const r = (x - lo) / h;
                int const i = std::clamp(static_cast<int>(r), 0, steps - 1);
                return F[i] + (F[i + 1] - F[i]) * (r - i);
            }) < 1.95);
        }
        wsl::random::von_mises const U({.mu = 0.0, .kappa = 0.0});                    // uniforme sur le cercle
        assert(ulps(pdf(U, 2.0), 1 / (2 * std::numbers::pi), 1) && pdf(U, 3.5) == 0 && U.circular_variance() == 1);
        assert(throws_domain_error([] { wsl::random::von_mises X({.mu = 0.0, .kappa = -1.0}); }));
        // p = 3 : t = μᵀX de répartition (e^{κt} − e^{−κ})/(e^κ − e^{−κ}), azimut uniforme.
        for (double kappa : {0.5, 10.0, 300.0})
        {
            wsl::random::von_mises_fisher<double, 3> const M({.mean_direction = {0, 0, 1}, .kappa = kappa});
            std::vector<double> t, phi;
            for (int i = 0; i < n; ++i)
            {
                auto const x = sample(M, g);
                assert(std::abs(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1) < 1e-15);
                t.push_back(x[2]);
                phi.push_back(std::atan2(x[1], x[0]));
            }
            assert(kolmogorov_smirnov(t, [&](double z) { return std::expm1(kappa * (z + 1)) / std::expm1(2 * kappa); }) < 1.95);
            assert(kolmogorov_smirnov(phi, [](double a) { return (a + std::numbers::pi) / (2 * std::numbers::pi); }) < 1.95);
        }
        wsl::random::von_mises_fisher<double> const D({.mean_direction = {1, 2, 2, 0, -1}, .kappa = 3.5});
        assert(ulps(D.mean_direction()[1], 2 / std::sqrt(10.0), 1) && ulps(D.entropy(), 2.447886519925983, 2) && ulps(D.mean_resultant_length(), 0.53928842320413421, 2));
        double m1 = 0;
        for (int i = 0; i < n; ++i)
            m1 += sample(D, g)[1];
        assert(std::abs(m1 / n - D.mean()[1]) < 0.01);
        wsl::random::von_mises_fisher<double, 3> const Z({.mean_direction = {0, -1, 0}, .kappa = 0.0});
        assert(ulps(pdf(Z, {0.0, 1.0, 0.0}), 1 / (4 * std::numbers::pi), 1));
        assert(throws_domain_error([] { wsl::random::von_mises_fisher<double, 2> X({.mean_direction = {0, 0}, .kappa = 1.0}); }));
    }

    // ---- Empreintes ----------------------------------------------------------------
    bool ok = true;
    {
        wsl::random::xoshiro256pp x(42);
        wsl::random::dirichlet<double, 3> const D({.alpha = {0.5, 2.0, 7.0}});
        wsl::random::multinomial<std::int64_t> const M({.trials = 50, .p = {0.1, 0.2, 0.3, 0.4}});
        wsl::random::multivariate_normal<double, 3> const X({.mean = {1, -1, 0}, .covariance = {{2, 0.3, 0.1}, {0.3, 1, -0.2}, {0.1, -0.2, 0.5}}});
        wsl::random::uniform_sphere<double, 4> const S;
        wsl::random::uniform_ball<double, 3> const B;
        ok &= check("dirichlet", fingerprint([&](fnv1a& h) { for (double v : sample(D, x)) h.add(bits(v)); }), 0xd64c83a93edccedb);
        ok &= check("multinomial", fingerprint([&](fnv1a& h) { for (auto v : sample(M, x)) h.add(static_cast<std::uint64_t>(v)); }), 0x9160e76ae8ecf5e9);
        ok &= check("multivariate_normal", fingerprint([&](fnv1a& h) { for (double v : sample(X, x)) h.add(bits(v)); }), 0x3b220932212ad366);
        ok &= check("mvn (fast)", fingerprint([&](fnv1a& h) { for (double v : sample(X, x, wsl::random::fast)) h.add(bits(v)); }), 0x885b367a151c3b43);
        ok &= check("uniform_sphere", fingerprint([&](fnv1a& h) { for (double v : sample(S, x)) h.add(bits(v)); }), 0xd472b57b4039a242);
        ok &= check("uniform_ball", fingerprint([&](fnv1a& h) { for (double v : sample(B, x)) h.add(bits(v)); }), 0x9a77e07637a1c221);
        wsl::random::multivariate_student_t<double, 2> const T2({.degrees_of_freedom = 3.5, .location = {1, 0}, .scale = {{2, -0.4}, {-0.4, 1}}});
        wsl::random::wishart<double> const W({.degrees_of_freedom = 4.5, .scale = {{1, 0.3, 0}, {0.3, 2, 0.1}, {0, 0.1, 0.5}}});
        wsl::random::inverse_wishart<double, 2> const IW({.degrees_of_freedom = 5, .scale = {{1, 0.2}, {0.2, 2}}});
        wsl::random::haar_orthogonal<double, 4> const HO;
        wsl::random::haar_unitary<double, 3> const HU;
        wsl::random::goe<double, 3> const GO;
        wsl::random::gue<double, 3> const GU;
        auto const all = [](fnv1a& h, auto const& m) { for (std::size_t i = 0; i < m.size(); ++i) for (std::size_t j = 0; j < m.size(); ++j) { if constexpr (std::is_same_v<std::remove_cvref_t<decltype(m(i, j))>, double>) h.add(bits(m(i, j))); else { h.add(bits(m(i, j).real())); h.add(bits(m(i, j).imag())); } } };
        ok &= check("multivariate_student", fingerprint([&](fnv1a& h) { for (double v : sample(T2, x)) h.add(bits(v)); }), 0x66189197e5f2a049);
        ok &= check("wishart", fingerprint([&](fnv1a& h) { all(h, sample(W, x)); }), 0x3526188a6af0a6bd);
        ok &= check("inverse_wishart", fingerprint([&](fnv1a& h) { all(h, sample(IW, x)); }), 0x621e8af1f4c88f70);
        ok &= check("haar_orthogonal", fingerprint([&](fnv1a& h) { all(h, sample(HO, x)); }), 0x57ac6de9964d4ba7);
        ok &= check("haar_unitary", fingerprint([&](fnv1a& h) { all(h, sample(HU, x)); }), 0xcf6b88432f8193b9);
        ok &= check("goe", fingerprint([&](fnv1a& h) { all(h, sample(GO, x)); }), 0xddb4938359d8b2e9);
        ok &= check("gue", fingerprint([&](fnv1a& h) { all(h, sample(GU, x)); }), 0xfd93e1b6acc6374f);
        wsl::random::von_mises const VM({.mu = -1.0, .kappa = 6.0});
        wsl::random::von_mises_fisher<double, 3> const VF({.mean_direction = {0.6, 0, -0.8}, .kappa = 20.0});
        wsl::random::von_mises_fisher<double> const VD({.mean_direction = {1, 1, 1, 1, 1, 1}, .kappa = 2.0});
        ok &= check("von_mises", fingerprint([&](fnv1a& h) { h.add(bits(sample(VM, x))); }), 0xc6ba70678de3be7f);
        ok &= check("von_mises_fisher", fingerprint([&](fnv1a& h) { for (double v : sample(VF, x)) h.add(bits(v)); }), 0x1aea2eb0da6a278b);
        ok &= check("vmf (dynamique)", fingerprint([&](fnv1a& h) { for (double v : sample(VD, x)) h.add(bits(v)); }), 0xb45be65456b07a76);
        ok &= check("functions 2", fingerprint([&](fnv1a& h) {
            auto const w = sample(W, x);
            h.add(bits(log_pdf(W, w)) ^ (bits(log_pdf(T2, {w(0, 0), w(1, 2)})) << 1) ^ (bits(log_pdf(GO, sample(GO, x))) << 2));
        }, 1 << 8), 0xcebc860185777128);
        ok &= check("functions", fingerprint([&](fnv1a& h) {
            namespace fp = wsl::special::detail;
            double const u = wsl::random::uniform_real(x), v = fp::fmul(0.5, u);
            h.add(bits(pdf(D, {v, 0.25, 0.75 - v})) ^ (bits(log_pdf(X, {u, v, -u})) << 1) ^ (bits(pmf(M, {5, 10, 15, 20})) << 2));
        }, 1 << 9), 0x5e2d79ac0279d5e6);
    }
    assert(ok);
    std::printf("test_random_multivariate : ok\n");
}
