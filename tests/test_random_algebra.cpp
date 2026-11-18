// Algèbre des lois de wsl::random : quantile de survie, image affine et identités,
// troncature, mélanges ; vérifications exactes, statistiques (graines fixes) et
// empreintes des tirages.

#include <wsl/random.hpp>

#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
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
    std::uint64_t fingerprint(F f, int count = 1 << 14)
    {
        fnv1a h;
        for (int i = 0; i < count; ++i)
            h.add(f());
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

    // Kolmogorov–Smirnov : √n · écart maximal, < 1,95 au seuil de 0,1 %.
    template <class L, class G>
    double kolmogorov_smirnov(L const& law, G& g, int n = 20000)
    {
        std::vector<double> v;
        for (int i = 0; i < n; ++i)
            v.push_back(double(sample(law, g)));
        std::ranges::sort(v);
        double worst = 0;
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            double const f = cdf(law, v[i]);
            worst = std::max({worst, std::abs(f - double(i) / n), std::abs(f - double(i + 1) / n)});
        }
        return worst * std::sqrt(double(n));
    }

    bool close(double a, double b, double relative) { return std::abs(a - b) <= relative * std::abs(b); }

    // ccdf(X, cquantile(X, q)) ≈ q, et cquantile(X, q) = quantile(X, 1 − q) pour q >= ½.
    template <class L>
    void check_cquantile(L const& X, double tolerance = 1e-12)
    {
        for (double q : {0.5, 0.75, 0.9375, 1 - 0x1p-40})
            assert(cquantile(X, q) == quantile(X, 1 - q));
        for (double q : {0.25, 1e-5, 1e-30, 1e-200})
        {
            double const x = cquantile(X, q);
            if (std::isfinite(x) && x != X.support().lower && x != X.support().upper && !close(double(ccdf(X, x)), q, tolerance))
            {
                std::printf("cquantile(%g) = %a, ccdf = %.17g\n", q, x, double(ccdf(X, x)));
                assert(false);
            }
        }
        assert(cquantile(X, 1.0) == quantile(X, 0.0) && cquantile(X, 0.0) == quantile(X, 1.0));
        assert(throws_domain_error([&] { (void)cquantile(X, 1.5); }));
    }
}

// Identités : les familles fermées gardent leur type.
using normal_t = wsl::normal<double>;
static_assert(std::is_same_v<decltype(2.0 * std::declval<normal_t>() + 1), normal_t>);
static_assert(std::is_same_v<decltype(-std::declval<wsl::cauchy<>>()), wsl::cauchy<>>);
static_assert(std::is_same_v<decltype(3 - std::declval<wsl::uniform<float>>()), wsl::uniform<float>>);
static_assert(std::is_same_v<decltype(std::declval<wsl::laplace<>>() / 4), wsl::laplace<>>);
static_assert(std::is_same_v<decltype(2 * std::declval<wsl::exponential<>>()), wsl::affine_image<wsl::exponential<>>>);
static_assert(std::is_same_v<decltype(2 * (2 * std::declval<wsl::student_t<>>()) + 1), wsl::affine_image<wsl::student_t<>>>);
static_assert(wsl::affine_closed<wsl::logistic<>> && wsl::scale_closed<wsl::gamma<>> && !wsl::scale_closed<wsl::beta<>>);
static_assert(!wsl::has_mean<decltype(2.0 * std::declval<wsl::cauchy<>>())>);
static_assert(wsl::has_mean<wsl::affine_image<wsl::student_t<>>> && wsl::has_cquantile<wsl::affine_image<wsl::beta<>>>);
static_assert(!wsl::has_mean<wsl::truncated<normal_t>> && wsl::has_cquantile<wsl::truncated<normal_t>>);
static_assert(wsl::probability_law<wsl::mixture<normal_t, wsl::gamma<>>> && wsl::has_pmf<wsl::mixture<wsl::poisson<>, wsl::binomial<>>>);
static_assert(wsl::mixture<normal_t, wsl::exponential<>>::sampling == wsl::exactness::exact);
static_assert(wsl::mixture<normal_t, wsl::gamma<>>::sampling == wsl::exactness::approximate);
static_assert(!wsl::has_pdf<wsl::mixture<wsl::poisson<>>> && wsl::has_cdf<wsl::mixture<wsl::poisson<>>>);
// Évaluables à la compilation.
static_assert((3.0 * wsl::normal({.mean = 1.0, .sd = 2.0}) - 1.0).mean() == 2.0);
static_assert(cdf(2.0 * wsl::exponential({.rate = 1.0}), 0.0) == 0.0);

int main()
{
    wsl::xoshiro256pp g(20261118);

    // ---- Quantile de survie --------------------------------------------------------
    check_cquantile(wsl::exponential({.rate = 2.0}));
    check_cquantile(wsl::normal({.mean = 1.0, .sd = 2.0}));
    assert(cquantile(wsl::uniform({.low = -1.0, .high = 3.0}), 0.75) == quantile(wsl::uniform({.low = -1.0, .high = 3.0}), 0.25));   // (queue : x arrondi près de b)
    check_cquantile(wsl::gamma({.shape = 2.5, .rate = 2.0}));
    check_cquantile(wsl::gamma({.shape = 1e-3, .scale = 1.0}), 1e-10);
    check_cquantile(wsl::beta({.alpha = 2.0, .beta = 0.5}), 1e-6);           // x près de 1 : 1 − x mal résolu
    check_cquantile(wsl::chi_squared({.degrees_of_freedom = 3.0}));
    check_cquantile(wsl::student_t({.degrees_of_freedom = 4.0}));
    check_cquantile(wsl::lognormal({.mu = 0.0, .sigma = 0.5}));
    check_cquantile(wsl::cauchy({.location = 1.0, .scale = 2.0}));
    check_cquantile(wsl::fisher_f({.d1 = 5.0, .d2 = 10.0}), 1e-10);
    check_cquantile(wsl::weibull({.shape = 1.5, .scale = 2.0}));
    check_cquantile(wsl::laplace({.location = 1.0, .scale = 2.0}));
    check_cquantile(wsl::logistic({.location = 1.0, .scale = 2.0}));
    check_cquantile(wsl::pareto({.scale = 1.0, .shape = 2.5}));
    assert(cquantile(wsl::uniform({.low = 1.0, .high = 3.0}), 0.1) == 2.8);
    assert(cquantile(wsl::exponential({.rate = 1.0}), 1e-300) > 690);     // 1 − q ne vaudrait rien

    // ---- Image affine : identités ---------------------------------------------------
    {
        auto const Y = 2.0 * wsl::normal({.mean = 1.0, .sd = 2.0}) + 3;
        assert(Y.mean() == 5 && Y.sd() == 4);
        auto const U = 2 - 2 * wsl::uniform({.low = 0.0, .high = 1.0});
        assert(U.low() == 0 && U.high() == 2);
        auto const C = -wsl::cauchy({.location = 1.0, .scale = 2.0});
        assert(C.location() == -1 && C.scale() == 2);
        // Familles fermées par changement d'échelle : mêmes bits que la loi de la famille.
        auto const E4 = 4 * wsl::exponential({.rate = 2.0});
        wsl::exponential const E({.rate = 0.5});
        auto const G3 = 3 * wsl::gamma({.shape = 2.5, .scale = 0.5});
        wsl::gamma const G({.shape = 2.5, .scale = 1.5});
        auto const W2 = 2 * wsl::weibull({.shape = 1.5, .scale = 2.0});
        wsl::weibull const W({.shape = 1.5, .scale = 4.0});
        auto const P5 = 0.5 * wsl::pareto({.scale = 2.0, .shape = 2.5});
        wsl::pareto const P({.scale = 1.0, .shape = 2.5});
        for (double x : {0.01, 0.5, 1.0, 3.0, 40.0})
        {
            assert(cdf(E4, x) == cdf(E, x) && pdf(E4, x) == pdf(E, x) && ccdf(G3, x) == ccdf(G, x) && pdf(G3, x) == pdf(G, x));
            assert(cdf(W2, x) == cdf(W, x) && log_pdf(W2, x) == log_pdf(W, x) && ccdf(P5, x) == ccdf(P, x));
        }
        for (double p : {1e-100, 0.25, 0.5, 0.99})
            assert(quantile(E4, p) == quantile(E, p) && cquantile(G3, p) == cquantile(G, p) && quantile(P5, p) == quantile(P, p));
        assert(E4.mean() == 2 && E4.variance() == E.variance() && G3.entropy() == G.entropy() && W2.median() == W.median());
        auto const L2 = 2 * wsl::lognormal({.mu = 0.0, .sigma = 0.5});
        assert(std::abs(L2.median() - 2) < 4e-16 && close(cdf(L2, 3.0), cdf(wsl::lognormal({.mu = std::log(2.0), .sigma = 0.5}), 3.0), 1e-15));
        // aX, a < 0 ou b ≠ 0 : composition.
        auto const S = 1 - 2 * wsl::exponential({.rate = 1.0});              // support (−∞, 1]
        assert(S.support().lower == -INFINITY && S.support().upper == 1 && cdf(S, 1.0) == 1 && cdf(S, -1.0) == ccdf(wsl::exponential({.rate = 1.0}), 1.0));
        assert(quantile(S, 0.25) == 1 - 2 * cquantile(wsl::exponential({.rate = 1.0}), 0.25));
        assert(throws_domain_error([] { (void)(0 * wsl::normal({.mean = 0.0, .sd = 1.0})); }));
        assert(throws_domain_error([] { (void)wsl::affine(wsl::student_t({.degrees_of_freedom = 3.0}), {.scale = 1, .shift = INFINITY}); }));
    }

    // ---- Image affine générique ----------------------------------------------------
    {
        wsl::student_t const T({.degrees_of_freedom = 4.0});
        auto const Z = wsl::affine(T, {.scale = -2, .shift = 5});
        for (double y : {-11.0, 1.0, 3.0, 5.0, 9.0, 1e6})                    // (y − 5)/(−2) exact
        {
            double const z = (y - 5) / -2;
            assert(cdf(Z, y) == ccdf(T, z) && ccdf(Z, y) == cdf(T, z) && pdf(Z, y) == pdf(T, z) / 2);
        }
        assert(Z.mean() == 5 && Z.variance() == 4 * T.variance() && Z.median() == 5 && Z.entropy() == T.entropy() + std::log(2.0));
        for (double p : {1e-12, 0.1, 0.5, 0.9})
            assert(quantile(Z, p) == 5 - 2 * cquantile(T, p) && cquantile(Z, p) == 5 - 2 * quantile(T, p));
        assert(kolmogorov_smirnov(Z, g) < 1.95);
        auto const B = 10 * wsl::beta({.alpha = 2.0, .beta = 3.0}) - 4;
        assert(B.support().lower == -4 && B.support().upper == 6 && std::abs(B.mean() - 0) < 1e-15 && B.skewness() == wsl::beta({.alpha = 2.0, .beta = 3.0}).skewness());
        assert(kolmogorov_smirnov(B, g) < 1.95);
        auto const N = -wsl::gamma({.shape = 3.0, .rate = 1.0});
        assert(N.skewness() == -wsl::gamma({.shape = 3.0, .rate = 1.0}).skewness() && kolmogorov_smirnov(N, g) < 1.95);
        auto const F = 2.0 * wsl::student_t<float>({.degrees_of_freedom = 3.0f}) + 1.0f;   // float
        static_assert(std::is_same_v<decltype(cdf(F, 1.0f)), float>);
        assert(cdf(F, 1.0f) == 0.5f);
    }

    // ---- Troncature ----------------------------------------------------------------
    {
        wsl::normal const N({.mean = 0.0, .sd = 1.0});
        auto const T = wsl::truncate(N, 10.0, INFINITY);                     // queue droite : par la survie
        assert(T.mass() == ccdf(N, 10.0) && T.support().lower == 10 && cdf(T, 9.0) == 0 && pdf(T, 9.0) == 0);
        assert(close(pdf(T, 10.0), pdf(N, 10.0) / ccdf(N, 10.0), 1e-15) && close(ccdf(T, 10.5), ccdf(N, 10.5) / ccdf(N, 10.0), 1e-13));
        // Aller-retour à 10⁻¹⁵ près en absolu : près d'une borne, (S(a) − S(x))/m se compense.
        for (double p : {1e-8, 0.3, 0.5, 0.999})
            assert(std::abs(cdf(T, quantile(T, p)) - p) < 1e-14 && std::abs(ccdf(T, cquantile(T, p)) - p) < 1e-14);
        assert(quantile(T, 0.0) == 10 && quantile(T, 1.0) == INFINITY);
        double s = 0;                                                        // E[N | N >= a] = φ(a)/Φ̄(a)
        for (int i = 0; i < 20000; ++i)
            s += sample(T, g);
        assert(std::abs(s / 20000 - pdf(N, 10.0) / ccdf(N, 10.0)) < 0.005);
        assert(kolmogorov_smirnov(T, g) < 1.95);
        auto const G = wsl::truncate(wsl::gamma({.shape = 2.0, .rate = 1.0}), 1.0, 3.0);   // masse 0,54 : rejet
        assert(close(G.mass(), cdf(wsl::gamma({.shape = 2.0, .rate = 1.0}), 3.0) - cdf(wsl::gamma({.shape = 2.0, .rate = 1.0}), 1.0), 1e-15));
        assert(kolmogorov_smirnov(G, g) < 1.95 && close(cdf(G, G.median()), 0.5, 1e-14));
        auto const L = wsl::truncate(wsl::cauchy({.location = 0.0, .scale = 1.0}), -INFINITY, -1e3);   // queue gauche
        assert(kolmogorov_smirnov(L, g) < 1.95 && L.support().upper == -1e3);
        auto const Z = wsl::truncate(2.0 * wsl::student_t({.degrees_of_freedom = 2.0}) + 1, -1.0, 0.5);   // tronquer une image
        assert(kolmogorov_smirnov(Z, g) < 1.95);
        assert(throws_domain_error([&] { (void)wsl::truncate(N, 40.0, 50.0); }));   // masse sous les sous-normaux
        assert(throws_domain_error([&] { (void)wsl::truncate(N, 1.0, 1.0); }));
        assert(throws_domain_error([] { (void)wsl::truncate(wsl::exponential({.rate = 1.0}), -5.0, -1.0); }));
    }

    // ---- Mélanges ------------------------------------------------------------------
    {
        wsl::normal const A({.mean = -2.0, .sd = 1.0});
        wsl::gamma const B({.shape = 3.0, .rate = 1.0});
        wsl::mixture const M(wsl::weighted(0.3, A), wsl::weighted(0.7, B));
        wsl::mixture const M10(wsl::weighted(3.0, A), wsl::weighted(7.0, B));    // poids non normalisés
        for (double x : {-3.0, 0.0, 0.5, 2.0, 10.0})
        {
            double const f = 0.3 * pdf(A, x) + 0.7 * pdf(B, x);
            assert(close(pdf(M, x), f, 1e-15) && close(pdf(M10, x), f, 1e-15) && close(std::exp(log_pdf(M, x)), f, 1e-14));
            assert(close(cdf(M, x) + ccdf(M, x), 1, 1e-15));
        }
        assert(close(M.mean(), 1.5, 1e-15) && close(M.variance(), 7.65, 1e-15) && std::abs(M.skewness()) < 1e-15);   // μ₃ = 0 ici
        for (double p : {1e-9, 0.1, 0.5, 0.9, 1 - 1e-9})
            assert(close(cdf(M, quantile(M, p)), p, 1e-13) && close(ccdf(M, cquantile(M, 1 - p)), 1 - p, 1e-13));
        assert(close(ccdf(M, cquantile(M, 1e-30)), 1e-30, 1e-12));
        assert(M.support().lower == -INFINITY && M.support().upper == INFINITY);
        assert(kolmogorov_smirnov(M, g) < 1.95);
        // Une composante : la loi elle-même.
        wsl::mixture const one(wsl::weighted(0.7, B));
        for (double x : {0.5, 2.0, 9.0})
            assert(pdf(one, x) == pdf(B, x) && cdf(one, x) == cdf(B, x));
        // Choix exact des composantes : poids 1, 1, 2 ; χ² à 2 ddl.
        wsl::mixture const C(wsl::weighted(1.0, wsl::uniform({.low = 0.0, .high = 1.0})), wsl::weighted(1.0, wsl::uniform({.low = 1.0, .high = 2.0})),
                             wsl::weighted(2.0, wsl::uniform({.low = 2.0, .high = 3.0})));
        int counts[3] = {};
        for (int i = 0; i < 40000; ++i)
            ++counts[static_cast<int>(sample(C, g))];
        double chi = 0;
        for (int i = 0; i < 3; ++i)
        {
            double const e = 40000.0 * (i == 2 ? 0.5 : 0.25);
            chi += (counts[i] - e) * (counts[i] - e) / e;
        }
        assert(chi < 13.8);                                                  // seuil de 0,1 %
        // Poids d'étendue immense : grands entiers ; la composante de poids 2^-1000 n'apparaît pas.
        wsl::mixture const H(wsl::weighted(0x1p1000, wsl::uniform({.low = 0.0, .high = 1.0})), wsl::weighted(0x1p-1000, wsl::uniform({.low = 5.0, .high = 6.0})),
                             wsl::weighted(0x1p999, wsl::uniform({.low = 1.0, .high = 2.0})));
        int above = 0;
        for (int i = 0; i < 30000; ++i)
        {
            double const x = sample(H, g);
            assert(x < 2);
            above += x >= 1;
        }
        assert(std::abs(above / 30000.0 - 1.0 / 3) < 0.012);
        // Discret : masse, répartition, quantile exact de la répartition du mélange.
        wsl::mixture const D(wsl::weighted(1.0, wsl::poisson({.mean = 2.0})), wsl::weighted(3.0, wsl::poisson({.mean = 20.0})));
        assert(close(pmf(D, 2), 0.25 * pmf(wsl::poisson({.mean = 2.0}), 2) + 0.75 * pmf(wsl::poisson({.mean = 20.0}), 2), 1e-15));
        for (double p : {0.01, 0.2, 0.5, 0.99})
        {
            auto const k = quantile(D, p);
            assert(cdf(D, k) >= p && (k == 0 || cdf(D, k - 1) < p));
        }
        assert(D.mean() == 15.5 && D.variance() == 76.25);
        assert(throws_domain_error([&] { (void)wsl::mixture(wsl::weighted(0.0, A), wsl::weighted(1.0, B)); }));
        assert(throws_domain_error([&] { (void)wsl::mixture(wsl::weighted(-1.0, A)); }));
    }

    // ---- Empreintes ----------------------------------------------------------------
    bool ok = true;
    {
        wsl::xoshiro256pp x(42);
        auto const Z = wsl::affine(wsl::student_t({.degrees_of_freedom = 5.0}), {.scale = -1.5, .shift = 2});
        auto const T = wsl::truncate(wsl::normal({.mean = 0.0, .sd = 1.0}), 3.0, INFINITY);     // inversion
        auto const R = wsl::truncate(wsl::laplace({.location = 0.0, .scale = 1.0}), -1.0, 2.0); // rejet, tirage de X exact
        wsl::mixture const M(wsl::weighted(0.3, wsl::normal({.mean = -2.0, .sd = 1.0})), wsl::weighted(0.7, wsl::gamma({.shape = 3.0, .rate = 1.0})));
        wsl::mixture const H(wsl::weighted(0x1p600, wsl::exponential({.rate = 1.0})), wsl::weighted(0x1p-600, wsl::normal({.mean = 0.0, .sd = 1.0})),
                             wsl::weighted(0x1p599, wsl::laplace({.location = 0.0, .scale = 1.0})));
        wsl::mixture const D(wsl::weighted(1.0, wsl::poisson({.mean = 2.0})), wsl::weighted(3.0, wsl::binomial({.trials = 40, .p = 0.5})));
        ok &= check("affine", fingerprint([&] { return bits(sample(Z, x)); }), 0x44509b09b529826c);
        ok &= check("truncated (inv.)", fingerprint([&] { return bits(sample(T, x)); }), 0x398b1148ffa1462f);
        ok &= check("truncated (rejet)", fingerprint([&] { return bits(sample(R, x)); }), 0x550c9bd3098ec373);
        ok &= check("mixture", fingerprint([&] { return bits(sample(M, x)); }), 0x1cd22859dc5c692f);
        ok &= check("mixture (grands)", fingerprint([&] { return bits(sample(H, x)); }), 0xa71e5cac21b536d6);
        ok &= check("mixture discrète", fingerprint([&] { return static_cast<std::uint64_t>(sample(D, x)); }), 0x2c54e55c5ac36aab);
        ok &= check("algebra functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::uniform_real(x), t = fp::fmul(20.0, u) - 10;
            return bits(cdf(Z, t)) ^ (bits(quantile(Z, u)) << 1) ^ (bits(cdf(T, t + 13)) << 2) ^ (bits(quantile(T, u)) << 3)
                 ^ (bits(cdf(M, t)) << 4) ^ (bits(quantile(M, u)) << 5) ^ (bits(cquantile(wsl::gamma({.shape = 0.5, .rate = 1.0}), u)) << 6);
        }, 1 << 9), 0xb06389a42d4e1276);
    }
    assert(ok);
    std::printf("test_random_algebra : ok\n");
}
