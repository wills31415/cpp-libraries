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
#include <tuple>
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
            v.push_back(double(wsl::random::sample(law, g)));
        std::ranges::sort(v);
        double worst = 0;
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            double const f = wsl::random::cdf(law, v[i]);
            worst = std::max({worst, std::abs(f - double(i) / n), std::abs(f - double(i + 1) / n)});
        }
        return worst * std::sqrt(double(n));
    }

    bool close(double a, double b, double relative) { return std::abs(a - b) <= relative * std::abs(b); }

    // wsl::random::ccdf(X, wsl::random::cquantile(X, q)) ≈ q, et wsl::random::cquantile(X, q) = wsl::random::quantile(X, 1 − q) pour q >= ½.
    template <class L>
    void check_cquantile(L const& X, double tolerance = 1e-12)
    {
        for (double q : {0.5, 0.75, 0.9375, 1 - 0x1p-40})
            assert(wsl::random::cquantile(X, q) == wsl::random::quantile(X, 1 - q));
        for (double q : {0.25, 1e-5, 1e-30, 1e-200})
        {
            double const x = wsl::random::cquantile(X, q);
            if (std::isfinite(x) && x != X.support().lower && x != X.support().upper && !close(double(wsl::random::ccdf(X, x)), q, tolerance))
            {
                std::printf("wsl::random::cquantile(%g) = %a, ccdf = %.17g\n", q, x, double(wsl::random::ccdf(X, x)));
                assert(false);
            }
        }
        assert(wsl::random::cquantile(X, 1.0) == wsl::random::quantile(X, 0.0) && wsl::random::cquantile(X, 0.0) == wsl::random::quantile(X, 1.0));
        assert(throws_domain_error([&] { (void)wsl::random::cquantile(X, 1.5); }));
    }
}

// Identités : les familles fermées gardent leur type.
using normal_t = wsl::random::normal<double>;
static_assert(std::is_same_v<decltype(wsl::random::affine(std::declval<normal_t>(), {.scale = 2.0, .shift = 1})), normal_t>);
static_assert(std::is_same_v<decltype(wsl::random::affine(std::declval<wsl::random::cauchy<>>(), {.scale = -1})), wsl::random::cauchy<>>);
static_assert(std::is_same_v<decltype(wsl::random::affine(std::declval<wsl::random::uniform<float>>(), {.scale = -1, .shift = 3})), wsl::random::uniform<float>>);
static_assert(std::is_same_v<decltype(wsl::random::affine(std::declval<wsl::random::laplace<>>(), {.scale = 0.25})), wsl::random::laplace<>>);
static_assert(std::is_same_v<decltype(wsl::random::affine(std::declval<wsl::random::exponential<>>(), {.scale = 2})), wsl::random::affine_image<wsl::random::exponential<>>>);
static_assert(std::is_same_v<decltype(wsl::random::affine(wsl::random::affine(std::declval<wsl::random::student_t<>>(), {.scale = 2}), {.scale = 2, .shift = 1})), wsl::random::affine_image<wsl::random::student_t<>>>);
static_assert(wsl::random::affine_closed<wsl::random::logistic<>> && wsl::random::scale_closed<wsl::random::gamma<>> && !wsl::random::scale_closed<wsl::random::beta<>>);
static_assert(!wsl::random::has_mean<decltype(wsl::random::affine(std::declval<wsl::random::cauchy<>>(), {.scale = 2.0}))>);
static_assert(wsl::random::has_mean<wsl::random::affine_image<wsl::random::student_t<>>> && wsl::random::has_cquantile<wsl::random::affine_image<wsl::random::beta<>>>);
static_assert(!wsl::random::has_mean<wsl::random::truncated<normal_t>> && wsl::random::has_cquantile<wsl::random::truncated<normal_t>>);
static_assert(wsl::random::probability_law<wsl::random::mixture<normal_t, wsl::random::gamma<>>> && wsl::random::has_pmf<wsl::random::mixture<wsl::random::poisson<>, wsl::random::binomial<>>>);
static_assert(wsl::random::mixture<normal_t, wsl::random::exponential<>>::sampling == wsl::random::exactness::exact);
static_assert(wsl::random::mixture<normal_t, wsl::random::gamma<>>::sampling == wsl::random::exactness::approximate);
static_assert(!wsl::random::has_pdf<wsl::random::mixture<wsl::random::poisson<>>> && wsl::random::has_cdf<wsl::random::mixture<wsl::random::poisson<>>>);
// Évaluables à la compilation.
static_assert(wsl::random::affine(wsl::random::normal({.mean = 1.0, .stddev = 2.0}), {.scale = 3.0, .shift = -1.0}).mean() == 2.0);
static_assert(wsl::random::cdf(wsl::random::affine(wsl::random::exponential({.rate = 1.0}), {.scale = 2.0}), 0.0) == 0.0);

int main()
{
    wsl::random::xoshiro256pp g(20261118);

    // ---- Quantile de survie --------------------------------------------------------
    check_cquantile(wsl::random::exponential({.rate = 2.0}));
    check_cquantile(wsl::random::normal({.mean = 1.0, .stddev = 2.0}));
    assert(wsl::random::cquantile(wsl::random::uniform({.lower = -1.0, .upper = 3.0}), 0.75) == wsl::random::quantile(wsl::random::uniform({.lower = -1.0, .upper = 3.0}), 0.25));   // (queue : x arrondi près de b)
    check_cquantile(wsl::random::gamma({.rate = 2.0, .shape = 2.5}));
    check_cquantile(wsl::random::gamma({.scale = 1.0, .shape = 1e-3}), 1e-10);
    check_cquantile(wsl::random::beta({.alpha = 2.0, .beta = 0.5}), 1e-6);           // x près de 1 : 1 − x mal résolu
    check_cquantile(wsl::random::chi_squared({.degrees_of_freedom = 3.0}));
    check_cquantile(wsl::random::student_t({.degrees_of_freedom = 4.0}));
    check_cquantile(wsl::random::lognormal({.log_mean = 0.0, .log_stddev = 0.5}));
    check_cquantile(wsl::random::cauchy({.location = 1.0, .scale = 2.0}));
    check_cquantile(wsl::random::fisher_f({.d1 = 5.0, .d2 = 10.0}), 1e-10);
    check_cquantile(wsl::random::weibull({.scale = 2.0, .shape = 1.5}));
    check_cquantile(wsl::random::laplace({.location = 1.0, .scale = 2.0}));
    check_cquantile(wsl::random::logistic({.location = 1.0, .scale = 2.0}));
    check_cquantile(wsl::random::pareto({.scale = 1.0, .shape = 2.5}));
    assert(wsl::random::cquantile(wsl::random::uniform({.lower = 1.0, .upper = 3.0}), 0.1) == 2.8);
    assert(wsl::random::cquantile(wsl::random::exponential({.rate = 1.0}), 1e-300) > 690);     // 1 − q ne vaudrait rien

    // ---- Image affine : identités ---------------------------------------------------
    {
        auto const Y = wsl::random::affine(wsl::random::normal({.mean = 1.0, .stddev = 2.0}), {.scale = 2.0, .shift = 3});
        assert(Y.mean() == 5 && Y.stddev() == 4);
        auto const U = wsl::random::affine(wsl::random::uniform({.lower = 0.0, .upper = 1.0}), {.scale = -2, .shift = 2});
        assert(U.lower() == 0 && U.upper() == 2);
        auto const C = wsl::random::affine(wsl::random::cauchy({.location = 1.0, .scale = 2.0}), {.scale = -1});
        assert(C.location() == -1 && C.scale() == 2);
        // Familles fermées par changement d'échelle : mêmes bits que la loi de la famille.
        auto const E4 = wsl::random::affine(wsl::random::exponential({.rate = 2.0}), {.scale = 4});
        wsl::random::exponential const E({.rate = 0.5});
        auto const G3 = wsl::random::affine(wsl::random::gamma({.scale = 0.5, .shape = 2.5}), {.scale = 3});
        wsl::random::gamma const G({.scale = 1.5, .shape = 2.5});
        auto const W2 = wsl::random::affine(wsl::random::weibull({.scale = 2.0, .shape = 1.5}), {.scale = 2});
        wsl::random::weibull const W({.scale = 4.0, .shape = 1.5});
        auto const P5 = wsl::random::affine(wsl::random::pareto({.scale = 2.0, .shape = 2.5}), {.scale = 0.5});
        wsl::random::pareto const P({.scale = 1.0, .shape = 2.5});
        for (double x : {0.01, 0.5, 1.0, 3.0, 40.0})
        {
            assert(wsl::random::cdf(E4, x) == wsl::random::cdf(E, x) && wsl::random::pdf(E4, x) == wsl::random::pdf(E, x) && wsl::random::ccdf(G3, x) == wsl::random::ccdf(G, x) && wsl::random::pdf(G3, x) == wsl::random::pdf(G, x));
            assert(wsl::random::cdf(W2, x) == wsl::random::cdf(W, x) && wsl::random::log_pdf(W2, x) == wsl::random::log_pdf(W, x) && wsl::random::ccdf(P5, x) == wsl::random::ccdf(P, x));
        }
        for (double p : {1e-100, 0.25, 0.5, 0.99})
            assert(wsl::random::quantile(E4, p) == wsl::random::quantile(E, p) && wsl::random::cquantile(G3, p) == wsl::random::cquantile(G, p) && wsl::random::quantile(P5, p) == wsl::random::quantile(P, p));
        assert(E4.mean() == 2 && E4.variance() == E.variance() && G3.entropy() == G.entropy() && W2.median() == W.median());
        auto const L2 = wsl::random::affine(wsl::random::lognormal({.log_mean = 0.0, .log_stddev = 0.5}), {.scale = 2});
        assert(std::abs(L2.median() - 2) < 4e-16 && close(wsl::random::cdf(L2, 3.0), wsl::random::cdf(wsl::random::lognormal({.log_mean = std::log(2.0), .log_stddev = 0.5}), 3.0), 1e-15));
        // aX, a < 0 ou b ≠ 0 : composition.
        auto const S = wsl::random::affine(wsl::random::exponential({.rate = 1.0}), {.scale = -2, .shift = 1});              // support (−∞, 1]
        assert(S.support().lower == -INFINITY && S.support().upper == 1 && wsl::random::cdf(S, 1.0) == 1 && wsl::random::cdf(S, -1.0) == wsl::random::ccdf(wsl::random::exponential({.rate = 1.0}), 1.0));
        assert(wsl::random::quantile(S, 0.25) == 1 - 2 * wsl::random::cquantile(wsl::random::exponential({.rate = 1.0}), 0.25));
        assert(throws_domain_error([] { (void)wsl::random::affine(wsl::random::normal({.mean = 0.0, .stddev = 1.0}), {.scale = 0}); }));
        assert(throws_domain_error([] { (void)wsl::random::affine(wsl::random::student_t({.degrees_of_freedom = 3.0}), {.scale = 1, .shift = INFINITY}); }));
    }

    // ---- Image affine générique ----------------------------------------------------
    {
        wsl::random::student_t const T({.degrees_of_freedom = 4.0});
        auto const Z = wsl::random::affine(T, {.scale = -2, .shift = 5});
        for (double y : {-11.0, 1.0, 3.0, 5.0, 9.0, 1e6})                    // (y − 5)/(−2) exact
        {
            double const z = (y - 5) / -2;
            assert(wsl::random::cdf(Z, y) == wsl::random::ccdf(T, z) && wsl::random::ccdf(Z, y) == wsl::random::cdf(T, z) && wsl::random::pdf(Z, y) == wsl::random::pdf(T, z) / 2);
        }
        assert(Z.mean() == 5 && Z.variance() == 4 * T.variance() && Z.median() == 5 && Z.entropy() == T.entropy() + std::log(2.0));
        for (double p : {1e-12, 0.1, 0.5, 0.9})
            assert(wsl::random::quantile(Z, p) == 5 - 2 * wsl::random::cquantile(T, p) && wsl::random::cquantile(Z, p) == 5 - 2 * wsl::random::quantile(T, p));
        assert(kolmogorov_smirnov(Z, g) < 1.95);
        auto const B = wsl::random::affine(wsl::random::beta({.alpha = 2.0, .beta = 3.0}), {.scale = 10, .shift = -4});
        assert(B.support().lower == -4 && B.support().upper == 6 && std::abs(B.mean() - 0) < 1e-15 && B.skewness() == wsl::random::beta({.alpha = 2.0, .beta = 3.0}).skewness());
        assert(kolmogorov_smirnov(B, g) < 1.95);
        auto const N = wsl::random::affine(wsl::random::gamma({.rate = 1.0, .shape = 3.0}), {.scale = -1});
        assert(N.skewness() == -wsl::random::gamma({.rate = 1.0, .shape = 3.0}).skewness() && kolmogorov_smirnov(N, g) < 1.95);
        auto const F = wsl::random::affine(wsl::random::student_t<float>({.degrees_of_freedom = 3.0f}), {.scale = 2.0, .shift = 1.0});   // float
        static_assert(std::is_same_v<decltype(wsl::random::cdf(F, 1.0f)), float>);
        assert(wsl::random::cdf(F, 1.0f) == 0.5f);
    }

    // ---- Troncature ----------------------------------------------------------------
    {
        wsl::random::normal const N({.mean = 0.0, .stddev = 1.0});
        auto const T = wsl::random::truncate(N, {.lower = 10.0});                     // queue droite : par la survie
        assert(T.mass() == wsl::random::ccdf(N, 10.0) && T.support().lower == 10 && wsl::random::cdf(T, 9.0) == 0 && wsl::random::pdf(T, 9.0) == 0);
        assert(close(wsl::random::pdf(T, 10.0), wsl::random::pdf(N, 10.0) / wsl::random::ccdf(N, 10.0), 1e-15) && close(wsl::random::ccdf(T, 10.5), wsl::random::ccdf(N, 10.5) / wsl::random::ccdf(N, 10.0), 1e-13));
        // Aller-retour à 10⁻¹⁵ près en absolu : près d'une borne, (S(a) − S(x))/m se compense.
        for (double p : {1e-8, 0.3, 0.5, 0.999})
            assert(std::abs(wsl::random::cdf(T, wsl::random::quantile(T, p)) - p) < 1e-14 && std::abs(wsl::random::ccdf(T, wsl::random::cquantile(T, p)) - p) < 1e-14);
        assert(wsl::random::quantile(T, 0.0) == 10 && wsl::random::quantile(T, 1.0) == INFINITY);
        double s = 0;                                                        // E[N | N >= a] = φ(a)/Φ̄(a)
        for (int i = 0; i < 20000; ++i)
            s += wsl::random::sample(T, g);
        assert(std::abs(s / 20000 - wsl::random::pdf(N, 10.0) / wsl::random::ccdf(N, 10.0)) < 0.005);
        assert(kolmogorov_smirnov(T, g) < 1.95);
        auto const G = wsl::random::truncate(wsl::random::gamma({.rate = 1.0, .shape = 2.0}), {.lower = 1.0, .upper = 3.0});   // masse 0,54 : rejet
        assert(close(G.mass(), wsl::random::cdf(wsl::random::gamma({.rate = 1.0, .shape = 2.0}), 3.0) - wsl::random::cdf(wsl::random::gamma({.rate = 1.0, .shape = 2.0}), 1.0), 1e-15));
        assert(kolmogorov_smirnov(G, g) < 1.95 && close(wsl::random::cdf(G, G.median()), 0.5, 1e-14));
        auto const L = wsl::random::truncate(wsl::random::cauchy({.location = 0.0, .scale = 1.0}), {.upper = -1e3});   // queue gauche
        assert(kolmogorov_smirnov(L, g) < 1.95 && L.support().upper == -1e3);
        auto const Z = wsl::random::truncate(wsl::random::affine(wsl::random::student_t({.degrees_of_freedom = 2.0}), {.scale = 2.0, .shift = 1}), {.lower = -1.0, .upper = 0.5});   // tronquer une image
        assert(kolmogorov_smirnov(Z, g) < 1.95);
        assert(throws_domain_error([&] { (void)wsl::random::truncate(N, {.lower = 40.0, .upper = 50.0}); }));   // masse sous les sous-normaux
        assert(throws_domain_error([&] { (void)wsl::random::truncate(N, {.lower = 1.0, .upper = 1.0}); }));
        assert(throws_domain_error([] { (void)wsl::random::truncate(wsl::random::exponential({.rate = 1.0}), {.lower = -5.0, .upper = -1.0}); }));
    }

    // ---- Mélanges ------------------------------------------------------------------
    {
        wsl::random::normal const A({.mean = -2.0, .stddev = 1.0});
        wsl::random::gamma const B({.rate = 1.0, .shape = 3.0});
        wsl::random::mixture const M(wsl::random::weighted(0.3, A), wsl::random::weighted(0.7, B));
        wsl::random::mixture const M10(wsl::random::weighted(3.0, A), wsl::random::weighted(7.0, B));    // poids non normalisés
        for (double x : {-3.0, 0.0, 0.5, 2.0, 10.0})
        {
            double const f = 0.3 * wsl::random::pdf(A, x) + 0.7 * wsl::random::pdf(B, x);
            assert(close(wsl::random::pdf(M, x), f, 1e-15) && close(wsl::random::pdf(M10, x), f, 1e-15) && close(std::exp(wsl::random::log_pdf(M, x)), f, 1e-14));
            assert(close(wsl::random::cdf(M, x) + wsl::random::ccdf(M, x), 1, 1e-15));
        }
        assert(close(M.mean(), 1.5, 1e-15) && close(M.variance(), 7.65, 1e-15) && std::abs(M.skewness()) < 1e-15);   // μ₃ = 0 ici
        for (double p : {1e-9, 0.1, 0.5, 0.9, 1 - 1e-9})
            assert(close(wsl::random::cdf(M, wsl::random::quantile(M, p)), p, 1e-13) && close(wsl::random::ccdf(M, wsl::random::cquantile(M, 1 - p)), 1 - p, 1e-13));
        assert(close(wsl::random::ccdf(M, wsl::random::cquantile(M, 1e-30)), 1e-30, 1e-12));
        assert(M.support().lower == -INFINITY && M.support().upper == INFINITY);
        assert(kolmogorov_smirnov(M, g) < 1.95);
        // Une composante : la loi elle-même.
        wsl::random::mixture const one(wsl::random::weighted(0.7, B));
        for (double x : {0.5, 2.0, 9.0})
            assert(wsl::random::pdf(one, x) == wsl::random::pdf(B, x) && wsl::random::cdf(one, x) == wsl::random::cdf(B, x));
        // Choix exact des composantes : poids 1, 1, 2 ; χ² à 2 ddl.
        wsl::random::mixture const C(wsl::random::weighted(1.0, wsl::random::uniform({.lower = 0.0, .upper = 1.0})), wsl::random::weighted(1.0, wsl::random::uniform({.lower = 1.0, .upper = 2.0})),
                             wsl::random::weighted(2.0, wsl::random::uniform({.lower = 2.0, .upper = 3.0})));
        int counts[3] = {};
        for (int i = 0; i < 40000; ++i)
            ++counts[static_cast<int>(wsl::random::sample(C, g))];
        double chi = 0;
        for (int i = 0; i < 3; ++i)
        {
            double const e = 40000.0 * (i == 2 ? 0.5 : 0.25);
            chi += (counts[i] - e) * (counts[i] - e) / e;
        }
        assert(chi < 13.8);                                                  // seuil de 0,1 %
        // Poids d'étendue immense : grands entiers ; la composante de poids 2^-1000 n'apparaît pas.
        wsl::random::mixture const H(wsl::random::weighted(0x1p1000, wsl::random::uniform({.lower = 0.0, .upper = 1.0})), wsl::random::weighted(0x1p-1000, wsl::random::uniform({.lower = 5.0, .upper = 6.0})),
                             wsl::random::weighted(0x1p999, wsl::random::uniform({.lower = 1.0, .upper = 2.0})));
        int above = 0;
        for (int i = 0; i < 30000; ++i)
        {
            double const x = wsl::random::sample(H, g);
            assert(x < 2);
            above += x >= 1;
        }
        assert(std::abs(above / 30000.0 - 1.0 / 3) < 0.012);
        // Discret : masse, répartition, quantile exact de la répartition du mélange.
        wsl::random::mixture const D(wsl::random::weighted(1.0, wsl::random::poisson({.mean = 2.0})), wsl::random::weighted(3.0, wsl::random::poisson({.mean = 20.0})));
        assert(close(wsl::random::pmf(D, 2), 0.25 * wsl::random::pmf(wsl::random::poisson({.mean = 2.0}), 2) + 0.75 * wsl::random::pmf(wsl::random::poisson({.mean = 20.0}), 2), 1e-15));
        for (double p : {0.01, 0.2, 0.5, 0.99})
        {
            auto const k = wsl::random::quantile(D, p);
            assert(wsl::random::cdf(D, k) >= p && (k == 0 || wsl::random::cdf(D, k - 1) < p));
        }
        assert(D.mean() == 15.5 && D.variance() == 76.25);
        assert(throws_domain_error([&] { (void)wsl::random::mixture(wsl::random::weighted(0.0, A), wsl::random::weighted(1.0, B)); }));
        assert(throws_domain_error([&] { (void)wsl::random::mixture(wsl::random::weighted(-1.0, A)); }));
    }

    // ---- Sommes indépendantes ------------------------------------------------------
    {
        static_assert(std::is_same_v<decltype(wsl::random::convolution(wsl::random::normal({.mean = 1.0, .stddev = 3.0}), wsl::random::normal({.mean = 2.0, .stddev = 4.0}))), wsl::random::normal<>>);
        static_assert(std::is_same_v<decltype(wsl::random::convolution(wsl::random::poisson({.mean = 1.0}), wsl::random::poisson({.mean = 2.0}))), wsl::random::poisson<>>);
        static_assert(std::is_same_v<decltype(wsl::random::convolution(wsl::random::gamma({.rate = 1.0, .shape = 1.0}), wsl::random::gamma({.rate = 1.0, .shape = 1.0}))), wsl::random::convolution_law<wsl::random::gamma<>, wsl::random::gamma<>>>);
        static_assert(wsl::random::has_cdf<wsl::random::convolution_law<wsl::random::gamma<>, wsl::random::exponential<>>> && !wsl::random::has_cdf<wsl::random::convolution_law<wsl::random::student_t<>, wsl::random::uniform<>>>);
        static_assert(wsl::random::convolution_law<wsl::random::poisson<>, wsl::random::geometric<>>::sampling == wsl::random::exactness::approximate);
        static_assert(wsl::random::convolution_law<wsl::random::discrete_uniform<>, wsl::random::discrete_uniform<>>::sampling == wsl::random::exactness::exact);
        auto const S = wsl::random::convolution(wsl::random::normal({.mean = 1.0, .stddev = 3.0}), wsl::random::normal({.mean = 2.0, .stddev = 4.0}));
        assert(S.mean() == 3 && S.stddev() == 5);
        auto const D = wsl::random::convolution(wsl::random::normal({.mean = 1.0, .stddev = 3.0}), wsl::random::affine(wsl::random::normal({.mean = 2.0, .stddev = 4.0}), {.scale = -1}));
        assert(D.mean() == -1 && D.stddev() == 5);
        auto const C = wsl::random::convolution(wsl::random::cauchy({.location = 1.0, .scale = 2.0}), wsl::random::cauchy({.location = -3.0, .scale = 0.5}));
        assert(C.location() == -2 && C.scale() == 2.5);
        assert(wsl::random::convolution(wsl::random::chi_squared({.degrees_of_freedom = 3.0}), wsl::random::chi_squared({.degrees_of_freedom = 4.5})).degrees_of_freedom() == 7.5);
        // Identités conditionnelles : mêmes bits que la loi de la famille.
        auto const G = wsl::random::convolution(wsl::random::gamma({.rate = 2.0, .shape = 2.0}), wsl::random::gamma({.scale = 0.5, .shape = 1.5}));   // paramétrisations mêlées
        wsl::random::gamma const G35({.rate = 2.0, .shape = 3.5});
        auto const E = wsl::random::convolution(wsl::random::exponential({.rate = 2.0}), wsl::random::gamma({.rate = 2.0, .shape = 2.5}));
        auto const B = wsl::random::convolution(wsl::random::binomial({.trials = 10, .p = 0.3}), wsl::random::binomial({.trials = 5, .p = 0.3}));
        wsl::random::binomial const B15({.trials = 15, .p = 0.3});
        assert(G.family() && E.family() && B.family());
        for (double x : {0.1, 1.0, 4.0})
            assert(wsl::random::cdf(G, x) == wsl::random::cdf(G35, x) && wsl::random::pdf(G, x) == wsl::random::pdf(G35, x) && wsl::random::ccdf(E, x) == wsl::random::ccdf(G35, x));
        assert(wsl::random::quantile(G, 0.3) == wsl::random::quantile(G35, 0.3) && G.variance() == G35.variance());
        for (int k = 0; k <= 15; ++k)
            assert(wsl::random::pmf(B, k) == wsl::random::pmf(B15, k) && wsl::random::cdf(B, k) == wsl::random::cdf(B15, k));
        // Sans identité : moments par les cumulants ; NaN pour une somme continue.
        auto const H = wsl::random::convolution(wsl::random::gamma({.rate = 2.0, .shape = 2.0}), wsl::random::gamma({.rate = 3.0, .shape = 1.5}));
        assert(!H.family() && std::isnan(wsl::random::cdf(H, 1.0)) && throws_domain_error([&] { (void)wsl::random::quantile(H, 0.5); }));
        assert(H.mean() == 1.5 && close(H.variance(), 0.5 + 1.5 / 9, 1e-15));
        double const k3 = 2.0 * 2 / 8 + 2.0 * 1.5 / 27, v = 0.5 + 1.5 / 9;                     // κ₃ = 2a/λ³
        assert(close(H.skewness(), k3 / (v * std::sqrt(v)), 1e-14));
        double m = 0;
        for (int i = 0; i < 40000; ++i)
            m += wsl::random::sample(H, g);
        assert(std::abs(m / 40000 - 1.5) < 0.02);
        // Convolution discrète : masse, répartition, survie, quantile exact.
        wsl::random::poisson const P({.mean = 3.5});
        wsl::random::binomial const Bi({.trials = 20, .p = 0.3});
        auto const PB = wsl::random::convolution(P, Bi);
        double total = 0;
        for (long long k = 0; k <= 60; ++k)
        {
            double direct = 0;
            for (long long j = 0; j <= std::min<long long>(k, 20); ++j)
                direct += wsl::random::pmf(Bi, j) * wsl::random::pmf(P, k - j);
            assert(close(wsl::random::pmf(PB, k), direct, 1e-13));
            assert(close(wsl::random::cdf(PB, k) + wsl::random::ccdf(PB, k), 1, 1e-15));
            total += wsl::random::pmf(PB, k);
        }
        assert(close(total, 1, 1e-14) && PB.mean() == 9.5 && close(PB.variance(), 3.5 + 4.2, 1e-15));
        for (double p : {1e-12, 0.05, 0.5, 0.99, 1 - 1e-10})
        {
            auto const k = wsl::random::quantile(PB, p);
            assert(wsl::random::cdf(PB, k) >= p && (k == 0 || wsl::random::cdf(PB, k - 1) < p));
        }
        assert(close(wsl::random::ccdf(PB, 80), wsl::random::ccdf(P, 60) * 0 + wsl::random::ccdf(PB, 80), 0) && wsl::random::ccdf(PB, 80) > 0 && wsl::random::ccdf(PB, 80) < 1e-30);   // queue lointaine, sans 1 − cdf
        auto const GU = wsl::random::convolution(wsl::random::geometric({.p = 0.2}), wsl::random::discrete_uniform({.lower = -3, .upper = 3}));        // supports négatifs
        assert(GU.support().lower == -3 && close(wsl::random::cdf(GU, 0), [&] { double s = 0; for (int j = -3; j <= 3; ++j) s += wsl::random::cdf(wsl::random::geometric({.p = 0.2}), -j) / 7; return s; }(), 1e-14));
    }

    // ---- Loi jointe et image par une fonction ----------------------------------------
    {
        namespace sp = wsl::special;
        wsl::random::normal const N({.mean = 0.0, .stddev = 1.0});
        wsl::random::exponential const E({.rate = 2.0});
        wsl::random::product_measure const J(N, E, wsl::random::poisson({.mean = 3.0}));
        static_assert(std::is_same_v<decltype(J)::value_type, std::tuple<double, double, std::int64_t>>);
        static_assert(decltype(J)::sampling == wsl::random::exactness::approximate && !wsl::random::has_pdf<decltype(J)>);
        wsl::random::product_measure const K(N, E);
        static_assert(decltype(K)::sampling == wsl::random::exactness::exact);
        assert(close(wsl::random::pdf(K, {0.5, 1.0}), wsl::random::pdf(N, 0.5) * wsl::random::pdf(E, 1.0), 1e-15) && close(wsl::random::cdf(K, {0.0, 1.0}), 0.5 * wsl::random::cdf(E, 1.0), 1e-15));
        assert(close(wsl::random::log_pdf(K, {0.5, 1.0}), wsl::random::log_pdf(N, 0.5) + wsl::random::log_pdf(E, 1.0), 1e-15) && close(K.entropy(), N.entropy() + E.entropy(), 1e-15));
        assert(wsl::random::pdf(K, {0.5, -1.0}) == 0 && wsl::random::log_pdf(K, {0.5, -1.0}) == -INFINITY && std::get<1>(K.mean()) == 0.5);
        wsl::random::xoshiro256pp a(9), b(9);                                       // tirage : composantes dans l'ordre
        auto const [x, y, k] = wsl::random::sample(J, a);
        assert(x == wsl::random::sample(N, b) && y == wsl::random::sample(E, b) && k == wsl::random::sample(wsl::random::poisson({.mean = 3.0}), b));
        // Image monotone : exp(N) est la log-normale ; 1/X, X gamma, la gamma inverse.
        auto const L = wsl::random::pushforward(N, wsl::random::increasing([](double t) { return sp::exp(t); }, [](double t) { return sp::log(t); }, [](double t) { return 1 / t; }));
        wsl::random::lognormal const LN({.log_mean = 0.0, .log_stddev = 1.0});
        for (double t : {0.1, 1.0, 2.0, 30.0})
            assert(close(wsl::random::cdf(L, t), wsl::random::cdf(LN, t), 1e-15) && close(wsl::random::pdf(L, t), wsl::random::pdf(LN, t), 1e-15) && close(wsl::random::log_pdf(L, t), wsl::random::log_pdf(LN, t), 1e-14));
        assert(close(wsl::random::quantile(L, 0.9), wsl::random::quantile(LN, 0.9), 1e-15) && L.median() == 1 && L.support().lower == 0 && kolmogorov_smirnov(L, g) < 1.95);
        auto const I = wsl::random::pushforward(wsl::random::gamma({.rate = 2.0, .shape = 3.0}), wsl::random::decreasing([](double t) { return 1 / t; }, [](double t) { return 1 / t; }, [](double t) { return -1 / (t * t); }));
        assert(close(wsl::random::cdf(I, 1.0), wsl::random::ccdf(wsl::random::gamma({.rate = 2.0, .shape = 3.0}), 1.0), 1e-15) && close(wsl::random::quantile(I, 0.25), 1 / wsl::random::cquantile(wsl::random::gamma({.rate = 2.0, .shape = 3.0}), 0.25), 1e-15));
        assert(kolmogorov_smirnov(I, g) < 1.95 && wsl::random::pdf(I, -1.0) == 0);
        // Fonction quelconque : tirage seul. N² est la loi du χ² à 1 degré.
        auto const Q = wsl::random::pushforward(N, [](double t) { return t * t; });
        static_assert(!wsl::random::has_cdf<decltype(Q)> && wsl::random::probability_law<decltype(Q)>);
        std::vector<double> v;
        for (int i = 0; i < 20000; ++i)
            v.push_back(wsl::random::sample(Q, g));
        std::ranges::sort(v);
        double worst = 0;
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            double const f = wsl::random::cdf(wsl::random::chi_squared({.degrees_of_freedom = 1.0}), v[i]);
            worst = std::max({worst, std::abs(f - double(i) / 20000), std::abs(f - double(i + 1) / 20000)});
        }
        assert(worst * std::sqrt(20000.0) < 1.95);
    }

    // ---- Empreintes ----------------------------------------------------------------
    bool ok = true;
    {
        wsl::random::xoshiro256pp x(42);
        auto const Z = wsl::random::affine(wsl::random::student_t({.degrees_of_freedom = 5.0}), {.scale = -1.5, .shift = 2});
        auto const T = wsl::random::truncate(wsl::random::normal({.mean = 0.0, .stddev = 1.0}), {.lower = 3.0});     // inversion
        auto const R = wsl::random::truncate(wsl::random::laplace({.location = 0.0, .scale = 1.0}), {.lower = -1.0, .upper = 2.0}); // rejet, tirage de X exact
        wsl::random::mixture const M(wsl::random::weighted(0.3, wsl::random::normal({.mean = -2.0, .stddev = 1.0})), wsl::random::weighted(0.7, wsl::random::gamma({.rate = 1.0, .shape = 3.0})));
        wsl::random::mixture const H(wsl::random::weighted(0x1p600, wsl::random::exponential({.rate = 1.0})), wsl::random::weighted(0x1p-600, wsl::random::normal({.mean = 0.0, .stddev = 1.0})),
                             wsl::random::weighted(0x1p599, wsl::random::laplace({.location = 0.0, .scale = 1.0})));
        wsl::random::mixture const D(wsl::random::weighted(1.0, wsl::random::poisson({.mean = 2.0})), wsl::random::weighted(3.0, wsl::random::binomial({.trials = 40, .p = 0.5})));
        ok &= check("affine", fingerprint([&] { return bits(wsl::random::sample(Z, x)); }), 0x44509b09b529826c);
        ok &= check("truncated (inv.)", fingerprint([&] { return bits(wsl::random::sample(T, x)); }), 0x398b1148ffa1462f);
        ok &= check("truncated (rejet)", fingerprint([&] { return bits(wsl::random::sample(R, x)); }), 0x550c9bd3098ec373);
        ok &= check("mixture", fingerprint([&] { return bits(wsl::random::sample(M, x)); }), 0x1cd22859dc5c692f);
        ok &= check("mixture (grands)", fingerprint([&] { return bits(wsl::random::sample(H, x)); }), 0xa71e5cac21b536d6);
        ok &= check("mixture discrète", fingerprint([&] { return static_cast<std::uint64_t>(wsl::random::sample(D, x)); }), 0x2c54e55c5ac36aab);
        auto const PB = wsl::random::convolution(wsl::random::poisson({.mean = 3.5}), wsl::random::binomial({.trials = 20, .p = 0.3}));
        auto const ST = wsl::random::convolution(wsl::random::student_t({.degrees_of_freedom = 3.0}), wsl::random::uniform({.lower = 0.0, .upper = 1.0}));
        wsl::random::product_measure const J(wsl::random::normal({.mean = 0.0, .stddev = 1.0}), wsl::random::gamma({.rate = 1.0, .shape = 2.0}));
        auto const LP = wsl::random::pushforward(wsl::random::normal({.mean = 0.0, .stddev = 1.0}), wsl::random::increasing([](double t) { return wsl::special::exp(t); }, [](double t) { return wsl::special::log(t); }));
        ok &= check("sum (discrète)", fingerprint([&] { return static_cast<std::uint64_t>(wsl::random::sample(PB, x)); }), 0xb95ba5ee227f675e);
        ok &= check("sum (continue)", fingerprint([&] { return bits(wsl::random::sample(ST, x)); }), 0xa5686f02baa4b6f2);
        ok &= check("independent", fingerprint([&] { auto const [u, w] = wsl::random::sample(J, x); return bits(u) ^ (bits(w) << 1); }), 0xf8b02c1c396edb9f);
        ok &= check("pushforward", fingerprint([&] { return bits(wsl::random::sample(LP, x)); }), 0x4bf00f5f9cc7fcfc);
        ok &= check("sum functions", fingerprint([&] {
            double const u = wsl::random::uniform_real(x);
            auto const k = static_cast<std::int64_t>(wsl::random::uniform_below(40, x));
            return bits(wsl::random::pmf(PB, k)) ^ (bits(wsl::random::cdf(PB, k)) << 1) ^ (bits(wsl::random::ccdf(PB, k)) << 2) ^ (static_cast<std::uint64_t>(wsl::random::quantile(PB, u)) << 3)
                 ^ (bits(wsl::random::cdf(LP, u + 0.5)) << 4) ^ (bits(wsl::random::quantile(LP, u)) << 5);
        }, 1 << 9), 0xfb59f13d32d4c878);
        ok &= check("algebra functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wsl::random::uniform_real(x), t = fp::fmul(20.0, u) - 10;
            return bits(wsl::random::cdf(Z, t)) ^ (bits(wsl::random::quantile(Z, u)) << 1) ^ (bits(wsl::random::cdf(T, t + 13)) << 2) ^ (bits(wsl::random::quantile(T, u)) << 3)
                 ^ (bits(wsl::random::cdf(M, t)) << 4) ^ (bits(wsl::random::quantile(M, u)) << 5) ^ (bits(wsl::random::cquantile(wsl::random::gamma({.rate = 1.0, .shape = 0.5}), u)) << 6);
        }, 1 << 9), 0x22bda422ad262729);
    }
    assert(ok);
    std::printf("test_random_algebra : ok\n");
}
