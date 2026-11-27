// Algèbre des lois de wsl::random : quantile de survie, image affine et identités,
// troncature, mélanges ; vérifications exactes, statistiques (graines fixes) et
// empreintes des tirages.

#include <wsl/probability.hpp>
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
            v.push_back(double(wr::sample(law, g)));
        std::ranges::sort(v);
        double worst = 0;
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            double const f = wr::cdf(law, v[i]);
            worst = std::max({worst, std::abs(f - double(i) / n), std::abs(f - double(i + 1) / n)});
        }
        return worst * std::sqrt(double(n));
    }

    bool close(double a, double b, double relative) { return std::abs(a - b) <= relative * std::abs(b); }

    // wr::ccdf(X, wr::cquantile(X, q)) ≈ q, et wr::cquantile(X, q) = wr::quantile(X, 1 − q) pour q >= ½.
    template <class L>
    void check_cquantile(L const& X, double tolerance = 1e-12)
    {
        for (double q : {0.5, 0.75, 0.9375, 1 - 0x1p-40})
            assert(wr::cquantile(X, q) == wr::quantile(X, 1 - q));
        for (double q : {0.25, 1e-5, 1e-30, 1e-200})
        {
            double const x = wr::cquantile(X, q);
            if (std::isfinite(x) && x != X.support().lower && x != X.support().upper && !close(double(wr::ccdf(X, x)), q, tolerance))
            {
                std::printf("wr::cquantile(%g) = %a, ccdf = %.17g\n", q, x, double(wr::ccdf(X, x)));
                assert(false);
            }
        }
        assert(wr::cquantile(X, 1.0) == wr::quantile(X, 0.0) && wr::cquantile(X, 0.0) == wr::quantile(X, 1.0));
        assert(throws_domain_error([&] { (void)wr::cquantile(X, 1.5); }));
    }
}

// Identités : les familles fermées gardent leur type.
using normal_t = wr::normal<double>;
static_assert(std::is_same_v<decltype(wr::affine(std::declval<normal_t>(), {.scale = 2.0, .shift = 1})), normal_t>);
static_assert(std::is_same_v<decltype(wr::affine(std::declval<wr::cauchy<>>(), {.scale = -1})), wr::cauchy<>>);
static_assert(std::is_same_v<decltype(wr::affine(std::declval<wr::uniform<float>>(), {.scale = -1, .shift = 3})), wr::uniform<float>>);
static_assert(std::is_same_v<decltype(wr::affine(std::declval<wr::laplace<>>(), {.scale = 0.25})), wr::laplace<>>);
static_assert(std::is_same_v<decltype(wr::affine(std::declval<wr::exponential<>>(), {.scale = 2})), wr::affine_image<wr::exponential<>>>);
static_assert(std::is_same_v<decltype(wr::affine(wr::affine(std::declval<wr::student_t<>>(), {.scale = 2}), {.scale = 2, .shift = 1})), wr::affine_image<wr::student_t<>>>);
static_assert(wr::affine_closed<wr::logistic<>> && wr::scale_closed<wr::gamma<>> && !wr::scale_closed<wr::beta<>>);
static_assert(!wr::has_mean<decltype(wr::affine(std::declval<wr::cauchy<>>(), {.scale = 2.0}))>);
static_assert(wr::has_mean<wr::affine_image<wr::student_t<>>> && wr::has_cquantile<wr::affine_image<wr::beta<>>>);
static_assert(!wr::has_mean<wr::truncated<normal_t>> && wr::has_cquantile<wr::truncated<normal_t>>);
static_assert(wr::probability_law<wr::mixture<normal_t, wr::gamma<>>> && wr::has_pmf<wr::mixture<wr::poisson<>, wr::binomial<>>>);
static_assert(wr::mixture<normal_t, wr::exponential<>>::sampling == wr::exactness::exact);
static_assert(wr::mixture<normal_t, wr::gamma<>>::sampling == wr::exactness::approximate);
static_assert(!wr::has_pdf<wr::mixture<wr::poisson<>>> && wr::has_cdf<wr::mixture<wr::poisson<>>>);
// Évaluables à la compilation.
static_assert(wr::affine(wr::normal({.mean = 1.0, .stddev = 2.0}), {.scale = 3.0, .shift = -1.0}).mean() == 2.0);
static_assert(wr::cdf(wr::affine(wr::exponential({.rate = 1.0}), {.scale = 2.0}), 0.0) == 0.0);

int main()
{
    wr::xoshiro256pp g(20261118);

    // ---- Quantile de survie --------------------------------------------------------
    check_cquantile(wr::exponential({.rate = 2.0}));
    check_cquantile(wr::normal({.mean = 1.0, .stddev = 2.0}));
    assert(wr::cquantile(wr::uniform({.lower = -1.0, .upper = 3.0}), 0.75) == wr::quantile(wr::uniform({.lower = -1.0, .upper = 3.0}), 0.25));   // (queue : x arrondi près de b)
    check_cquantile(wr::gamma({.rate = 2.0, .shape = 2.5}));
    check_cquantile(wr::gamma({.scale = 1.0, .shape = 1e-3}), 1e-10);
    check_cquantile(wr::beta({.alpha = 2.0, .beta = 0.5}), 1e-6);           // x près de 1 : 1 − x mal résolu
    check_cquantile(wr::chi_squared({.degrees_of_freedom = 3.0}));
    check_cquantile(wr::student_t({.degrees_of_freedom = 4.0}));
    check_cquantile(wr::lognormal({.log_mean = 0.0, .log_stddev = 0.5}));
    check_cquantile(wr::cauchy({.location = 1.0, .scale = 2.0}));
    check_cquantile(wr::fisher_f({.d1 = 5.0, .d2 = 10.0}), 1e-10);
    check_cquantile(wr::weibull({.scale = 2.0, .shape = 1.5}));
    check_cquantile(wr::laplace({.location = 1.0, .scale = 2.0}));
    check_cquantile(wr::logistic({.location = 1.0, .scale = 2.0}));
    check_cquantile(wr::pareto({.scale = 1.0, .shape = 2.5}));
    assert(wr::cquantile(wr::uniform({.lower = 1.0, .upper = 3.0}), 0.1) == 2.8);
    assert(wr::cquantile(wr::exponential({.rate = 1.0}), 1e-300) > 690);     // 1 − q ne vaudrait rien

    // ---- Image affine : identités ---------------------------------------------------
    {
        auto const Y = wr::affine(wr::normal({.mean = 1.0, .stddev = 2.0}), {.scale = 2.0, .shift = 3});
        assert(Y.mean() == 5 && Y.stddev() == 4);
        auto const U = wr::affine(wr::uniform({.lower = 0.0, .upper = 1.0}), {.scale = -2, .shift = 2});
        assert(U.lower() == 0 && U.upper() == 2);
        auto const C = wr::affine(wr::cauchy({.location = 1.0, .scale = 2.0}), {.scale = -1});
        assert(C.location() == -1 && C.scale() == 2);
        // Familles fermées par changement d'échelle : mêmes bits que la loi de la famille.
        auto const E4 = wr::affine(wr::exponential({.rate = 2.0}), {.scale = 4});
        wr::exponential const E({.rate = 0.5});
        auto const G3 = wr::affine(wr::gamma({.scale = 0.5, .shape = 2.5}), {.scale = 3});
        wr::gamma const G({.scale = 1.5, .shape = 2.5});
        auto const W2 = wr::affine(wr::weibull({.scale = 2.0, .shape = 1.5}), {.scale = 2});
        wr::weibull const W({.scale = 4.0, .shape = 1.5});
        auto const P5 = wr::affine(wr::pareto({.scale = 2.0, .shape = 2.5}), {.scale = 0.5});
        wr::pareto const P({.scale = 1.0, .shape = 2.5});
        for (double x : {0.01, 0.5, 1.0, 3.0, 40.0})
        {
            assert(wr::cdf(E4, x) == wr::cdf(E, x) && wr::pdf(E4, x) == wr::pdf(E, x) && wr::ccdf(G3, x) == wr::ccdf(G, x) && wr::pdf(G3, x) == wr::pdf(G, x));
            assert(wr::cdf(W2, x) == wr::cdf(W, x) && wr::log_pdf(W2, x) == wr::log_pdf(W, x) && wr::ccdf(P5, x) == wr::ccdf(P, x));
        }
        for (double p : {1e-100, 0.25, 0.5, 0.99})
            assert(wr::quantile(E4, p) == wr::quantile(E, p) && wr::cquantile(G3, p) == wr::cquantile(G, p) && wr::quantile(P5, p) == wr::quantile(P, p));
        assert(E4.mean() == 2 && E4.variance() == E.variance() && G3.entropy() == G.entropy() && W2.median() == W.median());
        auto const L2 = wr::affine(wr::lognormal({.log_mean = 0.0, .log_stddev = 0.5}), {.scale = 2});
        assert(std::abs(L2.median() - 2) < 4e-16 && close(wr::cdf(L2, 3.0), wr::cdf(wr::lognormal({.log_mean = std::log(2.0), .log_stddev = 0.5}), 3.0), 1e-15));
        // aX, a < 0 ou b ≠ 0 : composition.
        auto const S = wr::affine(wr::exponential({.rate = 1.0}), {.scale = -2, .shift = 1});              // support (−∞, 1]
        assert(S.support().lower == -INFINITY && S.support().upper == 1 && wr::cdf(S, 1.0) == 1 && wr::cdf(S, -1.0) == wr::ccdf(wr::exponential({.rate = 1.0}), 1.0));
        assert(wr::quantile(S, 0.25) == 1 - 2 * wr::cquantile(wr::exponential({.rate = 1.0}), 0.25));
        assert(throws_domain_error([] { (void)wr::affine(wr::normal({.mean = 0.0, .stddev = 1.0}), {.scale = 0}); }));
        assert(throws_domain_error([] { (void)wr::affine(wr::student_t({.degrees_of_freedom = 3.0}), {.scale = 1, .shift = INFINITY}); }));
    }

    // ---- Image affine générique ----------------------------------------------------
    {
        wr::student_t const T({.degrees_of_freedom = 4.0});
        auto const Z = wr::affine(T, {.scale = -2, .shift = 5});
        for (double y : {-11.0, 1.0, 3.0, 5.0, 9.0, 1e6})                    // (y − 5)/(−2) exact
        {
            double const z = (y - 5) / -2;
            assert(wr::cdf(Z, y) == wr::ccdf(T, z) && wr::ccdf(Z, y) == wr::cdf(T, z) && wr::pdf(Z, y) == wr::pdf(T, z) / 2);
        }
        assert(Z.mean() == 5 && Z.variance() == 4 * T.variance() && Z.median() == 5 && Z.entropy() == T.entropy() + std::log(2.0));
        for (double p : {1e-12, 0.1, 0.5, 0.9})
            assert(wr::quantile(Z, p) == 5 - 2 * wr::cquantile(T, p) && wr::cquantile(Z, p) == 5 - 2 * wr::quantile(T, p));
        assert(kolmogorov_smirnov(Z, g) < 1.95);
        auto const B = wr::affine(wr::beta({.alpha = 2.0, .beta = 3.0}), {.scale = 10, .shift = -4});
        assert(B.support().lower == -4 && B.support().upper == 6 && std::abs(B.mean() - 0) < 1e-15 && B.skewness() == wr::beta({.alpha = 2.0, .beta = 3.0}).skewness());
        assert(kolmogorov_smirnov(B, g) < 1.95);
        auto const N = wr::affine(wr::gamma({.rate = 1.0, .shape = 3.0}), {.scale = -1});
        assert(N.skewness() == -wr::gamma({.rate = 1.0, .shape = 3.0}).skewness() && kolmogorov_smirnov(N, g) < 1.95);
        auto const F = wr::affine(wr::student_t<float>({.degrees_of_freedom = 3.0f}), {.scale = 2.0, .shift = 1.0});   // float
        static_assert(std::is_same_v<decltype(wr::cdf(F, 1.0f)), float>);
        assert(wr::cdf(F, 1.0f) == 0.5f);
    }

    // ---- Troncature ----------------------------------------------------------------
    {
        wr::normal const N({.mean = 0.0, .stddev = 1.0});
        auto const T = wr::truncate(N, {.lower = 10.0});                     // queue droite : par la survie
        assert(T.mass() == wr::ccdf(N, 10.0) && T.support().lower == 10 && wr::cdf(T, 9.0) == 0 && wr::pdf(T, 9.0) == 0);
        assert(close(wr::pdf(T, 10.0), wr::pdf(N, 10.0) / wr::ccdf(N, 10.0), 1e-15) && close(wr::ccdf(T, 10.5), wr::ccdf(N, 10.5) / wr::ccdf(N, 10.0), 1e-13));
        // Aller-retour à 10⁻¹⁵ près en absolu : près d'une borne, (S(a) − S(x))/m se compense.
        for (double p : {1e-8, 0.3, 0.5, 0.999})
            assert(std::abs(wr::cdf(T, wr::quantile(T, p)) - p) < 1e-14 && std::abs(wr::ccdf(T, wr::cquantile(T, p)) - p) < 1e-14);
        assert(wr::quantile(T, 0.0) == 10 && wr::quantile(T, 1.0) == INFINITY);
        double s = 0;                                                        // E[N | N >= a] = φ(a)/Φ̄(a)
        for (int i = 0; i < 20000; ++i)
            s += wr::sample(T, g);
        assert(std::abs(s / 20000 - wr::pdf(N, 10.0) / wr::ccdf(N, 10.0)) < 0.005);
        assert(kolmogorov_smirnov(T, g) < 1.95);
        auto const G = wr::truncate(wr::gamma({.rate = 1.0, .shape = 2.0}), {.lower = 1.0, .upper = 3.0});   // masse 0,54 : rejet
        assert(close(G.mass(), wr::cdf(wr::gamma({.rate = 1.0, .shape = 2.0}), 3.0) - wr::cdf(wr::gamma({.rate = 1.0, .shape = 2.0}), 1.0), 1e-15));
        assert(kolmogorov_smirnov(G, g) < 1.95 && close(wr::cdf(G, G.median()), 0.5, 1e-14));
        auto const L = wr::truncate(wr::cauchy({.location = 0.0, .scale = 1.0}), {.upper = -1e3});   // queue gauche
        assert(kolmogorov_smirnov(L, g) < 1.95 && L.support().upper == -1e3);
        auto const Z = wr::truncate(wr::affine(wr::student_t({.degrees_of_freedom = 2.0}), {.scale = 2.0, .shift = 1}), {.lower = -1.0, .upper = 0.5});   // tronquer une image
        assert(kolmogorov_smirnov(Z, g) < 1.95);
        assert(throws_domain_error([&] { (void)wr::truncate(N, {.lower = 40.0, .upper = 50.0}); }));   // masse sous les sous-normaux
        assert(throws_domain_error([&] { (void)wr::truncate(N, {.lower = 1.0, .upper = 1.0}); }));
        assert(throws_domain_error([] { (void)wr::truncate(wr::exponential({.rate = 1.0}), {.lower = -5.0, .upper = -1.0}); }));
    }

    // ---- Mélanges ------------------------------------------------------------------
    {
        wr::normal const A({.mean = -2.0, .stddev = 1.0});
        wr::gamma const B({.rate = 1.0, .shape = 3.0});
        wr::mixture const M(wr::weighted(0.3, A), wr::weighted(0.7, B));
        wr::mixture const M10(wr::weighted(3.0, A), wr::weighted(7.0, B));    // poids non normalisés
        for (double x : {-3.0, 0.0, 0.5, 2.0, 10.0})
        {
            double const f = 0.3 * wr::pdf(A, x) + 0.7 * wr::pdf(B, x);
            assert(close(wr::pdf(M, x), f, 1e-15) && close(wr::pdf(M10, x), f, 1e-15) && close(std::exp(wr::log_pdf(M, x)), f, 1e-14));
            assert(close(wr::cdf(M, x) + wr::ccdf(M, x), 1, 1e-15));
        }
        assert(close(M.mean(), 1.5, 1e-15) && close(M.variance(), 7.65, 1e-15) && std::abs(M.skewness()) < 1e-15);   // μ₃ = 0 ici
        for (double p : {1e-9, 0.1, 0.5, 0.9, 1 - 1e-9})
            assert(close(wr::cdf(M, wr::quantile(M, p)), p, 1e-13) && close(wr::ccdf(M, wr::cquantile(M, 1 - p)), 1 - p, 1e-13));
        assert(close(wr::ccdf(M, wr::cquantile(M, 1e-30)), 1e-30, 1e-12));
        assert(M.support().lower == -INFINITY && M.support().upper == INFINITY);
        assert(kolmogorov_smirnov(M, g) < 1.95);
        // Une composante : la loi elle-même.
        wr::mixture const one(wr::weighted(0.7, B));
        for (double x : {0.5, 2.0, 9.0})
            assert(wr::pdf(one, x) == wr::pdf(B, x) && wr::cdf(one, x) == wr::cdf(B, x));
        // Choix exact des composantes : poids 1, 1, 2 ; χ² à 2 ddl.
        wr::mixture const C(wr::weighted(1.0, wr::uniform({.lower = 0.0, .upper = 1.0})), wr::weighted(1.0, wr::uniform({.lower = 1.0, .upper = 2.0})),
                             wr::weighted(2.0, wr::uniform({.lower = 2.0, .upper = 3.0})));
        int counts[3] = {};
        for (int i = 0; i < 40000; ++i)
            ++counts[static_cast<int>(wr::sample(C, g))];
        double chi = 0;
        for (int i = 0; i < 3; ++i)
        {
            double const e = 40000.0 * (i == 2 ? 0.5 : 0.25);
            chi += (counts[i] - e) * (counts[i] - e) / e;
        }
        assert(chi < 13.8);                                                  // seuil de 0,1 %
        // Poids d'étendue immense : grands entiers ; la composante de poids 2^-1000 n'apparaît pas.
        wr::mixture const H(wr::weighted(0x1p1000, wr::uniform({.lower = 0.0, .upper = 1.0})), wr::weighted(0x1p-1000, wr::uniform({.lower = 5.0, .upper = 6.0})),
                             wr::weighted(0x1p999, wr::uniform({.lower = 1.0, .upper = 2.0})));
        int above = 0;
        for (int i = 0; i < 30000; ++i)
        {
            double const x = wr::sample(H, g);
            assert(x < 2);
            above += x >= 1;
        }
        assert(std::abs(above / 30000.0 - 1.0 / 3) < 0.012);
        // Discret : masse, répartition, quantile exact de la répartition du mélange.
        wr::mixture const D(wr::weighted(1.0, wr::poisson({.mean = 2.0})), wr::weighted(3.0, wr::poisson({.mean = 20.0})));
        assert(close(wr::pmf(D, 2), 0.25 * wr::pmf(wr::poisson({.mean = 2.0}), 2) + 0.75 * wr::pmf(wr::poisson({.mean = 20.0}), 2), 1e-15));
        for (double p : {0.01, 0.2, 0.5, 0.99})
        {
            auto const k = wr::quantile(D, p);
            assert(wr::cdf(D, k) >= p && (k == 0 || wr::cdf(D, k - 1) < p));
        }
        assert(D.mean() == 15.5 && D.variance() == 76.25);
        assert(throws_domain_error([&] { (void)wr::mixture(wr::weighted(0.0, A), wr::weighted(1.0, B)); }));
        assert(throws_domain_error([&] { (void)wr::mixture(wr::weighted(-1.0, A)); }));
    }

    // ---- Sommes indépendantes ------------------------------------------------------
    {
        static_assert(std::is_same_v<decltype(wr::convolution(wr::normal({.mean = 1.0, .stddev = 3.0}), wr::normal({.mean = 2.0, .stddev = 4.0}))), wr::normal<>>);
        static_assert(std::is_same_v<decltype(wr::convolution(wr::poisson({.mean = 1.0}), wr::poisson({.mean = 2.0}))), wr::poisson<>>);
        static_assert(std::is_same_v<decltype(wr::convolution(wr::gamma({.rate = 1.0, .shape = 1.0}), wr::gamma({.rate = 1.0, .shape = 1.0}))), wr::convolution_law<wr::gamma<>, wr::gamma<>>>);
        static_assert(wr::has_cdf<wr::convolution_law<wr::gamma<>, wr::exponential<>>> && !wr::has_cdf<wr::convolution_law<wr::student_t<>, wr::uniform<>>>);
        static_assert(wr::convolution_law<wr::poisson<>, wr::geometric<>>::sampling == wr::exactness::approximate);
        static_assert(wr::convolution_law<wr::discrete_uniform<>, wr::discrete_uniform<>>::sampling == wr::exactness::exact);
        auto const S = wr::convolution(wr::normal({.mean = 1.0, .stddev = 3.0}), wr::normal({.mean = 2.0, .stddev = 4.0}));
        assert(S.mean() == 3 && S.stddev() == 5);
        auto const D = wr::convolution(wr::normal({.mean = 1.0, .stddev = 3.0}), wr::affine(wr::normal({.mean = 2.0, .stddev = 4.0}), {.scale = -1}));
        assert(D.mean() == -1 && D.stddev() == 5);
        auto const C = wr::convolution(wr::cauchy({.location = 1.0, .scale = 2.0}), wr::cauchy({.location = -3.0, .scale = 0.5}));
        assert(C.location() == -2 && C.scale() == 2.5);
        assert(wr::convolution(wr::chi_squared({.degrees_of_freedom = 3.0}), wr::chi_squared({.degrees_of_freedom = 4.5})).degrees_of_freedom() == 7.5);
        // Identités conditionnelles : mêmes bits que la loi de la famille.
        auto const G = wr::convolution(wr::gamma({.rate = 2.0, .shape = 2.0}), wr::gamma({.scale = 0.5, .shape = 1.5}));   // paramétrisations mêlées
        wr::gamma const G35({.rate = 2.0, .shape = 3.5});
        auto const E = wr::convolution(wr::exponential({.rate = 2.0}), wr::gamma({.rate = 2.0, .shape = 2.5}));
        auto const B = wr::convolution(wr::binomial({.trials = 10, .p = 0.3}), wr::binomial({.trials = 5, .p = 0.3}));
        wr::binomial const B15({.trials = 15, .p = 0.3});
        assert(G.family() && E.family() && B.family());
        for (double x : {0.1, 1.0, 4.0})
            assert(wr::cdf(G, x) == wr::cdf(G35, x) && wr::pdf(G, x) == wr::pdf(G35, x) && wr::ccdf(E, x) == wr::ccdf(G35, x));
        assert(wr::quantile(G, 0.3) == wr::quantile(G35, 0.3) && G.variance() == G35.variance());
        for (int k = 0; k <= 15; ++k)
            assert(wr::pmf(B, k) == wr::pmf(B15, k) && wr::cdf(B, k) == wr::cdf(B15, k));
        // Sans identité : moments par les cumulants ; NaN pour une somme continue.
        auto const H = wr::convolution(wr::gamma({.rate = 2.0, .shape = 2.0}), wr::gamma({.rate = 3.0, .shape = 1.5}));
        assert(!H.family() && std::isnan(wr::cdf(H, 1.0)) && throws_domain_error([&] { (void)wr::quantile(H, 0.5); }));
        assert(H.mean() == 1.5 && close(H.variance(), 0.5 + 1.5 / 9, 1e-15));
        double const k3 = 2.0 * 2 / 8 + 2.0 * 1.5 / 27, v = 0.5 + 1.5 / 9;                     // κ₃ = 2a/λ³
        assert(close(H.skewness(), k3 / (v * std::sqrt(v)), 1e-14));
        double m = 0;
        for (int i = 0; i < 40000; ++i)
            m += wr::sample(H, g);
        assert(std::abs(m / 40000 - 1.5) < 0.02);
        // Convolution discrète : masse, répartition, survie, quantile exact.
        wr::poisson const P({.mean = 3.5});
        wr::binomial const Bi({.trials = 20, .p = 0.3});
        auto const PB = wr::convolution(P, Bi);
        double total = 0;
        for (long long k = 0; k <= 60; ++k)
        {
            double direct = 0;
            for (long long j = 0; j <= std::min<long long>(k, 20); ++j)
                direct += wr::pmf(Bi, j) * wr::pmf(P, k - j);
            assert(close(wr::pmf(PB, k), direct, 1e-13));
            assert(close(wr::cdf(PB, k) + wr::ccdf(PB, k), 1, 1e-15));
            total += wr::pmf(PB, k);
        }
        assert(close(total, 1, 1e-14) && PB.mean() == 9.5 && close(PB.variance(), 3.5 + 4.2, 1e-15));
        for (double p : {1e-12, 0.05, 0.5, 0.99, 1 - 1e-10})
        {
            auto const k = wr::quantile(PB, p);
            assert(wr::cdf(PB, k) >= p && (k == 0 || wr::cdf(PB, k - 1) < p));
        }
        assert(close(wr::ccdf(PB, 80), wr::ccdf(P, 60) * 0 + wr::ccdf(PB, 80), 0) && wr::ccdf(PB, 80) > 0 && wr::ccdf(PB, 80) < 1e-30);   // queue lointaine, sans 1 − cdf
        auto const GU = wr::convolution(wr::geometric({.p = 0.2}), wr::discrete_uniform({.lower = -3, .upper = 3}));        // supports négatifs
        assert(GU.support().lower == -3 && close(wr::cdf(GU, 0), [&] { double s = 0; for (int j = -3; j <= 3; ++j) s += wr::cdf(wr::geometric({.p = 0.2}), -j) / 7; return s; }(), 1e-14));
    }

    // ---- Loi jointe et image par une fonction ----------------------------------------
    {
        namespace sp = wsl::special;
        wr::normal const N({.mean = 0.0, .stddev = 1.0});
        wr::exponential const E({.rate = 2.0});
        wr::product_measure const J(N, E, wr::poisson({.mean = 3.0}));
        static_assert(std::is_same_v<decltype(J)::value_type, std::tuple<double, double, std::int64_t>>);
        static_assert(decltype(J)::sampling == wr::exactness::approximate && !wr::has_pdf<decltype(J)>);
        wr::product_measure const K(N, E);
        static_assert(decltype(K)::sampling == wr::exactness::exact);
        assert(close(wr::pdf(K, {0.5, 1.0}), wr::pdf(N, 0.5) * wr::pdf(E, 1.0), 1e-15) && close(wr::cdf(K, {0.0, 1.0}), 0.5 * wr::cdf(E, 1.0), 1e-15));
        assert(close(wr::log_pdf(K, {0.5, 1.0}), wr::log_pdf(N, 0.5) + wr::log_pdf(E, 1.0), 1e-15) && close(K.entropy(), N.entropy() + E.entropy(), 1e-15));
        assert(wr::pdf(K, {0.5, -1.0}) == 0 && wr::log_pdf(K, {0.5, -1.0}) == -INFINITY && std::get<1>(K.mean()) == 0.5);
        wr::xoshiro256pp a(9), b(9);                                       // tirage : composantes dans l'ordre
        auto const [x, y, k] = wr::sample(J, a);
        assert(x == wr::sample(N, b) && y == wr::sample(E, b) && k == wr::sample(wr::poisson({.mean = 3.0}), b));
        // Image monotone : exp(N) est la log-normale ; 1/X, X gamma, la gamma inverse.
        auto const L = wr::pushforward(N, wr::increasing([](double t) { return sp::exp(t); }, [](double t) { return sp::log(t); }, [](double t) { return 1 / t; }));
        wr::lognormal const LN({.log_mean = 0.0, .log_stddev = 1.0});
        for (double t : {0.1, 1.0, 2.0, 30.0})
            assert(close(wr::cdf(L, t), wr::cdf(LN, t), 1e-15) && close(wr::pdf(L, t), wr::pdf(LN, t), 1e-15) && close(wr::log_pdf(L, t), wr::log_pdf(LN, t), 1e-14));
        assert(close(wr::quantile(L, 0.9), wr::quantile(LN, 0.9), 1e-15) && L.median() == 1 && L.support().lower == 0 && kolmogorov_smirnov(L, g) < 1.95);
        auto const I = wr::pushforward(wr::gamma({.rate = 2.0, .shape = 3.0}), wr::decreasing([](double t) { return 1 / t; }, [](double t) { return 1 / t; }, [](double t) { return -1 / (t * t); }));
        assert(close(wr::cdf(I, 1.0), wr::ccdf(wr::gamma({.rate = 2.0, .shape = 3.0}), 1.0), 1e-15) && close(wr::quantile(I, 0.25), 1 / wr::cquantile(wr::gamma({.rate = 2.0, .shape = 3.0}), 0.25), 1e-15));
        assert(kolmogorov_smirnov(I, g) < 1.95 && wr::pdf(I, -1.0) == 0);
        // Fonction quelconque : tirage seul. N² est la loi du χ² à 1 degré.
        auto const Q = wr::pushforward(N, [](double t) { return t * t; });
        static_assert(!wr::has_cdf<decltype(Q)> && wr::probability_law<decltype(Q)>);
        std::vector<double> v;
        for (int i = 0; i < 20000; ++i)
            v.push_back(wr::sample(Q, g));
        std::ranges::sort(v);
        double worst = 0;
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            double const f = wr::cdf(wr::chi_squared({.degrees_of_freedom = 1.0}), v[i]);
            worst = std::max({worst, std::abs(f - double(i) / 20000), std::abs(f - double(i + 1) / 20000)});
        }
        assert(worst * std::sqrt(20000.0) < 1.95);
    }

    // ---- Empreintes ----------------------------------------------------------------
    bool ok = true;
    {
        wr::xoshiro256pp x(42);
        auto const Z = wr::affine(wr::student_t({.degrees_of_freedom = 5.0}), {.scale = -1.5, .shift = 2});
        auto const T = wr::truncate(wr::normal({.mean = 0.0, .stddev = 1.0}), {.lower = 3.0});     // inversion
        auto const R = wr::truncate(wr::laplace({.location = 0.0, .scale = 1.0}), {.lower = -1.0, .upper = 2.0}); // rejet, tirage de X exact
        wr::mixture const M(wr::weighted(0.3, wr::normal({.mean = -2.0, .stddev = 1.0})), wr::weighted(0.7, wr::gamma({.rate = 1.0, .shape = 3.0})));
        wr::mixture const H(wr::weighted(0x1p600, wr::exponential({.rate = 1.0})), wr::weighted(0x1p-600, wr::normal({.mean = 0.0, .stddev = 1.0})),
                             wr::weighted(0x1p599, wr::laplace({.location = 0.0, .scale = 1.0})));
        wr::mixture const D(wr::weighted(1.0, wr::poisson({.mean = 2.0})), wr::weighted(3.0, wr::binomial({.trials = 40, .p = 0.5})));
        ok &= check("affine", fingerprint([&] { return bits(wr::sample(Z, x)); }), 0x44509b09b529826c);
        ok &= check("truncated (inv.)", fingerprint([&] { return bits(wr::sample(T, x)); }), 0x398b1148ffa1462f);
        ok &= check("truncated (rejet)", fingerprint([&] { return bits(wr::sample(R, x)); }), 0x550c9bd3098ec373);
        ok &= check("mixture", fingerprint([&] { return bits(wr::sample(M, x)); }), 0x1cd22859dc5c692f);
        ok &= check("mixture (grands)", fingerprint([&] { return bits(wr::sample(H, x)); }), 0xa71e5cac21b536d6);
        ok &= check("mixture discrète", fingerprint([&] { return static_cast<std::uint64_t>(wr::sample(D, x)); }), 0x2c54e55c5ac36aab);
        auto const PB = wr::convolution(wr::poisson({.mean = 3.5}), wr::binomial({.trials = 20, .p = 0.3}));
        auto const ST = wr::convolution(wr::student_t({.degrees_of_freedom = 3.0}), wr::uniform({.lower = 0.0, .upper = 1.0}));
        wr::product_measure const J(wr::normal({.mean = 0.0, .stddev = 1.0}), wr::gamma({.rate = 1.0, .shape = 2.0}));
        auto const LP = wr::pushforward(wr::normal({.mean = 0.0, .stddev = 1.0}), wr::increasing([](double t) { return wsl::special::exp(t); }, [](double t) { return wsl::special::log(t); }));
        ok &= check("sum (discrète)", fingerprint([&] { return static_cast<std::uint64_t>(wr::sample(PB, x)); }), 0xb95ba5ee227f675e);
        ok &= check("sum (continue)", fingerprint([&] { return bits(wr::sample(ST, x)); }), 0xa5686f02baa4b6f2);
        ok &= check("independent", fingerprint([&] { auto const [u, w] = wr::sample(J, x); return bits(u) ^ (bits(w) << 1); }), 0xf8b02c1c396edb9f);
        ok &= check("pushforward", fingerprint([&] { return bits(wr::sample(LP, x)); }), 0x4bf00f5f9cc7fcfc);
        ok &= check("sum functions", fingerprint([&] {
            double const u = wr::uniform_real(x);
            auto const k = static_cast<std::int64_t>(wr::uniform_below(40, x));
            return bits(wr::pmf(PB, k)) ^ (bits(wr::cdf(PB, k)) << 1) ^ (bits(wr::ccdf(PB, k)) << 2) ^ (static_cast<std::uint64_t>(wr::quantile(PB, u)) << 3)
                 ^ (bits(wr::cdf(LP, u + 0.5)) << 4) ^ (bits(wr::quantile(LP, u)) << 5);
        }, 1 << 9), 0xfb59f13d32d4c878);
        ok &= check("algebra functions", fingerprint([&] {
            namespace fp = wsl::special::detail;
            double const u = wr::uniform_real(x), t = fp::fmul(20.0, u) - 10;
            return bits(wr::cdf(Z, t)) ^ (bits(wr::quantile(Z, u)) << 1) ^ (bits(wr::cdf(T, t + 13)) << 2) ^ (bits(wr::quantile(T, u)) << 3)
                 ^ (bits(wr::cdf(M, t)) << 4) ^ (bits(wr::quantile(M, u)) << 5) ^ (bits(wr::cquantile(wr::gamma({.rate = 1.0, .shape = 0.5}), u)) << 6);
        }, 1 << 9), 0x22bda422ad262729);
    }
    // ---- Intensités rationnelles exactes (décision 5) ----------------------------------
    {
                // forme canonique : == est l'égalité en loi
        assert(wr::gamma({.rate = 2.0, .shape = 3.0}) == wr::gamma({.scale = 0.5, .shape = 3.0}));
        assert(wr::gamma({.rate = wr::ratio(1, 3), .shape = 3.0}) == wr::gamma({.scale = 3.0, .shape = 3.0}));
        assert(!(wr::gamma({.rate = 1.0 / 3, .shape = 3.0}) == wr::gamma({.scale = 3.0, .shape = 3.0})));   // 1/3 arrondi : autre loi
        assert(wr::exponential({.scale = 3.0}) == wr::exponential({.rate = wr::ratio(2, 6)}) && wr::exponential() == wr::exponential({.rate = 1.0}));
        assert(throws_domain_error([] { (void)wr::gamma({.rate = 1.0, .scale = 1.0, .shape = 2.0}); }));
        assert(throws_domain_error([] { (void)wr::exponential({.rate = wr::ratio(1, 0)}); }));
        // param() : forme lisible, aller-retour exact
        auto const G3 = wr::gamma({.scale = 3.0, .shape = 2.5});
        assert(!G3.param().rate && G3.param().scale->numerator() == 3 && wr::gamma(G3.param()) == G3);
        auto const H = wr::exponential({.rate = wr::ratio(2, 3)});                                 // échelle 1,5 : un double
        assert(!H.param().rate && H.param().scale->numerator() == 1.5 && wr::exponential(H.param()) == H);
        auto const R = wr::exponential({.rate = wr::ratio(3, 5)});                                 // ni λ ni 1/λ ne sont des doubles
        assert(R.param().rate && R.param().rate->numerator() == 3 && R.param().rate->denominator() == 5 && wr::exponential(R.param()) == R);
        // valeurs : celles de l'échelle exacte
        wr::exponential const E3({.scale = 3.0});
        assert(E3.mean() == 3 && E3.variance() == 9 && E3.scale() == 3 && E3.cdf(3.0) == wr::exponential({.rate = 1.0}).cdf(1.0));
        assert(kolmogorov_smirnov(E3, g) < 1.95 && kolmogorov_smirnov(wr::gamma({.rate = wr::ratio(3, 7), .shape = 0.6}), g) < 1.95);
        auto const K = wr::erlang({.scale = 3.0, .shape = 4});
        assert(K.mean() == 12 && K.scale() == 3 && kolmogorov_smirnov(K, g) < 1.95);
        // algèbre : c X et sommes, intensités exactes
        auto const S = wr::affine(wr::gamma({.rate = 2.0, .shape = 3.0}), {.scale = 3});            // gamma(3, 2/3) exactement
        assert(S.cdf(1.0) == wr::gamma({.rate = wr::ratio(2, 3), .shape = 3.0}).cdf(1.0));
        auto const C = wr::convolution(wr::exponential({.scale = 3.0}), wr::gamma({.rate = wr::ratio(1, 3), .shape = 2.0}));
        assert(C.family() && *C.family() == wr::gamma({.scale = 3.0, .shape = 3.0}));
        auto const C2 = wr::convolution(wr::exponential({.scale = 3.0}), wr::exponential({.rate = wr::ratio(1, 3)}));
        assert(C2.family() && *C2.family() == wr::gamma({.scale = 3.0, .shape = 2.0}));
    }
    assert(ok);
    std::printf("test_random_algebra : ok\n");
}
