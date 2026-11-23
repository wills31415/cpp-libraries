// Variables aléatoires (docs/variables.md) : identité, évaluation couplée, constantes.

#include <wsl/random.hpp>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <vector>
#include <type_traits>

namespace wr = wsl::random;

template <class Z>
concept has_expectation = requires(Z const& z) { wr::E(z); };

int main()
{
    wr::xoshiro256pp g(20261121);

    // ---- Identité et couplage (V1, V2) ------------------------------------------------------
    {
        wr::variable const X{wr::normal()};
        wr::variable const Y{wr::exponential()};
        auto const X2 = X;                                                  // copie : même variable
        auto const Xi = X.independent_copy();                               // nouvelle identité
        assert(X2.id() == X.id() && Xi.id() != X.id() && Xi.id() > Y.id());
        static_assert(wr::random_variable<decltype(X)> && wr::random_variable<decltype(2 * X + Y)> && !wr::random_variable<wr::normal<>>);
        for (int i = 0; i < 1000; ++i)
            assert((X - X)(g) == 0 && (X - X2)(g) == 0);
        // une variable seule : le tirage de sa loi, même consommation
        wr::xoshiro256pp a(7), b(7);
        assert(X(a) == wr::normal()(b) && a == b);
        // expression : chaque variable de base tirée une fois, par identifiant croissant
        wr::xoshiro256pp c(9), d(9);
        double const z = (2 * X + Y - X)(c);                                // X apparaît deux fois
        double const x = wr::normal()(d), y = wr::exponential()(d);
        assert(z == 2 * x + y - x && c == d);
        wr::xoshiro256pp e(9), f(9);
        double const w = (Y + X)(e);                                        // ordre d'écriture sans effet
        double const x2 = wr::normal()(f), y2 = wr::exponential()(f);
        assert(w == y2 + x2 && e == f);
        // même ω pour plusieurs expressions
        for (int i = 0; i < 100; ++i)
        {
            auto const [xs, zs, ws] = wr::evaluate(g, X, 2 * X + Y, Y);
            assert(std::abs(zs - (2 * xs + ws)) <= 4e-16 * (std::abs(zs) + 1));   // exact contre deux arrondis
        }
        // types : arithmétique usuelle
        wr::variable const K{wr::poisson({.mean = 3.0})};
        static_assert(std::is_same_v<decltype((K + K)(g)), std::int64_t> && std::is_same_v<decltype((K + X)(g)), double>);
        static_assert(std::is_same_v<decltype(X(3, g)), std::vector<double>> && std::is_same_v<decltype((X + 1)(3, g)), std::vector<double>>);
    }

    // ---- Constantes par linéarité ----------------------------------------------------------------
    {
        wr::variable const X{wr::normal({.mean = 1.0, .stddev = 2.0})};
        wr::variable const Y{wr::exponential({.rate = 0.5})};
        auto const Z = 2 * X + Y - 3;
        assert(wr::E(Z) == 2 * 1 + 2 - 3 && wr::Var(Z) == 4 * 4 + 4);
        assert(wr::Var(X + X) == 4 * wr::Var(X) && wr::Var(X - X) == 0 && wr::E(-X / 4) == -0.25);
        assert(wr::Var(X + X.independent_copy()) == 2 * wr::Var(X));
        assert(wr::Cov(X, Z) == 2 * 4 && wr::Cov(Y, X) == 0 && wr::Cov(Z, Z) == wr::Var(Z));
        // moyenne empirique
        double s = 0;
        int const n = 200000;
        for (int i = 0; i < n; ++i)
            s += Z(g);
        assert(std::abs(s / n - wr::E(Z)) < 0.05);
        static_assert(!has_expectation<decltype(X * Y)> && has_expectation<decltype(X * 2)>);   // non affine : pas de linéarité
    }
    // ---- Lot 1b : combinaisons affines exactes (V3) --------------------------------------------
    {
        wr::variable const X{wr::normal({.mean = 0.5, .stddev = 3.0})};
        wr::variable const U{wr::uniform()};
        wr::variable const E{wr::exponential({.scale = 3.0})};
        // identités vraies au bit près : l'arrondi correct de la même valeur exacte
        for (int i = 0; i < 2000; ++i)
        {
            auto const [x, a, b, c, d] = wr::evaluate(g, X, (X + X) / 2, 2 * X - X, (3 * X + 1 - 1) / 3, -(-X));
            assert(a == x && b == x && c == x && d == x);
            auto const [u, e, s] = wr::evaluate(g, U, E, U + E);
            assert(std::abs(s - (u + e)) <= 2.3e-16 * s);
        }
        // loi exacte d'une somme : N(0.5, 9) + N(0.5, 9) indépendantes = N(1, 18)
        auto const Z = X + X.independent_copy();
        wr::normal const N({.mean = 1.0, .variance = 18.0});
        std::vector<double> d;
        for (int i = 0; i < 20000; ++i)
            d.push_back(Z(g));
        std::ranges::sort(d);
        double ks = 0;
        for (std::size_t i = 0; i < d.size(); ++i)
            ks = std::max(ks, std::max(std::abs(N.cdf(d[i]) - double(i) / double(d.size())), std::abs(N.cdf(d[i]) - double(i + 1) / double(d.size()))));
        assert(ks * std::sqrt(double(d.size())) < 1.95);
        // repli flottant : constante non dyadique ((X + 1)/3), loi sans tirage paresseux (gamma)
        wr::variable const G{wr::gamma({.rate = 1.0, .shape = 2.0})};
        for (int i = 0; i < 100; ++i)
        {
            auto const [x, y, gx, h] = wr::evaluate(g, X, (X + 1) / 3, G, G + X);
            assert(std::abs(y - (x + 1) / 3) <= 4e-16 * (std::abs(y) + 1) && std::abs(h - (gx + x)) <= 4e-16 * (std::abs(h) + 1));
        }
        // voie rapide : calcul flottant sur les tirages rapides des lois, mêmes valeurs
        wr::xoshiro256pp a(5), b(5);
        wr::variable const Y{wr::exponential()};
        double const f = (2 * X + Y)(a, wr::fast);
        double const fx = wr::normal({.mean = 0.5, .stddev = 3.0})(b, wr::fast), fy = wr::exponential()(b, wr::fast);
        assert(f == 2 * fx + fy && a == b);
    }
    // ---- Lot 1c : loi d'une expression (V4) ----------------------------------------------------
    {
        wr::variable const X{wr::normal()};
        wr::variable const C{wr::cauchy({.location = 1.0, .scale = 2.0})};
        wr::variable const Y{wr::exponential()};
        static_assert(std::is_same_v<decltype(wr::law(X)), wr::normal<> const&>);       // une variable : sa loi
        // familles closes : la loi exacte, résolue à l'exécution
        auto const A = wr::law(2 * X + 1);
        assert(A.family() && *A.family() == wr::normal({.mean = 1.0, .stddev = 2.0}) && A.cdf(1.0) == 0.5);
        assert(*wr::law(X + X).family() == wr::normal({.variance = 4.0}));                 // 2X, pas P ∗ P
        assert(*wr::law(X + X.independent_copy()).family() == wr::normal({.variance = 2.0}));
        assert(*wr::law(C - 3 * C.independent_copy()).family() == wr::cauchy({.location = -2.0, .scale = 8.0}));
        wr::variable const G{wr::gamma({.rate = 2.0, .shape = 1.5})};
        assert(*wr::law(G + G.independent_copy()).family() == wr::gamma({.rate = 2.0, .shape = 3.0}));
        assert(*wr::law(3 * G).family() == wr::gamma({.rate = wr::ratio(2, 3), .shape = 1.5}));
        assert(!wr::law(G - G.independent_copy()).family() && std::isnan(wr::law(G - G.independent_copy()).cdf(0.0)));   // pas d'identité
        assert(!wr::law(X - X).family());                                                  // loi de Dirac en 0
        // sans famille : moments par les cumulants, tirage de l'expression
        auto const M = wr::law(2 * X + Y);
        static_assert(!decltype(M)::has_family && decltype(M)::sampling == wr::exactness::exact);
        assert(M.mean() == 1 && M.variance() == 5 && std::abs(M.skewness() - 2 / (5 * std::sqrt(5.0))) < 1e-15 && std::abs(M.excess_kurtosis() - 6.0 / 25) < 1e-15);
        assert(wr::law(X + X).excess_kurtosis() == 0 && wr::law(Y - Y.independent_copy()).skewness() == 0);
        wr::xoshiro256pp a(3), b(3);
        assert(M(a) == (2 * X + Y)(b) && a == b);                                          // tirage de l'expression
    }

    // ---- Lot 1d : fonctions de variables --------------------------------------------------------
    {
        namespace sp = wsl::special;
        wr::variable const X{wr::normal({.mean = 0.2, .stddev = 0.7})};
        wr::variable const U{wr::uniform()};
        for (int i = 0; i < 1000; ++i)
        {
            auto const [x, u, e, h, m, q, a, f] = wr::evaluate(g, X, U, exp(X), hypot(X, U), max(X, U), pow(X, 2.0), wr::abs(X), wr::apply([](double s, double t) { return s - t; }, X, U));
            assert(e == sp::exp(x) && h == sp::hypot(x, u) && m == std::max(x, u) && q == sp::pow(x, 2.0) && a == std::abs(x) && f == x - u);
        }
        // loi de f(X) : mesure image ; f monotone déclarée → répartition, quantiles, densité
        auto const LN = wr::law(exp(X));
        wr::lognormal const R({.log_mean = 0.2, .log_stddev = 0.7});
        for (double y : {0.3, 1.0, 2.5})
            assert(std::abs(LN.cdf(y) - R.cdf(y)) < 1e-15 && std::abs(LN.pdf(y) - R.pdf(y)) < 1e-14 * R.pdf(y));
        assert(std::abs(LN.quantile(0.9) - R.quantile(0.9)) < 1e-14 * R.quantile(0.9) && LN.median() == sp::exp(0.2));
        auto const S = wr::law(sin(X));                                     // fonction quelconque : tirage seul
        static_assert(!wr::has_cdf<decltype(S)> && wr::has_cdf<decltype(LN)>);
        auto const T = wr::law(sqrt(U));                                    // √U : répartition y² sur [0, 1]
        assert(T.cdf(0.5) == 0.25 && T.quantile(0.25) == 0.5);
        static_assert(!has_expectation<decltype(exp(X))>);                   // non affine : pas de linéarité
        // empreinte des fonctions
        wr::xoshiro256pp x(77);
        std::uint64_t hash = 0xcbf29ce484222325;
        for (int i = 0; i < 4096; ++i)
        {
            auto const [a1, a2, a3] = wr::evaluate(x, exp(X) * U + log1p(U), atan2(X, U) - erf(X), cbrt(X) / (1 + abs(X)));
            for (double v : {a1, a2, a3})
                hash = (hash ^ std::bit_cast<std::uint64_t>(v)) * 0x100000001b3;
        }
        std::printf("variables functions    %016llx\n", static_cast<unsigned long long>(hash));
        assert(hash == 0x93f4e7f5ae31e15a);
    }

    // ---- Lot 2 : événements, probabilités, conditionnement --------------------------------------
    {
        wr::variable const X{wr::normal()};
        wr::variable const Y{wr::normal({.mean = 1.0, .stddev = 2.0})};
        wr::variable const U{wr::uniform()};
        wr::variable const K{wr::poisson({.mean = 3.0})};
        // événements : indicatrices, évaluées sur le même ω
        for (int i = 0; i < 2000; ++i)
        {
            auto const [x, u, a, b, c, d] = wr::evaluate(g, X, U, X > 0, X - X > 0, X - X == 0, (X > 0) && !(U > 0.5));
            assert(a == (x > 0) && !b && c && d == (x > 0 && !(u > 0.5)));
        }
        static_assert(std::is_same_v<decltype((X > 0)(g)), bool> && wr::random_event<decltype(X > Y || U < 0.1)>);
        // probabilités exactes (E1)
        assert(wr::P(X > 0) == 0.5 && wr::P(0.0 < X) == 0.5 && wr::P(X == 1.0) == 0 && wr::P(X != 1.0) == 1);
        assert(wr::P(Y <= 1.0) == 0.5 && wr::P(X > Y) == wr::normal({.mean = -1.0, .variance = 5.0}).ccdf(0.0));   // X − Y ~ N(−1, 5)
        assert(wr::P(K >= 3) == wr::poisson({.mean = 3.0}).ccdf(2) && wr::P(K > 2.5) == wr::P(K >= 3) && wr::P(K == 2) == wr::poisson({.mean = 3.0}).pmf(2));
        assert(wr::P(K == 2.5) == 0 && wr::P(K < 0) == 0 && wr::P(K <= -1.0) == 0 && wr::P(K >= 0) == 1);
        assert(std::abs(wr::P((X > 0) && (U < 0.25)) - 0.125) < 1e-16 && std::abs(wr::P(!(X > 0) || (U < 0.25)) - 0.625) < 1e-16);
        assert(std::isnan(wr::P((X > 0) && (X < 1))) && std::isnan(wr::P(X * U > 0)));   // même variable, non affine : NaN
        // conditionnement (E2) : nouvelle variable, loi conditionnelle
        auto const T = wr::given(X, X > 1.0);
        assert(T.id() > U.id());
        wr::truncated const Tr = wr::truncate(wr::normal(), {.lower = 1.0});
        for (double t : {1.1, 1.5, 3.0})
            assert(wr::law(T).cdf(t) == Tr.cdf(t) && wr::law(T).pdf(t) == Tr.pdf(t));
        assert(wr::law(T).quantile(0.5) == Tr.quantile(0.5));
        for (int i = 0; i < 2000; ++i)
            assert(T(g) > 1.0);
        auto const I = wr::given(X, U > 0.5);                               // indépendant : la loi de X
        assert(wr::law(I).cdf(0.3) == wr::normal().cdf(0.3));
        auto const Mx = wr::given(X + U, U > 0.5);                          // sans identité : tirage seul
        double m = 0;
        for (int i = 0; i < 40000; ++i)
            m += Mx(g);
        assert(std::abs(m / 40000 - 0.75) < 0.03);                          // E[X + U | U > ½] = 0,75
        bool thrown = false;
        try { (void)wr::given(X, X > INFINITY); } catch (std::domain_error const&) { thrown = true; }
        assert(thrown);
    }

    // ---- Lot 3 : variables vectorielles, processus comme familles de variables --------------------
    {
        wr::variable const V{wr::multivariate_normal<double, 2>({.mean = {1.0, 2.0}, .covariance = {{2.0, 0.6}, {0.6, 1.0}}})};
        wr::variable const X{wr::normal()};
        // composantes : même tirage, couplées
        for (int i = 0; i < 1000; ++i)
        {
            auto const [v, a, b, c] = wr::evaluate(g, V, V[0], V[1], V[0] - 2 * V[1] + X);
            assert(a == v[0] && b == v[1] && std::abs(c - (v[0] - 2 * v[1])) < 20);
        }
        // linéarité avec covariances
        assert(wr::E(V[0] + 2 * V[1]) == 5 && std::abs(wr::Var(V[0] + 2 * V[1]) - 8.4) < 1e-15 && wr::Cov(V[0], V[1]) == 0.6);
        assert(wr::Cov(V[0], X) == 0 && wr::Var(V[1] - V[1]) == 0 && std::abs(wr::Cov(V[0] + V[1], V[1] + X) - 1.6) < 1e-15);
        // famille gaussienne : la loi d'une combinaison de composantes
        auto const D = wr::law(V[0] - V[1] + X);
        assert(D.family() && D.family()->mean() == -1 && std::abs(D.family()->variance() - 2.8) < 1e-15);
        assert(std::abs(wr::P(V[0] > V[1]) - wr::normal({.mean = -1.0, .variance = 1.8}).ccdf(0.0)) < 1e-15);   // V₀ − V₁ ~ N(−1, 1,8)
        // processus : le mouvement brownien sur une grille, ses accroissements
        wr::brownian_motion const B({.drift = 0.5, .volatility = 2.0});
        wr::variable const W{B.on_grid({1.0, 2.0, 4.0})};
        assert(wr::E(W[2] - W[0]) == 1.5 && wr::Var(W[2] - W[1]) == 8 && wr::Cov(W[0], W[2]) == 4);
        auto const I = wr::law(W[2] - W[0]);                               // accroissement : N(μΔt, σ²Δt)
        assert(I.family() && I.family()->mean() == 1.5 && I.family()->variance() == 12);
        double m = 0, q = 0;
        int const n = 40000;
        for (int i = 0; i < n; ++i)
        {
            auto const [w, d] = wr::evaluate(g, W, W[2] - W[0]);
            assert(d == w[2] - w[0]);
            m += d;
            q += d * d;
        }
        assert(std::abs(m / n - 1.5) < 0.06 && std::abs(q / n - (12 + 2.25)) < 0.4);
    }

    // ---- Empreinte (bit-exacte sur toutes les plateformes) -----------------------------------
    {
        wr::xoshiro256pp x(20261121);
        wr::variable const N{wr::normal({.mean = 1.0, .variance = 2.0})};
        wr::variable const L{wr::laplace({.location = -1.0, .scale = 0.5})};
        wr::variable const R{wr::exponential({.rate = wr::ratio(3, 5)})};
        wr::variable const P{wr::poisson({.mean = 4.0})};
        std::uint64_t h = 0xcbf29ce484222325;
        auto const add = [&](double v) { h = (h ^ std::bit_cast<std::uint64_t>(v)) * 0x100000001b3; };
        for (int i = 0; i < 4096; ++i)
        {
            auto const [a, b, c, d] = wr::evaluate(x, 2 * N - L / 4 + R, N + N.independent_copy(), N * L + P, (R - N) / 3);
            add(a);
            add(b);
            add(c);
            add(d);
        }
        std::printf("variables              %016llx\n", static_cast<unsigned long long>(h));
        assert(h == 0x0949c24e8827fb5f);
    }
    // ---- Empreinte des événements et du conditionnement ---------------------------------------
    {
        wr::xoshiro256pp x(4242);
        wr::variable const N{wr::normal({.mean = 0.5, .variance = 3.0})};
        wr::variable const E{wr::exponential({.scale = 2.0})};
        auto const C = wr::given(N, N > E);
        std::uint64_t h = 0xcbf29ce484222325;
        for (int i = 0; i < 4096; ++i)
        {
            auto const [a, b, c] = wr::evaluate(x, N - E > 1, (N > 0) || (E < 0.5), C);
            h = (h ^ (std::uint64_t{a} | std::uint64_t{b} << 1)) * 0x100000001b3;
            h = (h ^ std::bit_cast<std::uint64_t>(c)) * 0x100000001b3;
        }
        std::printf("variables events       %016llx\n", static_cast<unsigned long long>(h));
        assert(h == 0xd86c1f4992ee39c0);
    }
    // ---- Empreinte des variables vectorielles ---------------------------------------------------
    {
        wr::xoshiro256pp x(777);
        wr::variable const V{wr::multivariate_normal<double, 3>({.mean = {0.0, 1.0, -1.0}, .covariance = {{1.0, 0.2, 0.0}, {0.2, 2.0, 0.5}, {0.0, 0.5, 1.5}}})};
        wr::variable const W{wr::brownian_motion({.drift = -0.25, .volatility = 1.5}).on_grid({0.5, 1.0, 3.0})};
        std::uint64_t h = 0xcbf29ce484222325;
        for (int i = 0; i < 2048; ++i)
        {
            auto const [a, b, c] = wr::evaluate(x, V[0] + 2 * V[2], W[2] - W[1], V[1] * W[0] > 0);
            h = (h ^ std::bit_cast<std::uint64_t>(a)) * 0x100000001b3;
            h = (h ^ std::bit_cast<std::uint64_t>(b)) * 0x100000001b3;
            h = (h ^ std::uint64_t{c}) * 0x100000001b3;
        }
        std::printf("variables vectors      %016llx\n", static_cast<unsigned long long>(h));
        assert(h == 0x7d67eaa28ed5efff);
    }
    std::printf("test_random_variables : ok\n");
}
