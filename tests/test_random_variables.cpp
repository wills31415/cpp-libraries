// Variables aléatoires (docs/variables.md) : identité, évaluation couplée, constantes.

#include <wsl/random.hpp>

#include <cassert>
#include <cmath>
#include <cstdio>
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
    std::printf("test_random_variables : ok\n");
}
