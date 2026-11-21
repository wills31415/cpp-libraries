// Variables aléatoires (docs/variables.md) : identité, évaluation couplée, constantes.

#include <wsl/random.hpp>

#include <cassert>
#include <cmath>
#include <cstdio>
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
            assert(zs == 2 * xs + ws);
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
    std::printf("test_random_variables : ok\n");
}
