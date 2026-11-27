// wsl::mdarray, monde des tableaux entiers (design/mdarray.md, V8, I6, I11, I12, I13) :
// opérations terme à terme, diffusion des scalaires, fonctions, vectorize, masques,
// where, réductions, affectation et construction, aliasing, expressions temporaires.

#include <wsl/mdarray.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <type_traits>

namespace
{
    using D2 = wsl::dims<2>;
    using M = wsl::mdarray<double, D2>;
    using M23 = wsl::mdarray<double, wsl::extents<std::size_t, 2, 3>>;

    template <auto F>
    concept constant = requires { typename std::bool_constant<(F(), true)>; };

    template <class A, class B>
    concept addable = requires(A a, B b) { a + b; };

    // Refus (P2) : rangs ou étendues fixes différents (erreurs de compilation, avec
    // message : tests/mdarray/compile_fail/expression_array.cpp).
    static_assert(addable<M23, M> && addable<M, double> && addable<double, M>);
    static_assert(!addable<wsl::mdposition<2>, M>, "une position n'est pas un scalaire");
    static_assert(!std::is_convertible_v<decltype(std::declval<M>() < std::declval<M>()), bool>, "un masque n'est pas un bool");

    constexpr bool arithmetic()
    {
        M a{{1, 2, 3}, {4, 5, 6}};
        M23 b{{6, 5, 4}, {3, 2, 1}};
        M c = a + b;
        if (!equal(c, M{{7, 7, 7}, {7, 7, 7}}))
            return false;
        wsl::mdarray d = a * b - 1; // déduction : étendues fixes l'emportent
        static_assert(std::is_same_v<decltype(d)::extents_type, wsl::extents<std::size_t, 2, 3>>);
        if (d[0, 0] != 5 || d[1, 2] != 5)
            return false;
        M e(D2(2, 3));
        e = 2 * a + 1;
        if (!equal(e, M{{3, 5, 7}, {9, 11, 13}}))
            return false;
        e = -a / 2;
        if (e[1, 1] != -2.5)
            return false;
        // Un élément d'une expression.
        if ((a + b)[1, 2] != 7 || (a * a)[wsl::mdposition{1, 0}] != 16)
            return false;
        // Affectations composées.
        e = a;
        e += b;
        e *= 2;
        e -= a;
        if (!equal(e, M{{13, 12, 11}, {10, 9, 8}}))
            return false;
        // vectorize (les fonctions de <cmath> ne sont pas constexpr en C++23 : à l'exécution).
        M x{{0, 1}, {4, 9}};
        auto F = wsl::vectorize([](double u, double v) { return u * 10 + v; });
        M g = F(x, x);
        if (g[1, 1] != 99 || F(1.0, 2.0) != 12)
            return false;
        M h = max(x - 3, 0); // ReLU
        if (!equal(h, M{{0, 0}, {1, 6}}))
            return false;
        // Masques, where.
        auto m = x > 1;
        (void)m;
        if (count(x > 1) != 2 || !any(x > 8) || all(x > 0) || !all(x >= 0))
            return false;
        M w = where(x > 1, x, -1);
        if (!equal(w, M{{-1, -1}, {4, 9}}))
            return false;
        // Réductions de tableaux entiers.
        if (sum(a) != 21 || prod(a) != 720 || min(a) != 1 || max(b) != 6)
            return false;
        if (argmax(a) != wsl::mdposition{1, 2} || argmin(b) != wsl::mdposition{1, 2} || sum(a * a) != 91)
            return false;
        return true;
    }
    static_assert(arithmetic());
    static_assert(!constant<[] {
        M a(2, 3), b(3, 2);
        M c = a + b;
        return c[0, 0];
    }>, "formes différentes");

    void runtime()
    {
        assert(arithmetic());

        // Fonctions de <cmath>, élément par élément (I11).
        M x0{{0, 1}, {4, 9}};
        M r = sqrt(x0) + abs(-x0);
        assert((equal(r, M{{0, 2}, {6, 12}})));
        M p = pow(x0, 2.0) + atan2(x0, 1.0) * 0;
        assert((p[1, 1] == 81));

        // NaN propagé (I6), fmin l'ignore.
        double const nan = std::numeric_limits<double>::quiet_NaN();
        M x{{1, nan}, {3, 4}};
        assert(std::isnan(max(x)) && std::isnan(min(x)));
        assert((argmax(x) == wsl::mdposition{0, 1}));
        M f = fmin(x, 2.0);
        assert((f[0, 1] == 2));

        // Somme par paires (I9) : 1 + 2⁻⁵³ × n termes, ordre fixé.
        M big(1, 1000);
        big.fill(1.0);
        assert(sum(big) == 1000);

        // Aliasing (I12) : la vue transposée du même stockage oblige à un temporaire.
        M a{{1, 2}, {3, 4}};
        wsl::mdspan<double const, D2, wsl::layout_left> t(a.data(), 2, 2); // aᵀ
        a = a + t;                                                         // a + aᵀ
        assert((equal(a, M{{2, 5}, {5, 8}})));
        // Mêmes positions : sur place.
        a = a * 2;
        assert((a[1, 1] == 16));

        // std::move : l'échappatoire explicite (L2).
        auto e = a + 1;
        M y = std::move(e);
        assert((y[0, 0] == 5));
    }
}

int main()
{
    runtime();
    std::puts("mdarray.expression_array : OK");
}
