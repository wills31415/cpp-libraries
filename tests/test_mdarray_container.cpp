// wsl::mdarray, conteneur (design/mdarray.md, I13) : conteneur par défaut,
// construction (étendues, valeur, conteneur, littéral, vue, autre tableau),
// sémantique de valeur, accès, vue, parcours à plat, equal.

#include <wsl/mdarray.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace
{
    using E33 = wsl::extents<std::size_t, 3, 3>;
    using D2 = wsl::dims<2>;
    using M33 = wsl::mdarray<double, E33>;
    using M = wsl::mdarray<double, D2>;

    template <auto F>
    concept constant = requires { typename std::bool_constant<(F(), true)>; };

    // Conteneur par défaut : std::array pour des étendues fixes (aucune allocation).
    static_assert(std::is_same_v<M33::container_type, std::array<double, 9>>);
    static_assert(std::is_same_v<M::container_type, std::vector<double>>);
    static_assert(std::is_same_v<wsl::mdarray<int, E33, wsl::layout_left>::container_type, std::array<int, 9>>);
    static_assert(std::is_same_v<wsl::mdarray<int, E33, wsl::layout_stride>::container_type, std::vector<int>>,
                  "disposition non toujours exhaustive : taille inconnue à la compilation");
    static_assert(sizeof(M33) == 9 * sizeof(double), "étendues fixes : seulement les éléments");

    // Pas de conversion implicite d'entiers, de vue ni de littéral en rang 1.
    static_assert(!std::is_convertible_v<int, wsl::mdarray<double, wsl::dims<1>>>);
    static_assert(!std::is_constructible_v<wsl::mdarray<double, wsl::dims<1>>, std::initializer_list<double>>);
    static_assert(std::is_convertible_v<wsl::mdspan<double, D2>, M>, "copie depuis une vue de mêmes étendues : implicite (P1684)");
    static_assert(!std::is_convertible_v<wsl::mdspan<double, D2>, M33> && std::is_constructible_v<M33, wsl::mdspan<double, D2>>,
                  "étendues dynamiques vers fixes : explicite");

    constexpr bool construction()
    {
        M33 g; // zéros
        for (double x : g)
            if (x != 0)
                return false;
        M a(3, 4);
        if (a.extent(0) != 3 || a.extent(1) != 4 || a.size() != 12 || a[2, 3] != 0)
            return false;
        M b(D2(2, 2), 1.5);
        if (b[1, 1] != 1.5)
            return false;
        M c(D2(2, 3), std::vector<double>{1, 2, 3, 4, 5, 6});
        if (c[1, 0] != 4)
            return false;
        // Littéral : étendues dynamiques prises au littéral, fixes vérifiées.
        M lit{{1, 2, 3}, {4, 5, 6}};
        if (lit.extent(0) != 2 || lit.extent(1) != 3 || lit[1, 2] != 6 || !equal(lit, c))
            return false;
        M33 id{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        if (id[2, 2] != 1 || id[0, 1] != 0)
            return false;
        wsl::mdarray<int, wsl::dims<3>> cube{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};
        if (cube[1, 0, 1] != 6 || cube.extent(2) != 2)
            return false;
        // Colonne majeure : même littéral, autre stockage.
        wsl::mdarray<double, D2, wsl::layout_left> left{{1, 2, 3}, {4, 5, 6}};
        if (left[1, 2] != 6 || left.container()[1] != 4 || !equal(left, lit))
            return false;
        // Depuis une vue (copie), depuis un autre tableau.
        wsl::mdarray from_view(lit.to_mdspan());
        static_assert(std::is_same_v<decltype(from_view), M>);
        wsl::mdarray<int, D2> ints(wsl::mdarray<int, wsl::extents<std::size_t, 2, 2>>{{1, 2}, {3, 4}});
        if (!equal(from_view, lit) || ints[1, 0] != 3)
            return false;
        return true;
    }
    static_assert(construction());
    static_assert(!constant<[] { return M{{1, 2}, {3}}; }>, "littéral non rectangulaire");
    static_assert(!constant<[] { return M33{{1, 2}, {3, 4}}; }>, "étendue fixe différente");
    static_assert(!constant<[] { return M(D2(2, 3), std::vector<double>{1, 2}); }>, "conteneur trop petit");

    constexpr bool semantics()
    {
        M a{{1, 2}, {3, 4}};
        M b = a; // copie profonde
        b[0, 0] = 9;
        if (a[0, 0] != 1 || b[0, 0] != 9 || equal(a, b))
            return false;
        // Accès : entiers, std::array, mdposition, at.
        if (a[std::array{1, 1}] != 4 || a[wsl::mdposition{1, 0}] != 3 || a.at(0, 1) != 2)
            return false;
        // Vue sur les mêmes éléments, constance propagée.
        wsl::mdspan<double, D2> v = a;
        v[1, 1] = 7;
        M const& ca = a;
        wsl::mdspan<double const, D2> cv = ca;
        if (a[1, 1] != 7 || cv[1, 1] != 7 || !a.contains(1, 1) || a.contains(2, 0))
            return false;
        // Parcours à plat, remplissage, extraction.
        a.fill(2);
        if (std::accumulate(a.begin(), a.end(), 0.0) != 8)
            return false;
        std::vector<double> storage = std::move(a).extract_container();
        if (storage.size() != 4)
            return false;
        M x{{1, 2}, {3, 4}}, y{{5, 6}, {7, 8}};
        swap(x, y);
        if (x[0, 0] != 5 || y[1, 1] != 4)
            return false;
        // Rang 0, tableau vide.
        wsl::mdarray<int, wsl::extents<std::size_t>> scalar;
        scalar[] = 3;
        M none(0, 4);
        return scalar[] == 3 && none.empty() && none.size() == 0;
    }
    static_assert(semantics());

    void runtime()
    {
        assert(construction());
        assert(semantics());
        M a(3, 4);
        bool thrown = false;
        try
        {
            (void)a.at(3, 0);
        }
        catch (std::out_of_range const&)
        {
            thrown = true;
        }
        assert(thrown);
        // Disposition à pas, conteneur explicite (avec un trou entre les lignes).
        wsl::layout_stride::mapping<D2> m(D2(2, 2), std::array{3, 1});
        wsl::mdarray<int, D2, wsl::layout_stride> s(m, std::vector<int>(5, 1));
        s[1, 1] = 4;
        assert(s.container()[4] == 4 && !s.is_exhaustive());
    }
}

int main()
{
    runtime();
    std::puts("test_mdarray_container : OK");
}
