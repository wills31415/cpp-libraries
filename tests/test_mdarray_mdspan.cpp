// wsl::mdarray, étendues, dispositions et vue (design/mdarray.md, A5, I2, D19 à D23) :
// sémantique de <mdspan>, comparée à std::mdspan là où la bibliothèque standard
// l'offre (libc++) ; ajouts : contains, wrap, clamp, m[p], at.

#include <wsl/mdarray.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <algorithm>
#include <cstdio>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace
{
    using wsl::dynamic_extent;

    template <auto F>
    concept constant = requires { typename std::bool_constant<(F(), true)>; };

    template <class M, class I>
    concept indexable = requires(M m, I i) { m[i]; };
    template <class M, class I, class J>
    concept indexable2 = requires(M m, I i, J j) { m[i, j]; };

    using E34 = wsl::extents<std::size_t, 3, 4>;
    using E3d = wsl::extents<std::size_t, 3, dynamic_extent>;
    using D2 = wsl::dims<2>;

    // ---------------------------------------------------------------------
    // Étendues
    // ---------------------------------------------------------------------
    static_assert(E3d::rank() == 2 && E3d::rank_dynamic() == 1 && E3d::static_extent(0) == 3);
    static_assert(E3d::static_extent(1) == dynamic_extent);
    static_assert(std::is_same_v<wsl::dims<3>, wsl::extents<std::size_t, dynamic_extent, dynamic_extent, dynamic_extent>>);
    static_assert(std::is_same_v<wsl::dims<2, int>, wsl::dextents<int, 2>>);
    static_assert(std::is_same_v<decltype(wsl::extents(3, 4)), D2>);
    static_assert(E3d(4).extent(0) == 3 && E3d(4).extent(1) == 4 && E3d(3, 4).extent(1) == 4);
    static_assert(E3d(4) == E34{} && E34{} == D2(3, 4) && !(D2(3, 4) == D2(4, 3)));
    static_assert(!(wsl::extents<std::size_t, 3>{} == E34{}), "rangs différents");
    static_assert(wsl::extents<int>::rank() == 0 && wsl::extents<int>{} == wsl::extents<std::size_t>{});
    static_assert(!constant<[] { return E3d(2, 4); }>, "étendue fixe différente");
    static_assert(!constant<[] { return D2(-1, 4); }>, "étendue négative");
    static_assert(!constant<[] { return E34::static_extent(2); }>, "rang hors limites");

    // Appartenance (D19).
    static_assert(E34{}.contains(2, 3) && !E34{}.contains(3, 0) && !E34{}.contains(-1, 0));
    static_assert(E34{}.contains(wsl::mdposition{2, 3}) && !E34{}.contains(wsl::mdposition{0, 4}));
    static_assert(E34{}.contains(std::array{2, 3}) && !E34{}.contains(std::array{2, -1}));
    static_assert(E34{}.contains(wsl::mdposition{0, 0}, wsl::mddisplacement{2, 3}));
    static_assert(!E34{}.contains(wsl::mdposition{0, 0}, wsl::mddisplacement{-1, 0}), "sans former une position négative");
    static_assert(!E34{}.contains(wsl::mdposition{2, 3}, wsl::mddisplacement{1, 0}));
    static_assert(E34{}.contains(wsl::mdposition{2, 3}, wsl::mddisplacement{-2, -3}));
    static_assert(!D2(0, 4).contains(0, 0), "ensemble vide");

    // Bords (D22, D23).
    static_assert(wrap(wsl::mdposition{0, 0}, wsl::mddisplacement{-1, -1}, E34{}) == wsl::mdposition{2, 3});
    static_assert(wrap(wsl::mdposition{2, 3}, wsl::mddisplacement{1, 1}, E34{}) == wsl::mdposition{0, 0});
    static_assert(wrap(wsl::mdposition{1, 1}, wsl::mddisplacement{-7, 9}, E34{}) == wsl::mdposition{0, 2});
    static_assert(wrap(wsl::mdposition{5, 9}, E34{}) == wsl::mdposition{2, 1});
    static_assert(clamp(wsl::mdposition{0, 0}, wsl::mddisplacement{-1, -5}, E34{}) == wsl::mdposition{0, 0});
    static_assert(clamp(wsl::mdposition{2, 3}, wsl::mddisplacement{1, 7}, E34{}) == wsl::mdposition{2, 3});
    static_assert(clamp(wsl::mdposition{1, 1}, wsl::mddisplacement{1, -1}, E34{}) == wsl::mdposition{2, 0});
    static_assert(clamp(wsl::mdposition{7, 1}, E34{}) == wsl::mdposition{2, 1});
    static_assert(!constant<[] { return wrap(wsl::mdposition{0, 0}, D2(0, 4)); }>, "étendue nulle");

    // Toutes les valeurs de wrap et clamp, contre une définition directe.
    constexpr bool borders(std::ptrdiff_t range)
    {
        constexpr E34 e{};
        for (std::ptrdiff_t a = -range; a <= range; ++a)
            for (std::ptrdiff_t b = -range; b <= range; ++b)
                for (std::size_t i = 0; i < 3; ++i)
                    for (std::size_t j = 0; j < 4; ++j)
                    {
                        wsl::mdposition p{i, j};
                        wsl::mddisplacement d{a, b};
                        auto const wi = wsl::rem_euclid(static_cast<std::ptrdiff_t>(i) + a, std::ptrdiff_t{3});
                        auto const wj = wsl::rem_euclid(static_cast<std::ptrdiff_t>(j) + b, std::ptrdiff_t{4});
                        if (wrap(p, d, e) != wsl::mdposition{wi, wj})
                            return false;
                        auto const ci = std::clamp<std::ptrdiff_t>(static_cast<std::ptrdiff_t>(i) + a, 0, 2);
                        auto const cj = std::clamp<std::ptrdiff_t>(static_cast<std::ptrdiff_t>(j) + b, 0, 3);
                        if (clamp(p, d, e) != wsl::mdposition{ci, cj})
                            return false;
                        bool const inside = static_cast<std::ptrdiff_t>(i) + a >= 0 && static_cast<std::ptrdiff_t>(i) + a < 3
                                            && static_cast<std::ptrdiff_t>(j) + b >= 0 && static_cast<std::ptrdiff_t>(j) + b < 4;
                        if (e.contains(p, d) != inside)
                            return false;
                    }
        return true;
    }
    static_assert(borders(3));

    // ---------------------------------------------------------------------
    // Dispositions
    // ---------------------------------------------------------------------
    constexpr bool layouts()
    {
        wsl::layout_right::mapping<E34> r{};
        wsl::layout_left::mapping<E34> l{};
        if (r.required_span_size() != 12 || l.required_span_size() != 12)
            return false;
        if (r(1, 2) != 6 || l(1, 2) != 7 || r.stride(0) != 4 || r.stride(1) != 1 || l.stride(0) != 1 || l.stride(1) != 3)
            return false;
        wsl::layout_stride::mapping<E34> s(E34{}, std::array{1, 6});
        if (s(2, 3) != 20 || s.required_span_size() != 21 || s.is_exhaustive())
            return false;
        wsl::layout_stride::mapping<E34> t = r; // implicite : extents convertibles
        if (!(t == r) || !t.is_exhaustive() || t.strides() != std::array<std::size_t, 2>{4, 1})
            return false;
        wsl::layout_right::mapping<E34> back(t); // explicite depuis layout_stride
        if (!(back == r))
            return false;
        wsl::layout_right::mapping<D2> dyn = r; // fixe → dynamique, implicite
        if (dyn.extents().extent(1) != 4)
            return false;
        // Rang 0 : un seul élément.
        wsl::layout_right::mapping<wsl::extents<std::size_t>> z{};
        if (z.required_span_size() != 1 || z() != 0)
            return false;
        // Étendue nulle.
        wsl::layout_stride::mapping<D2> empty(D2(0, 4), std::array{4, 1});
        if (empty.required_span_size() != 0 || !empty.is_exhaustive())
            return false;
        // Longueur 1 et pas quelconque : exhaustif (aucune place ne manque).
        wsl::layout_stride::mapping<D2> one(D2(1, 3), std::array{100, 1});
        if (!one.is_exhaustive())
            return false;
        return true;
    }
    static_assert(layouts());
    static_assert(!std::is_convertible_v<wsl::layout_stride::mapping<E34>, wsl::layout_right::mapping<E34>>);
    static_assert(std::is_convertible_v<wsl::layout_right::mapping<E34>, wsl::layout_stride::mapping<E34>>);
    static_assert(!std::is_convertible_v<wsl::layout_right::mapping<D2>, wsl::layout_right::mapping<E34>>);
    static_assert(!constant<[] { return wsl::layout_stride::mapping<E34>(E34{}, std::array{1, 1}); }>, "pas non emboîtés");
    static_assert(!constant<[] { return wsl::layout_stride::mapping<E34>(E34{}, std::array{0, 1}); }>, "pas nul");
    static_assert(!constant<[] { return wsl::layout_right::mapping<E34>{}(3, 0); }>, "position hors de l'ensemble");

    // ---------------------------------------------------------------------
    // Vue
    // ---------------------------------------------------------------------
    static_assert(std::is_same_v<decltype(wsl::mdspan(static_cast<int*>(nullptr), 3, 4)), wsl::mdspan<int, D2>>);
    static_assert(!indexable<wsl::mdspan<int, D2>, int>, "rang différent");
    static_assert(!indexable<wsl::mdspan<int, D2>, wsl::mddisplacement<2>>, "un déplacement n'indexe pas (D21)");
    static_assert(indexable<wsl::mdspan<int, D2>, wsl::mdposition<2>> && indexable2<wsl::mdspan<int, D2>, int, long>);
    static_assert(std::is_convertible_v<wsl::mdspan<int, E34>, wsl::mdspan<int const, D2>>);
    static_assert(!std::is_convertible_v<wsl::mdspan<int const, D2>, wsl::mdspan<int, D2>>);
    static_assert(!std::is_convertible_v<wsl::mdspan<int, D2>, wsl::mdspan<int, E34>>);

    constexpr bool views()
    {
        std::array<int, 12> data{};
        for (int k = 0; k < 12; ++k)
            data[k] = k;
        wsl::mdspan m(data.data(), 3, 4);
        if (m.rank() != 2 || m.extent(0) != 3 || m.size() != 12 || m.empty())
            return false;
        if (m[1, 2] != 6 || m[std::array{2, 3}] != 11 || m[wsl::mdposition{2, 0}] != 8)
            return false;
        std::array<std::size_t, 2> idx{1, 1};
        if (m[std::span(idx)] != 5)
            return false;
        m[0, 0] = 42;
        if (data[0] != 42)
            return false;
        wsl::mdspan<int const, D2> c = m;
        if (c[0, 0] != 42 || c.data_handle() != data.data())
            return false;
        wsl::mdspan<int, E34, wsl::layout_left> left(data.data());
        if (left[1, 2] != data[7])
            return false;
        wsl::mdspan<int, D2, wsl::layout_stride> s(data.data(), wsl::layout_stride::mapping<D2>(D2(2, 2), std::array{4, 1}));
        if (s[1, 1] != 5 || s.is_exhaustive())
            return false;
        if (!m.contains(2, 3) || m.contains(3, 0) || !m.contains(wsl::mdposition{0, 3}))
            return false;
        int scalar = 7;
        wsl::mdspan<int, wsl::extents<std::size_t>> zero(&scalar);
        if (zero[] != 7 || zero.size() != 1)
            return false;
        wsl::mdspan<int, D2> none(data.data(), 0, 4);
        if (!none.empty() || none.size() != 0)
            return false;
        int carray[5] = {1, 2, 3, 4, 5};
        wsl::mdspan v(carray);
        static_assert(std::is_same_v<decltype(v), wsl::mdspan<int, wsl::extents<std::size_t, 5>>>);
        return v[4] == 5;
    }
    static_assert(views());
    static_assert(!constant<[] {
        std::array<int, 4> d{};
        return wsl::mdspan(d.data(), 2, 2).at(2, 0);
    }>, "at hors du tableau");

    void runtime()
    {
        assert(views());
        assert(borders(9));
        std::vector<double> data(12, 0.0);
        wsl::mdspan m(data.data(), 3, 4);
        bool thrown = false;
        try
        {
            (void)m.at(3, 0);
        }
        catch (std::out_of_range const&)
        {
            thrown = true;
        }
        assert(thrown);
        thrown = false;
        try
        {
            (void)m.at(-1, 0);
        }
        catch (std::out_of_range const&)
        {
            thrown = true;
        }
        assert(thrown);
        m.at(wsl::mdposition{2, 3}) = 1.5;
        assert(data[11] == 1.5);
    }

#if defined(__cpp_lib_mdspan)
    // ---------------------------------------------------------------------
    // Comparaison avec std::mdspan (même sémantique, A5)
    // ---------------------------------------------------------------------
    struct W
    {
        template <class I, std::size_t... E>
        using ext = wsl::extents<I, E...>;
        using right = wsl::layout_right;
        using left = wsl::layout_left;
        using stride = wsl::layout_stride;
        template <class T, class E, class L = right>
        using span = wsl::mdspan<T, E, L>;
    };
    struct S
    {
        template <class I, std::size_t... E>
        using ext = std::extents<I, E...>;
        using right = std::layout_right;
        using left = std::layout_left;
        using stride = std::layout_stride;
        template <class T, class E, class L = right>
        using span = std::mdspan<T, E, L>;
    };

    // Constructibilité et conversions implicites, mêmes réponses des deux côtés.
    template <class N>
    constexpr auto traits()
    {
        using e3 = typename N::template ext<std::size_t, 3>;
        using ed = typename N::template ext<std::size_t, dynamic_extent>;
        using ei = typename N::template ext<int, dynamic_extent>;
        using e34 = typename N::template ext<std::size_t, 3, 4>;
        using ed2 = typename N::template ext<std::size_t, dynamic_extent, dynamic_extent>;
        using r34 = typename N::right::template mapping<e34>;
        using l34 = typename N::left::template mapping<e34>;
        using s34 = typename N::stride::template mapping<e34>;
        using r3 = typename N::right::template mapping<e3>;
        using l3 = typename N::left::template mapping<e3>;
        return std::array{
            std::is_convertible_v<e3, ed>, std::is_convertible_v<ed, e3>, std::is_constructible_v<e3, ed>,
            std::is_convertible_v<ei, ed>, std::is_convertible_v<ed, ei>, std::is_constructible_v<ed, ei>,
            std::is_convertible_v<int, ed>, std::is_constructible_v<ed, int>, std::is_constructible_v<e34, int>,
            std::is_constructible_v<e34, int, int>, std::is_constructible_v<ed2, int>,
            std::is_convertible_v<std::array<int, 1>, ed>, std::is_convertible_v<std::array<int, 2>, e34>,
            std::is_constructible_v<e34, std::array<int, 2>>, std::is_convertible_v<r34, s34>, std::is_convertible_v<s34, r34>,
            std::is_constructible_v<r34, s34>, std::is_convertible_v<l34, r34>, std::is_constructible_v<l34, r34>,
            std::is_convertible_v<l3, r3>, std::is_convertible_v<r3, l3>,
            std::is_convertible_v<typename N::template span<int, e34>, typename N::template span<int const, ed2>>,
            std::is_convertible_v<typename N::template span<int, ed2>, typename N::template span<int, e34>>,
            std::is_constructible_v<typename N::template span<int, ed2>, int*, int, int>,
            std::is_convertible_v<typename N::template span<int, e34, typename N::left>, typename N::template span<int, e34>>,
            std::is_default_constructible_v<typename N::template span<int, e34>>,
            std::is_default_constructible_v<typename N::template span<int, ed2>>,
        };
    }
    static_assert(traits<W>() == traits<S>());

    // Mêmes places en mémoire, mêmes tailles, mêmes propriétés.
    template <class WMap, class SMap>
    bool same_mapping(WMap const& w, SMap const& s)
    {
        if (w.required_span_size() != s.required_span_size())
            return false;
        // Ensemble vide : exhaustif selon la norme ([mdspan.layout.stride.obs] 6.1) ;
        // libc++ suit une formulation antérieure, on ne compare pas.
        bool const empty = w.extents().extent(0) == 0 || w.extents().extent(1) == 0 || w.extents().extent(2) == 0;
        if (empty ? !w.is_exhaustive() : w.is_exhaustive() != s.is_exhaustive())
            return false;
        for (std::size_t i = 0; i < w.extents().extent(0); ++i)
            for (std::size_t j = 0; j < w.extents().extent(1); ++j)
                for (std::size_t k = 0; k < w.extents().extent(2); ++k)
                    if (w(i, j, k) != s(i, j, k))
                        return false;
        for (std::size_t r = 0; r < 3; ++r)
            if (w.stride(r) != s.stride(r))
                return false;
        return true;
    }

    void compare_with_std()
    {
        using WE = wsl::dims<3>;
        using SE = std::dextents<std::size_t, 3>;
        for (std::size_t a = 0; a <= 3; ++a)
            for (std::size_t b = 1; b <= 3; ++b)
                for (std::size_t c = 1; c <= 4; ++c)
                {
                    assert(same_mapping(wsl::layout_right::mapping<WE>(WE(a, b, c)), std::layout_right::mapping<SE>(SE(a, b, c))));
                    assert(same_mapping(wsl::layout_left::mapping<WE>(WE(a, b, c)), std::layout_left::mapping<SE>(SE(a, b, c))));
                    // Pas emboîtés, dans divers ordres et avec des trous.
                    std::array<std::array<std::size_t, 3>, 4> strides{{{b * c, c, 1},
                                                                       {1, a ? a : 1, (a ? a : 1) * b},
                                                                       {2 * b * c + 1, c + 1, 1},
                                                                       {b, 1, 3 * (a ? a : 1) * b}}};
                    for (auto const& s : strides)
                        assert(same_mapping(wsl::layout_stride::mapping<WE>(WE(a, b, c), s),
                                            std::layout_stride::mapping<SE>(SE(a, b, c), s)));
                }

        // Conversions dans les deux sens, sans copie.
        std::array<int, 12> data{};
        std::mdspan<int, std::dextents<std::size_t, 2>> s(data.data(), 3, 4);
        wsl::mdspan w = s;
        static_assert(std::is_same_v<decltype(w), wsl::mdspan<int, wsl::dims<2>>>);
        w[2, 1] = 9;
        assert((s[2, 1] == 9));
        std::mdspan<int, std::dextents<std::size_t, 2>> back = w;
        assert(back.data_handle() == data.data() && back.extent(1) == 4);
        std::mdspan<int, std::dextents<std::size_t, 2>, std::layout_stride> strided = w;
        assert(strided.stride(0) == 4 && strided.stride(1) == 1);
        wsl::mdspan<int, wsl::dims<2>, wsl::layout_left> l(data.data(), 3, 4);
        std::mdspan<int, std::dextents<std::size_t, 2>, std::layout_left> sl = l;
        assert((sl[1, 2] == l[1, 2]));
        static_assert(!std::is_convertible_v<wsl::mdspan<int, wsl::dims<2>, wsl::layout_left>, std::mdspan<int, std::dextents<std::size_t, 2>>>,
                      "pas de changement de disposition");
    }
#endif
}

int main()
{
    runtime();
#if defined(__cpp_lib_mdspan)
    compare_with_std();
    std::puts("test_mdarray_mdspan : OK (comparé à std::mdspan)");
#else
    std::puts("test_mdarray_mdspan : OK (std::mdspan absent : pas de comparaison)");
#endif
}
