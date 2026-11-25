// wsl::mdarray, positions et déplacements (docs/mdarray.md, §4) : construction,
// modèle affine, refus à la compilation (P1, P2), préconditions en évaluation
// constante (P3), ordres, assemblage, normes, division euclidienne, hachage,
// formatage.

#include <wsl/mdarray.hpp>

#include <array>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <format>
#include <span>
#include <sstream>
#include <type_traits>
#include <unordered_set>

namespace
{
    using wsl::mddisplacement;
    using wsl::mdposition;

    using P2 = mdposition<2>;
    using D2 = mddisplacement<2>;

    // Une expression est-elle une expression constante ? Une précondition violée
    // en évaluation constante doit la rendre invalide (P3).
    template <auto F>
    concept constant = requires { typename std::bool_constant<(F(), true)>; };

    // Opérations possibles (dans un contexte de modèle, une expression invalide rend
    // le concept faux au lieu d'être une erreur).
    template <class A, class B> concept addable = requires(A a, B b) { a + b; };
    template <class A, class B> concept subtractable = requires(A a, B b) { a - b; };
    template <class A, class B> concept multipliable = requires(A a, B b) { a * b; };
    template <class A, class B> concept divisible = requires(A a, B b) { a / b; };
    template <class A, class B> concept modulable = requires(A a, B b) { a % b; };
    template <class A, class B> concept comparable_eq = requires(A a, B b) { a == b; };
    template <class A, class B> concept comparable_less = requires(A a, B b) { a < b; };
    template <class A> concept negatable = requires(A a) { -a; };
    template <class A> concept incrementable = requires(A a) { ++a; };
    template <class A, class B> concept partially_ordered = requires(A a, B b) { all_less(a, b); };

    // ---------------------------------------------------------------------
    // Refus à la compilation (P1, P2, D2, D5, D11, D12, D13, D16, D21)
    // ---------------------------------------------------------------------
    static_assert(!addable<P2, P2>, "position + position");
    static_assert(!negatable<P2>, "−position");
    static_assert(!multipliable<int, P2> && !multipliable<P2, int>, "entier × position");
    static_assert(!divisible<P2, int> && !modulable<P2, int>, "position / entier, position % entier");
    static_assert(!incrementable<P2> && !incrementable<D2>, "++");
    static_assert(!addable<P2, int> && !addable<D2, int>, "position + entier, déplacement + entier");
    static_assert(!multipliable<D2, D2>, "déplacement × déplacement");
    static_assert(!comparable_eq<P2, D2> && !comparable_less<P2, D2>, "position et déplacement comparés");
    static_assert(!subtractable<P2, mdposition<3>> && !comparable_eq<P2, mdposition<3>>, "rangs différents");
    static_assert(!addable<P2, mddisplacement<3>>, "rangs différents");
    static_assert(addable<P2, D2> && subtractable<P2, P2> && multipliable<int, D2> && negatable<D2>, "opérations du modèle");
    static_assert(!std::is_constructible_v<P2, int>, "une seule composante pour un rang 2");
    static_assert(!std::is_constructible_v<P2, int, int, int>, "trois composantes pour un rang 2");
    static_assert(!std::is_constructible_v<P2, std::array<int, 3>>, "tableau de rang 3");
    static_assert(!std::is_convertible_v<int, mdposition<1>>, "pas de conversion implicite d'un entier");
    static_assert(!std::is_convertible_v<int, mddisplacement<1>>, "pas de conversion implicite d'un entier");
    static_assert(!std::is_constructible_v<P2, D2>, "pas de conversion déplacement → position");
    static_assert(!std::is_constructible_v<D2, P2>, "pas de conversion position → déplacement");
    static_assert(!std::is_constructible_v<P2, bool, bool>, "bool n'est pas un indice");
    static_assert(!partially_ordered<P2, D2> && partially_ordered<P2, P2>, "ordre partiel entre sortes différentes");

    // Conversions entre types de composantes (D28) : implicites sans perte.
    static_assert(std::is_convertible_v<mdposition<2, int>, mdposition<2, std::size_t>>);
    static_assert(!std::is_convertible_v<mdposition<2, std::size_t>, mdposition<2, int>>);
    static_assert(std::is_constructible_v<mdposition<2, int>, mdposition<2, std::size_t>>);
    static_assert(std::is_convertible_v<mddisplacement<2, int>, mddisplacement<2, long long>>);
    static_assert(!std::is_convertible_v<mddisplacement<2, long long>, mddisplacement<2, int>>);

    // Déduction (D4) et types des résultats (D28).
    static_assert(std::is_same_v<decltype(mdposition{2, 1}), mdposition<2, std::size_t>>);
    static_assert(std::is_same_v<decltype(mddisplacement{0, -1}), mddisplacement<2, std::ptrdiff_t>>);
    static_assert(std::is_same_v<decltype(P2{} - P2{}), mddisplacement<2, std::ptrdiff_t>>);
    static_assert(std::is_same_v<decltype(mdposition<2, int>{} - mdposition<2, int>{}), mddisplacement<2, int>>);
    static_assert(std::is_same_v<decltype(mdposition<2, int>{} + D2{}), mdposition<2, int>>);
    static_assert(std::is_same_v<decltype(mddisplacement<2, int>{} + D2{}), D2>);
    static_assert(std::is_same_v<decltype(concat(P2{}, mdposition<1>{})), mdposition<3>>);
    static_assert(P2::rank() == 2 && mdposition<0>::rank() == 0);

    // ---------------------------------------------------------------------
    // Préconditions en évaluation constante (P3)
    // ---------------------------------------------------------------------
    static_assert(constant<[] { return P2{2, 1}; }>);
    static_assert(!constant<[] { return P2{-1, 0}; }>, "composante négative");
    static_assert(!constant<[] { return mdposition<1, std::uint8_t>{300}; }>, "composante non représentable");
    static_assert(!constant<[] { return P2{0, 1} + D2{-1, 0}; }>, "position + déplacement négative");
    static_assert(!constant<[] { return P2{0, 1} - D2{1, 0}; }>, "position − déplacement négative");
    static_assert(!constant<[] { return mdposition<1, std::uint8_t>{250} + mddisplacement<1>{10}; }>, "dépassement");
    static_assert(!constant<[] { return mddisplacement<1, std::int8_t>{-128} - mddisplacement<1, std::int8_t>{1}; }>);
    static_assert(!constant<[] { return -mddisplacement<1, std::int8_t>{-128}; }>, "−min");
    static_assert(!constant<[] { return D2{1, 2} / 0; }>, "division par zéro");
    static_assert(!constant<[] { return D2::unit(2); }>, "direction hors du rang");
    static_assert(!constant<[] { return P2{1, 2}[2]; }>, "composante hors du rang");
    static_assert(!constant<[] { return wsl::rem_euclid(3, 0); }>, "diviseur nul");
    static_assert(!constant<[] { return mdposition<1, std::uint64_t>{std::uint64_t{1} << 63} - mdposition<1, std::uint64_t>{0}; }>,
                  "différence non représentable dans le type signé");
    static_assert(constant<[] { return mdposition<1, std::uint64_t>{0} - mdposition<1, std::uint64_t>{std::uint64_t{1} << 63}; }>,
                  "−2⁶³ est représentable");

    // ---------------------------------------------------------------------
    // Valeurs, à la compilation
    // ---------------------------------------------------------------------
    constexpr bool arithmetic()
    {
        P2 p{2, 1};
        D2 d{0, -1};
        P2 q = p + d;
        if (q != P2{2, 0} || q - p != d || d + p != q || q - d != p)
            return false;
        p += D2{1, 1};
        p -= D2{0, 2};
        if (p != P2{3, 0})
            return false;
        D2 e = D2{3, -4};
        if (e + d != D2{3, -5} || e - d != D2{3, -3} || -e != D2{-3, 4})
            return false;
        if (2 * e != D2{6, -8} || e * 2 != D2{6, -8})
            return false;
        // Troncature vers zéro, comme les entiers de C++ (D10).
        if (D2{-3, 5} / 2 != D2{-1, 2} || D2{-3, 5} % 2 != D2{-1, 1})
            return false;
        // Division euclidienne (D15).
        if (div_euclid(D2{-3, 5}, 2) != D2{-2, 2} || rem_euclid(D2{-3, 5}, 2) != D2{1, 1})
            return false;
        if (wsl::rem_euclid(-1, 4) != 3 || wsl::div_euclid(-1, 4) != -1 || wsl::rem_euclid(7u, 4u) != 3u)
            return false;
        // Normes (D26).
        if (norm1(e) != 7 || norm_inf(e) != 4 || norm2_squared(e) != 25)
            return false;
        // filled, unit (D12, D25).
        if (P2{1, 1} + D2::filled(1) != P2{2, 2} || P2::filled(3) != P2{3, 3})
            return false;
        if (P2{1, 1} + D2::unit(1) != P2{1, 2})
            return false;
        // Origine, déplacement nul (D2, D6).
        if (P2{} != P2{0, 0} || D2{} != D2{0, 0} || P2{4, 2} - P2{} != D2{4, 2})
            return false;
        // Types de composantes mêlés (D28) : valeurs mathématiques.
        if (mdposition<2, int>{2, 1} != P2{2, 1} || mddisplacement<2, int>{-1, 0} != D2{-1, 0})
            return false;
        if (mdposition<2, int>{0, 3} - P2{2, 1} != D2{-2, 2})
            return false;
        // Assemblage et extraction (D24).
        auto c = concat(P2{1, 2}, mdposition<1>{3});
        if (c != mdposition<3>{1, 2, 3} || c.first<2>() != P2{1, 2} || c.last<1>() != mdposition<1>{3})
            return false;
        if (concat(D2{-1, 0}, mddisplacement<1>{2}) != mddisplacement<3>{-1, 0, 2})
            return false;
        // Accès, tableau, décomposition (D4).
        auto [i, j] = P2{5, 7};
        if (i != 5 || j != 7 || P2{5, 7}[1] != 7 || P2{5, 7}.to_array() != std::array<std::size_t, 2>{5, 7})
            return false;
        std::array<int, 2> a{4, 6};
        P2 from_array = a;
        std::span<int const, 2> s(a);
        P2 from_span = s;
        P2 from_dynamic(std::span<int const>(a.data(), 2));
        if (from_array != P2{4, 6} || from_span != P2{4, 6} || from_dynamic != P2{4, 6})
            return false;
        return true;
    }
    static_assert(arithmetic());

    constexpr bool orders()
    {
        // Ordre lexicographique, celui du parcours en ligne majeure (D17).
        if (!(P2{0, 4} < P2{1, 0}) || !(P2{1, 0} > P2{0, 9}) || (P2{1, 2} <=> P2{1, 2}) != std::strong_ordering::equal)
            return false;
        if (!(D2{-1, 5} < D2{0, -5}))
            return false;
        // Pas de piège signé / non signé : valeurs mathématiques.
        if (!(mddisplacement<1, int>{-1} < mddisplacement<1, long long>{0}))
            return false;
        // Ordre composante par composante (D18).
        if (!all_less(P2{0, 1}, P2{1, 2}) || all_less(P2{0, 2}, P2{1, 2}) || !all_less_equal(P2{0, 2}, P2{1, 2}))
            return false;
        if (all_less_equal(P2{0, 7}, P2{4, 5}))
            return false;
        return true;
    }
    static_assert(orders());

    void runtime()
    {
        assert(arithmetic());
        assert(orders());

        // Hachage (D27).
        std::unordered_set<P2> seen{P2{1, 2}, P2{2, 1}, P2{1, 2}};
        assert(seen.size() == 2 && seen.contains(P2{2, 1}));
        std::unordered_set<D2> moves{D2{0, 1}, D2{0, -1}, D2{1, 0}, D2{-1, 0}};
        assert(moves.size() == 4);

        // Formatage (D27).
        assert(std::format("{}", P2{2, 1}) == "(2, 1)");
        assert(std::format("{}", D2{0, -1}) == "(0, -1)");
        assert(std::format("{:3}", P2{2, 1}) == "(  2,   1)");
        assert(std::format("{:x}", mdposition<3>{10, 11, 255}) == "(a, b, ff)");
        assert(std::format("{}", mdposition<0>{}) == "()");
        std::ostringstream out;
        out << P2{3, 4} << ' ' << D2{-1, 1};
        assert(out.str() == "(3, 4) (-1, 1)");

        // Valeurs à l'exécution (pas d'évaluation constante).
        volatile int k = 3;
        P2 p{k, k + 1};
        assert((p - P2{1, 1} == D2{2, 3}));
        assert(norm_inf(D2{-k, 2}) == 3);
    }
}

int main()
{
    runtime();
    std::puts("test_mdarray_position : OK");
}
