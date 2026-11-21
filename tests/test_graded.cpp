#include <wsl/graded.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <unordered_set>
#include <vector>

using namespace wsl;

// Un plafond généreux tient le rôle que jouait la variante largeamique.
using large = series<int, 8>;

// =======================================================================
// Construction et rétrécissement
// =======================================================================
static_assert(std::is_aggregate_v<jet<unsigned, 2>>);
static_assert(std::is_aggregate_v<series<unsigned, 2>>);
static_assert(std::is_aggregate_v<large>);

template <class T, auto... V>
concept buildable = requires { T{V...}; };

static_assert(buildable<jet<unsigned, 2>, 0, 3>);
static_assert(!buildable<jet<unsigned, 2>, -1, 3>);
static_assert(!buildable<jet<std::uint8_t, 1>, 300>);
static_assert(buildable<series<unsigned, 2>, 0, 3>);
static_assert(!buildable<series<unsigned, 2>, -1, 3>);
static_assert(buildable<index_jet<2>, 0, 3>);
static_assert(!buildable<index_jet<2>, -1, 3>);

// =======================================================================
// Le concept couvre les deux mondes
// =======================================================================
static_assert(lexicographic<jet<int, 3>>);
static_assert(lexicographic<series<int, 3>>);
static_assert(lexicographic<large>);
static_assert(!lexicographic<int>);
static_assert(!lexicographic<std::array<int, 3>>);

// =======================================================================
// Support : order, degree, sign -- identiques dans les deux mondes
// =======================================================================
static_assert(jet<int, 5>{0, 3, 0, 7, 0}.order() == 1);
static_assert(jet<int, 5>{0, 3, 0, 7, 0}.degree() == 3);
static_assert(series<int, 5>{0, 3, 0, 7, 0}.order() == 1);
static_assert(series<int, 5>{0, 3, 0, 7, 0}.degree() == 3);
static_assert(jet<int, 4>{}.order() == no_support);
static_assert(jet<int, 4>{}.is_zero());
static_assert(series<int, 4>{}.is_zero());
static_assert(jet<int, 0>{}.is_zero());

// taille contre plafond : deux mots pour deux notions
static_assert(jet<int, 4>{1, 2}.size() == 4);
static_assert(series<int, 4>{1, 2}.capacity() == 4);
static_assert(large::capacity() == 8);

// size() n'existe que là où la taille appartient à la valeur
template <class X>
concept has_size = requires { X::size(); };
static_assert(has_size<jet<int, 3>>);
static_assert(!has_size<series<int, 3>>);

// les alias du cas courant, un par arité
static_assert(std::same_as<index_jet<3>, jet<std::size_t, 3>>);
static_assert(std::same_as<index_series<3>, series<std::size_t, 3>>);
static_assert(has_size<index_jet<3>>);
static_assert(!has_size<index_series<3>>);
static_assert(index_jet<2>{1, 2} == index_series<5>{1, 2});

// accès total des deux côtés
static_assert(jet<int, 2>{1, 2}[9999] == 0);
static_assert(series<int, 2>{1, 2}[9999] == 0);

// order et degree bornent exactement le support
constexpr bool support()
{
    constexpr jet<int, 5> value{0, 3, 0, 7, 0};
    for (std::size_t coordinate = 0; coordinate < value.order(); ++coordinate)
        if (value[coordinate] != 0)
            return false;
    for (std::size_t coordinate = value.degree() + 1; coordinate < 20; ++coordinate)
        if (value[coordinate] != 0)
            return false;
    return value[value.order()] != 0 && value[value.degree()] != 0;
}
static_assert(support());

// order est une valuation, et gouverne l'ordre
constexpr bool valuation()
{
    constexpr jet<int, 4> x{0, 5, 1, 0};
    constexpr jet<int, 4> y{0, -5, 0, 2};
    constexpr jet<int, 4> a{1, 2, 9, 9};
    constexpr jet<int, 4> b{1, 3, 0, 0};
    constexpr auto difference = a - b;
    return (x + y).order() == 2 && (-x).order() == x.order() && (x - x).is_zero() &&
           difference.order() == 1 && difference[difference.order()] < 0 && a < b;
}
static_assert(valuation());

// signe
static_assert(sign(jet<int, 4>{}) == 0);
static_assert(sign(jet<int, 4>{0, 0, -3, 9}) == -1);
static_assert(sign(series<int, 4>{0, 0, 3, 0}) == 1);
static_assert(sign(large{0, -2}) == -1);
static_assert(sign(-jet<int, 3>{0, 4, 9}) == -sign(jet<int, 3>{0, 4, 9}));
static_assert(sign(jet<unsigned, 4>{0, 5}) == 1);

// =======================================================================
// Monde fixe : la taille fait partie de l'identité
// =======================================================================
template <class L, class R>
concept comparable = requires(L lhs, R rhs) { lhs < rhs; };
template <class L, class R>
concept equatable = requires(L lhs, R rhs) { lhs == rhs; };
template <class L, class R>
concept addable = requires(L lhs, R rhs) { lhs + rhs; };
template <class L, class R>
concept accumulable = requires(L lhs, R rhs) { lhs += rhs; };

static_assert(comparable<jet<int, 3>, jet<int, 3>>);
static_assert(addable<jet<int, 3>, jet<int, 3>>);
static_assert(addable<jet<int, 3>, jet<unsigned, 3>>);

// deux tailles différentes ne se rencontrent jamais
static_assert(!comparable<jet<int, 2>, jet<int, 3>>);
static_assert(!equatable<jet<int, 2>, jet<int, 3>>);
static_assert(!addable<jet<int, 2>, jet<int, 3>>);
static_assert(!accumulable<jet<int, 3>, jet<int, 2>>);
static_assert(!accumulable<jet<int, 2>, jet<int, 3>>);

// pas de conversion implicite entre tailles, ni vers le monde flexible
static_assert(!std::convertible_to<jet<int, 2>, jet<int, 3>>);
static_assert(!std::convertible_to<jet<int, 3>, jet<int, 2>>);
static_assert(!std::convertible_to<jet<int, 3>, series<int, 3>>);
static_assert(!std::convertible_to<jet<int, 3>, large>);

// arithmétique interne au monde fixe : le résultat y reste
static_assert(std::same_as<decltype(jet<int, 3>{} + jet<int, 3>{}),
                           jet<int, 3>>);
static_assert(std::same_as<decltype(jet<unsigned, 3>{} - jet<unsigned, 3>{}),
                           jet<unsigned, 3>>);
static_assert(std::same_as<decltype(-jet<unsigned, 3>{}), jet<unsigned, 3>>);
static_assert(std::same_as<decltype(+jet<short, 3>{}), jet<int, 3>>);
static_assert(std::same_as<decltype(jet<short, 2>{} + jet<short, 2>{}),
                           jet<int, 2>>);

// =======================================================================
// Monde flexible : l'equivalence
// =======================================================================
static_assert(series<int, 2>{1, 2} == series<int, 4>{1, 2, 0, 0});
static_assert(series<int, 2>{1, 2} < series<int, 3>{1, 2, 1});
static_assert(series<int, 0>{} == series<int, 3>{});
static_assert(large{1, 2, 0, 0} == series<int, 9>{1, 2});
static_assert(large{1, 2, 0, 0}.degree() == 1);

// élargissement implicite interne au monde flexible
static_assert(std::convertible_to<series<int, 2>, series<int, 4>>);
static_assert(std::convertible_to<series<int, 2>, large>);
static_assert(!std::convertible_to<series<int, 4>, series<int, 2>>);

// un préfixe n'est pas inférieur à son prolongement : {1} vaut {1, 0} et
// doit donc être supérieur à {1, -2}
constexpr bool prefixe_signe()
{
    return large{1} > large{1, -2} && large{1} < large{1, 2} &&
           series<int, 1>{1} > series<int, 2>{1, -2};
}
static_assert(prefixe_signe());

// ordre lexicographique non archimédien, compatible avec la loi de groupe
static_assert(jet<int, 2>{2, 0} > jet<int, 2>{1, 999});
constexpr bool translation_invariante()
{
    constexpr jet<int, 3> x{1, -2, 5};
    constexpr jet<int, 3> y{1, 0, -7};
    constexpr jet<int, 3> z{-4, 8, 3};
    return (x < y) == (x + z < y + z);
}
static_assert(translation_invariante());

// =======================================================================
// Mélange : à sens unique, le résultat est toujours flexible
// =======================================================================
static_assert(std::same_as<decltype(jet<int, 3>{} + series<int, 5>{}),
                           series<int, 5>>);
static_assert(std::same_as<decltype(series<int, 5>{} + jet<int, 3>{}),
                           series<int, 5>>);
static_assert(std::same_as<decltype(jet<int, 3>{} + large{}), large>);
static_assert(equatable<jet<int, 2>, series<int, 5>>);
static_assert(comparable<jet<int, 2>, large>);

constexpr bool melange()
{
    constexpr jet<int, 3> a{1, 2, 3};
    constexpr series<int, 5> b{0, 0, 0, 4};
    return (a + b) == series<int, 5>{1, 2, 3, 4} && a == series<int, 9>{1, 2, 3} &&
           (a + b).capacity() == 5;
}
static_assert(melange());

// affectations composées croisées
static_assert(accumulable<large, jet<int, 3>>);
static_assert(accumulable<series<int, 5>, jet<int, 3>>);
static_assert(accumulable<jet<int, 5>, series<int, 3>>);
static_assert(!accumulable<jet<int, 3>, series<int, 5>>);
static_assert(!accumulable<jet<int, 3>, large>);
static_assert(!accumulable<series<int, 3>, large>);

constexpr bool accumulation()
{
    jet<int, 4> fixe{4, 4, 4, 4};
    fixe += series<int, 2>{1, -1};
    fixe -= jet<int, 4>{2, 2, 2, 2};

    large souple{1, 2};
    souple += jet<int, 4>{0, 0, 0, 7};
    souple -= large{1, 2, 0, 7};

    return fixe == jet<int, 4>{3, 1, 2, 2} && souple.is_zero();
}
static_assert(accumulation());

// avancer une position non signée d'un déplacement signe
constexpr bool avance()
{
    jet<unsigned, 2> value{3, 7};
    value += jet<int, 2>{-3, 1};
    return value == jet<unsigned, 2>{0, 8};
}
static_assert(avance());

// =======================================================================
// Fidélité au type sous-jacent
// =======================================================================
static_assert((jet<unsigned, 2>{3, 7} - jet<unsigned, 2>{5, 2})[0] == 3u - 5u);
static_assert((jet<int, 2>{3, 7} - jet<int, 2>{5, 2})[0] == 3 - 5);
static_assert((-jet<unsigned, 2>{2, 0})[0] == -2u);
static_assert((jet<std::uint32_t, 1>{0} - jet<std::uint32_t, 1>{5})[0] == 0u - 5u);

// piège signe/non-signe reproduit tel quel
static_assert(jet<int, 1>{-1} > jet<unsigned, 1>{1u});
static_assert(jet<int, 1>{2} == jet<unsigned, 1>{2u});

// groupe : x + (y - x) == y
constexpr bool groupe()
{
    constexpr jet<unsigned, 2> x{3, 7};
    constexpr jet<unsigned, 2> y{1, 9};
    return x + (y - x) == y && y + (x - y) == x;
}
static_assert(groupe());

// a -= b et a = a - b coïncident
constexpr bool composee_coherente()
{
    jet<unsigned, 2> value{3, 7};
    value -= jet<unsigned, 2>{5, 2};
    return value == jet<unsigned, 2>{3, 7} - jet<unsigned, 2>{5, 2};
}
static_assert(composee_coherente());

// =======================================================================
// Le pont : to_flex et to_lex
// =======================================================================
static_assert(std::same_as<decltype(to_series(jet<int, 3>{})), series<int, 3>>);
static_assert(std::same_as<decltype(to_series(series<int, 3>{})), series<int, 3>>);
static_assert(std::same_as<decltype(to_series(large{})), large>);
static_assert(std::same_as<decltype(to_jet<2>(large{})), jet<int, 2>>);
static_assert(std::same_as<decltype(to_jet<2>(jet<int, 5>{})), jet<int, 2>>);

// to_flex s'applique uniformément aux deux mondes : le code générique
// écrit sur le concept n'a pas à les distinguer
template <lexicographic X>
constexpr auto souplifie(X const &value)
{
    return to_series(value);
}
static_assert(souplifie(jet<int, 3>{1, 2}) == souplifie(series<int, 3>{1, 2}));
static_assert(souplifie(large{1, 2}) == large{1, 2});

constexpr bool pont()
{
    constexpr jet<int, 3> a{1, 2, 0};
    return to_series(a) == series<int, 3>{1, 2} && to_jet<3>(to_series(a)) == a &&
           to_jet<5>(a) == jet<int, 5>{1, 2, 0, 0, 0} &&
           to_jet<2>(series<int, 5>{1, 2}) == jet<int, 2>{1, 2};
}
static_assert(pont());

// élargissement sans perte, aller-retour fidèle
static_assert(to_jet<2>(to_jet<7>(jet<int, 2>{1, 2})) == jet<int, 2>{1, 2});

// troncature : totale, bien définie, sans précondition
static_assert(to_jet<2>(jet<int, 3>{1, 2, 3}) == jet<int, 2>{1, 2});
static_assert(to_jet<1>(large{9, 8, 7}) == jet<int, 1>{9});
static_assert(to_series<2>(series<int, 4>{1, 2, 3, 4}) == series<int, 2>{1, 2});

// la troncature est un homomorphisme de groupes...
constexpr bool homomorphisme()
{
    constexpr jet<int, 4> x{5, 6, 7, 8};
    constexpr jet<int, 4> y{1, 2, 3, 4};
    return to_jet<2>(x + y) == to_jet<2>(x) + to_jet<2>(y) &&
           to_jet<2>(-x) == -to_jet<2>(x);
}
static_assert(homomorphisme());

// ... mais elle ne préserve pas l'ordre
static_assert(jet<int, 3>{1, 2, 3} > jet<int, 3>{1, 2, 0});
static_assert(to_jet<2>(jet<int, 3>{1, 2, 3}) == to_jet<2>(jet<int, 3>{1, 2, 0}));

// value_cast préserve le monde
static_assert(std::same_as<decltype(value_cast<long>(jet<int, 3>{})),
                           jet<long, 3>>);
static_assert(std::same_as<decltype(value_cast<long>(large{})), series<long, 8>>);
static_assert(value_cast<long>(jet<unsigned, 2>{1, 2}) == jet<long, 2>{1, 2});

// =======================================================================
// Concaténation
// =======================================================================
static_assert(std::same_as<decltype(concat(jet<int, 2>{}, jet<int, 3>{})),
                           jet<int, 5>>);
static_assert(std::same_as<decltype(concat(jet<int, 2>{}, series<int, 3>{})),
                           series<int, 5>>);
static_assert(std::same_as<decltype(concat(jet<int, 2>{}, large{})), series<int, 10>>);

static_assert(concat(jet<int, 2>{1, 2}, jet<int, 3>{3, 4, 5}) ==
              jet<int, 5>{1, 2, 3, 4, 5});
static_assert(concat(jet<int, 2>{1, 0}, jet<int, 2>{0, 5}) ==
              jet<int, 4>{1, 0, 0, 5});
static_assert(concat(jet<int, 2>{1, 2}, large{3}) == series<int, 3>{1, 2, 3});
static_assert(concat(jet<int, 0>{}, jet<int, 2>{1, 2}) == jet<int, 2>{1, 2});

// variadique, par pliage à gauche
static_assert(concat(jet<int, 1>{1}, jet<int, 1>{2}, jet<int, 2>{3, 4}) ==
              jet<int, 4>{1, 2, 3, 4});

// l'opérande gauche doit avoir une longueur : un series n'en a pas
template <class L, class R>
concept concatenable = requires(L lhs, R rhs) { concat(lhs, rhs); };
static_assert(concatenable<jet<int, 2>, jet<int, 2>>);
static_assert(concatenable<jet<int, 2>, large>);
static_assert(!concatenable<large, jet<int, 2>>);
static_assert(!concatenable<series<int, 2>, jet<int, 2>>);
// les value_type peuvent différer : le type résultat est hérite de l'addition
static_assert(concatenable<jet<int, 2>, jet<unsigned, 2>>);
static_assert(std::same_as<decltype(concat(jet<int, 2>{}, jet<unsigned, 2>{})),
                           jet<unsigned, 4>>);
static_assert(std::same_as<decltype(concat(jet<short, 2>{}, jet<short, 2>{})),
                           jet<int, 4>>);   // promotion entière

// la réduction littérale : concat(a, b) est la somme des deux, complétés de zéros
constexpr bool reduction()
{
    constexpr jet<int, 2> a{1, 2};
    constexpr jet<int, 3> b{3, 4, 5};
    return concat(a, b) ==
           concat(a, jet<int, 3>{}) + concat(jet<int, 2>{}, b);
}
static_assert(reduction());

// le piège signe/non-signe s'applique ici comme à l'addition, ni plus ni moins
static_assert(concat(jet<unsigned, 2>{1, 2}, jet<int, 2>{-3, 4})[2] ==
              static_cast<unsigned>(-3));
static_assert(concat(jet<int, 2>{1, 2}, jet<int, 2>{-3, 4}) ==
              jet<int, 4>{1, 2, -3, 4});

// c'est le produit lexicographique : l'ordre du couple, dans l'ordre
constexpr bool produit_lexicographique()
{
    constexpr jet<int, 1> a{1};
    constexpr jet<int, 1> b{2};
    constexpr jet<int, 1> x{7};
    constexpr jet<int, 1> y{9};
    return concat(a, y) < concat(b, x) &&   // a < b décide, quels que soient x et y
           concat(a, x) < concat(a, y) &&   // a == a, c'est x < y qui décide
           concat(a, x) == concat(a, x);
}
static_assert(produit_lexicographique());

// l'opérateur virgule est un alias, à employer entre parenthèses
static_assert((jet<int, 2>{1, 2}, jet<int, 2>{3, 4}) ==
              concat(jet<int, 2>{1, 2}, jet<int, 2>{3, 4}));
static_assert((jet<int, 1>{1}, jet<int, 1>{2}, jet<int, 1>{3}) ==
              jet<int, 3>{1, 2, 3});

// la garantie de séquencement gauche-avant-droite, rendue observable
constexpr bool sequencement()
{
    int compteur = 0;
    auto marque = [&compteur](int rang)
    { return jet<int, 1>{(compteur = compteur * 10 + rang, rang)}; };
    auto const resultat = (marque(1), marque(2));
    return compteur == 12 && resultat == jet<int, 2>{1, 2};
}
static_assert(sequencement());

// =======================================================================
// Le marqueur d'arité se propage
// =======================================================================
// exact + exact de même taille reste exact ; des qu'un plafond entre en jeu,
// le résultat est borné -- le mélange est à sens unique
static_assert(std::same_as<decltype(jet<int, 3>{} + jet<int, 3>{}), jet<int, 3>>);
static_assert(std::same_as<decltype(jet<int, 3>{} + series<int, 3>{}), series<int, 3>>);
static_assert(std::same_as<decltype(series<int, 3>{} + jet<int, 5>{}), series<int, 5>>);
static_assert(std::same_as<decltype(-jet<int, 3>{}), jet<int, 3>>);
static_assert(std::same_as<decltype(-series<int, 3>{}), series<int, 3>>);
static_assert(std::same_as<decltype(value_cast<long>(jet<int, 3>{})), jet<long, 3>>);

// la concaténation transmet l'arité de son opérande droit
static_assert(std::same_as<decltype(concat(jet<int, 2>{}, jet<int, 3>{})), jet<int, 5>>);
static_assert(std::same_as<decltype(concat(jet<int, 2>{}, series<int, 3>{})),
                           series<int, 5>>);

// deux tailles exactes différentes ne produisent aucun type résultat
template <class L, class R>
concept combinable = requires { typename wsl::detail::graded::combine<L, R, int>; };
static_assert(combinable<jet<int, 3>, jet<int, 3>>);
static_assert(!combinable<jet<int, 3>, jet<int, 2>>);
static_assert(combinable<jet<int, 3>, series<int, 2>>);
static_assert(combinable<series<int, 3>, series<int, 2>>);

// =======================================================================
// Composantes en lecture seule
// =======================================================================
template <class T>
concept subscript_writable = requires(T value) { value[0] = 1; };
static_assert(!subscript_writable<jet<int, 2>>);
static_assert(!subscript_writable<series<int, 2>>);
static_assert(!subscript_writable<large>);
static_assert(std::assignable_from<jet<int, 2> &, jet<int, 2>>);

int main()
{
    std::vector<jet<int, 2>> values{{2, 0}, {1, 9}, {1, 2}, {0, 5}};
    std::ranges::sort(values);
    assert((values.front() == jet<int, 2>{0, 5}));
    assert((values.back() == jet<int, 2>{2, 0}));

    std::vector<large> souples{{2, 0}, {1, 9}, {1, 2}, {0, 5}};
    std::ranges::sort(souples);
    assert((souples.front() == large{0, 5}));
    assert((souples.back() == large{2}));

    // hachage invariant par plafond et cohérent entre les deux arités
    std::hash<series<int, 2>> const court{};
    std::hash<series<int, 5>> const moyen{};
    std::hash<large> const genereux{};
    std::hash<jet<int, 7>> const exact{};
    assert(court(series<int, 2>{1, 2}) == moyen(series<int, 5>{1, 2, 0, 0, 0}));
    assert(court(series<int, 2>{1, 2}) == genereux(large{1, 2}));
    assert(court(series<int, 2>{1, 2}) == exact(jet<int, 7>{1, 2}));

    std::unordered_set<jet<int, 3>> ensemble{{1, 2, 0}};
    assert(ensemble.contains(jet<int, 3>{1, 2, 0}));

    return 0;
}
