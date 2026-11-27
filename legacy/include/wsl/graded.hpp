#pragma once

#include <array>
#include <compare>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <type_traits>

namespace wsl
{
    // =====================================================================
    // Modèle
    //
    // graded<T, N, A> est un N-uplet de coefficients ordonné
    // lexicographiquement, la coordonnée 0 étant la plus significative.
    //
    // C'est littéralement un jet : a0 + a1.e + a2.e^2 + ... où e est un
    // infinitésimal positif. Chaque coordonnée est infiniment moins
    // significative que la précédente, d'où un ordre total non archimédien
    // -- {0,1} ajouté à lui-même n'atteint jamais {1,0}. Tout le vocabulaire
    // vient de la : order() est la valuation, degree() le degré, sign() le
    // signe du coefficient dominant, et la concaténation s'écrit a + e^N.b.
    //
    // L'ordre lexicographique n'est pas un choix de commodité : c'est le
    // seul des deux ordres naturels sur Z^N qui soit total *et* compatible
    // avec la loi de groupe (x < y => x + z < y + z), l'ordre produit
    // n'étant que partiel.
    //
    // Le marqueur A décide si N appartient à l'identité de la valeur.
    //
    //   arity::exact   -- N est la taille. Deux tailles différentes ne se
    //                     rencontrent jamais : c'est une erreur de
    //                     compilation, pas un résultat faux. Les zéros de
    //                     queue sont de vraies coordonnées.
    //
    //   arity::bounded -- N n'est qu'un plafond de représentation. Les zéros
    //                     de queue sont insignifiants : {1,2} et {1,2,0}
    //                     sont la même valeur, à des plafonds différents.
    //
    // Sur une valeur donnée, les deux arités sont observationnellement
    // identiques -- mêmes coordonnées à tout indice, même order(), même
    // degree(), même sign(). La distinction est nominale, portée par le
    // système de types et par lui seul.
    //
    // Aucune convention de signe n'est introduite par la bibliothèque : un
    // graded se comporte sous ses opérateurs comme T sous les siens. Réserve
    // héritée du langage : unsigned est Z/2^n, qui n'admet aucun ordre
    // compatible avec l'addition.
    // =====================================================================

    enum class arity
    {
        exact,
        bounded
    };

    // Support vide. order() et degree() bornent le support ; l'élément neutre
    // n'en a pas, et aucun indice ne peut le représenter.
    inline constexpr std::size_t no_support = std::numeric_limits<std::size_t>::max();

    // ---------------------------------------------------------------------
    // Traits d'arithmétique.
    //
    // Principe unique : le type du résultat d'une opération est le type de
    // l'expression scalaire correspondante. Aucun arbitrage sur le signe.
    //
    // add_t et subtract_t coïncident pour les entiers, les conversions
    // arithmétiques usuelles ne dépendant pas de l'opérateur. Ils restent
    // distincts parce que le principe est énoncé opérateur par opérateur.
    // ---------------------------------------------------------------------
    template <std::integral T, std::integral U>
    using add_t = decltype(std::declval<T>() + std::declval<U>());

    template <std::integral T, std::integral U>
    using subtract_t = decltype(std::declval<T>() - std::declval<U>());

    template <std::integral T>
    using negate_t = decltype(-std::declval<T>());

    // Le plus unaire n'est pas l'identité : il promeut. +short vaut int.
    template <std::integral T>
    using promote_t = decltype(+std::declval<T>());

    template <std::integral T, std::size_t N, arity A>
    struct graded;

    // ---------------------------------------------------------------------
    // Les deux lectures, et le cas courant.
    // ---------------------------------------------------------------------
    template <std::integral T, std::size_t N>
    using jet = graded<T, N, arity::exact>;

    template <std::integral T, std::size_t N>
    using series = graded<T, N, arity::bounded>;

    // Le cas courant garde un nom court, le général reste disponible. Le
    // schéma <domaine>_<structure> se décline dans les deux arités.
    template <std::size_t N>
    using index_jet = jet<std::size_t, N>;

    template <std::size_t N>
    using index_series = series<std::size_t, N>;

    template <class>
    struct is_graded : std::false_type
    {
    };
    template <std::integral T, std::size_t N, arity A>
    struct is_graded<graded<T, N, A>> : std::true_type
    {
    };

    template <class X>
    concept lexicographic = is_graded<X>::value;

    namespace detail::graded
    {
        template <class>
        struct traits;
        template <std::integral T, std::size_t N, arity A>
        struct traits<wsl::graded<T, N, A>>
        {
            using value_type = T;
            static constexpr std::size_t capacity = N;
            static constexpr wsl::arity arity = A;
            template <std::integral To>
            using rebind = wsl::graded<To, N, A>;
        };

        template <lexicographic X>
        using value_of = typename traits<X>::value_type;

        template <lexicographic X>
        inline constexpr std::size_t capacity_of = traits<X>::capacity;

        template <lexicographic X>
        inline constexpr wsl::arity arity_of = traits<X>::arity;

        constexpr std::size_t larger(std::size_t lhs, std::size_t rhs) noexcept
        {
            return lhs < rhs ? rhs : lhs;
        }

        constexpr std::size_t smaller(std::size_t lhs, std::size_t rhs) noexcept
        {
            return rhs < lhs ? rhs : lhs;
        }

        // Fin du support : premier indice au-delà duquel tout est nul, donc
        // longueur minimale d'un représentant. Toutes les boucles s'y
        // arrêtent.
        template <lexicographic X>
        constexpr std::size_t support_end(X const &operand) noexcept
        {
            std::size_t const last = operand.degree();
            return last == no_support ? 0 : last + 1;
        }

        // Une affectation composée doit ranger son résultat dans l'opérande
        // de gauche sans changer sa taille. Entre deux arités exactes, la
        // taille étant l'identité, seule l'egalite est recevable ; des qu'un
        // plafond est en jeu, un plafond inférieur ou égal suffit.
        template <class Other, std::size_t N, wsl::arity A>
        concept assignable_into =
            lexicographic<Other> &&
            (arity_of<Other> == wsl::arity::exact && A == wsl::arity::exact
                 ? capacity_of<Other> == N
                 : capacity_of<Other> <= N);
    }

    // ---------------------------------------------------------------------
    // graded
    //
    // Agrégat, ce qui est le seul mécanisme appliquant la règle de
    // rétrécissement avec son exception pour les expressions constantes
    // représentables.
    //
    //     jet<unsigned, 2> p{0, 3};    // ok
    //     jet<unsigned, 2> q{-1, 3};   // mal forme
    //     jet<unsigned, 3> r{1, 2};    // ok, vaut {1, 2, 0}
    //
    // La troisième ligne remplit la queue de zéros. En arité exacte ce sont
    // de vraies coordonnées ; en arité bornée, c'est l'inclusion canonique
    // et {1,2} en est un représentant plus court.
    // ---------------------------------------------------------------------
    template <std::integral T, std::size_t N, arity A>
    struct graded
    {
        using value_type = T;
        using storage_type = std::array<T, N>;

        // Le plafond de représentation, toujours disponible.
        static constexpr std::size_t capacity() noexcept { return N; }

        // La taille n'existe que là où elle a un sens : en arité bornée,
        // elle n'appartient pas à la valeur et le mot serait un mensonge.
        static constexpr std::size_t size() noexcept
            requires(A == arity::exact)
        {
            return N;
        }

        storage_type elems;

        // Accès total : défini sur N entier, nul au-delà du plafond.
        // Renvoie par valeur, il n'existe pas de zéro adressable.
        constexpr value_type operator[](std::size_t coordinate) const noexcept
        {
            return coordinate < N ? elems[coordinate] : value_type{};
        }

        // order() : indice de la première coordonnée non nulle -- l'ordre de
        // la valeur au sens des séries formelles, c'est-à-dire sa valuation
        // dans le groupe ordonné. Elle vérifie
        // order(x + y) >= min(order(x), order(y)), elle est invariante par
        // passage à l'opposé, et elle gouverne l'ordre : x < y se décide
        // entièrement au vu du signe de (x - y) à sa coordonnée order.
        // Les zéros de tête sont significatifs, d'où un parcours en O(N).
        constexpr std::size_t order() const noexcept
        {
            for (std::size_t coordinate = 0; coordinate < N; ++coordinate)
                if (elems[coordinate] != value_type{})
                    return coordinate;
            return no_support;
        }

        // degree() : indice de la dernière coordonnée non nulle, par analogie
        // avec le degré d'un polynôme. En arité bornée, degree() + 1 est en
        // outre la longueur minimale d'un représentant.
        constexpr std::size_t degree() const noexcept
        {
            std::size_t coordinate = N;
            while (coordinate != 0)
                if (elems[--coordinate] != value_type{})
                    return coordinate;
            return no_support;
        }

        constexpr bool is_zero() const noexcept { return degree() == no_support; }

        // Élargissement du plafond : fonction de conversion et non
        // constructeur, ce qui préserve le statut d'agrégat. Réservé à
        // l'arité bornée, où il ne change pas la valeur -- d'où l'implicite.
        // En arité exacte il changerait l'identité, et passe donc par
        // to_jet, explicite.
        template <std::size_t Other>
            requires(A == arity::bounded && Other > N)
        constexpr operator graded<T, Other, arity::bounded>() const noexcept
        {
            graded<T, Other, arity::bounded> result{};
            for (std::size_t coordinate = 0; coordinate < N; ++coordinate)
                result.elems[coordinate] = elems[coordinate];
            return result;
        }

        // Affectations composées : valides exactement quand t op= u l'est
        // pour les scalaires sous-jacents.
        template <lexicographic Other>
            requires detail::graded::assignable_into<Other, N, A> &&
                     requires(T &left, detail::graded::value_of<Other> right) { left += right; }
        constexpr graded &operator+=(Other const &rhs) noexcept
        {
            std::size_t const bound = detail::graded::smaller(N, detail::graded::support_end(rhs));
            for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
                elems[coordinate] += rhs[coordinate];
            return *this;
        }

        template <lexicographic Other>
            requires detail::graded::assignable_into<Other, N, A> &&
                     requires(T &left, detail::graded::value_of<Other> right) { left -= right; }
        constexpr graded &operator-=(Other const &rhs) noexcept
        {
            std::size_t const bound = detail::graded::smaller(N, detail::graded::support_end(rhs));
            for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
                elems[coordinate] -= rhs[coordinate];
            return *this;
        }

        // Non-templates, donc préférés aux surcharges génériques quand les
        // deux opérandes ont exactement le même type. Licites parce que deux
        // graded de même plafond ont la même longueur : l'ordre
        // lexicographique sur les tableaux coïncide alors avec l'ordre sur
        // les suites étendues.
        constexpr bool operator==(graded const &) const = default;
        constexpr std::strong_ordering operator<=>(graded const &) const = default;
    };

    namespace detail::graded
    {
        // Type résultat d'une opération binaire.
        //
        // Deux arités exactes de tailles différentes ne satisfont pas la
        // contrainte : c'est cette absence qui produit l'erreur de
        // compilation faisant tout l'intérêt de l'arité exacte. Des qu'un
        // plafond est en jeu, le résultat est borné -- le mélange est à sens
        // unique.
        template <class A, class B, class V>
        struct combine_impl
        {
        };

        template <std::integral T, std::size_t N, wsl::arity AA, std::integral U,
                  std::size_t M, wsl::arity BB, class V>
            requires(AA == wsl::arity::bounded || BB == wsl::arity::bounded || N == M)
        struct combine_impl<wsl::graded<T, N, AA>, wsl::graded<U, M, BB>, V>
        {
            using type = wsl::graded<V, larger(N, M),
                                (AA == wsl::arity::exact && BB == wsl::arity::exact)
                                    ? wsl::arity::exact
                                    : wsl::arity::bounded>;
        };

        template <class A, class B, class V>
        using combine = typename combine_impl<A, B, V>::type;

        template <class A, class B>
        concept mixable = requires { typename combine_impl<A, B, int>::type; };

        // Une seule représentation, donc une seule façon de construire.
        template <class Result, class Generator>
        constexpr Result build(Generator generator) noexcept
        {
            Result result{};
            for (std::size_t coordinate = 0; coordinate < capacity_of<Result>; ++coordinate)
                result.elems[coordinate] = generator(coordinate);
            return result;
        }

        template <class Result, lexicographic A, lexicographic B, class Operation>
        constexpr Result zip(A const &lhs, B const &rhs, Operation operation) noexcept
        {
            return build<Result>([&](std::size_t coordinate)
                                 { return operation(lhs[coordinate], rhs[coordinate]); });
        }
    }

    // ---------------------------------------------------------------------
    // Comparaisons.
    //
    // Deux arités exactes de tailles différentes ne sont pas comparables :
    // c'est une erreur de compilation, pas un résultat faux. Tout le reste
    // l'est, les opérandes étant lus comme des suites étendues.
    //
    // Les conversions arithmétiques usuelles s'appliquent coordonnée par
    // coordonnée : le piège signe/non-signe est reproduit tel quel.
    //
    // Pour l'ordre, ce piège demande un detour. Le langage refuse int <=>
    // unsigned (conversion rétrécissante), mais accepte int < unsigned. On
    // procède donc comme std::pair, std::tuple et les conteneurs, par la
    // comparaison synthétisée de la norme (synth-three-way,
    // [expos.only.entity]) : <=> quand il existe, sinon < dans les deux sens.
    // Le résultat est alors un std::weak_ordering, comme pour std::pair.
    // ---------------------------------------------------------------------
    namespace detail::graded
    {
        template <class T, class U>
        constexpr auto synth_three_way(T const &lhs, U const &rhs) noexcept
        {
            if constexpr (std::three_way_comparable_with<T, U>)
                return lhs <=> rhs;
            else
            {
                if (lhs < rhs)
                    return std::weak_ordering::less;
                if (rhs < lhs)
                    return std::weak_ordering::greater;
                return std::weak_ordering::equivalent;
            }
        }

        template <class T, class U>
        using synth_three_way_result =
            decltype(synth_three_way(std::declval<T const &>(), std::declval<U const &>()));
    }

    template <lexicographic A, lexicographic B>
        requires detail::graded::mixable<A, B>
    constexpr bool operator==(A const &lhs, B const &rhs) noexcept
    {
        constexpr std::size_t bound =
            detail::graded::larger(detail::graded::capacity_of<A>, detail::graded::capacity_of<B>);
        for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
            if (!(lhs[coordinate] == rhs[coordinate]))
                return false;
        return true;
    }

    template <lexicographic A, lexicographic B>
        requires detail::graded::mixable<A, B>
    constexpr auto operator<=>(A const &lhs, B const &rhs) noexcept
        -> detail::graded::synth_three_way_result<detail::graded::value_of<A>, detail::graded::value_of<B>>
    {
        constexpr std::size_t bound =
            detail::graded::larger(detail::graded::capacity_of<A>, detail::graded::capacity_of<B>);
        for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
            if (auto const result = detail::graded::synth_three_way(lhs[coordinate], rhs[coordinate]);
                result != 0)
                return result;
        return std::strong_ordering::equal;
    }

    // ---------------------------------------------------------------------
    // Arithmétique.
    // ---------------------------------------------------------------------
    template <lexicographic A, lexicographic B>
        requires detail::graded::mixable<A, B>
    constexpr auto operator+(A const &lhs, B const &rhs) noexcept
        -> detail::graded::combine<A, B, add_t<detail::graded::value_of<A>, detail::graded::value_of<B>>>
    {
        using result_type =
            detail::graded::combine<A, B, add_t<detail::graded::value_of<A>, detail::graded::value_of<B>>>;
        return detail::graded::zip<result_type>(
            lhs, rhs, [](auto left, auto right) noexcept { return left + right; });
    }

    // Rien de particulier : left - right, comme pour deux scalaires. La
    // différence de deux valeurs non signées reste non signée et modulaire,
    // exactement comme 3u - 5u. Conséquence utile, a -= b et a = a - b
    // désignent la même chose.
    template <lexicographic A, lexicographic B>
        requires detail::graded::mixable<A, B>
    constexpr auto operator-(A const &lhs, B const &rhs) noexcept
        -> detail::graded::combine<A, B, subtract_t<detail::graded::value_of<A>, detail::graded::value_of<B>>>
    {
        using result_type =
            detail::graded::combine<A, B, subtract_t<detail::graded::value_of<A>, detail::graded::value_of<B>>>;
        return detail::graded::zip<result_type>(
            lhs, rhs, [](auto left, auto right) noexcept { return left - right; });
    }

    template <lexicographic A>
    constexpr auto operator-(A const &operand) noexcept ->
        typename detail::graded::traits<A>::template rebind<negate_t<detail::graded::value_of<A>>>
    {
        using result_type =
            typename detail::graded::traits<A>::template rebind<negate_t<detail::graded::value_of<A>>>;
        return detail::graded::build<result_type>(
            [&](std::size_t coordinate) { return -operand[coordinate]; });
    }

    template <lexicographic A>
    constexpr auto operator+(A const &operand) noexcept ->
        typename detail::graded::traits<A>::template rebind<promote_t<detail::graded::value_of<A>>>
    {
        using result_type =
            typename detail::graded::traits<A>::template rebind<promote_t<detail::graded::value_of<A>>>;
        return detail::graded::build<result_type>(
            [&](std::size_t coordinate) { return +operand[coordinate]; });
    }

    // ---------------------------------------------------------------------
    // Signe.
    //
    // Le signe d'une valeur non nulle est celui de sa coordonnée de rang
    // order() : les précédentes sont nulles par définition, les suivantes
    // sont infiniment moins significatives. Le neutre est de signe nul.
    // Équivalent à comparer la valeur au neutre.
    //
    // Attention pour T non signe : sign(x - y) ne décide pas la comparaison
    // de x et y, la soustraction y étant modulaire. Conséquence directe du
    // principe de fidélité au type sous-jacent.
    // ---------------------------------------------------------------------
    template <lexicographic A>
    constexpr int sign(A const &operand) noexcept
    {
        std::size_t const leading = operand.order();
        if (leading == no_support)
            return 0;
        if constexpr (std::is_signed_v<detail::graded::value_of<A>>)
            return operand[leading] < detail::graded::value_of<A>{} ? -1 : 1;
        else
            return 1;
    }

    // ---------------------------------------------------------------------
    // Concaténation.
    //
    // C'est le produit lexicographique : concat(a,b) < concat(a',b') si et
    // seulement si a < a', ou a == a' et b < b'. Autrement dit un
    // isomorphisme de groupes qui préserve l'ordre. En termes de jets,
    // concat(a, b) vaut a + e^N.b -- ce qui explique pourquoi l'opérande
    // gauche doit être d'arité exacte : sans longueur, pas de puissance de e.
    //
    // En arité bornée, l'opération ne passerait d'ailleurs pas au quotient :
    // si {1,0} et {1} sont la même valeur, les concaténer avec {5} devrait
    // rendre {1,0,5} dans un cas et {1,5} dans l'autre. L'opérande droit,
    // lui, n'a besoin que d'être décalé, et le décalage respecte
    // l'équivalence -- il peut donc être de l'une ou l'autre arité, qu'il
    // transmet au résultat.
    //
    // Le type résultat n'est pas inventé, il est hérité : concat(a, b) vaut
    // a complété de zéros à droite, plus b décalé et complété de zéros à
    // gauche. Chaque coordonnée subit donc exactement la conversion que
    // l'addition correspondante lui aurait fait subir.
    // ---------------------------------------------------------------------
    namespace detail::graded
    {
        // Réalise littéralement la somme des deux opérandes complétés de
        // zéros : chaque coordonnée est additionnée au zéro de l'autre type,
        // ce qui applique les conversions arithmétiques usuelles sans que la
        // bibliothèque ait à les reproduire.
        template <class Result, std::size_t N, lexicographic A, lexicographic B>
        constexpr Result concatenate(A const &lhs, B const &rhs) noexcept
        {
            return build<Result>(
                [&](std::size_t coordinate)
                {
                    return coordinate < N ? lhs[coordinate] + value_of<B>{}
                                          : value_of<A>{} + rhs[coordinate - N];
                });
        }
    }

    template <std::integral T, std::size_t N, std::integral U, std::size_t M, arity B>
    constexpr graded<add_t<T, U>, N + M, B> concat(jet<T, N> const &lhs,
                                                   graded<U, M, B> const &rhs) noexcept
    {
        return detail::graded::concatenate<graded<add_t<T, U>, N + M, B>, N>(lhs, rhs);
    }

    template <lexicographic A, lexicographic B, lexicographic... Rest>
        requires(sizeof...(Rest) != 0)
    constexpr auto concat(A const &first, B const &second, Rest const &...rest) noexcept
    {
        return concat(concat(first, second), rest...);
    }

    // Alias de concat. Depuis C++17 les opérandes d'un opérateur surcharge
    // sont séquencés dans l'ordre prescrit pour l'opérateur intégré, donc la
    // virgule surchargée conserve la garantie gauche-avant-droite.
    // A employer entre parenthèses : la virgule ayant la précédence la plus
    // basse du langage, x = a, b s'analyse (x = a), b.
    template <lexicographic A, lexicographic B>
        requires requires(A const &lhs, B const &rhs) { concat(lhs, rhs); }
    constexpr auto operator,(A const &lhs, B const &rhs) noexcept
    {
        return concat(lhs, rhs);
    }

    // ---------------------------------------------------------------------
    // Conversions, toutes explicites et toutes totales.
    //
    // Chacune nomme sa destination ou ce qu'elle change : value_cast le type
    // des coefficients, to_jet et to_series l'arité et le plafond.
    //
    // Aucune ne pose de précondition, aucune n'est un comportement indéfini.
    // Une troncature explicite est bien définie, sur le modèle de
    // static_cast entre entiers : l'appel écrit à la main est l'annonce.
    //
    // La troncature est un homomorphisme de groupes -- tronquer une somme
    // revient à sommer les troncatures -- mais elle ne préserve pas l'ordre,
    // {1,2,3} et {1,2,0} s'ecrasant sur la même valeur.
    // ---------------------------------------------------------------------
    template <std::integral To, lexicographic A>
    constexpr auto value_cast(A const &operand) noexcept ->
        typename detail::graded::traits<A>::template rebind<To>
    {
        using result_type = typename detail::graded::traits<A>::template rebind<To>;
        return detail::graded::build<result_type>([&](std::size_t coordinate)
                                          { return static_cast<To>(operand[coordinate]); });
    }

    // Entrée dans l'arité exacte, depuis n'importe quelle source. Élargit
    // sans perte quand la taille de destination suffit, tronque sinon.
    template <std::size_t N, lexicographic A>
    constexpr jet<detail::graded::value_of<A>, N> to_jet(A const &operand) noexcept
    {
        return detail::graded::build<jet<detail::graded::value_of<A>, N>>(
            [&](std::size_t coordinate) { return operand[coordinate]; });
    }

    // Entrée dans l'arité bornée. La forme sans plafond conserve celui de la
    // source ; elle rend l'injection canonique, qui ne perd jamais rien.
    template <std::size_t N, lexicographic A>
    constexpr series<detail::graded::value_of<A>, N> to_series(A const &operand) noexcept
    {
        return detail::graded::build<series<detail::graded::value_of<A>, N>>(
            [&](std::size_t coordinate) { return operand[coordinate]; });
    }

    template <lexicographic A>
    constexpr series<detail::graded::value_of<A>, detail::graded::capacity_of<A>>
    to_series(A const &operand) noexcept
    {
        return to_series<detail::graded::capacity_of<A>>(operand);
    }

    namespace detail::graded
    {
        // Le hachage porte sur le support, sans quoi des valeurs
        // équivalentes de plafonds différents auraient des empreintes
        // différentes. La cohérence est garantie à value_type fixe ; entre
        // types de signes différents, deux valeurs égales au sens des
        // conversions usuelles peuvent se hacher différemment.
        template <lexicographic A>
        constexpr std::size_t hash_of(A const &operand) noexcept
        {
            std::size_t result = 14695981039346656037ULL;
            std::size_t const bound = support_end(operand);
            for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
            {
                result ^= static_cast<std::size_t>(
                    static_cast<std::make_unsigned_t<value_of<A>>>(operand[coordinate]));
                result *= 1099511628211ULL;
            }
            return result;
        }
    }
}

template <std::integral T, std::size_t N, wsl::arity A>
struct std::hash<wsl::graded<T, N, A>>
{
    constexpr std::size_t operator()(wsl::graded<T, N, A> const &operand) const noexcept
    {
        return wsl::detail::graded::hash_of(operand);
    }
};

