#ifndef WSL_GRADED_HPP
#define WSL_GRADED_HPP

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
    // Modele
    //
    // graded<T, N, A> est un N-uplet de coefficients ordonne
    // lexicographiquement, la coordonnee 0 etant la plus significative.
    //
    // C'est litteralement un jet : a0 + a1.e + a2.e^2 + ... ou e est un
    // infinitesimal positif. Chaque coordonnee est infiniment moins
    // significative que la precedente, d'ou un ordre total non archimedien
    // -- {0,1} ajoute a lui-meme n'atteint jamais {1,0}. Tout le vocabulaire
    // vient de la : order() est la valuation, degree() le degre, sign() le
    // signe du coefficient dominant, et la concatenation s'ecrit a + e^N.b.
    //
    // L'ordre lexicographique n'est pas un choix de commodite : c'est le
    // seul des deux ordres naturels sur Z^N qui soit total *et* compatible
    // avec la loi de groupe (x < y => x + z < y + z), l'ordre produit
    // n'etant que partiel.
    //
    // Le marqueur A decide si N appartient a l'identite de la valeur.
    //
    //   arity::exact   -- N est la taille. Deux tailles differentes ne se
    //                     rencontrent jamais : c'est une erreur de
    //                     compilation, pas un resultat faux. Les zeros de
    //                     queue sont de vraies coordonnees.
    //
    //   arity::bounded -- N n'est qu'un plafond de representation. Les zeros
    //                     de queue sont insignifiants : {1,2} et {1,2,0}
    //                     sont la meme valeur, a des plafonds differents.
    //
    // Sur une valeur donnee, les deux arites sont observationnellement
    // identiques -- memes coordonnees a tout indice, meme order(), meme
    // degree(), meme sign(). La distinction est nominale, portee par le
    // systeme de types et par lui seul.
    //
    // Aucune convention de signe n'est introduite par la bibliotheque : un
    // graded se comporte sous ses operateurs comme T sous les siens. Reserve
    // heritee du langage : unsigned est Z/2^n, qui n'admet aucun ordre
    // compatible avec l'addition.
    // =====================================================================

    enum class arity
    {
        exact,
        bounded
    };

    // Support vide. order() et degree() bornent le support ; l'element neutre
    // n'en a pas, et aucun indice ne peut le representer.
    inline constexpr std::size_t no_support = std::numeric_limits<std::size_t>::max();

    // ---------------------------------------------------------------------
    // Traits d'arithmetique.
    //
    // Principe unique : le type du resultat d'une operation est le type de
    // l'expression scalaire correspondante. Aucun arbitrage sur le signe.
    //
    // add_t et subtract_t coincident pour les entiers, les conversions
    // arithmetiques usuelles ne dependant pas de l'operateur. Ils restent
    // distincts parce que le principe est enonce operateur par operateur.
    // ---------------------------------------------------------------------
    template <std::integral T, std::integral U>
    using add_t = decltype(std::declval<T>() + std::declval<U>());

    template <std::integral T, std::integral U>
    using subtract_t = decltype(std::declval<T>() - std::declval<U>());

    template <std::integral T>
    using negate_t = decltype(-std::declval<T>());

    // Le plus unaire n'est pas l'identite : il promeut. +short vaut int.
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

    // Le cas courant garde un nom court, le general reste disponible. Le
    // schema <domaine>_<structure> se decline dans les deux arites.
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

    namespace detail
    {
        template <class>
        struct traits;
        template <std::integral T, std::size_t N, arity A>
        struct traits<graded<T, N, A>>
        {
            using value_type = T;
            static constexpr std::size_t capacity = N;
            static constexpr wsl::arity arity = A;
            template <std::integral To>
            using rebind = graded<To, N, A>;
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

        // Fin du support : premier indice au-dela duquel tout est nul, donc
        // longueur minimale d'un representant. Toutes les boucles s'y
        // arretent.
        template <lexicographic X>
        constexpr std::size_t support_end(X const &operand) noexcept
        {
            std::size_t const last = operand.degree();
            return last == no_support ? 0 : last + 1;
        }

        // Une affectation composee doit ranger son resultat dans l'operande
        // de gauche sans changer sa taille. Entre deux arites exactes, la
        // taille etant l'identite, seule l'egalite est recevable ; des qu'un
        // plafond est en jeu, un plafond inferieur ou egal suffit.
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
    // Agregat, ce qui est le seul mecanisme appliquant la regle de
    // retrecissement avec son exception pour les expressions constantes
    // representables.
    //
    //     jet<unsigned, 2> p{0, 3};    // ok
    //     jet<unsigned, 2> q{-1, 3};   // mal forme
    //     jet<unsigned, 3> r{1, 2};    // ok, vaut {1, 2, 0}
    //
    // La troisieme ligne remplit la queue de zeros. En arite exacte ce sont
    // de vraies coordonnees ; en arite bornee, c'est l'inclusion canonique
    // et {1,2} en est un representant plus court.
    // ---------------------------------------------------------------------
    template <std::integral T, std::size_t N, arity A>
    struct graded
    {
        using value_type = T;
        using storage_type = std::array<T, N>;

        // Le plafond de representation, toujours disponible.
        static constexpr std::size_t capacity() noexcept { return N; }

        // La taille n'existe que la ou elle a un sens : en arite bornee,
        // elle n'appartient pas a la valeur et le mot serait un mensonge.
        static constexpr std::size_t size() noexcept
            requires(A == arity::exact)
        {
            return N;
        }

        storage_type elems;

        // Acces total : defini sur N entier, nul au-dela du plafond.
        // Renvoie par valeur, il n'existe pas de zero adressable.
        constexpr value_type operator[](std::size_t coordinate) const noexcept
        {
            return coordinate < N ? elems[coordinate] : value_type{};
        }

        // order() : indice de la premiere coordonnee non nulle -- l'ordre de
        // la valeur au sens des series formelles, c'est-a-dire sa valuation
        // dans le groupe ordonne. Elle verifie
        // order(x + y) >= min(order(x), order(y)), elle est invariante par
        // passage a l'oppose, et elle gouverne l'ordre : x < y se decide
        // entierement au vu du signe de (x - y) a sa coordonnee order.
        // Les zeros de tete sont significatifs, d'ou un parcours en O(N).
        constexpr std::size_t order() const noexcept
        {
            for (std::size_t coordinate = 0; coordinate < N; ++coordinate)
                if (elems[coordinate] != value_type{})
                    return coordinate;
            return no_support;
        }

        // degree() : indice de la derniere coordonnee non nulle, par analogie
        // avec le degre d'un polynome. En arite bornee, degree() + 1 est en
        // outre la longueur minimale d'un representant.
        constexpr std::size_t degree() const noexcept
        {
            std::size_t coordinate = N;
            while (coordinate != 0)
                if (elems[--coordinate] != value_type{})
                    return coordinate;
            return no_support;
        }

        constexpr bool is_zero() const noexcept { return degree() == no_support; }

        // Elargissement du plafond : fonction de conversion et non
        // constructeur, ce qui preserve le statut d'agregat. Reserve a
        // l'arite bornee, ou il ne change pas la valeur -- d'ou l'implicite.
        // En arite exacte il changerait l'identite, et passe donc par
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

        // Affectations composees : valides exactement quand t op= u l'est
        // pour les scalaires sous-jacents.
        template <lexicographic Other>
            requires detail::assignable_into<Other, N, A> &&
                     requires(T &left, detail::value_of<Other> right) { left += right; }
        constexpr graded &operator+=(Other const &rhs) noexcept
        {
            std::size_t const bound = detail::smaller(N, detail::support_end(rhs));
            for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
                elems[coordinate] += rhs[coordinate];
            return *this;
        }

        template <lexicographic Other>
            requires detail::assignable_into<Other, N, A> &&
                     requires(T &left, detail::value_of<Other> right) { left -= right; }
        constexpr graded &operator-=(Other const &rhs) noexcept
        {
            std::size_t const bound = detail::smaller(N, detail::support_end(rhs));
            for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
                elems[coordinate] -= rhs[coordinate];
            return *this;
        }

        // Non-templates, donc preferes aux surcharges generiques quand les
        // deux operandes ont exactement le meme type. Licites parce que deux
        // graded de meme plafond ont la meme longueur : l'ordre
        // lexicographique sur les tableaux coincide alors avec l'ordre sur
        // les suites etendues.
        constexpr bool operator==(graded const &) const = default;
        constexpr std::strong_ordering operator<=>(graded const &) const = default;
    };

    namespace detail
    {
        // Type resultat d'une operation binaire.
        //
        // Deux arites exactes de tailles differentes ne satisfont pas la
        // contrainte : c'est cette absence qui produit l'erreur de
        // compilation faisant tout l'interet de l'arite exacte. Des qu'un
        // plafond est en jeu, le resultat est borne -- le melange est a sens
        // unique.
        template <class A, class B, class V>
        struct combine_impl
        {
        };

        template <std::integral T, std::size_t N, wsl::arity AA, std::integral U,
                  std::size_t M, wsl::arity BB, class V>
            requires(AA == wsl::arity::bounded || BB == wsl::arity::bounded || N == M)
        struct combine_impl<graded<T, N, AA>, graded<U, M, BB>, V>
        {
            using type = graded<V, larger(N, M),
                                (AA == wsl::arity::exact && BB == wsl::arity::exact)
                                    ? wsl::arity::exact
                                    : wsl::arity::bounded>;
        };

        template <class A, class B, class V>
        using combine = typename combine_impl<A, B, V>::type;

        template <class A, class B>
        concept mixable = requires { typename combine_impl<A, B, int>::type; };

        // Une seule representation, donc une seule facon de construire.
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
    // Deux arites exactes de tailles differentes ne sont pas comparables :
    // c'est une erreur de compilation, pas un resultat faux. Tout le reste
    // l'est, les operandes etant lus comme des suites etendues.
    //
    // Les conversions arithmetiques usuelles s'appliquent coordonnee par
    // coordonnee : le piege signe/non-signe est reproduit tel quel.
    //
    // Pour l'ordre, ce piege demande un detour. Le langage refuse int <=>
    // unsigned (conversion retrecissante), mais accepte int < unsigned. On
    // procede donc comme std::pair, std::tuple et les conteneurs, par la
    // comparaison synthetisee de la norme (synth-three-way,
    // [expos.only.entity]) : <=> quand il existe, sinon < dans les deux sens.
    // Le resultat est alors un std::weak_ordering, comme pour std::pair.
    // ---------------------------------------------------------------------
    namespace detail
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
        requires detail::mixable<A, B>
    constexpr bool operator==(A const &lhs, B const &rhs) noexcept
    {
        constexpr std::size_t bound =
            detail::larger(detail::capacity_of<A>, detail::capacity_of<B>);
        for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
            if (!(lhs[coordinate] == rhs[coordinate]))
                return false;
        return true;
    }

    template <lexicographic A, lexicographic B>
        requires detail::mixable<A, B>
    constexpr auto operator<=>(A const &lhs, B const &rhs) noexcept
        -> detail::synth_three_way_result<detail::value_of<A>, detail::value_of<B>>
    {
        constexpr std::size_t bound =
            detail::larger(detail::capacity_of<A>, detail::capacity_of<B>);
        for (std::size_t coordinate = 0; coordinate < bound; ++coordinate)
            if (auto const result = detail::synth_three_way(lhs[coordinate], rhs[coordinate]);
                result != 0)
                return result;
        return std::strong_ordering::equal;
    }

    // ---------------------------------------------------------------------
    // Arithmetique.
    // ---------------------------------------------------------------------
    template <lexicographic A, lexicographic B>
        requires detail::mixable<A, B>
    constexpr auto operator+(A const &lhs, B const &rhs) noexcept
        -> detail::combine<A, B, add_t<detail::value_of<A>, detail::value_of<B>>>
    {
        using result_type =
            detail::combine<A, B, add_t<detail::value_of<A>, detail::value_of<B>>>;
        return detail::zip<result_type>(
            lhs, rhs, [](auto left, auto right) noexcept { return left + right; });
    }

    // Rien de particulier : left - right, comme pour deux scalaires. La
    // difference de deux valeurs non signees reste non signee et modulaire,
    // exactement comme 3u - 5u. Consequence utile, a -= b et a = a - b
    // designent la meme chose.
    template <lexicographic A, lexicographic B>
        requires detail::mixable<A, B>
    constexpr auto operator-(A const &lhs, B const &rhs) noexcept
        -> detail::combine<A, B, subtract_t<detail::value_of<A>, detail::value_of<B>>>
    {
        using result_type =
            detail::combine<A, B, subtract_t<detail::value_of<A>, detail::value_of<B>>>;
        return detail::zip<result_type>(
            lhs, rhs, [](auto left, auto right) noexcept { return left - right; });
    }

    template <lexicographic A>
    constexpr auto operator-(A const &operand) noexcept ->
        typename detail::traits<A>::template rebind<negate_t<detail::value_of<A>>>
    {
        using result_type =
            typename detail::traits<A>::template rebind<negate_t<detail::value_of<A>>>;
        return detail::build<result_type>(
            [&](std::size_t coordinate) { return -operand[coordinate]; });
    }

    template <lexicographic A>
    constexpr auto operator+(A const &operand) noexcept ->
        typename detail::traits<A>::template rebind<promote_t<detail::value_of<A>>>
    {
        using result_type =
            typename detail::traits<A>::template rebind<promote_t<detail::value_of<A>>>;
        return detail::build<result_type>(
            [&](std::size_t coordinate) { return +operand[coordinate]; });
    }

    // ---------------------------------------------------------------------
    // Signe.
    //
    // Le signe d'une valeur non nulle est celui de sa coordonnee de rang
    // order() : les precedentes sont nulles par definition, les suivantes
    // sont infiniment moins significatives. Le neutre est de signe nul.
    // Equivalent a comparer la valeur au neutre.
    //
    // Attention pour T non signe : sign(x - y) ne decide pas la comparaison
    // de x et y, la soustraction y etant modulaire. Consequence directe du
    // principe de fidelite au type sous-jacent.
    // ---------------------------------------------------------------------
    template <lexicographic A>
    constexpr int sign(A const &operand) noexcept
    {
        std::size_t const leading = operand.order();
        if (leading == no_support)
            return 0;
        if constexpr (std::is_signed_v<detail::value_of<A>>)
            return operand[leading] < detail::value_of<A>{} ? -1 : 1;
        else
            return 1;
    }

    // ---------------------------------------------------------------------
    // Concatenation.
    //
    // C'est le produit lexicographique : concat(a,b) < concat(a',b') si et
    // seulement si a < a', ou a == a' et b < b'. Autrement dit un
    // isomorphisme de groupes qui preserve l'ordre. En termes de jets,
    // concat(a, b) vaut a + e^N.b -- ce qui explique pourquoi l'operande
    // gauche doit etre d'arite exacte : sans longueur, pas de puissance de e.
    //
    // En arite bornee, l'operation ne passerait d'ailleurs pas au quotient :
    // si {1,0} et {1} sont la meme valeur, les concatener avec {5} devrait
    // rendre {1,0,5} dans un cas et {1,5} dans l'autre. L'operande droit,
    // lui, n'a besoin que d'etre decale, et le decalage respecte
    // l'equivalence -- il peut donc etre de l'une ou l'autre arite, qu'il
    // transmet au resultat.
    //
    // Le type resultat n'est pas invente, il est herite : concat(a, b) vaut
    // a complete de zeros a droite, plus b decale et complete de zeros a
    // gauche. Chaque coordonnee subit donc exactement la conversion que
    // l'addition correspondante lui aurait fait subir.
    // ---------------------------------------------------------------------
    namespace detail
    {
        // Realise litteralement la somme des deux operandes completes de
        // zeros : chaque coordonnee est additionnee au zero de l'autre type,
        // ce qui applique les conversions arithmetiques usuelles sans que la
        // bibliotheque ait a les reproduire.
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
        return detail::concatenate<graded<add_t<T, U>, N + M, B>, N>(lhs, rhs);
    }

    template <lexicographic A, lexicographic B, lexicographic... Rest>
        requires(sizeof...(Rest) != 0)
    constexpr auto concat(A const &first, B const &second, Rest const &...rest) noexcept
    {
        return concat(concat(first, second), rest...);
    }

    // Alias de concat. Depuis C++17 les operandes d'un operateur surcharge
    // sont sequences dans l'ordre prescrit pour l'operateur integre, donc la
    // virgule surchargee conserve la garantie gauche-avant-droite.
    // A employer entre parentheses : la virgule ayant la precedence la plus
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
    // des coefficients, to_jet et to_series l'arite et le plafond.
    //
    // Aucune ne pose de precondition, aucune n'est un comportement indefini.
    // Une troncature explicite est bien definie, sur le modele de
    // static_cast entre entiers : l'appel ecrit a la main est l'annonce.
    //
    // La troncature est un homomorphisme de groupes -- tronquer une somme
    // revient a sommer les troncatures -- mais elle ne preserve pas l'ordre,
    // {1,2,3} et {1,2,0} s'ecrasant sur la meme valeur.
    // ---------------------------------------------------------------------
    template <std::integral To, lexicographic A>
    constexpr auto value_cast(A const &operand) noexcept ->
        typename detail::traits<A>::template rebind<To>
    {
        using result_type = typename detail::traits<A>::template rebind<To>;
        return detail::build<result_type>([&](std::size_t coordinate)
                                          { return static_cast<To>(operand[coordinate]); });
    }

    // Entree dans l'arite exacte, depuis n'importe quelle source. Elargit
    // sans perte quand la taille de destination suffit, tronque sinon.
    template <std::size_t N, lexicographic A>
    constexpr jet<detail::value_of<A>, N> to_jet(A const &operand) noexcept
    {
        return detail::build<jet<detail::value_of<A>, N>>(
            [&](std::size_t coordinate) { return operand[coordinate]; });
    }

    // Entree dans l'arite bornee. La forme sans plafond conserve celui de la
    // source ; elle rend l'injection canonique, qui ne perd jamais rien.
    template <std::size_t N, lexicographic A>
    constexpr series<detail::value_of<A>, N> to_series(A const &operand) noexcept
    {
        return detail::build<series<detail::value_of<A>, N>>(
            [&](std::size_t coordinate) { return operand[coordinate]; });
    }

    template <lexicographic A>
    constexpr series<detail::value_of<A>, detail::capacity_of<A>>
    to_series(A const &operand) noexcept
    {
        return to_series<detail::capacity_of<A>>(operand);
    }

    namespace detail
    {
        // Le hachage porte sur le support, sans quoi des valeurs
        // equivalentes de plafonds differents auraient des empreintes
        // differentes. La coherence est garantie a value_type fixe ; entre
        // types de signes differents, deux valeurs egales au sens des
        // conversions usuelles peuvent se hacher differemment.
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
        return wsl::detail::hash_of(operand);
    }
};

#endif
