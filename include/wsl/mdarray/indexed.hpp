#pragma once

// Monde des axes (docs/mdarray.md, V2 à V4, V9, L1 à L14, I6, I10, I12) :
//
//   c[i, k] = a[i, j] * b[j, k];          // j absent à gauche : sommé (V3)
//   y[i] = tanh(W[i, j] * x[j] + b[i]);   double n = sqrt(x[i] * x[i]);
//   A[i, i] = 1;   y[i, j] *= w[i];       r[i] = sum(j)(a[i, j]);   h[i, j] = 1.0 / (i + j + 1);
//
// Le membre de gauche décide : un axe présent à gauche n'est jamais sommé ; un axe
// absent à gauche est sommé sur la plus petite sous-expression qui contient toutes
// ses occurrences. Garde-fous (V4), à la compilation : (1) un axe absent à gauche
// apparaît au moins deux fois ; (2) pas d'ambiguïté de portée autour d'une fonction ;
// (3) un axe lié par une réduction n'apparaît ni à gauche ni ailleurs. Un scalaire a
// un membre de gauche sans axe (`double s = …`, `wsl::value(…)`).

#include <wsl/mdarray/axis.hpp>
#include <wsl/mdarray/container.hpp>
#include <wsl/mdarray/detail/expression_base.hpp>
#include <wsl/mdarray/detail/function_list.hpp>
#include <wsl/mdarray/detail/loop.hpp>
#include <wsl/mdarray/detail/precondition.hpp>
#include <wsl/mdarray/detail/reduce.hpp>
#include <wsl/mdarray/detail/typelist.hpp>
#include <wsl/mdarray/expression.hpp>
#include <wsl/mdarray/mdspan.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// WSL_MDARRAY_INDEPENDENT_LOOP demande la vectorisation : Clang avertit quand elle est
// impossible (mode reproductible, fonctions), ce qui est attendu (detail/loop.hpp).
#if defined(__clang__)
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wpass-failed"
#endif

namespace wsl
{
    namespace detail::mdarray
    {
        // --- Environnement : position et longueur de chaque axe de l'instruction ---------

        template <class Axes>
        struct environment
        {
            std::array<std::ptrdiff_t, Axes::size> position{};
            std::array<std::ptrdiff_t, Axes::size> length{};
            std::array<bool, Axes::size> known{};

            template <class Id>
            constexpr std::ptrdiff_t at() const noexcept
            {
                return position[tindex_v<Id, Axes>];
            }

            // Longueur d'un axe donnée par une dimension (L5) : toutes doivent coïncider.
            template <class Id>
            constexpr void bind(std::ptrdiff_t n) noexcept
            {
                constexpr std::size_t k = tindex_v<Id, Axes>;
                if (known[k])
                    precondition(length[k] == n, "wsl::mdarray : un axe étiquette des dimensions de longueurs différentes (L5)");
                else
                {
                    length[k] = n;
                    known[k] = true;
                }
            }
        };

        // Contexte statique d'une instruction : axes, axes sommés, occurrences du membre de droite.
        template <class Axes, class Summed, class Root>
        struct context
        {
            using axes = Axes;
            using summed = Summed;
            using root = Root;
        };

        // Compteur à roues sur K axes de l'environnement : la k-ième position dans l'ordre
        // lexicographique, par incrément quand k suit le précédent (sans division), par
        // division sinon (saut de la somme par paires).
        template <std::size_t K>
        struct odometer
        {
            std::array<std::size_t, K> index{};
            std::array<std::ptrdiff_t, K> length{};
            std::size_t last = static_cast<std::size_t>(-1);

            template <class Env>
            constexpr odometer(std::array<std::size_t, K> const& indices, Env const& env) noexcept : index(indices)
            {
                for (std::size_t r = 0; r < K; ++r)
                    length[r] = env.length[index[r]];
            }

            constexpr std::size_t size() const noexcept
            {
                std::size_t n = 1;
                for (std::size_t r = 0; r < K; ++r)
                    n *= static_cast<std::size_t>(length[r]);
                return n;
            }

            template <class Env>
            constexpr void set(std::size_t k, Env& env) noexcept
            {
                step(k, [&](std::size_t r) -> std::ptrdiff_t& { return env.position[index[r]]; });
            }

            // Même parcours, dans des positions locales (accès à pas fixe).
            constexpr void set_local(std::size_t k) noexcept
            {
                step(k, [&](std::size_t r) -> std::ptrdiff_t& { return local[r]; });
            }

            std::array<std::ptrdiff_t, K> local{};

        private:
            template <class Position>
            constexpr void step(std::size_t k, Position position) noexcept
            {
                if constexpr (K == 1)
                    position(0) = static_cast<std::ptrdiff_t>(k);
                else
                {
                    if (k == last + 1 && last != static_cast<std::size_t>(-1))
                    {
                        std::size_t r = K;
                        while (r-- > 0)
                        {
                            if (++position(r) < length[r])
                                break;
                            position(r) = 0;
                        }
                    }
                    else
                    {
                        std::size_t q = k;
                        for (std::size_t r = K; r-- > 0;)
                        {
                            auto const len = static_cast<std::size_t>(length[r]);
                            position(r) = static_cast<std::ptrdiff_t>(q % len);
                            q /= len;
                        }
                    }
                    last = k;
                }
            }
        };

        // --- Emplacements d'un crochet ------------------------------------------------

        template <class Id>
        struct axis_slot
        {
            using id = Id;
        };
        struct fixed_slot
        {
            std::ptrdiff_t value;
        };

        template <class S>
        struct slot_occurrence
        {
            using type = tlist<>;
        };
        template <class Id>
        struct slot_occurrence<axis_slot<Id>>
        {
            using type = tlist<Id>;
        };

        template <class A>
        constexpr auto make_slot(A const& a) noexcept
        {
            if constexpr (axis_like<A>)
                return axis_slot<typename std::remove_cvref_t<A>::id_type>{};
            else
                return fixed_slot{static_cast<std::ptrdiff_t>(a)};
        }

        template <class Env>
        constexpr std::ptrdiff_t slot_value(fixed_slot const& s, Env const&) noexcept
        {
            return s.value;
        }
        template <class Id, class Env>
        constexpr std::ptrdiff_t slot_value(axis_slot<Id> const&, Env const& env) noexcept
        {
            return env.template at<Id>();
        }

        // --- Base des nœuds --------------------------------------------------------------

        template <class Node>
        struct no_children
        {
            using children = tlist<>;
        };

        template <class Node, class Ctx>
        struct home;

        template <class Ctx, class Node, class Env>
        constexpr auto evaluate(Node const& node, Env& env);

        // Le terme k d'une somme le long des axes Summed (à pas fixe si possible) ; Root :
        // node est le foyer lui-même (terme : node.compute), sinon un enfant (evaluate).
        template <class Ctx, class Summed, bool Root, class Node, class Env, std::size_t K>
        constexpr auto terms(Node const& node, Env& env, odometer<K>& wheels);

        // Évalue node avec un membre de gauche sans axe (scalaire, V9).
        template <class Node>
        constexpr auto evaluate_scalar(Node const& node);

        template <class Node>
        struct indexed_base : indexed_tag
        {
            // `double s = …` : membre de gauche sans axe ; temporaire seulement (L2).
            template <class T>
                requires(std::is_arithmetic_v<T> || std::is_class_v<T>) && (!indexed<T>)
            constexpr operator T() &&
            {
                return static_cast<T>(evaluate_scalar(default_mode(static_cast<Node&&>(*this))));
            }

            template <class T>
                requires(std::is_arithmetic_v<T> || std::is_class_v<T>) && (!indexed<T>)
            operator T() const&
            {
                static_assert(sizeof(T) == 0, "wsl::mdarray : une expression ne se garde pas (L2) ; écrire le calcul dans "
                                              "l'instruction, ou std::move(e)");
                return T{};
            }

            // max(i, j)(…), (a[i] + b[i])(…) : une expression ne s'appelle pas (I6).
            template <class... A>
            void operator()(A&&...) const
            {
                static_assert(sizeof...(A) < 0, "wsl::mdarray : une expression ne s'appelle pas ; pour réduire le long "
                                                "d'axes, écrire sum(j)(…), et pour max sur deux axes max(i)(max(j)(…))");
            }
        };

        // --- Feuilles ------------------------------------------------------------------

        template <class T>
        struct indexed_scalar : indexed_base<indexed_scalar<T>>, no_children<indexed_scalar<T>>
        {
            using occurrences = tlist<>;
            using bound = tlist<>;
            using provided = tlist<>;
            using static_lengths = tlist<>;
            using value_type = T;
            static constexpr bool is_function = false;
            T value;
            constexpr explicit indexed_scalar(T v) : value(std::move(v)) {}
            indexed_scalar(indexed_scalar const&) = delete;
            constexpr indexed_scalar(indexed_scalar&&) = default;
            template <class Ctx, class Env>
            constexpr T compute(Env&) const
            {
                return value;
            }
            template <class Env>
            constexpr void lengths(Env&) const noexcept
            {
            }
            template <class F>
            constexpr void for_each_leaf(F&&) const
            {
            }
        };

        // Un axe pris comme valeur : ses positions (L12).
        template <class Id, class T>
        struct indexed_axis : indexed_base<indexed_axis<Id, T>>, no_children<indexed_axis<Id, T>>
        {
            using occurrences = tlist<Id>;
            using bound = tlist<>;
            using provided = tlist<>;
            using static_lengths = tlist<>;
            using value_type = T;
            static constexpr bool is_function = false;
            constexpr indexed_axis() = default;
            indexed_axis(indexed_axis const&) = delete;
            constexpr indexed_axis(indexed_axis&&) = default;
            template <class Ctx, class Env>
            constexpr T compute(Env& env) const
            {
                return static_cast<T>(env.template at<Id>());
            }
            template <class Env>
            constexpr void lengths(Env&) const noexcept
            {
            }
            template <class F>
            constexpr void for_each_leaf(F&&) const
            {
            }
        };

        // δᵢⱼ (L13).
        template <class A, class B, class T>
        struct indexed_kronecker : indexed_base<indexed_kronecker<A, B, T>>, no_children<indexed_kronecker<A, B, T>>
        {
            using occurrences = tlist<A, B>;
            using bound = tlist<>;
            using provided = tlist<>;
            using static_lengths = tlist<>;
            using value_type = T;
            static constexpr bool is_function = false;
            constexpr indexed_kronecker() = default;
            indexed_kronecker(indexed_kronecker const&) = delete;
            constexpr indexed_kronecker(indexed_kronecker&&) = default;
            template <class Ctx, class Env>
            constexpr T compute(Env& env) const
            {
                return static_cast<T>(env.template at<A>() == env.template at<B>() ? 1 : 0);
            }
            template <class Env>
            constexpr void lengths(Env&) const noexcept
            {
            }
            template <class F>
            constexpr void for_each_leaf(F&&) const
            {
            }
        };

        template <class Id, std::size_t E>
        struct static_length
        {
            using id = Id;
            static constexpr std::size_t value = E;
        };

        template <class Source, class Seq, class... Slots>
        struct static_lengths_of;
        template <class Source, std::size_t... R, class... Slots>
        struct static_lengths_of<Source, std::index_sequence<R...>, Slots...>
        {
            template <std::size_t Rank, class S>
            struct one
            {
                using type = tlist<>;
            };
            template <std::size_t Rank, class Id>
            struct one<Rank, axis_slot<Id>>
            {
                static constexpr std::size_t e = Source::extents_type::static_extent(Rank);
                using type = std::conditional_t<e == dynamic_extent, tlist<>, tlist<static_length<Id, e>>>;
            };
            using type = tconcat_t<typename one<R, Slots>::type...>;
        };

        // Un tableau, une vue ou une expression de tableaux entiers, indicé : m[i, 0, k].
        template <class Source, class... Slots>
        struct indexed_leaf : indexed_base<indexed_leaf<Source, Slots...>>, no_children<indexed_leaf<Source, Slots...>>
        {
            using occurrences = tconcat_t<typename slot_occurrence<Slots>::type...>;
            using bound = tlist<>;
            using provided = tunique_t<occurrences>;
            using static_lengths = typename static_lengths_of<Source, std::index_sequence_for<Slots...>, Slots...>::type;
            using value_type = std::remove_cvref_t<decltype(std::declval<Source const&>()(std::declval<std::conditional_t<true, std::ptrdiff_t, Slots>>()...))>;
            static constexpr bool is_function = false;

            Source source;
            std::tuple<Slots...> slots;

            constexpr indexed_leaf(Source s, Slots... sl) : source(std::move(s)), slots(sl...) {}
            indexed_leaf(indexed_leaf const&) = delete;
            constexpr indexed_leaf(indexed_leaf&&) = default;

            template <class Ctx, class Env>
            constexpr value_type compute(Env& env) const
            {
                return std::apply([&](auto const&... s) { return source(slot_value(s, env)...); }, slots);
            }

            template <class Env>
            constexpr void lengths(Env& env) const noexcept
            {
                [&]<std::size_t... R>(std::index_sequence<R...>) {
                    (bind_slot<R>(env, std::get<R>(slots)), ...);
                }(std::index_sequence_for<Slots...>{});
            }

            template <class F>
            constexpr void for_each_leaf(F&& f) const
            {
                f(*this);
            }

        private:
            template <std::size_t R, class Env>
            constexpr void bind_slot(Env&, fixed_slot const& s) const noexcept
            {
                precondition(s.value >= 0 && s.value < static_cast<std::ptrdiff_t>(source.extents().extent(R)),
                             "wsl::mdarray : position fixe hors du tableau (L3)");
            }
            template <std::size_t R, class Env, class Id>
            constexpr void bind_slot(Env& env, axis_slot<Id> const&) const noexcept
            {
                env.template bind<Id>(static_cast<std::ptrdiff_t>(source.extents().extent(R)));
            }
        };

        // Source « vue » : m[i…] lu dans une vue.
        template <class View>
        struct view_source
        {
            View view;
            using extents_type = typename View::extents_type;
            constexpr extents_type const& extents() const noexcept { return view.extents(); }
            template <class... I>
            constexpr decltype(auto) operator()(I... i) const
            {
                return view[static_cast<typename View::index_type>(i)...];
            }
        };

        // Source « expression » : (a + b)[i, j].
        template <class Expr>
        struct expression_source
        {
            Expr expr;
            using extents_type = typename Expr::extents_type;
            constexpr extents_type const& extents() const noexcept { return expr.extents(); }
            template <class... I>
            constexpr auto operator()(I... i) const
            {
                return expr(static_cast<typename Expr::index_type>(i)...);
            }
        };

        // --- Nœuds ---------------------------------------------------------------------

        template <class Op>
        inline constexpr bool op_is_function = requires { Op::is_function; };

        template <class Op, class... Children>
        struct indexed_node : indexed_base<indexed_node<Op, Children...>>
        {
            using children = tlist<Children...>;
            using occurrences = tconcat_t<typename Children::occurrences...>;
            using bound = tconcat_t<typename Children::bound...>;
            using provided = tunique_t<tconcat_t<typename Children::provided...>>;
            using static_lengths = tconcat_t<typename Children::static_lengths...>;
            using value_type =
                std::remove_cvref_t<decltype(std::declval<Op const&>()(std::declval<typename Children::value_type>()...))>;
            static constexpr bool is_function = op_is_function<Op>;

            [[no_unique_address]] Op op;
            std::tuple<Children...> nodes;

            constexpr indexed_node(Op o, Children... c) : op(std::move(o)), nodes(std::move(c)...) {}
            indexed_node(indexed_node const&) = delete;
            constexpr indexed_node(indexed_node&&) = default;

            template <class Ctx, class Env>
            constexpr auto compute(Env& env) const
            {
                return std::apply([&](auto const&... c) { return op(evaluate<Ctx>(c, env)...); }, nodes);
            }

            template <class Env>
            constexpr void lengths(Env& env) const noexcept
            {
                std::apply([&](auto const&... c) { (c.lengths(env), ...); }, nodes);
            }

            template <class F>
            constexpr void for_each_leaf(F&& f) const
            {
                std::apply([&](auto const&... c) { (c.for_each_leaf(f), ...); }, nodes);
            }
        };

        // Réductions explicites le long d'axes liés (L8, I6).
        struct sum_kind
        {
        };
        struct prod_kind
        {
        };
        struct min_kind
        {
        };
        struct max_kind
        {
        };
        struct all_kind
        {
        };
        struct any_kind
        {
        };
        struct count_kind
        {
        };
        struct argmin_kind
        {
        };
        struct argmax_kind
        {
        };

        template <class Kind, class Bound, class Child>
        struct indexed_reduction : indexed_base<indexed_reduction<Kind, Bound, Child>>
        {
            using children = tlist<Child>;
            using occurrences = tminus_t<typename Child::occurrences, Bound>;
            using bound = tconcat_t<typename Child::bound, Bound>;
            using provided = typename Child::provided;
            using static_lengths = typename Child::static_lengths;
            using child_value = typename Child::value_type;
            using value_type = std::conditional_t<
                std::is_same_v<Kind, all_kind> || std::is_same_v<Kind, any_kind>, bool,
                std::conditional_t<std::is_same_v<Kind, count_kind> || std::is_same_v<Kind, argmin_kind>
                                       || std::is_same_v<Kind, argmax_kind>,
                                   std::ptrdiff_t, child_value>>;
            static constexpr bool is_function = false;

            Child child;

            constexpr explicit indexed_reduction(Child c) : child(std::move(c)) {}
            indexed_reduction(indexed_reduction const&) = delete;
            constexpr indexed_reduction(indexed_reduction&&) = default;

            template <class Ctx, class Env>
            constexpr auto compute(Env& env) const
            {
                constexpr std::size_t K = Bound::size;
                odometer<K> wheels(bound_indices<typename Ctx::axes>(Bound{}), env);
                std::size_t const n = wheels.size();
                auto const at = terms<Ctx, Bound, false>(child, env, wheels);
                if constexpr (std::is_same_v<Kind, sum_kind>)
                {
                    using T = decltype(at(0) + at(0));
                    return pairwise<T>(at, 0, n, T{}, detail::mdarray::add_op{});
                }
                else if constexpr (std::is_same_v<Kind, prod_kind>)
                {
                    using T = decltype(at(0) * at(0));
                    return pairwise<T>(at, 0, n, T{1}, detail::mdarray::multiply_op{});
                }
                else if constexpr (std::is_same_v<Kind, min_kind> || std::is_same_v<Kind, max_kind>)
                {
                    precondition(n > 0, "min, max : axe de longueur nulle");
                    auto m = at(0);
                    for (std::size_t k = 1; k < n; ++k)
                        m = std::conditional_t<std::is_same_v<Kind, min_kind>, min_op, max_op>{}(m, at(k));
                    return m;
                }
                else if constexpr (std::is_same_v<Kind, all_kind>)
                {
                    for (std::size_t k = 0; k < n; ++k)
                        if (!static_cast<bool>(at(k)))
                            return false;
                    return true;
                }
                else if constexpr (std::is_same_v<Kind, any_kind>)
                {
                    for (std::size_t k = 0; k < n; ++k)
                        if (static_cast<bool>(at(k)))
                            return true;
                    return false;
                }
                else if constexpr (std::is_same_v<Kind, count_kind>)
                {
                    std::ptrdiff_t c = 0;
                    for (std::size_t k = 0; k < n; ++k)
                        c += static_cast<bool>(at(k)) ? 1 : 0;
                    return c;
                }
                else
                {
                    static_assert(K == 1, "wsl::mdarray : argmin, argmax le long d'un seul axe");
                    precondition(n > 0, "argmin, argmax : axe de longueur nulle");
                    std::size_t best = 0;
                    auto v = at(0);
                    for (std::size_t k = 1; k < n && v == v; ++k)
                    {
                        auto const w = at(k);
                        bool const better = std::is_same_v<Kind, argmin_kind> ? w < v : v < w;
                        if (w != w || better)
                        {
                            best = k;
                            v = w;
                        }
                    }
                    return static_cast<std::ptrdiff_t>(best);
                }
            }

            template <class Env>
            constexpr void lengths(Env& env) const noexcept
            {
                child.lengths(env);
            }

            template <class F>
            constexpr void for_each_leaf(F&& f) const
            {
                child.for_each_leaf(f);
            }

        private:
            template <class Axes, class... Ids>
            static constexpr std::array<std::size_t, sizeof...(Ids)> bound_indices(tlist<Ids...>) noexcept
            {
                return {tindex_v<Ids, Axes>...};
            }
        };

        // --- Sommation aux nœuds « foyers » (V3) ------------------------------------------

        // Axes sommés dont node est le foyer : toutes leurs occurrences y sont, et aucun
        // enfant ne les contient toutes.
        template <class Id, class Node, class Ctx>
        struct is_home
        {
            static constexpr std::size_t total = tcount<Id, typename Ctx::root>;
            template <class C>
            struct below
            {
                static constexpr bool value = tcount<Id, typename C::occurrences> < total;
            };
            static constexpr bool value =
                tcount<Id, typename Node::occurrences> == total && tall<typename Node::children, below>;
        };

        template <class Node, class Ctx, class Summed>
        struct home_impl;
        template <class Node, class Ctx, class... S>
        struct home_impl<Node, Ctx, tlist<S...>>
        {
            using type = tconcat_t<std::conditional_t<is_home<S, Node, Ctx>::value, tlist<S>, tlist<>>...>;
        };

        template <class Node, class Ctx>
        struct home
        {
            using type = typename home_impl<Node, Ctx, typename Ctx::summed>::type;
        };

        template <class Axes, class... Ids>
        constexpr std::array<std::size_t, sizeof...(Ids)> axis_indices(tlist<Ids...>) noexcept
        {
            return {tindex_v<Ids, Axes>...};
        }

        // --- Accès à pas fixe -----------------------------------------------------------
        // Une somme le long d'axes Summed dont le sous-arbre ne lit que des vues à pas
        // (layout right, left, stride ; accesseur par défaut), des scalaires, des axes et
        // des δ, sans réduction explicite : chaque feuille calcule une fois son décalage et
        // ses pas le long des axes sommés, puis lit p[o + Σ posᵣ · pasᵣ], sans passer par
        // l'environnement. Un foyer imbriqué (g[p, i, j] v[p, i] v[p, j] : i sommé dans le
        // premier produit) somme de même ses propres axes, ajoutés aux positions. Mêmes
        // éléments, mêmes opérations, même ordre : résultats identiques à l'évaluation
        // générale.

        template <class View, class... Slots>
        struct indexed_target;

        template <class Ctx, class Node>
        struct strideable : std::false_type
        {
        };
        template <class Ctx, class V, class... S>
        struct strideable<Ctx, indexed_leaf<view_source<V>, S...>>
            : std::bool_constant<std::is_same_v<typename V::accessor_type, default_accessor<typename V::element_type>>
                                 && V::is_always_strided()>
        {
        };
        template <class Ctx, class V, class... S>
        struct strideable<Ctx, indexed_target<V, S...>> : strideable<Ctx, indexed_leaf<view_source<V>, S...>>
        {
        };
        template <class Ctx, class T>
        struct strideable<Ctx, indexed_scalar<T>> : std::true_type
        {
        };
        template <class Ctx, class Id, class T>
        struct strideable<Ctx, indexed_axis<Id, T>> : std::true_type
        {
        };
        template <class Ctx, class A, class B, class T>
        struct strideable<Ctx, indexed_kronecker<A, B, T>> : std::true_type
        {
        };
        template <class Ctx, class Op, class... C>
        struct strideable<Ctx, indexed_node<Op, C...>> : std::bool_constant<(strideable<Ctx, C>::value && ...)>
        {
        };

        template <std::size_t N>
        using positions = std::array<std::ptrdiff_t, N>;

        // Chaque plan a value(pos), sa valeur aux positions pos des axes de Summed, et
        // shift(pos), qui y ajoute une fois pour toutes la contribution de pos : une somme
        // imbriquée décale son corps à son entrée, puis le parcourt avec des positions
        // extérieures nulles (que le compilateur voit nulles).

        template <class Summed, class... Ids>
        constexpr std::array<std::size_t, sizeof...(Ids)> summed_indices(tlist<Ids...>) noexcept
        {
            return {tindex_v<Ids, Summed>...};
        }

        // Voies : W positions consécutives le long de l'axe d'indice I des positions,
        // calculées côte à côte (C[i, k…k+7] = A[i, j] B[j, k…k+7]) ; chaque voie suit
        // exactement les opérations et l'ordre de la voie scalaire : mêmes bits, mais des
        // lectures contiguës que le compilateur vectorise.
        // Somme voie à voie, convertie dans le type T de la somme.
        template <class T, std::size_t W>
        struct lane_add
        {
            template <class A, class B>
            constexpr lane_pack<T, W> operator()(lane_pack<A, W> const& a, lane_pack<B, W> const& b) const
            {
                lane_pack<T, W> r;
                for (std::size_t l = 0; l < W; ++l)
                    r.v[l] = detail::mdarray::add_op{}(a.v[l], b.v[l]);
                return r;
            }
        };

        // Une position d'axe : base + pos[r] le long d'un axe de Summed, base sinon.
        template <class Summed, class Id>
        struct axis_reader
        {
            std::ptrdiff_t base = 0;
            constexpr std::ptrdiff_t value(positions<Summed::size> const& pos) const noexcept
            {
                if constexpr (tcontains<Id, Summed>)
                    return base + pos[tindex_v<Id, Summed>];
                else
                    return base;
            }
            constexpr void shift(positions<Summed::size> const& pos) noexcept
            {
                if constexpr (tcontains<Id, Summed>)
                    base += pos[tindex_v<Id, Summed>];
            }
            // Avance d'une voie à la suivante.
            template <std::size_t I>
            static constexpr std::ptrdiff_t lane_step() noexcept
            {
                if constexpr (tcontains<Id, Summed>)
                    return tindex_v<Id, Summed> == I ? 1 : 0;
                else
                    return 0;
            }
        };

        template <class Summed, class Id, class Env>
        constexpr axis_reader<Summed, Id> make_reader(Env const& env) noexcept
        {
            if constexpr (tcontains<Id, Summed>)
                return {};
            else
                return {env.template at<Id>()};
        }

        // Accès d'une feuille : p[o + Σ posᵣ · pasᵣ], sur les axes de Summed présents
        // dans son crochet (Present, dans l'ordre de Summed). Unit[q] : pas connu égal à 1
        // à la compilation (axe seul à la dimension contiguë de layout_right ou
        // layout_left), pour que le compilateur voie les accès contigus et vectorise.
        template <class Summed, class Present, class Pointer, std::array<bool, Present::size> Unit>
        struct strided_access
        {
            static constexpr std::size_t P = Present::size;
            static constexpr std::array<std::size_t, P> where = summed_indices<Summed>(Present{});
            Pointer p;
            std::ptrdiff_t offset = 0;
            std::array<std::ptrdiff_t, P> stride{};

            constexpr std::ptrdiff_t index(positions<Summed::size> const& pos) const noexcept
            {
                return [&]<std::size_t... Q>(std::index_sequence<Q...>) {
                    return (offset + ... + (Unit[Q] ? pos[where[Q]] : pos[where[Q]] * stride[Q]));
                }(std::make_index_sequence<P>{});
            }
            constexpr decltype(auto) operator()(positions<Summed::size> const& pos) const noexcept
            {
                return p[index(pos)];
            }
            constexpr void shift(positions<Summed::size> const& pos) noexcept
            {
                offset = index(pos);
            }
            // Recopie (indexed_target::tiled) : l'accès ne dépend-il que des positions
            // d'indices I et J (et de J) ? emplacement du pas le long d'une position.
            static constexpr std::size_t positions_size = Summed::size;
            static constexpr bool uses(std::size_t i) noexcept
            {
                for (std::size_t q = 0; q < P; ++q)
                    if (where[q] == i)
                        return true;
                return false;
            }
            // Ne dépend que des positions i et [lo, hi), et de toutes celles de [lo, hi).
            static constexpr bool only(std::size_t i, std::size_t lo, std::size_t hi) noexcept
            {
                for (std::size_t q = 0; q < P; ++q)
                    if (where[q] != i && (where[q] < lo || where[q] >= hi))
                        return false;
                for (std::size_t r = lo; r < hi; ++r)
                    if (!uses(r))
                        return false;
                return true;
            }
            static constexpr std::array<bool, P> unit = Unit;
            static constexpr std::size_t slot(std::size_t i) noexcept
            {
                for (std::size_t q = 0; q < P; ++q)
                    if (where[q] == i)
                        return q;
                return P;
            }

            // Pas d'une voie à la suivante (axe d'indice I des positions).
            template <std::size_t I>
            constexpr std::ptrdiff_t lane_step() const noexcept
            {
                return [&]<std::size_t... Q>(std::index_sequence<Q...>) {
                    return (std::ptrdiff_t{0} + ... + (where[Q] == I ? (Unit[Q] ? 1 : stride[Q]) : 0));
                }(std::make_index_sequence<P>{});
            }
        };

        // Dimension contiguë d'une vue (pas 1 quelles que soient les longueurs), ou rang.
        template <class V>
        constexpr std::size_t contiguous_rank() noexcept
        {
            if constexpr (V::rank() == 0)
                return 0;
            else if constexpr (std::is_same_v<typename V::layout_type, layout_right>)
                return V::rank() - 1;
            else if constexpr (std::is_same_v<typename V::layout_type, layout_left>)
                return 0;
            else
                return V::rank();
        }

        // L'axe Id n'occupe qu'une dimension du crochet, la dimension contiguë.
        template <class V, class Id, class... S>
        constexpr bool unit_stride() noexcept
        {
            constexpr std::array<bool, sizeof...(S)> is_id{std::is_same_v<S, axis_slot<Id>>...};
            std::size_t count = 0, rank = 0;
            for (std::size_t r = 0; r < is_id.size(); ++r)
                if (is_id[r])
                {
                    ++count;
                    rank = r;
                }
            return count == 1 && rank == contiguous_rank<V>();
        }

        template <class V, class Present, class... S>
        constexpr std::array<bool, Present::size> unit_strides() noexcept
        {
            return []<class... Ids>(tlist<Ids...>) { return std::array<bool, sizeof...(Ids)>{unit_stride<V, Ids, S...>()...}; }(Present{});
        }

        template <class Summed, class V, class... S, class Env>
        constexpr auto make_access(indexed_leaf<view_source<V>, S...> const& leaf, Env const& env) noexcept
        {
            using present = tminus_t<Summed, tminus_t<Summed, typename indexed_leaf<view_source<V>, S...>::provided>>;
            auto const& view = leaf.source.view;
            strided_access<Summed, present, typename V::data_handle_type, unit_strides<V, present, S...>()> access{view.data_handle()};
            auto place = [&]<std::size_t R, class Slot>(Slot const& slot) {
                auto const sr = static_cast<std::ptrdiff_t>(view.mapping().stride(R));
                if constexpr (std::is_same_v<Slot, fixed_slot>)
                    access.offset += slot.value * sr;
                else if constexpr (tcontains<typename Slot::id, Summed>)
                    access.stride[tindex_v<typename Slot::id, present>] += sr;
                else
                    access.offset += env.template at<typename Slot::id>() * sr;
            };
            [&]<std::size_t... R>(std::index_sequence<R...>) {
                (place.template operator()<R>(std::get<R>(leaf.slots)), ...);
            }(std::index_sequence_for<S...>{});
            return access;
        }

        template <class Access, class T>
        struct strided_leaf
        {
            static constexpr bool reads_axes = false;
            Access access;
            template <class Pos>
            constexpr T value(Pos const& pos) const noexcept
            {
                return access(pos);
            }
            template <class Pos>
            constexpr void shift(Pos const& pos) noexcept
            {
                access.shift(pos);
            }
            template <class F>
            constexpr void for_each_access(F&& f)
            {
                f(access);
            }
            template <std::size_t I, std::size_t W, class Pos>
            constexpr lane_pack<T, W> lanes(Pos const& pos) const noexcept
            {
                std::ptrdiff_t const o = access.index(pos), step = access.template lane_step<I>();
                lane_pack<T, W> r;
                for (std::size_t l = 0; l < W; ++l)
                    r.v[l] = access.p[o + static_cast<std::ptrdiff_t>(l) * step];
                return r;
            }
        };

        template <class T>
        struct strided_scalar
        {
            static constexpr bool reads_axes = false;
            T v;
            template <class Pos>
            constexpr T value(Pos const&) const
            {
                return v;
            }
            template <class Pos>
            constexpr void shift(Pos const&) noexcept
            {
            }
            template <std::size_t I, std::size_t W, class Pos>
            constexpr lane_pack<T, W> lanes(Pos const&) const
            {
                lane_pack<T, W> r;
                for (std::size_t l = 0; l < W; ++l)
                    r.v[l] = v;
                return r;
            }
            template <class F>
            constexpr void for_each_access(F&&)
            {
            }
        };

        template <class T, class Reader>
        struct strided_axis
        {
            static constexpr bool reads_axes = true;
            Reader r;
            template <class Pos>
            constexpr T value(Pos const& pos) const noexcept
            {
                return static_cast<T>(r.value(pos));
            }
            template <class Pos>
            constexpr void shift(Pos const& pos) noexcept
            {
                r.shift(pos);
            }
            template <std::size_t I, std::size_t W, class Pos>
            constexpr lane_pack<T, W> lanes(Pos const& pos) const noexcept
            {
                lane_pack<T, W> out;
                for (std::size_t l = 0; l < W; ++l)
                    out.v[l] = static_cast<T>(r.value(pos) + static_cast<std::ptrdiff_t>(l) * Reader::template lane_step<I>());
                return out;
            }
            template <class F>
            constexpr void for_each_access(F&&)
            {
            }
        };

        template <class T, class RA, class RB>
        struct strided_kronecker
        {
            static constexpr bool reads_axes = true;
            RA a;
            RB b;
            template <class Pos>
            constexpr T value(Pos const& pos) const noexcept
            {
                return static_cast<T>(a.value(pos) == b.value(pos) ? 1 : 0);
            }
            template <class Pos>
            constexpr void shift(Pos const& pos) noexcept
            {
                a.shift(pos);
                b.shift(pos);
            }
            template <std::size_t I, std::size_t W, class Pos>
            constexpr lane_pack<T, W> lanes(Pos const& pos) const noexcept
            {
                lane_pack<T, W> out;
                for (std::size_t l = 0; l < W; ++l)
                {
                    auto const d = static_cast<std::ptrdiff_t>(l);
                    out.v[l] = static_cast<T>(a.value(pos) + d * RA::template lane_step<I>() == b.value(pos) + d * RB::template lane_step<I>()
                                                  ? 1
                                                  : 0);
                }
                return out;
            }
            template <class F>
            constexpr void for_each_access(F&&)
            {
            }
        };

        template <class Op, class... C>
        struct strided_node
        {
            static constexpr bool reads_axes = (C::reads_axes || ...);
            Op const* op;
            std::tuple<C...> children;
            template <class Pos>
            constexpr auto value(Pos const& pos) const
            {
                return std::apply([&](auto const&... c) { return (*op)(c.value(pos)...); }, children);
            }
            template <class Pos>
            constexpr void shift(Pos const& pos) noexcept
            {
                std::apply([&](auto&... c) { (c.shift(pos), ...); }, children);
            }
            template <class F>
            constexpr void for_each_access(F&& f)
            {
                std::apply([&](auto&... c) { (c.for_each_access(f), ...); }, children);
            }
            template <std::size_t I, std::size_t W, class Pos>
            constexpr auto lanes(Pos const& pos) const
            {
                return std::apply(
                    [&](auto const&... c) {
                        using R = decltype((*op)(c.value(pos)...));
                        return std::apply(
                            [&](auto const&... x) {
                                // Une opération qui traite un paquet entier (produit
                                // reproductible : barrière vectorielle), sinon voie par voie.
                                if constexpr (requires { op->lanes(x...); })
                                    return op->lanes(x...);
                                else
                                {
                                    lane_pack<R, W> r;
                                    for (std::size_t l = 0; l < W; ++l)
                                        r.v[l] = (*op)(x.v[l]...);
                                    return r;
                                }
                            },
                            std::tuple{c.template lanes<I, W>(pos)...});
                    },
                    children);
            }
        };

        template <class Ctx, class Summed, class Node, class Env>
        constexpr auto make_strided(Node const& node, Env const& env);

        template <class Ctx, class Summed, class V, class... S, class Env>
        constexpr auto strided_body(indexed_leaf<view_source<V>, S...> const& leaf, Env const& env) noexcept
        {
            using access = decltype(make_access<Summed>(leaf, env));
            return strided_leaf<access, typename indexed_leaf<view_source<V>, S...>::value_type>{make_access<Summed>(leaf, env)};
        }

        template <class Ctx, class Summed, class T, class Env>
        constexpr auto strided_body(indexed_scalar<T> const& s, Env const&)
        {
            return strided_scalar<T>{s.value};
        }

        template <class Ctx, class Summed, class Id, class T, class Env>
        constexpr auto strided_body(indexed_axis<Id, T> const&, Env const& env) noexcept
        {
            return strided_axis<T, axis_reader<Summed, Id>>{make_reader<Summed, Id>(env)};
        }

        template <class Ctx, class Summed, class A, class B, class T, class Env>
        constexpr auto strided_body(indexed_kronecker<A, B, T> const&, Env const& env) noexcept
        {
            return strided_kronecker<T, axis_reader<Summed, A>, axis_reader<Summed, B>>{make_reader<Summed, A>(env),
                                                                                         make_reader<Summed, B>(env)};
        }

        template <class Ctx, class Summed, class Op, class... C, class Env>
        constexpr auto strided_body(indexed_node<Op, C...> const& node, Env const& env)
        {
            auto children = std::apply([&](auto const&... c) { return std::tuple{make_strided<Ctx, Summed>(c, env)...}; }, node.nodes);
            return std::apply([&](auto&... c) { return strided_node<Op, std::remove_cvref_t<decltype(c)>...>{&node.op, children}; },
                              children);
        }

        template <class T>
        constexpr void check_not_bool_sum() noexcept
        {
            static_assert(!std::is_same_v<std::remove_cvref_t<T>, bool>,
                          "wsl::mdarray : somme de booléens (I10) ; écrire count(i)(…), any(i)(…) ou all(i)(…)");
        }

        // Longueur statique d'un axe dans une liste de static_length, ou dynamic_extent.
        template <class Id, class L>
        struct static_length_in : std::integral_constant<std::size_t, dynamic_extent>
        {
        };
        template <class Id, class X, class... Rest>
        struct static_length_in<Id, tlist<X, Rest...>>
            : std::conditional_t<std::is_same_v<typename X::id, Id>, std::integral_constant<std::size_t, X::value>,
                                 static_length_in<Id, tlist<Rest...>>>
        {
        };

        // Termes d'une somme le long d'un seul axe, d'indice S des positions : le k-ième, ou
        // 8 consécutifs à la fois (voies le long de l'axe sommé, pour pairwise).
        template <std::size_t S, class Body>
        struct axis_terms
        {
            Body body;
            constexpr auto operator()(std::size_t k) const
            {
                positions<S + 1> q{};
                q[S] = static_cast<std::ptrdiff_t>(k);
                return body.value(q);
            }
            template <std::size_t W>
            constexpr auto pack(std::size_t k) const
            {
                positions<S + 1> q{};
                q[S] = static_cast<std::ptrdiff_t>(k);
                return body.template lanes<S, W>(q);
            }
        };

        // Somme d'un corps le long des K derniers axes de ses positions (S extérieurs) ;
        // Static : longueurs connues à la compilation (dynamic_extent sinon), pour que les
        // petites sommes (3 × 3) soient déroulées.
        template <std::size_t S, std::size_t K, class Body, std::array<std::size_t, K> Static>
        struct strided_home
        {
            static constexpr bool reads_axes = Body::reads_axes;
            Body body;
            std::array<std::ptrdiff_t, K> length{};
            // Somme aplatie (flatten) : les termes ne se distinguent que par l'indice aplati,
            // porté par la dernière position.
            bool flat = false;

            [[gnu::always_inline]] constexpr auto value(positions<S> const& pos) const
            {
                using V = decltype(body.value(positions<S + K>{}));
                using T = decltype(std::declval<V>() + std::declval<V>());
                Body b = body;
                b.shift(extend(pos));
                std::size_t n = 1;
                for (std::size_t r = 0; r < K; ++r)
                    n *= len(r);
                if constexpr (K == 1)
                    return pairwise<T>(axis_terms<S, Body>{std::move(b)}, 0, n, T{}, detail::mdarray::add_op{});
                // Positions neuves à chaque terme : les S premières, nulles, restent connues
                // du compilateur (un tableau partagé s'échapperait dans pairwise_split).
                counter wheel;
                auto const at = [&](std::size_t k) {
                    positions<S + K> q{};
                    if constexpr (K == 1)
                        q[S] = static_cast<std::ptrdiff_t>(k);
                    else if (flat)
                        q[S + K - 1] = static_cast<std::ptrdiff_t>(k);
                    else
                    {
                        advance(wheel, k);
                        for (std::size_t r = 0; r < K; ++r)
                            q[S + r] = wheel.p[r];
                    }
                    return b.value(q);
                };
                return pairwise<T>(at, 0, n, T{}, detail::mdarray::add_op{});
            }
            constexpr void shift(positions<S> const& pos) noexcept
            {
                body.shift(extend(pos));
            }
            template <class F>
            constexpr void for_each_access(F&& f)
            {
                body.for_each_access(f);
            }
            // Somme le long de K ≥ 2 axes : si chaque accès qui lit un axe sommé les lit tous
            // avec des pas cohérents (pas_r = Mᵣ × pas_dernier, Mᵣ le produit des longueurs
            // suivantes : g[p, i, j], t[i, j, k] en ligne majeure), la somme ne dépend que de
            // l'indice aplati t = Σ posᵣ Mᵣ : pas des autres axes mis à 0, t porté par le
            // dernier ; mêmes éléments, mêmes termes, même ordre. Refusé si le corps lit des
            // valeurs d'axes ou des δ, ou si une dimension contiguë (pas 1 fixé) s'y oppose.
            constexpr void flatten() noexcept
            {
                if constexpr (K >= 2 && !Body::reads_axes)
                {
                    std::array<std::size_t, K> weight{};
                    for (std::size_t r = K; r-- > 0;)
                        weight[r] = r + 1 < K ? weight[r + 1] * len(r + 1) : 1;
                    bool ok = true;
                    body.for_each_access([&](auto& a) {
                        using A = std::remove_cvref_t<decltype(a)>;
                        bool any = false;
                        for (std::size_t r = 0; r < K; ++r)
                            any = any || A::uses(S + r);
                        if (!any)
                            return;
                        if (!A::uses(S + K - 1))
                        {
                            ok = false;
                            return;
                        }
                        auto const stride = [&](std::size_t r) -> std::ptrdiff_t {
                            std::size_t const q = A::slot(S + r);
                            return q == A::P ? 0 : (A::unit[q] ? 1 : a.stride[q]);
                        };
                        std::ptrdiff_t const last = stride(K - 1);
                        for (std::size_t r = 0; r + 1 < K; ++r)
                            if (stride(r) != static_cast<std::ptrdiff_t>(weight[r]) * last || A::unit[A::slot(S + r)])
                                ok = false;
                    });
                    if (!ok)
                        return;
                    body.for_each_access([&](auto& a) {
                        using A = std::remove_cvref_t<decltype(a)>;
                        for (std::size_t r = 0; r + 1 < K; ++r)
                            if (A::slot(S + r) != A::P)
                                a.stride[A::slot(S + r)] = 0;
                    });
                    flat = true;
                }
            }

            // Somme le long d'un seul axe : sa longueur, et la somme d'une tranche [begin, end)
            // de ses termes (au plus 128, pairwise_leaves) pour W positions à la fois.
            static constexpr std::size_t summed_axes = K;
            constexpr std::size_t sum_length() const noexcept
            {
                std::size_t n = 1;
                for (std::size_t r = 0; r < K; ++r)
                    n *= len(r);
                return n;
            }
            // Longueur du r-ième axe sommé.
            constexpr std::size_t axis_length(std::size_t r) const noexcept
            {
                return len(r);
            }
            template <std::size_t I, std::size_t W>
            constexpr auto lanes_partial(positions<S> const& pos, std::size_t begin, std::size_t end) const
            {
                using V = decltype(body.value(positions<S + K>{}));
                using U = decltype(std::declval<V>() + std::declval<V>());
                using T = lane_pack<U, W>;
                Body b = body;
                b.shift(extend(pos));
                counter wheel;
                auto const at = [&](std::size_t k) {
                    positions<S + K> q{};
                    if constexpr (K == 1)
                        q[S] = static_cast<std::ptrdiff_t>(k);
                    else if (flat)
                        q[S + K - 1] = static_cast<std::ptrdiff_t>(k);
                    else
                    {
                        advance(wheel, k);
                        for (std::size_t r = 0; r < K; ++r)
                            q[S + r] = wheel.p[r];
                    }
                    auto const x = b.template lanes<I, W>(q);
                    T t;
                    for (std::size_t l = 0; l < W; ++l)
                        t.v[l] = x.v[l];
                    return t;
                };
                return pairwise<T>(at, begin, end, T{}, lane_add<U, W>{});
            }

            // Mêmes termes, même ordre, pour W positions extérieures à la fois.
            template <std::size_t I, std::size_t W>
            [[gnu::always_inline]] constexpr auto lanes(positions<S> const& pos) const
            {
                using V = decltype(body.value(positions<S + K>{}));
                using T = lane_pack<decltype(std::declval<V>() + std::declval<V>()), W>;
                Body b = body;
                b.shift(extend(pos));
                std::size_t n = 1;
                for (std::size_t r = 0; r < K; ++r)
                    n *= len(r);
                counter wheel;
                auto const at = [&](std::size_t k) {
                    positions<S + K> q{};
                    if constexpr (K == 1)
                        q[S] = static_cast<std::ptrdiff_t>(k);
                    else if (flat)
                        q[S + K - 1] = static_cast<std::ptrdiff_t>(k);
                    else
                    {
                        advance(wheel, k);
                        for (std::size_t r = 0; r < K; ++r)
                            q[S + r] = wheel.p[r];
                    }
                    auto const x = b.template lanes<I, W>(q);
                    T t;
                    for (std::size_t l = 0; l < W; ++l)
                        t.v[l] = x.v[l];
                    return t;
                };
                return pairwise<T>(at, 0, n, T{}, lane_add<decltype(std::declval<V>() + std::declval<V>()), W>{});
            }

        private:
            constexpr std::size_t len(std::size_t r) const noexcept
            {
                return Static[r] != dynamic_extent ? Static[r] : static_cast<std::size_t>(length[r]);
            }

            // Positions le long des K axes sommés du k-ième terme : la somme par paires les
            // visite dans l'ordre au sein d'une tranche, d'où un incrément (sans division),
            // une division seulement aux sauts.
            struct counter
            {
                std::array<std::ptrdiff_t, K> p{};
                std::size_t last = static_cast<std::size_t>(-1);
            };
            constexpr void advance(counter& c, std::size_t k) const noexcept
            {
                if (k == c.last + 1 && c.last != static_cast<std::size_t>(-1))
                {
                    for (std::size_t r = K; r-- > 0;)
                    {
                        if (++c.p[r] < static_cast<std::ptrdiff_t>(len(r)))
                            break;
                        c.p[r] = 0;
                    }
                }
                else
                {
                    std::size_t q = k;
                    for (std::size_t r = K; r-- > 0;)
                    {
                        c.p[r] = static_cast<std::ptrdiff_t>(q % len(r));
                        q /= len(r);
                    }
                }
                c.last = k;
            }

            static constexpr positions<S + K> extend(positions<S> const& pos) noexcept
            {
                positions<S + K> p{};
                for (std::size_t r = 0; r < S; ++r)
                    p[r] = pos[r];
                return p;
            }
        };

        // Un nœud, et la somme le long de ses propres axes s'il en est le foyer.
        template <class Ctx, class Summed, class Node, class Env>
        constexpr auto make_strided(Node const& node, Env const& env)
        {
            using own = typename home<Node, Ctx>::type;
            if constexpr (own::size == 0)
                return strided_body<Ctx, Summed>(node, env);
            else
            {
                constexpr std::size_t S = Summed::size, K = own::size;
                auto body = strided_body<Ctx, tconcat_t<Summed, own>>(node, env);
                check_not_bool_sum<decltype(body.value(positions<S + K>{}))>();
                constexpr auto fixed = []<class... Ids>(tlist<Ids...>) {
                    return std::array<std::size_t, K>{static_length_in<Ids, typename Node::static_lengths>::value...};
                }(own{});
                strided_home<S, K, decltype(body), fixed> h{body};
                for (std::size_t r = 0; r < K; ++r)
                    h.length[r] = env.length[axis_indices<typename Ctx::axes>(own{})[r]];
                h.flatten();
                return h;
            }
        }

        // Toutes les lectures le long de l'axe Id sont-elles contiguës (pas unité) ? Sinon,
        // Clang vectorise la boucle de l'instruction par des lectures dispersées (gather
        // d'AVX2), bien plus lentes ici que la boucle scalaire (gₚᵢⱼ vₚᵢ vₚⱼ : 19 ns par point
        // contre 3,4 sur un i7-9750H ; la même boucle écrite à la main aussi) ; GCC ne le
        // fait pas.
        template <class Id, class Node>
        struct contiguous_along : std::true_type
        {
        };
        template <class Id, class V, class... S>
        struct contiguous_along<Id, indexed_leaf<view_source<V>, S...>>
            : std::bool_constant<!(std::is_same_v<S, axis_slot<Id>> || ...) || unit_stride<V, Id, S...>()>
        {
        };
        template <class Id, class V, class... S>
        struct contiguous_along<Id, indexed_target<V, S...>> : contiguous_along<Id, indexed_leaf<view_source<V>, S...>>
        {
        };
        template <class Id, class Op, class... C>
        struct contiguous_along<Id, indexed_node<Op, C...>> : std::bool_constant<(contiguous_along<Id, C>::value && ...)>
        {
        };

        template <class L>
        struct last_of;
        template <class... T>
        struct last_of<tlist<T...>>
        {
            using type = typename decltype((std::type_identity<T>{}, ...))::type;
        };

        template <class Ctx, class Summed, bool Root, class Node, class Env, std::size_t K>
        constexpr auto terms(Node const& node, Env& env, odometer<K>& wheels)
        {
            if constexpr (strideable<Ctx, Node>::value)
            {
                auto term = [&] {
                    if constexpr (Root)
                        return strided_body<Ctx, Summed>(node, env);
                    else
                        return make_strided<Ctx, Summed>(node, env);
                }();
                if constexpr (K == 1)
                    return axis_terms<0, decltype(term)>{std::move(term)};
                else
                    return [term, &wheels](std::size_t k) {
                        wheels.set_local(k);
                        return term.value(wheels.local);
                    };
            }
            else if constexpr (Root)
                return [&](std::size_t k) {
                    wheels.set(k, env);
                    return node.template compute<Ctx>(env);
                };
            else
                return [&](std::size_t k) {
                    wheels.set(k, env);
                    return evaluate<Ctx>(node, env);
                };
        }

        template <class Ctx, class Node, class Env>
        constexpr auto evaluate(Node const& node, Env& env)
        {
            using own = typename home<Node, Ctx>::type;
            if constexpr (own::size == 0)
                return node.template compute<Ctx>(env);
            else
            {
                check_not_bool_sum<decltype(node.template compute<Ctx>(env))>();
                constexpr std::size_t K = own::size;
                odometer<K> wheels(axis_indices<typename Ctx::axes>(own{}), env);
                std::size_t const n = wheels.size();
                auto const at = terms<Ctx, own, true>(node, env, wheels);
                using T = decltype(at(0) + at(0));
                return pairwise<T>(at, 0, n, T{}, detail::mdarray::add_op{});
            }
        }

        // --- Vérifications d'une instruction (V4, L5) --------------------------------------

        template <class Id, class Root>
        struct appears_twice
        {
            static constexpr bool value = tcount<Id, Root> >= 2;
        };

        // Garde-fou (2) : un axe sommé n'a pas à la fois ≥ 2 occurrences dans un argument
        // de fonction et d'autres ailleurs.
        template <class Node, class Root, class Summed>
        struct scope_ok;
        template <class Node, class Root, class... S>
        struct scope_ok<Node, Root, tlist<S...>>
        {
            template <class Arg>
            struct arg_ok
            {
                static constexpr bool value =
                    ((!(tcount<S, typename Arg::occurrences> >= 2 && tcount<S, Root> > tcount<S, typename Arg::occurrences>)) && ...);
            };
            template <class C>
            struct child_ok
            {
                static constexpr bool value = scope_ok<C, Root, tlist<S...>>::value;
            };
            static constexpr bool value = (!Node::is_function || tall<typename Node::children, arg_ok>)
                                          && tall<typename Node::children, child_ok>;
        };

        template <class A, class B>
        struct length_pair_ok
        {
            static constexpr bool value = !std::is_same_v<typename A::id, typename B::id> || A::value == B::value;
        };
        template <class L>
        struct static_lengths_ok;
        template <>
        struct static_lengths_ok<tlist<>>
        {
            static constexpr bool value = true;
        };
        template <class X, class... Rest>
        struct static_lengths_ok<tlist<X, Rest...>>
        {
            static constexpr bool value = (length_pair_ok<X, Rest>::value && ...) && static_lengths_ok<tlist<Rest...>>::value;
        };

        template <class Lhs, class LhsProvided, class LhsStatic, class Node>
        struct statement
        {
            using root = typename Node::occurrences;
            using bound = tunique_t<typename Node::bound>;
            using summed = tminus_t<tunique_t<root>, Lhs>;
            using axes = tunique_t<tconcat_t<Lhs, root, bound>>;
            using ctx = context<axes, summed, root>;
            using provided = tunique_t<tconcat_t<LhsProvided, typename Node::provided>>;

            template <class Id>
            struct is_provided
            {
                static constexpr bool value = tcontains<Id, provided>;
            };
            template <class Id>
            struct not_reused
            {
                static constexpr bool value = !tcontains<Id, Lhs> && !tcontains<Id, root>;
            };
            template <class Id>
            struct twice
            {
                static constexpr bool value = appears_twice<Id, root>::value;
            };

            static_assert(tall<bound, not_reused>, "wsl::mdarray : un axe lié par une réduction réapparaît à gauche ou hors "
                                                   "de la réduction (V4, garde-fou 3)");
            static_assert(!tall<bound, not_reused> || tall<summed, twice>, "wsl::mdarray : un axe absent à gauche n'apparaît qu'une fois à droite (V4, "
                                               "garde-fou 1) ; l'écrire à gauche, ou réduire explicitement : sum(j)(…)");
            static_assert(scope_ok<Node, root, summed>::value,
                          "wsl::mdarray : axe sommé à la fois deux fois dans un argument de fonction et ailleurs (V4, "
                          "garde-fou 2) : portée ambiguë ; renommer l'un des axes");
            static_assert(tall<axes, is_provided>, "wsl::mdarray : un axe n'étiquette aucune dimension d'un tableau : "
                                                   "longueur inconnue (L5)");
            static_assert(static_lengths_ok<tconcat_t<LhsStatic, typename Node::static_lengths>>::value,
                          "wsl::mdarray : un axe étiquette des dimensions fixes de longueurs différentes (L5)");
        };

        template <class Node>
        constexpr auto evaluate_scalar(Node const& node)
        {
            using st = statement<tlist<>, tlist<>, tlist<>, Node>;
            environment<typename st::axes> env;
            node.lengths(env);
            return evaluate<typename st::ctx>(node, env);
        }

        // --- Normalisation des opérandes ---------------------------------------------------

        template <class T>
        constexpr auto to_indexed(T&& t)
        {
            using U = std::remove_cvref_t<T>;
            if constexpr (indexed<U>)
            {
                static_assert(!std::is_lvalue_reference_v<T>,
                              "wsl::mdarray : une expression ne se garde pas (L2) ; écrire le calcul dans l'instruction, "
                              "ou std::move(e)");
                return U(std::move(t));
            }
            else if constexpr (axis_like<U>)
                return indexed_axis<typename U::id_type, std::ptrdiff_t>{};
            else
            {
                static_assert(!array_operand<U>, "wsl::mdarray : un tableau entier ne se mêle pas à une expression indicée ; "
                                                 "l'indexer : a[i, j]");
                return indexed_scalar<U>(std::forward<T>(t));
            }
        }

        template <class T>
        using indexed_t = decltype(to_indexed(std::declval<T>()));

        template <class T>
        concept indexed_operand = indexed<T> || axis_like<T>;

        template <class... T>
        concept with_indexed = (indexed_operand<T> || ...)
                               && ((indexed_operand<T> || scalar_operand<T> || array_operand<T>) && ...);

        template <class Op, class... T>
        constexpr auto make_node(Op op, T&&... t)
        {
            return indexed_node<Op, indexed_t<T>...>(std::move(op), to_indexed(std::forward<T>(t))...);
        }

        template <class Op, class... T>
        concept node_valid = with_indexed<T...>;

        // --- Cible d'une affectation : m[i, j] = … (V9) ---------------------------------------

        template <class View, class... Slots>
        struct indexed_target : indexed_leaf<view_source<View>, Slots...>
        {
            using leaf = indexed_leaf<view_source<View>, Slots...>;
            using leaf::leaf;
            indexed_target(indexed_target const&) = delete;
            constexpr indexed_target(indexed_target&&) = default;

            using lhs_axes = tunique_t<typename leaf::occurrences>;

            template <class Rhs>
            constexpr void assign(Rhs const& rhs) const
            {
                using st = statement<lhs_axes, lhs_axes, typename leaf::static_lengths, Rhs>;
                using ctx = typename st::ctx;
                using axes = typename st::axes;
                environment<axes> env;
                this->lengths(env);
                rhs.lengths(env);
                constexpr std::size_t K = lhs_axes::size;
                odometer<K> wheels(axis_indices<axes>(lhs_axes{}), env);
                std::size_t const n = wheels.size();
                auto place = [&](std::size_t k) { wheels.set(k, env); };
                if constexpr (strideable<ctx, Rhs>::value && strideable<ctx, leaf>::value)
                {
                    // Accès à pas fixe : le plan est construit une fois pour l'instruction.
                    auto const value = make_strided<ctx, lhs_axes>(rhs, env);
                    auto const target = make_access<lhs_axes>(static_cast<leaf const&>(*this), env);
                    using inner = typename last_of<lhs_axes>::type;
                    constexpr bool contiguous = contiguous_along<inner, Rhs>::value && contiguous_along<inner, leaf>::value;
                    auto loop = [&](auto&& step) {
                        if constexpr (contiguous)
                        {
                            WSL_MDARRAY_INDEPENDENT_LOOP
                            for (std::size_t k = 0; k < n; ++k)
                                step(k);
                        }
                        else
                        {
                            WSL_MDARRAY_SCALAR_LOOP
                            for (std::size_t k = 0; k < n; ++k)
                                step(k);
                        }
                    };
                    if (aliased<lhs_axes>(rhs))
                    {
                        using T = typename View::value_type;
                        std::vector<T> temporary(n);
                        loop([&](std::size_t k) {
                            wheels.set_local(k);
                            temporary[k] = static_cast<T>(value.value(wheels.local));
                        });
                        for (std::size_t k = 0; k < n; ++k)
                        {
                            wheels.set_local(k);
                            target(wheels.local) = temporary[k];
                        }
                    }
                    else if constexpr (contiguous && ctx::summed::size > 0)
                    {
                        // Une somme : W positions consécutives le long de l'axe intérieur à la
                        // fois (voies), le reste une à une ; mêmes bits.
                        constexpr std::size_t W = 8;
                        auto const inner_length = static_cast<std::size_t>(env.length[tindex_v<inner, axes>]);
                        // Paquet de colonnes à l'extérieur : les lectures d'un paquet (B[j, k…k+7]
                        // pour tout j) restent dans le cache d'une ligne à l'autre. Chaque
                        // élément est indépendant : l'ordre de parcours ne change pas les bits.
                        std::size_t const full = inner_length - inner_length % W;
                        if constexpr (requires { value.template lanes_partial<K - 1, W>(wheels.local, 0, 0); })
                            if (value.sum_length() > 128)
                            {
                                tiled<W>(value, target, wheels, n, inner_length, full);
                                for (std::size_t outer = 0; outer < n; outer += inner_length)
                                    for (std::size_t r = full; r < inner_length; ++r)
                                    {
                                        wheels.set_local(outer + r);
                                        target(wheels.local) = static_cast<typename View::value_type>(value.value(wheels.local));
                                    }
                                return;
                            }
                        for (std::size_t r = 0; r < full; r += W)
                            for (std::size_t outer = 0; outer < n; outer += inner_length)
                            {
                                wheels.set_local(outer + r);
                                auto const pack = value.template lanes<K - 1, W>(wheels.local);
                                std::ptrdiff_t const o = target.index(wheels.local), step = target.template lane_step<K - 1>();
                                for (std::size_t l = 0; l < W; ++l)
                                    target.p[o + static_cast<std::ptrdiff_t>(l) * step] =
                                        static_cast<typename View::value_type>(pack.v[l]);
                            }
                        for (std::size_t outer = 0; outer < n; outer += inner_length)
                            for (std::size_t r = full; r < inner_length; ++r)
                            {
                                wheels.set_local(outer + r);
                                target(wheels.local) = static_cast<typename View::value_type>(value.value(wheels.local));
                            }
                    }
                    else
                        loop([&](std::size_t k) {
                            wheels.set_local(k);
                            target(wheels.local) = static_cast<typename View::value_type>(value.value(wheels.local));
                        });
                    return;
                }
                auto& view = this->source.view;
                auto write = [&](auto const& value) {
                    std::apply([&](auto const&... s) { view[static_cast<typename View::index_type>(slot_value(s, env))...] = value; },
                               this->slots);
                };
                if (aliased<lhs_axes>(rhs))
                {
                    using T = typename View::value_type;
                    std::vector<T> temporary;
                    temporary.reserve(n);
                    for (std::size_t k = 0; k < n; ++k)
                    {
                        place(k);
                        temporary.push_back(static_cast<T>(evaluate<ctx>(rhs, env)));
                    }
                    for (std::size_t k = 0; k < n; ++k)
                    {
                        place(k);
                        write(temporary[k]);
                    }
                }
                else
                    for (std::size_t k = 0; k < n; ++k)
                    {
                        place(k);
                        write(static_cast<typename View::value_type>(evaluate<ctx>(rhs, env)));
                    }
            }

            template <class Rhs>
                requires(indexed_operand<Rhs> || scalar_operand<Rhs> || array_operand<Rhs>)
            constexpr void operator=(Rhs&& rhs) &&
            {
                assign(default_mode(to_indexed(std::forward<Rhs>(rhs))));
            }

            constexpr void operator=(indexed_target&& rhs) &&
            {
                assign(std::move(rhs));
            }

            template <class Rhs>
            void operator=(Rhs&&) &
            {
                static_assert(sizeof(Rhs) == 0, "wsl::mdarray : une cible m[i, j] ne se garde pas (L2) ; écrire m[i, j] = … "
                                                "directement");
            }

#define WSL_MDARRAY_TARGET_COMPOUND(op, base)                                                                          \
    template <class Rhs>                                                                                               \
        requires(indexed_operand<Rhs> || scalar_operand<Rhs>)                                                          \
    constexpr void operator op(Rhs&& rhs) &&                                                                           \
    {                                                                                                                  \
        assign(default_mode(make_node(base{}, reread(), std::forward<Rhs>(rhs))));                                     \
    }
            WSL_MDARRAY_TARGET_COMPOUND(+=, detail::mdarray::add_op)
            WSL_MDARRAY_TARGET_COMPOUND(-=, detail::mdarray::subtract_op)
            WSL_MDARRAY_TARGET_COMPOUND(*=, detail::mdarray::multiply_op)
            WSL_MDARRAY_TARGET_COMPOUND(/=, std::divides<>)
#undef WSL_MDARRAY_TARGET_COMPOUND

        private:
            // Somme le long d'un axe, plus de 128 termes : découpage du cache le long de
            // l'axe sommé. Pour un paquet de colonnes et un bloc de Rows positions
            // extérieures, chaque tranche de la somme par paires (au plus 128 termes : la
            // tranche de B tient dans le cache L1) est calculée pour tout le bloc, puis les
            // résultats partiels de chaque élément sont combinés selon l'arbre de la somme
            // par paires : mêmes bits que la somme d'un trait.
            template <std::size_t W, class Value, class Target, class Wheels>
            static void tiled(Value const& value, Target const& target, Wheels& wheels, std::size_t n, std::size_t inner_length,
                              std::size_t full)
            {
                constexpr std::size_t K = lhs_axes::size, Rows = 64;
                using pack = decltype(value.template lanes_partial<K - 1, W>(wheels.local, 0, 0));
                std::size_t const length = value.sum_length();
                std::vector<std::size_t> bounds;
                pairwise_leaves(0, length, [&](std::size_t b, std::size_t) { bounds.push_back(b); });
                bounds.push_back(length);
                std::size_t const chunks = bounds.size() - 1;
                std::size_t const outer_count = n / inner_length;
                std::vector<pack> partial(Rows * chunks);
                std::array<positions<K>, Rows> where{};

                // Recopie : une feuille du corps qui ne dépend que des axes sommés (positions
                // K à K + Ks − 1, tous) et de l'axe intérieur (K − 1), B[j, k], a les mêmes
                // valeurs pour tout le bloc ; ses valeurs d'une tranche et d'un paquet de
                // colonnes sont recopiées dans un tampon contigu (les lignes de B espacées
                // d'une puissance de 2 tomberaient dans le même ensemble du cache), et la
                // copie du plan lit ce tampon : mêmes valeurs, mêmes opérations, mêmes bits.
                // L'indice du tampon, (t − b) × largeur + voie, t = Σ posᵣ Mᵣ l'indice aplati
                // des termes, est linéaire : il suffit de changer décalage et pas. Une
                // dimension contiguë a un pas 1 fixé à la compilation : si le tampon en
                // demande un autre, la feuille reste lue en place.
                constexpr std::size_t Ks = Value::summed_axes;
                std::array<std::size_t, Ks> weight{};
                for (std::size_t r = Ks; r-- > 0;)
                    weight[r] = r + 1 < Ks ? weight[r + 1] * value.axis_length(r + 1) : 1;
                Value plan = value;
                std::vector<std::function<void(std::size_t, std::size_t, std::size_t)>> fills;
                plan.for_each_access([&](auto& a) {
                    using A = std::remove_cvref_t<decltype(a)>;
                    using T = std::remove_cv_t<std::remove_pointer_t<decltype(a.p)>>;
                    if constexpr (A::positions_size == K + Ks && A::only(K - 1, K, K + Ks) && std::is_trivially_copyable_v<T>)
                    {
                        constexpr bool wide = A::uses(K - 1);
                        constexpr std::size_t width = wide ? W : 1;
                        for (std::size_t r = 0; r < Ks; ++r)
                            if (A::unit[A::slot(K + r)] && weight[r] * width != 1)
                                return;
                        auto buffer = std::make_shared<std::vector<T>>(128 * width);
                        A const original = a;
                        fills.push_back([&a, original, buffer, weight, &value](std::size_t b, std::size_t e, std::size_t r) {
                            // Constantes redites par des types (GCC refuse ici les constantes
                            // locales non capturées).
                            constexpr std::size_t K = lhs_axes::size, Ks = Value::summed_axes;
                            constexpr bool wide = A::uses(K - 1);
                            constexpr std::size_t width = wide ? W : 1;
                            positions<K + Ks> q{};
                            for (std::size_t t = b; t < e; ++t)
                            {
                                if (value.flat)
                                    q[K + Ks - 1] = static_cast<std::ptrdiff_t>(t);
                                else
                                {
                                    std::size_t rest = t;
                                    for (std::size_t m = Ks; m-- > 0;)
                                    {
                                        q[K + m] = static_cast<std::ptrdiff_t>(rest % value.axis_length(m));
                                        rest /= value.axis_length(m);
                                    }
                                }
                                for (std::size_t l = 0; l < width; ++l)
                                {
                                    q[K - 1] = static_cast<std::ptrdiff_t>(r + l);
                                    (*buffer)[(t - b) * width + l] = original(q);
                                }
                            }
                            a.p = buffer->data();
                            a.offset = -static_cast<std::ptrdiff_t>(b * width) - (wide ? static_cast<std::ptrdiff_t>(r) : 0);
                            for (std::size_t m = 0; m < Ks; ++m)
                                a.stride[A::slot(K + m)] =
                                    value.flat ? (m + 1 == Ks ? static_cast<std::ptrdiff_t>(width) : 0)
                                               : static_cast<std::ptrdiff_t>(weight[m] * width);
                        });
                    }
                });

                for (std::size_t r = 0; r < full; r += W)
                    for (std::size_t o0 = 0; o0 < outer_count; o0 += Rows)
                    {
                        std::size_t const rows = std::min(Rows, outer_count - o0);
                        for (std::size_t o = 0; o < rows; ++o)
                        {
                            wheels.set_local((o0 + o) * inner_length + r);
                            where[o] = wheels.local;
                        }
                        for (std::size_t c = 0; c < chunks; ++c)
                        {
                            for (auto const& fill : fills)
                                fill(bounds[c], bounds[c + 1], r);
                            for (std::size_t o = 0; o < rows; ++o)
                                partial[o * chunks + c] = plan.template lanes_partial<K - 1, W>(where[o], bounds[c], bounds[c + 1]);
                        }
                        for (std::size_t o = 0; o < rows; ++o)
                        {
                            std::size_t next = 0;
                            auto const sum = pairwise_combine(0, length, partial.data() + o * chunks, next,
                                                              lane_add<std::remove_cvref_t<decltype(pack{}.v[0])>, W>{});
                            std::ptrdiff_t const t = target.index(where[o]), step = target.template lane_step<K - 1>();
                            for (std::size_t l = 0; l < W; ++l)
                                target.p[t + static_cast<std::ptrdiff_t>(l) * step] = static_cast<typename View::value_type>(sum.v[l]);
                        }
                    }
            }

            // La cible relue, aux mêmes positions (affectations composées, V9).
            constexpr leaf reread() const
            {
                return std::apply([&](auto const&... s) { return leaf(this->source, s...); }, this->slots);
            }

            // Une feuille lit-elle notre stockage ailleurs qu'aux mêmes positions ? (I12)
            template <class Lhs, class Rhs>
            constexpr bool aliased(Rhs const& rhs) const
            {
                if consteval
                {
                    return true;
                }
                else
                {
                    auto const& view = this->source.view;
                    void const* lo = view.data_handle();
                    void const* hi = view.data_handle() + view.mapping().required_span_size();
                    std::less<void const*> less;
                    bool result = false;
                    auto overlaps = [&](auto const& v) {
                        void const* p = v.data_handle();
                        void const* q = v.data_handle() + v.mapping().required_span_size();
                        return less(p, hi) && less(lo, q);
                    };
                    rhs.for_each_leaf([&](auto const& l) {
                        using L = std::remove_cvref_t<decltype(l)>;
                        if constexpr (requires { l.source.view; })
                        {
                            auto const& v = l.source.view;
                            if (!overlaps(v))
                                return;
                            if constexpr (std::is_same_v<L, leaf> && tminus_t<typename L::occurrences, Lhs>::size == 0
                                          && std::is_same_v<std::remove_cvref_t<decltype(v)>, View>)
                                if (v.data_handle() == view.data_handle() && v.mapping() == view.mapping()
                                    && same_fixed(l.slots, this->slots))
                                    return; // mêmes positions, hors somme : sûr
                            result = true;
                        }
                        else
                            l.source.expr.for_each_storage([&](auto const& v) {
                                if (overlaps(v))
                                    result = true;
                            });
                    });
                    return result;
                }
            }

            static constexpr bool same_fixed(std::tuple<Slots...> const& a, std::tuple<Slots...> const& b)
            {
                return [&]<std::size_t... R>(std::index_sequence<R...>) {
                    return (same_slot(std::get<R>(a), std::get<R>(b)) && ...);
                }(std::index_sequence_for<Slots...>{});
            }
            static constexpr bool same_slot(fixed_slot const& a, fixed_slot const& b) { return a.value == b.value; }
            template <class Id>
            static constexpr bool same_slot(axis_slot<Id> const&, axis_slot<Id> const&)
            {
                return true;
            }
        };

        template <class A>
        inline constexpr bool slot_argument = axis_like<A> || index_integer<std::remove_cvref_t<A>>;
    }

    // --- m[i, j] : make_indexed (déclarée dans expression_base.hpp) -------------------------

    template <class Source, class... Args>
    constexpr auto detail::mdarray::make_indexed(Source&& source, Args const&... args)
    {
        using S = std::remove_cvref_t<Source>;
        static_assert(sizeof...(Args) == S::rank(), "wsl::mdarray : nombre d'indices différent du rang");
        if constexpr (is_mdspan<S>)
        {
            if constexpr (std::is_const_v<typename S::element_type>)
                return indexed_leaf<view_source<S>, decltype(make_slot(args))...>(view_source<S>{source}, make_slot(args)...);
            else
                return indexed_target<S, decltype(make_slot(args))...>(view_source<S>{source}, make_slot(args)...);
        }
        else
        {
            static_assert(!std::is_lvalue_reference_v<Source>, "wsl::mdarray : une expression ne se garde pas (L2)");
            return indexed_leaf<expression_source<S>, decltype(make_slot(args))...>(expression_source<S>{std::move(source)},
                                                                                    make_slot(args)...);
        }
    }

    // --- Opérateurs du monde des axes (V3) -----------------------------------------------

#define WSL_MDARRAY_INDEXED_BINARY(op, functor)                                                                        \
    template <class A, class B>                                                                                        \
        requires detail::mdarray::node_valid<functor, A, B>                                                            \
    constexpr auto operator op(A&& a, B&& b)                                                                           \
    {                                                                                                                  \
        return detail::mdarray::make_node(functor{}, std::forward<A>(a), std::forward<B>(b));                          \
    }

    WSL_MDARRAY_INDEXED_BINARY(+, detail::mdarray::add_op)
    WSL_MDARRAY_INDEXED_BINARY(-, detail::mdarray::subtract_op)
    WSL_MDARRAY_INDEXED_BINARY(*, detail::mdarray::multiply_op)
    WSL_MDARRAY_INDEXED_BINARY(/, std::divides<>)
    WSL_MDARRAY_INDEXED_BINARY(%, std::modulus<>)
    WSL_MDARRAY_INDEXED_BINARY(==, std::equal_to<>)
    WSL_MDARRAY_INDEXED_BINARY(!=, std::not_equal_to<>)
    WSL_MDARRAY_INDEXED_BINARY(<, std::less<>)
    WSL_MDARRAY_INDEXED_BINARY(<=, std::less_equal<>)
    WSL_MDARRAY_INDEXED_BINARY(>, std::greater<>)
    WSL_MDARRAY_INDEXED_BINARY(>=, std::greater_equal<>)
    WSL_MDARRAY_INDEXED_BINARY(&&, detail::mdarray::logical_and_op)
    WSL_MDARRAY_INDEXED_BINARY(||, detail::mdarray::logical_or_op)
#undef WSL_MDARRAY_INDEXED_BINARY

    template <class A>
        requires detail::mdarray::with_indexed<A>
    constexpr auto operator-(A&& a)
    {
        return detail::mdarray::make_node(detail::mdarray::negate_op{}, std::forward<A>(a));
    }

    template <class A>
        requires detail::mdarray::with_indexed<A>
    constexpr auto operator!(A&& a)
    {
        return detail::mdarray::make_node(detail::mdarray::logical_not_op{}, std::forward<A>(a));
    }

    // +i : réservé (marques de variance des futurs tenseurs, L12).
    template <class Id>
    constexpr void operator+(mdaxis<Id> const&)
    {
        static_assert(sizeof(Id) == 0, "wsl::mdarray : +i est réservé (L12) ; un axe comme valeur s'écrit i");
    }

#define WSL_MDARRAY_INDEXED_UNARY_FUNCTION(name)                                                                       \
    template <class A>                                                                                                 \
        requires detail::mdarray::with_indexed<A>                                                                      \
    constexpr auto name(A&& a)                                                                                         \
    {                                                                                                                  \
        return detail::mdarray::make_node(detail::mdarray::name##_op{}, std::forward<A>(a));                           \
    }
#define WSL_MDARRAY_INDEXED_BINARY_FUNCTION(name)                                                                      \
    template <class A, class B>                                                                                        \
        requires detail::mdarray::with_indexed<A, B>                                                                   \
    constexpr auto name(A&& a, B&& b)                                                                                  \
    {                                                                                                                  \
        return detail::mdarray::make_node(detail::mdarray::name##_op{}, std::forward<A>(a), std::forward<B>(b));       \
    }
    WSL_MDARRAY_UNARY_FUNCTIONS(WSL_MDARRAY_INDEXED_UNARY_FUNCTION)
    WSL_MDARRAY_BINARY_FUNCTIONS(WSL_MDARRAY_INDEXED_BINARY_FUNCTION)
#undef WSL_MDARRAY_INDEXED_UNARY_FUNCTION
#undef WSL_MDARRAY_INDEXED_BINARY_FUNCTION

    template <class A, class B>
        requires detail::mdarray::with_indexed<A, B>
    constexpr auto min(A&& a, B&& b)
    {
        return detail::mdarray::make_node(detail::mdarray::min_op{}, std::forward<A>(a), std::forward<B>(b));
    }

    template <class A, class B>
        requires detail::mdarray::with_indexed<A, B>
    constexpr auto max(A&& a, B&& b)
    {
        return detail::mdarray::make_node(detail::mdarray::max_op{}, std::forward<A>(a), std::forward<B>(b));
    }

    template <class C, class A, class B>
        requires detail::mdarray::with_indexed<C, A, B>
    constexpr auto where(C&& c, A&& a, B&& b)
    {
        return detail::mdarray::make_node(detail::mdarray::where_op{}, std::forward<C>(c), std::forward<A>(a), std::forward<B>(b));
    }

    // Fonction de l'utilisateur dans le monde des axes (I11) : F(x[i], x[j]).
    template <class F>
    template <class... T>
    constexpr auto vectorized<F>::make_indexed_call(T&&... t) const
    {
        return detail::mdarray::make_node(detail::mdarray::call_op<F>{f_}, std::forward<T>(t)...);
    }

    // Axe comme valeur, dans un type choisi (L12) ; δᵢⱼ (L13).
    template <class T = std::ptrdiff_t, class Id>
    constexpr auto position(mdaxis<Id> const&)
    {
        return detail::mdarray::indexed_axis<Id, T>{};
    }

    template <class T = int, class A, class B>
    constexpr auto kronecker(mdaxis<A> const&, mdaxis<B> const&)
    {
        return detail::mdarray::indexed_kronecker<A, B, T>{};
    }

    // Une expression indicée, avec un membre de gauche sans axe (V9).
    template <class E>
        requires detail::mdarray::indexed<E>
    constexpr auto value(E&& e)
    {
        static_assert(!std::is_lvalue_reference_v<E>, "wsl::mdarray : une expression ne se garde pas (L2)");
        return detail::mdarray::evaluate_scalar(detail::mdarray::default_mode(std::move(e)));
    }

    // --- Réducteurs le long d'axes : sum(j)(…) (I6) ------------------------------------------

    template <class Kind, class... Ids>
    class reducer
    {
    public:
        template <class E>
        constexpr auto operator()(E&& e) const
        {
            auto child = detail::mdarray::to_indexed(std::forward<E>(e));
            using node = detail::mdarray::indexed_reduction<Kind, detail::mdarray::tlist<Ids...>, decltype(child)>;
            node r(std::move(child));
            if constexpr (node::occurrences::size == 0)
                return detail::mdarray::evaluate_scalar(detail::mdarray::default_mode(std::move(r))); // entièrement liée (V9)
            else
                return r;
        }
    };

#define WSL_MDARRAY_REDUCER(name, kind)                                                                                \
    template <class... Ids>                                                                                            \
    constexpr reducer<detail::mdarray::kind, Ids...> name(mdaxis<Ids> const&...) noexcept                             \
        requires(sizeof...(Ids) > 0)                                                                                   \
    {                                                                                                                  \
        return {};                                                                                                     \
    }
    WSL_MDARRAY_REDUCER(sum, sum_kind)
    WSL_MDARRAY_REDUCER(prod, prod_kind)
    WSL_MDARRAY_REDUCER(all, all_kind)
    WSL_MDARRAY_REDUCER(any, any_kind)
    WSL_MDARRAY_REDUCER(count, count_kind)
#undef WSL_MDARRAY_REDUCER

    // min, max, argmin, argmax : un seul axe (deux axes : max(i, j) est élément par élément).
    template <class Id>
    constexpr reducer<detail::mdarray::min_kind, Id> min(mdaxis<Id> const&) noexcept
    {
        return {};
    }
    template <class Id>
    constexpr reducer<detail::mdarray::max_kind, Id> max(mdaxis<Id> const&) noexcept
    {
        return {};
    }
    template <class Id>
    constexpr reducer<detail::mdarray::argmin_kind, Id> argmin(mdaxis<Id> const&) noexcept
    {
        return {};
    }
    template <class Id>
    constexpr reducer<detail::mdarray::argmax_kind, Id> argmax(mdaxis<Id> const&) noexcept
    {
        return {};
    }
}

#if defined(__clang__)
#    pragma clang diagnostic pop
#endif
