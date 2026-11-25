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
#include <wsl/mdarray/detail/precondition.hpp>
#include <wsl/mdarray/detail/reduce.hpp>
#include <wsl/mdarray/detail/typelist.hpp>
#include <wsl/mdarray/expression.hpp>
#include <wsl/mdarray/mdspan.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

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

        // --- Emplacements d'un crochet ------------------------------------------------

        template <class Id>
        struct axis_slot
        {
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
                auto const indices = bound_indices<typename Ctx::axes>(Bound{});
                std::size_t n = 1;
                for (std::size_t r = 0; r < K; ++r)
                    n *= static_cast<std::size_t>(env.length[indices[r]]);
                auto at = [&](std::size_t k) {
                    for (std::size_t r = K; r-- > 0;)
                    {
                        auto const len = static_cast<std::size_t>(env.length[indices[r]]);
                        env.position[indices[r]] = static_cast<std::ptrdiff_t>(k % len);
                        k /= len;
                    }
                    return evaluate<Ctx>(child, env);
                };
                if constexpr (std::is_same_v<Kind, sum_kind>)
                {
                    using T = decltype(at(0) + at(0));
                    return pairwise<T>(at, 0, n, T{}, std::plus<>{});
                }
                else if constexpr (std::is_same_v<Kind, prod_kind>)
                {
                    using T = decltype(at(0) * at(0));
                    return pairwise<T>(at, 0, n, T{1}, std::multiplies<>{});
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

        template <class Ctx, class Node, class Env>
        constexpr auto evaluate(Node const& node, Env& env)
        {
            using own = typename home<Node, Ctx>::type;
            if constexpr (own::size == 0)
                return node.template compute<Ctx>(env);
            else
            {
                static_assert(!std::is_same_v<std::remove_cvref_t<decltype(node.template compute<Ctx>(env))>, bool>,
                              "wsl::mdarray : somme de booléens (I10) ; écrire count(i)(…), any(i)(…) ou all(i)(…)");
                constexpr std::size_t K = own::size;
                auto const indices = axis_indices<typename Ctx::axes>(own{});
                std::size_t n = 1;
                for (std::size_t r = 0; r < K; ++r)
                    n *= static_cast<std::size_t>(env.length[indices[r]]);
                auto at = [&](std::size_t k) {
                    for (std::size_t r = K; r-- > 0;)
                    {
                        auto const len = static_cast<std::size_t>(env.length[indices[r]]);
                        env.position[indices[r]] = static_cast<std::ptrdiff_t>(k % len);
                        k /= len;
                    }
                    return node.template compute<Ctx>(env);
                };
                using T = decltype(at(0) + at(0));
                return pairwise<T>(at, 0, n, T{}, std::plus<>{});
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
                auto const indices = axis_indices<axes>(lhs_axes{});
                std::size_t n = 1;
                for (std::size_t r = 0; r < K; ++r)
                    n *= static_cast<std::size_t>(env.length[indices[r]]);
                auto place = [&](std::size_t k) {
                    for (std::size_t r = K; r-- > 0;)
                    {
                        auto const len = static_cast<std::size_t>(env.length[indices[r]]);
                        env.position[indices[r]] = static_cast<std::ptrdiff_t>(k % len);
                        k /= len;
                    }
                };
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
            WSL_MDARRAY_TARGET_COMPOUND(+=, std::plus<>)
            WSL_MDARRAY_TARGET_COMPOUND(-=, std::minus<>)
            WSL_MDARRAY_TARGET_COMPOUND(*=, std::multiplies<>)
            WSL_MDARRAY_TARGET_COMPOUND(/=, std::divides<>)
#undef WSL_MDARRAY_TARGET_COMPOUND

        private:
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

    WSL_MDARRAY_INDEXED_BINARY(+, std::plus<>)
    WSL_MDARRAY_INDEXED_BINARY(-, std::minus<>)
    WSL_MDARRAY_INDEXED_BINARY(*, std::multiplies<>)
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
