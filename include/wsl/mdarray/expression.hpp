#pragma once

#include <wsl/modes.hpp>

// Monde des tableaux entiers (design/mdarray.md, V8, I6, I11, I13) : opérations terme
// à terme, paresseuses, sur des opérandes de même forme ; un scalaire se diffuse.
//
//   y = x * x;   r = hypot(x, y);   mask = a < b;   z = where(a > 0, a, 0);
//   wsl::mdarray c = a * b + 1;     auto t = sum(a);   auto F = wsl::vectorize(f);
//
// Une expression n'est utilisable que comme temporaire (L2) : `auto e = a + b;`
// compile, mais `e` ne sert à rien (message clair) ; std::move(e) reste possible.
// Comparaisons : élément par élément (des masques, comme std::valarray) ;
// `if (masque)` est refusé ; equal(a, b) rend un bool.

#include <wsl/mdarray/axis.hpp>
#include <wsl/mdarray/container.hpp>
#include <wsl/mdarray/detail/expression_base.hpp>
#include <wsl/mdarray/detail/function_list.hpp>
#include <wsl/mdarray/detail/iterate.hpp>
#include <wsl/mdarray/detail/precondition.hpp>
#include <wsl/mdarray/detail/reduce.hpp>
#include <wsl/mdarray/extents.hpp>
#include <wsl/mdarray/mdspan.hpp>
#include <wsl/mdarray/position.hpp>

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace wsl::inline WSL_DEFAULTS
{
    namespace detail::mdarray
    {
        // Opérande « tableau » : un tableau, une vue ou une expression.
        template <class T>
        concept array_operand = array_like<T> || expression<T>;

        // Opérande scalaire : tout ce qui n'est pas l'un de nos objets d'indices.
        template <class T>
        concept scalar_operand = !array_operand<T> && !axis_like<T> && !position_like<T> && !displacement_like<T>
                                 && !is_extents<std::remove_cvref_t<T>> && !indexed<T>;

        template <class T>
        concept operand = array_operand<T> || scalar_operand<T>;

        template <class... T>
        concept with_array = (operand<T> && ...) && (array_operand<T> || ...);

        // --- Feuilles ----------------------------------------------------------

        template <class View>
        struct view_leaf
        {
            View view;
            using extents_type = typename View::extents_type;
            constexpr extents_type const& extents() const noexcept { return view.extents(); }
            template <class... I>
            constexpr decltype(auto) operator()(I... i) const
            {
                return view[i...];
            }
            constexpr void const* storage() const noexcept { return view.data_handle(); }

            // Accès direct : la k-ième valeur en ordre lexicographique est la k-ième du
            // stockage (ligne majeure, accesseur par défaut).
            static constexpr bool flat_ok = std::is_same_v<typename View::layout_type, wsl::layout_right>
                                            && std::is_same_v<typename View::accessor_type,
                                                              wsl::default_accessor<typename View::element_type>>;
            constexpr decltype(auto) flat(std::size_t k) const { return view.data_handle()[k]; }
            // W valeurs consécutives du stockage.
            using lane_value = std::remove_cv_t<typename View::element_type>;
            static constexpr bool lane_barrier = false;
            template <std::size_t W>
            constexpr auto flat_lanes(std::size_t k) const
            {
                lane_pack<std::remove_cv_t<typename View::element_type>, W> r;
                for (std::size_t l = 0; l < W; ++l)
                    r.v[l] = view.data_handle()[k + l];
                return r;
            }
        };

        template <class T>
        struct scalar_leaf
        {
            static constexpr bool flat_ok = true;
            constexpr T const& flat(std::size_t) const noexcept { return value; }
            using lane_value = T;
            static constexpr bool lane_barrier = false;
            template <std::size_t W>
            constexpr lane_pack<T, W> flat_lanes(std::size_t) const
            {
                lane_pack<T, W> r;
                for (std::size_t l = 0; l < W; ++l)
                    r.v[l] = value;
                return r;
            }
            T value;
            template <class... I>
            constexpr T const& operator()(I...) const noexcept
            {
                return value;
            }
        };

        template <class T>
        inline constexpr bool is_scalar_leaf = false;
        template <class T>
        inline constexpr bool is_scalar_leaf<scalar_leaf<T>> = true;

        template <class T>
        constexpr auto make_operand(T&& t)
        {
            using U = std::remove_cvref_t<T>;
            if constexpr (expression<U>)
            {
                static_assert(!std::is_lvalue_reference_v<T>,
                              "wsl::mdarray : une expression ne se garde pas (L2) ; écrire le calcul dans l'affectation "
                              "(y = a + b), ou std::move(e) si c'est voulu");
                return U(std::move(t));
            }
            else if constexpr (is_mdarray<U>)
                return view_leaf<typename U::const_mdspan_type>{std::as_const(t).to_mdspan()};
            else if constexpr (is_mdspan<U>)
                return view_leaf<U>{t};
            else
                return scalar_leaf<U>{std::forward<T>(t)};
        }

        template <class T>
        using operand_t = decltype(make_operand(std::declval<T>()));

        // --- Étendues communes ------------------------------------------------------

        // Étendues de deux opérandes : même rang, étendues fixes égales ; une fixe
        // l'emporte sur une dynamique.
        template <class E1, class E2, class Seq = std::make_index_sequence<E1::rank()>>
        struct combine_extents;

        template <class E1, class E2, std::size_t... R>
        struct combine_extents<E1, E2, std::index_sequence<R...>>
        {
            static_assert(E1::rank() == E2::rank(), "wsl::mdarray : opérandes de rangs différents (V8 : même forme)");
            static_assert(((E1::static_extent(R) == dynamic_extent || E2::static_extent(R) == dynamic_extent
                            || E1::static_extent(R) == E2::static_extent(R))
                           && ...),
                          "wsl::mdarray : opérandes d'étendues fixes différentes (V8 : même forme)");
            using type = wsl::extents<std::common_type_t<typename E1::index_type, typename E2::index_type>,
                                      (E1::static_extent(R) != dynamic_extent ? E1::static_extent(R) : E2::static_extent(R))...>;
        };

        // Étendues communes des feuilles (les scalaires n'en ont pas).
        template <class E, class... Leaves>
        struct fold_extents
        {
            using type = E;
        };
        template <class E, class L, class... Rest>
            requires is_scalar_leaf<L>
        struct fold_extents<E, L, Rest...>
        {
            using type = typename fold_extents<E, Rest...>::type;
        };
        template <class L, class... Rest>
            requires(!is_scalar_leaf<L>)
        struct fold_extents<void, L, Rest...>
        {
            using type = typename fold_extents<typename L::extents_type, Rest...>::type;
        };
        template <class E, class L, class... Rest>
            requires(!is_scalar_leaf<L> && !std::is_void_v<E>)
        struct fold_extents<E, L, Rest...>
        {
            using type = typename fold_extents<typename combine_extents<E, typename L::extents_type>::type, Rest...>::type;
        };

        // Type des valeurs : op appliquée aux valeurs des feuilles en une position.
        template <class Op, class Index, class Seq, class... Leaves>
        struct value_type_of;
        template <class Op, class Index, std::size_t... R, class... Leaves>
        struct value_type_of<Op, Index, std::index_sequence<R...>, Leaves...>
        {
            using type = std::remove_cvref_t<decltype(std::declval<Op const&>()(
                std::declval<Leaves const&>()(std::declval<std::conditional_t<true, Index, std::integral_constant<std::size_t, R>>>()...)...))>;
        };

        template <class Leaf>
        constexpr auto const* extents_of(Leaf const& leaf) noexcept
        {
            if constexpr (is_scalar_leaf<Leaf>)
                return static_cast<void const*>(nullptr);
            else
                return &leaf.extents();
        }

        // --- Nœud élément par élément --------------------------------------------------

        template <class Op, class... Leaves>
        class elementwise : public expression_tag
        {
        public:
            using extents_type = typename fold_extents<void, Leaves...>::type;
            using index_type = typename extents_type::index_type;
            using rank_type = std::size_t;
            using value_type =
                typename value_type_of<Op, index_type, std::make_index_sequence<extents_type::rank()>, Leaves...>::type;
            using op_type = Op;

            static constexpr rank_type rank() noexcept { return extents_type::rank(); }

            constexpr elementwise(Op op, Leaves... leaves) : op_(std::move(op)), leaves_(std::move(leaves)...)
            {
                std::array<index_type, rank()> v{};
                bool first = true;
                std::apply(
                    [&](auto const&... l) {
                        (
                            [&](auto const& leaf) {
                                if constexpr (!is_scalar_leaf<std::remove_cvref_t<decltype(leaf)>>)
                                {
                                    if (first)
                                    {
                                        for (rank_type r = 0; r < rank(); ++r)
                                            v[r] = static_cast<index_type>(leaf.extents().extent(r));
                                        first = false;
                                    }
                                    else
                                        for (rank_type r = 0; r < rank(); ++r)
                                            precondition(std::cmp_equal(leaf.extents().extent(r), v[r]),
                                                         "wsl::mdarray : opérandes de formes différentes (V8)");
                                }
                            }(l),
                            ...);
                    },
                    leaves_);
                if constexpr (rank() > 0)
                    extents_ = extents_type(v);
            }

            elementwise(elementwise const&) = delete;
            constexpr elementwise(elementwise&&) = default;

            constexpr extents_type const& extents() const noexcept { return extents_; }

            template <class... I>
            constexpr value_type operator()(I... i) const
            {
                return std::apply([&](auto const&... l) { return op_(l(i...)...); }, leaves_);
            }

            // Accès direct, quand toutes les feuilles le permettent (formes égales).
            static constexpr bool flat_ok = (Leaves::flat_ok && ...);
            constexpr value_type flat(std::size_t k) const
            {
                return std::apply([&](auto const&... l) { return op_(l.flat(k)...); }, leaves_);
            }

            // Paquets de W valeurs consécutives (voies) : une opération qui traite un paquet
            // entier (produit reproductible : une barrière par registre vectoriel au lieu
            // d'une par produit, mode.hpp), sinon voie par voie ; mêmes bits que flat.
            // lane_barrier : l'expression contient une telle opération (sinon, la boucle
            // simple se vectorise déjà).
            using lane_value = value_type;
            static constexpr bool lane_barrier =
                requires(Op const& o) { o.lanes(std::declval<lane_pack<typename Leaves::lane_value, 8>>()...); }
                || (Leaves::lane_barrier || ...);
            template <std::size_t W>
            constexpr lane_pack<value_type, W> flat_lanes(std::size_t k) const
            {
                return std::apply(
                    [&](auto const&... l) {
                        return [&](auto const&... x) {
                            if constexpr (requires { op_.lanes(x...); })
                            {
                                auto const p = op_.lanes(x...);
                                lane_pack<value_type, W> r;
                                for (std::size_t q = 0; q < W; ++q)
                                    r.v[q] = static_cast<value_type>(p.v[q]);
                                return r;
                            }
                            else
                            {
                                lane_pack<value_type, W> r;
                                for (std::size_t q = 0; q < W; ++q)
                                    r.v[q] = op_(x.v[q]...);
                                return r;
                            }
                        }(l.template flat_lanes<W>(k)...);
                    },
                    leaves_);
            }

            // Un élément d'une expression : (a + b)[1, 2].
            template <class... I>
                requires(sizeof...(I) == rank() && (index_integer<I> && ...))
            constexpr value_type operator[](I... i) const
            {
                precondition(extents_.contains(i...), "wsl::mdarray : position hors de l'expression");
                return (*this)(static_cast<index_type>(i)...);
            }

            // (a + b)[i, j] : une expression indicée (monde des axes).
            template <class... A>
                requires(sizeof...(A) == rank() && (axis_like<A> || ...)
                         && ((axis_like<A> || index_integer<std::remove_cvref_t<A>>) && ...))
            constexpr auto operator[](A const&... a) &&
            {
                return make_indexed(std::move(*this), a...);
            }

            template <class I>
            constexpr value_type operator[](mdposition<rank(), I> const& p) const
            {
                return [&]<std::size_t... K>(std::index_sequence<K...>) { return (*this)[p[K]...]; }(std::make_index_sequence<rank()>{});
            }

            // Le même arbre, dans un mode (mode.hpp).
            template <class Mode>
            constexpr auto rebuild() &&
            {
                return std::apply(
                    [&](auto&&... l) {
                        return elementwise<decltype(with_mode_op<Mode>(std::move(op_))), decltype(with_mode<Mode>(std::move(l)))...>(
                            with_mode_op<Mode>(std::move(op_)), with_mode<Mode>(std::move(l))...);
                    },
                    std::move(leaves_));
            }

            // Plages mémoire lues (aliasing, I12).
            template <class F>
            constexpr void for_each_storage(F&& f) const
            {
                std::apply([&](auto const&... l) { (visit_storage(l, f), ...); }, leaves_);
            }

            // `if (masque)` : refusé (I13).
            template <class T = void>
            explicit operator bool() const
            {
                static_assert(sizeof(T) == 0, "wsl::mdarray : un masque n'est pas un booléen ; écrire all(m), any(m) ou count(m)");
                return false;
            }

        private:
            template <class L, class F>
            static constexpr void visit_storage(L const& l, F& f)
            {
                if constexpr (requires { l.storage(); })
                    f(l.view);
                else if constexpr (expression<L>)
                    l.for_each_storage(f);
            }

            [[no_unique_address]] Op op_;
            std::tuple<Leaves...> leaves_;
            [[no_unique_address]] extents_type extents_{};
        };

        template <class Op, class... T>
        constexpr auto make_elementwise(Op op, T&&... t)
        {
            return elementwise<Op, operand_t<T>...>(std::move(op), make_operand(std::forward<T>(t))...);
        }

        // L'opération est-elle définie sur les éléments ?
        template <class T>
        struct element_of
        {
            using type = std::remove_cvref_t<T>;
        };
        template <class T>
            requires array_operand<T>
        struct element_of<T>
        {
            using type = typename std::remove_cvref_t<T>::value_type;
        };
        template <class T>
        using element_t = typename element_of<T>::type;

        template <class Op, class... T>
        concept elementwise_valid = with_array<T...> && requires(Op op, element_t<T>... e) { op(e...); };

        // Valeur d'un tableau ou d'une expression en une position.
        template <class A, class... I>
        constexpr decltype(auto) value_at(A const& a, I... i)
        {
            if constexpr (expression<A>)
                return a(i...);
            else
                return as_view(a)[i...];
        }

        // Toutes les valeurs, dans l'ordre lexicographique des positions.
        template <class A>
        constexpr auto const& extents_of_array(A const& a) noexcept
        {
            if constexpr (expression<A>)
                return a.extents();
            else
                return a.extents();
        }

        struct negate_op
        {
            template <class T>
            constexpr auto operator()(T const& a) const -> decltype(-a)
            {
                return -a;
            }
        };
        struct identity_op
        {
            template <class T>
            constexpr auto operator()(T const& a) const -> decltype(+a)
            {
                return +a;
            }
        };
        struct logical_not_op
        {
            template <class T>
            constexpr auto operator()(T const& a) const -> decltype(!a)
            {
                return !a;
            }
        };
        struct logical_and_op
        {
            template <class A, class B>
            constexpr auto operator()(A const& a, B const& b) const -> decltype(a && b)
            {
                return a && b; // pas de court-circuit : les deux opérandes sont déjà calculés
            }
        };
        struct logical_or_op
        {
            template <class A, class B>
            constexpr auto operator()(A const& a, B const& b) const -> decltype(a || b)
            {
                return a || b;
            }
        };

        // min, max élément par élément, propageant NaN (I6).
        struct min_op
        {
            static constexpr bool is_function = true;
            template <class A, class B>
            constexpr auto operator()(A const& a, B const& b) const
            {
                using C = std::common_type_t<A, B>;
                if constexpr (std::is_floating_point_v<C>)
                {
                    if (a != a)
                        return static_cast<C>(a);
                    if (b != b)
                        return static_cast<C>(b);
                }
                return static_cast<C>(b < a ? b : a);
            }
        };
        struct max_op
        {
            static constexpr bool is_function = true;
            template <class A, class B>
            constexpr auto operator()(A const& a, B const& b) const
            {
                using C = std::common_type_t<A, B>;
                if constexpr (std::is_floating_point_v<C>)
                {
                    if (a != a)
                        return static_cast<C>(a);
                    if (b != b)
                        return static_cast<C>(b);
                }
                return static_cast<C>(a < b ? b : a);
            }
        };
        struct where_op
        {
            static constexpr bool is_function = true;
            template <class C, class A, class B>
            constexpr auto operator()(C const& c, A const& a, B const& b) const
            {
                using R = std::common_type_t<A, B>;
                return c ? static_cast<R>(a) : static_cast<R>(b);
            }
        };

        template <class F>
        struct call_op
        {
            static constexpr bool is_function = true;
            F f;
            template <class... T>
            constexpr auto operator()(T const&... t) const -> decltype(std::invoke(f, t...))
            {
                return std::invoke(f, t...);
            }
        };
    }

    // --- Opérateurs (V8), trouvés par ADL --------------------------------------------

#define WSL_MDARRAY_BINARY_OPERATOR(op, functor)                                                                       \
    template <class A, class B>                                                                                        \
        requires detail::mdarray::elementwise_valid<functor, A, B>                                                     \
    constexpr auto operator op(A&& a, B&& b)                                                                           \
    {                                                                                                                  \
        return detail::mdarray::make_elementwise(functor{}, std::forward<A>(a), std::forward<B>(b));                   \
    }

    WSL_MDARRAY_BINARY_OPERATOR(+, detail::mdarray::add_op)
    WSL_MDARRAY_BINARY_OPERATOR(-, detail::mdarray::subtract_op)
    WSL_MDARRAY_BINARY_OPERATOR(*, detail::mdarray::multiply_op)
    WSL_MDARRAY_BINARY_OPERATOR(/, std::divides<>)
    WSL_MDARRAY_BINARY_OPERATOR(%, std::modulus<>)
    WSL_MDARRAY_BINARY_OPERATOR(==, std::equal_to<>)
    WSL_MDARRAY_BINARY_OPERATOR(!=, std::not_equal_to<>)
    WSL_MDARRAY_BINARY_OPERATOR(<, std::less<>)
    WSL_MDARRAY_BINARY_OPERATOR(<=, std::less_equal<>)
    WSL_MDARRAY_BINARY_OPERATOR(>, std::greater<>)
    WSL_MDARRAY_BINARY_OPERATOR(>=, std::greater_equal<>)
    WSL_MDARRAY_BINARY_OPERATOR(&&, detail::mdarray::logical_and_op)
    WSL_MDARRAY_BINARY_OPERATOR(||, detail::mdarray::logical_or_op)
#undef WSL_MDARRAY_BINARY_OPERATOR

    template <class A>
        requires detail::mdarray::elementwise_valid<detail::mdarray::negate_op, A>
    constexpr auto operator-(A&& a)
    {
        return detail::mdarray::make_elementwise(detail::mdarray::negate_op{}, std::forward<A>(a));
    }

    template <class A>
        requires detail::mdarray::elementwise_valid<detail::mdarray::identity_op, A>
    constexpr auto operator+(A&& a)
    {
        return detail::mdarray::make_elementwise(detail::mdarray::identity_op{}, std::forward<A>(a));
    }

    template <class A>
        requires detail::mdarray::elementwise_valid<detail::mdarray::logical_not_op, A>
    constexpr auto operator!(A&& a)
    {
        return detail::mdarray::make_elementwise(detail::mdarray::logical_not_op{}, std::forward<A>(a));
    }

    // Affectations composées d'un tableau : y += x, y *= 2 (mêmes positions : sur place).
#define WSL_MDARRAY_COMPOUND(op, base)                                                                                  \
    template <class Element, class Extents, class Layout, class Container, class B>                                    \
        requires(detail::mdarray::array_operand<B> || detail::mdarray::scalar_operand<B>)                              \
                && requires(Element& x, detail::mdarray::element_t<B> y) { x op y; }                                     \
    constexpr mdarray<Element, Extents, Layout, Container>& operator op(mdarray<Element, Extents, Layout, Container>& a, B&& b) \
    {                                                                                                                  \
        return a = (std::as_const(a) base std::forward<B>(b));                                                         \
    }
    WSL_MDARRAY_COMPOUND(+=, +)
    WSL_MDARRAY_COMPOUND(-=, -)
    WSL_MDARRAY_COMPOUND(*=, *)
    WSL_MDARRAY_COMPOUND(/=, /)
#undef WSL_MDARRAY_COMPOUND

    // --- Fonctions élément par élément (I11) -----------------------------------------

#define WSL_MDARRAY_UNARY_FUNCTION(name)                                                                               \
    namespace detail::mdarray                                                                                          \
    {                                                                                                                  \
        struct name##_op                                                                                               \
        {                                                                                                              \
            static constexpr bool is_function = true;                                                                  \
            template <class T>                                                                                         \
            constexpr auto operator()(T const& x) const -> decltype(std::name(x))                                      \
            {                                                                                                          \
                return std::name(x);                                                                                   \
            }                                                                                                          \
        };                                                                                                             \
    }                                                                                                                  \
    template <class A>                                                                                                 \
        requires detail::mdarray::elementwise_valid<detail::mdarray::name##_op, A>                                     \
    constexpr auto name(A&& a)                                                                                         \
    {                                                                                                                  \
        return detail::mdarray::make_elementwise(detail::mdarray::name##_op{}, std::forward<A>(a));                    \
    }

#define WSL_MDARRAY_BINARY_FUNCTION(name)                                                                              \
    namespace detail::mdarray                                                                                          \
    {                                                                                                                  \
        struct name##_op                                                                                               \
        {                                                                                                              \
            static constexpr bool is_function = true;                                                                  \
            template <class T, class U>                                                                                \
            constexpr auto operator()(T const& x, U const& y) const -> decltype(std::name(x, y))                       \
            {                                                                                                          \
                return std::name(x, y);                                                                                \
            }                                                                                                          \
        };                                                                                                             \
    }                                                                                                                  \
    template <class A, class B>                                                                                        \
        requires detail::mdarray::elementwise_valid<detail::mdarray::name##_op, A, B>                                  \
    constexpr auto name(A&& a, B&& b)                                                                                  \
    {                                                                                                                  \
        return detail::mdarray::make_elementwise(detail::mdarray::name##_op{}, std::forward<A>(a), std::forward<B>(b)); \
    }

    WSL_MDARRAY_UNARY_FUNCTIONS(WSL_MDARRAY_UNARY_FUNCTION)
    WSL_MDARRAY_BINARY_FUNCTIONS(WSL_MDARRAY_BINARY_FUNCTION)
#undef WSL_MDARRAY_UNARY_FUNCTION
#undef WSL_MDARRAY_BINARY_FUNCTION

    // min, max à deux arguments : élément par élément, NaN propagé (I6).
    template <class A, class B>
        requires detail::mdarray::elementwise_valid<detail::mdarray::min_op, A, B>
    constexpr auto min(A&& a, B&& b)
    {
        return detail::mdarray::make_elementwise(detail::mdarray::min_op{}, std::forward<A>(a), std::forward<B>(b));
    }

    template <class A, class B>
        requires detail::mdarray::elementwise_valid<detail::mdarray::max_op, A, B>
    constexpr auto max(A&& a, B&& b)
    {
        return detail::mdarray::make_elementwise(detail::mdarray::max_op{}, std::forward<A>(a), std::forward<B>(b));
    }

    // where(c, x, y) : x là où c, y ailleurs (I10).
    template <class C, class A, class B>
        requires detail::mdarray::elementwise_valid<detail::mdarray::where_op, C, A, B>
    constexpr auto where(C&& c, A&& a, B&& b)
    {
        return detail::mdarray::make_elementwise(detail::mdarray::where_op{}, std::forward<C>(c), std::forward<A>(a),
                                                 std::forward<B>(b));
    }

    // Fonctions de l'utilisateur (I11) : F(x, z) élément par élément ; F(1.0, 2.0)
    // reste un appel ordinaire.
    template <class F>
    class vectorized
    {
    public:
        constexpr explicit vectorized(F f) : f_(std::move(f)) {}

        template <class... T>
        constexpr auto operator()(T&&... t) const
        {
            if constexpr ((detail::mdarray::indexed<T> || ...) || (detail::mdarray::axis_like<T> || ...))
                return make_indexed_call(std::forward<T>(t)...);
            else if constexpr (detail::mdarray::with_array<T...>)
                return detail::mdarray::make_elementwise(detail::mdarray::call_op<F>{f_}, std::forward<T>(t)...);
            else
                return std::invoke(f_, std::forward<T>(t)...);
        }

    private:
        // Défini dans indexed.hpp (monde des axes) : F(x[i], x[j]).
        template <class... T>
        constexpr auto make_indexed_call(T&&... t) const;

        F f_;
    };

    template <class F>
    constexpr vectorized<std::decay_t<F>> vectorize(F&& f)
    {
        return vectorized<std::decay_t<F>>(std::forward<F>(f));
    }

    // --- Réductions de tableaux entiers (I6), calculées immédiatement ------------------

    namespace detail::mdarray
    {
        template <class A>
        concept reducible = array_operand<A>;

        // Opérande d'une réduction : une expression prend le mode par défaut.
        template <class A>
        constexpr decltype(auto) reduction_operand(A&& a)
        {
            if constexpr (expression<A> && !std::is_lvalue_reference_v<A>)
                return default_mode(std::move(a));
            else
                return (a);
        }

        // Les valeurs dans l'ordre lexicographique, accédées par rang k.
        template <class A>
        struct flat
        {
            A const& a;
            std::array<std::size_t, std::remove_cvref_t<A>::rank()> strides{};

            constexpr explicit flat(A const& array) : a(array)
            {
                constexpr std::size_t R = std::remove_cvref_t<A>::rank();
                if constexpr (R > 0)
                {
                    strides[R - 1] = 1;
                    for (std::size_t r = R - 1; r-- > 0;)
                        strides[r] = strides[r + 1] * static_cast<std::size_t>(extents_of_array(a).extent(r + 1));
                }
            }

            constexpr std::size_t size() const
            {
                std::size_t n = 1;
                for (std::size_t r = 0; r < std::remove_cvref_t<A>::rank(); ++r)
                    n *= static_cast<std::size_t>(extents_of_array(a).extent(r));
                return n;
            }

            static constexpr bool direct = [] {
                using U = std::remove_cvref_t<A>;
                if constexpr (expression<U>)
                    return U::flat_ok;
                else if constexpr (is_mdarray<U> || is_mdspan<U>)
                    return std::is_same_v<typename U::layout_type, wsl::layout_right>;
                else
                    return false;
            }();

            constexpr auto operator()(std::size_t k) const
            {
                if constexpr (direct)
                {
                    if constexpr (expression<std::remove_cvref_t<A>>)
                        return a.flat(k);
                    else
                        return as_view(a).data_handle()[k];
                }
                else
                    return unravelled(k);
            }

            // 8 valeurs consécutives à la fois, pour la somme par paires (reduce.hpp) :
            // seulement quand l'expression a une barrière à poser par paquets.
            template <std::size_t W>
                requires direct && expression<std::remove_cvref_t<A>> && std::remove_cvref_t<A>::lane_barrier
            constexpr auto pack(std::size_t k) const
            {
                return a.template flat_lanes<W>(k);
            }

            constexpr auto unravelled(std::size_t k) const
            {
                using I = typename std::remove_cvref_t<decltype(extents_of_array(a))>::index_type;
                return [&]<std::size_t... R>(std::index_sequence<R...>) {
                    return value_at(a, static_cast<I>((k / strides[R]) % static_cast<std::size_t>(extents_of_array(a).extent(R)))...);
                }(std::make_index_sequence<std::remove_cvref_t<A>::rank()>{});
            }
        };

        template <class A>
        using value_of = std::remove_cvref_t<decltype(std::declval<flat<A>>()(0))>;
    }

    // Somme de tous les éléments, dans l'ordre fixé d'I9.
    template <detail::mdarray::reducible A>
    constexpr auto sum(A&& source)
    {
        auto&& a = detail::mdarray::reduction_operand(std::forward<A>(source));
        detail::mdarray::flat<std::remove_cvref_t<decltype(a)>> x(a);
        using V = detail::mdarray::value_of<std::remove_cvref_t<decltype(a)>>;
        using T = decltype(std::declval<V>() + std::declval<V>());
        return detail::mdarray::pairwise<T>(x, 0, x.size(), T{}, detail::mdarray::add_op{});
    }

    template <detail::mdarray::reducible A>
    constexpr auto prod(A&& source)
    {
        auto&& a = detail::mdarray::reduction_operand(std::forward<A>(source));
        detail::mdarray::flat<std::remove_cvref_t<decltype(a)>> x(a);
        using V = detail::mdarray::value_of<std::remove_cvref_t<decltype(a)>>;
        using T = decltype(std::declval<V>() * std::declval<V>());
        return detail::mdarray::pairwise<T>(x, 0, x.size(), T{1}, detail::mdarray::multiply_op{});
    }

    // min, max d'un tableau entier (précondition : non vide), NaN propagé.
    template <detail::mdarray::reducible A>
    constexpr auto min(A&& source)
    {
        auto&& a = detail::mdarray::reduction_operand(std::forward<A>(source));
        detail::mdarray::flat<std::remove_cvref_t<decltype(a)>> x(a);
        std::size_t const n = x.size();
        detail::mdarray::precondition(n > 0, "min : tableau vide");
        auto m = x(0);
        for (std::size_t k = 1; k < n; ++k)
            m = detail::mdarray::min_op{}(m, x(k));
        return m;
    }

    template <detail::mdarray::reducible A>
    constexpr auto max(A&& source)
    {
        auto&& a = detail::mdarray::reduction_operand(std::forward<A>(source));
        detail::mdarray::flat<std::remove_cvref_t<decltype(a)>> x(a);
        std::size_t const n = x.size();
        detail::mdarray::precondition(n > 0, "max : tableau vide");
        auto m = x(0);
        for (std::size_t k = 1; k < n; ++k)
            m = detail::mdarray::max_op{}(m, x(k));
        return m;
    }

    namespace detail::mdarray
    {
        template <class A, class Better>
        constexpr auto arg_extremum(A const& a, Better better)
        {
            flat<A> x(a);
            std::size_t const n = x.size();
            precondition(n > 0, "argmin, argmax : tableau vide");
            std::size_t best = 0;
            auto v = x(0);
            for (std::size_t k = 1; k < n && v == v; ++k)
            {
                auto const w = x(k);
                if (w != w || better(w, v))
                {
                    best = k;
                    v = w;
                }
            }
            constexpr std::size_t R = A::rank();
            using I = typename std::remove_cvref_t<decltype(extents_of_array(a))>::index_type;
            std::array<I, R> p{};
            for (std::size_t r = 0; r < R; ++r)
                p[r] = static_cast<I>((best / x.strides[r]) % static_cast<std::size_t>(extents_of_array(a).extent(r)));
            return wsl::mdposition<R, I>(p);
        }
    }

    // Position du minimum, du maximum : la première en ordre lexicographique ; un NaN
    // l'emporte (comme NumPy).
    template <detail::mdarray::reducible A>
    constexpr auto argmin(A&& source)
    {
        auto&& a = detail::mdarray::reduction_operand(std::forward<A>(source));
        return detail::mdarray::arg_extremum(a, [](auto const& w, auto const& v) { return w < v; });
    }

    template <detail::mdarray::reducible A>
    constexpr auto argmax(A&& source)
    {
        auto&& a = detail::mdarray::reduction_operand(std::forward<A>(source));
        return detail::mdarray::arg_extremum(a, [](auto const& w, auto const& v) { return v < w; });
    }

    // Masques (I10, I13).
    template <detail::mdarray::reducible A>
    constexpr bool all(A&& source)
    {
        auto&& a = detail::mdarray::reduction_operand(std::forward<A>(source));
        return detail::mdarray::for_each_index(detail::mdarray::extents_of_array(a),
                                               [&](auto... i) { return static_cast<bool>(detail::mdarray::value_at(a, i...)); });
    }

    template <detail::mdarray::reducible A>
    constexpr bool any(A&& source)
    {
        auto&& a = detail::mdarray::reduction_operand(std::forward<A>(source));
        return !detail::mdarray::for_each_index(detail::mdarray::extents_of_array(a),
                                                [&](auto... i) { return !static_cast<bool>(detail::mdarray::value_at(a, i...)); });
    }

    template <detail::mdarray::reducible A>
    constexpr auto count(A&& source)
    {
        auto&& a = detail::mdarray::reduction_operand(std::forward<A>(source));
        using I = typename std::remove_cvref_t<decltype(detail::mdarray::extents_of_array(a))>::index_type;
        std::make_signed_t<I> n = 0;
        detail::mdarray::for_each_index(detail::mdarray::extents_of_array(a), [&](auto... i) {
            if (static_cast<bool>(detail::mdarray::value_at(a, i...)))
                ++n;
        });
        return n;
    }

    // Construction d'un tableau depuis une expression : mdarray c = a * b + 1; (V8)
    template <class Op, class... Leaves>
    mdarray(detail::mdarray::elementwise<Op, Leaves...>&&)
        -> mdarray<typename detail::mdarray::elementwise<Op, Leaves...>::value_type,
                   typename detail::mdarray::elementwise<Op, Leaves...>::extents_type>;
}
