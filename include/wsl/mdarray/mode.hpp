#pragma once

#include <wsl/modes.hpp>

// Modes de calcul flottant (design/mdarray.md, I9 ; design/interfaces.md, décisions 24 à 27) :
//
//   y[i] = wsl::portable(A[i, j] * x[j] + b[i]);       // mêmes bits partout
//   y = wsl::deterministic(exp(x) * 2);                 // fusion en fma permise, <cmath>
//
// Déterministe (défaut) : fusion en fma permise, fonctions de <cmath> ; mêmes résultats pour
// une même configuration (compilateur, options, machine). Portable : aucun produit ne se
// fusionne (barrière de wsl::special), fonctions de wsl::special (bit-exactes, constexpr) ;
// une fonction que wsl::special ne couvre pas est une erreur de compilation, comme une
// fonction de l'utilisateur non signée wsl::assume_portable(f). L'ordre des sommes est fixé
// dans les deux modes. Mode des instructions non marquées : wsl::default_reproducibility
// (option CMake WSL_PORTABLE, la même dans tout le projet).

#include <wsl/mdarray/expression.hpp>
#include <wsl/mdarray/indexed.hpp>
#include <wsl/special.hpp>

#include <type_traits>
#include <utility>


namespace wsl::inline WSL_DEFAULTS
{
    namespace detail::mdarray
    {
        using fast_mode = deterministic_t;
        using reproducible_mode = portable_t;

        // Une opération, dans un mode.
        template <class Mode, class Op>
        struct moded_op;

        template <class Op>
        struct moded_op<fast_mode, Op> : Op
        {
            constexpr moded_op() = default;
            constexpr explicit moded_op(Op op) : Op(std::move(op)) {}
        };

        template <class T>
        constexpr auto as_floating(T const& x) noexcept
        {
            if constexpr (std::is_integral_v<T>)
                return static_cast<double>(x);
            else
                return x;
        }

        // Portable : par défaut, l'opération telle quelle (exacte ou sans risque de
        // fusion : +, −, ÷, comparaisons, min, max, where, abs, floor…).
        template <class Op>
        struct moded_op<reproducible_mode, Op> : Op
        {
            constexpr moded_op() = default;
            constexpr explicit moded_op(Op op) : Op(std::move(op)) {}
        };

        // Fonction de l'utilisateur : la bibliothèque ne peut pas inspecter son corps, qui
        // peut fusionner des produits ; l'utilisateur la signe (décision 25).
        template <class F>
        struct moded_op<reproducible_mode, call_op<F>> : call_op<F>
        {
            static_assert(is_assumed_portable<F>,
                          "wsl::portable : une fonction de l'utilisateur doit être signée wsl::assume_portable(f), "
                          "la garantie reposant sur son corps (design/interfaces.md, décision 25)");
            constexpr moded_op() = default;
            constexpr explicit moded_op(call_op<F> op) : call_op<F>(std::move(op)) {}
        };

        // Barrière sur un paquet de voies : comme wsl::special::detail::opaque, mais un
        // registre vectoriel à la fois, pour que les produits du paquet restent vectorisés.
        template <class T, std::size_t W>
        constexpr void opaque_lanes(T (&v)[W]) noexcept
        {
            if consteval
            {
                return;
            }
            if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>)
            {
#if defined(__AVX__)
                constexpr std::size_t bytes = 32;
#else
                constexpr std::size_t bytes = 16;
#endif
                constexpr std::size_t chunk = bytes / sizeof(T);
                using vector [[gnu::vector_size(bytes)]] = T;
                std::size_t l = 0;
                for (; l + chunk <= W; l += chunk)
                {
                    vector x;
                    __builtin_memcpy(&x, v + l, bytes);
#if defined(__x86_64__)
                    asm("" : "+x"(x));
#elif defined(__aarch64__)
                    asm("" : "+w"(x));
#else
                    asm("" : "+m"(x));
#endif
                    __builtin_memcpy(v + l, &x, bytes);
                }
                for (; l < W; ++l)
                    v[l] = wsl::special::detail::opaque(v[l]);
            }
            else
                for (std::size_t l = 0; l < W; ++l)
                    v[l] = wsl::special::detail::opaque(v[l]);
        }

        // Produit : jamais fusionné avec une addition.
        template <>
        struct moded_op<reproducible_mode, detail::mdarray::multiply_op>
        {
            constexpr moded_op() = default;
            constexpr explicit moded_op(detail::mdarray::multiply_op) {}
            template <class A, class B>
            constexpr auto operator()(A const& a, B const& b) const
            {
                using C = decltype(a * b);
                if constexpr (std::is_floating_point_v<C>)
                    return wsl::special::detail::fmul(static_cast<C>(a), static_cast<C>(b));
                else
                    return a * b;
            }
            // Un paquet de voies : les W produits, puis une barrière par registre vectoriel.
            template <class A, class B, std::size_t W>
            constexpr auto lanes(lane_pack<A, W> const& a, lane_pack<B, W> const& b) const
            {
                using C = decltype(a.v[0] * b.v[0]);
                lane_pack<C, W> r;
                for (std::size_t l = 0; l < W; ++l)
                    r.v[l] = static_cast<C>(a.v[l]) * static_cast<C>(b.v[l]);
                if constexpr (std::is_floating_point_v<C>)
                    opaque_lanes(r.v);
                return r;
            }
        };

#define WSL_MDARRAY_REPRODUCIBLE_UNARY(name)                                                                           \
    template <>                                                                                                        \
    struct moded_op<reproducible_mode, name##_op>                                                                      \
    {                                                                                                                  \
        static constexpr bool is_function = true;                                                                      \
        constexpr moded_op() = default;                                                                                \
        constexpr explicit moded_op(name##_op) {}                                                                      \
        template <class T>                                                                                             \
        constexpr auto operator()(T const& x) const                                                                    \
        {                                                                                                              \
            return wsl::special::name(as_floating(x));                                                                 \
        }                                                                                                              \
    };
#define WSL_MDARRAY_REPRODUCIBLE_BINARY(name)                                                                          \
    template <>                                                                                                        \
    struct moded_op<reproducible_mode, name##_op>                                                                      \
    {                                                                                                                  \
        static constexpr bool is_function = true;                                                                      \
        constexpr moded_op() = default;                                                                                \
        constexpr explicit moded_op(name##_op) {}                                                                      \
        template <class T, class U>                                                                                    \
        constexpr auto operator()(T const& x, U const& y) const                                                        \
        {                                                                                                              \
            using C = std::common_type_t<decltype(as_floating(x)), decltype(as_floating(y))>;                          \
            return wsl::special::name(static_cast<C>(x), static_cast<C>(y));                                           \
        }                                                                                                              \
    };
#define WSL_MDARRAY_REPRODUCIBLE_UNCOVERED(name)                                                                       \
    template <>                                                                                                        \
    struct moded_op<reproducible_mode, name##_op>                                                                      \
    {                                                                                                                  \
        static constexpr bool is_function = true;                                                                      \
        constexpr moded_op() = default;                                                                                \
        constexpr explicit moded_op(name##_op) {}                                                                      \
        template <class T>                                                                                             \
        constexpr auto operator()(T const& x) const                                                                    \
        {                                                                                                              \
            static_assert(sizeof(T) == 0, "wsl::mdarray : " #name " n'est pas couverte par wsl::special : pas de mode " \
                                          "portable pour elle (I9) ; wsl::deterministic(…) pour la calculer");             \
            return x;                                                                                                  \
        }                                                                                                              \
    };

        WSL_MDARRAY_REPRODUCIBLE_UNARY(sqrt)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(cbrt)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(exp)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(exp2)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(expm1)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(log)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(log2)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(log10)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(log1p)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(sin)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(cos)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(tan)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(atan)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(erf)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(erfc)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(tgamma)
        WSL_MDARRAY_REPRODUCIBLE_UNARY(lgamma)
        WSL_MDARRAY_REPRODUCIBLE_BINARY(pow)
        WSL_MDARRAY_REPRODUCIBLE_BINARY(atan2)
        WSL_MDARRAY_REPRODUCIBLE_BINARY(hypot)
        WSL_MDARRAY_REPRODUCIBLE_UNCOVERED(asin)
        WSL_MDARRAY_REPRODUCIBLE_UNCOVERED(acos)
        WSL_MDARRAY_REPRODUCIBLE_UNCOVERED(sinh)
        WSL_MDARRAY_REPRODUCIBLE_UNCOVERED(cosh)
        WSL_MDARRAY_REPRODUCIBLE_UNCOVERED(tanh)
        WSL_MDARRAY_REPRODUCIBLE_UNCOVERED(asinh)
        WSL_MDARRAY_REPRODUCIBLE_UNCOVERED(acosh)
        WSL_MDARRAY_REPRODUCIBLE_UNCOVERED(atanh)
#undef WSL_MDARRAY_REPRODUCIBLE_UNARY
#undef WSL_MDARRAY_REPRODUCIBLE_BINARY
#undef WSL_MDARRAY_REPRODUCIBLE_UNCOVERED

        // Une opération déjà dans un mode le garde (wsl::deterministic dans une instruction
        // portable, et inversement).
        template <class Mode, class Op>
        struct mode_of_op
        {
            using type = moded_op<Mode, Op>;
        };
        template <class Mode, class Inner, class Op>
        struct mode_of_op<Mode, moded_op<Inner, Op>>
        {
            using type = moded_op<Inner, Op>;
        };

        template <class Mode, class Op>
        constexpr auto with_mode_op(Op op)
        {
            // (déclarée dans expression_base.hpp)
            using M = typename mode_of_op<Mode, Op>::type;
            if constexpr (std::is_same_v<M, Op>)
                return op;
            else
                return M(std::move(op));
        }

        template <class T>
        inline constexpr bool is_moded_op = false;
        template <class Mode, class Op>
        inline constexpr bool is_moded_op<moded_op<Mode, Op>> = true;

        // --- Réécriture d'un arbre dans un mode -------------------------------------------

        template <class Mode, class T>
        constexpr auto with_mode(T&& t)
        {
            using U = std::remove_cvref_t<T>;
            if constexpr (requires { std::move(t).template rebuild<Mode>(); })
                return std::move(t).template rebuild<Mode>();
            else if constexpr (requires { t.op; t.nodes; })
                return std::apply(
                    [&](auto&&... c) {
                        return indexed_node<decltype(with_mode_op<Mode>(std::move(t.op))),
                                            decltype(with_mode<Mode>(std::move(c)))...>(with_mode_op<Mode>(std::move(t.op)),
                                                                                       with_mode<Mode>(std::move(c))...);
                    },
                    std::move(t.nodes));
            else if constexpr (requires { t.child; })
                return reduction_with_mode<Mode>(std::move(t));
            else if constexpr (requires { t.source.expr; })
                return leaf_with_mode<Mode>(std::move(t));
            else
                return U(std::move(t)); // feuilles : rien à réécrire
        }

        template <class Mode, class Kind, class Bound, class Child>
        constexpr auto reduction_with_mode(indexed_reduction<Kind, Bound, Child>&& r)
        {
            auto child = with_mode<Mode>(std::move(r.child));
            return indexed_reduction<Kind, Bound, decltype(child)>(std::move(child));
        }

        template <class Mode, class Expr, class... Slots>
        constexpr auto leaf_with_mode(indexed_leaf<expression_source<Expr>, Slots...>&& l)
        {
            auto expr = with_mode<Mode>(std::move(l.source.expr));
            return std::apply(
                [&](auto const&... s) {
                    return indexed_leaf<expression_source<decltype(expr)>, Slots...>(expression_source<decltype(expr)>{std::move(expr)},
                                                                                    s...);
                },
                l.slots);
        }

        // Une racine déjà mise dans un mode (son opération l'indique) ; une feuille n'a
        // rien à réécrire.
        template <class T>
        constexpr bool has_mode() noexcept
        {
            using U = std::remove_cvref_t<T>;
            if constexpr (requires { typename U::op_type; })
                return is_moded_op<typename U::op_type>;
            else if constexpr (requires(U u) { u.op; })
                return is_moded_op<std::remove_cvref_t<decltype(std::declval<U>().op)>>;
            else if constexpr (requires(U u) { u.child; })
                return has_mode<decltype(std::declval<U>().child)>();
            else
                return true;
        }

        using default_mode_t = default_reproducibility;

        // Mode par défaut d'une instruction non marquée.
        template <class T>
        constexpr decltype(auto) default_mode(T&& t)
        {
            if constexpr (std::is_same_v<default_mode_t, fast_mode> || has_mode<T>())
                return std::forward<T>(t);
            else
                return with_mode<default_mode_t>(std::move(t));
        }
    }

    // wsl::deterministic(e), wsl::portable(e) : le mode d'une expression (I9), par les
    // étiquettes de <wsl/modes>, qui appellent apply_mode.
    template <reproducibility_mode M, class E>
        requires(detail::mdarray::expression<E> || detail::mdarray::indexed<E>)
    constexpr auto apply_mode(M, E&& e)
    {
        static_assert(!std::is_lvalue_reference_v<E>, "wsl::mdarray : une expression ne se garde pas (L2)");
        return detail::mdarray::with_mode<M>(std::move(e));
    }
}
