#pragma once

// Variables aléatoires (voir docs/variables.md) : applications mesurables X : Ω → E sur un
// espace probabilisé commun. wsl::random::variable(P) ajoute à Ω un facteur indépendant de
// loi P ; les expressions (X + Y, 2 * X, …) sont des fonctions des variables de base,
// évaluées paresseusement.
//
//   namespace wr = wsl::random;
//   wr::variable X{wr::normal()};
//   wr::variable Y{wr::exponential()};
//   auto Z = 2 * X + Y;
//   Z(g);                                  // un ω : X et Y tirées une fois chacune
//   (X - X)(g);                            // 0
//   auto [x, z] = wr::evaluate(g, X, Z);   // même ω
//   wr::E(Z); wr::Var(Z); wr::Cov(X, Z);   // par linéarité : 1, 5, 2
//
// Identité (V1) : une variable de base reçoit un identifiant unique à sa création ; ses
// copies désignent la même fonction de ω. X.independent_copy() : même loi, nouvelle identité.
//
// Évaluation (V2) : à chaque ω, chaque variable de base de l'expression est tirée une fois,
// par l'algorithme de sa loi, dans l'ordre croissant des identifiants (ordre de création) ;
// ses autres occurrences reçoivent la même valeur. L'ordre et le nombre de valeurs
// consommées font partie du contrat.
//
// Arrondi (V3, exact quand c'est possible) : une combinaison affine Σ aᵢ Xᵢ + b de variables
// à tirage paresseux (normale, exponentielle, Laplace, uniforme) est l'arrondi correct de sa
// valeur exacte (réels paresseux partagés par toutes les expressions d'un même ω : `lazy_sum`),
// coefficients et diviseurs exacts ; ailleurs (produits de variables, autres lois, coefficient
// non représentable), calcul flottant sur les valeurs des variables de base, elles-mêmes
// tirées exactement. Z(g, wsl::random::fast) : calcul flottant sur des tirages rapides.
//
// Loi (V4) : wsl::random::law(X) est la loi de X ; law(Z) d'une expression, une
// expression_law : tirage de Z, moments par les cumulants (κₙ(Σ aᵢ Xᵢ) = Σ aᵢⁿ κₙ(Xᵢ),
// variables de base distinctes indépendantes) ; pour une combinaison affine de variables
// d'une même famille close (normale, Cauchy ; gamma à coefficients positifs et sans
// constante, de même intensité), la loi de la famille, résolue à l'exécution (family()),
// dont densité, répartition et quantiles sont ceux de l'expression ; NaN sans identité.
//
// Fonctions : wsl::random::apply(f, X, Y…) applique f (ou une bijection monotone déclarée,
// increasing / decreasing de pushforward.hpp) ; exp, log, sqrt, cbrt, expm1, log1p, atan,
// erf (monotones, déclarées comme telles), sin, cos, tan, abs, pow, hypot, atan2, min, max,
// toutes de wsl::special (bit-exactes). Calcul flottant ; law(f(X)), X variable de base :
// pushforward(loi de X, f) (répartition, quantiles, densité si f est monotone).
//
// Constantes : E, Var et Cov par linéarité sur la forme affine Σ aᵢ Xᵢ + b, les coefficients
// d'une même variable regroupés (Var(X + X) = 4 Var X ; variables de base distinctes
// indépendantes).

#include <wsl/random/affine.hpp>
#include <wsl/random/convolution.hpp>
#include <wsl/random/detail/lazy_sum.hpp>
#include <wsl/random/engines.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/pushforward.hpp>
#include <wsl/random/truncated.hpp>
#include <wsl/special.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <atomic>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace wsl::random
{
    template <class L>
    class variable;
    template <class Op, class... A>
    class expression;

    namespace detail
    {
        inline std::uint64_t new_variable_id() noexcept
        {
            static std::atomic<std::uint64_t> next{1};
            return next.fetch_add(1, std::memory_order_relaxed);
        }

        // Constante d'une expression (scalaire arithmétique).
        template <class S>
        struct constant
        {
            using value_type = S;
            S value;
        };

        template <class T>
        struct is_variable_node : std::false_type {};
        template <class L>
        struct is_variable_node<variable<L>> : std::true_type {};
        template <class Op, class... A>
        struct is_variable_node<expression<Op, A...>> : std::true_type {};

        template <class T>
        concept expression_scalar = std::is_arithmetic_v<T> && !std::same_as<T, bool>;
    }

    // Une variable de base, ou une expression de variables.
    template <class T>
    concept random_variable = detail::is_variable_node<std::remove_cvref_t<T>>::value;

    namespace detail
    {
        // ---- Feuilles : variables de base, dans l'ordre d'apparition ----------------------
        template <class N>
        inline constexpr std::size_t leaf_count = 0;
        template <class L>
        inline constexpr std::size_t leaf_count<variable<L>> = 1;
        template <class Op, class... A>
        inline constexpr std::size_t leaf_count<expression<Op, A...>> = (std::size_t{0} + ... + leaf_count<A>);

        // Rang de la première feuille du K-ième argument.
        template <std::size_t K, class... A>
        constexpr std::size_t leaf_offset() noexcept
        {
            constexpr std::size_t counts[] = {leaf_count<A>..., 0};
            std::size_t s = 0;
            for (std::size_t i = 0; i < K; ++i)
                s += counts[i];
            return s;
        }

        template <class L>
        constexpr auto leaves(variable<L> const& v) { return std::tuple<variable<L> const*>{&v}; }
        template <class S>
        constexpr auto leaves(constant<S> const&) { return std::tuple<>{}; }
        template <class Op, class... A>
        constexpr auto leaves(expression<Op, A...> const& e)
        {
            return std::apply([](auto const&... a) { return std::tuple_cat(leaves(a)...); }, e.arguments());
        }

        // Appelle f.template operator()<J>() pour J == j (j < N, connu à l'exécution).
        template <std::size_t N, class F>
        constexpr void with_index(std::size_t j, F&& f)
        {
            [&]<std::size_t... J>(std::index_sequence<J...>) { (void)((j == J ? (f.template operator()<J>(), true) : false) || ...); }(std::make_index_sequence<N>{});
        }

        // Tirage des feuilles : chaque variable de base une fois, par identifiant croissant.
        template <class G, class... L>
        constexpr std::tuple<typename L::value_type...> draw_leaves(G& g, std::tuple<variable<L> const*...> const& leaf, bool fast = false)
        {
            constexpr std::size_t n = sizeof...(L);
            std::tuple<std::optional<typename L::value_type>...> slots;
            if constexpr (n > 0)
            {
                std::array<std::uint64_t, n> ids{};
                [&]<std::size_t... J>(std::index_sequence<J...>) { ((ids[J] = std::get<J>(leaf)->id()), ...); }(std::make_index_sequence<n>{});
                std::array<std::size_t, n> order{};
                for (std::size_t i = 0; i < n; ++i)
                    order[i] = i;
                std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return ids[a] != ids[b] ? ids[a] < ids[b] : a < b; });
                std::size_t first = 0;
                for (std::size_t p = 0; p < n; ++p)
                {
                    std::size_t const j = order[p];
                    if (p == 0 || ids[j] != ids[order[p - 1]])
                    {
                        first = j;
                        with_index<n>(j, [&]<std::size_t J>() {
                            std::get<J>(slots) = fast ? std::get<J>(leaf)->law()(g, wsl::random::fast) : std::get<J>(leaf)->law()(g);
                        });
                        continue;
                    }
                    with_index<n>(j, [&]<std::size_t J>() {
                        with_index<n>(first, [&]<std::size_t F>() {
                            using vj = std::tuple_element_t<J, std::tuple<typename L::value_type...>>;
                            using vf = std::tuple_element_t<F, std::tuple<typename L::value_type...>>;
                            if constexpr (std::is_same_v<vj, vf>)                // même identifiant : même variable
                                std::get<J>(slots) = std::get<F>(slots);
                        });
                    });
                }
            }
            return std::apply([](auto&... s) { return std::tuple<typename L::value_type...>{std::move(*s)...}; }, slots);
        }

        // ---- Valeur d'un nœud, les feuilles tirées (I : rang de sa première feuille) --------
        template <std::size_t I, class L, class V>
        constexpr auto const& value_at(variable<L> const&, V const& values) { return std::get<I>(values); }
        template <std::size_t I, class S, class V>
        constexpr S value_at(constant<S> const& c, V const&) { return c.value; }
        template <std::size_t I, class Op, class... A, class V>
        constexpr auto value_at(expression<Op, A...> const& e, V const& values)
        {
            return [&]<std::size_t... K>(std::index_sequence<K...>) {
                return e.operation()(value_at<I + leaf_offset<K, A...>()>(std::get<K>(e.arguments()), values)...);
            }(std::index_sequence_for<A...>{});
        }

        template <class N>
        struct node_value;
        template <class L>
        struct node_value<variable<L>> { using type = typename L::value_type; };
        template <class S>
        struct node_value<constant<S>> { using type = S; };
        template <class Op, class... A>
        struct node_value<expression<Op, A...>> { using type = std::remove_cvref_t<std::invoke_result_t<Op const&, typename node_value<A>::type...>>; };

        // ---- Opérations ------------------------------------------------------------------
        struct plus_op { constexpr auto operator()(auto const& a, auto const& b) const { return a + b; } };
        struct minus_op { constexpr auto operator()(auto const& a, auto const& b) const { return a - b; } };
        // Produit protégé de la contraction en fma (-ffp-contract=fast chez l'appelant) : un
        // produit flottant suivi d'une addition ne doit pas changer de résultat selon les options.
        struct multiplies_op
        {
            constexpr auto operator()(auto const& a, auto const& b) const
            {
                using R = decltype(a * b);
                if constexpr (std::floating_point<R>)
                    return static_cast<R>(sd::fmul(double(a), double(b)));
                else
                    return a * b;
            }
        };
        struct divides_op { constexpr auto operator()(auto const& a, auto const& b) const { return a / b; } };
        struct negate_op { constexpr auto operator()(auto const& a) const { return -a; } };
        // Événements : comparaisons (indicatrices), et, ou, non.
        enum class relation
        {
            less,
            less_equal,
            greater,
            greater_equal,
            equal,
            not_equal
        };
        constexpr bool decide(relation r, int sign) noexcept                // sign : signe de a − b
        {
            switch (r)
            {
            case relation::less: return sign < 0;
            case relation::less_equal: return sign <= 0;
            case relation::greater: return sign > 0;
            case relation::greater_equal: return sign >= 0;
            case relation::equal: return sign == 0;
            default: return sign != 0;
            }
        }
        template <relation R>
        struct compare_op
        {
            static constexpr relation kind = R;
            constexpr bool operator()(auto const& a, auto const& b) const
            {
                if constexpr (R == relation::less) return a < b;
                else if constexpr (R == relation::less_equal) return a <= b;
                else if constexpr (R == relation::greater) return a > b;
                else if constexpr (R == relation::greater_equal) return a >= b;
                else if constexpr (R == relation::equal) return a == b;
                else return a != b;
            }
        };
        struct and_op { constexpr bool operator()(bool a, bool b) const { return a && b; } };
        struct or_op { constexpr bool operator()(bool a, bool b) const { return a || b; } };
        struct not_op { constexpr bool operator()(bool a) const { return !a; } };

        // f appliquée aux valeurs ; pour une bijection monotone déclarée, sa partie directe.
        template <class F>
        struct apply_op
        {
            F f;
            constexpr auto operator()(auto const&... a) const
            {
                if constexpr (is_monotone_map<F>::value)
                    return f.forward(a...);
                else
                    return std::invoke(f, a...);
            }
        };

        template <class T>
        constexpr auto as_node(T const& x)
        {
            if constexpr (random_variable<T>)
                return x;
            else
                return constant<T>{x};
        }

        // ---- Forme affine Σ aᵢ Xᵢ + b (constantes par linéarité) -----------------------------
        template <class N>
        inline constexpr bool is_affine = false;
        template <class L>
        inline constexpr bool is_affine<variable<L>> = true;
        template <class S>
        inline constexpr bool is_affine<constant<S>> = true;
        template <class A, class B>
        inline constexpr bool is_affine<expression<plus_op, A, B>> = is_affine<A> && is_affine<B>;
        template <class A, class B>
        inline constexpr bool is_affine<expression<minus_op, A, B>> = is_affine<A> && is_affine<B>;
        template <class A>
        inline constexpr bool is_affine<expression<negate_op, A>> = is_affine<A>;
        template <class A, class S>
        inline constexpr bool is_affine<expression<multiplies_op, A, constant<S>>> = is_affine<A>;
        template <class S, class A>
        inline constexpr bool is_affine<expression<multiplies_op, constant<S>, A>> = is_affine<A>;
        template <class A, class S>
        inline constexpr bool is_affine<expression<divides_op, A, constant<S>>> = is_affine<A>;

        // Constantes connues de chaque feuille.
        template <class N>
        inline constexpr bool leaves_have_mean = true;
        template <class L>
        inline constexpr bool leaves_have_mean<variable<L>> = has_mean<L>;
        template <class Op, class... A>
        inline constexpr bool leaves_have_mean<expression<Op, A...>> = (leaves_have_mean<A> && ...);
        template <class N>
        inline constexpr bool leaves_have_variance = true;
        template <class L>
        inline constexpr bool leaves_have_variance<variable<L>> = has_variance<L>;
        template <class Op, class... A>
        inline constexpr bool leaves_have_variance<expression<Op, A...>> = (leaves_have_variance<A> && ...);

        struct affine_term
        {
            std::uint64_t id;
            double coefficient, mean, variance;
        };
        struct affine_form
        {
            std::vector<affine_term> terms;                                 // regroupés par identifiant, croissants
            double constant = 0;
        };

        template <class N>
        constexpr void collect_affine(N const& node, double scale, affine_form& out)
        {
            if constexpr (requires { typename N::law_type; })
            {
                double m = 0, v = 0;
                if constexpr (has_mean<typename N::law_type>)
                    m = double(node.law().mean());
                if constexpr (has_variance<typename N::law_type>)
                    v = double(node.law().variance());
                out.terms.push_back({node.id(), scale, m, v});
            }
            else if constexpr (requires { node.value; })
                out.constant += sd::fmul(scale, double(node.value));
            else
            {
                using op = typename N::operation_type;
                auto const& args = node.arguments();
                if constexpr (std::is_same_v<op, plus_op>)
                {
                    collect_affine(std::get<0>(args), scale, out);
                    collect_affine(std::get<1>(args), scale, out);
                }
                else if constexpr (std::is_same_v<op, minus_op>)
                {
                    collect_affine(std::get<0>(args), scale, out);
                    collect_affine(std::get<1>(args), -scale, out);
                }
                else if constexpr (std::is_same_v<op, negate_op>)
                    collect_affine(std::get<0>(args), -scale, out);
                else if constexpr (std::is_same_v<op, multiplies_op>)
                {
                    if constexpr (requires { std::get<1>(args).value; })
                        collect_affine(std::get<0>(args), sd::fmul(scale, double(std::get<1>(args).value)), out);
                    else
                        collect_affine(std::get<1>(args), sd::fmul(scale, double(std::get<0>(args).value)), out);
                }
                else                                                        // divides_op, diviseur constant
                    collect_affine(std::get<0>(args), scale / double(std::get<1>(args).value), out);
            }
        }

        template <class N>
        constexpr affine_form affine_of(N const& node)
        {
            affine_form f;
            collect_affine(node, 1.0, f);
            std::ranges::sort(f.terms, [](affine_term const& a, affine_term const& b) { return a.id < b.id; });
            std::vector<affine_term> merged;
            for (auto const& t : f.terms)
                if (!merged.empty() && merged.back().id == t.id)
                    merged.back().coefficient += t.coefficient;
                else
                    merged.push_back(t);
            f.terms = std::move(merged);
            return f;
        }

        template <class N>
        concept affine_with_mean = is_affine<std::remove_cvref_t<N>> && leaves_have_mean<std::remove_cvref_t<N>>;
        template <class N>
        concept affine_with_variance = is_affine<std::remove_cvref_t<N>> && leaves_have_variance<std::remove_cvref_t<N>>;
    }

    namespace detail
    {
        // ---- Voie exacte : combinaisons affines de variables à tirage paresseux ----------------
        template <class N>
        inline constexpr bool lazy_leaves = true;                           // constantes
        template <class L>
        inline constexpr bool lazy_leaves<variable<L>> = lazy_samplable<L, splitmix64>;
        template <class Op, class... A>
        inline constexpr bool lazy_leaves<expression<Op, A...>> = (lazy_leaves<A> && ...);

        template <class N>
        concept exact_affine = is_affine<N> && lazy_leaves<N> && std::floating_point<typename node_value<N>::type>;

        template <class N>
        inline constexpr bool is_exact_comparison = false;
        template <relation R, class A, class B>
        inline constexpr bool is_exact_comparison<expression<compare_op<R>, A, B>> = is_affine<A> && is_affine<B> && lazy_leaves<A> && lazy_leaves<B>
                                                                                   && (leaf_count<A> + leaf_count<B> > 0);
        template <class N>
        concept exact_comparison = is_exact_comparison<N>;
        template <class N>
        inline constexpr bool is_logical = false;
        template <class A, class B>
        inline constexpr bool is_logical<expression<and_op, A, B>> = true;
        template <class A, class B>
        inline constexpr bool is_logical<expression<or_op, A, B>> = true;
        template <class A>
        inline constexpr bool is_logical<expression<not_op, A>> = true;
        template <class N>
        concept logical = is_logical<N>;

        constexpr bool is_power_of_two(double d) noexcept
        {
            return d > 0 && is_finite(d) && (std::bit_cast<std::uint64_t>(d) & ((std::uint64_t{1} << 52) - 1)) == 0
                   && (std::bit_cast<std::uint64_t>(d) >> 52) != 0;
        }
        // a · b exactement en double ?
        constexpr bool exact_product(double a, double b, double& out) noexcept
        {
            auto const [p, e] = sd::two_product(a, b);
            out = p;
            return e == 0 && is_finite(p) && p != 0;
        }

        struct exact_term
        {
            std::uint64_t id;
            std::size_t leaf;                                               // rang de la feuille (dans l'expression)
            dyadic_value numerator;
            double denominator;                                             // coefficient : numerator / denominator
        };
        // Constante numerator / denominator, regroupée par dénominateur avant d'être réduite.
        struct exact_constant
        {
            dyadic_value numerator;
            double denominator;
        };
        struct exact_form
        {
            std::vector<exact_term> terms;
            std::vector<exact_constant> constants;
            bool ok = true;
        };

        // Σ des constantes, dyadique ; vide si l'une des sommes par dénominateur ne l'est pas
        // (n / d dyadique si la mantisse de n est divisible par la partie impaire de d).
        constexpr std::optional<dyadic_value> exact_constant_sum(std::vector<exact_constant> c)
        {
            std::ranges::sort(c, [](exact_constant const& a, exact_constant const& b) { return a.denominator < b.denominator; });
            dyadic_value total;
            for (std::size_t i = 0; i < c.size();)
            {
                dyadic_value n = c[i].numerator;
                double const den = c[i].denominator;
                std::size_t k = i + 1;
                for (; k < c.size() && c[k].denominator == den; ++k)
                    n = n + c[k].numerator;
                i = k;
                if (n.is_zero())
                    continue;
                auto const d = decompose(den);                              // den = m 2^e, m impair
                std::uint64_t const m = d.mantissa >> std::countr_zero(d.mantissa);
                int const e = d.exponent + std::countr_zero(d.mantissa);
                auto const [q, r] = divide(n.m.magnitude, m);
                if (r != 0)
                    return std::nullopt;
                total = total + dyadic_value{{n.m.negative, q}, n.e - e};
            }
            return total;
        }

        template <class N>
        constexpr void collect_exact(N const& node, dyadic_value const& num, double den, exact_form& out, std::size_t& leaf)
        {
            if (!out.ok)
                return;
            if constexpr (requires { typename N::law_type; })
                out.terms.push_back({node.id(), leaf++, num, den});
            else if constexpr (requires { node.value; })
            {
                if (double(node.value) != 0)
                    out.constants.push_back({num * to_dyadic(double(node.value)), den});
            }
            else
            {
                using op = typename N::operation_type;
                auto const& args = node.arguments();
                if constexpr (std::is_same_v<op, plus_op>)
                {
                    collect_exact(std::get<0>(args), num, den, out, leaf);
                    collect_exact(std::get<1>(args), num, den, out, leaf);
                }
                else if constexpr (std::is_same_v<op, minus_op>)
                {
                    collect_exact(std::get<0>(args), num, den, out, leaf);
                    collect_exact(std::get<1>(args), -num, den, out, leaf);
                }
                else if constexpr (std::is_same_v<op, negate_op>)
                    collect_exact(std::get<0>(args), -num, den, out, leaf);
                else if constexpr (std::is_same_v<op, multiplies_op>)
                {
                    if constexpr (requires { std::get<1>(args).value; })
                        collect_exact(std::get<0>(args), num * to_dyadic(double(std::get<1>(args).value)), den, out, leaf);
                    else
                        collect_exact(std::get<1>(args), num * to_dyadic(double(std::get<0>(args).value)), den, out, leaf);
                }
                else                                                        // divides_op, diviseur constant
                {
                    double const d = double(std::get<1>(args).value);
                    double q = 0;
                    if (!(d != 0) || !is_finite(d) || !exact_product(den, d < 0 ? -d : d, q))
                    {
                        out.ok = false;
                        return;
                    }
                    collect_exact(std::get<0>(args), d < 0 ? -num : num, q, out, leaf);
                }
            }
        }

        // Ce que ω a tiré pour une feuille à tirage paresseux : offset + c V / q (ou c √v V).
        struct lazy_leaf
        {
            lazy_real* v = nullptr;                                         // dans le ω (adresse stable)
            dyadic_value offset, coefficient;
            double divisor = 1;
            bool root = false;
            scaled variance{1.0, 0.0, 0};
        };

        // Un ω : les feuilles (de toutes les expressions évaluées ensemble) tirées une fois par
        // identifiant, paresseusement si leur loi le permet, sinon par leur loi.
        template <class G, class... L>
        class omega
        {
        public:
            using values = std::tuple<typename L::value_type...>;
            static constexpr std::size_t n = sizeof...(L);

            omega(omega const&) = delete;                                   // les sommes pointent sur reals_
            omega& operator=(omega const&) = delete;

            constexpr omega(G& g, std::tuple<variable<L> const*...> const& leaf) : g_(g), leaf_(leaf)
            {
                if constexpr (n > 0)
                {
                    std::array<std::uint64_t, n> ids{};
                    [&]<std::size_t... J>(std::index_sequence<J...>) { ((ids[J] = std::get<J>(leaf_)->id()), ...); }(std::make_index_sequence<n>{});
                    std::array<std::size_t, n> order{};
                    for (std::size_t i = 0; i < n; ++i)
                        order[i] = i;
                    std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return ids[a] != ids[b] ? ids[a] < ids[b] : a < b; });
                    for (std::size_t p = 0; p < n; ++p)
                    {
                        std::size_t const j = order[p];
                        if (p == 0 || ids[j] != ids[order[p - 1]])
                        {
                            first_ = j;
                            with_index<n>(j, [&]<std::size_t J>() { draw<J>(); });
                        }
                        else
                            with_index<n>(j, [&]<std::size_t J>() { with_index<n>(first_, [&]<std::size_t F>() { share<J, F>(); }); });
                    }
                }
            }

            // Valeur de la feuille J (arrondie une fois, exactement, puis mémorisée).
            template <std::size_t J>
            constexpr auto const& leaf_value()
            {
                auto& slot = std::get<J>(slots_);
                if (!slot)
                {
                    using T = std::tuple_element_t<J, values>;
                    if constexpr (std::floating_point<T> && lazy_samplable<std::remove_cvref_t<decltype(std::get<J>(leaf_)->law())>, splitmix64>)
                    {
                        lazy_leaf const& d = lazy_[J];
                        lazy_sum s;
                        s.add_constant(d.offset);
                        if (d.root)
                            s.use_root(d.variance);
                        s.add_shared(d.coefficient, s.divisor_slot(d.divisor), d.v, g_, d.root);
                        slot = s.template rounded<T>(g_, {}, {}, 0, to_dyadic(1.0));
                    }
                }
                return *slot;
            }
            constexpr values leaf_values()
            {
                return [&]<std::size_t... J>(std::index_sequence<J...>) { return values{leaf_value<J>()...}; }(std::make_index_sequence<n>{});
            }

            // Combinaison affine exacte (feuilles first.. de l'expression) ; vide si un
            // coefficient, un diviseur ou une constante n'est pas représentable exactement.
            template <class T, class N>
            constexpr std::optional<T> exact(N const& node, std::size_t first)
            {
                exact_form f;
                std::size_t leaf = first;
                collect_exact(node, to_dyadic(1.0), 1.0, f, leaf);
                auto s = exact_sum(std::move(f));
                if (!s)
                    return std::nullopt;
                return s->sum.template rounded<T>(g_, {}, {}, 0, to_dyadic(1.0));
            }

            // Signe exact de A − B (comparaison affine) ; vide hors de la voie exacte.
            template <class A, class B>
            constexpr std::optional<int> exact_sign(A const& a, B const& b, std::size_t first)
            {
                exact_form f;
                std::size_t leaf = first;
                collect_exact(a, to_dyadic(1.0), 1.0, f, leaf);
                collect_exact(b, to_dyadic(-1.0), 1.0, f, leaf);
                auto s = exact_sum(std::move(f));
                if (!s)
                    return std::nullopt;
                if (s->sum.size() == 0)                                     // constante seule (X − X…)
                    return s->constant.is_zero() ? 0 : s->constant.m.negative ? -1 : 1;
                return s->sum.compare(g_, {}, {}, 0, to_dyadic(1.0), {});
            }

            // Somme paresseuse exacte d'une forme affine (et sa constante) ; vide si un
            // coefficient, un diviseur ou une constante n'est pas représentable exactement.
            struct exact_result
            {
                lazy_sum sum;
                dyadic_value constant;
            };
            constexpr std::optional<exact_result> exact_sum(exact_form f)
            {
                if (!f.ok)
                    return std::nullopt;
                std::ranges::sort(f.terms, [](exact_term const& a, exact_term const& b) { return a.id != b.id ? a.id < b.id : a.leaf < b.leaf; });
                lazy_sum s;
                for (std::size_t i = 0; i < f.terms.size();)
                {
                    exact_term t = f.terms[i];                               // regroupe une même variable
                    std::size_t k = i + 1;
                    for (; k < f.terms.size() && f.terms[k].id == t.id; ++k)
                    {
                        if (f.terms[k].denominator != t.denominator)
                            return std::nullopt;
                        t.numerator = t.numerator + f.terms[k].numerator;
                    }
                    i = k;
                    if (t.numerator.is_zero())
                        continue;
                    lazy_leaf const& d = lazy_[t.leaf];
                    double q = 0;
                    if (!exact_product(t.denominator, d.divisor, q))
                        return std::nullopt;
                    if (!d.offset.is_zero())                                 // offset · numerator / denominator
                        f.constants.push_back({d.offset * t.numerator, t.denominator});
                    if (d.root)
                        s.use_root(d.variance);
                    s.add_shared(t.numerator * d.coefficient, s.divisor_slot(q), d.v, g_, d.root);
                }
                auto const constant = exact_constant_sum(std::move(f.constants));
                if (!constant)
                    return std::nullopt;
                s.add_constant(*constant);
                return exact_result{std::move(s), *constant};
            }

            // Valeur d'une expression dont la première feuille est first : combinaison affine et
            // comparaison affine exactes ; et, ou, non sur les valeurs exactes de leurs opérandes.
            template <std::size_t First, class N>
            constexpr auto value_of(N const& node)
            {
                using T = typename node_value<N>::type;
                if constexpr (exact_affine<N>)
                {
                    if (auto const r = exact<T>(node, First))
                        return *r;
                }
                else if constexpr (exact_comparison<N>)
                {
                    auto const& args = node.arguments();
                    if (auto const sign = exact_sign(std::get<0>(args), std::get<1>(args), First))
                        return decide(N::operation_type::kind, *sign);
                }
                else if constexpr (logical<N>)
                    return children<First>(node);
                return value_at<First>(node, leaf_values());
            }
            template <std::size_t First, class Op, class... A>
            constexpr auto children(expression<Op, A...> const& e)
            {
                return [&]<std::size_t... K>(std::index_sequence<K...>) {
                    return e.operation()(value_of<First + leaf_offset<K, A...>()>(std::get<K>(e.arguments()))...);
                }(std::index_sequence_for<A...>{});
            }

        private:
            template <std::size_t J>
            constexpr void draw()
            {
                auto const& law = std::get<J>(leaf_)->law();
                using law_type = std::remove_cvref_t<decltype(law)>;
                if constexpr (lazy_samplable<law_type, G>)
                {
                    lazy_draw d = lazy_sample(law, g_);
                    reals_[J] = std::move(d.v);
                    lazy_[J] = {&reals_[J], std::move(d.offset), std::move(d.coefficient), d.divisor, d.root, d.variance};
                }
                else
                    std::get<J>(slots_) = law(g_);
            }
            template <std::size_t J, std::size_t F>
            constexpr void share()
            {
                if constexpr (std::is_same_v<std::tuple_element_t<J, values>, std::tuple_element_t<F, values>>)   // même variable
                {
                    lazy_[J] = lazy_[F];
                    std::get<J>(slots_) = std::get<F>(slots_);
                }
            }

            G& g_;
            std::tuple<variable<L> const*...> leaf_;
            std::array<lazy_real, n> reals_{};
            std::array<lazy_leaf, n> lazy_{};
            std::tuple<std::optional<typename L::value_type>...> slots_;
            std::size_t first_ = 0;
        };

        template <class G, class... L>
        omega(G&, std::tuple<variable<L> const*...> const&) -> omega<G, L...>;
    }

    // -----------------------------------------------------------------------
    // Variable de base : un facteur de Ω, de loi L.
    // -----------------------------------------------------------------------
    template <class L>
    class variable
    {
    public:
        using law_type = L;
        using value_type = typename L::value_type;

        explicit variable(L law) : law_(std::move(law)), id_(detail::new_variable_id()) {}

        constexpr L const& law() const noexcept { return law_; }
        constexpr std::uint64_t id() const noexcept { return id_; }
        // Même loi, nouvelle identité : une variable indépendante.
        variable independent_copy() const { return variable(law_); }

        template <std::uniform_random_bit_generator G>
        constexpr value_type operator()(G& g) const { return law_(g); }
        template <std::uniform_random_bit_generator G>
        constexpr value_type operator()(G& g, fast_t) const { return law_(g, wsl::random::fast); }
        template <std::uniform_random_bit_generator G>
        constexpr std::vector<value_type> operator()(std::size_t n, G& g) const { return law_(n, g); }

    private:
        L law_;
        std::uint64_t id_;
    };

    template <class L>
    variable(L) -> variable<L>;

    // -----------------------------------------------------------------------
    // Expression : Op appliquée aux arguments (variables, expressions, constantes).
    // -----------------------------------------------------------------------
    template <class Op, class... A>
    class expression
    {
    public:
        using operation_type = Op;
        using value_type = typename detail::node_value<expression>::type;

        constexpr explicit expression(Op op, A... args) : op_(std::move(op)), args_(std::move(args)...) {}

        constexpr Op const& operation() const noexcept { return op_; }
        constexpr std::tuple<A...> const& arguments() const noexcept { return args_; }

        // Un ω neuf : chaque variable de base tirée une fois ; exact pour une combinaison affine.
        template <std::uniform_random_bit_generator G>
        constexpr value_type operator()(G& g) const
        {
            auto const leaf = detail::leaves(*this);
            detail::omega w(g, leaf);
            return w.template value_of<0>(*this);
        }
        // Calcul flottant sur des tirages rapides.
        template <std::uniform_random_bit_generator G>
        constexpr value_type operator()(G& g, fast_t) const
        {
            auto const values = detail::draw_leaves(g, detail::leaves(*this), true);
            return detail::value_at<0>(*this, values);
        }
        template <std::uniform_random_bit_generator G>
        constexpr std::vector<value_type> operator()(std::size_t n, G& g) const
        {
            std::vector<value_type> out;
            out.reserve(n);
            for (std::size_t i = 0; i < n; ++i)
                out.push_back((*this)(g));
            return out;
        }

    private:
        Op op_;
        std::tuple<A...> args_;
    };

    // Plusieurs expressions sur le même ω.
    template <std::uniform_random_bit_generator G, random_variable... E>
    constexpr std::tuple<typename detail::node_value<std::remove_cvref_t<E>>::type...> evaluate(G& g, E const&... e)
    {
        auto const leaf = std::tuple_cat(detail::leaves(e)...);
        detail::omega w(g, leaf);
        constexpr std::size_t counts[] = {detail::leaf_count<std::remove_cvref_t<E>>..., 0};
        return [&]<std::size_t... K>(std::index_sequence<K...>) {
            constexpr auto offset = [](std::size_t k) {
                std::size_t s = 0;
                for (std::size_t i = 0; i < k; ++i)
                    s += counts[i];
                return s;
            };
            return std::tuple<typename detail::node_value<std::remove_cvref_t<E>>::type...>{w.template value_of<offset(K)>(e)...};
        }(std::index_sequence_for<E...>{});
    }

    // -----------------------------------------------------------------------
    // Opérateurs : variables, expressions et scalaires arithmétiques (bool exclu).
    // -----------------------------------------------------------------------
    template <class A, class B>
    concept expression_operands = (random_variable<A> || random_variable<B>) && (random_variable<A> || detail::expression_scalar<A>)
                                  && (random_variable<B> || detail::expression_scalar<B>);

    template <class A, class B>
        requires expression_operands<A, B>
    constexpr auto operator+(A const& a, B const& b)
    {
        return expression<detail::plus_op, decltype(detail::as_node(a)), decltype(detail::as_node(b))>({}, detail::as_node(a), detail::as_node(b));
    }
    template <class A, class B>
        requires expression_operands<A, B>
    constexpr auto operator-(A const& a, B const& b)
    {
        return expression<detail::minus_op, decltype(detail::as_node(a)), decltype(detail::as_node(b))>({}, detail::as_node(a), detail::as_node(b));
    }
    template <class A, class B>
        requires expression_operands<A, B>
    constexpr auto operator*(A const& a, B const& b)
    {
        return expression<detail::multiplies_op, decltype(detail::as_node(a)), decltype(detail::as_node(b))>({}, detail::as_node(a), detail::as_node(b));
    }
    template <class A, class B>
        requires expression_operands<A, B>
    constexpr auto operator/(A const& a, B const& b)
    {
        return expression<detail::divides_op, decltype(detail::as_node(a)), decltype(detail::as_node(b))>({}, detail::as_node(a), detail::as_node(b));
    }
    template <random_variable A>
    constexpr auto operator-(A const& a)
    {
        return expression<detail::negate_op, A>({}, a);
    }

    // Événements : comparaisons (variables indicatrices, valeurs bool), et, ou, non.
#define WSL_RANDOM_COMPARISON(symbol, kind)                                                                                                   \
    template <class A, class B>                                                                                                               \
        requires expression_operands<A, B>                                                                                                    \
    constexpr auto operator symbol(A const& a, B const& b)                                                                                    \
    {                                                                                                                                         \
        return expression<detail::compare_op<detail::relation::kind>, decltype(detail::as_node(a)), decltype(detail::as_node(b))>(           \
            {}, detail::as_node(a), detail::as_node(b));                                                                                      \
    }
    WSL_RANDOM_COMPARISON(<, less)
    WSL_RANDOM_COMPARISON(<=, less_equal)
    WSL_RANDOM_COMPARISON(>, greater)
    WSL_RANDOM_COMPARISON(>=, greater_equal)
    WSL_RANDOM_COMPARISON(==, equal)
    WSL_RANDOM_COMPARISON(!=, not_equal)
#undef WSL_RANDOM_COMPARISON

    template <class A>
    concept random_event = random_variable<A> && std::same_as<typename detail::node_value<std::remove_cvref_t<A>>::type, bool>;

    template <random_event A, random_event B>
    constexpr auto operator&&(A const& a, B const& b) { return expression<detail::and_op, A, B>({}, a, b); }
    template <random_event A, random_event B>
    constexpr auto operator||(A const& a, B const& b) { return expression<detail::or_op, A, B>({}, a, b); }
    template <random_event A>
    constexpr auto operator!(A const& a) { return expression<detail::not_op, A>({}, a); }

    // -----------------------------------------------------------------------
    // Constantes, par linéarité : E, Var, Cov (notation des manuels).
    // -----------------------------------------------------------------------
    namespace detail::cpo
    {
        struct expectation_fn
        {
            template <random_variable Z>
                requires affine_with_mean<Z>
            constexpr double operator()(Z const& z) const
            {
                auto const f = affine_of(z);
                double s = f.constant;
                for (auto const& t : f.terms)
                    s += sd::fmul(t.coefficient, t.mean);
                return s;
            }
        };
        struct variance_of_fn
        {
            template <random_variable Z>
                requires affine_with_variance<Z>
            constexpr double operator()(Z const& z) const
            {
                double s = 0;
                for (auto const& t : affine_of(z).terms)
                    s += sd::fmul(sd::fmul(t.coefficient, t.coefficient), t.variance);
                return s;
            }
        };
        struct covariance_fn
        {
            template <random_variable Z1, random_variable Z2>
                requires affine_with_variance<Z1> && affine_with_variance<Z2>
            constexpr double operator()(Z1 const& a, Z2 const& b) const
            {
                auto const fa = affine_of(a), fb = affine_of(b);
                double s = 0;
                std::size_t j = 0;
                for (auto const& t : fa.terms)                              // termes triés : fusion
                {
                    while (j < fb.terms.size() && fb.terms[j].id < t.id)
                        ++j;
                    if (j < fb.terms.size() && fb.terms[j].id == t.id)
                        s += sd::fmul(sd::fmul(t.coefficient, fb.terms[j].coefficient), t.variance);
                }
                return s;
            }
        };
    }

    inline namespace operations
    {
        inline constexpr detail::cpo::expectation_fn E{};
        inline constexpr detail::cpo::variance_of_fn Var{};
        inline constexpr detail::cpo::covariance_fn Cov{};
    }
    namespace detail
    {
        // Loi commune des feuilles : void sans feuille, mixed_laws si elles diffèrent.
        struct mixed_laws {};
        template <class A, class B>
        using merge_laws = std::conditional_t<std::is_void_v<A>, B, std::conditional_t<std::is_void_v<B> || std::is_same_v<A, B>, A, mixed_laws>>;
        template <class N>
        struct common_leaf_law { using type = void; };
        template <class L>
        struct common_leaf_law<variable<L>> { using type = L; };
        template <class Op, class... A>
        struct common_leaf_law<expression<Op, A...>>
        {
            template <class... T>
            struct fold { using type = void; };
            template <class T, class... R>
            struct fold<T, R...> { using type = merge_laws<T, typename fold<R...>::type>; };
            using type = typename fold<typename common_leaf_law<A>::type...>::type;
        };

        // Famille close : toute image affine et toute somme y restent (normale, Cauchy), ou
        // seulement l'échelle positive et les sommes de même paramètre (gamma).
        template <class L>
        concept closed_family = affine_closed<L> && requires(L const& a) { { sum_closure(a, a) } -> std::same_as<L>; };
        template <class L>
        concept scale_family = scale_closed<L> && requires(L const& a) { { sum_closure(a, a) } -> std::same_as<std::optional<L>>; };

        template <class N>
        using family_of_expression = typename common_leaf_law<N>::type;
        template <class N>
        concept has_family = !std::is_void_v<family_of_expression<N>> && !std::is_same_v<family_of_expression<N>, mixed_laws> && is_affine<N>
                             && (closed_family<family_of_expression<N>> || scale_family<family_of_expression<N>>);
        template <class N>
        struct family_storage { using type = char; };
        template <class N>
            requires has_family<N>
        struct family_storage<N> { using type = std::optional<family_of_expression<N>>; };

        // Termes a·Pᵢ (variables regroupées) et constante, pour une famille L.
        template <class L, class N>
        constexpr void collect_family(N const& node, double scale, std::vector<std::tuple<std::uint64_t, double, L>>& terms, double& constant)
        {
            if constexpr (requires { typename N::law_type; })
                terms.emplace_back(node.id(), scale, node.law());
            else if constexpr (requires { node.value; })
                constant += sd::fmul(scale, double(node.value));
            else
            {
                using op = typename N::operation_type;
                auto const& args = node.arguments();
                if constexpr (std::is_same_v<op, plus_op>)
                {
                    collect_family<L>(std::get<0>(args), scale, terms, constant);
                    collect_family<L>(std::get<1>(args), scale, terms, constant);
                }
                else if constexpr (std::is_same_v<op, minus_op>)
                {
                    collect_family<L>(std::get<0>(args), scale, terms, constant);
                    collect_family<L>(std::get<1>(args), -scale, terms, constant);
                }
                else if constexpr (std::is_same_v<op, negate_op>)
                    collect_family<L>(std::get<0>(args), -scale, terms, constant);
                else if constexpr (std::is_same_v<op, multiplies_op>)
                {
                    if constexpr (requires { std::get<1>(args).value; })
                        collect_family<L>(std::get<0>(args), sd::fmul(scale, double(std::get<1>(args).value)), terms, constant);
                    else
                        collect_family<L>(std::get<1>(args), sd::fmul(scale, double(std::get<0>(args).value)), terms, constant);
                }
                else
                    collect_family<L>(std::get<0>(args), scale / double(std::get<1>(args).value), terms, constant);
            }
        }

        template <class N>
        std::optional<family_of_expression<N>> family_law(N const& node)
        {
            using L = family_of_expression<N>;
            std::vector<std::tuple<std::uint64_t, double, L>> raw;
            double constant = 0;
            collect_family<L>(node, 1.0, raw, constant);
            std::ranges::sort(raw, [](auto const& a, auto const& b) { return std::get<0>(a) < std::get<0>(b); });
            std::vector<std::pair<double, L>> terms;                        // variables regroupées
            for (std::size_t i = 0; i < raw.size(); ++i)
                if (i > 0 && std::get<0>(raw[i]) == std::get<0>(raw[i - 1]))
                    terms.back().first += std::get<1>(raw[i]);
                else
                    terms.emplace_back(std::get<1>(raw[i]), std::get<2>(raw[i]));
            std::erase_if(terms, [](auto const& t) { return t.first == 0; });
            if (terms.empty())
                return std::nullopt;                                        // loi de Dirac : pas de famille
            if constexpr (closed_family<L>)
            {
                std::optional<L> r;
                for (std::size_t i = 0; i < terms.size(); ++i)
                {
                    L const image = affine_closure(terms[i].second, terms[i].first, i == 0 ? constant : 0.0);
                    r = r ? sum_closure(*r, image) : image;
                }
                return r;
            }
            else
            {
                if (constant != 0)
                    return std::nullopt;
                std::optional<L> r;
                for (auto const& [a, law] : terms)
                {
                    if (!(a > 0))
                        return std::nullopt;
                    auto const image = scale_closure(law, a);
                    if (!image)
                        return std::nullopt;
                    r = r ? sum_closure(*r, *image) : image;
                    if (!r)
                        return std::nullopt;
                }
                return r;
            }
        }

        // Cumulants de la forme affine : κ₂, κ₃, κ₄ (feuilles munies de ces constantes).
        template <class N>
        inline constexpr bool leaves_have_shape = true;
        template <class L>
        inline constexpr bool leaves_have_shape<variable<L>> = has_skewness<L> && has_excess_kurtosis<L> && has_variance<L>;
        template <class Op, class... A>
        inline constexpr bool leaves_have_shape<expression<Op, A...>> = (leaves_have_shape<A> && ...);

        // Feuilles d'une forme affine : identifiant, coefficient, variance, asymétrie, kurtosis.
        using shape_term = std::tuple<std::uint64_t, double, double, double, double>;
        template <class M>
        void collect_shape(M const& n, double scale, std::vector<shape_term>& raw)
        {
            if constexpr (requires { typename M::law_type; })
                raw.emplace_back(n.id(), scale, double(n.law().variance()), double(n.law().skewness()), double(n.law().excess_kurtosis()));
            else if constexpr (requires { n.value; })
                return;
            else
            {
                using op = typename M::operation_type;
                auto const& args = n.arguments();
                if constexpr (std::is_same_v<op, plus_op>)
                {
                    collect_shape(std::get<0>(args), scale, raw);
                    collect_shape(std::get<1>(args), scale, raw);
                }
                else if constexpr (std::is_same_v<op, minus_op>)
                {
                    collect_shape(std::get<0>(args), scale, raw);
                    collect_shape(std::get<1>(args), -scale, raw);
                }
                else if constexpr (std::is_same_v<op, negate_op>)
                    collect_shape(std::get<0>(args), -scale, raw);
                else if constexpr (std::is_same_v<op, multiplies_op>)
                {
                    if constexpr (requires { std::get<1>(args).value; })
                        collect_shape(std::get<0>(args), sd::fmul(scale, double(std::get<1>(args).value)), raw);
                    else
                        collect_shape(std::get<1>(args), sd::fmul(scale, double(std::get<0>(args).value)), raw);
                }
                else
                    collect_shape(std::get<0>(args), scale / double(std::get<1>(args).value), raw);
            }
        }

        template <class N>
        std::array<double, 3> cumulants(N const& node)                      // κ₂, κ₃, κ₄
        {
            std::vector<shape_term> raw;
            collect_shape(node, 1.0, raw);
            std::ranges::sort(raw, [](auto const& a, auto const& b) { return std::get<0>(a) < std::get<0>(b); });
            std::array<double, 3> k{};
            for (std::size_t i = 0; i < raw.size();)
            {
                auto [id, a, v, g1, g2] = raw[i];
                std::size_t j = i + 1;
                for (; j < raw.size() && std::get<0>(raw[j]) == id; ++j)
                    a += std::get<1>(raw[j]);
                i = j;
                double const a2 = sd::fmul(a, a), s = wsl::special::sqrt(v);
                k[0] += sd::fmul(a2, v);
                k[1] += sd::fmul(sd::fmul(sd::fmul(a2, a), g1), sd::fmul(v, s));
                k[2] += sd::fmul(sd::fmul(a2, a2), sd::fmul(g2, sd::fmul(v, v)));
            }
            return k;
        }
    }

    // -----------------------------------------------------------------------
    // Loi d'une expression.
    // -----------------------------------------------------------------------
    template <class N>
    class expression_law : public law_interface
    {
    public:
        using value_type = typename detail::node_value<N>::type;
        static constexpr exactness sampling = detail::exact_affine<N> ? exactness::exact : exactness::approximate;
        static constexpr bool has_family = detail::has_family<N>;

        explicit expression_law(N z) : z_(std::move(z))
        {
            if constexpr (has_family)
                family_ = detail::family_law(z_);
        }

        constexpr N const& expression() const noexcept { return z_; }
        // La loi de la famille, quand une identité la donne (sinon vide).
        auto const& family() const noexcept requires has_family { return family_; }

        // Moments : par les cumulants (feuilles indépendantes).
        double mean() const requires detail::affine_with_mean<N> { return E(z_); }
        double variance() const requires detail::affine_with_variance<N> { return Var(z_); }
        double skewness() const requires detail::is_affine<N> && detail::leaves_have_shape<N>
        {
            auto const k = detail::cumulants(z_);
            return k[0] > 0 ? k[1] / (k[0] * wsl::special::sqrt(k[0])) : std::numeric_limits<double>::quiet_NaN();
        }
        double excess_kurtosis() const requires detail::is_affine<N> && detail::leaves_have_shape<N>
        {
            auto const k = detail::cumulants(z_);
            return k[0] > 0 ? k[2] / (k[0] * k[0]) : std::numeric_limits<double>::quiet_NaN();
        }

        // Fonctions de la loi : celles de la famille ; NaN sans identité (règle 3b).
        value_type pdf(value_type x) const requires has_family { return family_ ? family_->pdf(x) : nan(); }
        value_type log_pdf(value_type x) const requires has_family { return family_ ? family_->log_pdf(x) : nan(); }
        value_type cdf(value_type x) const requires has_family { return family_ ? family_->cdf(x) : nan(); }
        value_type ccdf(value_type x) const requires has_family { return family_ ? family_->ccdf(x) : nan(); }
        value_type quantile(double p) const requires has_family { return family_ ? family_->quantile(p) : nan(); }
        value_type cquantile(double q) const requires has_family { return family_ ? family_->cquantile(q) : nan(); }

        template <std::uniform_random_bit_generator G>
        value_type sample(G& g) const { return z_(g); }
        template <std::uniform_random_bit_generator G>
        value_type sample(G& g, fast_t) const { return z_(g, wsl::random::fast); }

    private:
        static constexpr value_type nan() noexcept { return std::numeric_limits<value_type>::quiet_NaN(); }

        N z_;
        typename detail::family_storage<N>::type family_{};
    };

    namespace detail::cpo
    {
        struct law_fn
        {
            template <class L>
            constexpr L const& operator()(variable<L> const& x) const noexcept { return x.law(); }
            // f(X), X variable de base : la mesure image de la loi de X.
            template <class F, class L>
            auto operator()(expression<apply_op<F>, variable<L>> const& z) const
            {
                return pushforward(std::get<0>(z.arguments()).law(), z.operation().f);
            }
            template <class Op, class... A>
            expression_law<expression<Op, A...>> operator()(expression<Op, A...> const& z) const { return expression_law<expression<Op, A...>>(z); }
        };
    }
    inline namespace operations
    {
        inline constexpr detail::cpo::law_fn law{};
    }

    // -----------------------------------------------------------------------
    // Fonctions de variables (lot 1d).
    // -----------------------------------------------------------------------
    template <class F, class... A>
        requires(sizeof...(A) > 0) && (random_variable<A> || ...) && ((random_variable<A> || detail::expression_scalar<A>) && ...)
    constexpr auto apply(F f, A const&... a)
    {
        return expression<detail::apply_op<F>, decltype(detail::as_node(a))...>(detail::apply_op<F>{std::move(f)}, detail::as_node(a)...);
    }

    namespace detail::functions
    {
        namespace sp = wsl::special;
        // Bijections croissantes : f, f⁻¹, (f⁻¹)' (produits protégés de la contraction).
        inline constexpr auto exp = increasing([](auto x) { return sp::exp(x); }, [](auto y) { return sp::log(y); }, [](auto y) { return 1 / y; });
        inline constexpr auto log = increasing([](auto x) { return sp::log(x); }, [](auto y) { return sp::exp(y); }, [](auto y) { return sp::exp(y); });
        inline constexpr auto sqrt = increasing([](auto x) { return sp::sqrt(x); }, [](auto y) { return decltype(y)(sd::fmul(double(y), double(y))); },
                                                [](auto y) { return y + y; });
        inline constexpr auto cbrt = increasing([](auto x) { return sp::cbrt(x); },
                                                [](auto y) { return decltype(y)(sd::fmul(sd::fmul(double(y), double(y)), double(y))); },
                                                [](auto y) { return decltype(y)(sd::fmul(3.0, sd::fmul(double(y), double(y)))); });
        inline constexpr auto expm1 = increasing([](auto x) { return sp::expm1(x); }, [](auto y) { return sp::log1p(y); }, [](auto y) { return 1 / (1 + y); });
        inline constexpr auto log1p = increasing([](auto x) { return sp::log1p(x); }, [](auto y) { return sp::expm1(y); }, [](auto y) { return sp::exp(y); });
        inline constexpr auto atan = increasing([](auto x) { return sp::atan(x); }, [](auto y) { return sp::tan(y); },
                                                [](auto y) { auto const t = sp::tan(y); return decltype(y)(1 + sd::fmul(double(t), double(t))); });
        // (erf⁻¹)'(y) = (√π/2) e^{erf⁻¹(y)²}
        inline constexpr auto erf = increasing([](auto x) { return sp::erf(x); }, [](auto y) { return sp::erf_inv(y); },
                                               [](auto y) {
                                                   double const u = double(sp::erf_inv(y));
                                                   return decltype(y)(sd::fmul(0x1.c5bf891b4ef6bp-1, sp::exp(sd::fmul(u, u))));
                                               });
        // Fonctions ordinaires.
        inline constexpr auto sin = [](auto x) { return sp::sin(x); };
        inline constexpr auto cos = [](auto x) { return sp::cos(x); };
        inline constexpr auto tan = [](auto x) { return sp::tan(x); };
        inline constexpr auto abs = [](auto x) { return x < 0 ? -x : x; };
        inline constexpr auto pow = [](auto x, auto y) { return sp::pow(x, y); };
        inline constexpr auto hypot = [](auto x, auto y) { return sp::hypot(x, y); };
        inline constexpr auto atan2 = [](auto y, auto x) { return sp::atan2(y, x); };
        inline constexpr auto min = [](auto x, auto y) { return y < x ? y : x; };
        inline constexpr auto max = [](auto x, auto y) { return x < y ? y : x; };
    }

    template <random_variable A> constexpr auto exp(A const& a) { return apply(detail::functions::exp, a); }
    template <random_variable A> constexpr auto log(A const& a) { return apply(detail::functions::log, a); }
    template <random_variable A> constexpr auto sqrt(A const& a) { return apply(detail::functions::sqrt, a); }
    template <random_variable A> constexpr auto cbrt(A const& a) { return apply(detail::functions::cbrt, a); }
    template <random_variable A> constexpr auto expm1(A const& a) { return apply(detail::functions::expm1, a); }
    template <random_variable A> constexpr auto log1p(A const& a) { return apply(detail::functions::log1p, a); }
    template <random_variable A> constexpr auto atan(A const& a) { return apply(detail::functions::atan, a); }
    template <random_variable A> constexpr auto erf(A const& a) { return apply(detail::functions::erf, a); }
    template <random_variable A> constexpr auto sin(A const& a) { return apply(detail::functions::sin, a); }
    template <random_variable A> constexpr auto cos(A const& a) { return apply(detail::functions::cos, a); }
    template <random_variable A> constexpr auto tan(A const& a) { return apply(detail::functions::tan, a); }
    template <random_variable A> constexpr auto abs(A const& a) { return apply(detail::functions::abs, a); }
    template <class A, class B> requires expression_operands<A, B> constexpr auto pow(A const& a, B const& b) { return apply(detail::functions::pow, a, b); }
    template <class A, class B> requires expression_operands<A, B> constexpr auto hypot(A const& a, B const& b) { return apply(detail::functions::hypot, a, b); }
    template <class A, class B> requires expression_operands<A, B> constexpr auto atan2(A const& y, B const& x) { return apply(detail::functions::atan2, y, x); }
    template <class A, class B> requires expression_operands<A, B> constexpr auto min(A const& a, B const& b) { return apply(detail::functions::min, a, b); }
    template <class A, class B> requires expression_operands<A, B> constexpr auto max(A const& a, B const& b) { return apply(detail::functions::max, a, b); }

    // -----------------------------------------------------------------------
    // Probabilités (lot 2, E1) : exactes quand la loi les donne, NaN sinon.
    // -----------------------------------------------------------------------
    namespace detail
    {
        constexpr relation flipped(relation r) noexcept                     // c R X  ⇔  X R' c
        {
            switch (r)
            {
            case relation::less: return relation::greater;
            case relation::less_equal: return relation::greater_equal;
            case relation::greater: return relation::less;
            case relation::greater_equal: return relation::less_equal;
            default: return r;
            }
        }

        // P(X > k), k entier.
        template <class L>
        double k_survival(L const& law, double k, double lo, double hi)
        {
            using T = typename L::value_type;
            return k < lo ? 1.0 : k >= hi ? 0.0 : double(law.ccdf(static_cast<T>(k)));
        }

        // P(X R c) par la loi de X (masse comprise pour une loi discrète).
        template <class L>
        double threshold_probability(L const& law, relation r, double c)
        {
            constexpr double nan = std::numeric_limits<double>::quiet_NaN();
            using T = typename L::value_type;
            if constexpr (std::integral<T>)
            {
                if constexpr (has_cdf<L> && has_ccdf<L> && has_pmf<L>)
                {
                    if (c != c)
                        return nan;
                    double const lo = double(law.support().lower), hi = double(law.support().upper);
                    auto const at = [&](double k) { return static_cast<T>(k < lo ? lo - 1 : k > hi ? hi : k); };
                    auto const cdf_at = [&](double k) { return k < lo ? 0.0 : k >= hi ? 1.0 : double(law.cdf(at(k))); };     // P(X <= k)
                    double const fl = std::floor(c), ce = std::ceil(c);
                    switch (r)
                    {
                    case relation::less_equal: return cdf_at(fl);
                    case relation::less: return cdf_at(ce - 1);
                    case relation::greater: return k_survival(law, fl, lo, hi);
                    case relation::greater_equal: return k_survival(law, ce - 1, lo, hi);
                    case relation::equal: return fl == c && c >= lo && c <= hi ? double(law.pmf(static_cast<T>(c))) : 0.0;
                    default: return fl == c && c >= lo && c <= hi ? 1.0 - double(law.pmf(static_cast<T>(c))) : 1.0;
                    }
                }
                else
                    return nan;
            }
            else if constexpr (has_cdf<L> && has_ccdf<L>)                    // loi continue : pas d'atome
            {
                switch (r)
                {
                case relation::less:
                case relation::less_equal: return double(law.cdf(static_cast<T>(c)));
                case relation::greater:
                case relation::greater_equal: return double(law.ccdf(static_cast<T>(c)));
                case relation::equal: return c == c ? 0.0 : nan;
                default: return c == c ? 1.0 : nan;
                }
            }
            else
                return nan;
        }

        // Identifiants des variables de base d'une expression.
        template <class N>
        std::vector<std::uint64_t> ids_of(N const& node)
        {
            std::vector<std::uint64_t> ids;
            std::apply([&](auto const*... v) { (ids.push_back(v->id()), ...); }, leaves(node));
            std::ranges::sort(ids);
            return ids;
        }
        template <class A, class B>
        bool disjoint(A const& a, B const& b)
        {
            auto const x = ids_of(a), y = ids_of(b);
            std::size_t i = 0, j = 0;
            while (i < x.size() && j < y.size())
            {
                if (x[i] == y[j])
                    return false;
                x[i] < y[j] ? ++i : ++j;
            }
            return true;
        }

        template <class N>
        double probability(N const& node)
        {
            constexpr double nan = std::numeric_limits<double>::quiet_NaN();
            if constexpr (requires { typename N::law_type; })                // variable indicatrice de base
            {
                if constexpr (has_pmf<typename N::law_type>)
                    return double(node.law().pmf(true));
                else
                    return nan;
            }
            else
            {
                using op = typename N::operation_type;
                auto const& args = node.arguments();
                if constexpr (std::is_same_v<op, not_op>)
                    return 1 - probability(std::get<0>(args));
                else if constexpr (std::is_same_v<op, and_op> || std::is_same_v<op, or_op>)
                {
                    if (!disjoint(std::get<0>(args), std::get<1>(args)))     // même variable : pas de produit
                        return nan;
                    double const a = probability(std::get<0>(args)), b = probability(std::get<1>(args));
                    return std::is_same_v<op, and_op> ? sd::fmul(a, b) : a + b - sd::fmul(a, b);
                }
                else if constexpr (requires { op::kind; })
                {
                    using A = std::remove_cvref_t<decltype(std::get<0>(args))>;
                    using B = std::remove_cvref_t<decltype(std::get<1>(args))>;
                    if constexpr (requires { typename A::law_type; } && requires { std::get<1>(args).value; })
                        return threshold_probability(std::get<0>(args).law(), op::kind, double(std::get<1>(args).value));
                    else if constexpr (requires { std::get<0>(args).value; } && requires { typename B::law_type; })
                        return threshold_probability(std::get<1>(args).law(), flipped(op::kind), double(std::get<0>(args).value));
                    else if constexpr (is_affine<A> && is_affine<B> && has_family<expression<minus_op, A, B>>)
                    {
                        auto const difference = expression<minus_op, A, B>({}, std::get<0>(args), std::get<1>(args));
                        auto const law = family_law(difference);           // A − B, de loi connue
                        return law ? threshold_probability(*law, op::kind, 0.0) : nan;
                    }
                    else
                        return nan;
                }
                else
                    return nan;
            }
        }
    }

    namespace detail::cpo
    {
        struct probability_fn
        {
            template <random_event A>
            double operator()(A const& a) const { return detail::probability(a); }
        };
    }
    inline namespace operations
    {
        inline constexpr detail::cpo::probability_fn P{};
    }

    // -----------------------------------------------------------------------
    // Conditionnement (lot 2, E2) : given(Z, A), variable de base nouvelle, de loi celle de
    // Z sachant A ; tirage par rejet (ω tirés jusqu'à ce que A soit vraie).
    // -----------------------------------------------------------------------
    template <random_variable Z, random_event A>
    class conditional_law : public law_interface
    {
        template <class N>
        static constexpr bool base_variable = requires { typename N::law_type; };
        static constexpr bool truncatable = [] {
            if constexpr (base_variable<Z>)
                return real_law<typename Z::law_type> && has_cdf<typename Z::law_type> && has_ccdf<typename Z::law_type>
                       && std::floating_point<typename Z::law_type::value_type>;
            else
                return false;
        }();
        static constexpr bool known = base_variable<Z> && truncatable;

    public:
        using value_type = typename detail::node_value<Z>::type;
        static constexpr exactness sampling =
            (detail::exact_affine<Z> || requires { typename Z::law_type; }) && (detail::exact_comparison<A> || detail::logical<A>) ? exactness::exact
                                                                                                                            : exactness::approximate;

        conditional_law(Z z, A a) : z_(std::move(z)), a_(std::move(a))
        {
            if (detail::probability(a_) == 0)
                throw std::domain_error("wsl::random::given : événement de probabilité nulle");
            if constexpr (requires { typename Z::law_type; })
            {
                if (detail::disjoint(z_, a_))                               // A indépendant de Z : la loi de Z
                    independent_ = true;
                else if constexpr (truncatable)
                    truncated_ = truncation();
            }
        }

        constexpr Z const& variable_of() const noexcept { return z_; }
        constexpr A const& event() const noexcept { return a_; }

        // Fonctions de la loi, quand une identité la donne : indépendance (loi de Z), seuil sur
        // Z (loi tronquée) ; NaN sinon.
        value_type cdf(value_type x) const requires known { return call([&](auto const& l) { return l.cdf(x); }); }
        value_type ccdf(value_type x) const requires known { return call([&](auto const& l) { return l.ccdf(x); }); }
        value_type pdf(value_type x) const requires known { return call([&](auto const& l) { return l.pdf(x); }); }
        value_type quantile(double p) const requires known { return call([&](auto const& l) { return l.quantile(p); }); }

        template <std::uniform_random_bit_generator G>
        value_type sample(G& g) const
        {
            for (;;)
            {
                auto const [z, ok] = evaluate(g, z_, a_);
                if (ok)
                    return z;
            }
        }

    private:

        // A = seuil sur Z lui-même (Z R c ou c R Z) : la loi tronquée.
        auto truncation() const
        {
            using L = typename Z::law_type;
            using T = typename L::value_type;
            std::optional<truncated<L>> out;
            if constexpr (requires { A::operation_type::kind; })
            {
                auto const& args = a_.arguments();
                using X = std::remove_cvref_t<decltype(std::get<0>(args))>;
                using Y = std::remove_cvref_t<decltype(std::get<1>(args))>;
                auto make = [&](detail::relation r, double c) {
                    constexpr double inf = std::numeric_limits<double>::infinity();
                    if (r == detail::relation::greater || r == detail::relation::greater_equal)
                        out = truncate(z_.law(), {.lower = c, .upper = inf});
                    else if (r == detail::relation::less || r == detail::relation::less_equal)
                        out = truncate(z_.law(), {.lower = -inf, .upper = c});
                };
                if constexpr (std::is_same_v<X, Z> && requires { std::get<1>(args).value; })
                {
                    if (std::get<0>(args).id() == z_.id())
                        make(A::operation_type::kind, double(std::get<1>(args).value));
                }
                else if constexpr (std::is_same_v<Y, Z> && requires { std::get<0>(args).value; })
                {
                    if (std::get<1>(args).id() == z_.id())
                        make(detail::flipped(A::operation_type::kind), double(std::get<0>(args).value));
                }
            }
            (void)sizeof(T);
            return out;
        }

        template <class F>
        value_type call(F const& f) const
        {
            if (independent_)
                return static_cast<value_type>(f(z_.law()));
            if constexpr (truncatable)
                if (truncated_)
                    return static_cast<value_type>(f(*truncated_));
            return std::numeric_limits<value_type>::quiet_NaN();
        }

        Z z_;
        A a_;
        bool independent_ = false;
        std::conditional_t<truncatable, std::optional<truncated<typename std::conditional_t<base_variable<Z>, Z, variable<normal<>>>::law_type>>, char> truncated_{};
    };

    // Variable aléatoire de loi celle de Z sachant A (nouvelle identité).
    template <random_variable Z, random_event A>
    auto given(Z const& z, A const& a)
    {
        return variable<conditional_law<Z, A>>{conditional_law<Z, A>(z, a)};
    }
}
