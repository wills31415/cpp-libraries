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
// consommées font partie du contrat. L'expression est ensuite calculée en flottant.
//
// Constantes : E, Var et Cov par linéarité sur la forme affine Σ aᵢ Xᵢ + b, les coefficients
// d'une même variable regroupés (Var(X + X) = 4 Var X ; variables de base distinctes
// indépendantes).

#include <wsl/random/law.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <concepts>
#include <cstddef>
#include <cstdint>
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
        constexpr std::tuple<typename L::value_type...> draw_leaves(G& g, std::tuple<variable<L> const*...> const& leaf)
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
                std::ranges::stable_sort(order, [&](std::size_t a, std::size_t b) { return ids[a] < ids[b]; });
                std::size_t first = 0;
                for (std::size_t p = 0; p < n; ++p)
                {
                    std::size_t const j = order[p];
                    if (p == 0 || ids[j] != ids[order[p - 1]])
                    {
                        first = j;
                        with_index<n>(j, [&]<std::size_t J>() { std::get<J>(slots) = std::get<J>(leaf)->law()(g); });
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
        struct multiplies_op { constexpr auto operator()(auto const& a, auto const& b) const { return a * b; } };
        struct divides_op { constexpr auto operator()(auto const& a, auto const& b) const { return a / b; } };
        struct negate_op { constexpr auto operator()(auto const& a) const { return -a; } };

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
                out.constant += scale * double(node.value);
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
                        collect_affine(std::get<0>(args), scale * double(std::get<1>(args).value), out);
                    else
                        collect_affine(std::get<1>(args), scale * double(std::get<0>(args).value), out);
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
            std::ranges::stable_sort(f.terms, [](affine_term const& a, affine_term const& b) { return a.id < b.id; });
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

        // Un ω neuf : chaque variable de base tirée une fois, puis l'expression en flottant.
        template <std::uniform_random_bit_generator G>
        constexpr value_type operator()(G& g) const
        {
            auto const values = detail::draw_leaves(g, detail::leaves(*this));
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
        auto const values = detail::draw_leaves(g, leaf);
        constexpr std::size_t counts[] = {detail::leaf_count<std::remove_cvref_t<E>>..., 0};
        return [&]<std::size_t... K>(std::index_sequence<K...>) {
            constexpr auto offset = [](std::size_t k) {
                std::size_t s = 0;
                for (std::size_t i = 0; i < k; ++i)
                    s += counts[i];
                return s;
            };
            return std::tuple<typename detail::node_value<std::remove_cvref_t<E>>::type...>{detail::value_at<offset(K)>(e, values)...};
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
                    s += t.coefficient * t.mean;
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
                    s += t.coefficient * t.coefficient * t.variance;
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
                        s += t.coefficient * fb.terms[j].coefficient * t.variance;
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
}
