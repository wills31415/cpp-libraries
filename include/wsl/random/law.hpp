#pragma once

// Le concept de loi de probabilité (voir docs/random.md §4.3 et docs/interfaces.md, décisions 8).
//
// Une loi est un type qui déclare :
//   - value_type : le type des valeurs tirées ;
//   - sampling   : l'exactitude de son tirage par défaut (exact ou approximate) ;
//   - support()  : l'intervalle [lower, upper] des valeurs possibles, quand il a un
//                  sens (has_support : pas pour les objets combinatoires) ;
//   - sa mesure de référence (reference_measure : lebesgue, counting…), par rapport à
//     laquelle density est la dérivée de Radon–Nikodym dP/dμ (déduite de pdf ou pmf à
//     défaut de déclaration) ;
//   - ses constantes, quand elles existent : mean(), variance(), skewness(),
//     excess_kurtosis(), entropy(), median(), mode() ;
//   - ses opérations : sample, density / log_density (pdf / log_pdf pour Lebesgue,
//     pmf / log_pmf pour le comptage), cdf, ccdf (survie, 1 − cdf calculée sans
//     annulation), quantile, cquantile (quantile de survie : x tel que P(X > x) = q).
//
// Toutes s'appellent par les objets-fonctions de wsl::random — X.pdf(x),
// wsl::random::mean(X), X(g) —, qui prennent le membre de la loi
// (X.pdf(x), X.mean(), X.sample(g)) ou, à défaut, une fonction libre trouvée par ADL (lois
// tierces). Les lois qui héritent de law_interface ont en plus X(g), X(g, wsl::random::fast),
// X(n, g), X.density(x).
//
// Ce qui n'existe pas mathématiquement — l'espérance d'une loi de Cauchy — ou
// n'est pas implémenté n'est pas déclaré : l'appel ne compile pas, et les
// concepts has_mean, has_cdf… permettent de le tester.

#include <wsl/random/uniform.hpp>

#include <concepts>
#include <cstddef>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

namespace wsl::random
{
    enum class exactness
    {
        exact,          // valeur rendue : l'arrondi (règle documentée) d'un réel tiré exactement selon la loi
        approximate     // algorithme classique, erreur documentée ; toujours bit-exact
    };

    // Demande explicite du tirage rapide : sample(X, g, wsl::random::fast).
    struct fast_t
    {
        explicit fast_t() = default;
    };
    inline constexpr fast_t fast{};

    // Intervalle fermé [lower, upper] (bornes éventuellement infinies).
    template <class T>
    struct interval
    {
        T lower, upper;
        friend constexpr bool operator==(interval const&, interval const&) = default;
    };

    template <class L>
    concept probability_law = requires {
        typename L::value_type;
        { L::sampling } -> std::convertible_to<exactness>;
    };

    // -----------------------------------------------------------------------
    // Mesures de référence : la densité est la dérivée de Radon–Nikodym dP/dμ.
    // -----------------------------------------------------------------------
    struct lebesgue
    {
    };
    struct counting
    {
    };
    // Mesure de surface (de Hausdorff) d'une sphère.
    struct surface
    {
    };

    namespace detail::cpo
    {
        // Pilules empoisonnées : l'appel non qualifié ne trouve que les fonctions de la loi (ADL).
        void pdf() = delete;
        void log_pdf() = delete;
        void pmf() = delete;
        void log_pmf() = delete;
        void cdf() = delete;
        void ccdf() = delete;
        void quantile() = delete;
        void cquantile() = delete;
        void sample() = delete;
        void mean() = delete;
        void variance() = delete;
        void skewness() = delete;
        void excess_kurtosis() = delete;
        void entropy() = delete;
        void median() = delete;
        void mode() = delete;
        void support() = delete;

// Le membre de la loi s'il existe, sinon la fonction libre trouvée par ADL.
#define WSL_RANDOM_OPERATION(name)                                                                   \
        struct name##_fn                                                                             \
        {                                                                                            \
            template <class L, class... A>                                                           \
                requires requires(L const& law, A&&... a) { law.name(std::forward<A>(a)...); }       \
                         || requires(L const& law, A&&... a) { name(law, std::forward<A>(a)...); }   \
            constexpr decltype(auto) operator()(L const& law, A&&... a) const                        \
            {                                                                                        \
                if constexpr (requires { law.name(std::forward<A>(a)...); })                          \
                    return law.name(std::forward<A>(a)...);                                          \
                else                                                                                 \
                    return name(law, std::forward<A>(a)...);                                         \
            }                                                                                        \
        };
// Opérations en un point de la loi : l'argument est un value_type (listes entre
// accolades et conversions acceptées, comme pour une fonction ordinaire).
#define WSL_RANDOM_POINT_OPERATION(name)                                                             \
        struct name##_fn                                                                             \
        {                                                                                            \
            template <class L>                                                                       \
                requires requires(L const& law, typename L::value_type const& x) { law.name(x); }    \
                         || requires(L const& law, typename L::value_type const& x) { name(law, x); } \
            constexpr decltype(auto) operator()(L const& law, typename L::value_type const& x) const \
            {                                                                                        \
                if constexpr (requires { law.name(x); })                                             \
                    return law.name(x);                                                              \
                else                                                                                 \
                    return name(law, x);                                                             \
            }                                                                                        \
        };
        WSL_RANDOM_POINT_OPERATION(pdf)
        WSL_RANDOM_POINT_OPERATION(log_pdf)
        WSL_RANDOM_POINT_OPERATION(pmf)
        WSL_RANDOM_POINT_OPERATION(log_pmf)
        WSL_RANDOM_POINT_OPERATION(cdf)
        WSL_RANDOM_POINT_OPERATION(ccdf)
#undef WSL_RANDOM_POINT_OPERATION
        WSL_RANDOM_OPERATION(quantile)
        WSL_RANDOM_OPERATION(cquantile)
        WSL_RANDOM_OPERATION(mean)
        WSL_RANDOM_OPERATION(variance)
        WSL_RANDOM_OPERATION(skewness)
        WSL_RANDOM_OPERATION(excess_kurtosis)
        WSL_RANDOM_OPERATION(entropy)
        WSL_RANDOM_OPERATION(median)
        WSL_RANDOM_OPERATION(mode)
        WSL_RANDOM_OPERATION(support)
#undef WSL_RANDOM_OPERATION

        // sample(X, g) ; sample(X, g, fast) : la variante rapide de la loi, sinon le tirage par défaut.
        struct sample_fn
        {
            template <class L, std::uniform_random_bit_generator G>
                requires requires(L const& law, G& g) { law.sample(g); } || requires(L const& law, G& g) { sample(law, g); }
            constexpr auto operator()(L const& law, G& g) const
            {
                if constexpr (requires { law.sample(g); })
                    return law.sample(g);
                else
                    return sample(law, g);
            }
            template <class L, std::uniform_random_bit_generator G>
                requires requires(sample_fn f, L const& law, G& g) { f(law, g); }
            constexpr auto operator()(L const& law, G& g, fast_t f) const
            {
                if constexpr (requires { law.sample(g, f); })
                    return law.sample(g, f);
                else if constexpr (requires { sample(law, g, f); })
                    return sample(law, g, f);
                else
                    return (*this)(law, g);
            }
        };

        // densité dP/dμ : pdf (Lebesgue) ou pmf (comptage), ou le membre density de la loi.
        struct density_fn
        {
            template <class L>
                requires requires(L const& law, typename L::value_type const& x) { law.density(x); }
                         || requires(L const& law, typename L::value_type const& x) { pdf_fn{}(law, x); }
                         || requires(L const& law, typename L::value_type const& x) { pmf_fn{}(law, x); }
            constexpr auto operator()(L const& law, typename L::value_type const& x) const
            {
                if constexpr (requires { law.density(x); })
                    return law.density(x);
                else if constexpr (requires { pdf_fn{}(law, x); })
                    return pdf_fn{}(law, x);
                else
                    return pmf_fn{}(law, x);
            }
        };
        struct log_density_fn
        {
            template <class L>
                requires requires(L const& law, typename L::value_type const& x) { law.log_density(x); }
                         || requires(L const& law, typename L::value_type const& x) { log_pdf_fn{}(law, x); }
                         || requires(L const& law, typename L::value_type const& x) { log_pmf_fn{}(law, x); }
            constexpr auto operator()(L const& law, typename L::value_type const& x) const
            {
                if constexpr (requires { law.log_density(x); })
                    return law.log_density(x);
                else if constexpr (requires { log_pdf_fn{}(law, x); })
                    return log_pdf_fn{}(law, x);
                else
                    return log_pmf_fn{}(law, x);
            }
        };
    }

    // Objets-fonctions : X.pdf(x), wsl::random::mean(X)… (espace inline :
    // pas de conflit avec les fonctions amies des lois).
    inline namespace operations
    {
        inline constexpr detail::cpo::pdf_fn pdf{};
        inline constexpr detail::cpo::log_pdf_fn log_pdf{};
        inline constexpr detail::cpo::pmf_fn pmf{};
        inline constexpr detail::cpo::log_pmf_fn log_pmf{};
        inline constexpr detail::cpo::density_fn density{};
        inline constexpr detail::cpo::log_density_fn log_density{};
        inline constexpr detail::cpo::cdf_fn cdf{};
        inline constexpr detail::cpo::ccdf_fn ccdf{};
        inline constexpr detail::cpo::quantile_fn quantile{};
        inline constexpr detail::cpo::cquantile_fn cquantile{};
        inline constexpr detail::cpo::sample_fn sample{};
        inline constexpr detail::cpo::mean_fn mean{};
        inline constexpr detail::cpo::variance_fn variance{};
        inline constexpr detail::cpo::skewness_fn skewness{};
        inline constexpr detail::cpo::excess_kurtosis_fn excess_kurtosis{};
        inline constexpr detail::cpo::entropy_fn entropy{};
        inline constexpr detail::cpo::median_fn median{};
        inline constexpr detail::cpo::mode_fn mode{};
        inline constexpr detail::cpo::support_fn support{};
    }

    // Support, quand un intervalle (ou un pavé) a un sens : lois numériques et
    // vectorielles ; les objets combinatoires (permutations, graphes…) n'en ont pas.
    template <class L>
    concept has_support = probability_law<L> && requires(L const& law) {
        { wsl::random::support(law) } -> std::same_as<interval<typename L::value_type>>;
    };

    // Constantes.
    template <class L> concept has_mean = requires(L const& law) { wsl::random::mean(law); };
    template <class L> concept has_variance = requires(L const& law) { wsl::random::variance(law); };
    template <class L> concept has_skewness = requires(L const& law) { wsl::random::skewness(law); };
    template <class L> concept has_excess_kurtosis = requires(L const& law) { wsl::random::excess_kurtosis(law); };
    template <class L> concept has_entropy = requires(L const& law) { wsl::random::entropy(law); };
    template <class L> concept has_median = requires(L const& law) { wsl::random::median(law); };
    template <class L> concept has_mode = requires(L const& law) { wsl::random::mode(law); };

    // Opérations.
    template <class L> concept has_pdf = requires(L const& law, typename L::value_type x) { wsl::random::pdf(law, x); };
    template <class L> concept has_pmf = requires(L const& law, typename L::value_type x) { wsl::random::pmf(law, x); };
    template <class L> concept has_density = requires(L const& law, typename L::value_type x) { wsl::random::density(law, x); };
    template <class L> concept has_cdf = requires(L const& law, typename L::value_type x) { wsl::random::cdf(law, x); };
    template <class L> concept has_ccdf = requires(L const& law, typename L::value_type x) { wsl::random::ccdf(law, x); };
    template <class L> concept has_quantile = requires(L const& law) { wsl::random::quantile(law, 0.5); };
    template <class L> concept has_cquantile = requires(L const& law) { wsl::random::cquantile(law, 0.5); };

    template <class L, class G>
    concept samplable = std::uniform_random_bit_generator<G> && requires(L const& law, G& g) { wsl::random::sample(law, g); };

    // Mesure de référence : déclarée par la loi, sinon déduite (pdf : Lebesgue, pmf : comptage).
    namespace detail
    {
        template <class L>
        consteval auto reference_measure_of()
        {
            if constexpr (requires { typename L::reference_measure; })
                return typename L::reference_measure{};
            else if constexpr (has_pmf<L>)
                return counting{};
            else
                return lebesgue{};
        }
    }
    template <class L>
    using reference_measure_t = decltype(detail::reference_measure_of<L>());

    // -----------------------------------------------------------------------
    // Base commune des lois (paramètre this explicite) : X(g), X(g, fast), X(n, g),
    // X(n, g, fast) et X.density(x), X.log_density(x), à partir de X.sample(g) et de
    // X.pdf / X.pmf. X(n, g) rend exactement les n valeurs, et consomme exactement les
    // mots du moteur, de n appels X(g) successifs.
    // -----------------------------------------------------------------------
    struct law_interface
    {
        template <class Self, std::uniform_random_bit_generator G>
        constexpr auto operator()(this Self const& self, G& g)
        {
            return wsl::random::sample(self, g);
        }
        template <class Self, std::uniform_random_bit_generator G>
        constexpr auto operator()(this Self const& self, G& g, fast_t f)
        {
            return wsl::random::sample(self, g, f);
        }
        template <class Self, std::uniform_random_bit_generator G>
        constexpr auto operator()(this Self const& self, std::size_t n, G& g)
        {
            std::vector<std::remove_cvref_t<decltype(wsl::random::sample(self, g))>> out;
            out.reserve(n);
            for (std::size_t i = 0; i < n; ++i)
                out.push_back(wsl::random::sample(self, g));
            return out;
        }
        template <class Self, std::uniform_random_bit_generator G>
        constexpr auto operator()(this Self const& self, std::size_t n, G& g, fast_t f)
        {
            std::vector<std::remove_cvref_t<decltype(wsl::random::sample(self, g, f))>> out;
            out.reserve(n);
            for (std::size_t i = 0; i < n; ++i)
                out.push_back(wsl::random::sample(self, g, f));
            return out;
        }
        template <class Self, class X>
            requires requires(Self const& self, X const& x) { self.pdf(x); } || requires(Self const& self, X const& x) { self.pmf(x); }
        constexpr auto density(this Self const& self, X const& x)
        {
            if constexpr (requires { self.pdf(x); })
                return self.pdf(x);
            else
                return self.pmf(x);
        }
        template <class Self, class X>
            requires requires(Self const& self, X const& x) { self.log_pdf(x); } || requires(Self const& self, X const& x) { self.log_pmf(x); }
        constexpr auto log_density(this Self const& self, X const& x)
        {
            if constexpr (requires { self.log_pdf(x); })
                return self.log_pdf(x);
            else
                return self.log_pmf(x);
        }
        friend constexpr bool operator==(law_interface const&, law_interface const&) noexcept = default;
    };

    // -----------------------------------------------------------------------
    // Adaptateur : une loi scalaire vue comme une distribution de <random>
    // (exigence RandomNumberDistribution, sans les flux << et >>).
    // -----------------------------------------------------------------------
    template <probability_law L>
        requires std::is_arithmetic_v<typename L::value_type> && has_support<L>
    class distribution
    {
    public:
        using result_type = typename L::value_type;

        struct param_type
        {
            using distribution_type = distribution;
            L law;
            friend constexpr bool operator==(param_type const&, param_type const&) = default;
        };

        constexpr explicit distribution(L law) : param_{law} {}
        constexpr explicit distribution(param_type const& p) : param_(p) {}

        constexpr void reset() noexcept {}

        template <std::uniform_random_bit_generator G>
        constexpr result_type operator()(G& g) const
        {
            return sample(param_.law, g);
        }

        template <std::uniform_random_bit_generator G>
        constexpr result_type operator()(G& g, param_type const& p) const
        {
            return sample(p.law, g);
        }

        constexpr param_type param() const { return param_; }
        constexpr void param(param_type const& p) { param_ = p; }
        constexpr L const& law() const noexcept { return param_.law; }

        constexpr result_type min() const { return param_.law.support().lower; }
        constexpr result_type max() const { return param_.law.support().upper; }

        friend constexpr bool operator==(distribution const&, distribution const&) = default;

    private:
        param_type param_;
    };

    template <probability_law L>
        requires std::is_arithmetic_v<typename L::value_type> && has_support<L>
    constexpr distribution<L> as_distribution(L law)
    {
        return distribution<L>(law);
    }
}
