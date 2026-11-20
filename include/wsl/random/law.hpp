#pragma once

// Le concept de loi de probabilité (voir docs/random.md §4.3).
//
// Une loi est un type qui déclare :
//   - value_type : le type des valeurs tirées ;
//   - sampling   : l'exactitude de son tirage par défaut (exact ou approximate) ;
//   - support()  : l'intervalle [lower, upper] des valeurs possibles, quand il a un
//                  sens (has_support : pas pour les objets combinatoires) ;
//   - ses constantes, en membres, quand elles existent : mean(), variance(),
//     skewness(), excess_kurtosis(), entropy(), median(), mode() ;
//   - ses opérations, en fonctions libres trouvées par ADL (appel sans
//     qualification, comme swap) : sample(X, g), sample(X, g, wsl::random::fast),
//     pdf / pmf (densité, fonction de masse), log_pdf / log_pmf, cdf, ccdf
//     (survie, 1 − cdf calculée sans annulation), quantile, cquantile (quantile de
//     survie : x tel que P(X > x) = q, sans passer par 1 − q).
//
// Ce qui n'existe pas mathématiquement — l'espérance d'une loi de Cauchy — ou
// n'est pas implémenté n'est pas déclaré : l'appel ne compile pas, et les
// concepts has_mean, has_cdf… permettent de le tester.

#include <wsl/random/uniform.hpp>

#include <concepts>
#include <random>
#include <type_traits>

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

    // Support, quand un intervalle (ou un pavé) a un sens : lois numériques et
    // vectorielles ; les objets combinatoires (permutations, graphes…) n'en ont pas.
    template <class L>
    concept has_support = probability_law<L> && requires(L const& law) {
        { law.support() } -> std::same_as<interval<typename L::value_type>>;
    };

    // Constantes.
    template <class L> concept has_mean = requires(L const& law) { law.mean(); };
    template <class L> concept has_variance = requires(L const& law) { law.variance(); };
    template <class L> concept has_skewness = requires(L const& law) { law.skewness(); };
    template <class L> concept has_excess_kurtosis = requires(L const& law) { law.excess_kurtosis(); };
    template <class L> concept has_entropy = requires(L const& law) { law.entropy(); };
    template <class L> concept has_median = requires(L const& law) { law.median(); };
    template <class L> concept has_mode = requires(L const& law) { law.mode(); };

    // Opérations (ADL).
    template <class L> concept has_pdf = requires(L const& law, typename L::value_type x) { pdf(law, x); };
    template <class L> concept has_pmf = requires(L const& law, typename L::value_type x) { pmf(law, x); };
    template <class L> concept has_cdf = requires(L const& law, typename L::value_type x) { cdf(law, x); };
    template <class L> concept has_ccdf = requires(L const& law, typename L::value_type x) { ccdf(law, x); };
    template <class L> concept has_quantile = requires(L const& law) { quantile(law, 0.5); };
    template <class L> concept has_cquantile = requires(L const& law) { cquantile(law, 0.5); };

    template <class L, class G>
    concept samplable = std::uniform_random_bit_generator<G> && requires(L const& law, G& g) { sample(law, g); };

    // Tirage rapide demandé à une loi qui n'en a pas de propre : le tirage par
    // défaut. (Une loi qui en a un le déclare ; il est préféré, étant plus
    // spécialisé.)
    template <probability_law L, std::uniform_random_bit_generator G>
        requires samplable<L, G>
    constexpr auto sample(L const& law, G& g, fast_t)
    {
        return sample(law, g);
    }

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
