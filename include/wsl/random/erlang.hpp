#pragma once

// Loi d'Erlang : somme de k exponentielles indépendantes de taux λ, gamma(k, λ) à forme
// entière.
//
//   wsl::random::erlang E({.rate = 2, .shape = 3});   E.cdf(1.0); E.quantile(0.5); E(g);
//
// Fonctions : celles de wsl::random::gamma (sous 0,5 ulp). Tirage **exact** : somme paresseuse de k
// exponentielles exactes, arrondie correctement (detail/lazy_sum.hpp) ; coût O(k).

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/lazy_sum.hpp>
#include <wsl/random/gamma.hpp>
#include <wsl/random/law.hpp>

#include <cstdint>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    template <law_real T = double>
    class erlang : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::exact;

        // Intensité .rate ou échelle .scale, comme gamma (exactes toutes deux).
        struct parameters
        {
            std::optional<ratio> rate, scale;
            std::int64_t shape;
        };
        using param_type = parameters;
        constexpr param_type param() const
        {
            auto const [r, s] = detail::parameters_of(g_.exact_rate());
            return {.rate = r, .scale = s, .shape = k_};
        }

        constexpr erlang(parameters p) : k_(p.shape), g_(check(p))
        {
        }

        constexpr std::int64_t shape() const noexcept { return k_; }
        constexpr T rate() const noexcept { return g_.rate(); }
        constexpr T scale() const noexcept { return g_.scale(); }
        constexpr gamma<T> const& as_gamma() const noexcept { return g_; }
        constexpr interval<T> support() const noexcept { return g_.support(); }
        constexpr T mean() const noexcept { return g_.mean(); }
        constexpr T variance() const noexcept { return g_.variance(); }
        constexpr T mode() const noexcept { return g_.mode(); }
        constexpr T skewness() const noexcept { return g_.skewness(); }
        constexpr T excess_kurtosis() const noexcept { return g_.excess_kurtosis(); }
        constexpr T entropy() const noexcept { return g_.entropy(); }

        constexpr T pdf(T x) const noexcept { auto const& E = *this; return wsl::random::pdf(E.g_, x); }
        constexpr T log_pdf(T x) const noexcept { auto const& E = *this; return wsl::random::log_pdf(E.g_, x); }
        constexpr T cdf(T x) const noexcept { auto const& E = *this; return wsl::random::cdf(E.g_, x); }
        constexpr T ccdf(T x) const noexcept { auto const& E = *this; return wsl::random::ccdf(E.g_, x); }
        constexpr T quantile(double p) const { auto const& E = *this; return wsl::random::quantile(E.g_, p); }
        constexpr T cquantile(double q) const { auto const& E = *this; return wsl::random::cquantile(E.g_, q); }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& E = *this;
            auto const [r, sc] = detail::parameters_of(E.g_.exact_rate());
            exponential<double> const e({.rate = r, .scale = sc});
            double const r0 = e.exact_rate().as_double();
            detail::lazy_sum s(std::vector<double>{r0 != 0 ? r0 : e.exact_rate().sampling_factors().divisor});
            for (std::int64_t i = 0; i < E.k_; ++i)
                s.add(detail::lazy_sample(e, g), g);
            return s.template rounded<T>(g, {}, {}, 0, detail::to_dyadic(1.0));
        }
        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g, fast_t) const
        {
            auto const& E = *this;
            return wsl::random::sample(E.g_, g);
        }

        friend constexpr bool operator==(erlang const& a, erlang const& b) noexcept { return a.k_ == b.k_ && a.g_ == b.g_; }

    private:
        static constexpr gamma<T> check(parameters const& p)
        {
            if (p.shape <= 0 || p.shape > (std::int64_t{1} << 53))
                throw std::domain_error("wsl::random::erlang : il faut 1 <= shape <= 2^53");
            return gamma<T>({.rate = p.rate, .scale = p.scale, .shape = static_cast<T>(p.shape)});
        }

        std::int64_t k_;
        gamma<T> g_;
    };
}
