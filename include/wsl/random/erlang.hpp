#pragma once

// Loi d'Erlang : somme de k exponentielles indépendantes de taux λ, gamma(k, λ) à forme
// entière.
//
//   wsl::random::erlang E({.shape = 3, .rate = 2});   cdf(E, 1.0); quantile(E, 0.5); sample(E, g);
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
    class erlang
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::int64_t shape;
            T rate;
        };

        constexpr erlang(parameters p) : k_(p.shape), g_(check(p))
        {
        }

        constexpr std::int64_t shape() const noexcept { return k_; }
        constexpr T rate() const noexcept { return g_.rate(); }
        constexpr gamma<T> const& as_gamma() const noexcept { return g_; }
        constexpr interval<T> support() const noexcept { return g_.support(); }
        constexpr T mean() const noexcept { return g_.mean(); }
        constexpr T variance() const noexcept { return g_.variance(); }
        constexpr T mode() const noexcept { return g_.mode(); }
        constexpr T skewness() const noexcept { return g_.skewness(); }
        constexpr T excess_kurtosis() const noexcept { return g_.excess_kurtosis(); }
        constexpr T entropy() const noexcept { return g_.entropy(); }

        friend constexpr T pdf(erlang const& E, T x) noexcept { return pdf(E.g_, x); }
        friend constexpr T log_pdf(erlang const& E, T x) noexcept { return log_pdf(E.g_, x); }
        friend constexpr T cdf(erlang const& E, T x) noexcept { return cdf(E.g_, x); }
        friend constexpr T ccdf(erlang const& E, T x) noexcept { return ccdf(E.g_, x); }
        friend constexpr T quantile(erlang const& E, double p) { return quantile(E.g_, p); }
        friend constexpr T cquantile(erlang const& E, double q) { return cquantile(E.g_, q); }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(erlang const& E, G& g)
        {
            detail::lazy_sum s(std::vector<double>{double(E.g_.rate())});
            for (std::int64_t i = 0; i < E.k_; ++i)
                s.add(detail::lazy_sample(exponential<double>({.rate = double(E.g_.rate())}), g), g);
            return s.template rounded<T>(g, {}, {}, 0, detail::to_dyadic(1.0));
        }
        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(erlang const& E, G& g, fast_t)
        {
            return sample(E.g_, g);
        }

        friend constexpr bool operator==(erlang const& a, erlang const& b) noexcept { return a.k_ == b.k_ && a.g_ == b.g_; }

    private:
        static constexpr gamma<T> check(parameters const& p)
        {
            if (p.shape <= 0 || p.shape > (std::int64_t{1} << 53))
                throw std::domain_error("wsl::random::erlang : il faut 1 <= shape <= 2^53");
            return gamma<T>(typename gamma<T>::shape_rate{static_cast<T>(p.shape), p.rate});
        }

        std::int64_t k_;
        gamma<T> g_;
    };
}
