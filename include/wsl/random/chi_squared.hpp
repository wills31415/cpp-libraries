#pragma once

// Loi du χ² à k > 0 degrés de liberté : la loi gamma de forme k/2 et d'échelle 2
// (somme des carrés de k normales centrées réduites, pour k entier).
//
//   wsl::chi_squared X({.degrees_of_freedom = 3});
//   X.mean(); cdf(X, 7.81); quantile(X, 0.95); sample(X, g);
//
// Tout est délégué à wsl::gamma (k/2 et l'échelle 2 sont exacts) : mêmes méthodes,
// même précision.

#include <wsl/random/gamma.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl
{
    template <law_real T = double>
    class chi_squared
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T degrees_of_freedom;
        };

        constexpr chi_squared(parameters p) : k_(p.degrees_of_freedom), g_(make(p.degrees_of_freedom)) {}

        constexpr T degrees_of_freedom() const noexcept { return k_; }
        constexpr interval<T> support() const noexcept { return g_.support(); }
        constexpr T mean() const noexcept { return k_; }
        constexpr T variance() const noexcept { return static_cast<T>(detail::sd::fmul(2.0, double(k_))); }
        constexpr T mode() const noexcept { return g_.mode(); }
        constexpr T skewness() const noexcept { return g_.skewness(); }
        constexpr T excess_kurtosis() const noexcept { return g_.excess_kurtosis(); }
        constexpr T entropy() const noexcept { return g_.entropy(); }

        friend constexpr T pdf(chi_squared const& X, T x) noexcept { return pdf(X.g_, x); }
        friend constexpr T log_pdf(chi_squared const& X, T x) noexcept { return log_pdf(X.g_, x); }
        friend constexpr T cdf(chi_squared const& X, T x) noexcept { return cdf(X.g_, x); }
        friend constexpr T ccdf(chi_squared const& X, T x) noexcept { return ccdf(X.g_, x); }
        friend constexpr T quantile(chi_squared const& X, double p) { return quantile(X.g_, p); }
        friend constexpr T cquantile(chi_squared const& X, double q) { return cquantile(X.g_, q); }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(chi_squared const& X, G& g)
        {
            return sample(X.g_, g);
        }

        friend constexpr bool operator==(chi_squared const& a, chi_squared const& b) noexcept { return a.k_ == b.k_; }

    private:
        static constexpr gamma<T> make(T k)
        {
            if (!(k > 0) || !detail::is_finite(k))
                throw std::domain_error("wsl::chi_squared : il faut degrees_of_freedom > 0, fini");
            return gamma<T>({.shape = static_cast<T>(double(k) / 2), .scale = T(2)});
        }

        T k_;
        gamma<T> g_;
    };
}
