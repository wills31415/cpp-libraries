#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Loi du χ² à k > 0 degrés de liberté : la loi gamma de forme k/2 et d'échelle 2
// (somme des carrés de k normales centrées réduites, pour k entier).
//
//   wsl::probability::laws::chi_squared X({.degrees_of_freedom = 3});
//   X.mean(); X.cdf(7.81); X.quantile(0.95); X(g);
//
// Tout est délégué à wsl::probability::laws::gamma (k/2 et l'échelle 2 sont exacts) : mêmes méthodes,
// même précision.

#include <wsl/probability/laws/gamma.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_real T = double>
    class chi_squared : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T degrees_of_freedom;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.degrees_of_freedom = degrees_of_freedom()}; }

        constexpr chi_squared(parameters p) : k_(p.degrees_of_freedom), g_(make(p.degrees_of_freedom)) {}

        constexpr T degrees_of_freedom() const noexcept { return k_; }
        constexpr interval<T> support() const noexcept { return g_.support(); }
        constexpr T mean() const noexcept { return k_; }
        constexpr T variance() const noexcept { return static_cast<T>(detail::sd::fmul(2.0, double(k_))); }
        constexpr T mode() const noexcept { return g_.mode(); }
        constexpr T skewness() const noexcept { return g_.skewness(); }
        constexpr T excess_kurtosis() const noexcept { return g_.excess_kurtosis(); }
        constexpr T entropy() const noexcept { return g_.entropy(); }

        constexpr T pdf(T x) const noexcept { auto const& X = *this; return wsl::probability::pdf(X.g_, x); }
        constexpr T log_pdf(T x) const noexcept { auto const& X = *this; return wsl::probability::log_pdf(X.g_, x); }
        constexpr T cdf(T x) const noexcept { auto const& X = *this; return wsl::probability::cdf(X.g_, x); }
        constexpr T ccdf(T x) const noexcept { auto const& X = *this; return wsl::probability::ccdf(X.g_, x); }
        constexpr T quantile(double p) const { auto const& X = *this; return wsl::probability::quantile(X.g_, p); }
        constexpr T cquantile(double q) const { auto const& X = *this; return wsl::probability::cquantile(X.g_, q); }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            return detail::native_sample(X.g_, g);
        }

        friend constexpr bool operator==(chi_squared const& a, chi_squared const& b) noexcept { return a.k_ == b.k_; }

    private:
        static constexpr gamma<T> make(T k)
        {
            if (!(k > 0) || !detail::is_finite(k))
                throw std::domain_error("wsl::probability::laws::chi_squared : il faut degrees_of_freedom > 0, fini");
            return gamma<T>({.scale = T(2), .shape = static_cast<T>(double(k) / 2)});
        }

        T k_;
        gamma<T> g_;
    };
}
