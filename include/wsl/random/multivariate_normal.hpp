#pragma once

// Loi normale multivariée N(μ, Σ), Σ symétrique définie positive.
//
//   wsl::random::multivariate_normal<double, 2> X({.mean = {1, 2}, .covariance = {{1, 0.5}, {0.5, 2}}});
//   X.pdf({1.5, 1.0}); X.entropy(); X.marginal(1); X(g);
//
// Σ = L Lᵀ (Cholesky, à la construction ; std::domain_error si Σ n'est pas symétrique
// définie positive). Log-densité −½ |z|² − Σ log Lᵢᵢ − (n/2) log 2π, z = L⁻¹(x − μ) par
// substitution, sommes en double-double ; densité : un arrondi de l'exponentielle.
// Entropie (n/2)(1 + log 2π) + Σ log Lᵢᵢ ; marginale i : N(μᵢ, √Σᵢᵢ).
//
// Tirage (approché, bit-exact) : μ + L z, z de composantes normales centrées réduites
// tirées exactement (wsl::random::fast : ziggourat), dans l'ordre ; chaque composante arrondie
// une fois.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/special.hpp>

#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    template <law_real T = double, std::size_t N = dynamic>
    class multivariate_normal : public law_interface
    {
    public:
        using value_type = vector_of<T, N>;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            vector_of<T, N> mean;
            square_matrix<T, N> covariance;
        };

        constexpr multivariate_normal(parameters p) : mu_(p.mean), sigma_(p.covariance), l_(make_matrix(p.mean.size()))
        {
            if (mu_.size() < 1 || sigma_.size() != mu_.size())
                throw std::domain_error("wsl::random::multivariate_normal : moyenne et covariance de dimensions différentes");
            for (T const m : mu_)
                if (!detail::is_finite(m))
                    throw std::domain_error("wsl::random::multivariate_normal : moyenne non finie");
            if (!detail::cholesky(sigma_, l_))
                throw std::domain_error("wsl::random::multivariate_normal : covariance non symétrique définie positive");
            detail::sd::double_double s{0.0, 0.0};
            for (std::size_t i = 0; i < dimension(); ++i)
                s = detail::sum(s, detail::log_precise(l_(i, i)));
            log_det_half_ = s;
        }

        constexpr std::size_t dimension() const noexcept { return mu_.size(); }
        constexpr value_type const& mean() const noexcept { return mu_; }
        constexpr square_matrix<T, N> const& covariance() const noexcept { return sigma_; }
        constexpr square_matrix<double, N> const& cholesky_factor() const noexcept { return l_; }   // L, Σ = L Lᵀ
        constexpr value_type variance() const noexcept
        {
            value_type v = detail::make_vector<T, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                v[i] = sigma_(i, i);
            return v;
        }
        constexpr interval<value_type> support() const noexcept
        {
            return {detail::make_vector<T, N>(dimension(), -std::numeric_limits<T>::infinity()), detail::make_vector<T, N>(dimension(), std::numeric_limits<T>::infinity())};
        }
        // (n/2)(1 + log 2π) + ½ log det Σ
        constexpr double entropy() const noexcept
        {
            auto const h = detail::sum(detail::mul({double(dimension()), 0.0}, {detail::normal_entropy_hi, detail::normal_entropy_lo}), log_det_half_);
            return h.hi + h.lo;
        }
        constexpr normal<T> marginal(std::size_t i) const
        {
            return normal<T>({.mean = mu_[i], .sd = static_cast<T>(wsl::special::sqrt(double(sigma_(i, i))))});
        }
        // (x − μ)ᵀ Σ⁻¹ (x − μ)
        constexpr double mahalanobis_squared(value_type const& x) const noexcept
        {
            auto const q = quadratic(x);
            return q.hi + q.lo;
        }

        constexpr double log_pdf(value_type const& x) const noexcept
        {
            auto const& X = *this;
            auto const l = X.log_density_dd(x);
            return l.hi + l.lo;
        }
        constexpr double pdf(value_type const& x) const noexcept
        {
            auto const& X = *this;
            auto const l = X.log_density_dd(x);
            if (!(l.hi > -1500))
                return l.hi == l.hi ? 0.0 : l.hi;
            if (l.hi > 1400)
                return std::numeric_limits<double>::infinity();
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& X = *this;
            auto z = detail::make_vector<double, N>(X.dimension());
            normal<double> const standard({.mean = 0.0, .sd = 1.0});
            for (auto& v : z)
                v = wsl::random::sample(standard, g);
            return X.transform(z);
        }
        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g, fast_t) const
        {
            auto const& X = *this;
            auto z = detail::make_vector<double, N>(X.dimension());
            for (auto& v : z)
                v = detail::ziggurat_normal(g);
            return X.transform(z);
        }

        friend constexpr bool operator==(multivariate_normal const& a, multivariate_normal const& b) noexcept { return a.mu_ == b.mu_ && a.sigma_ == b.sigma_; }

    private:
        static constexpr square_matrix<double, N> make_matrix(std::size_t n)
        {
            if constexpr (N == dynamic)
                return square_matrix<double, N>(n);
            else
                return square_matrix<double, N>();
        }

        // μ + L z, chaque composante arrondie une fois.
        constexpr value_type transform(vector_of<double, N> const& z) const noexcept
        {
            value_type x = detail::make_vector<T, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
            {
                auto const e = detail::dot(i + 1, [&](std::size_t k) { return std::pair{l_(i, k), z[k]}; });
                auto const t = detail::sum(e, {double(mu_[i]), 0.0});
                x[i] = static_cast<T>(t.hi + t.lo);
            }
            return x;
        }

        // |L⁻¹(x − μ)|², en double-double.
        constexpr detail::sd::double_double quadratic(value_type const& x) const noexcept
        {
            auto d = detail::make_vector<double, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                d[i] = double(x[i]) - double(mu_[i]);
            auto const z = detail::forward_substitution(l_, d);
            return detail::dot(dimension(), [&](std::size_t k) { return std::pair{z[k], z[k]}; });
        }

        constexpr detail::sd::double_double log_density_dd(value_type const& x) const noexcept
        {
            if (x.size() != dimension())
                return {-std::numeric_limits<double>::infinity(), 0.0};
            for (T const v : x)
                if (!detail::is_finite(v))
                    return {v == v ? -std::numeric_limits<double>::infinity() : double(v), 0.0};
            auto const q = quadratic(x);
            detail::sd::double_double l{detail::sd::fmul(q.hi, -0.5), detail::sd::fmul(q.lo, -0.5)};
            l = detail::sum(l, {-log_det_half_.hi, -log_det_half_.lo});
            auto const c = detail::mul({double(dimension()), 0.0}, {detail::half_log_2pi_hi, detail::half_log_2pi_lo});
            return detail::sum(l, {-c.hi, -c.lo});
        }

        vector_of<T, N> mu_;
        square_matrix<T, N> sigma_;
        square_matrix<double, N> l_;
        detail::sd::double_double log_det_half_{};
    };
}
