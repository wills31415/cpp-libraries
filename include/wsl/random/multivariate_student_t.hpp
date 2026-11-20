#pragma once

// Loi de Student multivariée t_ν(μ, Σ) : μ + L z √(ν/V), Σ = L Lᵀ, z normale centrée
// réduite de dimension n, V de loi χ²(ν) indépendante.
//
//   wsl::random::multivariate_student_t<double, 2> X({.location = {0, 1}, .scale = {{1, 0.3}, {0.3, 2}},
//                                                     .degrees_of_freedom = 4});
//   X.pdf({0.5, 0.5}); X.covariance(); X.marginal(0); X(g);
//
// Log-densité lgamma((ν + n)/2) − lgamma(ν/2) − (n/2) log(νπ) − ½ log det Σ
// − ((ν + n)/2) log1p(q/ν), q = (x − μ)ᵀ Σ⁻¹ (x − μ), en double-double. Moyenne μ (ν > 1),
// covariance ν/(ν − 2) Σ (ν > 2) ; NaN, ou +∞, sinon. Entropie ½ log det Σ
// + log((νπ)^{n/2} Γ(ν/2)/Γ((ν + n)/2)) + ((ν + n)/2)(ψ((ν + n)/2) − ψ(ν/2)) ;
// marginale i : √Σᵢᵢ T + μᵢ, T de Student (wsl::random::affine_image).
//
// Tirage (approché, bit-exact) : z (normales exactes), puis V (gamma), dans cet ordre.

#include <wsl/random/affine.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/gamma.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/random/student_t.hpp>
#include <wsl/special.hpp>

#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    template <law_real T = double, std::size_t N = dynamic>
    class multivariate_student_t : public law_interface
    {
    public:
        using value_type = vector_of<T, N>;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            vector_of<T, N> location;
            square_matrix<T, N> scale;
            T degrees_of_freedom;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.location = location(), .scale = scale(), .degrees_of_freedom = degrees_of_freedom()}; }

        constexpr multivariate_student_t(parameters p)
            : nu_(p.degrees_of_freedom), mu_(p.location), sigma_(p.scale), l_(detail::make_matrix<double, N>(p.location.size())),
              chi_(typename gamma<double>::shape_scale{.scale = 2.0, .shape = detail::sd::fmul(double(p.degrees_of_freedom), 0.5)})
        {
            if (!(nu_ > 0) || !detail::is_finite(nu_))
                throw std::domain_error("wsl::random::multivariate_student_t : il faut degrees_of_freedom > 0, fini");
            if (mu_.size() < 1 || sigma_.size() != mu_.size())
                throw std::domain_error("wsl::random::multivariate_student_t : position et échelle de dimensions différentes");
            for (T const m : mu_)
                if (!detail::is_finite(m))
                    throw std::domain_error("wsl::random::multivariate_student_t : position non finie");
            if (!detail::cholesky(sigma_, l_))
                throw std::domain_error("wsl::random::multivariate_student_t : échelle non symétrique définie positive");
            namespace fp = detail::sd;
            double const n = double(dimension()), nu = nu_;
            double const h = fp::fmul(nu, 0.5);
            auto const [s_hi, s_lo] = fp::two_sum(h, fp::fmul(n, 0.5));     // (ν + n)/2
            half_sum_ = {s_hi, s_lo};
            // lgamma((ν + n)/2) − lgamma(ν/2) − (n/2) log(νπ) − ½ log det Σ
            auto c = detail::log_gamma(s_hi, s_lo);
            auto const lh = detail::log_gamma(h);
            c = detail::sum(c, {-lh.hi, -lh.lo});
            auto const lnp = detail::sum(detail::log_precise(nu), {fp::log_pi_hi, fp::log_pi_lo});
            c = detail::sum(c, {-fp::fmul(lnp.hi, fp::fmul(n, 0.5)), -fp::fmul(lnp.lo, fp::fmul(n, 0.5))});
            log_det_half_ = detail::half_log_determinant(l_);
            constant_ = detail::sum(c, {-log_det_half_.hi, -log_det_half_.lo});
        }

        constexpr std::size_t dimension() const noexcept { return mu_.size(); }
        constexpr T degrees_of_freedom() const noexcept { return nu_; }
        constexpr value_type const& location() const noexcept { return mu_; }
        constexpr square_matrix<T, N> const& scale() const noexcept { return sigma_; }
        constexpr interval<value_type> support() const noexcept
        {
            return {detail::make_vector<T, N>(dimension(), -std::numeric_limits<T>::infinity()), detail::make_vector<T, N>(dimension(), std::numeric_limits<T>::infinity())};
        }

        constexpr value_type mean() const noexcept
        {
            return nu_ > 1 ? mu_ : detail::make_vector<T, N>(dimension(), std::numeric_limits<T>::quiet_NaN());
        }
        // ν/(ν − 2) Σ
        constexpr square_matrix<T, N> covariance() const
        {
            auto c = detail::make_matrix<T, N>(dimension());
            double const nu = nu_;
            for (std::size_t i = 0; i < dimension(); ++i)
                for (std::size_t j = 0; j < dimension(); ++j)
                {
                    double v = std::numeric_limits<double>::quiet_NaN();
                    if (nu > 2)
                        v = double(sigma_(i, j)) == 0 ? 0.0 : (double(sigma_(i, j)) < 0 ? -1 : 1) *
                                detail::rounded(detail::normalized(nu, 0.0) * detail::normalized(detail::magnitude(double(sigma_(i, j))), 0.0) / (nu - 2));
                    else if (nu > 1)
                        v = i == j ? std::numeric_limits<double>::infinity() : std::numeric_limits<double>::quiet_NaN();
                    c(i, j) = static_cast<T>(v);
                }
            return c;
        }
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            double const nu = nu_;
            // −constante + ((ν + n)/2)(ψ((ν + n)/2) − ψ(ν/2)) (ψ par wsl::special)
            double const d = wsl::special::digamma(half_sum_.hi + half_sum_.lo) - wsl::special::digamma(fp::fmul(nu, 0.5));
            auto h = detail::sum({-constant_.hi, -constant_.lo}, detail::mul(half_sum_, {d, 0.0}));
            return static_cast<T>(h.hi + h.lo);
        }
        constexpr affine_image<student_t<T>> marginal(std::size_t i) const
        {
            return affine_image<student_t<T>>(student_t<T>({.degrees_of_freedom = nu_}), {.scale = wsl::special::sqrt(double(sigma_(i, i))), .shift = double(mu_[i])});
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
            normal<double> const standard({.mean = 0.0, .stddev = 1.0});
            for (auto& v : z)
                v = wsl::random::sample(standard, g);
            double const w = wsl::random::sample(X.chi_, g);
            double const f = wsl::special::sqrt(double(X.nu_) / w);         // √(ν/V)
            value_type x = detail::make_vector<T, N>(X.dimension());
            for (std::size_t i = 0; i < X.dimension(); ++i)
            {
                auto const e = detail::dot(i + 1, [&](std::size_t k) { return std::pair{X.l_(i, k), z[k]}; });
                auto const t = detail::sum(detail::mul(e, {f, 0.0}), {double(X.mu_[i]), 0.0});
                x[i] = static_cast<T>(t.hi + t.lo);
            }
            return x;
        }

        friend constexpr bool operator==(multivariate_student_t const& a, multivariate_student_t const& b) noexcept
        {
            return a.nu_ == b.nu_ && a.mu_ == b.mu_ && a.sigma_ == b.sigma_;
        }

    private:
        constexpr detail::sd::double_double log_density_dd(value_type const& x) const noexcept
        {
            if (x.size() != dimension())
                return {-std::numeric_limits<double>::infinity(), 0.0};
            for (T const v : x)
                if (!detail::is_finite(v))
                    return {v == v ? -std::numeric_limits<double>::infinity() : double(v), 0.0};
            auto d = detail::make_vector<double, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                d[i] = double(x[i]) - double(mu_[i]);
            auto const z = detail::forward_substitution(l_, d);
            auto const q = detail::dot(dimension(), [&](std::size_t k) { return std::pair{z[k], z[k]}; });
            auto const r = detail::sd::quotient_dd(q, {double(nu_), 0.0});   // q/ν
            auto const lp = detail::log1p_positive(r.hi, r.lo);
            auto const t = detail::mul(half_sum_, lp);
            return detail::sum(constant_, {-t.hi, -t.lo});
        }

        T nu_;
        vector_of<T, N> mu_;
        square_matrix<T, N> sigma_;
        square_matrix<double, N> l_;
        gamma<double> chi_;                                                 // χ²(ν) = gamma(ν/2, échelle 2)
        detail::sd::double_double half_sum_{}, log_det_half_{}, constant_{};
    };
}
