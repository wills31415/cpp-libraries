#pragma once

// Lois de Wishart W_p(V, n) et inverse-Wishart W⁻¹_p(Ψ, ν) sur les matrices p × p
// symétriques définies positives (n, ν > p − 1).
//
//   wsl::random::wishart<double, 2> W({.degrees_of_freedom = 5, .scale = {{1, 0.2}, {0.2, 2}}});
//   wsl::random::inverse_wishart<double, 2> I({.degrees_of_freedom = 6, .scale = {{1, 0.2}, {0.2, 2}}});
//   W.mean(); W.log_pdf(X); W.entropy(); W(g);
//
// Wishart : log-densité ((n − p − 1)/2) log det X − ½ tr(V⁻¹X) − (np/2) log 2
// − (n/2) log det V − log Γ_p(n/2), par Cholesky de X (−∞ si X n'est pas symétrique
// définie positive) et tr(V⁻¹X) = ‖L⁻¹ L_X‖²_F, en double-double ; moyenne nV,
// variances n(Vᵢⱼ² + VᵢᵢVⱼⱼ), mode (n − p − 1)V (n >= p + 1), entropie
// ((p + 1)/2) log det V + (p(p + 1)/2) log 2 + log Γ_p(n/2) − ((n − p − 1)/2) ψ_p(n/2)
// + np/2 ; marginale Wᵢᵢ : gamma(n/2, échelle 2Vᵢᵢ).
// Inverse-Wishart (X⁻¹ ~ W(Ψ⁻¹, ν)) : log-densité (ν/2) log det Ψ − (νp/2) log 2
// − log Γ_p(ν/2) − ((ν + p + 1)/2) log det X − ½ tr(Ψ X⁻¹) ; moyenne Ψ/(ν − p − 1)
// (ν > p + 1), mode Ψ/(ν + p + 1), entropie celle de W(Ψ⁻¹, ν) − (p + 1) E[log det W].
//
// Tirage (approché, bit-exact) : décomposition de Bartlett, A triangulaire inférieure,
// Aᵢᵢ = √χ²(n − i), Aᵢⱼ normales exactes (i > j), tirés ligne par ligne (Aᵢᵢ puis
// Aᵢ₀, …, Aᵢ,ᵢ₋₁) ; W = (LA)(LA)ᵀ, V = LLᵀ ; inverse-Wishart X = (LA⁻ᵀ)(LA⁻ᵀ)ᵀ, Ψ = LLᵀ.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/gamma.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/special.hpp>

#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace wsl::random
{
    namespace detail
    {
        // Bartlett : A triangulaire inférieure, Aᵢᵢ² ~ χ²(n − i), Aᵢⱼ ~ N(0, 1).
        template <std::size_t N, class G>
        constexpr square_matrix<double, N> bartlett(std::vector<gamma<double>> const& chi, G& g)
        {
            std::size_t const p = chi.size();
            auto a = make_matrix<double, N>(p);
            normal<double> const standard({.mean = 0.0, .sd = 1.0});
            for (std::size_t i = 0; i < p; ++i)
            {
                a(i, i) = wsl::special::sqrt(wsl::random::sample(chi[i], g));
                for (std::size_t j = 0; j < i; ++j)
                    a(i, j) = wsl::random::sample(standard, g);
            }
            return a;
        }

        // Parties communes : Cholesky de l'échelle, lois χ²(n − i).
        template <class T, std::size_t N>
        struct wishart_base
        {
            T n;
            square_matrix<T, N> v;
            square_matrix<double, N> l;
            std::vector<gamma<double>> chi;
            sd::double_double log_det_half{}, log_gamma_p{};

            constexpr wishart_base(T degrees, square_matrix<T, N> const& scale, char const* name)
                : n(degrees), v(scale), l(make_matrix<double, N>(scale.size()))
            {
                std::size_t const p = v.size();
                if (p < 1)
                    throw std::domain_error(name);
                if (!(double(n) > double(p) - 1) || !is_finite(n))
                    throw std::domain_error("wsl::random::wishart : il faut degrees_of_freedom > dimension − 1, fini");
                if (!cholesky(v, l))
                    throw std::domain_error("wsl::random::wishart : échelle non symétrique définie positive");
                for (std::size_t i = 0; i < p; ++i)
                    chi.push_back(gamma<double>(typename gamma<double>::shape_scale{sd::fmul(double(n) - double(i), 0.5), 2.0}));
                log_det_half = half_log_determinant(l);
                log_gamma_p = log_multivariate_gamma(p, sd::fmul(double(n), 0.5));
            }

            constexpr std::size_t dimension() const noexcept { return v.size(); }

            // Entropie de W_p(LLᵀ, n), ½ log det = s (s = ± log_det_half).
            constexpr sd::double_double wishart_entropy(sd::double_double s) const noexcept
            {
                double const p = double(dimension()), nn = n;
                auto h = mul({p + 1, 0.0}, s);                              // ((p + 1)/2) log det V
                h = sum(h, mul({sd::fmul(sd::fmul(p, p + 1), 0.5), 0.0}, {sd::ln2_hi, sd::ln2_lo}));
                h = sum(h, log_gamma_p);
                auto const psi = product(sd::fmul(nn - p - 1, 0.5), multivariate_digamma(dimension(), sd::fmul(nn, 0.5)));
                h = sum(h, {-psi.hi, -psi.lo});
                return sum(h, product(nn, sd::fmul(p, 0.5)));
            }

            constexpr interval<square_matrix<T, N>> matrix_support() const
            {
                auto lo = make_matrix<T, N>(dimension()), up = make_matrix<T, N>(dimension());
                for (std::size_t i = 0; i < dimension(); ++i)
                    for (std::size_t j = 0; j < dimension(); ++j)
                    {
                        lo(i, j) = i == j ? T(0) : -std::numeric_limits<T>::infinity();
                        up(i, j) = std::numeric_limits<T>::infinity();
                    }
                return {lo, up};
            }

            constexpr square_matrix<T, N> scaled_scale(double c) const      // c V, un arrondi par coefficient
            {
                auto m = make_matrix<T, N>(dimension());
                for (std::size_t i = 0; i < dimension(); ++i)
                    for (std::size_t j = 0; j < dimension(); ++j)
                        m(i, j) = static_cast<T>(sd::fmul(c, double(v(i, j))));
                return m;
            }
        };
    }

    template <law_real T = double, std::size_t N = dynamic>
    class wishart : public law_interface
    {
    public:
        using value_type = square_matrix<T, N>;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T degrees_of_freedom;
            square_matrix<T, N> scale;
        };

        constexpr wishart(parameters p) : b_(p.degrees_of_freedom, p.scale, "wsl::random::wishart : dimension nulle") {}

        constexpr std::size_t dimension() const noexcept { return b_.dimension(); }
        constexpr T degrees_of_freedom() const noexcept { return b_.n; }
        constexpr square_matrix<T, N> const& scale() const noexcept { return b_.v; }
        constexpr interval<value_type> support() const { return b_.matrix_support(); }

        constexpr value_type mean() const { return b_.scaled_scale(double(b_.n)); }
        // n (Vᵢⱼ² + Vᵢᵢ Vⱼⱼ)
        constexpr value_type variance() const
        {
            auto m = detail::make_matrix<T, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                for (std::size_t j = 0; j < dimension(); ++j)
                {
                    double const vij = b_.v(i, j), vii = b_.v(i, i), vjj = b_.v(j, j);
                    auto const s = detail::sum(detail::product(vij, vij), detail::product(vii, vjj));
                    m(i, j) = static_cast<T>(detail::rounded(detail::normalized(double(b_.n), 0.0) * detail::normalized(s.hi, s.lo)));
                }
            return m;
        }
        constexpr value_type mode() const
        {
            double const c = double(b_.n) - double(dimension()) - 1;
            return c >= 0 ? b_.scaled_scale(c) : b_.scaled_scale(std::numeric_limits<double>::quiet_NaN());
        }
        constexpr double entropy() const noexcept
        {
            auto const h = b_.wishart_entropy(b_.log_det_half);
            return h.hi + h.lo;
        }
        // Wᵢᵢ ~ Vᵢᵢ χ²(n) = gamma(n/2, échelle 2Vᵢᵢ)
        constexpr gamma<T> marginal(std::size_t i) const
        {
            return gamma<T>(typename gamma<T>::shape_scale{static_cast<T>(detail::sd::fmul(double(b_.n), 0.5)), static_cast<T>(detail::sd::fmul(2.0, double(b_.v(i, i))))});
        }

        constexpr double log_pdf(value_type const& x) const
        {
            auto const& W = *this;
            auto const l = W.log_density(x);
            return l.hi + l.lo;
        }
        constexpr double pdf(value_type const& x) const
        {
            auto const& W = *this;
            auto const l = W.log_density(x);
            if (!(l.hi > -1500))
                return 0;
            if (l.hi > 1400)
                return std::numeric_limits<double>::infinity();
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& W = *this;
            auto const a = detail::bartlett<N>(W.b_.chi, g);
            return detail::gram<T>(detail::multiply(W.b_.l, a));
        }

        friend constexpr bool operator==(wishart const& a, wishart const& b) noexcept { return a.b_.n == b.b_.n && a.b_.v == b.b_.v; }

    private:
        constexpr detail::sd::double_double log_density(value_type const& x) const
        {
            namespace fp = detail::sd;
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            std::size_t const p = dimension();
            if (x.size() != p)
                return none;
            auto lx = detail::make_matrix<double, N>(p);
            if (!detail::cholesky(x, lx))
                return none;
            double const nn = b_.n;
            auto l = detail::mul({nn - double(p) - 1, 0.0}, detail::half_log_determinant(lx));   // ((n − p − 1)/2) log det X
            auto const t = detail::trace_solve(b_.l, lx);                   // tr(V⁻¹ X)
            l = detail::sum(l, {fp::fmul(t.hi, -0.5), fp::fmul(t.lo, -0.5)});
            auto const c = detail::mul({fp::fmul(nn, fp::fmul(double(p), 0.5)), 0.0}, {fp::ln2_hi, fp::ln2_lo});
            l = detail::sum(l, {-c.hi, -c.lo});
            auto const d = detail::mul({nn, 0.0}, b_.log_det_half);         // (n/2) log det V
            l = detail::sum(l, {-d.hi, -d.lo});
            return detail::sum(l, {-b_.log_gamma_p.hi, -b_.log_gamma_p.lo});
        }

        detail::wishart_base<T, N> b_;
    };

    template <law_real T = double, std::size_t N = dynamic>
    class inverse_wishart : public law_interface
    {
    public:
        using value_type = square_matrix<T, N>;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T degrees_of_freedom;
            square_matrix<T, N> scale;
        };

        constexpr inverse_wishart(parameters p) : b_(p.degrees_of_freedom, p.scale, "wsl::random::inverse_wishart : dimension nulle") {}

        constexpr std::size_t dimension() const noexcept { return b_.dimension(); }
        constexpr T degrees_of_freedom() const noexcept { return b_.n; }
        constexpr square_matrix<T, N> const& scale() const noexcept { return b_.v; }
        constexpr interval<value_type> support() const { return b_.matrix_support(); }

        // Ψ/(ν − p − 1), ν > p + 1
        constexpr value_type mean() const
        {
            double const c = double(b_.n) - double(dimension()) - 1;
            return b_.scaled_scale(c > 0 ? 1 / c : std::numeric_limits<double>::quiet_NaN());
        }
        constexpr value_type mode() const { return b_.scaled_scale(1 / (double(b_.n) + double(dimension()) + 1)); }
        // H(W(Ψ⁻¹, ν)) − (p + 1)(ψ_p(ν/2) + p log 2 − log det Ψ)
        constexpr double entropy() const noexcept
        {
            namespace fp = detail::sd;
            double const p = double(dimension());
            auto h = b_.wishart_entropy({-b_.log_det_half.hi, -b_.log_det_half.lo});
            auto e = detail::sum({multivariate_digamma_half(), 0.0}, detail::mul({p, 0.0}, {fp::ln2_hi, fp::ln2_lo}));
            e = detail::sum(e, {fp::fmul(b_.log_det_half.hi, -2.0), fp::fmul(b_.log_det_half.lo, -2.0)});
            auto const t = detail::mul({p + 1, 0.0}, e);
            h = detail::sum(h, {-t.hi, -t.lo});
            return h.hi + h.lo;
        }

        constexpr double log_pdf(value_type const& x) const
        {
            auto const& W = *this;
            auto const l = W.log_density(x);
            return l.hi + l.lo;
        }
        constexpr double pdf(value_type const& x) const
        {
            auto const& W = *this;
            auto const l = W.log_density(x);
            if (!(l.hi > -1500))
                return 0;
            if (l.hi > 1400)
                return std::numeric_limits<double>::infinity();
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& W = *this;
            auto const a = detail::bartlett<N>(W.b_.chi, g);
            auto const inv = detail::lower_inverse(a);                      // A⁻¹
            std::size_t const p = W.dimension();
            auto t = detail::make_matrix<double, N>(p);                     // A⁻ᵀ
            for (std::size_t i = 0; i < p; ++i)
                for (std::size_t j = 0; j < p; ++j)
                    t(i, j) = inv(j, i);
            return detail::gram<T>(detail::multiply(W.b_.l, t));
        }

        friend constexpr bool operator==(inverse_wishart const& a, inverse_wishart const& b) noexcept { return a.b_.n == b.b_.n && a.b_.v == b.b_.v; }

    private:
        constexpr double multivariate_digamma_half() const noexcept { return detail::multivariate_digamma(dimension(), detail::sd::fmul(double(b_.n), 0.5)); }

        constexpr detail::sd::double_double log_density(value_type const& x) const
        {
            namespace fp = detail::sd;
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            std::size_t const p = dimension();
            if (x.size() != p)
                return none;
            auto lx = detail::make_matrix<double, N>(p);
            if (!detail::cholesky(x, lx))
                return none;
            double const nu = b_.n;
            auto l = detail::mul({nu, 0.0}, b_.log_det_half);              // (ν/2) log det Ψ
            auto const c = detail::mul({fp::fmul(nu, fp::fmul(double(p), 0.5)), 0.0}, {fp::ln2_hi, fp::ln2_lo});
            l = detail::sum(l, {-c.hi, -c.lo});
            l = detail::sum(l, {-b_.log_gamma_p.hi, -b_.log_gamma_p.lo});
            auto const d = detail::mul({nu + double(p) + 1, 0.0}, detail::half_log_determinant(lx));
            l = detail::sum(l, {-d.hi, -d.lo});
            auto const t = detail::trace_solve(lx, b_.l);                   // tr(Ψ X⁻¹) = ‖L_X⁻¹ L‖²
            return detail::sum(l, {fp::fmul(t.hi, -0.5), fp::fmul(t.lo, -0.5)});
        }

        detail::wishart_base<T, N> b_;
    };
}
