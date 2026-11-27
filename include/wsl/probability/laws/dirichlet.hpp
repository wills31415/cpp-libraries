#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Loi de Dirichlet de paramètres α₁, …, α_K > 0 (K >= 2), sur le simplexe
// {x : xᵢ > 0, Σ xᵢ = 1}.
//
//   wsl::probability::dirichlet<double, 3> D({.alpha = {1.0, 2.0, 3.0}});
//   wsl::probability::dirichlet<double, wsl::probability::dynamic> E({.alpha = std::vector<double>(10, 0.5)});
//   D.mean(); D.covariance(); D.pdf({0.2, 0.3, 0.5}); D.marginal(0); D(g);
//
// Densité Γ(α₀)/Π Γ(αᵢ) Π xᵢ^{αᵢ − 1}, α₀ = Σ αᵢ, calculée en logarithme double-double
// (log B(α) = Σ lgamma(αᵢ) − lgamma(α₀)) puis un arrondi ; x doit être dans le
// simplexe à K·2^-50 près (Σ xᵢ = 1 ne tombe pas juste en flottant), sinon 0.
// Moyenne αᵢ/α₀, variances et covariances, mode (αᵢ > 1), entropie
// log B(α) + (α₀ − K) ψ(α₀) − Σ (αᵢ − 1) ψ(αᵢ) ; marginale i : bêta(αᵢ, α₀ − αᵢ).
//
// Tirage (approché, bit-exact) : xᵢ = Gᵢ/Σ G, Gᵢ de lois gamma(αᵢ), tirées dans l'ordre.

#include <wsl/probability/laws/beta.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/loader.hpp>
#include <wsl/probability/laws/gamma.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/multivariate.hpp>
#include <wsl/special.hpp>

#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_real T = double, std::size_t N = dynamic>
    class dirichlet : public law_interface
    {
    public:
        using value_type = vector_of<T, N>;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            vector_of<T, N> alpha;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.alpha = alpha()}; }

        constexpr dirichlet(parameters p) : alpha_(p.alpha)
        {
            if (alpha_.size() < 2)
                throw std::domain_error("wsl::probability::dirichlet : il faut au moins deux paramètres");
            for (T const a : alpha_)
                if (!(a > 0) || !detail::is_finite(a))
                    throw std::domain_error("wsl::probability::dirichlet : il faut alpha > 0, fini");
            detail::sd::double_double s{0.0, 0.0}, lb{0.0, 0.0};
            for (T const a : alpha_)
            {
                s = detail::sum(s, {double(a), 0.0});
                lb = detail::sum(lb, detail::log_gamma(double(a)));
                gammas_.push_back(gamma<double>(typename gamma<double>::parameters{.rate = 1.0, .shape = double(a)}));
            }
            total_ = s.hi + s.lo;
            auto const l0 = detail::log_gamma(total_);
            log_beta_ = detail::sum(lb, {-l0.hi, -l0.lo});
        }

        constexpr std::size_t dimension() const noexcept { return alpha_.size(); }
        constexpr value_type const& alpha() const noexcept { return alpha_; }
        constexpr interval<value_type> support() const noexcept
        {
            return {detail::make_vector<T, N>(dimension(), T(0)), detail::make_vector<T, N>(dimension(), T(1))};
        }

        // αᵢ/α₀
        constexpr value_type mean() const noexcept
        {
            value_type m = detail::make_vector<T, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                m[i] = static_cast<T>(double(alpha_[i]) / total_);
            return m;
        }
        // αᵢ (α₀ − αᵢ) / (α₀² (α₀ + 1))
        constexpr value_type variance() const noexcept
        {
            value_type v = detail::make_vector<T, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                v[i] = static_cast<T>(covariance_entry(i, i));
            return v;
        }
        constexpr square_matrix<T, N> covariance() const
        {
            square_matrix<T, N> c = make_matrix();
            for (std::size_t i = 0; i < dimension(); ++i)
                for (std::size_t j = 0; j < dimension(); ++j)
                    c[i, j] = static_cast<T>(covariance_entry(i, j));
            return c;
        }
        // (αᵢ − 1)/(α₀ − K), tous les αᵢ > 1 ; NaN sinon.
        constexpr value_type mode() const noexcept
        {
            value_type m = detail::make_vector<T, N>(dimension());
            double const d = total_ - double(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                m[i] = alpha_[i] > 1 ? static_cast<T>((double(alpha_[i]) - 1) / d) : std::numeric_limits<T>::quiet_NaN();
            return m;
        }
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            double const k = double(dimension());
            auto h = detail::sum(log_beta_, detail::product(total_ - k, wsl::special::digamma(total_)));
            for (T const a : alpha_)
            {
                auto const t = detail::product(double(a) - 1, wsl::special::digamma(double(a)));
                h = detail::sum(h, {-t.hi, -t.lo});
            }
            return static_cast<T>(h.hi + h.lo);
        }
        // Loi de la composante i : bêta(αᵢ, α₀ − αᵢ).
        constexpr beta<T> marginal(std::size_t i) const
        {
            auto const r = detail::sum({total_, 0.0}, {-double(alpha_[i]), 0.0});
            return beta<T>({.alpha = alpha_[i], .beta = static_cast<T>(r.hi + r.lo)});
        }

        constexpr double log_pdf(value_type const& x) const noexcept
        {
            auto const& D = *this;
            if (x.size() != D.dimension() || !D.in_simplex(x))
                return -std::numeric_limits<double>::infinity();
            detail::sd::double_double l{-D.log_beta_.hi, -D.log_beta_.lo};
            for (std::size_t i = 0; i < D.dimension(); ++i)
                if (D.alpha_[i] != 1)
                    l = detail::sum(l, detail::mul({double(D.alpha_[i]) - 1, 0.0}, detail::log_precise(double(x[i]))));
            return l.hi + l.lo;
        }
        constexpr double pdf(value_type const& x) const noexcept
        {
            auto const& D = *this;
            if (x.size() != D.dimension() || !D.in_simplex(x))
                return 0;
            detail::sd::double_double l{-D.log_beta_.hi, -D.log_beta_.lo};
            for (std::size_t i = 0; i < D.dimension(); ++i)
                if (D.alpha_[i] != 1)
                    l = detail::sum(l, detail::mul({double(D.alpha_[i]) - 1, 0.0}, detail::log_precise(double(x[i]))));
            if (l.hi > 1400)
                return std::numeric_limits<double>::infinity();
            if (l.hi < -1500)
                return 0;
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& D = *this;
            value_type x = detail::make_vector<T, N>(D.dimension());
            detail::sd::double_double s{0.0, 0.0};
            for (std::size_t i = 0; i < D.dimension(); ++i)
            {
                double const v = wsl::probability::sample(D.gammas_[i], g);
                x[i] = static_cast<T>(v);
                s = detail::sum(s, {v, 0.0});
            }
            double const total = s.hi + s.lo;
            if (!(total > 0))                                               // αᵢ minuscules : tout s'arrondit à 0
            {
                for (std::size_t i = 0; i < D.dimension(); ++i)
                    x[i] = T(0);
                x[0] = T(1);
                return x;
            }
            for (std::size_t i = 0; i < D.dimension(); ++i)
                x[i] = static_cast<T>(double(x[i]) / total);
            return x;
        }

        friend constexpr bool operator==(dirichlet const& a, dirichlet const& b) noexcept { return a.alpha_ == b.alpha_; }

    private:
        constexpr square_matrix<T, N> make_matrix() const
        {
            if constexpr (N == dynamic)
                return square_matrix<T, N>(dimension());
            else
                return square_matrix<T, N>();
        }

        // Cov(Xᵢ, Xⱼ) = (δᵢⱼ αᵢ α₀ − αᵢ αⱼ) / (α₀² (α₀ + 1)), numérateur exact en double-double.
        constexpr double covariance_entry(std::size_t i, std::size_t j) const noexcept
        {
            double const ai = alpha_[i], aj = alpha_[j];
            detail::scaled const den = detail::normalized(total_, 0.0) * detail::normalized(total_, 0.0) * detail::normalized(total_ + 1, 0.0);
            if (i == j)                                                     // αᵢ (α₀ − αᵢ)
            {
                auto const r = detail::sum({total_, 0.0}, {-ai, 0.0});
                if (!(r.hi > 0))
                    return 0;
                return detail::rounded(detail::normalized(ai, 0.0) * detail::normalized(r.hi, r.lo) / den);
            }
            return -detail::rounded(detail::normalized(ai, 0.0) * detail::normalized(aj, 0.0) / den);
        }

        // xᵢ > 0 et |Σ xᵢ − 1| <= K 2^-50.
        constexpr bool in_simplex(value_type const& x) const noexcept
        {
            detail::sd::double_double s{0.0, 0.0};
            for (T const v : x)
            {
                if (!(v > 0))
                    return false;
                s = detail::sum(s, {double(v), 0.0});
            }
            auto const d = detail::sum(s, {-1.0, 0.0});
            double const e = d.hi + d.lo;
            return (e < 0 ? -e : e) <= double(dimension()) * 0x1p-50;
        }

        vector_of<T, N> alpha_;
        std::vector<gamma<double>> gammas_;
        double total_ = 0;
        detail::sd::double_double log_beta_{};
    };
}
