#pragma once

// Lois de von Mises (sur le cercle) et de von Mises–Fisher (sur la sphère Sᵖ⁻¹).
//
//   wsl::von_mises V({.mu = 0.5, .kappa = 4});                   // angles dans [μ − π, μ + π]
//   wsl::von_mises_fisher<double, 3> F({.mean_direction = {0, 0, 1}, .kappa = 10});
//   pdf(V, 1.0); V.circular_variance(); sample(V, g); pdf(F, x); F.mean(); sample(F, g);
//
// von Mises : densité e^{κ cos(x − μ)}/(2π I₀(κ)) sur [μ − π, μ + π], calculée comme
// exp(−2κ sin²((x − μ)/2) − log 2π − log(e^{−κ} I₀(κ))) (sans compensation près du mode ;
// x − μ et le sinus en double-double ; I₀ par wsl::special, detail::log_cyl_bessel_i),
// un arrondi. Longueur résultante moyenne A = I₁(κ)/I₀(κ), variance circulaire 1 − A,
// entropie log(2π I₀(κ)) − κA. (Pas de répartition : pas de forme close.)
//
// von Mises–Fisher : densité C_p(κ) e^{κ μᵀx} par rapport à la mesure de surface,
// C_p(κ) = κ^{p/2−1}/((2π)^{p/2} I_{p/2−1}(κ)), calculée comme
// exp(κ (μᵀx − 1) + log C_p(κ) + κ), μᵀx − 1 en double-double (sans compensation près
// du mode, ni débordement pour κ grand) ; x doit être de norme 1 à p 2^-50 près (sinon 0). μ est normalisé à la
// construction. Moyenne A_p(κ) μ, A_p = I_{p/2}/I_{p/2−1} ; entropie −log C_p − κ A_p ;
// κ = 0 : loi uniforme sur la sphère.
//
// Tirages (approchés, bit-exacts) : von Mises, Best et Fisher (1979) ; von Mises–Fisher,
// Wood (1994) : la composante w selon μ par rejet (loi bêta, uniforme), une direction
// uniforme orthogonale (sphère Sᵖ⁻²), puis une réflexion de Householder qui envoie e₁ sur μ.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/beta.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/sphere.hpp>
#include <wsl/special.hpp>

#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl
{
    namespace detail
    {
        inline constexpr double pi_hi = 0x1.921fb54442d18p+1, pi_lo = 0x1.1a62633145c07p-53;

        // sin(h), |h| <= π/2 + ε, h en double-double (à 2^-70 près).
        constexpr sd::double_double sine_dd(double h_hi, double h_lo) noexcept
        {
            namespace fp = sd;
            double const a = magnitude(h_hi);
            if (a == 0)
                return {h_hi, h_lo};
            if (a <= 0x1.921fb54442d18p-1)                                  // |h| <= π/4
                return fp::sin_of(fp::sincos_prepare(h_hi, h_lo));
            bool const negative = h_hi < 0;                                 // sin|h| = cos(π/2 − |h|)
            auto const r = sum({fp::fmul(pi_hi, 0.5), fp::fmul(pi_lo, 0.5)}, {-a, negative ? h_lo : -h_lo});
            auto const c = fp::cos_of(fp::sincos_prepare(r.hi, r.lo));
            return negative ? fp::double_double{-c.hi, -c.lo} : c;
        }

        // log(e^{−κ} I_ν(κ)), κ > 0, en double-double.
        constexpr sd::double_double log_scaled_bessel(double nu, double kappa) noexcept
        {
            return sum(wsl::special::detail::log_cyl_bessel_i(nu, kappa), {-kappa, 0.0});
        }

        // A_ν = I_{ν+1}(κ)/I_ν(κ), et 1 − A_ν, en double-double.
        struct bessel_ratio
        {
            sd::double_double a, one_minus_a;
        };
        constexpr bessel_ratio bessel_ratio_of(double nu, double kappa) noexcept
        {
            auto const l0 = wsl::special::detail::log_cyl_bessel_i(nu, kappa);
            auto const l = sum(wsl::special::detail::log_cyl_bessel_i(nu + 1, kappa), {-l0.hi, -l0.lo});
            sd::double_double a{0.0, 0.0};
            if (l.hi > -1500)
                a = unscaled(exp_scaled(l.hi, l.lo));
            return {a, sum({1.0, 0.0}, {-a.hi, -a.lo})};
        }
    }

    template <law_real T = double>
    class von_mises
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T mu, kappa;
        };

        constexpr von_mises(parameters p) : mu_(p.mu), kappa_(p.kappa)
        {
            if (!detail::is_finite(mu_) || !(kappa_ >= 0) || !detail::is_finite(kappa_))
                throw std::domain_error("wsl::von_mises : il faut mu fini et kappa >= 0, fini");
            namespace fp = detail::sd;
            double const k = kappa_;
            log_norm_ = {fp::fmul(detail::half_log_2pi_hi, 2.0), fp::fmul(detail::half_log_2pi_lo, 2.0)};   // log 2π
            if (k > 0)
                log_norm_ = detail::sum(log_norm_, detail::log_scaled_bessel(0.0, k));
            // Best–Fisher : τ = 1 + √(1 + 4κ²), ρ = (τ − √(2τ))/(2κ), r = (1 + ρ²)/(2ρ).
            if (k >= 0x1p-30)
            {
                double const tau = 1 + wsl::special::sqrt(1 + fp::fmul(4.0, fp::fmul(k, k)));
                double const rho = (tau - wsl::special::sqrt(fp::fmul(2.0, tau))) / fp::fmul(2.0, k);
                r_ = (1 + fp::fmul(rho, rho)) / fp::fmul(2.0, rho);
            }
        }

        constexpr T mu() const noexcept { return mu_; }
        constexpr T kappa() const noexcept { return kappa_; }
        constexpr interval<T> support() const noexcept
        {
            return {static_cast<T>(double(mu_) - detail::pi_hi), static_cast<T>(double(mu_) + detail::pi_hi)};
        }
        constexpr T mean() const noexcept { return mu_; }                   // direction moyenne
        constexpr T median() const noexcept { return mu_; }
        constexpr T mode() const noexcept { return mu_; }
        // A = I₁(κ)/I₀(κ)
        constexpr T mean_resultant_length() const noexcept
        {
            if (kappa_ == 0)
                return 0;
            auto const r = detail::bessel_ratio_of(0.0, double(kappa_));
            return static_cast<T>(r.a.hi + r.a.lo);
        }
        constexpr T circular_variance() const noexcept
        {
            if (kappa_ == 0)
                return 1;
            auto const r = detail::bessel_ratio_of(0.0, double(kappa_));
            return static_cast<T>(r.one_minus_a.hi + r.one_minus_a.lo);
        }
        // log 2π + log(e^{−κ} I₀) + κ (1 − A)
        constexpr double entropy() const noexcept
        {
            if (kappa_ == 0)
                return log_norm_.hi + log_norm_.lo;
            auto const r = detail::bessel_ratio_of(0.0, double(kappa_));
            auto const h = detail::sum(log_norm_, detail::mul({double(kappa_), 0.0}, r.one_minus_a));
            return h.hi + h.lo;
        }

        friend constexpr T log_pdf(von_mises const& V, T x) noexcept
        {
            auto const l = V.log_density(x);
            return static_cast<T>(l.hi + l.lo);
        }
        friend constexpr T pdf(von_mises const& V, T x) noexcept
        {
            auto const l = V.log_density(x);
            if (!(l.hi > -1500))
                return l.hi == l.hi ? T(0) : static_cast<T>(l.hi);
            return static_cast<T>(detail::rounded(detail::exp_scaled(l.hi, l.lo)));
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(von_mises const& V, G& g)
        {
            namespace fp = detail::sd;
            namespace sp = wsl::special;
            double const mu = V.mu_, k = V.kappa_;
            if (V.r_ == 0)                                                  // κ minuscule : uniforme
                return static_cast<T>(detail::interpolate(mu - detail::pi_hi, mu + detail::pi_hi, uniform_real<double>(g)));
            for (;;)
            {
                double const u1 = uniform_real<double>(g), u2 = detail::uniform_positive(g), u3 = uniform_real<double>(g);
                double const z = sp::cos(fp::fmul(detail::pi_hi, u1));
                double const f = (1 + fp::fmul(V.r_, z)) / (V.r_ + z);
                double const c = fp::fmul(k, V.r_ - f);
                if (fp::fmul(c, 2 - c) - u2 > 0 || sp::log(c / u2) + 1 - c >= 0)
                {
                    double const ff = f > 1 ? 1.0 : f < -1 ? -1.0 : f;
                    double const angle = sp::atan2(sp::sqrt(fp::fmul(1 - ff, 1 + ff)), ff);   // arccos f
                    return static_cast<T>(u3 < 0.5 ? mu - angle : mu + angle);
                }
            }
        }

        friend constexpr bool operator==(von_mises const& a, von_mises const& b) noexcept { return a.mu_ == b.mu_ && a.kappa_ == b.kappa_; }

    private:
        // −2κ sin²((x − μ)/2) − log 2π − log(e^{−κ} I₀(κ))
        constexpr detail::sd::double_double log_density(T x) const noexcept
        {
            namespace fp = detail::sd;
            double const v = x;
            if (!(v == v))
                return {v, 0.0};
            auto const [d, d_error] = fp::two_sum(v, -double(mu_));
            auto const lim = detail::sum({detail::pi_hi, detail::pi_lo}, {-detail::magnitude(d), d < 0 ? d_error : -d_error});
            if (!detail::is_finite(d) || lim.hi < 0)                        // hors de [μ − π, μ + π]
                return {-std::numeric_limits<double>::infinity(), 0.0};
            fp::double_double l{-log_norm_.hi, -log_norm_.lo};
            if (kappa_ == 0 || d == 0)
                return l;
            auto const s = detail::sine_dd(fp::fmul(d, 0.5), fp::fmul(d_error, 0.5));
            auto const s2 = detail::mul(s, s);
            auto const e = detail::mul({fp::fmul(double(kappa_), -2.0), 0.0}, s2);
            return detail::sum(l, e);
        }

        T mu_, kappa_;
        double r_ = 0;
        detail::sd::double_double log_norm_{};                              // log 2π + log(e^{−κ} I₀(κ))
    };

    template <law_real T = double, std::size_t N = dynamic>
    class von_mises_fisher
    {
    public:
        using value_type = vector_of<T, N>;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            vector_of<T, N> mean_direction;
            T kappa;
        };

        constexpr von_mises_fisher(parameters q) : mu_(q.mean_direction), kappa_(q.kappa)
        {
            namespace fp = detail::sd;
            std::size_t const n = mu_.size();
            if (n < 2)
                throw std::domain_error("wsl::von_mises_fisher : il faut une dimension >= 2");
            if (!(kappa_ >= 0) || !detail::is_finite(kappa_))
                throw std::domain_error("wsl::von_mises_fisher : il faut kappa >= 0, fini");
            auto const r2 = detail::dot(n, [&](std::size_t k) { return std::pair{double(mu_[k]), double(mu_[k])}; });
            if (!(r2.hi > 0) || !detail::is_finite(r2.hi))
                throw std::domain_error("wsl::von_mises_fisher : direction moyenne nulle ou non finie");
            detail::scaled const r = detail::sqrt(detail::normalized(r2.hi, r2.lo));
            for (auto& v : mu_)                                             // μ/|μ|
            {
                double const a = detail::magnitude(double(v));
                double const q = a == 0 ? 0.0 : detail::rounded(detail::normalized(a, 0.0) / r);
                v = static_cast<T>(v < 0 ? -q : q);
            }
            double const p = double(n), nu = fp::fmul(p, 0.5) - 1;
            double const k = kappa_;
            if (k == 0)
                log_c_plus_kappa_ = {-detail::log_sphere_area(n).hi, -detail::log_sphere_area(n).lo};
            else
            {
                // (p/2 − 1) log κ − (p/2) log 2π − log(e^{−κ} I_ν(κ))
                auto l = detail::mul({nu, 0.0}, detail::log_precise(k));
                l = detail::sum(l, detail::mul({-p, 0.0}, {detail::half_log_2pi_hi, detail::half_log_2pi_lo}));
                auto const b = detail::log_scaled_bessel(nu, k);
                log_c_plus_kappa_ = detail::sum(l, {-b.hi, -b.lo});
                // Wood : b = (p − 1)/(2κ + √(4κ² + (p − 1)²)), x₀ = (1 − b)/(1 + b), c = κx₀ + (p − 1) log(1 − x₀²).
                double const m = p - 1;
                b_ = m / (fp::fmul(2.0, k) + wsl::special::sqrt(fp::fmul(fp::fmul(4.0, k), k) + fp::fmul(m, m)));
                x0_ = (1 - b_) / (1 + b_);
                c_ = fp::fmul(k, x0_) + fp::fmul(m, wsl::special::log(fp::fmul(1 - x0_, 1 + x0_)));
            }
        }

        constexpr std::size_t dimension() const noexcept { return mu_.size(); }
        constexpr value_type const& mean_direction() const noexcept { return mu_; }
        constexpr T kappa() const noexcept { return kappa_; }
        constexpr interval<value_type> support() const noexcept
        {
            return {detail::make_vector<T, N>(dimension(), T(-1)), detail::make_vector<T, N>(dimension(), T(1))};
        }
        // A_p(κ) = I_{p/2}(κ)/I_{p/2−1}(κ)
        constexpr T mean_resultant_length() const noexcept
        {
            if (kappa_ == 0)
                return 0;
            auto const r = detail::bessel_ratio_of(nu(), double(kappa_));
            return static_cast<T>(r.a.hi + r.a.lo);
        }
        constexpr value_type mean() const noexcept
        {
            double const a = mean_resultant_length();
            value_type m = detail::make_vector<T, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                m[i] = static_cast<T>(detail::sd::fmul(a, double(mu_[i])));
            return m;
        }
        // −(log C_p + κ) + κ (1 − A_p)
        constexpr double entropy() const noexcept
        {
            detail::sd::double_double h{-log_c_plus_kappa_.hi, -log_c_plus_kappa_.lo};
            if (kappa_ > 0)
            {
                auto const r = detail::bessel_ratio_of(nu(), double(kappa_));
                h = detail::sum(h, detail::mul({double(kappa_), 0.0}, r.one_minus_a));
            }
            return h.hi + h.lo;
        }

        friend constexpr double log_pdf(von_mises_fisher const& F, value_type const& x) noexcept
        {
            auto const l = F.log_density(x);
            return l.hi + l.lo;
        }
        friend constexpr double pdf(von_mises_fisher const& F, value_type const& x) noexcept
        {
            auto const l = F.log_density(x);
            if (!(l.hi > -1500))
                return 0;
            if (l.hi > 1400)
                return std::numeric_limits<double>::infinity();
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr value_type sample(von_mises_fisher const& F, G& g)
        {
            namespace fp = detail::sd;
            namespace sp = wsl::special;
            std::size_t const n = F.dimension();
            if (F.kappa_ == 0)
                return detail::uniform_direction<T, N>(n, g);
            double const m = double(n) - 1, k = F.kappa_;
            beta<double> const z_law({.alpha = fp::fmul(m, 0.5), .beta = fp::fmul(m, 0.5)});
            double w = 0;
            for (;;)
            {
                double const z = sample(z_law, g);
                double const u = detail::uniform_positive(g);
                w = (1 - fp::fmul(1 + F.b_, z)) / (1 - fp::fmul(1 - F.b_, z));
                double const t = 1 - fp::fmul(F.x0_, w);
                if (t > 0 && (fp::fmul(k, w) + fp::fmul(m, sp::log(t))) - F.c_ >= sp::log(u))
                    break;
            }
            // y = (w, √(1 − w²) v), v uniforme sur Sⁿ⁻², puis x = ±H y, H = I − 2uuᵀ/|u|², u = e₁ ∓ μ.
            auto const v = detail::uniform_direction<double, dynamic>(n - 1, g);
            double const s = sp::sqrt(fp::fmul(1 - w, 1 + w));
            auto y = detail::make_vector<double, N>(n);
            y[0] = w;
            for (std::size_t i = 1; i < n; ++i)
                y[i] = fp::fmul(s, v[i - 1]);
            bool const flip = F.mu_[0] < 0;                                 // u = e₁ + μ, x = −H y : |u|² >= 2
            auto u = detail::make_vector<double, N>(n);
            for (std::size_t i = 0; i < n; ++i)
                u[i] = flip ? double(F.mu_[i]) : -double(F.mu_[i]);
            u[0] = u[0] + 1;
            auto const uu = detail::dot(n, [&](std::size_t i) { return std::pair{u[i], u[i]}; });
            auto const uy = detail::dot(n, [&](std::size_t i) { return std::pair{u[i], y[i]}; });
            value_type x = detail::make_vector<T, N>(n);
            if (!(uu.hi > 0))                                               // μ = e₁
            {
                for (std::size_t i = 0; i < n; ++i)
                    x[i] = static_cast<T>(y[i]);
                return x;
            }
            double const f = fp::fmul(2.0, (uy.hi + uy.lo) / (uu.hi + uu.lo));
            for (std::size_t i = 0; i < n; ++i)
            {
                double const hy = y[i] - fp::fmul(f, u[i]);
                x[i] = static_cast<T>(flip ? -hy : hy);
            }
            return x;
        }

        friend constexpr bool operator==(von_mises_fisher const& a, von_mises_fisher const& b) noexcept { return a.mu_ == b.mu_ && a.kappa_ == b.kappa_; }

    private:
        constexpr double nu() const noexcept { return detail::sd::fmul(double(dimension()), 0.5) - 1; }

        // κ (μᵀx − 1) + log C_p + κ ; −∞ hors de la sphère.
        constexpr detail::sd::double_double log_density(value_type const& x) const noexcept
        {
            namespace fp = detail::sd;
            constexpr fp::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            std::size_t const n = dimension();
            if (x.size() != n)
                return none;
            double const off = detail::squared_norm_minus_one(x);
            if (!((off < 0 ? -off : off) <= double(n) * 0x1p-50))
                return none;
            if (kappa_ == 0)
                return log_c_plus_kappa_;
            auto const c = detail::dot(n, [&](std::size_t i) { return std::pair{double(mu_[i]), double(x[i])}; });
            auto const t = detail::mul({double(kappa_), 0.0}, detail::sum(c, {-1.0, 0.0}));   // κ (μᵀx − 1), en double-double
            return detail::sum(log_c_plus_kappa_, t);
        }

        vector_of<T, N> mu_;
        T kappa_;
        detail::sd::double_double log_c_plus_kappa_{};                      // log C_p(κ) + κ
        double b_ = 0, x0_ = 0, c_ = 0;                                     // Wood
    };
}
