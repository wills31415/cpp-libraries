#pragma once

// Loi log-normale : X = e^Y, Y normale de moyenne μ et d'écart type σ > 0 (paramètres
// de log X) ; densité φ((log x − μ)/σ)/(σx) sur (0, +∞).
//
//   wsl::random::lognormal X({.mu = 0, .sigma = 0.5});
//   X.median(); X.cdf(2.0); X.quantile(0.9); X(g);
//
// Les fonctions de la loi passent par celles de la normale : z = (log x − μ)/σ avec
// log x en double-double, erfc corrigée au premier ordre pour la répartition ;
// quantile exp(μ + σΦ⁻¹(p)), l'exposant en double-double, un seul arrondi. Moments
// e^{μ + σ²/2}, (e^{σ²} − 1) e^{2μ + σ²}… de même.
//
// Tirage (approché, bit-exact) : e^Y, Y tirée exactement (wsl::random::normal) ; wsl::random::fast :
// Y par la ziggourat.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    template <law_real T = double>
    class lognormal : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T mu, sigma;
        };

        constexpr lognormal(parameters p) : mu_(p.mu), sigma_(p.sigma)
        {
            if (!detail::is_finite(mu_) || !(sigma_ > 0) || !detail::is_finite(sigma_))
                throw std::domain_error("wsl::random::lognormal : il faut mu fini et sigma > 0, fini");
        }

        constexpr T mu() const noexcept { return mu_; }
        constexpr T sigma() const noexcept { return sigma_; }
        constexpr interval<T> support() const noexcept { return {T(0), std::numeric_limits<T>::infinity()}; }

        constexpr T median() const noexcept { return exponential_of(mu_, 0.0); }
        constexpr T mean() const noexcept { return exponential_of(mu_, 0.5); }       // e^{μ + σ²/2}
        constexpr T mode() const noexcept { return exponential_of(mu_, -1.0); }      // e^{μ − σ²}
        // (e^{σ²} − 1) e^{2μ + σ²}
        constexpr T variance() const noexcept
        {
            namespace fp = detail::sd;
            auto const [s_hi, s_lo] = fp::two_product(double(sigma_), double(sigma_));
            fp::double_double m1{};
            if (s_hi <= 1)
                m1 = detail::expm1_double_double(s_hi, s_lo);
            else
            {
                auto const [e_hi, e_lo] = detail::unscaled(detail::exp_scaled(s_hi > 1400 ? 1400 : s_hi, s_lo));
                m1 = detail::sum({e_hi, e_lo}, {-1.0, 0.0});
            }
            auto const [t_hi, t_lo] = detail::sum({fp::fmul(double(mu_), 2.0), 0.0}, {s_hi, s_lo});
            if (t_hi > 1400 || s_hi > 1400)
                return std::numeric_limits<T>::infinity();
            return static_cast<T>(detail::rounded(detail::normalized(m1.hi, m1.lo) * detail::exp_scaled(t_hi, t_lo)));
        }
        constexpr T skewness() const noexcept
        {
            double const e = wsl::special::exp(detail::sd::fmul(double(sigma_), double(sigma_)));
            return static_cast<T>(detail::sd::fmul(e + 2, wsl::special::sqrt(e - 1)));
        }
        constexpr T excess_kurtosis() const noexcept
        {
            namespace fp = detail::sd;
            double const e = wsl::special::exp(fp::fmul(double(sigma_), double(sigma_)));
            double const e2 = fp::fmul(e, e);
            return static_cast<T>(((fp::fmul(e2, e2) + fp::fmul(2.0, fp::fmul(e2, e))) + fp::fmul(3.0, e2)) - 6.0);
        }
        // μ + log σ + (1 + log 2π)/2
        constexpr T entropy() const noexcept
        {
            auto const [hi, lo] = detail::log_precise(double(sigma_));
            auto e = detail::sum({double(mu_), 0.0}, {hi, lo});
            e = detail::sum(e, {detail::normal_entropy_hi, detail::normal_entropy_lo});
            return static_cast<T>(e.hi + e.lo);
        }

        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            double const v = x;
            if (!(v > 0))
                return v == v ? T(0) : x;
            if (v == std::numeric_limits<double>::infinity())
                return 0;
            auto const [z_hi, z_lo] = X.standard(v);
            if (!(detail::magnitude(z_hi) < 54.6))
                return 0;
            detail::scaled e{1.0, 0.0, 0};
            if (detail::magnitude(z_hi) >= 0x1p-480)
            {
                auto const [s, s_error] = fp::two_product(z_hi, z_hi);
                double const s_lo = s_error + fp::fmul(2.0, fp::fmul(z_hi, z_lo));
                e = detail::exp_scaled(fp::fmul(s, -0.5), fp::fmul(s_lo, -0.5));
            }
            detail::scaled const c = detail::normalized(detail::inv_sqrt_2pi_hi, detail::inv_sqrt_2pi_lo);
            return static_cast<T>(detail::rounded(e * c / double(X.sigma_) / v));
        }

        // −z²/2 − log σ − log x − log(2π)/2
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            double const v = x;
            if (!(v > 0))
                return v == v ? -std::numeric_limits<T>::infinity() : x;
            if (v == std::numeric_limits<double>::infinity())
                return -std::numeric_limits<T>::infinity();
            auto const [z_hi, z_lo] = X.standard(v);
            auto const [s, s_error] = fp::two_product(z_hi, z_hi);
            fp::double_double e{fp::fmul(s, -0.5), fp::fmul(s_error + fp::fmul(2.0, fp::fmul(z_hi, z_lo)), -0.5)};
            auto const ls = detail::log_precise(double(X.sigma_)), lx = detail::log_precise(v);
            e = detail::sum(e, {-ls.hi, -ls.lo});
            e = detail::sum(e, {-lx.hi, -lx.lo});
            e = detail::sum(e, {-detail::half_log_2pi_hi, -detail::half_log_2pi_lo});
            return static_cast<T>(e.hi + e.lo);
        }

        constexpr T cdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, true); }

        // e^{μ + σ Φ⁻¹(p)}
        constexpr T quantile(double p) const
        {
            auto const& X = *this;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::random::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return 0;
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            auto const [x_hi, x_lo] = detail::standard_quantile(p);
            auto const y = detail::sum({double(X.mu_), 0.0}, detail::mul({double(X.sigma_), 0.0}, {x_hi, x_lo}));
            if (y.hi > 1400)
                return std::numeric_limits<T>::infinity();
            if (y.hi < -1500)
                return 0;
            return static_cast<T>(detail::rounded(detail::exp_scaled(y.hi, y.lo)));
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        constexpr T cquantile(double q) const
        {
            auto const& X = *this;
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::random::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return wsl::random::quantile(X, 1.0 - q);                                // 1 − q exact
            if (q == 0)
                return std::numeric_limits<T>::infinity();
            auto const [x_hi, x_lo] = detail::standard_quantile(q);         // e^{μ − σ Φ⁻¹(q)}
            auto const y = detail::sum({double(X.mu_), 0.0}, detail::mul({double(X.sigma_), 0.0}, {-x_hi, -x_lo}));
            if (y.hi > 1400)
                return std::numeric_limits<T>::infinity();
            if (y.hi < -1500)
                return 0;
            return static_cast<T>(detail::rounded(detail::exp_scaled(y.hi, y.lo)));
        }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            return static_cast<T>(wsl::special::exp(double(wsl::random::sample(normal<double>({.mean = X.mu_, .sd = X.sigma_}), g))));
        }
        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g, fast_t) const
        {
            auto const& X = *this;
            return static_cast<T>(wsl::special::exp(double(X.mu_) + detail::sd::fmul(double(X.sigma_), detail::ziggurat_normal(g))));
        }

        friend constexpr bool operator==(lognormal const&, lognormal const&) = default;

    private:
        // e^{μ + c σ²}, c dans {0, ½, −1}, l'exposant en double-double.
        constexpr T exponential_of(double mu, double c) const noexcept
        {
            namespace fp = detail::sd;
            auto const [s_hi, s_lo] = fp::two_product(double(sigma_), double(sigma_));
            auto const y = detail::sum({mu, 0.0}, {fp::fmul(s_hi, c), fp::fmul(s_lo, c)});   // c σ² exact (c = ½, −1)
            if (y.hi > 1400)
                return std::numeric_limits<T>::infinity();
            if (y.hi < -1500)
                return 0;
            return static_cast<T>(detail::rounded(detail::exp_scaled(y.hi, y.lo)));
        }

        // z = (log x − μ)/σ en double-double.
        constexpr detail::sd::double_double standard(double x) const noexcept
        {
            namespace fp = detail::sd;
            auto const [l_hi, l_lo] = detail::log_precise(x);
            auto const [d, d_error] = fp::two_sum(l_hi, -double(mu_));
            auto const [n_hi, n_lo] = fp::fast_two_sum(d, d_error + l_lo);
            if (n_hi == 0)
                return {0.0, 0.0};
            bool const negative = n_hi < 0;
            detail::scaled const q = detail::normalized(negative ? -n_hi : n_hi, negative ? -n_lo : n_lo) / double(sigma_);
            auto const [hi, lo] = detail::unscaled(q);
            return negative ? fp::double_double{-hi, -lo} : fp::double_double{hi, lo};
        }

        // Φ(z) (ou 1 − Φ(z)) = erfc(∓z/√2)/2.
        constexpr T distribution(T x, bool upper) const noexcept
        {
            namespace fp = detail::sd;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0))
                return upper ? T(1) : T(0);
            if (v == std::numeric_limits<double>::infinity())
                return upper ? T(0) : T(1);
            auto const [z_hi, z_lo] = standard(v);
            double w_hi = z_hi, w_lo = 0;                                   // |z| immense : seul le signe compte
            if (z_hi != 0 && detail::magnitude(z_hi) < 0x1p500)
            {
                auto const [m, m_error] = fp::two_product(z_hi, detail::inv_sqrt2_hi);
                auto const [h, l] = fp::fast_two_sum(m, m_error + (fp::fmul(z_hi, detail::inv_sqrt2_lo) + fp::fmul(z_lo, detail::inv_sqrt2_hi)));
                w_hi = h;
                w_lo = l;
            }
            return static_cast<T>(upper ? detail::half_erfc(w_hi, w_lo) : detail::half_erfc(-w_hi, -w_lo));
        }

        T mu_, sigma_;
    };
}
