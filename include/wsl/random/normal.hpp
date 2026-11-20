#pragma once

// Loi normale (gaussienne) de moyenne μ et d'écart type σ > 0 :
// densité e^{−z²/2} / (σ √(2π)), z = (x − μ)/σ.
//
//   wsl::random::normal X({.mean = 0, .sd = 1});
//   X.variance(); X.pdf(1.0); X.cdf(-3.0); X.quantile(0.975); X(g);
//
// Tirage exact : |Z| tirée exactement (rejet exponentiel décidé en arithmétique
// d'intervalles sur des réels paresseux, voir detail/lazy.hpp), signe uniforme,
// puis l'arrondi correct de μ + σ Z. Tirage rapide (wsl::random::fast) : ziggourat.
//
// Densité, répartition, survie, quantile : z en double-double, erfc et son
// développement au premier ordre pour l'argument double-double, un seul arrondi
// final ; quantile par erfc_inv puis un pas de Newton en double-double. Erreurs
// mesurées dans tests/test_random_precision.cpp.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/gaussian.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/detail/ziggurat_tables.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    namespace detail
    {
        // Ziggourat de Marsaglia et Tsang (2000), 256 couches : un mot de 64 bits
        // donne la couche (8 bits), le signe (1 bit) et l'abscisse (53 bits) ; la
        // queue au-delà de R par la méthode de Marsaglia (1964). Z normale centrée
        // réduite, approchée (abscisses sur 53 bits), bit-exacte.
        template <class G>
        constexpr double ziggurat_normal(G& g)
        {
            for (;;)
            {
                std::uint64_t const b = random_bits(g);
                auto const i = static_cast<std::size_t>(b & 0xff);
                bool const negative = ((b >> 8) & 1) != 0;
                double const x = sd::fmul(static_cast<double>(b >> 11) * 0x1p-53, ziggurat_x[i]);
                if (x < ziggurat_x[i + 1])                                  // dans le rectangle : 99 % des cas
                    return negative ? -x : x;
                if (i == 0)                                                 // queue : R + a, a = −log(U₁)/R,
                {                                                           // si −2 log(U₂) > a²
                    double a = 0, t = 0;
                    do
                    {
                        a = sd::fmul(-wsl::special::log(uniform_positive(g)), ziggurat_inv_r);
                        t = -wsl::special::log(uniform_positive(g));
                    } while (t + t < sd::fmul(a, a));
                    double const y = ziggurat_r + a;
                    return negative ? -y : y;
                }
                double const y = ziggurat_f[i] + sd::fmul(uniform_real<double>(g), ziggurat_f[i + 1] - ziggurat_f[i]);
                if (y < wsl::special::exp(sd::fmul(sd::fmul(x, x), -0.5)))  // dans le coin, sous la courbe
                    return negative ? -x : x;
            }
        }
    }

    template <law_real T = double>
    class normal : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            T mean, sd;
        };

        constexpr normal(parameters p) : mean_(p.mean), sd_(p.sd)
        {
            if (!detail::is_finite(mean_) || !(sd_ > 0) || !detail::is_finite(sd_))
                throw std::domain_error("wsl::random::normal : il faut mean fini et sd > 0, fini");
        }

        constexpr T sd() const noexcept { return sd_; }
        constexpr interval<T> support() const noexcept
        {
            return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()};
        }

        constexpr T mean() const noexcept { return mean_; }
        constexpr T median() const noexcept { return mean_; }
        constexpr T mode() const noexcept { return mean_; }
        constexpr T variance() const noexcept { return static_cast<T>(detail::sd::fmul(double(sd_), double(sd_))); }
        constexpr T skewness() const noexcept { return 0; }
        constexpr T excess_kurtosis() const noexcept { return 0; }
        // log σ + (1 + log 2π)/2.
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            auto const [hi, lo] = detail::log_precise(double(sd_));
            auto const [s, s_error] = fp::two_sum(hi, detail::normal_entropy_hi);
            return static_cast<T>(s + (s_error + (lo + detail::normal_entropy_lo)));
        }

        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            auto const [z_hi, z_lo] = detail::standardize(x, X.mean_, X.sd_);
            double const a = detail::magnitude(z_hi);
            if (!(a < 54.6))                                                // e^{−z²/2}/(σ√2π) < 2^-1075 ; NaN
                return a == a ? T(0) : T(z_hi);
            detail::scaled e{1.0, 0.0, 0};
            if (a >= 0x1p-480)                                              // sinon e^{−z²/2} s'arrondit à 1
            {
                auto const [s, s_error] = fp::two_product(z_hi, z_hi);
                double const s_lo = s_error + fp::fmul(2.0, fp::fmul(z_hi, z_lo));
                e = detail::exp_scaled(fp::fmul(s, -0.5), fp::fmul(s_lo, -0.5));
            }
            detail::scaled const c = detail::normalized(detail::inv_sqrt_2pi_hi, detail::inv_sqrt_2pi_lo);
            return static_cast<T>(detail::rounded(e * c / double(X.sd_)));
        }

        // −z²/2 − log σ − log(2π)/2.
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            auto const [z_hi, z_lo] = detail::standardize(x, X.mean_, X.sd_);
            double const a = detail::magnitude(z_hi);
            if (!(a < 0x1p500))                                             // z²/2 domine (−∞ compris) ; NaN
                return static_cast<T>(-fp::fmul(fp::fmul(z_hi, z_hi), 0.5));
            double y_hi = 0, y_lo = 0;
            if (a >= 0x1p-480)
            {
                auto const [s, s_error] = fp::two_product(z_hi, z_hi);
                y_hi = fp::fmul(s, 0.5);
                y_lo = fp::fmul(s_error + fp::fmul(2.0, fp::fmul(z_hi, z_lo)), 0.5);
            }
            auto const [l_hi, l_lo] = detail::log_precise(double(X.sd_));
            auto const [s1, e1] = fp::two_sum(-y_hi, -l_hi);
            auto const [s2, e2] = fp::two_sum(s1, -detail::half_log_2pi_hi);
            return static_cast<T>(s2 + ((e1 + e2) - ((y_lo + l_lo) + detail::half_log_2pi_lo)));
        }

        // Φ(z) = erfc(−z/√2)/2.
        constexpr T cdf(T x) const noexcept
        {
            auto const& X = *this;
            auto const [w_hi, w_lo] = X.scaled_argument(x);
            return static_cast<T>(detail::half_erfc(-w_hi, -w_lo));
        }

        // 1 − Φ(z) = erfc(z/√2)/2, sans annulation.
        constexpr T ccdf(T x) const noexcept
        {
            auto const& X = *this;
            auto const [w_hi, w_lo] = X.scaled_argument(x);
            return static_cast<T>(detail::half_erfc(w_hi, w_lo));
        }

        // μ + σ Φ^{-1}(p).
        constexpr T quantile(double p) const
        {
            auto const& X = *this;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::random::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return -std::numeric_limits<T>::infinity();
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            auto const [x_hi, x_lo] = detail::standard_quantile(p);
            if (x_hi == 0)
                return X.mean_;
            return static_cast<T>(detail::affine(X.mean_, X.sd_, x_hi, x_lo));
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
            auto const [x_hi, x_lo] = detail::standard_quantile(q);         // μ − σ Φ⁻¹(q)
            return static_cast<T>(detail::affine(X.mean_, X.sd_, -x_hi, -x_lo));
        }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            detail::lazy_real v = detail::exact_half_normal(g);
            bool const negative = (random_bits(g) >> 63) != 0;
            return detail::round_lazy<T>(X.mean_, X.sd_, 1.0, negative, v, g);
        }

        // Ziggourat (approché, bit-exact) : μ + σ Z, deux arrondis.
        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g, fast_t) const
        {
            auto const& X = *this;
            return static_cast<T>(double(X.mean_) + detail::sd::fmul(double(X.sd_), detail::ziggurat_normal(g)));
        }

        friend constexpr bool operator==(normal const&, normal const&) = default;

    private:
        // z/√2 en double-double.
        constexpr detail::sd::double_double scaled_argument(T x) const noexcept
        {
            namespace fp = detail::sd;
            auto const [z_hi, z_lo] = detail::standardize(x, mean_, sd_);
            if (z_hi == 0 || !detail::is_finite(z_hi))
                return {z_hi, 0.0};
            auto const [m, m_error] = fp::two_product(z_hi, detail::inv_sqrt2_hi);
            auto const [hi, lo] = fp::fast_two_sum(m, m_error + (fp::fmul(z_hi, detail::inv_sqrt2_lo) + fp::fmul(z_lo, detail::inv_sqrt2_hi)));
            return {hi, lo};
        }

        T mean_, sd_;
    };
}
