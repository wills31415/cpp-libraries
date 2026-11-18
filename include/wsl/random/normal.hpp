#pragma once

// Loi normale (gaussienne) de moyenne μ et d'écart type σ > 0 :
// densité e^{−z²/2} / (σ √(2π)), z = (x − μ)/σ.
//
//   wsl::normal X({.mean = 0, .sd = 1});
//   X.variance(); pdf(X, 1.0); cdf(X, -3.0); quantile(X, 0.975); sample(X, g);
//
// Tirage exact : |Z| tirée exactement (rejet exponentiel décidé en arithmétique
// d'intervalles sur des réels paresseux, voir detail/lazy.hpp), signe uniforme,
// puis l'arrondi correct de μ + σ Z. Tirage rapide (wsl::fast) : ziggourat.
//
// Densité, répartition, survie, quantile : z en double-double, erfc et son
// développement au premier ordre pour l'argument double-double, un seul arrondi
// final ; quantile par erfc_inv puis un pas de Newton en double-double. Erreurs
// mesurées dans tests/test_random_precision.cpp.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
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

namespace wsl
{
    namespace detail
    {
        // z = (x − μ)/σ en double-double. ±∞ si |z| > 2^600, 0 si |z| < 2^-600
        // (sans effet sur aucune des fonctions de la loi).
        constexpr sd::double_double standardize(double x, double mu, double sigma) noexcept
        {
            if (!is_finite(x))
                return {x, 0.0};                                            // ±∞, NaN
            auto [n, n_error] = sd::two_sum(x, -mu);
            int doubled = 0;
            if (!is_finite(n))                                              // sur les moitiés
            {
                auto const [h, h_error] = sd::two_sum(sd::fmul(x, 0.5), -sd::fmul(mu, 0.5));
                n = h;
                n_error = h_error;
                doubled = 1;
            }
            if (n == 0)
                return {0.0, 0.0};
            bool const negative = n < 0;
            double const a = negative ? -n : n, a_error = negative ? -n_error : n_error;
            int const e = sd::exponent_of(a) + doubled - sd::exponent_of(sigma);
            if (e > 600)
                return {negative ? -sd::infinity : sd::infinity, 0.0};
            if (e < -600)
                return {0.0, 0.0};
            scaled q = normalized(a, a_error) / sigma;
            q.e += doubled;
            auto const [hi, lo] = unscaled(q);
            return negative ? sd::double_double{-hi, -lo} : sd::double_double{hi, lo};
        }

        // erfc(w)/2, w = w_hi + w_lo, arrondi une fois. erfc(w_hi + w_lo) =
        // erfc(w_hi) (1 − w_lo (2/√π) / F(w_hi)) au premier ordre, F(x) = e^{x²}
        // erfc(x) (le second ordre est sous 2^-90 relativement).
        constexpr double half_erfc(double w_hi, double w_lo) noexcept
        {
            if (!(w_hi == w_hi))
                return w_hi;
            if (w_hi >= 28.0)                                               // erfc(28) < 2^-1100
                return 0.0;
            if (w_hi <= -6.0)                                               // 1 − erfc(6)/2 s'arrondit à 1
                return 1.0;
            double const aw = magnitude(w_hi);
            if (aw < 0.25)
            {
                if (aw < 0x1p-900)
                    return 0.5;
                auto const [e_hi, e_lo] = sd::erf_small(w_hi);
                double const derivative = sd::fmul(sd::fmul(w_lo, sd::two_over_sqrt_pi_hi), 1.0 - sd::fmul(w_hi, w_hi));
                auto const [s, s_error] = sd::two_sum(0.5, -sd::fmul(e_hi, 0.5));
                return s + (s_error - sd::fmul(e_lo + derivative, 0.5));
            }
            double const shift = w_hi > 0 ? w_lo : -w_lo;                   // erfc(aw + shift)
            auto const [hi, lo, e] = sd::erfc_scaled(aw);
            auto const [f_hi, f_lo] = sd::erfc_factor(aw);
            double const correction = sd::fmul(shift, sd::two_over_sqrt_pi_hi) / f_hi;
            scaled const half{hi, lo - sd::fmul(hi + lo, correction), e - 1};   // hi + lo : la paire n'est pas normalisée
            if (w_hi > 0)
                return rounded(half);
            auto const [c_hi, c_lo] = unscaled(half);                       // >= 2^-57 : normal
            auto const [s, s_error] = sd::two_sum(1.0, -c_hi);
            return s + (s_error - c_lo);
        }

        // Φ^{-1}(p), p dans (0, 1), en double-double : w = erfc_inv(2q), q = min(p,
        // 1 − p) (1 − p exact pour p >= 1/2), corrigé d'un pas de Newton sur
        // erfc(w)/2 = q avec résidu précis, puis Φ^{-1}(q) = −√2 w.
        constexpr sd::double_double standard_quantile(double p) noexcept
        {
            bool const upper = p > 0.5;
            double const q = upper ? 1.0 - p : p;
            if (q == 0.5)
                return {0.0, 0.0};
            double const w = wsl::special::erfc_inv(2.0 * q);              // > 0
            double delta = 0;
            if (w < 0.25)                                                   // q > 0,36 : résidu (1/2 − q) − erf(w)/2
            {
                auto const [e_hi, e_lo] = sd::erf_small(w);
                double const r = ((0.5 - q) - sd::fmul(e_hi, 0.5)) - sd::fmul(e_lo, 0.5);   // la première différence est exacte
                delta = sd::fmul(sd::fmul(r, 2.0 * half_sqrt_pi), wsl::special::exp(sd::fmul(w, w)));
            }
            else                                                            // ρ = erfc(w)/(2q) − 1, δ = (√π/2) F ρ/(1 + ρ)
            {
                auto const [hi, lo, e] = sd::erfc_scaled(w);
                auto const [u_hi, u_lo] = unscaled(scaled{hi, lo, e - 1} / q);
                auto const [s, s_error] = sd::two_sum(u_hi, -1.0);
                double const rho = s + (s_error + u_lo);
                auto const [f_hi, f_lo] = sd::erfc_factor(w);
                delta = sd::fmul(sd::fmul(half_sqrt_pi, f_hi), rho) / (1.0 + rho);
            }
            auto const [w_hi, w_lo] = sd::fast_two_sum(w, delta);
            auto const [m, m_error] = sd::two_product(w_hi, sqrt2_hi);
            auto const [x_hi, x_lo] = sd::fast_two_sum(m, m_error + (sd::fmul(w_hi, sqrt2_lo) + sd::fmul(w_lo, sqrt2_hi)));
            return upper ? sd::double_double{x_hi, x_lo} : sd::double_double{-x_hi, -x_lo};
        }

        // μ ± t arrondi une fois, t > 0 à l'échelle (au-delà des doubles s'il le faut).
        constexpr double affine_scaled(double mu, bool negative, scaled t) noexcept
        {
            if (mu == 0)
                return negative ? -rounded(t) : rounded(t);
            // Près du débordement, calcul à l'échelle 2^-k (exacte dans les deux sens :
            // le résultat est arrondi une fois, à grande échelle).
            int const k = std::max(t.e, sd::exponent_of(magnitude(mu))) - 1020;
            if (k > 0)
            {
                t.e -= k;
                mu = sd::scale_by_power_of_two(mu, -k);
            }
            auto const [t_hi, t_lo] = unscaled(t);
            auto const [s, s_error] = sd::two_sum(mu, negative ? -t_hi : t_hi);
            double const y = s + (s_error + (negative ? -t_lo : t_lo));
            return k > 0 ? sd::scale_by_power_of_two(y, k) : y;
        }

        // μ + σ x arrondi une fois, x en double-double non nul.
        constexpr double affine(double mu, double sigma, double x_hi, double x_lo) noexcept
        {
            bool const negative = x_hi < 0;
            return affine_scaled(mu, negative, normalized(negative ? -x_hi : x_hi, negative ? -x_lo : x_lo) * normalized(sigma, 0.0));
        }
    }

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
    class normal
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
                throw std::domain_error("wsl::normal : il faut mean fini et sd > 0, fini");
        }

        constexpr T sd() const noexcept { return sd_; }
        constexpr interval<T> support() const noexcept
        {
            return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()};
        }

        constexpr T mean() const noexcept { return mean_; }
        constexpr T median() const noexcept { return mean_; }
        constexpr T mode() const noexcept { return mean_; }
        constexpr T variance() const noexcept { return static_cast<T>(double(sd_) * double(sd_)); }
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

        friend constexpr T pdf(normal const& X, T x) noexcept
        {
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
        friend constexpr T log_pdf(normal const& X, T x) noexcept
        {
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
        friend constexpr T cdf(normal const& X, T x) noexcept
        {
            auto const [w_hi, w_lo] = X.scaled_argument(x);
            return static_cast<T>(detail::half_erfc(-w_hi, -w_lo));
        }

        // 1 − Φ(z) = erfc(z/√2)/2, sans annulation.
        friend constexpr T ccdf(normal const& X, T x) noexcept
        {
            auto const [w_hi, w_lo] = X.scaled_argument(x);
            return static_cast<T>(detail::half_erfc(w_hi, w_lo));
        }

        // μ + σ Φ^{-1}(p).
        friend constexpr T quantile(normal const& X, double p)
        {
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return -std::numeric_limits<T>::infinity();
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            auto const [x_hi, x_lo] = detail::standard_quantile(p);
            if (x_hi == 0)
                return X.mean_;
            return static_cast<T>(detail::affine(X.mean_, X.sd_, x_hi, x_lo));
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(normal const& X, G& g)
        {
            detail::lazy_real v = detail::exact_half_normal(g);
            bool const negative = (random_bits(g) >> 63) != 0;
            return detail::round_lazy<T>(X.mean_, X.sd_, 1.0, negative, v, g);
        }

        // Ziggourat (approché, bit-exact) : μ + σ Z, deux arrondis.
        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(normal const& X, G& g, fast_t)
        {
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
