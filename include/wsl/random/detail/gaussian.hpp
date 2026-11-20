#pragma once

// Aides gaussiennes des lois : forme réduite (x − μ)/σ en double-double ou à l'échelle,
// erfc(w)/2, Φ⁻¹, μ ± t arrondi une fois, Φ(w) à l'échelle.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/special.hpp>
#include <algorithm>

namespace wsl::detail
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

    // |x − μ|/σ à l'échelle, sans saturation, et son signe ; x fini.
    struct standardized_scaled
    {
        bool negative;
        scaled z;                                                       // z.hi == 0 si x = μ
    };

    constexpr standardized_scaled standardize_scaled(double x, double mu, double sigma) noexcept
    {
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
            return {false, {0.0, 0.0, 0}};
        bool const negative = n < 0;
        scaled q = normalized(negative ? -n : n, negative ? -n_error : n_error) / sigma;
        q.e += doubled;
        return {negative, q};
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

    // Φ(w), w en double-double, à l'échelle (0 sous 2^-1100 environ).
    constexpr scaled normal_lower_scaled(sd::double_double w) noexcept
    {
        auto const [x_hi, x_error] = sd::two_product(w.hi, inv_sqrt2_hi);   // erfc(−w/√2)/2
        double const x_lo = x_error + (sd::fmul(w.lo, inv_sqrt2_hi) + sd::fmul(w.hi, inv_sqrt2_lo));
        auto const [y_hi, y_lo] = sd::fast_two_sum(-x_hi, -x_lo);      // y = −w/√2
        if (y_hi >= 27.0)                                              // erfc(27) < 2^-1060
            return {0.0, 0.0, 0};
        double const ay = magnitude(y_hi);
        if (ay < 0.25)
        {
            if (ay < 0x1p-900)
                return {1.0, 0.0, -1};
            auto const [e_hi, e_lo] = sd::erf_small(y_hi);
            double const derivative = sd::fmul(sd::fmul(y_lo, sd::two_over_sqrt_pi_hi), 1.0 - sd::fmul(y_hi, y_hi));
            auto const [s, s_error] = sd::two_sum(0.5, -sd::fmul(e_hi, 0.5));
            auto const [h, l] = sd::fast_two_sum(s, s_error - sd::fmul(e_lo + derivative, 0.5));
            return normalized(h, l);
        }
        if (y_hi <= -6.5)                                              // 1 − erfc(6,5)/2, erfc(6,5) < 2^-64
        {
            if (y_hi <= -27.0)
                return {1.0, 0.0, 0};
            auto const [hi, lo, e] = sd::erfc_scaled(ay);
            scaled const half{hi, lo, e - 1};
            auto const [c_hi, c_lo] = unscaled(half);
            auto const [s, s_error] = sd::two_sum(1.0, -c_hi);
            auto const [h, l] = sd::fast_two_sum(s, s_error - c_lo);
            return normalized(h, l);
        }
        double const shift = y_hi > 0 ? y_lo : -y_lo;
        auto const [hi, lo, e] = sd::erfc_scaled(ay);
        auto const [f_hi, f_lo] = sd::erfc_factor(ay);
        (void)f_lo;
        double const correction = sd::fmul(shift, sd::two_over_sqrt_pi_hi) / f_hi;
        scaled const half = normalized(hi, lo - sd::fmul(hi + lo, correction), e - 1);
        if (y_hi > 0)
            return half;
        auto const [c_hi, c_lo] = unscaled(half);
        auto const [s, s_error] = sd::two_sum(1.0, -c_hi);
        auto const [h, l] = sd::fast_two_sum(s, s_error - c_lo);
        return normalized(h, l);
    }
}
