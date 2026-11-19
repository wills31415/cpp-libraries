#pragma once

// Moments des lois de Weibull et de Fréchet : lgamma(1 + ε) et c² [Γ(1 + 2ε) − Γ(1 + ε)²]
// (ε = 1/k, ou −1/α), la différence sans compensation : D = lgamma(1 + 2ε) − 2 lgamma(1 + ε)
// par la série de lgamma(1 + a) pour |ε| <= 1/8 (les termes linéaires s'annulent).

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/special.hpp>

#include <limits>

namespace wsl::detail
{
    // lgamma(1 + ε), ε = e + e_lo > −1, en double-double.
    constexpr sd::double_double lgamma_one_plus(double e, double e_lo) noexcept
    {
        namespace fp = sd;
        if (e <= 0.25 && e >= -0.25)
        {
            auto const [g_hi, g_lo] = fp::lgamma1p_small(e);               // + ψ(1 + e) e_lo
            return sum({g_hi, g_lo}, {fp::fmul(digamma_approximate(1.0 + e), e_lo), 0.0});
        }
        auto const [s_hi, s_lo0] = fp::two_sum(1.0, e);
        return log_gamma(s_hi, s_lo0 + e_lo);
    }

    // c² [Γ(1 + 2ε) − Γ(1 + ε)²], log c = log_c, arrondi une fois.
    constexpr double gamma_variance(double e, double e_lo, sd::double_double log_c) noexcept
    {
        namespace fp = sd;
        auto const g1 = lgamma_one_plus(e, e_lo);
        auto const base = sum({fp::fmul(g1.hi, 2.0), fp::fmul(g1.lo, 2.0)}, {fp::fmul(log_c.hi, 2.0), fp::fmul(log_c.lo, 2.0)});   // log(c² Γ(1+ε)²)
        fp::double_double d{};                                              // D = lgamma(1 + 2ε) − 2 lgamma(1 + ε)
        if (e <= 0.125 && e >= -0.125)
        {
            // Σ_{j>=2} l_j (2^j − 2) ε^j, l_j les coefficients de lgamma(1 + a) (Horner en double-double)
            constexpr int terms = sizeof(fp::lgamma1p_series) / sizeof(fp::lgamma1p_series[0]);
            fp::double_double s{0.0, 0.0};
            for (int j = terms; j >= 2; --j)
            {
                double const w = fp::scale_by_power_of_two(1.0, j) - 2.0;  // 2^j − 2, exact
                fp::double_double const c = mul({fp::lgamma1p_series[j - 1][0], fp::lgamma1p_series[j - 1][1]}, {w, 0.0});
                s = sum(c, mul(s, {e, e_lo}));
            }
            d = mul(mul(s, {e, e_lo}), {e, e_lo});                          // s ε² (s porte déjà ε^{j−2})
        }
        else
        {
            auto const g2 = lgamma_one_plus(fp::fmul(e, 2.0), fp::fmul(e_lo, 2.0));
            d = sum(g2, {fp::fmul(g1.hi, -2.0), fp::fmul(g1.lo, -2.0)});
        }
        fp::double_double m{};
        if (d.hi <= 1)
            m = expm1_double_double(d.hi, d.lo);
        else
            m = sum(unscaled(exp_scaled(d.hi > 1400 ? 1400 : d.hi, d.lo)), {-1.0, 0.0});
        if (base.hi > 1400)
            return std::numeric_limits<double>::infinity();
        return rounded(normalized(m.hi, m.lo) * exp_scaled(base.hi, base.lo));
    }
}
