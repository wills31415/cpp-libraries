#pragma once

// Décomposition de Loader (2000), « Fast and accurate computation of binomial
// probabilities » : les masses de Poisson et binomiales (et les densités gamma
// et bêta) s'écrivent avec deux termes petits et bien conditionnés,
//
//   δ(y)      = lgamma(y + 1) − [(y + 1/2) log y − y + log(2π)/2]   (erreur de Stirling)
//   bd0(x, M) = x log(x/M) + M − x >= 0                              (déviance)
//
// par exemple P(N = k) = e^{−δ(k) − bd0(k, μ)} / √(2πk) pour N de Poisson. Le
// calcul direct, exp(k log μ − μ − lgamma(k + 1)), retranche des termes qui
// peuvent valoir des milliers de fois le résultat : l'erreur relative du
// résultat est l'erreur absolue de l'exposant. Ici, δ et bd0 sont calculés en
// double-double à 2^-60 près relativement (ou mieux).

#include <wsl/special/detail/extended.hpp>

namespace wsl::special::detail::extended
{
    // δ(y) = lgamma(y) − [(y − 1/2) log y − y + log(2π)/2], y > 0 (la même
    // quantité que ci-dessus, écrite avec lgamma(y) = lgamma(y + 1) − log y).
    // - y >= 16 : série de Stirling Σ B_2k / (2k (2k − 1) y^(2k−1)), le terme de
    //   tête 1/(12 y) en double-double, le reste (moins de 2^-20 de lui) en double ;
    // - sinon : par récurrence (ci-dessous).
    constexpr sd::double_double stirling_error(double s_hi, double s_lo) noexcept;

    constexpr sd::double_double stirling_error(double y) noexcept
    {
        if (y >= 16.0)
        {
            double const z = 1.0 / y;
            auto const [p, p_error] = sd::two_product(z, y);
            double const z_lo = ((1.0 - p) - p_error) / y;                  // 1/y = z + z_lo
            double const z2 = z * z;
            constexpr int terms = sizeof(sd::stirling) / sizeof(sd::stirling[0]);
            double rest = sd::stirling[terms - 1];
            for (int n = terms - 2; n >= 1; --n)
                rest = sd::stirling[n] + sd::fmul(z2, rest);
            auto const [m, m_error] = sd::two_product(z, one_twelfth_hi);
            double const tail = (m_error + (sd::fmul(z, one_twelfth_lo) + sd::fmul(z_lo, one_twelfth_hi)))
                              + sd::fmul(sd::fmul(z, z2), rest);
            auto const [hi, lo] = sd::fast_two_sum(m, tail);
            return {hi, lo};
        }
        // y < 16 : récurrence jusqu'à s = y + n dans [16, 17), lgamma(y) = lgamma(s) − log P,
        // P = y (y + 1) … (y + n − 1) à l'échelle (exact à 2^-100 près, quel que soit y) :
        //     δ(y) = δ(s) + (s − 1/2) log s − (y − 1/2) log y − log P − n.
        // (lgamma_positive de wsl::special, précise à 2^-58 près, ne suffit pas ici.)
        int const n = 16 - static_cast<int>(y);
        auto const [s_hi, s_lo] = sd::two_sum(y, static_cast<double>(n));
        scaled product = normalized(y, 0.0);
        for (int i = 1; i < n; ++i)
        {
            auto const [f_hi, f_lo] = sd::two_sum(y, static_cast<double>(i));
            product = product * normalized(f_hi, f_lo);
        }
        sd::double_double result = stirling_error(s_hi, s_lo);
        auto const [sh_hi, sh_lo] = sd::two_sum(s_hi, -0.5);
        result = sum(result, mul({sh_hi, sh_lo + s_lo}, log(normalized(s_hi, s_lo))));
        auto const [yh_hi, yh_lo] = sd::two_sum(y, -0.5);
        auto const ly = mul({yh_hi, yh_lo}, log_precise(y));
        result = sum(result, {-ly.hi, -ly.lo});
        auto const lp = log(product);
        result = sum(result, {-lp.hi, -lp.lo});
        return sum(result, {-static_cast<double>(n), 0.0});
    }

    // δ(s_hi + s_lo), s_hi >= 16 : δ(s_hi) + δ'(s_hi) s_lo, δ'(s) ≈ −1/(12 s²) (le terme
    // suivant, 1/(120 s⁴), et le second ordre sont négligeables).
    constexpr sd::double_double stirling_error(double s_hi, double s_lo) noexcept
    {
        auto const [hi, lo] = stirling_error(s_hi);
        if (s_lo == 0)
            return {hi, lo};
        auto const [h, l] = sd::fast_two_sum(hi, lo - s_lo / sd::fmul(12.0, sd::fmul(s_hi, s_hi)));
        return {h, l};
    }

    // bd0(x, M) = x log(x/M) + M − x, x >= 0, M = m_hi + m_lo > 0, en
    // double-double. Homogène de degré 1 : mis à l'échelle s'il le faut.
    // - |x − M| < (x + M)/10 : v = (x − M)/(x + M), bd0 = (x − M) v + 2x Σ v^(2j+1)/(2j+1)
    //   (sans annulation ; le terme de tête en double-double, la série en double,
    //   moins de 1/30 du tout) ;
    // - sinon : x log(x/M), log du rapport en double-double (|log| > 0,2 : précis
    //   relativement), puis M − x ; l'annulation coûte au plus un facteur 10.
    constexpr sd::double_double deviance(double x, double m_hi, double m_lo) noexcept
    {
        if (x == 0)
            return {m_hi, m_lo};
        int shift = 0;
        double const largest = x > m_hi ? x : m_hi;
        if (largest > 0x1p900)
            shift = -200;
        else if (largest < 0x1p-900)
            shift = 200;
        if (shift != 0)
        {
            x = sd::scale_by_power_of_two(x, shift);
            m_hi = sd::scale_by_power_of_two(m_hi, shift);
            m_lo = sd::scale_by_power_of_two(m_lo, shift);
        }
        double hi = 0, lo = 0;
        if (magnitude(x - m_hi) < sd::fmul(x + m_hi, 0.1))
        {
            auto const [d1, e1] = sd::two_sum(x, -m_hi);
            auto const [d_hi, d_lo] = sd::fast_two_sum(d1, e1 - m_lo);       // x − M
            auto const [s1, f1] = sd::two_sum(x, m_hi);
            auto const [s_hi, s_lo] = sd::fast_two_sum(s1, f1 + m_lo);       // x + M
            double const v = d_hi / s_hi;                                   // v = (x − M)/(x + M) en double-double
            auto const [q, q_error] = sd::two_product(v, s_hi);
            double const v_lo = ((((d_hi - q) - q_error) + d_lo) - sd::fmul(v, s_lo)) / s_hi;
            if (d_hi == 0)
                return {0.0, 0.0};
            sd::double_double const vv{v, v_lo};
            sd::double_double const head = mul({d_hi, d_lo}, vv);           // (x − M) v
            if (magnitude(v) < 0x1p-400)                                    // série sous 2^-800 de la tête
            {
                hi = head.hi;
                lo = head.lo;
            }
            else
            {
                // 2x v w S(w), w = v², S(w) = Σ_{j>=1} w^(j−1)/(2j + 1), en double-double :
                // la série pèse jusqu'à 1/15 de bd0, qui peut valoir des centaines.
                sd::double_double const w = mul(vv, vv);
                sd::double_double series = reciprocal(29.0);
                for (int j = 13; j >= 1; --j)                               // w^14 < 2^-93
                    series = sum(reciprocal(static_cast<double>(2 * j + 1)), mul(w, series));
                sd::double_double const tail = mul(mul({sd::fmul(2.0, x), 0.0}, vv), mul(w, series));
                auto const [h, l] = sum(head, tail);
                hi = h;
                lo = l;
            }
        }
        else
        {
            scaled const ratio = normalized(x, 0.0) / normalized(m_hi, m_lo);
            auto const [l_hi, l_lo] = log(ratio);
            double xl_hi = 0, xl_lo = 0;
            double const approx = magnitude(sd::fmul(x, l_hi));
            if (approx >= 0x1p-900)
            {
                auto const [p, p_error] = balanced_product(x, l_hi);
                xl_hi = p;
                xl_lo = p_error + sd::fmul(x, l_lo);
            }
            else
                xl_hi = sd::fmul(x, l_hi);                                  // négligeable devant M
            auto const [s1, e1] = sd::two_sum(xl_hi, m_hi);
            auto const [s2, e2] = sd::two_sum(s1, -x);
            auto const [h, l] = sd::fast_two_sum(s2, (e1 + e2) + (xl_lo + m_lo));
            hi = h;
            lo = l;
        }
        if (shift != 0)
        {
            hi = sd::scale_by_power_of_two(hi, -shift);
            lo = sd::scale_by_power_of_two(lo, -shift);
        }
        return {hi, lo};
    }

    // lgamma(y) = (y − 1/2) log y − y + log(2π)/2 + δ(y), en double-double, y > 0.
    constexpr sd::double_double log_gamma(double y) noexcept
    {
        auto const [h_hi, h_lo] = sd::two_sum(y, -0.5);
        sd::double_double result = mul({h_hi, h_lo}, log_precise(y));
        result = sum(result, {-y, 0.0});
        result = sum(result, {sd::half_log_two_pi_hi, sd::half_log_two_pi_lo});
        return sum(result, stirling_error(y));
    }

    // ψ(s), s > 0, à 2^-30 près relativement : pour corriger lgamma d'un argument
    // double-double au premier ordre (récurrence jusqu'à s >= 6, puis
    // développement asymptotique).
    constexpr double digamma_approximate(double s) noexcept
    {
        double shift = 0;
        while (s < 6)
        {
            shift -= 1 / s;
            s += 1;
        }
        double const z = 1 / sd::fmul(s, s);
        return shift + ((wsl::special::log(s) - 0.5 / s)
                        - sd::fmul(z, 1.0 / 12 - sd::fmul(z, 1.0 / 120 - sd::fmul(z, 1.0 / 252))));
    }

    // lgamma(s_hi + s_lo) = lgamma(s_hi) + ψ(s_hi) s_lo au premier ordre.
    constexpr sd::double_double log_gamma(double s_hi, double s_lo) noexcept
    {
        auto const [g_hi, g_lo] = log_gamma(s_hi);
        if (s_lo == 0)
            return {g_hi, g_lo};
        auto const [hi, lo] = sd::fast_two_sum(g_hi, g_lo + sd::fmul(digamma_approximate(s_hi), s_lo));
        return {hi, lo};
    }

    // log B(a, b) = lgamma(a) + lgamma(b) − lgamma(a + b), en double-double, a, b > 0.
    // Si le plus grand q dépasse 16, lgamma(q) − lgamma(p + q) se compense : on
    // l'écrit δ(q) − δ(p + q) − (q − 1/2) log1p(p/q) − p log(p + q) + p.
    constexpr sd::double_double log_beta(double a, double b) noexcept
    {
        double const p = a < b ? a : b, q = a < b ? b : a;
        auto const [s_hi, s_lo] = sd::two_sum(p, q);
        sd::double_double result = log_gamma(p);
        if (q < 16)
        {
            result = sum(result, log_gamma(q));
            auto const g = log_gamma(s_hi, s_lo);
            return sum(result, {-g.hi, -g.lo});
        }
        scaled const ratio = normalized(p, 0.0) / q;
        auto const [r_hi, r_lo] = unscaled(ratio);
        auto const l1 = log1p_positive(r_hi, r_lo);
        auto const ls = log(normalized(s_hi, s_lo));
        auto const [h_hi, h_lo] = sd::two_sum(q, -0.5);
        result = sum(result, stirling_error(q));
        auto const ds = stirling_error(s_hi, s_lo);
        result = sum(result, {-ds.hi, -ds.lo});
        auto const t1 = mul({h_hi, h_lo}, l1);
        result = sum(result, {-t1.hi, -t1.lo});
        auto const t2 = mul({p, 0.0}, ls);
        result = sum(result, {-t2.hi, -t2.lo});
        return sum(result, {p, 0.0});
    }
}
