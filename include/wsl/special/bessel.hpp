#pragma once

// Fonction de Bessel modifiée de première espèce I_ν(x), bit-exacte (voir
// detail/fp.hpp pour les conditions), constexpr, nommée comme std::cyl_bessel_i :
//
//   cyl_bessel_i(ν, x) = Σ_k (x/2)^{2k + ν} / (k! Γ(k + ν + 1))
//
// ν >= 0 et x >= 0 ; ν entier : aussi ν < 0 (I_{−n} = I_n) et x < 0 (I_n(−x) = (−1)ⁿ I_n(x)) ;
// NaN sinon. I_0(0) = 1, I_ν(0) = 0 (ν > 0), I_ν(+∞) = +∞, I_∞(x) = 0.
//
// Régimes (tout en double-double, un arrondi final) :
// - ν >= 50 : développement uniforme de Debye (DLMF 10.41.3),
//   I_ν(νz) = e^{νη} Σ_k u_k(t)/ν^k / (√(2πν) (1 + z²)^{1/4}), t = 1/√(1 + z²),
//   η = √(1 + z²) + log(z/(1 + √(1 + z²))), 28 polynômes u_k engendrés exactement ;
// - sinon, x >= max(60, ν²/2) : développement de Hankel,
//   e^x/√(2πx) Σ_k (−1)^k a_k(ν)/x^k (au plus 36 termes ; le terme exponentiellement
//   petit omis est sous e^{−2x} relativement) ;
// - sinon : série entière, termes positifs (pas de compensation), la somme à l'échelle
//   (elle dépasse 2^1024 pour x proche de 1250).
// detail::log_cyl_bessel_i rend log I_ν(x) en double-double (lois de von Mises).

#include <wsl/special/detail/bessel_tables.hpp>
#include <wsl/special/detail/extended.hpp>
#include <wsl/special/detail/fp.hpp>
#include <wsl/special/detail/loader.hpp>
#include <wsl/special/elementary.hpp>
#include <wsl/special/gamma.hpp>
#include <wsl/special/incomplete_gamma.hpp>

#include <limits>

namespace wsl::special
{
    namespace detail
    {
        namespace xd = extended;

        // I_ν(x) = e^{log_prefactor} · sum.
        struct bessel_parts
        {
            double_double log_prefactor;
            xd::scaled sum;
        };

        // log Γ(1 + ν), ν >= 0, en double-double.
        constexpr double_double lgamma_one_plus(double nu) noexcept
        {
            if (nu == 0)
                return {0.0, 0.0};
            if (nu <= 0.25)
                return lgamma1p_small(nu);
            auto const [s_hi, s_lo] = two_sum(1.0, nu);
            return xd::log_gamma(s_hi, s_lo);
        }

        // Série : (x/2)^ν/Γ(ν + 1) · Σ_k (x²/4)^k / (k! (ν + 1)_k).
        constexpr bessel_parts bessel_series(double nu, double x) noexcept
        {
            auto const lx = xd::sum(xd::log_precise(x), {-ln2_hi, -ln2_lo});
            double_double pref = nu == 0 ? double_double{0.0, 0.0} : xd::mul({nu, 0.0}, lx);
            auto const lg = lgamma_one_plus(nu);
            pref = xd::sum(pref, {-lg.hi, -lg.lo});
            if (x < 0x1p-400)                                               // x²/4 (ν + 1) < 2^-800 : un terme
                return {pref, xd::normalized(1.0, 0.0)};
            auto const [p, p_error] = two_product(x, x);
            double_double const q{fmul(p, 0.25), fmul(p_error, 0.25)};      // x²/4, exact
            double_double s{1.0, 0.0}, t{1.0, 0.0};
            int e = 0;                                                      // s, t à l'échelle 2^e
            for (int k = 1;; ++k)
            {
                auto const [a, a_error] = two_sum(nu, double(k));           // k (k + ν)
                double_double const d = xd::mul({a, a_error}, {double(k), 0.0});
                t = quotient_dd(xd::mul(t, q), d);
                s = xd::sum(s, t);
                if (s.hi > 0x1p600)
                {
                    s = {scale_by_power_of_two(s.hi, -600), scale_by_power_of_two(s.lo, -600)};
                    t = {scale_by_power_of_two(t.hi, -600), scale_by_power_of_two(t.lo, -600)};
                    e += 600;
                }
                if (t.hi <= fmul(s.hi, 0x1p-112) && double(k) * double(k) >= q.hi)
                    break;
            }
            xd::scaled sum = xd::normalized(s.hi, s.lo);
            sum.e += e;
            return {pref, sum};
        }

        // Hankel : e^x/√(2πx) Σ_k (−1)^k a_k(ν)/x^k, a_k = a_{k−1} (4ν² − (2k − 1)²)/(8k).
        constexpr bessel_parts bessel_hankel(double nu, double x) noexcept
        {
            auto const lx = xd::log_precise(x);
            double_double pref = xd::sum({x, 0.0}, {-xd::half_log_2pi_hi, -xd::half_log_2pi_lo});
            pref = xd::sum(pref, {fmul(lx.hi, -0.5), fmul(lx.lo, -0.5)});
            auto const [m, m_error] = two_product(nu, nu);
            double_double const mu4{fmul(m, 4.0), fmul(m_error, 4.0)};
            double_double s{1.0, 0.0}, t{1.0, 0.0};
            for (int k = 1; k < 200; ++k)
            {
                double const odd = double(2 * k - 1);
                auto const f = xd::sum(mu4, {-fmul(odd, odd), 0.0});
                auto const [d, d_error] = two_product(double(8 * k), x);
                t = quotient_dd(xd::mul(t, f), {d, d_error});
                t = {-t.hi, -t.lo};
                s = xd::sum(s, t);
                if (xd::magnitude(t.hi) <= 0x1p-112)
                    break;
            }
            return {pref, xd::normalized(s.hi, s.lo)};
        }

        // Debye : I_ν(νz) = e^{νη} √t / √(2πν) · Σ_k u_k(t) ν^{−k}.
        constexpr bessel_parts bessel_uniform(double nu, double x) noexcept
        {
            double_double const z = quotient_dd({x, 0.0}, {nu, 0.0});
            double_double const w = xd::sum({1.0, 0.0}, xd::mul(z, z));     // 1 + z²
            xd::scaled const ws = xd::normalized(w.hi, w.lo);
            double_double const sw = xd::unscaled(xd::sqrt(ws));            // √(1 + z²)
            double_double const r = quotient_dd(z, xd::sum({1.0, 0.0}, sw));
            double_double const eta = xd::sum(sw, xd::log(xd::normalized(r.hi, r.lo)));
            double_double pref = xd::mul({nu, 0.0}, eta);
            pref = xd::sum(pref, {-xd::half_log_2pi_hi, -xd::half_log_2pi_lo});
            auto const ln = xd::log_precise(nu);
            pref = xd::sum(pref, {fmul(ln.hi, -0.5), fmul(ln.lo, -0.5)});
            auto const lw = xd::log(ws);                                    // ½ log t = −¼ log(1 + z²)
            pref = xd::sum(pref, {fmul(lw.hi, -0.25), fmul(lw.lo, -0.25)});
            double_double const t = quotient_dd({1.0, 0.0}, sw);
            double_double const tau = quotient_dd({1.0, 0.0}, w);           // t²
            double_double const step = quotient_dd(t, {nu, 0.0});           // t/ν
            // Σ_k (t/ν)^k P_k(t²), P_k(τ) = Σ_j c_kj τ^j ; Horner en k, puis en τ.
            double_double s{0.0, 0.0};
            for (int k = bessel_debye_terms - 1; k >= 0; --k)
            {
                int const base = k * (k + 1) / 2;
                double_double pk{bessel_debye[base + k][0], bessel_debye[base + k][1]};
                for (int j = k - 1; j >= 0; --j)
                    pk = xd::sum(xd::mul(pk, tau), {bessel_debye[base + j][0], bessel_debye[base + j][1]});
                s = xd::sum(xd::mul(s, step), pk);
            }
            return {pref, xd::normalized(s.hi, s.lo)};
        }

        // ν >= 0, x > 0 finis.
        constexpr bessel_parts bessel_parts_of(double nu, double x) noexcept
        {
            if (nu >= bessel_debye_from)
                return bessel_uniform(nu, x);
            double const half_square = fmul(fmul(nu, nu), 0.5);
            if (x >= bessel_hankel_from && x >= half_square)
                return bessel_hankel(nu, x);
            return bessel_series(nu, x);
        }

        // I_ν(x), ν >= 0, x > 0 finis, arrondi une fois.
        constexpr double bessel_positive(double nu, double x) noexcept
        {
            if (nu >= bessel_debye_from && x < 0x1p-900)                    // (x/2)^ν/Γ(ν + 1) < 2^-45000
                return 0;
            auto const [pref, sum] = bessel_parts_of(nu, x);
            if (pref.hi > 1400)
                return std::numeric_limits<double>::infinity();
            if (pref.hi < -1500 && sum.e < 600)
                return 0;
            return xd::rounded(xd::exp_scaled(pref.hi, pref.lo) * sum);
        }

        // log I_ν(x), ν >= 0, x > 0 finis, en double-double.
        constexpr double_double log_cyl_bessel_i(double nu, double x) noexcept
        {
            auto const [pref, sum] = bessel_parts_of(nu, x);
            return xd::sum(pref, xd::log(sum));
        }

        // ν entier impair ? (pair au-delà de 2^53.)
        constexpr bool is_odd_integer(double v) noexcept
        {
            double const a = xd::magnitude(v);
            if (a >= 0x1p53)
                return false;
            double const h = fmul(a, 0.5);
            return xd::floor_nonnegative(h) != h;
        }
    }

    constexpr double cyl_bessel_i(double nu, double x) noexcept
    {
        if (!(x == x))
            return detail::opaque(x + x);
        if (!(nu == nu))
            return detail::opaque(nu + nu);
        bool const integer = detail::is_integer(nu);
        if (nu < 0)
        {
            if (!integer)
                return std::numeric_limits<double>::quiet_NaN();
            nu = -nu;
        }
        bool negative = false;
        if (x < 0)
        {
            if (!integer)
                return std::numeric_limits<double>::quiet_NaN();
            negative = detail::is_odd_integer(nu);
            x = -x;
        }
        double r = 0;
        if (x == 0)
            r = nu == 0 ? 1.0 : 0.0;
        else if (nu == std::numeric_limits<double>::infinity())
            r = 0;
        else if (x == std::numeric_limits<double>::infinity())
            r = std::numeric_limits<double>::infinity();
        else
            r = detail::bessel_positive(nu, x);
        return negative ? -r : r;
    }

    constexpr float cyl_bessel_i(float nu, float x) noexcept
    {
        return detail::to_float(cyl_bessel_i(static_cast<double>(nu), static_cast<double>(x)));
    }
}
