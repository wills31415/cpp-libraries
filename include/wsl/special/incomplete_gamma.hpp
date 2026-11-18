#pragma once

// Fonctions gamma incomplètes régularisées, bit-exactes (voir detail/fp.hpp pour
// les conditions), constexpr :
//
//   gamma_p(a, x) = γ(a, x)/Γ(a) = (1/Γ(a)) ∫₀ˣ t^{a−1} e^{−t} dt
//   gamma_q(a, x) = Γ(a, x)/Γ(a) = 1 − gamma_p(a, x)
//
// a > 0, x >= 0 (sinon NaN). gamma_p(a, 0) = 0, gamma_p(a, +∞) = 1 ; a = +∞ : P = 0.
// Chacune précise relativement, même minuscule : le calcul rend la plus petite
// des deux, l'autre vaut 1 − elle.
//
// Régimes (tout en double-double, un seul arrondi final) :
// - a >= 50 et |η| <= 1 (x/a dans [0,30 ; 2,36]) : développement uniforme de
//   Temme (DLMF 8.12), Q = erfc(w)/2 + R, w = ±√bd0(a, x), le facteur e^{−bd0}
//   mis à part (pas de sous-débordement), coefficients engendrés exactement ;
// - sinon, x < a + 1 : série P = D Σ xⁿ/((a + 1)…(a + n)), D = xᵃ e^{−x}/Γ(a + 1)
//   par la décomposition de Loader (detail/loader.hpp) ;
//   a < 1/4 de plus : Q = −expm1(a log x − lgamma(1 + a)) − (xᵃ/Γ(a + 1)) a Σ (−x)ⁿ/(n!(a + n)),
//   lgamma(1 + a) en série (Q tend vers 0 avec a : 1 − P ne suffirait pas) ;
// - sinon : fraction continue de Legendre pour Q (Lentz).

#include <wsl/special/detail/extended.hpp>
#include <wsl/special/detail/fp.hpp>
#include <wsl/special/detail/incomplete_gamma_tables.hpp>
#include <wsl/special/detail/loader.hpp>
#include <wsl/special/elementary.hpp>
#include <wsl/special/error_function.hpp>


namespace wsl::special
{
    namespace detail
    {
        namespace xd = extended;

        // P ou Q, la plus petite (valeur nulle : hi = 0).
        struct incomplete_gamma_result
        {
            xd::scaled value;
            bool upper;                                                     // value est Q (sinon P)
        };

        inline constexpr xd::scaled scaled_zero{0.0, 0.0, 0};

        constexpr xd::scaled to_scaled(double_double v) noexcept
        {
            return v.hi > 0 ? xd::normalized(v.hi, v.lo) : scaled_zero;
        }

        // a/b en double-double.
        constexpr double_double quotient_dd(double_double a, double_double b) noexcept
        {
            double const q = a.hi / b.hi;
            auto const [p, p_error] = two_product(q, b.hi);
            auto const [hi, lo] = fast_two_sum(q, ((((a.hi - p) - p_error) + a.lo) - fmul(q, b.lo)) / b.hi);
            return {hi, lo};
        }

        // D = xᵃ e^{−x}/Γ(a + 1) = e^{−δ(a) − bd0(a, x)} / √(2πa), à l'échelle.
        constexpr xd::scaled gamma_prefix(double a, double x) noexcept
        {
            auto const [d_hi, d_lo] = xd::stirling_error(a);
            auto const [b_hi, b_lo] = xd::deviance(a, x, 0.0);
            auto const [e, e_error] = two_sum(-d_hi, -b_hi);
            if (!(e > -1500))                                               // < 2^-1080 même divisé par √(2πa)
                return scaled_zero;
            auto const [hi, lo] = fast_two_sum(e, e_error - (d_lo + b_lo));
            return xd::exp_scaled(hi, lo) * xd::normalized(xd::inv_sqrt_2pi_hi, xd::inv_sqrt_2pi_lo)
                 / xd::sqrt(xd::normalized(a, 0.0));
        }

        // P = D Σ_{n>=0} tₙ, t₀ = 1, tₙ = tₙ₋₁ x/(a + n).
        constexpr xd::scaled gamma_p_series(double a, double x) noexcept
        {
            xd::scaled const prefix = gamma_prefix(a, x);
            if (prefix.hi == 0)
                return prefix;
            double_double term{1.0, 0.0}, sum{1.0, 0.0};
            for (int n = 1; n < 1000000; ++n)
            {
                auto const [s_hi, s_lo] = two_sum(a, static_cast<double>(n));
                term = xd::mul(term, quotient_dd({x, 0.0}, {s_hi, s_lo}));
                sum = xd::sum(sum, term);
                if (term.hi < fmul(sum.hi, 0x1p-80))
                    break;
            }
            return prefix * xd::normalized(sum.hi, sum.lo);
        }

        // Q = a D / f, f = b₀ + a₁/(b₁ + a₂/(b₂ + …)), b_n = x + 2n + 1 − a,
        // a_n = −n (n − a) (Legendre), par l'algorithme de Lentz.
        constexpr xd::scaled gamma_q_fraction(double a, double x) noexcept
        {
            xd::scaled const prefix = gamma_prefix(a, x);
            if (prefix.hi == 0)
                return prefix;
            auto const [xa_hi, xa_lo] = two_sum(x, -a);
            double_double const x_minus_a{xa_hi, xa_lo};
            double_double f = xd::sum(x_minus_a, {1.0, 0.0});
            double_double c = f, d{0.0, 0.0};
            for (int n = 1; n < 1000000; ++n)
            {
                auto const [na_hi, na_lo] = two_sum(static_cast<double>(n), -a);
                double_double const an = xd::mul({-static_cast<double>(n), 0.0}, {na_hi, na_lo});
                double_double const bn = xd::sum(x_minus_a, {static_cast<double>(2 * n + 1), 0.0});
                d = xd::sum(bn, xd::mul(an, d));
                c = xd::sum(bn, quotient_dd(an, c));
                d = quotient_dd({1.0, 0.0}, d);
                double_double const delta = xd::mul(c, d);
                f = xd::mul(f, delta);
                double const change = (delta.hi - 1.0) + delta.lo;
                if (xd::magnitude(change) < 0x1p-80)
                    break;
            }
            return prefix * xd::normalized(a, 0.0) / xd::normalized(f.hi, f.lo);
        }

        // lgamma(1 + a), 0 < a <= 1/4, en double-double (série en a).
        constexpr double_double lgamma1p_small(double a) noexcept
        {
            constexpr int terms = sizeof(lgamma1p_series) / sizeof(lgamma1p_series[0]);
            double_double s{lgamma1p_series[terms - 1][0], lgamma1p_series[terms - 1][1]};
            for (int k = terms - 2; k >= 0; --k)
                s = xd::sum({lgamma1p_series[k][0], lgamma1p_series[k][1]}, xd::mul(s, {a, 0.0}));
            return xd::mul(s, {a, 0.0});
        }

        // Q pour a < 1/4, 0 < x < a + 1 :
        //     Q = −expm1(L) − u a S₁,  L = a log x − lgamma(1 + a),  u = e^L,
        //     S₁ = Σ_{n>=1} (−x)ⁿ / (n! (a + n)).
        constexpr xd::scaled gamma_q_small_a(double a, double x) noexcept
        {
            double_double const lg = lgamma1p_small(a);
            double_double const lx = xd::log_precise(x);
            double_double const l = xd::sum(xd::mul({a, 0.0}, lx), {-lg.hi, -lg.lo});
            double_double first{1.0, 0.0}, u{0.0, 0.0};
            if (l.hi < -40)                                                 // e^L < 2^-57 : 1 − e^L
            {
                u = xd::unscaled(xd::exp_scaled(l.hi, l.lo));
                first = xd::sum({1.0, 0.0}, {-u.hi, -u.lo});
            }
            else
            {
                auto const [m_hi, m_lo] = xd::expm1_double_double(l.hi, l.lo);
                first = {-m_hi, -m_lo};
                u = xd::sum({m_hi, m_lo}, {1.0, 0.0});
            }
            double_double power{1.0, 0.0}, s1{0.0, 0.0};
            for (int n = 1; n < 200; ++n)
            {
                power = xd::mul(power, quotient_dd({-x, 0.0}, {static_cast<double>(n), 0.0}));   // (−x)ⁿ/n!
                auto const [s_hi, s_lo] = two_sum(a, static_cast<double>(n));
                double_double const term = quotient_dd(power, {s_hi, s_lo});
                s1 = xd::sum(s1, term);
                if (xd::magnitude(term.hi) < fmul(xd::magnitude(s1.hi), 0x1p-80))
                    break;
            }
            double_double const second = xd::mul(xd::mul(u, {a, 0.0}), s1);
            return to_scaled(xd::sum(first, {-second.hi, -second.lo}));
        }

        // F(w) = e^{w²} erfc(w), w = w_hi + w_lo >= 1/4, en double-double :
        // F(w_hi) + F'(w_hi) w_lo, F' = 2wF − 2/√π ; au-delà de 28, développement
        // asymptotique (1/(√π w)) Σ (−1)ⁿ (2n − 1)!! / (2w²)ⁿ.
        constexpr double_double erfc_factor_dd(double w_hi, double w_lo) noexcept
        {
            double_double f{};
            if (w_hi < 28.0)
            {
                auto const [h, l] = erfc_factor(w_hi);
                f = {h, l};
            }
            else
            {
                double const z = 1.0 / fmul(2.0, fmul(w_hi, w_hi));
                double series = 0;
                for (int n = 9; n >= 1; --n)                                // −z + 3z² − 15z³ …
                    series = fmul(-static_cast<double>(2 * n - 1), fmul(z, 1.0 + series));
                double_double const lead = quotient_dd({xd::inv_sqrt_pi, 0x1.1ae3a914fed80p-57}, {w_hi, 0.0});
                f = xd::mul(lead, {1.0, series});
            }
            if (w_lo == 0)
                return f;
            double const slope = fmul(fmul(2.0, w_hi), f.hi) - two_over_sqrt_pi_hi;
            return xd::sum(f, {fmul(slope, w_lo), 0.0});
        }

        // Temme : a >= 50, bd0(a, x) = b <= a/2.
        constexpr incomplete_gamma_result gamma_temme(double a, double x, double_double b) noexcept
        {
            bool const above = x > a;
            double_double w{0.0, 0.0}, eta{0.0, 0.0};
            if (b.hi > 0)
            {
                xd::scaled const bs = xd::normalized(b.hi, b.lo);
                w = xd::unscaled(xd::sqrt(bs));
                eta = xd::unscaled(xd::sqrt(bs * xd::normalized(2.0, 0.0) / a));
                if (!above)
                {
                    w = {-w.hi, -w.lo};
                    eta = {-eta.hi, -eta.lo};
                }
            }
            // T = Σ_k c_k(η) a^{−k} : c₀ en double-double, le reste en double.
            double_double c0{temme_c0[temme_degree[0]][0], temme_c0[temme_degree[0]][1]};
            for (int n = temme_degree[0] - 1; n >= 0; --n)
                c0 = xd::sum({temme_c0[n][0], temme_c0[n][1]}, xd::mul(c0, eta));
            double const z = 1.0 / a;
            double rest = 0;
            for (int k = temme_terms - 1; k >= 1; --k)
            {
                auto const& row = temme_c[k - 1];
                double ck = row[temme_degree[k]];
                for (int n = temme_degree[k] - 1; n >= 0; --n)
                    ck = row[n] + fmul(eta.hi, ck);
                rest = fmul(z, ck + rest);
            }
            double_double const t = xd::sum(c0, {rest, 0.0});
            double_double const scale = xd::unscaled(xd::normalized(xd::inv_sqrt_2pi_hi, xd::inv_sqrt_2pi_lo)
                                                     / xd::sqrt(xd::normalized(a, 0.0)));
            double_double const r = xd::mul(t, scale);                      // R = e^{−b} r
            if (xd::magnitude(w.hi) < 0.25)                                 // Q = 1/2 − erf(w)/2 + R, près de 1/2
            {
                double_double e{0.0, 0.0};
                if (xd::magnitude(w.hi) >= 0x1p-900)
                {
                    auto const [e_hi, e_lo] = erf_small(w.hi);
                    double const slope = fmul(two_over_sqrt_pi_hi, 1.0 - fmul(w.hi, w.hi));
                    e = {e_hi, e_lo + fmul(slope, w.lo)};
                }
                double_double const ex = xd::unscaled(xd::exp_scaled(-b.hi, -b.lo));
                double_double q = xd::sum({0.5, 0.0}, {fmul(e.hi, -0.5), fmul(e.lo, -0.5)});
                q = xd::sum(q, xd::mul(ex, r));
                return {to_scaled(q), true};
            }
            if (b.hi > 1500)                                                // e^{−bd0} F/2 < 2^-2000 : nul (et hors
                return {scaled_zero, above};                                // du domaine de exp_scaled)
            double_double const aw = above ? w : double_double{-w.hi, -w.lo};
            double_double const f = erfc_factor_dd(aw.hi, aw.lo);
            double_double const bracket = xd::sum({fmul(f.hi, 0.5), fmul(f.lo, 0.5)}, above ? r : double_double{-r.hi, -r.lo});
            return {xd::exp_scaled(-b.hi, -b.lo) * xd::normalized(bracket.hi, bracket.lo), above};
        }

        // a > 0, x > 0 finis.
        constexpr incomplete_gamma_result incomplete_gamma(double a, double x) noexcept
        {
            if (a >= temme_from)
            {
                auto const [b_hi, b_lo] = xd::deviance(a, x, 0.0);
                if (b_hi <= fmul(a, 0.5))
                    return gamma_temme(a, x, {b_hi, b_lo});
            }
            if (x < a + 1)
            {
                xd::scaled const p = gamma_p_series(a, x);
                if (a < 0.25)
                {
                    xd::scaled const q = gamma_q_small_a(a, x);
                    bool const q_smaller = q.e < p.e || (q.e == p.e && q.hi < p.hi);
                    return q_smaller ? incomplete_gamma_result{q, true} : incomplete_gamma_result{p, false};
                }
                return {p, false};
            }
            return {gamma_q_fraction(a, x), true};
        }

        // 1 − v, v dans [0, 1].
        constexpr double one_minus(xd::scaled const& v) noexcept
        {
            if (v.hi == 0)
                return 1.0;
            auto const [hi, lo] = xd::unscaled(v);
            auto const [s, s_error] = two_sum(1.0, -hi);
            return s + (s_error - lo);
        }

        // NaN, domaine, bords : vrai si result est fixé.
        constexpr bool incomplete_gamma_special(double a, double x, bool upper, double& result) noexcept
        {
            if (!(a == a) || !(x == x))
                result = opaque(a + x);
            else if (!(a > 0) || !(x >= 0) || (a == infinity && x == infinity))
                result = quiet_nan;
            else if (x == 0 || a == infinity)
                result = upper ? 1.0 : 0.0;
            else if (x == infinity)
                result = upper ? 0.0 : 1.0;
            else
                return false;
            return true;
        }
    }

    namespace detail
    {
        // ---------------------------------------------------------------------
        // Inverses : x tel que P(a, x) = p (ou Q(a, x) = q).
        //
        // On résout sur le côté minoritaire (P si la cible t <= 1/2, sinon Q avec
        // t = 1 − p, exact), en log : y = log x, L(y) = log F(e^y), pente
        // dL/dy = x f(x)/F(x) (f la densité, x f = a D). Newton, encadré (bissection
        // en y si le pas sort de l'encadrement), depuis Wilson–Hilferty (a >= 1) ou
        // les approximations des queues. Le résidu L − log t est exact à 2^-64 près
        // en absolu ; l'erreur sur x en est le quotient par la pente, à 2^-62 près
        // sauf pour a petit dans la queue basse (pente a) : là, log P est calculé
        // par a log x − lgamma(1 + a) + log1p(a S₁), d'erreur absolue en a 2^-64.
        // Le dernier pas, x (1 − Δ), en double-double, arrondi une fois.
        // ---------------------------------------------------------------------

        constexpr xd::scaled one_minus_scaled(xd::scaled const& v) noexcept
        {
            if (v.hi == 0)
                return xd::normalized(1.0, 0.0);
            auto const [hi, lo] = xd::unscaled(v);
            return to_scaled(xd::sum({1.0, 0.0}, {-hi, -lo}));
        }

        // log P(a, x), a < 1/4, 0 < x < a + 1 : a log x − lgamma(1 + a) + log1p(a S₁),
        // S₁ = Σ_{n>=1} (−x)ⁿ/(n! (a + n)) ; erreur absolue proportionnelle à a.
        constexpr double_double log_p_small_a(double a, double x) noexcept
        {
            double_double power{1.0, 0.0}, s1{0.0, 0.0};
            for (int n = 1; n < 200; ++n)
            {
                power = xd::mul(power, quotient_dd({-x, 0.0}, {static_cast<double>(n), 0.0}));
                auto const [s_hi, s_lo] = two_sum(a, static_cast<double>(n));
                double_double const term = quotient_dd(power, {s_hi, s_lo});
                s1 = xd::sum(s1, term);
                if (xd::magnitude(term.hi) < fmul(xd::magnitude(s1.hi), 0x1p-80))
                    break;
            }
            double_double const u = xd::mul({a, 0.0}, s1);                  // a S₁, dans (−1, 0)
            auto const [l1_hi, l1_lo] = xd::log1p_dd(u.hi);
            double_double const l1{l1_hi, l1_lo + u.lo / (1.0 + u.hi)};
            double_double const lg = lgamma1p_small(a);
            double_double result = xd::mul({a, 0.0}, xd::log_precise(x));
            result = xd::sum(result, {-lg.hi, -lg.lo});
            return xd::sum(result, l1);
        }

        struct side_log
        {
            double_double log_value;                                        // log F(x), F = P ou Q
            double slope;                                                   // x f(x) / F(x) > 0
        };

        constexpr side_log gamma_side(double a, double x, bool upper) noexcept
        {
            xd::scaled const density = xd::normalized(a, 0.0) * gamma_prefix(a, x);   // x f(x) = a D
            if (!upper && a < 0.25 && x < a + 1)
            {
                double_double const l = log_p_small_a(a, x);
                xd::scaled const value = xd::exp_scaled(l.hi, l.lo);
                return {l, xd::rounded(density / value)};
            }
            auto const [v, v_upper] = incomplete_gamma(a, x);
            xd::scaled const value = v_upper == upper ? v : one_minus_scaled(v);
            if (value.hi == 0)
                return {{-infinity, 0.0}, infinity};
            return {xd::log(value), density.hi == 0 ? 0.0 : xd::rounded(density / value)};
        }

        // Point de départ (double).
        constexpr double gamma_inverse_start(double a, double t, bool upper) noexcept
        {
            double const z = fmul(xd::sqrt2_hi, erfc_inv(fmul(2.0, t)));      // Φ⁻¹ de la probabilité côté P, au signe près
            double const zp = upper ? z : -z;
            if (a >= 1)                                                     // Wilson–Hilferty
            {
                double const c = 1.0 / fmul(9.0, a);
                double const w = (1.0 - c) + fmul(zp, sqrt(c));
                if (w > 0)
                {
                    double const x = fmul(a, fmul(fmul(w, w), w));
                    if (x > 0 && x < 0x1p1023)
                        return x;
                }
            }
            if (!upper)                                                     // P ≈ xᵃ/Γ(a + 1)
            {
                double const l = (log(t) + lgamma(a + 1.0)) / a;
                return l < -744 ? 0x1p-1074 : l > 709 ? 0x1p1023 : exp(l);
            }
            double x = a > 1 ? a : 1.0;                                     // Q ≈ x^{a−1} e^{−x}/Γ(a)
            for (int i = 0; i < 6; ++i)
            {
                double const next = (-log(t) - lgamma(a)) + fmul(a - 1.0, log(x));
                x = next > 0x1p-1074 ? next : 0x1p-1074;
            }
            return x;
        }

        // En double-double (le dernier pas) : pour qui doit encore diviser le résultat.
        constexpr double_double gamma_inverse_dd(double a, double t, bool upper) noexcept
        {
            auto const [lt_hi, lt_lo] = xd::log_precise(t);
            if (upper && a < 1)                                             // Q(a, x) <= Q(a, 2^-1074) < t :
            {                                                               // la solution est sous le plus petit double
                auto const [l, slope] = gamma_side(a, 0x1p-1074, true);
                if ((l.hi - lt_hi) + (l.lo - lt_lo) < 0)
                    return {0.0, 0.0};
            }
            double x = gamma_inverse_start(a, t, upper);
            double low = 0, high = infinity;                                // encadrement de x
            for (int iteration = 0; iteration < 200; ++iteration)
            {
                auto const [l, slope] = gamma_side(a, x, upper);
                double const g = (l.hi - lt_hi) + (l.lo - lt_lo);           // L − log t
                double const s = upper ? -slope : slope;                    // dL/d(log x)
                if (g == 0)
                    return {x, 0.0};
                bool const too_far = (g > 0) == (s > 0);                    // x trop grand
                if (too_far)
                    high = x;
                else
                    low = x;
                double const delta = s != 0 && xd::is_finite(s) ? g / s : (too_far ? 4.0 : -4.0);
                // Le pas se fait sur x (log x n'a pas assez de résolution) : x e^{−Δ}.
                double next = xd::magnitude(delta) < 0x1p-20 ? x - fmul(x, delta) : fmul(x, exp(-delta));
                if (xd::magnitude(delta) < 0x1p-56 || next == x)            // dernier pas, en double-double
                {
                    auto const [m, m_error] = two_product(x, delta);
                    auto const [r, r_error] = two_sum(x, -m);
                    auto const [hi, lo] = fast_two_sum(r, r_error - m_error);
                    return {hi, lo};
                }
                if (!(next > low && next < high))                           // bissection, en log
                {
                    if (low > 0 && high < infinity)
                        next = exp(fmul(log(low) + log(high), 0.5));
                    else
                        next = fmul(x, too_far ? 0x1.2c1599b2a5cb3p-6 : 0x1.b4c902e273a58p+5);   // e^∓4
                }
                if (!(next > 0))                                            // sous le plus petit sous-normal
                {
                    if (high <= 0x1p-1074)
                        return {0.0, 0.0};
                    next = 0x1p-1074;
                }
                if (!(next < infinity))
                {
                    if (low >= 0x1.fffffffffffffp+1023)
                        return {infinity, 0.0};
                    next = 0x1.fffffffffffffp+1023;
                }
                x = next;
            }
            return {x, 0.0};
        }

        constexpr double gamma_inverse(double a, double t, bool upper) noexcept
        {
            auto const [hi, lo] = gamma_inverse_dd(a, t, upper);
            return hi + lo;
        }

        constexpr bool gamma_inverse_special(double a, double p, double& result) noexcept
        {
            if (!(a == a) || !(p == p))
                result = opaque(a + p);
            else if (!(a > 0) || a == infinity || !(p >= 0 && p <= 1))
                result = quiet_nan;
            else if (p == 0)
                result = 0.0;
            else if (p == 1)
                result = infinity;
            else
                return false;
            return true;
        }
    }

    constexpr double gamma_p(double a, double x) noexcept
    {
        double special = 0;
        if (detail::incomplete_gamma_special(a, x, false, special))
            return special;
        auto const [value, upper] = detail::incomplete_gamma(a, x);
        return upper ? detail::one_minus(value) : detail::extended::rounded(value);
    }

    constexpr double gamma_q(double a, double x) noexcept
    {
        double special = 0;
        if (detail::incomplete_gamma_special(a, x, true, special))
            return special;
        auto const [value, upper] = detail::incomplete_gamma(a, x);
        return upper ? detail::extended::rounded(value) : detail::one_minus(value);
    }

    constexpr float gamma_p(float a, float x) noexcept { return detail::to_float(gamma_p(static_cast<double>(a), static_cast<double>(x))); }
    constexpr float gamma_q(float a, float x) noexcept { return detail::to_float(gamma_q(static_cast<double>(a), static_cast<double>(x))); }

    // Inverses : gamma_p_inv(a, p) = x tel que P(a, x) = p ; gamma_q_inv(a, q) = x tel que
    // Q(a, x) = q. p (q) dans [0, 1] : 0 en 0 (en 1), +∞ en 1 (en 0) ; sinon NaN.
    constexpr double gamma_p_inv(double a, double p) noexcept
    {
        double special = 0;
        if (detail::gamma_inverse_special(a, p, special))
            return special;
        return p <= 0.5 ? detail::gamma_inverse(a, p, false) : detail::gamma_inverse(a, 1.0 - p, true);
    }

    constexpr double gamma_q_inv(double a, double q) noexcept
    {
        if (!(a == a) || !(q == q))
            return detail::opaque(a + q);
        if (!(a > 0) || a == detail::infinity || !(q >= 0 && q <= 1))
            return detail::quiet_nan;
        if (q == 1)
            return 0.0;
        if (q == 0)
            return detail::infinity;
        return q <= 0.5 ? detail::gamma_inverse(a, q, true) : detail::gamma_inverse(a, 1.0 - q, false);
    }

    constexpr float gamma_p_inv(float a, float p) noexcept { return detail::to_float(gamma_p_inv(static_cast<double>(a), static_cast<double>(p))); }
    constexpr float gamma_q_inv(float a, float q) noexcept { return detail::to_float(gamma_q_inv(static_cast<double>(a), static_cast<double>(q))); }
}
