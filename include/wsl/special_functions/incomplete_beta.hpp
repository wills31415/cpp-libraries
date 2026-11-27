#pragma once

#include <wsl/modes.hpp>

// Fonction bêta et bêta incomplète régularisée, bit-exactes (voir detail/fp.hpp
// pour les conditions), constexpr :
//
//   beta(a, b)  = B(a, b) = Γ(a) Γ(b) / Γ(a + b)      lbeta(a, b) = log B(a, b)
//   (lbeta, comme lgamma près de ses zéros, est précise en valeur absolue — 2^-70 —
//   près de la courbe B(a, b) = 1 ; exacte pour a = 1 ou b = 1)
//   ibeta(a, b, x)  = Iₓ(a, b) = B(a, b)⁻¹ ∫₀ˣ t^{a−1} (1 − t)^{b−1} dt
//   ibetac(a, b, x) = 1 − Iₓ(a, b) = I_{1−x}(b, a)
//
// a, b > 0 et x dans [0, 1] (sinon NaN) ; Iₓ = 0 en x = 0, 1 en x = 1.
// Chacune précise relativement, même minuscule : le calcul rend la plus petite
// des deux, l'autre vaut 1 − elle.
//
// Préfixe xᵃ (1 − x)ᵇ / B(a, b) : par la décomposition de Loader si a, b > 2
// (binomiale généralisée), sinon directement, log B en double-double sans
// compensation ; puis la fraction continue de Iₓ(a, b) (Lentz, double-double)
// pour x < (a + 1)/(a + b + 2), celle de I_{1−x}(b, a) sinon. Si la valeur obtenue
// dépasse 1/2 avec un paramètre petit, l'autre côté lui est proportionnel (1 − v
// perdrait ses chiffres) : formule propre, comme pour la gamma incomplète.

#include <wsl/special_functions/detail/extended.hpp>
#include <wsl/floating_point_arithmetic/detail/primitives.hpp>
#include <wsl/special_functions/detail/loader.hpp>
#include <wsl/special_functions/incomplete_gamma.hpp>

#include <wsl/elementary_functions.hpp>
#include <wsl/elementary_functions/detail/pinned_portable.hpp>

// Couches du dessous : primitives flottantes et fonctions élémentaires.
namespace wsl::inline WSL_DEFAULTS::special_functions::detail
{
    using namespace floating_point_arithmetic::detail;
    using namespace elementary_functions::detail;
}
// Transition (refonte, étape 4a) : les fonctions élémentaires, trouvées avant les surcharges
// de mdarray pour les tableaux ; l'étape 4b les appellera qualifiées, avec leur mode.
namespace wsl::inline WSL_DEFAULTS::special_functions
{
    using elementary_functions::detail::pinned_portable::sqrt;
    using elementary_functions::detail::pinned_portable::cbrt;
    using elementary_functions::detail::pinned_portable::hypot;
    using elementary_functions::detail::pinned_portable::exp;
    using elementary_functions::detail::pinned_portable::exp2;
    using elementary_functions::detail::pinned_portable::expm1;
    using elementary_functions::detail::pinned_portable::log;
    using elementary_functions::detail::pinned_portable::log2;
    using elementary_functions::detail::pinned_portable::log10;
    using elementary_functions::detail::pinned_portable::log1p;
    using elementary_functions::detail::pinned_portable::pow;
    using elementary_functions::detail::pinned_portable::sin;
    using elementary_functions::detail::pinned_portable::cos;
    using elementary_functions::detail::pinned_portable::tan;
    using elementary_functions::detail::pinned_portable::atan;
    using elementary_functions::detail::pinned_portable::atan2;
}

namespace wsl::inline WSL_DEFAULTS::special_functions
{
    namespace detail
    {
        // xᵃ yᵇ / B(a, b), y = 1 − x (en double-double), à l'échelle.
        constexpr xd::scaled beta_prefix(double a, double b, double x, double_double y) noexcept
        {
            if (a > 2 && b > 2)                                             // Loader : √(ab/(2πn)) e^{δ(n) − δ(a) − δ(b) − bd0(a, nx) − bd0(b, ny)}
            {
                auto const [n_hi, n_lo] = two_sum(a, b);
                double_double const nx = xd::mul({n_hi, n_lo}, {x, 0.0});
                double_double const ny = xd::mul({n_hi, n_lo}, y);
                double_double e = xd::stirling_error(n_hi, n_lo);
                auto const da = xd::stirling_error(a), db = xd::stirling_error(b);
                auto const ba = xd::deviance(a, nx.hi, nx.lo), bb = xd::deviance(b, ny.hi, ny.lo);
                e = xd::sum(e, {-da.hi, -da.lo});
                e = xd::sum(e, {-db.hi, -db.lo});
                e = xd::sum(e, {-ba.hi, -ba.lo});
                e = xd::sum(e, {-bb.hi, -bb.lo});
                if (!(e.hi > -1500))
                    return scaled_zero;
                xd::scaled const factor = xd::sqrt(xd::normalized(a, 0.0) * xd::normalized(b, 0.0)
                                                   / xd::normalized(n_hi, n_lo))
                                        * xd::normalized(xd::inv_sqrt_2pi_hi, xd::inv_sqrt_2pi_lo);
                return xd::exp_scaled(e.hi, e.lo) * factor;
            }
            // a log x + b log y − log B(a, b)
            double_double l = xd::mul({a, 0.0}, xd::log_precise(x));
            double_double const ly = xd::log(xd::normalized(y.hi, y.lo));
            l = xd::sum(l, xd::mul({b, 0.0}, ly));
            double_double const lb = xd::log_beta(a, b);
            l = xd::sum(l, {-lb.hi, -lb.lo});
            if (!(l.hi > -1500))
                return scaled_zero;
            return xd::exp_scaled(l.hi, l.lo);
        }

        // log(xᵃ yᵇ / B(a, b)), x et y = 1 − x en double-double tous deux : Loader si a, b > 2
        // (nx et ny en double-double), sinon a log x + b log y − log B(a, b).
        constexpr double_double beta_log_prefix(double a, double b, double_double x, double_double y) noexcept
        {
            if (a > 2 && b > 2)
            {
                auto const [n_hi, n_lo] = two_sum(a, b);
                double_double const nx = xd::mul({n_hi, n_lo}, x), ny = xd::mul({n_hi, n_lo}, y);
                double_double e = xd::stirling_error(n_hi, n_lo);
                auto const da = xd::stirling_error(a), db = xd::stirling_error(b);
                auto const ba = xd::deviance(a, nx.hi, nx.lo), bb = xd::deviance(b, ny.hi, ny.lo);
                e = xd::sum(e, {-da.hi, -da.lo});
                e = xd::sum(e, {-db.hi, -db.lo});
                e = xd::sum(e, {-ba.hi, -ba.lo});
                e = xd::sum(e, {-bb.hi, -bb.lo});
                // + ½ (log a + log b − log n) − ½ log 2π
                auto const la = xd::log_precise(a), lb = xd::log_precise(b);
                auto const ln = xd::log(xd::normalized(n_hi, n_lo));
                double_double f = xd::sum(la, lb);
                f = xd::sum(f, {-ln.hi, -ln.lo});
                e = xd::sum(e, {fmul(f.hi, 0.5), fmul(f.lo, 0.5)});
                return xd::sum(e, {-half_log_two_pi_hi, -half_log_two_pi_lo});
            }
            double_double l = xd::mul({a, 0.0}, xd::log(xd::normalized(x.hi, x.lo)));
            l = xd::sum(l, xd::mul({b, 0.0}, xd::log(xd::normalized(y.hi, y.lo))));
            double_double const lb = xd::log_beta(a, b);
            return xd::sum(l, {-lb.hi, -lb.lo});
        }

        // h tel que Iₓ(a, b) = préfixe · h / a : fraction continue
        //     1/(1 + d₁/(1 + d₂/(1 + …))),  d_{2m+1} = −(a + m)(a + b + m) x / ((a + 2m)(a + 2m + 1)),
        //                                   d_{2m}   = m (b − m) x / ((a + 2m − 1)(a + 2m)),
        // par l'algorithme de Lentz, en double-double (x = x_hi + x_lo).
        constexpr double_double beta_fraction(double a, double b, double_double x) noexcept
        {
            constexpr double tiny = 0x1p-900;
            auto const fix = [](double_double v) { return xd::magnitude(v.hi) < tiny ? double_double{tiny, 0.0} : v; };
            auto const [ab_hi, ab_lo] = two_sum(a, b);
            double_double const ab{ab_hi, ab_lo};
            auto const [a1_hi, a1_lo] = two_sum(a, 1.0);
            double_double c{1.0, 0.0};
            double_double d = fix(xd::sum({1.0, 0.0}, quotient_dd(xd::mul(ab, x), {-a1_hi, -a1_lo})));
            d = quotient_dd({1.0, 0.0}, d);
            double_double h = d;
            for (int m = 1; m < 10000000; ++m)
            {
                double const md = static_cast<double>(m);
                auto const [p_hi, p_lo] = two_sum(a, static_cast<double>(2 * m - 1));           // a + 2m − 1
                auto const [q_hi, q_lo] = two_sum(a, static_cast<double>(2 * m));             // a + 2m
                auto const [r_hi, r_lo] = two_sum(a, static_cast<double>(2 * m + 1));           // a + 2m + 1
                auto const [bm_hi, bm_lo] = two_sum(b, -md);                // b − m
                double_double const even = quotient_dd(xd::mul(xd::mul({md, 0.0}, {bm_hi, bm_lo}), x),
                                                       xd::mul({p_hi, p_lo}, {q_hi, q_lo}));
                d = quotient_dd({1.0, 0.0}, fix(xd::sum({1.0, 0.0}, xd::mul(even, d))));
                c = fix(xd::sum({1.0, 0.0}, quotient_dd(even, c)));
                h = xd::mul(h, xd::mul(d, c));
                auto const [am_hi, am_lo] = two_sum(a, md);                 // a + m
                double_double const abm = xd::sum(ab, {md, 0.0});           // a + b + m
                double_double const odd = quotient_dd(xd::mul(xd::mul({-am_hi, -am_lo}, abm), x),
                                                      xd::mul({q_hi, q_lo}, {r_hi, r_lo}));
                d = quotient_dd({1.0, 0.0}, fix(xd::sum({1.0, 0.0}, xd::mul(odd, d))));
                c = fix(xd::sum({1.0, 0.0}, quotient_dd(odd, c)));
                double_double const delta = xd::mul(d, c);
                h = xd::mul(h, delta);
                if (xd::magnitude((delta.hi - 1.0) + delta.lo) < 0x1p-80)
                    break;
            }
            return h;
        }

        // lgamma(t) − lgamma(t + s), 0 < s <= 1/4, t > 0, en double-double, précis
        // relativement à s : récurrence jusqu'à u = t + n >= 16 (termes log1p(s/(t + i))),
        // puis Stirling, lgamma(u) − lgamma(u + s) =
        //     −(u − 1/2) log1p(s/u) − s (log(u + s) − 1) − (δ(u + s) − δ(u)),
        // δ(u + s) − δ(u) terme à terme (le premier exactement : −s / (12 u (u + s))).
        constexpr double_double log_gamma_ratio(double s, double t) noexcept
        {
            double_double result{0.0, 0.0};
            int n = 0;
            if (t < 16)
            {
                n = 16 - static_cast<int>(t);
                for (int i = 0; i < n; ++i)                                 // + log1p(s/(t + i))
                {
                    auto const [f_hi, f_lo] = two_sum(t, static_cast<double>(i));
                    double_double const r = quotient_dd({s, 0.0}, {f_hi, f_lo});
                    result = xd::sum(result, xd::log1p_positive(r.hi, r.lo));
                }
            }
            auto const [u_hi, u_lo] = two_sum(t, static_cast<double>(n));
            double_double const u{u_hi, u_lo};
            double_double const r = quotient_dd({s, 0.0}, u);
            double_double const l1 = xd::log1p_positive(r.hi, r.lo);
            double_double const uh = xd::sum(u, {-0.5, 0.0});
            double_double term = xd::mul(uh, l1);
            result = xd::sum(result, {-term.hi, -term.lo});
            double_double const us = xd::sum(u, {s, 0.0});
            double_double const lus = xd::log(xd::normalized(us.hi, us.lo));
            term = xd::mul({s, 0.0}, xd::sum(lus, {-1.0, 0.0}));
            result = xd::sum(result, {-term.hi, -term.lo});
            // δ(u + s) − δ(u) = Σ_k c_k (u + s)^{1−2k} − c_k u^{1−2k}
            double_double const first = quotient_dd(xd::mul({s, 0.0}, {xd::one_twelfth_hi, xd::one_twelfth_lo}), xd::mul(u, us));
            double rest = 0;
            constexpr int terms = sizeof(stirling) / sizeof(stirling[0]);
            double const w = 1.0 / fmul(u.hi, u.hi);
            double power = 1.0 / u.hi;                                      // u^{1−2k}
            for (int k = 2; k <= terms; ++k)
            {
                power = fmul(power, w);
                double const m = static_cast<double>(1 - 2 * k);
                rest = rest + fmul(fmul(stirling[k - 1], power), expm1(fmul(m, log1p(r.hi))));
            }
            // −(δ(u + s) − δ(u)) = first − rest
            result = xd::sum(result, first);
            return xd::sum(result, {-rest, 0.0});
        }

        // 1 − I_z(s, t), s <= 1/4, z <= 3/4, t z <= 8 (z = z_hi + z_lo) : comme pour la gamma,
        //     1 − I = −expm1(L) − u s T,   L = s log z − log(s B(s, t)),  u = e^L,
        //     T = Σ_{j>=1} (1 − t)_j z^j / (j! (s + j)),
        // log(s B(s, t)) = lgamma(1 + s) + lgamma(t) − lgamma(s + t) : tout est proportionnel à s.
        constexpr xd::scaled beta_small_complement(double s, double t, double_double z) noexcept
        {
            double_double const lsb = xd::sum(lgamma1p_small(s), log_gamma_ratio(s, t));
            double_double const lz = xd::log(xd::normalized(z.hi, z.lo));
            double_double const l = xd::sum(xd::mul({s, 0.0}, lz), {-lsb.hi, -lsb.lo});
            double_double first{1.0, 0.0}, u{0.0, 0.0};
            if (l.hi < -40)
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
            double_double term{1.0, 0.0}, total{0.0, 0.0};
            double const peak = fmul(t, z.hi);
            for (int j = 1; j < 100000; ++j)
            {
                double const jd = static_cast<double>(j);
                auto const [c_hi, c_lo] = two_sum(jd, -t);                  // j − t
                term = xd::mul(term, quotient_dd(xd::mul({c_hi, c_lo}, z), {jd, 0.0}));
                auto const [d_hi, d_lo] = two_sum(s, jd);
                double_double const add = quotient_dd(term, {d_hi, d_lo});
                total = xd::sum(total, add);
                if (jd > peak && xd::magnitude(add.hi) <= fmul(xd::magnitude(total.hi), 0x1p-80))
                    break;
            }
            double_double const second = xd::mul(xd::mul(u, {s, 0.0}), total);
            return to_scaled(xd::sum(first, {-second.hi, -second.lo}));
        }

        struct incomplete_beta_result
        {
            xd::scaled value;
            bool complement;                                                // value est 1 − Iₓ (sinon Iₓ)
        };

        // a, b > 0 finis, 0 < x < 1.
        constexpr incomplete_beta_result incomplete_beta(double a, double b, double x) noexcept
        {
            auto const [y_hi, y_lo] = two_sum(1.0, -x);                     // 1 − x, exact
            xd::scaled const prefix = beta_prefix(a, b, x, {y_hi, y_lo});
            if (prefix.hi == 0)
            {
                // Préfixe nul : la plus petite des deux est nulle ; laquelle ? celle du côté de x.
                bool const lower = x < a / (a + b);
                return {scaled_zero, !lower};
            }
            // Fraction de Iₓ(a, b) en deçà du point de bascule, de I_{1−x}(b, a) au-delà ; si
            // la valeur obtenue dépasse 1/2 (un paramètre minuscule : l'autre côté lui est
            // proportionnel, et 1 − v perdrait ses chiffres), l'autre côté directement —
            // sa fraction converge encore, plus lentement.
            auto const lower = [&] { double_double const h = beta_fraction(a, b, {x, 0.0}); return prefix * xd::normalized(h.hi, h.lo) / a; };
            auto const upper = [&] { double_double const h = beta_fraction(b, a, {y_hi, y_lo}); return prefix * xd::normalized(h.hi, h.lo) / b; };
            bool const below = x < (a + 1.0) / ((a + b) + 2.0);
            xd::scaled const v = below ? lower() : upper();
            if (v.hi == 0 || v.e < -1)                                     // v < 1/2 (hi dans [1, 2))
                return {v, !below};
            // L'autre côté est 1 − I_z(s, t), (s, t, z) = (a, b, x) si v = Iₓ(a, b), (b, a, 1 − x) sinon.
            double const s_ = below ? a : b, t_ = below ? b : a;
            double_double const z = below ? double_double{x, 0.0} : double_double{y_hi, y_lo};
            if (s_ <= 0.25 && z.hi <= 0.75 && fmul(t_, z.hi) <= 8)
                return {beta_small_complement(s_, t_, z), below};
            return {v, !below};                                             // 1 − v : assez précis ici
        }

        // Iₓ(a, b) du côté minoritaire en x = x_hi + x_lo, y = 1 − x = y_hi + y_lo (en
        // double-double tous deux) : évaluée en la plus petite des parties hautes, puis
        // corrigée au premier ordre (dI/dx = xᵃ⁻¹yᵇ⁻¹/B, le second ordre est négligeable).
        constexpr incomplete_beta_result incomplete_beta_dd(double a, double b, double_double x, double_double y) noexcept
        {
            bool const in_x = x.hi <= y.hi;
            double const s = in_x ? a : b, t = in_x ? b : a;
            double_double const u = in_x ? x : y;
            auto [value, complement] = incomplete_beta(s, t, u.hi);         // I_u(s, t), au point u_hi
            if (u.lo != 0 && value.hi != 0)
            {
                auto const [w_hi, w_lo] = two_sum(1.0, -u.hi);
                xd::scaled const density = beta_prefix(s, t, u.hi, {w_hi, w_lo}) / xd::normalized(u.hi, 0.0)
                                         / xd::normalized(w_hi, w_lo);
                double const r = xd::rounded(density * xd::normalized(xd::magnitude(u.lo), 0.0) / value);
                double const signed_r = (u.lo > 0) != complement ? r : -r;  // dI/du = f, d(1 − I)/du = −f
                auto const [c_hi, c_lo] = two_sum(1.0, signed_r);
                value = value * xd::normalized(c_hi, c_lo);
            }
            return {value, in_x ? complement : !complement};              // Iₓ(a, b) = 1 − I_y(b, a)
        }

        constexpr bool incomplete_beta_special(double a, double b, double x, bool complement, double& result) noexcept
        {
            if (!(a == a) || !(b == b) || !(x == x))
                result = opaque(a + b + x);
            else if (!(a > 0) || !(b > 0) || !(x >= 0 && x <= 1) || (a == infinity && b == infinity))
                result = quiet_nan;
            else if (x == 0 || a == infinity)
                result = complement ? 1.0 : 0.0;
            else if (x == 1 || b == infinity)
                result = complement ? 0.0 : 1.0;
            else
                return false;
            return true;
        }
    }

    namespace detail
    {
        // ---------------------------------------------------------------------
        // Inverses : x tel que Iₓ(a, b) = p (ou 1 − Iₓ(a, b) = q).
        //
        // Variable u = x si la solution est <= 1/2 (p <= I_{1/2}), sinon u = 1 − x avec
        // les paramètres échangés (Iₓ(a, b) = 1 − I_{1−x}(b, a)) : u reste petit, donc
        // précis. Puis, comme pour la gamma : Newton encadré sur log F(u) du côté
        // minoritaire, pente u f(u)/F(u) = préfixe/((1 − u) F), pas sur u lui-même,
        // dernier pas en double-double ; queue basse d'un petit paramètre s :
        // log I = s log u − log(sB) + log1p(s T), d'erreur absolue proportionnelle à s.
        // ---------------------------------------------------------------------

        constexpr double lbeta_double(double a, double b) noexcept
        {
            auto const [hi, lo] = xd::log_beta(a, b);
            return hi + lo;
        }

        // I_{1/2}(a, b) et 1 − I_{1/2}(a, b), chacune calculée de son côté (pour situer la
        // solution sans former 1 − p, qui vaut 1 pour p minuscule).
        struct beta_halves
        {
            double lower, upper;
        };

        constexpr beta_halves ibeta_halves(double a, double b) noexcept
        {
            auto const [v, complement] = incomplete_beta(a, b, 0.5);
            double const direct = xd::rounded(v), other = one_minus(v);
            return complement ? beta_halves{other, direct} : beta_halves{direct, other};
        }

        // log Iᵤ(s, t), s <= 1/4, u <= 3/4, t u <= 8.
        constexpr double_double log_ibeta_small(double s, double t, double u) noexcept
        {
            double_double const lsb = xd::sum(lgamma1p_small(s), log_gamma_ratio(s, t));
            double_double term{1.0, 0.0}, total{0.0, 0.0};
            double const peak = fmul(t, u);
            for (int j = 1; j < 100000; ++j)
            {
                double const jd = static_cast<double>(j);
                auto const [c_hi, c_lo] = two_sum(jd, -t);
                term = xd::mul(term, quotient_dd(xd::mul({c_hi, c_lo}, {u, 0.0}), {jd, 0.0}));
                auto const [d_hi, d_lo] = two_sum(s, jd);
                double_double const add = quotient_dd(term, {d_hi, d_lo});
                total = xd::sum(total, add);
                if (jd > peak && xd::magnitude(add.hi) <= fmul(xd::magnitude(total.hi), 0x1p-80))
                    break;
            }
            double_double const st = xd::mul({s, 0.0}, total);              // dans (−1, ∞)
            double_double l1{};
            if (st.hi >= 0)
                l1 = xd::log1p_positive(st.hi, st.lo);
            else
            {
                auto const [h, lo] = xd::log1p_dd(st.hi);
                l1 = {h, lo + st.lo / (1.0 + st.hi)};
            }
            double_double result = xd::mul({s, 0.0}, xd::log_precise(u));
            result = xd::sum(result, {-lsb.hi, -lsb.lo});
            return xd::sum(result, l1);
        }

        // log F(u) et pente u f(u)/F(u), F = Iᵤ(a, b) (ou 1 − Iᵤ si upper).
        constexpr side_log beta_side(double a, double b, double u, bool upper) noexcept
        {
            auto const [y_hi, y_lo] = two_sum(1.0, -u);
            xd::scaled const density = beta_prefix(a, b, u, {y_hi, y_lo}) / xd::normalized(y_hi, y_lo);   // u f(u)
            if (!upper && a <= 0.25 && u <= 0.75 && fmul(b, u) <= 8)
            {
                double_double const l = log_ibeta_small(a, b, u);
                return {l, xd::rounded(density / xd::exp_scaled(l.hi, l.lo))};
            }
            auto const [v, complement] = incomplete_beta(a, b, u);
            xd::scaled const value = complement == upper ? v : one_minus_scaled(v);
            if (value.hi == 0)
                return {{-infinity, 0.0}, infinity};
            return {xd::log(value), density.hi == 0 ? 0.0 : xd::rounded(density / value)};
        }

        // u dans (0, 1/2] tel que F(u) = t, F = Iᵤ(a, b) ou 1 − Iᵤ(a, b) (upper) ; en
        // double-double.
        constexpr double_double beta_inverse_u(double a, double b, double t, bool upper) noexcept
        {
            auto const [lt_hi, lt_lo] = xd::log_precise(t);
            // Départ : queue basse Iᵤ ≈ uᵃ/(a B) si elle donne u < 1/2, sinon 1/4.
            double u = 0.25;
            if (!upper)
            {
                double const l = (log(t) + log(a) + lbeta_double(a, b)) / a;
                if (l < log(0.5))
                    u = l < -744 ? 0x1p-1074 : exp(l);
            }
            double low = 0, high = 0.5;
            for (int iteration = 0; iteration < 300; ++iteration)
            {
                auto const [l, slope] = beta_side(a, b, u, upper);
                double const g = (l.hi - lt_hi) + (l.lo - lt_lo);
                double const s = upper ? -slope : slope;
                if (g == 0)
                    return {u, 0.0};
                bool const too_far = (g > 0) == (s > 0);
                if (too_far)
                    high = u;
                else
                    low = u;
                double const delta = s != 0 && xd::is_finite(s) ? g / s : (too_far ? 4.0 : -4.0);
                double next = xd::magnitude(delta) < 0x1p-20 ? u - fmul(u, delta) : fmul(u, exp(-delta));
                if (xd::magnitude(delta) < 0x1p-56 || next == u)
                {
                    auto const [m, m_error] = two_product(u, delta);
                    auto const [r, r_error] = two_sum(u, -m);
                    auto const [hi, lo] = fast_two_sum(r, r_error - m_error);
                    return {hi, lo};
                }
                if (!(next > low && (next < high || (next == 0.5 && high == 0.5))))   // 1/2 elle-même est permise
                    next = low > 0 ? exp(fmul(log(low) + log(high), 0.5)) : fmul(high, 0x1.2c1599b2a5cb3p-6);   // e^-4
                if (!(next > 0))
                {
                    if (high <= 0x1p-1074)
                        return {0.0, 0.0};
                    next = 0x1p-1074;
                }
                u = next;
            }
            return {u, 0.0};
        }

        // x et y = 1 − x, en double-double tous deux, tels que Iₓ(a, b) = p (complement :
        // 1 − Iₓ(a, b) = p), p dans (0, 1) : le plus petit est exact à l'arrondi près,
        // l'autre vaut 1 − lui.
        struct beta_solution
        {
            double_double x, y;
        };

        constexpr beta_solution beta_inverse_xy(double a, double b, double p, bool complement) noexcept
        {
            auto const [half_lower, half_upper] = ibeta_halves(a, b);
            bool const in_x = complement ? p >= half_upper : p <= half_lower;
            bool const minority_upper = complement ? p <= 0.5 : p > 0.5;
            double const t = p <= 0.5 ? p : 1.0 - p;
            double_double const u = in_x ? beta_inverse_u(a, b, t, minority_upper) : beta_inverse_u(b, a, t, !minority_upper);
            auto const [r, r_error] = two_sum(1.0, -u.hi);
            auto const [w_hi, w_lo] = fast_two_sum(r, r_error - u.lo);
            double_double const w{w_hi, w_lo};
            return in_x ? beta_solution{u, w} : beta_solution{w, u};
        }

        // x tel que Iₓ(a, b) = p (complement : 1 − Iₓ(a, b) = p), p dans (0, 1).
        constexpr double beta_inverse(double a, double b, double p, bool complement) noexcept
        {
            auto const [half_lower, half_upper] = ibeta_halves(a, b);
            bool const in_x = complement ? p >= half_upper : p <= half_lower;   // solution x <= 1/2
            bool const minority_upper = complement ? p <= 0.5 : p > 0.5;    // côté 1 − I minoritaire
            double const t = p <= 0.5 ? p : 1.0 - p;                        // exact
            if (in_x)
            {
                auto const [hi, lo] = beta_inverse_u(a, b, t, minority_upper);
                return hi + lo;
            }
            // u = 1 − x, paramètres échangés : Iₓ(a, b) = 1 − Iᵤ(b, a), les côtés s'échangent.
            auto const [hi, lo] = beta_inverse_u(b, a, t, !minority_upper);
            auto const [r, r_error] = two_sum(1.0, -hi);
            return r + (r_error - lo);
        }
    }

    constexpr double lbeta(double a, double b) noexcept
    {
        if (!(a == a) || !(b == b))
            return detail::opaque(a + b);
        if (!(a > 0) || !(b > 0))
            return detail::quiet_nan;
        if (a == detail::infinity || b == detail::infinity)
            return -detail::infinity;
        if (a == 1 || b == 1)                                               // B(1, b) = 1/b
            return -log(a == 1 ? b : a);
        auto const [hi, lo] = detail::extended::log_beta(a, b);
        return hi + lo;
    }

    constexpr double beta(double a, double b) noexcept
    {
        if (!(a == a) || !(b == b))
            return detail::opaque(a + b);
        if (!(a > 0) || !(b > 0))
            return detail::quiet_nan;
        if (a == detail::infinity || b == detail::infinity)
            return 0.0;
        if (a == 1 || b == 1)
            return 1.0 / (a == 1 ? b : a);
        auto const [hi, lo] = detail::extended::log_beta(a, b);
        if (hi > 1400)
            return detail::infinity;
        if (hi < -1500)
            return 0.0;
        return detail::extended::rounded(detail::extended::exp_scaled(hi, lo));
    }

    constexpr double ibeta(double a, double b, double x) noexcept
    {
        double special = 0;
        if (detail::incomplete_beta_special(a, b, x, false, special))
            return special;
        auto const [value, complement] = detail::incomplete_beta(a, b, x);
        return complement ? detail::one_minus(value) : detail::extended::rounded(value);
    }

    constexpr double ibetac(double a, double b, double x) noexcept
    {
        double special = 0;
        if (detail::incomplete_beta_special(a, b, x, true, special))
            return special;
        auto const [value, complement] = detail::incomplete_beta(a, b, x);
        return complement ? detail::extended::rounded(value) : detail::one_minus(value);
    }

    constexpr float lbeta(float a, float b) noexcept { return detail::to_float(lbeta(static_cast<double>(a), static_cast<double>(b))); }
    constexpr float beta(float a, float b) noexcept { return detail::to_float(beta(static_cast<double>(a), static_cast<double>(b))); }
    constexpr float ibeta(float a, float b, float x) noexcept
    {
        return detail::to_float(ibeta(static_cast<double>(a), static_cast<double>(b), static_cast<double>(x)));
    }
    constexpr float ibetac(float a, float b, float x) noexcept
    {
        return detail::to_float(ibetac(static_cast<double>(a), static_cast<double>(b), static_cast<double>(x)));
    }

    // Inverses : ibeta_inv(a, b, p) = x tel que Iₓ(a, b) = p ; ibetac_inv(a, b, q) = x tel que
    // 1 − Iₓ(a, b) = q. p (q) dans [0, 1] : bords 0 et 1 ; sinon NaN.
    constexpr double ibeta_inv(double a, double b, double p) noexcept
    {
        if (!(a == a) || !(b == b) || !(p == p))
            return detail::opaque(a + b + p);
        if (!(a > 0) || !(b > 0) || a == detail::infinity || b == detail::infinity || !(p >= 0 && p <= 1))
            return detail::quiet_nan;
        if (p == 0 || p == 1)
            return p;
        return detail::beta_inverse(a, b, p, false);
    }

    constexpr double ibetac_inv(double a, double b, double q) noexcept
    {
        if (!(a == a) || !(b == b) || !(q == q))
            return detail::opaque(a + b + q);
        if (!(a > 0) || !(b > 0) || a == detail::infinity || b == detail::infinity || !(q >= 0 && q <= 1))
            return detail::quiet_nan;
        if (q == 0 || q == 1)
            return 1.0 - q;
        return detail::beta_inverse(a, b, q, true);
    }

    constexpr float ibeta_inv(float a, float b, float p) noexcept
    {
        return detail::to_float(ibeta_inv(static_cast<double>(a), static_cast<double>(b), static_cast<double>(p)));
    }
    constexpr float ibetac_inv(float a, float b, float q) noexcept
    {
        return detail::to_float(ibetac_inv(static_cast<double>(a), static_cast<double>(b), static_cast<double>(q)));
    }

    // Arguments mêlés (entiers, float, double), en double, comme <cmath>.
    template <class A, class B> requires detail::mixed_arithmetic<A, B>
    constexpr double beta(A a, B b) noexcept { return beta(static_cast<double>(a), static_cast<double>(b)); }
    template <class A, class B> requires detail::mixed_arithmetic<A, B>
    constexpr double lbeta(A a, B b) noexcept { return lbeta(static_cast<double>(a), static_cast<double>(b)); }
    template <class A, class B, class C> requires detail::mixed_arithmetic<A, B, C>
    constexpr double ibeta(A a, B b, C x) noexcept
    {
        return ibeta(static_cast<double>(a), static_cast<double>(b), static_cast<double>(x));
    }
    template <class A, class B, class C> requires detail::mixed_arithmetic<A, B, C>
    constexpr double ibetac(A a, B b, C x) noexcept
    {
        return ibetac(static_cast<double>(a), static_cast<double>(b), static_cast<double>(x));
    }
    template <class A, class B, class C> requires detail::mixed_arithmetic<A, B, C>
    constexpr double ibeta_inv(A a, B b, C p) noexcept
    {
        return ibeta_inv(static_cast<double>(a), static_cast<double>(b), static_cast<double>(p));
    }
    template <class A, class B, class C> requires detail::mixed_arithmetic<A, B, C>
    constexpr double ibetac_inv(A a, B b, C q) noexcept
    {
        return ibetac_inv(static_cast<double>(a), static_cast<double>(b), static_cast<double>(q));
    }

    long double beta(long double, long double) = delete;
    long double lbeta(long double, long double) = delete;
    long double ibeta(long double, long double, long double) = delete;
    long double ibetac(long double, long double, long double) = delete;
    long double ibeta_inv(long double, long double, long double) = delete;
    long double ibetac_inv(long double, long double, long double) = delete;
}
