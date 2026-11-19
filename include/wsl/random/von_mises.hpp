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
// entropie log(2π I₀(κ)) − κA. Répartition : la part la plus petite,
// ∫_{|θ|}^π e^{κ(cos t − 1)} dt/(2π e^{−κ} I₀(κ)), par morceaux de Taylor en double-double
// (detail::von_mises_tail), l'autre par complément ; quantiles par Newton sur log G.
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
#include <wsl/random/detail/discrete.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/normal.hpp>
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

        constexpr sd::double_double neg_dd(sd::double_double v) noexcept { return {-v.hi, -v.lo}; }

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

        // a + b, a, b >= 0 à l'échelle (zéro : hi = 0).
        constexpr scaled add_scaled(scaled const& a, scaled const& b) noexcept
        {
            if (a.hi == 0)
                return b;
            if (b.hi == 0)
                return a;
            scaled const& big = a.e >= b.e ? a : b;
            scaled const& small = a.e >= b.e ? b : a;
            if (big.e - small.e > 120)
                return big;
            int const shift = small.e - big.e;
            auto const t = sum({big.hi, big.lo}, {sd::scale_by_power_of_two(small.hi, shift), sd::scale_by_power_of_two(small.lo, shift)});
            return normalized(t.hi, t.lo, big.e);
        }

        // Gauss–Legendre à 12 nœuds sur [−1, 1] (6 nœuds positifs et leurs poids), en
        // double-double : Newton sur P₁₂ (récurrence de Bonnet), à la compilation.
        inline constexpr int legendre_order = 12;
        struct legendre_rule
        {
            sd::double_double node[legendre_order / 2], weight[legendre_order / 2];
        };
        constexpr legendre_rule make_legendre_rule() noexcept
        {
            legendre_rule r{};
            constexpr int m = legendre_order;
            for (int i = 0; i < m / 2; ++i)
            {
                double const guess = wsl::special::cos(pi_hi * (static_cast<double>(i) + 0.75) / (static_cast<double>(m) + 0.5));
                sd::double_double x{guess, 0.0}, derivative{1.0, 0.0};
                for (int iteration = 0; iteration < 8; ++iteration)
                {
                    sd::double_double p0{1.0, 0.0}, p1 = x;                 // P_{k−1}, P_k
                    for (int k = 1; k < m; ++k)
                    {
                        auto const next = sd::quotient_dd(sum(mul({static_cast<double>(2 * k + 1), 0.0}, mul(x, p1)), neg_dd(mul({static_cast<double>(k), 0.0}, p0))),
                                                          {static_cast<double>(k + 1), 0.0});
                        p0 = p1;
                        p1 = next;
                    }
                    auto const one_minus_x2 = sum({1.0, 0.0}, neg_dd(mul(x, x)));
                    derivative = sd::quotient_dd(mul({static_cast<double>(m), 0.0}, sum(p0, neg_dd(mul(x, p1)))), one_minus_x2);   // P′ = m (P_{m−1} − x P_m)/(1 − x²)
                    x = sum(x, neg_dd(sd::quotient_dd(p1, derivative)));
                }
                r.node[i] = x;
                auto const one_minus_x2 = sum({1.0, 0.0}, neg_dd(mul(x, x)));
                r.weight[i] = sd::quotient_dd({2.0, 0.0}, mul(one_minus_x2, mul(derivative, derivative)));
            }
            return r;
        }
        inline constexpr legendre_rule legendre = make_legendre_rule();

        // κ (cos t − 1) = −2κ sin²(t/2), 0 <= t <= π en double-double.
        constexpr sd::double_double von_mises_exponent(double kappa, sd::double_double t) noexcept
        {
            auto const s = sine_dd(sd::fmul(t.hi, 0.5), sd::fmul(t.lo, 0.5));
            return mul({sd::fmul(kappa, -2.0), 0.0}, mul(s, s));
        }

        // ∫_a^b e^{κ (cos t − 1)} dt, 0 <= a <= b <= π en double-double, à l'échelle : par morceaux
        // [t₀, t₀ + H], Gauss–Legendre à 12 nœuds sur e^{κ (cos t − cos t₀)} (erreur relative
        // sous 2^-100 pour κ H sin t₀ <= 2 et κ H² |cos t₀| <= ½ ; nœuds à 2^-70, ceux de
        // exp_scaled), le facteur e^{κ (cos t₀ − 1)} à l'échelle. L'intégrande décroît sur
        // [0, π] : arrêt quand e^{κ (cos t₀ − 1)} (b − t₀) devient négligeable.
        constexpr scaled von_mises_integral(double kappa, sd::double_double a, sd::double_double b) noexcept
        {
            namespace fp = sd;
            scaled total{0.0, 0.0, 0};
            sd::double_double t0 = a;
            for (int piece = 0; piece < 1000000; ++piece)
            {
                auto const rem = sum(b, {-t0.hi, -t0.lo});
                if (!(rem.hi > 0))
                    break;
                auto const s = sine_dd(fp::fmul(t0.hi, 0.5), fp::fmul(t0.lo, 0.5));                  // sin(t₀/2)
                auto const half_rem = sum({fp::fmul(pi_hi, 0.5), fp::fmul(pi_lo, 0.5)}, {-fp::fmul(t0.hi, 0.5), -fp::fmul(t0.lo, 0.5)});
                auto const c = sine_dd(half_rem.hi, half_rem.lo);                                    // cos(t₀/2)
                auto const s2 = mul(s, s);
                double const sin_t = fp::fmul(2.0, fp::fmul(s.hi, c.hi));
                double const cos_t = fp::fmul(c.hi, c.hi) - fp::fmul(s.hi, s.hi);
                auto const h0 = mul({fp::fmul(kappa, -2.0), 0.0}, s2);                              // κ (cos t₀ − 1)
                if (!(h0.hi > -1500))
                    break;
                scaled const level = exp_scaled(h0.hi, h0.lo);
                if (total.hi != 0)                                          // reste <= e^{h₀} (b − t₀)
                {
                    scaled const bound = level * normalized(rem.hi, rem.lo);
                    if (bound.e < total.e - 80)
                        break;
                }
                double width = rem.hi > 0.5 ? 0.5 : rem.hi;
                double const ks = fp::fmul(kappa, magnitude(sin_t)), kc = fp::fmul(kappa, magnitude(cos_t));
                if (fp::fmul(ks, width) > 2.0)
                    width = 2.0 / ks;
                if (fp::fmul(kc, fp::fmul(width, width)) > 0.5)
                    width = wsl::special::sqrt(0.5 / kc);
                bool const last = width >= rem.hi;
                sd::double_double const h = last ? rem : sd::double_double{width, 0.0};
                sd::double_double const half{fp::fmul(h.hi, 0.5), fp::fmul(h.lo, 0.5)};
                auto const middle = sum(t0, half);
                sd::double_double integral{0.0, 0.0};
                for (int i = 0; i < legendre_order / 2; ++i)
                {
                    auto const offset = mul(half, legendre.node[i]);
                    sd::double_double value{0.0, 0.0};
                    for (int side = 0; side < 2; ++side)
                    {
                        auto const t = sum(middle, side == 0 ? offset : neg_dd(offset));
                        auto const e = sum(von_mises_exponent(kappa, t), neg_dd(h0));     // κ (cos t − cos t₀)
                        if (e.hi > -1500)
                            value = sum(value, unscaled(exp_scaled(e.hi, e.lo)));
                    }
                    integral = sum(integral, mul(legendre.weight[i], value));
                }
                integral = mul(integral, half);
                if (integral.hi > 0)
                    total = add_scaled(total, level * normalized(integral.hi, integral.lo));
                if (last)
                    break;
                t0 = sum(t0, h);
            }
            return total;
        }

        constexpr scaled von_mises_tail(double kappa, sd::double_double phi) noexcept
        {
            return von_mises_integral(kappa, phi, {pi_hi, pi_lo});
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

        // P(X <= x) sur [μ − π, μ + π] : la part la plus petite, ∫_{|θ|}^π e^{κ(cos t − 1)} dt/(2π e^{−κ} I₀(κ))
        // (θ = x − μ en double-double), l'autre par complément ; un arrondi.
        friend constexpr T cdf(von_mises const& V, T x) noexcept { return V.distribution(x, false); }
        friend constexpr T ccdf(von_mises const& V, T x) noexcept { return V.distribution(x, true); }
        friend constexpr T quantile(von_mises const& V, double p) { return V.inverse(p, false); }
        friend constexpr T cquantile(von_mises const& V, double q) { return V.inverse(q, true); }

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

        constexpr T distribution(T x, bool upper) const noexcept
        {
            namespace fp = detail::sd;
            double const v = x;
            if (!(v == v))
                return x;
            auto const [d, d_error] = fp::two_sum(v, -double(mu_));
            if (!detail::is_finite(d))
                return (d > 0) != upper ? T(1) : T(0);
            auto const lim = detail::sum({detail::pi_hi, detail::pi_lo}, {-detail::magnitude(d), d < 0 ? d_error : -d_error});
            if (!(lim.hi > 0))                                              // hors de (μ − π, μ + π)
                return (d > 0) != upper ? T(1) : T(0);
            bool const negative = d < 0;
            fp::double_double const phi = negative ? fp::double_double{-d, -d_error} : fp::double_double{d, d_error};
            detail::scaled g = detail::von_mises_tail(double(kappa_), phi); // F si θ < 0, S sinon
            if (g.hi != 0)
                g = g * detail::exp_scaled(-log_norm_.hi, -log_norm_.lo);
            bool const direct = negative != upper;
            if (direct)
                return static_cast<T>(detail::rounded(g));
            return static_cast<T>(fp::one_minus(g));
        }

        // φ tel que G(φ) = t, t <= ½ : pour t <= ¼, Newton sur log G (G(φ) = ∫_φ^π g,
        // décroissante) ; au-delà, sur log D, D(φ) = ∫_0^φ g = ½ − t (exact), croissante, qui
        // fixe φ petit sans compensation. Encadré dans [0, π], puis deux pas à résidu précis.
        constexpr detail::sd::double_double solve(double t) const
        {
            namespace fp = detail::sd;
            bool const head = t > 0.25;
            double const target = head ? 0.5 - t : t;                       // exact (Sterbenz)
            if (target == 0)
                return {0.0, 0.0};
            auto const log_t = detail::log_precise(target);
            detail::scaled const inv_norm = detail::exp_scaled(-log_norm_.hi, -log_norm_.lo);
            auto const value_at = [&](fp::double_double phi) {
                detail::scaled v = head ? detail::von_mises_integral(double(kappa_), {0.0, 0.0}, phi) : detail::von_mises_tail(double(kappa_), phi);
                return v.hi == 0 ? v : v * inv_norm;
            };
            struct step
            {
                int side;                                                   // > 0 : φ trop petit
                fp::double_double dphi;
                double residual;                                            // |log(V/t)|
                bool finite;
            };
            auto const evaluate = [&](fp::double_double phi) {
                detail::scaled const value = value_at(phi);
                int const c = value.hi == 0 ? -1 : detail::compare(value, {target, 0.0});
                step s{head ? -c : c, {0.0, 0.0}, 0.0, false};
                auto const sh = detail::sine_dd(fp::fmul(phi.hi, 0.5), fp::fmul(phi.lo, 0.5));
                auto const e = detail::mul({fp::fmul(double(kappa_), -2.0), 0.0}, detail::mul(sh, sh));
                if (value.hi == 0 || !(e.hi > -1500))
                    return s;
                detail::scaled const density = detail::exp_scaled(e.hi, e.lo) * inv_norm;
                auto const h = detail::sum(detail::log(value), {-log_t.hi, -log_t.lo});   // log(V/t)
                auto const r = detail::unscaled(value / density);                         // V/g
                s.dphi = detail::mul(h, r);                                                 // G : φ + h V/g ; D : φ − h V/g
                if (head)
                    s.dphi = detail::neg_dd(s.dphi);
                s.residual = detail::magnitude(h.hi);
                s.finite = detail::is_finite(s.dphi.hi);
                return s;
            };
            if (!head)                                                      // solution au-delà de fl(π) : G(φ) = (π − φ) g(π)
            {
                detail::scaled const at_pi = value_at({detail::pi_hi, 0.0});
                if (at_pi.hi != 0 && detail::compare(at_pi, {target, 0.0}) > 0)
                {
                    detail::scaled const g_pi = detail::exp_scaled(fp::fmul(double(kappa_), -2.0), 0.0) * inv_norm;
                    auto const d = detail::unscaled(detail::normalized(target, 0.0) / g_pi);
                    return detail::sum({detail::pi_hi, detail::pi_lo}, detail::neg_dd(d));
                }
            }
            double lo = 0, hi = detail::pi_hi;
            double phi = 0;
            if (kappa_ > 0)                                                 // départ : loi normale d'écart type 1/√κ
            {
                auto const z = detail::standard_quantile(t);               // < 0
                phi = -fp::fmul(z.hi, 1.0 / wsl::special::sqrt(double(kappa_)));
                if (!(phi < 3.0))
                    phi = 3.0;
            }
            else
                phi = fp::fmul(detail::pi_hi, 1.0 - fp::fmul(2.0, t));
            if (!(phi > 0))
                phi = 0x1p-30;
            for (int iteration = 0; iteration < 200; ++iteration)
            {
                step const s = evaluate({phi, 0.0});
                if (s.side == 0)
                    return {phi, 0.0};
                if (s.side > 0)
                    lo = phi > lo ? phi : lo;
                else
                    hi = phi < hi ? phi : hi;
                double next = phi + s.dphi.hi;
                // résidu relatif sous 2^-36, ou pas plus fin que les doubles (près de π) : deux pas en double-double
                if (s.finite && (0x1p-36 > s.residual || next == phi))
                {
                    fp::double_double r = detail::sum({phi, 0.0}, s.dphi);
                    step const s2 = evaluate(r);
                    if (s2.finite)
                        r = detail::sum(r, s2.dphi);
                    return r;
                }
                if (!s.finite || !(next > lo && next < hi))
                    next = fp::fmul(lo, 0.5) + fp::fmul(hi, 0.5);
                if (next == phi)
                    return {phi, 0.0};
                phi = next;
            }
            return {phi, 0.0};
        }

        constexpr T inverse(double t, bool upper) const
        {
            namespace fp = detail::sd;
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::quantile : probabilité hors de [0, 1]");
            bool from_low = !upper;                                         // x = μ − φ (côté bas) ou μ + φ
            double small = t;
            if (t > 0.5)
            {
                small = 1.0 - t;                                            // exact
                from_low = !from_low;
            }
            fp::double_double phi{detail::pi_hi, detail::pi_lo};
            if (small > 0)
                phi = solve(small);
            fp::double_double const x = detail::sum({double(mu_), 0.0}, from_low ? detail::neg_dd(phi) : phi);
            return static_cast<T>(x.hi + x.lo);
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
