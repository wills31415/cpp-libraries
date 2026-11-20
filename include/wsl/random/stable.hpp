#pragma once

// Lois stables S(α, β, γ, δ), paramétrisation S1 (celle de Chambers–Mallows–Stuck) :
// E e^{itX} = exp(−γ^α |t|^α (1 − iβ sign(t) tan(πα/2)) + iδt) (α ≠ 1),
//             exp(−γ |t| (1 + iβ (2/π) sign(t) log|t|) + iδt)        (α = 1),
// 0 < α <= 2, −1 <= β <= 1.
//
//   wsl::random::stable X({.location = 0, .scale = 1, .alpha = 1.5, .beta = 0.5});
//   X.pdf(0.3); X.cdf(-2); X.quantile(0.99); X(g);
//
// Forme réduite Z : X = γZ + δ (α ≠ 1), γZ + (2/π)βγ log γ + δ (α = 1). Densité et
// répartition de Z par les intégrales de Nolan (1997), dans la variable φ = θ + θ₀ (sans
// compensation près de −θ₀), θ₀ = arctan(β tan(πα/2))/α, en double-double :
//     g(φ) = y^{α/(α−1)} V(θ),  log V = log cos(αθ₀)/(α − 1) + (α/(α−1)) (log cos θ − log sin αφ)
//                                       + log cos(θ₀ + (α − 1)φ) − log cos θ,
// pour y > 0 (y < 0 : f(y; β) = f(−y; −β)). Avec A = (1/π) ∫ e^{−g}, B = (1/π) ∫ (1 − e^{−g})
// (expm1) et c = (π/2 − θ₀)/π : α < 1 : F = c + A, S = B ; α > 1 : F = c + B, S = A — toujours
// des sommes de termes positifs ; densité α/(π |α − 1| y) ∫ g e^{−g}. α = 1 (β > 0) :
// g = e^{−πy/(2β)} V, log V = log(2/π) + log(π/2 + βθ) − log cos θ + (π/2 + βθ) tan θ/β,
// F = A, S = B, densité (1/(2β)) ∫ g e^{−g}. Intégration : coupure au pic (g = 1, g
// monotone), puis Gauss–Legendre adaptatif (12 nœuds, comparaison à deux demi-intervalles),
// intégrandes à l'échelle. Cas particuliers en forme close : α = 2 (normale de variance
// 2γ²), α = 1 et β = 0 (Cauchy), y = 0.
// Quantiles : Newton sur le logarithme de la part visée, encadré.
//
// Tirage (approché, bit-exact) : Chambers, Mallows et Stuck (1976), forme de Weron (1996).

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/trigonometry.hpp>
#include <wsl/random/detail/quadrature.hpp>
#include <wsl/random/detail/root.hpp>
#include <wsl/random/detail/gaussian.hpp>
#include <wsl/random/detail/discrete.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/exponential.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    namespace detail
    {

        // Noyau de Nolan pour la forme réduite, β fixé (déjà réfléchi si y < 0).
        struct stable_kernel
        {
            double alpha = 1, beta = 0;
            sd::double_double theta0{}, top{};                              // θ₀, Φ = π/2 + θ₀ (α ≠ 1) ou π (α = 1)
            sd::double_double c{};                                          // (π/2 − θ₀)/π
            sd::double_double log_cos_a0{};                                 // log cos(αθ₀)
            sd::double_double a1{}, a2{}, am1{};                            // α/(α − 1), 1/(α − 1), α − 1
            bool one = false;

            constexpr stable_kernel() = default;
            constexpr stable_kernel(double a, double b) noexcept : alpha(a), beta(b), one(a == 1)
            {
                if (one)
                {
                    top = {pi_hi, pi_lo};
                    return;
                }
                auto const u = mul({sd::fmul(pi_hi, 0.5), sd::fmul(pi_lo, 0.5)}, {a, 0.0});   // πα/2 dans (0, π)
                auto const s = sin_any(u);
                auto const cu = sum(half_pi_dd, negated(u));
                auto const co = sine_dd(cu.hi, cu.lo);                     // cos(πα/2)
                auto const t = mul({b, 0.0}, sd::quotient_dd(s, co));      // β tan(πα/2)
                auto const at = atan_dd(t);                                 // αθ₀
                theta0 = sd::quotient_dd(at, {a, 0.0});
                top = sum(half_pi_dd, theta0);
                c = mul(sum(half_pi_dd, negated(theta0)), {inv_pi_hi, inv_pi_lo});
                // log cos(αθ₀) = −½ log(1 + t²)
                double const mt = magnitude(t.hi);
                sd::double_double l{};
                if (mt == 0)
                    l = {0.0, 0.0};
                else if (mt < 0x1p500)
                    l = log_positive(sum({1.0, 0.0}, mul(t, t)));
                else
                {
                    auto const lt = log_positive(t.hi > 0 ? t : negated(t));
                    l = {sd::fmul(lt.hi, 2.0), sd::fmul(lt.lo, 2.0)};
                }
                log_cos_a0 = {sd::fmul(l.hi, -0.5), sd::fmul(l.lo, -0.5)};
                am1 = sum({a, 0.0}, {-1.0, 0.0});
                a1 = sd::quotient_dd({a, 0.0}, am1);
                a2 = sd::quotient_dd({1.0, 0.0}, am1);
            }

            // log g(φ), shift = (α/(α−1)) log y (α ≠ 1) ou −πy/(2β) (α = 1).
            constexpr sd::double_double log_g(sd::double_double phi, sd::double_double shift) const noexcept
            {
                if (one)
                {
                    auto const sp = sin_any(phi);                           // cos θ = sin φ
                    auto const hc = sum(half_pi_dd, negated(phi));
                    auto const cp = sine_dd(hc.hi, hc.lo);                  // cos φ ; tan θ = −cos φ/sin φ
                    auto const lin = sum(mul(sum({1.0, 0.0}, {-beta, 0.0}), half_pi_dd), mul({beta, 0.0}, phi));   // π/2 + βθ = (1 − β)π/2 + βφ
                    auto l = sum(sum({sd::ln2_hi, sd::ln2_lo}, {-sd::log_pi_hi, -sd::log_pi_lo}), log_positive(lin));   // log(2/π) + log(π/2 + βθ)
                    l = sum(l, negated(log_positive(sp)));
                    auto const tan_t = negated(sd::quotient_dd(cp, sp));
                    l = sum(l, sd::quotient_dd(mul(lin, tan_t), {beta, 0.0}));
                    return sum(l, shift);
                }
                auto const sa = sin_any(mul({alpha, 0.0}, phi));           // sin αφ
                auto const ct = sin_any(sum(top, negated(phi)));            // cos θ = sin(Φ − φ)
                auto const arg = sum(sum(half_pi_dd, negated(theta0)), negated(mul(am1, phi)));   // π/2 − θ₀ − (α − 1)φ
                auto const cm = sin_any(arg);                               // cos(θ₀ + (α − 1)φ)
                auto const lct = log_positive(ct);
                auto l = mul(a2, log_cos_a0);
                l = sum(l, mul(a1, sum(lct, negated(log_positive(sa)))));
                l = sum(l, log_positive(cm));
                l = sum(l, negated(lct));
                return sum(l, shift);
            }
        };

        enum class stable_integrand
        {
            survival_like,                                                  // e^{−g}
            complement_like,                                                // 1 − e^{−g}
            density                                                         // g e^{−g}
        };

        constexpr scaled stable_value(stable_integrand kind, sd::double_double lg) noexcept
        {
            if (!(lg.hi == lg.hi))
                return {0.0, 0.0, 0};
            if (kind == stable_integrand::complement_like)
            {
                if (lg.hi < -700)
                    return exp_scaled(lg.hi, lg.lo);                        // 1 − e^{−g} = g (1 − g/2 …)
                if (lg.hi > 4)                                              // g > 54 : e^{−g} < 2^-77
                    return {1.0, 0.0, 0};
                auto const g = unscaled(exp_scaled(lg.hi, lg.lo));
                if (g.hi <= 1)
                {
                    auto const e = expm1_double_double(-g.hi, -g.lo);       // < 0
                    return normalized(-e.hi, -e.lo);
                }
                auto const e = unscaled(exp_scaled(-g.hi, -g.lo));
                auto const r = sum({1.0, 0.0}, negated(e));
                return normalized(r.hi, r.lo);
            }
            if (lg.hi > 7.3)                                                // g > 1480
                return {0.0, 0.0, 0};
            if (lg.hi < -1400)
                return kind == stable_integrand::survival_like ? scaled{1.0, 0.0, 0} : scaled{0.0, 0.0, 0};
            scaled const gs = exp_scaled(lg.hi, lg.lo);
            sd::double_double g{0.0, 0.0};
            if (gs.e > -1000)
                g = unscaled(gs);
            scaled const e = exp_scaled(-g.hi, -g.lo);
            return kind == stable_integrand::survival_like ? e : gs * e;
        }

    }

    template <law_real T = double>
    class stable : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T location = 0, scale = 1, alpha, beta = 0;
        };

        constexpr stable(parameters p) : alpha_(p.alpha), beta_(p.beta), gamma_(p.scale), delta_(p.location)
        {
            if (!(alpha_ > 0 && alpha_ <= 2) || !(beta_ >= -1 && beta_ <= 1) || !(gamma_ > 0) || !detail::is_finite(gamma_) || !detail::is_finite(delta_))
                throw std::domain_error("wsl::random::stable : il faut 0 < alpha <= 2, −1 <= beta <= 1, scale > 0 fini, location fini");
            double const a = alpha_, b = beta_;
            if (a < 2 && !(a == 1 && b == 0))
            {
                positive_ = detail::stable_kernel(a, a == 1 ? detail::magnitude(b) : b);
                negative_ = detail::stable_kernel(a, a == 1 ? detail::magnitude(b) : -b);
            }
            if (a == 1 && b != 0)                                           // (2/π) β log γ
            {
                auto const lg = detail::log_precise(double(gamma_));
                log_shift_ = detail::mul(detail::mul({b, 0.0}, lg), {detail::sd::fmul(detail::inv_pi_hi, 2.0), detail::sd::fmul(detail::inv_pi_lo, 2.0)});
            }
        }

        constexpr T alpha() const noexcept { return alpha_; }
        constexpr T beta() const noexcept { return beta_; }
        constexpr T scale() const noexcept { return gamma_; }
        constexpr T location() const noexcept { return delta_; }
        constexpr interval<T> support() const noexcept
        {
            constexpr T inf = std::numeric_limits<T>::infinity();
            if (alpha_ < 1 && beta_ == 1)
                return {delta_, inf};
            if (alpha_ < 1 && beta_ == -1)
                return {-inf, delta_};
            return {-inf, inf};
        }
        // δ pour α > 1, indéfinie sinon
        constexpr T mean() const noexcept { return alpha_ > 1 ? delta_ : std::numeric_limits<T>::quiet_NaN(); }
        // 2γ² pour α = 2, infinie sinon
        constexpr T variance() const noexcept
        {
            if (alpha_ < 2)
                return std::numeric_limits<T>::infinity();
            return static_cast<T>(detail::rounded(detail::normalized(double(gamma_), 0.0) * detail::normalized(double(gamma_), 0.0) * detail::normalized(2.0, 0.0)));
        }
        constexpr T median() const { return wsl::random::quantile(*this, 0.5); }

        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (v == std::numeric_limits<double>::infinity() || v == -std::numeric_limits<double>::infinity())
                return 0;
            auto const y = X.reduced(v);
            detail::scaled const f = X.reduced_density(y);
            if (f.hi == 0)
                return 0;
            return static_cast<T>(detail::rounded(f / double(X.gamma_)));
        }
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v == v))
                return x;
            auto const y = X.reduced(v);
            detail::scaled const f = X.reduced_density(y);
            if (f.hi == 0)
                return -std::numeric_limits<T>::infinity();
            auto const l = detail::sum(detail::log(f), detail::negated(detail::log_precise(double(X.gamma_))));
            return static_cast<T>(l.hi + l.lo);
        }
        constexpr T cdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, true); }
        constexpr T quantile(double p) const { auto const& X = *this; return X.inverse(p, false); }
        constexpr T cquantile(double q) const { auto const& X = *this; return X.inverse(q, true); }

        // Chambers–Mallows–Stuck : V uniforme sur (−π/2, π/2), W exponentielle.
        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            namespace sp = wsl::special;
            double const a = X.alpha_, b = X.beta_;
            double u = 0;
            do
                u = uniform_real<double>(g);
            while (u == 0);
            double const v = sd::fmul(detail::pi_hi, u - 0.5);
            double const w = wsl::random::sample(exponential<double>({.rate = 1.0}), g);
            double z = 0;
            if (a == 1)
            {
                double const h = sd::fmul(detail::pi_hi, 0.5);
                double const lin = h + sd::fmul(b, v);
                z = sd::fmul(2.0 / detail::pi_hi, sd::fmul(lin, sp::tan(v)) - sd::fmul(b, sp::log(sd::fmul(sd::fmul(h, w), sp::cos(v)) / lin)));
                return static_cast<T>((double(X.delta_) + sd::fmul(double(X.gamma_), z)) + sd::fmul(sd::fmul(2.0 / detail::pi_hi, b), sd::fmul(double(X.gamma_), sp::log(double(X.gamma_)))));
            }
            double const t = sd::fmul(b, sp::tan(sd::fmul(sd::fmul(detail::pi_hi, 0.5), a)));
            double const bb = sp::atan(t) / a;                              // θ₀
            double const s = sp::exp(sd::fmul(0.5 / a, sp::log1p(sd::fmul(t, t))));
            double const av = sd::fmul(a, v + bb);
            double const num = sp::sin(av);
            double const l = sd::fmul((1 - a) / a, sp::log(sp::cos(v - av)) - sp::log(w)) - sp::log(sp::cos(v)) / a;
            z = sd::fmul(sd::fmul(s, num), sp::exp(l));
            return static_cast<T>(double(X.delta_) + sd::fmul(double(X.gamma_), z));
        }

        friend constexpr bool operator==(stable const& a, stable const& b) noexcept
        {
            return a.alpha_ == b.alpha_ && a.beta_ == b.beta_ && a.gamma_ == b.gamma_ && a.delta_ == b.delta_;
        }

    private:
        // forme réduite y en double-double
        constexpr detail::sd::double_double reduced(double x) const noexcept
        {
            auto y = detail::standardize(x, double(delta_), double(gamma_));
            if (alpha_ == 1 && beta_ != 0)
                y = detail::sum(y, detail::negated(log_shift_));
            return y;
        }

        // F(y) (upper = false) ou S(y) de la forme réduite, à l'échelle : une seule intégrale.
        constexpr detail::scaled reduced_part(detail::sd::double_double y, bool upper) const noexcept
        {
            namespace sd = detail::sd;
            double const a = alpha_, b = beta_;
            if (!(y.hi == y.hi))
                return {0.0, 0.0, 0};
            detail::scaled const inv_pi = detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo);
            if (a == 2)                                                     // normale de variance 2 : Φ(±y/√2)
            {
                auto const w = detail::mul(y, {detail::inv_sqrt2_hi, detail::inv_sqrt2_lo});
                return detail::normal_lower_scaled(upper ? detail::negated(w) : w);
            }
            if (a == 1 && b == 0)                                           // Cauchy : (1/π) atan(1/|y|) pour la queue
            {
                bool const negative = y.hi < 0;
                if (y.hi == 0)
                    return {1.0, 0.0, -1};
                auto const ay = negative ? detail::negated(y) : y;
                auto const tail = detail::mul(detail::atan_dd(sd::quotient_dd({1.0, 0.0}, ay)), {detail::inv_pi_hi, detail::inv_pi_lo});
                detail::scaled const t = tail.hi > 0 ? detail::normalized(tail.hi, tail.lo) : detail::scaled{0.0, 0.0, 0};
                return negative != upper ? t : detail::scaled_complement(t);
            }
            if (a == 1)                                                     // β > 0 après réflexion : F = A, S = B
            {
                bool const reflect = b < 0;
                auto const yy = reflect ? detail::negated(y) : y;
                bool const want_upper = upper != reflect;
                auto const shift = sd::quotient_dd(detail::mul(yy, {detail::pi_hi, detail::pi_lo}), {-sd::fmul(2.0, positive_.beta), 0.0});   // −πy/(2β)
                detail::scaled const I = integral(positive_, want_upper ? detail::stable_integrand::complement_like : detail::stable_integrand::survival_like, shift);
                return I.hi == 0 ? I : I * inv_pi;
            }
            bool const negative = y.hi < 0;
            auto const& k = negative ? negative_ : positive_;
            auto const yy = negative ? detail::negated(y) : y;
            bool const want_upper = upper != negative;                      // part voulue pour yy >= 0
            detail::scaled const c = k.c.hi > 0 ? detail::normalized(k.c.hi, k.c.lo) : detail::scaled{0.0, 0.0, 0};
            if (yy.hi == 0)
                return want_upper ? detail::scaled_complement(c) : c;
            auto const shift = detail::mul(k.a1, detail::log_positive(yy));
            // α < 1 : F = c + A, S = B ; α > 1 : F = c + B, S = A
            bool const use_a = (a < 1) != want_upper;
            detail::scaled I = integral(k, use_a ? detail::stable_integrand::survival_like : detail::stable_integrand::complement_like, shift);
            if (I.hi != 0)
                I = I * inv_pi;
            return want_upper ? I : detail::add_scaled(c, I);
        }

        // densité de la forme réduite en y, à l'échelle
        constexpr detail::scaled reduced_density(detail::sd::double_double y) const noexcept
        {
            namespace sd = detail::sd;
            double const a = alpha_, b = beta_;
            if (a == 2)                                                     // e^{−y²/4}/(2√π)
            {
                auto const y2 = detail::mul(y, y);
                if (!(y2.hi < 6000))
                    return {0.0, 0.0, 0};
                detail::scaled f = detail::exp_scaled(sd::fmul(y2.hi, -0.25), sd::fmul(y2.lo, -0.25)) * detail::normalized(detail::inv_sqrt_pi, 0x1.1ae3a914fed80p-57);
                f.e -= 1;
                return f;
            }
            if (a == 1 && b == 0)                                           // 1/(π(1 + y²))
            {
                if (detail::magnitude(y.hi) > 0x1p500)
                {
                    detail::scaled const ay = detail::normalized(detail::magnitude(y.hi), y.hi < 0 ? -y.lo : y.lo);
                    return detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo) / (ay * ay);
                }
                auto const d = detail::sum({1.0, 0.0}, detail::mul(y, y));
                return detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo) / detail::normalized(d.hi, d.lo);
            }
            if (a == 1)
            {
                bool const reflect = b < 0;
                auto const yy = reflect ? detail::negated(y) : y;
                double const bk = positive_.beta;
                auto const shift = sd::quotient_dd(detail::mul(yy, {detail::pi_hi, detail::pi_lo}), {-sd::fmul(2.0, bk), 0.0});
                detail::scaled const I = integral(positive_, detail::stable_integrand::density, shift);
                return I.hi == 0 ? I : I / sd::fmul(2.0, bk);
            }
            bool const negative = y.hi < 0;
            auto const& k = negative ? negative_ : positive_;
            auto const yy = negative ? detail::negated(y) : y;
            if (yy.hi == 0)                                                 // Γ(1 + 1/α) cos θ₀ cos(αθ₀)^{1/α}/π
            {
                if (!(k.top.hi > 0) || !(k.c.hi > 0))                       // α < 1, |β| = 1 : bord du support, densité nulle
                    return {0.0, 0.0, 0};
                auto const arg = detail::sum({1.0, 0.0}, sd::quotient_dd({1.0, 0.0}, {a, 0.0}));   // 1 + 1/α en double-double
                auto const lg = detail::log_gamma(arg.hi, arg.lo);
                auto const ct = detail::sin_any(detail::sum(detail::half_pi_dd, detail::negated(k.theta0)));
                auto l = detail::sum(lg, detail::log_positive(ct));
                l = detail::sum(l, sd::quotient_dd(k.log_cos_a0, {a, 0.0}));
                return detail::exp_scaled(l.hi, l.lo) * detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo);
            }
            auto const shift = detail::mul(k.a1, detail::log_positive(yy));
            detail::scaled const I = integral(k, detail::stable_integrand::density, shift);
            if (I.hi == 0)
                return I;
            // α/(π |α − 1| y)
            detail::scaled const factor = detail::normalized(a, 0.0) * detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo)
                                        / (detail::normalized(detail::magnitude(k.am1.hi), k.am1.hi < 0 ? -k.am1.lo : k.am1.lo) * detail::normalized(yy.hi, yy.lo));
            return I * factor;
        }

        // ∫_0^Φ de l'intégrande, coupé au pic g = 1 (log g = 0, g monotone en φ).
        static constexpr detail::scaled integral(detail::stable_kernel const& k, detail::stable_integrand kind, detail::sd::double_double shift) noexcept
        {
            namespace sd = detail::sd;
            if (!(k.top.hi > 0))
                return {0.0, 0.0, 0};
            auto const f = [&](sd::double_double phi) { return detail::stable_value(kind, k.log_g(phi, shift)); };
            // bornes ouvertes : on reste à 2^-60 Φ des extrémités pour la recherche du pic
            sd::double_double const eps{sd::fmul(k.top.hi, 0x1p-60), 0.0};
            auto const lo_end = eps, hi_end = detail::sum(k.top, detail::negated(eps));
            double const l0 = k.log_g(lo_end, shift).hi, l1 = k.log_g(hi_end, shift).hi;
            if ((l0 > 0) == (l1 > 0) || !(l0 == l0) || !(l1 == l1))
                return detail::adaptive_integral(f, {0.0, 0.0}, k.top, {0.0, 0.0, 0});
            sd::double_double a = lo_end, b = hi_end;                     // log g change de signe dans (a, b)
            bool const increasing = l1 > l0;
            for (int i = 0; i < 80; ++i)
            {
                auto const w = detail::sum(b, detail::negated(a));
                auto const m = detail::sum(a, {sd::fmul(w.hi, 0.5), sd::fmul(w.lo, 0.5)});
                double const lm = k.log_g(m, shift).hi;
                if ((lm > 0) == increasing)
                    b = m;
                else
                    a = m;
                if (!(w.hi > sd::fmul(k.top.hi, 0x1p-50)))
                    break;
            }
            // de part et d'autre du pic, morceaux géométriques depuis le pic (premier : 1/8 de sa
            // largeur, puis ×4) : la masse se concentre près de lui (largeur ~ y^{−α} dans les
            // queues), où des nœuds grossiers ne verraient qu'un intégrande nul ; arrêt quand
            // l'intégrande décroît et que le morceau devient négligeable.
            // largeur du pic : 1/|∂ log g/∂φ| au point de coupure (différence centrée, pas resserré au besoin)
            double width = 0;
            {
                double h = sd::fmul(k.top.hi, 0x1p-24);
                for (int i = 0; i < 4; ++i)
                {
                    double const up = k.log_g(detail::sum(a, {h, 0.0}), shift).hi, down = k.log_g(detail::sum(a, {-h, 0.0}), shift).hi;
                    double const slope = detail::magnitude(up - down) / sd::fmul(2.0, h);
                    width = slope > 0 ? 1.0 / slope : k.top.hi;
                    if (width >= sd::fmul(64.0, h) || !(width > 0))
                        break;
                    h = width / 256;
                }
            }
            detail::scaled outer{0.0, 0.0, 0};                             // masse déjà intégrée (référence des tolérances)
            auto const side = [&](sd::double_double end) {
                bool const forward = detail::sum(end, detail::negated(a)).hi > 0;
                // largeurs : 1/8 de celle du pic, puis bornes à 4× la distance au pic
                auto const piece_width = [&](sd::double_double t) {
                    double const d = detail::magnitude(detail::sum(t, detail::negated(a)).hi);
                    return d == 0 ? width / 8 : sd::fmul(3.0, d);
                };
                detail::scaled previous = f(a), last_piece{0.0, 0.0, 0};
                bool first = true;
                auto const piece = [&](sd::double_double lo, sd::double_double hi, detail::scaled const& total) {
                    last_piece = detail::adaptive_integral(f, lo, hi, detail::add_scaled(outer, total));
                    return last_piece;
                };
                // arrêt quand l'intégrande décroît et que le dernier morceau et f(t) sont négligeables
                auto const stop = [&](sd::double_double t, detail::scaled const& total, sd::double_double) {
                    if (first)
                    {
                        first = false;
                        return false;
                    }
                    detail::scaled const at = f(t);
                    bool const decreasing = at.hi == 0 || (previous.hi != 0 && (at.e < previous.e || (at.e == previous.e && at.hi < previous.hi)));
                    previous = at;
                    return decreasing && total.hi != 0 && (last_piece.hi == 0 || last_piece.e < total.e - 100) && (at.hi == 0 || at.e < total.e - 100);
                };
                detail::scaled const total = detail::march_integral(a, end, forward, piece_width, piece, stop, 80);
                outer = detail::add_scaled(outer, total);
                return total;
            };
            detail::scaled const right = side(k.top);
            return detail::add_scaled(side({0.0, 0.0}), right);
        }

        constexpr T distribution(T x, bool upper) const noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (v == std::numeric_limits<double>::infinity() || v == -std::numeric_limits<double>::infinity())
                return (v > 0) != upper ? T(1) : T(0);
            return static_cast<T>(detail::rounded(reduced_part(reduced(v), upper)));
        }

        // y tel que F(y) = t (ou S(y) = t), t <= ½ : solveur commun (detail/root.hpp) ; hors
        // encadrement, y ± 2 max(1, |y|) ; bissection géométrique sur un encadrement large d'un
        // seul signe (queues lourdes) ; demi-droite pour α < 1, |β| = 1.
        constexpr detail::sd::double_double solve(double t, bool upper) const
        {
            namespace sd = detail::sd;
            auto const evaluate = [&](sd::double_double y) {
                detail::part_value r{reduced_part(y, upper), {0.0, 0.0, 0}, y};
                if (r.value.hi != 0)
                    r.slope = reduced_density(y);
                return r;
            };
            detail::root_bounds b{.start = 0};
            if (alpha_ < 1 && beta_ == 1)
                b = {.start = 1, .lo = 0};
            else if (alpha_ < 1 && beta_ == -1)
                b = {.start = -1, .hi = 0};
            auto const outward = [](double y, int side) {
                double const mag = detail::magnitude(y) < 1 ? 1.0 : detail::magnitude(y);
                return side > 0 ? y - sd::fmul(2.0, mag) : y + sd::fmul(2.0, mag);
            };
            auto const midpoint = [](double lo, double hi) {
                if (lo == 0 || hi == 0 || (lo < 0) != (hi < 0) || hi - lo < 1)
                    return sd::fmul(lo, 0.5) + sd::fmul(hi, 0.5);
                return lo > 0 ? wsl::special::sqrt(lo) * wsl::special::sqrt(hi) : -(wsl::special::sqrt(-lo) * wsl::special::sqrt(-hi));
            };
            return detail::solve_part(evaluate, !upper, t, b, outward, midpoint);
        }

        constexpr T inverse(double t, bool upper) const
        {
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::random::quantile : probabilité hors de [0, 1]");
            auto const sup = support();
            if (t == 0)
                return upper ? sup.upper : sup.lower;
            if (t == 1)
                return upper ? sup.lower : sup.upper;
            auto y = t <= 0.5 ? solve(t, upper) : solve(1.0 - t, !upper);
            if (alpha_ == 1 && beta_ != 0)
                y = detail::sum(y, log_shift_);
            if (y.hi == 0)
                return delta_;
            return static_cast<T>(detail::affine(double(delta_), double(gamma_), y.hi, y.lo));
        }

        T alpha_, beta_, gamma_, delta_;
        detail::stable_kernel positive_{}, negative_{};
        detail::sd::double_double log_shift_{0.0, 0.0};
    };
}
