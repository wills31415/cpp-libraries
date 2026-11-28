#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Loi de Cauchy de position x₀ et d'échelle γ > 0 : densité 1/(πγ (1 + z²)),
// z = (x − x₀)/γ. Ni moyenne ni variance (elles n'existent pas : les membres
// correspondants sont absents, has_mean<cauchy> est faux).
//
//   wsl::probability::laws::cauchy C({.location = 0, .scale = 1});
//   C.cdf(3.0); C.quantile(0.99); C.entropy(); C(g);
//
// z en double-double (comme pour la normale). Répartition : P(X > x) = atan(1/z)/π
// pour z >= 1, 1/2 − atan(z)/π sinon (noyau atan de special_functions à argument
// double-double) ; quantile x₀ + γ tan(π(p − ½)) = x₀ − γ cot(πp) pour p <= ½,
// cos(πp) et sin(πp) en double-double, un seul arrondi.
//
// Tirage (approché, bit-exact) : x₀ ± γ cot(πW/2), W uniforme dans (0, 1), signe
// uniforme.

#include <wsl/probability/laws/basic.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/trigonometry.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/laws/normal.hpp>
#include <wsl/random/detail/portable_math.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
}

// Outils internes : dans wsl::probability::detail (un seul « detail » visible depuis probability).
namespace wsl::inline WSL_DEFAULTS::probability::detail
{
        // cot(π u), 0 < u <= ½, à l'échelle : cos(πu)/sin(πu) en double-double ; pour
        // u < 2^-30, (1/(πu)) (1 − (πu)²/3) (le terme suivant est sous 2^-120, et le
        // produit exact de sin(πu) ne vaut plus sous 2^-970).
        constexpr scaled cot_pi(double u) noexcept
        {
            if (u < 0x1p-30)
            {
                scaled const r = normalized(inv_pi_hi, inv_pi_lo) / u;
                double const pu = sd::fmul(sd::pi_hi, u);
                double const c = sd::fmul(pu, pu) / 3;
                auto const [c_hi, c_lo] = sd::two_sum(1.0, -c);
                return r * normalized(c_hi, c_lo);
            }
            auto const [s_hi, s_lo] = sd::sin_pi_abs(u);
            auto const [c_hi, c_lo] = sd::sin_pi_abs(0.5 - u);              // cos(πu), ½ − u exact
            return normalized(c_hi, c_lo) / normalized(s_hi, s_lo);
        }
}

namespace wsl::inline WSL_DEFAULTS::probability::laws
{

    template <law_real T = double>
    class cauchy : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T location = 0, scale = 1;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.location = location(), .scale = scale()}; }

        constexpr cauchy() : cauchy(parameters{}) {}

        constexpr cauchy(parameters p) : x0_(p.location), gamma_(p.scale)
        {
            if (!detail::is_finite(x0_) || !(gamma_ > 0) || !detail::is_finite(gamma_))
                throw std::domain_error("wsl::probability::laws::cauchy : il faut location fini et scale > 0, fini");
        }

        constexpr T location() const noexcept { return x0_; }
        constexpr T scale() const noexcept { return gamma_; }
        constexpr interval<T> support() const noexcept
        {
            return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()};
        }
        constexpr T median() const noexcept { return x0_; }
        constexpr T mode() const noexcept { return x0_; }
        // log(4πγ)
        constexpr T entropy() const noexcept
        {
            auto const l = detail::log_precise(double(gamma_));
            auto const e = detail::sum(l, {detail::log_4pi_hi, detail::log_4pi_lo});
            return static_cast<T>(e.hi + e.lo);
        }

        // 1/(πγ(1 + z²))
        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            if (!(x == x))
                return x;
            if (!detail::is_finite(x))
                return 0;
            detail::scaled const c = detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo) / double(X.gamma_);
            return static_cast<T>(detail::rounded(c / one_plus_square(X.standard(x).z)));
        }

        // −log π − log γ − log(1 + z²)
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            if (!(x == x))
                return x;
            if (!detail::is_finite(x))
                return -std::numeric_limits<T>::infinity();
            auto e = detail::log(one_plus_square(X.standard(x).z));
            auto const lg = detail::log_precise(double(X.gamma_));
            e = detail::sum(e, lg);
            e = detail::sum(e, {fp::log_pi_hi, fp::log_pi_lo});
            return static_cast<T>(-(e.hi + e.lo));
        }

        constexpr T cdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, true); }

        // x₀ + γ tan(π(p − ½)) : −cot(πp) pour p < ½, cot(π(1 − p)) sinon.
        constexpr T quantile(double p) const
        {
            auto const& X = *this;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return -std::numeric_limits<T>::infinity();
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            if (p == 0.5)
                return X.x0_;
            bool const lower = p < 0.5;
            detail::scaled const t = detail::cot_pi(lower ? p : 1.0 - p) * detail::normalized(double(X.gamma_), 0.0);
            return static_cast<T>(detail::affine_scaled(X.x0_, lower, t));
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        constexpr T cquantile(double q) const
        {
            auto const& X = *this;
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::probability::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return wsl::probability::quantile(X, 1.0 - q);                                // 1 − q exact
            if (q == 0)
                return std::numeric_limits<T>::infinity();
            detail::scaled const t = detail::cot_pi(q) * detail::normalized(double(X.gamma_), 0.0);   // x₀ + γ cot(πq)
            return static_cast<T>(detail::affine_scaled(X.x0_, false, t));
        }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            double const w = detail::uniform_positive(g);                   // |Z| = cot(πW/2)
            bool const negative = (random::random_bits(g) >> 63) != 0;
            detail::scaled const t = detail::cot_pi(detail::sd::fmul(w, 0.5)) * detail::normalized(double(X.gamma_), 0.0);
            return static_cast<T>(detail::affine_scaled(X.x0_, negative, t));
        }

        friend constexpr bool operator==(cauchy const&, cauchy const&) = default;

    private:
        // 1 + z², à l'échelle.
        // |x − x₀|/γ à l'échelle, sans saturation : la queue de Cauchy reste représentable bien
        // au-delà de |z| = 2^600 (detail::standardize sature).
        constexpr detail::standardized_scaled standard(double x) const noexcept { return detail::standardize_scaled(x, double(x0_), double(gamma_)); }

        static constexpr detail::scaled one_plus_square(detail::scaled const& z) noexcept
        {
            if (z.hi == 0)
                return detail::normalized(1.0, 0.0);
            detail::scaled const z2 = z * z;
            if (z2.e > 120)                                                  // 1 négligeable
                return z2;
            auto const [s_hi, s_lo] = detail::unscaled(z2);
            auto const [u_hi, u_lo] = detail::sum({1.0, 0.0}, {s_hi, s_lo});
            return detail::normalized(u_hi, u_lo);
        }

        // P(X <= x) (ou P(X > x)) : queue P(Z > |z|) = atan(1/|z|)/π (|z| >= 1), ½ − atan(|z|)/π sinon.
        constexpr T distribution(T x, bool upper) const noexcept
        {
            namespace fp = detail::sd;
            if (!(x == x))
                return x;
            if (!detail::is_finite(x))
                return (x < 0) != upper ? T(0) : T(1);
            auto const [negative, z] = standard(x);
            if (z.hi == 0 || z.e < -600)                                    // ½ ± z/π s'arrondit à ½
                return T(0.5);
            bool const tail_side = negative != upper;
            detail::scaled tail{};
            if (z.e >= 512)                                                 // atan(1/z)/π = (1/z)/π à 2^-1024 près ; et
                                                                            // 1/z par Dekker sortirait de son domaine (2^995)
                tail = detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo) / z;
            else if (z.e >= 0)
            {
                auto const [a, a_lo] = detail::unscaled(z);
                auto const [v_hi, v_lo] = fp::quotient_dd({1.0, 0.0}, {a, a_lo});   // 1/|z|
                fp::double_double const at = v_hi < 0x1p-27                    // atan v = v − v³/3 (le terme suivant : 2^-108)
                    ? fp::double_double{v_hi, v_lo - fp::fmul(fp::fmul(v_hi, v_hi), v_hi) / 3} : fp::atan_kernel(v_hi, v_lo);
                tail = detail::normalized(at.hi, at.lo) * detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo);
            }
            else
            {
                auto const [a, a_lo] = detail::unscaled(z);
                auto const at = fp::atan_kernel(a, a_lo);
                auto const q = detail::unscaled(detail::normalized(at.hi, at.lo) * detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo));
                auto const [t_hi, t_lo] = detail::sum({0.5, 0.0}, {-q.hi, -q.lo});
                tail = detail::normalized(t_hi, t_lo);
            }
            return static_cast<T>(tail_side ? detail::rounded(tail) : fp::one_minus(tail));
        }

        T x0_, gamma_;
    };
}
