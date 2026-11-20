#pragma once

// Loi de Fisher (F de Snedecor) à d₁ > 0 et d₂ > 0 degrés de liberté : quotient
// (χ²(d₁)/d₁)/(χ²(d₂)/d₂).
//
//   wsl::random::fisher_f F({.d1 = 5, .d2 = 10});
//   F.cdf(3.3); F.quantile(0.95); F(g);
//
// Par la loi bêta (a = d₁/2, b = d₂/2) : u = d₁x/(d₁x + d₂), 1 − u = d₂/(d₁x + d₂) en
// double-double ; répartition I_u(a, b) (bêta incomplète à argument double-double),
// densité préfixe(u)/x = uᵃ(1 − u)ᵇ/(B(a, b) x) en logarithme double-double, un
// seul arrondi ; quantile x = d₂ u/(d₁ (1 − u)), u et 1 − u par l'inverse de la
// bêta. Entropie : celle de la bêta + log(b/a) + 2(ψ(a + b) − ψ(b)), ψ réécrite avec
// δ' (sans compensation). Moments au sens usuel (NaN là où ils n'existent pas).
//
// Tirage (approché, bit-exact) : (G₁/a)/(G₂/b), G₁, G₂ de lois gamma(a), gamma(b).

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/beta.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/gamma.hpp>
#include <wsl/random/law.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    template <law_real T = double>
    class fisher_f : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T d1, d2;
        };

        constexpr fisher_f(parameters p) : d1_(p.d1), d2_(p.d2)
        {
            namespace fp = detail::sd;
            if (!(d1_ > 0) || !detail::is_finite(d1_) || !(d2_ > 0) || !detail::is_finite(d2_))
                throw std::domain_error("wsl::random::fisher_f : il faut d1 > 0 et d2 > 0, finis");
            double const a = half(d1_), b = half(d2_);
            double const sa = a >= 1 ? a : a + 1, sb = b >= 1 ? b : b + 1;
            d_a_ = sa - 1.0 / 3;
            c_a_ = 1 / wsl::special::sqrt(fp::fmul(9.0, d_a_));
            d_b_ = sb - 1.0 / 3;
            c_b_ = 1 / wsl::special::sqrt(fp::fmul(9.0, d_b_));
        }

        constexpr T d1() const noexcept { return d1_; }
        constexpr T d2() const noexcept { return d2_; }
        constexpr interval<T> support() const noexcept { return {T(0), std::numeric_limits<T>::infinity()}; }

        // d₂/(d₂ − 2), d₂ > 2.
        constexpr T mean() const noexcept
        {
            double const d2 = d2_;
            if (d2 <= 2)
                return std::numeric_limits<T>::quiet_NaN();
            auto const [m_hi, m_lo] = detail::sd::two_sum(d2, -2.0);
            return static_cast<T>(detail::rounded(detail::normalized(d2, 0.0) / detail::normalized(m_hi, m_lo)));
        }
        // 2 d₂² (d₁ + d₂ − 2) / (d₁ (d₂ − 2)² (d₂ − 4)), d₂ > 4.
        constexpr T variance() const noexcept
        {
            namespace fp = detail::sd;
            double const d1 = d1_, d2 = d2_;
            if (d2 <= 2)
                return std::numeric_limits<T>::quiet_NaN();
            if (d2 <= 4)
                return std::numeric_limits<T>::infinity();
            auto const [s1, e1] = fp::two_sum(d1, d2);
            auto const [n_hi, n_lo] = detail::sum({s1, e1}, {-2.0, 0.0});
            auto const [m2_hi, m2_lo] = fp::two_sum(d2, -2.0);
            auto const [m4_hi, m4_lo] = fp::two_sum(d2, -4.0);
            detail::scaled const m2 = detail::normalized(m2_hi, m2_lo);
            detail::scaled const num = detail::normalized(2.0, 0.0) * detail::normalized(d2, 0.0) * detail::normalized(d2, 0.0) * detail::normalized(n_hi, n_lo);
            detail::scaled const den = detail::normalized(d1, 0.0) * m2 * m2 * detail::normalized(m4_hi, m4_lo);
            return static_cast<T>(detail::rounded(num / den));
        }
        // (d₁ − 2)/d₁ · d₂/(d₂ + 2), d₁ > 2 ; 0 sinon.
        constexpr T mode() const noexcept
        {
            namespace fp = detail::sd;
            double const d1 = d1_, d2 = d2_;
            if (d1 <= 2)
                return 0;
            auto const [a_hi, a_lo] = fp::two_sum(d1, -2.0);
            auto const [b_hi, b_lo] = fp::two_sum(d2, 2.0);
            return static_cast<T>(detail::rounded(detail::normalized(a_hi, a_lo) * detail::normalized(d2, 0.0)
                                                  / (detail::normalized(d1, 0.0) * detail::normalized(b_hi, b_lo))));
        }
        constexpr T skewness() const noexcept
        {
            double const d1 = d1_, d2 = d2_;
            if (d2 <= 6)
                return std::numeric_limits<T>::quiet_NaN();
            return static_cast<T>((detail::sd::fmul(2.0, d1) + d2 - 2) * wsl::special::sqrt(8 * (d2 - 4)) / ((d2 - 6) * wsl::special::sqrt(d1 * (d1 + d2 - 2))));
        }
        constexpr T excess_kurtosis() const noexcept
        {
            namespace fp = detail::sd;
            double const d1 = d1_, d2 = d2_;
            if (d2 <= 8)
                return std::numeric_limits<T>::quiet_NaN();
            double const num = fp::fmul(fp::fmul(d1, fp::fmul(5.0, d2) - 22), d1 + d2 - 2) + fp::fmul(fp::fmul(d2 - 4, d2 - 2), d2 - 2);
            return static_cast<T>(12 * num / (d1 * (d2 - 6) * (d2 - 8) * (d1 + d2 - 2)));
        }
        // H(bêta(a, b)) + log(b/a) + 2(ψ(a + b) − ψ(b)), ψ(n) − ψ(b) =
        // log1p(a/b) − 1/(2n) + 1/(2b) + δ'(n) − δ'(b).
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            double const a = half(d1_), b = half(d2_);
            auto h = detail::beta_entropy(a, b);
            auto const la = detail::log_precise(a), lb = detail::log_precise(b);
            h = detail::sum(h, lb);
            h = detail::sum(h, {-la.hi, -la.lo});
            auto const [n_hi, n_lo] = fp::two_sum(a, b);
            fp::double_double const r = fp::quotient_dd({a, 0.0}, {b, 0.0});
            fp::double_double d = detail::log1p_positive(r.hi, r.lo);
            fp::double_double const inv2n = fp::quotient_dd({0.5, 0.0}, {n_hi, n_lo});
            fp::double_double const inv2b = fp::quotient_dd({0.5, 0.0}, {b, 0.0});
            d = detail::sum(d, {-inv2n.hi, -inv2n.lo});
            d = detail::sum(d, inv2b);
            auto const dn = detail::stirling_error_derivative(n_hi), dbb = detail::stirling_error_derivative(b);
            d = detail::sum(d, dn);
            d = detail::sum(d, {-dbb.hi, -dbb.lo});
            h = detail::sum(h, {fp::fmul(d.hi, 2.0), fp::fmul(d.lo, 2.0)});
            return static_cast<T>(h.hi + h.lo);
        }

        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0) || v == std::numeric_limits<double>::infinity())
                return X.at_zero(v, false);
            auto const l = X.log_density_dd(v);
            if (l.hi < -1500)
                return 0;
            if (l.hi > 1400)
                return std::numeric_limits<T>::infinity();
            return static_cast<T>(detail::rounded(detail::exp_scaled(l.hi, l.lo)));
        }
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0) || v == std::numeric_limits<double>::infinity())
                return X.at_zero(v, true);
            auto const l = X.log_density_dd(v);
            return static_cast<T>(l.hi + l.lo);
        }

        constexpr T cdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, true); }

        // x = d₂ u / (d₁ (1 − u)), I_u(a, b) = p.
        constexpr T quantile(double p) const
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::random::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return 0;
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            return X.inverse(p, false);
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        constexpr T cquantile(double q) const
        {
            auto const& X = *this;
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::random::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return wsl::random::quantile(X, 1.0 - q);                                // 1 − q exact
            if (q == 0)
                return std::numeric_limits<T>::infinity();
            return X.inverse(q, true);
        }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            double const a = half(X.d1_), b = half(X.d2_);
            double ga = detail::marsaglia_tsang(X.d_a_, X.c_a_, g);
            double gb = detail::marsaglia_tsang(X.d_b_, X.c_b_, g);
            if (a < 1)
                ga = fp::fmul(ga, wsl::special::exp(wsl::special::log(detail::uniform_positive(g)) / a));
            if (b < 1)
                gb = fp::fmul(gb, wsl::special::exp(wsl::special::log(detail::uniform_positive(g)) / b));
            return static_cast<T>((ga / a) / (gb / b));
        }

        friend constexpr bool operator==(fisher_f const& p, fisher_f const& q) noexcept { return p.d1_ == q.d1_ && p.d2_ == q.d2_; }

    private:
        // x tel que I_u(a, b) = t (1 − I_u(a, b) = t si complement), t dans (0, 1).
        constexpr T inverse(double t, bool complement) const noexcept
        {
            namespace fp = detail::sd;
            double const a = half(d1_), b = half(d2_);
            auto const [u, w] = fp::beta_inverse_xy(a, b, t, complement);
            // u (ou 1 − u) sous 2^-960 n'a plus tous ses bits : I_u(a, b) = uᵃ/(a B(a, b)) (1 + O(b u)),
            // d'où log u = (log p + log a + log B)/a, et x = (d₂/d₁) u/(1 − u) en un arrondi.
            bool const tiny_u = u.hi < 0x1p-960 && fp::fmul(u.hi, b) < 0x1p-900;
            bool const tiny_w = !tiny_u && w.hi < 0x1p-960 && fp::fmul(w.hi, a) < 0x1p-900;
            if (tiny_u || tiny_w)
            {
                double const s = tiny_u ? a : b;
                fp::double_double const lt = detail::log_precise(t), l1 = detail::log1p_dd(-t);   // log t, log(1 − t)
                auto l = detail::sum(tiny_u == complement ? l1 : lt, detail::log_precise(s));      // log de la queue de u (ou de 1 − u)
                l = detail::sum(l, detail::log_beta(a, b));
                l = fp::quotient_dd(l, {s, 0.0});                           // log u (ou log(1 − u))
                if (!tiny_u)
                    l = {-l.hi, -l.lo};
                l = detail::sum(l, detail::log_precise(double(d2_)));
                auto const ld1 = detail::log_precise(double(d1_));
                l = detail::sum(l, {-ld1.hi, -ld1.lo});
                if (l.hi < -1500)
                    return 0;
                if (l.hi > 1400)
                    return std::numeric_limits<T>::infinity();
                return static_cast<T>(detail::rounded(detail::exp_scaled(l.hi, l.lo)));
            }
            if (u.hi == 0)
                return 0;
            if (w.hi == 0)
                return std::numeric_limits<T>::infinity();
            detail::scaled const x = detail::normalized(double(d2_), 0.0) * detail::normalized(u.hi, u.lo)
                                   / (detail::normalized(double(d1_), 0.0) * detail::normalized(w.hi, w.lo));
            return static_cast<T>(detail::rounded(x));
        }

        static constexpr double half(T d) noexcept { return detail::sd::fmul(double(d), 0.5); }

        // u = d₁x/(d₁x + d₂), 1 − u = d₂/(d₁x + d₂), en double-double.
        struct uv
        {
            detail::sd::double_double u, w;
        };
        constexpr uv arguments(double x) const noexcept
        {
            detail::scaled const p = detail::normalized(double(d1_), 0.0) * detail::normalized(x, 0.0);   // d₁x
            detail::scaled const q = detail::normalized(double(d2_), 0.0);
            if (p.e > q.e + 1100)
                return {{1.0, 0.0}, detail::unscaled(q / p)};
            if (q.e > p.e + 1100)
                return {detail::unscaled(p / q), {1.0, 0.0}};
            auto const [p_hi, p_lo] = detail::unscaled(p);
            auto const [s_hi, s_lo] = detail::sum({p_hi, p_lo}, {double(d2_), 0.0});
            detail::scaled const s = detail::normalized(s_hi, s_lo);
            return {detail::unscaled(p / s), detail::unscaled(q / s)};
        }

        // log f(x) = log(uᵃ (1 − u)ᵇ / B(a, b)) − log x.
        constexpr detail::sd::double_double log_density_dd(double x) const noexcept
        {
            namespace fp = detail::sd;
            double const a = half(d1_), b = half(d2_);
            if (a == 1)                                                     // d₁ = 2 : f = (1 + r)^{−(b+1)}, r = 2x/d₂
            {                                                               // (log f proche de 0 près de 0 : sans compensation)
                detail::scaled const r = detail::normalized(double(d1_), 0.0) * detail::normalized(x, 0.0) / double(d2_);
                auto const [b1_hi, b1_lo] = fp::two_sum(b, 1.0);
                fp::double_double l{};
                if (r.e < -1000)                                            // log1p(r) = r
                    l = {detail::rounded(detail::normalized(b1_hi, b1_lo) * r), 0.0};
                else if (r.e > 1000)
                    l = detail::mul({b1_hi, b1_lo}, detail::log(r));
                else
                {
                    auto const [r_hi, r_lo] = detail::unscaled(r);
                    l = detail::mul({b1_hi, b1_lo}, detail::log1p_positive(r_hi, r_lo));
                }
                return {-l.hi, -l.lo};
            }
            auto const [u, w] = arguments(x);
            if (u.hi == 0 || w.hi == 0)
                return {-std::numeric_limits<double>::infinity(), 0.0};
            auto l = fp::beta_log_prefix(a, b, u, w);
            auto const lx = detail::log_precise(x);
            return detail::sum(l, {-lx.hi, -lx.lo});
        }

        // Densité (ou log-densité) en x <= 0 ou x = +∞.
        constexpr T at_zero(double v, bool logarithm) const noexcept
        {
            double const d1 = d1_;
            if (v == 0 && d1 < 2)
                return std::numeric_limits<T>::infinity();
            if (v == 0 && d1 == 2)                                          // f(0) = 1 (d₁ = 2)
                return logarithm ? T(0) : T(1);
            return logarithm ? -std::numeric_limits<T>::infinity() : T(0);
        }

        constexpr T distribution(T x, bool upper) const noexcept
        {
            namespace fp = detail::sd;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0))
                return upper ? T(1) : T(0);
            if (v == std::numeric_limits<double>::infinity())
                return upper ? T(0) : T(1);
            auto const [u, w] = arguments(v);
            if (u.hi == 0)
                return upper ? T(1) : T(0);
            if (w.hi == 0)
                return upper ? T(0) : T(1);
            auto const [value, complement] = fp::incomplete_beta_dd(half(d1_), half(d2_), u, w);
            return static_cast<T>(complement == upper ? detail::rounded(value) : fp::one_minus(value));
        }

        T d1_, d2_;
        double d_a_ = 0, c_a_ = 0, d_b_ = 0, c_b_ = 0;
    };
}
