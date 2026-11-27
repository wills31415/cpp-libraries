#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Loi de Student à ν > 0 degrés de liberté : densité
//   Γ((ν + 1)/2) / (√(νπ) Γ(ν/2)) · (1 + t²/ν)^{−(ν+1)/2}.
//
//   wsl::probability::student_t T({.degrees_of_freedom = 5});
//   T.cdf(2.0); T.quantile(0.975); T(g);
//
// Répartition : P(T > |t|) = I_x(ν/2, 1/2)/2, x = ν/(ν + t²) et 1 − x = t²/(ν + t²)
// en double-double, la bêta incomplète prise à argument double-double (premier
// ordre) ; quantile : l'inverse de la bêta en (x, 1 − x) double-double, puis
// t² = ν (1 − x)/x, une racine, un arrondi. Densité : Γ((ν+1)/2)/Γ(ν/2) sans
// compensation (detail::log_gamma_ratio), (1 + t²/ν)^{−(ν+1)/2} par log1p en
// double-double. Moments au sens usuel (NaN là où ils n'existent pas : moyenne
// pour ν <= 1, variance infinie pour 1 < ν <= 2…).
//
// Tirage (approché, bit-exact) : Z / √(V/ν), Z normale (ziggourat), V de loi χ²(ν)
// (gamma, Marsaglia–Tsang).

#include <wsl/probability/laws/basic.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/loader.hpp>
#include <wsl/probability/laws/gamma.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/laws/normal.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_real T = double>
    class student_t : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T degrees_of_freedom;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.degrees_of_freedom = degrees_of_freedom()}; }

        constexpr student_t(parameters p) : nu_(p.degrees_of_freedom)
        {
            namespace fp = detail::sd;
            if (!(nu_ > 0) || !detail::is_finite(nu_))
                throw std::domain_error("wsl::probability::student_t : il faut degrees_of_freedom > 0, fini");
            double const h = detail::sd::fmul(double(nu_), 0.5);
            double const a = h >= 1 ? h : h + 1;                            // Marsaglia–Tsang pour V/2 ~ Γ(ν/2)
            d_ = a - 1.0 / 3;
            c_ = 1 / wsl::special::sqrt(fp::fmul(9.0, d_));
        }

        constexpr T degrees_of_freedom() const noexcept { return nu_; }
        constexpr interval<T> support() const noexcept
        {
            return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()};
        }
        constexpr T mean() const noexcept { return nu_ > 1 ? T(0) : std::numeric_limits<T>::quiet_NaN(); }
        constexpr T median() const noexcept { return 0; }
        constexpr T mode() const noexcept { return 0; }
        // ν/(ν − 2) (ν > 2), +∞ (1 < ν <= 2).
        constexpr T variance() const noexcept
        {
            double const v = nu_;
            if (v <= 1)
                return std::numeric_limits<T>::quiet_NaN();
            if (v <= 2)
                return std::numeric_limits<T>::infinity();
            auto const [d_hi, d_lo] = detail::sd::two_sum(v, -2.0);
            return static_cast<T>(detail::rounded(detail::normalized(v, 0.0) / detail::normalized(d_hi, d_lo)));
        }
        constexpr T skewness() const noexcept { return nu_ > 3 ? T(0) : std::numeric_limits<T>::quiet_NaN(); }
        constexpr T excess_kurtosis() const noexcept
        {
            double const v = nu_;
            if (v <= 2)
                return std::numeric_limits<T>::quiet_NaN();
            if (v <= 4)
                return std::numeric_limits<T>::infinity();
            return static_cast<T>(6.0 / (v - 4.0));
        }
        // (ν + 1)/2 (ψ((ν + 1)/2) − ψ(ν/2)) + ½ log ν + log B(ν/2, ½), avec ψ(y) = log y − 1/(2y) + δ'(y)
        // et log B(h, ½) = log π/2 + lgamma(h) − lgamma(h + ½) : sans compensation.
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            double const h = detail::sd::fmul(double(nu_), 0.5);
            auto const [hp_hi, hp_lo] = fp::two_sum(h, 0.5);                // h + ½
            // ψ(h + ½) − ψ(h) = log1p(1/(2h)) − 1/(2(h + ½)) + 1/(2h) + δ'(h + ½) − δ'(h)
            double const r = 0.5 / h;
            auto const [m, m_error] = fp::two_product(r, h);
            fp::double_double const inv2h{r, ((0.5 - m) - m_error) / h};    // 1/(2h)
            fp::double_double diff = detail::log1p_positive(inv2h.hi, inv2h.lo);
            fp::double_double const inv2hp = fp::quotient_dd({0.5, 0.0}, {hp_hi, hp_lo});
            diff = detail::sum(diff, {-inv2hp.hi, -inv2hp.lo});
            diff = detail::sum(diff, inv2h);
            auto const dp1 = detail::stirling_error_derivative(hp_hi), dp0 = detail::stirling_error_derivative(h);
            diff = detail::sum(diff, dp1);
            diff = detail::sum(diff, {-dp0.hi, -dp0.lo});
            fp::double_double e = detail::mul({hp_hi, hp_lo}, diff);
            auto const ln = detail::log_precise(double(nu_));
            e = detail::sum(e, {fp::fmul(ln.hi, 0.5), fp::fmul(ln.lo, 0.5)});
            e = detail::sum(e, {fp::fmul(fp::log_pi_hi, 0.5), fp::fmul(fp::log_pi_lo, 0.5)});
            e = detail::sum(e, fp::log_gamma_ratio(0.5, h));                 // lgamma(h) − lgamma(h + ½)
            return static_cast<T>(e.hi + e.lo);
        }

        constexpr T pdf(T t) const noexcept
        {
            auto const& X = *this;
            if (!(t == t))
                return t;
            auto const [e_hi, e_lo] = X.log_density_dd(double(t));
            if (e_hi < -1500)
                return 0;
            return static_cast<T>(detail::rounded(detail::exp_scaled(e_hi, e_lo)));
        }
        constexpr T log_pdf(T t) const noexcept
        {
            auto const& X = *this;
            if (!(t == t))
                return t;
            auto const [e_hi, e_lo] = X.log_density_dd(double(t));
            return static_cast<T>(e_hi + e_lo);
        }

        constexpr T cdf(T t) const noexcept { auto const& X = *this; return X.distribution(double(t), false); }
        constexpr T ccdf(T t) const noexcept { auto const& X = *this; return X.distribution(double(t), true); }

        constexpr T quantile(double p) const
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return -std::numeric_limits<T>::infinity();
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            if (p == 0.5)
                return 0;
            double const tail = p < 0.5 ? p : 1.0 - p;                      // exact
            double const nu = X.nu_;
            auto const [x, y] = fp::beta_inverse_xy(detail::sd::fmul(nu, 0.5), 0.5, 2.0 * tail, false);   // I_x(ν/2, ½) = 2 tail
            if (x.hi < 0x1p-900)                                            // x perdu sous les sous-normaux, t encore
            {                                                               // représentable : I_x(a, ½) = xᵃ/(a B(a, ½))
                double const a = detail::sd::fmul(nu, 0.5);                                    // à 2^-900 près, d'où log x, puis
                auto l = detail::log_precise(2.0 * tail);                   // log t = (log ν − log x)/2
                l = detail::sum(l, detail::log_precise(a));
                l = detail::sum(l, detail::log_beta(a, 0.5));
                fp::double_double const lx = fp::quotient_dd(l, {a, 0.0});
                auto const ln = detail::log_precise(nu);
                auto const d = detail::sum(ln, {-lx.hi, -lx.lo});
                double const lt_hi = fp::fmul(d.hi, 0.5), lt_lo = fp::fmul(d.lo, 0.5);
                double const t = lt_hi > 1400 ? std::numeric_limits<double>::infinity() : detail::rounded(detail::exp_scaled(lt_hi, lt_lo));
                return static_cast<T>(p < 0.5 ? -t : t);
            }
            // t² = ν y / x
            detail::scaled const t2 = detail::normalized(nu, 0.0) * detail::normalized(y.hi, y.lo) / detail::normalized(x.hi, x.lo);
            double const t = detail::rounded(detail::sqrt(t2));
            return static_cast<T>(p < 0.5 ? -t : t);
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        constexpr T cquantile(double q) const
        {
            auto const& X = *this;
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::probability::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return wsl::probability::quantile(X, 1.0 - q);                                // 1 − q exact
            return -wsl::probability::quantile(X, q);                                         // loi symétrique
        }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            double const z = detail::ziggurat_normal(g);
            double const h = double(X.nu_) / 2;
            double w = detail::marsaglia_tsang(X.d_, X.c_, g);              // Γ(h) (ou Γ(h + 1))
            if (h < 1)
                w = fp::fmul(w, wsl::special::exp(wsl::special::log(detail::uniform_positive(g)) / h));
            return static_cast<T>(z / wsl::special::sqrt(w / h));           // V/ν = 2W/(2h) = W/h
        }

        friend constexpr bool operator==(student_t const& a, student_t const& b) noexcept { return a.nu_ == b.nu_; }

    private:
        // log f(t) = lgamma(h + ½) − lgamma(h) − ½ log(νπ) − (h + ½) log1p(t²/ν), h = ν/2.
        constexpr detail::sd::double_double log_density_dd(double t) const noexcept
        {
            namespace fp = detail::sd;
            double const nu = nu_, h = detail::sd::fmul(nu, 0.5);
            fp::double_double const ratio = fp::log_gamma_ratio(0.5, h);    // lgamma(h) − lgamma(h + ½)
            fp::double_double e{-ratio.hi, -ratio.lo};
            auto const ln = detail::log_precise(nu);
            e = detail::sum(e, {fp::fmul(ln.hi, -0.5), fp::fmul(ln.lo, -0.5)});
            e = detail::sum(e, {fp::fmul(fp::log_pi_hi, -0.5), fp::fmul(fp::log_pi_lo, -0.5)});
            double const a = detail::magnitude(t);
            if (a == 0)
                return e;
            if (a == std::numeric_limits<double>::infinity())
                return {-std::numeric_limits<double>::infinity(), 0.0};
            // w = t²/ν, à l'échelle ; log1p(w) = log w + log1p(1/w) pour w immense.
            detail::scaled const w = detail::normalized(a, 0.0) * detail::normalized(a, 0.0) / nu;
            fp::double_double l{};
            if (w.e > 60)
            {
                l = detail::log(w);
                auto const [iw_hi, iw_lo] = detail::unscaled(detail::normalized(1.0, 0.0) / w);
                l = detail::sum(l, {iw_hi, iw_lo});
            }
            else
            {
                auto const [w_hi, w_lo] = detail::unscaled(w);
                l = detail::log1p_positive(w_hi, w_lo);
            }
            auto const [hp_hi, hp_lo] = fp::two_sum(h, 0.5);
            fp::double_double const term = detail::mul({hp_hi, hp_lo}, l);
            return detail::sum(e, {-term.hi, -term.lo});
        }

        // P(T <= t) (ou P(T > t)) : queue P(T > |t|) = I_x(ν/2, ½)/2.
        constexpr T distribution(double t, bool upper) const noexcept
        {
            namespace fp = detail::sd;
            if (!(t == t))
                return static_cast<T>(t);
            if (t == 0)
                return T(0.5);
            double const nu = nu_, a = detail::magnitude(t);
            bool const tail_side = (t < 0) != upper;                        // la valeur cherchée est la queue
            if (a == std::numeric_limits<double>::infinity())
                return tail_side ? T(0) : T(1);
            // x = ν/(ν + t²), y = t²/(ν + t²), à l'échelle puis en double-double.
            detail::scaled const t2 = detail::normalized(a, 0.0) * detail::normalized(a, 0.0);
            detail::scaled const nus = detail::normalized(nu, 0.0);
            detail::scaled x{}, y{};
            if (t2.e > nus.e + 1100)                                        // t² immense : x = ν/t²
            {
                x = nus / t2;
                y = detail::normalized(1.0, 0.0);
            }
            else if (nus.e > t2.e + 1100)                                    // t² minuscule devant ν
            {
                x = detail::normalized(1.0, 0.0);
                y = t2 / nus;
            }
            else
            {
                auto const [th, tl] = detail::unscaled(t2);
                auto const [s_hi, s_lo] = detail::sum({nu, 0.0}, {th, tl});
                detail::scaled const s = detail::normalized(s_hi, s_lo);
                x = nus / s;
                y = t2 / s;
            }
            auto const xv = detail::unscaled(x), yv = detail::unscaled(y);
            if (xv.hi == 0)                                                  // queue < 2^-1075
                return tail_side ? T(0) : T(1);
            if (yv.hi == 0)
                return T(0.5);
            auto [value, complement] = fp::incomplete_beta_dd(detail::sd::fmul(nu, 0.5), 0.5, xv, yv);
            // queue = I_x/2 : value = I_x si !complement (précise), sinon value = 1 − I_x.
            value.e -= 1;                                                    // ½ ·
            if (!complement)
                return static_cast<T>(tail_side ? detail::rounded(value) : fp::one_minus(value));
            // queue = ½ − value/2 (près de ½)
            auto const [v_hi, v_lo] = detail::unscaled(value);
            auto const [q, q_error] = fp::two_sum(0.5, -v_hi);
            double const tail = q + (q_error - v_lo);
            if (tail_side)
                return static_cast<T>(tail);
            auto const [r, r_error] = fp::two_sum(0.5, v_hi);
            return static_cast<T>(r + (r_error + v_lo));
        }

        T nu_;
        double d_ = 0, c_ = 0;
    };
}
