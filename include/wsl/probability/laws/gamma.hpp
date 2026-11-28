#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Loi gamma de forme a > 0, de taux λ > 0 ou d'échelle θ = 1/λ :
// densité λ^a x^{a−1} e^{−λx} / Γ(a) sur [0, +∞).
//
//   wsl::probability::laws::gamma X({.rate = 2.0, .shape = 2.5});
//   wsl::probability::laws::gamma Y({.scale = 0.5, .shape = 2.5});        // la même loi
//   X.mean(); X.pdf(1.0); X.log_pdf(1.0); X(g);
//
// Le paramètre donné (taux ou échelle) est gardé tel quel : les calculs
// multiplient par λ ou divisent par θ, sans arrondir 1/θ.
//
// Densité : comme dgamma de R, par la masse de Poisson généralisée,
//   a >= 1 : λ p(a − 1, λx),   a < 1 : (a/x) p(a, λx),   p(ν, t) = t^ν e^{−t}/Γ(ν + 1),
// p(ν, t) = e^{−δ(ν) − bd0(ν, t)} / √(2πν) par la décomposition de Loader (voir
// detail/loader.hpp), un seul arrondi. Répartition et survie : P(a, λx) et Q(a, λx)
// de special_functions, λx en double-double (correction du premier ordre f(t) t_lo) ;
// quantile : l'inverse de special_functions en double-double, divisé par λ (ou multiplié
// par θ) avec un seul arrondi ; entropie réécrite avec l'erreur de Stirling δ :
//     a − log λ + lgamma(a) + (1 − a)ψ(a) = ½ log a + ½ log 2π + ½ − 1/(2a) + δ(a) + (1 − a)δ'(a) − log λ,
// sans les compensations de la forme usuelle.
//
// Tirage (approché, bit-exact) : Marsaglia et Tsang (2000), la normale par la
// ziggourat ; a < 1 : G(a + 1) U^{1/a}.

#include <wsl/probability/laws/basic.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/detail/rational.hpp>
#include <wsl/probability/detail/loader.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/laws/normal.hpp>
#include <wsl/random/detail/portable_math.hpp>

#include <limits>
#include <optional>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
}

// Outils internes : dans wsl::probability::detail (un seul « detail » visible depuis probability).
namespace wsl::inline WSL_DEFAULTS::probability::detail
{
        // −δ(ν) − bd0(ν, t), ν > 0, t > 0 à l'échelle (au-delà des doubles s'il le
        // faut : bd0 est homogène de degré 1).
        constexpr sd::double_double poisson_exponent(double nu, scaled t) noexcept
        {
            int shift = 0;
            if (t.e >= 1000)
                shift = -100;
            else if (t.e < -1000)
                shift = 100;
            double const nu_scaled = sd::scale_by_power_of_two(nu, shift);
            if (!is_finite(nu_scaled))                                      // t minuscule, ν immense : bd0 immense
                return {-std::numeric_limits<double>::infinity(), 0.0};
            t.e += shift;
            auto const [t_hi, t_lo] = unscaled(t);
            auto [b_hi, b_lo] = deviance(nu_scaled, t_hi, t_lo);
            b_hi = sd::scale_by_power_of_two(b_hi, -shift);
            b_lo = sd::scale_by_power_of_two(b_lo, -shift);
            auto const [d_hi, d_lo] = stirling_error(nu);
            auto const [s, s_error] = sd::two_sum(-d_hi, -b_hi);
            if (!is_finite(s))
                return {-std::numeric_limits<double>::infinity(), 0.0};
            auto const [hi, lo] = sd::fast_two_sum(s, s_error - (d_lo + b_lo));
            return {hi, lo};
        }

        // Gamma(a, 1), a >= 1 : Marsaglia et Tsang, d = a − 1/3, c = 1/√(9d).
        template <class G>
        constexpr double marsaglia_tsang(double d, double c, G& g)
        {
            for (;;)
            {
                double z = 0, v = 0;
                do
                {
                    z = ziggurat_normal(g);
                    v = 1.0 + sd::fmul(c, z);
                } while (v <= 0);
                v = sd::fmul(sd::fmul(v, v), v);
                double const u = random::uniform_real<double>(g, wsl::exact);
                double const z2 = sd::fmul(z, z);
                if (u < 1.0 - sd::fmul(0.0331, sd::fmul(z2, z2)))
                    return sd::fmul(d, v);
                if (wsl::random::detail::pm::log(u) < sd::fmul(0.5, z2) + sd::fmul(d, (1.0 - v) + wsl::random::detail::pm::log(v)))
                    return sd::fmul(d, v);
            }
        }
}

namespace wsl::inline WSL_DEFAULTS::probability::laws
{

    template <law_real T = double>
    class gamma : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        // Intensité .rate ou échelle .scale (l'une des deux au plus ; aucune : λ = 1), double ou
        // rnd::ratio(a, b) : exactes toutes deux (λ rationnel, voir detail/rational.hpp).
        struct parameters
        {
            std::optional<random::ratio> rate = std::nullopt, scale = std::nullopt;
            T shape;
        };
        using param_type = parameters;
        constexpr param_type param() const
        {
            auto const [r, s] = detail::parameters_of(lambda_);
            return {.rate = r, .scale = s, .shape = shape_};
        }

        constexpr gamma(parameters p)
            : gamma(p.shape, detail::rate_of(p.rate, p.scale, "wsl::probability::laws::gamma : il faut .rate ou .scale (pas les deux), > 0 et finis"))
        {
        }

        constexpr T shape() const noexcept { return shape_; }
        // λ et 1/λ arrondis (exacts quand ils sont des doubles).
        constexpr T rate() const noexcept { return static_cast<T>(rate_ != 0 ? rate_ : detail::rounded(lambda_.value())); }
        constexpr T scale() const noexcept
        {
            return static_cast<T>(scale_ != 0 ? scale_ : detail::rounded(lambda_.per(detail::normalized(1.0, 0.0))));
        }
        constexpr detail::exact_rate const& exact_rate() const noexcept { return lambda_; }
        constexpr interval<T> support() const noexcept { return {T(0), std::numeric_limits<T>::infinity()}; }
        constexpr T mean() const noexcept { return static_cast<T>(detail::rounded(per_rate(detail::normalized(shape_, 0.0)))); }
        constexpr T variance() const noexcept
        {
            return static_cast<T>(detail::rounded(per_rate(per_rate(detail::normalized(shape_, 0.0)))));
        }
        constexpr T mode() const noexcept
        {
            double const nu = double(shape_) - 1;
            return nu <= 0 ? T(0) : static_cast<T>(detail::rounded(per_rate(detail::normalized(nu, 0.0))));
        }
        constexpr T skewness() const noexcept { return static_cast<T>(2 / wsl::random::detail::pm::sqrt(double(shape_))); }
        constexpr T excess_kurtosis() const noexcept { return static_cast<T>(6 / double(shape_)); }

        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x, a = X.shape_;
            if (!(v >= 0))
                return v == v ? T(0) : x;
            if (v == 0)
                return a < 1 ? std::numeric_limits<T>::infinity() : a == 1 ? static_cast<T>(detail::rounded(X.per_scale(detail::normalized(1.0, 0.0)))) : T(0);
            if (v == std::numeric_limits<double>::infinity())
                return 0;
            detail::scaled const t = X.argument(v);
            detail::scaled const c = detail::normalized(detail::inv_sqrt_2pi_hi, detail::inv_sqrt_2pi_lo);
            if (a < 1)                                                      // e^E c √a / x
            {
                auto const [e_hi, e_lo] = detail::poisson_exponent(a, t);
                if (e_hi < -1500)
                    return 0;
                if (e_hi > 1400)
                    return std::numeric_limits<T>::infinity();
                return static_cast<T>(detail::rounded(detail::exp_scaled(e_hi, e_lo) * c * detail::sqrt(detail::normalized(a, 0.0)) / v));
            }
            double const nu = a - 1;                                        // exact sous 2^53
            if (nu == 0)                                                    // λ e^{−λx}
            {
                if (t.e > 11)
                    return 0;
                auto const [t_hi, t_lo] = detail::unscaled(t);
                if (t_hi > 1500)
                    return 0;
                return static_cast<T>(detail::rounded(X.per_scale(detail::exp_scaled(-t_hi, -t_lo))));
            }
            auto const [e_hi, e_lo] = detail::poisson_exponent(nu, t);      // λ e^E c / √ν
            if (e_hi < -1500)
                return 0;
            return static_cast<T>(detail::rounded(X.per_scale(detail::exp_scaled(e_hi, e_lo) * c / detail::sqrt(detail::normalized(nu, 0.0)))));
        }

        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            double const v = x, a = X.shape_;
            if (!(v >= 0))
                return v == v ? -std::numeric_limits<T>::infinity() : x;
            if (v == 0)
                return a < 1 ? std::numeric_limits<T>::infinity() : a == 1 ? static_cast<T>(X.log_rate().hi + X.log_rate().lo) : -std::numeric_limits<T>::infinity();
            if (v == std::numeric_limits<double>::infinity())
                return -std::numeric_limits<T>::infinity();
            detail::scaled const t = X.argument(v);
            fp::double_double const minus_half_log_2pi{-detail::half_log_2pi_hi, -detail::half_log_2pi_lo};
            if (a < 1)                                                      // E − log(2π)/2 + (log a)/2 − log x
            {
                auto e = detail::poisson_exponent(a, t);
                if (!detail::is_finite(e.hi))
                    return static_cast<T>(e.hi);
                auto const la = detail::log_precise(a);
                auto const lx = detail::log_precise(v);
                e = detail::sum(e, minus_half_log_2pi);
                e = detail::sum(e, {fp::fmul(la.hi, 0.5), fp::fmul(la.lo, 0.5)});
                e = detail::sum(e, {-lx.hi, -lx.lo});
                return static_cast<T>(e.hi + e.lo);
            }
            double const nu = a - 1;
            if (nu == 0)                                                    // log λ − λx
            {
                auto const [t_hi, t_lo] = detail::unscaled(t);
                if (!detail::is_finite(t_hi))
                    return -std::numeric_limits<T>::infinity();
                auto const e = detail::sum(X.log_rate(), {-t_hi, -t_lo});
                return static_cast<T>(e.hi + e.lo);
            }
            auto e = detail::poisson_exponent(nu, t);                       // E − log(2π)/2 − (log ν)/2 + log λ
            if (!detail::is_finite(e.hi))
                return static_cast<T>(e.hi);
            auto const ln = detail::log_precise(nu);
            e = detail::sum(e, minus_half_log_2pi);
            e = detail::sum(e, {fp::fmul(ln.hi, -0.5), fp::fmul(ln.lo, -0.5)});
            e = detail::sum(e, X.log_rate());
            return static_cast<T>(e.hi + e.lo);
        }

        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            double const a = shape_;
            auto const la = detail::log_precise(a);
            fp::double_double h{fp::fmul(la.hi, 0.5), fp::fmul(la.lo, 0.5)};
            h = detail::sum(h, {fp::half_log_two_pi_hi, fp::half_log_two_pi_lo});
            h = detail::sum(h, {0.5, 0.0});
            double const r = 0.5 / a;
            auto const [m, m_error] = fp::two_product(r, a);
            h = detail::sum(h, {-r, -(((0.5 - m) - m_error) / a)});         // −1/(2a)
            h = detail::sum(h, detail::stirling_error(a));
            auto const [c_hi, c_lo] = fp::two_sum(1.0, -a);
            h = detail::sum(h, detail::mul({c_hi, c_lo}, detail::stirling_error_derivative(a)));
            auto const lr = log_rate();
            h = detail::sum(h, {-lr.hi, -lr.lo});
            return static_cast<T>(h.hi + h.lo);
        }

        constexpr T cdf(T x) const noexcept { auto const& X = *this; return X.tail(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& X = *this; return X.tail(x, true); }

        constexpr T quantile(double p) const
        {
            auto const& X = *this;
            namespace fp = detail::sd;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return 0;
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            auto const [t_hi, t_lo] = p <= 0.5 ? fp::gamma_inverse_dd(double(X.shape_), p, false)
                                              : fp::gamma_inverse_dd(double(X.shape_), 1.0 - p, true);
            if (!detail::is_finite(t_hi))
                return std::numeric_limits<T>::infinity();
            if (t_hi < 0x1p-900)                                            // t sous-normal ou presque : perdu à
            {                                                               // l'échelle 2^-1074 ; or P(a, t) = tᵃ/Γ(a + 1)
                double const a = X.shape_;                                  // à 2^-900 près, d'où
                auto const lp = detail::log_precise(p);                     // log t = (log p + lgamma(1 + a))/a
                fp::double_double lg{};
                if (a <= 0.25)
                    lg = fp::lgamma1p_small(a);
                else
                {
                    auto const [s_hi, s_lo] = fp::two_sum(1.0, a);
                    lg = detail::log_gamma(s_hi, s_lo);
                }
                fp::double_double l = fp::quotient_dd(detail::sum(lp, lg), {a, 0.0});
                auto const lr = X.log_rate();
                l = detail::sum(l, {-lr.hi, -lr.lo});                       // x = t/λ = t θ
                if (l.hi < -1500)
                    return 0;
                return static_cast<T>(detail::rounded(detail::exp_scaled(l.hi, l.lo)));
            }
            return static_cast<T>(detail::rounded(X.per_rate(detail::normalized(t_hi, t_lo))));
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        constexpr T cquantile(double q) const
        {
            auto const& X = *this;
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::probability::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return wsl::probability::quantile(X, 1.0 - q);                                // 1 − q exact
            namespace fp = detail::sd;
            if (q == 0)
                return std::numeric_limits<T>::infinity();
            auto const [t_hi, t_lo] = fp::gamma_inverse_dd(double(X.shape_), q, true);
            if (!detail::is_finite(t_hi))
                return std::numeric_limits<T>::infinity();
            if (t_hi < 0x1p-900)                                            // (petite forme) comme pour quantile,
            {                                                               // log p = log1p(−q)
                double const a = X.shape_;
                auto const lp = detail::log1p_dd(-q);
                fp::double_double lg{};
                if (a <= 0.25)
                    lg = fp::lgamma1p_small(a);
                else
                {
                    auto const [s_hi, s_lo] = fp::two_sum(1.0, a);
                    lg = detail::log_gamma(s_hi, s_lo);
                }
                fp::double_double l = fp::quotient_dd(detail::sum(lp, lg), {a, 0.0});
                auto const lr = X.log_rate();
                l = detail::sum(l, {-lr.hi, -lr.lo});
                if (l.hi < -1500)
                    return 0;
                return static_cast<T>(detail::rounded(detail::exp_scaled(l.hi, l.lo)));
            }
            return static_cast<T>(detail::rounded(X.per_rate(detail::normalized(t_hi, t_lo))));
        }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            double y = 0;
            if (X.shape_ >= 1)
                y = detail::marsaglia_tsang(X.d_, X.c_, g);
            else                                                            // G(a + 1) U^{1/a}
            {
                double const h = detail::marsaglia_tsang(X.d_, X.c_, g);
                double const u = detail::uniform_positive(g);
                y = fmul_exp(h, wsl::random::detail::pm::log(u) / double(X.shape_));
            }
            if (X.rate_ != 0)
                return static_cast<T>(y / X.rate_);
            if (X.scale_ != 0)
                return static_cast<T>(detail::sd::fmul(y, X.scale_));
            return static_cast<T>(y > 0 ? detail::rounded(X.lambda_.per(detail::normalized(y, 0.0))) : 0.0);
        }

        // Identité de l'algèbre (voir sum.hpp) : gamma(a, λ) + gamma(b, λ) = gamma(a + b, λ), dans
        // la paramétrisation de X ; vide si les intensités diffèrent.
        friend constexpr std::optional<gamma> sum_closure(gamma const& X, gamma const& Y) noexcept
        {
            auto const a = static_cast<T>(double(X.shape_) + double(Y.shape_));
            if (!(X.lambda_ == Y.lambda_) || !detail::is_finite(a))
                return std::nullopt;
            return gamma(a, X.lambda_);
        }

        // Identité de l'algèbre (voir affine.hpp) : c X, c > 0, est gamma(a, λ/c), λ/c exact ;
        // vide si son dénominateur dépasse 2^53.
        friend constexpr std::optional<gamma> scale_closure(gamma const& X, double c) noexcept
        {
            auto const l = X.lambda_.divided(c);
            if (!l)
                return std::nullopt;
            return gamma(X.shape_, *l);
        }

        friend constexpr bool operator==(gamma const& a, gamma const& b) noexcept
        {
            return a.shape_ == b.shape_ && a.lambda_ == b.lambda_;
        }

    private:
        // Forme canonique : λ réduit ; voies rapides (calculs inchangés) quand λ ou 1/λ est un double.
        constexpr gamma(T shape, detail::exact_rate l) : shape_(shape), lambda_(l), rate_(l.as_double()), scale_(rate_ != 0 ? 0.0 : l.inverse().as_double())
        {
            validate();
        }

        constexpr void validate()
        {
            namespace fp = detail::sd;
            if (!(shape_ > 0) || !detail::is_finite(shape_))
                throw std::domain_error("wsl::probability::laws::gamma : il faut shape > 0, fini");
            double const a = shape_ >= 1 ? double(shape_) : double(shape_) + 1;
            d_ = a - 1.0 / 3;
            c_ = 1 / wsl::random::detail::pm::sqrt(fp::fmul(9.0, d_));
        }

        // P(X <= x) (ou P(X > x)) : P(a, t) de special_functions en t = t_hi, corrigé de
        // f(t_hi) t_lo (f la densité de la gamma standard), t = λx (ou x/θ) en double-double.
        constexpr T tail(T x, bool upper) const noexcept
        {
            namespace fp = detail::sd;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0))
                return upper ? T(1) : T(0);
            if (v == std::numeric_limits<double>::infinity())
                return upper ? T(0) : T(1);
            double t_hi = 0, t_lo = 0;
            if (rate_ != 0)
            {
                auto const [h, l] = detail::product(v, rate_);
                t_hi = h;
                t_lo = l;
            }
            else
            {
                auto const [h, l] = detail::unscaled(scale_ != 0 ? detail::normalized(v, 0.0) / scale_ : lambda_.times(detail::normalized(v, 0.0)));
                t_hi = h;
                t_lo = l;
            }
            double special = 0;
            if (fp::incomplete_gamma_special(double(shape_), t_hi, upper, special))
                return static_cast<T>(special);
            auto [value, value_upper] = fp::incomplete_gamma(double(shape_), t_hi);
            if (t_lo != 0 && value.hi != 0)
            {
                detail::scaled const density = detail::normalized(double(shape_), 0.0) * fp::gamma_prefix(double(shape_), t_hi) / t_hi;
                double const r = detail::rounded(density * detail::normalized(detail::magnitude(t_lo), 0.0) / value);
                double const signed_r = (t_lo > 0) != value_upper ? r : -r;  // dP/dt = f, dQ/dt = −f
                auto const [c_hi, c_lo] = fp::two_sum(1.0, signed_r);
                value = value * detail::normalized(c_hi, c_lo);
            }
            return static_cast<T>(value_upper == upper ? detail::rounded(value) : fp::one_minus(value));
        }

        static constexpr double fmul_exp(double h, double l) noexcept
        {
            return detail::sd::fmul(h, wsl::random::detail::pm::exp(l));
        }

        // λx (ou x/θ), à l'échelle.
        constexpr detail::scaled argument(double x) const noexcept
        {
            return per_scale(detail::normalized(x, 0.0));
        }
        // v · λ (ou v / θ), et v / λ (ou v · θ).
        constexpr detail::scaled per_scale(detail::scaled const& v) const noexcept
        {
            return rate_ != 0 ? v * detail::normalized(rate_, 0.0) : scale_ != 0 ? v / scale_ : lambda_.times(v);
        }
        constexpr detail::scaled per_rate(detail::scaled const& v) const noexcept
        {
            return rate_ != 0 ? v / rate_ : scale_ != 0 ? v * detail::normalized(scale_, 0.0) : lambda_.per(v);
        }
        constexpr detail::sd::double_double log_rate() const noexcept
        {
            if (rate_ == 0 && scale_ == 0)
                return lambda_.log();
            auto const [hi, lo] = detail::log_precise(rate_ != 0 ? rate_ : scale_);
            return rate_ != 0 ? detail::sd::double_double{hi, lo} : detail::sd::double_double{-hi, -lo};
        }

        T shape_;
        detail::exact_rate lambda_;
        double rate_ = 0, scale_ = 0;                                       // λ, 1/λ s'ils sont des doubles (sinon 0)
        double d_ = 0, c_ = 0;                                              // Marsaglia–Tsang
    };
}
