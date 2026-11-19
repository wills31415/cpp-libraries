#pragma once

// Loi de Weibull de forme k > 0 et d'échelle λ > 0 : P(X > x) = e^{−(x/λ)^k}, x >= 0.
//
//   wsl::weibull W({.shape = 1.5, .scale = 2});
//   cdf(W, 1.0); quantile(W, 0.5); W.mean(); sample(W, g);
//
// t = (x/λ)^k = e^{k (log x − log λ)} en double-double ; répartition −expm1(−t),
// survie e^{−t}, densité k t e^{−t}/x ; quantile λ (−log1p(−p))^{1/k}, un arrondi.
// Moments : Γ(1 + 1/k) par lgamma en double-double ; la variance
// λ² [Γ(1 + 2/k) − Γ(1 + 1/k)²], qui se compense pour k grand, par
// Γ(1 + ε)² (e^D − 1), D = lgamma(1 + 2ε) − 2 lgamma(1 + ε) en série (les termes
// linéaires s'annulent exactement).
//
// Tirage (approché, bit-exact) : λ E^{1/k}, E exponentielle tirée exactement.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/gamma_moments.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/exponential.hpp>
#include <wsl/random/law.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl
{
    template <law_real T = double>
    class weibull
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T shape, scale;
        };

        constexpr weibull(parameters p) : k_(p.shape), lambda_(p.scale)
        {
            if (!(k_ > 0) || !detail::is_finite(k_) || !(lambda_ > 0) || !detail::is_finite(lambda_))
                throw std::domain_error("wsl::weibull : il faut shape > 0 et scale > 0, finis");
        }

        constexpr T shape() const noexcept { return k_; }
        constexpr T scale() const noexcept { return lambda_; }
        constexpr interval<T> support() const noexcept { return {T(0), std::numeric_limits<T>::infinity()}; }

        // λ Γ(1 + 1/k)
        constexpr T mean() const noexcept
        {
            auto const lg = lgamma_1p(1.0 / double(k_), inv_k_lo());
            return static_cast<T>(exp_of(detail::sum(lg, detail::log_precise(double(lambda_)))));
        }
        // λ² [Γ(1 + 2/k) − Γ(1 + 1/k)²]
        constexpr T variance() const noexcept
        {
            return static_cast<T>(detail::gamma_variance(1.0 / double(k_), inv_k_lo(), detail::log_precise(double(lambda_))));
        }
        // λ (log 2)^{1/k}
        constexpr T median() const noexcept
        {
            namespace fp = detail::sd;
            auto const l2 = detail::log(detail::normalized(fp::ln2_hi, fp::ln2_lo));   // log log 2
            return static_cast<T>(exp_of(detail::sum(fp::quotient_dd(l2, {double(k_), 0.0}), detail::log_precise(double(lambda_)))));
        }
        // λ ((k − 1)/k)^{1/k} si k > 1, 0 sinon
        constexpr T mode() const noexcept
        {
            namespace fp = detail::sd;
            double const k = k_;
            if (k <= 1)
                return 0;
            double const r = 1.0 / k;
            auto const l = detail::log1p_dd(-r);                            // log((k − 1)/k), au premier ordre en 1/k
            return static_cast<T>(exp_of(detail::sum(fp::quotient_dd({l.hi, l.lo - inv_k_lo() / (1 - r)}, {k, 0.0}), detail::log_precise(double(lambda_)))));
        }
        constexpr T skewness() const noexcept
        {
            double const k = k_;
            namespace fp = detail::sd;                                       // (approché : la compensation n'est pas traitée)
            double const g1 = wsl::special::tgamma(1 + 1 / k), g2 = wsl::special::tgamma(1 + 2 / k), g3 = wsl::special::tgamma(1 + 3 / k);
            double const v = g2 - fp::fmul(g1, g1);
            return static_cast<T>(((g3 - fp::fmul(fp::fmul(3.0, g1), v)) - fp::fmul(fp::fmul(g1, g1), g1)) / (v * wsl::special::sqrt(v)));
        }
        // γ (1 − 1/k) + log(λ/k) + 1
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            double const k = k_;
            auto const [c_hi, c_lo] = detail::sum({1.0, 0.0}, {-1.0 / k, -inv_k_lo()});
            fp::double_double h = detail::mul({fp::euler_gamma_hi, fp::euler_gamma_lo}, {c_hi, c_lo});
            h = detail::sum(h, detail::log_precise(double(lambda_)));
            auto const lk = detail::log_precise(k);
            h = detail::sum(h, {-lk.hi, -lk.lo});
            h = detail::sum(h, {1.0, 0.0});
            return static_cast<T>(h.hi + h.lo);
        }

        friend constexpr T pdf(weibull const& X, T x) noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v >= 0) || v == std::numeric_limits<double>::infinity())
                return 0;
            if (v == 0)
                return double(X.k_) < 1 ? std::numeric_limits<T>::infinity() : double(X.k_) == 1 ? static_cast<T>(1.0 / double(X.lambda_)) : T(0);
            auto const l = X.log_density(v);
            if (l.hi < -1500)
                return 0;
            if (l.hi > 1400)
                return std::numeric_limits<T>::infinity();
            return static_cast<T>(detail::rounded(detail::exp_scaled(l.hi, l.lo)));
        }
        friend constexpr T log_pdf(weibull const& X, T x) noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0) || v == std::numeric_limits<double>::infinity())
                return v == 0 && double(X.k_) <= 1 ? static_cast<T>(double(X.k_) < 1 ? std::numeric_limits<double>::infinity() : -wsl::special::log(double(X.lambda_)))
                                                   : -std::numeric_limits<T>::infinity();
            auto const l = X.log_density(v);
            return static_cast<T>(l.hi + l.lo);
        }

        friend constexpr T cdf(weibull const& X, T x) noexcept { return X.distribution(x, false); }
        friend constexpr T ccdf(weibull const& X, T x) noexcept { return X.distribution(x, true); }

        // λ (−log1p(−p))^{1/k}
        friend constexpr T quantile(weibull const& X, double p)
        {
            namespace fp = detail::sd;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return 0;
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            auto const [l_hi, l_lo] = detail::log1p_dd(-p);                  // < 0
            auto const ll = detail::log(detail::normalized(-l_hi, -l_lo));   // log(−log(1 − p))
            return static_cast<T>(exp_of(detail::sum(fp::quotient_dd(ll, {double(X.k_), 0.0}), detail::log_precise(double(X.lambda_)))));
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        friend constexpr T cquantile(weibull const& X, double q)
        {
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return quantile(X, 1.0 - q);                                // 1 − q exact
            namespace fp = detail::sd;
            if (q == 0)
                return std::numeric_limits<T>::infinity();
            auto const [l_hi, l_lo] = detail::log_precise(q);               // λ (−log q)^{1/k}
            auto const ll = detail::log(detail::normalized(-l_hi, -l_lo));
            return static_cast<T>(exp_of(detail::sum(fp::quotient_dd(ll, {double(X.k_), 0.0}), detail::log_precise(double(X.lambda_)))));
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(weibull const& X, G& g)
        {
            double const e = sample(exponential<double>({.rate = 1.0}), g);
            if (e == 0)
                return 0;
            return static_cast<T>(detail::sd::fmul(double(X.lambda_), wsl::special::exp(wsl::special::log(e) / double(X.k_))));
        }

        friend constexpr bool operator==(weibull const&, weibull const&) = default;

    private:
        // 1/k − fl(1/k).
        constexpr double inv_k_lo() const noexcept
        {
            double const k = k_, r = 1.0 / k;
            auto const [m, m_error] = detail::sd::two_product(r, k);
            return ((1.0 - m) - m_error) / k;
        }

        // lgamma(1 + ε), ε = e + e_lo, en double-double.
        static constexpr detail::sd::double_double lgamma_1p(double e, double e_lo) noexcept { return detail::lgamma_one_plus(e, e_lo); }

        static constexpr double exp_of(detail::sd::double_double l) noexcept
        {
            if (l.hi > 1400)
                return std::numeric_limits<double>::infinity();
            if (l.hi < -1500)
                return 0;
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }

        // y = k log(x/λ) = log t, en double-double : log du rapport r = x/λ (en double-double),
        // par log1p(r − 1) près de 1 — précis relativement, car k l'amplifie.
        constexpr detail::sd::double_double log_t(double x) const noexcept
        {
            namespace fp = detail::sd;
            detail::scaled const r = detail::normalized(x, 0.0) / double(lambda_);
            fp::double_double l{};
            if (r.e == 0 || (r.e == -1 && r.hi > 1.0))                     // r dans (0,5 ; 2)
            {
                auto const [r_hi, r_lo] = detail::unscaled(r);
                auto const [m_hi, m_lo] = fp::fast_two_sum(r_hi - 1.0, r_lo);    // r − 1 : r_hi − 1 exact
                if (m_hi >= 0)
                    l = detail::log1p_positive(m_hi, m_lo);
                else
                {
                    auto const [h, lo] = detail::log1p_dd(m_hi);
                    l = {h, lo + m_lo / (1.0 + m_hi)};
                }
            }
            else
                l = detail::log(r);
            return detail::mul({double(k_), 0.0}, l);
        }

        // log f = log k + y − t − log x.
        constexpr detail::sd::double_double log_density(double x) const noexcept
        {
            auto const y = log_t(x);
            if (y.hi > 710)
                return {-std::numeric_limits<double>::infinity(), 0.0};
            auto const t = detail::unscaled(detail::exp_scaled(y.hi, y.lo));
            auto l = detail::sum(detail::log_precise(double(k_)), y);
            l = detail::sum(l, {-t.hi, -t.lo});
            auto const lx = detail::log_precise(x);
            return detail::sum(l, {-lx.hi, -lx.lo});
        }

        constexpr T distribution(T x, bool upper) const noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0))
                return upper ? T(1) : T(0);
            if (v == std::numeric_limits<double>::infinity())
                return upper ? T(0) : T(1);
            auto const y = log_t(v);
            if (y.hi > 7)                                                   // t > 1096 : e^{−t} < 2^-1500
                return upper ? T(0) : T(1);
            if (y.hi < -1500)                                               // t sous 2^-2000
                return upper ? T(1) : T(0);
            detail::scaled const ts = detail::exp_scaled(y.hi, y.lo);
            if (ts.e < -900)                                                // 1 − e^{−t} = t (1 − t/2 …)
                return upper ? T(1) : static_cast<T>(detail::rounded(ts));
            auto const t = detail::unscaled(ts);
            if (upper)
                return static_cast<T>(detail::rounded(detail::exp_scaled(-t.hi, -t.lo)));
            if (t.hi > 40)
                return 1;
            return static_cast<T>(-detail::expm1_dd(-t.hi, -t.lo));
        }

        T k_, lambda_;
    };
}
