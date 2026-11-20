#pragma once

// Loi de Pareto (type I) d'échelle x_m > 0 et de forme α > 0 : P(X > x) = (x_m/x)^α,
// x >= x_m.
//
//   wsl::random::pareto P({.scale = 1, .shape = 2.5});
//   P.mean(); P.cdf(3.0); P.quantile(0.99); P(g);
//
// L = α log(x/x_m) en double-double (log du rapport, par log1p près de 1) ; survie
// e^{−L}, répartition −expm1(−L), densité (α/x) e^{−L} ; quantile x_m e^{−log1p(−p)/α} ;
// moments au sens usuel (moyenne pour α > 1, variance pour α > 2, sinon +∞).
//
// Tirage (approché, bit-exact) : x_m e^{E/α}, E exponentielle tirée exactement.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/exponential.hpp>
#include <wsl/random/law.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    template <law_real T = double>
    class pareto : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T scale = 1, shape;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.scale = scale(), .shape = shape()}; }

        constexpr pareto(parameters p) : xm_(p.scale), alpha_(p.shape)
        {
            if (!(xm_ > 0) || !detail::is_finite(xm_) || !(alpha_ > 0) || !detail::is_finite(alpha_))
                throw std::domain_error("wsl::random::pareto : il faut scale > 0 et shape > 0, finis");
        }

        constexpr T scale() const noexcept { return xm_; }
        constexpr T shape() const noexcept { return alpha_; }
        constexpr interval<T> support() const noexcept { return {xm_, std::numeric_limits<T>::infinity()}; }

        // α x_m/(α − 1), α > 1
        constexpr T mean() const noexcept
        {
            double const a = alpha_;
            if (a <= 1)
                return std::numeric_limits<T>::infinity();
            auto const [m_hi, m_lo] = detail::sd::two_sum(a, -1.0);
            return static_cast<T>(detail::rounded(detail::normalized(a, 0.0) * detail::normalized(double(xm_), 0.0) / detail::normalized(m_hi, m_lo)));
        }
        // x_m² α / ((α − 1)² (α − 2)), α > 2
        constexpr T variance() const noexcept
        {
            double const a = alpha_;
            if (a <= 2)
                return std::numeric_limits<T>::infinity();
            auto const [m1_hi, m1_lo] = detail::sd::two_sum(a, -1.0);
            auto const [m2_hi, m2_lo] = detail::sd::two_sum(a, -2.0);
            detail::scaled const x = detail::normalized(double(xm_), 0.0), m1 = detail::normalized(m1_hi, m1_lo);
            return static_cast<T>(detail::rounded(x * x * detail::normalized(a, 0.0) / (m1 * m1 * detail::normalized(m2_hi, m2_lo))));
        }
        // x_m 2^{1/α}
        constexpr T median() const noexcept
        {
            auto const l = detail::sd::quotient_dd({detail::sd::ln2_hi, detail::sd::ln2_lo}, {double(alpha_), 0.0});
            return static_cast<T>(scaled_exp(l));
        }
        constexpr T mode() const noexcept { return xm_; }
        constexpr T skewness() const noexcept
        {
            double const a = alpha_;
            if (a <= 3)
                return std::numeric_limits<T>::quiet_NaN();
            return static_cast<T>(detail::sd::fmul(2 * (1 + a) / (a - 3), wsl::special::sqrt((a - 2) / a)));
        }
        constexpr T excess_kurtosis() const noexcept
        {
            namespace fp = detail::sd;
            double const a = alpha_;
            if (a <= 4)
                return std::numeric_limits<T>::quiet_NaN();
            double const num = ((fp::fmul(fp::fmul(a, a), a) + fp::fmul(a, a)) - fp::fmul(6.0, a)) - 2.0;
            return static_cast<T>(6 * num / (a * (a - 3) * (a - 4)));
        }
        // log(x_m/α) + 1 + 1/α
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            auto h = detail::sum(detail::log_precise(double(xm_)), {1.0, 0.0});
            auto const la = detail::log_precise(double(alpha_));
            h = detail::sum(h, {-la.hi, -la.lo});
            h = detail::sum(h, fp::quotient_dd({1.0, 0.0}, {double(alpha_), 0.0}));
            return static_cast<T>(h.hi + h.lo);
        }

        // (α/x) (x_m/x)^α
        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v >= double(X.xm_)) || v == std::numeric_limits<double>::infinity())
                return 0;
            auto const l = X.log_ratio(v);
            if (l.hi > 1500)
                return 0;
            return static_cast<T>(detail::rounded(detail::exp_scaled(-l.hi, -l.lo) * detail::normalized(double(X.alpha_), 0.0) / v));
        }
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v >= double(X.xm_)) || v == std::numeric_limits<double>::infinity())
                return -std::numeric_limits<T>::infinity();
            auto const r = X.log_ratio(v);
            auto l = detail::sum(detail::log_precise(double(X.alpha_)), {-r.hi, -r.lo});
            auto const lx = detail::log_precise(v);
            l = detail::sum(l, {-lx.hi, -lx.lo});
            return static_cast<T>(l.hi + l.lo);
        }

        constexpr T cdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, true); }

        // x_m e^{−log1p(−p)/α}
        constexpr T quantile(double p) const
        {
            auto const& X = *this;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::random::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return X.xm_;
            if (p == 1)
                return std::numeric_limits<T>::infinity();
            auto const [l_hi, l_lo] = detail::log1p_dd(-p);
            return static_cast<T>(X.scaled_exp(detail::sd::quotient_dd({-l_hi, -l_lo}, {double(X.alpha_), 0.0})));
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
            auto const [l_hi, l_lo] = detail::log_precise(q);               // x_m e^{−log q/α}
            return static_cast<T>(X.scaled_exp(detail::sd::quotient_dd({-l_hi, -l_lo}, {double(X.alpha_), 0.0})));
        }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            double const e = wsl::random::sample(exponential<double>({.rate = 1.0}), g);
            return static_cast<T>(detail::sd::fmul(double(X.xm_), wsl::special::exp(e / double(X.alpha_))));
        }

        friend constexpr bool operator==(pareto const&, pareto const&) = default;

    private:
        // x_m e^{l}
        constexpr double scaled_exp(detail::sd::double_double l) const noexcept
        {
            auto const e = detail::sum(l, detail::log_precise(double(xm_)));
            if (e.hi > 1400)
                return std::numeric_limits<double>::infinity();
            return detail::rounded(detail::exp_scaled(e.hi, e.lo));
        }

        // α log(x/x_m), x >= x_m, en double-double (log1p du rapport près de 1).
        constexpr detail::sd::double_double log_ratio(double x) const noexcept
        {
            namespace fp = detail::sd;
            detail::scaled const r = detail::normalized(x, 0.0) / double(xm_);
            fp::double_double l{};
            if (r.e == 0)                                                   // r dans [1, 2)
            {
                auto const [r_hi, r_lo] = detail::unscaled(r);
                auto const [m_hi, m_lo] = fp::fast_two_sum(r_hi - 1.0, r_lo);
                l = m_hi >= 0 ? detail::log1p_positive(m_hi, m_lo) : fp::double_double{0.0, 0.0};
            }
            else
                l = detail::log(r);
            return detail::mul({double(alpha_), 0.0}, l);
        }

        constexpr T distribution(T x, bool upper) const noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > double(xm_)))
                return upper ? T(1) : T(0);
            if (v == std::numeric_limits<double>::infinity())
                return upper ? T(0) : T(1);
            auto const l = log_ratio(v);
            if (upper)
                return l.hi > 1500 ? T(0) : static_cast<T>(detail::rounded(detail::exp_scaled(-l.hi, -l.lo)));
            if (l.hi > 40)
                return 1;
            if (l.hi < 0x1p-900)
                return static_cast<T>(l.hi);
            return static_cast<T>(-detail::expm1_dd(-l.hi, -l.lo));
        }

        T xm_, alpha_;
    };
}
