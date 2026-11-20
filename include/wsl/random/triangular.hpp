#pragma once

// Loi triangulaire sur [a, b] de mode c (a <= c <= b, a < b).
//
//   wsl::random::triangular T({.lower = 0, .mode = 1, .upper = 3});
//   T.cdf(2.0); T.quantile(0.9); T.variance(); T(g);
//
// Répartition (x − a)²/((b − a)(c − a)) sous c, 1 − (b − x)²/((b − a)(b − c)) au-dessus
// (survie symétrique), différences exactes, un arrondi ; quantile a + √(p (b − a)(c − a))
// sous F(c) = (c − a)/(b − a), b − √((1 − p)(b − a)(b − c)) au-dessus (1 − p en
// double-double). Variance ((a − b)² + (a − c)² + (b − c)²)/36, sans compensation.
//
// Tirage (approché, bit-exact) : inversion d'un uniforme.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    template <law_real T = double>
    class triangular : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T lower, mode, upper;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.lower = lower(), .mode = mode(), .upper = upper()}; }

        constexpr triangular(parameters p) : a_(p.lower), c_(p.mode), b_(p.upper)
        {
            if (!detail::is_finite(a_) || !detail::is_finite(b_) || !(a_ < b_) || !(a_ <= c_ && c_ <= b_))
                throw std::domain_error("wsl::random::triangular : il faut lower <= mode <= upper, lower < upper, finis");
        }

        constexpr T lower() const noexcept { return a_; }
        constexpr T mode() const noexcept { return c_; }
        constexpr T upper() const noexcept { return b_; }
        constexpr interval<T> support() const noexcept { return {a_, b_}; }

        // (a + b + c)/3
        constexpr T mean() const noexcept
        {
            auto const s = detail::sum(detail::sum({double(a_), 0.0}, {double(b_), 0.0}), {double(c_), 0.0});
            auto const q = detail::sd::quotient_dd(s, {3.0, 0.0});
            return static_cast<T>(q.hi + q.lo);
        }
        constexpr T variance() const noexcept
        {
            auto const ab = diff(a_, b_), ac = diff(a_, c_), bc = diff(b_, c_);
            auto s = detail::sum(detail::mul(ab, ab), detail::mul(ac, ac));
            s = detail::sum(s, detail::mul(bc, bc));
            auto const q = detail::sd::quotient_dd(s, {36.0, 0.0});
            return static_cast<T>(q.hi + q.lo);
        }
        // a + √((b − a)(c − a)/2) si c >= (a + b)/2, sinon b − √((b − a)(b − c)/2)
        constexpr T median() const noexcept
        {
            auto const ba = diff(b_, a_);
            auto const cm = detail::sum(diff(c_, a_), {-detail::sd::fmul(ba.hi, 0.5), -detail::sd::fmul(ba.lo, 0.5)});
            bool const right = cm.hi >= 0;
            auto const other = right ? diff(c_, a_) : diff(b_, c_);
            detail::scaled w = detail::normalized(ba.hi, ba.lo) * detail::normalized(other.hi, other.lo);
            w.e -= 1;
            return static_cast<T>(detail::affine_scaled(right ? double(a_) : double(b_), !right, detail::sqrt(w)));
        }
        // ½ + log((b − a)/2)
        constexpr T entropy() const noexcept
        {
            auto const ba = diff(b_, a_);
            auto h = detail::log(detail::normalized(ba.hi, ba.lo));
            h = detail::sum(h, {0.5 - detail::sd::ln2_hi, -detail::sd::ln2_lo});
            return static_cast<T>(h.hi + h.lo);
        }
        // √2 (a + b − 2c)(2a − b − c)(a − 2b + c) / (5 Q^{3/2}), Q = a² + b² + c² − ab − ac − bc
        constexpr T skewness() const noexcept
        {
            namespace fp = detail::sd;
            double const a = a_, b = b_, c = c_;
            double const q = fp::fmul((fp::fmul(a - b, a - b) + fp::fmul(a - c, a - c)) + fp::fmul(b - c, b - c), 0.5);
            double const n = fp::fmul(fp::fmul((a + b) - fp::fmul(2.0, c), (fp::fmul(2.0, a) - b) - c), (a - fp::fmul(2.0, b)) + c);
            return static_cast<T>(fp::fmul(detail::sqrt2_hi, n) / fp::fmul(5.0, fp::fmul(q, wsl::special::sqrt(q))));
        }
        constexpr T excess_kurtosis() const noexcept { return static_cast<T>(-0.6); }

        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v >= double(X.a_) && v <= double(X.b_)))
                return 0;
            auto const ba = diff(X.b_, X.a_);
            bool const left = v < double(X.c_) || (v == double(X.c_) && X.c_ == X.b_);
            auto const num = left ? diff(T(v), X.a_) : diff(X.b_, T(v));
            auto const den = left ? diff(X.c_, X.a_) : diff(X.b_, X.c_);
            if (num.hi == 0)
                return 0;
            detail::scaled const r = detail::normalized(num.hi, num.lo) / (detail::normalized(ba.hi, ba.lo) * detail::normalized(den.hi, den.lo));
            return static_cast<T>(detail::rounded(r * detail::normalized(2.0, 0.0)));
        }
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            T const f = wsl::random::pdf(X, x);
            if (!(f > 0))
                return f == f ? -std::numeric_limits<T>::infinity() : f;
            double const v = x;
            auto const ba = diff(X.b_, X.a_);
            bool const left = v < double(X.c_) || (v == double(X.c_) && X.c_ == X.b_);
            auto const num = left ? diff(T(v), X.a_) : diff(X.b_, T(v));
            auto const den = left ? diff(X.c_, X.a_) : diff(X.b_, X.c_);
            detail::scaled const r = detail::normalized(num.hi, num.lo) / (detail::normalized(ba.hi, ba.lo) * detail::normalized(den.hi, den.lo));
            auto l = detail::log(r);
            l = detail::sum(l, {detail::sd::ln2_hi, detail::sd::ln2_lo});
            return static_cast<T>(l.hi + l.lo);
        }
        constexpr T cdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, true); }
        constexpr T quantile(double p) const { auto const& X = *this; return X.inverse(p, false); }
        constexpr T cquantile(double q) const { auto const& X = *this; return X.inverse(q, true); }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            return X.inverse(uniform_real<double>(g), false);
        }

        friend constexpr bool operator==(triangular const&, triangular const&) = default;

    private:
        static constexpr detail::sd::double_double diff(T x, T y) noexcept
        {
            auto const [h, l] = detail::sd::two_sum(double(x), -double(y));
            return {h, l};
        }

        // P(X <= x) (ou P(X > x)) : (x − a)²/((b − a)(c − a)) d'un côté, son complément de l'autre.
        constexpr T distribution(T x, bool upper) const noexcept
        {
            namespace fp = detail::sd;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > double(a_)))
                return upper ? T(1) : T(0);
            if (!(v < double(b_)))
                return upper ? T(0) : T(1);
            auto const ba = diff(b_, a_);
            bool const left = v <= double(c_);
            auto const num = left ? diff(x, a_) : diff(b_, x);
            auto const den = left ? diff(c_, a_) : diff(b_, c_);
            detail::scaled const n = detail::normalized(num.hi, num.lo);
            detail::scaled const part = n * n / (detail::normalized(ba.hi, ba.lo) * detail::normalized(den.hi, den.lo));
            bool const direct = left != upper;                              // la part calculée est la valeur demandée
            return static_cast<T>(direct ? detail::rounded(part) : fp::one_minus(part));
        }

        constexpr T inverse(double t, bool upper) const
        {
            namespace fp = detail::sd;
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::random::quantile : probabilité hors de [0, 1]");
            if (t == 0)
                return upper ? b_ : a_;
            if (t == 1)
                return upper ? a_ : b_;
            // p = P(X <= x) et 1 − p, en double-double
            auto const [o_hi, o_lo] = fp::two_sum(1.0, -t);
            fp::double_double const p = upper ? fp::double_double{o_hi, o_lo} : fp::double_double{t, 0.0};
            fp::double_double const q = upper ? fp::double_double{t, 0.0} : fp::double_double{o_hi, o_lo};
            auto const ba = diff(b_, a_), ca = diff(c_, a_), bc = diff(b_, c_);
            // p (b − a) <= c − a ?
            auto const lhs = detail::mul(p, ba);
            auto const gap = detail::sum(lhs, {-ca.hi, -ca.lo});
            bool const left = gap.hi < 0 || (gap.hi == 0 && gap.lo <= 0);
            if (left)
            {
                if (ca.hi == 0)
                    return a_;
                detail::scaled const w = detail::normalized(p.hi, p.lo) * detail::normalized(ba.hi, ba.lo) * detail::normalized(ca.hi, ca.lo);
                return static_cast<T>(detail::affine_scaled(double(a_), false, detail::sqrt(w)));
            }
            if (bc.hi == 0 || q.hi == 0)
                return b_;
            detail::scaled const w = detail::normalized(q.hi, q.lo) * detail::normalized(ba.hi, ba.lo) * detail::normalized(bc.hi, bc.lo);
            return static_cast<T>(detail::affine_scaled(double(b_), true, detail::sqrt(w)));
        }

        T a_, c_, b_;
    };
}
