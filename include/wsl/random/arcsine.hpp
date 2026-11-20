#pragma once

// Loi de l'arcsinus sur [a, b] : densité 1/(π √((x − a)(b − x))).
//
//   wsl::arcsine A({.low = 0, .high = 1});   cdf(A, 0.1); quantile(A, 0.9); sample(A, g);
//
// u = (x − a)/(b − a) et 1 − u = (b − x)/(b − a) en double-double ; répartition
// (2/π) atan √(u/(1 − u)) (u <= ½), 1 − (2/π) atan √((1 − u)/u) sinon, noyau atan à argument
// double-double ; quantile a + (b − a) sin²(πp/2) (p <= ½), b − (b − a) sin²(π(1 − p)/2) ;
// un arrondi. Moyenne et médiane (a + b)/2, variance (b − a)²/8, entropie log(π(b − a)/4).
//
// Tirage (approché, bit-exact) : quantile d'un uniforme.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/trigonometry.hpp>
#include <wsl/random/detail/discrete.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl
{
    template <law_real T = double>
    class arcsine
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T low = 0, high = 1;
        };

        constexpr arcsine(parameters p) : a_(p.low), b_(p.high)
        {
            if (!detail::is_finite(a_) || !detail::is_finite(b_) || !(a_ < b_))
                throw std::domain_error("wsl::arcsine : il faut low < high, finis");
        }

        constexpr T low() const noexcept { return a_; }
        constexpr T high() const noexcept { return b_; }
        constexpr interval<T> support() const noexcept { return {a_, b_}; }
        constexpr T mean() const noexcept { return static_cast<T>(detail::sd::fmul(double(a_), 0.5) + detail::sd::fmul(double(b_), 0.5)); }
        constexpr T median() const noexcept { return mean(); }
        // (b − a)²/8
        constexpr T variance() const noexcept
        {
            auto const w = width();
            detail::scaled const s = detail::normalized(w.hi, w.lo);
            detail::scaled v = s * s;
            v.e -= 3;
            return static_cast<T>(detail::rounded(v));
        }
        constexpr T skewness() const noexcept { return 0; }
        constexpr T excess_kurtosis() const noexcept { return static_cast<T>(-1.5); }
        // log π + log(b − a) − 2 log 2
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            auto const w = width();
            auto h = detail::sum(detail::log(detail::normalized(w.hi, w.lo)), {fp::log_pi_hi, fp::log_pi_lo});
            h = detail::sum(h, {fp::fmul(fp::ln2_hi, -2.0), fp::fmul(fp::ln2_lo, -2.0)});
            return static_cast<T>(h.hi + h.lo);
        }

        // 1/(π (b − a) √(u (1 − u)))
        friend constexpr T pdf(arcsine const& A, T x) noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > double(A.a_) && v < double(A.b_)))
                return v == double(A.a_) || v == double(A.b_) ? std::numeric_limits<T>::infinity() : T(0);
            auto const [u, w] = A.parts(v);
            auto const ba = A.width();
            detail::scaled const d = detail::sqrt(u * w) * detail::normalized(ba.hi, ba.lo);
            return static_cast<T>(detail::rounded(detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo) / d));
        }
        friend constexpr T log_pdf(arcsine const& A, T x) noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > double(A.a_) && v < double(A.b_)))
                return v == double(A.a_) || v == double(A.b_) ? std::numeric_limits<T>::infinity() : -std::numeric_limits<T>::infinity();
            auto const [u, w] = A.parts(v);
            auto const ba = A.width();
            auto l = detail::log(u * w);
            l = {detail::sd::fmul(l.hi, -0.5), detail::sd::fmul(l.lo, -0.5)};
            auto const lb = detail::log(detail::normalized(ba.hi, ba.lo));
            l = detail::sum(l, {-lb.hi, -lb.lo});
            l = detail::sum(l, {-detail::sd::log_pi_hi, -detail::sd::log_pi_lo});
            return static_cast<T>(l.hi + l.lo);
        }
        friend constexpr T cdf(arcsine const& A, T x) noexcept { return A.distribution(x, false); }
        friend constexpr T ccdf(arcsine const& A, T x) noexcept { return A.distribution(x, true); }
        friend constexpr T quantile(arcsine const& A, double p) { return A.inverse(p, false); }
        friend constexpr T cquantile(arcsine const& A, double q) { return A.inverse(q, true); }

        template <std::uniform_random_bit_generator G>
        friend constexpr T sample(arcsine const& A, G& g)
        {
            return A.inverse(uniform_real<double>(g), false);
        }

        friend constexpr bool operator==(arcsine const&, arcsine const&) = default;

    private:
        constexpr detail::sd::double_double width() const noexcept
        {
            auto const [h, l] = detail::sd::two_sum(double(b_), -double(a_));
            return {h, l};
        }

        // u = (x − a)/(b − a), 1 − u = (b − x)/(b − a), à l'échelle.
        struct uw
        {
            detail::scaled u, w;
        };
        constexpr uw parts(double x) const noexcept
        {
            auto const ba = width();
            auto const [p_hi, p_lo] = detail::sd::two_sum(x, -double(a_));
            auto const [q_hi, q_lo] = detail::sd::two_sum(double(b_), -x);
            detail::scaled const d = detail::normalized(ba.hi, ba.lo);
            return {detail::normalized(p_hi, p_lo) / d, detail::normalized(q_hi, q_lo) / d};
        }

        // (2/π) atan √(s/l), s <= l (le petit côté), à l'échelle.
        static constexpr detail::scaled side(detail::scaled const& s, detail::scaled const& l) noexcept
        {
            namespace fp = detail::sd;
            auto const [v_hi, v_lo] = detail::unscaled(detail::sqrt(s / l));
            fp::double_double const at = v_hi < 0x1p-27 ? fp::double_double{v_hi, v_lo - fp::fmul(fp::fmul(v_hi, v_hi), v_hi) / 3} : fp::atan_kernel(v_hi, v_lo);
            detail::scaled r = detail::normalized(at.hi, at.lo) * detail::normalized(detail::inv_pi_hi, detail::inv_pi_lo);
            r.e += 1;
            return r;
        }

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
            auto const [u, w] = parts(v);
            bool const low_side = detail::compare(u, detail::unscaled(w)) <= 0;   // u <= 1 − u
            detail::scaled const part = low_side ? side(u, w) : side(w, u);       // F si low_side, sinon 1 − F
            bool const direct = low_side != upper;
            return static_cast<T>(direct ? detail::rounded(part) : fp::one_minus(part));
        }

        constexpr T inverse(double t, bool upper) const
        {
            namespace fp = detail::sd;
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::quantile : probabilité hors de [0, 1]");
            if (t == 0)
                return upper ? b_ : a_;
            if (t == 1)
                return upper ? a_ : b_;
            // part r <= ½ mesurée depuis la borne la plus proche : x = a + (b − a) sin²(π r/2) ou b − …
            bool const from_low = upper ? t > 0.5 : t <= 0.5;
            double const r = (upper ? t > 0.5 : t > 0.5) ? 1.0 - t : t;     // exact (Sterbenz)
            auto const h = detail::mul({fp::fmul(detail::pi_hi, 0.5), fp::fmul(detail::pi_lo, 0.5)}, {r, 0.0});
            auto const s = detail::sine_dd(h.hi, h.lo);
            if (s.hi == 0)
                return from_low ? a_ : b_;
            detail::scaled const sn = detail::normalized(s.hi, s.lo);
            auto const ba = width();
            detail::scaled const d = sn * sn * detail::normalized(ba.hi, ba.lo);
            return static_cast<T>(detail::affine_scaled(from_low ? double(a_) : double(b_), !from_low, d));
        }

        T a_, b_;
    };
}
