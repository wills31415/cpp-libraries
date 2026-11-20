#pragma once

// Loi de Rayleigh de paramètre σ > 0 : densité (x/σ²) e^{−x²/(2σ²)}, x >= 0.
//
//   wsl::random::rayleigh R({.scale = 2});   R.cdf(1.0); R.quantile(0.5); R.mean(); R(g);
//
// h = z²/2, z = x/σ en double-double ; répartition −expm1(−h), survie e^{−h}, densité
// (z/σ) e^{−h} ; quantile σ √(−2 log1p(−p)) (cquantile : σ √(−2 log q)) ; un arrondi.
// Moyenne σ √(π/2), variance σ² (4 − π)/2, médiane σ √(2 log 2), entropie
// 1 + log(σ/√2) + γ/2.
//
// Tirage (approché, bit-exact) : σ √(2E), E exponentielle exacte.

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
    class rayleigh : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T scale = 1;
        };

        constexpr rayleigh() : rayleigh(parameters{}) {}

        constexpr rayleigh(parameters p) : sigma_(p.scale)
        {
            if (!(sigma_ > 0) || !detail::is_finite(sigma_))
                throw std::domain_error("wsl::random::rayleigh : il faut scale > 0, fini");
        }

        constexpr T scale() const noexcept { return sigma_; }
        constexpr interval<T> support() const noexcept { return {T(0), std::numeric_limits<T>::infinity()}; }
        constexpr T mean() const noexcept { return times(detail::sqrt_half_pi_hi, detail::sqrt_half_pi_lo); }
        constexpr T median() const noexcept { return times(detail::sqrt_two_ln2_hi, detail::sqrt_two_ln2_lo); }
        constexpr T mode() const noexcept { return sigma_; }
        constexpr T variance() const noexcept
        {
            detail::scaled const s = detail::normalized(double(sigma_), 0.0);
            return static_cast<T>(detail::rounded(s * s * detail::normalized(detail::four_minus_pi_half_hi, detail::four_minus_pi_half_lo)));
        }
        constexpr T skewness() const noexcept { return static_cast<T>(detail::rayleigh_skewness_hi); }
        constexpr T excess_kurtosis() const noexcept { return static_cast<T>(detail::rayleigh_kurtosis_hi); }
        // 1 + log σ − log(2)/2 + γ/2
        constexpr T entropy() const noexcept
        {
            namespace fp = detail::sd;
            auto h = detail::sum({1.0, 0.0}, detail::log_precise(double(sigma_)));
            h = detail::sum(h, {fp::fmul(fp::ln2_hi, -0.5), fp::fmul(fp::ln2_lo, -0.5)});
            h = detail::sum(h, {fp::fmul(fp::euler_gamma_hi, 0.5), fp::fmul(fp::euler_gamma_lo, 0.5)});
            return static_cast<T>(h.hi + h.lo);
        }

        constexpr T pdf(T x) const noexcept
        {
            auto const& R = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0) || v == std::numeric_limits<double>::infinity())
                return 0;
            auto const [z, h] = R.reduced(v);
            if (h.hi > 1500)
                return 0;
            return static_cast<T>(detail::rounded(z * detail::exp_scaled(-h.hi, -h.lo) / double(R.sigma_)));
        }
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& R = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0) || v == std::numeric_limits<double>::infinity())
                return -std::numeric_limits<T>::infinity();
            auto const [z, h] = R.reduced(v);
            auto l = detail::log(z);
            auto const ls = detail::log_precise(double(R.sigma_));
            l = detail::sum(l, {-ls.hi, -ls.lo});
            l = detail::sum(l, {-h.hi, -h.lo});
            return static_cast<T>(l.hi + l.lo);
        }
        constexpr T cdf(T x) const noexcept { auto const& R = *this; return R.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& R = *this; return R.distribution(x, true); }

        // σ √(−2 log(1 − p)) ; cquantile : σ √(−2 log q)
        constexpr T quantile(double p) const { auto const& R = *this; return R.inverse(p, false); }
        constexpr T cquantile(double q) const { auto const& R = *this; return R.inverse(q, true); }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& R = *this;
            double const e = wsl::random::sample(exponential<double>({.rate = 1.0}), g);
            return static_cast<T>(detail::sd::fmul(double(R.sigma_), wsl::special::sqrt(detail::sd::fmul(2.0, e))));
        }

        friend constexpr bool operator==(rayleigh const&, rayleigh const&) = default;

    private:
        constexpr T times(double c_hi, double c_lo) const noexcept
        {
            return static_cast<T>(detail::rounded(detail::normalized(double(sigma_), 0.0) * detail::normalized(c_hi, c_lo)));
        }

        // z = x/σ à l'échelle, h = z²/2 en double-double (+∞ si hors de portée).
        struct parts
        {
            detail::scaled z;
            detail::sd::double_double h;
        };
        constexpr parts reduced(double x) const noexcept
        {
            detail::scaled const z = detail::normalized(x, 0.0) / double(sigma_);
            detail::scaled h = z * z;
            h.e -= 1;
            if (h.e > 1020)
                return {z, {std::numeric_limits<double>::infinity(), 0.0}};
            if (h.e < -1020)
                return {z, {0.0, 0.0}};
            return {z, detail::unscaled(h)};
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
            detail::scaled const z = detail::normalized(v, 0.0) / double(sigma_);
            detail::scaled h = z * z;
            h.e -= 1;
            if (h.e < -900)                                                 // 1 − e^{−h} = h (1 − h/2 …)
                return upper ? T(1) : static_cast<T>(detail::rounded(h));
            if (h.e > 11)                                                   // h > 4096
                return upper ? T(0) : T(1);
            auto const [h_hi, h_lo] = detail::unscaled(h);
            if (upper)
                return static_cast<T>(detail::rounded(detail::exp_scaled(-h_hi, -h_lo)));
            return h_hi > 40 ? T(1) : static_cast<T>(-detail::expm1_dd(-h_hi, -h_lo));
        }

        constexpr T inverse(double t, bool upper) const
        {
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::random::quantile : probabilité hors de [0, 1]");
            if (t == 0)
                return upper ? std::numeric_limits<T>::infinity() : T(0);
            if (t == 1)
                return upper ? T(0) : std::numeric_limits<T>::infinity();
            auto const l = upper ? detail::log_precise(t) : detail::log1p_dd(-t);   // < 0
            detail::scaled w = detail::normalized(-l.hi, -l.lo);
            w.e += 1;                                                       // −2 log(…)
            return static_cast<T>(detail::rounded(detail::sqrt(w) * detail::normalized(double(sigma_), 0.0)));
        }

        T sigma_;
    };
}
