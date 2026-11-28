#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Algèbre des lois : loi réelle tronquée (conditionnée) à [a, b] (voir design/random.md §4.6).
//
//   wsl::probability::normal N({.mean = 0, .stddev = 1});
//   auto T = wsl::probability::truncate(N, {.lower = 10});   // N | N >= 10, masse ~ 7,6·10⁻²⁴
//   T.cdf(10.5); T.quantile(0.5); T(g);
//
// Masse m = P(a <= X <= b), F_T(x) = P(a <= X <= x)/m et S_T(x) = P(x < X <= b)/m sont
// des différences de répartitions F ou de survies S de X, choisies point par point pour
// que les deux termes soient <= ½ quand c'est possible (pas de compensation : la
// troncature loin dans une queue garde sa précision). Quantiles : F(x) = F(a) + p m ou
// S(x) = S(b) + (1 − p) m, sommes de positifs, par quantile ou cquantile de X selon que
// la cible est sous ½ ou non. Composition, bit-exacte, erreur documentée (design/random.md) :
// elle croît près des bornes (conditionnement |x|/(x − a)) et quand une valeur de F ou de
// S passe sous les sous-normaux. Pas de moments (pas de forme générale).
//
// Tirage : par rejet sur X si m >= 1/4 (au plus 4 tirages en moyenne, exact si celui
// de X l'est) ; par inversion sinon (approché). Le type est donc marqué approximate.

#include <wsl/probability/affine.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability
{
    template <real_law L>
        requires has_cdf<L> && has_ccdf<L>
    class truncated : public law_interface
    {
    public:
        using value_type = typename L::value_type;
        using law_type = L;
        static constexpr exactness sampling = exactness::approximate;

        constexpr truncated(L const& law, value_type lower, value_type upper) : law_(law)
        {
            if (!(lower < upper))
                throw std::domain_error("wsl::probability::truncate : il faut lower < upper");
            auto const [lo, up] = law_.support();
            low_ = lower > lo ? lower : lo;
            high_ = upper < up ? upper : up;
            fa_ = wsl::probability::cdf(law_, low_);
            sa_ = wsl::probability::ccdf(law_, low_);
            fb_ = wsl::probability::cdf(law_, high_);
            sb_ = wsl::probability::ccdf(law_, high_);
            if (fb_ <= 0.5)
                mass_ = fb_ - fa_;
            else if (sa_ <= 0.5)
                mass_ = sa_ - sb_;
            else
                mass_ = 1.0 - (fa_ + sb_);                                  // à cheval sur la médiane : la plus grande part
            if (!(low_ <= high_) || !(mass_ > 0))
                throw std::domain_error("wsl::probability::truncate : l'intervalle doit avoir une probabilité non nulle");
        }

        constexpr L const& law() const noexcept { return law_; }
        constexpr value_type lower() const noexcept { return low_; }
        constexpr value_type upper() const noexcept { return high_; }
        constexpr double mass() const noexcept { return mass_; }               // P(a <= X <= b)
        constexpr interval<value_type> support() const noexcept { return {low_, high_}; }

        constexpr value_type median() const requires has_quantile<L> && has_cquantile<L> { return wsl::probability::quantile(*this, 0.5); }

        constexpr value_type pdf(value_type x) const noexcept requires has_pdf<L>
        {
            auto const& Y = *this;
            if (!(x == x))
                return static_cast<value_type>(detail::sd::opaque(double(x)));   // NaN opaque : GCC fusionnerait x + x chez l'appelant
            if (!(x >= Y.low_ && x <= Y.high_))
                return 0;
            return static_cast<value_type>(double(wsl::probability::pdf(Y.law_, x)) / Y.mass_);
        }
        constexpr value_type log_pdf(value_type x) const noexcept requires requires(L const& X, value_type v) { wsl::probability::log_pdf(X, v); }
        {
            auto const& Y = *this;
            if (!(x == x))
                return static_cast<value_type>(detail::sd::opaque(double(x)));   // NaN opaque : GCC fusionnerait x + x chez l'appelant
            if (!(x >= Y.low_ && x <= Y.high_))
                return -std::numeric_limits<value_type>::infinity();
            auto const lm = detail::log_precise(Y.mass_);
            auto const l = detail::sum({double(wsl::probability::log_pdf(Y.law_, x)), 0.0}, {-lm.hi, -lm.lo});
            return static_cast<value_type>(l.hi + l.lo);
        }

        constexpr value_type cdf(value_type x) const noexcept
        {
            auto const& Y = *this;
            if (!(x == x))
                return static_cast<value_type>(detail::sd::opaque(double(x)));   // NaN opaque : GCC fusionnerait x + x chez l'appelant
            if (x < Y.low_)
                return 0;
            if (x >= Y.high_)
                return 1;
            double const f = wsl::probability::cdf(Y.law_, x);                                 // P(a <= X <= x)
            double const d = f <= 0.5 || Y.sa_ > 0.5 ? f - Y.fa_ : Y.sa_ - double(wsl::probability::ccdf(Y.law_, x));
            return clamp01(d / Y.mass_);
        }
        constexpr value_type ccdf(value_type x) const noexcept
        {
            auto const& Y = *this;
            if (!(x == x))
                return static_cast<value_type>(detail::sd::opaque(double(x)));   // NaN opaque : GCC fusionnerait x + x chez l'appelant
            if (x < Y.low_)
                return 1;
            if (x >= Y.high_)
                return 0;
            double const t = wsl::probability::ccdf(Y.law_, x);                                // P(x < X <= b)
            double const d = t <= 0.5 || Y.fb_ > 0.5 ? t - Y.sb_ : Y.fb_ - double(wsl::probability::cdf(Y.law_, x));
            return clamp01(d / Y.mass_);
        }

        constexpr value_type quantile(double p) const requires has_quantile<L> && has_cquantile<L>
        {
            auto const& Y = *this;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            return Y.solve(p, false);
        }
        constexpr value_type cquantile(double q) const requires has_quantile<L> && has_cquantile<L>
        {
            auto const& Y = *this;
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::probability::cquantile : q doit être dans [0, 1]");
            return Y.solve(q, true);
        }

        template <std::uniform_random_bit_generator G>
            requires samplable<L, G>
        constexpr value_type sample(G& g) const
        {
            auto const& Y = *this;
            if constexpr (has_quantile<L> && has_cquantile<L>)
                if (Y.mass_ < 0.25)
                    return Y.solve(uniform_real<double>(g, wsl::exact), false);
            for (;;)
            {
                value_type const x = detail::native_sample(Y.law_, g);
                if (x >= Y.low_ && x <= Y.high_)
                    return x;
            }
        }
        template <std::uniform_random_bit_generator G>
            requires samplable<L, G>
        constexpr value_type sample(G& g, approximate_t) const
        {
            auto const& Y = *this;
            if constexpr (has_quantile<L> && has_cquantile<L>)
                if (Y.mass_ < 0.25)
                    return Y.solve(uniform_real<double>(g, wsl::exact), false);
            for (;;)
            {
                value_type const x = wsl::probability::sample(Y.law_, g, wsl::approximate);
                if (x >= Y.low_ && x <= Y.high_)
                    return x;
            }
        }

        friend constexpr bool operator==(truncated const& p, truncated const& q) noexcept
        {
            return p.law_ == q.law_ && p.low_ == q.low_ && p.high_ == q.high_;
        }

    private:
        static constexpr value_type clamp01(double v) noexcept { return static_cast<value_type>(v < 0 ? 0.0 : v > 1 ? 1.0 : v); }

        // x tel que F_T(x) = t (S_T(x) = t si upper) : F(x) = F(a) + p m ou S(x) = S(b) + (1 − p) m
        // (p = F_T(x)), sommes de positifs, la cible prise sous ½.
        constexpr value_type solve(double t, bool upper) const
        {
            namespace fp = detail::sd;
            if (t == 0)
                return upper ? high_ : low_;
            if (t == 1)
                return upper ? low_ : high_;
            double const p = upper ? 1.0 - t : t, q = upper ? t : 1.0 - t;   // F_T, S_T (l'un des deux exact)
            double const target_f = fa_ + fp::fmul(p, mass_);
            value_type const x = target_f <= 0.5 ? wsl::probability::quantile(law_, target_f) : wsl::probability::cquantile(law_, sb_ + fp::fmul(q, mass_));
            return x < low_ ? low_ : x > high_ ? high_ : x;
        }

        L law_;
        value_type low_{}, high_{};
        double fa_ = 0, sa_ = 1, fb_ = 1, sb_ = 0, mass_ = 1;                // F(a), S(a), F(b), S(b), m
    };

    // Bornes d'une troncature, chacune facultative : truncate(P, {.lower = 0}).
    struct truncation_bounds
    {
        double lower = -std::numeric_limits<double>::infinity(), upper = std::numeric_limits<double>::infinity();
    };

    // Loi de X sachant lower <= X <= upper : P conditionnée à [lower, upper].
    template <real_law L>
        requires has_cdf<L> && has_ccdf<L>
    constexpr truncated<L> truncate(L const& P, truncation_bounds b)
    {
        using T = typename L::value_type;
        return truncated<L>(P, static_cast<T>(b.lower), static_cast<T>(b.upper));
    }
}
