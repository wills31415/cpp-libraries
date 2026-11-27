#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Loi logarithmique (série logarithmique) de paramètre p dans (0, 1) :
// P(X = k) = −p^k/(k log(1 − p)), k >= 1.
//
//   wsl::probability::logarithmic X({.p = 0.9});   X.pmf(3); X.cdf(3); X.quantile(0.99); X(g);
//
// Masse : e^{k log p}/(k L), L = −log1p(−p), en double-double, un arrondi. Répartition :
// F(k) = Σ_{j<=k} p^j/(j L) si elle est la plus petite part, sinon la survie
// Σ_{j>k} p^j/(j L) (termes décroissants) ; coût O(k) ou O(1/(1 − p)). Quantile exact par
// sommes partielles croissantes. Moments : moyenne p/((1 − p) L), variance
// p (L − p)/((1 − p)² L²) (L − p par sa série pour p petit).
//
// Tirage (approché, bit-exact) : algorithme LK de Kemp (1981), deux uniformes.

#include <wsl/probability/laws/basic.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/discrete.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_integer I = std::int64_t>
    class logarithmic : public law_interface
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            double p;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.p = p()}; }

        constexpr logarithmic(parameters q) : p_(q.p)
        {
            if (!(p_ > 0 && p_ < 1))
                throw std::domain_error("wsl::probability::logarithmic : il faut p dans (0, 1)");
            auto const l = detail::log1p_dd(-p_);
            L_ = {-l.hi, -l.lo};
            log_p_ = detail::log_precise(p_);
        }

        constexpr double p() const noexcept { return p_; }
        constexpr interval<I> support() const noexcept { return {I(1), std::numeric_limits<I>::max()}; }
        constexpr I mode() const noexcept { return 1; }
        // p/((1 − p) L)
        constexpr double mean() const noexcept
        {
            return detail::rounded(detail::normalized(p_, 0.0) / (one_minus_p() * detail::normalized(L_.hi, L_.lo)));
        }
        // p (L − p)/((1 − p)² L²)
        constexpr double variance() const noexcept
        {
            detail::scaled const q = one_minus_p(), l = detail::normalized(L_.hi, L_.lo);
            auto const d = l_minus_p();
            return detail::rounded(detail::normalized(p_, 0.0) * detail::normalized(d.hi, d.lo) / (q * q * l * l));
        }

        constexpr double pmf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 1)
                return 0;
            auto const l = detail::mul({double(k), 0.0}, X.log_p_);
            if (!(l.hi > -1500))
                return 0;
            return detail::rounded(detail::exp_scaled(l.hi, l.lo) / (detail::normalized(double(k), 0.0) * detail::normalized(X.L_.hi, X.L_.lo)));
        }
        constexpr double log_pmf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 1)
                return -std::numeric_limits<double>::infinity();
            auto l = detail::mul({double(k), 0.0}, X.log_p_);
            auto const lk = detail::log_precise(double(k));
            auto const lL = detail::log(detail::normalized(X.L_.hi, X.L_.lo));
            l = detail::sum(l, {-lk.hi, -lk.lo});
            l = detail::sum(l, {-lL.hi, -lL.lo});
            return l.hi + l.lo;
        }
        constexpr double cdf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 1)
                return 0;
            auto const [v, is_f] = X.distribution(k);
            return is_f ? detail::rounded(v) : detail::sd::one_minus(v);
        }
        constexpr double ccdf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 1)
                return 1;
            auto const [v, is_f] = X.distribution(k);
            return is_f ? detail::sd::one_minus(v) : detail::rounded(v);
        }

        // Plus petit k tel que P(X <= k) >= u : sommes partielles Σ p^j/j comparées à u L.
        constexpr I quantile(double u) const
        {
            auto const& X = *this;
            if (!(u >= 0 && u <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            if (u == 0)
                return 1;
            constexpr I top = std::numeric_limits<I>::max() < (I(1) << 62) ? std::numeric_limits<I>::max() : (I(1) << 62);
            if (u == 1)
                return top;
            if (u > 0.5)                                                    // par la survie : plus petit k avec S(k) <= 1 − u
            {
                auto const ge = [&](I k) {
                    auto const [v, is_f] = X.distribution(k);
                    return detail::at_least(v, is_f, u);
                };
                double const guess = wsl::special::log1p(-u) / (X.log_p_.hi);    // p^k ≈ 1 − u
                return detail::smallest_true(static_cast<I>(guess > 1 ? (guess < 0x1p61 ? guess : 0x1p61) : 1.0), I(1), top, ge);
            }
            detail::sd::double_double const target = detail::mul({u, 0.0}, X.L_);
            detail::sd::double_double s{0.0, 0.0}, power{1.0, 0.0};
            for (I k = 1; k < top; ++k)
            {
                power = detail::mul(power, {X.p_, 0.0});
                s = detail::sum(s, detail::sd::quotient_dd(power, {double(k), 0.0}));
                auto const d = detail::sum(s, {-target.hi, -target.lo});
                if (d.hi > 0 || (d.hi == 0 && d.lo >= 0))
                    return k;
            }
            return top;
        }

        // Kemp (1981), LK : V = U₁ ; V >= p : 1 ; sinon q = 1 − e^{r U₂} (r = log(1 − p)),
        // V <= q² : ⌊1 + log V/log q⌋ ; V >= q : 1 ; sinon 2.
        template <std::uniform_random_bit_generator G>
        constexpr I sample(G& g) const
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            double const v = uniform_real<double>(g);
            if (v >= X.p_)
                return 1;
            double const u = uniform_real<double>(g);
            double const q = -wsl::special::expm1(sd::fmul(-(X.L_.hi), u));
            if (v <= sd::fmul(q, q))
            {
                double const k = 1 + wsl::special::log(v) / wsl::special::log(q);
                return k >= 0x1p62 ? (I(1) << 62) : static_cast<I>(k);
            }
            return v >= q ? I(1) : I(2);
        }

        friend constexpr bool operator==(logarithmic const& a, logarithmic const& b) noexcept { return a.p_ == b.p_; }

    private:
        constexpr detail::scaled one_minus_p() const noexcept
        {
            auto const [q, q_error] = detail::sd::two_sum(1.0, -p_);
            return detail::normalized(q, q_error);
        }

        // L − p = Σ_{j>=2} p^j/j (série pour p <= ¼, sinon différence en double-double)
        constexpr detail::sd::double_double l_minus_p() const noexcept
        {
            if (p_ > 0.25)
                return detail::sum(L_, {-p_, 0.0});
            detail::sd::double_double s{0.0, 0.0}, power{p_, 0.0};
            for (int j = 2; j < 200; ++j)
            {
                power = detail::mul(power, {p_, 0.0});
                auto const term = detail::sd::quotient_dd(power, {double(j), 0.0});
                s = detail::sum(s, term);
                if (term.hi < detail::sd::fmul(s.hi, 0x1p-112))
                    break;
            }
            return s;
        }

        // Σ_{j<=k} p^j/j / L (is_f) si elle est la plus petite part, sinon Σ_{j>k} p^j/j / L.
        struct part
        {
            detail::scaled value;
            bool is_f;
        };
        constexpr part distribution(I k) const noexcept
        {
            // tête si k <= 64 et sa somme ne dépasse pas L/2 ; sinon la queue
            if (k <= 64)
            {
                detail::sd::double_double s{0.0, 0.0}, power{1.0, 0.0};
                for (I j = 1; j <= k; ++j)
                {
                    power = detail::mul(power, {p_, 0.0});
                    s = detail::sum(s, detail::sd::quotient_dd(power, {double(j), 0.0}));
                }
                auto const half = detail::sum(s, {-detail::sd::fmul(L_.hi, 0.5), -detail::sd::fmul(L_.lo, 0.5)});
                if (half.hi <= 0)
                    return {detail::normalized(s.hi, s.lo) / detail::normalized(L_.hi, L_.lo), true};
            }
            // queue : p^{k+1} Σ_{m>=0} p^m/(k + 1 + m)
            auto const lead = detail::mul({double(k) + 1, 0.0}, log_p_);
            if (!(lead.hi > -1500))
                return {{0.0, 0.0, 0}, false};
            detail::sd::double_double s{0.0, 0.0}, power{1.0, 0.0};
            for (std::int64_t m = 0; m < (std::int64_t{1} << 40); ++m)
            {
                auto const term = detail::sd::quotient_dd(power, {double(k) + 1 + double(m), 0.0});
                s = detail::sum(s, term);
                if (term.hi < detail::sd::fmul(s.hi, 0x1p-112))
                    break;
                power = detail::mul(power, {p_, 0.0});
            }
            return {detail::exp_scaled(lead.hi, lead.lo) * detail::normalized(s.hi, s.lo) / detail::normalized(L_.hi, L_.lo), false};
        }

        double p_;
        detail::sd::double_double L_{}, log_p_{};
    };
}
