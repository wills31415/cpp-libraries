#pragma once

// Loi binomiale négative : nombre d'échecs avant le r-ième succès d'épreuves de Bernoulli
// de paramètre p, r > 0 réel : P(X = k) = Γ(k + r)/(k! Γ(r)) p^r (1 − p)^k, k >= 0.
//
//   wsl::negative_binomial X({.successes = 3.5, .p = 0.4});
//   X.mean(); pmf(X, 6); cdf(X, 6); quantile(X, 0.99); sample(X, g);
//
// Masse : r/(k + r) fois la masse binomiale de r succès en n = k + r épreuves (Loader,
// n réel en double-double) :
//   √(r/(2π k (k + r))) e^{δ(k + r) − δ(r) − δ(k) − bd0(r, np) − bd0(k, nq)}, un arrondi ;
// k = 0 : p^r. Répartition I_p(r, k + 1), survie 1 − I_p(r, k + 1) (bêta incomplète de
// wsl::special) ; quantile exact (detail/discrete.hpp).
//
// Tirage (approché, bit-exact) : mélange gamma–Poisson, λ ~ Γ(r, (1 − p)/p), X ~ P(λ).

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/binomial.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/discrete.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/gamma.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/poisson.hpp>
#include <wsl/special.hpp>

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl
{
    template <law_integer I = std::int64_t>
    class negative_binomial
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            double successes, p;
        };

        constexpr negative_binomial(parameters q) : r_(q.successes), p_(q.p)
        {
            if (!(r_ > 0) || !detail::is_finite(r_) || !(p_ > 0 && p_ <= 1))
                throw std::domain_error("wsl::negative_binomial : il faut successes > 0 fini et p dans (0, 1]");
        }

        constexpr double successes() const noexcept { return r_; }
        constexpr double p() const noexcept { return p_; }
        constexpr interval<I> support() const noexcept { return {I(0), std::numeric_limits<I>::max()}; }
        // r q/p
        constexpr double mean() const noexcept { return detail::rounded(detail::normalized(r_, 0.0) * q_scaled() / p_); }
        // r q/p²
        constexpr double variance() const noexcept { return detail::rounded(detail::normalized(r_, 0.0) * q_scaled() / p_ / p_); }
        // (2 − p)/√(r q)
        constexpr double skewness() const noexcept
        {
            if (p_ == 1)
                return std::numeric_limits<double>::infinity();
            return detail::rounded(detail::normalized(2.0 - p_, 0.0) / detail::sqrt(detail::normalized(r_, 0.0) * q_scaled()));
        }
        // 6/r + p²/(r q)
        constexpr double excess_kurtosis() const noexcept
        {
            if (p_ == 1)
                return std::numeric_limits<double>::infinity();
            namespace sd = detail::sd;
            double const q = 1 - p_;
            return 6 / r_ + sd::fmul(p_, p_) / sd::fmul(r_, q);
        }
        // ⌊(r − 1) q/p⌋ pour r > 1, 0 sinon
        constexpr I mode() const noexcept
        {
            if (!(r_ > 1) || p_ == 1)
                return 0;
            return static_cast<I>(detail::sd::fmul(r_ - 1, 1 - p_) / p_);
        }

        friend constexpr double pmf(negative_binomial const& X, I k) noexcept
        {
            if (k < 0)
                return 0;
            if (X.p_ == 1)
                return k == 0 ? 1.0 : 0.0;
            auto const [e_hi, e_lo] = X.log_mass(k);
            if (!(e_hi > -1500))
                return 0;
            detail::scaled value = detail::exp_scaled(e_hi, e_lo);
            if (k != 0)
                value = value * X.prefactor(k);
            return detail::rounded(value);
        }

        friend constexpr double log_pmf(negative_binomial const& X, I k) noexcept
        {
            if (k < 0)
                return -std::numeric_limits<double>::infinity();
            if (X.p_ == 1)
                return k == 0 ? 0.0 : -std::numeric_limits<double>::infinity();
            auto e = X.log_mass(k);
            if (k != 0)
                e = detail::sum(e, detail::log(X.prefactor(k)));
            return e.hi + e.lo;
        }

        // P(X <= k) = I_p(r, k + 1).
        friend constexpr double cdf(negative_binomial const& X, I k) noexcept
        {
            if (k < 0)
                return 0;
            if (X.p_ == 1)
                return 1;
            return wsl::special::ibeta(X.r_, static_cast<double>(k) + 1, X.p_);
        }

        // P(X > k) = 1 − I_p(r, k + 1).
        friend constexpr double ccdf(negative_binomial const& X, I k) noexcept
        {
            if (k < 0)
                return 1;
            if (X.p_ == 1)
                return 0;
            return wsl::special::ibetac(X.r_, static_cast<double>(k) + 1, X.p_);
        }

        // Plus petit k tel que P(X <= k) >= u.
        friend constexpr I quantile(negative_binomial const& X, double u)
        {
            namespace sd = detail::sd;
            if (!(u >= 0 && u <= 1))
                throw std::domain_error("wsl::quantile : p doit être dans [0, 1]");
            if (u == 0 || X.p_ == 1)
                return 0;
            constexpr I top = std::numeric_limits<I>::max() < (I(1) << 62) ? std::numeric_limits<I>::max() : (I(1) << 62);
            if (u == 1)
                return top;
            double const q = 1 - X.p_;
            double const mean = sd::fmul(X.r_, q) / X.p_;
            double const sd_ = wsl::special::sqrt(sd::fmul(X.r_, q)) / X.p_;
            double const z = detail::normal_quantile_guess(u);
            double const skew = (2 - X.p_) / wsl::special::sqrt(sd::fmul(X.r_, q));
            double guess = (mean + sd::fmul(z, sd_)) + sd::fmul(sd::fmul(sd::fmul(z, z) - 1.0, skew), sd_) / 6;   // Cornish–Fisher
            guess = guess < 0 ? 0 : guess > 0x1p61 ? 0x1p61 : guess;
            auto const ge = [&](I k) {                                      // I_p(r, k + 1) >= u
                auto const [value, complement] = sd::incomplete_beta(X.r_, static_cast<double>(k) + 1, X.p_);
                return detail::at_least(value, !complement, u);
            };
            return detail::smallest_true(static_cast<I>(guess), I(0), top, ge);
        }

        // Mélange gamma–Poisson.
        template <std::uniform_random_bit_generator G>
        friend constexpr I sample(negative_binomial const& X, G& g)
        {
            if (X.p_ == 1)
                return 0;
            double const scale = (1 - X.p_) / X.p_;
            double const lambda = sample(gamma<double>(typename gamma<double>::shape_scale{X.r_, scale}), g);
            if (!(lambda > 0))
                return 0;
            return sample(poisson<I>({.mean = lambda}), g);
        }

        friend constexpr bool operator==(negative_binomial const& a, negative_binomial const& b) noexcept { return a.r_ == b.r_ && a.p_ == b.p_; }

    private:
        // 1 − p à l'échelle (exact)
        constexpr detail::scaled q_scaled() const noexcept
        {
            auto const [q, q_error] = detail::sd::two_sum(1.0, -p_);
            return q == 0 ? detail::scaled{0.0, 0.0, 0} : detail::normalized(q, q_error);
        }

        // δ(k + r) − δ(r) − δ(k) − bd0(r, n p) − bd0(k, n q), n = k + r ; k = 0 : r log p.
        constexpr detail::sd::double_double log_mass(I k) const noexcept
        {
            namespace sd = detail::sd;
            if (k == 0)
                return detail::mul({r_, 0.0}, detail::log_precise(p_));
            double const x = static_cast<double>(k);
            auto const [n_hi, n_lo] = sd::two_sum(x, r_);
            auto const [q, q_error] = sd::two_sum(1.0, -p_);
            auto const np = detail::mul({n_hi, n_lo}, {p_, 0.0});
            auto const nq = detail::mul({n_hi, n_lo}, {q, q_error});
            auto e = detail::stirling_error(n_hi, n_lo);
            e = detail::sum(e, detail::negated(detail::stirling_error(r_)));
            e = detail::sum(e, detail::negated(detail::stirling_error(x)));
            e = detail::sum(e, detail::negated(detail::deviance(r_, np.hi, np.lo)));
            e = detail::sum(e, detail::negated(detail::deviance(x, nq.hi, nq.lo)));
            return e;
        }

        // √(r/(2π k (k + r)))
        constexpr detail::scaled prefactor(I k) const noexcept
        {
            double const x = static_cast<double>(k);
            auto const [n_hi, n_lo] = detail::sd::two_sum(x, r_);
            detail::scaled const d = detail::normalized(x, 0.0) * detail::normalized(n_hi, n_lo);
            return detail::sqrt(detail::normalized(r_, 0.0) / d) * detail::normalized(detail::inv_sqrt_2pi_hi, detail::inv_sqrt_2pi_lo);
        }

        double r_, p_;
    };
}
