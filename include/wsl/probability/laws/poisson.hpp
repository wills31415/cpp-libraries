#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Loi de Poisson de moyenne μ > 0 : P(N = k) = e^{−μ} μ^k / k!, k >= 0.
//
//   wsl::probability::poisson N({.mean = 3.5});
//   N.variance(); N.pmf(2); N.log_pmf(40); N(g);
//
// Masse : décomposition de Loader, e^{−δ(k) − bd0(k, μ)} / √(2πk) (voir
// detail/loader.hpp), un seul arrondi. Répartition P(N <= k) = Q(k + 1, μ) et survie
// P(k + 1, μ) par la gamma incomplète de wsl::special ; quantile exact (recherche
// depuis Cornish–Fisher, comparaisons en double-double, voir detail/discrete.hpp).
//
// Tirage (approché, bit-exact) : inversion séquentielle pour μ < 10 ; au-delà,
// PTRS, rejet transformé de Hörmann (1993), en O(1).

#include <wsl/probability/laws/basic.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/discrete.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/probability/detail/loader.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/special.hpp>

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_integer I = std::int64_t>
    class poisson : public law_interface
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            double mean;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.mean = mean()}; }

        constexpr poisson(parameters p) : mu_(p.mean)
        {
            namespace sd = detail::sd;
            if (!(mu_ > 0) || !detail::is_finite(mu_))
                throw std::domain_error("wsl::probability::poisson : il faut mean > 0, fini");
            if (mu_ < 10)
                exp_minus_mu_ = wsl::special::exp(-mu_);
            else                                                            // constantes de PTRS
            {
                double const smu = wsl::special::sqrt(mu_);
                b_ = 0.931 + sd::fmul(2.53, smu);
                a_ = -0.059 + sd::fmul(0.02483, b_);
                log_inv_alpha_ = wsl::special::log(1.1239 + 1.1328 / (b_ - 3.4));
                vr_ = 0.9277 - 3.6224 / (b_ - 2);
                log_mu_ = wsl::special::log(mu_);
            }
        }

        constexpr interval<I> support() const noexcept { return {I(0), std::numeric_limits<I>::max()}; }
        constexpr double mean() const noexcept { return mu_; }
        constexpr double variance() const noexcept { return mu_; }
        constexpr I mode() const noexcept
        {
            return mu_ >= static_cast<double>(std::numeric_limits<I>::max()) ? std::numeric_limits<I>::max() : static_cast<I>(mu_);
        }
        constexpr double skewness() const noexcept { return 1 / wsl::special::sqrt(mu_); }
        constexpr double excess_kurtosis() const noexcept { return 1 / mu_; }

        constexpr double pmf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 0)
                return 0;
            if (k == 0)
                return X.mu_ > 746 ? 0.0 : detail::rounded(detail::exp_scaled(-X.mu_, 0.0));
            double const x = static_cast<double>(k);
            auto const [e_hi, e_lo] = exponent(x, X.mu_);
            if (e_hi < -746)                                                // √(2πk) > 1
                return 0;
            detail::scaled const c = detail::normalized(detail::inv_sqrt_2pi_hi, detail::inv_sqrt_2pi_lo);
            return detail::rounded(detail::exp_scaled(e_hi, e_lo) * c / detail::sqrt(detail::normalized(x, 0.0)));
        }

        // −δ(k) − bd0(k, μ) − log(2πk)/2.
        constexpr double log_pmf(I k) const noexcept
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            if (k < 0)
                return -std::numeric_limits<double>::infinity();
            if (k == 0)
                return -X.mu_;
            double const x = static_cast<double>(k);
            auto const [e_hi, e_lo] = exponent(x, X.mu_);
            auto const [l_hi, l_lo] = detail::log_precise(x);
            auto const [s1, f1] = sd::two_sum(e_hi, -sd::half_log_two_pi_hi);
            auto const [s2, f2] = sd::two_sum(s1, sd::fmul(l_hi, -0.5));
            return s2 + ((f1 + f2) + (e_lo - (sd::half_log_two_pi_lo + sd::fmul(l_lo, 0.5))));
        }

        // P(N <= k) = Q(k + 1, μ).
        constexpr double cdf(I k) const noexcept
        {
            auto const& X = *this;
            return k < 0 ? 0.0 : wsl::special::gamma_q(static_cast<double>(k) + 1, X.mu_);
        }

        // P(N > k) = P(k + 1, μ).
        constexpr double ccdf(I k) const noexcept
        {
            auto const& X = *this;
            return k < 0 ? 1.0 : wsl::special::gamma_p(static_cast<double>(k) + 1, X.mu_);
        }

        // Plus petit k tel que P(N <= k) >= p.
        constexpr I quantile(double p) const
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            constexpr I top = std::numeric_limits<I>::max();
            if (p == 0)
                return 0;
            if (p == 1)
                return top;
            double const z = detail::normal_quantile_guess(p);
            double const smu = wsl::special::sqrt(X.mu_);
            double guess = (X.mu_ + sd::fmul(z, smu)) + (sd::fmul(z, z) - 1.0) / 6;   // Cornish–Fisher
            guess = guess < 0 ? 0 : guess >= static_cast<double>(top) / 2 ? static_cast<double>(top) / 2 : guess;
            auto const ge = [&](I k) {                                      // Q(k + 1, μ) >= p
                auto const [value, upper] = sd::incomplete_gamma(static_cast<double>(k) + 1, X.mu_);
                return detail::at_least(value, upper, p);
            };
            return detail::smallest_true(static_cast<I>(guess), I(0), top, ge);
        }

        template <std::uniform_random_bit_generator G>
        constexpr I sample(G& g) const
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            constexpr auto top = std::numeric_limits<I>::max();
            if (X.mu_ < 10)                                                 // inversion : P(N <= k) >= U
            {
                double const u = uniform_real<double>(g, wsl::exact);
                double probability = X.exp_minus_mu_, sum = probability;
                I k = 0;
                while (u > sum && probability > 0 && k < top)
                {
                    ++k;
                    probability = sd::fmul(probability, X.mu_) / static_cast<double>(k);
                    sum = sum + probability;
                }
                return k;
            }
            for (;;)                                                        // PTRS
            {
                double const u = uniform_real<double>(g, wsl::exact) - 0.5;
                double const v = uniform_real<double>(g, wsl::exact);
                double const us = 0.5 - (u < 0 ? -u : u);
                if (us == 0)
                    continue;
                double const kd = sd::fmul(sd::fmul(2.0, X.a_) / us + X.b_, u) + X.mu_ + 0.43;
                if (!(kd >= 0))
                    continue;
                double const k = detail::floor_nonnegative(kd);
                if (k >= static_cast<double>(top))
                    return top;
                if (us >= 0.07 && v <= X.vr_)
                    return static_cast<I>(k);
                if (us < 0.013 && v > us)
                    continue;
                double const lhs = wsl::special::log(v) + X.log_inv_alpha_ - wsl::special::log(X.a_ / sd::fmul(us, us) + X.b_);
                double const rhs = -X.mu_ + sd::fmul(k, X.log_mu_) - wsl::special::lgamma(k + 1);
                if (lhs <= rhs)
                    return static_cast<I>(k);
            }
        }

        friend constexpr bool operator==(poisson const& a, poisson const& b) noexcept { return a.mu_ == b.mu_; }

    private:
        // −δ(x) − bd0(x, μ), x >= 1.
        static constexpr detail::sd::double_double exponent(double x, double mu) noexcept
        {
            auto const [d_hi, d_lo] = detail::stirling_error(x);
            auto const [b_hi, b_lo] = detail::deviance(x, mu, 0.0);
            auto const [s, s_error] = detail::sd::two_sum(-d_hi, -b_hi);
            auto const [hi, lo] = detail::sd::fast_two_sum(s, s_error - (d_lo + b_lo));
            return {hi, lo};
        }

        double mu_;
        double exp_minus_mu_ = 0;                                           // μ < 10
        double a_ = 0, b_ = 0, log_inv_alpha_ = 0, vr_ = 0, log_mu_ = 0;    // μ >= 10
    };
}
