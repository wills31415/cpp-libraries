#pragma once

// Loi binomiale : nombre de succès en n épreuves de Bernoulli indépendantes de
// paramètre p, P(X = k) = C(n, k) p^k (1 − p)^(n−k), 0 <= k <= n.
//
//   wsl::random::binomial X({.trials = 20, .p = 0.3});
//   X.mean(); pmf(X, 6); log_pmf(X, 6); sample(X, g);
//
// Masse : décomposition de Loader (voir detail/loader.hpp),
//   e^{δ(n) − δ(k) − δ(n−k) − bd0(k, np) − bd0(n−k, n(1−p))} √(n / (2π k (n−k))),
// np et n(1 − p) en double-double, un seul arrondi. Répartition P(X <= k) =
// 1 − I_p(k + 1, n − k) et survie I_p(k + 1, n − k) par la bêta incomplète de
// wsl::special (p tel quel) ; quantile exact (voir detail/discrete.hpp).
//
// Tirage (approché, bit-exact), avec p' = min(p, 1 − p) puis symétrie :
// inversion séquentielle si n p' < 10, BTRS (rejet transformé, Hörmann 1993)
// sinon, en O(1).

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/discrete.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/law.hpp>
#include <wsl/special.hpp>

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    namespace detail
    {
    }

    template <law_integer I = std::int64_t>
    class binomial
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            I trials;
            double p;
        };

        constexpr binomial(parameters q) : n_(q.trials), p_(q.p)
        {
            namespace sd = detail::sd;
            if (!(n_ >= 0) || !(p_ >= 0 && p_ <= 1))
                throw std::domain_error("wsl::random::binomial : il faut trials >= 0 et 0 <= p <= 1");
            // Tirage : p' = min(p, 1 − p).
            mirrored_ = p_ > 0.5;
            double const pp = mirrored_ ? 1 - p_ : p_, qq = 1 - pp, n = static_cast<double>(n_);
            pp_ = pp;
            if (n_ == 0 || pp == 0)
                return;
            if (sd::fmul(n, pp) < 10)
            {
                start_ = wsl::special::exp(sd::fmul(n, wsl::special::log1p(-pp)));
                ratio_ = pp / qq;
                return;
            }
            double const spq = wsl::special::sqrt(sd::fmul(sd::fmul(n, pp), qq));
            b_ = 1.15 + sd::fmul(2.53, spq);
            a_ = (-0.0873 + sd::fmul(0.0248, b_)) + sd::fmul(0.01, pp);
            c_ = sd::fmul(n, pp) + 0.5;
            vr_ = 0.92 - 4.2 / b_;
            alpha_ = sd::fmul(2.83 + 5.1 / b_, spq);
            lpq_ = wsl::special::log(pp / qq);
            m_ = detail::floor_nonnegative(sd::fmul(n + 1, pp));
            h_ = wsl::special::lgamma(m_ + 1) + wsl::special::lgamma((n - m_) + 1);
        }

        constexpr I trials() const noexcept { return n_; }
        constexpr double p() const noexcept { return p_; }
        constexpr interval<I> support() const noexcept { return {I(0), n_}; }

        constexpr double mean() const noexcept { return detail::sd::fmul(static_cast<double>(n_), p_); }
        // n p (1 − p).
        constexpr double variance() const noexcept
        {
            if (n_ == 0 || p_ == 0 || p_ == 1)
                return 0;
            auto const [q, q_error] = detail::sd::two_sum(1.0, -p_);
            return detail::rounded(detail::normalized(static_cast<double>(n_), 0.0) * detail::normalized(p_, 0.0) * detail::normalized(q, q_error));
        }
        constexpr I mode() const noexcept
        {
            double const m = detail::floor_nonnegative(detail::sd::fmul(static_cast<double>(n_) + 1, p_));
            return m > static_cast<double>(n_) ? n_ : static_cast<I>(m);
        }
        constexpr double skewness() const noexcept
        {
            double const q = 1 - p_;
            return (q - p_) / wsl::special::sqrt(detail::sd::fmul(detail::sd::fmul(static_cast<double>(n_), p_), q));
        }
        constexpr double excess_kurtosis() const noexcept
        {
            double const pq = detail::sd::fmul(p_, 1 - p_);
            return (1 - detail::sd::fmul(6.0, pq)) / detail::sd::fmul(static_cast<double>(n_), pq);
        }

        friend constexpr double pmf(binomial const& X, I k) noexcept
        {
            if (k < 0 || k > X.n_)
                return 0;
            if (X.p_ == 0 || X.p_ == 1 || X.n_ == 0)
                return k == (X.p_ == 1 ? X.n_ : I(0)) ? 1.0 : 0.0;
            auto const [e_hi, e_lo] = X.log_mass(k);
            if (e_hi < -746)
                return 0;
            detail::scaled value = detail::exp_scaled(e_hi, e_lo);
            if (k != 0 && k != X.n_)
                value = value * X.prefactor(k);
            return detail::rounded(value);
        }

        friend constexpr double log_pmf(binomial const& X, I k) noexcept
        {
            if (k < 0 || k > X.n_)
                return -std::numeric_limits<double>::infinity();
            if (X.p_ == 0 || X.p_ == 1 || X.n_ == 0)
                return k == (X.p_ == 1 ? X.n_ : I(0)) ? 0.0 : -std::numeric_limits<double>::infinity();
            auto e = X.log_mass(k);
            if (k != 0 && k != X.n_)
                e = detail::sum(e, detail::log(X.prefactor(k)));
            return e.hi + e.lo;
        }

        // P(X <= k) = 1 − I_p(k + 1, n − k).
        friend constexpr double cdf(binomial const& X, I k) noexcept
        {
            if (k < 0)
                return 0;
            if (k >= X.n_ || X.p_ == 0)
                return 1;
            if (X.p_ == 1)
                return 0;
            return wsl::special::ibetac(static_cast<double>(k) + 1, static_cast<double>(X.n_ - k), X.p_);
        }

        // P(X > k) = I_p(k + 1, n − k).
        friend constexpr double ccdf(binomial const& X, I k) noexcept
        {
            if (k < 0)
                return 1;
            if (k >= X.n_ || X.p_ == 0)
                return 0;
            if (X.p_ == 1)
                return 1;
            return wsl::special::ibeta(static_cast<double>(k) + 1, static_cast<double>(X.n_ - k), X.p_);
        }

        // Plus petit k tel que P(X <= k) >= u.
        friend constexpr I quantile(binomial const& X, double u)
        {
            namespace sd = detail::sd;
            if (!(u >= 0 && u <= 1))
                throw std::domain_error("wsl::random::quantile : p doit être dans [0, 1]");
            if (u == 0 || X.p_ == 0 || X.n_ == 0)
                return 0;
            if (u == 1 || X.p_ == 1)
                return X.n_;
            double const n = static_cast<double>(X.n_), q = 1 - X.p_;
            double const z = detail::normal_quantile_guess(u);
            double const mean = sd::fmul(n, X.p_);
            double guess = (mean + sd::fmul(z, wsl::special::sqrt(sd::fmul(mean, q))))
                         + sd::fmul(sd::fmul(z, z) - 1.0, q - X.p_) / 6;     // Cornish–Fisher
            guess = guess < 0 ? 0 : guess > n ? n : guess;
            auto const ge = [&](I k) {                                      // P(X <= k) >= u
                if (k >= X.n_)
                    return true;
                auto const [value, complement] = sd::incomplete_beta(static_cast<double>(k) + 1, static_cast<double>(X.n_ - k), X.p_);
                return detail::at_least(value, complement, u);              // complement : value = 1 − I = F
            };
            return detail::smallest_true(static_cast<I>(guess), I(0), X.n_, ge);
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr I sample(binomial const& X, G& g)
        {
            namespace sd = detail::sd;
            if (X.n_ == 0 || X.pp_ == 0)
                return X.mirrored_ ? X.n_ : I(0);
            double const n = static_cast<double>(X.n_);
            I k = 0;
            if (X.b_ == 0)                                                  // inversion : P(X <= k) >= U
            {
                double const u = uniform_real<double>(g);
                double probability = X.start_, sum = probability;
                while (u > sum && k < X.n_ && probability > 0)
                {
                    probability = sd::fmul(sd::fmul(probability, X.ratio_), static_cast<double>(X.n_ - k)) / static_cast<double>(k + 1);
                    ++k;
                    sum = sum + probability;
                }
            }
            else                                                            // BTRS
            {
                for (;;)
                {
                    double const u = uniform_real<double>(g) - 0.5;
                    double v = uniform_real<double>(g);
                    double const us = 0.5 - (u < 0 ? -u : u);
                    if (us == 0)
                        continue;
                    double const kd = sd::fmul(sd::fmul(2.0, X.a_) / us + X.b_, u) + X.c_;
                    if (!(kd >= 0) || kd >= n + 1)
                        continue;
                    double const kf = detail::floor_nonnegative(kd);
                    if (us >= 0.07 && v <= X.vr_)
                    {
                        k = static_cast<I>(kf);
                        break;
                    }
                    v = wsl::special::log(sd::fmul(v, X.alpha_) / (X.a_ / sd::fmul(us, us) + X.b_));
                    double const bound = (X.h_ - wsl::special::lgamma(kf + 1)) - wsl::special::lgamma((n - kf) + 1)
                                       + sd::fmul(kf - X.m_, X.lpq_);
                    if (v <= bound)
                    {
                        k = static_cast<I>(kf);
                        break;
                    }
                }
            }
            return X.mirrored_ ? static_cast<I>(X.n_ - k) : k;
        }

        friend constexpr bool operator==(binomial const& a, binomial const& b) noexcept { return a.n_ == b.n_ && a.p_ == b.p_; }

    private:
        // Exposant de la masse, en double-double : n log(1 − p) (k = 0), n log p
        // (k = n), sinon δ(n) − δ(k) − δ(n−k) − bd0(k, np) − bd0(n−k, n(1−p)).
        constexpr detail::sd::double_double log_mass(I k) const noexcept
        {
            namespace sd = detail::sd;
            auto const nn = detail::exact_split(static_cast<std::uint64_t>(n_));
            if (k == 0)
                return detail::times(nn, detail::log1p_dd(-p_));
            if (k == n_)
                return detail::times(nn, detail::log_precise(p_));
            double const n = static_cast<double>(n_), x = static_cast<double>(k), y = static_cast<double>(n_ - k);
            auto const [q, q_error] = sd::two_sum(1.0, -p_);
            auto const np = detail::product(n, p_);
            auto const [nq1, nq1_error] = detail::product(n, q);
            auto const [nq_hi, nq_lo] = sd::fast_two_sum(nq1, nq1_error + sd::fmul(n, q_error));
            auto e = detail::stirling_error(n);
            e = detail::sum(e, detail::negated(detail::stirling_error(x)));
            e = detail::sum(e, detail::negated(detail::stirling_error(y)));
            e = detail::sum(e, detail::negated(detail::deviance(x, np.hi, np.lo)));
            e = detail::sum(e, detail::negated(detail::deviance(y, nq_hi, nq_lo)));
            return e;
        }

        // √(n / (2π k (n − k))), 0 < k < n.
        constexpr detail::scaled prefactor(I k) const noexcept
        {
            detail::scaled const r = detail::normalized(static_cast<double>(n_), 0.0)
                                   / (detail::normalized(static_cast<double>(k), 0.0) * detail::normalized(static_cast<double>(n_ - k), 0.0));
            return detail::sqrt(r) * detail::normalized(detail::inv_sqrt_2pi_hi, detail::inv_sqrt_2pi_lo);
        }

        I n_;
        double p_;
        bool mirrored_ = false;                                             // tirage : p > 1/2, X = n − X'
        double pp_ = 0, start_ = 0, ratio_ = 0;                             // p' ; inversion
        double a_ = 0, b_ = 0, c_ = 0, vr_ = 0, alpha_ = 0, lpq_ = 0, m_ = 0, h_ = 0;   // BTRS
    };
}
