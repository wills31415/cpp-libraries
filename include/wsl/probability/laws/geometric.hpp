#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Loi géométrique de paramètre p dans (0, 1] : nombre d'échecs avant le premier
// succès d'épreuves de Bernoulli indépendantes, P(X = k) = p (1 − p)^k, k >= 0
// (la convention de std::geometric_distribution et de R ; le nombre d'essais
// est X + 1).
//
//   wsl::probability::laws::geometric X({.p = 0.2});
//   X.mean(); X.pmf(3); X.cdf(3); X.quantile(0.5); X(g);
//
// Masse, répartition, survie : (1 − p)^k = e^{k log1p(−p)}, k log1p(−p) en
// double-double (k exact), un seul arrondi. Quantile exact de la vraie loi,
// ⌈log1p(−u)/log1p(−p)⌉ − 1 avec le quotient en double-double (sauf si ce
// quotient tombe à 2^-100 près, relativement, d'un entier).
//
// Tirage par inversion, ⌊log(U)/log1p(−p)⌋ : approché (deux arrondis), bit-exact.

#include <wsl/probability/laws/basic.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/random/detail/portable_math.hpp>

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_integer I = std::int64_t>
    class geometric : public law_interface
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::approximate;

        struct probability
        {
            double p;
        };
        using param_type = probability;
        constexpr param_type param() const { return {.p = p()}; }

        constexpr geometric(probability q) : p_(q.p)
        {
            if (!(p_ > 0 && p_ <= 1))
                throw std::domain_error("wsl::probability::laws::geometric : il faut 0 < p <= 1");
            if (p_ < 1)
                log_q_ = detail::log1p_dd(-p_);
        }

        constexpr double p() const noexcept { return p_; }
        constexpr interval<I> support() const noexcept { return {I(0), std::numeric_limits<I>::max()}; }

        // (1 − p)/p.
        constexpr double mean() const noexcept
        {
            if (p_ == 1)
                return 0;
            return detail::rounded(failure() / p_);
        }
        // (1 − p)/p².
        constexpr double variance() const noexcept
        {
            if (p_ == 1)
                return 0;
            return detail::rounded(failure() / p_ / p_);
        }
        constexpr I mode() const noexcept { return 0; }
        constexpr I median() const noexcept { return wsl::probability::quantile(*this, 0.5); }
        constexpr double skewness() const noexcept
        {
            if (p_ == 1)                                                    // loi dégénérée : 0/0
                return std::numeric_limits<double>::quiet_NaN();
            return (2 - p_) / wsl::random::detail::pm::sqrt(1 - p_);
        }
        constexpr double excess_kurtosis() const noexcept
        {
            if (p_ == 1)
                return std::numeric_limits<double>::quiet_NaN();
            return 6 + detail::sd::fmul(p_, p_) / (1 - p_);
        }
        // −[(1 − p) log(1 − p) + p log p]/p = −((1 − p)/p) log(1 − p) − log p.
        constexpr double entropy() const noexcept
        {
            namespace sd = detail::sd;
            auto const [l_hi, l_lo] = detail::log_precise(p_);
            if (p_ == 1)
                return 0;
            auto const [r_hi, r_lo] = detail::unscaled(failure() / p_);    // (1 − p)/p, normal ou nul
            auto const [m, m_error] = detail::product(r_hi, -log_q_.hi);
            double const m_lo = m_error + (sd::fmul(r_hi, -log_q_.lo) - sd::fmul(r_lo, log_q_.hi));
            auto const [s, s_error] = sd::two_sum(m, -l_hi);
            return s + (s_error + (m_lo - l_lo));
        }

        constexpr double pmf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 0)
                return 0;
            if (X.p_ == 1)
                return k == 0 ? 1.0 : 0.0;
            auto const [m_hi, m_lo] = X.log_survival(static_cast<std::uint64_t>(k));   // k log(1 − p)
            if (m_hi < -746)
                return 0;
            if (m_hi > -0x1p-900)                                            // p (1 + k log(1 − p)) s'arrondit à p
                return X.p_;
            return detail::rounded(detail::exp_scaled(m_hi, m_lo) * detail::normalized(X.p_, 0.0));
        }

        constexpr double log_pmf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 0)
                return -std::numeric_limits<double>::infinity();
            if (X.p_ == 1)
                return k == 0 ? 0.0 : -std::numeric_limits<double>::infinity();
            auto const [m_hi, m_lo] = X.log_survival(static_cast<std::uint64_t>(k));
            auto const [l_hi, l_lo] = detail::log_precise(X.p_);
            auto const [s, s_error] = detail::sd::two_sum(l_hi, m_hi);
            return s + (s_error + (l_lo + m_lo));
        }

        // 1 − (1 − p)^(k+1) = −expm1((k + 1) log(1 − p)).
        constexpr double cdf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 0)
                return 0;
            if (X.p_ == 1)
                return 1;
            auto const [m_hi, m_lo] = X.log_survival(static_cast<std::uint64_t>(k) + 1);
            if (m_hi < -40)                                                  // e^{−40} < 2^-57
                return 1;
            if (m_hi > -0x1p-900)
                return -m_hi;
            return -detail::expm1_dd(m_hi, m_lo);
        }

        // (1 − p)^(k+1).
        constexpr double ccdf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 0)
                return 1;
            if (X.p_ == 1)
                return 0;
            auto const [m_hi, m_lo] = X.log_survival(static_cast<std::uint64_t>(k) + 1);
            if (m_hi < -746)
                return 0;
            if (m_hi > -0x1p-900)
                return 1;
            return detail::rounded(detail::exp_scaled(m_hi, m_lo));
        }

        // Plus petit k tel que P(X <= k) >= u : ⌈log1p(−u)/log1p(−p)⌉ − 1.
        constexpr I quantile(double u) const
        {
            auto const& X = *this;
            if (!(u >= 0 && u <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            if (u == 0 || X.p_ == 1)
                return 0;
            if (u == 1)
                return std::numeric_limits<I>::max();
            auto const [a_hi, a_lo] = detail::log1p_dd(-u);
            detail::scaled const t = detail::normalized(-a_hi, -a_lo) / detail::normalized(-X.log_q_.hi, -X.log_q_.lo);
            if (t.e >= 62)
                return std::numeric_limits<I>::max();
            if (t.e < 0)                                                    // 0 < t < 1 : ⌈t⌉ = 1 (t peut sous-dépasser)
                return 0;
            auto const [t_hi, t_lo] = detail::unscaled(t);
            std::uint64_t const k = detail::ceiling(t_hi, t_lo) - 1;        // t > 0 : ⌈t⌉ >= 1
            return k >= static_cast<std::uint64_t>(std::numeric_limits<I>::max()) ? std::numeric_limits<I>::max() : static_cast<I>(k);
        }

        template <std::uniform_random_bit_generator G>
        constexpr I sample(G& g) const
        {
            auto const& X = *this;
            if (X.p_ == 1)
                return 0;
            double const r = wsl::random::detail::pm::log(detail::uniform_positive(g)) / X.log_q_.hi;   // >= 0
            if (!(r < 0x1p62) || r >= static_cast<double>(std::numeric_limits<I>::max()))
                return std::numeric_limits<I>::max();
            return static_cast<I>(r);                                       // ⌊r⌋
        }

        friend constexpr bool operator==(geometric const& a, geometric const& b) noexcept { return a.p_ == b.p_; }   // log(1 − p) s'en déduit

    private:
        // 1 − p, à l'échelle (exact : deux doubles).
        constexpr detail::scaled failure() const noexcept
        {
            auto const [q, q_error] = detail::sd::two_sum(1.0, -p_);
            return detail::normalized(q, q_error);
        }

        // n log(1 − p) en double-double, n entier.
        constexpr detail::sd::double_double log_survival(std::uint64_t n) const noexcept
        {
            if (n == 0)
                return {0.0, 0.0};
            return detail::times(detail::exact_split(n), log_q_);
        }

        double p_;
        detail::sd::double_double log_q_{0.0, 0.0};                         // log(1 − p) ; 0 pour p = 1
    };
}
