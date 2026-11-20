#pragma once

// Intensités rationnelles exactes (docs/interfaces.md, décision 5) : λ = P/Q · 2^e, P et Q
// entiers impairs premiers entre eux (< 2^53, exacts en double). Une intensité r et une
// échelle s données en double, ou un rapport wr::ratio(a, b), s'y écrivent exactement :
// r = P·2^k, 1/s = 1/(Q·2^k). La forme réduite est unique : deux lois de même λ ont la
// même représentation, donc les mêmes calculs, au bit près.

#include <wsl/random/detail/arithmetic.hpp>

#include <bit>
#include <cstdint>
#include <numeric>
#include <optional>
#include <stdexcept>

namespace wsl::random
{
    // Rapport a/b de deux doubles, pour donner exactement une intensité ou une échelle :
    // gamma({.rate = wr::ratio(1, 3), .shape = 2}). Un double v s'y convertit (v/1).
    class ratio
    {
    public:
        constexpr ratio(double value) noexcept : num_(value), den_(1) {}
        constexpr ratio(double numerator, double denominator) noexcept : num_(numerator), den_(denominator) {}

        constexpr double numerator() const noexcept { return num_; }
        constexpr double denominator() const noexcept { return den_; }

    private:
        double num_, den_;
    };

    namespace detail
    {
        struct exact_rate
        {
            double p = 1, q = 1;                                            // impairs, premiers entre eux
            int e = 0;

            friend constexpr bool operator==(exact_rate const&, exact_rate const&) = default;

            // x > 0 fini : x = m · 2^k, m impair.
            static constexpr std::uint64_t odd_part(double x, int& k) noexcept
            {
                std::uint64_t const b = std::bit_cast<std::uint64_t>(x);
                std::uint64_t m = b & ((std::uint64_t{1} << 52) - 1);
                int const biased = static_cast<int>(b >> 52);
                if (biased == 0)
                    k = -1074;                                              // sous-normal
                else
                {
                    m |= std::uint64_t{1} << 52;
                    k = biased - 1075;
                }
                int const z = std::countr_zero(m);
                k += z;
                return m >> z;
            }

            // a/b, a et b > 0 finis.
            static constexpr exact_rate of(double a, double b) noexcept
            {
                int ka = 0, kb = 0;
                std::uint64_t const ma = odd_part(a, ka), mb = odd_part(b, kb);
                std::uint64_t const g = std::gcd(ma, mb);
                return {static_cast<double>(ma / g), static_cast<double>(mb / g), ka - kb};
            }

            constexpr exact_rate inverse() const noexcept { return {q, p, -e}; }

            // λ en double s'il l'est exactement (Q = 1, sans débordement ni sous-normal inexact), 0 sinon.
            constexpr double as_double() const noexcept
            {
                if (q != 1)
                    return 0;
                double const x = sd::scale_by_power_of_two(p, e);
                if (!(x > 0) || !is_finite(x))
                    return 0;
                int k = 0;
                return odd_part(x, k) == static_cast<std::uint64_t>(p) && k == e ? x : 0.0;
            }

            // v · λ et v / λ, à l'échelle (v > 0).
            constexpr scaled times(scaled const& v) const noexcept
            {
                scaled const r = v * normalized(p, 0.0, e);
                return q == 1 ? r : r / q;
            }
            constexpr scaled per(scaled const& v) const noexcept
            {
                scaled const r = v * normalized(q, 0.0, -e);
                return p == 1 ? r : r / p;
            }
            constexpr scaled value() const noexcept { return times(normalized(1.0, 0.0)); }

            // log λ = log P − log Q + e log 2, en double-double.
            constexpr sd::double_double log() const noexcept
            {
                auto const [m, m_error] = sd::two_product(static_cast<double>(e), sd::ln2_hi);
                sd::double_double s{m, m_error + sd::fmul(static_cast<double>(e), sd::ln2_lo)};
                if (p != 1)
                    s = sum(s, log_precise(p));
                if (q != 1)
                    s = sum(s, negated(log_precise(q)));
                return s;
            }

            // λ/c exactement (c > 0 fini) ; vide si Q dépasse 2^53.
            constexpr std::optional<exact_rate> divided(double c) const noexcept
            {
                int kc = 0;
                std::uint64_t const mc = odd_part(c, kc);
                auto const pi = static_cast<std::uint64_t>(p), qi = static_cast<std::uint64_t>(q);
                std::uint64_t const g = std::gcd(pi, mc);
                std::uint64_t const m = mc / g;
                if (qi > ((std::uint64_t{1} << 53) - 1) / m)
                    return std::nullopt;
                return exact_rate{static_cast<double>(pi / g), static_cast<double>(qi * m), e - kc};
            }

            // X = V/λ = V · σ / d avec σ et d doubles (tirage exact) : σ = Q·2^-e, d = P, ou
            // σ = Q, d = P·2^e si 2^-e déborde.
            struct factors
            {
                double sigma, divisor;
            };
            constexpr factors sampling_factors() const noexcept
            {
                double const s = sd::scale_by_power_of_two(q, -e);
                int k = 0;
                if (s > 0 && is_finite(s) && odd_part(s, k) == static_cast<std::uint64_t>(q) && k == -e)
                    return {s, p};
                return {q, sd::scale_by_power_of_two(p, e)};
            }
        };

        // Intensité donnée par .rate ou par .scale (l'un des deux au plus ; aucun : λ = 1).
        constexpr exact_rate rate_of(std::optional<wsl::random::ratio> const& rate, std::optional<wsl::random::ratio> const& scale, char const* message)
        {
            if (rate && scale)
                throw std::domain_error(message);
            if (!rate && !scale)
                return {};
            wsl::random::ratio const r = rate ? *rate : *scale;
            double const a = r.numerator(), b = r.denominator();
            if (!(a > 0) || !is_finite(a) || !(b > 0) || !is_finite(b))
                throw std::domain_error(message);
            return rate ? exact_rate::of(a, b) : exact_rate::of(b, a);
        }

        // Paramètres équivalents, pour param() : .rate (double) si λ en est un, sinon .scale
        // (double) si 1/λ en est un, sinon .rate = wsl::random::ratio(a, b).
        struct rate_or_scale
        {
            std::optional<wsl::random::ratio> rate, scale;
        };
        constexpr rate_or_scale parameters_of(exact_rate const& l) noexcept
        {
            if (double const r = l.as_double(); r != 0)
                return {r, std::nullopt};
            if (double const s = l.inverse().as_double(); s != 0)
                return {std::nullopt, s};
            double const a = sd::scale_by_power_of_two(l.p, l.e);
            if (a > 0 && is_finite(a) && exact_rate::of(a, l.q) == l)
                return {wsl::random::ratio(a, l.q), std::nullopt};
            return {wsl::random::ratio(l.p, sd::scale_by_power_of_two(l.q, -l.e)), std::nullopt};
        }
    }
}
