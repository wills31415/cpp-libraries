#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Loi catégorielle : P(X = i) = pᵢ / Σ p, i dans {0, …, K − 1} (poids pᵢ >= 0 finis,
// normalisés exactement par leur somme ; des probabilités dont la somme flottante n'est
// pas exactement 1 conviennent donc).
//
//   wsl::probability::laws::categorical C({.p = {0.2, 0.5, 0.3}});
//   C(g); C.pmf(1); C.cdf(1); C.quantile(0.9); C.mean(); C.entropy();
//
// Tirage exact en O(1) : les poids sont des entiers Nᵢ à une même puissance de 2 près ;
// table d'alias de Walker–Vose construite sur des entiers (colonne de capacité N = Σ Nᵢ,
// poids K Nᵢ), une colonne uniforme puis un entier uniforme sous N : P(i) = Nᵢ/N
// exactement. Si K N ne tient pas sur 64 bits : recherche dans les cumuls en grands entiers
// (O(log K)). Quantile exact (cumuls entiers comparés à p N exactement) ; masse,
// répartition, moments, entropie en double-double.

#include <wsl/probability/laws/basic.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/big.hpp>
#include <wsl/random/detail/choice.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/random/detail/portable_math.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_integer I = std::size_t>
    class categorical : public law_interface
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::vector<double> p;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.p = w_}; }

        constexpr explicit categorical(parameters q) : w_(std::move(q.p))
        {
            std::size_t const k = w_.size();
            if (k < 1)
                throw std::domain_error("wsl::probability::laws::categorical : il faut au moins une issue");
            for (double const v : w_)
                if (!(v >= 0) || !detail::is_finite(v))
                    throw std::domain_error("wsl::probability::laws::categorical : il faut des poids >= 0, finis");
            for (double const v : w_)
                total_ = detail::sum(total_, {v, 0.0});
            if (!(total_.hi > 0))
                throw std::domain_error("wsl::probability::laws::categorical : poids tous nuls");
            integer_weights();
        }

        constexpr std::size_t size() const noexcept { return w_.size(); }
        constexpr std::vector<double> const& weights() const noexcept { return w_; }
        constexpr interval<I> support() const noexcept { return {I(0), static_cast<I>(w_.size() - 1)}; }

        constexpr double mean() const noexcept
        {
            auto const s = average([](std::size_t i) { return double(i); });
            return s.hi + s.lo;
        }
        constexpr double variance() const noexcept
        {
            double const m = mean();
            auto const s = average([m](std::size_t i) { double const d = double(i) - m; return detail::sd::fmul(d, d); });
            return s.hi + s.lo;
        }
        constexpr I mode() const noexcept
        {
            std::size_t best = 0;
            for (std::size_t i = 1; i < w_.size(); ++i)
                if (w_[i] > w_[best])
                    best = i;
            return static_cast<I>(best);
        }
        // −Σ pᵢ log pᵢ, pᵢ = wᵢ/W : log W − Σ (wᵢ/W) log wᵢ
        constexpr double entropy() const noexcept
        {
            auto const lw = detail::log(detail::normalized(total_.hi, total_.lo));
            detail::sd::double_double s = lw;
            for (double const v : w_)
                if (v > 0)
                {
                    auto const t = detail::mul(detail::sd::quotient_dd({v, 0.0}, total_), detail::log_precise(v));
                    s = detail::sum(s, {-t.hi, -t.lo});
                }
            return s.hi + s.lo;
        }

        constexpr double pmf(I k) const noexcept
        {
            auto const& C = *this;
            if (k < 0 || static_cast<std::uint64_t>(k) >= C.w_.size())
                return 0;
            auto const q = detail::sd::quotient_dd({C.w_[static_cast<std::size_t>(k)], 0.0}, C.total_);
            return q.hi + q.lo;
        }
        constexpr double log_pmf(I k) const noexcept
        {
            auto const& C = *this;
            double const p = wsl::probability::pmf(C, k);
            if (!(p > 0))
                return -std::numeric_limits<double>::infinity();
            auto const l = detail::sum(detail::log_precise(C.w_[static_cast<std::size_t>(k)]), [&] { auto const t = detail::log(detail::normalized(C.total_.hi, C.total_.lo)); return detail::sd::double_double{-t.hi, -t.lo}; }());
            return l.hi + l.lo;
        }
        // P(X <= k), P(X > k) : sommes de poids, en double-double.
        constexpr double cdf(I k) const noexcept
        {
            auto const& C = *this;
            if (k < 0)
                return 0;
            if (static_cast<std::uint64_t>(k) >= C.w_.size() - 1)
                return 1;
            return C.partial(0, static_cast<std::size_t>(k) + 1);
        }
        constexpr double ccdf(I k) const noexcept
        {
            auto const& C = *this;
            if (k < 0)
                return 1;
            if (static_cast<std::uint64_t>(k) >= C.w_.size() - 1)
                return 0;
            return C.partial(static_cast<std::size_t>(k) + 1, C.w_.size());
        }
        // Plus petit k tel que P(X <= k) >= p, exactement : C_k >= p N en entiers.
        constexpr I quantile(double p) const
        {
            auto const& C = *this;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            if (p == 0)
                return static_cast<I>(C.first_positive());
            auto const [m, e] = detail::dyadic_of(p);                       // p = m 2^e
            detail::big const pm = detail::big(m) * C.cumulative_.back();   // C_k 2^{−e} >= m N
            std::size_t lo = 0, hi = C.w_.size() - 1;
            while (lo < hi)
            {
                std::size_t const mid = lo + (hi - lo) / 2;
                detail::big const c = e < 0 ? C.cumulative_[mid].shifted_left(-e) : C.cumulative_[mid];
                bool const ge = e < 0 ? !(c < pm) : !(c < pm.shifted_left(e));
                if (ge)
                    hi = mid;
                else
                    lo = mid + 1;
            }
            return static_cast<I>(lo);
        }

        template <std::uniform_random_bit_generator G>
        constexpr I sample(G& g) const
        {
            auto const& C = *this;
            std::size_t const k = C.w_.size();
            if (k == 1)
                return 0;
            if (!C.prob_.empty())                                           // alias, 64 bits
            {
                auto const c = static_cast<std::size_t>(random::uniform_below(k, g));
                std::uint64_t const u = random::uniform_below(C.capacity_, g);
                return static_cast<I>(u < C.prob_[c] ? c : C.alias_[c]);
            }
            detail::big const z = detail::uniform_below_big(g, C.cumulative_.back());
            auto const it = std::ranges::upper_bound(C.cumulative_, z);     // premier C_i > z
            return static_cast<I>(it - C.cumulative_.begin());
        }

        friend constexpr bool operator==(categorical const& a, categorical const& b) noexcept { return a.w_ == b.w_; }

    private:
        template <class F>
        constexpr detail::sd::double_double average(F f) const noexcept
        {
            detail::sd::double_double s{0.0, 0.0};
            for (std::size_t i = 0; i < w_.size(); ++i)
                if (w_[i] > 0)
                    s = detail::sum(s, detail::product(w_[i], f(i)));
            return detail::sd::quotient_dd(s, total_);
        }

        constexpr double partial(std::size_t from, std::size_t to) const noexcept
        {
            detail::sd::double_double s{0.0, 0.0};
            for (std::size_t i = from; i < to; ++i)
                s = detail::sum(s, {w_[i], 0.0});
            auto const q = detail::sd::quotient_dd(s, total_);
            double const v = q.hi + q.lo;
            return v > 1 ? 1.0 : v;
        }

        constexpr std::size_t first_positive() const noexcept
        {
            std::size_t i = 0;
            while (w_[i] == 0)
                ++i;
            return i;
        }

        // Nᵢ = mᵢ 2^{eᵢ − e_min} ; cumuls en grands entiers ; alias si K N < 2^63.
        constexpr void integer_weights()
        {
            std::size_t const k = w_.size();
            std::vector<detail::dyadic> d(k);
            int low = std::numeric_limits<int>::max();
            for (std::size_t i = 0; i < k; ++i)
                if (w_[i] > 0)
                {
                    d[i] = detail::dyadic_of(w_[i]);
                    low = d[i].e < low ? d[i].e : low;
                }
            std::vector<detail::big> n(k);
            detail::big c;
            for (std::size_t i = 0; i < k; ++i)
            {
                if (w_[i] > 0)
                    n[i] = detail::big(d[i].m).shifted_left(d[i].e - low);
                c = c + n[i];
                cumulative_.push_back(c);
            }
            if (c.bit_length() + std::bit_width(k) > 63)
                return;
            capacity_ = c.bits_at(0);                                       // N
            std::vector<std::uint64_t> scaled(k);
            for (std::size_t i = 0; i < k; ++i)
                scaled[i] = (n[i].is_zero() ? 0 : n[i].bits_at(0)) * k;      // K Nᵢ < 2^63
            prob_.assign(k, capacity_);
            alias_.resize(k);
            for (std::size_t i = 0; i < k; ++i)
                alias_[i] = i;
            std::vector<std::size_t> small, large;
            for (std::size_t i = 0; i < k; ++i)
                (scaled[i] < capacity_ ? small : large).push_back(i);
            while (!small.empty() && !large.empty())
            {
                std::size_t const s = small.back(), l = large.back();
                small.pop_back();
                prob_[s] = scaled[s];
                alias_[s] = l;
                scaled[l] -= capacity_ - scaled[s];
                if (scaled[l] < capacity_)
                {
                    large.pop_back();
                    small.push_back(l);
                }
            }
        }

        std::vector<double> w_;
        detail::sd::double_double total_{0.0, 0.0};
        std::vector<detail::big> cumulative_;                               // C_i = N_0 + … + N_i
        std::uint64_t capacity_ = 0;
        std::vector<std::uint64_t> prob_;
        std::vector<std::size_t> alias_;
    };
}
