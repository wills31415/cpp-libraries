#pragma once

// Lois sur les sous-ensembles, combinaisons, mots et partitions d'entiers. Valeurs :
// std::vector<std::size_t> (éléments croissants d'un sous-ensemble de {0, …, n − 1},
// lettres d'un mot, parts décroissantes d'une partition).
//
//   wsl::random::uniform_subset S({.size = 20});                    // chaque sous-ensemble : 2^-n
//   wsl::random::bernoulli_subset B({.size = 20, .p = 0.1});        // chaque élément avec probabilité p
//   wsl::random::uniform_combination C({.size = 50, .count = 6});   // k-sous-ensembles : 1/C(n, k)
//   wsl::random::uniform_word W({.length = 8, .alphabet = 4});      // k^-ℓ
//   wsl::random::uniform_integer_partition P({.total = 100});           // 1/p(n)
//   wsl::random::uniform_set_partition Q({.size = 10});              // 1/B_n (Bell)
//   wsl::random::uniform_composition K({.total = 10});                   // 2^{1−n} ; {.n, .parts = k} : 1/C(n − 1, k − 1)
//
// Tirages exacts : bits uniformes ; épreuves de Bernoulli exactes ; algorithme de Floyd
// (k entiers uniformes) ; lettres uniformes ; partitions par l'algorithme de Nijenhuis et
// Wilf (1978) — (d, j) de poids d p(m − jd), m p(m) au total, p(·) en grands entiers
// (récurrence pentagonale d'Euler), puis j parts égales à d ; partitions d'ensemble :
// l'élément suivant rejoint chaque bloc existant avec poids T(r − 1, k), ou en ouvre un
// avec poids T(r − 1, k + 1), T(r, k) = k T(r − 1, k) + T(r − 1, k + 1) (nombre de façons
// de finir, r éléments restants, k blocs), en grands entiers (table en O(n²)) ;
// compositions : points de coupure uniformes (sous-ensemble ou combinaison).
// Partition d'ensemble : mot à croissance restreinte a (a₀ = 0, aᵢ <= max(a₀…aᵢ₋₁) + 1,
// aᵢ le bloc de i). Masses en logarithme
// double-double ; 0 hors du support (éléments non croissants, hors de {0, …, n − 1}…).

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/big.hpp>
#include <wsl/random/detail/choice.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/permutations.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace wsl::random
{
    namespace detail
    {
        // Éléments strictement croissants de {0, …, n − 1}.
        constexpr bool is_subset(std::vector<std::size_t> const& s, std::size_t n) noexcept
        {
            for (std::size_t i = 0; i < s.size(); ++i)
                if (s[i] >= n || (i > 0 && s[i] <= s[i - 1]))
                    return false;
            return true;
        }

        // log C(n, k)
        constexpr sd::double_double log_binomial(std::uint64_t n, std::uint64_t k) noexcept
        {
            auto l = log_factorial(n);
            auto const a = log_factorial(k), b = log_factorial(n - k);
            l = sum(l, {-a.hi, -a.lo});
            return sum(l, {-b.hi, -b.lo});
        }

        // log x, x grand entier > 0, en double-double (64 bits de tête).
        constexpr sd::double_double log_big(big const& x) noexcept
        {
            int const b = x.bit_length();
            int const shift = b > 64 ? b - 64 : 0;
            std::uint64_t const top = x.bits_at(shift);
            auto const [hi, lo] = exact_split(top);
            auto l = log(normalized(hi, lo));
            return sum(l, mul({double(shift), 0.0}, {sd::ln2_hi, sd::ln2_lo}));
        }
    }

    class uniform_subset : public law_interface
    {
    public:
        using value_type = std::vector<std::size_t>;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t size;
        };

        constexpr explicit uniform_subset(parameters p) : n_(p.size) {}

        constexpr std::size_t size() const noexcept { return n_; }
        constexpr double mean_size() const noexcept { return detail::sd::fmul(double(n_), 0.5); }
        constexpr double entropy() const noexcept { return detail::sd::fmul(double(n_), detail::sd::ln2_hi); }

        constexpr double pmf(value_type const& s) const noexcept
        {
            auto const& S = *this;
            return detail::is_subset(s, S.n_) ? detail::sd::scale_by_power_of_two(1.0, -static_cast<int>(S.n_ > 2000 ? 2000 : S.n_)) : 0.0;
        }
        constexpr double log_pmf(value_type const& s) const noexcept
        {
            auto const& S = *this;
            if (!detail::is_subset(s, S.n_))
                return -std::numeric_limits<double>::infinity();
            auto const l = detail::mul({-double(S.n_), 0.0}, {detail::sd::ln2_hi, detail::sd::ln2_lo});
            return l.hi + l.lo;
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& S = *this;
            value_type s;
            std::uint64_t word = 0;
            for (std::size_t i = 0; i < S.n_; ++i)
            {
                if (i % 64 == 0)
                    word = random_bits(g);
                if ((word >> (i % 64)) & 1)
                    s.push_back(i);
            }
            return s;
        }

        friend constexpr bool operator==(uniform_subset const&, uniform_subset const&) = default;

    private:
        std::size_t n_;
    };

    class bernoulli_subset : public law_interface
    {
    public:
        using value_type = std::vector<std::size_t>;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t size;
            double p;
        };

        constexpr explicit bernoulli_subset(parameters q) : n_(q.size), p_(q.p)
        {
            if (!(p_ >= 0 && p_ <= 1))
                throw std::domain_error("wsl::random::bernoulli_subset : il faut p dans [0, 1]");
        }

        constexpr std::size_t size() const noexcept { return n_; }
        constexpr double p() const noexcept { return p_; }
        constexpr double mean_size() const noexcept { return detail::sd::fmul(double(n_), p_); }

        // p^k (1 − p)^{n − k}
        constexpr double log_pmf(value_type const& s) const noexcept
        {
            auto const& B = *this;
            auto const l = B.log_mass(s);
            return l.hi + l.lo;
        }
        constexpr double pmf(value_type const& s) const noexcept { auto const& B = *this; return detail::exp_of_dd(B.log_mass(s)); }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& B = *this;
            value_type s;
            for (std::size_t i = 0; i < B.n_; ++i)
                if (bernoulli_trial(B.p_, g))
                    s.push_back(i);
            return s;
        }

        friend constexpr bool operator==(bernoulli_subset const&, bernoulli_subset const&) = default;

    private:
        constexpr detail::sd::double_double log_mass(value_type const& s) const noexcept
        {
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            if (!detail::is_subset(s, n_))
                return none;
            double const k = double(s.size()), m = double(n_) - k;
            detail::sd::double_double l{0.0, 0.0};
            if (k > 0)
            {
                if (p_ == 0)
                    return none;
                l = detail::mul({k, 0.0}, detail::log_precise(p_));
            }
            if (m > 0)
            {
                if (p_ == 1)
                    return none;
                l = detail::sum(l, detail::mul({m, 0.0}, detail::log1p_dd(-p_)));
            }
            return l;
        }

        std::size_t n_;
        double p_;
    };

    class uniform_combination : public law_interface
    {
    public:
        using value_type = std::vector<std::size_t>;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t size, count;
        };

        constexpr explicit uniform_combination(parameters p) : n_(p.size), k_(p.count)
        {
            if (k_ > n_)
                throw std::domain_error("wsl::random::uniform_combination : il faut count <= size");
            log_count_ = detail::log_binomial(n_, k_);
        }

        constexpr std::size_t size() const noexcept { return n_; }
        constexpr std::size_t count() const noexcept { return k_; }
        constexpr double entropy() const noexcept { return log_count_.hi + log_count_.lo; }

        constexpr double log_pmf(value_type const& s) const noexcept
        {
            auto const& C = *this;
            return s.size() == C.k_ && detail::is_subset(s, C.n_) ? -C.entropy() : -std::numeric_limits<double>::infinity();
        }
        constexpr double pmf(value_type const& s) const noexcept
        {
            auto const& C = *this;
            return s.size() == C.k_ && detail::is_subset(s, C.n_) ? detail::exp_of_dd({-C.log_count_.hi, -C.log_count_.lo}) : 0.0;
        }

        // Floyd : pour j = n − k, …, n − 1, t uniforme sur [0, j] ; t s'il est neuf, sinon j.
        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& C = *this;
            value_type s;
            s.reserve(C.k_);
            for (std::size_t j = C.n_ - C.k_; j < C.n_; ++j)
            {
                auto const t = static_cast<std::size_t>(uniform_below(std::uint64_t(j) + 1, g));
                auto const at = std::ranges::lower_bound(s, t);
                if (at != s.end() && *at == t)
                    s.insert(std::ranges::lower_bound(s, j), j);
                else
                    s.insert(at, t);
            }
            return s;
        }

        friend constexpr bool operator==(uniform_combination const& a, uniform_combination const& b) noexcept { return a.n_ == b.n_ && a.k_ == b.k_; }

    private:
        std::size_t n_, k_;
        detail::sd::double_double log_count_{};
    };

    class uniform_word : public law_interface
    {
    public:
        using value_type = std::vector<std::size_t>;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t length, alphabet;
        };

        constexpr explicit uniform_word(parameters p) : length_(p.length), k_(p.alphabet)
        {
            if (k_ < 1)
                throw std::domain_error("wsl::random::uniform_word : il faut alphabet >= 1");
        }

        constexpr std::size_t length() const noexcept { return length_; }
        constexpr std::size_t alphabet() const noexcept { return k_; }
        constexpr double entropy() const noexcept
        {
            auto const l = detail::mul({double(length_), 0.0}, detail::log_precise(double(k_)));
            return l.hi + l.lo;
        }

        constexpr double log_pmf(value_type const& w) const noexcept
        {
            auto const& W = *this;
            return W.valid(w) ? -W.entropy() : -std::numeric_limits<double>::infinity();
        }
        constexpr double pmf(value_type const& w) const noexcept
        {
            auto const& W = *this;
            if (!W.valid(w))
                return 0;
            auto const l = detail::mul({-double(W.length_), 0.0}, detail::log_precise(double(W.k_)));
            return detail::exp_of_dd(l);
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& W = *this;
            value_type w(W.length_);
            for (auto& c : w)
                c = static_cast<std::size_t>(uniform_below(W.k_, g));
            return w;
        }

        friend constexpr bool operator==(uniform_word const&, uniform_word const&) = default;

    private:
        constexpr bool valid(value_type const& w) const noexcept
        {
            if (w.size() != length_)
                return false;
            for (std::size_t c : w)
                if (c >= k_)
                    return false;
            return true;
        }

        std::size_t length_, k_;
    };

    class uniform_set_partition : public law_interface
    {
    public:
        using value_type = std::vector<std::size_t>;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t size;
        };

        constexpr explicit uniform_set_partition(parameters q) : n_(q.size)
        {
            // T(r, k), r = 0 … n, k = 0 … n − r + 1 : T(0, k) = 1, T(r, k) = k T(r − 1, k) + T(r − 1, k + 1).
            table_.resize(n_ + 1);
            for (std::size_t r = 0; r <= n_; ++r)
            {
                std::size_t const width = n_ - r + 2;
                table_[r].assign(width, detail::big(1));
                if (r == 0)
                    continue;
                for (std::size_t k = 0; k < width; ++k)
                    table_[r][k] = detail::big(k) * table_[r - 1][k] + table_[r - 1][k + 1];
            }
            log_count_ = n_ == 0 ? detail::sd::double_double{0.0, 0.0} : detail::log_big(table_[n_][0]);   // B_n = T(n, 0)
        }

        constexpr std::size_t size() const noexcept { return n_; }
        // B_n (en double ; +∞ au-delà).
        constexpr double count() const noexcept { return detail::exp_of_dd(log_count_); }
        constexpr double entropy() const noexcept { return log_count_.hi + log_count_.lo; }

        constexpr double log_pmf(value_type const& a) const noexcept
        {
            auto const& P = *this;
            return P.valid(a) ? -P.entropy() : -std::numeric_limits<double>::infinity();
        }
        constexpr double pmf(value_type const& a) const noexcept
        {
            auto const& P = *this;
            return P.valid(a) ? detail::exp_of_dd({-P.log_count_.hi, -P.log_count_.lo}) : 0.0;
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& P = *this;
            value_type a(P.n_);
            std::size_t blocks = 0;
            for (std::size_t i = 0; i < P.n_; ++i)
            {
                std::size_t const rest = P.n_ - i - 1;                      // éléments après i
                auto const& next = P.table_[rest];
                // z sous T(rest + 1, blocks) = blocks T(rest, blocks) + T(rest, blocks + 1)
                detail::big z = detail::uniform_below_big(g, P.table_[rest + 1][blocks]);
                if (z < next[blocks + 1])
                    a[i] = blocks++;                                        // nouveau bloc
                else
                {
                    auto const [quotient, remainder] = detail::divide(z - next[blocks + 1], next[blocks]);
                    (void)remainder;
                    a[i] = static_cast<std::size_t>(quotient.bits_at(0));
                }
            }
            return a;
        }

        friend constexpr bool operator==(uniform_set_partition const& a, uniform_set_partition const& b) noexcept { return a.n_ == b.n_; }

    private:
        constexpr bool valid(value_type const& a) const noexcept
        {
            if (a.size() != n_)
                return false;
            std::size_t blocks = 0;
            for (std::size_t v : a)
            {
                if (v > blocks)
                    return false;
                if (v == blocks)
                    ++blocks;
            }
            return true;
        }

        std::size_t n_;
        std::vector<std::vector<detail::big>> table_;
        detail::sd::double_double log_count_{};
    };

    class uniform_composition : public law_interface
    {
    public:
        using value_type = std::vector<std::size_t>;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t total;
            std::size_t parts = 0;                                          // 0 : nombre de parts libre
        };

        constexpr explicit uniform_composition(parameters q) : n_(q.total), k_(q.parts)
        {
            if (n_ < 1 || k_ > n_)
                throw std::domain_error("wsl::random::uniform_composition : il faut total >= 1 et parts <= total");
            log_count_ = k_ == 0 ? detail::mul({double(n_ - 1), 0.0}, {detail::sd::ln2_hi, detail::sd::ln2_lo}) : detail::log_binomial(n_ - 1, k_ - 1);
        }

        constexpr std::size_t total() const noexcept { return n_; }
        constexpr std::size_t parts() const noexcept { return k_; }
        constexpr double entropy() const noexcept { return log_count_.hi + log_count_.lo; }

        constexpr double log_pmf(value_type const& c) const noexcept
        {
            auto const& C = *this;
            return C.valid(c) ? -C.entropy() : -std::numeric_limits<double>::infinity();
        }
        constexpr double pmf(value_type const& c) const noexcept
        {
            auto const& C = *this;
            return C.valid(c) ? detail::exp_of_dd({-C.log_count_.hi, -C.log_count_.lo}) : 0.0;
        }

        // Points de coupure dans {1, …, n − 1} : sous-ensemble uniforme, ou combinaison de k − 1.
        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& C = *this;
            auto const cuts = C.k_ == 0 ? wsl::random::sample(uniform_subset({.size = C.n_ - 1}), g) : wsl::random::sample(uniform_combination({.size = C.n_ - 1, .count = C.k_ - 1}), g);
            value_type c;
            std::size_t last = 0;
            for (std::size_t x : cuts)
            {
                c.push_back(x + 1 - last);
                last = x + 1;
            }
            c.push_back(C.n_ - last);
            return c;
        }

        friend constexpr bool operator==(uniform_composition const& a, uniform_composition const& b) noexcept { return a.n_ == b.n_ && a.k_ == b.k_; }

    private:
        constexpr bool valid(value_type const& c) const noexcept
        {
            if (k_ != 0 && c.size() != k_)
                return false;
            std::size_t total = 0;
            for (std::size_t v : c)
            {
                if (v == 0 || v > n_ - total)
                    return false;
                total += v;
            }
            return total == n_;
        }

        std::size_t n_, k_;
        detail::sd::double_double log_count_{};
    };

    class uniform_integer_partition : public law_interface
    {
    public:
        using value_type = std::vector<std::size_t>;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t total;
        };

        constexpr explicit uniform_integer_partition(parameters q) : n_(q.total)
        {
            // p(m) = Σ_{k>=1} (−1)^{k+1} [p(m − k(3k − 1)/2) + p(m − k(3k + 1)/2)]
            p_.assign(n_ + 1, detail::big());
            p_[0] = detail::big(1);
            for (std::size_t m = 1; m <= n_; ++m)
            {
                detail::big plus, minus;
                for (std::size_t k = 1;; ++k)
                {
                    std::size_t const g1 = k * (3 * k - 1) / 2, g2 = k * (3 * k + 1) / 2;
                    if (g1 > m)
                        break;
                    detail::big& target = k % 2 == 1 ? plus : minus;
                    target = target + p_[m - g1];
                    if (g2 <= m)
                        target = target + p_[m - g2];
                }
                p_[m] = plus - minus;
            }
            log_count_ = detail::log_big(p_[n_]);
        }

        constexpr std::size_t total() const noexcept { return n_; }
        // p(n) (arrondi en double ; +∞ au-delà de 1,8·10³⁰⁸).
        constexpr double count() const noexcept { return detail::exp_of_dd(log_count_); }
        constexpr double entropy() const noexcept { return log_count_.hi + log_count_.lo; }

        constexpr double log_pmf(value_type const& parts) const noexcept
        {
            auto const& P = *this;
            return P.valid(parts) ? -P.entropy() : -std::numeric_limits<double>::infinity();
        }
        constexpr double pmf(value_type const& parts) const noexcept
        {
            auto const& P = *this;
            return P.valid(parts) ? detail::exp_of_dd({-P.log_count_.hi, -P.log_count_.lo}) : 0.0;
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& P = *this;
            value_type parts;
            std::size_t m = P.n_;
            while (m > 0)
            {
                // z uniforme sous m p(m) = Σ_{d, j : jd <= m} d p(m − jd)
                detail::big z = detail::uniform_below_big(g, detail::big(m) * P.p_[m]);
                std::size_t chosen_d = 1, chosen_j = 1;
                bool found = false;
                for (std::size_t d = 1; d <= m && !found; ++d)
                    for (std::size_t j = 1; j * d <= m; ++j)
                    {
                        detail::big const w = detail::big(d) * P.p_[m - j * d];
                        if (z < w)
                        {
                            chosen_d = d;
                            chosen_j = j;
                            found = true;
                            break;
                        }
                        z = z - w;
                    }
                for (std::size_t r = 0; r < chosen_j; ++r)
                    parts.push_back(chosen_d);
                m -= chosen_j * chosen_d;
            }
            std::ranges::sort(parts, std::greater<>{});
            return parts;
        }

        friend constexpr bool operator==(uniform_integer_partition const& a, uniform_integer_partition const& b) noexcept { return a.n_ == b.n_; }

    private:
        constexpr bool valid(value_type const& parts) const noexcept
        {
            std::size_t total = 0;
            for (std::size_t i = 0; i < parts.size(); ++i)
            {
                if (parts[i] == 0 || (i > 0 && parts[i] > parts[i - 1]) || parts[i] > n_ - total)
                    return false;
                total += parts[i];
            }
            return total == n_;
        }

        std::size_t n_;
        std::vector<detail::big> p_;
        detail::sd::double_double log_count_{};
    };
}
