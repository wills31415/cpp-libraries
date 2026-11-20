#pragma once

// Grands entiers non signés, minimaux, pour les tirages exacts : addition,
// soustraction, multiplication, décalages, comparaison, division euclidienne.
// Mots de 64 bits, poids faible d'abord, sans mot nul de tête. constexpr.
//
// Ils ne servent qu'à décider des comparaisons et des arrondis de réels
// paresseux : quelques centaines de bits au plus en pratique, l'efficacité n'est
// pas le but.

#include <wsl/random/engines.hpp>

#include <algorithm>
#include <bit>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace wsl::random::detail
{
    class big
    {
    public:
        constexpr big() = default;
        constexpr explicit big(std::uint64_t v)
        {
            if (v != 0)
                words_.push_back(v);
        }

        // Entier dont les chiffres en base 2^64 sont donnés poids fort d'abord.
        static constexpr big from_digits(std::vector<std::uint64_t> const& most_significant_first)
        {
            big r;
            r.words_.assign(most_significant_first.rbegin(), most_significant_first.rend());
            r.trim();
            return r;
        }

        constexpr bool is_zero() const noexcept { return words_.empty(); }

        constexpr int bit_length() const noexcept
        {
            if (words_.empty())
                return 0;
            return static_cast<int>(64 * (words_.size() - 1)) + std::bit_width(words_.back());
        }

        constexpr bool bit(int i) const noexcept
        {
            auto const word = static_cast<std::size_t>(i / 64);
            return word < words_.size() && ((words_[word] >> (i % 64)) & 1) != 0;
        }

        // Un bit à 1 strictement sous le rang i ?
        constexpr bool any_below(int i) const noexcept
        {
            if (i <= 0)
                return false;
            auto const full = static_cast<std::size_t>(i / 64);
            for (std::size_t w = 0; w < std::min(full, words_.size()); ++w)
                if (words_[w] != 0)
                    return true;
            if (full < words_.size() && i % 64 != 0)
                return (words_[full] & ((std::uint64_t{1} << (i % 64)) - 1)) != 0;
            return false;
        }

        // Les 64 bits de rang [i, i + 64).
        constexpr std::uint64_t bits_at(int i) const noexcept
        {
            auto const word = static_cast<std::size_t>(i / 64);
            int const shift = i % 64;
            std::uint64_t const low = word < words_.size() ? words_[word] >> shift : 0;
            std::uint64_t const high = (shift != 0 && word + 1 < words_.size()) ? words_[word + 1] << (64 - shift) : 0;
            return low | high;
        }

        constexpr big shifted_left(int n) const
        {
            if (is_zero() || n == 0)
                return *this;
            big r;
            auto const whole = static_cast<std::size_t>(n / 64);
            int const shift = n % 64;
            r.words_.assign(whole, 0);
            std::uint64_t carry = 0;
            for (std::uint64_t w : words_)
            {
                r.words_.push_back(shift == 0 ? w : (w << shift) | carry);
                carry = shift == 0 ? 0 : w >> (64 - shift);
            }
            if (carry != 0)
                r.words_.push_back(carry);
            return r;
        }

        constexpr big shifted_right(int n) const
        {
            big r;
            auto const whole = static_cast<std::size_t>(n / 64);
            if (whole >= words_.size())
                return r;
            int const shift = n % 64;
            for (std::size_t i = whole; i < words_.size(); ++i)
            {
                std::uint64_t w = words_[i] >> shift;
                if (shift != 0 && i + 1 < words_.size())
                    w |= words_[i + 1] << (64 - shift);
                r.words_.push_back(w);
            }
            r.trim();
            return r;
        }

        friend constexpr big operator+(big const& a, big const& b)
        {
            big r;
            std::size_t const n = std::max(a.words_.size(), b.words_.size());
            r.words_.resize(n + 1, 0);
            std::uint64_t carry = 0;
            for (std::size_t i = 0; i < n; ++i)
            {
                uint128 const s = uint128{i < a.words_.size() ? a.words_[i] : 0} + (i < b.words_.size() ? b.words_[i] : 0) + carry;
                r.words_[i] = static_cast<std::uint64_t>(s);
                carry = static_cast<std::uint64_t>(s >> 64);
            }
            r.words_[n] = carry;
            r.trim();
            return r;
        }

        // a − b, pour a >= b.
        friend constexpr big operator-(big const& a, big const& b)
        {
            big r;
            r.words_.resize(a.words_.size(), 0);
            std::uint64_t borrow = 0;
            for (std::size_t i = 0; i < a.words_.size(); ++i)
            {
                std::uint64_t const bi = i < b.words_.size() ? b.words_[i] : 0;
                std::uint64_t const d = a.words_[i] - bi - borrow;
                borrow = (a.words_[i] < bi || (a.words_[i] == bi && borrow != 0)) ? 1 : 0;
                r.words_[i] = d;
            }
            r.trim();
            return r;
        }

        friend constexpr big operator*(big const& a, big const& b)
        {
            big r;
            if (a.is_zero() || b.is_zero())
                return r;
            r.words_.assign(a.words_.size() + b.words_.size(), 0);
            for (std::size_t i = 0; i < a.words_.size(); ++i)
            {
                std::uint64_t carry = 0;
                for (std::size_t j = 0; j < b.words_.size(); ++j)
                {
                    uint128 const t = uint128{a.words_[i]} * b.words_[j] + r.words_[i + j] + carry;
                    r.words_[i + j] = static_cast<std::uint64_t>(t);
                    carry = static_cast<std::uint64_t>(t >> 64);
                }
                r.words_[i + b.words_.size()] = carry;
            }
            r.trim();
            return r;
        }

        friend constexpr std::strong_ordering operator<=>(big const& a, big const& b) noexcept
        {
            if (a.words_.size() != b.words_.size())
                return a.words_.size() <=> b.words_.size();
            for (std::size_t i = a.words_.size(); i-- > 0;)
                if (a.words_[i] != b.words_[i])
                    return a.words_[i] <=> b.words_[i];
            return std::strong_ordering::equal;
        }
        friend constexpr bool operator==(big const&, big const&) = default;

    private:
        constexpr void trim() noexcept
        {
            while (!words_.empty() && words_.back() == 0)
                words_.pop_back();
        }

        std::vector<std::uint64_t> words_;
    };

    // Division euclidienne bit à bit : quotient et reste (b > 0).
    struct big_division
    {
        big quotient, remainder;
    };

    constexpr big_division divide(big const& a, big const& b)
    {
        big_division d;
        for (int i = a.bit_length() - 1; i >= 0; --i)
        {
            d.remainder = d.remainder.shifted_left(1);
            if (a.bit(i))
                d.remainder = d.remainder + big(1);
            if (!(d.remainder < b))
            {
                d.remainder = d.remainder - b;
                d.quotient = d.quotient + big(1).shifted_left(i);
            }
        }
        return d;
    }

    // Division par un entier d'un mot (d > 0), mot par mot sur 128 bits.
    struct big_by_word
    {
        big quotient;
        std::uint64_t remainder;
    };

    constexpr big_by_word divide(big const& a, std::uint64_t d)
    {
        int const n = (a.bit_length() + 63) / 64;
        std::vector<std::uint64_t> digits(static_cast<std::size_t>(n), 0);  // poids fort d'abord
        uint128 r = 0;
        for (int i = n - 1; i >= 0; --i)
        {
            uint128 const cur = (r << 64) | a.bits_at(64 * i);
            digits[static_cast<std::size_t>(n - 1 - i)] = static_cast<std::uint64_t>(cur / d);
            r = cur % d;
        }
        return {big::from_digits(digits), static_cast<std::uint64_t>(r)};
    }

    // Grand entier signé (signe et valeur absolue), pour les bornes d'intervalles.
    struct signed_big
    {
        bool negative = false;
        big magnitude;

        friend constexpr signed_big operator+(signed_big const& a, signed_big const& b)
        {
            if (a.negative == b.negative)
                return {a.negative, a.magnitude + b.magnitude};
            if (!(a.magnitude < b.magnitude))
                return {a.negative && !(a.magnitude == b.magnitude), a.magnitude - b.magnitude};
            return {b.negative, b.magnitude - a.magnitude};
        }
        friend constexpr signed_big operator*(signed_big const& a, signed_big const& b)
        {
            big m = a.magnitude * b.magnitude;
            bool const negative = (a.negative != b.negative) && !m.is_zero();
            return {negative, std::move(m)};
        }
    };
}
