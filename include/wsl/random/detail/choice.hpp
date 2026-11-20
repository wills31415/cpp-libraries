#pragma once

// Choix exacts à poids flottants : un flottant est un entier à une puissance de 2
// près (w = m 2^e), deux poids le sont donc à une même puissance près ; un entier
// uniforme sous leur somme décide exactement. Sur 128 bits quand les entiers y
// tiennent, en grands entiers sinon.
//
//   choose(g, a, b)     vrai avec probabilité a/(a + b) exactement, a, b >= 0, a + b > 0 ;
//   fraction(g, a, b)   vrai avec probabilité a/b exactement, 0 <= a <= b, b > 0.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/big.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <bit>
#include <cstdint>
#include <vector>

namespace wsl::random::detail
{
    // Entier uniforme sur [0, n), n > 0 en grand entier : k = bit_length(n) bits
    // tirés (mot de tête d'abord), rejetés s'ils dépassent n (moins de 2 essais en moyenne).
    template <class G>
    constexpr big uniform_below_big(G& g, big const& n)
    {
        int const k = n.bit_length();
        int const words = (k + 63) / 64;
        int const top = k - 64 * (words - 1);
        std::vector<std::uint64_t> digits(static_cast<std::size_t>(words));
        for (;;)
        {
            for (auto& d : digits)
                d = random_bits(g);
            if (top < 64)
                digits[0] &= (std::uint64_t(1) << top) - 1;
            big r = big::from_digits(digits);
            if (r < n)
                return r;
        }
    }

    // Entier uniforme sur [0, n), n > 0 sur 128 bits : Lemire sous 2^64, rejet sinon.
    template <class G>
    constexpr sd::uint128 uniform_below_128(G& g, sd::uint128 n)
    {
        if ((n >> 64) == 0)
            return uniform_below(static_cast<std::uint64_t>(n), g);
        int const k = 64 + std::bit_width(static_cast<std::uint64_t>(n >> 64));
        std::uint64_t const mask = k - 64 == 64 ? ~std::uint64_t{0} : (std::uint64_t{1} << (k - 64)) - 1;
        for (;;)
        {
            std::uint64_t const hi = random_bits(g) & mask;
            std::uint64_t const lo = random_bits(g);
            sd::uint128 const r = (sd::uint128{hi} << 64) | lo;
            if (r < n)
                return r;
        }
    }

    // w = m 2^e, m impair (w > 0 fini).
    struct dyadic
    {
        std::uint64_t m;
        int e;
    };
    constexpr dyadic dyadic_of(double w) noexcept
    {
        std::uint64_t const b = std::bit_cast<std::uint64_t>(w);
        int const biased = static_cast<int>(b >> 52);
        std::uint64_t m = biased == 0 ? (b & 0x000fffffffffffff) : ((b & 0x000fffffffffffff) | 0x0010000000000000);
        int e = biased == 0 ? -1074 : biased - 1075;
        int const z = std::countr_zero(m);
        return {m >> z, e + z};
    }

    // A = m_a 2^{e_a − e}, B = m_b 2^{e_b − e}, e = min ; u < A, u uniforme sous total(A, B).
    template <bool Fraction, class G>
    constexpr bool integer_choice(G& g, double a, double b)
    {
        auto const da = dyadic_of(a), db = dyadic_of(b);
        int const e = da.e < db.e ? da.e : db.e;
        int const wa = std::bit_width(da.m) + (da.e - e), wb = std::bit_width(db.m) + (db.e - e);
        if (wa <= 126 && wb <= 126)
        {
            sd::uint128 const ia = sd::uint128{da.m} << (da.e - e), ib = sd::uint128{db.m} << (db.e - e);
            return uniform_below_128(g, Fraction ? ib : ia + ib) < ia;
        }
        big const ia = big(da.m).shifted_left(da.e - e), ib = big(db.m).shifted_left(db.e - e);
        return uniform_below_big(g, Fraction ? ib : ia + ib) < ia;
    }

    // Vrai avec probabilité a/(a + b), a, b >= 0 finis, a + b > 0.
    template <class G>
    constexpr bool choose(G& g, double a, double b)
    {
        if (a == 0)
            return false;
        if (b == 0)
            return true;
        return integer_choice<false>(g, a, b);
    }

    // Vrai avec probabilité a/b, 0 <= a <= b, b > 0 finis.
    template <class G>
    constexpr bool fraction(G& g, double a, double b)
    {
        if (a == 0)
            return false;
        if (a >= b)
            return true;
        return integer_choice<true>(g, a, b);
    }

    // Indice i avec probabilité wᵢ/Σ w, poids entiers (grands entiers), total = Σ w > 0.
    template <class G>
    constexpr std::size_t pick(G& g, std::vector<big> const& weights, big const& total)
    {
        big z = uniform_below_big(g, total);
        for (std::size_t i = 0; i + 1 < weights.size(); ++i)
        {
            if (z < weights[i])
                return i;
            z = z - weights[i];
        }
        return weights.size() - 1;
    }

    // w = M 2^e exactement, M entier en grand entier, pour un exposant e <= celui de w.
    constexpr big integer_at(double w, int e)
    {
        if (w == 0)
            return big();
        auto const d = dyadic_of(w);
        return big(d.m).shifted_left(d.e - e);
    }

    // Poids entiers qui croissent (urnes, restaurants) : sur 64 bits tant que le total
    // reste sous 2^63, en grands entiers ensuite (le passage ne dépend que des poids :
    // tirages reproductibles ; exacts dans les deux modes).
    class growing_weights
    {
    public:
        constexpr growing_weights() = default;
        constexpr explicit growing_weights(std::vector<big> const& w)
        {
            for (auto const& v : w)
                push(v);
        }

        constexpr std::size_t size() const noexcept { return small_ ? w64_.size() : big_.size(); }

        constexpr void push(big const& v)
        {
            if (small_ && fits(v))
            {
                w64_.push_back(v.is_zero() ? 0 : v.bits_at(0));
                total64_ += w64_.back();
                if (total64_ < (std::uint64_t{1} << 63))
                    return;
                to_big();
                return;
            }
            if (small_)
                to_big();
            big_.push_back(v);
            total_ = total_ + v;
        }

        // w_i += v
        constexpr void add(std::size_t i, big const& v)
        {
            if (small_ && fits(v))
            {
                std::uint64_t const x = v.is_zero() ? 0 : v.bits_at(0);
                if (total64_ + x < (std::uint64_t{1} << 63))
                {
                    w64_[i] += x;
                    total64_ += x;
                    return;
                }
            }
            if (small_)
                to_big();
            big_[i] = big_[i] + v;
            total_ = total_ + v;
        }

        // Indice i avec probabilité wᵢ/Σ w.
        template <class G>
        constexpr std::size_t pick(G& g) const
        {
            if (!small_)
                return wsl::random::detail::pick(g, big_, total_);
            std::uint64_t z = uniform_below(total64_, g);
            for (std::size_t i = 0; i + 1 < w64_.size(); ++i)
            {
                if (z < w64_[i])
                    return i;
                z -= w64_[i];
            }
            return w64_.size() - 1;
        }

    private:
        static constexpr bool fits(big const& v) noexcept { return v.bit_length() <= 62; }

        constexpr void to_big()
        {
            for (std::uint64_t x : w64_)
                big_.push_back(big(x));
            total_ = big(total64_);
            w64_.clear();
            small_ = false;
        }

        bool small_ = true;
        std::vector<std::uint64_t> w64_;
        std::uint64_t total64_ = 0;
        std::vector<big> big_;
        big total_;
    };
}
