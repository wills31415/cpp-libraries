#pragma once

// Moteurs aléatoires : sources de bits uniformes, entièrement spécifiées, donc
// reproductibles au bit près sur toute plateforme. Tous satisfont
// std::uniform_random_bit_generator (utilisables avec <random> et avec
// wsl::random), et s'emploient en évaluation constante.
//
//   splitmix64     suite de Weyl brouillée ; sert surtout à initialiser les autres
//   xoshiro256pp   rapide, excellente qualité ; sauts de 2^128 et 2^192 (flux disjoints)
//   pcg64          PCG XSL RR 128/64 ; flux sélectionnables, avance en O(log n)
//   philox4x64     Philox 4×64-10, à compteur : la n-ième valeur s'obtient directement,
//                  flux parallèles indépendants par construction
//
// Les sorties sont vérifiées contre des sources indépendantes de ce code
// (NumPy pour pcg64 et philox4x64, sauts de xoshiro par matrices sur GF(2)) :
// tools/random/engine_vectors.py.

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>

namespace wsl::random
{
    namespace detail
    {
        // Entier non signé de 128 bits (extension de GCC et Clang).
        __extension__ typedef unsigned __int128 uint128;
    }

    // -----------------------------------------------------------------------
    // splitmix64 (Steele, Lea, Flood ; Vigna)
    // -----------------------------------------------------------------------
    class splitmix64
    {
    public:
        using result_type = std::uint64_t;

        static constexpr result_type min() noexcept { return 0; }
        static constexpr result_type max() noexcept { return ~result_type{0}; }

        constexpr explicit splitmix64(std::uint64_t seed = 0) noexcept : state_(seed) {}

        constexpr void seed(std::uint64_t seed) noexcept { state_ = seed; }

        constexpr result_type operator()() noexcept
        {
            std::uint64_t z = (state_ += 0x9e3779b97f4a7c15);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
            z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
            return z ^ (z >> 31);
        }

        // L'état est une suite de Weyl : avancer de n coûte une multiplication.
        constexpr void discard(unsigned long long n) noexcept { state_ += n * 0x9e3779b97f4a7c15; }

        constexpr std::uint64_t state() const noexcept { return state_; }

        friend constexpr bool operator==(splitmix64 const&, splitmix64 const&) = default;

    private:
        std::uint64_t state_;
    };

    // -----------------------------------------------------------------------
    // xoshiro256++ (Blackman, Vigna)
    // -----------------------------------------------------------------------
    class xoshiro256pp
    {
    public:
        using result_type = std::uint64_t;
        using state_type = std::array<std::uint64_t, 4>;

        static constexpr result_type min() noexcept { return 0; }
        static constexpr result_type max() noexcept { return ~result_type{0}; }

        // État initialisé par quatre sorties de splitmix64(seed), comme le
        // recommandent les auteurs : jamais nul.
        constexpr explicit xoshiro256pp(std::uint64_t seed = 0) noexcept { this->seed(seed); }

        // État brut : ne doit pas être entièrement nul.
        constexpr explicit xoshiro256pp(state_type const& state) noexcept : s_(state) {}

        constexpr void seed(std::uint64_t seed) noexcept
        {
            splitmix64 g(seed);
            for (auto& word : s_)
                word = g();
        }

        constexpr result_type operator()() noexcept
        {
            std::uint64_t const result = std::rotl(s_[0] + s_[3], 23) + s_[0];
            std::uint64_t const t = s_[1] << 17;
            s_[2] ^= s_[0];
            s_[3] ^= s_[1];
            s_[1] ^= s_[2];
            s_[0] ^= s_[3];
            s_[2] ^= t;
            s_[3] = std::rotl(s_[3], 45);
            return result;
        }

        constexpr void discard(unsigned long long n) noexcept
        {
            while (n-- != 0)
                (*this)();
        }

        // Avance de 2^128 pas : 2^128 flux disjoints de 2^128 valeurs chacun.
        constexpr void jump() noexcept
        {
            constexpr std::uint64_t table[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c};
            jump_by(table);
        }

        // Avance de 2^192 pas : 2^64 groupes de flux, à répartir entre machines.
        constexpr void long_jump() noexcept
        {
            constexpr std::uint64_t table[] = {0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635};
            jump_by(table);
        }

        constexpr state_type const& state() const noexcept { return s_; }

        friend constexpr bool operator==(xoshiro256pp const&, xoshiro256pp const&) = default;

    private:
        constexpr void jump_by(std::uint64_t const (&table)[4]) noexcept
        {
            state_type acc{};
            for (std::uint64_t word : table)
                for (int b = 0; b < 64; ++b)
                {
                    if ((word >> b) & 1)
                        for (int i = 0; i < 4; ++i)
                            acc[i] ^= s_[i];
                    (*this)();
                }
            s_ = acc;
        }

        state_type s_{};
    };

    // -----------------------------------------------------------------------
    // pcg64 : PCG XSL RR 128/64 (O'Neill), celui de pcg-cpp et de NumPy (PCG64)
    // -----------------------------------------------------------------------
    class pcg64
    {
    public:
        using result_type = std::uint64_t;

        static constexpr result_type min() noexcept { return 0; }
        static constexpr result_type max() noexcept { return ~result_type{0}; }

        // Flux par défaut : celui de pcg-cpp.
        static constexpr detail::uint128 default_stream =
            ((detail::uint128{0x5851f42d4c957f2d} << 64) | 0x14057b7ef767814f) >> 1;

        // Initialisation de pcg-cpp (setseq) : incrément 2 stream + 1 (impair), état
        // nul, un pas, + seed, un pas.
        constexpr explicit pcg64(detail::uint128 seed = 0, detail::uint128 stream = default_stream) noexcept
        {
            this->seed(seed, stream);
        }

        constexpr void seed(detail::uint128 seed, detail::uint128 stream = default_stream) noexcept
        {
            inc_ = (stream << 1) | 1;
            state_ = 0;
            step();
            state_ += seed;
            step();
        }

        // État et incrément bruts (l'incrément doit être impair).
        struct raw_state
        {
            detail::uint128 state, increment;
        };

        constexpr explicit pcg64(raw_state raw) noexcept : state_(raw.state), inc_(raw.increment | 1) {}

        constexpr result_type operator()() noexcept
        {
            step();
            auto const rotation = static_cast<int>(state_ >> 122);
            auto const folded = static_cast<std::uint64_t>(state_ ^ (state_ >> 64));
            return std::rotr(folded, rotation);
        }

        // Avance de delta pas en O(log delta) (Brown, 1994).
        constexpr void advance(detail::uint128 delta) noexcept
        {
            detail::uint128 multiplier = 1, increment = 0;
            detail::uint128 current_multiplier = lcg_multiplier, current_increment = inc_;
            while (delta != 0)
            {
                if (delta & 1)
                {
                    multiplier *= current_multiplier;
                    increment = increment * current_multiplier + current_increment;
                }
                current_increment = (current_multiplier + 1) * current_increment;
                current_multiplier *= current_multiplier;
                delta >>= 1;
            }
            state_ = multiplier * state_ + increment;
        }

        constexpr void discard(unsigned long long n) noexcept { advance(n); }

        constexpr raw_state state() const noexcept { return {state_, inc_}; }

        friend constexpr bool operator==(pcg64 const& a, pcg64 const& b) noexcept
        {
            return a.state_ == b.state_ && a.inc_ == b.inc_;
        }

    private:
        static constexpr detail::uint128 lcg_multiplier =
            (detail::uint128{0x2360ed051fc65da4} << 64) | 0x4385df649fccf645;

        constexpr void step() noexcept { state_ = state_ * lcg_multiplier + inc_; }

        detail::uint128 state_ = 0, inc_ = 1;
    };

    // -----------------------------------------------------------------------
    // philox4x64 : Philox 4×64, 10 tours (Salmon et al., Random123)
    //
    // Fonction de (clé, compteur) : bloc(c) = 4 valeurs ; le moteur rend bloc(c)
    // dans l'ordre, puis passe à c + 1. Deux clés différentes donnent deux flux
    // indépendants ; le compteur se fixe directement (accès aléatoire).
    // -----------------------------------------------------------------------
    class philox4x64
    {
    public:
        using result_type = std::uint64_t;
        using key_type = std::array<std::uint64_t, 2>;
        using counter_type = std::array<std::uint64_t, 4>;

        static constexpr result_type min() noexcept { return 0; }
        static constexpr result_type max() noexcept { return ~result_type{0}; }

        constexpr explicit philox4x64(std::uint64_t seed = 0) noexcept : key_{seed, 0} {}

        constexpr philox4x64(key_type key, counter_type counter) noexcept : key_(key), counter_(counter) {}

        constexpr void seed(std::uint64_t seed) noexcept { *this = philox4x64(seed); }

        constexpr result_type operator()() noexcept
        {
            if (index_ == 4)
            {
                buffer_ = block(key_, counter_);
                increment(counter_);
                index_ = 0;
            }
            return buffer_[static_cast<std::size_t>(index_++)];
        }

        // Le bloc de clé et de compteur donnés, sans état.
        static constexpr counter_type block(key_type key, counter_type counter) noexcept
        {
            for (int round = 0; round < 10; ++round)
            {
                if (round != 0)
                {
                    key[0] += 0x9e3779b97f4a7c15;
                    key[1] += 0xbb67ae8584caa73b;
                }
                auto const p0 = detail::uint128{0xd2e7470ee14c6c93} * counter[0];
                auto const p1 = detail::uint128{0xca5a826395121157} * counter[2];
                auto const hi0 = static_cast<std::uint64_t>(p0 >> 64), lo0 = static_cast<std::uint64_t>(p0);
                auto const hi1 = static_cast<std::uint64_t>(p1 >> 64), lo1 = static_cast<std::uint64_t>(p1);
                counter = {hi1 ^ counter[1] ^ key[0], lo1, hi0 ^ counter[3] ^ key[1], lo0};
            }
            return counter;
        }

        // Avance de n valeurs : n / 4 blocs sautés d'un coup.
        constexpr void discard(unsigned long long n) noexcept
        {
            while (n != 0 && index_ != 4)
            {
                (*this)();
                --n;
            }
            for (unsigned long long blocks = n / 4; blocks != 0; --blocks)
                increment(counter_);
            for (n %= 4; n != 0; --n)
                (*this)();
        }

        constexpr key_type const& key() const noexcept { return key_; }
        // Compteur du prochain bloc à calculer.
        constexpr counter_type const& counter() const noexcept { return counter_; }

        friend constexpr bool operator==(philox4x64 const& a, philox4x64 const& b) noexcept
        {
            // Même position dans le même flux : même clé, et même prochaine valeur.
            auto next_counter = [](philox4x64 const& e) {
                counter_type c = e.counter_;
                if (e.index_ != 4)                                         // bloc en cours : c'est c − 1
                    decrement(c);
                return c;
            };
            return a.key_ == b.key_ && next_counter(a) == next_counter(b) && (a.index_ % 4) == (b.index_ % 4);
        }

    private:
        static constexpr void increment(counter_type& c) noexcept
        {
            for (auto& word : c)
                if (++word != 0)
                    break;
        }

        static constexpr void decrement(counter_type& c) noexcept
        {
            for (auto& word : c)
                if (word-- != 0)
                    break;
        }

        key_type key_{};
        counter_type counter_{};
        counter_type buffer_{};
        int index_ = 4;                                                     // 4 : tampon vide
    };
}
