#pragma once

#include <wsl/modes.hpp>

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
// Moteurs complets au sens de RandomNumberEngine : E(), E(s), E(q) et seed(q) depuis une
// SeedSequence (mots de 32 bits assemblés deux à deux, le premier en poids faible, comme
// mt19937_64), seed(), discard(n), ==, os << e et is >> e (entiers décimaux séparés par
// des espaces ; une lecture invalide laisse le moteur intact et lève failbit).
//
// Les sorties sont vérifiées contre des sources indépendantes de ce code
// (NumPy pour pcg64 et philox4x64, sauts de xoshiro par matrices sur GF(2)) :
// tools/random/engine_vectors.py.

#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <type_traits>

namespace wsl::inline WSL_DEFAULTS::random::engines
{
    namespace detail
    {
        // Entier non signé de 128 bits (extension de GCC et Clang).
        __extension__ typedef unsigned __int128 uint128;

        // SeedSequence (seed_seq…) pour le moteur E : generate(début, fin) en mots de 32 bits ;
        // ni E ni un entier (E(E&) reste la copie, E(s) l'initialisation par une graine).
        template <class Q, class E>
        concept seed_sequence_for = !std::is_convertible_v<Q, typename E::result_type> && !std::same_as<std::remove_cv_t<Q>, E>
                                    && requires(Q& q, std::uint32_t* p) { q.generate(p, p); };

        // N mots de 64 bits tirés de q : 2N mots de 32 bits, le premier de chaque paire en poids faible.
        template <std::size_t N, class Q>
        constexpr std::array<std::uint64_t, N> seed_words(Q& q)
        {
            std::array<std::uint32_t, 2 * N> w{};
            q.generate(w.begin(), w.end());
            std::array<std::uint64_t, N> out{};
            for (std::size_t i = 0; i < N; ++i)
                out[i] = std::uint64_t{w[2 * i]} | (std::uint64_t{w[2 * i + 1]} << 32);
            return out;
        }

        // Écriture : entiers décimaux séparés par une espace, format du flux rétabli ensuite.
        template <class C, class Tr, std::size_t N>
        std::basic_ostream<C, Tr>& write_words(std::basic_ostream<C, Tr>& os, std::array<std::uint64_t, N> const& w)
        {
            using ios = std::basic_ios<C, Tr>;
            auto const flags = os.flags();
            auto const fill = os.fill();
            os.flags(ios::dec | ios::left);
            os.fill(os.widen(' '));
            for (std::size_t i = 0; i < N; ++i)
            {
                if (i != 0)
                    os << os.widen(' ');
                os << w[i];
            }
            os.flags(flags);
            os.fill(fill);
            return os;
        }

        // Lecture : rend faux (et lève failbit) si le flux échoue ou si valid refuse les mots.
        template <class C, class Tr, std::size_t N, class Valid>
        bool read_words(std::basic_istream<C, Tr>& is, std::array<std::uint64_t, N>& w, Valid const& valid)
        {
            using ios = std::basic_ios<C, Tr>;
            auto const flags = is.flags();
            is.flags(ios::dec | ios::skipws);
            for (auto& word : w)
                is >> word;
            is.flags(flags);
            if (is.fail())
                return false;
            if (!valid(w))
            {
                is.setstate(ios::failbit);
                return false;
            }
            return true;
        }

        constexpr uint128 join(std::uint64_t hi, std::uint64_t lo) noexcept { return (uint128{hi} << 64) | lo; }
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

        constexpr splitmix64() noexcept : splitmix64(0) {}
        constexpr explicit splitmix64(std::uint64_t seed) noexcept : state_(seed) {}
        template <detail::seed_sequence_for<splitmix64> Q>
        constexpr explicit splitmix64(Q& q) { seed(q); }

        constexpr void seed() noexcept { seed(0); }
        constexpr void seed(std::uint64_t seed) noexcept { state_ = seed; }
        template <detail::seed_sequence_for<splitmix64> Q>
        constexpr void seed(Q& q) { state_ = detail::seed_words<1>(q)[0]; }

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

        template <class C, class Tr>
        friend std::basic_ostream<C, Tr>& operator<<(std::basic_ostream<C, Tr>& os, splitmix64 const& e)
        {
            return detail::write_words(os, std::array<std::uint64_t, 1>{e.state_});
        }
        template <class C, class Tr>
        friend std::basic_istream<C, Tr>& operator>>(std::basic_istream<C, Tr>& is, splitmix64& e)
        {
            std::array<std::uint64_t, 1> w{};
            if (detail::read_words(is, w, [](auto const&) { return true; }))
                e.state_ = w[0];
            return is;
        }

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
        constexpr xoshiro256pp() noexcept : xoshiro256pp(0) {}
        constexpr explicit xoshiro256pp(std::uint64_t seed) noexcept { this->seed(seed); }
        template <detail::seed_sequence_for<xoshiro256pp> Q>
        constexpr explicit xoshiro256pp(Q& q) { seed(q); }

        // État brut : ne doit pas être entièrement nul.
        constexpr explicit xoshiro256pp(state_type const& state) noexcept : s_(state) {}

        constexpr void seed() noexcept { seed(0); }
        constexpr void seed(std::uint64_t seed) noexcept
        {
            splitmix64 g(seed);
            for (auto& word : s_)
                word = g();
        }
        // Huit mots de 32 bits ; un état entièrement nul (interdit) devient celui de seed(0).
        template <detail::seed_sequence_for<xoshiro256pp> Q>
        constexpr void seed(Q& q)
        {
            s_ = detail::seed_words<4>(q);
            if (s_ == state_type{})
                seed(0);
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

        template <class C, class Tr>
        friend std::basic_ostream<C, Tr>& operator<<(std::basic_ostream<C, Tr>& os, xoshiro256pp const& e)
        {
            return detail::write_words(os, e.s_);
        }
        template <class C, class Tr>
        friend std::basic_istream<C, Tr>& operator>>(std::basic_istream<C, Tr>& is, xoshiro256pp& e)
        {
            state_type w{};
            if (detail::read_words(is, w, [](state_type const& v) { return v != state_type{}; }))
                e.s_ = w;
            return is;
        }

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
        constexpr pcg64() noexcept : pcg64(0) {}
        constexpr explicit pcg64(detail::uint128 seed, detail::uint128 stream = default_stream) noexcept
        {
            this->seed(seed, stream);
        }
        template <detail::seed_sequence_for<pcg64> Q>
        constexpr explicit pcg64(Q& q) { seed(q); }

        constexpr void seed() noexcept { seed(0); }
        // Huit mots de 32 bits : graine (128 bits) puis flux (128 bits), poids faible d'abord.
        template <detail::seed_sequence_for<pcg64> Q>
        constexpr void seed(Q& q)
        {
            auto const w = detail::seed_words<4>(q);
            seed(detail::join(w[1], w[0]), detail::join(w[3], w[2]));
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

        // État (poids fort, poids faible) puis incrément (poids fort, poids faible) ; incrément impair.
        template <class C, class Tr>
        friend std::basic_ostream<C, Tr>& operator<<(std::basic_ostream<C, Tr>& os, pcg64 const& e)
        {
            auto const hi = [](detail::uint128 v) { return static_cast<std::uint64_t>(v >> 64); };
            auto const lo = [](detail::uint128 v) { return static_cast<std::uint64_t>(v); };
            return detail::write_words(os, std::array<std::uint64_t, 4>{hi(e.state_), lo(e.state_), hi(e.inc_), lo(e.inc_)});
        }
        template <class C, class Tr>
        friend std::basic_istream<C, Tr>& operator>>(std::basic_istream<C, Tr>& is, pcg64& e)
        {
            std::array<std::uint64_t, 4> w{};
            if (detail::read_words(is, w, [](auto const& v) { return (v[3] & 1) != 0; }))
            {
                e.state_ = detail::join(w[0], w[1]);
                e.inc_ = detail::join(w[2], w[3]);
            }
            return is;
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

        constexpr philox4x64() noexcept : philox4x64(0) {}
        constexpr explicit philox4x64(std::uint64_t seed) noexcept : key_{seed, 0} {}
        template <detail::seed_sequence_for<philox4x64> Q>
        constexpr explicit philox4x64(Q& q) { seed(q); }

        constexpr philox4x64(key_type key, counter_type counter) noexcept : key_(key), counter_(counter) {}

        constexpr void seed() noexcept { seed(0); }
        constexpr void seed(std::uint64_t seed) noexcept { *this = philox4x64(seed); }
        // Quatre mots de 32 bits : la clé ; compteur nul.
        template <detail::seed_sequence_for<philox4x64> Q>
        constexpr void seed(Q& q)
        {
            auto const w = detail::seed_words<2>(q);
            *this = philox4x64(key_type{w[0], w[1]}, counter_type{});
        }

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

        // Avance de n valeurs : n / 4 blocs sautés d'un coup (compteur + n / 4, O(1)).
        constexpr void discard(unsigned long long n) noexcept
        {
            while (n != 0 && index_ != 4)
            {
                (*this)();
                --n;
            }
            add(counter_, n / 4);
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

        // Clé (2 mots), compteur du prochain bloc (4 mots), position dans le bloc en cours (0 à 4,
        // 4 : tampon vide) ; le tampon se recalcule.
        template <class C, class Tr>
        friend std::basic_ostream<C, Tr>& operator<<(std::basic_ostream<C, Tr>& os, philox4x64 const& e)
        {
            auto const& k = e.key_;
            auto const& c = e.counter_;
            return detail::write_words(os, std::array<std::uint64_t, 7>{k[0], k[1], c[0], c[1], c[2], c[3], std::uint64_t(e.index_)});
        }
        template <class C, class Tr>
        friend std::basic_istream<C, Tr>& operator>>(std::basic_istream<C, Tr>& is, philox4x64& e)
        {
            std::array<std::uint64_t, 7> w{};
            if (detail::read_words(is, w, [](auto const& v) { return v[6] <= 4; }))
            {
                philox4x64 r(key_type{w[0], w[1]}, counter_type{w[2], w[3], w[4], w[5]});
                r.index_ = static_cast<int>(w[6]);
                if (r.index_ != 4)
                {
                    counter_type previous = r.counter_;
                    decrement(previous);
                    r.buffer_ = block(r.key_, previous);
                }
                e = r;
            }
            return is;
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

        // c + n, retenue propagée sur les quatre mots (modulo 2^256).
        static constexpr void add(counter_type& c, std::uint64_t n) noexcept
        {
            for (auto& word : c)
            {
                word += n;
                if (word >= n)                                              // pas de retenue
                    break;
                n = 1;
            }
        }

        key_type key_{};
        counter_type counter_{};
        counter_type buffer_{};
        int index_ = 4;                                                     // 4 : tampon vide
    };
}

// Transition (refonte, étape 5a) : les outils internes des moteurs restent visibles comme
// wsl::random::detail, où les emploient les primitives et les tirages exacts.
namespace wsl::inline WSL_DEFAULTS::random::detail
{
    using namespace engines::detail;
}
