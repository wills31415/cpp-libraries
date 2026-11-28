#pragma once

#include <wsl/modes.hpp>

// Primitives uniformes exactes : la fondation de wsl::random.
//
//   random_bits(g)              64 bits uniformes, de n'importe quel générateur
//   uniform_below(n, g)         entier uniforme sur [0, n)        (Lemire)
//   uniform_int(a, b, g)        entier uniforme sur [a, b]
//   uniform_real<T>(g)          réel uniforme sur [0, 1), exact : chaque flottant x
//                               est rendu avec exactement la probabilité de [x, x⁺),
//                               sous-normaux compris
//   uniform_real_fast<T>(g)     réel uniforme sur [0, 1), rapide : 2^53 (2^24) valeurs
//                               également espacées
//   bernoulli_trial(p, g)             vrai avec probabilité exactement p (p flottant)
//   bernoulli_trial(num, den, g)      vrai avec probabilité exactement num/den
//   shuffle(range, g)           mélange de Fisher–Yates : les n! ordres équiprobables
//
// « Exact » : si le générateur rend des bits uniformes et indépendants, les
// probabilités sont exactement celles annoncées, sans biais d'arrondi ni de
// modulo. Algorithmes fixés et documentés, indépendants de la bibliothèque
// standard : mêmes résultats au bit près partout. Le nombre de valeurs consommées
// sur le générateur fait aussi partie du contrat.

#include <wsl/random/engines.hpp>

#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <random>
#include <ranges>
#include <stdexcept>
#include <type_traits>

namespace wsl::inline WSL_DEFAULTS::random
{
    // -----------------------------------------------------------------------
    // 64 bits uniformes. Générateur de plage 2^64 : une valeur. Plage 2^w :
    // valeurs concaténées, la première en poids fort (les bits de tête en excès
    // sont écartés). Autre plage (minstd_rand…) : rejet jusqu'à une valeur sous
    // 2^w, w = ⌊log2(plage)⌋, puis concaténation.
    // -----------------------------------------------------------------------
    template <std::uniform_random_bit_generator G>
    constexpr std::uint64_t random_bits(G& g)
    {
        constexpr auto low = static_cast<std::uint64_t>(G::min());
        constexpr std::uint64_t span = static_cast<std::uint64_t>(G::max()) - low;   // plage − 1
        static_assert(span != 0, "wsl::random::random_bits : générateur à une seule valeur");
        if constexpr (span == ~std::uint64_t{0})
            return static_cast<std::uint64_t>(g()) - low;
        else
        {
            constexpr bool power_of_two = (span & (span + 1)) == 0;
            constexpr int width = std::bit_width(span + 1) - 1;              // ⌊log2(plage)⌋
            std::uint64_t result = 0;
            for (int filled = 0; filled < 64;)
            {
                std::uint64_t const v = static_cast<std::uint64_t>(g()) - low;
                if constexpr (!power_of_two)
                    if (v >> width)                                         // hors de [0, 2^w) : rejet
                        continue;
                result = (result << width) | v;
                filled += width;
            }
            return result;
        }
    }

    // -----------------------------------------------------------------------
    // Entier uniforme sur [0, n), n > 0 : méthode de Lemire (2019). Produit
    // 128 bits m = x n ; sa partie haute est uniforme sur [0, n) une fois rejetées
    // les valeurs dont la partie basse tombe sous 2^64 mod n. Une seule valeur du
    // générateur presque toujours, et une division seulement dans le cas rare.
    // -----------------------------------------------------------------------
    template <std::uniform_random_bit_generator G>
    constexpr std::uint64_t uniform_below(std::uint64_t n, G& g)
    {
        if (n == 0)
            throw std::domain_error("wsl::random::uniform_below : n doit être > 0");
        detail::uint128 m = detail::uint128{random_bits(g)} * n;
        auto low = static_cast<std::uint64_t>(m);
        if (low < n)
        {
            std::uint64_t const threshold = (0 - n) % n;                     // 2^64 mod n
            while (low < threshold)
            {
                m = detail::uint128{random_bits(g)} * n;
                low = static_cast<std::uint64_t>(m);
            }
        }
        return static_cast<std::uint64_t>(m >> 64);
    }

    // Entier uniforme sur [a, b] (types entiers jusqu'à 64 bits, sauf bool).
    template <std::integral T, std::uniform_random_bit_generator G>
        requires(!std::same_as<T, bool> && sizeof(T) <= 8)
    constexpr T uniform_int(T a, T b, G& g)
    {
        if (b < a)
            throw std::domain_error("wsl::random::uniform_int : a doit être <= b");
        using U = std::make_unsigned_t<T>;
        auto const span = static_cast<std::uint64_t>(static_cast<U>(static_cast<U>(b) - static_cast<U>(a)));
        std::uint64_t const offset = span == ~std::uint64_t{0} ? random_bits(g) : uniform_below(span + 1, g);
        return static_cast<T>(static_cast<U>(static_cast<U>(a) + static_cast<U>(offset)));
    }

    // -----------------------------------------------------------------------
    // Réel uniforme exact sur [0, 1).
    //
    // U = 0,b₁b₂b₃… (bits indépendants) ; si le premier 1 est en position k,
    // U ∈ [2^-k, 2^(1-k)) et les bits suivants donnent la mantisse. Le résultat
    // est U arrondi vers le bas : chaque flottant x l'est avec probabilité
    // x⁺ − x exactement, et 1 n'est jamais atteint. k suit une loi géométrique,
    // lue sur les zéros de tête ; la mantisse vient des bits restants du même mot
    // quand il en reste assez (presque toujours : une valeur consommée), d'un mot
    // neuf sinon.
    // -----------------------------------------------------------------------
    namespace detail
    {
    template <std::floating_point T, std::uniform_random_bit_generator G>
        requires(std::same_as<T, float> || std::same_as<T, double>)
    constexpr T uniform_real_exact(G& g)
    {
        constexpr int mantissa = std::numeric_limits<T>::digits - 1;        // 52 ou 23
        constexpr int last_normal = -std::numeric_limits<T>::min_exponent + 1;   // 1022 ou 126 : k maximal normal
        constexpr int last = last_normal + mantissa;                         // 1074 ou 149 : plus petit sous-normal

        int k = 0;
        std::uint64_t word = random_bits(g);
        while (word == 0)
        {
            k += 64;
            if (k >= last)                                                   // U < 2^-last : 0
                return T(0);
            word = random_bits(g);
        }
        int const zeros = std::countl_zero(word);
        k += zeros + 1;
        if (k > last)
            return T(0);
        std::uint64_t const fraction = zeros <= 63 - mantissa ? word << (zeros + 1) : random_bits(g);

        using Bits = std::conditional_t<std::same_as<T, double>, std::uint64_t, std::uint32_t>;
        if (k <= last_normal)
        {
            auto const biased = static_cast<Bits>(last_normal + 1 - k);       // exposant −k, biaisé
            auto const m = static_cast<Bits>(fraction >> (64 - mantissa));
            return std::bit_cast<T>(static_cast<Bits>((biased << mantissa) | m));
        }
        // Sous-normal : last − k bits après le premier 1, valeur N · 2^-last.
        int const kept = last - k;                                           // dans [0, mantissa − 1]
        std::uint64_t const n = (std::uint64_t{1} << kept) | (kept == 0 ? 0 : fraction >> (64 - kept));
        return std::bit_cast<T>(static_cast<Bits>(n));
    }

    // Réel uniforme approché sur [0, 1) : les 53 (24) bits de tête d'une valeur, sur une grille
    // de pas 2^-53 (2^-24).
    template <std::floating_point T, std::uniform_random_bit_generator G>
        requires(std::same_as<T, float> || std::same_as<T, double>)
    constexpr T uniform_real_approximate(G& g)
    {
        if constexpr (std::same_as<T, double>)
            return static_cast<double>(random_bits(g) >> 11) * 0x1p-53;
        else
            return static_cast<float>(random_bits(g) >> 40) * 0x1p-24f;
    }
    }

    // Réel uniforme sur [0, 1), dans le mode des étiquettes (design/interfaces.md, décisions 24
    // à 27) : approché par défaut (une valeur, grille de pas 2^-53), exact sur wsl::exact (chaque
    // flottant x avec probabilité x⁺ − x, presque toujours une valeur). Les deux sont portables.
    template <std::floating_point T = double, std::uniform_random_bit_generator G, class... M>
        requires(std::same_as<T, float> || std::same_as<T, double>) && mode_set<M...>
    constexpr T uniform_real(G& g, M...)
    {
        if constexpr (std::same_as<algorithm_of<M...>, exact_t>)
            return detail::uniform_real_exact<T>(g);
        else
            return detail::uniform_real_approximate<T>(g);
    }

    // -----------------------------------------------------------------------
    // Bernoulli exacte : vrai si U < p, U uniforme sur [0, 1), en comparant les
    // bits de U à ceux de p (un flottant : développement binaire fini) mot de 64
    // bits par mot de 64 bits. Le premier mot qui diffère décide ; si p n'a plus
    // de bits, U >= p. Une valeur consommée presque toujours (2^-64 de chances
    // d'en lire une autre).
    // -----------------------------------------------------------------------
    template <std::uniform_random_bit_generator G>
    constexpr bool bernoulli_trial(double p, G& g)
    {
        if (!(p == p) || p < 0 || p > 1)
            throw std::domain_error("wsl::random::bernoulli_trial : p doit être dans [0, 1]");
        if (p == 0)
            return false;
        if (p == 1)
            return true;
        std::uint64_t const b = std::bit_cast<std::uint64_t>(p);
        auto const biased = static_cast<int>(b >> 52);
        std::uint64_t const m = biased == 0 ? (b & 0x000fffffffffffff) : ((b & 0x000fffffffffffff) | 0x0010000000000000);
        int const exponent = biased == 0 ? -1074 : biased - 1075;            // p = m 2^exponent
        // Mot j de p : ⌊p 2^(64(j+1))⌋ mod 2^64 = ⌊m 2^s⌋ mod 2^64, s = exponent + 64(j+1).
        for (int j = 0;; ++j)
        {
            int const s = exponent + 64 * (j + 1);
            std::uint64_t const p_word = s >= 64 ? 0 : (s >= 0 ? m << s : (s > -64 ? m >> -s : 0));
            std::uint64_t const u_word = random_bits(g);
            if (u_word != p_word)
                return u_word < p_word;
            // Plus aucun bit de p au-delà de ce mot : U >= p.
            if (s >= 0 || (s > -64 && (m & ((std::uint64_t{1} << -s) - 1)) == 0))
                return false;
        }
    }

    // Bernoulli rationnelle : vrai avec probabilité num/den exactement (0 <= num <= den, den > 0).
    template <std::uniform_random_bit_generator G>
    constexpr bool bernoulli_trial(std::uint64_t num, std::uint64_t den, G& g)
    {
        if (den == 0 || num > den)
            throw std::domain_error("wsl::random::bernoulli_trial : il faut 0 <= num <= den, den > 0");
        return uniform_below(den, g) < num;
    }

    // -----------------------------------------------------------------------
    // Mélange de Fisher–Yates : pour i de n − 1 à 1, échange de l'élément i avec
    // un élément uniforme de [0, i]. Les n! ordres sont équiprobables ; n − 1
    // tirages uniform_below.
    // -----------------------------------------------------------------------
    template <std::ranges::random_access_range R, std::uniform_random_bit_generator G>
        requires std::permutable<std::ranges::iterator_t<R>> && std::ranges::sized_range<R>
    constexpr void shuffle(R&& range, G& g)
    {
        auto const first = std::ranges::begin(range);
        auto const n = static_cast<std::uint64_t>(std::ranges::size(range));
        for (std::uint64_t i = n; i > 1; --i)
        {
            auto const j = uniform_below(i, g);
            std::ranges::iter_swap(first + static_cast<std::ptrdiff_t>(i - 1), first + static_cast<std::ptrdiff_t>(j));
        }
    }
}
