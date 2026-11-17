#pragma once

// Réels paresseux et tirages exacts (exponentielle, normale).
//
// Un réel uniforme paresseux U dans [0, 1) n'existe que par les mots de 64 bits
// de son développement binaire déjà tirés : U ∈ [A, A + 1] · 2^-64n. On n'en
// tire de nouveaux que lorsqu'une comparaison ou un arrondi l'exige ; deux réels
// uniformes indépendants sont égaux avec probabilité nulle, si bien que toute
// comparaison se décide en un nombre fini (presque toujours un) de mots.
//
// Exponentielle exacte (von Neumann, 1951) : x uniforme est accepté si, dans la
// suite x > U₁ > U₂ > …, la première remontée survient à un rang impair, ce qui
// arrive avec probabilité e^-x exactement ; chaque rejet ajoute 1 à la partie
// entière, qui suit ainsi la loi géométrique de raison e^-1. Résultat : k + x.
//
// Demi-normale exacte par rejet : la loi exponentielle E₁ est acceptée si une
// autre exponentielle E₂ dépasse (E₁ − 1)²/2, ce qui arrive avec probabilité
// e^{−(E₁−1)²/2} : le rapport des densités demi-normale et exponentielle. La
// comparaison se fait par arithmétique d'intervalles exacte (grands entiers) sur
// les bits tirés, raffinés jusqu'à décision.
//
// Arrondi : μ + σ V / Q au plus proche (pair), en double ou en float, V réel
// paresseux, en raffinant V jusqu'à ce que les deux bornes de son intervalle
// donnent le même flottant : c'est l'arrondi correct d'une valeur tirée
// exactement selon la loi. Une voie rapide en double-double décide presque
// toujours ; la voie exacte, en grands entiers, tranche les cas serrés.

#include <wsl/random/detail/big.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special/detail/fp.hpp>

#include <bit>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace wsl::detail
{
    // Le premier mot est rangé sur place (sans allocation) : presque tous les
    // réels paresseux s'en contentent.
    class lazy_uniform
    {
    public:
        template <class G>
        constexpr void refine(G& g)
        {
            std::uint64_t const w = random_bits(g);
            if (count_ == 0)
                first_ = w;
            else
                rest_.push_back(w);
            ++count_;
        }
        constexpr int size() const noexcept { return count_; }
        constexpr int bits() const noexcept { return 64 * count_; }
        constexpr std::uint64_t word(int i) const noexcept
        {
            return i == 0 ? first_ : rest_[static_cast<std::size_t>(i - 1)];
        }
        constexpr big numerator() const
        {
            std::vector<std::uint64_t> digits;
            digits.reserve(static_cast<std::size_t>(count_));
            for (int i = 0; i < count_; ++i)
                digits.push_back(word(i));
            return big::from_digits(digits);
        }

    private:
        std::uint64_t first_ = 0;
        std::vector<std::uint64_t> rest_;
        int count_ = 0;
    };

    // U uniforme dans (0, 1) : uniform_real, 0 rejeté (probabilité 2^-1074…).
    template <class G>
    constexpr double uniform_positive(G& g)
    {
        for (;;)
            if (double const u = uniform_real<double>(g); u != 0)
                return u;
    }

    // a < b, en tirant les mots nécessaires (égalité : probabilité nulle).
    template <class G>
    constexpr bool less(lazy_uniform& a, lazy_uniform& b, G& g)
    {
        for (int i = 0;; ++i)
        {
            if (i == a.size())
                a.refine(g);
            if (i == b.size())
                b.refine(g);
            if (a.word(i) != b.word(i))
                return a.word(i) < b.word(i);
        }
    }

    // k + u, k entier, u uniforme paresseux.
    struct lazy_real
    {
        std::uint64_t integer = 0;
        lazy_uniform fraction;
    };

    template <class G>
    constexpr lazy_real exact_exponential(G& g)
    {
        for (std::uint64_t k = 0;; ++k)
        {
            lazy_uniform x;
            lazy_uniform previous;
            if (!less(previous, x, g))                                      // U₁ > x : remontée au rang 1
                return {k, std::move(x)};
            int n = 2;
            for (;; ++n)
            {
                lazy_uniform next;
                if (!less(next, previous, g))                               // remontée au rang n
                    break;
                previous = std::move(next);
            }
            if (n % 2 == 1)
                return {k, std::move(x)};
        }
    }

    // E₂ > (E₁ − 1)²/2 ? Bornes exactes de chaque côté sur les bits tirés ; on
    // raffine le plus imprécis des deux jusqu'à séparation.
    template <class G>
    constexpr bool exceeds_shifted_half_square(lazy_real& e2, lazy_real& e1, G& g)
    {
        if (e1.fraction.size() == 0)
            e1.fraction.refine(g);
        if (e2.fraction.size() == 0)
            e2.fraction.refine(g);
        // Raccourci : avec un mot chacun, bornes en double à 2^-50 près
        // relativement au plus ; une séparation de plus de 2^-40 (relativement,
        // plus 2^-60 en absolu) décide sûrement. Sinon, le calcul exact.
        if (e1.fraction.size() == 1 && e2.fraction.size() == 1 && e1.integer < (1u << 20) && e2.integer < (1u << 20))
        {
            namespace sd = wsl::special::detail;                             // produits protégés de la contraction :
            double const u1 = sd::fmul(static_cast<double>(e1.fraction.word(0) >> 11), 0x1p-53);   // la décision, donc
            double const u2 = sd::fmul(static_cast<double>(e2.fraction.word(0) >> 11), 0x1p-53);   // le nombre de mots tirés,
            // en dépend (u₁ par défaut, à 2^-53 près)
            double const t_low = static_cast<double>(e1.integer) - 1.0 + u1;
            double const t_high = t_low + 0x1p-52;
            double const square_low = e1.integer >= 1 ? sd::fmul(t_low, t_low) : sd::fmul(t_high, t_high);
            double const square_high = e1.integer >= 1 ? sd::fmul(t_high, t_high) : sd::fmul(t_low, t_low);
            double const g_low = sd::fmul(0.5, square_low), g_high = sd::fmul(0.5, square_high);
            double const e_low = static_cast<double>(e2.integer) + u2;
            double const e_high = e_low + 0x1p-52;
            double const margin = 0x1p-40;
            if (e_low > sd::fmul(g_high, 1 + margin) + 0x1p-60)
                return true;
            if (sd::fmul(e_high, 1 + margin) + 0x1p-60 < g_low)
                return false;
        }
        for (;;)
        {
            int const b1 = e1.fraction.bits(), b2 = e2.fraction.bits();
            big const a1 = e1.fraction.numerator(), a2 = e2.fraction.numerator();
            // t = E₁ − 1 dans [t_low, t_high] · 2^-b1 ; (t²/2) en unités de 2^-(2 b1 + 1).
            big square_low, square_high;
            if (e1.integer >= 1)
            {
                big const t_low = big(e1.integer - 1).shifted_left(b1) + a1;
                big const t_high = t_low + big(1);
                square_low = t_low * t_low;
                square_high = t_high * t_high;
            }
            else                                                            // t = u₁ − 1 dans [−1, 0] : |t| décroît
            {
                big const one = big(1).shifted_left(b1);
                big const abs_high = one - a1;                              // |t| <= (2^b1 − A₁) · 2^-b1
                big const abs_low = abs_high - big(1);
                square_low = abs_low * abs_low;
                square_high = abs_high * abs_high;
            }
            big const e2_low = big(e2.integer).shifted_left(b2) + a2;
            big const e2_high = e2_low + big(1);
            // E₂ > t²/2  ⇔  e2 · 2^(2 b1 + 1) > square · 2^b2
            if (e2_low.shifted_left(2 * b1 + 1) > square_high.shifted_left(b2))
                return true;
            if (e2_high.shifted_left(2 * b1 + 1) < square_low.shifted_left(b2))
                return false;
            if (b1 <= b2)
                e1.fraction.refine(g);
            else
                e2.fraction.refine(g);
        }
    }

    // |Z| exacte, Z normale centrée réduite.
    template <class G>
    constexpr lazy_real exact_half_normal(G& g)
    {
        for (;;)
        {
            lazy_real e1 = exact_exponential(g);
            lazy_real e2 = exact_exponential(g);
            if (exceeds_shifted_half_square(e2, e1, g))
                return e1;
        }
    }

    // ---------------------------------------------------------------------
    // Arrondi au plus proche (pair) d'une valeur dyadique ou rationnelle, en
    // double ou en float.
    // ---------------------------------------------------------------------
    struct decomposed
    {
        bool negative;
        std::uint64_t mantissa;                                             // x = ± mantissa · 2^exponent
        int exponent;
    };

    // Mantisse impaire (zéros de queue retirés) : une puissance de 2 a pour
    // mantisse 1, ce qui évite toute division par elle.
    constexpr decomposed decompose(double x) noexcept
    {
        std::uint64_t const b = std::bit_cast<std::uint64_t>(x);
        auto const biased = static_cast<int>((b >> 52) & 0x7ff);
        std::uint64_t m = b & 0x000fffffffffffff;
        int e = -1074;
        if (biased != 0)
        {
            m |= 0x0010000000000000;
            e = biased - 1075;
        }
        if (m != 0)
        {
            int const zeros = std::countr_zero(m);
            m >>= zeros;
            e += zeros;
        }
        return {(b >> 63) != 0, m, e};
    }

    // Format cible : précision, plus petit exposant normal, débordement.
    template <class T>
    struct format
    {
        static constexpr int digits = std::numeric_limits<T>::digits;      // 53 ou 24
        static constexpr int min_exponent = std::numeric_limits<T>::min_exponent - 1;   // −1022 ou −126
        static constexpr int max_exponent = std::numeric_limits<T>::max_exponent;       // 1024 ou 128
        static constexpr int denorm_exponent = min_exponent - digits + 1;  // −1074 ou −149
    };

    // ± (q + sticky·ε) · 2^g, q entier, au plus proche (pair) en T ; sticky : la
    // valeur dépasse strictement q · 2^g. Quand q tient sans arrondi, l'appelant
    // garantit sticky faux (assez de bits de garde).
    template <class T>
    constexpr T round_to(bool negative, big const& q, bool sticky, int g)
    {
        using F = format<T>;
        if (q.is_zero())
            return negative ? -T(0) : T(0);
        int const length = q.bit_length();
        int const top = length - 1 + g;                                     // exposant du bit de tête
        if (top >= F::max_exponent)
            return negative ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity();
        int const drop = top >= F::min_exponent ? length - F::digits : F::denorm_exponent - g;
        std::uint64_t mantissa = 0;
        int scale = g;
        if (drop <= 0)
            mantissa = q.bits_at(0);
        else
        {
            if (drop > length + 1)                                          // sous la moitié du plus petit sous-normal
                return negative ? -T(0) : T(0);
            mantissa = q.shifted_right(drop).bits_at(0);
            bool const round_bit = q.bit(drop - 1);
            bool const rest = sticky || q.any_below(drop - 1);
            if (round_bit && (rest || (mantissa & 1)))
                ++mantissa;
            scale = g + drop;
        }
        // mantissa · 2^scale est représentable en T (ou déborde vers ±∞) : exact en
        // double, puis converti sans arrondi.
        auto const magnitude = static_cast<T>(wsl::special::detail::scale_by_power_of_two(static_cast<double>(mantissa), scale));
        return negative ? -magnitude : magnitude;
    }

    // Arrondi de μ + σ · s · W / (2^b · Q), s = ±1, W entier : une borne de
    // l'intervalle du réel paresseux. σ > 0, Q > 0 finis.
    template <class T>
    constexpr T round_affine(double mu, double sigma, double divisor, bool v_negative, big const& w, int b)
    {
        decomposed const m = decompose(mu), s = decompose(sigma), d = decompose(divisor);
        // y = [mμ mQ 2^(eμ + b + eQ) ± mσ W 2^eσ] / (mQ 2^(b + eQ)).
        int const e1 = m.exponent + b + d.exponent, e2 = s.exponent;
        int const base = m.mantissa == 0 ? e2 : std::min(e1, e2);
        signed_big const t1{m.negative, (big(m.mantissa) * big(d.mantissa)).shifted_left(m.mantissa == 0 ? 0 : e1 - base)};
        signed_big const t2{v_negative, (big(s.mantissa) * w).shifted_left(e2 - base)};
        signed_big const n = t1 + t2;
        int const g = base - b - d.exponent;                                // y = (n / mQ) · 2^g
        if (d.mantissa == 1)
            return round_to<T>(n.negative, n.magnitude, false, g);
        // Quotient par mQ avec au moins 56 bits de tête, reste en « sticky ».
        int const extra = std::max(0, 56 + 64 - n.magnitude.bit_length());
        auto const [q, r] = divide(n.magnitude.shifted_left(extra), d.mantissa);
        return round_to<T>(n.negative, q, r != 0, g - extra);
    }

    // Voie rapide de l'arrondi, quand le réel paresseux n'a qu'un mot : valeur du
    // milieu de l'intervalle en double-double, demi-largeur et erreur de calcul
    // majorées ; si tout l'intervalle tombe strictement dans la cellule d'arrondi
    // d'un même T, c'est lui. Sinon (rare), rien : la voie exacte décide, avec les
    // mêmes bits. Le résultat et le nombre de mots tirés sont donc ceux de la voie
    // exacte, quelle que soit la voie prise.
    template <class T>
    constexpr bool round_quickly(double mu, double sigma, double divisor, bool negative, lazy_real const& v, T& result) noexcept
    {
        namespace sd = wsl::special::detail;
        using F = format<T>;
        if (v.fraction.size() != 1 || v.integer >= (std::uint64_t{1} << 20))
            return false;
        auto const in_range = [](double x) { return x >= 0x1p-800 && x <= 0x1p800; };
        if (!in_range(sigma) || !in_range(divisor) || (mu != 0 && !(mu >= -0x1p800 && mu <= 0x1p800)))
            return false;

        // V = k + c · 2^-65, c = 2A + 1 (milieu), c = c₁ 2^12 + c₀ : somme exacte, puis reste.
        sd::uint128 const c = (sd::uint128{v.fraction.word(0)} << 1) | 1;
        double const c1 = sd::fmul(static_cast<double>(static_cast<std::uint64_t>(c >> 12)), 0x1p-53);
        double const c0 = sd::fmul(static_cast<double>(static_cast<std::uint64_t>(c & 0xfff)), 0x1p-65);
        auto const [v_hi, v_error] = sd::two_sum(static_cast<double>(v.integer), c1);
        double const v_lo = v_error + c0;

        // t = σ V / Q en double-double.
        auto const [p, p_error] = sd::two_product(sigma, v_hi);
        auto const [t_hi, t_lo] = sd::fast_two_sum(p, p_error + sd::fmul(sigma, v_lo));
        double q_hi = t_hi, q_lo = t_lo;
        if (divisor != 1.0)
        {
            q_hi = t_hi / divisor;
            if (!in_range(q_hi < 0 ? -q_hi : q_hi))
                return false;
            auto const [m, m_error] = sd::two_product(q_hi, divisor);
            q_lo = (((t_hi - m) - m_error) + t_lo) / divisor;
        }
        if (negative)
        {
            q_hi = -q_hi;
            q_lo = -q_lo;
        }

        // y = μ + q.
        auto const [s, s_error] = sd::two_sum(mu, q_hi);
        auto const [y_hi, y_lo] = sd::fast_two_sum(s, s_error + q_lo);
        double const magnitude_bound = (mu < 0 ? -mu : mu) + (q_hi < 0 ? -q_hi : q_hi);
        double const half_width = sd::fmul(sd::fmul(sigma, 0x1p-65) / divisor, 1 + 0x1p-20);
        double const bound = half_width + sd::fmul(magnitude_bound, 0x1p-95);

        // Cellule d'arrondi de r = T(y) : [r − bas/2, r + haut/2], à l'écart d'une marge.
        auto const r = static_cast<T>(y_hi);
        double const a = r < 0 ? -double(r) : double(r);
        if (!(a >= std::numeric_limits<T>::min() * 2 && a < std::numeric_limits<T>::max() / 2))
            return false;
        int const e = sd::exponent_of(a);
        double const up = sd::power_of_two(e - F::digits + 1);
        double const down = a == sd::power_of_two(e) ? up * 0.5 : up;
        double deviation = (y_hi - double(r)) + y_lo;                       // y_hi − r exact (Sterbenz)
        if (r < 0)
            deviation = -deviation;                                         // écart vers l'extérieur
        double const margin = 1 - 0x1p-30;
        if (deviation + bound < sd::fmul(up * 0.5, margin) && bound - deviation < sd::fmul(down * 0.5, margin))
        {
            result = r;
            return true;
        }
        return false;
    }

    // Arrondi correct (en T) de μ + σ · V / Q, V = ±(k + u) paresseux : raffine u
    // jusqu'à ce que les deux bornes de son intervalle s'arrondissent au même T.
    template <class T, class G>
    constexpr T round_lazy(double mu, double sigma, double divisor, bool negative, lazy_real& v, G& g)
    {
        if (v.fraction.size() == 0)                                          // sans aucun bit, l'arrondi
            v.fraction.refine(g);                                            // n'est presque jamais décidé
        T quick{};
        if (round_quickly(mu, sigma, divisor, negative, v, quick))
            return quick;
        for (;;)
        {
            int const b = v.fraction.bits();
            big const low = big(v.integer).shifted_left(b) + v.fraction.numerator();
            T const a = round_affine<T>(mu, sigma, divisor, negative, low, b);
            T const c = round_affine<T>(mu, sigma, divisor, negative, low + big(1), b);
            using bits_type = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;
            if (std::bit_cast<bits_type>(a) == std::bit_cast<bits_type>(c))
                return a;
            v.fraction.refine(g);
        }
    }
}
