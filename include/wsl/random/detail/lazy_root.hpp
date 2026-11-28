#pragma once

#include <wsl/modes.hpp>

// Tirage exact de μ ± √v · V (normale donnée par sa variance, √v irrationnel) : arrondi
// correct sans jamais calculer √v. Pour un point dyadique x, le signe de μ ± √A − x
// (A = v t², t rationnel) se décide exactement en comparant A au carré de x − μ, en
// entiers. Un candidat approché est corrigé jusqu'à ce que les deux milieux qui
// l'encadrent soient de part et d'autre de la valeur : c'est l'arrondi au plus proche
// (jamais d'égalité : la valeur est irrationnelle, ou vaut μ si t = 0).

#include <wsl/random/detail/portable_math.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/big.hpp>
#include <wsl/random/detail/lazy.hpp>

#include <algorithm>
#include <bit>
#include <compare>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace wsl::inline WSL_DEFAULTS::random::detail
{
    // m · 2^e, exact.
    struct exact_dyadic
    {
        signed_big m;
        int e = 0;
    };

    constexpr exact_dyadic exact_of(double x)
    {
        auto const d = decompose(x);
        return {{d.negative, big(d.mantissa)}, d.exponent};
    }

    constexpr exact_dyadic exact_sum(exact_dyadic a, exact_dyadic b)
    {
        if (a.m.magnitude.is_zero())
            return b;
        if (b.m.magnitude.is_zero())
            return a;
        int const base = std::min(a.e, b.e);
        a.m.magnitude = a.m.magnitude.shifted_left(a.e - base);
        b.m.magnitude = b.m.magnitude.shifted_left(b.e - base);
        return {a.m + b.m, base};
    }

    constexpr exact_dyadic exact_negated(exact_dyadic a)
    {
        if (!a.m.magnitude.is_zero())
            a.m.negative = !a.m.negative;
        return a;
    }

    // v = (hi + lo) · 2^e > 0 (à l'échelle) en entier : v = R · 2^e.
    struct exact_square
    {
        big r;
        int e = 0;
    };

    constexpr exact_square exact_square_of(scaled const& v)
    {
        exact_dyadic const t = exact_sum(exact_of(v.hi), exact_of(v.lo));
        return {t.m.magnitude, t.e + v.e};
    }

    // Signe de P + √v B (v = R 2^ev > 0, P et B dyadiques) : immédiat si P et B ne sont pas de
    // signes opposés, sinon par comparaison de P² et v B².
    constexpr int sign_plus_root(exact_dyadic const& p, exact_square const& v, exact_dyadic const& b)
    {
        bool const p_zero = p.m.magnitude.is_zero(), b_zero = b.m.magnitude.is_zero();
        int const sp = p_zero ? 0 : p.m.negative ? -1 : 1, sb = b_zero ? 0 : b.m.negative ? -1 : 1;
        if (sb == 0 || sp == sb)
            return sp;
        if (sp == 0)
            return sb;
        big const p2 = p.m.magnitude * p.m.magnitude, vb2 = v.r * (b.m.magnitude * b.m.magnitude);
        int const ep = 2 * p.e, eb = v.e + 2 * b.e;
        int const base = std::min(ep, eb);
        auto const c = p2.shifted_left(ep - base) <=> vb2.shifted_left(eb - base);   // P² face à v B²
        int const s = c > 0 ? 1 : c < 0 ? -1 : 0;
        return sp > 0 ? s : -s;
    }

    // Cellule d'arrondi de r (T fini) : ses deux milieux, dyadiques exacts.
    template <class T>
    constexpr void rounding_cell(T r, exact_dyadic& below, exact_dyadic& above)
    {
        using F = format<T>;
        double const x = r;
        double const magnitude = x < 0 ? -x : x;
        int ulp_e = F::denorm_exponent;
        bool power = false;
        if (magnitude != 0)
        {
            int const e = wsl::random::detail::sd::exponent_of(magnitude);
            if (e >= F::min_exponent)
            {
                ulp_e = e - F::digits + 1;
                power = magnitude == wsl::random::detail::sd::power_of_two(e) && e > F::min_exponent;
            }
        }
        exact_dyadic const up{{false, big(1)}, ulp_e - 1}, down{{false, big(1)}, power ? ulp_e - 2 : ulp_e - 1};
        exact_dyadic const c = exact_of(x);
        if (x < 0)                                                          // vers zéro : l'écart « down » est au-dessus
        {
            above = exact_sum(c, down);
            below = exact_sum(c, exact_negated(up));
        }
        else
        {
            above = exact_sum(c, up);
            below = exact_sum(c, exact_negated(down));
        }
    }

    // Arrondi au plus proche de y, connu par side(x) = signe de y − x (x dyadique ; y n'est
    // jamais un milieu) : pas d'un T à son voisin depuis le candidat, puis bissection sur
    // l'ordre des T si le candidat était loin (compensation).
    template <class T, class Side>
    constexpr T round_by_side(Side const& side, T candidate)
    {
        using bits_type = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;
        constexpr bits_type sign_bit = bits_type{1} << (8 * sizeof(T) - 1);
        // clé croissante avec la valeur (−∞ … +∞), et sa réciproque
        auto const key = [](T x) {
            bits_type const b = std::bit_cast<bits_type>(x);
            return (b & sign_bit) != 0 ? bits_type(~b) : bits_type(b | sign_bit);
        };
        auto const of_key = [](bits_type k) { return std::bit_cast<T>((k & sign_bit) != 0 ? bits_type(k & ~sign_bit) : bits_type(~k)); };
        constexpr T inf = std::numeric_limits<T>::infinity();
        // y au-dessus de la cellule de r ? (cellule de +∞ : au-delà du milieu haut du plus grand fini)
        auto const above_cell = [&](T r) {
            if (r == inf)
                return false;
            exact_dyadic below, above;
            rounding_cell(r == -inf ? -std::numeric_limits<T>::max() : r, below, above);
            return side(r == -inf ? below : above) > 0;
        };
        if (candidate != candidate)
            candidate = 0;
        T y = candidate == 0 ? T(0) : candidate;
        for (int step = 0; step < 64; ++step)
        {
            if (above_cell(y))
            {
                y = of_key(bits_type(key(y) + 1));
                continue;
            }
            T const previous = y == -inf ? y : of_key(bits_type(key(y) - 1));
            if (y != -inf && above_cell(previous) == false && !(previous == T(0) && y == T(0)))
            {
                y = previous;
                continue;
            }
            return y == 0 ? (side(exact_dyadic{}) < 0 ? -T(0) : T(0)) : y;   // zéro du signe de la valeur
        }
        // bissection : plus petit r (en clé) dont la cellule n'est pas sous y
        bits_type lo = key(-inf), hi = key(inf);
        while (lo < hi)
        {
            bits_type const mid = lo + (hi - lo) / 2;
            if (above_cell(of_key(mid)))
                lo = mid + 1;
            else
                hi = mid;
        }
        T const r = of_key(lo);
        return r == 0 ? (side(exact_dyadic{}) < 0 ? -T(0) : T(0)) : r;
    }

    // Arrondi correct en T de μ ± √v · w 2^-b (w >= 0 entier), σ ≈ √v pour le candidat.
    template <class T>
    constexpr T round_root_affine(double mu, exact_square const& v, double sigma, bool negative, big const& w, int b)
    {
        if (w.is_zero())
            return static_cast<T>(mu);
        // y − x = (μ − x) ± √v · w 2^-b
        exact_dyadic const coefficient{{negative, w}, -b};
        exact_dyadic const m = exact_of(mu);
        auto const side = [&](exact_dyadic const& x) { return sign_plus_root(exact_sum(m, exact_negated(x)), v, coefficient); };
        // Candidat : t ≈ w 2^-b par ses 64 bits de tête.
        int const length = w.bit_length();
        int const shift = std::max(0, length - 64);
        double const t = wsl::random::detail::sd::scale_by_power_of_two(static_cast<double>(w.shifted_right(shift).bits_at(0)), shift - b);
        double const step = wsl::random::detail::sd::fmul(sigma, t);
        return round_by_side<T>(side, static_cast<T>(negative ? mu - step : mu + step));
    }

    // Arrondi correct en T de μ ± √v · V, V = k + u paresseux : raffine u jusqu'à ce que les
    // deux bornes de son intervalle s'arrondissent au même T (même consommation que round_lazy ;
    // la voie rapide ne conclut que si le résultat est certain, donc mêmes bits et même consommation).
    template <class T, class G>
    constexpr T round_lazy_root(double mu, scaled const& variance, double sigma, bool negative, lazy_real& v, G& g)
    {
        if (v.fraction.size() == 0)
            v.fraction.refine(g);
        T quick{};                                                          // voie rapide : seulement si certaine
        auto const root = unscaled(sqrt(variance));
        if (round_quickly(mu, root.hi, 1.0, negative, v, quick, root.lo))
            return quick;
        exact_square const s = exact_square_of(variance);
        using bits_type = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;
        for (;;)
        {
            int const b = v.fraction.bits();
            big const low = big(v.integer).shifted_left(b) + v.fraction.numerator();
            T const a = round_root_affine<T>(mu, s, sigma, negative, low, b);
            T const c = round_root_affine<T>(mu, s, sigma, negative, low + big(1), b);
            if (std::bit_cast<bits_type>(a) == std::bit_cast<bits_type>(c))
                return a;
            v.fraction.refine(g);
        }
    }
}
