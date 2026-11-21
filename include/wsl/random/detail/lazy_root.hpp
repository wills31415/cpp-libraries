#pragma once

// Tirage exact de μ ± √v · V (normale donnée par sa variance, √v irrationnel) : arrondi
// correct sans jamais calculer √v. Pour un point dyadique x, le signe de μ ± √A − x
// (A = v t², t rationnel) se décide exactement en comparant A au carré de x − μ, en
// entiers. Un candidat approché est corrigé jusqu'à ce que les deux milieux qui
// l'encadrent soient de part et d'autre de la valeur : c'est l'arrondi au plus proche
// (jamais d'égalité : la valeur est irrationnelle, ou vaut μ si t = 0).

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/big.hpp>
#include <wsl/random/detail/lazy.hpp>

#include <algorithm>
#include <bit>
#include <compare>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace wsl::random::detail
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

    // signe de √A − d (A = R 2^ea >= 0) : comparaison de A et d² si d >= 0.
    constexpr int root_minus_sign(big const& r, int ea, exact_dyadic const& d)
    {
        if (d.m.magnitude.is_zero())
            return r.is_zero() ? 0 : 1;
        if (d.m.negative)
            return 1;
        big const d2 = d.m.magnitude * d.m.magnitude;
        int const ed = 2 * d.e;
        int const base = std::min(ea, ed);
        auto const c = r.shifted_left(ea - base) <=> d2.shifted_left(ed - base);
        return c > 0 ? 1 : c < 0 ? -1 : 0;
    }

    // Arrondi correct en T de μ ± √v · w 2^-b (w >= 0 entier), σ ≈ √v pour le candidat.
    template <class T>
    constexpr T round_root_affine(double mu, exact_square const& v, double sigma, bool negative, big const& w, int b)
    {
        using F = format<T>;
        using bits_type = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;
        if (w.is_zero())
            return static_cast<T>(mu);
        big const r = v.r * (w * w);                                        // A = v w² 2^-2b
        int const ea = v.e - 2 * b;
        exact_dyadic const m = exact_negated(exact_of(mu));
        // signe de y − x, y = μ ± √A.
        auto const side = [&](exact_dyadic const& x) {
            exact_dyadic const d = exact_sum(x, m);                         // x − μ
            return negative ? -root_minus_sign(r, ea, exact_negated(d)) : root_minus_sign(r, ea, d);
        };
        // Candidat : t ≈ w 2^-b par ses 64 bits de tête.
        int const length = w.bit_length();
        int const shift = std::max(0, length - 64);
        double const t = wsl::special::detail::scale_by_power_of_two(static_cast<double>(w.shifted_right(shift).bits_at(0)), shift - b);
        double const step = wsl::special::detail::fmul(sigma, t);
        T y = static_cast<T>(negative ? mu - step : mu + step);
        // Demi-écarts aux voisins de y (dyadiques exacts) : milieux de la cellule d'arrondi.
        auto const half_gaps = [](T a, exact_dyadic& below, exact_dyadic& above) {
            double const x = a;
            double const magnitude = x < 0 ? -x : x;
            int ulp_e = F::denorm_exponent;
            bool power = false;
            if (magnitude != 0)
            {
                int const e = wsl::special::detail::exponent_of(magnitude);
                if (e >= F::min_exponent)
                {
                    ulp_e = e - F::digits + 1;
                    power = magnitude == wsl::special::detail::power_of_two(e) && e > F::min_exponent;
                }
            }
            exact_dyadic const up{{false, big(1)}, ulp_e - 1}, down{{false, big(1)}, power ? ulp_e - 2 : ulp_e - 1};
            exact_dyadic const c = exact_of(x);
            if (x < 0)                                                      // vers zéro : l'écart « down » est au-dessus
            {
                above = exact_sum(c, down);
                below = exact_sum(c, exact_negated(up));
            }
            else
            {
                above = exact_sum(c, up);
                below = exact_sum(c, exact_negated(down));
            }
        };
        auto const next = [](T a, bool up) {
            if (a == 0)
                return up ? std::numeric_limits<T>::denorm_min() : -std::numeric_limits<T>::denorm_min();
            bits_type const bits = std::bit_cast<bits_type>(a);
            return std::bit_cast<T>((a > 0) == up ? bits_type(bits + 1) : bits_type(bits - 1));
        };
        for (int guard = 0; guard < 4096; ++guard)
        {
            if (y == std::numeric_limits<T>::infinity() || y == -std::numeric_limits<T>::infinity())
            {
                T const edge = y > 0 ? std::numeric_limits<T>::max() : -std::numeric_limits<T>::max();
                exact_dyadic below, above;
                half_gaps(edge, below, above);
                if (y > 0 ? side(above) < 0 : side(below) > 0)               // dans la cellule du plus grand fini
                    y = edge;
                else
                    return y;
                continue;
            }
            exact_dyadic below, above;
            half_gaps(y, below, above);
            if (side(above) > 0)
                y = next(y, true);
            else if (side(below) < 0)
                y = next(y, false);
            else
                return y == 0 ? (side(exact_dyadic{}) < 0 ? -T(0) : T(0)) : y;   // zéro du signe de la valeur
        }
        return y;
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
