#pragma once

// Quadrature des lois : Gauss–Legendre à 12 nœuds (nœuds et poids en double-double,
// calculés à la compilation), sur [a, b] et adaptative (intégrandes à l'échelle).

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/trigonometry.hpp>
#include <wsl/special.hpp>

namespace wsl::random::detail
{
    // Gauss–Legendre à 12 nœuds sur [−1, 1] (6 nœuds positifs et leurs poids), en
    // double-double : Newton sur P₁₂ (récurrence de Bonnet), à la compilation.
    inline constexpr int legendre_order = 12;

    struct legendre_rule
    {
        sd::double_double node[legendre_order / 2], weight[legendre_order / 2];
    };

    constexpr legendre_rule make_legendre_rule() noexcept
    {
        legendre_rule r{};
        constexpr int m = legendre_order;
        for (int i = 0; i < m / 2; ++i)
        {
            double const guess = wsl::special::cos(pi_hi * (static_cast<double>(i) + 0.75) / (static_cast<double>(m) + 0.5));
            sd::double_double x{guess, 0.0}, derivative{1.0, 0.0};
            for (int iteration = 0; iteration < 8; ++iteration)
            {
                sd::double_double p0{1.0, 0.0}, p1 = x;                 // P_{k−1}, P_k
                for (int k = 1; k < m; ++k)
                {
                    auto const next = sd::quotient_dd(sum(mul({static_cast<double>(2 * k + 1), 0.0}, mul(x, p1)), negated(mul({static_cast<double>(k), 0.0}, p0))),
                                                      {static_cast<double>(k + 1), 0.0});
                    p0 = p1;
                    p1 = next;
                }
                auto const one_minus_x2 = sum({1.0, 0.0}, negated(mul(x, x)));
                derivative = sd::quotient_dd(mul({static_cast<double>(m), 0.0}, sum(p0, negated(mul(x, p1)))), one_minus_x2);   // P′ = m (P_{m−1} − x P_m)/(1 − x²)
                x = sum(x, negated(sd::quotient_dd(p1, derivative)));
            }
            r.node[i] = x;
            auto const one_minus_x2 = sum({1.0, 0.0}, negated(mul(x, x)));
            r.weight[i] = sd::quotient_dd({2.0, 0.0}, mul(one_minus_x2, mul(derivative, derivative)));
        }
        return r;
    }

    inline constexpr legendre_rule legendre = make_legendre_rule();

    // ∫_a^b, Gauss–Legendre à 12 nœuds sur [a, b] (à l'échelle).
    template <class F>
    constexpr scaled legendre12(F const& f, sd::double_double a, sd::double_double b) noexcept
    {
        auto const width = sum(b, negated(a));
        sd::double_double const half{sd::fmul(width.hi, 0.5), sd::fmul(width.lo, 0.5)};
        if (!(half.hi > 0))
            return {0.0, 0.0, 0};
        auto const middle = sum(a, half);
        scaled total{0.0, 0.0, 0};
        for (int i = 0; i < legendre_order / 2; ++i)
        {
            auto const offset = mul(half, legendre.node[i]);
            scaled const w = normalized(legendre.weight[i].hi, legendre.weight[i].lo);
            scaled const v1 = f(sum(middle, offset)), v2 = f(sum(middle, negated(offset)));
            if (v1.hi != 0)
                total = add_scaled(total, v1 * w);
            if (v2.hi != 0)
                total = add_scaled(total, v2 * w);
        }
        return total.hi == 0 ? total : total * normalized(half.hi, half.lo);
    }

    // |a − b| <= 2^-k max(|a|, |b|) (à l'échelle, positifs) ?
    constexpr bool scaled_close(scaled const& a, scaled const& b, scaled const& reference, int k) noexcept
    {
        if (a.hi == 0 && b.hi == 0)
            return true;
        scaled const& big = a.hi == 0 ? b : b.hi == 0 ? a : (a.e > b.e || (a.e == b.e && a.hi >= b.hi) ? a : b);
        scaled const& small = &big == &a ? b : a;
        if (small.hi == 0)
            return reference.hi != 0 && big.e < reference.e - k;
        scaled const d = scaled_difference(big, small);
        if (d.hi == 0)
            return true;
        scaled const& ref = reference.hi == 0 ? big : reference;
        return d.e < ref.e - k;
    }

    // ∫_a^b f adaptatif : un intervalle est accepté quand Q(a, b) et Q(a, m) + Q(m, b) concordent à
    // 2^-82 près de max(leur valeur, reference) (reference : la masse déjà intégrée ailleurs).
    template <class F>
    constexpr scaled adaptive_integral(F const& f, sd::double_double a, sd::double_double b, scaled const& reference) noexcept
    {
        struct interval
        {
            sd::double_double a, b;
            scaled q;
            int depth;
        };
        interval stack[128]{};
        int top = 0;
        stack[top++] = {a, b, legendre12(f, a, b), 0};
        scaled total{0.0, 0.0, 0};
        int evaluations = 0;
        while (top > 0)
        {
            interval const it = stack[--top];
            auto const width = sum(it.b, negated(it.a));
            auto const m = sum(it.a, {sd::fmul(width.hi, 0.5), sd::fmul(width.lo, 0.5)});
            scaled const left = legendre12(f, it.a, m), right = legendre12(f, m, it.b);
            scaled const both = add_scaled(left, right);
            ++evaluations;
            scaled const ref = add_scaled(reference, add_scaled(total, both));
            if (scaled_close(it.q, both, ref, 82) || it.depth >= 50 || top + 2 > 128 || evaluations > 2000)
            {
                total = add_scaled(total, both);
                continue;
            }
            stack[top++] = {m, it.b, right, it.depth + 1};
            stack[top++] = {it.a, m, left, it.depth + 1};
        }
        return total;
    }

    // ∫ par morceaux depuis start vers end (forward : end > start) : avant chaque morceau,
    // stop(t, total, reste) peut arrêter la marche ; sinon le morceau [t, t ± width(t)], tronqué
    // à end, est intégré par piece(bas, haut, total). Intégrandes et résultat à l'échelle.
    template <class Width, class Piece, class Stop>
    constexpr scaled march_integral(sd::double_double start, sd::double_double end, bool forward, Width const& width,
                                    Piece const& piece, Stop const& stop, int max_pieces = 1000000) noexcept
    {
        scaled total{0.0, 0.0, 0};
        sd::double_double t = start;
        for (int k = 0; k < max_pieces; ++k)
        {
            auto const rest = forward ? sum(end, negated(t)) : sum(t, negated(end));
            if (!(rest.hi > 0) || stop(t, total, rest))
                break;
            double const w = width(t);
            bool const last = !(w < rest.hi);
            sd::double_double const h = last ? rest : sd::double_double{w, 0.0};
            auto const next = last ? end : forward ? sum(t, h) : sum(t, negated(h));
            scaled const p = forward ? piece(t, next, total) : piece(next, t, total);
            if (p.hi != 0)
                total = add_scaled(total, p);
            if (last)
                break;
            t = next;
        }
        return total;
    }
}
