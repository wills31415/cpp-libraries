#pragma once

// Ordre fixé des sommes et produits flottants (docs/mdarray.md, I9) : par paires,
// comme NumPy (pairwise_sum) — en dessous de 8 termes, dans l'ordre ; jusqu'à 128,
// 8 accumulateurs entrelacés combinés ((r₀ + r₁) + (r₂ + r₃)) + ((r₄ + r₅) + (r₆ + r₇)),
// puis le reste dans l'ordre ; au-delà, deux moitiés (la première d'une longueur
// multiple de 8) réduites de même. Reproductible, vectorisable, erreur en O(log n ε).

#include <cstddef>

namespace wsl::detail::mdarray
{
    // x(k) pour k dans [begin, end) ; op associatif « en principe » (+ ou ×).
    template <class T, class X, class Op>
    constexpr T pairwise(X const& x, std::size_t begin, std::size_t end, T identity, Op op)
    {
        std::size_t const n = end - begin;
        if (n < 8)
        {
            T r = identity;
            for (std::size_t k = begin; k < end; ++k)
                r = op(r, x(k));
            return r;
        }
        if (n <= 128)
        {
            T r[8] = {x(begin), x(begin + 1), x(begin + 2), x(begin + 3),
                      x(begin + 4), x(begin + 5), x(begin + 6), x(begin + 7)};
            std::size_t k = begin + 8;
            for (; k + 8 <= end; k += 8)
                for (std::size_t j = 0; j < 8; ++j)
                    r[j] = op(r[j], x(k + j));
            T s = op(op(op(r[0], r[1]), op(r[2], r[3])), op(op(r[4], r[5]), op(r[6], r[7])));
            for (; k < end; ++k)
                s = op(s, x(k));
            return s;
        }
        std::size_t half = n / 2;
        half -= half % 8;
        return op(pairwise(x, begin, begin + half, identity, op), pairwise(x, begin + half, end, identity, op));
    }
}
