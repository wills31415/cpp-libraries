#pragma once

#include <wsl/modes.hpp>

// Ordre fixé des sommes et produits flottants (design/mdarray.md, I9) : par paires,
// comme NumPy (pairwise_sum) — en dessous de 8 termes, dans l'ordre ; jusqu'à 128,
// 8 accumulateurs entrelacés combinés ((r₀ + r₁) + (r₂ + r₃)) + ((r₄ + r₅) + (r₆ + r₇)),
// puis le reste dans l'ordre ; au-delà, deux moitiés (la première d'une longueur
// multiple de 8) réduites de même. Reproductible, vectorisable, erreur en O(log n ε).

#include <cstddef>

#include <wsl/mdarray/detail/expression_base.hpp>

namespace wsl::inline WSL_DEFAULTS::detail::mdarray
{
    template <class T, class X, class Op>
    constexpr T pairwise_split(X const& x, std::size_t begin, std::size_t end, T identity, Op op);

    // x(k) pour k dans [begin, end) ; op associatif « en principe » (+ ou ×). Non
    // récursive, pour être intégrée aux petites sommes ; le partage en deux moitiés
    // est dans pairwise_split.
    template <class T, class X, class Op>
    [[gnu::always_inline]] constexpr T pairwise(X const& x, std::size_t begin, std::size_t end, T identity, Op op)
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
            T r[8];
            std::size_t k = begin + 8;
            if constexpr (requires { x.template pack<8>(begin).v[0]; })
            {
                // Termes par paquets de 8 (voies, indexed.hpp) : mêmes termes, même ordre.
                auto const first = x.template pack<8>(begin);
                for (std::size_t j = 0; j < 8; ++j)
                    r[j] = first.v[j];
                for (; k + 8 <= end; k += 8)
                {
                    auto const p = x.template pack<8>(k);
                    for (std::size_t j = 0; j < 8; ++j)
                        r[j] = op(r[j], p.v[j]);
                }
            }
            else
            {
                for (std::size_t j = 0; j < 8; ++j)
                    r[j] = x(begin + j);
                for (; k + 8 <= end; k += 8)
                    for (std::size_t j = 0; j < 8; ++j)
                        r[j] = op(r[j], x(k + j));
            }
            T s = op(op(op(r[0], r[1]), op(r[2], r[3])), op(op(r[4], r[5]), op(r[6], r[7])));
            for (; k < end; ++k)
                s = op(s, x(k));
            return s;
        }
        return pairwise_split(x, begin, end, identity, op);
    }

    // Hors ligne : sinon le compilateur la fusionne dans pairwise, qui redevient récursive
    // et cesse d'être intégrée.
    template <class T, class X, class Op>
    [[gnu::noinline]] constexpr T pairwise_split(X const& x, std::size_t begin, std::size_t end, T identity, Op op)
    {
        std::size_t half = (end - begin) / 2;
        half -= half % 8;
        return op(pairwise(x, begin, begin + half, identity, op), pairwise(x, begin + half, end, identity, op));
    }

    // Tranches de la somme par paires : pairwise(x, begin, end) somme d'un bloc chaque
    // tranche [b, e) de pairwise_leaves (au plus 128 termes, dans l'ordre), puis combine
    // leurs résultats selon l'arbre de pairwise_combine. Calculer les tranches séparément
    // (découpage du cache) puis les combiner ainsi donne les mêmes bits.
    template <class F>
    constexpr void pairwise_leaves(std::size_t begin, std::size_t end, F&& f)
    {
        if (end - begin <= 128)
        {
            f(begin, end);
            return;
        }
        std::size_t half = (end - begin) / 2;
        half -= half % 8;
        pairwise_leaves(begin, begin + half, f);
        pairwise_leaves(begin + half, end, f);
    }

    template <class T, class Op>
    constexpr T pairwise_combine(std::size_t begin, std::size_t end, T const* leaves, std::size_t& next, Op op)
    {
        if (end - begin <= 128)
            return leaves[next++];
        std::size_t half = (end - begin) / 2;
        half -= half % 8;
        T const left = pairwise_combine(begin, begin + half, leaves, next, op);
        T const right = pairwise_combine(begin + half, end, leaves, next, op);
        return op(left, right);
    }
}
