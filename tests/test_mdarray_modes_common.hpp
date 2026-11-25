// Calculs communs aux tests des modes (test_mdarray_modes*.cpp) : leur empreinte en
// mode reproductible doit être la même sur toute plateforme, fma comprise (I9).
#pragma once

#include <wsl/mdarray.hpp>

#include <bit>
#include <cstddef>
#include <cstdint>

namespace modes
{
    using V = wsl::mdarray<double, wsl::dims<1>>;
    using M = wsl::mdarray<double, wsl::dims<2>>;

    struct fingerprint
    {
        std::uint64_t h = 0xcbf29ce484222325ULL;
        void add(double x)
        {
            auto const b = std::bit_cast<std::uint64_t>(x);
            for (int k = 0; k < 8; ++k)
            {
                h ^= (b >> (8 * k)) & 0xff;
                h *= 0x100000001b3ULL;
            }
        }
        template <class A>
        void add_all(A const& a)
        {
            for (double x : a)
                add(x);
        }
    };

    // R : une fonction qui met une expression dans le mode à éprouver.
    template <class R>
    std::uint64_t compute(R mode)
    {
        auto [i, j, k] = wsl::axis<3>();
        M A(5, 6), B(6, 4), C(5, 4);
        V x(6), b(5), y(5), s6(6);
        for (std::size_t p = 0; p < 5; ++p)
        {
            b[p] = (double(p) + 0.5) / 3.0;
            for (std::size_t q = 0; q < 6; ++q)
                A[p, q] = (double(p) + 1) * (double(q) + 2) / 13.0;
        }
        for (std::size_t p = 0; p < 6; ++p)
        {
            x[p] = (double(p) + 1) / 7.0;
            for (std::size_t q = 0; q < 4; ++q)
                B[p, q] = (double(q) - double(p)) / 11.0;
        }
        fingerprint f;
        C[i, k] = mode(A[i, j] * B[j, k]);
        f.add_all(C);
        y[i] = mode(A[i, j] * x[j] + b[i]);
        f.add_all(y);
        double const n = mode(sqrt(x[i] * x[i]));
        f.add(n);
        s6[i] = mode(exp(x[i]) / sum(j)(exp(x[j])));
        f.add_all(s6);
        V z = mode(exp(x) * log1p(x) + hypot(x, x * 2) / pow(x, 1.5));
        f.add_all(z);
        V big(1000);
        for (std::size_t p = 0; p < 1000; ++p)
            big[p] = 1.0 / (double(p) + 3.0);
        f.add(sum(mode(big * big + 0.1 * big)));
        double const q = mode(x[i] * A[k, i] * b[k]);
        f.add(q);
        return f.h;
    }

    inline auto reproducible = [](auto&& e) { return wsl::reproducible(std::move(e)); };
    inline auto fast = [](auto&& e) { return wsl::fast(std::move(e)); };
    inline auto unmarked = [](auto&& e) { return std::move(e); };
}
