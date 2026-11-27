// Étiquettes : modes
// wsl::mdarray, modes de calcul (design/mdarray.md, I9) : empreinte du mode
// reproductible, identique sur toute plateforme et avec toutes les options (fma
// comprise) ; calcul constexpr avec des fonctions ; mode rapide proche.

#include "modes_common.hpp"

#include <bit>
#include <cassert>
#include <cinttypes>
#include <cstdint>
#include <cstdio>

namespace
{
    // Empreinte attendue, la même partout (un changement voulu d'algorithme la change :
    // rupture, au CHANGELOG ; un changement sur une seule plateforme est un bug).
    constexpr std::uint64_t expected = 0xf33b15d552bca1bcULL;

    // Reproductible : calcul à la compilation possible (fonctions de wsl::special).
    constexpr bool compile_time()
    {
        auto [i, j] = wsl::axis<2>();
        wsl::mdarray<double, wsl::extents<std::size_t, 3>> x{}, y{};
        x[0] = 3;
        x[1] = 4;
        double const n = wsl::portable(sqrt(x[i] * x[i]));
        y[i] = wsl::portable(exp(x[i] * 0) + x[j] * x[j] * 0);
        return n == 5 && y[2] == 1;
    }
    static_assert(compile_time());

    // Tableaux entiers par paquets de 8 (barrière vectorielle) : mêmes bits que la voie
    // scalaire, prise par un tableau en layout_left (pas d'accès direct).
    bool same(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }

    void lanes_match_scalar()
    {
        using E = wsl::dims<2>;
        constexpr std::size_t n = 13, m = 7; // 91 éléments : paquets de 8 et un reste
        wsl::mdarray<double, E> a(n, m), b(n, m), c(n, m);
        wsl::mdarray<double, E, wsl::layout_left> al(n, m), bl(n, m), cl(n, m);
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t q = 0; q < m; ++q)
            {
                a[r, q] = al[r, q] = 1 / (1.0 + r + 3 * q);
                b[r, q] = bl[r, q] = 0.1 * r - 0.7 * q + 1 / 3.0;
            }
        c = wsl::portable(a * b + 1 - a * a * b);
        cl = wsl::portable(al * bl + 1 - al * al * bl);
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t q = 0; q < m; ++q)
                assert(same(c[r, q], cl[r, q]));
        assert(same(wsl::sum(wsl::portable(a * b + a)), wsl::sum(wsl::portable(al * bl + al))));
    }
}

int main()
{
    lanes_match_scalar();
    std::uint64_t const r = modes::compute(modes::reproducible);
    std::printf("mdarray.modes : empreinte reproductible 0x%016" PRIx64 "\n", r);
    // Sans marque : mode rapide par défaut ; résultats proches du mode reproductible.
    (void)modes::compute(modes::unmarked);
    (void)modes::compute(modes::fast);
    if (r != expected)
    {
        std::printf("empreinte attendue 0x%016" PRIx64 "\n", expected);
        return 1;
    }
    std::puts("mdarray.modes : OK");
}
