// wsl::mdarray, modes de calcul (docs/mdarray.md, I9) : empreinte du mode
// reproductible, identique sur toute plateforme et avec toutes les options (fma
// comprise) ; calcul constexpr avec des fonctions ; mode rapide proche.

#include "test_mdarray_modes_common.hpp"

#include <cassert>
#include <cinttypes>
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
        double const n = wsl::reproducible(sqrt(x[i] * x[i]));
        y[i] = wsl::reproducible(exp(x[i] * 0) + x[j] * x[j] * 0);
        return n == 5 && y[2] == 1;
    }
    static_assert(compile_time());
}

int main()
{
    std::uint64_t const r = modes::compute(modes::reproducible);
    std::printf("test_mdarray_modes : empreinte reproductible 0x%016" PRIx64 "\n", r);
    // Sans marque : mode rapide par défaut ; résultats proches du mode reproductible.
    (void)modes::compute(modes::unmarked);
    (void)modes::compute(modes::fast);
    if (r != expected)
    {
        std::printf("empreinte attendue 0x%016" PRIx64 "\n", expected);
        return 1;
    }
    std::puts("test_mdarray_modes : OK");
}
