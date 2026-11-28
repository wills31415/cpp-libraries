// Permutations de {0, …, n − 1} : composition, inverse, cycles, signature, constexpr.

#include <wsl/permutation>

#include <cstdio>

int main()
{
    // Évaluable à la compilation (la mémoire allouée y est rendue avant la fin) :
    static_assert((wsl::permutation{2, 0, 1, 4, 3} * wsl::permutation{2, 0, 1, 4, 3}.inverse()).is_identity());
    wsl::permutation const s{2, 0, 1, 4, 3};
    std::printf("s a %zu cycles :", s.cycle_count());
    for (auto const& c : s.cycles())
    {
        std::printf(" (");
        for (auto x : c)
            std::printf(" %zu", static_cast<std::size_t>(x));
        std::printf(" )");
    }
    std::printf(", %llu inversions\n", static_cast<unsigned long long>(s.inversions()));
    return 0;
}
