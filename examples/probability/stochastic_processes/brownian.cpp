// Processus stochastiques : une trajectoire brownienne observée sur une grille de temps.

#include <wsl/probability/stochastic_processes>
#include <wsl/random>

#include <cstdio>

namespace stochproc = wsl::probability::stochastic_processes;

int main()
{
    wsl::random::engines::pcg64 g(2026);
    auto const W = stochproc::processes::brownian_motion({.drift = 0.0, .volatility = 1.0}).on_grid({0.25, 0.5, 1.0});
    auto const w = W(g);   // un tirage : les valeurs aux trois instants
    std::printf("W(0.25) = %.6f, W(0.5) = %.6f, W(1) = %.6f\n", w[0], w[1], w[2]);
    return 0;
}
