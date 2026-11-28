// Modes des tirages : refus à la compilation (design/interfaces.md, décisions 24 à 27, 41).
#if CASE == 3 || CASE == 4
#    define WSL_EXACT 1
#endif
#include <wsl/probability.hpp>
#include <wsl/random.hpp>

// ERREUR 1 : n'a pas de tirage exact
// ERREUR 2 : n'a pas de tirage exact
// ERREUR 3 : n'a pas de tirage exact
// ERREUR 4 : n'a pas de tirage exact

namespace proba = wsl::probability;

double f()
{
    wsl::random::engines::pcg64 g(1);
    proba::laws::gamma const G({.rate = 1.0, .shape = 2.5});      // tirage approché seulement
    proba::laws::exponential const E({.rate = 2.0});               // tirage exact
#if CASE == 1
    return G(g, wsl::exact);
#elif CASE == 2
    return proba::sample(G, g, wsl::exact);
#elif CASE == 3
    return G(g);                                                    // WSL_EXACT : le défaut est exact
#elif CASE == 4
    return proba::laws::normal<double>({.mean = 0.0, .stddev = 1.0})(g) + G(g, wsl::portable);
#else
    return E(g) + E(g, wsl::exact) + G(g) + G(g, wsl::approximate, wsl::portable);
#endif
}
