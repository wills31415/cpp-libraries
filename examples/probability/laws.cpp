// Lois de probabilité : des objets mathématiques, paramètres nommés, opérations en membres,
// tirage X(g). Algèbre des lois aux noms de la théorie de la mesure.

#include <wsl/probability>
#include <wsl/random>

#include <cstdio>

namespace proba = wsl::probability;
namespace rnd = wsl::random;

int main()
{
    proba::laws::normal const N({.mean = 1.0, .stddev = 2.0});
    std::printf("N : P(X <= 3) = %.17g, quantile 0.975 = %.17g, moyenne %g\n", N.cdf(3.0), N.quantile(0.975), N.mean());

    rnd::engines::pcg64 g(7);
    double const x = N(g);                  // tirage approché (ziggourat), le défaut
    double const y = N(g, wsl::exact);      // tirage exact : arrondi correct d'un tirage réel exact
    std::printf("tirages : %.17g (approché), %.17g (exact)\n", x, y);

    // Somme de deux lois indépendantes : Poisson(2) + Poisson(3) = Poisson(5), reconnue.
    auto const S = proba::convolution(proba::laws::poisson({.mean = 2.0}), proba::laws::poisson({.mean = 3.0}));
    std::printf("Poisson(2) * Poisson(3) : P(S = 4) = %.17g\n", S.pmf(4));
    return N.quantile(N.cdf(3.0)) > 2.999 ? 0 : 1;
}
