// Fonctions élémentaires : les deux modes de reproductibilité (design/interfaces.md, décisions 24
// à 27). Par défaut, le mode déterministe : même résultat d'une exécution à l'autre sur une même
// machine, le compilateur pouvant fusionner les produits en fma. wsl::portable : mêmes bits sur
// toute plateforme, au prix de quelques pour cent.

#include <wsl/elementary_functions>

#include <cstdio>

namespace ef = wsl::elementary_functions;

int main()
{
    double const x = 0.75;
    double const fast = ef::exp(x);                    // déterministe (le défaut)
    double const same = ef::exp(x, wsl::portable);     // mêmes bits sur toute plateforme
    std::printf("exp(%g) = %.17g (déterministe), %.17g (portable)\n", x, fast, same);

    // Calcul à la compilation : c'est la voie portable, la même qu'à l'exécution.
    constexpr double c = ef::log(2.0, wsl::portable);
    std::printf("log(2) = %.17g, constexpr\n", c);
    return c == ef::log(2.0, wsl::portable) ? 0 : 1;
}
