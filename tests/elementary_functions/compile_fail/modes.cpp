// Modes des fonctions élémentaires : refus à la compilation (design/interfaces.md, décision 25).
#include <wsl/elementary_functions.hpp>

// ERREUR 1 : pas de version exacte
// ERREUR 2 : pas de version exacte
// ERREUR 3 : pas de version exacte

namespace ef = wsl::elementary_functions;

double f(double x)
{
#if CASE == 1
    return ef::exp(x, wsl::exact);                     // erreur < 0,51 ulp : pas l'arrondi correct
#elif CASE == 2
    return ef::sin(x, wsl::portable, wsl::exact);
#elif CASE == 3
    return ef::pow(x, 2.0, wsl::exact);
#else
    return ef::sqrt(x, wsl::exact) + ef::exp(x, wsl::approximate);   // sqrt est exacte ; approximate permis
#endif
}
