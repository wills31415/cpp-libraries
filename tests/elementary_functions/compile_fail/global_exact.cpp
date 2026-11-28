// WSL_EXACT, mode strict provisoire (design/interfaces.md, décision 41) : le défaut du projet est
// exact, et une fonction sans version exacte est refusée, même appelée sans étiquette.
#define WSL_EXACT 1
#include <wsl/elementary_functions.hpp>
#include <wsl/special_functions.hpp>

// ERREUR 1 : pas de version exacte
// ERREUR 2 : pas de version exacte
// ERREUR 3 : pas de version exacte

double f(double x)
{
#if CASE == 1
    return wsl::elementary_functions::exp(x);                 // le défaut est exact
#elif CASE == 2
    return wsl::special_functions::erf(x);
#elif CASE == 3
    return wsl::elementary_functions::log(x, wsl::portable);  // l'étiquette de l'autre axe n'y change rien
#else
    return wsl::elementary_functions::sqrt(x) + wsl::elementary_functions::exp(x, wsl::approximate);
#endif
}
