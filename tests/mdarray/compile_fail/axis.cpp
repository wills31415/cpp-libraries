// Axes : refus à la compilation (design/mdarray.md, L1), avec leur message.
#include <wsl/mdarray.hpp>

// ERREUR 1 : étiquettes en double
// ERREUR 2 : au moins un axe
// ERREUR 3 : deleted
// ERREUR 4 : auto
// ERREUR 5 : no matching function

void f()
{
#if CASE == 1
    auto [a, b] = wsl::axis<"x", "x">();
#elif CASE == 2
    auto t = wsl::axis<0>();
#elif CASE == 3
    auto i = wsl::axis();
    auto j = i; // copie : deux noms pour un même axe
#elif CASE == 4
    auto i = wsl::axis(), j = wsl::axis(); // déclaration groupée : types différents
#elif CASE == 5
    auto t = wsl::axis<2, "x">(); // nombre et étiquette mêlés
#endif
}
