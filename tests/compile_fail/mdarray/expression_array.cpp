// Monde des tableaux entiers : refus à la compilation (docs/mdarray.md, V8, L2, I13).
#include <wsl/mdarray.hpp>

// ERREUR 1 : rangs différents
// ERREUR 2 : étendues fixes différentes
// ERREUR 3 : une expression ne se garde pas
// ERREUR 4 : un masque n'est pas un booléen
// ERREUR 5 : une expression ne se garde pas

void f()
{
    wsl::mdarray<double, wsl::dims<2>> a(2, 2);
    wsl::mdarray<double, wsl::extents<std::size_t, 2, 3>> b;
    (void)a;
    (void)b;
#if CASE == 1
    wsl::mdarray<double, wsl::dims<3>> c(2, 2, 2);
    auto e = a + c;
#elif CASE == 2
    wsl::mdarray<double, wsl::extents<std::size_t, 3, 2>> c;
    auto e = b + c;
#elif CASE == 3
    auto e = a + 1;
    auto g = e * 2; // expression nommée réutilisée
#elif CASE == 4
    if (a < a)
        return;
#elif CASE == 5
    auto e = a + 1;
    wsl::mdarray<double, wsl::dims<2>> y = e;
#endif
}
