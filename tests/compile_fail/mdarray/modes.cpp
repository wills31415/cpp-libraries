// Modes de calcul : refus à la compilation (docs/mdarray.md, I9, L2).
#include <wsl/mdarray.hpp>

// ERREUR 1 : n'est pas couverte par wsl::special
// ERREUR 2 : ne se garde pas
// ERREUR 3 : n'est pas couverte par wsl::special

void f()
{
    auto [i] = wsl::axis<1>();
    wsl::mdarray<double, wsl::dims<1>> x(3), y(3);
    (void)i;
#if CASE == 1
    y = wsl::reproducible(tanh(x));
#elif CASE == 2
    auto e = x * 2;
    y = wsl::reproducible(e);
#elif CASE == 3
    y[i] = wsl::reproducible(asin(x[i]));
#endif
    (void)y;
}
