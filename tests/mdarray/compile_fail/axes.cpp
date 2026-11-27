// Monde des axes : refus à la compilation (design/mdarray.md, V2, V4, V9, L2, L5, L12, I6, I10).
#include <wsl/mdarray.hpp>

// ERREUR 1 : n'apparaît qu'une fois
// ERREUR 2 : portée ambiguë
// ERREUR 3 : lié par une réduction
// ERREUR 4 : lié par une réduction
// ERREUR 5 : somme de booléens
// ERREUR 6 : longueur inconnue
// ERREUR 7 : longueurs différentes
// ERREUR 8 : est une valeur, pas un axe
// ERREUR 9 : est réservé
// ERREUR 10 : ne s'appelle pas
// ERREUR 11 : ne se garde pas
// ERREUR 12 : ne se mêle pas
// ERREUR 13 : ne se garde pas
// ERREUR 14 : operator[]

void f()
{
    auto [i, j] = wsl::axis<2>();
    wsl::mdarray<double, wsl::dims<1>> x(3), y(3);
    wsl::mdarray<double, wsl::dims<2>> a(3, 3);
    (void)i;
    (void)j;
#if CASE == 1
    y[i] = a[i, j]; // j une seule fois, absent à gauche
#elif CASE == 2
    double s = sqrt(x[j] * x[j]) + y[j] * y[j];
#elif CASE == 3
    y[j] = sum(j)(a[i, j]) * x[i];
#elif CASE == 4
    y[i] = sum(j)(a[i, j]) * x[j];
#elif CASE == 5
    double s = x[i] < y[i];
#elif CASE == 6
    double s = i * i;
#elif CASE == 7
    wsl::mdarray<double, wsl::extents<std::size_t, 2>> u;
    wsl::mdarray<double, wsl::extents<std::size_t, 3>> v;
    double s = u[i] * v[i];
#elif CASE == 8
    y[i] = x[i + 1];
#elif CASE == 9
    y[i] = x[i] * +i;
#elif CASE == 10
    double s = max(i, j)(a[i, j]);
#elif CASE == 11
    auto e = x[i] * x[i];
    double s = e;
#elif CASE == 12
    y[i] = x + x[i];
#elif CASE == 13
    auto t = y[i];
    t = 1;
#elif CASE == 14
    double s = a[i] * a[i];
#endif
}
