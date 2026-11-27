// Modes de calcul : refus à la compilation (design/mdarray.md, I9, L2 ; design/interfaces.md, décision 25).
#include <wsl/mdarray.hpp>

// ERREUR 1 : n'a pas de version portable
// ERREUR 2 : ne se garde pas
// ERREUR 3 : n'a pas de version portable
// ERREUR 4 : doit être signée wsl::assume_portable
// ERREUR 5 : doit être signée wsl::assume_portable

void f()
{
    auto [i] = wsl::axis<1>();
    wsl::mdarray<double, wsl::dims<1>> x(3), y(3);
    (void)i;
#if CASE == 1
    y = wsl::portable(tanh(x));
#elif CASE == 2
    auto e = x * 2;
    y = wsl::portable(e);
#elif CASE == 3
    y[i] = wsl::portable(asin(x[i]));
#elif CASE == 4
    auto square = wsl::vectorize([](double v) { return v * v + 1; });
    y = wsl::portable(square(x));
#elif CASE == 5
    auto square = wsl::vectorize([](double v) { return v * v + 1; });
    y[i] = wsl::portable(square(x[i]));
#endif
    // Signée, la même fonction est acceptée.
    auto signed_square = wsl::vectorize(wsl::assume_portable([](double v) { return v * v + 1; }));
    y = wsl::portable(signed_square(x));
    y[i] = wsl::portable(signed_square(x[i]));
    (void)y;
}
