// Tableaux multidimensionnels : le calcul par les axes, écrit comme au tableau. Le membre de
// gauche décide des axes libres ; un axe répété à droite seulement est sommé.

#include <wsl/mdarray>

#include <cstdio>

int main()
{
    using matrix = wsl::mdarray<double, wsl::dims<2>>;
    using vector = wsl::mdarray<double, wsl::dims<1>>;
    matrix const a{{1, 2}, {3, 4}}, b{{0, 1}, {1, 0}};
    matrix c(2, 2);
    vector r(2);
    auto [i, j, k] = wsl::axis<3>();
    c[i, k] = a[i, j] * b[j, k];    // produit matriciel
    r[i] = sum(j)(a[i, j]);         // sommes des lignes
    std::printf("c = [[%g, %g], [%g, %g]], r = [%g, %g]\n", c[0, 0], c[0, 1], c[1, 0], c[1, 1], r[0], r[1]);
    return c[0, 0] == 2 && r[1] == 7 ? 0 : 1;
}
