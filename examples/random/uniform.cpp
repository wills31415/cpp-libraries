// Tirages simples : moteurs reproductibles et tirages uniformes. Le tirage réel est approché par
// défaut (le plus rapide) ; wsl::exact donne le réel correctement arrondi d'un uniforme réel.

#include <wsl/random>

#include <array>
#include <cstdio>

namespace rnd = wsl::random;

int main()
{
    rnd::engines::xoshiro256pp g(42);                          // même suite partout
    double const u = rnd::uniform_real(g);                     // approché (le défaut)
    double const v = rnd::uniform_real(g, wsl::exact);         // exact
    auto const d = rnd::uniform_int(1, 6, g);                  // dé équilibré, sans biais
    std::array<int, 5> cards{1, 2, 3, 4, 5};
    rnd::shuffle(cards, g);
    std::printf("u = %.17g, v = %.17g, dé = %lld, mélange :", u, v, static_cast<long long>(d));
    for (int c : cards)
        std::printf(" %d", c);
    std::printf("\n");
    return d >= 1 && d <= 6 ? 0 : 1;
}
