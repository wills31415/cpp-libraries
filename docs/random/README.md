# `wsl::random` (`rnd`)

`#include <wsl/random>` — engendrer le hasard : moteurs et fonctions de tirage simples. Les
lois sont dans [`probability`](../probability/README.md).

**Reproductibilité** : à graine égale, les mêmes valeurs au bit près sur toute plateforme,
et le même nombre de valeurs consommées sur le moteur (cela fait partie du contrat). Aucune
loi de `<random>` n'est employée (leurs algorithmes changent d'une bibliothèque standard à
l'autre) ; les moteurs standard sont acceptés partout où un moteur l'est
(`std::uniform_random_bit_generator`).

| Notion | Noms |
|---|---|
| Moteurs (`rnd::engines`) | `splitmix64`, `xoshiro256pp` (sauts compris), `pcg64`, `philox4x64` (à compteur) |
| Tirages | `random_bits`, `uniform_below(n, g)`, `uniform_int(a, b, g)`, `uniform_real(g, modes…)`, `bernoulli_trial`, `shuffle` |

```cpp
namespace rnd = wsl::random;
rnd::engines::pcg64 g(42);
double u = rnd::uniform_real(g);                 // approché (défaut) : le plus rapide
double v = rnd::uniform_real(g, wsl::exact);     // arrondi correct d'un uniforme réel
auto   d = rnd::uniform_int(1, 6, g);            // exact, sans biais
```

Les tirages entiers sont toujours exacts (sans biais) ; `uniform_real` est approché par
défaut, exact sur demande ([modes](../README.md#modes-de-calcul)). Conception :
[design/random.md](../../design/random.md). Exemple :
[`examples/random/uniform.cpp`](../../examples/random/uniform.cpp).
