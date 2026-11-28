# Documentation

Une bibliothèque scientifique C++23 en en-têtes seuls, sous l'espace de noms `wsl` (nom
provisoire, remplacé avant la première version stable). Cette documentation suit l'arbre des
espaces de noms : un dossier par domaine, un `README.md` par dossier. Les exemples cités sont
dans [`examples/`](../examples), compilés et exécutés par la CI ; la conception de chaque
domaine (choix, algorithmes, mesures) est dans [`design/`](../design).

| Domaine | Alias conseillé | En-tête | Cible CMake | Contenu |
|---|---|---|---|---|
| [`wsl::floating_point_arithmetic`](floating_point_arithmetic/README.md) | `fpa` | `<wsl/floating_point_arithmetic>` | `wsl::floating_point_arithmetic` | primitives flottantes (interface à venir) |
| [`wsl::elementary_functions`](elementary_functions/README.md) | `ef` | `<wsl/elementary_functions>` | `wsl::elementary_functions` | `exp`, `log`, `pow`, `sin`, `sqrt`… reproductibles |
| [`wsl::special_functions`](special_functions/README.md) | `sf` | `<wsl/special_functions>` | `wsl::special_functions` | gamma, erreur, bêta et gamma incomplètes, Bessel |
| [`wsl::random`](random/README.md) | `rnd` | `<wsl/random>` | `wsl::random` | moteurs, tirages uniformes exacts, mélange |
| [`wsl::probability`](probability/README.md) | `proba` | `<wsl/probability>` | `wsl::probability` | lois, algèbre des lois, variables aléatoires |
| [`wsl::probability::stochastic_processes`](probability/stochastic_processes/README.md) | `stochproc` | `<wsl/probability/stochastic_processes>` | `wsl::probability::stochastic_processes` | processus et trajectoires |
| [`wsl::measure_theory`](measure_theory/README.md) | `meas` | `<wsl/measure_theory>` | `wsl::measure_theory` | mesures de référence |
| [`wsl::mdarray`](mdarray/README.md) | — | `<wsl/mdarray>` | `wsl::mdarray` | tableaux multidimensionnels, calcul par les axes |
| `wsl::permutation` (ci-dessous) | — | `<wsl/permutation>` | `wsl::permutation` | permutations de {0, …, n − 1} |

`wsl::wsl` lie tout. Un alias est défini chez l'utilisateur
(`namespace ef = wsl::elementary_functions;`), jamais par la bibliothèque.

## Conventions communes

- **En-têtes sans extension** : `#include <wsl/elementary_functions>` (le chemin est l'espace
  de noms) ; `<wsl/elementary_functions.hpp>` reste valable. Un dossier a son en-tête
  parapluie, un fichier par notion en dessous (`<wsl/special_functions/gamma>`).
- **Paramètres nommés** des objets mathématiques : `normal({.mean = 0.0, .stddev = 1.0})`.
- **`constexpr`** partout où le standard le permet ; `noexcept` quand rien ne peut échouer.
- **C++23**, transmis au projet qui lie une cible.

## Modes de calcul

Deux axes, les mêmes dans toutes les bibliothèques (`<wsl/modes>`, décisions 24 à 27 de
[design/interfaces.md](../design/interfaces.md)). Le défaut est toujours la **permission**, la
voie la plus rapide ; une étiquette **exige** davantage.

| Axe | Permission (défaut) | Exigence |
|---|---|---|
| Reproductibilité | `wsl::deterministic` : mêmes bits d'une exécution à l'autre sur une même configuration (le compilateur peut fusionner en fma) | `wsl::portable` : mêmes bits sur toute configuration IEEE 754 |
| Algorithme | `wsl::approximate` : erreur bornée et documentée | `wsl::exact` : le résultat de l'algorithme idéal (tirage exact, arrondi correct) |

```cpp
double y = ef::exp(x);                             // déterministe, approché
double z = ef::exp(x, wsl::portable);              // mêmes bits partout
double u = N(g, wsl::exact);                       // tirage exact d'une loi
y_[i] = wsl::portable(A[i, j] * x[j] + b[i]);      // enveloppe : toute l'expression
```

- Une exigence qu'une fonction ne sait pas honorer est **refusée à la compilation**
  (`N(g, wsl::exact)` pour une loi sans tirage exact) ; une permission est toujours acceptée,
  la fonction rend ce qu'elle a de mieux. `wsl::provides<F, Mode>` dit si `F` a une vraie
  variante pour une permission.
- Dans une expression portable, une fonction de l'utilisateur se signe
  `wsl::assume_portable(f)` : la bibliothèque ne peut pas inspecter son corps.
- **Défaut du projet** : les options CMake `WSL_PORTABLE` et `WSL_EXACT` (propagées à tout
  consommateur) renforcent le défaut, jamais ne l'affaiblissent. `WSL_EXACT` est aujourd'hui
  strict (décision 41, provisoire) : une fonction sans version exacte est refusée.
- Les **lois** de `probability` n'ont pour l'instant que leur version portable (décision 42) :
  `wsl::deterministic` y est accepté et rend le résultat portable.

## Permutations

`wsl::permutation` représente une permutation de {0, …, n − 1} : composition `s * t`,
`inverse()`, puissances, `cycles()`, `cycle_count()`, `sign()`, `inversions()` en
O(n log n), rang lexicographique, tout `constexpr`. Exemple :
[`examples/permutation/cycles.cpp`](../examples/permutation/cycles.cpp).
