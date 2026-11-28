# cpp-libraries

[![CI](https://github.com/wills31415/cpp-libraries/actions/workflows/ci.yml/badge.svg)](https://github.com/wills31415/cpp-libraries/actions/workflows/ci.yml)

Une bibliothèque scientifique C++23 en en-têtes seuls : mathématiques, physique, calcul
scientifique, et leurs fondations. Espace de noms `wsl` (nom provisoire, remplacé avant la
première version stable). Objectif : des objets mathématiques au plus près de la
bibliothèque standard (nommage, garanties, `constexpr`, concepts), les options les plus
rapides par défaut, la **reproductibilité** au bit près sur demande, et des compatibilités
vérifiées à la compilation.

Testée à chaque commit avec GCC 14 et Clang 18 (libstdc++, Linux) et Apple Clang (libc++,
macOS), en Debug (ASan + UBSan) et en Release, avec et sans fma, avertissements traités
comme des erreurs.

```cpp
#include <wsl/elementary_functions>
#include <wsl/probability>
#include <wsl/random>

namespace ef = wsl::elementary_functions;
namespace proba = wsl::probability;

double y = ef::exp(0.75, wsl::portable);                  // mêmes bits sur toute plateforme
proba::laws::normal const N({.mean = 1.0, .stddev = 2.0});
wsl::random::engines::pcg64 g(7);
double p = N.cdf(3.0), x = N(g), z = N(g, wsl::exact);   // tirage approché, puis exact
```

## Domaines

Documentation : [docs/](docs/README.md), un guide par domaine ; exemples compilés et
exécutés par la CI : [examples/](examples) ; conception et charte : [design/](design).

| Domaine (alias conseillé) | Contenu |
|---|---|
| [`wsl::elementary_functions`](docs/elementary_functions/README.md) (`ef`) | `exp`, `log`, `pow`, `sin`, `atan2`, `cbrt`, `hypot`… reproductibles, `constexpr`, erreur près de 0,5 ulp |
| [`wsl::special_functions`](docs/special_functions/README.md) (`sf`) | gamma, digamma, erreur et réciproques, gamma et bêta incomplètes et réciproques, Bessel `I_ν` |
| [`wsl::random`](docs/random/README.md) (`rnd`) | moteurs `splitmix64`, `xoshiro256pp`, `pcg64`, `philox4x64` ; tirages uniformes exacts, mélange |
| [`wsl::probability`](docs/probability/README.md) (`proba`) | plus de 80 lois (continues, discrètes, multivariées, combinatoires, matrices aléatoires), algèbre des lois, variables aléatoires |
| [`wsl::probability::stochastic_processes`](docs/probability/stochastic_processes/README.md) (`stochproc`) | chaînes de Markov, marches, urnes, mouvement brownien, Ornstein–Uhlenbeck, Poisson… à trajectoires paresseuses |
| [`wsl::measure_theory`](docs/measure_theory/README.md) (`meas`) | mesures de référence (à étendre) |
| [`wsl::mdarray`](docs/mdarray/README.md) | tableaux multidimensionnels, calcul par les axes (`c[i, k] = a[i, j] * b[j, k]`), `mdspan` |
| [`wsl::permutation`](docs/README.md#permutations) | permutations : composition, cycles, signature, rang |
| [`wsl::floating_point_arithmetic`](docs/floating_point_arithmetic/README.md) (`fpa`) | primitives flottantes (interface à venir) |

Les **modes de calcul** sont les mêmes partout : par défaut la voie la plus rapide
(déterministe, approchée) ; `wsl::portable` exige les mêmes bits sur toute plateforme,
`wsl::exact` le résultat exact ; une exigence impossible est refusée à la compilation
([détails](docs/README.md#modes-de-calcul)).

Chaque domaine a sa cible CMake du même nom (`wsl::random`,
`wsl::probability::stochastic_processes`…) et son en-tête sans extension
(`<wsl/random>`) ; `wsl::wsl` lie tout. Toutes exigent C++23, transmis au projet
consommateur. L'historique des changements est dans [CHANGELOG.md](CHANGELOG.md). Les
anciennes bibliothèques (`array2d`, `graded`, `graded_io`, `elo`, `progress`) attendent leur
refonte dans `legacy/`, figées.

## Utiliser les bibliothèques dans un projet

### FetchContent, figé sur une version (recommandé)

```cmake
include(FetchContent)
FetchContent_Declare(wsl
    GIT_REPOSITORY https://github.com/wills31415/cpp-libraries.git
    GIT_TAG        v0.1.0
    GIT_SHALLOW    TRUE
    EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(wsl)

target_link_libraries(mon_programme PRIVATE wsl::mdarray wsl::random)
```

Ces noms (`wsl`, `wsl::wsl`) valent à partir de la version qui suivra `v0.1.0`, qui
emploie encore `wsl_libraries` et `wsl::libraries`.

Chaque projet garde la version qu'il a choisie : une nouvelle version de
cpp-libraries ne change rien tant que `GIT_TAG` n'est pas modifié. Les tests et
exemples de ce dépôt ne sont pas construits chez le consommateur.


### Modifier une bibliothèque en même temps qu'un projet

Pointer FetchContent sur un clone local, sans toucher au `CMakeLists.txt` du
projet, par un preset dans `CMakeUserPresets.json` (fichier propre à la machine,
à ne pas versionner) :

```json
{
    "version": 6,
    "configurePresets": [{
        "name": "debug-wsl-local", "inherits": "debug",
        "cacheVariables": { "FETCHCONTENT_SOURCE_DIR_WSL": "/chemin/vers/cpp-libraries" }
    }]
}
```

Les modifications du clone sont alors visibles immédiatement, et « aller à la
définition » ouvre le vrai fichier source.

### Paquet installé

```sh
cmake --preset release && cmake --build --preset release && ctest --preset release
cmake --install build/release --prefix ~/.local
```

```cmake
find_package(wsl 0.1 REQUIRED)
target_link_libraries(mon_programme PRIVATE wsl::mdarray)
```

Pratique pour un essai rapide, mais une seule version existe alors pour toute la
machine : une version incompatible installée plus tard casse les anciens projets.

## Développement

| Action | VSCode | Terminal |
|---|---|---|
| Compiler | ⇧⌘B | `cmake --preset debug && cmake --build --preset debug -j` |
| Tests (ASan + UBSan, sans les tests longs) | tâche « Tests » | `ctest --preset debug` |
| Tests en Release (tous) | tâche « Tests (Release) » | `ctest --preset release` |
| Vérification rapide (bibliothèques touchées, Release et fma) | — | `tools/check.sh` |
| Déboguer le test ouvert | F5 → « Déboguer le test courant » | — |

Les tests sont rangés comme les espaces de noms : `tests/probability/laws.cpp` donne le test
`probability.laws`, `tests/mdarray/axis.cpp` le test `mdarray.axis`. Une ligne
`// Étiquettes : precision long` en tête du fichier lui donne ses étiquettes
(`ctest -L precision`, `ctest -LE long`). Les exemples, `examples/<domaine>/<tâche>.cpp`,
sont lancés comme tests `example.<domaine>.<tâche>` (étiquette `example`). Un test réussit
s'il renvoie 0. `tests/package/` vérifie, en CI, l'installation et FetchContent du point de
vue d'un consommateur. La CI ne compile que ce qu'un push peut avoir changé, et ce qui en
dépend.

### Ajouter une bibliothèque

1. Ses en-têtes dans `include/wsl/<chemin>`, le chemin étant l'espace de noms
   (`include/wsl/linear_algebra.hpp`, `include/wsl/linear_algebra/…`), dans
   `namespace wsl::inline WSL_DEFAULTS::<chemin>` si elle dépend des modes.
2. Dans `CMakeLists.txt` : `wsl_add_library(<nom avec ::> STANDARD 23 [DEPENDS ...])` ;
   les relais sans extension sont engendrés.
3. Ses tests dans `tests/<chemin>.cpp` (ou `tests/<chemin>/<batterie>.cpp`), sa documentation
   dans `docs/<chemin>/README.md`, un exemple dans `examples/<domaine>/`.
4. Une entrée dans `CHANGELOG.md` ; les règles de nommage sont dans
   [design/interfaces.md](design/interfaces.md) (§6).

### Publier une version

Numérotation `MAJEUR.MINEUR.CORRECTIF` ; tant que `MAJEUR` vaut 0, une version
mineure peut casser l'interface, un correctif jamais.

1. Dans un commit poussé sur `main` : `VERSION` dans `project()` de
   `CMakeLists.txt`, et la section `[Non publié]` du `CHANGELOG.md` renommée en
   `[X.Y.Z] — date`.
2. `tools/release.sh X.Y.Z` : vérifie version, journal, tests locaux et CI verte,
   puis crée le tag `vX.Y.Z` et la release GitHub.

## Outils

- **`tools/release.sh`** : publication d'une version (voir plus haut).
- **`tools/check.sh`** : vérification rapide des bibliothèques touchées.
- **`tools/elementary_functions/generate.py`** : tables et valeurs de référence (mpmath) des
  fonctions ; **`tools/elementary_functions/benchmark.sh`** : mesures de performance.
- **`tools/random/`** : vecteurs des moteurs, tables de la ziggourat ;
  **`tools/probability/`** : références des lois (mpmath) ;
  **`tools/mdarray/benchmark.sh`** : performances de mdarray.

## Licence

[MIT](LICENSE).
