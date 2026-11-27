# cpp-libraries

[![CI](https://github.com/wills31415/cpp-libraries/actions/workflows/ci.yml/badge.svg)](https://github.com/wills31415/cpp-libraries/actions/workflows/ci.yml)

Collection de bibliothèques C++ génériques, en en-têtes seuls, sous le namespace
`wsl`. Objectif : du code réutilisable, optimisé et aussi proche du style et des
garanties de la bibliothèque standard que possible.

Testées à chaque commit avec GCC 14 et Clang 18 (libstdc++, Linux) et Apple Clang
(libc++, macOS), en Debug (ASan + UBSan) et en Release, en C++23,
avertissements traités comme des erreurs.

## Catalogue

| Cible CMake | En-tête | Standard | Contenu |
|---|---|---|---|
| `wsl::permutation` | `<wsl/permutation.hpp>` | C++23 | Permutations de {0, …, n − 1} : composition, inverse, puissances, cycles, signature, inversions (O(n log n)), rang lexicographique, `constexpr` |
| `wsl::special` | `<wsl/special.hpp>` | C++23 | Fonctions mathématiques **bit-exactes** (même résultat au bit près sur toute plateforme, `constexpr`) : `exp`, `exp2`, `expm1`, `log`, `log2`, `log10`, `log1p`, `pow`, `sin`, `cos`, `tan`, `atan`, `atan2`, `sqrt`, `cbrt`, `hypot`, `erf`, `erfc`, `erf_inv`, `erfc_inv`, `lgamma`, `tgamma`, `digamma`, `gamma_p`, `gamma_q`, `gamma_p_inv`, `gamma_q_inv`, `beta`, `lbeta`, `ibeta`, `ibetac`, `ibeta_inv`, `ibetac_inv`, `cyl_bessel_i` ([conception](design/special.md)) |
| `wsl::random` | `<wsl/random.hpp>` | C++23 | Aléatoire reproductible au bit près sur toute plateforme : moteurs `splitmix64`, `xoshiro256pp`, `pcg64`, `philox4x64` ; entiers et réels uniformes **exacts**, épreuve de Bernoulli exacte, mélange ; lois de probabilité en objets mathématiques (opérations en membres, `X.cdf(x)`, tirage `X(g)` ; paramètres nommés et exacts ; `uniform`, `discrete_uniform`, `bernoulli`, `exponential`, `normal`, `geometric`, `poisson`, `binomial`, `gamma`, `beta`, `chi_squared`, `student_t`, `lognormal`, `cauchy`, `fisher_f`, `weibull`, `laplace`, `logistic`, `pareto`, `gumbel`, `frechet`, `rayleigh`, `triangular`, `arcsine`, `erlang`, `inverse_gaussian`, `negative_binomial`, `hypergeometric`, `beta_binomial`, `logarithmic`, `discrete_gaussian`, `skellam`, `zipf`, `rice`, `skew_normal`, `stable` : moments, densité ou masse, répartition, quantile et quantile de survie sous 0,51 ulp ; tirage exact pour les cinq premières, `laplace` et `erlang`) ; algèbre des lois aux noms de la théorie de la mesure : `affine`, `convolution` (identités de familles), `truncate`, `mixture`, `product_measure`, `pushforward` ; lois multivariées (`dirichlet`, `multinomial`, `multivariate_normal`, `multivariate_student_t`, `wishart`, `inverse_wishart`, `uniform_sphere`, `uniform_ball`, `von_mises`, `von_mises_fisher`, dimension fixe ou `wsl::random::dynamic`) ; objets combinatoires à tirage exact (`uniform_permutation`, `ewens`, `mallows`, `uniform_subset`, `bernoulli_subset`, `uniform_combination`, `uniform_word`, `uniform_integer_partition`, `uniform_set_partition`, `uniform_composition`, `erdos_renyi`, `uniform_graph`, `barabasi_albert`, `random_regular_graph`, `uniform_tree`, `uniform_binary_tree`) ; processus à trajectoires paresseuses et exactes (`markov_chain`, `random_walk`, `polya_urn`, `chinese_restaurant`, `brownian_motion`, `brownian_bridge`, `ornstein_uhlenbeck`, `poisson_process`, `compound_poisson`, `continuous_markov_chain`…) et `categorical` exacte ; matrices aléatoires (`haar_orthogonal`, `haar_unitary`, `goe`, `gue`) ([conception](design/random.md)) |
| `wsl::mdarray` | `<wsl/mdarray.hpp>` | C++23 | Tableaux multidimensionnels (lot 1) : modes de calcul rapide et reproductible (`wsl::reproducible(…)` : mêmes bits sur toute plateforme, fma comprise, fonctions de `wsl::special`) ; calcul par les axes, le membre de gauche décide (`c[i, k] = a[i, j] * b[j, k]`, `y[i] = tanh(W[i, j] * x[j] + b[i])`, `double n = sqrt(x[i] * x[i])`, `A[i, i] = 1`, `r[i] = sum(j)(a[i, j])`, `h[i, j] = 1.0 / (i + j + 1)`, garde-fous à la compilation) ; calcul sur tableaux entiers (`y = x * x`, `r = hypot(x, y)`, `mdarray c = a * b + 1`, fonctions de `<cmath>`, `vectorize`, masques et `where`, réductions `sum`, `prod`, `min`, `max`, `argmin`, `argmax`, `all`, `any`, `count`) ; axes (`wsl::axis()`, `wsl::axis<"mu", "nu">()`, `wsl::axes::i`…) ; conteneur `mdarray` (P1684 : sémantique de valeur, `std::array` sans allocation pour des étendues fixes, littéral `{{1, 2}, {3, 4}}`, `equal`) ; `extents`, `dextents`, `dims`, `layout_right`, `layout_left`, `layout_stride`, `mdspan` (noms et sémantique de `<mdspan>`, sans GCC 16 ni libc++, conversions avec `std::mdspan` là où il existe ; ajouts : `m[p]`, `at`, `contains`, `wrap`, `clamp`) ; positions `mdposition` et déplacements `mddisplacement` (modèle affine, ordre lexicographique et composante par composante, normes, hachage, formatage) ; division euclidienne `div_euclid`, `rem_euclid` ([conception](design/mdarray.md)) |
| `wsl::libraries` | — | C++23 | Toutes les bibliothèques ci-dessus |

Les anciennes bibliothèques (`array2d`, `graded`, `graded_io`, `elo`, `progress`) attendent
leur refonte dans `legacy/`, figées ; `dlist` a quitté le dépôt. Refonte en cours :
voir design/interfaces.md (décisions 10 à 40).

Toutes les bibliothèques exigent C++23, transmis au projet consommateur : lier une
cible `wsl::…` à un projet en C++20 le fait passer en C++23. L'historique des
changements, bibliothèque par bibliothèque, est dans [CHANGELOG.md](CHANGELOG.md).

En préparation dans `wsl::random` : d'autres lois (Rice, zêta…)
([conception](design/random.md)).

## Utiliser les bibliothèques dans un projet

### FetchContent, figé sur une version (recommandé)

```cmake
include(FetchContent)
FetchContent_Declare(wsl_libraries
    GIT_REPOSITORY https://github.com/wills31415/cpp-libraries.git
    GIT_TAG        v0.1.0
    GIT_SHALLOW    TRUE
    EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(wsl_libraries)

target_link_libraries(mon_programme PRIVATE wsl::mdarray wsl::random)
```

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
        "cacheVariables": { "FETCHCONTENT_SOURCE_DIR_WSL_LIBRARIES": "/chemin/vers/cpp-libraries" }
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
find_package(wsl_libraries 0.1 REQUIRED)
target_link_libraries(mon_programme PRIVATE wsl::mdarray)
```

Pratique pour un essai rapide, mais une seule version existe alors pour toute la
machine : une version incompatible installée plus tard casse les anciens projets.

## Développement

| Action | VSCode | Terminal |
|---|---|---|
| Compiler | ⇧⌘B | `cmake --preset debug && cmake --build --preset debug -j` |
| Tests (ASan + UBSan) | tâche « Tests » | `ctest --preset debug` |
| Tests en Release | tâche « Tests (Release) » | `ctest --preset release` |
| Déboguer le test ouvert | F5 → « Déboguer le test courant » | — |

Chaque `tests/test_<nom>.cpp` donne `test_<nom>`, compilé en C++23 (et, si une
bibliothèque redescendait un jour plus bas, une fois par standard supporté :
`test_<nom>_cxx20`…). Un test réussit s'il renvoie 0. `tests/package/` vérifie, en
CI, l'installation et FetchContent du point de vue d'un consommateur.

### Ajouter une bibliothèque

1. `include/wsl/<nom>.hpp`, dans le namespace `wsl`.
2. `tests/test_<nom>.cpp`.
3. Dans `CMakeLists.txt` : `wsl_add_library(<nom> STANDARD 23 [DEPENDS ...])`
   et `wsl_add_test(<nom>)`.
4. Une ligne dans le catalogue ci-dessus, une entrée dans `CHANGELOG.md`.

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

## Licence

[MIT](LICENSE).
