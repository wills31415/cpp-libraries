# cpp-libraries

Collection centrale des bibliothèques C++ génériques de l'utilisateur, en en-têtes
seuls, namespace `wsl`, paquet CMake `wsl_libraries` (cibles `wsl::<nom>`).

## Règles

- **Claude implémente** (demande de l'utilisateur, 2026-11-16) : écrire, corriger et
  faire évoluer le code des bibliothèques et de leurs tests est attendu. Exception à
  la règle par défaut de `cpp-init`.
- **Objectif** : centraliser le maximum de code générique, aussi optimisé, générique
  et proche de la bibliothèque standard que possible (nommage, garanties,
  `constexpr`/`noexcept`, concepts, comportement aligné sur les conteneurs et
  algorithmes standard).
- **Langue** : identifiants, noms de fichiers et de cibles en anglais ; commentaires
  et messages en français.
- Standard minimal par bibliothèque : C++20 si possible, C++23 seulement si une
  fonctionnalité l'exige. Le déclarer dans `wsl_add_library(... STANDARD ...)`.
- Aucune dépendance externe dans une bibliothèque, sauf décision explicite.

## Dépôt, versions, CI

- Dépôt **public** `github.com/wills31415/cpp-libraries`, licence MIT : rien de
  personnel dans les fichiers suivis (pas de chemin absolu `/Users/...`, pas de
  détail sur le homelab ou les projets privés au-delà de leur nom).
- Toute modification d'une bibliothèque s'accompagne d'une entrée dans la section
  `[Non publié]` de `CHANGELOG.md`, sous le nom de la bibliothèque.
- Publication : `tools/release.sh X.Y.Z` (voir README). Ne publier, pousser ou
  taguer qu'à la demande de l'utilisateur.
- CI (`.github/workflows/ci.yml`) : GCC 14 et Clang 18 (Ubuntu 24.04) et Apple Clang
  (macOS), Debug + Release, `-DWSL_LIBRARIES_WERROR=ON`, plus `tests/package/`.
  Reproduire la partie Linux en local : conteneur `ubuntu:24.04` avec
  `cmake make g++-14 clang-18 libclang-rt-18-dev` (Docker est installé).
  Y compiler les tests **optimisés** (`-O2`, comme les jobs release) : certains
  avertissements de GCC (`-Warray-bounds`…) n'apparaissent qu'avec l'optimisation ;
  `-fsyntax-only` ne suffit pas.
- Les projets consomment le dépôt par FetchContent figé sur un tag ; cpp-init le
  configure par défaut.

## Commandes

- Compiler : `cmake --preset debug && cmake --build --preset debug -j`
- Tests : `ctest --preset debug` (ASan + UBSan) ; `ctest --preset release`
- Objectif : 0 avertissement avec `-Wall -Wextra -Wpedantic`, en C++20 et C++23,
  sur les trois compilateurs de la CI.

## Outils

- `tools/cpp-init` : **source unique** du script de création de projets ;
  `~/.local/bin/cpp-init` est un lien vers lui, et le skill `cpp-init`
  (`~/.claude/skills/cpp-init/`) le lance. Le modifier ici.
- `tools/release.sh` : publication d'une version.

## wsl::special (fonctions bit-exactes)

- Conception et mesures : `docs/special.md` ; lois de probabilité à venir : `docs/random.md`.
- **Règle d'écriture** : uniquement `+ − × ÷` (et `sqrt`, `fma` explicites), jamais de
  fonction transcendante de `<cmath>`. **Tout produit qui entre dans une addition ou
  une soustraction passe par `detail::fmul`** (barrière contre la contraction en fma,
  qui changerait les résultats selon la plateforme et les options du consommateur).
- Tables engendrées et versionnées : `uv run --no-project --with mpmath python
  tools/special/generate.py --tables` (ne jamais modifier `detail/*tables.hpp` à la
  main). Valeurs de référence des tests : **non versionnées**, engendrées par CMake
  dans `build/<preset>/generated/` (uv, ou python3 + mpmath ; sinon le test de
  précision est sauté), et seulement quand `generate.py` change.
- `tests/test_special_bitexact.cpp` contient les empreintes attendues : un changement
  voulu d'algorithme les change partout (rupture, à noter au CHANGELOG) ; un
  changement sur une seule plateforme est un bug.
- Les entrées des tests doivent elles aussi être protégées de la contraction
  (`detail::fmul`), sinon ce sont les entrées qui changent sous `-ffp-contract=fast`.
- Le test `special_no_contraction` (cmake/check_no_contraction.cmake) refuse toute fma
  dans le code engendré ; il ne voit pas un produit non protégé à usages multiples.
  Il essaie aussi `1 + f(x)` pour chaque fonction : GCC fusionne, après intégration, un
  produit **rendu** (`return y * 2^e`) ou un `x + x` (chemins NaN) avec l'addition de
  l'appelant ; d'où `fmul` sur les produits rendus et `opaque(x + x)`. Le reproduire
  avec GCC aussi (conteneur), Apple Clang ne le fait pas.
- `two_product` : fma(a, b, −p) quand la cible a une fma (`detail::hardware_fma`),
  Dekker sinon et en évaluation constante ; mêmes bits **tant que Dekker est exact** :
  facteurs sous 2^995 (sinon le découpage déborde en NaN) et somme des exposants
  >= −970 (pas de produit partiel sous-normal). Vérifier ces bornes à chaque nouvel
  appel (garder sinon) ; les empreintes comparent les deux voies (CI : ARM64 et jobs
  FMA contre x86-64 par défaut ; en local : un build `-march=native`).
- Tests : ne jamais exiger l'égalité avec la libm du système (glibc n'est pas
  correctement arrondie partout) ; comparer à des valeurs exactes (mpmath) ou à 1 ulp.
- Générateur : toute constante qui demande plus de 256 bits (2/π de Payne–Hanek…)
  doit être calculée sous `mp.workprec(...)` ; les racines complexes de mpmath
  (`mp.cbrt` d'un négatif) ne sont pas les racines réelles.
- Mesurer les performances : `tools/special/benchmark.sh` (un programme par fonction ;
  ne pas mesurer plusieurs fonctions dans un même programme, les chiffres s'y faussent).
- Reproduire la CI Linux, configurations FMA comprises, en local : conteneur
  `ubuntu:24.04` (voir plus haut), avec `-DCMAKE_CXX_FLAGS="-march=x86-64-v3 -ffp-contract=fast"`.

## wsl::random

- Conception : `docs/random.md`. Mêmes exigences que wsl::special : algorithmes fixés,
  tirages bit-exacts partout ; le nombre de valeurs consommées sur le moteur fait partie
  du contrat. Jamais de loi de `<random>` (non reproductibles d'une bibliothèque à
  l'autre) ; les moteurs standard sont acceptés.
- Vecteurs des moteurs : `uv run --no-project --with numpy python
  tools/random/engine_vectors.py` écrit `tests/random_engine_vectors.hpp` (versionné,
  petit ; NumPy pour pcg64 et philox4x64, sauts de xoshiro vérifiés par GF(2)).
- Empreintes des suites tirées : `tests/test_random_uniform.cpp`, `tests/test_random_laws.cpp`.
- Précision des lois : références mpmath non versionnées, engendrées par CMake
  (`tools/random/law_references.py`, qui vérifie aussi les constantes de
  `wsl/special/detail/extended.hpp`) ; tables de la ziggourat versionnées
  (`uv run --no-project --with mpmath python tools/random/ziggurat.py`).
- Tirages exacts (`detail/lazy.hpp`) : toute voie rapide ne doit rendre un résultat que
  s'il est certain, pour que résultat et nombre de mots tirés soient ceux de la voie
  exacte ; ses décisions dépendent de flottants : produits par `fmul` là aussi.

## Structure

- `include/wsl/<nom>.hpp` : une bibliothèque par en-tête (dépendances internes via
  `DEPENDS`, p. ex. `graded_io` → `graded`).
- `tests/test_<nom>.cpp` : compilé une fois par standard supporté (`test_<nom>`
  en C++23, `test_<nom>_cxx20`).
- `examples/*.cpp` : un exécutable chacun, lié à `wsl::libraries`.
- Compilateur imposé : `/usr/bin/clang++` (miniconda est en tête du PATH).

## Provenance

Les originaux sont restés en place, non modifiés : `~/projets/cpp-array2d`,
`~/projets/cpp-container-distance-list`, `cpp-sudoku/src/exact_random.hpp` (retirée
d'ici depuis, remplacée par `wsl::random`),
`cpp-traversal/src/util/{array2d,progress}.hpp`, `cpp-traversal/src/agents/elo.hpp`.
cpp-traversal utilise encore sa propre copie d'`array2d.hpp`.

Le paquet `wsl` de `~/projets/cpp-tensors` (cibles `wsl::wsl`, `wsl::formula`,
en-têtes `include/wsl/`) partage le namespace et le dossier d'en-têtes : éviter les
noms d'en-têtes qui y existent déjà (`index.hpp`, `space.hpp`, `field.hpp`…).
