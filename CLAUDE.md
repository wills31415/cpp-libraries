# cpp-libraries

Grande bibliothèque scientifique C++23 en en-têtes seuls (mathématiques, physique, calcul
scientifique, et leurs fondations : mdarray, hasard, fonctions spéciales…), espace de noms
`wsl` (**nom provisoire**, à remplacer avant la première version stable), paquet CMake
`wsl` (cibles `wsl::<espace de noms>`, `wsl::wsl` pour tout ; `v0.1.0` : `wsl_libraries`). Identité : design/interfaces.md,
décision 34 ; ce qui n'en relève pas sort du dépôt ou y devient un domaine.

## Règles

- **Claude implémente** (demande de l'utilisateur, 2026-11-16) : écrire, corriger et
  faire évoluer le code des bibliothèques et de leurs tests est attendu.
- **Objectif** : une bibliothèque cohérente, aussi optimisée, générique et proche de la
  bibliothèque standard que possible (nommage, garanties,
  `constexpr`/`noexcept`, concepts, comportement aligné sur les conteneurs et
  algorithmes standard).
- **Langue** : identifiants, noms de fichiers et de cibles en anglais ; commentaires
  et messages en français.
- Standard : **C++23 partout** (décision 2026-11-20, design/interfaces.md) ; C++20 seulement
  sur décision explicite, pour un besoin de compatibilité identifié. Le déclarer dans
  `wsl_add_library(... STANDARD ...)`.
- Interfaces : les décisions de la revue (espaces de noms, paramètres, types, erreurs…) sont
  dans `design/interfaces.md` ; tout nouveau code les suit.
- **Charte** (design/interfaces.md §6, décisions 10 à 36) : reproductibilité (10), espaces de
  noms, bibliothèques, abréviations et alias (11 à 23), modes `deterministic`/`portable` et
  `approximate`/`exact`, le défaut étant la permission (24 à 27), registre provisoire de `wsl::`
  (28), chemins, en-têtes sans extension, cibles, organisation du dépôt, paquet (29 à 33),
  identité du dépôt et sort des anciennes bibliothèques (34 à 36), découpage, barrières,
  mesures, processus, exactitude et lois (37 à 42). Tout nouveau code la suit.
- **Priorités de conception** : simple pour l'utilisateur, compatibilités vérifiées à la
  compilation, options les plus rapides, reproductibilité.
- Objets mathématiques (algèbre linéaire, tenseurs, analyse) : vision, méthode et feuille
  de route dans `design/algebra.md`. Ne citer que des ouvrages publiés, jamais de supports de
  cours (ni formation, ni école, ni enseignants).
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
  (macOS), Debug (-O1, ASan + UBSan) + Release, configurations FMA en Release seulement, `-DWSL_WERROR=ON`, plus `tests/package/`.
  Reproduire la partie Linux en local : conteneur `ubuntu:24.04` avec
  `cmake make g++-14 clang-18 libclang-rt-18-dev` (Docker est installé).
  Y compiler les tests **optimisés** (`-O2`, comme les jobs release) : certains
  avertissements de GCC (`-Warray-bounds`…) n'apparaissent qu'avec l'optimisation ;
  `-fsyntax-only` ne suffit pas.
- La CI ne compile et ne teste que ce qu'un push peut avoir changé : chemins modifiés
  depuis le dernier passage réussi (`tools/ci/changed_files.sh`), bibliothèques touchées
  et celles qui en dépendent (`cmake/affected.cmake`, option `WSL_AFFECTED_BY`) ;
  documentation seule (`cmake/ci_ignored_paths.txt`) : rien ; infrastructure, tag,
  lancement manuel, passage hebdomadaire : tout. Une nouvelle bibliothèque est reconnue
  par ses chemins (`include/wsl/<nom>…`, `tools/<nom>/`, `tests/test_<nom>…`).
- Les projets consomment le dépôt par FetchContent figé sur un tag ; l'outil personnel
  cpp-init (dépôt privé `wills31415/cpp-init`, `~/projets/cpp-init`) le configure par défaut.

## Commandes

- Compiler : `cmake --preset debug && cmake --build --preset debug -j`
- Tests : `ctest --preset debug` (ASan + UBSan, sans les tests étiquetés `long`) ;
  `ctest --preset release` (8 tests en parallèle, tous) ; `ctest -L precision`, `-R mdarray`…
  Noms et étiquettes : tests rangés comme les espaces de noms (`tests/mdarray/axis.cpp` →
  `mdarray.axis`), étiquettes en tête du fichier (`// Étiquettes : precision long`).
- Vérification rapide : `tools/check.sh [release|fma|debug]…` (par défaut release et fma)
  ne compile et ne teste que les bibliothèques touchées depuis le dernier push (même règle
  que la CI, `cmake/affected.cmake`), dans `build/check-<configuration>` ; `--all` pour
  tout. La suite debug complète est laissée à la CI.
- Objectif : 0 avertissement avec `-Wall -Wextra -Wpedantic`, en C++23,
  sur les trois compilateurs de la CI.

## Outils

- `tools/release.sh` : publication d'une version ; `tools/check.sh` : vérification rapide.
- `tools/elementary_functions/` (tables, références, mesures), `tools/random/` (vecteurs des
  moteurs, ziggourat), `tools/probability/` (références des lois), `tools/mdarray/` (mesures).

## Fonctions : floating_point_arithmetic, elementary_functions, special_functions

- Conception et mesures : `design/special.md` ; guides : `docs/<domaine>/README.md`.
- Les lois et le socle double-double (`special_functions/detail/{extended,loader}.hpp`)
  n'ont qu'une version, portable (décision 42) : ils appellent les fonctions figées en mode
  portable (`detail::pinned_portable` de ef et sf, `random::detail::pm`), jamais `ef::exp(x)`
  sans étiquette (bascule sans bruit en déterministe, attrapée par
  `special_functions.no_contraction`).
- **Règle d'écriture** : uniquement `+ − × ÷` (et `sqrt`, `fma` explicites), jamais de
  fonction transcendante de `<cmath>`. **Tout produit qui entre dans une addition ou une
  soustraction passe par `fpa::mul<R>`** (barrière en mode portable, fusion permise en mode
  déterministe) **ou, si l'algorithme s'appuie sur sa valeur arrondie** (argument direct d'une
  somme exacte, d'un découpage de Veltkamp, d'une astuce d'arrondi), **par
  `fpa::unfused_mul`** (décision 38). Chaque fonction s'écrit une fois, le mode `R` en
  paramètre ; tout appel interne porte son `<R>` (un appel qui l'oublie ne compile pas,
  plutôt que de basculer sans bruit en déterministe). `fmul` reste l'ancien nom de
  `unfused_mul`, dans le socle double-double et random, pas encore convertis.
- Tables engendrées et versionnées : `uv run --no-project --with mpmath python
  tools/elementary_functions/generate.py --tables` (ne jamais modifier `elementary_functions/detail/*tables.hpp` ni
  `special_functions/detail/*tables.hpp` à la
  main). Valeurs de référence des tests : **non versionnées**, engendrées par CMake
  dans `build/<preset>/generated/` (uv, ou python3 + mpmath ; sinon le test de
  précision est sauté), et seulement quand `generate.py` change.
- `tests/special_functions/fingerprints.cpp` contient les empreintes attendues : un changement
  voulu d'algorithme les change partout (rupture, à noter au CHANGELOG) ; un
  changement sur une seule plateforme est un bug.
- Les entrées des tests doivent elles aussi être protégées de la contraction
  (`detail::fmul`), sinon ce sont les entrées qui changent sous `-ffp-contract=fast`.
- Le test `special_functions.no_contraction` (cmake/check_no_contraction.cmake) refuse toute fma
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
- Mesurer les performances : `tools/elementary_functions/benchmark.sh` (un programme par fonction ;
  ne pas mesurer plusieurs fonctions dans un même programme, les chiffres s'y faussent).
- Reproduire la CI Linux, configurations FMA comprises, en local : conteneur
  `ubuntu:24.04` (voir plus haut), avec `-DCMAKE_CXX_FLAGS="-march=x86-64-v3 -ffp-contract=fast"`.

## random, probability

- Conception : `design/random.md`. Mêmes exigences que les fonctions : algorithmes fixés,
  tirages bit-exacts partout ; le nombre de valeurs consommées sur le moteur fait partie
  du contrat. Jamais de loi de `<random>` (non reproductibles d'une bibliothèque à
  l'autre) ; les moteurs standard sont acceptés.
- Vecteurs des moteurs : `uv run --no-project --with numpy python
  tools/random/engine_vectors.py` écrit `tests/random/engine_vectors.hpp` (versionné,
  petit ; NumPy pour pcg64 et philox4x64, sauts de xoshiro vérifiés par GF(2)).
- Empreintes des suites tirées : `tests/random/uniform.cpp`, `tests/probability/laws.cpp`.
- Précision des lois : références mpmath non versionnées, engendrées par CMake, **une
  famille par module** de `tools/probability/references/` (`common.py` : outils et écriture ;
  la famille `base` vérifie aussi les constantes de `wsl/special_functions/detail/extended.hpp`).
  Seules les familles modifiées sont recalculées : une nouvelle loi a son propre module
  (déclaré dans `law_families` du CMakeLists). Mesurer le temps d'un module seul
  (`law_references.py FAMILLE FICHIER`) avant la génération complète ; éviter les
  fonctions spéciales de mpmath dans les boucles (récurrences, sommes en précision fixée) ;
  tables de la ziggourat versionnées
  (`uv run --no-project --with mpmath python tools/random/ziggurat.py`).
- Tirages exacts (`detail/lazy.hpp`) : toute voie rapide ne doit rendre un résultat que
  s'il est certain, pour que résultat et nombre de mots tirés soient ceux de la voie
  exacte ; ses décisions dépendent de flottants : produits par `fmul` là aussi.
- Tirage par défaut **approché** (`X(g)` : variante approchée si la loi en a une) ;
  `X(g, wsl::exact)` refusé sans tirage exact (`L::sampling`). **Tout tirage interne d'une loi
  ou d'une variable** passe par `detail::native_sample(X, g)` (le tirage natif), sinon un
  tirage exact tire ses composantes en approché sans bruit. Les tests relèvent les empreintes
  sur le tirage natif (`wr::native`).

## Structure

- `include/wsl/<chemin>.hpp` : le chemin est l'espace de noms (décision 29) ; une
  bibliothèque par espace de noms de domaine, `wsl_add_library(<nom avec ::> … DEPENDS …)`.
- `tests/<chemin>.cpp` : test « chemin.avec.points », en C++23 (décision 32).
- `docs/<chemin>/README.md` : guide de chaque domaine (en français, langue à décider) ;
  `docs/README.md` : index, conventions, modes.
- `examples/<domaine>/<tâche>.cpp` : un exécutable chacun, lié à `wsl::wsl`, lancé comme
  test `example.<domaine>.<tâche>` (étiquette `example`) ; cité par la documentation.
- Compilateur imposé : `/usr/bin/clang++` (miniconda est en tête du PATH).

## Provenance

Les anciennes bibliothèques sont dans `legacy/` (figées, décision 36) ; `dlist` a été rendue à
`~/projets/cpp-container-distance-list/from-cpp-libraries/` (2026-11-28).
Les originaux sont restés en place, non modifiés : `~/projets/cpp-array2d`,
`~/projets/cpp-container-distance-list`, `cpp-sudoku/src/exact_random.hpp` (retirée
d'ici depuis, remplacée par `wsl::random`),
`cpp-traversal/src/util/{array2d,progress}.hpp`, `cpp-traversal/src/agents/elo.hpp`.
cpp-traversal utilise encore sa propre copie d'`array2d.hpp`.

Les autres dépôts de l'utilisateur sont des projets qui peuvent consommer la bibliothèque.
Les noms s'y choisissent librement : le paquet `wsl` de `~/projets/cpp-tensors` (même espace, même
dossier d'en-têtes) est un premier essai, qui s'adaptera en étant refondu sur les objets
mathématiques d'ici (design/algebra.md). Aucun projet ne combine aujourd'hui les deux.

## Refonte (faite, 2026-11-28) et feuille de route

La refonte décidée par la charte est exécutée : rangement (`legacy/`, `design/`, cpp-init à
part), infrastructure (paquet `wsl`, cibles nommées par l'espace de noms, relais sans
extension, tests rangés et étiquetés, contrôle des inclusions), modes (`<wsl/modes>`, ef et
sf dans les deux modes, `WSL_EXACT` strict, tirage approché par défaut), découpage de
l'ancien `wsl::special` et de l'ancien `wsl::random` en domaines (décision 37), transitions
retirées, documentation (`docs/`, `examples/`). cpp-init garde `wsl_libraries` tant que le
dernier tag est `v0.1.0` (à changer au prochain tag).

Ensuite, dans l'ordre :
1. Fonctions de `ef` correctement arrondies (`wsl::exact`, méthode de Ziv, décision 41) :
   exp, exp2, expm1, puis logarithmes, trigonométrie, cbrt et hypot, pow ; puis revoir
   `WSL_EXACT` (strict, ou séparé en tirages exacts et arrondi correct).
2. Lois et socle double-double dans les deux modes (règle D, décision 38), algèbre des lois
   dans `measure_theory` avec sa généralisation aux mesures (décisions 39 et 42).
3. Interface publique de `floating_point_arithmetic` ; référence engendrée depuis les
   en-têtes (outil à choisir) ; langue de la documentation (décision 32).
4. Lot 2 de mdarray (vues), puis les objets mathématiques de `design/algebra.md`.
