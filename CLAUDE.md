# cpp-libraries

Grande bibliothèque scientifique C++23 en en-têtes seuls (mathématiques, physique, calcul
scientifique, et leurs fondations : mdarray, hasard, fonctions spéciales…), espace de noms
`wsl` (**nom provisoire**, à remplacer avant la première version stable), paquet CMake
`wsl_libraries` (cibles `wsl::<nom>`) jusqu'à la refonte. Identité : docs/interfaces.md,
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
- Standard : **C++23 partout** (décision 2026-11-20, docs/interfaces.md) ; C++20 seulement
  sur décision explicite, pour un besoin de compatibilité identifié. Le déclarer dans
  `wsl_add_library(... STANDARD ...)`.
- Interfaces : les décisions de la revue (espaces de noms, paramètres, types, erreurs…) sont
  dans `docs/interfaces.md` ; tout nouveau code les suit.
- **Charte** (docs/interfaces.md §6, décisions 10 à 36) : reproductibilité (10), espaces de
  noms, bibliothèques, abréviations et alias (11 à 23), modes `deterministic`/`portable` et
  `approximate`/`exact`, le défaut étant la permission (24 à 27), registre provisoire de `wsl::`
  (28), chemins, en-têtes sans extension, cibles, organisation du dépôt, paquet (29 à 33),
  identité du dépôt et sort des anciennes bibliothèques (34 à 36). Tout nouveau code la suit ;
  les sections ci-dessous marquées « jusqu'à la refonte » décrivent l'état actuel.
- **Priorités de conception** : simple pour l'utilisateur, compatibilités vérifiées à la
  compilation, options les plus rapides, reproductibilité.
- Objets mathématiques (algèbre linéaire, tenseurs, analyse) : vision, méthode et feuille
  de route dans `docs/algebra.md`. Ne citer que des ouvrages publiés, jamais de supports de
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
  (macOS), Debug (-O1, ASan + UBSan) + Release, configurations FMA en Release seulement, `-DWSL_LIBRARIES_WERROR=ON`, plus `tests/package/`.
  Reproduire la partie Linux en local : conteneur `ubuntu:24.04` avec
  `cmake make g++-14 clang-18 libclang-rt-18-dev` (Docker est installé).
  Y compiler les tests **optimisés** (`-O2`, comme les jobs release) : certains
  avertissements de GCC (`-Warray-bounds`…) n'apparaissent qu'avec l'optimisation ;
  `-fsyntax-only` ne suffit pas.
- La CI ne compile et ne teste que ce qu'un push peut avoir changé : chemins modifiés
  depuis le dernier passage réussi (`tools/ci/changed_files.sh`), bibliothèques touchées
  et celles qui en dépendent (`cmake/affected.cmake`, option `WSL_LIBRARIES_AFFECTED_BY`) ;
  documentation seule (`cmake/ci_ignored_paths.txt`) : rien ; infrastructure, tag,
  lancement manuel, passage hebdomadaire : tout. Une nouvelle bibliothèque est reconnue
  par ses chemins (`include/wsl/<nom>…`, `tools/<nom>/`, `tests/test_<nom>…`).
- Les projets consomment le dépôt par FetchContent figé sur un tag ; cpp-init le
  configure par défaut.

## Commandes

- Compiler : `cmake --preset debug && cmake --build --preset debug -j`
- Tests : `ctest --preset debug` (ASan + UBSan) ; `ctest --preset release` (8 tests en
  parallèle ; la suite complète en séquentiel : près de 3 min en release, bien plus en
  debug, où les tests de `wsl::random` prennent 30 à 90 s chacun)
- Vérification rapide : `tools/check.sh [release|fma|debug]…` (par défaut release et fma)
  ne compile et ne teste que les bibliothèques touchées depuis le dernier push (même règle
  que la CI, `cmake/affected.cmake`), dans `build/check-<configuration>` ; `--all` pour
  tout. La suite debug complète est laissée à la CI.
- Objectif : 0 avertissement avec `-Wall -Wextra -Wpedantic`, en C++23,
  sur les trois compilateurs de la CI.

## Outils (jusqu'à la refonte)

- `tools/cpp-init` : **source unique** du script de création de projets ;
  `~/.local/bin/cpp-init` est un lien vers lui, et le skill `cpp-init`
  (`~/.claude/skills/cpp-init/`) le lance. Le modifier ici.
- `tools/release.sh` : publication d'une version.

## wsl::special (fonctions bit-exactes ; jusqu'à la refonte)

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

## wsl::random (jusqu'à la refonte)

- Conception : `docs/random.md`. Mêmes exigences que wsl::special : algorithmes fixés,
  tirages bit-exacts partout ; le nombre de valeurs consommées sur le moteur fait partie
  du contrat. Jamais de loi de `<random>` (non reproductibles d'une bibliothèque à
  l'autre) ; les moteurs standard sont acceptés.
- Vecteurs des moteurs : `uv run --no-project --with numpy python
  tools/random/engine_vectors.py` écrit `tests/random_engine_vectors.hpp` (versionné,
  petit ; NumPy pour pcg64 et philox4x64, sauts de xoshiro vérifiés par GF(2)).
- Empreintes des suites tirées : `tests/test_random_uniform.cpp`, `tests/test_random_laws.cpp`.
- Précision des lois : références mpmath non versionnées, engendrées par CMake, **une
  famille par module** de `tools/random/references/` (`common.py` : outils et écriture ;
  la famille `base` vérifie aussi les constantes de `wsl/special/detail/extended.hpp`).
  Seules les familles modifiées sont recalculées : une nouvelle loi a son propre module
  (déclaré dans `law_families` du CMakeLists). Mesurer le temps d'un module seul
  (`law_references.py FAMILLE FICHIER`) avant la génération complète ; éviter les
  fonctions spéciales de mpmath dans les boucles (récurrences, sommes en précision fixée) ;
  tables de la ziggourat versionnées
  (`uv run --no-project --with mpmath python tools/random/ziggurat.py`).
- Tirages exacts (`detail/lazy.hpp`) : toute voie rapide ne doit rendre un résultat que
  s'il est certain, pour que résultat et nombre de mots tirés soient ceux de la voie
  exacte ; ses décisions dépendent de flottants : produits par `fmul` là aussi.

## Structure (jusqu'à la refonte)

- `include/wsl/<nom>.hpp` : une bibliothèque par en-tête (dépendances internes via
  `DEPENDS`, p. ex. `graded_io` → `graded`).
- `tests/test_<nom>.cpp` : `test_<nom>`, en C++23 (une variante `_cxx20` par standard
  inférieur supporté, aucune aujourd'hui).
- `examples/*.cpp` : un exécutable chacun, lié à `wsl::libraries`.
- Compilateur imposé : `/usr/bin/clang++` (miniconda est en tête du PATH).

## Provenance

Les originaux sont restés en place, non modifiés : `~/projets/cpp-array2d`,
`~/projets/cpp-container-distance-list`, `cpp-sudoku/src/exact_random.hpp` (retirée
d'ici depuis, remplacée par `wsl::random`),
`cpp-traversal/src/util/{array2d,progress}.hpp`, `cpp-traversal/src/agents/elo.hpp`.
cpp-traversal utilise encore sa propre copie d'`array2d.hpp`.

Les autres dépôts de l'utilisateur sont des projets qui peuvent consommer la bibliothèque.
Les noms s'y choisissent librement : le paquet `wsl` de `~/projets/cpp-tensors` (même espace, même
dossier d'en-têtes) est un premier essai, qui s'adaptera en étant refondu sur les objets
mathématiques d'ici (docs/algebra.md). Aucun projet ne combine aujourd'hui les deux.

## Refonte à faire (décidée, pas encore exécutée)

- Organisation (décision 32) : `include/wsl/`, `docs/` (guides, `README.md` par dossier),
  `examples/<domaine>/`, `tests/` (sans préfixe, nom dérivé du chemin, dossier dès deux
  batteries, étiquettes de nature et `long`), `tools/`, `design/` (conception et charte).
- Chemins d'après les espaces de noms (29) ; relais sans extension engendrés par CMake, après
  un prototype qui étend `tests/package/` (30) ; cibles nommées par l'espace de noms,
  `wsl::wsl` pour tout, contrôle des inclusions en CI (31).
- Paquet `find_package(wsl)`, options et macros `WSL_*` (`WSL_PORTABLE`, `WSL_EXACT`,
  `WSL_WERROR`, `WSL_AFFECTED_BY`) (33).
- Modes (24 à 27) : `wsl::fast` et `wsl::random::fast` disparaissent, `wsl::reproducible`
  devient `wsl::portable`, `approximate` par défaut, `wsl::provides`, `wsl::assume_portable`,
  espace de noms en ligne du défaut.
- Découpage de special et random selon les domaines (`floating_point_arithmetic`,
  `elementary_functions`, `special_functions`, `random`, `probability`, `measure_theory`) ; carte : décision 37 ;
  audit de chaque `fmul` en `fpa::mul` (selon le mode) ou `fpa::unfused_mul` (justesse) : décision 38.
- Anciennes bibliothèques dans `legacy/` (array2d, graded, graded_io, elo, progress) ; `dlist`
  rendue à `~/projets/cpp-container-distance-list` (dossier simple, pas un dépôt git :
  regarder avant de remplacer) (36).
- cpp-init vers un dépôt privé à lui, sur confirmation (35) ; README réécrit à la refonte.
