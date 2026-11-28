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
  identité du dépôt et sort des anciennes bibliothèques (34 à 36). Tout nouveau code la suit ;
  les sections ci-dessous marquées « jusqu'à la refonte » décrivent l'état actuel.
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

## Outils (jusqu'à la refonte)

- `tools/release.sh` : publication d'une version.

## wsl::special (fonctions bit-exactes ; jusqu'à la refonte)

- Conception et mesures : `design/special.md` ; lois de probabilité à venir : `design/random.md`.
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

## wsl::random (jusqu'à la refonte)

- Conception : `design/random.md`. Mêmes exigences que wsl::special : algorithmes fixés,
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

## Structure

- `include/wsl/<chemin>.hpp` : le chemin est l'espace de noms (décision 29) ; une
  bibliothèque par espace de noms de domaine, `wsl_add_library(<nom avec ::> … DEPENDS …)`.
- `tests/<chemin>.cpp` : test « chemin.avec.points », en C++23 (décision 32).
- `examples/*.cpp` : un exécutable chacun, lié à `wsl::wsl`.
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

## Refonte à faire (décidée, pas encore exécutée)

Ordre (validé le 2026-11-28, exécuté en autonomie) : 0 prototype (relais sans extension,
cibles ; arrêt s'il échoue) ; 1 rangement (`legacy/`, `dlist`, `docs/` → `design/`, cpp-init :
confirmation) ; 2 infrastructure (cibles, paquet, options, tests, CI : accord) ; 3 modes ;
4 special ; 5 random ; 6 mdarray ; 7 documentation, README, CLAUDE.md ; puis chantier des
fonctions de `ef` correctement arrondies (décision 41) ; conversion des lois aux deux modes et
algèbre des lois dans `measure_theory` (décision 42). Exécution des étapes
4 à 6 en deux passes (2026-11-28) : d'abord le mécanique, sans changer un bit (4a special
découpée, 5a random découpée, 6a mdarray sur fpa/ef/sf, puis retrait de `<wsl/special>`, la
transition) ; ensuite le délicat (4b audit des `fmul` et modes de ef et sf, 5b modes de random
et `approximate` par défaut), fonction par fonction, bornes vérifiées dans les deux modes.

- Organisation (décision 32) : `include/wsl/`, `docs/` (guides, `README.md` par dossier),
  `examples/<domaine>/`, `tests/` (sans préfixe, nom dérivé du chemin, dossier dès deux
  batteries, étiquettes de nature et `long`), `tools/`, `design/` (conception et charte).
- Chemins d'après les espaces de noms (29) ; relais sans extension engendrés par CMake, après
  un prototype qui étend `tests/package/` (30) ; cibles nommées par l'espace de noms,
  `wsl::wsl` pour tout, contrôle des inclusions en CI (31).
- ~~Paquet `find_package(wsl)`, options `WSL_*`, cibles, tests rangés (permutation, mdarray),
  étiquettes, contrôle des inclusions~~ : fait (étape 2) ; `WSL_PORTABLE` et `WSL_EXACT` viennent
  avec l'étape 3. cpp-init garde `wsl_libraries` tant que le dernier tag est `v0.1.0`.
- ~~Modes (24 à 27)~~ : faits (étapes 3, 4b, 5b partie 1 : `approximate` par défaut pour les
  tirages, `wsl::random::fast` retiré). Les lois restent portables (axe de la
  reproductibilité accepté sans effet) en attendant leur conversion (5b partie 2, à décider).
- ~~Découpage mécanique de special (4a) et random (5a)~~ : fait. random : `wsl::random`
  (moteurs `random::engines`, primitives, tirages exacts), `wsl::probability` (lois
  `probability::laws`, variables, algèbre des lois provisoirement), `wsl::probability::
  stochastic_processes` (`processes`, `trajectory`), `wsl::measure_theory` (mesures de
  référence) ; transitions `<wsl/special>` et `probability/detail/transition.hpp` : retirées
  (les lois appellent `random::detail::pm`, fonctions figées en mode portable, décision 42).
- Découpage de special et random selon les domaines (`floating_point_arithmetic`,
  `elementary_functions`, `special_functions`, `random`, `probability`, `measure_theory`) ; carte : décision 37 ;
  audit de chaque `fmul` en `fpa::mul` (selon le mode) ou `fpa::unfused_mul` (justesse) : décision 38.
- ~~Anciennes bibliothèques dans `legacy/`, `dlist` rendue~~ : fait (étape 1).
- ~~cpp-init vers un dépôt privé à lui~~ : fait (étape 1) ; il faudra y changer le nom du
  paquet et la cible quand l'étape 2 les renommera. README réécrit à la refonte.
