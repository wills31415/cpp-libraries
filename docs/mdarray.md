# wsl::mdarray : tableaux multidimensionnels

Chantier ouvert le 2026-11-24, préalable aux objets mathématiques (docs/algebra.md), qui
s'en serviront pour stocker des composantes. La bibliothèque est une infrastructure
générale, indépendante de ces objets : elle sert aussi bien des images, des grilles ou des
plateaux de jeu. `array2d` en deviendra un alias de rang 2.

Livraison par lots (en-têtes séparés, une seule cible), ordonnés selon le but : une syntaxe
simple et lisible, et des opérations vectorisées (réordonnés le 2026-11-24) :

1. étendues (I2), tableau et vue, axes (§5), expressions (Einstein, diffusion, scalaires),
   calcul en boucles vectorisables : `c[i, k] = a[i, j] * b[j, k]`, `y[i] = 2 * x[i] + 1` ;
   positions et déplacements (§4) ;
2. vues de découpage et de décalage, parcours (I3) ;
3. exécution parallèle et GPU (I8).

État des décisions : ✅ validée ; ◻︎ proposée, en attente de validation.

## 1. Objet

| Objet | En mathématiques | Type |
|---|---|---|
| indice entier | une position i ∈ {0, …, n − 1} d'une dimension | `std::size_t` (le type d'indice des étendues) |
| position | un point p ∈ ℕᴺ de l'espace affine ℤᴺ | `wsl::mdsize<N>` |
| déplacement | un vecteur d ∈ ℤᴺ (ℤ-module), différence de deux positions | `wsl::mddiff<N>` |
| ensemble d'indices | I = [0, n₁) × ⋯ × [0, n_N), décrit par ses étendues | `wsl::extents`, `wsl::dextents` (I2) |
| axe | un nom pour toutes les positions d'une dimension | `wsl::mdaxis` |
| vue, conteneur | un tableau de composantes indexé par I | `wsl::mdspan`, `wsl::mdarray` |

Hors du périmètre :

- **multi-indices de l'analyse** (∂^α, x^α, Taylor, |α|, α!, binomiaux) : objet différent,
  le monoïde (ℕᴺ, +), dans une bibliothèque à part, `wsl::mdindex` (pan 3 de
  docs/algebra.md), avec le devenir de `graded` ;
- **structures algébriques** (`wsl::algebraic`, docs/algebra.md) : voir V6 ;
- **variance** des emplacements d'un tenseur (bibliothèque des tenseurs : un tenseur de
  type (1, 1) est un élément de V ⊗ V* ; les axes restent neutres).

## 2. Décisions validées

| # | Décision |
|---|---|
| ✅ V1 | **Modèle affine**, celui de `std::chrono` : `mdsize<N>` est une **position** dans un tableau (un point, comme `time_point`), `mddiff<N>` un **déplacement** (un vecteur de ℤᴺ, comme `duration`). Position − position = déplacement ; position ± déplacement = position. On n'ajoute, ne multiplie ni ne divise des positions, comme on ne le fait pas de dates. `mdsize` ne sert **qu'aux positions** : les multi-indices de l'analyse sont un autre type (`wsl::mdindex`, plus tard). |
| ✅ V2 | **Un axe désigne toutes les positions d'une dimension** (et non une variable qui les parcourt). Dans des crochets, il étiquette une dimension ; comme valeur, il vaut ses positions (`1.0 / (i + j + 1)`, L12). Aucune arithmétique ne produit un axe : `u[i + 1]` est refusé (`i + 1` est une valeur) ; décalages et découpages se font sur les tableaux (vues, §6). Précisé le 2026-11-25. |
| ✅ V3 | **Sommes : le membre de gauche décide** (convention du mode explicite de `numpy.einsum`, choisie le 2026-11-25 après comparaison avec Einstein, voir §5bis). Un axe présent à gauche n'est jamais sommé : c'est un axe de la boucle, élément par élément à travers toutes les opérations, quel que soit le nombre de ses occurrences. Un axe absent à gauche est sommé sur la plus petite sous-expression qui contient toutes ses occurrences. Un scalaire a un membre de gauche sans axe (`double s = …`, `if (…)`, `f(…)`). Raison décisive : les calculs par point d'un champ (`n2[p] = g[p, i, j] * v[p, i] * v[p, j]`), impossibles avec Einstein. |
| ✅ V4 | **Diffusion** : un axe de gauche absent à droite se diffuse (`y[i] = 1`, `c[i, j] = a[i]`, `m[i, 0] = 0`, `c[i, j] = a[i] + b[j]` ; un tableau n'a pas de base) ; un axe répété à gauche désigne la diagonale (`A[i, i] = 1`). Refusée dans la future bibliothèque des tenseurs, qui aura sa propre règle (variances). Garde-fous : (1) un axe absent à gauche apparaît au moins deux fois à droite (chaque position d'un crochet compte, `position` et `kronecker` aussi), sinon `sum(j, …)` : rattrape l'oubli d'axe (`c[i] = a[i, j]`) et l'ambiguïté de `exp(a[i, j])` ; (2) un axe sommé n'a pas à la fois au moins deux occurrences dans un argument de fonction et d'autres ailleurs (`sqrt(x[j] * x[j]) + y[j] * y[j]` refusé ; `exp(a[i, j]) * b[j]` et `huber(p[i], t[i])` acceptés : Σ) ; (3) un axe lié par `sum(j, …)` n'apparaît ni à gauche ni hors de cette somme. |
| ✅ V5 | **Nommage par préfixe**, comme la norme : `wsl::mdsize`, `wsl::mddiff`, `wsl::mdaxis`, `wsl::mdspan`, `wsl::mdarray` ; sans préfixe ce que la norme ne préfixe pas (`wsl::extents`, `wsl::layout_right`…). **Types préfixés, fonctions surchargées** : `wsl::sum`, `wsl::kronecker`, `wsl::position`… Le nom `mdindex` est réservé aux multi-indices de l'analyse. **Espaces de noms** (2026-11-25) : tout dans `wsl` ; toute fonction publique est contrainte à nos types (pas de modèle non contraint dans `wsl`), et se trouve par ADL (`wrap(p, d, E)`, `sum(j, …)`) ; les fonctions sur des entiers ordinaires (`div_euclid`, `rem_euclid`) sont de l'arithmétique générale, dans `wsl` ; internes dans `wsl::detail::mdarray` ; axes prédéfinis dans `wsl::axes`. |
| ✅ V6 | **Indépendance de `wsl::algebraic`** : mdarray ne dépend d'aucune bibliothèque mathématique ; ses exigences sur les éléments sont syntaxiques, comme dans la norme (`+`, `*`…). La déclaration des structures (`mddiff` est un ℤ-module, `mdsize` un espace affine) viendra plus tard dans un en-tête de liaison, cible qui dépend des deux ; ni l'une ni l'autre ne dépend de l'autre. |
| ✅ V7 | **`array2d` devient un alias** de tableau de rang 2. |
| ✅ V8 | **Deux mondes**. Tableaux entiers : `a + b`, `a * b`, `exp(a)`, `hypot(x, y)`, `2 * a + 1` terme à terme, opérandes de même forme (erreur de compilation entre étendues fixes, précondition sinon), scalaires diffusés, aucune autre diffusion ; paresseux et temporaires seulement (L2) ; indexables (`hypot(x, y)[i]` est un facteur ordinaire) et constructibles (`mdarray c = a * b + 1;`). Axes : V3, V4. |
| ✅ V9 | **Affectations** : `=` parcourt les axes de gauche ; `*=`, `/=`, `+=`, `−=` sont de vrais raccourcis de la forme développée (`y[i, j] *= w[i]` vaut `y[i, j] = y[i, j] * w[i]` : i est à gauche, jamais sommé). **Scalaires** : une expression dont les sommes implicites dépendent du membre de gauche se range dans un type explicite (`double s = x[i] * x[i];`) ou passe par `wsl::value(…)` ; raison : l'expression est construite avant que C++ ne voie la cible, et seule la conversion en nombre dit « membre de gauche sans axe » (`auto` garderait l'expression, refusée par L2). Les réductions explicites entièrement liées (`sum(i, …)`, `count(i, …)`, `all(i, …)`, `any(i, …)`) et les réductions de tableaux entiers sont calculées immédiatement (`auto` permis), sans détection de la réutilisation de leur axe dans l'instruction (`y[i] = count(i, …) * x[i]` passe : déconseillé, documenté). Précisé le 2026-11-25. |

## 3. Principes

| # | Principe | État |
|---|---|---|
| P1 | Des opérateurs seulement pour les opérations du modèle (V1) : `p − q`, `p ± d` ; sur les déplacements `d ± e`, `−d`, et `*`, `/`, `%` par un entier, avec la sémantique des entiers de C++ (troncature vers zéro, comme `chrono::duration`) ; `==` et `<=>`. Le reste prend une fonction nommée | ✅ |
| P2 | Ce qui est absurde est refusé à la compilation (rangs différents, déplacement pris pour une position, garde-fous de V4, étendues fixes incompatibles). Le rang est toujours connu à la compilation, comme dans `std::mdspan`. Un refus est l'absence de l'opération, avec un `static_assert` explicite là où le message du compilateur serait obscur | ✅ |
| P3 | Ce qui est valide sous condition est une précondition : aucune vérification en Release, une assertion (arrêt avec message) quand `NDEBUG` n'est pas défini, forçable par `WSL_ASSERTIONS` (`0` ou `1`) ; toujours vérifiée en évaluation constante (erreur de compilation). Seuls `at()` (`std::out_of_range`, décision 3c de docs/interfaces.md) et `contains()` vérifient toujours. Complète 3c : celle-ci couvre les erreurs toujours détectées (construction, `at()`), P3 les conditions des opérations rapides | ✅ |
| P4 | Dans mdarray, la norme d'abord : ses noms, ses ordres, ses types, ses conventions (indices à partir de 0, intervalles [a, b), `rank()` comme `mdspan`, `layout_right` par défaut) ; puis C++26 et les propositions autour de `mdspan` (`submdspan`, P1673, P1684) ; puis une pratique établie, citée dans la décision qui l'emprunte. Les conventions mathématiques appartiennent aux bibliothèques mathématiques (docs/algebra.md, A1) | ✅ |
| P5 | Un axe n'a pas de valeur, seulement une identité (son type) : le compilateur connaît toute la structure d'une expression | ✅ (V2) |
| P6 | Le domaine d'un calcul est l'axe entier ; pour une partie, on découpe les tableaux | ✅ (V2) |

## 4. Positions et déplacements

Toutes les décisions sont validées (2026-11-25). Écrits au lot 1 (petits, autonomes) ; leur maintien est réexaminé avant la première version
publiée qui les contient (les retirer ensuite serait une rupture).

Exemples : `wsl::mdsize<2> p{2, 1}, q{1, 3};`, `wsl::mddiff<2> d{0, 1}, e{-3, 0};`, un
ensemble d'indices E de forme (4, 5).

**Types et construction**

| # | Cas | Décision | État |
|---|---|---|---|
| D1 | type des composantes | `mdsize<N, IndexType = std::size_t>` (tout entier, comme `std::extents` ; signé : composantes ≥ 0 en précondition), `mddiff<N, DiffType = std::ptrdiff_t>` (signé) ; `p - q` donne `mddiff<N, make_signed_t<IndexType>>`. Rang d'abord : il est toujours donné, le type rarement. Alias `mdsize_for<E>` pour des étendues E (I2) | ✅ |
| D2 | nombre de composantes | exactement N (`mdsize<2>{2}` refusé) ; `mdsize<N>{}` est l'origine, `mddiff<N>{}` le déplacement nul ; pas d'agrégat : constructeurs, comme `std::extents` | ✅ |
| D3 | `mdsize<2>{-1, 0}` | constructeurs acceptant tout type entier, comme `std::extents` ; précondition : composante d'un `mdsize` positive et représentable dans `IndexType`, d'un `mddiff` représentable dans `DiffType` (P3) | ✅ |
| D4 | construction, décomposition, accès | déduction (`mdsize{2, 1}` → `mdsize<2>` en `size_t`, `mddiff{0, -1}` → `mddiff<2>` en `ptrdiff_t`), rang 0, depuis `std::array<I, N>` et `std::span<I, N>` (implicites), depuis `std::span<I>` (explicite, précondition `size() == N`) ; décomposition par le protocole des tuples ; composantes en lecture seule (`p[k]`), `rank()` statique ; `to_array()` vers `std::array<IndexType, N>` | ✅ |
| D5 | rangs différents | refusés à la compilation dans toutes les opérations et conversions (P2) ; étendue dynamique de taille différente : précondition (D4) ; assemblage par `concat` (D24) | ✅ |
| D6 | position ↔ déplacement | aucune conversion : `p - mdsize<2>{}` (déplacement depuis l'origine, comme `time_since_epoch()`), `mdsize<2>{} + d` (précondition : composantes ≥ 0) | ✅ (V1) |

**Opérations**

| # | Cas | Résultat | Décision | État |
|---|---|---|---|---|
| D7 | `p - q` | `mddiff` | ✓ toujours défini | ✅ (P1) |
| D8 | `p ± d`, `d + p`, `+=`, `-=` | `mdsize` | ✓ précondition : résultat ≥ 0 | ✅ (P1) |
| D9 | `d ± e`, `-d`, `+=`, `-=` | `mddiff` | ✓ | ✅ (P1) |
| D10 | `k * d`, `d * k`, `d / k`, `d % k`, `*=`, `/=`, `%=` | `mddiff` | ✓ troncature vers zéro ; k ≠ 0 : précondition | ✅ (P1) |
| D11 | `p + q`, `-p`, `k * p`, `p / k`, `p % k` | — | ✗ | ✅ (V1) |
| D12 | `p + 1`, `d + 1`, `d * e` | — | refusés (P1) ; `mddiff<N>::filled(k)`, `mdsize<N>::filled(k)` : toutes les composantes à k (`p + mddiff<2>::filled(1)`) ; nom du vocabulaire de la norme (`fill`) | ✅ |
| D13 | `++p`, `--p`, `++d` | — | refusés : la position suivante dépend d'un ensemble d'indices ; parcours (I3), où une éventuelle `next(p, E)` sera décidée | ✅ |
| D14 | dépassement du type des composantes (`p + d`, `p - q`, `d + e`, `k * d`) | — | précondition (P3) ; variantes saturées éventuelles en fonctions nommées (I6) | ✅ |
| D15 | division euclidienne | `mddiff` | `div_euclid(d, k)`, `rem_euclid(d, k)` (noms de Rust), composante par composante, reste dans {0, …, k − 1}, précondition k > 0 ; aussi pour les entiers (`wrap`) | ✅ |

**Comparaisons**

| # | Cas | Décision | État |
|---|---|---|---|
| D16 | `==`, `!=` | ✓ entre positions, entre déplacements ; ✗ entre une position et un déplacement | ✅ (V1) |
| D17 | `<`, `<=`, `<=>` | ordre lexicographique total (`std::strong_ordering`), comme `std::array` : ordre du parcours en ligne majeure, algorithmes et conteneurs ordonnés standard ; appartenance au tableau : `E.contains(p)`, pas `p < shape` | ✅ |
| D18 | q ≤ p composante par composante | `all_less(q, p)`, `all_less_equal(q, p)` (ordre partiel), `mdsize` et `mddiff` de même rang ; boîte [lo, hi) : I2 | ✅ |

**Ensembles d'indices et bords**

| # | Cas | Décision | État |
|---|---|---|---|
| D19 | `E.contains(p)` | `E.contains(p)`, `E.contains(i, j, …)` (un entier négatif est hors de E), `E.contains(p, d)` (p + d ∈ E sans former p + d) ; toujours vérifié ; aussi `m.contains(…)` sur les tableaux et les vues | ✅ |
| D20 | `m[p]` hors de E | précondition (P3) ; `m.at(p)`, `m.at(i, j)` : toujours vérifié, `std::out_of_range` (3c), message en français ; mêmes formes de positions que `m[…]`, sans axes | ✅ |
| D21 | `m[d]` | ✗ | ✅ (V1) |
| D22 | `wrap(p, d, E)`, `wrap(p, E)` | ✓ position ramenée dans E par `rem_euclid` composante par composante (bord périodique, tore) ; précondition : étendues non nulles ; nom du mode de bord de NumPy ; périodicité sur certaines dimensions : lot 2, avec les vues | ✅ |
| D23 | `clamp(p, d, E)`, `clamp(p, E)` | ✓ position ramenée au bord le plus proche, composantes bornées à [0, nₖ − 1] (bord recopié) ; précondition : étendues non nulles ; pas de conflit avec `std::clamp` ; autres modes (miroir…) : lot 2, avec les vues | ✅ |

**Divers**

| # | Cas | Décision | État |
|---|---|---|---|
| D24 | `concat(p, q, …)` | ✓ `mdsize<N + M + …>`, de même pour `mddiff` ; position et déplacement mêlés : refusés ; inverse : `p.first<K>()`, `p.last<K>()`, comme `std::span` | ✅ |
| D25 | `mddiff<N>::unit(k)` | ✓ vecteur de base eₖ ; précondition k < N ; pas de `mdsize::unit` ; voisinages (von Neumann, Moore) : lot 2 | ✅ |
| D26 | normes | ✓ `norm1(d)`, `norm_inf(d)`, `norm2_squared(d)`, de type `DiffType` (dépassement : précondition) ; pas de norme d'une position (`norm1(p - q)`) | ✅ |
| D27 | `std::hash`, `std::format`, `<<` | ✓ affichage `(2, 1)`, options appliquées à chaque composante, `<<` par `std::format`, pour `mdsize` et `mddiff` ; pas de lecture pour l'instant | ✅ |
| D28 | types de composantes différents | `p ± d` garde le type de p ; `p - q`, `d ± e` : type commun (signé pour `mddiff`) ; comparaisons sur les valeurs mathématiques ; conversions implicites sans perte, explicites sinon (précondition), comme `std::extents` | ✅ |

## 5. Axes

| # | Écriture | Décision |
|---|---|---|
| L1 | `auto i = wsl::axis();`, `auto mu = wsl::axis<"mu">();`, `auto [mu, nu] = wsl::axis<"mu", "nu">();`, `auto [i, j] = wsl::axis<2>();` | ✅ Un axe est un objet de type `wsl::mdaxis<Id>`, non copiable, dont l'identité est le type ; chaque appel de `wsl::axis` crée des axes distincts (lambda en argument par défaut d'un modèle de fonction). Refusés : étiquettes en double, groupe vide, mélange nombre et étiquettes, déclaration groupée (`auto i = axis(), j = axis();`, règle de `auto`), copie (`auto j = i;`). `wsl::axes` fournit des axes prédéfinis étiquetés (`i, j, k…`). Au niveau d'un espace de noms dans un en-tête : étiquette de type explicite (`inline constexpr wsl::mdaxis<tag> i;`), le lambda n'y étant fiable sur aucun compilateur. Les opérations prennent les axes par référence. Avec Clang 18, une décomposition non `static` ne sert pas directement dans un `static_assert` (P2280). Mesures (2026-11-24, Apple Clang, GCC 14, Clang 18) : `source_location` en argument de modèle donne la ligne du modèle ; `mdaxis i;` par guide de déduction échoue dans les modèles avec Clang ; `mdaxis<> a, b;` fusionne les axes ; réflexion absente (GCC 16 seulement, sans auto-nommage) |
| L2 | `m[1, 2]` contre `m[i, j]` | ✅ entiers : l'élément (référence, comme `mdspan`) ; axes : une expression paresseuse, **utilisable seulement comme temporaire** (opérations sur `Expr&&`) : calculée à l'affectation (`c[i, k] = …`) en une boucle fusionnée ; une expression nommée (`auto e = …; y[i] = e;`) est refusée, avec un message clair ; `std::move` reste une échappatoire explicite ; on nomme une sous-expression par un lambda ou une fonction. Un résultat scalaire se range dans un type explicite (V9). Raisons : pas de référence pendante (les temporaires vivent jusqu'à la fin de l'instruction), pas de surprise si un tableau change entre création et calcul, pas de calcul caché répété. Vérifié sur Apple Clang, GCC 14, Clang 18 (2026-11-24) |
| L3 | `m[1, j]`, `m[i, k]` (k entier, constant ou variable) | ✅ une ligne, une colonne, et en rang N toute section par des positions fixées ; une expression au sens de L2 (temporaire seulement), en lecture comme en affectation (`m[1, j] = v[j]`). Entiers : ceux qu'accepte `mdspan` ; hors du tableau : précondition (P3). Une section qu'on garde est une vue (lot 2) |
| L4 | `m[p]` | ✅ p : `mdsize`, `std::array` ou `std::span` de rang `m.rank()` (comme `mdspan` : compatibilité, positions venues d'ailleurs) : l'élément, comme `m[1, 2]` ; nécessaire au code générique en rang. Position partielle mêlée d'axes (`m[q, j]`) et rang différent : refusés (ajout possible plus tard) |
| L5 | longueur d'un axe | ✅ celle de chaque dimension qu'il étiquette dans l'instruction, membre de gauche compris ; toutes coïncident (axes sommés compris : `a[i, j] * b[j, k]` exige `a.extent(1) == b.extent(0)`) : erreur de compilation entre étendues fixes, précondition (vérifiée une fois par instruction) sinon ; une étendue fixe donne la longueur à une dynamique. Un axe qui n'étiquette aucune dimension d'un tableau est refusé (`double s = i * i;` : longueur inconnue) |
| L6 | `u[i + 1]`, `u[2 * i]`, `b[k - j]`… | ✗ ✅ aucune arithmétique sur un axe (V2) |
| L7 | axe absent à gauche | ✅ sommé sur la plus petite sous-expression qui contient ses occurrences (V3) |
| L8 | `wsl::sum(j, a[i, j])` | ✅ somme explicite, nécessaire quand l'axe n'apparaît qu'une fois (`r[i] = sum(j, a[i, j])`, `exp(x[i]) / sum(j, exp(x[j]))`) ; même forme pour les autres réductions (noms : I6) |
| L9 | diagonale, produit de Hadamard | ✅ `d[i] = A[i, i]`, `A[i, i] = 1` (V3, V4) ; Hadamard : `a * b` (V8) |
| L10 | `c[i, j] = a[i] + b[j]` | ✅ diffusion (V4) |
| L11 | axe multiple (`r[I] = a[I, j]`) | ✅ reporté jusqu'à un besoin réel (précédent : l'ellipse `...` d'`einsum`) ; le calcul élément par élément de tout rang est couvert par V8 ; ajout possible sans rupture |
| L12 | `h[i, j] = 1.0 / (i + j + 1)`, `x[i] = x0 + h * i`, `double m = i * p[i]` | ✅ un axe comme valeur : ses positions 0, …, n − 1, de type signé (`make_signed_t<index_type>`, pour que `i - j` soit juste) ; `wsl::position<T>(i)` pour un autre type ; une occurrence au sens de V4 ; longueur donnée par les tableaux (L5) ; `-i` permis, `+i` refusé. Comparaisons d'axes (masques, `where(i >= j, a[i, j], 0)`) : I10 |
| L13 | `wsl::kronecker(i, j)` | ✅ δᵢⱼ, valeurs 0 ou 1 de type `int` (`kronecker<T>(i, j)` pour un autre type) ; une occurrence de i et de j au sens de V4 ; équivalent lisible de `(i == j)` ; nom sans équivoque (`delta` : Dirac, différences finies) |
| L14 | `2 * m[i, j]`, `m[i, j] + 1` | ✅ un scalaire se diffuse |

## 5bis. Convention des sommes : exemples et comparaison

| Calcul | Écriture |
|---|---|
| produit matriciel, application affine | `c[i, k] = a[i, j] * b[j, k];` `y[i] = A[i, j] * x[j] + b[i];` |
| norme, trace, forme bilinéaire | `double n = sqrt(x[i] * x[i]);` `double t = A[i, i];` `double s = x[p] * A[p, q] * y[q];` |
| diagonale, lignes multipliées, normes de lignes | `d[i] = A[i, i];` `y[i, j] = d[i] * A[i, j];` `n[i] = sqrt(A[i, j] * A[i, j]);` |
| champ, en chaque point p | `n2[p] = g[p, i, j] * v[p, i] * v[p, j];` `a[p, i] = u[p, j] * du[p, i, j];` |
| matrice de masse, Hooke | `M[i, j] = w[q] * phi[i, q] * phi[j, q];` `s[i, j] = lambda * e[k, k] * kronecker(i, j) + 2 * mu * e[i, j];` |
| softmax, distances, perte | `y[i] = exp(x[i]) / sum(j, exp(x[j]));` `D[i, j] = hypot(x[i] - x[j], y[i] - y[j]);` `double L = huber(p[i], t[i]);` |
| élément par élément | `y = x * x;` `r = hypot(x, y);` `y = x / sqrt(x * x);` |

Conventions écartées (2026-11-24 et 25) : Einstein (un axe répété est sommé ; étendue aux
quotients et fonctions, avec axes ponctuels, `keep`, marques `+i`, axes libres et muets
déclarés) : sens local et fidélité aux livres, mais calculs par point d'un champ impossibles
avec des axes ; membre de gauche décide sans garde-fou : oubli d'axe sommé en silence,
ambiguïté de `exp(a[i, j])` ; sommes toujours explicites (Halide) : lourd. Coût accepté : le
sens d'une sous-expression dépend du membre de gauche de son instruction ; `+i`, `-i` restent
réservés aux tenseurs (variances).

## 6. Stockage, vues, calcul (à concevoir : stockage et calcul au lot 1, vues au lot 2)

- **Infrastructure** (✅ A5, validée le 2026-11-25) : `wsl::extents`, `wsl::dextents`,
  `layout_right`, `layout_left`, `layout_stride`, `wsl::mdspan`, aux noms et à la sémantique
  de `<mdspan>`, dont ils ne s'écartent que par des **ajouts** (axes, `mdsize`, `contains`,
  `at`, `concat`… ; ajouts de C++26 repris sous leur nom) : tout code valide avec
  `std::mdspan` garde le même sens. Conversions avec `std::mdspan` là où il existe
  (`__cpp_lib_mdspan`). Raisons : `std::mdspan` manque à GCC 14 et à Clang 18 avec
  libstdc++ ; surtout, `m[i, j]` avec des axes exige notre propre vue (on n'ajoute pas de
  surcharge à un type de la norme) ; pas de dépendance externe (Kokkos).
- **Ensembles d'indices** (✅ I2, 2026-11-25) : `wsl::extents`, `wsl::dextents` de la norme ;
  un `extents` est l'ensemble d'indices E (ajouts : `contains`, parcours, `wrap`, `clamp`) ;
  pas de conversion vers `mdsize` (une forme est une taille, pas une position) ;
  `wsl::dims<N>` (C++26, P2389) ; étendues nulles et rang 0 permis ; boîtes [lo, hi) : I3 et
  lot 2.
- **Conteneur** (✅ I13, 2026-11-25) : `wsl::mdarray<T, Extents, Layout = layout_right,
  Container>` (P1684), possédant, sémantique de valeur (copie profonde). Conteneur par défaut :
  `std::array<T, n₁ ⋯ n_N>` si toutes les étendues sont fixes (aucune allocation, `constexpr`
  complet : tenseurs 3 × 3, 4 × 4), `std::vector<T>` sinon (écart assumé avec P1684).
  Construction : `(extents)`, `(n₁, n₂…)` (valeurs initialisées), `(extents, valeur)`, depuis
  un `mdspan` (copie), depuis une expression de tableaux entiers (V8) ; littérale pour des
  étendues fixes (`{{1, 2}, {3, 4}}`, ajout). Affectation d'une expression : formes égales
  (précondition), pas de redimensionnement implicite. `m.fill(v)` ; `to_mdspan()` et
  conversion implicite en `wsl::mdspan` (P1684), en `std::mdspan` là où il existe ; accès de
  `mdspan` plus `data()`, `container()`, `extract_container()`, `size()`, `extents()`,
  `mapping()` ; parcours à plat `begin()`/`end()` sur le stockage (ajout) ; `constexpr`
  partout. **Comparaisons des tableaux entiers** (`==`, `!=`, `<`, `<=`, `>`, `>=`) :
  élément par élément, des masques (cohérent avec V8 ; précédent : `std::valarray`) ;
  `equal(a, b)` rend un `bool` (nom de `std::ranges::equal`) ; `all(m)`, `any(m)` réduisent
  un masque ; `if (masque)` refusé avec un message. Coût accepté : le tableau n'est pas
  `std::equality_comparable`. `array2d` (V7) : décidé à la fin du lot 2, quand les vues de
  lignes et de colonnes existeront ; inchangé d'ici là.
- **Algèbre des expressions** (✅ I5, 2026-11-25) : toutes les vérifications à la compilation
  (garde-fous de V4, longueurs fixes, mélanges interdits), messages en français nommant la
  règle ; seules les longueurs dynamiques sont vérifiées à l'exécution, une fois par
  instruction. Une boucle fusionnée par instruction, sans tableau intermédiaire caché ; ordre
  des boucles de gauche choisi par la bibliothèque (ordre du stockage du tableau de gauche),
  sans effet sur le résultat ; axes sommés parcourus dans l'ordre de l'écriture (première
  apparition), indépendamment des dispositions ; sommes imbriquées calculées à leur niveau ;
  représentation interne libre.
- **Calcul flottant** (✅ I9 et fusion, 2026-11-25) : ordre des sommes et réductions fixé par
  la bibliothèque, par paires avec des feuilles en K accumulateurs (K fixé, comme NumPy) :
  reproductible, vectorisable, erreur en O(log n ε). Deux modes : **rapide** (fusion en `fma`
  permise, fonctions de `<cmath>` ; mêmes résultats pour un même programme, une même machine,
  les mêmes options) et **reproductible** (sans fusion : barrières placées par le moteur sur
  des vecteurs, évaluation par paquets ; fonctions de `wsl::special` ; une fonction non
  couverte est une erreur de compilation ; mêmes bits partout, `constexpr` = exécution).
  Choix par instruction (`wsl::fast(…)`, `wsl::reproducible(…)`) ; mode par défaut fixé par
  une option CMake de la cible (rapide par défaut), cohérente dans tout le projet (une macro
  libre violerait la règle de définition unique). Mesure (Apple Clang, AVX2) : la barrière
  scalaire de `wsl::special` empêche la vectorisation, la barrière sur un vecteur la garde.
- **Comparaisons avec des axes** (✅ I10, 2026-11-25) : élément par élément selon la convention
  C ; une somme de booléens est refusée (message vers `any`, `all`, `count`), y compris
  `if (a[i] < b[i])` et `double s = a[i] < b[i];` ; une comparaison dans un produit numérique
  est un crochet d'Iverson (`(a[i] < b[i]) * x[i]`) ; réductions `all`, `any`, `count` sur
  les masques entiers et par les axes (`count(i, …)`) ; `count` rend
  `make_signed_t<index_type>` (comme `std::ranges::count`) ; `&&`, `||`, `!` élément par
  élément sans court-circuit (comme `std::valarray`) ; `where(c, x, y)`.
- **Vues de découpage et de décalage**, noms provisoires : `u.slice(a, b)` (positions
  [a, b), comme `submdspan` de C++26), `roll(u, k)` (décalage circulaire, `numpy.roll`),
  `shift(u, k, wsl::clamp)` ou `shift(u, k, wsl::fill(0))` ; décalages de rang N par
  `mddiff`. Stencil : `v.slice(1, n - 1)[i] = u.slice(2, n)[i] - 2 * u.slice(1, n - 1)[i]
  + u.slice(0, n - 2)[i]` ; périodique : `v[i] = roll(u, 1)[i] - 2 * u[i] + roll(u, -1)[i]`.
- **Calcul des expressions** : évaluation paresseuse, aliasing (`m[i, j] = m[j, i]` sur
  place).

## 7. Décisions restantes

| # | Sujet |
|---|---|
| I3 | Parcours d'un ensemble d'indices (ordre lexicographique, ordre de stockage) |
| I6 | Noms des réductions (produit, minimum, maximum, tous, au moins un ; conflits avec `std::max`) ; `filled`, comparaison composante par composante, division euclidienne et normes sont réglés (D12, D15, D18, D26) |
| I8 | Exécution parallèle et GPU (appel depuis le code des noyaux, absence d'exceptions, moteur d'exécution : dépendance externe à décider) ; d'ici là, types `constexpr`, trivialement copiables, sans exception hors `at()` |
| I11 | Fonctions de l'utilisateur sur des expressions et des tableaux entiers (`wsl::map(f, x, y)`) ; surcharges de `<cmath>` pour les expressions |
| I12 | Aliasing : `y[i] = A[i, j] * y[j]`, `m[i, j] = m[j, i]` sur place |
