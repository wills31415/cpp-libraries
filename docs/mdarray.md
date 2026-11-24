# wsl::mdarray : tableaux multidimensionnels

Chantier ouvert le 2026-11-24, préalable aux objets mathématiques (docs/algebra.md), qui
s'en serviront pour stocker des composantes. La bibliothèque est une infrastructure
générale, indépendante de ces objets : elle sert aussi bien des images, des grilles ou des
plateaux de jeu. `array2d` en deviendra un alias de rang 2.

Livraison par lots (en-têtes séparés, une seule cible) :

1. positions, déplacements et axes (§4, §5) ;
2. ensembles d'indices, dispositions et vue (§6) ;
3. stockage et calcul des expressions (§6).

État des décisions : ✅ validée ; ◻︎ proposée, en attente de validation.

## 1. Objet

| Objet | En mathématiques | Type |
|---|---|---|
| indice entier | une position i ∈ {0, …, n − 1} d'une dimension | `std::size_t` (le type d'indice des étendues) |
| position | un point p ∈ ℕᴺ de l'espace affine ℤᴺ | `wsl::mdsize<N>` |
| déplacement | un vecteur d ∈ ℤᴺ (ℤ-module), différence de deux positions | `wsl::mddiff<N>` |
| ensemble d'indices | I = [0, n₁) × ⋯ × [0, n_N), décrit par ses étendues | `wsl::extents`, `wsl::dextents` (I2) |
| axe | un nom pour toutes les positions d'une dimension | `wsl::mdaxis` |
| vue, conteneur | un tableau de composantes indexé par I | `wsl::mdspan`, `wsl::mdarray` (lots 2 et 3) |

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
| ✅ V2 | **Un axe désigne toutes les positions d'une dimension** (et non une variable qui les parcourt) : aucune arithmétique sur un axe ; décalages et découpages se font sur les tableaux (vues, §6). |
| ✅ V3 | **Sommes selon la convention d'Einstein** : un axe répété deux fois dans un produit est sommé au niveau de ce produit (`sqrt(v[i] * v[i])` est la norme) ; une troisième occurrence est refusée. |
| ✅ V4 | **Diffusion acceptée** (`c[i, j] = a[i] + b[j]` : un tableau n'a pas de base), refusée dans la future bibliothèque des tenseurs (Einstein strict : a_i + b_j n'est pas une opération tensorielle). Deux garde-fous : (1) les axes libres du membre de droite, réunis sur ses termes, sont exactement ceux du membre de gauche ; (2) tout axe du membre de gauche apparaît à droite. |
| ✅ V5 | **Nommage par préfixe**, comme la norme : `wsl::mdsize`, `wsl::mddiff`, `wsl::mdaxis`, `wsl::mdspan`, `wsl::mdarray` ; sans préfixe ce que la norme ne préfixe pas (`wsl::extents`, `wsl::layout_right`…). **Types préfixés, fonctions surchargées** : `wsl::sum`, `wsl::kronecker`, `wsl::position`… Le nom `mdindex` est réservé aux multi-indices de l'analyse. |
| ✅ V6 | **Indépendance de `wsl::algebraic`** : mdarray ne dépend d'aucune bibliothèque mathématique ; ses exigences sur les éléments sont syntaxiques, comme dans la norme (`+`, `*`…). La déclaration des structures (`mddiff` est un ℤ-module, `mdsize` un espace affine) viendra plus tard dans un en-tête de liaison, cible qui dépend des deux ; ni l'une ni l'autre ne dépend de l'autre. |
| ✅ V7 | **`array2d` devient un alias** de tableau de rang 2. |

## 3. Principes

| # | Principe | État |
|---|---|---|
| P1 | Des opérateurs seulement pour les opérations du modèle (V1) : `p − q`, `p ± d` ; sur les déplacements `d ± e`, `−d`, et `*`, `/`, `%` par un entier, avec la sémantique des entiers de C++ (troncature vers zéro, comme `chrono::duration`) ; `==` et `<=>`. Le reste prend une fonction nommée | ✅ |
| P2 | Ce qui est absurde est refusé à la compilation (rangs différents, déplacement pris pour une position…) | ◻︎ |
| P3 | Ce qui est valide sous condition est une précondition, comme dans la norme : aucun coût en Release, une assertion en débogage, une erreur de compilation en évaluation constante ; seuls `at()` et `contains()` vérifient toujours | ◻︎ |
| P4 | La norme d'abord : ses noms, ses ordres, ses types | ◻︎ |
| P5 | Un axe n'a pas de valeur, seulement une identité (son type) : le compilateur connaît toute la structure d'une expression | ✅ (V2) |
| P6 | Le domaine d'un calcul est l'axe entier ; pour une partie, on découpe les tableaux | ✅ (V2) |

## 4. Positions et déplacements

Exemples : `wsl::mdsize<2> p{2, 1}, q{1, 3};`, `wsl::mddiff<2> d{0, 1}, e{-3, 0};`, un
ensemble d'indices E de forme (4, 5).

**Types et construction**

| # | Cas | Décision | État |
|---|---|---|---|
| D1 | type des composantes | `std::size_t` par défaut, paramétrable comme `std::extents<IndexType, …>` ; `mddiff` prend le type signé correspondant | ◻︎ |
| D2 | `mdsize<2>{2}` | ✗ exactement N composantes ; `mdsize<2>{}` est l'origine | ◻︎ |
| D3 | `mdsize<2>{-1, 0}` | ✗ pour une constante ; précondition sinon | ◻︎ |
| D4 | `mdsize{2, 1}`, `mdsize<0>{}`, depuis `std::array`, `std::span`, `auto [i, j] = p` | ✓ | ◻︎ |
| D5 | rangs différents | ✗ | ◻︎ |
| D6 | position ↔ déplacement | aucune conversion : `p - mdsize<2>{}` (déplacement depuis l'origine, comme `time_since_epoch()`), `mdsize<2>{} + d` (précondition : composantes ≥ 0) | ✅ (V1) |

**Opérations**

| # | Cas | Résultat | Décision | État |
|---|---|---|---|---|
| D7 | `p - q` | `mddiff` | ✓ toujours défini | ✅ (P1) |
| D8 | `p ± d`, `d + p`, `+=`, `-=` | `mdsize` | ✓ précondition : résultat ≥ 0 | ✅ (P1) |
| D9 | `d ± e`, `-d`, `+=`, `-=` | `mddiff` | ✓ | ✅ (P1) |
| D10 | `k * d`, `d * k`, `d / k`, `d % k`, `*=`, `/=`, `%=` | `mddiff` | ✓ troncature vers zéro ; k ≠ 0 : précondition | ✅ (P1) |
| D11 | `p + q`, `-p`, `k * p`, `p / k`, `p % k` | — | ✗ | ✅ (V1) |
| D12 | `p + 1`, `d + 1`, `d * e` | — | ✗ ; `p + mddiff<2>::filled(1)` (nom en I6) | ◻︎ |
| D13 | `++p` | — | ✗ ; itérateurs du parcours (I3) | ◻︎ |
| D14 | dépassement du type des composantes | — | précondition | ◻︎ |
| D15 | division euclidienne de `d` par k (reste dans {0, …, k − 1}) | `mddiff` | fonction nommée, comme `std::chrono::floor` (nom en I6) | ◻︎ |

**Comparaisons**

| # | Cas | Décision | État |
|---|---|---|---|
| D16 | `==`, `!=` | ✓ entre positions, entre déplacements ; ✗ entre une position et un déplacement | ✅ (V1) |
| D17 | `<`, `<=`, `<=>` | ✓ ordre lexicographique, total (`std::strong_ordering`) : ordre du parcours en ligne majeure, tris, `std::map` | ◻︎ |
| D18 | q ≤ p composante par composante | fonction nommée (nom en I6) | ◻︎ |

**Ensembles d'indices et bords**

| # | Cas | Décision | État |
|---|---|---|---|
| D19 | `E.contains(p)` | ✓ | ◻︎ |
| D20 | `m[p]` hors de E | précondition ; `m.at(p)` lève `std::out_of_range` | ◻︎ |
| D21 | `m[d]` | ✗ | ✅ (V1) |
| D22 | `wrap(p, d, E)` | ✓ position p + d ramenée dans E, modulo euclidien ((−1) mod 4 = 3) : bord périodique, tore ℤᴺ/(n₁ℤ × ⋯ × n_Nℤ) | ◻︎ |
| D23 | `clamp(p, d, E)` | ✓ position p + d ramenée au bord le plus proche : bord recopié | ◻︎ |

**Divers**

| # | Cas | Décision | État |
|---|---|---|---|
| D24 | `concat(p, q)` | ✓ `mdsize<4>` : position dans un produit de tableaux | ◻︎ |
| D25 | `mddiff<N>::unit(k)` | ✓ base canonique e_k de ℤᴺ (voisins d'un stencil) | ◻︎ |
| D26 | `norm1(d)`, `norm_inf(d)` | ✓ fonctions nommées (rayon d'un stencil, voisinages) (I6) | ◻︎ |
| D27 | `std::hash`, `std::format`, `<<` | ✓ (I7) | ◻︎ |
| D28 | types de composantes différents | ✓ type commun, comme l'expression scalaire | ◻︎ |

## 5. Axes (◻︎ sauf mention)

| # | Écriture | Décision |
|---|---|---|
| L1 | `constexpr wsl::mdaxis i, j, k;` | ✓ chacun de son propre type, réutilisable |
| L2 | `m[1, 2]` contre `m[i, j]` | ✓ entiers : un élément ; axes : une expression paresseuse |
| L3 | `m[1, j]`, `m[i, 0]` | ✓ une ligne, une colonne |
| L4 | `m[p]` | ✓ une position fixe |
| L5 | longueur d'un axe | déduite des tableaux ; toutes les dimensions qu'il étiquette ont la même longueur (✗ à la compilation pour des étendues fixes, précondition sinon) |
| L6 | `u[i + 1]`, `u[2 * i]`, `b[k - j]`… | ✗ ✅ aucune arithmétique sur un axe (V2) |
| L7 | axe répété dans un produit | ✅ sommé (V3) |
| L8 | `wsl::sum(j, a[i, j])` | ✓ somme explicite (sommes de lignes) |
| L9 | diagonale, produit de Hadamard | ✓ opérations nommées sur les tableaux entiers (`diagonal(m)`, `a * b`) |
| L10 | `c[i, j] = a[i] + b[j]` | ✅ diffusion, avec ses deux garde-fous (V4) |
| L11 | `wsl::mdaxis<2> I;` puis `c[I] = a[I] + b[I]` | ✓ plusieurs axes à la fois, pour les expressions élément par élément de tout rang (coexistence avec L1 : I5) |
| L12 | `h[i, j] = 1.0 / (wsl::position(i) + wsl::position(j) + 1)` | ✓ `position(i)` : les positions 0, …, n − 1 le long de l'axe i (`arange` de NumPy) ; l'arithmétique porte sur des valeurs, jamais sur l'axe |
| L13 | `wsl::kronecker(i, j)` | ✓ δ_ij sur les axes i et j |
| L14 | `2 * m[i, j]`, `m[i, j] + 1` | ✓ un scalaire se diffuse |

## 6. Stockage, vues, calcul (lots 2 et 3, à concevoir)

- **Infrastructure** (◻︎, proposée le 2026-11-23) : étendues, dispositions et vue propres, aux
  noms et à la sémantique de `<mdspan>`, dont elles ne s'écartent que par des **ajouts** :
  tout code valide avec `std::mdspan` garde le même sens. Raison : `std::mdspan` manque à
  GCC 14 et 15 (GCC 16 ou libc++ requis), et la vue gagnera des extensions (positions
  `mdsize`, axes, indexation partielle). Conversion vers `std::mdspan` là où il existe.
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
| I2 | Ensembles d'indices : `extents` fixes ou dynamiques, alignement sur `std::extents`, lien avec `mdsize` (forme d'un tableau) |
| I3 | Parcours d'un ensemble d'indices (ordre lexicographique, ordre de stockage) |
| I5 | Algèbre des expressions d'axes à la compilation : axes libres, répétés, garde-fous de V4 ; coexistence de `mdaxis` et `mdaxis<N>` |
| I6 | Fonctions nommées : `filled`, comparaison composante par composante, division euclidienne, normes |
| I7 | Affichage et lecture des positions et déplacements |
