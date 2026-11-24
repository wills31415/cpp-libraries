# wsl::mdindex : indices, multi-indices, ensembles d'indices, axes

Chantier ouvert le 2026-11-24, premier des trois chantiers des objets mathématiques :
`wsl::algebraic` et `wsl::mdindex`, puis `wsl::mdarray`, puis docs/algebra.md. Il remplace
`graded` et `graded_io`, dont l'utile est refondu ici avant leur suppression.

État des décisions : ✅ validée ; ◻︎ proposée, en attente de validation.

## 1. Objet

Tout le vocabulaire des indices, sans aucun tableau : `mdindex` **dit** ce que signifie une
expression indicielle, `mdarray` (chantier suivant) la **calcule**.

| Objet | En mathématiques | Type |
|---|---|---|
| indice entier | une position i ∈ {0, …, n − 1} | `std::size_t` (le type d'indice des étendues) |
| multi-indice | α ∈ ℕᴺ : position dans un tableau, forme, ordre de dérivation | `wsl::mdsize<N>` |
| décalage | δ ∈ ℤᴺ, différence de deux multi-indices | `wsl::mddiff<N>` |
| ensemble d'indices | I = [0, n₁) × ⋯ × [0, n_N), décrit par ses étendues | `wsl::extents`, `wsl::dextents` |
| axe | un nom pour toutes les positions d'une dimension | `wsl::mdaxis`, `wsl::mdaxis<N>` |
| opérations d'analyse | \|α\|, α!, binomiaux, ordre composante par composante | fonctions (I6) |

Hors du périmètre : dispositions, vues et conteneurs (mdarray) ; variance des emplacements
d'un tenseur (bibliothèque des tenseurs : un tenseur de type (1, 1) est un élément de
V ⊗ V* ; les axes restent neutres).

## 2. Décisions validées

| # | Décision |
|---|---|
| ✅ V1 | **Structure B** : (ℕᴺ, +) est un monoïde commutatif simplifiable (positions et formes), (ℤᴺ, +) son groupe de Grothendieck (décalages) et un ℤ-module ; ℤᴺ agit sur ℕᴺ (position + décalage). Comme `size_t` et `ptrdiff_t`. |
| ✅ V2 | **Un axe désigne toutes les positions d'une dimension** (et non une variable qui les parcourt) : aucune arithmétique sur un axe ; décalages et découpages se font sur les tableaux (vues de mdarray, §6). |
| ✅ V3 | **Sommes selon la convention d'Einstein** : un axe répété deux fois dans un produit est sommé au niveau de ce produit (`sqrt(v[i] * v[i])` est la norme) ; une troisième occurrence est refusée. |
| ✅ V4 | **Diffusion acceptée** dans mdindex et mdarray (`c[i, j] = a[i] + b[j]` : un tableau n'a pas de base), refusée dans la future bibliothèque des tenseurs (Einstein strict : a_i + b_j n'est pas une opération tensorielle). Deux garde-fous : (1) les axes libres du membre de droite, réunis sur ses termes, sont exactement ceux du membre de gauche ; (2) tout axe du membre de gauche apparaît à droite. |
| ✅ V5 | **Nommage par préfixe**, comme la norme : `wsl::mdsize`, `wsl::mddiff`, `wsl::mdaxis`, puis `wsl::mdspan`, `wsl::mdarray` ; sans préfixe ce que la norme ne préfixe pas (`wsl::extents`, `wsl::layout_right`…). **Types préfixés, fonctions surchargées** : `wsl::sum`, `wsl::kronecker`, `wsl::position`… |
| ✅ V6 | **Structures algébriques par déclaration explicite** : un type n'est un monoïde, un groupe… que s'il le déclare (trait à spécialiser, qui donne aussi l'élément neutre) ; le concept vérifie la syntaxe en plus. Raison : la syntaxe ne dit pas quelle structure un type modélise (`size_t` avec `+` est ℤ/2⁶⁴ℤ pour la machine, ℕ pour l'intention). Dans une **bibliothèque séparée** (◻︎ nom : `wsl::algebraic`), écrite dans ce chantier : magma, demi-groupe, monoïde, groupe et leurs variantes commutatives ; étendue plus tard aux anneaux, corps, modules, espaces vectoriels. |

## 3. Principes

| # | Principe | État |
|---|---|---|
| P1 | Des opérateurs seulement pour les opérations de la structure (V1) ; le reste prend une fonction nommée | ◻︎ |
| P2 | Ce qui est absurde est refusé à la compilation (rangs différents, décalage pris pour une position…) | ◻︎ |
| P3 | Ce qui est valide sous condition est une précondition, comme dans la norme : aucun coût en Release, une assertion en débogage, une erreur de compilation en évaluation constante ; seuls `at()` et `contains()` vérifient toujours | ◻︎ |
| P4 | La norme d'abord : ses noms, ses ordres, ses types | ◻︎ |
| P5 | Un axe n'a pas de valeur, seulement une identité (son type) : le compilateur connaît toute la structure d'une expression | ✅ (V2) |
| P6 | Le domaine d'un calcul est l'axe entier ; pour une partie, on découpe les tableaux | ✅ (V2) |

## 4. Arithmétique des multi-indices et des décalages (◻︎)

Exemples : `wsl::mdsize<2> alpha{2, 1}, beta{1, 3}, gamma{1, 0};`,
`wsl::mddiff<2> d{0, 1}, e{-3, 0};`, un ensemble d'indices E de forme (4, 5).

**Types et construction**

| # | Cas | Décision |
|---|---|---|
| D1 | type des composantes | `std::size_t` par défaut, paramétrable comme `std::extents<IndexType, …>` ; `mddiff` prend le type signé correspondant |
| D2 | `mdsize<2>{2}` | ✗ exactement N composantes ; `mdsize<2>{}` est l'origine |
| D3 | `mdsize<2>{-1, 0}` | ✗ pour une constante ; précondition sinon |
| D4 | `mdsize{2, 1}`, `mdsize<0>{}`, depuis `std::array`, `std::span`, `auto [i, j] = alpha` | ✓ |
| D5 | rangs différents | ✗ |
| D6 | `mddiff d = alpha` | ✓ conversion implicite (ℕᴺ ↪ ℤᴺ) |
| D7 | `mdsize(d)` | ✓ explicite, précondition : composantes ≥ 0 |

**Opérations**

| # | Cas | Résultat | Décision |
|---|---|---|---|
| D8 | `alpha + beta`, `+=` | `mdsize` | ✓ ; dépassement : précondition |
| D9 | `alpha - beta` | `mddiff` | ✓ toujours défini |
| D10 | `alpha -= beta` | `mdsize` | ✓ précondition β ≤ α |
| D11 | `alpha ± d`, `d + alpha`, `+=`, `-=` | `mdsize` | ✓ précondition : résultat ≥ 0 |
| D12 | `d ± e`, `-d`, `-alpha` | `mddiff` | ✓ |
| D13 | `k * alpha`, `alpha * k` | `mdsize` | ✓ précondition k ≥ 0 |
| D14 | `k * d`, `d * k` | `mddiff` | ✓ |
| D15 | `alpha + 1` | — | ✗ ; `alpha + mdsize<2>::filled(1)` (nom à fixer) |
| D16 | `++alpha` | — | ✗ ; itérateurs du parcours (I3) |
| D17 | `alpha * beta`, `alpha / k`, `alpha % k` | — | ✗ ; fonctions nommées plus tard si besoin |
| D18 | `sub_sat(alpha, beta)` | `mdsize` | ✓ max(α − β, 0), au nom de C++26 (`std::sub_sat`) |

**Comparaisons**

| # | Cas | Décision |
|---|---|---|
| D19 | `==`, `!=` | ✓, y compris entre `mdsize` et `mddiff` (valeurs mathématiques) |
| D20 | `<`, `<=`, `<=>` | ✓ ordre lexicographique, total (`std::strong_ordering`) : tris, `std::map`, parcours en ligne majeure |
| D21 | β ≤ α composante par composante | ✓ fonction nommée (nom en I6) |

**Ensembles d'indices et bords**

| # | Cas | Décision |
|---|---|---|
| D22 | `E.contains(alpha)` | ✓ |
| D23 | `m[alpha]` hors de I | précondition ; `m.at(alpha)` lève `std::out_of_range` |
| D24 | `m[d]` | ✗ |
| D25 | `wrap(x, E)` | ✓ modulo euclidien ((−1) mod 4 = 3) : bord périodique, tore ℤᴺ/(n₁ℤ × ⋯ × n_Nℤ) |
| D26 | `clamp(x, E)` | ✓ bord recopié |

**Divers**

| # | Cas | Décision |
|---|---|---|
| D27 | `concat(alpha, beta)` | ✓ `mdsize<4>` (indice d'un produit tensoriel) |
| D28 | `mddiff<N>::unit(k)` | ✓ base canonique e_k de ℤᴺ |
| D29 | \|α\|, α!, binomiaux, `norm1(d)` | ✓ fonctions nommées (I6) |
| D30 | `std::hash`, `std::format`, `<<` | ✓ (I7) |
| D31 | types de composantes différents | ✓ type commun, comme l'expression scalaire |

## 5. Axes (◻︎ sauf mention)

| # | Écriture | Décision |
|---|---|---|
| L1 | `constexpr wsl::mdaxis i, j, k;` | ✓ chacun de son propre type, réutilisable |
| L2 | `m[1, 2]` contre `m[i, j]` | ✓ entiers : un élément ; axes : une expression paresseuse |
| L3 | `m[1, j]`, `m[i, 0]` | ✓ une ligne, une colonne |
| L4 | `m[alpha]` | ✓ une position fixe |
| L5 | longueur d'un axe | déduite des tableaux ; toutes les dimensions qu'il étiquette ont la même longueur (✗ à la compilation pour des étendues fixes, précondition sinon) |
| L6 | `u[i + 1]`, `u[2 * i]`, `b[k - j]`… | ✗ ✅ aucune arithmétique sur un axe (V2) |
| L7 | axe répété dans un produit | ✅ sommé (V3) |
| L8 | `wsl::sum(j, a[i, j])` | ✓ somme explicite (sommes de lignes) |
| L9 | diagonale, produit de Hadamard | ✓ opérations nommées sur les tableaux entiers (`diagonal(m)`, `a * b`), dans mdarray |
| L10 | `c[i, j] = a[i] + b[j]` | ✅ diffusion, avec ses deux garde-fous (V4) |
| L11 | `wsl::mdaxis<2> I;` puis `c[I] = a[I] + b[I]` | ✓ plusieurs axes à la fois, pour les expressions élément par élément de tout rang |
| L12 | `h[i, j] = 1.0 / (wsl::position(i) + wsl::position(j) + 1)` | ✓ `position(i)` : les positions 0, …, n − 1 le long de l'axe i (`arange` de NumPy) ; l'arithmétique porte sur des valeurs, jamais sur l'axe |
| L13 | `wsl::kronecker(i, j)` | ✓ δ_ij sur les axes i et j |
| L14 | `2 * m[i, j]`, `m[i, j] + 1` | ✓ un scalaire se diffuse |

## 6. Transmis à mdarray

- Calcul effectif des expressions, aliasing (`m[i, j] = m[j, i]` sur place).
- Vues de découpage et de décalage, noms provisoires : `u.slice(a, b)` (positions [a, b),
  comme `submdspan` de C++26), `roll(u, k)` (décalage circulaire, `numpy.roll`),
  `shift(u, k, wsl::clamp)` ou `shift(u, k, wsl::fill(0))` ; décalages de rang N par
  `mddiff`. Stencil : `v.slice(1, n - 1)[i] = u.slice(2, n)[i] - 2 * u.slice(1, n - 1)[i]
  + u.slice(0, n - 2)[i]` ; périodique : `v[i] = roll(u, 1)[i] - 2 * u[i] + roll(u, -1)[i]`.

## 7. Décisions restantes

| # | Sujet |
|---|---|
| I2 | Ensembles d'indices : `extents` fixes ou dynamiques, alignement sur `std::extents` |
| I3 | Parcours d'un ensemble d'indices (ordre lexicographique, ordre de stockage) |
| I5 | Algèbre des expressions d'axes à la compilation : axes libres, répétés, garde-fous de V4 |
| I6 | Opérations d'analyse : \|α\|, α!, binomiaux, ordre partiel, x^α ; noms `filled`, comparaison composante par composante |
| I7 | Affichage et lecture (repris de `graded_io`) ; suppression de `graded` et `graded_io` (rupture, au CHANGELOG) |
