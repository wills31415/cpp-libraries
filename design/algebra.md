# Objets mathématiques : algèbre linéaire, tenseurs, analyse

Chantier ouvert le 2026-11-23, en remplacement du chantier `mdarray` de la revue des
interfaces (design/interfaces.md, décision 6). Il s'appuiera sur un chantier préalable,
mené à part : `wsl::mdarray` (design/mdarray.md), le stockage des composantes, indépendant des
objets mathématiques. Ce document fixe la vision, la méthode et la
feuille de route ; les notions y seront ensuite conçues une à une, pan par pan.

## 1. Vision

Écrire et manipuler les objets des mathématiques appliquées comme dans un cours ou un livre :
espaces vectoriels munis de leurs structures, bases, applications linéaires, formes,
tenseurs, fonctions entre espaces, différentielles, dérivées partielles ∂^α, puis espaces
de fonctions, opérateurs et variétés. Les tableaux ne sont que la **représentation en
coordonnées** de ces objets : une matrice est une application linéaire lue dans deux
bases, un tableau multidimensionnel stocke des composantes.

C'est la démarche de wsl::random (des lois qui sont des mesures, des variables aléatoires
qui sont des applications mesurables), appliquée à l'algèbre linéaire et à l'analyse.

Consommateurs visés, hors de ce dépôt : un projet de relativité générale (équations,
visualisations ; le dépôt cpp-tensors en est un premier essai, à refondre sur ces
bibliothèques) et un futur projet de mécanique des milieux continus et des fluides.

## 2. Décisions de méthode

| # | Décision |
|---|---|
| A1 | **Définitions tirées de la littérature publiée**, citée dans chaque section, avec la correspondance entre symboles et noms ; ne rien inventer quand un usage établi existe. Les conventions qui divergent entre communautés (signature (−+++) ou (+−−−), linéarité du produit hermitien, place des indices, convention d'Einstein) sont choisies une fois et écrites ici. |
| A2 | **Une définition = un concept C++ ; les hypothèses d'un théorème = les contraintes du code.** Exemple : `grad(f)` n'existe que si l'espace a un produit scalaire (isomorphisme de Riesz) ; sans lui, seule la différentielle Df(x) ∈ E* est définie. |
| A3 | **Un espace vectoriel est une structure, pas une dimension** : corps des scalaires, dimension, puis structures ajoutées (norme, produit scalaire, forme bilinéaire non dégénérée de signature donnée, orientation, structure affine). |
| A4 | **Tout à la compilation par défaut** (vérifications, absence de coût à l'exécution, vectorisation, `constexpr`) ; dimension dynamique (`wsl::dynamic`) en exception assumée, là où elle n'est connue qu'à l'exécution (discrétisations). |
| A5 | **Stockage des composantes** : `wsl::mdarray` (design/mdarray.md), bibliothèque générale, qui ne dépend d'aucune bibliothèque mathématique ; ses positions (`mdposition`) et déplacements (`mddisplacement`) ne servent qu'à indexer des tableaux. |
| A6 | **Multi-indices de l'analyse** (∂^α, x^α, Taylor, \|α\|, α!, binomiaux) : le monoïde (ℕᴺ, +), dans une bibliothèque à part, `wsl::mdindex` (pan 3), distincte des positions de tableaux ; elle reprendra ou remplacera `graded` (jets, développements de Taylor tronqués). |
| A7 | **Plusieurs bibliothèques, pas une** : chaque pan devient une bibliothèque du dépôt (en-tête, cible CMake, dépendances déclarées) ; noms décidés pan par pan. |
| A8 | **Des lots qui livrent du code utilisable**, testé et poussé ; la conception ne devance le code que d'un pan. |
| A9 | **Structures algébriques par déclaration explicite** : un type n'est un monoïde, un groupe… que s'il le déclare (trait à spécialiser, qui donne aussi l'élément neutre) ; le concept vérifie la syntaxe en plus. Raison : la syntaxe ne dit pas quelle structure un type modélise (`size_t` avec `+` est ℤ/2⁶⁴ℤ pour la machine, ℕ pour l'intention). Dans une bibliothèque séparée (nom proposé : `wsl::algebraic`) : magma, demi-groupe, monoïde, groupe, actions de groupe, espaces affines, puis anneaux, corps, modules, espaces vectoriels. Les types de `wsl::mdarray` y sont déclarés par un en-tête de liaison, sans dépendance entre les deux bibliothèques. |

## 3. Feuille de route

| Pan | Contenu | Débloque |
|---|---|---|
| 1. Algèbre linéaire | corps des scalaires ; espaces vectoriels (Kⁿ et sa base canonique, espaces généraux, sous-espaces, familles, bases, dimension) ; applications linéaires L(E, F), matrices dans des bases, noyau, image, rang, changement de bases ; dual ; réduction ; formes bilinéaires et produits scalaires, orthogonalité ; espaces normés de dimension finie, normes subordonnées ; espace affine ; en support : structures algébriques (`wsl::algebraic`) ; stockage des composantes par `wsl::mdarray` | calcul en coordonnées, systèmes linéaires, valeurs propres ; points et vecteurs en mécanique |
| 2. Tenseurs | produit tensoriel, composantes, variances, contractions, notation indicielle (convention d'Einstein), symétries | tenseurs des contraintes et des déformations ; métrique, Riemann |
| 3. Calcul différentiel | fonctions entre espaces, différentielle, dérivées d'ordre supérieur, ∂^α, Taylor, multi-indices (`wsl::mdindex`), différentiation automatique (les jets de `graded` sont des développements de Taylor tronqués) | gradients, jacobiennes, hessiennes ; opérateurs grad, div, rot, Δ |
| 4. Au-delà (esquisse) | espaces de fonctions et opérateurs, discrétisations (différences finies, éléments finis), variétés, cartes, champs, connexions | équations de la relativité générale, de la MMC et des fluides |

Premier lot visé : scalaires, Kⁿ et espaces vectoriels, applications linéaires et leurs
matrices, en dimension fixe (A(u), B ∘ A, déterminant, inverse).

## 4. Questions ouvertes

- Syntaxe des espaces : Kⁿ (base canonique, vecteurs = n-uplets) et espaces nommés
  (distincts même de même dimension, coordonnées seulement dans une base explicite).
- Application et composition : `A(u)` et `B * A` (proposé), ou style matriciel `A * u`.
- Dérivées : opérateurs sur les fonctions (`D(f)`, `grad(f)`, `partial(f, α)`), et/ou en un
  point.
- Composantes : notation indicielle nommée (`C[i, k] = A[i, j] * B[j, k]`) et fonctions
  nommées (`contract`, `trace`).
- Familles de fonctions indexées (`m[i, j](x, y)` et `m(x, y)[i, j]`, curryfications de
  Y^(I×X)) ; diffusion entre rangs (aucune implicite hors scalaires, proposé).
- Place de `square_matrix` (wsl::random) ; `array2d` devient un alias de rang 2 de `wsl::mdarray`.
