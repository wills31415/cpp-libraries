# wsl::random — variables aléatoires (proposition)

> Document de conception, 2026-11-21. ✅ : décidé (V1 à V6, le 2026-11-21). Prolonge [random.md](random.md) (les lois) et [interfaces.md](interfaces.md)
> (décisions 8b, 8c : `X(g)` évalue une variable en ω, `X + Y` est réservé aux variables).
>
> Écriture : `wr::variable X{wr::normal()};` (accolades : `wr::variable X(wr::normal());`
> déclarerait une fonction).

## 1. Le modèle mathématique

Une **loi** est une mesure de probabilité P sur un espace E (ce que wsl::random fournit
aujourd'hui). Une **variable aléatoire** est une application mesurable X : Ω → E sur un
espace probabilisé (Ω, 𝓕, ℙ) commun ; sa loi est la mesure image ℙ ∘ X⁻¹.

Ce que les variables ajoutent aux lois, c'est la **dépendance** : X + X n'est pas X + Y
pour Y indépendante de même loi (2X contre la convolution P ∗ P), et X − X vaut 0 en tout
ω. Les lois ne peuvent pas l'exprimer, puisqu'elles ne voient pas Ω.

Réalisation : Ω est un produit d'espaces, un facteur par **variable de base** ;
`wr::variable{P}` ajoute un facteur de loi P (donc une variable indépendante de toutes
les autres) ; les expressions sont des fonctions des variables de base ; évaluer en ω,
c'est tirer chaque variable de base une fois, puis calculer l'expression.

```cpp
namespace wr = wsl::random;
wr::variable X{wr::normal()};                     // X ~ N(0, 1)
wr::variable Y{wr::exponential()};                // Y ~ E(1), indépendante de X
auto Z = 2 * X + Y;                               // expression paresseuse
Z(g);                                             // un ω : X et Y tirées une fois chacune
(X - X)(g);                                       // 0, en tout ω
auto [x, z] = wr::evaluate(g, X, Z);              // même ω : z, arrondi exact de 2X + Y
wr::law(2 * X + 1);                               // normal({.mean = 1, .stddev = 2}), exacte
wr::law(X + X);                                   // normal({.variance = 4}) : 2X, pas P ∗ P
wr::law(X + wr::variable{wr::normal()});          // normal({.variance = 2}) : convolution
wr::E(Z);                                         // 1 (linéarité, sans tirage)
```

## 2. Décisions à prendre

### V1. Identité d'une variable ✅ (A, 2026-11-21)

Une variable de base a une **identité** : deux copies du même objet désignent la même
fonction de ω (copier `X` ne crée pas une variable indépendante).

| Option | Écriture | Pour | Contre |
|---|---|---|---|
| A. Identifiant à l'exécution (recommandé) | `variable(P)` reçoit un numéro unique, croissant | naturel, copies partagées, dynamique (variables en boucle) | l'identité ne se voit pas dans le type |
| B. Identité dans le type | `variable(P, [] {})` (lambda : type unique) | tout résolu à la compilation | écriture lourde, pas de variables créées en boucle |

Copie indépendante explicite : `X.independent_copy()` (même loi, nouvelle identité).

### V2. Évaluation et ω ✅ (2026-11-21)

- `Z(g)` : un ω neuf à chaque appel ; chaque variable de base de Z est tirée **une fois**
  (mémorisée pendant l'évaluation), dans l'ordre croissant des identifiants (ordre de
  création) : le nombre de valeurs consommées et l'ordre font partie du contrat.
- `wr::evaluate(g, X, Y, …)` : plusieurs expressions sur le **même** ω (tuple).
- `Z(n, g)` : n évaluations indépendantes (comme pour les lois, décision 8b).
- Variante : tirer **toutes** les variables de base existantes à chaque ω (plus simple,
  mais coûteux et dépendant de variables étrangères à l'expression) : déconseillé.

### V3. Arrondi des expressions ✅ (B, 2026-11-21)

Que rend `(X + Y)(g)` ?

| Option | Sens | Pour | Contre |
|---|---|---|---|
| A. Calcul flottant | chaque variable de base arrondie, puis l'expression en `T` (comme du code ordinaire) | simple, rapide, toute expression | `(X + Y)(g)` n'a pas exactement la loi `law(X + Y)` (deux arrondis) |
| B. Exact quand c'est possible (recommandé) | combinaisons affines de variables à tirage exact : arrondi correct de la valeur exacte (réels paresseux, `lazy_sum`) ; ailleurs, calcul flottant documenté | cohérent avec « exact par défaut » : `(2X + Y)(g)` est un tirage exact de sa loi | deux régimes à documenter ; plus lent pour les combinaisons |
| C. Toujours flottant, exact sur demande | A par défaut, `Z(g, wr::exact)` pour B | rapide par défaut | contraire à la règle des lois (exact par défaut, `wr::fast` sur demande) |

Avec B, `Z(g, wr::fast)` donne le calcul flottant (A).

### V4. Loi d'une expression ✅ (2026-11-21)

`wr::law(Z)` rend la **loi exacte** quand une identité la donne, sinon une loi générique :

- combinaison affine Σ aᵢ Xᵢ + b de variables de base **distinctes** (donc
  indépendantes) : `convolution` des images `affine`, résolue par les identités des
  familles (normales, Cauchy, gamma de même intensité…) ; X + X se réduit d'abord en 2X ;
- f(X) pour une seule variable et f monotone déclarée (`wr::increasing`) : `pushforward` ;
- sinon : loi générique `pushforward(product_measure(P₁, …, Pₙ), f)` — tirage seulement.

Les constantes suivent la linéarité quand elle suffit, sans connaître la loi :
𝔼[Σ aᵢ Xᵢ + b] = Σ aᵢ 𝔼Xᵢ + b ; Var(Σ aᵢ Xᵢ) = Σ aᵢ² Var Xᵢ (variables de base
indépendantes, coefficients regroupés : Var(X + X) = 4 Var X) ; Cov(Z₁, Z₂) de même.

### V5. Opérations ✅ (2026-11-21)

| Opération | Écriture | Lot |
|---|---|---|
| Arithmétique | `+ − × ÷` entre variables et scalaires, `-X` | 1 |
| Fonctions | `wr::apply(f, X, Y…)` ; `exp`, `log`, `sqrt`, `abs`… de wsl::special | 1 |
| Constantes | `wr::E(Z)`, `wr::Var(Z)`, `wr::Cov(Z₁, Z₂)` (linéarité ; sinon par la loi si connue) | 1 |
| Loi | `wr::law(Z)` | 1 |
| Événements | `X > 0`, `X <= Y`… : variables indicatrices (loi de Bernoulli), `&&`, `||`, `!` | 2 |
| Probabilité | `wr::P(X > 0)` : exacte par la loi quand elle est connue | 2 |
| Conditionnement | `wr::given(Z, X > 0)` : évaluation par rejet sur ω ; loi tronquée quand c'est une identité (`truncate`) | 2 |
| Composantes | `V[i]` pour une variable vectorielle (normale multivariée…) | 3 |

### V6. Noms ✅ (2026-11-21 : notation des manuels)

`wr::variable`, `wr::evaluate`, `wr::law`, `wr::apply`, `X.independent_copy()` ; pour les
constantes et les probabilités, la notation des manuels : `wr::E(Z)`, `wr::Var(Z)`,
`wr::Cov(Z₁, Z₂)`, `wr::P(X > 0)` (objets-fonctions ; exception assumée au `snake_case`,
ce sont des symboles mathématiques ; seul risque : une macro `E` ou `P` chez l'utilisateur).
`wr::given(Z, X > 0)` pour le conditionnement.

## 3. Ce qui ne change pas

Les lois restent des mesures, sans opérateurs arithmétiques (décision 8c) ; `affine`,
`convolution`… restent les opérations sur les lois. Une variable n'est qu'une couche
au-dessus : elle ne modifie aucun tirage existant.

## 4. Réalisation

- Lot 1a ✅ : `wr::variable X{P}` (écriture à accolades : `wr::variable X(P());` déclarerait
  une fonction), expressions, évaluation couplée, `wr::evaluate`, `wr::E`, `wr::Var`,
  `wr::Cov` par linéarité.
- Lot 1b ✅ : combinaisons affines exactes. Dans une expression, chaque variable de base à
  tirage paresseux est tirée comme réel paresseux, partagé par toutes les expressions du
  même ω ; une combinaison affine en est l'arrondi correct (`lazy_sum`), les coefficients
  d'une même variable regroupés, les divisions par des constantes gardées exactes comme
  diviseurs, les constantes réduites par dénominateur. Repli sur le calcul flottant (sur
  des valeurs tirées exactement) dès qu'un coefficient n'est pas représentable exactement.
  Une variable de base seule, `X(g)`, est le tirage de sa loi.

## 5. Lots proposés

1. Variables de base, expressions arithmétiques et fonctions, évaluation couplée,
   `law` des combinaisons affines, constantes par linéarité.
2. Événements, probabilités, conditionnement.
3. Variables vectorielles, processus vus comme familles de variables.
