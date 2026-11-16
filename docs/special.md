# wsl::special — fonctions mathématiques bit-exactes

> Document de conception. Décisions ✅ prises avec l'utilisateur le 2026-11-16 ;
> les propositions ◻︎ ont été validées en bloc le même jour.
>
> **État** : phase 1 en cours — `exp`, `exp2`, `expm1`, `log`, `log2`, `log10` et
> `log1p` (float, double) implémentés et vérifiés (§6).

## 1. Pourquoi cette bibliothèque

La bibliothèque d'aléatoire ([random.md](random.md)) doit donner **les mêmes tirages au
bit près sur toutes les plateformes**. Or :

- les fonctions de `<cmath>` ne sont pas identiques d'une bibliothèque standard à
  l'autre. Mesuré sur 1 000 000 d'entrées identiques, macOS (libc++, libm d'Apple)
  contre Linux (glibc) :

  | Fonction | Résultats différents |
  |---|---|
  | `std::log` | 0,14 % |
  | `std::exp` | 0,22 % |
  | `std::pow` | 0,16 % |
  | `std::sin` | 15,7 % |
  | `std::lgamma` | 41,8 % |

- libc++ ne fournit pas les fonctions spéciales de C++17 (`std::beta`, `std::expint`,
  `std::riemann_zeta`, `std::cyl_bessel_*`, `std::assoc_legendre`… : vérifié absentes),
  et aucune bibliothèque standard ne fournit les fonctions gamma et bêta **incomplètes**
  ni leurs inverses, indispensables aux fonctions de répartition et aux quantiles des
  lois gamma, χ², bêta, Student et Fisher.

✅ Ces fonctions forment une bibliothèque à part entière, utilisable sans la
bibliothèque d'aléatoire (analyse numérique, physique, cpp-tensors…).

## 2. Garanties

### 2.1 Reproductibilité ✅

Pour `float` et `double`, chaque fonction rend **le même résultat au bit près** sur
toute plateforme conforme (voir 2.3), quels que soient le compilateur, la bibliothèque
standard et le niveau d'optimisation.

Moyen : les fonctions sont écrites uniquement avec les opérations dont IEEE 754 impose
l'arrondi correct — `+ − × ÷`, `sqrt`, `fma` — et des constantes exactes. Jamais
d'appel à `<cmath>` pour une fonction transcendante.

Les résultats font partie du contrat : changer l'algorithme d'une fonction (même pour
la rendre plus précise) change des résultats, c'est donc une rupture, annoncée dans le
`CHANGELOG.md` par une version mineure tant que la version majeure vaut 0.

### 2.2 Précision ◻︎

- Cible des fonctions élémentaires : erreur **< 1 ulp** (arrondi fidèle) sur tout leur
  domaine, avec l'**arrondi correct** visé à terme pour `exp` et `log`.
- Fonctions spéciales : erreur maximale mesurée et documentée, fonction par fonction
  et par intervalle, contre des valeurs de référence en précision arbitraire.
- Les erreurs mesurées sont publiées dans la documentation de chaque fonction et
  vérifiées par les tests (une régression de précision fait échouer la CI).

### 2.3 Conditions de la garantie

Ces conditions sont vérifiées à la compilation ou neutralisées par les en-têtes :

| Risque | Vérifié | Parade |
|---|---|---|
| `-ffast-math` (réassociation, `x * (1/y)`…) | `__FAST_MATH__`, `__FINITE_MATH_ONLY__` | `#error` |
| Calcul en précision étendue (x87, 32 bits) | `FLT_EVAL_METHOD != 0` | `#error` |
| **Contraction `a*b + c` → `fma`** | — | Barrière sur chaque produit (voir ci-dessous) |
| `long double` (64, 80 ou 128 bits selon la plateforme) | — | Refusé à la compilation |

La contraction est le piège le plus sournois, vérifié sur la même ligne de code :

| Compilation de `a * b + c` | Instructions |
|---|---|
| Clang, x86-64 | multiplication puis addition |
| Clang, ARM64 (Mac Apple Silicon, macOS de la CI) | `fmadd` |
| GCC 14, `-march=haswell`, même en `-std=c++20` strict | `vfmadd` |

Les bibliothèques étant en en-têtes seuls, elles sont compilées avec les options **du
projet consommateur** : la parade doit être dans les en-têtes eux-mêmes. Les directives
de compilation ne suffisent pas (vérifié) :

| Parade | Clang, options par défaut | Clang, `-ffp-contract=fast` | GCC 14 |
|---|---|---|---|
| `#pragma clang fp contract(off)` | tient | **ne tient pas** (le générateur de code de LLVM fusionne quand même) | — |
| `#pragma GCC optimize("fp-contract=off")` | — | — | tient, mais **empêche l'inlining** dans le code consommateur |
| **Barrière** `asm("" : "+x"(p))` sur le produit | tient | tient | tient |

La barrière retenue (`detail::fmul`, dans `detail/fp.hpp`) rend le produit opaque pour
l'optimiseur : rien ne peut être fusionné à travers. Elle ne coûte aucune instruction,
n'empêche pas l'inlining, et s'efface en évaluation constante (`constexpr`), qui ne
fusionne jamais. Règle d'écriture : **tout produit dont le résultat entre dans une
addition ou une soustraction passe par `fmul`**. Un oubli serait détecté par le test
d'empreintes dans les configurations FMA de la CI.

Là où une `fma` serait voulue, elle s'écrirait explicitement (`std::fma`) ; la phase 1
n'en utilise aucune.

## 3. Interface ◻︎

- Cible CMake `wsl::special`, en-tête parapluie `<wsl/special.hpp>`, en-têtes par
  famille sous `<wsl/special/…>` (`elementary.hpp`, `gamma.hpp`, `error.hpp`…).
- Namespace **`wsl::special`** : les noms de fonctions (`beta`, `gamma`, `zeta`…) sont
  aussi des noms de lois dans `wsl` ; les séparer évite toute collision. Employer des
  noms qualifiés ou un alias (`namespace sp = wsl::special;`) : un
  `using namespace wsl::special;` rendrait `exp(2.0)` ambigu avec le `::exp` de
  `<math.h>`.
- Noms et sémantique **de `<cmath>`** quand la fonction y existe (`exp`, `log1p`,
  `lgamma`, `tgamma`, `erf`, `beta`, `expint`, `riemann_zeta`, `cyl_bessel_i`…),
  y compris le traitement des cas limites (±0, ±∞, NaN, valeurs sous-normales) tel
  qu'Annexe F de C le décrit. Noms des manuels sinon (`gamma_p`, `gamma_q`, `beta_inc`,
  `erf_inv`…).
- Domaine : hors domaine, NaN comme `<cmath>` (pas d'exception : ce sont des fonctions
  numériques, pas des constructeurs de lois).
- `constexpr` partout où c'est possible : écrites en arithmétique pure, ces fonctions
  peuvent être évaluées à la compilation, ce que `<cmath>` ne permet pas en C++23.
- Surcharges `float` et `double`, plus les entiers (traités en `double`, comme
  `<cmath>`) ; `long double` est supprimé (`= delete`).
- En phase 1, la version `float` est calculée en `double` puis arrondie une fois :
  déterministe, et presque toujours l'arrondi correct.

## 4. Tests ◻︎

1. **Précision** : tables de valeurs de référence calculées par `mpmath` (précision
   arbitraire), générées par un script de `tools/reference/` lancé avec
   `uv run --with mpmath` et versionnées. Points choisis aux bords de domaine, autour
   des zéros, des pôles, des changements d'algorithme, plus des points aléatoires.
   Le test vérifie l'erreur en ulp de chaque point contre la borne documentée.
2. **Bit-exactitude** : empreinte (hachage) des résultats sur un grand ensemble fixe
   d'entrées, comparée à une valeur versionnée. Identique sur les trois environnements
   de la CI (x86-64 Linux GCC et Clang, ARM64 macOS) : c'est la preuve de la garantie
   2.1, et non une supposition.
3. **Cas limites** : ±0, ±∞, NaN, sous-normaux, plus grand et plus petit fini.
4. **Performance** : banc d'essai comparé à `<cmath>`, informatif (n'échoue pas).

## 5. Contenu, par phases

| Phase | Fonctions | Utilisées par |
|---|---|---|
| **1** | `exp`, `exp2`, `expm1`, `log`, `log2`, `log1p`, `pow`, `sin`, `cos`, `tan`, `atan`, `atan2`, `cbrt`, `hypot` | Toutes les lois continues : exponentielle, normale (Box–Muller), Cauchy, sphère, von Mises… |
| **1** | `lgamma`, `tgamma`, `erf`, `erfc`, `erf_inv`, `erfc_inv` | Densités de gamma, bêta, Poisson ; normale (répartition, quantile) |
| **4** | `gamma_p`, `gamma_q` (gamma incomplètes régularisées) et inverses ; `beta`, `lbeta`, `beta_inc` et inverse ; `digamma` | Répartition et quantiles : gamma, χ², Poisson, bêta, binomiale, Student, Fisher |
| plus tard | Bessel (`cyl_bessel_i`…), `riemann_zeta` et zeta de Hurwitz, `expint`, Owen T, fonctions de Lambert, polylogarithme… | von Mises, Zipf, Skellam, normale asymétrique, Rice… |

L'ordre dans la phase 1 suit les dépendances : `exp`/`log` d'abord (tout le reste
s'appuie dessus), puis `pow`, puis les fonctions trigonométriques (réduction d'argument
exacte), puis `lgamma`/`erf`.

## 6. État et mesures

Phase 1, au 2026-11-16 : famille exponentielle et logarithme, en `float` et en `double`.

**Précision**, contre `mpmath` (2 500 à 5 600 points en double, 1 200 à 3 100 en float,
choisis sur tout le domaine, aux bords des sous-intervalles et des tables, près des
zones délicates — petits arguments de `expm1` et `log1p`, voisinage de −1 pour `log1p`
—, sous-normaux compris) :

| Fonction | Erreur max (double) | Arrondis corrects (double) | Erreur max (float) |
|---|---|---|---|
| `exp` | 0,504 ulp | 99,89 % | 0,500 ulp (100 %) |
| `exp2` | 0,502 ulp | 99,92 % | 0,500 ulp (100 %) |
| `expm1` | 0,500 ulp | 99,92 % | 0,500 ulp (100 %) |
| `log` | 0,500 ulp | 100 % | 0,500 ulp (100 %) |
| `log2` | 0,500 ulp | 100 % | 0,500 ulp (100 %) |
| `log10` | 0,500 ulp | 100 % | 0,500 ulp (100 %) |
| `log1p` | 0,500 ulp | 100 % | 0,499 ulp (100 %) |

Exactitudes garanties par les tests : `exp2(n) = 2ⁿ` et `log2(2ⁿ) = n` pour tout
entier n de [−1074, 1023], `log10(10ⁿ) = n` pour n de 0 à 22 (les puissances de 10
exactes en double).

**Bit-exactitude** : les 14 empreintes (7 fonctions × 2 types) sont identiques sur macOS
x86-64 et ARM64 (Apple Clang), Linux x86-64 (GCC 14 et Clang 18, Debug et Release) et,
configurations hostiles, avec `-march=x86-64-v3 -ffp-contract=fast` (GCC, Clang) ou
`-ffp-contract=fast` (ARM64).

Organisation du code : un noyau par famille (`reduce_exp`, `exp_finish`, `scale_exp` ;
`log_kernel`, qui rend le logarithme en double-double), que les fonctions prolongent.
`expm1` garde la réduction en double-double et ajoute exactement les termes qui
s'annulent ; `log1p` évalue directement le polynôme pour |x| < 2⁻¹⁰ (passer par
1 + x y coûtait jusqu'à 0,57 ulp), et sinon corrige log(1 + x) de l'erreur exacte de
1 + x ; `log2` et `log10` multiplient le logarithme double-double par 1/ln 2 ou
1/ln 10 en double-double.

**Coût**, en nanosecondes par appel, contre la libm d'Apple (Intel i7-9750H,
`tools/special/benchmark.sh`, minimum de 9 mesures, un programme par fonction ; les
chiffres varient de 10 à 20 % d'une exécution à l'autre) :

| Fonction | wsl::special | `<cmath>` | Rapport |
|---|---|---|---|
| `exp` double | 6,5 | 5,3 à 6,0 | ×1,1 à 1,2 |
| `log` double | 8,8 à 10,7 | 4,7 à 5,3 | ×1,6 à 2,3 |
| `exp` float | 8,4 | 4,4 | ×1,9 |
| `log` float | 9,8 | 4,5 à 4,9 | ×2,0 à 2,3 |

`log` a été accéléré d'environ 30 % (13 → 9 à 10 ns), sans changer sa précision :
table deux fois plus fine (polynôme de degré 7 au lieu de 8), évaluation d'Estrin au
lieu de Horner, `log(c)` tabulé sur une grille de 2⁻⁴² pour que `e·ln2 + log(c)` soit
exact (une somme exacte de moins), réduction et indice sans branche ni conversion.

Les mesures sont faites par des programmes séparés : dans un même programme, de petites
variations de la boucle de mesure ont faussé les chiffres d'un facteur 3 à 5, de façon
reproductible et pour `<cmath>` aussi (constaté sur Intel).

Piste restante : pour `float`, des algorithmes propres (encore évalués en `double`, mais
avec moins de termes). Toute optimisation qui change des résultats est une rupture
(§2.1) : à faire avant la version 1.0.

**Vérification de la règle `fmul`** : le test `special_no_contraction` compile les
fonctions avec fma disponible et `-ffp-contract=fast`, et refuse toute instruction fma
dans le code engendré (contre-épreuve à l'appui). Limite : un compilateur ne fusionne
qu'un produit à usage unique ; un oubli sur un produit utilisé plusieurs fois lui
échappe (constaté : c'était le cas de `r_hi` dans la première version de `log`) — il
ne change alors aucun résultat, et les empreintes le verraient si un compilateur futur
le fusionnait.
