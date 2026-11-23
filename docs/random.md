# wsl::random — les lois de probabilité, du point de vue mathématique

> Document de conception. Décisions ✅ prises avec l'utilisateur à partir du 2026-11-16 ;
> ◻︎ : proposition ou piste non réalisée.
>
> **État** (2026-11-21) : phases 1 à 8 faites (fonctions spéciales, socle, lois univariées,
> algèbre des lois, lois multivariées, objets combinatoires, processus) ; la **revue des
> interfaces** ([interfaces.md](interfaces.md), §6) est appliquée : les §1 à §7 décrivent
> l'interface actuelle, les §8 à §21 sont le journal des phases (écritures mises à jour).

## 1. Vision ✅

Pouvoir tirer selon **toutes** les lois de probabilité — de l'uniforme à Dirichlet et
au-delà — en se plaçant du point de vue du probabiliste :

1. **Une loi est un objet mathématique complet**, pas seulement un générateur : elle a
   un support, une densité ou une fonction de masse, une fonction de répartition, des
   quantiles, des moments, une entropie… Le tirage n'est qu'une opération parmi d'autres.
2. **Une algèbre des lois** : construire des lois à partir d'autres (image par une
   transformation affine ou une fonction, troncature et conditionnement, mélanges,
   produits de lois indépendantes), avec les identités connues quand elles existent.
3. **Exactitude** : des tirages exacts quand c'est possible, sans approximation
   flottante cachée ; une précision documentée sinon.
4. **Le vocabulaire des manuels** : lois et paramétrisations nommées comme dans la
   littérature, paramètres vérifiés.

« Toutes les lois » est une liste sans fin : la bibliothèque n'est donc pas une
collection de lois écrites une à une, mais une **structure** où ajouter une loi coûte
peu. Une loi déclare ce que l'on sait d'elle ; la bibliothèque en déduit le reste.

## 2. Décisions

| Sujet | Décision | |
|---|---|---|
| Périmètre | Lois univariées (discrètes et continues), multivariées, objets combinatoires, processus | ✅ |
| Reproductibilité | **Bit-exacte sur toutes les plateformes**, lois continues comprises (voir §3) | ✅ |
| Fonctions spéciales | Bibliothèque à part, [`wsl::special`](special.md) | ✅ |
| Interface | Objets mathématiques **et** adaptateur vers le modèle de `<random>` | ✅ |
| Espace de noms | `wsl::random::` (internes `wsl::random::detail`) | ✅ |
| Paramètres | Nommés et vérifiés : `gamma({.rate = 2.0, .shape = 2.5})` ; `std::domain_error` hors domaine ; ordre position → échelle → formes, défauts position 0, échelle 1 | ✅ |
| Opérations | Membres (`X.cdf(x)`, `X(g)`) et objets-fonctions qualifiés (`wsl::random::cdf(X, x)`) | ✅ |
| Exactitude | **Exact par défaut** quand un algorithme exact est connu ; variante rapide sur demande explicite | ✅ |
| Types flottants | `float` et `double` ; `long double` refusé à la compilation | ✅ |
| Moteurs | Tout moteur standard, plus xoshiro256++, PCG64, Philox | ✅ |
| Dimension (multivarié) | Fixe ou dynamique, sur le modèle de `std::span` (`dirichlet<double, 3>`, `dirichlet<double, wsl::random::dynamic>`) ; vecteurs `std::array` ou `std::vector`, matrices `wsl::random::square_matrix` (§12) | ✅ |
| Processus | Vues paresseuses (`std::ranges`) : `wsl::random::trajectory(P, g)` | ✅ |

## 3. Reproductibilité ✅

**Garantie** : même moteur, même graine ⇒ mêmes tirages au bit près, sur toute
plateforme conforme (conditions dans [special.md §2.3](special.md)).

Pourquoi ne pas s'appuyer sur `<random>` : la norme spécifie complètement les moteurs
(`mt19937`…), mais pas les algorithmes des lois. Vérifié, même moteur `mt19937_64` de
graine 42 :

| | `normal_distribution` | `gamma_distribution(2.5)` | `uniform_int_distribution(1, 6)` |
|---|---|---|---|
| macOS (libc++) | 1.293820423273 | 3.041677468034 | 1 |
| Linux (libstdc++) | 0.704988266421 | 3.378867072039 | 5 |

Deux conséquences :

- **Aucune loi de `<random>` n'est utilisée**, et aucune fonction transcendante de
  `<cmath>` : tout passe par `wsl::special` (une différence d'un bit dans un `log` suffit
  à faire basculer un test d'acceptation, donc à désynchroniser le flux entier).
- **Un algorithme de tirage fait partie du contrat.** Le changer — même pour
  l'accélérer — change les tirages : c'est une rupture, annoncée dans le `CHANGELOG.md`.
  Le nombre de valeurs consommées sur le moteur fait aussi partie du contrat.

## 4. Architecture

```
wsl::special   fonctions bit-exactes (voir special.md)
     ↑
wsl::random
  ├─ engines      moteurs : ceux de <random> + xoshiro256++, PCG64, Philox, SplitMix64
  ├─ uniform      bits → entiers uniformes exacts, réels uniformes exacts, Bernoulli exacte
  ├─ laws         les lois : chacune déclare support, paramètres, et ce qu'on sait d'elle
  ├─ algebra      affine, truncate, mixture, convolution, product_measure, pushforward
  ├─ processes    marches, chaînes de Markov, Poisson, brownien, urnes (vues paresseuses)
  ├─ adaptor      as_distribution : une loi vue comme une distribution de <random>
  └─ (generic)    ◻︎ tirage déduit de ce qui est connu (inversion, rejet…), non fait
```

### 4.1 Moteurs ✅

Tout type satisfaisant `std::uniform_random_bit_generator` est accepté. La bibliothèque
fournit en plus :

| Moteur | Intérêt |
|---|---|
| `xoshiro256pp` | Très rapide, excellente qualité statistique, saut de 2¹²⁸ pour des flux disjoints |
| `pcg64` | Qualité statistique, flux sélectionnables |
| `philox4x64` | À compteur : la *n*-ième valeur s'obtient directement ; flux parallèles indépendants par construction (calcul réparti, GPU) |
| `splitmix64` | Initialisation des autres à partir d'une graine de 64 bits |

Tous ont un état et une initialisation entièrement spécifiés, donc reproductibles, et sont
complets au sens de *RandomNumberEngine* : `E()`, `E(s)`, `E(seed_seq)` et `seed(…)`
(mots de 32 bits assemblés deux à deux, poids faible d'abord, comme `mt19937_64`),
`discard(n)`, `==`, `os << e` et `is >> e` (décimal ; lecture invalide : moteur intact,
`failbit`).

### 4.2 Uniformes, la fondation ✅ (phase 2, §8)

Tout le reste en dépend ; c'est là que l'exactitude commence. Le moteur est toujours le
dernier argument, comme dans la norme.

- **Entier uniforme** : `uniform_below(n, g)` sur [0, n) (Lemire : multiplication et rejet,
  presque sans division), `uniform_int(a, b, g)` sur [a, b] ; exactement uniformes.
- **Réel uniforme sur [0, 1)**, exact : `uniform_real<T>(g)`. Chaque `double` de [0, 1) est
  atteint avec exactement la probabilité de son intervalle d'arrondi, y compris les plus
  petits (exposant tiré selon une loi géométrique, puis mantisse). Le classique
  `(x >> 11) * 2⁻⁵³` n'atteint que 2⁵³ valeurs régulièrement espacées.
- **Bernoulli exacte** : `bernoulli_trial(p, g)` (p flottant, lu bit à bit : 2 bits en
  moyenne) et `bernoulli_trial(num, den, g)` (p rationnel).
- **Mélange de Fisher–Yates** exact : `shuffle(range, g)`.

### 4.3 Le concept de loi ✅

Une loi est un type qui hérite de `law_interface` et déclare :

- `value_type` : `double`, `std::int64_t`, `std::array<double, 3>`, une permutation… ;
- ses **paramètres**, dans une structure nommée `parameters` (avec `param()` et
  `param_type` : `L(X.param()) == X`), vérifiés à la construction ;
- son **support** (`X.support()`) et sa **mesure de référence** (`lebesgue`, `counting`,
  `surface`…) ;
- un **niveau d'exactitude** du tirage (`sampling` : `exact`, ou `approximate` avec la
  borne d'erreur documentée) ;
- et, **quand ils sont connus**, ses autres attributs mathématiques, chacun détecté par
  un concept (`has_cdf<L>`…) :

| Attribut | Membre | Objet-fonction | Exemple |
|---|---|---|---|
| Densité (dérivée de Radon–Nikodym dP/dμ) | `X.density(x)`, `X.log_density(x)` | `wsl::random::density(X, x)` | toutes les lois |
| Densité / fonction de masse (noms usuels) | `X.pdf(x)`, `X.log_pdf(x)` / `X.pmf(k)`, `X.log_pmf(k)` | `wsl::random::pdf(X, x)`… | lois continues / discrètes |
| Fonction de répartition, survie | `X.cdf(x)`, `X.ccdf(x)` | `wsl::random::cdf(X, x)`… | normale, gamma… |
| Quantile, quantile de survie | `X.quantile(p)`, `X.cquantile(q)` (x tel que P(X > x) = q, sans former 1 − q) | `wsl::random::quantile(X, p)`… | lois continues |
| Tirage | `X(g)`, `X(g, wsl::random::fast)`, `X(n, g)` (vecteur, mêmes valeurs et même consommation que n appels) | `wsl::random::sample(X, g)`, vue `wsl::random::samples(X, g)` | toutes |
| Moments | `X.mean()`, `X.variance()`, `X.skewness()`, `X.excess_kurtosis()` | `wsl::random::mean(X)`… | absents pour Cauchy |
| Mode, médiane, entropie | `X.mode()`, `X.median()`, `X.entropy()` | | |

Constantes (décisions 3a, 3b) : en `T` pour les lois continues ; pour les lois discrètes,
médiane, mode, quantiles et support en `I`, moments, masse et répartition en `double`. Un
attribut **jamais défini** (l'espérance d'une Cauchy) est absent : l'appel ne compile pas.
Un attribut défini selon les paramètres suit IEEE, sans exception : +∞ quand l'intégrale
diverge vers +∞ (moyenne de Pareto α <= 1), NaN quand elle n'a pas de sens (moyenne de
Student ν <= 1, asymétrie d'une loi dégénérée).

Paramètres exacts (décision 5) : une seule forme canonique par loi, que `==` compare
(égalité en loi exacte ; deux objets égaux donnent les mêmes bits, tirages compris).
Gamma, exponentielle et Erlang stockent une intensité rationnelle λ = P/Q·2ᵉ, donnée par
`.rate` ou `.scale` (double ou `wsl::random::ratio(a, b)`) ; la normale stocke sa variance
exacte, donnée par `.stddev` ou `.variance` (√v irrationnel : σ en double-double pour les
fonctions, tirage exact par comparaison des carrés).

Exemple d'utilisation :

```cpp
namespace wr = wsl::random;
wr::gamma X({.rate = 2.0, .shape = 2.5});
wr::gamma Y({.scale = 0.5, .shape = 2.5});        // la même loi : X == Y
wr::normal Z({.mean = 0, .variance = 4});         // N(μ, σ²) des manuels, == {.stddev = 2}
wr::xoshiro256pp g(42);

X.mean();                // 1.25
X.cdf(1.0);              // P(X <= 1)
Z.quantile(0.975);       // 3.919…
double x = X(g);         // tirage (exact si la loi le permet)
auto xs = X(1000, g);    // 1000 tirages

wr::gamma({.rate = 1.0, .scale = 1.0, .shape = 2.0});   // std::domain_error : .rate ou .scale, pas les deux
wr::gamma({.rate = 1.0, .shape = -1.0});                 // std::domain_error : shape doit être > 0
```

Pourquoi des paramètres nommés : `gamma(2.5, 2.0)` a pour espérance 5 en C++ et en
NumPy (second paramètre = échelle), 1,25 en R (second paramètre = intensité). En
positionnel, se tromper de convention donne une loi fausse **sans aucune erreur**.

### 4.4 Tirage générique ◻︎ (non fait)

Une loi peut fournir son propre algorithme de tirage ; sinon, la bibliothèque en déduit
un de ce qui est connu :

| Ce qui est connu | Tirage déduit |
|---|---|
| Quantile | Inversion : `quantile(X, U)` |
| Densité bornée sur un support borné | Rejet sous un rectangle |
| Densité log-concave | Rejet adaptatif (ARS) |
| Fonction de masse sur un support fini | Méthode de l'alias (Walker–Vose), exacte avec des poids rationnels |

Les tirages déduits seraient marqués `approximate` (le quantile ou la densité est calculé
en flottant), avec la précision qui en découle. Aujourd'hui, chaque loi fournit son
algorithme ; l'alias de Walker–Vose sert `categorical`.

### 4.5 Exactitude ✅

Un tirage est **exact** quand la valeur rendue est l'arrondi, selon une règle
documentée, d'un réel tiré *exactement* selon la loi : la seule approximation est la
représentation finale en `float` ou `double`.

Exemples d'algorithmes exacts connus : uniforme entière et réelle, Bernoulli
rationnelle, géométrique, exponentielle (von Neumann), normale (Karney 2016),
gaussienne discrète (Canonne–Kamath–Steinke 2020), binomiale et Poisson par des
méthodes exactes, mélange de Fisher–Yates.

Quand un algorithme exact existe, il est utilisé par défaut ; la variante rapide se
demande explicitement (`X(g, wsl::random::fast)`) et porte sa borne d'erreur.
Quand aucun n'est connu, le meilleur algorithme est utilisé et marqué `approximate`.

### 4.6 Algèbre des lois ✅ (phase 5)

| Construction | Écriture | Ce qui est hérité |
|---|---|---|
| Image affine (mesure image par x ↦ ax + b) | `wsl::random::affine(P, {.scale, .shift})` | densité, répartition, quantiles, moments, entropie transformés |
| Troncature, conditionnement | `wsl::random::truncate(P, {.lower, .upper})` (bornes facultatives) | densité, répartition, quantiles renormalisés |
| Mélange | `wsl::random::mixture(wsl::random::weighted(w1, X1), wsl::random::weighted(w2, X2), …)` | densité ou masse, répartition, quantiles, moments |
| Convolution (loi de X + Y, X ⫫ Y) | `wsl::random::convolution(P, Q)` ; différence : `convolution(P, affine(Q, {.scale = -1}))` | identités (normale, Cauchy, Poisson, χ² ; gamma de même intensité, binomiales de même p) ; sinon moments par les cumulants, et pour les lois discrètes masse, répartition, quantile par convolution exacte |
| Mesure produit (loi jointe de X ⫫ Y) | `wsl::random::product_measure(P, Q, …)` | densité ou masse, répartition et survie jointes, entropie, moyennes |
| Image par une fonction | `wsl::random::pushforward(X, f)` | tirage ; avec `wsl::random::increasing(f, f⁻¹ [, (f⁻¹)'])` ou `decreasing` : répartition, quantiles, médiane, densité |

Les identités connues sont résolues **à la compilation** : `affine(normal, {.scale = 2,
.shift = 1})` a pour type `normal`, avec ses paramètres calculés (voir §10). Sans identité,
le résultat est une loi générique qui compose. Les lois n'ont pas d'opérateurs
arithmétiques : `X + Y` est réservé aux futures variables aléatoires (décision 8c).

### 4.7 Adaptateur vers `<random>` ✅

`wsl::random::as_distribution(X)` satisfait l'exigence *RandomNumberDistribution* de la norme
(`result_type`, `param_type`, `operator()(g)`, `operator()(g, p)`, `min`, `max`,
`reset`, `==`, `<<`, `>>`) : utilisable partout où la norme attend une distribution.
La norme exige un `result_type` arithmétique : l'adaptateur est donc réservé aux lois
scalaires.

## 5. Catalogue visé (non exhaustif)

Tout est fait, sauf les lois marquées ◻︎.

| Famille | Lois |
|---|---|
| Discrètes univariées | Bernoulli, binomiale, géométrique, binomiale négative, Poisson, hypergéométrique, uniforme discrète, catégorielle, Zipf, zeta ◻︎, logarithmique, Skellam, bêta-binomiale, gaussienne discrète |
| Continues univariées | uniforme, normale, log-normale, exponentielle, gamma, χ², Erlang, bêta, Student, Fisher, Cauchy, Laplace, logistique, Weibull, Pareto, Gumbel, Fréchet, Rayleigh, Rice ◻︎, von Mises, triangulaire, arcsinus, inverse-gaussienne, normale asymétrique, stables (α, β) |
| Multivariées | Dirichlet, multinomiale, normale multivariée, Student multivariée, Wishart et inverse-Wishart, uniformes sur la sphère, la boule et le simplexe (Dirichlet(1, …, 1)), von Mises–Fisher, matrices aléatoires (Haar orthogonale et unitaire, GOE, GUE) |
| Combinatoires | permutations (uniformes, Ewens, Mallows), sous-ensembles, combinaisons, partitions d'entiers et d'ensembles, compositions, graphes (Erdős–Rényi, Barabási–Albert, graphes réguliers), arbres (Cayley, binaires de Catalan), mots |
| Processus | marches aléatoires, chaînes de Markov (discrètes et continues), processus de Poisson (homogène, inhomogène, composé), mouvement brownien et pont brownien, processus de Dirichlet (par le restaurant chinois, Pitman–Yor), restaurant chinois, urnes de Pólya |

## 6. Tests ✅

Les tirages étant reproductibles, **les tests le sont aussi** : aucun test statistique
« qui échoue une fois sur cent », la graine est fixée.

1. **Empreintes au bit près** : pour chaque loi, l'empreinte (FNV-1a) d'une longue suite
   de tirages et de valeurs de fonctions, sur graine fixée, est versionnée dans les tests.
   Les environnements de la CI (GCC 14 et Clang 18 sur x86-64 Linux, Apple Clang sur ARM64
   macOS, chacun avec et sans FMA ni contraction) doivent les reproduire.
2. **Adéquation** : χ² et Kolmogorov–Smirnov contre la répartition de la loi, sur de grands
   échantillons ; tirages exacts confirmés par raffinement des réels paresseux.
3. **Précision** : densités, répartitions et quantiles comparés à des références mpmath
   (non versionnées, engendrées par CMake), en ulps.
4. **Cohérence mathématique** : `cdf(quantile(p)) ≈ p`, intégrale de la densité ≈ 1,
   moments empiriques ≈ moments théoriques, identités de l'algèbre vérifiées par tirage.
5. **Paramètres** : bords de domaine, valeurs refusées, paramétrisations équivalentes
   (`L(X.param()) == X` pour toutes les lois).

## 7. Feuille de route

| Phase | Contenu |
|---|---|
| 1 | ✅ `wsl::special`, socle : fonctions élémentaires bit-exactes, `lgamma`, `erf`, `erf_inv` ([special.md](special.md)) |
| 2 | ✅ `wsl::random`, socle : moteurs, uniformes exacts, Bernoulli exacte ; `exact_random` est absorbé |
| 3 | ✅ Concept de loi, adaptateur (tirage générique ◻︎) ; premières lois : uniformes, Bernoulli, binomiale, géométrique, Poisson, exponentielle, normale, gamma, bêta |
| 4 | ✅ `wsl::special`, suite : gamma et bêta incomplètes et inverses (répartition et quantiles de gamma, χ², bêta, Student, Fisher) |
| 5 | Algèbre des lois ✅ : image affine, troncature, mélanges, quantile de survie ; sommes indépendantes, loi jointe, image par une fonction |
| 6 | Lois multivariées ✅ (Dirichlet, multinomiale, normales et Student multivariées, Wishart, inverse-Wishart, sphère, boule, von Mises, von Mises–Fisher, Haar, GOE, GUE) |
| 7 | Objets combinatoires ✅ (avec `wsl::permutation`) : permutations, sous-ensembles, combinaisons, mots, partitions d'entiers et d'ensembles, compositions, graphes, arbres |
| 8 | Processus ✅ (temps discret et continu, exacts) |
| Revue | ✅ Revue des interfaces de toutes les bibliothèques du dépôt (2026-11-20 et 21) : décisions dans [interfaces.md](interfaces.md) §6 ; les incohérences du §15 sont résolues |
| (chantiers futurs) | Variables aléatoires (expressions paresseuses, `X + Y`, tirages couplés) ; `mdarray` (avec `graded` pour indice) ; Zipf, Skellam, Rice ; tirage générique |

Ensuite : migration des projets qui utilisent les lois de `<random>` (14 fichiers dans
cpp-traversal, cpp-sudoku et puissance-4, dont `random_agent`, dont la promesse « même
graine, mêmes coups » n'est aujourd'hui tenue que sur une même plateforme). `wsl::elo`
garde `<random>` pour son rééchantillonnage (décision 9).

## 8. Phase 2 : ce qui a été fixé

**Moteurs** (`<wsl/random/engines.hpp>`), tous `std::uniform_random_bit_generator`,
`constexpr`, comparables (`==`), avec `seed` et `discard` :

| Moteur | Initialisation par une graine | Particularités |
|---|---|---|
| `splitmix64` | état = graine | `discard(n)` en O(1) (suite de Weyl) |
| `xoshiro256pp` | 4 sorties de `splitmix64(graine)` | `jump()` (2¹²⁸), `long_jump()` (2¹⁹²) |
| `pcg64` | (graine, flux) à la manière de pcg-cpp | `advance(n)` en O(log n) ; état brut accepté |
| `philox4x64` | clé {graine, 0}, compteur 0 | `block(clé, compteur)` sans état ; rend bloc(c), puis c + 1 |

Vérifiés contre des sources indépendantes de ce code : NumPy (PCG64, Philox, état brut
fixé) ; xoshiro256++ par la première sortie calculable à la main et ses sauts par
élévation de la matrice de transition dans GF(2).

**Primitives** (`<wsl/random/uniform.hpp>`), pour tout générateur standard :

- `random_bits(g)` : 64 bits uniformes ; plage 2^w : concaténation, premier appel en
  poids fort ; autre plage (`minstd_rand`) : rejet sous 2^⌊log2 plage⌋.
- `uniform_below(n, g)` : méthode de Lemire, exacte ; une valeur presque toujours.
- `uniform_int(a, b, g)` : [a, b], jusqu'à la pleine plage de 64 bits.
- `uniform_real<T>(g)` : **U arrondi vers le bas**, chaque flottant x de [0, 1) rendu avec
  probabilité x⁺ − x exactement (sous-normaux compris, 1 jamais) ; exposant géométrique
  lu sur les zéros de tête, mantisse prise dans le même mot quand il en reste assez :
  une valeur consommée presque toujours. `uniform_real_fast<T>(g)` : 53 (24) bits de
  tête, une valeur.
- `bernoulli_trial(p, g)` : U < p par comparaison des bits, mot par mot (une valeur presque
  toujours) ; `bernoulli_trial(num, den, g)` : `uniform_below(den, g) < num`.
- `shuffle(range, g)` : Fisher–Yates, n − 1 tirages.

`wsl::exact_random` est retirée : `wsl::random` la remplace.

## 9. Phase 3 : architecture et premières lois

Décisions du 2026-11-17 : constantes en membres et opérations en fonctions libres ;
exactitude progressive (exactes d'abord : uniformes, Bernoulli, exponentielle, normale ;
les autres en algorithmes classiques bit-exacts marqués `approximate`) ; densités,
répartitions et quantiles à moins d'un ulp, mesurés ; variantes rapides là où l'écart
le justifie.

**Concept** (`<wsl/random/law.hpp>`) : `value_type`, `sampling` (`exact` ou
`approximate`), `support()` ; constantes en membres quand elles existent (`mean`,
`variance`, `skewness`, `excess_kurtosis`, `entropy`, `median`, `mode`) ; opérations en
fonctions libres trouvées par ADL — à appeler **sans qualification**, comme `swap` :
`sample(X, g)`, `sample(X, g, wsl::random::fast)`, `pdf`/`pmf`, `log_pdf`/`log_pmf`, `cdf`,
`ccdf`, `quantile`. Concepts `has_mean`, `has_cdf`… pour tester leur présence. Lois
templates sur le type des valeurs, `double` par défaut (`wsl::random::normal X({…})`,
`wsl::random::normal<float>`). `wsl::random::as_distribution(X)` : adaptateur *RandomNumberDistribution*.

> Remplacé par la revue des interfaces (décisions 8a, 8b) : opérations en membres
> (`X.cdf(x)`, `X(g)`) et objets-fonctions qualifiés (`wsl::random::cdf(X, x)`), plus
> d'appel non qualifié par ADL ; voir §4.3.

**Lois de l'étape A** (`<wsl/random/basic_laws.hpp>`) :

| Loi | Paramètres | Tirage | Attributs |
|---|---|---|---|
| `uniform<T>` | `{.lower = 0, .upper = 1}`, [lower, upper) | approximate : low + U (high − low), un arrondi | moyenne, variance, médiane, asymétrie, aplatissement, entropie ; pdf, log_pdf, cdf, ccdf, quantile |
| `discrete_uniform<I>` | `{.lower, .upper}`, bornes comprises | exact (`uniform_int`) | moyenne, variance, médiane inférieure, asymétrie, entropie ; pmf, log_pmf, cdf, ccdf, quantile exact |
| `bernoulli` | `{.p}` | exact (`bernoulli_trial`) | moyenne, variance, entropie ; pmf, cdf, ccdf |

Précision : calculs en double-double arrondis une fois (`detail/arithmetic.hpp`), mis
à l'échelle pour que les produits exacts le restent ; l'intervalle `[−DBL_MAX, DBL_MAX)`
fonctionne (b − a, qui déborde, est pris sur les moitiés).

Quantile d'une loi discrète : le **quantile exact** de la vraie loi, min{k : F(k) >= p}.
La cdf calculée étant arrondie, `quantile(cdf(k))` vaut k ou k + 1. (Rendre le quantile
cohérent avec la cdf arrondie a été essayé : intenable au-delà de 2⁵³ valeurs, où des
milliers de k voisins ont la même cdf arrondie.)

La primitive `bernoulli(g, p)` de la phase 2 est devenue `bernoulli_trial(p, g)`, pour
laisser le nom à la loi.

### Étape B : exponentielle et normale

| Loi | Paramètres | Tirage (exact) | `wsl::random::fast` |
|---|---|---|---|
| `exponential<T>` | `{.rate}` ou `{.scale}` (exacts, double ou `ratio`) | exponentielle de von Neumann sur réels paresseux, arrondi correct de E/λ | inversion −log(U)/λ |
| `normal<T>` | `{.mean, .stddev}` ou `{.mean, .variance}` (exacts) | demi-normale par rejet exponentiel exact, signe uniforme, arrondi correct de μ + σZ | ziggourat (256 couches) |

Les deux : moyenne, variance, médiane, mode, asymétrie, aplatissement, entropie ;
`pdf`, `log_pdf`, `cdf`, `ccdf`, `quantile`.

**Tirage exact** (`detail/lazy.hpp`). Un réel uniforme paresseux n'existe que par les
mots de 64 bits déjà tirés de son développement binaire ; on en tire d'autres quand
une comparaison ou un arrondi l'exige.
- Exponentielle (von Neumann, 1951) : x est accepté si la suite x > U₁ > U₂ > … remonte
  d'abord à un rang impair (probabilité e^{−x}) ; chaque rejet ajoute 1 à la partie
  entière.
- Demi-normale : E₁ acceptée si une autre exponentielle E₂ dépasse (E₁ − 1)²/2
  (probabilité e^{−(E₁−1)²/2}, le rapport des densités) ; comparaison décidée en
  arithmétique d'intervalles exacte (grands entiers, `detail/big.hpp`), raffinée
  jusqu'à décision. (Plus simple que l'algorithme de Karney, et exact aussi.)
- Arrondi : on raffine jusqu'à ce que les deux bornes de l'intervalle s'arrondissent
  au même flottant (double ou float) : la valeur rendue est l'arrondi correct d'un
  tirage exact de la loi. Une voie rapide en double-double (bornes majorées, marge de
  2⁻³⁰) décide 99,5 % des cas ; elle ne rend un résultat que s'il est certain, si bien
  que résultat **et nombre de mots tirés** sont ceux de la voie exacte.
- Coût mesuré (i7-9750H, xoshiro256++) : exponentielle 79 ns, normale 136 ns ; en
  `wsl::random::fast` 13 ns et 8 ns (`std::exponential_distribution` 15 ns,
  `std::normal_distribution` 18 ns, avec mt19937_64).

**Fonctions de la loi** : z ou λx en double-double (produits exacts, facteurs
rééquilibrés), `exp`, `expm1`, `log1p` à argument double-double, double-double à
exposant séparé (`detail::scaled`) pour ne jamais déborder ni passer par les
sous-normaux avant l'arrondi final, unique. `cdf` de la normale : erfc(w_hi) corrigé
au premier ordre pour w_lo ; quantile : `erfc_inv` puis un pas de Newton sur un
résidu précis.

**Précision mesurée** (contre mpmath, `tests/test_random_precision.cpp`, références
engendrées par `tools/random/law_references.py`) : en double, au plus 0,5012 ulp pour
`pdf`, `cdf`, `ccdf`, `quantile` des deux lois (arrondi correct dans 99,9 % des cas ou
plus) ; `log_pdf` aussi sur les points mesurés, mais sa précision relative se dégrade
quand ses termes se compensent (résultat proche de 0). En float, au plus 0,5000 ulp.
Paramètres extrêmes compris (σ = 10⁻³⁰⁰ ou 10²⁰⁰, λ = 2^±1000).

Tables de la ziggourat : `tools/random/ziggurat.py` (mpmath) écrit
`detail/ziggurat_tables.hpp`, versionné.

### Étape C1 : géométrique, Poisson, binomiale

| Loi | Paramètres | Fonctions | Tirage (approché, bit-exact) |
|---|---|---|---|
| `geometric<I>` | `{.p}`, échecs avant le premier succès (k >= 0) | pmf, log_pmf, cdf, ccdf, quantile **exact** ; moyenne, variance, mode, médiane, asymétrie, aplatissement, entropie | inversion ⌊log U / log1p(−p)⌋ |
| `poisson<I>` | `{.mean}` | pmf, log_pmf ; moyenne, variance, mode, asymétrie, aplatissement | inversion séquentielle (μ < 10), PTRS (Hörmann 1993) |
| `binomial<I>` | `{.trials, .p}` | pmf, log_pmf ; moyenne, variance, mode, asymétrie, aplatissement | inversion (n min(p, 1−p) < 10), BTRS (Hörmann 1993), symétrie |

Convention de la géométrique : celle de `std::geometric_distribution` et de R (support
{0, 1, …}) ; le nombre d'essais est X + 1.

Masses de Poisson et binomiale : décomposition de Loader (2000), `detail/loader.hpp` —
l'erreur de Stirling δ et la déviance bd0 en double-double (série de bd0 près de
x = M en double-double aussi : elle pèse jusqu'à 1/15 d'un exposant qui peut valoir des
centaines). Répartitions et quantiles de Poisson et binomiale : avec la gamma et la
bêta incomplètes (phase 4 de wsl::special).

Précision mesurée : au plus 0,5009 ulp pour toutes les masses, répartitions et survies,
de μ = 10⁻⁵ à 10¹², n jusqu'à 10⁹, p jusqu'à 2⁻⁶⁰ ; quantile géométrique exact sur
1 034 points (hors quotients à 2⁻⁹⁰ près d'un entier).

Coût du tirage (même machine) : géométrique 18 ns (`std` 87), Poisson(4) 34 ns (53),
Poisson(500) 40 ns (71), binomiale(20 ; 0,2) 40 ns (33), binomiale(10⁴ ; 0,3) 49 ns
(165).

### Étape C2 : gamma et bêta

| Loi | Paramètres | Fonctions | Tirage (approché, bit-exact) |
|---|---|---|---|
| `gamma<T>` | `{.rate, .shape}` ou `{.scale, .shape}` (exacts, double ou `ratio`) ou `{.shape, .scale}` (le paramètre donné est gardé tel quel) | pdf, log_pdf ; moyenne, variance, mode, asymétrie, aplatissement | Marsaglia–Tsang (normale par la ziggourat) ; a < 1 : G(a + 1) U^{1/a} |
| `beta<T>` | `{.alpha, .beta}` (lus par `a()`, `b()`) | pdf, log_pdf ; moyenne, variance, mode, asymétrie, aplatissement | G₁/(G₁ + G₂), en logarithmes si α < 1 ou β < 1 |

Densités comme `dgamma` et `dbeta` de R : masse de Poisson généralisée pour la gamma,
binomiale généralisée pour la bêta (α, β > 2), forme directe sinon avec log B(α, β) en
double-double sans compensation (Stirling quand le plus grand paramètre dépasse 16 ;
lgamma(a + b) corrigé au premier ordre par ψ). Répartitions, quantiles et entropies :
phase 4 (gamma et bêta incomplètes, digamma).

Deux leçons, qui ont amélioré au passage les lois des étapes B et C1 (empreintes des
fonctions changées en conséquence) :
- l'erreur de Stirling δ(y) pour y < 16 se calcule par récurrence jusqu'à y + n >= 16,
  le produit y (y + 1) … tenu à l'échelle, et non par `lgamma_positive` de wsl::special,
  précise à 2⁻⁵⁸ seulement ; δ(s) d'un argument double-double prend sa correction du
  premier ordre −s_lo/(12 s²) ;
- l'exponentielle à argument double-double a son propre noyau (Taylor de degré 7,
  produit de tête exact) : 2⁻⁷⁰ près au lieu de 2⁻⁶⁰.

Précision mesurée : densités et log-densités gamma et bêta au plus 0,5000 ulp
(paramètres de 10⁻³⁰⁰ à 10¹², taux et échelles extrêmes) ; sur l'ensemble des lois,
le maximum est maintenant 0,5032 ulp (`cdf` normale, limitée par erfc de wsl::special,
à 2⁻⁶⁰ près).

Coût du tirage : gamma(2,5) 23 ns (`std` 74), gamma(0,3) 62 ns (62), bêta(2, 5) 45 ns.

### Phase 4 côté lois : répartitions, quantiles, entropies

Avec les fonctions incomplètes de wsl::special (phase 4) :

| Loi | Nouveau | Méthode |
|---|---|---|
| `poisson` | cdf, ccdf, quantile exact | Q(k + 1, μ), P(k + 1, μ) ; quantile : plus petit k tel que F(k) >= p, recherche depuis Cornish–Fisher, F(k) comparée à p en double-double du côté minoritaire (`detail/discrete.hpp`) |
| `binomial` | cdf, ccdf, quantile exact | 1 − I_p(k + 1, n − k), I_p(k + 1, n − k), p tel quel ; quantile de même |
| `gamma` | cdf, ccdf, quantile, entropie | P(a, t), Q(a, t), t = λx (ou x/θ) en double-double par la correction du premier ordre f(t) t_lo ; quantile : l'inverse en double-double divisé par λ (multiplié par θ) une fois — et, pour t < 2⁻⁹⁰⁰ (résultat près des sous-normaux), x = exp((log p + lgamma(1 + a))/a + log θ) directement |
| `beta` | cdf, ccdf, quantile, entropie | `ibeta`, `ibetac`, `ibeta_inv` |

Entropies réécrites avec l'erreur de Stirling δ et sa dérivée δ' (la forme usuelle
compense des termes en a log a) :
- gamma : ½ log a + ½ log 2π + ½ − 1/(2a) + δ(a) + (1 − a) δ'(a) − log λ ;
- bêta (n = α + β) : ½ log(αβ/n³) + ½ log 2π + ½(1 − 1/α − 1/β + 2/n) + δ(α) + δ(β) − δ(n)
  − (α − 1)δ'(α) − (β − 1)δ'(β) + (n − 2)δ'(n) (attention : (3/2) log n en double-double,
  le produit par 1,5 n'est pas exact).

Précision mesurée (références à 320 et 400 bits) : toutes ces fonctions sous 0,5 ulp,
arrondi correct à 100 % (double et float) ; entropies sous 0,48 ulp ; quantiles de
Poisson et binomiaux exacts sur 1 300 points.

Ce qui reste hors du catalogue en route : entropies de Poisson et binomiale (sans forme
close ; séries à venir), médianes de gamma, bêta, Poisson, binomiale (quantile(0,5)).

### Lois dérivées, premier lot : χ², Student, log-normale, Cauchy

| Loi | Paramètres | Méthode | Tirage (approché, bit-exact) |
|---|---|---|---|
| `chi_squared<T>` | `{.degrees_of_freedom}` | délégué à `gamma` (forme k/2, échelle 2, exactes) | gamma |
| `student_t<T>` | `{.degrees_of_freedom}` | queue I_x(ν/2, ½)/2, x et 1 − x en double-double, bêta incomplète à argument double-double (premier ordre) ; quantile par l'inverse de la bêta en (x, 1 − x) double-double, t² = ν(1 − x)/x ; densité sans compensation (lgamma((ν+1)/2) − lgamma(ν/2) par récurrence et Stirling) ; moments NaN ou ∞ là où ils n'existent pas | Z/√(V/ν) |
| `lognormal<T>` | `{.log_mean, .log_stddev}` (de log X) | z = (log x − μ)/σ, log x en double-double, puis la normale ; quantile e^{μ+σΦ⁻¹(p)}, un arrondi ; moments en double-double | e^Y, Y exacte ; `wsl::random::fast` : ziggourat |
| `cauchy<T>` | `{.location, .scale}` | queue atan(1/z)/π (z >= 1) ou ½ − atan(z)/π, noyau atan à argument double-double ; quantile x₀ − γ cot(πp) ; pas de moyenne ni de variance (membres absents) | x₀ ± γ cot(πW/2) |

Précision mesurée : tout sous 0,5 ulp, arrondi correct à 100 % (double et float).
Défauts trouvés en route, tous dans les queues extrêmes :
- quantile de Student pour p minuscule : x = I⁻¹(2p) sous les sous-normaux alors que
  t ~ √(ν/x) reste représentable (rendait ∞) ; x y est donné par log x =
  (log 2p + log a + log B(a, ½))/a, d'où t directement ;
- cot(πp) pour p < 2⁻⁹⁷⁰ (produit exact hors domaine) : (1/(πp))(1 − (πp)²/3) à
  l'échelle, sans débordement intermédiaire ;
- atan(v) ≈ v pour v < 2⁻²⁷ oubliait −v³/3 (0,55 ulp).

### Lois dérivées, second lot : Fisher, Weibull, Laplace, logistique, Pareto

| Loi | Paramètres | Méthode | Tirage |
|---|---|---|---|
| `fisher_f<T>` | `{.d1, .d2}` | u = d₁x/(d₁x + d₂) et 1 − u en double-double ; répartition I_u(d₁/2, d₂/2) (bêta incomplète à argument double-double) ; densité par le préfixe bêta logarithmique (`beta_log_prefix`, Loader si a, b > 2), exacte près de 0 pour d₁ = 2 (−(b + 1) log1p(2x/d₂)) ; quantile x = d₂u/(d₁(1 − u)), en logarithme quand u ou 1 − u tombe sous 2⁻⁹⁶⁰ ; entropie par celle de la bêta et ψ réécrite avec δ′ | (G₁/a)/(G₂/b), approché |
| `weibull<T>` | `{.scale = 1, .shape}` | t = (x/λ)ᵏ = e^{k log(x/λ)}, log du rapport par log1p près de 1 ; répartition −expm1(−t), survie e^{−t} ; variance λ²Γ(1 + ε)²(e^D − 1), D = lgamma(1 + 2ε) − 2 lgamma(1 + ε) en série pour ε = 1/k <= 1/8 (sans compensation pour k grand) | λE^{1/k}, approché |
| `laplace<T>` | `{.location, .scale}` | z en double-double, queue e^{−|z|}/2 ; quantile μ ∓ b log(2 min(p, 1 − p)), un arrondi | **exact** : μ ± bE, E exponentielle exacte |
| `logistic<T>` | `{.location, .scale}` | queue e^{−|z|}/(1 + e^{−|z|}) ; quantile μ + s logit(p), logit par log1p((2p − 1)/(1 − p)) près de ½ | μ + s logit(U), approché |
| `pareto<T>` | `{.scale, .shape}` | L = α log(x/x_m) en double-double ; survie e^{−L}, répartition −expm1(−L) ; quantile x_m e^{−log1p(−p)/α} ; moyenne et variance +∞ quand elles divergent | x_m e^{E/α}, approché |

Précision mesurée : tout sous 0,5 ulp (double et float), arrondi correct sur
99,7 à 100 % des points. Défauts trouvés en route :
- `exp_scaled` (wsl::special, détail) rendait une paire non normalisée (bas ~ 2⁻²²
  relativement au haut) : réinjectée dans une autre exponentielle (survie de Weibull),
  elle perdait 20 bits. Elle est désormais renormalisée (aucune empreinte ne change) ;
- log-densité de Fisher pour d₁ = 2 : f(0) = 1, log f ≈ 0 près de 0 se compensait
  (120 ulp) ;
- quantile de Fisher sous-normal : u sous-normal n'a plus tous ses bits (2 ulp).

## 10. Phase 5 : algèbre des lois, lot 1

Décision (2026-11-18) : **identités de familles, puis composition**. Une famille fermée
donne une loi de la famille, sous 0,5 ulp comme toute loi ; sinon la loi construite
compose les fonctions de ses lois de base, bit-exacte partout, avec une erreur mesurée
et bornée par le conditionnement.

### Quantile de survie

`cquantile(X, q)`, x tel que P(X > x) = q, pour les quatorze lois continues : pour
q >= ½, `quantile(X, 1 − q)` (1 − q exact) ; sinon une branche propre à chaque loi, q
exact (−log q/λ, μ − σΦ⁻¹(q), Q⁻¹(a, q), I⁻¹ complémentaire, x₀ + γ cot(πq)…). Tout sous
0,5 ulp, arrondi correct à 99,4–100 % (14 lois × double et float, références mpmath
écrites sur la queue). Sans lui, l'image par a < 0 et la troncature dans une queue
perdraient tout pour q minuscule.

### Image affine

| Famille | aX + b |
|---|---|
| normale, Cauchy, Laplace, logistique, uniforme | même famille (`affine_closure`), type conservé ; paramètres arrondis une fois par appel de `affine` |
| exponentielle, gamma, Weibull, Pareto, log-normale | `affine_image<L>`, mais pour a > 0, b = 0 toutes les fonctions sont celles de la loi de la famille (`scale_closure`) : mêmes bits |
| autres | `affine_image<L>` : z = (y − b)/a arrondi une fois, cdf/ccdf échangées si a < 0, quantile a·Q(p) + b ou a·cQ(p) + b en un arrondi, moments a·m + b, a²v, ±asymétrie, H + log a |

Les points de personnalisation `affine_closure(X, a, b)` et `scale_closure(X, c)` sont
des fonctions libres trouvées par ADL : une loi de l'utilisateur peut les fournir.

### Troncature

`truncate(P, {.lower = a, .upper = b})`, loi de X sachant a <= X <= b (bornes absentes : ∓∞). Masse, répartition et survie sont des
différences de F ou de S de X choisies point par point (deux termes <= ½ dès que
possible) ; les cibles des quantiles sont des sommes de positifs, F(a) + p·m ou
S(b) + (1 − p)·m, résolues par `quantile` ou `cquantile` de X : N | N >= 10 (masse
7,6·10⁻²⁴) garde sa précision. Tirage par rejet sur X si m >= ¼ (exact si celui de X
l'est), par inversion sinon.

### Mélange

`mixture(weighted(w₁, X₁), …)`, composantes de types quelconques, de même `value_type`
(réel ou entier). Densité, masse, répartition : Σ wᵢ fᵢ / W en double-double ;
variance Σ wᵢ (vᵢ + (mᵢ − m)²) / W, sans compensation ; asymétrie et kurtosis par les
moments centrés ; quantiles par dichotomie sur les flottants (ou les entiers) entre
les quantiles extrêmes des composantes. **Tirage exact** si chaque composante l'est :
les poids flottants sont des entiers à un 2^e commun près, la composante est tirée
uniformément sous leur somme (64 bits, ou grands entiers quand les poids s'étendent
de 2^-1000 à 2^1000).

### Erreur de composition mesurée

| Loi construite | pdf | cdf | ccdf | quantile |
|---|---|---|---|---|
| −2T + 5, T de Student (ν = 1,5 ; 4 ; 30) | 0,50 | 0,50 | 0,50 | 3,7 (quand 5 + 2Q se compense) |
| N(0, 1) tronquée (cinq intervalles, queues comprises) | — | 2,2 | 5,2 | 61 (près du zéro du quantile) |
| 0,3 N(−2, 1) + 0,7 gamma(3, 1) | 1,1 | 0,8 | — | 1,5 |

(ulp, double ; points loin des bornes de troncature.) Limites inhérentes à la
composition, documentées : près d'une borne de troncature, l'erreur suit le
conditionnement |x|/(x − a) ; près d'un zéro du résultat (quantile tronqué, image
affine), seule l'erreur absolue reste petite ; une valeur de F ou de S de X sous les
sous-normaux (N tronquée au-delà de |x| = 38) perd ses bits.

### Défauts trouvés en route

- `ccdf` de Cauchy nulle au-delà de |z| = 2⁶⁰⁰ (`detail::standardize` sature, juste pour
  les lois à queue exponentielle) : z gardé à l'échelle, sans saturation ; et 1/z par
  Dekker sortait de son domaine au-delà de 2⁹⁹⁵. Log-densités de Laplace et de la
  logistique de même (−|z|, et non −∞). Références étendues à |z| = 2¹⁰⁰⁰.
- Produits rendus sans `fmul` (variance de la normale et du χ², **tirage gamma** par
  l'échelle) : GCC les fusionnait avec l'addition de l'appelant.

## 11. Phase 5 : algèbre des lois, lot 2

### Sommes indépendantes

`convolution(P, Q)` est la loi de X + Y, X ~ P et Y ~ Q **indépendantes** (`convolution(P, P)` :
deux copies indépendantes). Différence : `convolution(P, affine(Q, {.scale = -1}))`. Les
lois n'ont pas d'opérateurs arithmétiques : `X + Y` est réservé aux variables aléatoires
(chantier à venir), où `X + X` vaut 2X. Identités par
`sum_closure(X, Y)` (ADL) :

| Somme | Résultat |
|---|---|
| normale + normale | `normal`, μ₁ + μ₂, √(σ₁² + σ₂²) par `hypot` (un arrondi chacun) |
| Cauchy + Cauchy | `cauchy`, x₁ + x₂, γ₁ + γ₂ |
| Poisson + Poisson, χ² + χ² | `poisson`, `chi_squared` |
| gamma + gamma (même intensité, paramétrisations mêlées comprises), exponentielle + exponentielle ou gamma | `convolution_law`, fonctions de la gamma de forme a + b |
| binomiale + binomiale (même p) | `convolution_law`, fonctions de la binomiale n₁ + n₂ |

L'identité conditionnelle est testée à la construction (`S.family()`) : intensités
**exactement** égales (θλ = 1 exactement entre échelle et intensité). Sans identité :
moyenne, variance, asymétrie et kurtosis par les cumulants, qui s'ajoutent ; lois
discrètes : masse Σⱼ P(X = j) P(Y = k − j), répartition et survie par des sommes finies
de termes positifs (la survie lointaine n'est jamais 1 − cdf), en double-double, en
O(longueur du support utile) ; quantile exact de cette répartition. Une somme continue
sans identité n'a pas de densité calculée (pas de convolution numérique) : quand
l'identité conditionnelle échoue, ses fonctions rendent NaN et ses quantiles lèvent
`std::domain_error` ; tirage et moments restent disponibles.

Tirage : x puis y, x + y (exact pour deux lois discrètes exactes ; un arrondi pour des
réels), quelle que soit l'identité (le nombre de valeurs consommées ne dépend pas des
paramètres).

### Loi jointe

`wsl::random::product_measure(P₁, …, Pₙ)` (type et CTAD), `value_type` = `std::tuple` : densité (toutes continues)
ou masse (toutes discrètes) produit en un arrondi, log-densité somme, répartition et
survie jointes produits, entropie somme ; tirage des composantes dans l'ordre, exact si
chacune l'est.

### Image par une fonction

`wsl::random::pushforward(X, f)` : tirage f(x) seul pour une fonction quelconque (à valeurs
arithmétiques). Une bijection monotone se déclare avec sa réciproque (et la dérivée de
celle-ci pour la densité) : `wsl::random::increasing(f, g[, g'])`, `wsl::random::decreasing(f, g[, g'])`.
La précision et la reproductibilité au bit près sont celles de f et g : avec
`wsl::special::exp` et `log`, exp(N) redonne les bits de la log-normale.

### Erreur mesurée

| Loi construite | pmf | cdf | ccdf |
|---|---|---|---|
| Poisson(μ) + binomiale(20, p), quatre cas, queues comprises | 1,4 | 0,84 | 0,95 |

(Mélange : quantile jusqu'à 1,5 ulp sur le nouveau tirage des points, loin du zéro
du quantile.)

## 12. Phase 6 : lois multivariées, lot 1

### Socle (`multivariate.hpp`)

Dimension N fixe ou `wsl::random::dynamic` (= `std::dynamic_extent`) : valeurs `std::array<T, N>`
ou `std::vector<T>`, la dimension dynamique donnée par les paramètres. Matrices :
`wsl::random::square_matrix<T, N>`, ligne par ligne, construite par accolades imbriquées
(`.covariance = {{1, 0.5}, {0.5, 2}}`, dimension vérifiée). Le support d'une loi
multivariée est le pavé qui contient ses valeurs (le concept de loi est inchangé).
Calcul : produits scalaires en double-double, produits par `fmul`, Cholesky et
substitution triangulaire bit-exactes.

| Loi | Paramètres | Fonctions | Tirage (approché, bit-exact) |
|---|---|---|---|
| `dirichlet<T, N>` | `{.alpha}` | densité en logarithme double-double (log B(α) = Σ lgamma − lgamma(α₀)), x dans le simplexe à K·2⁻⁵⁰ près ; moyenne, variances, covariance, mode, entropie, marginale bêta | Gᵢ/Σ G, gammas tirées dans l'ordre |
| `multinomial<I, N>` | `{.trials, .p}` (Σ p = 1 à 2⁻⁴⁰ près) | masse n!/Π kᵢ! Π pᵢ^kᵢ en logarithme double-double, un arrondi ; moyenne, covariance, marginale binomiale | binomiales conditionnelles, restes en double-double |
| `multivariate_normal<T, N>` | `{.mean, .covariance}` | Cholesky à la construction (`std::domain_error` sinon) ; log-densité −½\|L⁻¹(x − μ)\|² − Σ log Lᵢᵢ − (n/2) log 2π ; entropie, Mahalanobis, marginale normale | μ + L z, z normales exactes ; `wsl::random::fast` : ziggourat |
| `uniform_sphere<T, N>`, `uniform_ball<T, N>` | dimension (N, ou `{.dimension}`) | densité 1/aire (mesure de surface), 1/volume ; aire, volume, entropie | z/\|z\| ; boule : × U^{1/n} |

Vérifications : valeurs contre mpmath à 2 ulp (densités, masses, entropies),
moments empiriques, marginales par Kolmogorov–Smirnov, x₃ uniforme sur [−1, 1] pour la
sphère S² (Archimède), \|x\|ⁿ uniforme pour la boule.

### Lot 2

| Loi | Paramètres | Fonctions | Tirage (approché, bit-exact) |
|---|---|---|---|
| `multivariate_student_t<T, N>` | `{.location, .scale, .degrees_of_freedom}` | log-densité en double-double (log1p(q/ν)), moyenne, covariance ν/(ν − 2) Σ, entropie, marginale `affine_image<student_t>` | μ + L z √(ν/V) |
| `wishart<T, N>` | `{.scale, .degrees_of_freedom}` | log-densité par Cholesky de X, tr(V⁻¹X) = ‖L⁻¹L_X‖²_F ; moyenne, variances, mode, entropie (ψ_p), marginale gamma | Bartlett : (LA)(LA)ᵀ |
| `inverse_wishart<T, N>` | `{.scale, .degrees_of_freedom}` | log-densité, moyenne, mode, entropie | (LA⁻ᵀ)(LA⁻ᵀ)ᵀ |
| `haar_orthogonal<T, N>`, `haar_unitary<T, N>` | dimension | moyenne | Q de la QR d'une matrice gaussienne (Gram–Schmidt réorthogonalisé, diagonale de R positive : exactement Haar) |
| `goe<T, N>`, `gue<T, N>` | dimension | densité ∝ e^{−tr H²/2} (diagonale N(0, 1), hors diagonale N(0, ½)), log en double-double | triangle supérieur ligne par ligne |

Valeurs complexes : `std::complex<T>`, mais les produits sont faits ici avec `fmul`,
jamais par l'opérateur `*` de `std::complex` (contractable en fma dans l'unité de
compilation de l'utilisateur). Vérifications : densités et entropies contre mpmath
(2 ulp ; entropie de l'inverse-Wishart à 3 ulp, ψ en double), cas n = 1 (Student,
gamma), moyennes empiriques, marginales par Kolmogorov–Smirnov, orthogonalité et
unitarité à 2·10⁻¹⁵, colonne de Haar uniforme sur S², |u₀₀|² de loi bêta(1, n − 1),
déterminant ±1 équiprobable.

### von Mises et von Mises–Fisher

Avec `wsl::special::cyl_bessel_i` (docs/special.md §8) :

| Loi | Paramètres | Fonctions | Tirage (approché, bit-exact) |
|---|---|---|---|
| `von_mises<T>` | `{.location = 0, .concentration}`, support [μ − π, μ + π] | densité exp(−2κ sin²((x − μ)/2) − log 2π − log(e^{−κ}I₀(κ))), x − μ et le sinus en double-double ; longueur résultante A = I₁/I₀, variance circulaire 1 − A, entropie | Best et Fisher (1979) |
| `von_mises_fisher<T, N>` | `{.mean_direction, .concentration}` (μ normalisé) | densité exp(κ(μᵀx − 1) + log C_p + κ), μᵀx − 1 en double-double ; moyenne A_p μ, entropie ; κ = 0 : uniforme | Wood (1994), réflexion de Householder e₁ → μ |

Mesure contre mpmath : densités, log-densités, A, variance circulaire, entropies sous
0,5 ulp, arrondi correct sur 100 % des 2 900 points (κ de 0,01 à 10⁵, dimensions 3 à
120). Tirages : Kolmogorov–Smirnov contre la répartition intégrée (von Mises) et contre
la loi exacte de μᵀX pour p = 3, azimut uniforme, E[X] = A_p μ. Pas de répartition de
von Mises (pas de forme close ; série de Bessel possible plus tard).

## 13. Phase 7 : objets combinatoires, lot 1

Bibliothèque à part, **`wsl::permutation`** (`<wsl/permutation.hpp>`, C++23, dont
`wsl::random` dépend) : permutations de {0, …, n − 1}, `(s * t)(i) = s(t(i))`, inverse,
puissances, cycles canoniques, signature, inversions en O(n log n), rang
lexicographique. Le concept de loi n'exige plus de support (`has_support` le teste) :
un intervalle n'a pas de sens pour une permutation.

Tous les tirages sont **exacts** : l'objet tiré a exactement la probabilité annoncée.
Les probabilités qui ne sont pas des flottants (θ/(θ + i), 1/q) sont décidées par
`detail::choose` et `detail::fraction` : deux flottants sont des entiers à une même
puissance de 2 près, un entier uniforme sous leur somme décide (128 bits, ou grands
entiers).

| Loi | Paramètres | Masse | Tirage |
|---|---|---|---|
| `uniform_permutation` | `{.size}` | 1/n! | Fisher–Yates |
| `ewens` | `{.size, .theta}` | θ^{cycles}/θ(θ + 1)…(θ + n − 1) ; `mean_cycles()` | restaurant chinois : nouveau cycle avec probabilité θ/(θ + i) (exacte), sinon insertion après un élément uniforme |
| `mallows` | `{.size, .q}` (q > 0) | q^{inversions}/Π [i]_q, [i]_q = (1 − qⁱ)/(1 − q) en double-double ; `mean_inversions()` | code d'insertion, P(V_i = v) ∝ q^v par épreuves de Bernoulli exactes : géométrique tronquée si i(1 − q) >= 1, uniforme acceptée avec probabilité q^v sinon ; q > 1 par symétrie (rapport 1/q exact) |
| `uniform_subset` | `{.size}` | 2⁻ⁿ | un bit par élément |
| `bernoulli_subset` | `{.size, .p}` | p^k (1 − p)^{n−k} | épreuves exactes |
| `uniform_combination` | `{.size, .count}` | 1/C(n, k) | Floyd (k entiers uniformes) |
| `uniform_word` | `{.length, .alphabet}` | k^{−ℓ} | lettres uniformes |
| `uniform_integer_partition` | `{.total}` | 1/p(n) ; `count()` | Nijenhuis–Wilf, p(·) en grands entiers (récurrence pentagonale) |

Valeurs : `wsl::permutation` ; `std::vector<std::size_t>` (éléments croissants,
lettres, parts décroissantes). Vérifications : χ² contre la masse exacte sur tous les
objets d'une petite taille (S₄, 2⁴ sous-ensembles, C(6, 3), 3³ mots, 22 partitions de 8),
masses sommant à 1 sur S₆, espérances (cycles, inversions) exactes, p(100) et p(1000).

### Lot 2 : graphes, arbres, partitions d'ensemble, compositions

Valeurs : `wsl::random::edge_list` (arêtes (u, v), u < v, triées) ; `wsl::random::binary_tree` (enfants
{gauche, droit} numérotés en préordre, racine 0, `binary_tree::none`) ; partition
d'ensemble : mot à croissance restreinte ; composition : liste des parts.

| Loi | Paramètres | Masse | Tirage exact |
|---|---|---|---|
| `erdos_renyi` | `{.vertices, .p}` | p^m (1 − p)^{N−m} | une épreuve par paire (O(n²)) |
| `uniform_graph` | `{.vertices, .edges}` | 1/C(N, m) | Floyd sur les indices des paires |
| `barabasi_albert` | `{.vertices, .edges}` | — | K_{m+1}, puis m cibles distinctes ∝ degré (extrémité uniforme d'arête, doublons rejetés) |
| `random_regular_graph` | `{.vertices, .degree}` | — (uniforme) | modèle de configuration rejeté jusqu'à un graphe simple (e^{(d² − 1)/4} essais : petits d) |
| `uniform_tree` | `{.vertices}` | n^{−(n−2)} | suite de Prüfer, décodage O(n) |
| `uniform_binary_tree` | `{.vertices}` | 1/C_n | Rémy (arbres complets), feuilles retirées |
| `uniform_set_partition` | `{.size}` | 1/B_n ; `count()` | élément par élément, poids T(r, k) des complétions (grands entiers, table O(n²)) |
| `uniform_composition` | `{.total, .parts = 0}` ou `{.n, .parts}` | 2^{1−n}, 1/C(n − 1, k − 1) | points de coupure uniformes |

Vérifications : χ² contre la loi exacte sur tous les objets d'une petite taille (64
graphes à 4 sommets, C(10, 3), les 70 graphes 2-réguliers à 6 sommets, les 6 issues de
Barabási–Albert (4, 1), 16 arbres de Cayley, 14 arbres binaires, 52 partitions, 32 et 10
compositions) ; grandes tailles : degrés, arbres valides, B₂₀.

## 14. Phase 8 : processus ✅ (exemples : interface réalisée)

### Principe

Un processus est, comme une loi, un **objet mathématique** : il décrit la loi d'une
trajectoire (X_t). Trois façons de s'en servir, du plus fini au plus paresseux :

| Usage | Écriture | Valeur |
|---|---|---|
| Loi marginale | `P.at(t)` | une loi de la bibliothèque, quand elle est connue (`normal`, `poisson`, `categorical`…) |
| Restriction finie | `P.steps(n)`, `P.on(T)`, `P.on_grid(ts)` | une **loi** ordinaire sur un objet fini (`sample`, `pmf`/`pdf` de la trajectoire quand ils existent) |
| Trajectoire paresseuse | `trajectory(P, g)` (`events(P, g)` ◻︎ non fait : `trajectory` rend déjà les instants des processus ponctuels) | une vue `std::ranges` (input_range, infinie ou jusqu'à absorption), tirée au fur et à mesure |

**Contrat de reproductibilité** : la k-ième valeur d'une trajectoire ne dépend que de
l'état initial du moteur et de k ; tronquer (`views::take(k)`, `take_while`) ne change
pas les valeurs gardées, et le moteur a alors consommé exactement ce qu'il faut pour
elles. La vue tient le moteur **par référence** (comme `uniform_below(n, g)`) : elle
l'avance, et il doit lui survivre. ✅ (2026-11-18)

Ordre ✅ : temps discret d'abord (catégorielle exacte, marches, chaînes de Markov, urnes
de Pólya, restaurant chinois), puis temps continu.

### Temps discret

```cpp
namespace wr = wsl::random;
wr::random_walk W(wr::normal({.mean = 0, .stddev = 1}), {.start = 0});
W.variance(10);                             // 10 : moyenne et variance de la position au pas n
for (double x : wr::trajectory(W, g) | std::views::take(100)) …

wr::markov_chain C({.transition = {{0.9, 0.1}, {0.5, 0.5}}, .start = 0});
C.at(5);                                    // wr::categorical, π₀ P⁵
C.stationary();                             // wr::categorical
auto path = C.steps(100)(g);                // 101 états
C.steps(100).pmf(path);                     // produit des transitions
```

Brique nécessaire, en tête de phase : `wsl::random::categorical({.p = {…}})` (tirer i avec
probabilité pᵢ), **exacte** (poids entiers à une puissance de 2 près, table d'alias de
Walker–Vose pour l'O(1)).

### Temps continu

```cpp
wr::poisson_process N({.rate = 2.0});
N.at(3.0);                                  // poisson({.mean = 6})
for (double t : wr::trajectory(N, g) | std::views::take_while([](double t) { return t < 10; })) …
auto times = N.on(10.0)(g);                 // instants dans [0, 10]

wr::inhomogeneous_poisson_process M(f, 5.0);                          // intensité f <= 5 (amincissement)
wr::compound_poisson S(2.0, wr::exponential({.rate = 1.0}));           // taux, loi des sauts

wr::brownian_motion B({.drift = 0, .volatility = 1});
B.at(2.0);                                  // normal(0, √2)
auto x = B.on_grid({0.5, 1.0, 2.0})(g);     // accroissements N(μΔt, σ²Δt) exacts, cumulés
for (auto [t, v] : wr::trajectory(B, g, {.step = 0.01})) …
wr::brownian_bridge D({.from = {0, 0}, .to = {1, 0}});                 // grilles seulement (regard en avant)
wr::continuous_markov_chain Q({.generator = {{-1, 1}, {2, -2}}, .start = 0});   // instants et états des sauts
// aussi : wr::geometric_brownian_motion, wr::ornstein_uhlenbeck (transitions exactes)
```

### Modèles bayésiens non paramétriques, urnes

```cpp
wr::chinese_restaurant R({.theta = 1.5});
for (std::size_t table : wr::trajectory(R, g) | std::views::take(100)) …   // table de chaque client
R.at(100);                                  // loi de la partition (formule d'Ewens)
wr::polya_urn U({.balls = {1, 2}, .added = 1});                        // couleurs tirées ; U.at(n) : bêta-binomiale…
// processus de Dirichlet ◻︎ (bâton brisé) : non fait ; le restaurant en donne la partition
```

### Exactitude ✅ (2026-11-18 : option B)

Les choix discrets sont exacts (chaînes de Markov, restaurant, urnes : poids flottants
décidés exactement). En temps continu et pour les marches à pas réels, **exact point par
point par défaut** : chaque position est l'arrondi correct d'une trajectoire tirée
exactement (somme de réels paresseux, affinés à la demande), la loi étant celle des
paramètres arrondis une fois (√Δt, e^{−θΔt}…) ; `wsl::random::fast` : accroissements approchés
(ziggourat), sommes en double-double. Brique nouvelle : réels paresseux additionnables
et multipliables par un double (lot 2).

### Concepts

`stochastic_process<P>` : `state_type`, `trajectory(P, g)` ; `has_marginal<P>` :
`P.at(t)` ; les restrictions finies sont des `probability_law`.

### Lot 1 fait : temps discret

| Objet | Paramètres | Fonctions | Tirage exact |
|---|---|---|---|
| `categorical<I>` | `{.p}` (poids normalisés exactement) | masse, répartition, quantile exact, moments, entropie | table d'alias de Walker–Vose sur des entiers exacts, O(1) (grands entiers sinon) |
| `rademacher` | `{.p}` | ±1 | Bernoulli exacte |
| `random_walk<L>` | `(pas, {.start})`, pas entiers | `mean(n)`, `variance(n)`, `steps(n)` (masse de la trajectoire) | sommes entières |
| `markov_chain` | `{.transition, .start}` ou `{.transition, .initial}` | `at(n)` (π₀ Pⁿ), `stationary()`, `probability(i, j)`, `steps(n)` | lignes catégorielles |
| `polya_urn` | `{.balls, .added}` (réels >= 0) | `at(n)` : `dirichlet_multinomial` ; `steps(n)` | poids entiers exacts (64 bits, grands entiers au-delà) |
| `chinese_restaurant` | `{.theta, .discount}` (Pitman–Yor) | `at(n)` : `ewens_partition` (EPPF, `mean_blocks()`) | d = 0 : θ/(θ + n) exact puis un client uniforme ; d > 0 : poids entiers tenus à jour |

Vue paresseuse générique (`process.hpp`) : un processus fournit un marcheur
(`first(g)`, `next(g)`) ; une valeur n'est tirée qu'à sa lecture, de sorte que
`views::take(k)` ne consomme rien au-delà des k valeurs (vérifié : mêmes valeurs et même
état du moteur que `P.steps(k − 1)(g)`). Les marches à pas réels viennent au lot 2,
avec les réels paresseux additionnables (exactitude point par point).

Vérifications : χ² contre la loi exacte (catégorielle, position d'une marche contre la
binomiale, X₅ d'une chaîne contre π₀P⁵, les 81 trajectoires de longueur 4, couleurs
d'une urne contre la Dirichlet-multinomiale, partitions du restaurant contre l'EPPF pour
quatre (θ, d), poids non entiers compris) ; π P = π ; masses de trajectoires sommant à 1 ;
échangeabilité de l'urne.

### Lot 2 fait : temps continu, exact point par point (option B)

Brique : `detail::lazy_sum`, somme de termes c V / q (c dyadique exact, V réel
paresseux, q diviseur flottant), arrondie correctement. Voie rapide en double-double avec
une borne rigoureuse séparée en largeurs (réductibles en raffinant) et erreurs d'arrondi
(fixes) ; si les largeurs dominent, tous les termes gagnent un mot (amorti) ; sinon, voie
exacte en grands entiers (bornes rationnelles, dénominateur commun des diviseurs). Les
verdicts de la voie rapide décident du nombre de mots tirés : tous leurs calculs passent
par `fmul` (vérifié par le scan de contraction, qui tire désormais des trajectoires).
Contrôle : 1 500 sommes aléatoires (compensations, poids α, β, diviseurs) et 300 valeurs
visées au milieu de deux doubles, toutes confirmées par raffinement de deux mots de plus.
Coût mesuré : environ 2 µs par point de brownien, constant jusqu'à 10⁵ points (voie rapide
`wsl::random::fast` : quelques dizaines de ns) ; mémoire linéaire en la longueur.

| Processus | Paramètres | Usage | Construction exacte |
|---|---|---|---|
| `brownian_motion` | `{.drift, .volatility, .start}` | `at(t)`, `on_grid(ts)`, `trajectory(B, g, {.step})`, `wsl::random::fast` | x₀ + μt (exact) + σ Σ √Δt Z, √Δt arrondi une fois |
| `brownian_bridge` | `{.from, .to, .volatility}` | `at(t)`, `on_grid(ts)` | a + r(b − a) + σ(W(t) − r W(T)), r arrondi ; extrémités exactes |
| `ornstein_uhlenbeck` | `{.theta, .mean, .volatility, .start}` | `at(t)`, `stationary()`, `on_grid`, `trajectory(…, {.step})` | changement de temps μ + a(x₀ − μ) + r W(u), a, r, u arrondis par instant ; θt <= 350 |
| `geometric_brownian_motion` | `{.drift, .volatility, .start}` | `at(t)` (log-normale), `on_grid` | logarithme exact, exponentielle approchée (wsl::special::exp) |
| `poisson_process` | `{.rate}` | `at(t)`, `on(T)`, `trajectory` (instants), `wsl::random::fast` | Tₖ = Σ Eⱼ/λ |
| `inhomogeneous_poisson_process` | `(f, borne)` | `on(T)`, `trajectory` | amincissement, f(t)/Λ exact, f à l'instant arrondi (approché) |
| `compound_poisson` | `(taux, loi des sauts)` | `mean(t)`, `variance(t)`, `trajectory` (instant, valeur) | sauts paresseux (normale, exponentielle, Laplace, uniforme) ou arrondis sommés exactement |
| `continuous_markov_chain` | `{.generator, .start}` | `at(t)` (π₀ e^{Qt}), `trajectory` (instant, état) | séjours E/q_x exacts, sauts catégoriels ; +∞ après absorption |
| `random_walk` (pas réels) | `(pas, {.start})` | `trajectory`, `steps(n)` | somme paresseuse exacte |

Vérifications : Kolmogorov–Smirnov contre les lois exactes (B(1), accroissements, pont au
milieu, Ornstein–Uhlenbeck à deux instants et stationnaire, brownien géométrique, marches
normales et exponentielles contre normale et gamma, premier instant de Poisson), χ² des
comptages de Poisson et de l'état d'une chaîne en temps continu contre π₀ e^{Qt}, moyennes
du Poisson inhomogène (∫ f) et composé ; contrat de préfixe (trajectoire = grille, même
consommation). Défauts trouvés en route : boucle sans fin près de 0 (les erreurs d'arrondi
ne diminuent pas en raffinant : désormais séparées des largeurs), et poids qui s'annulent
(pont en T : largeurs pondérées par α + β et β séparément).

## 15. Incohérences notées pour la revue des interfaces

Liste relevée avant la revue des interfaces (relevé complet : [interfaces.md](interfaces.md)).

> **État après la revue** (décisions de interfaces.md §6) : réglés — noms des paramètres
> (`.stddev`, `.log_mean`/`.log_stddev`, `.location`/`.concentration`, `.lower`/`.upper`,
> `.total`, `.vertices`), types des moments (3a), erreurs et constantes indéfinies (3b,
> 3c), noms des constructions (`product_measure`, `convolution`/`convolution_law`),
> défauts de position et d'échelle (Gumbel, arcsinus…). Inchangés : types des valeurs
> combinatoires, constantes propres à certaines lois, support facultatif, normalisation des
> probabilités, conventions des processus, forme entière d'`erlang`, `as_gamma()`.

- types des valeurs : `wsl::permutation`, `wsl::random::edge_list`, `wsl::random::binary_tree` (types
  propres) mais `std::vector<std::size_t>` pour sous-ensembles, mots, partitions et
  compositions ; `wsl::random::square_matrix` pour les matrices, `std::array`/`std::vector` pour
  les vecteurs ;
- noms des paramètres de taille : `.size`, `.n`, `.vertices`, `.nodes`, `.length`,
  `.dimension` ; `.mu`/`.sigma` (log-normale, von Mises) contre `.mean`/`.sd` (normale) ;
- types des moments : `T` pour les lois continues, `double` pour les discrètes, tuples
  pour `independent`, vecteurs (variances) et matrices (covariances) pour les lois
  multivariées ;
- constantes propres à certaines lois hors du vocabulaire commun : `mean_cycles()`,
  `mean_inversions()`, `mean_size()`, `count()`, `mean_resultant_length()` ;
- noms des constructions : `wsl::random::affine` (fonction) et `affine_image` (type),
  `independent` (type, CTAD) et `independent_sum`, `truncate`/`truncated`,
  `pushforward`/`image`/`monotone_image` ;
- erreurs : `std::domain_error` à la construction et pour `quantile` hors de [0, 1],
  mais NaN (sans exception) pour une somme continue sans identité ;
- le support, facultatif depuis la phase 7 (`has_support`) ;
- duplication interne : poids entiers d'un mélange (`mixture`) et `detail::choose`.
- probabilités : `categorical` normalise ses poids exactement, `multinomial` exige une
  somme égale à 1 (à 2⁻⁴⁰ près) ; `polya_urn::at` divise les boules par `added`
  (paramètres arrondis) ;
- processus : pas régulier par `trajectory(P, g, {.step})` pour le brownien et
  Ornstein–Uhlenbeck, mais `trajectory(P, g)` ailleurs ; valeurs (instant, valeur) ou
  instants seuls ; constructeurs positionnels (`compound_poisson(taux, loi)`,
  `inhomogeneous_poisson_process(f, borne)`) ; `at(n)` rend des lois de types variés ; `steps(n)` est une loi de
  trajectoire, sauf pour le restaurant (`steps(n)` = `at(n + 1)`, partition) ;
  `random_walk` se construit par position (`random_walk(pas, {.start})`), les autres par
  paramètres nommés.
- lois des compléments : `erlang` prend une forme entière (`std::int64_t`) quand la gamma
  prend `T` ; `arcsine` a des paramètres par défaut ([0, 1]), pas `gumbel` ; la position
  s'appelle `.location` (Gumbel, Fréchet, Laplace) mais `.low`/`.high` pour les lois à
  support borné ; `as_gamma()` n'a pas d'équivalent pour les autres lois spéciales d'une
  famille (χ², exponentielle).

## 16. Compléments au catalogue, lot 1

| Loi | Paramètres | Méthode | Tirage |
|---|---|---|---|
| `gumbel<T>` | `{.location, .scale}` | y = (x − μ)/β en double-double, t = e^{−y} ; répartition e^{−t}, survie −expm1(−t) (t à l'échelle quand il est minuscule) ; quantile μ − β log(−log p), cquantile μ − β log(−log1p(−q)) | μ − β log E, approché |
| `frechet<T>` | `{.location = 0, .scale = 1, .shape}` | même noyau, y = α log((x − m)/s) (log1p près de 1) ; quantile m + s e^{y/α} ; moyenne et variance par Γ(1 − 1/α) et Γ(1 − 2/α), sans compensation (ε = −1/α) | m + sE^{−1/α}, approché |
| `rayleigh<T>` | `{.scale}` | h = z²/2 à l'échelle ; répartition −expm1(−h), survie e^{−h} ; quantile σ√(−2 log1p(−p)) | σ√(2E), approché |
| `triangular<T>` | `{.lower, .mode, .upper}` | différences exactes (two_sum), (x − a)²/((b − a)(c − a)) d'un côté du mode, complément de l'autre ; quantile a + √(p(b − a)(c − a)) ou b − √((1 − p)(b − a)(b − c)), 1 − p en double-double | inversion, approché |
| `arcsine<T>` | `{.lower = 0, .upper = 1}` | u et 1 − u en double-double ; (2/π) atan √(u/(1 − u)) du côté le plus proche, complément sinon ; quantile a + (b − a) sin²(πp/2) depuis la borne la plus proche | inversion, approché |
| `erlang<T>` | `{.rate, .shape}` ou `{.scale, .shape}` (shape entier) entier, `{.rate}` | fonctions de la gamma (`as_gamma()`) | **exact** : somme paresseuse de k exponentielles exactes, O(k) |

Moments partagés avec Weibull (`detail/gamma_moments.hpp`) : lgamma(1 + ε) en série pour
|ε| <= ¼ et variance c²[Γ(1 + 2ε) − Γ(1 + ε)²] en série pour |ε| <= ⅛, ε positif
(Weibull, 1/k) ou négatif (Fréchet, −1/α) ; Weibull garde ses bits. Précision mesurée :
tout sous 0,5 ulp en double et float, arrondi correct sur 100 % des points
(≈ 9 000 points, références mpmath ; constantes double-double nouvelles vérifiées par
`law_references.py`).

## 17. Compléments au catalogue, lot 2 : inverse-gaussienne

`inverse_gaussian<T>({.mean = μ, .shape = λ})`. Avec c = √(λ/x), a = c(x − μ)/μ, δ = 2c
et le rapport de Mills R(s) = e^{s²/2} Φ(−s), les deux termes de la répartition usuelle
(Φ(a) + e^{2λ/μ} Φ(−b)) partagent le facteur e^{−a²/2} :

    P(X > x) = e^{−a²/2} (R(a) − R(a + δ)),   P(X <= x) = e^{−a²/2} (R(|a|) + R(a + δ))  (a < 0),

sans e^{2λ/μ} (qui déborderait) ni compensation, sauf dans R(a) − R(a + δ) quand δ est
petit devant 1 + a⁺ : queue droite (la perte vaut environ (x − μ)/(2μ)) ou λ/x petit.
On y somme alors −Σ r_k(a) δᵏ, coefficients de Taylor de R en a, termes décroissants d'un
facteur 4 au moins ; au-dessus du seuil δ = (1 + a⁺)/4, la différence directe perd au plus
un facteur 6.

R et ses coefficients, à 2⁻⁸⁵ près (le noyau `erfc_factor` de wsl::special, à 2⁻⁶¹,
ne suffisait pas sous une compensation de 6) :

| s | Méthode |
|---|---|
| < 4,5 | ½ Σ (s²/2)ⁿ/n! − (1/√(2π)) Σ s^{2n+1}/(2n+1)!!, termes positifs (perte de 2²² au plus) ; coefficients par la récurrence (k + 1) r_{k+1} = a r_k + r_{k−1}, r_1 = a r_0 − 1/√(2π), avancée |
| [4,5, 2048) | récurrence rétrograde de Miller depuis k = 90 (solution minimale), renormalisée par 2⁻³⁰⁰ pour rester dans le domaine de Dekker, normalisée par a r_0 − r_1 = 1/√(2π) |
| >= 2048 | (1/(√(2π) s)) (1 − v + 3v² − … − 945v⁵), v = 1/s² |

Quantiles : Newton en log x sur log F (ou log S), encadré (bissection si le pas sort de
l'encadrement), 3 à 13 évaluations ; puis deux pas de Newton à résidu précis depuis le
flottant trouvé. Moments : variance μ³/λ, mode μ/(√(1 + t²) + t) (t = 3μ/(2λ), sans
compensation), asymétrie 3√(μ/λ), aplatissement 15μ/λ ; médiane par le quantile. Pas
d'entropie (elle demande E₁(2λ/μ)). Tirage approché, bit-exact : Michael, Schucany et
Haas (normale exacte, Bernoulli exacte ; la petite racine écrite μ²/x₂, sans compensation).

Coût : répartition et survie 2 à 5 µs, quantile 15 à 60 µs.

### Répartition et quantiles de von Mises

F(x) sur [μ − π, μ + π], θ = x − μ en double-double : la part la plus petite,
G(|θ|) = ∫_{|θ|}^π e^{κ(cos t − 1)} dt/(2π e^{−κ} I₀(κ)), l'autre par complément. La série
de Fourier ½ + θ/2π + Σ I_n(κ) sin(nθ)/(nπ I₀) se compense dans les queues d'un facteur
e^{2κ} (2⁵⁸ dès κ = 20) : l'intégrale est calculée directement, par morceaux [t₀, t₀ + H] et
Gauss–Legendre à 12 nœuds (nœuds et poids en double-double, calculés à la compilation par
Newton sur P₁₂) appliqué à e^{κ(cos t − cos t₀)}, le facteur e^{κ(cos t₀ − 1)} à l'échelle.
H <= ½, 2/(κ sin t₀), √(½/(κ |cos t₀|)) : l'erreur de la règle reste sous 2⁻¹⁰⁰ ; la
précision est celle de exp_scaled (2⁻⁷⁰). L'intégrande décroît sur [0, π] : on s'arrête
dès que e^{κ(cos t₀ − 1)}(π − t₀) devient négligeable. Un premier essai par séries de Taylor
de l'exponentielle (récurrence d_n = (1/n) Σ k c_k d_{n−k}) donnait les mêmes bits, 14 fois
plus lentement.

Quantiles : Newton sur log G pour t <= ¼ ; au-delà, sur log D, D(φ) = ∫_0^φ = ½ − t (exact),
sans quoi φ petit (t près de ½) n'était connu qu'à 2⁻⁷¹/(φ g) près ; départ au quantile
normal d'écart type 1/√κ, encadrement dans [0, π], deux pas à résidu précis. Quand la
solution tombe entre fl(π) et π (t < G(fl(π)), κ petit et t minuscule), φ = π − t/g(π).
Près de π (κ petit, t minuscule), π − φ ~ 10⁻¹⁴ : le pas de Newton tombe sous l'espacement
des doubles et n'était pas appliqué (jusqu'à 1,3 ulp) ; il l'est désormais en double-double.
Précision mesurée : répartition, survie et quantiles sous 0,5 ulp (double et float), arrondi
correct sur 100 % des points. Coût : répartition 1 à 15 µs, quantile 20 à 250 µs.

Références (mpmath) : mp.quad en t ne convergeait pas pour κ grand (erreur relative 2·10⁻¹¹
pour κ = 1000) ; en w = 2κ(sin²(t/2) − sin²(φ/2)), avec u = √(w(π) − w) sur la seconde
moitié (singularité en t = π levée) et des coupures géométriques depuis 2κ sin²(φ/2), elles
concordent avec la série de Fourier à haute précision.

## 18. Lois discrètes, lot 3 : binomiale négative, hypergéométrique

| Loi | Paramètres | Méthode | Tirage |
|---|---|---|---|
| `negative_binomial<I>` | `{.successes = r (réel > 0), .p}` | masse r/(k + r) b(r; k + r, p) par Loader (n = k + r réel en double-double) ; répartition I_p(r, k + 1), survie par complément de la bêta incomplète ; quantile exact | gamma–Poisson, approché |
| `hypergeometric<I>` | `{.population = N, .successes = K, .draws = n}` (N <= 2⁵³) | masse b(k; K, π) b(n − k; N − K, π)/b(n; N, π), π = n/N en double-double (dhyper de R) ; répartition : la queue la plus petite, sommée par le rapport de masses consécutives (O(σ) termes) ; quantile exact | inversion depuis le mode, approché |

Précision mesurée : masses, répartitions et moments sous 0,5 ulp (arrondi correct sur 100 %
des points) ; quantiles exacts sur 1096 points. Générateur : la recherche des quantiles
discrets de référence avançait pas à pas depuis une estimation normale (des milliers de bêtas
incomplètes à 400 bits pour σ ~ 10⁴) ; recherche exponentielle puis dichotomique.

## 19. Lois discrètes, lot 4 : bêta-binomiale, logarithmique, gaussienne discrète

| Loi | Paramètres | Méthode | Tirage |
|---|---|---|---|
| `beta_binomial<I>` | `{.trials, .alpha, .beta}` | masse par huit log Γ en double-double (les −y de Stirling s'annulent ; erreur ~ 2⁻¹⁰⁶ n log n, 0,5 ulp mesuré jusqu'à n = 10⁶) ; répartition : queue du côté de k sommée par le rapport de masses (arrêt anticipé si α, β >= 1, sinon O(n)) ; quantile exact | bêta puis binomiale, approché |
| `logarithmic<I>` | `{.p}`, 0 < p < 1 | masse e^{k log p}/(k L), L = −log1p(−p) ; répartition : tête (k <= 64) ou queue p^{k+1} Σ p^m/(k + 1 + m), O(1/(1 − p)) ; variance par L − p en série pour p petit ; quantile exact par sommes partielles | Kemp (LK), approché |
| `discrete_gaussian<I>` | `{.location = 0, .scale = 1}`, μ réel | Z, moyenne et variance par sommation de Poisson (σ >= 1, deux termes) ou somme directe ; répartition : queue sommée depuis k, O(σ) ; quantile exact | rejet depuis Laplace discrète (Canonne, Kamath, Steinke), approché |

Précision mesurée : tout sous 0,5 ulp (arrondi correct sur 100 % des points), quantiles
exacts sur 1299 points. Le tirage de Kemp avait sa dernière branche inversée (V >= q rend 1),
vu par le test du χ². Références : sommes cumulées gardées par processus (logarithmique,
bêta-binomiale, gaussienne discrète), sans quoi la recherche des quantiles les recalculait à
chaque essai ; la gaussienne discrète tabulée sur μ ± (40σ + 40) (les quantiles de
probabilités sous-normales tombent jusque vers 39σ).

## 20. Normale asymétrique

`skew_normal<T>({.location = ξ, .scale = ω, .shape = α})`, densité (2/ω) φ(z) Φ(αz). La forme
Φ(z) − 2T(z, α) (T d'Owen) se compense dans la queue courte (α = 1 : F = Φ² calculée
Φ − Φ(1 − Φ)) : la plus petite queue, du côté de z par rapport à la moyenne, est intégrée
directement, par morceaux et Gauss–Legendre à 12 nœuds, de largeur <= 1, 4/|(log f)′|,
1/√|(log f)″| et max(1/|α|, |t|/2) — Φ(αt) bascule sur une échelle 1/|α| autour de 0, qu'un
morceau réglé sur les seules dérivées locales enjambait (erreurs de 10¹⁴ ulp pour α = 50).
Densité log-concave : la queue restante est bornée par f/|(log f)′| au-delà du mode. Quantiles
par Newton sur le logarithme de la part visée. Coût : répartition ~ 50 µs, quantile ~ 300 µs.
Précision : densité, répartition, survie, quantiles et moments sous 0,5 ulp (arrondi correct
sur 100 % des points). Défaut trouvé en route : w/√2 arrondi sans son erreur (Φ(αz) à 10⁻¹³).

Références mpmath : mp.quad juge sa convergence en erreur absolue (~ 2^-prec) ; sur un
intégrande de 10⁻⁵² elle se trompait de 10⁻⁹ relativement. L'intégrande est désormais
normalisé par sa valeur en z ; contrôles par Φ (α = 0), Φ² (α = 1) et T d'Owen.

## 21. Lois stables

`stable<T>({.location = δ, .scale = γ, .alpha, .beta})`, paramétrisation S1 (celle de
Chambers–Mallows–Stuck et la plus courante) ; forme réduite Z, X = γZ + δ (α ≠ 1),
γZ + (2/π)βγ log γ + δ (α = 1).

Densité et répartition : intégrales de Nolan (1997) dans la variable φ = θ + θ₀ (sans
compensation près de −θ₀), noyau log g en double-double. Avec A = (1/π) ∫ e^{−g} et
B = (1/π) ∫ (1 − e^{−g}) (expm1), A + B = (π/2 + θ₀)/π et c = (π/2 − θ₀)/π, la répartition et
la survie sont toujours c + A, c + B, A ou B : sommes de termes positifs, y compris dans les
queues (où la forme usuelle 1 − (1/π) ∫ e^{−g} se compense). Une seule intégrale par appel.
Intégration : coupure au pic (g = 1, g monotone en φ), morceaux géométriques de part et
d'autre (premier : 1/8 de la largeur du pic, 1/|∂ log g/∂φ|, puis ×4), Gauss–Legendre
adaptatif sur chacun (tolérance 2⁻⁸² de la masse déjà intégrée), arrêt quand l'intégrande
décroît et devient négligeable. Formes closes : α = 2 (normale de variance 2γ²), α = 1 et β = 0
(Cauchy), y = 0 (Γ(1 + 1/α) cos θ₀ cos(αθ₀)^{1/α}/π ; 0 au bord du support pour α < 1,
|β| = 1). Coût : densité ou répartition ~ 0,5 ms, quantile ~ 6 ms.

Précision mesurée (9 paramétrages, α de 0,3 à 2, β de −1 à 1, y jusqu'à 4096) : densité,
répartition et survie sous 0,5 ulp, arrondi correct sur 100 % des points. Quantiles sans
référence mpmath (Newton dans les queues lourdes y coûtait 5 à 20 s par point) : vérifiés
par encadrement F(x⁻) <= p <= F(x⁺) à un ulp. Défauts trouvés en route :
- des nœuds grossiers manquaient la masse concentrée près du pic dans les queues (y = 100 :
  facteur 3 sur la densité) : morceaux géométriques depuis le pic ;
- π/2 + βθ formé avec (1 − β)π/2 arrondi, α − 1 inexact pour α < ½, 1 + 1/α arrondi dans
  Γ(1 + 1/α) : jusqu'à 5,6 ulp ; en double-double désormais.

Références : mêmes intégrales en mpmath (coupure au pic, morceaux géométriques, intégrande
normalisé), contrôlées par l'inversion de Fourier de la fonction caractéristique et les cas
Lévy, Cauchy et normal. Tirage : Chambers, Mallows et Stuck (forme de Weron), approché,
bit-exact.

## 22. Skellam, Zipf, Rice

| Loi | Paramètres | Calcul | Tirage |
|---|---|---|---|
| `skellam<I>` | `{.mean1, .mean2}` (> 0) | masse e^{−(μ₁+μ₂)} (μ₁/μ₂)^{k/2} I_|k|(2√(μ₁μ₂)) en logarithme double-double (`log_cyl_bessel_i`, corrigée au premier ordre de l'arrondi de z) ; répartition : la queue du côté de k sommée depuis k, l'autre par complément ; quantiles exacts ; moyenne, variance, asymétrie, kurtosis exactes | N₁ − N₂ (deux Poisson) |
| `zipf<I>` | `{.size, .exponent}` (1 <= N <= 2⁵³, s >= 0) | sommes Σ_{j=a}^{b} j^p (p = −s, 1 − s, 2 − s) en double-double : directes sous a₀ = max(64, 8\|p\|), Euler–Maclaurin au-delà (12 termes de Bernoulli, relativement à a₀^p ; intégrale par expm1, stable près de s = 1) ; répartition H_k/H_N, survie par la queue directe ; moyenne 1 + A/H et variance (BH − A²)/H² sans compensation | rejet par inversion (Hörmann et Derflinger) ; uniforme si s = 0 |

Mesures (références mpmath à 200 bits, table étendue jusqu'à des masses sous 2⁻¹¹³⁰) :
masse, répartition, survie, asymétrie, kurtosis sous 0,5 ulp, arrondi correct sur tous les
points ; 600 quantiles exacts, probabilités sous-normales comprises. Le générateur de
références sommait la survie par différence avec le total quand F > ½ (perte de précision dans
les queues) : corrigé, pour la gaussienne discrète aussi.
Zipf : masse, répartition, survie, moyenne, variance sous 0,5 ulp, arrondi correct sur tous les
points (N jusqu'à 2⁵⁰, s de 0 à 30), 700 quantiles exacts. Le ζ de Hurwitz de mpmath n'étant
précis qu'en absolu, les références des petites queues sont des sommes directes. La loi zêta
(N infini, s > 1) attend `riemann_zeta` dans wsl::special.

