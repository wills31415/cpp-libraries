# wsl::random — les lois de probabilité, du point de vue mathématique

> Document de conception. Décisions ✅ prises avec l'utilisateur le 2026-11-16 ;
> les propositions ◻︎ ont été validées en bloc le même jour.
>
> **État** : phase 1 (`wsl::special`) terminée ; phase 2 (socle) faite le 2026-11-17 :
> moteurs et primitives uniformes exactes (§8) ; **phase 3 en cours** : architecture et
> premières lois (§9).

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
| Paramètres | Nommés et vérifiés : `gamma({.shape = 2.5, .rate = 2.0})` ; `std::domain_error` hors domaine | ✅ |
| Exactitude | **Exact par défaut** quand un algorithme exact est connu ; variante rapide sur demande explicite | ✅ |
| Types flottants | `float` et `double` ; `long double` refusé à la compilation | ✅ |
| Moteurs | Tout moteur standard, plus xoshiro256++, PCG64, Philox | ✅ |
| Dimension (multivarié) | Fixe ou dynamique, sur le modèle de `std::span` (`dirichlet<double, 3>`, `dirichlet<double, dynamic>`) | ◻︎ |
| Processus | Vues paresseuses (`std::ranges`) | ◻︎ |

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
  ├─ algebra      image, troncature, mélange, produit ; identités connues
  ├─ generic      tirage déduit de ce qui est connu (inversion, rejet…)
  ├─ adaptor      as_distribution : une loi vue comme une distribution de <random>
  └─ tests        tests statistiques (χ², Kolmogorov–Smirnov, Anderson–Darling)
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

Tous ont un état et une initialisation entièrement spécifiés, donc reproductibles.

### 4.2 Uniformes, la fondation ◻︎

Tout le reste en dépend ; c'est là que l'exactitude commence.

- **Entier uniforme sur [0, n)** : méthode de Lemire (multiplication et rejet, presque
  sans division), exactement uniforme, sur 64 et 128 bits. Remplace `exact_random`.
- **Réel uniforme sur [0, 1)**, exact : chaque `double` de [0, 1) est atteint avec
  exactement la probabilité de son intervalle d'arrondi, y compris les plus petits
  (exposant tiré selon une loi géométrique, puis mantisse). Le classique
  `(x >> 11) * 2⁻⁵³` n'atteint que 2⁵³ valeurs régulièrement espacées : ce sera la
  variante rapide.
- **Bernoulli(p) exacte** pour p rationnel, ou p donné par ses bits (on ne consomme que
  les bits nécessaires : 2 en moyenne).
- **Mélange de Fisher–Yates** exact, à partir de l'entier uniforme.

### 4.3 Le concept de loi ◻︎

Une loi est un type qui déclare :

- `value_type` : `double`, `std::int64_t`, `std::array<double, 3>`, une permutation… ;
- ses **paramètres**, dans une ou plusieurs structures nommées (une par
  paramétrisation), vérifiés à la construction ;
- son **support** ;
- un **niveau d'exactitude** du tirage (`exact`, `approximate` avec la borne d'erreur
  documentée) ;
- et, **quand ils sont connus**, ses autres attributs mathématiques, chacun détecté par
  un concept :

| Attribut | Forme | Exemple |
|---|---|---|
| Densité / fonction de masse | `pdf(X, x)`, `log_pdf(X, x)` | toutes les lois classiques |
| Fonction de répartition, survie | `cdf(X, x)`, `ccdf(X, x)` | normale, gamma… |
| Quantile | `quantile(X, p)` | exponentielle, Cauchy… |
| Moments | `X.mean()`, `X.variance()`, `X.skewness()`, `X.kurtosis()` | absents pour Cauchy : l'appel ne compile pas |
| Mode, médiane, entropie | `X.mode()`, `X.median()`, `X.entropy()` | |
| Fonctions caractéristique et génératrice | `characteristic(X, t)`, `mgf(X, t)` | |

Un attribut **inexistant** mathématiquement (l'espérance d'une Cauchy) n'est pas une
valeur spéciale : c'est une **erreur de compilation**. Un attribut existant mais **non
implémenté** est absent de la même façon, et la documentation de la loi le dit.

Exemple d'utilisation :

```cpp
wsl::gamma X({.shape = 2.5, .rate = 2.0});
wsl::gamma Y({.shape = 2.5, .scale = 0.5});   // la même loi, autre paramétrisation
wsl::normal Z({.mean = 0, .variance = 4});    // N(μ, σ²) des manuels

X.mean();                // 1.25
cdf(X, 1.0);             // P(X ≤ 1)
quantile(Z, 0.975);      // 3.919…
auto x = sample(X, gen); // tirage exact si la loi le permet

wsl::gamma({.shape = 2.5});              // erreur de compilation : paramétrisation ambiguë
wsl::gamma({.shape = -1, .rate = 1});    // std::domain_error : shape doit être > 0
```

Pourquoi des paramètres nommés : `gamma(2.5, 2.0)` a pour espérance 5 en C++ et en
NumPy (second paramètre = échelle), 1,25 en R (second paramètre = intensité). En
positionnel, se tromper de convention donne une loi fausse **sans aucune erreur**. Le
mécanisme fonctionne avec GCC et Clang (vérifié).

### 4.4 Tirage générique ◻︎

Une loi peut fournir son propre algorithme de tirage ; sinon, la bibliothèque en déduit
un de ce qui est connu :

| Ce qui est connu | Tirage déduit |
|---|---|
| Quantile | Inversion : `quantile(X, U)` |
| Densité bornée sur un support borné | Rejet sous un rectangle |
| Densité log-concave | Rejet adaptatif (ARS) |
| Fonction de masse sur un support fini | Méthode de l'alias (Walker–Vose), exacte avec des poids rationnels |

Les tirages déduits sont marqués `approximate` (le quantile ou la densité est calculé
en flottant), avec la précision qui en découle.

### 4.5 Exactitude ✅ (définition ◻︎)

Un tirage est **exact** quand la valeur rendue est l'arrondi, selon une règle
documentée, d'un réel tiré *exactement* selon la loi : la seule approximation est la
représentation finale en `float` ou `double`.

Exemples d'algorithmes exacts connus : uniforme entière et réelle, Bernoulli
rationnelle, géométrique, exponentielle (von Neumann), normale (Karney 2016),
gaussienne discrète (Canonne–Kamath–Steinke 2020), binomiale et Poisson par des
méthodes exactes, mélange de Fisher–Yates.

Quand un algorithme exact existe, il est utilisé par défaut ; la variante rapide se
demande explicitement (`sample(X, gen, wsl::fast)`) et porte sa borne d'erreur.
Quand aucun n'est connu, le meilleur algorithme est utilisé et marqué `approximate`.

### 4.6 Algèbre des lois ◻︎

| Construction | Écriture | Ce qui est hérité |
|---|---|---|
| Image affine | `2 * X + 1` | densité, répartition, quantile, moments transformés exactement |
| Image par une fonction | `pushforward(X, f)` | tirage ; le reste seulement si `f` est déclarée monotone et dérivable |
| Troncature, conditionnement | `truncate(X, a, b)`, `condition(X, event)` | densité, répartition, quantile renormalisés |
| Mélange | `mixture({w1, X1}, {w2, X2})` | densité, répartition, moments |
| Produit (indépendance) | `product(X, Y)` | loi jointe ; densité produit |
| Somme indépendante | `X + Y` (lois de même famille) | identités : normale + normale, gamma + gamma de même intensité, Poisson + Poisson… |

Les identités connues sont résolues **à la compilation** : `normal + normal` a pour type
`normal`, avec ses paramètres calculés. Sans identité, le résultat est une loi
générique (tirage par composition ; densité par convolution seulement si elle est
déclarée).

### 4.7 Adaptateur vers `<random>` ✅

`wsl::as_distribution(X)` satisfait l'exigence *RandomNumberDistribution* de la norme
(`result_type`, `param_type`, `operator()(g)`, `operator()(g, p)`, `min`, `max`,
`reset`, `==`, `<<`, `>>`) : utilisable partout où la norme attend une distribution.
La norme exige un `result_type` arithmétique : l'adaptateur est donc réservé aux lois
scalaires.

## 5. Catalogue visé (non exhaustif)

| Famille | Lois |
|---|---|
| Discrètes univariées | Bernoulli, binomiale, géométrique, binomiale négative, Poisson, hypergéométrique, uniforme discrète, catégorielle, Zipf, zeta, logarithmique, Skellam, bêta-binomiale, gaussienne discrète |
| Continues univariées | uniforme, normale, log-normale, exponentielle, gamma, χ², Erlang, bêta, Student, Fisher, Cauchy, Laplace, logistique, Weibull, Pareto, Gumbel, Fréchet, Rayleigh, Rice, von Mises, triangulaire, arcsinus, inverse-gaussienne, normale asymétrique, stables (α, β) |
| Multivariées | Dirichlet, multinomiale, normale multivariée, Student multivariée, Wishart et inverse-Wishart, uniformes sur la sphère, la boule et le simplexe, von Mises–Fisher, matrices aléatoires (Haar orthogonale et unitaire, GOE, GUE) |
| Combinatoires | permutations (uniformes, Ewens, Mallows), sous-ensembles, combinaisons, partitions d'entiers et d'ensembles, compositions, graphes (Erdős–Rényi, Barabási–Albert, graphes réguliers), arbres (Cayley, binaires de Catalan), mots |
| Processus | marches aléatoires, chaînes de Markov (discrètes et continues), processus de Poisson (homogène, inhomogène, composé), mouvement brownien et pont brownien, processus de Dirichlet, restaurant chinois, urnes de Pólya |

## 6. Tests ◻︎

Les tirages étant reproductibles, **les tests le sont aussi** : aucun test statistique
« qui échoue une fois sur cent », la graine est fixée.

1. **Valeurs de référence au bit près** : pour chaque loi, les premiers tirages de
   graines fixées sont versionnés. Les trois environnements de la CI (x86-64 Linux GCC
   et Clang, ARM64 macOS), plus un job FMA (`-march=x86-64-v3`), doivent les reproduire.
2. **Adéquation** : χ², Kolmogorov–Smirnov, Anderson–Darling contre la répartition de
   la loi, sur de grands échantillons.
3. **Cohérence mathématique** : `cdf(quantile(p)) ≈ p`, intégrale de la densité ≈ 1,
   moments empiriques ≈ moments théoriques, identités de l'algèbre vérifiées par tirage.
4. **Paramètres** : bords de domaine, valeurs refusées, paramétrisations équivalentes.

## 7. Feuille de route

| Phase | Contenu |
|---|---|
| 1 | `wsl::special`, socle : fonctions élémentaires bit-exactes, `lgamma`, `erf`, `erf_inv` ([special.md](special.md)) |
| 2 | `wsl::random`, socle : moteurs, uniformes exacts, Bernoulli exacte ; `exact_random` est absorbé |
| 3 | Concept de loi, tirage générique, adaptateur ; premières lois : uniformes, Bernoulli, binomiale, géométrique, Poisson, exponentielle, normale, gamma, bêta |
| 4 | `wsl::special`, suite : gamma et bêta incomplètes et inverses (répartition et quantiles de gamma, χ², bêta, Student, Fisher) |
| 5 | Algèbre des lois |
| 6 | Lois multivariées |
| 7 | Objets combinatoires (avec `wsl::permutation`) |
| 8 | Processus |

Ensuite : migration des projets qui utilisent les lois de `<random>` (14 fichiers dans
cpp-traversal, cpp-sudoku et puissance-4, dont `random_agent`, dont la promesse « même
graine, mêmes coups » n'est aujourd'hui tenue que sur une même plateforme) et de
`wsl::elo` (rééchantillonnage par `std::discrete_distribution`).

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
- `uniform_below(g, n)` : méthode de Lemire, exacte ; une valeur presque toujours.
- `uniform_int(g, a, b)` : [a, b], jusqu'à la pleine plage de 64 bits.
- `uniform_real<T>(g)` : **U arrondi vers le bas**, chaque flottant x de [0, 1) rendu avec
  probabilité x⁺ − x exactement (sous-normaux compris, 1 jamais) ; exposant géométrique
  lu sur les zéros de tête, mantisse prise dans le même mot quand il en reste assez :
  une valeur consommée presque toujours. `uniform_real_fast<T>(g)` : 53 (24) bits de
  tête, une valeur.
- `bernoulli_trial(g, p)` : U < p par comparaison des bits, mot par mot (une valeur presque
  toujours) ; `bernoulli_trial(g, num, den)` : `uniform_below(g, den) < num`.
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
`sample(X, g)`, `sample(X, g, wsl::fast)`, `pdf`/`pmf`, `log_pdf`/`log_pmf`, `cdf`,
`ccdf`, `quantile`. Concepts `has_mean`, `has_cdf`… pour tester leur présence. Lois
templates sur le type des valeurs, `double` par défaut (`wsl::normal X({…})`,
`wsl::normal<float>`). `wsl::as_distribution(X)` : adaptateur *RandomNumberDistribution*.

**Lois de l'étape A** (`<wsl/random/basic_laws.hpp>`) :

| Loi | Paramètres | Tirage | Attributs |
|---|---|---|---|
| `uniform<T>` | `{.low, .high}`, [low, high) | approximate : low + U (high − low), un arrondi | moyenne, variance, médiane, asymétrie, aplatissement, entropie ; pdf, log_pdf, cdf, ccdf, quantile |
| `discrete_uniform<I>` | `{.low, .high}`, bornes comprises | exact (`uniform_int`) | moyenne, variance, médiane inférieure, asymétrie, entropie ; pmf, log_pmf, cdf, ccdf, quantile exact |
| `bernoulli` | `{.p}` | exact (`bernoulli_trial`) | moyenne, variance, entropie ; pmf, cdf, ccdf |

Précision : calculs en double-double arrondis une fois (`detail/arithmetic.hpp`), mis
à l'échelle pour que les produits exacts le restent ; l'intervalle `[−DBL_MAX, DBL_MAX)`
fonctionne (b − a, qui déborde, est pris sur les moitiés).

Quantile d'une loi discrète : le **quantile exact** de la vraie loi, min{k : F(k) >= p}.
La cdf calculée étant arrondie, `quantile(cdf(k))` vaut k ou k + 1. (Rendre le quantile
cohérent avec la cdf arrondie a été essayé : intenable au-delà de 2⁵³ valeurs, où des
milliers de k voisins ont la même cdf arrondie.)

La primitive `bernoulli(g, p)` de la phase 2 est devenue `bernoulli_trial(g, p)`, pour
laisser le nom à la loi.

### Étape B : exponentielle et normale

| Loi | Paramètres | Tirage (exact) | `wsl::fast` |
|---|---|---|---|
| `exponential<T>` | `{.rate}` | exponentielle de von Neumann sur réels paresseux, arrondi correct de E/λ | inversion −log(U)/λ |
| `normal<T>` | `{.mean, .sd}` | demi-normale par rejet exponentiel exact, signe uniforme, arrondi correct de μ + σZ | ziggourat (256 couches) |

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
  `wsl::fast` 13 ns et 8 ns (`std::exponential_distribution` 15 ns,
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
| `gamma<T>` | `{.shape, .rate}` ou `{.shape, .scale}` (le paramètre donné est gardé tel quel) | pdf, log_pdf ; moyenne, variance, mode, asymétrie, aplatissement | Marsaglia–Tsang (normale par la ziggourat) ; a < 1 : G(a + 1) U^{1/a} |
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
