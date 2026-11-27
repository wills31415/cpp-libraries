# Revue des interfaces : constats

Document de travail du chantier d'harmonisation (rupture attendue, version mineure 0.x).
Relevé du dépôt au 2026-11-19 ; les incohérences notées en route sont dans
[random.md §15](random.md). Chaque piste est à trancher avant d'être appliquée.

## 1. Structure transversale

### 1.1 Espaces de noms

- Les lois, processus et objets combinatoires de wsl::random sont directement dans `wsl::`
  (une centaine de noms : `normal`, `gamma`, `beta`, `image`, `independent`, `sample`,
  `pdf`…), alors que wsl::special a son espace `wsl::special`. Le paquet `wsl` de
  cpp-tensors partage le même espace et le même dossier d'en-têtes : le risque de collision
  croît avec chaque loi.
- `wsl::detail` est partagé par graded (traits, `larger`, `smaller`, `zip`…) et par
  wsl::random (des centaines d'aides) : rien ne les sépare.
- Pistes : `wsl::random::…` (avec, au besoin, des `using` dans `wsl`), et un `detail`
  par bibliothèque (`wsl::random::detail`, `wsl::graded_detail` ou `wsl::detail::graded`).

### 1.2 Organisation des en-têtes

- Des aides internes vivent dans des en-têtes de lois : `sine_dd`, `legendre` (Gauss–Legendre)
  et `add_scaled` dans von_mises.hpp ; `normal_lower_scaled` dans skew_normal.hpp ;
  `standardize`, `affine`, `standard_quantile` dans normal.hpp ; `negated` dans binomial.hpp.
  stable.hpp inclut ainsi skew_normal.hpp et von_mises.hpp pour leurs seuls détails.
- Piste : des en-têtes `detail/` dédiés (quadrature, trigonométrie double-double, queue
  normale, recherche de racine), les en-têtes de lois n'incluant que ce qu'ils utilisent.

### 1.3 Doublons internes

| Rôle | Exemplaires |
|---|---|
| opposé d'un double-double | `neg` (inverse_gaussian), `neg_dd` (von_mises), `negated` (binomial) |
| somme de positifs à l'échelle | `add_scaled` (von_mises), `scaled_sum` (inverse_gaussian) |
| différence à l'échelle | `scaled_difference` (inverse_gaussian), `scaled_difference_any` (stable) |
| complément 1 − v à l'échelle | `complement` (inverse_gaussian), lambdas locales (stable, skew_normal, hypergeometric) |
| e^l arrondi | `exp_of` ×3 (gumbel, weibull, random_matrix) |
| quantile par Newton sur log de la plus petite part, encadré, deux pas à résidu précis | inverse_gaussian, von_mises, skew_normal, stable (et truncated) |
| intégration par morceaux de Gauss–Legendre | von_mises, skew_normal, stable (adaptatif) |
| poids entiers et choix exacts | `mixture` et `detail::choose` (§15) |

Pistes : une arithmétique à l'échelle complète dans `detail/arithmetic.hpp`, un solveur
de quantiles et un intégrateur communs.

### 1.4 Style du code

- Membres `x_` partout, sauf dlist (`m_x`) ; paramètre de modèle `IsConst` (dlist) contre
  `snake_case` ailleurs.
- Concept `wsl::streamable` (dlist) exporté sous un nom très générique.
- `constexpr` partout, sauf dlist, elo et progress.

### 1.5 Erreurs

- Construction invalide : `std::domain_error` (random), `std::invalid_argument`
  (permutation), `static_assert` (array2d, dimensions nulles).
- Probabilité hors de [0, 1] dans `quantile` : `std::domain_error`.
- Valeur indéfinie : NaN (moyenne d'une stable α <= 1, somme continue sans identité), +∞
  (variance d'une stable, moyenne de Fréchet α <= 1), membre absent (moyenne de Cauchy).
- `noexcept` irrégulier : `median()` des lois dont le quantile peut lever ne l'est pas.
- Piste : une règle écrite (exception à la construction seulement ? NaN ou absence pour
  l'indéfini ? +∞ pour le divergent ?).

### 1.6 Reproductibilité hors de wsl::random

- elo rééchantillonne avec `std::mt19937` et `std::discrete_distribution` : résultats
  différents d'une bibliothèque standard à l'autre (ce que wsl::random corrige) ; `std::pow`
  pour les conversions Elo.
- progress dépend de `isatty`/`fileno` (POSIX) et écrit toujours sur `stderr`.

## 2. wsl::random

### 2.1 Paramètres nommés

| Rôle | Noms en usage |
|---|---|
| position | `.location` (cauchy, laplace, gumbel…), `.mean` (normal, inverse_gaussian), `.mu` (lognormal, von_mises, discrete_gaussian), `.low`/`.high` (uniform, triangular, arcsine), `.start` (processus) |
| échelle | `.scale`, `.sd` (normal), `.sigma` (lognormal, discrete_gaussian), `.volatility` (processus), `.rate` (exponential, gamma, erlang, poisson_process) |
| forme | `.shape`, `.alpha`/`.beta`, `.kappa`, `.degrees_of_freedom`, `.successes`, `.theta`, `.q`, `.discount` |
| taille | `.size`, `.n`, `.trials`, `.draws`, `.vertices`, `.nodes`, `.edges`, `.dimension`, `.population` |

- Défauts irréguliers : `skew_normal` a tout par défaut, `gumbel` rien, `arcsine` [0, 1],
  `frechet` `.location = 0`, `stable` `.beta = 0, .scale = 1, .location = 0`.
- Types : `T` pour les lois continues, mais `double` pour `poisson::mean`,
  `negative_binomial::successes`, `discrete_gaussian::mu/sigma`, `logarithmic::p` ;
  `std::int64_t` pour `erlang::shape` quand `gamma::shape` est `T`.
- Piste : un vocabulaire canonique (location / scale / shape, sauf noms consacrés à
  arrêter un par un), des défauts réguliers (location 0, scale 1).

### 2.2 Types des moments

- `T` pour les continues, `double` pour les discrètes, `value_type` (multivariées), `auto`
  (independent) ; `entropy()` rend `double` dans 24 lois et `T` dans 20, y compris parmi
  les continues (von Mises : `double`).
- `quantile(law, double p)` prend toujours un `double`, même pour une loi en `float`.
- Piste : moments en `T` (continues) et `double` (discrètes), ou `double` partout ; règle
  pour l'indéfini (§1.5).

### 2.3 Fonctions présentes

Matrice des trous les plus visibles (le reste est complet) :

| Manque | Lois |
|---|---|
| `entropy` | inverse_gaussian, skew_normal, stable, negative_binomial, hypergeometric, beta_binomial, logarithmic, discrete_gaussian, poisson, binomial |
| `cquantile` | toutes les discrètes |
| `skewness`, `excess_kurtosis` | beta_binomial, hypergeometric, logarithmic, discrete_gaussian, categorical, stable (indéfinis) |
| `mode` | arcsine (bimodale), skew_normal, inverse_gaussian a `mode`, stable non |
| `median` | la plupart des discrètes |

Piste : compléter ce qui a une forme close ou un calcul sûr ; sinon l'absence est la règle
(concepts `has_*`).

### 2.4 Tirage

- `sampling` (exact / approché) est déclaré par chaque loi ; `sample(law, g, wsl::random::fast)`
  existe pour une partie seulement (le défaut sinon).
- Primitives : `uniform_real(g)`, `uniform_below(g, n)`, `bernoulli_trial(g, p)`,
  `shuffle(range, g)` (générateur en premier ou en dernier selon la fonction) ; lois
  `uniform`, `bernoulli` à côté.
- Processus et objets combinatoires : §15 de random.md (pas, valeurs, constructeurs
  positionnels, types de résultats).
- Le concept `stochastic_process` prend `std::mt19937_64` comme générateur témoin.

### 2.5 Moteurs

- Graines : `std::uint64_t` (splitmix64, xoshiro256pp, philox4x64), `uint128` + flux (pcg64),
  clé + compteur (philox4x64) ; états `state()` de types variés, `key()`.
- Sauts : `jump()`/`long_jump()` (xoshiro), `advance(δ)` (pcg64), `discard(n)` partout.
- Ni `seed(seed_seq)`, ni `operator<<`/`>>` : ce sont des `uniform_random_bit_generator`,
  pas des *RandomNumberEngine* au sens du standard.
- Piste : décider si l'on vise les exigences complètes du standard.

## 3. wsl::special

- Noms mêlés : ceux du standard (`lgamma`, `tgamma`, `erf`, `beta`, `cyl_bessel_i`) et ceux
  de Boost (`ibeta`, `ibetac`, `gamma_p`, `gamma_q`, `erf_inv`, `digamma`) ; `ibeta_inv` /
  `gamma_p_inv`.
- Surcharges entières pour sqrt, cbrt, exp…, sin…, erf, lgamma, mais pas pour hypot, atan2,
  beta, ibeta, gamma_p, cyl_bessel_i ; pas de `long double`.
- Manques de `<cmath>` : asin, acos, sinh, cosh, tanh et réciproques, floor/round (wsl::random
  en a recodé un : `detail::floor_small`), fma, ldexp/frexp.
- Les noyaux double-double et à l'échelle (`special/detail/extended.hpp`) servent de fait de
  bibliothèque interne à wsl::random : à rendre explicite (nom, place, garanties).

## 4. Petites bibliothèques

| Bibliothèque | Constats |
|---|---|
| array2d | agrégat C++23 (deducing this, mdspan) ; `row_count()`/`column_count()` plutôt que `extent(0)`/`extent(1)` ; étiquettes `row_major`/`column_major` ; `strided_view` public ; dimensions nulles refusées |
| dlist | seule bibliothèque non `constexpr` ; `dump(std::ostream&)` dans l'en-tête principal (graded sépare l'affichage dans graded_io) ; style `m_`, `IsConst` ; insertion et effacement par position et par itérateur |
| graded / graded_io | `jet`, `series`, `index_jet` ; concept `lexicographic` ; `hash_of` et `std::hash` ; `value_cast`, `to_series`, `concat` |
| progress | POSIX ; `std::print` (C++23) ; `stderr` imposé ; appel `meter(done, total)` |
| elo | `<random>` non reproductible ; `tally`, `rate`, `bradley_terry`, `estimate_elo` |
| permutation | cohérente ; rang lexicographique limité à n <= 20 (`std::uint64_t`) |

Un seul exemple (`examples/dlist_demo.cpp`) pour neuf bibliothèques.

## 5. Lecture approfondie (sources lues en entier)

### 5.1 Défauts

- **dlist** : `erase(const_iterator)` et `erase(first, last)` ne compilent pas
  (`basic_iterator<true>` rendu là où `iterator` est attendu) ; aucun test ne les appelle.
- **dlist** : aucune position n'est accessible publiquement (les décalages cumulés sont
  privés) ; les tests ne vérifient les positions qu'à travers le texte de `dump`. Chaque
  opération parcourt la liste depuis le début (`node_offset`), d'où `extract` et
  `erase(first, last)` en O(n²).
- **progress** : `done * 20` déborde pour `done` > 2⁵⁹ (tranches de 5 %).
- **gamma** : construite par l'échelle, `rate()` est 1/scale arrondi (et inversement) ; deux
  paramétrisations d'une même loi ne sont pas égales (`==` compare la forme stockée).

### 5.2 Écarts de nommage

- **beta** : paramètres `.alpha`, `.beta`, accesseurs `a()`, `b()`.
- Accesseurs : chaque loi expose ses paramètres sous leurs noms (`sd()`, `kappa()`,
  `trials()`…), parfois pas du tout (`categorical`) ; aucun accès commun (`param()` de
  `<random>`).
- Constructeurs : `explicit` pour les objets combinatoires et les processus, implicites pour
  les lois (une fonction qui prend une `normal<>` accepte `{.mean = …, .sd = …}`).
- Constructions positionnelles au milieu de paramètres nommés : `truncate(X, a, b)`,
  `weighted(w, X)`, `random_walk(pas, {.start})`, `compound_poisson(taux, loi)`,
  `inhomogeneous_poisson_process(f, borne)`.
- Noms génériques dans `wsl` : `interval`, `distribution`, `dynamic`, `fast`, `exactness`,
  `image`, `weighted`, `tally`, `rate`, `streamable`, `to_string`, `sign`.

### 5.3 Recouvrements entre bibliothèques

- Matrices : `array2d<T, R, C>` (taille fixe, deux dispositions, C++23) et
  `square_matrix<T, N>` (wsl::random, taille fixe ou dynamique) ; elo passe ses tableaux
  en `std::vector<std::vector<tally>>`. Piste : un type matriciel commun.
- Vecteurs : `vector_of<T, N>` (std::array ou std::vector) dans wsl::random seulement.
- Arithmétique double-double : `wsl::special::detail::extended` importée par
  `wsl::detail` (random), plus des compléments éparpillés dans les lois.
- Quatre `struct {elo, low, high}` en deux noms dans elo (`elo_estimate`, `rating`).

### 5.4 Code et documentation divergents

| Document | Dit | Code |
|---|---|---|
| random.md §4.3 | `wsl::random::normal Z({.mean = 0, .variance = 4})` | `.sd` |
| random.md §4.3 | `X.kurtosis()` | `excess_kurtosis()` |
| random.md §4.3 | `characteristic(X, t)`, `mgf(X, t)` | absents |
| random.md §4.4 | tirage générique (inversion, rejet, ARS, alias) | absent (alias seulement dans `categorical`) |
| random.md §4.3 | attribut inexistant = erreur de compilation | moyenne de Fréchet (α <= 1) : +∞ ; d'une stable (α <= 1) : NaN |
| random.md §4.7 | adaptateur avec `<<` et `>>` | sans les flux (law.hpp le dit) |
| random.md §4.2–4.4 | marqués ◻︎ | faits |
| special.md §3 | `beta_inc` | `ibeta`, `ibetac` |

### 5.5 Style

- Gardes : `#pragma once` partout, sauf dlist, graded, graded_io (`#ifndef`).
- Commentaires : accentués partout, en ASCII dans graded et graded_io ; dlist et array2d
  sans commentaire d'en-tête ; messages d'exception en anglais dans dlist, en français
  ailleurs.
- Accolades et espaces : trois styles (`if(` dans array2d/elo/progress, `if (` ailleurs ;
  accolades sur une ligne ou non).
- graded surcharge la virgule (alias de `concat`) ; `elems` public dans graded, `elems_`
  public dans array2d (agrégats).
- Hachage de graded : constantes FNV de 64 bits dans un `std::size_t` (32 bits ailleurs).

## 6. Décisions prises

| # | Point | Décision (date) |
|---|---|---|
| 1 | Espaces de noms | `wsl::random::…` pour tout wsl::random (lois, processus, objets combinatoires, moteurs, primitives ; appels non qualifiés par ADL) ; internes `wsl::random::detail`, `wsl::special::detail`, `wsl::detail::<bibliothèque>` pour les petites (2026-11-20) |
| 2a | Construction des lois | Structure de paramètres nommés conservée (arguments nommés, plusieurs paramétrisations sans ambiguïté, défauts libres) ; ajout d'un `param()` commun (2026-11-20) |
| 2b | Vocabulaire des paramètres | Compromis : moyenne, écart type, variance en `.mean`, `.stddev` (nom de `<random>`), `.variance`, plusieurs paramétrisations au besoin (normale par l'écart type ou la variance) ; position, échelle, intensité en `.location`, `.scale`, `.rate` ; formes par leur lettre usuelle quand elle est stable (`.alpha`, `.beta`, `.kappa`, `.nu`), sinon `.shape` ; effectifs en mots (`.trials`, `.size`, `.vertices`…). Chaque loi documente sa formule et la table symbole ↔ champ ; accesseurs nommés comme les champs. Défauts : location 0, scale 1 quand une forme standard existe, aucun pour les formes ; champs dans l'ordre position, échelle, formes (2026-11-20) |
| 3a | Types des constantes | Lois continues (valeurs `T`) : tout en `T`. Lois discrètes (valeurs `I`) : médiane, mode, quantiles, support en `I` (valeurs de la loi) ; moyenne, variance, asymétrie, aplatissement, entropie, masse, répartition, survie en `double` fixe (réelles par nature ; pas de paramètre réel supplémentaire). Entropie à régulariser (2026-11-20) |
| 3b | Constantes indéfinies | Jamais définie : membre absent (erreur de compilation, Cauchy). Selon les paramètres, comme IEEE et `<cmath>`, sans exception, `noexcept` : +∞ quand elle diverge vers +∞ (variance de Student ν <= 2, moyenne de Pareto α <= 1), NaN quand elle n'est pas définie (moyenne de Student ν <= 1) (2026-11-20) |
| 3c | Erreurs | Exceptions de la bibliothèque standard selon la nature : `std::domain_error` (paramètre ou probabilité hors domaine), `std::invalid_argument` (structure invalide), `std::out_of_range` (indice) ; NaN dans wsl::special, comme `<cmath>` (2026-11-20) |
| 4 | Standard | Moteurs complets au sens de *RandomNumberEngine* (`seed(seed_seq)`, `operator<<`/`>>`) ; noms de wsl::special : ceux de `<cmath>` (fonctions spéciales de C++17 comprises), sinon ceux de Boost ; special.md corrigé (`beta_inc` → `ibeta`) (2026-11-20) |
| 5 | Lois à plusieurs paramétrisations | Une forme canonique exacte par loi ; `==` la compare : égalité en loi exacte, et deux objets égaux identiques au bit près (tirages compris). Gamma, exponentielle, Erlang : intensité rationnelle λ = P/Q·2ᵉ réduite (P, Q impairs, premiers entre eux), donnée par `.rate` ou `.scale` (double, ou `wr::ratio(p, q)`) : les deux formes exactes, un seul chemin de calcul, voies rapides inchangées quand λ ou 1/λ est un double. Normale : variance exacte (`.stddev` : σ² exact en double-double ; `.variance` : v), σ en double-double pour les fonctions, tirage exact par comparaison des carrés ; voie rapide inchangée quand σ est un double (révisée le 2026-11-21 : remplace la forme en double-double, qui rendait inexactes les lois données par leur échelle) (2026-11-20) |
| 6 | Matrices, arithmétique étendue | Pas de fusion pendant ce chantier : `square_matrix` (random) et `array2d` alignés sur des conventions communes (`extent(i)` au lieu de `size()`, `m[i, j]`) ; elo inchangé. Chantier futur `mdarray` (sémantique de P1684, étendues fixes ou dynamiques par dimension, couche d'expressions paresseuses, intégration de cpp-tensors), dont `graded` devient le type d'indice (arité exacte : rang fixe ; bornée : broadcasting). Arithmétique étendue rassemblée en un composant interne documenté (2026-11-20). Le chantier `mdarray` est remplacé (2026-11-23) par celui des objets mathématiques ([algebra.md](algebra.md)), dont le stockage multidimensionnel n'est qu'un support ; ce stockage redevient un chantier à part (2026-11-24, [mdarray.md](mdarray.md)), indépendant des objets mathématiques, et `graded` n'en est plus le type d'indice |
| 7 | Standard du dépôt | C++23 partout ; C++20 seulement sur décision explicite, pour un besoin de compatibilité identifié (exception possible, non appliquée) ; tests compilés une seule fois (2026-11-20) |
| 8a | Appel des opérations | Opérations en membres (`X.pdf(x)`, `X.cdf(x)`, `X.quantile(p)`…) et objets-fonctions qualifiés dans `wsl::random` (`wr::pdf(X, x)`, passables en argument ; second recours par ADL pour les lois tierces) ; plus d'amies cachées. Constantes en membres, avec `wr::mean(X)`… pour le code générique (2026-11-20) |
| 8b | Tirages | `X(g)` (la variable aléatoire évaluée en ω, comme `d(g)` de `<random>`) et `wr::sample(X, g)` ; `X(n, g)` → `std::vector`, mêmes valeurs et même consommation que n appels `X(g)` ; vue paresseuse `wr::samples(X, g)` ; `wr::fast` en dernier argument ; moteur en dernier partout, primitives comprises (`uniform_below(n, g)`, `bernoulli_trial(p, g)`) ; pas de `X()` ni de moteur par défaut (2026-11-20) |
| 8c | Opérations sur les lois | Plus d'opérateurs arithmétiques sur les lois (réservés aux variables aléatoires) ; noms de la théorie de la mesure : `wr::affine(P, {.scale, .shift})`, `wr::convolution(P, Q)` (type `convolution_law<L1, L2>` : un type et une fonction ne peuvent partager un nom ; différence par `affine`), `wr::truncate(P, {.lower, .upper})` (bornes facultatives), `wr::mixture(wr::weighted(w, P), …)`, `wr::pushforward(P, f)` (avec `increasing`/`decreasing`), `wr::product_measure(P, Q, …)` (2026-11-20) |
| 8d | Densité | Une seule notion, la dérivée de Radon–Nikodym dP/dμ : `X.density(x)`, `X.log_density(x)`, `wr::density(X, x)` pour toutes les lois ; chaque loi déclare sa mesure de référence (`reference_measure` : `wr::lebesgue`, `wr::counting`, `wr::surface` pour les sphères…, combinaisons pour `product_measure` et `mixture`) ; `pdf` (Lebesgue) et `pmf` (comptage) restent comme noms usuels (2026-11-20) |
| 9 | elo | Espace `wsl::elo` (noms trop génériques dans `wsl::`) : `tally`, `from_score`, `to_score`, `pairwise`, `bradley_terry`, `rate`, une seule structure `rating {elo, low, high}`, `format_signed` ; rééchantillonnage par `<random>` gardé (non reproductible d'une bibliothèque standard à l'autre, choix assumé) (2026-11-21) |
| 10 | Reproductibilité (règle du dépôt) | Trois garanties, nommées de la même façon dans toutes les bibliothèques. **Déterministe** (défaut, partout) : mêmes bits pour une même configuration (compilateur et version, options, jeu d'instructions, version de la bibliothèque), quels que soient le lancement, le nombre de fils, l'ordonnancement et l'alignement en mémoire ; arbres de réduction fixés par les données, jamais par les fils. **Précision** (toutes les fonctions numériques) : bornes d'erreur documentées et testées (ulp, contre mpmath), indépendantes du niveau de reproductibilité. **Portable** (sur décision explicite, par bibliothèque ou par fonction) : mêmes bits sur toute configuration conforme à IEEE 754, hors `-ffast-math`, mise à zéro des sous-normaux et précision étendue x87 (à refuser quand c'est détectable) ; obtenue d'office par l'arrondi correct (résultat unique) ; ailleurs par une suite d'opérations fixée, des barrières contre la contraction, un ordre indépendant de la machine et sans libm du système ; se transmet par les dépendances (une partie portable ne s'appuie que sur des parties portables). L'exactitude bit à bit entre toutes les configurations n'est plus exigée partout (2026-11-26) |
| 11 | Charte de nommage : espaces de noms | `wsl` n'est qu'un préfixe d'auteur (comme `std`, `boost`) ; modèle de la bibliothèque standard : types de vocabulaire (servant à travers les domaines) et étiquettes communes (modes) directement dans `wsl::`, un espace par domaine aux noms nombreux ou génériques, au nom de sa bibliothèque (`wsl::random`, `wsl::special`…), `detail` par bibliothèque (décision 1). Liste des noms admis dans `wsl::` à fixer plus tard (pistes : `mdarray` et ses vues, `permutation`, `ragged`, étiquettes de modes, scalaires étendus) ; vérifier alors que les étiquettes dans `wsl::` ne bloquent rien. Principes : un nom exprime ce qu'il représente ; l'usage par défaut s'écrit sans étiquette, un autre mode est un cas particulier marqué (2026-11-26) |
| 12 | Charte de nommage : bibliothèques | Le nom dit ce que la bibliothèque contient : un type principal → ce type au singulier (`mdarray`, `permutation`) ; une collection de choses de même nature → pluriel (`preconditioners`, `engines`), comme `std::ranges`, `std::numbers` ; un domaine → le nom du champ (`linalg`, `geometry`, `statistics`). L'espace de noms suit le domaine, la bibliothèque suit les dépendances : un domaine peut se scinder en plusieurs bibliothèques qui partagent `wsl::<domaine>`, celle qui porte le nom du domaine les regroupe toutes (`<wsl/random>` et `wsl::random` : tout ; `<wsl/random/engines>` et `wsl::random_engines` : la partie). Renommages libres jusqu'à la première version stable (aucun utilisateur aujourd'hui). À suivre : imbrication domaine / sous-domaine (arborescence), abréviations et noms intérieurs (2026-11-26) |
| 13 | Charte de nommage : abréviations et niveaux | Les noms de domaines et de bibliothèques décrivent (on doit les trouver) : mots entiers par défaut ; abréviation seulement sur une liste fermée (nom dominant dans la littérature et un seul sens en calcul scientifique) ; un alias canonique par domaine long, recommandé par la charte et employé dans toute la documentation, défini chez l'utilisateur (`namespace sf = wsl::special_functions;`, comme `np` pour NumPy). Les noms intérieurs suivent la notation : opérations par la notation du domaine (`lgamma`, `erfc`, `lu`, `svd`, `grad`), sinon mots entiers (`matrix_product`) ; types en mots entiers descriptifs (`lu_decomposition` ; pour les décompositions, voir la décision 18) ; étiquettes par la garantie (`portable`). Un glossaire commun fixe chaque abréviation et chaque mot partagé. Liste fermée : `linalg`, `simd`, `meta` (standard), `mdarray`, `ode`, `pde`, `sde`, `io` (`fft` retiré, décision 20 ; `sde` ajouté, décision 23) ; refusées : `fp`, `uq`, `fem`, `autodiff`, `stats` (en mots entiers). Alias courts sans préfixe (`sf`, `la`), libres de l'ambiguïté interdite aux noms officiels ; table des domaines et alias : décision 20 (2026-11-26) |
| 14 | `random` et `probability` | Deux domaines : `wsl::random` (alias `rnd`) engendre le hasard : moteurs et **fonctions** de tirage simples (`random_bits`, `uniform_below`, `uniform_int`, `uniform_real`, `bernoulli_trial`, `shuffle`) ; `wsl::probability` (alias `proba`) contient les **objets** mathématiques : lois (mesures : densité, répartition, quantile, moments, tirage `X(g)` appuyé sur les primitives, mêmes bits), variables aléatoires, processus. Comme Boost (random / math::distributions), Julia (Random / Distributions.jl). Frontière provisoire, à améliorer au besoin ; vérifier à la refonte que les primitives exactes ne font pas dépendre `random` des couches du dessus (2026-11-26) |
| 15 | Charte de nommage : placement et appel | (1) Un en-tête par notion, aucun fichier fourre-tout ; en-têtes d'ensemble qui ne font qu'inclure. (2) Un type va dans `wsl::` s'il sert à travers les domaines, sinon dans son domaine. (3) Ses opérations vont avec lui : propriétés en membres (`p.degree()`, `A.extent(0)`), opérations mathématiques en fonctions libres (`wsl::roots(p)`). (4) Une fonction sur des types qui ne sont pas les nôtres va toujours dans un domaine (`wsl::special_functions::gamma(x)`). (5) Rien de générique non contraint dans un espace partagé. (6) Registre des noms de `wsl::` dans la charte, complété et vérifié à chaque nouvelle bibliothèque. (7) Appels qualifiés par défaut, code et documentation (`la::lu(A)`). (8) Opérations publiques en objets-fonctions, passables en argument (`std::ranges::transform(xs, out, sf::gamma)`), second recours par ADL (8a). (9) Écriture courte au cas par cas par `using wsl::dot, wsl::grad;` dans une fonction. Conséquence : pas de domaines `polynomials` ni `tensors` (types `wsl::polynomial`, `wsl::tensor`, placements provisoires ; familles orthogonales dans les fonctions spéciales) ; `units` reste un domaine pour son catalogue (2026-11-26) |
| 16 | Structures algébriques ; formation des alias | Domaine `wsl::algebraic_structures` (alias `algst`) : concepts des structures (A9) et leur déclaration ; pas dans `wsl::` (`field` y heurterait le champ physique, `module` le mot de C++20, `ring`, `group`). Formation des alias : (1) domaine d'un mot → troncature reconnaissable (`proba`, `quad`) ; (2) plusieurs mots → initiales (`sf`, `la`, `dg`) ; (3) si c'est un mot courant, un nom de variable fréquent ou un autre nom du dépôt → troncature des mots (`algst`) ; (4) vérification au registre (noms du dépôt, unités, liste fermée) ; (5) pas de fragment au sens établi en programmation (`str`, `ptr`, `fn`, `int`, `num`, `idx`). À appliquer à toute la table à la fin de la discussion (`cm`, `cg`, `fd`, `opt` à revoir) (2026-11-26) |
| 17 | Différences, volumes et éléments finis | Placés sous `wsl::pde` par la décision 22. Domaines nommés par leur champ, la méthode, comme dans la littérature (« *finite element method* », adjectif au singulier) : `wsl::finite_difference_method` (alias `fdm`), `wsl::finite_volume_method` (`fvm`), `wsl::finite_element_method` (`fem`) ; pas de répétition dans `fem::elements::lagrange<2>`. Règle des alias, deux points de plus : (6) un alias peut être le sigle établi du domaine dans la littérature, même s'il ne reprend pas exactement ses initiales ; (7) des domaines frères ont des alias de même forme (2026-11-27) |
| 18 | Charte de nommage : imbrication et catalogues | (1) Un domaine est plat par défaut ; chaque niveau d'imbrication se justifie (pas de profondeur chiffrée) ; les options d'un algorithme sont des membres de son type (énumérations, types imbriqués), jamais un espace de noms : `la::preconditioners::amg{.coarsening = ...}`. (2) Un **catalogue** (choix interchangeables qui modélisent un même concept) a toujours son sous-espace au pluriel, le concept au singulier dans le domaine : `rnd::engine` et `rnd::engines::pcg64`, `la::preconditioner` et `la::preconditioners::ilu0`, `proba::laws::normal`, `fem::elements::lagrange<2>` ; on découvre ainsi le catalogue (auto-complétion, documentation), la règle est uniforme et stable, et le mot de plus ne s'écrit qu'à la construction. Les opérations (non interchangeables) restent à plat : `sf::gamma`, `la::det`, `la::solve`. (3) Pas de sous-espace sans catalogue : le creux est une propriété du stockage (`la::csr_matrix`, opérations surchargées), pas de `la::sparse::` ; pas d'alias de sous-espace. (4) `linalg`, catalogues par rôle comme Ginkgo, AMGCL, PETSc, SciPy : `la::solvers` (richardson, chebyshev, cg, gmres, bicgstab, minres), `la::preconditioners` (jacobi, gauss_seidel, sor, ssor, ilu0, ilut, ic0, amg), `la::eigensolvers` (lanczos, arnoldi, lobpcg, qr, jacobi), `la::decompositions` (lu, cholesky, ldlt, qr, svd, schur, hessenberg, polar) ; chaque objet mathématique une seule fois, quel que soit son rôle : la méthode de Jacobi est `solvers::richardson{.preconditioner = preconditioners::jacobi{}}` (découpage A = M − N, Saad ; Barrett et coll.), un lisseur est une option de `amg`. (5) Une décomposition est un objet du catalogue, construit par déduction des arguments (`auto F = la::decompositions::lu(A);`, `auto [U, S, V] = la::decompositions::svd(A);`, forme courte locale par `using`) ; pas de fonction à plat pour l'instant, ajout possible plus tard sans rupture (`la::lu(A)` rendant l'objet). Restent ouverts : l'en-tête d'un sous-espace (`<wsl/linalg/solvers>`, question de l'arborescence) et la place de la matrice creuse (registre) (2026-11-27) |
| 19 | Mémoire et exécution | Deux domaines, modèle **séparé** (comme std, Kokkos, RAJA et Umpire, Thrust ; Ginkgo et SYCL couplent) : `wsl::memory` (alias `mem`) dit où résident les données (concept `mem::space`, catalogue `mem::spaces`, allocation, copies explicites) ; `wsl::execution` (alias `exec`) dit où s'exécute le calcul et comment il se découpe (concept `exec::space`, catalogue `exec::spaces`, boucles et réductions parallèles déterministes, concept d'accessibilité vérifié à la compilation) ; `execution` dépend de `memory`, mdarray seulement de `memory`. Raisons : durées de vie différentes (la mémoire suit la donnée, l'exécution se choisit par calcul), mémoire unifiée (une mémoire, plusieurs espaces d'exécution), dépendances plus légères. Défaut sans étiquette : l'hôte (`wsl::mdarray<double, 2>`). Écartés : `parallel` (partiel), `gpu` (un appareil), `hardware`/`backends`/`platform` (fourre-tout), un domaine couplé. Les alias gardent le sens établi de leur fragment (point 5 respecté). Exigences : simple pour l'utilisateur, compatibilités vérifiées à la compilation, options les plus rapides, reproductibilité (décision 10). Contenu, mot `device` et programmation du GPU : plus tard, avant le lot 2 (2026-11-27) |
| 20 | Table des domaines et des alias ; domaines parents | **Règle des parents** : on regroupe sous un domaine parent les domaines qui s'utilisent séparément et partagent un vocabulaire (jamais par le seul mot) ; les domaines mathématiques, combinés dans presque tout programme, restent à plat. La profondeur ne se paie qu'à la déclaration de l'alias : **tout domaine utilisé seul a son alias canonique, quelle que soit sa profondeur** (complète la décision 18, où les sous-domaines n'avaient pas d'alias). **Physique** : `wsl::physics` (`phys` ; constantes physiques, lien avec `units` à voir), puis les grandes théories (Landau et Lifchitz ; archives d'arXiv) : `classical_mechanics` (`classmech`), `continuum_mechanics` (`contmech`) avec les sous-domaines `fluids` et `solids` (alias `fluids`, `solids`), `electromagnetism` (`emag`), `statistical_physics` (`statphys`), `quantum::mechanics` (`qm`), `quantum::field_theory` (`qft`), `quantum::information` (`qinfo`), `relativity::special` (`srel`), `relativity::general` (`grel`), `relativity::numerical` (`numrel`) ; les domaines à plusieurs parents directement sous `physics` : `plasma_physics` (`plasma`), `cosmology` (`cosmo`), `condensed_matter` (`condmat`), `molecular_dynamics` (`moldyn`), `acoustics`. **Autres domaines** (règle des alias appliquée) : `floating_point_arithmetic` (`fpa` ; `fp` heurte `FILE* fp`), `elementary_functions` (`ef`), `special_functions` (`sf`), `algebraic_structures` (`algst`), `quadrature` (`quadr` ; `quad` désigne la précision quadruple), `interpolation` (`interp`), `fourier_analysis` (`fourier`, opération `fourier::fft(x)` ; remplace `fft`, nom d'un algorithme ; alias modifié par la décision 23), `automatic_differentiation` (`ad`), `symbolic` (`sym`), `random` (`rnd`), `probability` (`proba`), `statistics` (`stats`, l'usage de R et SciPy), `uncertainty_quantification` (`uq`), `linalg` (`la`), `lie_groups` (`lie`), `geometric_algebra` (`ga`), `units`, `computational_geometry` (`cgeom` ; `cg` est `la::solvers::cg`), `calculus` (`calc`), `optimization` (`optim` ; `opt` désigne une option), `differential_geometry` (`dg`), `meshes` (`msh` ; `mesh` est le nom de variable), `finite_*_method` (`fdm`, `fvm`, `fem`), `memory` (`mem`), `execution` (`exec`), `ode`, `pde`, `io`. Domaines mathématiques et numériques probables : décision 23 (2026-11-27) |
| 21 | Équations non linéaires | Un seul domaine `wsl::nonlinear_equations` (alias `nle`) pour F(x) = 0, scalaire et systèmes, comme la MSC 65H (65H05 équation seule, 65H10 systèmes, 65H20 homotopie), Ortega et Rheinboldt, Kelley, Deuflhard, NonlinearSolve.jl ; il absorbe `root_finding`. Un seul point d'entrée, symétrique de `la::solve` : `nle::solve(f, nle::bracket{a, b}, nle::solvers::brent{})`, `nle::solve(F, x0, nle::solvers::newton{.linear = la::solvers::gmres{}})` ; les exigences de chaque solveur (scalaire encadré, dérivée ou jacobienne) sont des concepts vérifiés à la compilation. Les racines d'un polynôme restent avec son type (`wsl::roots(p)`, décision 15). Optimisation : domaine frère, problème distinct malgré le vocabulaire commun (Dennis et Schnabel) (2026-11-27) |
| 22 | EDO, EDP, intégration en temps | Pas de domaine `time_integration` : par la méthode des lignes (MSC 65M20 ; Hundsdorfer et Verwer ; LeVeque), une EDP d'évolution discrétisée en espace est un système d'EDO u′ = F(u), et avancer en temps, c'est résoudre cette EDO (IMEX, SSP, découpage de Strang compris) ; la MSC n'a pas de rubrique « intégration en temps » (65L EDO, 65M et 65N EDP). `wsl::ode` : point d'entrée `ode::solve`, catalogue `ode::integrators` (concept `ode::integrator` ; `rk4`, `dopri5`, `bdf`, `radau5`, `imex_ark4`, `ssprk3`, `stormer_verlet`) pour les problèmes à valeur initiale, équations algébro-différentielles comprises ; un catalogue se nomme par ce que sont ses membres (« *integrator* » : Hairer, Lubich et Wanner ; SUNDIALS), pas par symétrie avec `la::solvers` ; l'uniformité porte sur le verbe (`X::solve`) et la structure (concept, catalogue, vérifications à la compilation). Les problèmes aux limites auront leur propre catalogue. `wsl::pde` : vocabulaire commun (conditions aux limites et initiales) et parent, par la règle des parents (décision 20), des méthodes de discrétisation : `wsl::pde::finite_difference_method`, `finite_volume_method`, `finite_element_method`, alias inchangés (`fdm`, `fvm`, `fem`) ; le pont est `fem::semi_discretize(...)`, qui rend un problème pour `ode::solve`. `ode` et `pde` restent frères, puisqu'ils se combinent (2026-11-27) |
| 23 | Domaines mathématiques et numériques probables | Classés par la MSC2020, règle des alias appliquée : `number_theory` (`nt`), `combinatorics` (`combin`), `complex_analysis` (`ca`), `dynamical_systems` (`dynsys`), `topology` (`topo`), `algebraic_topology` (`algtop`), `integral_equations` (`inteq`), `calculus_of_variations` (`calcvar`), `control_theory` (`ctrl`), `information_theory` (`infoth`), `signal_processing` (`sigproc`), `wavelets`, `sde` (liste fermée, frère de `ode` et `pde`), `monte_carlo` (`mc`), `interval_arithmetic` (`ia`), `inverse_problems` (`invprob`), `data_assimilation` (`da`), `model_order_reduction` (`mor`), `machine_learning` (`mlearn`) ; sous `pde` : `spectral_method` (`spm`), `discontinuous_galerkin_method` (`dgm`), `boundary_element_method` (`bem`). La famille « … *analysis* » n'est pas une fratrie (même mot, pas même vocabulaire) : chacun prend l'alias le plus lisible ou d'usage établi : `functional_analysis` (`fa`, le sigle des mathématiciens), `fourier_analysis` (`fourier`, le nom propre suffit), `complex_analysis` (`ca`), `convex_analysis` (`cvx`, comme CVX et CVXPY). En discussion : `differential_geometry` comme parent, `approximation_theory`, `measure_theory`, `graph_theory` (2026-11-27) |

### Pistes en réserve (pour une refonte ultérieure)

- **Variables aléatoires** (chantier futur, souhaité) : `wr::variable(P)`, feuilles munies
  d'une identité, expressions paresseuses (`2 * X + Y`), tirage couplé (`(X - X)(g) == 0`),
  loi exacte quand l'algèbre des lois la connaît (formes affines de feuilles indépendantes,
  identités), moments par linéarité.
- Opérations sur les lois à venir : `maximum`/`minimum` (F_P F_Q), sommes aléatoires
  `compound(N, P)`, marges `marginal<k>(J)`, conditionnement sur un événement.

- Construction : constructeurs positionnels (comme `<random>`), fabriques nommées
  (`gamma::with_rate(2, 1)`), types forts (`gamma(shape{2}, rate{1})`, ordre libre), ou
  variante mixte (positionnel pour les lois à un paramètre sans ambiguïté, `bernoulli(0.3)`,
  `poisson(4)`). Écartées pour l'instant : conventions silencieuses ou surcoût de types ;
  la structure impose l'ordre de déclaration des champs (C++20).

## 7. Décisions à prendre d'abord

1. Espaces de noms : `wsl::random::…` ou `wsl::…` ; un `detail` par bibliothèque.
2. Vocabulaire des paramètres et leurs défauts.
3. Types des moments et règle de l'indéfini et des erreurs.
4. Niveau d'ambition envers le standard (moteurs complets, noms de `<cmath>` ou de Boost).

5. Défauts à corriger sans attendre la refonte (dlist, progress) ou avec elle.
6. Recouvrements : type matriciel commun, arithmétique étendue partagée et son nom.

Ensuite : bibliothèque par bibliothèque, une entrée de rupture au CHANGELOG chacune.
