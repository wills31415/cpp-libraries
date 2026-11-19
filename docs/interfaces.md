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

- `sampling` (exact / approché) est déclaré par chaque loi ; `sample(law, g, wsl::fast)`
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
| random.md §4.3 | `wsl::normal Z({.mean = 0, .variance = 4})` | `.sd` |
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

## 6. Décisions à prendre d'abord

1. Espaces de noms : `wsl::random::…` ou `wsl::…` ; un `detail` par bibliothèque.
2. Vocabulaire des paramètres et leurs défauts.
3. Types des moments et règle de l'indéfini et des erreurs.
4. Niveau d'ambition envers le standard (moteurs complets, noms de `<cmath>` ou de Boost).

5. Défauts à corriger sans attendre la refonte (dlist, progress) ou avec elle.
6. Recouvrements : type matriciel commun, arithmétique étendue partagée et son nom.

Ensuite : bibliothèque par bibliothèque, une entrée de rupture au CHANGELOG chacune.
