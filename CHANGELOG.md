# Journal des versions

Numérotation `MAJEUR.MINEUR.CORRECTIF`. Tant que `MAJEUR` vaut 0, une version
mineure peut casser l'interface ; un correctif, jamais. Chaque entrée est rangée
par bibliothèque, pour savoir quoi vérifier en montant de version dans un projet.

## [Non publié]

### Toutes les bibliothèques
- Rupture : C++23 exigé partout (dlist, graded, graded_io, elo, special, permutation et
  random descendaient jusqu'à C++20) ; tests compilés une seule fois, en C++23. Décision de
  la revue des interfaces (docs/interfaces.md).

### mdarray

- Nouvelle bibliothèque (conception : docs/mdarray.md), lot 1 en cours. Positions
  `mdposition<N, IndexType = std::size_t>` et déplacements `mddisplacement<N, DiffType =
  std::ptrdiff_t>` : modèle affine (position − position = déplacement, position ±
  déplacement = position ; `*`, `/`, `%` par un entier sur les déplacements seulement,
  troncature vers zéro), construction à exactement N composantes de tout type entier,
  déduction, conversions depuis `std::array` et `std::span`, décomposition, `first`, `last`,
  `concat`, `filled`, `unit`, ordre lexicographique (`<=>`) et composante par composante
  (`all_less`, `all_less_equal`), normes (`norm1`, `norm_inf`, `norm2_squared`),
  `std::hash`, `std::formatter` (« (2, 1) ») et `operator<<`. Arithmétique entière :
  `div_euclid`, `rem_euclid`. Préconditions vérifiées en débogage (`WSL_ASSERTIONS`) et
  toujours en évaluation constante.

### elo

- Rupture : espace `wsl::elo` (noms trop génériques dans `wsl::`) ; `elo_from_score` →
  `from_score`, `score_from_elo` → `to_score`, `signed_elo` → `format_signed`,
  `estimate_elo` → `pairwise`, une seule structure `rating {elo, low, high}` (au lieu de
  `elo_estimate` et `rating`). Bilans incohérents et `samples < 40` : `std::invalid_argument`
  (au lieu d'`assert`, sans effet en Release). Rééchantillonnage par `<random>` inchangé.

### array2d

- Rupture : `extent(i)` et `rank()` (comme `std::mdspan`) remplacent `row_count()` et
  `column_count()` (décision 6 : conventions communes avec `square_matrix`). Commentaire
  d'en-tête (exemple vérifié), message d'exception préfixé `wsl::array2d`, espaces après
  `if`/`for`.

### graded, graded_io

- Rupture (interne) : aides internes dans `wsl::detail::graded` (au lieu de `wsl::detail`,
  partagé). Gardes `#pragma once` ; commentaires accentués.

### permutation (nouvelle)
- `wsl::permutation` : construction vérifiée (images, cycles), composition `s * t` (s∘t),
  inverse, puissances, cycles canoniques, type de cycles, points fixes, signature,
  inversions en O(n log n), rang lexicographique (n <= 20), `permute`, `constexpr`.

- Style : espaces après `if`/`for` (aligné sur le dépôt).

### dlist

- Rupture (mineure) : concept `wsl::streamable` déplacé dans `wsl::detail::dlist` ;
  messages d'exception en français. Commentaire d'en-tête (exemple vérifié), style des
  membres (`x_`) et des paramètres de modèle (`is_const`) alignés sur le reste du dépôt.
- Correction : `erase(const_iterator)` et `erase(first, last)` ne compilaient pas.
- `position(it)` : position absolue d'un élément (jusqu'ici lisible seulement par `dump`).
- Recherche de nœud, effacement et `extract` en un seul parcours (O(n) au lieu de O(n²)).

### progress
- Style : espaces après `if`/`for` (aligné sur le dépôt).
- Correction : le calcul des tranches de 5 % débordait pour `done` > 2⁵⁹.

### random (nouvelle)
- Rupture : tout wsl::random passe dans l'espace de noms `wsl::random` (lois, processus,
  objets combinatoires, moteurs, primitives : `wsl::random::normal`, `wsl::random::xoshiro256pp`…),
  ses internes dans `wsl::random::detail` (plus de `wsl::detail` partagé avec graded). Les
  appels non qualifiés (`pdf(X, x)`, `sample(X, g)`) restent trouvés par ADL. Alias conseillé :
  `namespace wr = wsl::random;`.
- Objets-fonctions `wsl::random::pdf`, `log_pdf`, `pmf`, `log_pmf`, `density`, `log_density`,
  `cdf`, `ccdf`, `quantile`, `cquantile`, `sample`, `mean`, `variance`, `skewness`,
  `excess_kurtosis`, `entropy`, `median`, `mode`, `support` : appel qualifié et passage en
  argument (`std::bind_front(wr::pdf, X)`) ; ils prennent le membre de la loi, sinon la
  fonction trouvée par ADL. Densité unique (dérivée de Radon–Nikodym) et mesures de
  référence `lebesgue`, `counting`, `surface` (`reference_measure_t<L>`). Vue paresseuse
  `samples(X, g)`. Rupture : `sample(X, g, fast)` non qualifié pour une loi sans variante
  rapide ne compile plus (`wr::sample(X, g, wr::fast)`).
- Rupture : opérations des lois en membres (`X.pdf(x)`, `X.cdf(x)`, `X.quantile(p)`,
  `X.sample(g)`…) au lieu de fonctions amies : `pdf(X, x)` non qualifié ne compile plus,
  `wsl::random::pdf(X, x)` (objet-fonction) le remplace. Toutes les lois héritent de
  `law_interface` : `X(g)`, `X(g, wr::fast)`, `X(n, g)` (vecteur de n tirages, mêmes valeurs
  et même consommation que n appels `X(g)`), `X.density(x)`, `X.log_density(x)`.
- Rupture : moteur en dernier argument pour les primitives, comme dans la norme :
  `uniform_below(n, g)`, `uniform_int(a, b, g)`, `bernoulli_trial(p, g)`,
  `bernoulli_trial(num, den, g)` (`random_bits(g)`, `uniform_real<T>(g)` et `shuffle(range, g)`
  inchangés). Mêmes tirages.
- Rupture : plus d'opérateurs arithmétiques sur les lois (réservés aux futures variables
  aléatoires) ; noms de la théorie de la mesure : `affine(P, {.scale, .shift})` remplace
  `a * X + b`, `-X`, `X / a` ; `convolution(P, Q)` remplace `X + Y` (type
  `convolution_law`, ex-`independent_sum` ; différence par `convolution(P, affine(Q,
  {.scale = -1}))`) ; `product_measure(P, Q, …)` remplace `independent` ;
  `truncate(P, {.lower, .upper})`, bornes nommées et facultatives, remplace
  `truncate(X, a, b)`. En-têtes `convolution.hpp` et `product_measure.hpp` (ex-`sum.hpp`,
  `independent.hpp`). Mêmes résultats.
- Rupture : vocabulaire des paramètres (docs/interfaces.md, décision 2b). Renommés :
  normale `.sd` → `.stddev` ; log-normale `.mu`, `.sigma` → `.log_mean`, `.log_stddev` ;
  gaussienne discrète `.mu`, `.sigma` → `.location`, `.scale` ; von Mises `.mu`, `.kappa` →
  `.location`, `.concentration` (von Mises–Fisher : `.concentration`) ; uniformes, arc
  sinus, triangulaire `.low`, `.high` → `.lower`, `.upper` ; composition et partition
  d'entier `.n` → `.total` ; arbre binaire `.nodes` → `.vertices`. Accesseurs renommés de
  même (`stddev()`, `lower()`, `upper()`, `concentration()`… ; hypergéométrique `lower()`,
  `upper()`). Ordre des champs position → échelle → formes : gamma `{.rate | .scale,
  .shape}`, Erlang `{.rate, .shape}`, Weibull `{.scale, .shape}`, Fréchet `{.location,
  .scale, .shape}`, stable `{.location, .scale, .alpha, .beta}`, Student multivariée
  `{.location, .scale, .degrees_of_freedom}`, Wishart `{.scale, .degrees_of_freedom}`.
  Valeurs par défaut position 0, échelle 1 (Cauchy, Gumbel, Laplace, logistique, Rayleigh,
  Weibull, Pareto, Fréchet, stable, gaussienne discrète, von Mises, uniforme sur [0, 1]) ;
  constructeurs par défaut pour les lois standard (`normal()`, `exponential()`, `uniform()`,
  `cauchy()`…). Mêmes tirages et mêmes valeurs.
- `param()` et `param_type` sur les 61 lois à structure de paramètres : `L(X.param()) == X`
  (gamma, à deux paramétrisations, attend la forme canonique). Loi bêta : pas d'accesseur
  `beta()` (nom de la classe), `B.param().alpha` et `B.param().beta`.
- Correction : `operator==` de `geometric` était implicitement supprimé (membre
  double-double sans égalité) ; il compare désormais `p`.
- Constantes indéfinies (décision 3b) : +∞ quand l'intégrale diverge vers +∞, NaN quand elle
  n'a pas de sens. Corrigés : moyenne de Fisher pour d₂ <= 2 (+∞, NaN avant) ; asymétrie et
  kurtosis +∞ quand la variance est finie et le moment d'ordre 3 ou 4 infini (Pareto,
  Fisher ; NaN avant) ; variance NaN quand la moyenne n'est pas finie (Pareto et Fréchet
  α <= 1, stable α <= 1 ; +∞ avant) ; moyenne stable ±∞ pour α < 1, β = ±1 ; moyenne de
  Wishart inverse pour n <= p + 1 : +∞ sur la diagonale, NaN ailleurs ; asymétrie et
  kurtosis NaN pour les lois dégénérées (binomiale p ∈ {0, 1} ou n = 0, géométrique et
  binomiale négative p = 1, uniforme discrète a = b ; ±∞ ou 0 avant).
- Types des constantes (décision 3a) : entropie en T pour von Mises et les lois continues
  multivariées (double avant) ; `product_measure` : type commun des composantes. Médianes
  `noexcept` (géométrique, inverse gaussienne, stable).
- Moteurs complets au sens de *RandomNumberEngine* : `E()` (non explicite), `E(q)` et
  `seed(q)` depuis une *SeedSequence* (mots de 32 bits assemblés deux à deux, poids faible
  d'abord, comme `mt19937_64` ; xoshiro : état nul remplacé par celui de `seed(0)`),
  `seed()`, `os << e` et `is >> e` (décimal ; lecture invalide : moteur intact, `failbit`).
  `philox4x64::discard` en O(1) (compteur + n/4 ; il avançait bloc par bloc).
- Rupture : `square_matrix` suit les conventions d'`array2d` (décision 6) : accès `m[i, j]`
  (au lieu de `m(i, j)`), dimensions `extent(i)` et `rank()` ; `size()`, qui rendait la
  dimension n et non le nombre d'éléments, est retiré (plutôt que de changer de sens).
  Dans `assert`, `m[i, j]` demande des parenthèses doubles (virgule de macro).
- Variables aléatoires, lot 1a (docs/variables.md, décisions V1 à V6) : `wr::variable X{P}`
  (identité unique, copies partagées, `independent_copy()`), expressions `+ − × ÷` avec
  variables et scalaires, évaluation couplée (`(X − X)(g) == 0` ; chaque variable de base
  tirée une fois par ω, par identifiant croissant), `wr::evaluate(g, X, Y, …)` sur un même
  ω, `wr::E`, `wr::Var`, `wr::Cov` par linéarité (coefficients regroupés : Var(X + X) =
  4 Var X).
- Variables aléatoires, lot 1b (V3) : une combinaison affine Σ aᵢ Xᵢ + b de variables à
  tirage paresseux (normale, exponentielle, Laplace, uniforme) est l'arrondi correct de sa
  valeur exacte, coefficients regroupés par identifiant, divisions par des constantes
  exactes (`(X + X) / 2`, `2 * X − X`, `(3 * X + 1 − 1) / 3` rendent X au bit près) ; réels
  paresseux partagés par toutes les expressions d'un même ω (`lazy_sum::add_shared`).
  Ailleurs (produits, autres lois, coefficient non représentable), calcul flottant sur des
  valeurs tirées exactement ; `Z(g, wr::fast)` : calcul flottant sur tirages rapides.
- Variables aléatoires, lot 1c (V4) : `wr::law(X)` (la loi de X) et `wr::law(Z)`
  (`expression_law`) : tirage de Z, moyenne, variance, asymétrie et kurtosis par les
  cumulants ; `family()` : loi exacte pour une famille close (normale, Cauchy ; gamma de
  même intensité à coefficients positifs), résolue à l'exécution ; ses fonctions sont
  celles de l'expression, NaN sans identité (règle 3b).
- Correction (lot 1b) : le produit de deux variables et les constantes (E, Var, Cov,
  cumulants) passent par `fmul` : sous `-ffp-contract=fast`, la contraction en fma changeait
  les valeurs (CI FMA).
- Variables aléatoires, lot 1d : `wr::apply(f, X, …)` ; fonctions `exp`, `log`, `sqrt`,
  `cbrt`, `expm1`, `log1p`, `atan`, `erf` (bijections monotones déclarées), `sin`, `cos`,
  `tan`, `abs`, `pow`, `hypot`, `atan2`, `min`, `max` (wsl::special) ; `wr::law(f(X))` :
  mesure image de la loi de X (`pushforward`).
- Variables aléatoires, lot 2 : événements (comparaisons, `&&`, `||`, `!`), décidés
  exactement pour les comparaisons affines ; `wr::P(A)` exacte (seuils, comparaisons de loi
  connue, événements indépendants), NaN ailleurs ; `wr::given(Z, A)` (loi conditionnelle, par
  rejet ; loi tronquée ou loi de Z quand c'est une identité).
- Skellam : `skellam({.mean1, .mean2})` (différence de deux Poisson indépendantes) : masse,
  log-masse, répartition, survie, asymétrie, kurtosis sous 0,51 ulp (masse par
  `log_cyl_bessel_i`, corrigée au premier ordre de l'arrondi de z = 2√(μ₁μ₂) ; queues sommées
  en double-double), quantiles exacts, tirage N₁ − N₂.
- Rice : `rice({.scale = 1, .noncentrality})` : densité, log-densité, répartition, survie
  (séries de Marcum à termes positifs), quantiles, moyenne, variance sous 0,51 ulp ;
  ν = 0 délégué à `rayleigh` ; tirage hypot(ν + σZ₁, σZ₂).
- Correction : survies sous-normales (erfc(y)/2 pour 27 <= y < 27,25) rendues nulles, pour
  `stable` (α = 2) et `skew_normal`.
- Correction : quantiles de `inverse_gaussian` faux de quelques ulp dans les queues des lois
  très concentrées (λ/μ grand) : finition entre deux doubles voisins.
- Références de précision des lois : une famille par module de `tools/random/references/`,
  engendrée séparément par CMake ; seules les familles modifiées sont recalculées.
- Zipf (finie) : `zipf({.size, .exponent})` sur {1, …, N} : masse, log-masse, répartition,
  survie, moyenne, variance sous 0,51 ulp (Euler–Maclaurin en double-double au-delà de
  max(64, 8|p|)), quantiles exacts, tirage par rejet par inversion (Hörmann et Derflinger).
- Correction : quantile de `geometric` pour u sous-normal (t = log(1 − u)/log(1 − p) sous-dépassait,
  ⌈0⌉ − 1 rendait le plus grand entier au lieu de 0).
- Variables aléatoires, lot 3 : composantes `V[i]` (couplées), E, Var, Cov avec la matrice de
  covariance, famille gaussienne (normales et composantes de vecteurs gaussiens) ;
  `multivariate_normal` et la grille de `brownian_motion` se déclarent gaussiennes, la grille
  gagne `mean()` et `covariance()`.
- Documentation (docs/random.md) réalignée sur l'interface : §1 à §7 décrivent l'interface
  actuelle (concept de loi, uniformes, exactitude, paramètres exacts), exemples et noms de
  paramètres des tableaux à jour, §14 sur l'interface réelle des processus (exemples
  compilés), état des points du §15 ; README à jour.
- Intensités rationnelles exactes (décision 5 révisée) : gamma, exponentielle et Erlang
  stockent λ = P/Q·2ᵉ réduite ; `.rate` ou `.scale` (double, ou `wr::ratio(a, b)`), exactes
  toutes deux, aucune des deux : λ = 1. `==` est l'égalité en loi
  (`gamma({.rate = 2})` == `gamma({.scale = 0.5})`), `param()` rend la forme lisible
  (`.rate` ou `.scale` en double si possible, sinon `ratio`). Tirage exact pour toute
  intensité rationnelle (V·Q·2⁻ᵉ/P). Rupture : `gamma::shape_rate` et `shape_scale`
  remplacés par `gamma::parameters` ; `exponential({.scale})` ajouté ; `c·X` d'une gamma
  ou d'une exponentielle garde λ/c exact (arrondi avant). Voies rapides inchangées quand
  λ ou 1/λ est un double : mêmes valeurs et mêmes tirages.
- Normale à variance exacte (décision 5, N4) : `{.mean, .stddev}` ou `{.mean, .variance}`
  (aucun : σ = 1) ; forme canonique la variance exacte (σ² en double-double), `==` égalité en
  loi (`{.variance = 4}` == `{.stddev = 2}`). Si √v n'est pas un double : fonctions avec σ = √v
  en double-double, tirage exact de μ ± √v·V par comparaison des carrés (sans calculer √v),
  avec voie rapide certaine (176 ns contre 136 ns pour σ double). Identités de l'algèbre :
  variances additionnées (convolution) et multipliées par a² (image affine), à 2⁻¹⁰⁶ près
  (σ arrondis avant). Sommes exactes (marche aléatoire, Poisson composé) : un groupe de
  termes c √v V par somme, décidé par comparaison de P² et v B² ; une somme mêlant deux
  variances à racines irrationnelles différentes lève `std::domain_error`. Voie rapide
  inchangée quand σ est un double : mêmes valeurs et mêmes tirages.
- Moteurs `splitmix64`, `xoshiro256pp` (sauts de 2^128 et 2^192), `pcg64` (PCG XSL RR
  128/64, avance en O(log n)), `philox4x64` (à compteur) : `std::uniform_random_bit_generator`,
  `constexpr`, sorties vérifiées contre NumPy et, pour les sauts de xoshiro, par
  matrices sur GF(2).
- `random_bits`, `uniform_below` (Lemire), `uniform_int`, `uniform_real` (exact : chaque
  flottant de [0, 1) avec exactement sa probabilité, sous-normaux compris),
  `uniform_real_fast`, `bernoulli_trial` (p flottant ou rationnel, exacte), `shuffle`
  (Fisher–Yates), pour tout générateur standard.

- Concept de loi (`probability_law`, `has_mean`, `has_cdf`…), tirage `sample(X, g)` et
  `sample(X, g, wsl::fast)`, adaptateur `as_distribution` vers `<random>`.
- Lois `uniform<T>`, `discrete_uniform<I>`, `bernoulli` : constantes, densité ou masse,
  répartition, survie, quantile (exact pour la loi discrète), tirage.
- La primitive `bernoulli(g, p)` devient `bernoulli_trial(g, p)`.
- Lois `exponential<T>` (`{.rate}`) et `normal<T>` (`{.mean, .sd}`) : tirage **exact**
  (valeur tirée exactement selon la loi puis arrondie au plus proche, réels paresseux),
  tirage rapide `wsl::fast` (inversion ; ziggourat), densité, log-densité, répartition,
  survie, quantile mesurés sous 0,502 ulp contre mpmath.
- Lois `geometric<I>` (`{.p}`, échecs avant le premier succès), `poisson<I>` (`{.mean}`),
  `binomial<I>` (`{.trials, .p}`) : masses et log-masses sous 0,501 ulp (décomposition de
  Loader), répartition, survie et quantile exact de la géométrique ; tirages bit-exacts
  (inversion, PTRS, BTRS).
- Lois `gamma<T>` (`{.shape, .rate}` ou `{.shape, .scale}`) et `beta<T>` (`{.alpha, .beta}`) :
  densités et log-densités sous 0,501 ulp, moments, tirages bit-exacts (Marsaglia–Tsang).
- Exponentielle à argument double-double et erreur de Stirling plus précises : quelques
  résultats des fonctions des lois changent d'un ulp (empreintes des fonctions mises à jour).
- Répartitions, survies et quantiles de `poisson`, `binomial` (quantiles exacts), `gamma`,
  `beta` ; entropies de `gamma` et `beta` : tout sous 0,5 ulp, arrondi correct à 100 % des
  points mesurés.
- Lois `chi_squared`, `student_t`, `lognormal`, `cauchy` : densité, log-densité, répartition,
  survie, quantile, moments (là où ils existent), entropie, tirage ; tout sous 0,5 ulp.
- Lois `fisher_f`, `weibull`, `laplace` (tirage exact), `logistic`, `pareto` : densité,
  log-densité, répartition, survie, quantile, moments, entropie, tirage ; tout sous 0,5 ulp.
- `scale_by_power_of_two` (wsl::special, détail) accepte tout exposant.
- `cquantile(X, q)` (quantile de survie, sans former 1 − q) pour les quatorze lois
  continues, sous 0,5 ulp ; concept `has_cquantile`.
- Algèbre des lois, lot 1 : image affine (`a * X + b`, `-X`, `X / a`, `wsl::affine`), les
  familles fermées gardant leur type (normale, Cauchy, Laplace, logistique, uniforme) ou
  leurs fonctions (exponentielle, gamma, Weibull, Pareto, log-normale) ;
  `wsl::truncate(X, a, b)` ; `wsl::mixture(wsl::weighted(w, X)…)`, tirage exact si les
  composantes le sont. Erreurs de composition mesurées (docs/random.md §10).
- Algèbre des lois, lot 2 : sommes indépendantes `X + Y`, `X - Y` (identités normale,
  Cauchy, Poisson, χ² ; gamma de même intensité, exponentielles, binomiales de même p ;
  cumulants ; convolution exacte des lois discrètes), loi jointe `wsl::independent`,
  image `wsl::pushforward(X, f)` (bijections monotones `wsl::increasing`,
  `wsl::decreasing`) ; `gamma::has_rate`.
- Lois multivariées, lot 1 : `dirichlet`, `multinomial`, `multivariate_normal`,
  `uniform_sphere`, `uniform_ball`, en dimension fixe ou `wsl::dynamic` ;
  `wsl::square_matrix`. Densités et masses en logarithme double-double, moments,
  covariances, entropies, marginales.
- Objets combinatoires, lot 1, tirages exacts : `uniform_permutation`, `ewens`,
  `mallows` (valeurs `wsl::permutation`), `uniform_subset`, `bernoulli_subset`,
  `uniform_combination`, `uniform_word`, `uniform_integer_partition` ; masses en
  logarithme double-double. Choix exacts à poids flottants (`detail::choose`,
  `detail::fraction`). `wsl::random` dépend désormais de `wsl::permutation`.
- Le support devient facultatif dans le concept de loi (`has_support`).
- Processus, lot 2 (temps continu, exacts point par point) : `detail::lazy_sum` (sommes
  de réels paresseux arrondies correctement), `brownian_motion`, `brownian_bridge`,
  `ornstein_uhlenbeck`, `geometric_brownian_motion`, `poisson_process`,
  `inhomogeneous_poisson_process`, `compound_poisson`, `continuous_markov_chain` ;
  `random_walk` à pas réels ; variantes `wsl::fast`.
- Processus, lot 1 (temps discret, tirages exacts) : vue paresseuse `trajectory(P, g)`
  (`std::ranges`, valeurs tirées à la lecture), concept `stochastic_process` ;
  `categorical` (alias de Walker–Vose exact), `rademacher`, `random_walk`,
  `markov_chain` (`at`, `stationary`, `steps`), `polya_urn`, `dirichlet_multinomial`,
  `chinese_restaurant` (Pitman–Yor), `ewens_partition`.
- Objets combinatoires, lot 2, tirages exacts : `erdos_renyi`, `uniform_graph`,
  `barabasi_albert`, `random_regular_graph`, `uniform_tree` (Cayley, Prüfer),
  `uniform_binary_tree` (Rémy), `uniform_set_partition`, `uniform_composition` ;
  `wsl::edge_list`, `wsl::binary_tree`.
- Lois `von_mises` (cercle) et `von_mises_fisher` (sphère, dimension fixe ou dynamique) :
  densités, longueurs résultantes, entropies sous 0,5 ulp ; tirages de Best–Fisher et Wood.
- Lois multivariées, lot 2 : `multivariate_student_t`, `wishart`, `inverse_wishart`,
  `haar_orthogonal`, `haar_unitary`, `goe`, `gue` (valeurs `wsl::square_matrix`, réelles
  ou `std::complex`, produits complexes protégés de la contraction).
- `ccdf`, `pdf`, `log_pdf` de Cauchy justes au-delà de |z| = 2⁶⁰⁰ (rendaient 0 ou −∞) ;
  log-densités de Laplace et de la logistique de même.
- Variance de la normale et du χ², tirage gamma par l'échelle : produits rendus protégés
  de la contraction (GCC les fusionnait avec l'addition de l'appelant ; mêmes bits sans
  contraction).
- `exp_scaled` (wsl::special, détail) rend une paire normalisée (précision des
  exponentielles enchaînées).
- Compléments au catalogue : `gumbel<T>` (`{.location, .scale}`), `frechet<T>`
  (`{.shape, .scale, .location}`), `rayleigh<T>` (`{.scale}`), `triangular<T>`
  (`{.low, .mode, .high}`), `arcsine<T>` (`{.low, .high}`) : densité, log-densité,
  répartition, survie, quantile et quantile de survie sous 0,5 ulp contre mpmath,
  moments et entropie ; tirages bit-exacts. `erlang<T>` (`{.shape}` entier, `{.rate}`) :
  fonctions de la gamma, tirage **exact** (somme paresseuse de k exponentielles exactes).
- `inverse_gaussian<T>` (`{.mean, .shape}`) : densité, log-densité, répartition, survie
  (rapport de Mills à 2⁻⁸⁵, série de Taylor quand la survie se compense), quantiles par
  Newton à résidu précis, moments ; tirage bit-exact (Michael, Schucany et Haas).
- `von_mises` : répartition, survie, quantile et quantile de survie (intégrale de la queue
  par Gauss–Legendre par morceaux, quantiles par Newton sur la plus petite part).
- `negative_binomial<I>` (`{.successes, .p}`, r réel) et `hypergeometric<I>`
  (`{.population, .successes, .draws}`) : masses, répartitions sous 0,5 ulp, quantiles
  exacts, moments ; tirages bit-exacts.
- `beta_binomial<I>` (`{.trials, .alpha, .beta}`), `logarithmic<I>` (`{.p}`),
  `discrete_gaussian<I>` (`{.mu, .sigma}`) : masses, répartitions et moments sous 0,5 ulp,
  quantiles exacts ; tirages bit-exacts.
- `skew_normal<T>` (`{.location, .scale, .shape}`) : densité, répartition, survie, quantiles
  et moments sous 0,5 ulp (queue courte intégrée directement) ; tirage bit-exact.
- `stable<T>` (`{.alpha, .beta, .scale, .location}`, paramétrisation S1) : densité,
  répartition et survie sous 0,5 ulp (intégrales de Nolan en double-double, toujours en
  sommes positives), quantiles par Newton ; tirage de Chambers–Mallows–Stuck, bit-exact.
- Interne : aides regroupées dans `detail/arithmetic.hpp` (compléments à l'échelle et en
  double-double), `detail/trigonometry.hpp`, `detail/quadrature.hpp` (Gauss–Legendre) et
  `detail/gaussian.hpp` ; doublons supprimés (trois opposés, deux sommes et deux différences
  à l'échelle, trois `exp_of`…) ; stable.hpp n'inclut plus skew_normal.hpp ni von_mises.hpp.
  Aucun résultat ne change (empreintes identiques).
- Interne : un seul solveur de quantiles (`detail/root.hpp` : Newton sur log de la plus
  petite part, encadré, deux pas à résidu précis) pour inverse_gaussian, von_mises,
  skew_normal et stable, et une seule marche d'intégration par morceaux
  (`march_integral`, detail/quadrature.hpp) pour von_mises, skew_normal et stable ; environ
  210 lignes de moins. Précision inchangée (sous 0,5 ulp), empreintes identiques.
- Variance de Weibull et lgamma(1 + ε) mis en commun (`detail/gamma_moments.hpp`, ε
  négatif accepté pour Fréchet) : mêmes bits.

### exact_random (retirée)
- Remplacée par `wsl::random` : `uniform_below(g, n)` (exact, tout générateur, plus
  rapide), `shuffle(range, g)`, et les moteurs initialisés à volonté. Rupture : l'ordre
  des arguments et les suites tirées changent.

### cpp-init (outil)
- Preset release des nouveaux projets en `-march=native` sur x86-64 (fma, AVX2…).

### special (nouvelle)
- Arguments entiers ou mêlés (entiers, `float`, `double`) pour toutes les fonctions à
  plusieurs arguments (`hypot`, `atan2`, `beta`, `lbeta`, `ibeta`, `ibetac`, `ibeta_inv`,
  `ibetac_inv`, `gamma_p`, `gamma_q` et inverses, `cyl_bessel_i`) et entiers pour `erf_inv`,
  `erfc_inv`, calculés en `double` comme `<cmath>` ; `long double` supprimé partout
  (`digamma` et ces fonctions compris). Règle de nommage écrite (special.md §3).
- `cyl_bessel_i(ν, x)` (Bessel modifiée I_ν, float et double) : série, Hankel, Debye
  (polynômes engendrés exactement) ; 0,4993 ulp mesurée, arrondi correct partout.
- `wsl::special::exp` et `wsl::special::log`, en `float` et en `double` : même
  résultat au bit près sur toute plateforme, quelles que soient les options de
  compilation du projet (fma, `-ffp-contract=fast`), et `constexpr`. Erreur
  maximale mesurée 0,504 ulp (`exp` double), 0,500 ulp (les autres).
- `log` en double environ 30 % plus rapide qu'à sa première écriture (même précision).
- `exp2`, `expm1`, `log2`, `log10`, `log1p` (float, double), mêmes garanties ; erreur
  maximale mesurée de 0,499 à 0,502 ulp. `exp2(n)`, `log2(2ⁿ)` et `log10(10ⁿ)` exacts.
- `pow` (float, double, arguments mixtes comme `<cmath>`), mêmes garanties ; erreur
  maximale mesurée 0,502 ulp ; cas à part identiques à `<cmath>` (Annexe F de C).
- `sin`, `cos`, `tan`, `atan`, `atan2` (float, double), mêmes garanties, arguments
  immenses compris (réduction de Payne–Hanek exacte) ; erreur maximale mesurée 0,500 ulp,
  100 % d'arrondis corrects sur les échantillons. Encore 3,5 à 5 fois plus lents que la
  libm d'Apple.
- `sqrt` (arrondi correct, évaluable à la compilation), `cbrt`, `hypot` (float,
  double) ; erreur maximale mesurée 0,500 ulp, 100 % d'arrondis corrects.
- `erf`, `erfc` (float, double), mêmes garanties, `erfc` sous-normaux compris ; erreur
  maximale mesurée 0,500 ulp (`erf`), 0,503 ulp (`erfc`).
- `lgamma`, `tgamma` (float, double), mêmes garanties ; erreur maximale mesurée
  0,500 ulp (`lgamma`, x > 0) et 0,507 ulp (`tgamma`). Pour x < 0 près des zéros de
  lgamma, erreur bornée à 2⁻⁵⁸ en valeur absolue seulement (voir docs/special.md).
- `erf_inv`, `erfc_inv` (float, double), mêmes garanties, `erfc_inv` jusqu'à 2⁻¹⁰⁷⁴ ;
  erreur maximale mesurée 0,502 et 0,506 ulp. `erf` et `erfc` affinées pour leur servir
  de résidu (mêmes précisions mesurées, bits changés pour quelques arguments).
- Produits exacts par fma quand la cible en a une (ARM64, x86-64 avec `-mfma` ou
  `-march=x86-64-v3`…) : 10 à 28 % plus rapide selon la fonction, **mêmes bits** que sans
  fma. `WSL_SPECIAL_NO_FMA` pour s'en passer.
- `digamma` (float, double) : erreur maximale mesurée 0,500 ulp pour x > 0 ; pour x < 0,
  0,51 ulp ou 2⁻⁶² en valeur absolue près des zéros (réflexion).
- `gamma_p`, `gamma_q` (gamma incomplètes régularisées, float et double) : erreur maximale
  mesurée 0,500 ulp, arrondi correct sur 100 % des points ; chacune précise même
  minuscule (Temme pour a >= 50 près de x = a, série, fraction continue, série en a pour
  les petits a). Domaine complet : a et x jusqu'à DBL_MAX (P + Q = 1 et croissance
  vérifiées au-delà des références).
- `gamma_p_inv`, `gamma_q_inv` (float, double) : erreur maximale mesurée 0,500 ulp,
  arrondi correct sur 100 % des points (Newton encadré en double-double sur le côté
  minoritaire ; queue basse des petits a en logarithme).
- `beta`, `lbeta`, `ibeta`, `ibetac` (bêta incomplète régularisée et complément, float et
  double) : erreur maximale mesurée 0,500 ulp, arrondi correct sur 100 % des points, chacune
  précise même minuscule (formule propre aux petits paramètres).
- `ibeta_inv`, `ibetac_inv` (float, double) : erreur maximale mesurée 0,500 ulp, arrondi
  correct sur 100 % des points. Phase 4 terminée.
- Garantie bit-exacte renforcée face à GCC : les produits rendus (`return y · 2^e`) et les
  `x + x` des chemins NaN ne peuvent plus être fusionnés avec une addition de l'appelant
  après intégration (`1 + exp(x)`…). Aucun résultat ne change ; le contrôle anti-fma
  essaie désormais `1 + f(x)` pour chaque fonction.
- Constantes double-double internes (`detail/extended.hpp`) pour les lois de Gumbel et de
  Rayleigh (√(π/2), (4 − π)/2, √(2 log 2), asymétries), vérifiées contre mpmath.

## [0.1.0] — 2026-11-16

Première version : regroupement des bibliothèques existantes.

### array2d
- Reprise de cpp-array2d, inchangée.

### dlist
- Reprise de cpp-container-distance-list.
- Retrait des alias `offset` et `position`, qui renvoyaient à un `index`
  supprimé (remplacé par `graded` / `index_jet`) et empêchaient la compilation.

### graded, graded_io
- Reprise de cpp-container-distance-list ; `graded-io.hpp` devient `graded_io.hpp`.
- `operator<=>` entre arités ou types différents passe par la comparaison
  synthétisée de la norme (comme `std::pair`) : `jet<int> <=> jet<unsigned>`
  compile désormais et rend un `std::weak_ordering`.

### exact_random
- Reprise de cpp-sudoku, rendue autonome (type `wsl::int128`) et en en-tête seul.

### elo, progress
- Reprise de cpp-traversal, namespace `wsl`.

[Non publié]: https://github.com/wills31415/cpp-libraries/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/wills31415/cpp-libraries/releases/tag/v0.1.0
