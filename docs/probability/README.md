# `wsl::probability` (`proba`)

`#include <wsl/probability>` — les lois de probabilité comme **objets mathématiques** :
paramètres nommés, opérations en membres, tirage `X(g)`.

```cpp
namespace proba = wsl::probability;
proba::laws::normal const N({.mean = 1.0, .stddev = 2.0});
N.cdf(3.0); N.quantile(0.975); N.pdf(0.5); N.mean(); N.variance();
double x = N(g);                  // tirage approché (le défaut : ziggourat)
double y = N(g, wsl::exact);      // tirage exact
std::vector<double> v = N(1000, g);
```

## Lois (`proba::laws`)

- **Continues** : `uniform`, `exponential`, `normal`, `gamma`, `beta`, `chi_squared`,
  `student_t`, `lognormal`, `cauchy`, `fisher_f`, `weibull`, `laplace`, `logistic`, `pareto`,
  `gumbel`, `frechet`, `rayleigh`, `triangular`, `arcsine`, `erlang`, `inverse_gaussian`,
  `rice`, `skew_normal`, `stable`.
- **Discrètes** : `discrete_uniform`, `bernoulli`, `geometric`, `poisson`, `binomial`,
  `negative_binomial`, `hypergeometric`, `beta_binomial`, `logarithmic`,
  `discrete_gaussian`, `skellam`, `zipf`, `categorical`.
- **Multivariées** : `dirichlet`, `multinomial`, `multivariate_normal`,
  `multivariate_student_t`, `wishart`, `inverse_wishart`, `uniform_sphere`, `uniform_ball`,
  `von_mises`, `von_mises_fisher` (dimension fixe ou dynamique).
- **Combinatoires** (tirage exact) : permutations (`uniform_permutation`, `ewens`,
  `mallows`), parties et mots, partitions, compositions, graphes (`erdos_renyi`,
  `barabasi_albert`, `random_regular_graph`…), arbres ; matrices aléatoires
  (`haar_orthogonal`, `haar_unitary`, `goe`, `gue`).

Moments, densité ou masse, répartition, survie, quantile et quantile de survie : erreur
mesurée sous 0,51 ulp, contre des références mpmath.

## Tirages et modes

`X(g)` prend la variante approchée de la loi si elle en a une, son seul tirage sinon ;
`X(g, wsl::exact)` demande le tirage exact (arrondi correct d'un tirage réel exact, nombre de
bits consommés fixé), **refusé à la compilation** pour une loi qui n'en a pas
(`L::sampling` le dit). Les lois sont portables au bit près dans tous les modes (décision 42 :
leur version déterministe viendra plus tard).

## Algèbre des lois et variables aléatoires

Aux noms de la théorie de la mesure : `affine`, `convolution` (qui reconnaît les familles :
Poisson + Poisson = Poisson), `truncate`, `mixture`, `product_measure`, `pushforward`.
Variables aléatoires sur un espace probabilisé commun : `proba::variable X{law}`, expressions
`2 * X + Y` évaluées paresseusement, `proba::evaluate(g, X, Z)` sur le même ω, espérance et
covariance par linéarité ([design/variables.md](../../design/variables.md)).

Sous-domaine : [processus stochastiques](stochastic_processes/README.md). Conception :
[design/random.md](../../design/random.md). Exemple :
[`examples/probability/laws.cpp`](../../examples/probability/laws.cpp).
