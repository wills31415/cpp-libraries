# Journal des versions

Numérotation `MAJEUR.MINEUR.CORRECTIF`. Tant que `MAJEUR` vaut 0, une version
mineure peut casser l'interface ; un correctif, jamais. Chaque entrée est rangée
par bibliothèque, pour savoir quoi vérifier en montant de version dans un projet.

## [Non publié]

### random (nouvelle)
- Moteurs `splitmix64`, `xoshiro256pp` (sauts de 2^128 et 2^192), `pcg64` (PCG XSL RR
  128/64, avance en O(log n)), `philox4x64` (à compteur) : `std::uniform_random_bit_generator`,
  `constexpr`, sorties vérifiées contre NumPy et, pour les sauts de xoshiro, par
  matrices sur GF(2).
- `random_bits`, `uniform_below` (Lemire), `uniform_int`, `uniform_real` (exact : chaque
  flottant de [0, 1) avec exactement sa probabilité, sous-normaux compris),
  `uniform_real_fast`, `bernoulli` (p flottant ou rationnel, exacte), `shuffle`
  (Fisher–Yates), pour tout générateur standard.

### exact_random (retirée)
- Remplacée par `wsl::random` : `uniform_below(g, n)` (exact, tout générateur, plus
  rapide), `shuffle(range, g)`, et les moteurs initialisés à volonté. Rupture : l'ordre
  des arguments et les suites tirées changent.

### cpp-init (outil)
- Preset release des nouveaux projets en `-march=native` sur x86-64 (fma, AVX2…).

### special (nouvelle)
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
