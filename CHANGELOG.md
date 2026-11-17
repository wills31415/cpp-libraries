# Journal des versions

Numérotation `MAJEUR.MINEUR.CORRECTIF`. Tant que `MAJEUR` vaut 0, une version
mineure peut casser l'interface ; un correctif, jamais. Chaque entrée est rangée
par bibliothèque, pour savoir quoi vérifier en montant de version dans un projet.

## [Non publié]

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
