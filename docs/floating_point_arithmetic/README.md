# `wsl::floating_point_arithmetic` (`fpa`)

Les primitives flottantes sur lesquelles reposent `elementary_functions` et
`special_functions` : transformations sans erreur (`two_sum`, `fast_two_sum`,
`two_product`), découpage de Veltkamp, arithmétique double-double, et les barrières contre
la fusion d'un produit en fma (`mul<R>` selon le mode, `unfused_mul` quand l'algorithme
s'appuie sur le produit arrondi : décision 38 de [design/interfaces.md](../../design/interfaces.md)).

**Pas encore d'interface publique** : tout est dans `fpa::detail`, à exposer après
conception. Conditions de la garantie bit-exacte (types IEEE 754, pas de `-ffast-math`,
évaluation en double, GCC ou Clang) : [design/special.md](../../design/special.md), §2.3 ;
elles sont vérifiées à la compilation.
