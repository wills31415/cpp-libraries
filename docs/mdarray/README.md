# `wsl::mdarray`

`#include <wsl/mdarray>` — tableaux multidimensionnels et calcul **par les axes**, écrit
comme au tableau : le membre de gauche décide des axes libres, un axe répété à droite
seulement est sommé.

```cpp
auto [i, j, k] = wsl::axis<3>();
c[i, k] = a[i, j] * b[j, k];                  // produit matriciel
y[i] = tanh(W[i, j] * x[j] + b[i]);
r[i] = sum(j)(a[i, j]);
A[i, i] = 1;
y[i] = wsl::portable(W[i, j] * x[j]);         // mêmes bits partout, fma comprise
```

- **Conteneur** `mdarray` (P1684 : sémantique de valeur, sans allocation pour des étendues
  fixes, littéral `{{1, 2}, {3, 4}}`) ; `extents`, `dextents`, `dims`, `layout_right`,
  `layout_left`, `layout_stride`, `mdspan` (noms et sémantique de `<mdspan>`).
- **Calcul sur tableaux entiers** : `y = x * x`, `hypot(x, y)`, fonctions de `<cmath>`,
  `vectorize`, masques et `where`, réductions `sum`, `prod`, `min`, `max`, `argmin`,
  `argmax`, `all`, `any`, `count`.
- **Positions** `mdposition` et déplacements `mddisplacement` ; `div_euclid`, `rem_euclid`.
- **Garde-fous à la compilation** : axes incompatibles, axe libre absent à gauche…
- **Modes** : déterministe par défaut ; `wsl::portable(…)` passe toute l'expression en
  mode portable, fonctions comprises (celles de `ef` et `sf`).

Feuille de route (vues, algèbre linéaire, tenseurs) : [design/mdarray.md](../../design/mdarray.md)
et [design/algebra.md](../../design/algebra.md). Exemple :
[`examples/mdarray/axes.cpp`](../../examples/mdarray/axes.cpp).
