# `wsl::elementary_functions` (`ef`)

`#include <wsl/elementary_functions>` — les fonctions élémentaires de `<cmath>`, en `float`
et `double`, **reproductibles** et `constexpr`, sans aucun appel à la libm du système.

| En-tête | Fonctions |
|---|---|
| `<wsl/elementary_functions/exponential>` | `exp`, `exp2`, `expm1`, `log`, `log2`, `log10`, `log1p`, `pow` |
| `<wsl/elementary_functions/trigonometric>` | `sin`, `cos`, `tan`, `atan`, `atan2` |
| `<wsl/elementary_functions/algebraic>` | `sqrt`, `cbrt`, `hypot` |

```cpp
namespace ef = wsl::elementary_functions;
double a = ef::exp(x);                    // mode déterministe (défaut) : le plus rapide
double b = ef::exp(x, wsl::portable);     // mêmes bits sur toute plateforme
constexpr double c = ef::log(2.0, wsl::portable);
```

- **Deux modes** de reproductibilité ([modes](../README.md#modes-de-calcul)) : la voie
  déterministe laisse le compilateur fusionner les produits en fma (de 10 à 50 % plus
  rapide selon la fonction) ; la voie portable les en empêche. Le calcul à la compilation
  suit la voie portable.
- **Précision** : erreur mesurée de 0,50 ulp environ (arrondi fidèle, moins de 1 ulp
  garanti) ; table : [design/special.md](../../design/special.md), §6. `wsl::exact`
  (arrondi correct) est **refusé** pour l'instant : c'est le prochain chantier (décision 41).
- Valeurs particulières, signe des zéros, infinis et NaN : ceux de l'annexe F du standard C.

Exemple : [`examples/elementary_functions/modes.cpp`](../../examples/elementary_functions/modes.cpp).
