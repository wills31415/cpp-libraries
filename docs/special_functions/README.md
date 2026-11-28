# `wsl::special_functions` (`sf`)

`#include <wsl/special_functions>` — fonctions spéciales en `float` et `double`,
reproductibles et `constexpr`, dans les deux modes comme
[`elementary_functions`](../elementary_functions/README.md).

| En-tête | Fonctions |
|---|---|
| `<wsl/special_functions/gamma>` | `tgamma`, `lgamma`, `digamma` |
| `<wsl/special_functions/error_functions>` | `erf`, `erfc`, `erf_inv`, `erfc_inv` |
| `<wsl/special_functions/incomplete_gamma>` | `gamma_p`, `gamma_q`, `gamma_p_inv`, `gamma_q_inv` |
| `<wsl/special_functions/incomplete_beta>` | `beta`, `lbeta`, `ibeta`, `ibetac`, `ibeta_inv`, `ibetac_inv` |
| `<wsl/special_functions/bessel>` | `cyl_bessel_i` |

```cpp
namespace sf = wsl::special_functions;
double g = sf::tgamma(4.5);
double q = sf::erfc(3.0, wsl::portable);
double x = sf::gamma_p_inv(2.5, 0.15);
```

- Noms et conventions de `<cmath>` (C++17, fonctions spéciales mathématiques) quand ils
  existent ; ceux de la littérature sinon (`gamma_p`, `ibeta`).
- **Précision** mesurée près de 0,50 ulp (absolue près des zéros de `lgamma` et `digamma`) :
  [design/special.md](../../design/special.md), §6 à 8.
- Les fonctions élémentaires sont dans `ef`, pas ici : `ef::exp`, non `sf::exp`.

Exemple : [`examples/special_functions/gamma.cpp`](../../examples/special_functions/gamma.cpp).
