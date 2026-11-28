# `wsl::probability::stochastic_processes` (`stochproc`)

`#include <wsl/probability/stochastic_processes>` — processus stochastiques
(`stochproc::processes`) et leurs trajectoires, paresseuses et exactes : on n'engendre que les
instants observés, et la valeur en un instant ne dépend pas de la grille demandée.

- **Temps discret** : `markov_chain`, `random_walk`, `polya_urn`, `chinese_restaurant`.
- **Temps continu** : `brownian_motion`, `brownian_bridge`, `ornstein_uhlenbeck`,
  `poisson_process`, `compound_poisson`, `continuous_markov_chain`.

```cpp
namespace stochproc = wsl::probability::stochastic_processes;
auto W = stochproc::processes::brownian_motion({.drift = 0.0, .volatility = 1.0}).on_grid({0.25, 0.5, 1.0});
auto w = W(g);      // les valeurs aux trois instants
```

Exemple : [`examples/probability/stochastic_processes/brownian.cpp`](../../../examples/probability/stochastic_processes/brownian.cpp).
