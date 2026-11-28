# `wsl::measure_theory` (`meas`)

`#include <wsl/measure_theory>` — pour l'instant les **mesures de référence** `lebesgue`,
`counting` et `surface` : la densité d'une loi est sa dérivée de Radon–Nikodym par rapport à
l'une d'elles (`proba::reference_measure_t<L>`). L'intégrale, les mesures à densité, les produits et
les images (où rejoindra l'algèbre des lois de `probability`) sont à concevoir (décisions 39
et 42 de [design/interfaces.md](../../design/interfaces.md)).
