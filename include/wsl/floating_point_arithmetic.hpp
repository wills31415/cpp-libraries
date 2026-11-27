#pragma once

// wsl::floating_point_arithmetic (alias conseillé : namespace fpa = wsl::floating_point_arithmetic;)
// — primitives flottantes : transformations sans erreur (two_sum, two_product), découpage de
// Veltkamp, barrières contre la fusion en fma (design/interfaces.md, décisions 37 et 38).
// Interface publique à concevoir : tout est encore dans detail (« découper maintenant,
// exposer plus tard »).

#include <wsl/floating_point_arithmetic/detail/primitives.hpp>
