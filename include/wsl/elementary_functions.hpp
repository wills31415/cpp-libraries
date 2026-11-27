#pragma once

// wsl::elementary_functions (alias conseillé : namespace ef = wsl::elementary_functions;) —
// fonctions élémentaires bit-exactes (design/special.md) : même résultat au bit près sur toute
// plateforme, en float et en double, évaluables à la compilation.
//   algébriques     sqrt, cbrt, hypot
//   exponentielles  exp, exp2, expm1, log, log2, log10, log1p, pow
//   trigonométriques sin, cos, tan, atan, atan2
// Noms qualifiés ou alias plutôt qu'un using namespace, qui rendrait exp(2.0) ambigu avec le
// ::exp de <math.h>.

#include <wsl/elementary_functions/algebraic.hpp>
#include <wsl/elementary_functions/exponential.hpp>
#include <wsl/elementary_functions/trigonometric.hpp>
