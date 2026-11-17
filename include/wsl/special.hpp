#pragma once

// wsl::special — fonctions mathématiques bit-exactes (voir docs/special.md).
//
// Même résultat au bit près sur toute plateforme, en float et en double, et
// évaluables à la compilation. Namespace wsl::special : employer des noms
// qualifiés ou un alias (namespace sp = wsl::special;) plutôt qu'un
// using namespace, qui rendrait exp(2.0) ambigu avec le ::exp de <math.h>.

#include <wsl/special/algebraic.hpp>
#include <wsl/special/elementary.hpp>
#include <wsl/special/error_function.hpp>
#include <wsl/special/gamma.hpp>
#include <wsl/special/incomplete_gamma.hpp>
#include <wsl/special/trigonometric.hpp>
