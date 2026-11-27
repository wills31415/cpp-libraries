#pragma once

// Fonctions de <cmath> élément par élément (design/mdarray.md, I11), une seule liste
// pour les deux mondes (tableaux entiers, axes) : X(nom).

#define WSL_MDARRAY_UNARY_FUNCTIONS(X)                                                                                 \
    X(abs) X(sqrt) X(cbrt) X(exp) X(exp2) X(expm1) X(log) X(log2) X(log10) X(log1p) X(sin) X(cos) X(tan) X(asin)      \
    X(acos) X(atan) X(sinh) X(cosh) X(tanh) X(asinh) X(acosh) X(atanh) X(erf) X(erfc) X(tgamma) X(lgamma) X(floor)    \
    X(ceil) X(trunc) X(round) X(isnan) X(isinf) X(isfinite)

#define WSL_MDARRAY_BINARY_FUNCTIONS(X) X(pow) X(atan2) X(hypot) X(fmod) X(fmin) X(fmax) X(copysign)
