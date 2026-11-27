#pragma once

// wsl::special_functions (alias conseillé : namespace sf = wsl::special_functions;) — fonctions
// spéciales bit-exactes (design/special.md), une notion par chapitre du DLMF :
//   gamma             lgamma, tgamma, digamma, beta, lbeta (DLMF 5)
//   error_functions   erf, erfc, erf_inv, erfc_inv (DLMF 7)
//   incomplete_gamma  gamma_p, gamma_q et leurs inverses (DLMF 8)
//   incomplete_beta   ibeta, ibetac et leurs inverses (DLMF 8.17)
//   bessel            cyl_bessel_i (DLMF 10)

#include <wsl/special_functions/bessel.hpp>
#include <wsl/special_functions/error_functions.hpp>
#include <wsl/special_functions/gamma.hpp>
#include <wsl/special_functions/incomplete_beta.hpp>
#include <wsl/special_functions/incomplete_gamma.hpp>
