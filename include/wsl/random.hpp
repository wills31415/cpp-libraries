#pragma once

// wsl::random (alias conseillé : namespace rnd = wsl::random;) — le hasard (design/random.md ;
// design/interfaces.md, décision 14) : moteurs (catalogue wsl::random::engines : splitmix64,
// xoshiro256pp, pcg64, philox4x64) et fonctions de tirage (random_bits, uniform_below,
// uniform_int, uniform_real, bernoulli_trial, shuffle). Les lois sont dans
// wsl::probability (<wsl/probability>).

#include <wsl/random/engines.hpp>
#include <wsl/random/uniform.hpp>
