#pragma once

// wsl::random — les lois de probabilité, du point de vue mathématique (voir
// docs/random.md). Phase 2 : le socle.
//
//   moteurs      splitmix64, xoshiro256pp, pcg64, philox4x64
//   uniformes    random_bits, uniform_below, uniform_int, uniform_real(_fast),
//                bernoulli, shuffle
//
// Tirages reproductibles au bit près sur toute plateforme : algorithmes fixés,
// indépendants de la bibliothèque standard. Tout générateur standard (mt19937…)
// est accepté.

#include <wsl/random/engines.hpp>
#include <wsl/random/uniform.hpp>
