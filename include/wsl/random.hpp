#pragma once

// wsl::random — les lois de probabilité, du point de vue mathématique (voir
// docs/random.md).
//
//   moteurs      splitmix64, xoshiro256pp, pcg64, philox4x64
//   uniformes    random_bits, uniform_below, uniform_int, uniform_real(_fast),
//                bernoulli_trial, shuffle
//   lois         concept et opérations (law.hpp) ; uniform, discrete_uniform,
//                bernoulli, exponential, normal, geometric, poisson, binomial
//
// Tirages reproductibles au bit près sur toute plateforme : algorithmes fixés,
// indépendants de la bibliothèque standard. Tout générateur standard (mt19937…)
// est accepté.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/binomial.hpp>
#include <wsl/random/engines.hpp>
#include <wsl/random/exponential.hpp>
#include <wsl/random/geometric.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/random/poisson.hpp>
#include <wsl/random/uniform.hpp>
