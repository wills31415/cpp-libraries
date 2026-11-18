#pragma once

// wsl::random — les lois de probabilité, du point de vue mathématique (voir
// docs/random.md).
//
//   moteurs      splitmix64, xoshiro256pp, pcg64, philox4x64
//   uniformes    random_bits, uniform_below, uniform_int, uniform_real(_fast),
//                bernoulli_trial, shuffle
//   lois         concept et opérations (law.hpp) ; uniform, discrete_uniform,
//                bernoulli, exponential, normal, geometric, poisson, binomial,
//                gamma, beta, chi_squared, student_t, lognormal, cauchy, fisher_f,
//                weibull, laplace, logistic, pareto
//
// Tirages reproductibles au bit près sur toute plateforme : algorithmes fixés,
// indépendants de la bibliothèque standard. Tout générateur standard (mt19937…)
// est accepté.

#include <wsl/random/affine.hpp>
#include <wsl/random/basic_laws.hpp>
#include <wsl/random/beta.hpp>
#include <wsl/random/binomial.hpp>
#include <wsl/random/cauchy.hpp>
#include <wsl/random/chi_squared.hpp>
#include <wsl/random/engines.hpp>
#include <wsl/random/exponential.hpp>
#include <wsl/random/fisher_f.hpp>
#include <wsl/random/gamma.hpp>
#include <wsl/random/independent.hpp>
#include <wsl/random/geometric.hpp>
#include <wsl/random/laplace.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/logistic.hpp>
#include <wsl/random/mixture.hpp>
#include <wsl/random/lognormal.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/random/pareto.hpp>
#include <wsl/random/poisson.hpp>
#include <wsl/random/pushforward.hpp>
#include <wsl/random/student_t.hpp>
#include <wsl/random/sum.hpp>
#include <wsl/random/truncated.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/random/weibull.hpp>
#include <wsl/random/dirichlet.hpp>
#include <wsl/random/multinomial.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/multivariate_normal.hpp>
#include <wsl/random/sphere.hpp>
#include <wsl/random/multivariate_student_t.hpp>
#include <wsl/random/wishart.hpp>
#include <wsl/random/random_matrix.hpp>
#include <wsl/random/von_mises.hpp>
#include <wsl/random/subsets.hpp>
#include <wsl/random/graphs.hpp>
#include <wsl/random/discrete_processes.hpp>
#include <wsl/random/continuous_processes.hpp>
#include <wsl/random/process.hpp>
#include <wsl/random/categorical.hpp>
#include <wsl/random/permutations.hpp>
