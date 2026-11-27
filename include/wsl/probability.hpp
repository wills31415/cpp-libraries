#pragma once

// wsl::probability (alias conseillé : namespace proba = wsl::probability;) — les lois de
// probabilité comme objets mathématiques (design/random.md ; design/interfaces.md, décisions
// 14, 39, 40) : concept de loi, catalogue wsl::probability::laws, variables aléatoires,
// algèbre des lois (provisoirement ici : elle passera dans wsl::measure_theory), et le
// sous-domaine wsl::probability::stochastic_processes.

#include <wsl/probability/law.hpp>
#include <wsl/probability/affine.hpp>
#include <wsl/probability/convolution.hpp>
#include <wsl/probability/laws/arcsine.hpp>
#include <wsl/probability/laws/basic.hpp>
#include <wsl/probability/laws/beta.hpp>
#include <wsl/probability/laws/beta_binomial.hpp>
#include <wsl/probability/laws/binomial.hpp>
#include <wsl/probability/laws/categorical.hpp>
#include <wsl/probability/laws/cauchy.hpp>
#include <wsl/probability/laws/chi_squared.hpp>
#include <wsl/probability/laws/dirichlet.hpp>
#include <wsl/probability/laws/discrete_gaussian.hpp>
#include <wsl/probability/laws/erlang.hpp>
#include <wsl/probability/laws/exponential.hpp>
#include <wsl/probability/laws/fisher_f.hpp>
#include <wsl/probability/laws/gamma.hpp>
#include <wsl/probability/laws/geometric.hpp>
#include <wsl/probability/laws/graphs.hpp>
#include <wsl/probability/laws/gumbel.hpp>
#include <wsl/probability/laws/hypergeometric.hpp>
#include <wsl/probability/laws/inverse_gaussian.hpp>
#include <wsl/probability/laws/laplace.hpp>
#include <wsl/probability/laws/logarithmic.hpp>
#include <wsl/probability/laws/logistic.hpp>
#include <wsl/probability/laws/lognormal.hpp>
#include <wsl/probability/laws/multinomial.hpp>
#include <wsl/probability/laws/multivariate_normal.hpp>
#include <wsl/probability/laws/multivariate_student_t.hpp>
#include <wsl/probability/laws/negative_binomial.hpp>
#include <wsl/probability/laws/normal.hpp>
#include <wsl/probability/laws/pareto.hpp>
#include <wsl/probability/laws/permutations.hpp>
#include <wsl/probability/laws/poisson.hpp>
#include <wsl/probability/laws/random_matrix.hpp>
#include <wsl/probability/laws/rayleigh.hpp>
#include <wsl/probability/laws/rice.hpp>
#include <wsl/probability/laws/skellam.hpp>
#include <wsl/probability/laws/skew_normal.hpp>
#include <wsl/probability/laws/sphere.hpp>
#include <wsl/probability/laws/stable.hpp>
#include <wsl/probability/laws/student_t.hpp>
#include <wsl/probability/laws/subsets.hpp>
#include <wsl/probability/laws/triangular.hpp>
#include <wsl/probability/laws/von_mises.hpp>
#include <wsl/probability/laws/weibull.hpp>
#include <wsl/probability/laws/wishart.hpp>
#include <wsl/probability/laws/zipf.hpp>
#include <wsl/probability/mixture.hpp>
#include <wsl/probability/multivariate.hpp>
#include <wsl/probability/product_measure.hpp>
#include <wsl/probability/pushforward.hpp>
#include <wsl/probability/truncated.hpp>
#include <wsl/probability/variable.hpp>
#include <wsl/probability/stochastic_processes.hpp>
