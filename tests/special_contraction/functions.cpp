// Compilé par cmake/check_no_contraction.cmake dans la configuration la plus
// hostile (fma disponible, -ffp-contract=fast) : le code engendré pour ces
// fonctions ne doit contenir aucune fma. Une fma ici veut dire qu'un produit
// entrant dans une addition n'est pas passé par detail::fmul.

#include <wsl/random.hpp>
#include <wsl/special.hpp>

extern "C"
{
    double check_exp(double x) { return wsl::special::exp(x); }
    double check_log(double x) { return wsl::special::log(x); }
    float check_expf(float x) { return wsl::special::exp(x); }
    float check_logf(float x) { return wsl::special::log(x); }
    double check_exp2(double x) { return wsl::special::exp2(x); }
    double check_expm1(double x) { return wsl::special::expm1(x); }
    double check_log2(double x) { return wsl::special::log2(x); }
    double check_log10(double x) { return wsl::special::log10(x); }
    double check_log1p(double x) { return wsl::special::log1p(x); }
    double check_pow(double x, double y) { return wsl::special::pow(x, y); }
    double check_sin(double x) { return wsl::special::sin(x); }
    double check_cos(double x) { return wsl::special::cos(x); }
    double check_tan(double x) { return wsl::special::tan(x); }
    double check_atan(double x) { return wsl::special::atan(x); }
    double check_atan2(double y, double x) { return wsl::special::atan2(y, x); }
    double check_cbrt(double x) { return wsl::special::cbrt(x); }
    double check_hypot(double x, double y) { return wsl::special::hypot(x, y); }
    double check_erf(double x) { return wsl::special::erf(x); }
    double check_erfc(double x) { return wsl::special::erfc(x); }
    double check_lgamma(double x) { return wsl::special::lgamma(x); }
    double check_tgamma(double x) { return wsl::special::tgamma(x); }
    double check_digamma(double x) { return wsl::special::digamma(x); }
    double check_gamma_p(double a, double x) { return wsl::special::gamma_p(a, x); }
    double check_gamma_q(double a, double x) { return wsl::special::gamma_q(a, x); }
    double check_bessel_i(double nu, double x) { return wsl::special::cyl_bessel_i(nu, x); }
    float check_bessel_if(float nu, float x) { return wsl::special::cyl_bessel_i(nu, x); }
    double check_sum_bessel_i(double nu, double x) { return 1.0 + wsl::special::cyl_bessel_i(nu, x); }
    double check_sum_gamma_p(double a, double x) { return 1.0 + wsl::special::gamma_p(a, x); }
    double check_gamma_p_inv(double a, double p) { return wsl::special::gamma_p_inv(a, p); }
    double check_ibeta(double a, double b, double x) { return wsl::special::ibeta(a, b, x); }
    double check_ibetac(double a, double b, double x) { return wsl::special::ibetac(a, b, x); }
    double check_sum_ibeta(double a, double b, double x) { return 1.0 + wsl::special::ibeta(a, b, x); }
    double check_beta(double a, double b) { return wsl::special::beta(a, b); }
    double check_ibeta_inv(double a, double b, double p) { return wsl::special::ibeta_inv(a, b, p); }
    double check_ibetac_inv(double a, double b, double q) { return wsl::special::ibetac_inv(a, b, q); }
    double check_sum_ibeta_inv(double a, double b, double p) { return 1.0 + wsl::special::ibeta_inv(a, b, p); }
    double check_sum_lbeta(double a, double b) { return 1.0 + wsl::special::lbeta(a, b); }
    double check_gamma_q_inv(double a, double q) { return wsl::special::gamma_q_inv(a, q); }
    double check_sum_gamma_q_inv(double a, double q) { return 1.0 + wsl::special::gamma_q_inv(a, q); }
    double check_sum_digamma(double x) { return 1.0 + wsl::special::digamma(x); }
    double check_erf_inv(double y) { return wsl::special::erf_inv(y); }
    double check_erfc_inv(double z) { return wsl::special::erfc_inv(z); }
    double check_uniform_cdf(double a, double b, double x) { return cdf(wsl::uniform({.low = a, .high = b}), x); }
    double check_uniform_quantile(double a, double b, double p) { return quantile(wsl::uniform({.low = a, .high = b}), p); }
    double check_uniform_variance(double a, double b) { return wsl::uniform({.low = a, .high = b}).variance(); }
    double check_discrete_variance(long long a, long long b) { return wsl::discrete_uniform({.low = a, .high = b}).variance(); }
    double check_exponential_pdf(double l, double x) { return pdf(wsl::exponential({.rate = l}), x); }
    double check_exponential_cdf(double l, double x) { return cdf(wsl::exponential({.rate = l}), x); }
    double check_exponential_ccdf(double l, double x) { return ccdf(wsl::exponential({.rate = l}), x); }
    double check_exponential_log_pdf(double l, double x) { return log_pdf(wsl::exponential({.rate = l}), x); }
    double check_exponential_quantile(double l, double p) { return quantile(wsl::exponential({.rate = l}), p); }
    double check_exponential_variance(double l) { return wsl::exponential({.rate = l}).variance(); }
    double check_normal_pdf(double m, double s, double x) { return pdf(wsl::normal({.mean = m, .sd = s}), x); }
    double check_normal_cdf(double m, double s, double x) { return cdf(wsl::normal({.mean = m, .sd = s}), x); }
    double check_normal_log_pdf(double m, double s, double x) { return log_pdf(wsl::normal({.mean = m, .sd = s}), x); }
    double check_normal_quantile(double m, double s, double p) { return quantile(wsl::normal({.mean = m, .sd = s}), p); }
    double check_normal_entropy(double m, double s) { return wsl::normal({.mean = m, .sd = s}).entropy(); }
    double check_exponential_sample(double l, wsl::xoshiro256pp& g) { return sample(wsl::exponential({.rate = l}), g); }
    double check_exponential_fast(double l, wsl::xoshiro256pp& g) { return sample(wsl::exponential({.rate = l}), g, wsl::fast); }
    double check_normal_sample(double m, double s, wsl::xoshiro256pp& g) { return sample(wsl::normal({.mean = m, .sd = s}), g); }
    float check_normal_sample_float(float m, float s, wsl::xoshiro256pp& g) { return sample(wsl::normal<float>({.mean = m, .sd = s}), g); }
    double check_normal_fast(double m, double s, wsl::xoshiro256pp& g) { return sample(wsl::normal({.mean = m, .sd = s}), g, wsl::fast); }
    double check_geometric_pmf(double p, long long k) { return pmf(wsl::geometric({.p = p}), k); }
    double check_geometric_cdf(double p, long long k) { return cdf(wsl::geometric({.p = p}), k); }
    long long check_geometric_quantile(double p, double u) { return quantile(wsl::geometric({.p = p}), u); }
    double check_geometric_entropy(double p) { return wsl::geometric({.p = p}).entropy(); }
    long long check_geometric_sample(double p, wsl::xoshiro256pp& g) { return sample(wsl::geometric({.p = p}), g); }
    double check_poisson_pmf(double m, long long k) { return pmf(wsl::poisson({.mean = m}), k); }
    double check_poisson_log_pmf(double m, long long k) { return log_pmf(wsl::poisson({.mean = m}), k); }
    long long check_poisson_sample(wsl::poisson<> const& P, wsl::xoshiro256pp& g) { return sample(P, g); }
    void check_poisson_make(double m, wsl::poisson<>* out) { *out = wsl::poisson({.mean = m}); }
    double check_binomial_pmf(long long n, double p, long long k) { return pmf(wsl::binomial({.trials = n, .p = p}), k); }
    double check_binomial_log_pmf(long long n, double p, long long k) { return log_pmf(wsl::binomial({.trials = n, .p = p}), k); }
    double check_binomial_variance(long long n, double p) { return wsl::binomial({.trials = n, .p = p}).variance(); }
    long long check_binomial_sample(wsl::binomial<> const& B, wsl::xoshiro256pp& g) { return sample(B, g); }
    void check_binomial_make(long long n, double p, wsl::binomial<>* out) { *out = wsl::binomial({.trials = n, .p = p}); }
    double check_gamma_pdf(double a, double r, double x) { return pdf(wsl::gamma({.shape = a, .rate = r}), x); }
    double check_gamma_cdf(double a, double r, double x) { return cdf(wsl::gamma({.shape = a, .rate = r}), x); }
    double check_student_cdf(double n, double t) { return cdf(wsl::student_t({.degrees_of_freedom = n}), t); }
    double check_fisher_cdf(double a, double b, double x) { return cdf(wsl::fisher_f({.d1 = a, .d2 = b}), x); }
    double check_fisher_pdf(double a, double b, double x) { return pdf(wsl::fisher_f({.d1 = a, .d2 = b}), x); }
    double check_fisher_quantile(double a, double b, double p) { return quantile(wsl::fisher_f({.d1 = a, .d2 = b}), p); }
    double check_fisher_entropy(double a, double b) { return wsl::fisher_f({.d1 = a, .d2 = b}).entropy(); }
    double check_fisher_moments(double a, double b) { wsl::fisher_f const F({.d1 = a, .d2 = b}); return F.variance() + F.skewness() + F.excess_kurtosis() + F.mode(); }
    double check_weibull_cdf(double k, double l, double x) { return cdf(wsl::weibull({.shape = k, .scale = l}), x); }
    double check_weibull_pdf(double k, double l, double x) { return pdf(wsl::weibull({.shape = k, .scale = l}), x); }
    double check_weibull_quantile(double k, double l, double p) { return quantile(wsl::weibull({.shape = k, .scale = l}), p); }
    double check_weibull_moments(double k, double l) { wsl::weibull const W({.shape = k, .scale = l}); return W.mean() + W.variance() + W.median() + W.mode() + W.skewness() + W.entropy(); }
    double check_laplace_cdf(double m, double b, double x) { return cdf(wsl::laplace({.location = m, .scale = b}), x); }
    double check_laplace_quantile(double m, double b, double p) { return quantile(wsl::laplace({.location = m, .scale = b}), p); }
    double check_logistic_pdf(double m, double s, double x) { return pdf(wsl::logistic({.location = m, .scale = s}), x); }
    double check_logistic_quantile(double m, double s, double p) { return quantile(wsl::logistic({.location = m, .scale = s}), p); }
    double check_laplace_variance(double m, double b) { return 1 + wsl::laplace({.location = m, .scale = b}).variance(); }
    double check_pareto_cdf(double a, double b, double x) { return cdf(wsl::pareto({.scale = a, .shape = b}), x); }
    // Quantiles de survie et algèbre des lois.
    double check_cquantile_gamma(double a, double q) { return cquantile(wsl::gamma({.shape = a, .rate = 1.0}), q); }
    double check_cquantile_fisher(double a, double b, double q) { return cquantile(wsl::fisher_f({.d1 = a, .d2 = b}), q); }
    double check_cquantile_sum(double a, double q) { return 1 + cquantile(wsl::weibull({.shape = a, .scale = 2.0}), q) + cquantile(wsl::logistic({.location = a, .scale = 2.0}), q)
                                                            + cquantile(wsl::lognormal({.mu = a, .sigma = 2.0}), q) + cquantile(wsl::uniform({.low = a, .high = a + 2}), q); }
    double check_affine_closure(double a, double b, double x) { return cdf(a * wsl::normal({.mean = 1.0, .sd = 2.0}) + b, x) + (a * wsl::laplace({.location = 1.0, .scale = 2.0}) + b).location(); }
    double check_affine_cdf(double a, double b, double x) { return cdf(wsl::affine(wsl::student_t({.degrees_of_freedom = 3.0}), {.scale = a, .shift = b}), x); }
    double check_affine_quantile(double a, double b, double p) { return 1 + quantile(wsl::affine(wsl::student_t({.degrees_of_freedom = 3.0}), {.scale = a, .shift = b}), p); }
    double check_affine_moments(double a, double b) { auto const Z = wsl::affine(wsl::gamma({.shape = 3.0, .rate = 1.0}), {.scale = a, .shift = b}); return Z.mean() + Z.variance() + Z.entropy() + Z.support().upper; }
    double check_truncated_cdf(double a, double b, double x) { auto const T = wsl::truncate(wsl::normal({.mean = 0.0, .sd = 1.0}), a, b); return cdf(T, x) + ccdf(T, x) + log_pdf(T, x); }
    double check_truncated_quantile(double a, double b, double p) { auto const T = wsl::truncate(wsl::normal({.mean = 0.0, .sd = 1.0}), a, b); return 1 + quantile(T, p) + cquantile(T, p); }
    double check_mixture_functions(double w, double x) { wsl::mixture const M(wsl::weighted(w, wsl::normal({.mean = -2.0, .sd = 1.0})), wsl::weighted(1.0, wsl::gamma({.shape = 3.0, .rate = 1.0})));
                                                         return pdf(M, x) + cdf(M, x) + log_pdf(M, x) + quantile(M, 0.25) + M.mean() + M.variance() + M.skewness() + M.excess_kurtosis(); }
    double check_sum_normal(double a, double b) { auto const S = wsl::normal({.mean = a, .sd = b}) + wsl::normal({.mean = 1.0, .sd = 2.0}); return 1 + S.mean() + S.sd(); }
    double check_sum_moments(double a, double b) { auto const S = wsl::gamma({.shape = a, .rate = b}) + wsl::student_t({.degrees_of_freedom = 9.0}); return S.mean() + S.variance() + S.skewness() + S.excess_kurtosis(); }
    double check_sum_convolution(double a, double b, long long k) { auto const S = wsl::poisson({.mean = a}) + wsl::binomial({.trials = 20, .p = b}); return 1 + pmf(S, k) + cdf(S, k) + ccdf(S, k) + double(quantile(S, b)); }
    double check_sum_gamma_family(double a, double b, double x) { auto const S = wsl::gamma({.shape = a, .rate = b}) + wsl::exponential({.rate = b}); return 1 + cdf(S, x) + S.mean(); }
    double check_independent(double a, double b) { wsl::independent const J(wsl::normal({.mean = a, .sd = 1.0}), wsl::gamma({.shape = b, .rate = 1.0})); return 1 + pdf(J, {a, b}) + log_pdf(J, {a, b}) + cdf(J, {a, b}) + J.entropy(); }
    double check_pushforward(double a, double y) { auto const L = wsl::pushforward(wsl::normal({.mean = a, .sd = 1.0}), wsl::increasing([](double t) { return wsl::special::exp(t); }, [](double t) { return wsl::special::log(t); }, [](double t) { return 1 / t; }));
                                                   return 1 + cdf(L, y) + pdf(L, y) + quantile(L, 0.3); }
    // Lois multivariées.
    double check_dirichlet(double a, double x) { wsl::dirichlet<double, 3> const D({.alpha = {a, 2.0, 3.0}}); auto const m = D.mean(); auto const c = D.covariance();
                                                 return 1 + pdf(D, {x, 0.25, 0.75 - x}) + log_pdf(D, {x, 0.25, 0.75 - x}) + D.entropy() + m[0] + c(0, 1) + D.variance()[1]; }
    double check_multinomial(double p, long long k) { wsl::multinomial<std::int64_t, 3> const M({.trials = 20, .p = {p, 0.5 - p, 0.5}});
                                                       return 1 + pmf(M, {k, 5, 15 - k}) + M.mean()[0] + M.covariance()(0, 1) + M.variance()[2]; }
    double check_mvn(double a, double x) { wsl::multivariate_normal<double, 2> const X({.mean = {a, 0}, .covariance = {{2, 0.5}, {0.5, 1}}});
                                           return 1 + pdf(X, {x, a}) + log_pdf(X, {x, a}) + X.entropy() + X.mahalanobis_squared({x, x}); }
    double check_sphere(double x) { wsl::uniform_sphere<double, 3> const S; wsl::uniform_ball<double, 5> const B;
                                    return 1 + pdf(S, {x, 0.0, 0.0}) + S.surface_area() + S.entropy() + B.volume() + pdf(B, {x, 0, 0, 0, 0}); }
    double check_mvt(double nu, double x) { wsl::multivariate_student_t<double, 2> const X({.degrees_of_freedom = nu, .location = {x, 0}, .scale = {{2, 0.5}, {0.5, 1}}});
                                            return 1 + pdf(X, {x, 1}) + log_pdf(X, {x, x}) + X.entropy() + X.covariance()(0, 1); }
    double check_wishart(double nu, double x) { wsl::wishart<double, 2> const W({.degrees_of_freedom = nu, .scale = {{2, 0.5}, {0.5, 1}}}); wsl::square_matrix<double, 2> const A{{x, 0.1}, {0.1, 1}};
                                                wsl::inverse_wishart<double, 2> const I({.degrees_of_freedom = nu, .scale = {{2, 0.5}, {0.5, 1}}});
                                                return 1 + pdf(W, A) + log_pdf(W, A) + W.entropy() + W.mean()(0, 1) + W.variance()(0, 0) + W.mode()(1, 1) + log_pdf(I, A) + I.entropy() + I.mean()(0, 0); }
    double check_random_matrix(double x) { wsl::goe<double, 2> const H; wsl::gue<double, 2> const K;
                                           return 1 + pdf(H, wsl::square_matrix<double, 2>{{x, 0.5}, {0.5, x}}) + log_pdf(K, wsl::square_matrix<std::complex<double>, 2>{{x, {0.5, x}}, {{0.5, -x}, 1}}); }
    double check_von_mises(double m, double k, double x) { wsl::von_mises const V({.mu = m, .kappa = k}); return 1 + pdf(V, x) + log_pdf(V, x) + V.entropy() + V.circular_variance() + V.mean_resultant_length(); }
    double check_vmf(double k, double t) { wsl::von_mises_fisher<double, 3> const F({.mean_direction = {0, 0.6, 0.8}, .kappa = k});
                                           return 1 + pdf(F, {t, 0.0, 0.0}) + log_pdf(F, {0.0, t, t}) + F.entropy() + F.mean()[1]; }
    double check_permutation_laws(double theta, double q) { wsl::permutation const p{2, 0, 1, 4, 3}; wsl::ewens const E({.size = 5, .theta = theta}); wsl::mallows const M({.size = 5, .q = q});
                                                            return 1 + pmf(E, p) + log_pmf(M, p) + E.mean_cycles() + M.mean_inversions() + pmf(wsl::uniform_permutation({.size = 5}), p); }
    double check_subset_laws(double p) { wsl::bernoulli_subset const B({.size = 10, .p = p}); wsl::uniform_combination const C({.size = 10, .count = 3});
                                         return 1 + pmf(B, {1, 4}) + log_pmf(C, {1, 4, 7}) + B.mean_size() + C.entropy() + wsl::uniform_word({.length = 3, .alphabet = 5}).entropy(); }
    double check_graph_laws(double p) { wsl::erdos_renyi const G({.vertices = 6, .p = p}); wsl::uniform_graph const H({.vertices = 6, .edges = 3});
                                        return 1 + pmf(G, {{0, 1}, {2, 3}}) + log_pmf(G, {{0, 1}}) + G.mean_edges() + H.entropy() + wsl::uniform_tree({.vertices = 9}).entropy()
                                             + wsl::uniform_binary_tree({.nodes = 7}).entropy() + wsl::uniform_set_partition({.size = 9}).entropy() + wsl::uniform_composition({.n = 9, .parts = 3}).entropy(); }
    double check_categorical(double a, double b) { wsl::categorical const C({.p = {a, b, 0.5}}); return 1 + pmf(C, 1) + log_pmf(C, 2) + cdf(C, 1) + ccdf(C, 0) + C.mean() + C.variance() + C.entropy(); }
    double check_markov(double a) { wsl::markov_chain const M({.transition = {{a, 1 - a}, {0.5, 0.5}}, .start = 0}); return 1 + pmf(M.at(7), 1) + pmf(M.stationary(), 0) + pmf(M.steps(2), {0, 1, 0}); }
    double check_urn_laws(double a, double t) { wsl::dirichlet_multinomial const D({.trials = 5, .alpha = {a, 1.5}}); wsl::ewens_partition const E({.size = 6, .theta = t, .discount = 0.3});
                                                return 1 + pmf(D, {2, 3}) + D.mean()[0] + D.variance()[1] + pmf(E, {0, 1, 0, 2, 1, 0}) + E.mean_blocks() + pmf(wsl::polya_urn({.balls = {a, 2}, .added = t}).steps(2), {0, 1, 1}); }
    double check_pareto_moments(double a, double b) { wsl::pareto const P({.scale = a, .shape = b}); return P.mean() + P.variance() + P.skewness() + P.excess_kurtosis() + P.entropy(); }
    double check_student_quantile(double n, double p) { return quantile(wsl::student_t({.degrees_of_freedom = n}), p); }
    double check_student_pdf(double n, double t) { return pdf(wsl::student_t({.degrees_of_freedom = n}), t); }
    double check_student_entropy(double n) { return wsl::student_t({.degrees_of_freedom = n}).entropy(); }
    double check_student_sample(double n, wsl::xoshiro256pp& g) { return sample(wsl::student_t({.degrees_of_freedom = n}), g); }
    double check_lognormal_cdf(double m, double s, double x) { return cdf(wsl::lognormal({.mu = m, .sigma = s}), x); }
    double check_lognormal_quantile(double m, double s, double p) { return quantile(wsl::lognormal({.mu = m, .sigma = s}), p); }
    double check_lognormal_variance(double m, double s) { return wsl::lognormal({.mu = m, .sigma = s}).variance(); }
    double check_cauchy_cdf(double a, double b, double x) { return cdf(wsl::cauchy({.location = a, .scale = b}), x); }
    double check_cauchy_quantile(double a, double b, double p) { return quantile(wsl::cauchy({.location = a, .scale = b}), p); }
    double check_cauchy_sample(double a, double b, wsl::xoshiro256pp& g) { return sample(wsl::cauchy({.location = a, .scale = b}), g); }
    double check_gamma_quantile(double a, double t, double p) { return quantile(wsl::gamma({.shape = a, .scale = t}), p); }
    double check_gamma_entropy(double a, double r) { return wsl::gamma({.shape = a, .rate = r}).entropy(); }
    double check_beta_entropy(double a, double b) { return wsl::beta({.alpha = a, .beta = b}).entropy(); }
    double check_beta_cdf(double a, double b, double x) { return cdf(wsl::beta({.alpha = a, .beta = b}), x); }
    double check_poisson_cdf(double m, long long k) { return cdf(wsl::poisson({.mean = m}), k); }
    long long check_poisson_quantile(double m, double u) { return quantile(wsl::poisson({.mean = m}), u); }
    long long check_binomial_quantile(long long n, double p, double u) { return quantile(wsl::binomial({.trials = n, .p = p}), u); }
    double check_gamma_log_pdf(double a, double t, double x) { return log_pdf(wsl::gamma({.shape = a, .scale = t}), x); }
    double check_gamma_sample(wsl::gamma<> const& G, wsl::xoshiro256pp& g) { return sample(G, g); }
    void check_gamma_make(double a, double r, wsl::gamma<>* out) { *out = wsl::gamma({.shape = a, .rate = r}); }
    double check_gamma_variance(double a, double t) { return wsl::gamma({.shape = a, .scale = t}).variance(); }
    double check_beta_pdf(double a, double b, double x) { return pdf(wsl::beta({.alpha = a, .beta = b}), x); }
    double check_beta_log_pdf(double a, double b, double x) { return log_pdf(wsl::beta({.alpha = a, .beta = b}), x); }
    double check_beta_sample(wsl::beta<> const& B, wsl::xoshiro256pp& g) { return sample(B, g); }
    void check_beta_make(double a, double b, wsl::beta<>* out) { *out = wsl::beta({.alpha = a, .beta = b}); }
    double check_beta_moments(double a, double b) { wsl::beta const B({.alpha = a, .beta = b}); return B.variance() + B.excess_kurtosis() + B.skewness() + B.mode(); }
    // Produits rendus puis additionnés par l'appelant (GCC les fusionne après intégration).
    double check_sum_exp(double x) { return 1.0 + wsl::special::exp(x); }
    double check_sum_log(double x) { return 1.0 + wsl::special::log(x); }
    double check_sum_exp2(double x) { return 1.0 + wsl::special::exp2(x); }
    double check_sum_expm1(double x) { return 1.0 + wsl::special::expm1(x); }
    double check_sum_log2(double x) { return 1.0 + wsl::special::log2(x); }
    double check_sum_log10(double x) { return 1.0 + wsl::special::log10(x); }
    double check_sum_log1p(double x) { return 1.0 + wsl::special::log1p(x); }
    double check_sum_sin(double x) { return 1.0 + wsl::special::sin(x); }
    double check_sum_cos(double x) { return 1.0 + wsl::special::cos(x); }
    double check_sum_tan(double x) { return 1.0 + wsl::special::tan(x); }
    double check_sum_atan(double x) { return 1.0 + wsl::special::atan(x); }
    double check_sum_cbrt(double x) { return 1.0 + wsl::special::cbrt(x); }
    double check_sum_erf(double x) { return 1.0 + wsl::special::erf(x); }
    double check_sum_erfc(double x) { return 1.0 + wsl::special::erfc(x); }
    double check_sum_lgamma(double x) { return 1.0 + wsl::special::lgamma(x); }
    double check_sum_tgamma(double x) { return 1.0 + wsl::special::tgamma(x); }
    double check_sum_pow(double x, double y) { return 1.0 + wsl::special::pow(x, y); }
    double check_sum_atan2(double y, double x) { return 1.0 + wsl::special::atan2(y, x); }
    double check_sum_hypot(double x, double y) { return 1.0 + wsl::special::hypot(x, y); }
}
