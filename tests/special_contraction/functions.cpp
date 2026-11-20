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
    double check_uniform_cdf(double a, double b, double x) { return wsl::random::cdf(wsl::random::uniform({.low = a, .high = b}), x); }
    double check_uniform_quantile(double a, double b, double p) { return wsl::random::quantile(wsl::random::uniform({.low = a, .high = b}), p); }
    double check_uniform_variance(double a, double b) { return wsl::random::uniform({.low = a, .high = b}).variance(); }
    double check_discrete_variance(long long a, long long b) { return wsl::random::discrete_uniform({.low = a, .high = b}).variance(); }
    double check_exponential_pdf(double l, double x) { return wsl::random::pdf(wsl::random::exponential({.rate = l}), x); }
    double check_exponential_cdf(double l, double x) { return wsl::random::cdf(wsl::random::exponential({.rate = l}), x); }
    double check_exponential_ccdf(double l, double x) { return wsl::random::ccdf(wsl::random::exponential({.rate = l}), x); }
    double check_exponential_log_pdf(double l, double x) { return wsl::random::log_pdf(wsl::random::exponential({.rate = l}), x); }
    double check_exponential_quantile(double l, double p) { return wsl::random::quantile(wsl::random::exponential({.rate = l}), p); }
    double check_exponential_variance(double l) { return wsl::random::exponential({.rate = l}).variance(); }
    double check_normal_pdf(double m, double s, double x) { return wsl::random::pdf(wsl::random::normal({.mean = m, .sd = s}), x); }
    double check_normal_cdf(double m, double s, double x) { return wsl::random::cdf(wsl::random::normal({.mean = m, .sd = s}), x); }
    double check_normal_log_pdf(double m, double s, double x) { return wsl::random::log_pdf(wsl::random::normal({.mean = m, .sd = s}), x); }
    double check_normal_quantile(double m, double s, double p) { return wsl::random::quantile(wsl::random::normal({.mean = m, .sd = s}), p); }
    double check_normal_entropy(double m, double s) { return wsl::random::normal({.mean = m, .sd = s}).entropy(); }
    double check_exponential_sample(double l, wsl::random::xoshiro256pp& g) { return wsl::random::sample(wsl::random::exponential({.rate = l}), g); }
    double check_exponential_fast(double l, wsl::random::xoshiro256pp& g) { return wsl::random::sample(wsl::random::exponential({.rate = l}), g, wsl::random::fast); }
    double check_normal_sample(double m, double s, wsl::random::xoshiro256pp& g) { return wsl::random::sample(wsl::random::normal({.mean = m, .sd = s}), g); }
    float check_normal_sample_float(float m, float s, wsl::random::xoshiro256pp& g) { return wsl::random::sample(wsl::random::normal<float>({.mean = m, .sd = s}), g); }
    double check_normal_fast(double m, double s, wsl::random::xoshiro256pp& g) { return wsl::random::sample(wsl::random::normal({.mean = m, .sd = s}), g, wsl::random::fast); }
    double check_geometric_pmf(double p, long long k) { return wsl::random::pmf(wsl::random::geometric({.p = p}), k); }
    double check_geometric_cdf(double p, long long k) { return wsl::random::cdf(wsl::random::geometric({.p = p}), k); }
    long long check_geometric_quantile(double p, double u) { return wsl::random::quantile(wsl::random::geometric({.p = p}), u); }
    double check_geometric_entropy(double p) { return wsl::random::geometric({.p = p}).entropy(); }
    long long check_geometric_sample(double p, wsl::random::xoshiro256pp& g) { return wsl::random::sample(wsl::random::geometric({.p = p}), g); }
    double check_poisson_pmf(double m, long long k) { return wsl::random::pmf(wsl::random::poisson({.mean = m}), k); }
    double check_poisson_log_pmf(double m, long long k) { return wsl::random::log_pmf(wsl::random::poisson({.mean = m}), k); }
    long long check_poisson_sample(wsl::random::poisson<> const& P, wsl::random::xoshiro256pp& g) { return wsl::random::sample(P, g); }
    void check_poisson_make(double m, wsl::random::poisson<>* out) { *out = wsl::random::poisson({.mean = m}); }
    double check_binomial_pmf(long long n, double p, long long k) { return wsl::random::pmf(wsl::random::binomial({.trials = n, .p = p}), k); }
    double check_binomial_log_pmf(long long n, double p, long long k) { return wsl::random::log_pmf(wsl::random::binomial({.trials = n, .p = p}), k); }
    double check_binomial_variance(long long n, double p) { return wsl::random::binomial({.trials = n, .p = p}).variance(); }
    long long check_binomial_sample(wsl::random::binomial<> const& B, wsl::random::xoshiro256pp& g) { return wsl::random::sample(B, g); }
    void check_binomial_make(long long n, double p, wsl::random::binomial<>* out) { *out = wsl::random::binomial({.trials = n, .p = p}); }
    double check_gamma_pdf(double a, double r, double x) { return wsl::random::pdf(wsl::random::gamma({.shape = a, .rate = r}), x); }
    double check_gamma_cdf(double a, double r, double x) { return wsl::random::cdf(wsl::random::gamma({.shape = a, .rate = r}), x); }
    double check_student_cdf(double n, double t) { return wsl::random::cdf(wsl::random::student_t({.degrees_of_freedom = n}), t); }
    double check_fisher_cdf(double a, double b, double x) { return wsl::random::cdf(wsl::random::fisher_f({.d1 = a, .d2 = b}), x); }
    double check_fisher_pdf(double a, double b, double x) { return wsl::random::pdf(wsl::random::fisher_f({.d1 = a, .d2 = b}), x); }
    double check_fisher_quantile(double a, double b, double p) { return wsl::random::quantile(wsl::random::fisher_f({.d1 = a, .d2 = b}), p); }
    double check_fisher_entropy(double a, double b) { return wsl::random::fisher_f({.d1 = a, .d2 = b}).entropy(); }
    double check_fisher_moments(double a, double b) { wsl::random::fisher_f const F({.d1 = a, .d2 = b}); return F.variance() + F.skewness() + F.excess_kurtosis() + F.mode(); }
    double check_weibull_cdf(double k, double l, double x) { return wsl::random::cdf(wsl::random::weibull({.shape = k, .scale = l}), x); }
    double check_weibull_pdf(double k, double l, double x) { return wsl::random::pdf(wsl::random::weibull({.shape = k, .scale = l}), x); }
    double check_weibull_quantile(double k, double l, double p) { return wsl::random::quantile(wsl::random::weibull({.shape = k, .scale = l}), p); }
    double check_weibull_moments(double k, double l) { wsl::random::weibull const W({.shape = k, .scale = l}); return W.mean() + W.variance() + W.median() + W.mode() + W.skewness() + W.entropy(); }
    double check_laplace_cdf(double m, double b, double x) { return wsl::random::cdf(wsl::random::laplace({.location = m, .scale = b}), x); }
    double check_laplace_quantile(double m, double b, double p) { return wsl::random::quantile(wsl::random::laplace({.location = m, .scale = b}), p); }
    double check_logistic_pdf(double m, double s, double x) { return wsl::random::pdf(wsl::random::logistic({.location = m, .scale = s}), x); }
    double check_logistic_quantile(double m, double s, double p) { return wsl::random::quantile(wsl::random::logistic({.location = m, .scale = s}), p); }
    double check_laplace_variance(double m, double b) { return 1 + wsl::random::laplace({.location = m, .scale = b}).variance(); }
    double check_pareto_cdf(double a, double b, double x) { return wsl::random::cdf(wsl::random::pareto({.scale = a, .shape = b}), x); }
    // Quantiles de survie et algèbre des lois.
    double check_cquantile_gamma(double a, double q) { return wsl::random::cquantile(wsl::random::gamma({.shape = a, .rate = 1.0}), q); }
    double check_cquantile_fisher(double a, double b, double q) { return wsl::random::cquantile(wsl::random::fisher_f({.d1 = a, .d2 = b}), q); }
    double check_cquantile_sum(double a, double q) { return 1 + wsl::random::cquantile(wsl::random::weibull({.shape = a, .scale = 2.0}), q) + wsl::random::cquantile(wsl::random::logistic({.location = a, .scale = 2.0}), q)
                                                            + wsl::random::cquantile(wsl::random::lognormal({.mu = a, .sigma = 2.0}), q) + wsl::random::cquantile(wsl::random::uniform({.low = a, .high = a + 2}), q); }
    double check_affine_closure(double a, double b, double x) { return wsl::random::cdf(a * wsl::random::normal({.mean = 1.0, .sd = 2.0}) + b, x) + (a * wsl::random::laplace({.location = 1.0, .scale = 2.0}) + b).location(); }
    double check_affine_cdf(double a, double b, double x) { return wsl::random::cdf(wsl::random::affine(wsl::random::student_t({.degrees_of_freedom = 3.0}), {.scale = a, .shift = b}), x); }
    double check_affine_quantile(double a, double b, double p) { return 1 + wsl::random::quantile(wsl::random::affine(wsl::random::student_t({.degrees_of_freedom = 3.0}), {.scale = a, .shift = b}), p); }
    double check_affine_moments(double a, double b) { auto const Z = wsl::random::affine(wsl::random::gamma({.shape = 3.0, .rate = 1.0}), {.scale = a, .shift = b}); return Z.mean() + Z.variance() + Z.entropy() + Z.support().upper; }
    double check_truncated_cdf(double a, double b, double x) { auto const T = wsl::random::truncate(wsl::random::normal({.mean = 0.0, .sd = 1.0}), a, b); return wsl::random::cdf(T, x) + wsl::random::ccdf(T, x) + wsl::random::log_pdf(T, x); }
    double check_truncated_quantile(double a, double b, double p) { auto const T = wsl::random::truncate(wsl::random::normal({.mean = 0.0, .sd = 1.0}), a, b); return 1 + wsl::random::quantile(T, p) + wsl::random::cquantile(T, p); }
    double check_mixture_functions(double w, double x) { wsl::random::mixture const M(wsl::random::weighted(w, wsl::random::normal({.mean = -2.0, .sd = 1.0})), wsl::random::weighted(1.0, wsl::random::gamma({.shape = 3.0, .rate = 1.0})));
                                                         return wsl::random::pdf(M, x) + wsl::random::cdf(M, x) + wsl::random::log_pdf(M, x) + wsl::random::quantile(M, 0.25) + M.mean() + M.variance() + M.skewness() + M.excess_kurtosis(); }
    double check_sum_normal(double a, double b) { auto const S = wsl::random::normal({.mean = a, .sd = b}) + wsl::random::normal({.mean = 1.0, .sd = 2.0}); return 1 + S.mean() + S.sd(); }
    double check_sum_moments(double a, double b) { auto const S = wsl::random::gamma({.shape = a, .rate = b}) + wsl::random::student_t({.degrees_of_freedom = 9.0}); return S.mean() + S.variance() + S.skewness() + S.excess_kurtosis(); }
    double check_sum_convolution(double a, double b, long long k) { auto const S = wsl::random::poisson({.mean = a}) + wsl::random::binomial({.trials = 20, .p = b}); return 1 + wsl::random::pmf(S, k) + wsl::random::cdf(S, k) + wsl::random::ccdf(S, k) + double(wsl::random::quantile(S, b)); }
    double check_sum_gamma_family(double a, double b, double x) { auto const S = wsl::random::gamma({.shape = a, .rate = b}) + wsl::random::exponential({.rate = b}); return 1 + wsl::random::cdf(S, x) + S.mean(); }
    double check_independent(double a, double b) { wsl::random::independent const J(wsl::random::normal({.mean = a, .sd = 1.0}), wsl::random::gamma({.shape = b, .rate = 1.0})); return 1 + wsl::random::pdf(J, {a, b}) + wsl::random::log_pdf(J, {a, b}) + wsl::random::cdf(J, {a, b}) + J.entropy(); }
    double check_pushforward(double a, double y) { auto const L = wsl::random::pushforward(wsl::random::normal({.mean = a, .sd = 1.0}), wsl::random::increasing([](double t) { return wsl::special::exp(t); }, [](double t) { return wsl::special::log(t); }, [](double t) { return 1 / t; }));
                                                   return 1 + wsl::random::cdf(L, y) + wsl::random::pdf(L, y) + wsl::random::quantile(L, 0.3); }
    // Lois multivariées.
    double check_dirichlet(double a, double x) { wsl::random::dirichlet<double, 3> const D({.alpha = {a, 2.0, 3.0}}); auto const m = D.mean(); auto const c = D.covariance();
                                                 return 1 + wsl::random::pdf(D, {x, 0.25, 0.75 - x}) + wsl::random::log_pdf(D, {x, 0.25, 0.75 - x}) + D.entropy() + m[0] + c(0, 1) + D.variance()[1]; }
    double check_multinomial(double p, long long k) { wsl::random::multinomial<std::int64_t, 3> const M({.trials = 20, .p = {p, 0.5 - p, 0.5}});
                                                       return 1 + wsl::random::pmf(M, {k, 5, 15 - k}) + M.mean()[0] + M.covariance()(0, 1) + M.variance()[2]; }
    double check_mvn(double a, double x) { wsl::random::multivariate_normal<double, 2> const X({.mean = {a, 0}, .covariance = {{2, 0.5}, {0.5, 1}}});
                                           return 1 + wsl::random::pdf(X, {x, a}) + wsl::random::log_pdf(X, {x, a}) + X.entropy() + X.mahalanobis_squared({x, x}); }
    double check_sphere(double x) { wsl::random::uniform_sphere<double, 3> const S; wsl::random::uniform_ball<double, 5> const B;
                                    return 1 + wsl::random::pdf(S, {x, 0.0, 0.0}) + S.surface_area() + S.entropy() + B.volume() + wsl::random::pdf(B, {x, 0, 0, 0, 0}); }
    double check_mvt(double nu, double x) { wsl::random::multivariate_student_t<double, 2> const X({.degrees_of_freedom = nu, .location = {x, 0}, .scale = {{2, 0.5}, {0.5, 1}}});
                                            return 1 + wsl::random::pdf(X, {x, 1}) + wsl::random::log_pdf(X, {x, x}) + X.entropy() + X.covariance()(0, 1); }
    double check_wishart(double nu, double x) { wsl::random::wishart<double, 2> const W({.degrees_of_freedom = nu, .scale = {{2, 0.5}, {0.5, 1}}}); wsl::random::square_matrix<double, 2> const A{{x, 0.1}, {0.1, 1}};
                                                wsl::random::inverse_wishart<double, 2> const I({.degrees_of_freedom = nu, .scale = {{2, 0.5}, {0.5, 1}}});
                                                return 1 + wsl::random::pdf(W, A) + wsl::random::log_pdf(W, A) + W.entropy() + W.mean()(0, 1) + W.variance()(0, 0) + W.mode()(1, 1) + wsl::random::log_pdf(I, A) + I.entropy() + I.mean()(0, 0); }
    double check_random_matrix(double x) { wsl::random::goe<double, 2> const H; wsl::random::gue<double, 2> const K;
                                           return 1 + wsl::random::pdf(H, wsl::random::square_matrix<double, 2>{{x, 0.5}, {0.5, x}}) + wsl::random::log_pdf(K, wsl::random::square_matrix<std::complex<double>, 2>{{x, {0.5, x}}, {{0.5, -x}, 1}}); }
    double check_von_mises(double m, double k, double x) { wsl::random::von_mises const V({.mu = m, .kappa = k}); return 1 + wsl::random::pdf(V, x) + wsl::random::log_pdf(V, x) + V.entropy() + V.circular_variance() + V.mean_resultant_length(); }
    double check_vmf(double k, double t) { wsl::random::von_mises_fisher<double, 3> const F({.mean_direction = {0, 0.6, 0.8}, .kappa = k});
                                           return 1 + wsl::random::pdf(F, {t, 0.0, 0.0}) + wsl::random::log_pdf(F, {0.0, t, t}) + F.entropy() + F.mean()[1]; }
    double check_permutation_laws(double theta, double q) { wsl::permutation const p{2, 0, 1, 4, 3}; wsl::random::ewens const E({.size = 5, .theta = theta}); wsl::random::mallows const M({.size = 5, .q = q});
                                                            return 1 + wsl::random::pmf(E, p) + wsl::random::log_pmf(M, p) + E.mean_cycles() + M.mean_inversions() + wsl::random::pmf(wsl::random::uniform_permutation({.size = 5}), p); }
    double check_subset_laws(double p) { wsl::random::bernoulli_subset const B({.size = 10, .p = p}); wsl::random::uniform_combination const C({.size = 10, .count = 3});
                                         return 1 + wsl::random::pmf(B, {1, 4}) + wsl::random::log_pmf(C, {1, 4, 7}) + B.mean_size() + C.entropy() + wsl::random::uniform_word({.length = 3, .alphabet = 5}).entropy(); }
    double check_graph_laws(double p) { wsl::random::erdos_renyi const G({.vertices = 6, .p = p}); wsl::random::uniform_graph const H({.vertices = 6, .edges = 3});
                                        return 1 + wsl::random::pmf(G, {{0, 1}, {2, 3}}) + wsl::random::log_pmf(G, {{0, 1}}) + G.mean_edges() + H.entropy() + wsl::random::uniform_tree({.vertices = 9}).entropy()
                                             + wsl::random::uniform_binary_tree({.nodes = 7}).entropy() + wsl::random::uniform_set_partition({.size = 9}).entropy() + wsl::random::uniform_composition({.n = 9, .parts = 3}).entropy(); }
    double check_categorical(double a, double b) { wsl::random::categorical const C({.p = {a, b, 0.5}}); return 1 + wsl::random::pmf(C, 1) + wsl::random::log_pmf(C, 2) + wsl::random::cdf(C, 1) + wsl::random::ccdf(C, 0) + C.mean() + C.variance() + C.entropy(); }
    double check_markov(double a) { wsl::random::markov_chain const M({.transition = {{a, 1 - a}, {0.5, 0.5}}, .start = 0}); return 1 + wsl::random::pmf(M.at(7), 1) + wsl::random::pmf(M.stationary(), 0) + wsl::random::pmf(M.steps(2), {0, 1, 0}); }
    double check_urn_laws(double a, double t) { wsl::random::dirichlet_multinomial const D({.trials = 5, .alpha = {a, 1.5}}); wsl::random::ewens_partition const E({.size = 6, .theta = t, .discount = 0.3});
                                                return 1 + wsl::random::pmf(D, {2, 3}) + D.mean()[0] + D.variance()[1] + wsl::random::pmf(E, {0, 1, 0, 2, 1, 0}) + E.mean_blocks() + wsl::random::pmf(wsl::random::polya_urn({.balls = {a, 2}, .added = t}).steps(2), {0, 1, 1}); }
    // Tirages exacts des processus (sommes paresseuses : verdicts de la voie rapide bit-exacts).
    double check_process_sampling(double a, unsigned long long seed)
    {
        wsl::random::xoshiro256pp g(seed);
        double s = 1;
        for (double v : wsl::random::sample(wsl::random::brownian_motion({.drift = a, .volatility = 1.5}).on_grid({0.5, 1.0}), g)) s += v;
        for (double v : wsl::random::sample(wsl::random::brownian_bridge({.from = {0, a}, .to = {1, 0}}).on_grid({0.5}), g)) s += v;
        for (double v : wsl::random::sample(wsl::random::ornstein_uhlenbeck({.theta = a, .mean = 1.0, .volatility = 0.5, .start = 0.0}).on_grid({0.5, 2.0}), g)) s += v;
        for (double v : wsl::random::sample(wsl::random::poisson_process({.rate = a}).on(2.0), g)) s += v;
        for (double v : wsl::random::sample(wsl::random::random_walk(wsl::random::laplace({.location = 0.0, .scale = a})).steps(3), g)) s += v;
        return s;
    }
    double check_pareto_moments(double a, double b) { wsl::random::pareto const P({.scale = a, .shape = b}); return P.mean() + P.variance() + P.skewness() + P.excess_kurtosis() + P.entropy(); }
    double check_student_quantile(double n, double p) { return wsl::random::quantile(wsl::random::student_t({.degrees_of_freedom = n}), p); }
    double check_student_pdf(double n, double t) { return wsl::random::pdf(wsl::random::student_t({.degrees_of_freedom = n}), t); }
    double check_student_entropy(double n) { return wsl::random::student_t({.degrees_of_freedom = n}).entropy(); }
    double check_student_sample(double n, wsl::random::xoshiro256pp& g) { return wsl::random::sample(wsl::random::student_t({.degrees_of_freedom = n}), g); }
    double check_lognormal_cdf(double m, double s, double x) { return wsl::random::cdf(wsl::random::lognormal({.mu = m, .sigma = s}), x); }
    double check_lognormal_quantile(double m, double s, double p) { return wsl::random::quantile(wsl::random::lognormal({.mu = m, .sigma = s}), p); }
    double check_lognormal_variance(double m, double s) { return wsl::random::lognormal({.mu = m, .sigma = s}).variance(); }
    double check_cauchy_cdf(double a, double b, double x) { return wsl::random::cdf(wsl::random::cauchy({.location = a, .scale = b}), x); }
    double check_cauchy_quantile(double a, double b, double p) { return wsl::random::quantile(wsl::random::cauchy({.location = a, .scale = b}), p); }
    double check_cauchy_sample(double a, double b, wsl::random::xoshiro256pp& g) { return wsl::random::sample(wsl::random::cauchy({.location = a, .scale = b}), g); }
    double check_gamma_quantile(double a, double t, double p) { return wsl::random::quantile(wsl::random::gamma({.shape = a, .scale = t}), p); }
    double check_gamma_entropy(double a, double r) { return wsl::random::gamma({.shape = a, .rate = r}).entropy(); }
    double check_beta_entropy(double a, double b) { return wsl::random::beta({.alpha = a, .beta = b}).entropy(); }
    double check_beta_cdf(double a, double b, double x) { return wsl::random::cdf(wsl::random::beta({.alpha = a, .beta = b}), x); }
    double check_poisson_cdf(double m, long long k) { return wsl::random::cdf(wsl::random::poisson({.mean = m}), k); }
    long long check_poisson_quantile(double m, double u) { return wsl::random::quantile(wsl::random::poisson({.mean = m}), u); }
    long long check_binomial_quantile(long long n, double p, double u) { return wsl::random::quantile(wsl::random::binomial({.trials = n, .p = p}), u); }
    double check_gamma_log_pdf(double a, double t, double x) { return wsl::random::log_pdf(wsl::random::gamma({.shape = a, .scale = t}), x); }
    double check_gamma_sample(wsl::random::gamma<> const& G, wsl::random::xoshiro256pp& g) { return wsl::random::sample(G, g); }
    void check_gamma_make(double a, double r, wsl::random::gamma<>* out) { *out = wsl::random::gamma({.shape = a, .rate = r}); }
    double check_gamma_variance(double a, double t) { return wsl::random::gamma({.shape = a, .scale = t}).variance(); }
    double check_beta_pdf(double a, double b, double x) { return wsl::random::pdf(wsl::random::beta({.alpha = a, .beta = b}), x); }
    double check_beta_log_pdf(double a, double b, double x) { return wsl::random::log_pdf(wsl::random::beta({.alpha = a, .beta = b}), x); }
    double check_beta_sample(wsl::random::beta<> const& B, wsl::random::xoshiro256pp& g) { return wsl::random::sample(B, g); }
    void check_beta_make(double a, double b, wsl::random::beta<>* out) { *out = wsl::random::beta({.alpha = a, .beta = b}); }
    double check_beta_moments(double a, double b) { wsl::random::beta const B({.alpha = a, .beta = b}); return B.variance() + B.excess_kurtosis() + B.skewness() + B.mode(); }
    // Compléments : 1 + f, pour voir aussi les produits rendus.
    double check_gumbel(double a, double b, double x, double p) { wsl::random::gumbel const G({.location = a, .scale = b}); return 1 + wsl::random::pdf(G, x) + wsl::random::log_pdf(G, x) + wsl::random::cdf(G, x) + wsl::random::ccdf(G, x) + wsl::random::quantile(G, p) + wsl::random::cquantile(G, p); }
    double check_gumbel_moments(double a, double b) { wsl::random::gumbel const G({.location = a, .scale = b}); return 1 + G.mean() + G.variance() + G.median() + G.entropy(); }
    double check_frechet(double a, double b, double c, double x, double p) { wsl::random::frechet const F({.shape = a, .scale = b, .location = c}); return 1 + wsl::random::pdf(F, x) + wsl::random::log_pdf(F, x) + wsl::random::cdf(F, x) + wsl::random::ccdf(F, x) + wsl::random::quantile(F, p) + wsl::random::cquantile(F, p); }
    double check_frechet_moments(double a, double b, double c) { wsl::random::frechet const F({.shape = a, .scale = b, .location = c}); return 1 + F.mean() + F.variance() + F.median() + F.mode() + F.entropy(); }
    double check_rayleigh(double a, double x, double p) { wsl::random::rayleigh const R({.scale = a}); return 1 + wsl::random::pdf(R, x) + wsl::random::log_pdf(R, x) + wsl::random::cdf(R, x) + wsl::random::ccdf(R, x) + wsl::random::quantile(R, p) + wsl::random::cquantile(R, p) + R.mean() + R.variance() + R.median() + R.entropy(); }
    double check_triangular(double a, double c, double b, double x, double p) { wsl::random::triangular const T({.low = a, .mode = c, .high = b}); return 1 + wsl::random::pdf(T, x) + wsl::random::log_pdf(T, x) + wsl::random::cdf(T, x) + wsl::random::ccdf(T, x) + wsl::random::quantile(T, p) + wsl::random::cquantile(T, p) + T.mean() + T.variance() + T.median() + T.skewness() + T.entropy(); }
    double check_arcsine(double a, double b, double x, double p) { wsl::random::arcsine const A({.low = a, .high = b}); return 1 + wsl::random::pdf(A, x) + wsl::random::log_pdf(A, x) + wsl::random::cdf(A, x) + wsl::random::ccdf(A, x) + wsl::random::quantile(A, p) + wsl::random::cquantile(A, p) + A.mean() + A.variance() + A.entropy(); }
    double check_inverse_gaussian(double a, double b, double x, double p) { wsl::random::inverse_gaussian const W({.mean = a, .shape = b}); return 1 + wsl::random::pdf(W, x) + wsl::random::log_pdf(W, x) + wsl::random::cdf(W, x) + wsl::random::ccdf(W, x) + wsl::random::quantile(W, p) + wsl::random::cquantile(W, p) + W.variance() + W.mode() + W.skewness() + W.excess_kurtosis(); }
    double check_inverse_gaussian_sample(double a, double b, wsl::random::xoshiro256pp& g) { return 1 + wsl::random::sample(wsl::random::inverse_gaussian({.mean = a, .shape = b}), g); }
    double check_von_mises_distribution(double a, double b, double x, double p) { wsl::random::von_mises const V({.mu = a, .kappa = b}); return 1 + wsl::random::cdf(V, x) + wsl::random::ccdf(V, x) + wsl::random::quantile(V, p) + wsl::random::cquantile(V, p); }
    double check_negative_binomial(double r, double p, long long k, double u) { wsl::random::negative_binomial const X({.successes = r, .p = p}); return 1 + wsl::random::pmf(X, k) + wsl::random::log_pmf(X, k) + wsl::random::cdf(X, k) + wsl::random::ccdf(X, k) + double(wsl::random::quantile(X, u)) + X.mean() + X.variance() + X.skewness() + X.excess_kurtosis(); }
    double check_hypergeometric(long long n, long long k, double u, wsl::random::xoshiro256pp& g) { wsl::random::hypergeometric const X({.population = 1000, .successes = 300, .draws = n}); return 1 + wsl::random::pmf(X, k) + wsl::random::log_pmf(X, k) + wsl::random::cdf(X, k) + wsl::random::ccdf(X, k) + double(wsl::random::quantile(X, u)) + X.mean() + X.variance() + double(wsl::random::sample(X, g)); }
    double check_beta_binomial(double a, double b, long long k, double u, wsl::random::xoshiro256pp& g) { wsl::random::beta_binomial const X({.trials = 30, .alpha = a, .beta = b}); return 1 + wsl::random::pmf(X, k) + wsl::random::log_pmf(X, k) + wsl::random::cdf(X, k) + wsl::random::ccdf(X, k) + double(wsl::random::quantile(X, u)) + X.mean() + X.variance() + double(wsl::random::sample(X, g)); }
    double check_logarithmic(double p, long long k, double u, wsl::random::xoshiro256pp& g) { wsl::random::logarithmic const X({.p = p}); return 1 + wsl::random::pmf(X, k) + wsl::random::log_pmf(X, k) + wsl::random::cdf(X, k) + wsl::random::ccdf(X, k) + double(wsl::random::quantile(X, u)) + X.mean() + X.variance() + double(wsl::random::sample(X, g)); }
    double check_discrete_gaussian(double m, double s, long long k, double u, wsl::random::xoshiro256pp& g) { wsl::random::discrete_gaussian const X({.mu = m, .sigma = s}); return 1 + wsl::random::pmf(X, k) + wsl::random::log_pmf(X, k) + wsl::random::cdf(X, k) + wsl::random::ccdf(X, k) + double(wsl::random::quantile(X, u)) + X.mean() + X.variance() + double(wsl::random::sample(X, g)); }
    double check_skew_normal(double a, double b, double c, double x, double p, wsl::random::xoshiro256pp& g) { wsl::random::skew_normal const X({.location = a, .scale = b, .shape = c}); return 1 + wsl::random::pdf(X, x) + wsl::random::log_pdf(X, x) + wsl::random::cdf(X, x) + wsl::random::ccdf(X, x) + wsl::random::quantile(X, p) + wsl::random::cquantile(X, p) + X.mean() + X.variance() + X.skewness() + X.excess_kurtosis() + wsl::random::sample(X, g); }
    double check_stable(double a, double b, double x, double p, wsl::random::xoshiro256pp& g) { wsl::random::stable const X({.alpha = a, .beta = b, .scale = 2.0, .location = 0.5}); return 1 + wsl::random::pdf(X, x) + wsl::random::log_pdf(X, x) + wsl::random::cdf(X, x) + wsl::random::ccdf(X, x) + wsl::random::quantile(X, p) + wsl::random::sample(X, g); }
    double check_complement_samples(double a, wsl::random::xoshiro256pp& g)
    {
        return 1 + wsl::random::sample(wsl::random::gumbel({.location = a, .scale = 2.0}), g) + wsl::random::sample(wsl::random::frechet({.shape = a, .scale = 2.0}), g) + wsl::random::sample(wsl::random::rayleigh({.scale = a}), g)
             + wsl::random::sample(wsl::random::triangular({.low = 0.0, .mode = a, .high = 5.0}), g) + wsl::random::sample(wsl::random::arcsine({.low = 0.0, .high = a}), g) + wsl::random::sample(wsl::random::erlang({.shape = 3, .rate = a}), g);
    }
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
