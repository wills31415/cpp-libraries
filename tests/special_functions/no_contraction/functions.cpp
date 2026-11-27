// Compilé par cmake/check_no_contraction.cmake dans la configuration la plus
// hostile (fma disponible, -ffp-contract=fast) : le code engendré pour ces
// fonctions ne doit contenir aucune fma. Une fma ici veut dire qu'un produit
// entrant dans une addition n'est pas passé par detail::fmul.

#include <wsl/probability.hpp>
#include <wsl/random.hpp>
#include <wsl/special.hpp>

// Les noms de l'ancien wsl::random, réunis (refonte, étape 5a).
namespace wr
{
    using namespace wsl::random;
    using namespace wsl::random::engines;
    using namespace wsl::probability;
    using namespace wsl::probability::laws;
    using namespace wsl::probability::stochastic_processes;
    using namespace wsl::probability::stochastic_processes::processes;
}

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
    double check_uniform_cdf(double a, double b, double x) { return wr::cdf(wr::uniform({.lower = a, .upper = b}), x); }
    double check_uniform_quantile(double a, double b, double p) { return wr::quantile(wr::uniform({.lower = a, .upper = b}), p); }
    double check_uniform_variance(double a, double b) { return wr::uniform({.lower = a, .upper = b}).variance(); }
    double check_discrete_variance(long long a, long long b) { return wr::discrete_uniform({.lower = a, .upper = b}).variance(); }
    double check_exponential_pdf(double l, double x) { return wr::pdf(wr::exponential({.rate = l}), x); }
    double check_exponential_cdf(double l, double x) { return wr::cdf(wr::exponential({.rate = l}), x); }
    double check_exponential_ccdf(double l, double x) { return wr::ccdf(wr::exponential({.rate = l}), x); }
    double check_exponential_log_pdf(double l, double x) { return wr::log_pdf(wr::exponential({.rate = l}), x); }
    double check_exponential_quantile(double l, double p) { return wr::quantile(wr::exponential({.rate = l}), p); }
    double check_exponential_variance(double l) { return wr::exponential({.rate = l}).variance(); }
    double check_normal_pdf(double m, double s, double x) { return wr::pdf(wr::normal({.mean = m, .stddev = s}), x); }
    double check_normal_cdf(double m, double s, double x) { return wr::cdf(wr::normal({.mean = m, .stddev = s}), x); }
    double check_normal_log_pdf(double m, double s, double x) { return wr::log_pdf(wr::normal({.mean = m, .stddev = s}), x); }
    double check_normal_quantile(double m, double s, double p) { return wr::quantile(wr::normal({.mean = m, .stddev = s}), p); }
    double check_normal_entropy(double m, double s) { return wr::normal({.mean = m, .stddev = s}).entropy(); }
    double check_exponential_sample(double l, wr::xoshiro256pp& g) { return wr::sample(wr::exponential({.rate = l}), g); }
    double check_exponential_fast(double l, wr::xoshiro256pp& g) { return wr::sample(wr::exponential({.rate = l}), g, wr::fast); }
    double check_normal_sample(double m, double s, wr::xoshiro256pp& g) { return wr::sample(wr::normal({.mean = m, .stddev = s}), g); }
    float check_normal_sample_float(float m, float s, wr::xoshiro256pp& g) { return wr::sample(wr::normal<float>({.mean = m, .stddev = s}), g); }
    double check_normal_fast(double m, double s, wr::xoshiro256pp& g) { return wr::sample(wr::normal({.mean = m, .stddev = s}), g, wr::fast); }
    double check_geometric_pmf(double p, long long k) { return wr::pmf(wr::geometric({.p = p}), k); }
    double check_geometric_cdf(double p, long long k) { return wr::cdf(wr::geometric({.p = p}), k); }
    long long check_geometric_quantile(double p, double u) { return wr::quantile(wr::geometric({.p = p}), u); }
    double check_geometric_entropy(double p) { return wr::geometric({.p = p}).entropy(); }
    long long check_geometric_sample(double p, wr::xoshiro256pp& g) { return wr::sample(wr::geometric({.p = p}), g); }
    double check_poisson_pmf(double m, long long k) { return wr::pmf(wr::poisson({.mean = m}), k); }
    double check_poisson_log_pmf(double m, long long k) { return wr::log_pmf(wr::poisson({.mean = m}), k); }
    long long check_poisson_sample(wr::poisson<> const& P, wr::xoshiro256pp& g) { return wr::sample(P, g); }
    void check_poisson_make(double m, wr::poisson<>* out) { *out = wr::poisson({.mean = m}); }
    double check_binomial_pmf(long long n, double p, long long k) { return wr::pmf(wr::binomial({.trials = n, .p = p}), k); }
    double check_binomial_log_pmf(long long n, double p, long long k) { return wr::log_pmf(wr::binomial({.trials = n, .p = p}), k); }
    double check_binomial_variance(long long n, double p) { return wr::binomial({.trials = n, .p = p}).variance(); }
    long long check_binomial_sample(wr::binomial<> const& B, wr::xoshiro256pp& g) { return wr::sample(B, g); }
    void check_binomial_make(long long n, double p, wr::binomial<>* out) { *out = wr::binomial({.trials = n, .p = p}); }
    double check_gamma_pdf(double a, double r, double x) { return wr::pdf(wr::gamma({.rate = r, .shape = a}), x); }
    double check_gamma_cdf(double a, double r, double x) { return wr::cdf(wr::gamma({.rate = r, .shape = a}), x); }
    double check_student_cdf(double n, double t) { return wr::cdf(wr::student_t({.degrees_of_freedom = n}), t); }
    double check_fisher_cdf(double a, double b, double x) { return wr::cdf(wr::fisher_f({.d1 = a, .d2 = b}), x); }
    double check_fisher_pdf(double a, double b, double x) { return wr::pdf(wr::fisher_f({.d1 = a, .d2 = b}), x); }
    double check_fisher_quantile(double a, double b, double p) { return wr::quantile(wr::fisher_f({.d1 = a, .d2 = b}), p); }
    double check_fisher_entropy(double a, double b) { return wr::fisher_f({.d1 = a, .d2 = b}).entropy(); }
    double check_fisher_moments(double a, double b) { wr::fisher_f const F({.d1 = a, .d2 = b}); return F.variance() + F.skewness() + F.excess_kurtosis() + F.mode(); }
    double check_weibull_cdf(double k, double l, double x) { return wr::cdf(wr::weibull({.scale = l, .shape = k}), x); }
    double check_weibull_pdf(double k, double l, double x) { return wr::pdf(wr::weibull({.scale = l, .shape = k}), x); }
    double check_weibull_quantile(double k, double l, double p) { return wr::quantile(wr::weibull({.scale = l, .shape = k}), p); }
    double check_weibull_moments(double k, double l) { wr::weibull const W({.scale = l, .shape = k}); return W.mean() + W.variance() + W.median() + W.mode() + W.skewness() + W.entropy(); }
    double check_laplace_cdf(double m, double b, double x) { return wr::cdf(wr::laplace({.location = m, .scale = b}), x); }
    double check_laplace_quantile(double m, double b, double p) { return wr::quantile(wr::laplace({.location = m, .scale = b}), p); }
    double check_logistic_pdf(double m, double s, double x) { return wr::pdf(wr::logistic({.location = m, .scale = s}), x); }
    double check_logistic_quantile(double m, double s, double p) { return wr::quantile(wr::logistic({.location = m, .scale = s}), p); }
    double check_laplace_variance(double m, double b) { return 1 + wr::laplace({.location = m, .scale = b}).variance(); }
    double check_pareto_cdf(double a, double b, double x) { return wr::cdf(wr::pareto({.scale = a, .shape = b}), x); }
    // Quantiles de survie et algèbre des lois.
    double check_cquantile_gamma(double a, double q) { return wr::cquantile(wr::gamma({.rate = 1.0, .shape = a}), q); }
    double check_cquantile_fisher(double a, double b, double q) { return wr::cquantile(wr::fisher_f({.d1 = a, .d2 = b}), q); }
    double check_cquantile_sum(double a, double q) { return 1 + wr::cquantile(wr::weibull({.scale = 2.0, .shape = a}), q) + wr::cquantile(wr::logistic({.location = a, .scale = 2.0}), q)
                                                            + wr::cquantile(wr::lognormal({.log_mean = a, .log_stddev = 2.0}), q) + wr::cquantile(wr::uniform({.lower = a, .upper = a + 2}), q); }
    double check_affine_closure(double a, double b, double x) { return wr::cdf(wr::affine(wr::normal({.mean = 1.0, .stddev = 2.0}), {.scale = a, .shift = b}), x) + wr::affine(wr::laplace({.location = 1.0, .scale = 2.0}), {.scale = a, .shift = b}).location(); }
    double check_affine_cdf(double a, double b, double x) { return wr::cdf(wr::affine(wr::student_t({.degrees_of_freedom = 3.0}), {.scale = a, .shift = b}), x); }
    double check_affine_quantile(double a, double b, double p) { return 1 + wr::quantile(wr::affine(wr::student_t({.degrees_of_freedom = 3.0}), {.scale = a, .shift = b}), p); }
    double check_affine_moments(double a, double b) { auto const Z = wr::affine(wr::gamma({.rate = 1.0, .shape = 3.0}), {.scale = a, .shift = b}); return Z.mean() + Z.variance() + Z.entropy() + Z.support().upper; }
    double check_truncated_cdf(double a, double b, double x) { auto const T = wr::truncate(wr::normal({.mean = 0.0, .stddev = 1.0}), {.lower = a, .upper = b}); return wr::cdf(T, x) + wr::ccdf(T, x) + wr::log_pdf(T, x); }
    double check_truncated_quantile(double a, double b, double p) { auto const T = wr::truncate(wr::normal({.mean = 0.0, .stddev = 1.0}), {.lower = a, .upper = b}); return 1 + wr::quantile(T, p) + wr::cquantile(T, p); }
    double check_mixture_functions(double w, double x) { wr::mixture const M(wr::weighted(w, wr::normal({.mean = -2.0, .stddev = 1.0})), wr::weighted(1.0, wr::gamma({.rate = 1.0, .shape = 3.0})));
                                                         return wr::pdf(M, x) + wr::cdf(M, x) + wr::log_pdf(M, x) + wr::quantile(M, 0.25) + M.mean() + M.variance() + M.skewness() + M.excess_kurtosis(); }
    double check_sum_normal(double a, double b) { auto const S = wr::convolution(wr::normal({.mean = a, .stddev = b}), wr::normal({.mean = 1.0, .stddev = 2.0})); return 1 + S.mean() + S.stddev(); }
    double check_sum_moments(double a, double b) { auto const S = wr::convolution(wr::gamma({.rate = b, .shape = a}), wr::student_t({.degrees_of_freedom = 9.0})); return S.mean() + S.variance() + S.skewness() + S.excess_kurtosis(); }
    double check_sum_convolution(double a, double b, long long k) { auto const S = wr::convolution(wr::poisson({.mean = a}), wr::binomial({.trials = 20, .p = b})); return 1 + wr::pmf(S, k) + wr::cdf(S, k) + wr::ccdf(S, k) + double(wr::quantile(S, b)); }
    double check_sum_gamma_family(double a, double b, double x) { auto const S = wr::convolution(wr::gamma({.rate = b, .shape = a}), wr::exponential({.rate = b})); return 1 + wr::cdf(S, x) + S.mean(); }
    double check_independent(double a, double b) { wr::product_measure const J(wr::normal({.mean = a, .stddev = 1.0}), wr::gamma({.rate = 1.0, .shape = b})); return 1 + wr::pdf(J, {a, b}) + wr::log_pdf(J, {a, b}) + wr::cdf(J, {a, b}) + J.entropy(); }
    double check_pushforward(double a, double y) { auto const L = wr::pushforward(wr::normal({.mean = a, .stddev = 1.0}), wr::increasing([](double t) { return wsl::special::exp(t); }, [](double t) { return wsl::special::log(t); }, [](double t) { return 1 / t; }));
                                                   return 1 + wr::cdf(L, y) + wr::pdf(L, y) + wr::quantile(L, 0.3); }
    // Lois multivariées.
    double check_dirichlet(double a, double x) { wr::dirichlet<double, 3> const D({.alpha = {a, 2.0, 3.0}}); auto const m = D.mean(); auto const c = D.covariance();
                                                 return 1 + wr::pdf(D, {x, 0.25, 0.75 - x}) + wr::log_pdf(D, {x, 0.25, 0.75 - x}) + D.entropy() + m[0] + c[0, 1] + D.variance()[1]; }
    double check_multinomial(double p, long long k) { wr::multinomial<std::int64_t, 3> const M({.trials = 20, .p = {p, 0.5 - p, 0.5}});
                                                       return 1 + wr::pmf(M, {k, 5, 15 - k}) + M.mean()[0] + M.covariance()[0, 1] + M.variance()[2]; }
    double check_mvn(double a, double x) { wr::multivariate_normal<double, 2> const X({.mean = {a, 0}, .covariance = {{2, 0.5}, {0.5, 1}}});
                                           return 1 + wr::pdf(X, {x, a}) + wr::log_pdf(X, {x, a}) + X.entropy() + X.mahalanobis_squared({x, x}); }
    double check_sphere(double x) { wr::uniform_sphere<double, 3> const S; wr::uniform_ball<double, 5> const B;
                                    return 1 + wr::pdf(S, {x, 0.0, 0.0}) + S.surface_area() + S.entropy() + B.volume() + wr::pdf(B, {x, 0, 0, 0, 0}); }
    double check_mvt(double nu, double x) { wr::multivariate_student_t<double, 2> const X({.location = {x, 0}, .scale = {{2, 0.5}, {0.5, 1}}, .degrees_of_freedom = nu});
                                            return 1 + wr::pdf(X, {x, 1}) + wr::log_pdf(X, {x, x}) + X.entropy() + X.covariance()[0, 1]; }
    double check_wishart(double nu, double x) { wr::wishart<double, 2> const W({.scale = {{2, 0.5}, {0.5, 1}}, .degrees_of_freedom = nu}); wr::square_matrix<double, 2> const A{{x, 0.1}, {0.1, 1}};
                                                wr::inverse_wishart<double, 2> const I({.scale = {{2, 0.5}, {0.5, 1}}, .degrees_of_freedom = nu});
                                                return 1 + wr::pdf(W, A) + wr::log_pdf(W, A) + W.entropy() + W.mean()[0, 1] + W.variance()[0, 0] + W.mode()[1, 1] + wr::log_pdf(I, A) + I.entropy() + I.mean()[0, 0]; }
    double check_random_matrix(double x) { wr::goe<double, 2> const H; wr::gue<double, 2> const K;
                                           return 1 + wr::pdf(H, wr::square_matrix<double, 2>{{x, 0.5}, {0.5, x}}) + wr::log_pdf(K, wr::square_matrix<std::complex<double>, 2>{{x, {0.5, x}}, {{0.5, -x}, 1}}); }
    double check_von_mises(double m, double k, double x) { wr::von_mises const V({.location = m, .concentration = k}); return 1 + wr::pdf(V, x) + wr::log_pdf(V, x) + V.entropy() + V.circular_variance() + V.mean_resultant_length(); }
    double check_vmf(double k, double t) { wr::von_mises_fisher<double, 3> const F({.mean_direction = {0, 0.6, 0.8}, .concentration = k});
                                           return 1 + wr::pdf(F, {t, 0.0, 0.0}) + wr::log_pdf(F, {0.0, t, t}) + F.entropy() + F.mean()[1]; }
    double check_permutation_laws(double theta, double q) { wsl::permutation const p{2, 0, 1, 4, 3}; wr::ewens const E({.size = 5, .theta = theta}); wr::mallows const M({.size = 5, .q = q});
                                                            return 1 + wr::pmf(E, p) + wr::log_pmf(M, p) + E.mean_cycles() + M.mean_inversions() + wr::pmf(wr::uniform_permutation({.size = 5}), p); }
    double check_subset_laws(double p) { wr::bernoulli_subset const B({.size = 10, .p = p}); wr::uniform_combination const C({.size = 10, .count = 3});
                                         return 1 + wr::pmf(B, {1, 4}) + wr::log_pmf(C, {1, 4, 7}) + B.mean_size() + C.entropy() + wr::uniform_word({.length = 3, .alphabet = 5}).entropy(); }
    double check_graph_laws(double p) { wr::erdos_renyi const G({.vertices = 6, .p = p}); wr::uniform_graph const H({.vertices = 6, .edges = 3});
                                        return 1 + wr::pmf(G, {{0, 1}, {2, 3}}) + wr::log_pmf(G, {{0, 1}}) + G.mean_edges() + H.entropy() + wr::uniform_tree({.vertices = 9}).entropy()
                                             + wr::uniform_binary_tree({.vertices = 7}).entropy() + wr::uniform_set_partition({.size = 9}).entropy() + wr::uniform_composition({.total = 9, .parts = 3}).entropy(); }
    double check_categorical(double a, double b) { wr::categorical const C({.p = {a, b, 0.5}}); return 1 + wr::pmf(C, 1) + wr::log_pmf(C, 2) + wr::cdf(C, 1) + wr::ccdf(C, 0) + C.mean() + C.variance() + C.entropy(); }
    double check_markov(double a) { wr::markov_chain const M({.transition = {{a, 1 - a}, {0.5, 0.5}}, .start = 0}); return 1 + wr::pmf(M.at(7), 1) + wr::pmf(M.stationary(), 0) + wr::pmf(M.steps(2), {0, 1, 0}); }
    double check_urn_laws(double a, double t) { wr::dirichlet_multinomial const D({.trials = 5, .alpha = {a, 1.5}}); wr::ewens_partition const E({.size = 6, .theta = t, .discount = 0.3});
                                                return 1 + wr::pmf(D, {2, 3}) + D.mean()[0] + D.variance()[1] + wr::pmf(E, {0, 1, 0, 2, 1, 0}) + E.mean_blocks() + wr::pmf(wr::polya_urn({.balls = {a, 2}, .added = t}).steps(2), {0, 1, 1}); }
    // Tirages exacts des processus (sommes paresseuses : verdicts de la voie rapide bit-exacts).
    double check_process_sampling(double a, unsigned long long seed)
    {
        wr::xoshiro256pp g(seed);
        double s = 1;
        for (double v : wr::sample(wr::brownian_motion({.drift = a, .volatility = 1.5}).on_grid({0.5, 1.0}), g)) s += v;
        for (double v : wr::sample(wr::brownian_bridge({.from = {0, a}, .to = {1, 0}}).on_grid({0.5}), g)) s += v;
        for (double v : wr::sample(wr::ornstein_uhlenbeck({.theta = a, .mean = 1.0, .volatility = 0.5, .start = 0.0}).on_grid({0.5, 2.0}), g)) s += v;
        for (double v : wr::sample(wr::poisson_process({.rate = a}).on(2.0), g)) s += v;
        for (double v : wr::sample(wr::random_walk(wr::laplace({.location = 0.0, .scale = a})).steps(3), g)) s += v;
        return s;
    }
    double check_pareto_moments(double a, double b) { wr::pareto const P({.scale = a, .shape = b}); return P.mean() + P.variance() + P.skewness() + P.excess_kurtosis() + P.entropy(); }
    double check_student_quantile(double n, double p) { return wr::quantile(wr::student_t({.degrees_of_freedom = n}), p); }
    double check_student_pdf(double n, double t) { return wr::pdf(wr::student_t({.degrees_of_freedom = n}), t); }
    double check_student_entropy(double n) { return wr::student_t({.degrees_of_freedom = n}).entropy(); }
    double check_student_sample(double n, wr::xoshiro256pp& g) { return wr::sample(wr::student_t({.degrees_of_freedom = n}), g); }
    double check_lognormal_cdf(double m, double s, double x) { return wr::cdf(wr::lognormal({.log_mean = m, .log_stddev = s}), x); }
    double check_lognormal_quantile(double m, double s, double p) { return wr::quantile(wr::lognormal({.log_mean = m, .log_stddev = s}), p); }
    double check_lognormal_variance(double m, double s) { return wr::lognormal({.log_mean = m, .log_stddev = s}).variance(); }
    double check_cauchy_cdf(double a, double b, double x) { return wr::cdf(wr::cauchy({.location = a, .scale = b}), x); }
    double check_cauchy_quantile(double a, double b, double p) { return wr::quantile(wr::cauchy({.location = a, .scale = b}), p); }
    double check_cauchy_sample(double a, double b, wr::xoshiro256pp& g) { return wr::sample(wr::cauchy({.location = a, .scale = b}), g); }
    double check_gamma_quantile(double a, double t, double p) { return wr::quantile(wr::gamma({.scale = t, .shape = a}), p); }
    double check_gamma_entropy(double a, double r) { return wr::gamma({.rate = r, .shape = a}).entropy(); }
    double check_beta_entropy(double a, double b) { return wr::beta({.alpha = a, .beta = b}).entropy(); }
    double check_beta_cdf(double a, double b, double x) { return wr::cdf(wr::beta({.alpha = a, .beta = b}), x); }
    double check_poisson_cdf(double m, long long k) { return wr::cdf(wr::poisson({.mean = m}), k); }
    long long check_poisson_quantile(double m, double u) { return wr::quantile(wr::poisson({.mean = m}), u); }
    long long check_binomial_quantile(long long n, double p, double u) { return wr::quantile(wr::binomial({.trials = n, .p = p}), u); }
    double check_gamma_log_pdf(double a, double t, double x) { return wr::log_pdf(wr::gamma({.scale = t, .shape = a}), x); }
    double check_gamma_sample(wr::gamma<> const& G, wr::xoshiro256pp& g) { return wr::sample(G, g); }
    void check_gamma_make(double a, double r, wr::gamma<>* out) { *out = wr::gamma({.rate = r, .shape = a}); }
    double check_gamma_variance(double a, double t) { return wr::gamma({.scale = t, .shape = a}).variance(); }
    double check_beta_pdf(double a, double b, double x) { return wr::pdf(wr::beta({.alpha = a, .beta = b}), x); }
    double check_beta_log_pdf(double a, double b, double x) { return wr::log_pdf(wr::beta({.alpha = a, .beta = b}), x); }
    double check_beta_sample(wr::beta<> const& B, wr::xoshiro256pp& g) { return wr::sample(B, g); }
    void check_beta_make(double a, double b, wr::beta<>* out) { *out = wr::beta({.alpha = a, .beta = b}); }
    double check_beta_moments(double a, double b) { wr::beta const B({.alpha = a, .beta = b}); return B.variance() + B.excess_kurtosis() + B.skewness() + B.mode(); }
    // Compléments : 1 + f, pour voir aussi les produits rendus.
    double check_gumbel(double a, double b, double x, double p) { wr::gumbel const G({.location = a, .scale = b}); return 1 + wr::pdf(G, x) + wr::log_pdf(G, x) + wr::cdf(G, x) + wr::ccdf(G, x) + wr::quantile(G, p) + wr::cquantile(G, p); }
    double check_gumbel_moments(double a, double b) { wr::gumbel const G({.location = a, .scale = b}); return 1 + G.mean() + G.variance() + G.median() + G.entropy(); }
    double check_frechet(double a, double b, double c, double x, double p) { wr::frechet const F({.location = c, .scale = b, .shape = a}); return 1 + wr::pdf(F, x) + wr::log_pdf(F, x) + wr::cdf(F, x) + wr::ccdf(F, x) + wr::quantile(F, p) + wr::cquantile(F, p); }
    double check_frechet_moments(double a, double b, double c) { wr::frechet const F({.location = c, .scale = b, .shape = a}); return 1 + F.mean() + F.variance() + F.median() + F.mode() + F.entropy(); }
    double check_rayleigh(double a, double x, double p) { wr::rayleigh const R({.scale = a}); return 1 + wr::pdf(R, x) + wr::log_pdf(R, x) + wr::cdf(R, x) + wr::ccdf(R, x) + wr::quantile(R, p) + wr::cquantile(R, p) + R.mean() + R.variance() + R.median() + R.entropy(); }
    double check_triangular(double a, double c, double b, double x, double p) { wr::triangular const T({.lower = a, .mode = c, .upper = b}); return 1 + wr::pdf(T, x) + wr::log_pdf(T, x) + wr::cdf(T, x) + wr::ccdf(T, x) + wr::quantile(T, p) + wr::cquantile(T, p) + T.mean() + T.variance() + T.median() + T.skewness() + T.entropy(); }
    double check_arcsine(double a, double b, double x, double p) { wr::arcsine const A({.lower = a, .upper = b}); return 1 + wr::pdf(A, x) + wr::log_pdf(A, x) + wr::cdf(A, x) + wr::ccdf(A, x) + wr::quantile(A, p) + wr::cquantile(A, p) + A.mean() + A.variance() + A.entropy(); }
    double check_inverse_gaussian(double a, double b, double x, double p) { wr::inverse_gaussian const W({.mean = a, .shape = b}); return 1 + wr::pdf(W, x) + wr::log_pdf(W, x) + wr::cdf(W, x) + wr::ccdf(W, x) + wr::quantile(W, p) + wr::cquantile(W, p) + W.variance() + W.mode() + W.skewness() + W.excess_kurtosis(); }
    double check_inverse_gaussian_sample(double a, double b, wr::xoshiro256pp& g) { return 1 + wr::sample(wr::inverse_gaussian({.mean = a, .shape = b}), g); }
    double check_von_mises_distribution(double a, double b, double x, double p) { wr::von_mises const V({.location = a, .concentration = b}); return 1 + wr::cdf(V, x) + wr::ccdf(V, x) + wr::quantile(V, p) + wr::cquantile(V, p); }
    double check_negative_binomial(double r, double p, long long k, double u) { wr::negative_binomial const X({.successes = r, .p = p}); return 1 + wr::pmf(X, k) + wr::log_pmf(X, k) + wr::cdf(X, k) + wr::ccdf(X, k) + double(wr::quantile(X, u)) + X.mean() + X.variance() + X.skewness() + X.excess_kurtosis(); }
    double check_hypergeometric(long long n, long long k, double u, wr::xoshiro256pp& g) { wr::hypergeometric const X({.population = 1000, .successes = 300, .draws = n}); return 1 + wr::pmf(X, k) + wr::log_pmf(X, k) + wr::cdf(X, k) + wr::ccdf(X, k) + double(wr::quantile(X, u)) + X.mean() + X.variance() + double(wr::sample(X, g)); }
    double check_beta_binomial(double a, double b, long long k, double u, wr::xoshiro256pp& g) { wr::beta_binomial const X({.trials = 30, .alpha = a, .beta = b}); return 1 + wr::pmf(X, k) + wr::log_pmf(X, k) + wr::cdf(X, k) + wr::ccdf(X, k) + double(wr::quantile(X, u)) + X.mean() + X.variance() + double(wr::sample(X, g)); }
    double check_logarithmic(double p, long long k, double u, wr::xoshiro256pp& g) { wr::logarithmic const X({.p = p}); return 1 + wr::pmf(X, k) + wr::log_pmf(X, k) + wr::cdf(X, k) + wr::ccdf(X, k) + double(wr::quantile(X, u)) + X.mean() + X.variance() + double(wr::sample(X, g)); }
    double check_discrete_gaussian(double m, double s, long long k, double u, wr::xoshiro256pp& g) { wr::discrete_gaussian const X({.location = m, .scale = s}); return 1 + wr::pmf(X, k) + wr::log_pmf(X, k) + wr::cdf(X, k) + wr::ccdf(X, k) + double(wr::quantile(X, u)) + X.mean() + X.variance() + double(wr::sample(X, g)); }
    double check_skellam(double a, double b, long long k, double u, wr::xoshiro256pp& g) { wr::skellam const X({.mean1 = a, .mean2 = b}); return 1 + wr::pmf(X, k) + wr::log_pmf(X, k) + wr::cdf(X, k) + wr::ccdf(X, k) + double(wr::quantile(X, u)) + X.skewness() + X.excess_kurtosis() + double(wr::sample(X, g)); }
    double check_zipf(double s, long long k, double u, wr::xoshiro256pp& g) { wr::zipf const X({.size = 1000, .exponent = s}); return 1 + wr::pmf(X, k) + wr::log_pmf(X, k) + wr::cdf(X, k) + wr::ccdf(X, k) + double(wr::quantile(X, u)) + X.mean() + X.variance() + double(wr::sample(X, g)); }
    double check_rice(double a, double b, double x, double p, wr::xoshiro256pp& g) { wr::rice const X({.scale = a, .noncentrality = b}); return 1 + wr::pdf(X, x) + wr::log_pdf(X, x) + wr::cdf(X, x) + wr::ccdf(X, x) + wr::quantile(X, p) + wr::cquantile(X, p) + X.mean() + X.variance() + wr::sample(X, g); }
    double check_skew_normal(double a, double b, double c, double x, double p, wr::xoshiro256pp& g) { wr::skew_normal const X({.location = a, .scale = b, .shape = c}); return 1 + wr::pdf(X, x) + wr::log_pdf(X, x) + wr::cdf(X, x) + wr::ccdf(X, x) + wr::quantile(X, p) + wr::cquantile(X, p) + X.mean() + X.variance() + X.skewness() + X.excess_kurtosis() + wr::sample(X, g); }
    double check_stable(double a, double b, double x, double p, wr::xoshiro256pp& g) { wr::stable const X({.location = 0.5, .scale = 2.0, .alpha = a, .beta = b}); return 1 + wr::pdf(X, x) + wr::log_pdf(X, x) + wr::cdf(X, x) + wr::ccdf(X, x) + wr::quantile(X, p) + wr::sample(X, g); }
    double check_complement_samples(double a, wr::xoshiro256pp& g)
    {
        return 1 + wr::sample(wr::gumbel({.location = a, .scale = 2.0}), g) + wr::sample(wr::frechet({.scale = 2.0, .shape = a}), g) + wr::sample(wr::rayleigh({.scale = a}), g)
             + wr::sample(wr::triangular({.lower = 0.0, .mode = a, .upper = 5.0}), g) + wr::sample(wr::arcsine({.lower = 0.0, .upper = a}), g) + wr::sample(wr::erlang({.rate = a, .shape = 3}), g);
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
