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
