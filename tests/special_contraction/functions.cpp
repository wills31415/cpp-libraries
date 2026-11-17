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
    wsl::poisson<> check_poisson_make(double m) { return wsl::poisson({.mean = m}); }
    double check_binomial_pmf(long long n, double p, long long k) { return pmf(wsl::binomial({.trials = n, .p = p}), k); }
    double check_binomial_log_pmf(long long n, double p, long long k) { return log_pmf(wsl::binomial({.trials = n, .p = p}), k); }
    double check_binomial_variance(long long n, double p) { return wsl::binomial({.trials = n, .p = p}).variance(); }
    long long check_binomial_sample(wsl::binomial<> const& B, wsl::xoshiro256pp& g) { return sample(B, g); }
    wsl::binomial<> check_binomial_make(long long n, double p) { return wsl::binomial({.trials = n, .p = p}); }
}
