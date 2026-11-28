// Étiquettes : precision long
// Précision des lois de wsl::random : erreur en ulp contre les valeurs de
// référence de mpmath (tools/probability/law_references.py, engendrées par CMake).

#include <wsl/probability.hpp>
#include <wsl/random.hpp>

#include "data/random_law_reference.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <iterator>
#include <span>

// Les noms de l'ancien wsl::random, réunis (refonte, étape 5a) : moteurs, primitives, lois,
// processus, mesures de référence.
namespace wr
{
    using namespace wsl::random;
    using namespace wsl::random::engines;
    using namespace wsl::probability;
    using namespace wsl::probability::laws;
    using namespace wsl::probability::stochastic_processes;
    using namespace wsl::probability::stochastic_processes::processes;
    using namespace wsl::measure_theory;
    namespace detail = wsl::probability::detail;
}


namespace
{
    using law_reference_data::point;

    template <class T>
    T ulp(T x)
    {
        T const a = std::abs(x);
        T const next = std::nextafter(a, std::numeric_limits<T>::infinity());
        return next - a;
    }

    bool failed = false;

    // Erreur de f(a, b, x) en ulp de l'arrondi correct : |f − valeur exacte| / ulp.
    template <class T, class F>
    void check(char const* name, std::span<point<T> const> points, F f, double bound)
    {
        double worst = 0, at = 0;
        std::size_t correct = 0;
        for (auto const& r : points)
        {
            T const y = f(r.a, r.b, r.x);
            double const error = std::abs(double(y - r.rounded) / double(ulp(r.rounded)) - r.offset);
            if (!(error <= worst))
            {
                worst = error;
                at = r.x;
            }
            correct += y == r.rounded;
        }
        bool const ok = worst <= bound;
        std::printf("%-28s %5zu points : erreur max %.4f ulp (x = %a), arrondi correct %.2f %%%s\n", name, points.size(),
                    worst, at, 100.0 * double(correct) / double(points.size()), ok ? "" : "  ÉCHEC");
        std::fflush(stdout);
        failed |= !ok;
    }

    template <class T>
    wr::exponential<T> exponential(double rate, double) { return wr::exponential<T>({.rate = T(rate)}); }
    template <class T>
    wr::normal<T> normal(double mean, double sd) { return wr::normal<T>({.mean = T(mean), .stddev = T(sd)}); }
}

// Normale asymétrique : (α, ξ) portés par (a, b), ω retrouvé dans la table des références.
template <class T>
wr::skew_normal<T> skew_of(double a, double xi)
{
    struct { double a, xi, omega; } const table[] = {{4.0, 0.5, 2.0}, {1.0, 0.0, 1.0}, {-3.0, -1.0, 0.5}, {50.0, 0.0, 1.0}, {0.3, 2.0, 3.0}};
    for (auto const& t : table)
        if (float(t.a) == float(a) && float(t.xi) == float(xi))
            return wr::skew_normal<T>({.location = T(t.xi), .scale = T(t.omega), .shape = T(t.a)});
    return wr::skew_normal<T>({.location = T(xi), .scale = T(1), .shape = T(a)});
}

// Bêta-binomiale : (n, α) portés par (a, b), β retrouvé dans la table des références.
template <class Int>
wr::beta_binomial<Int> beta_binomial_of(double n, double alpha)
{
    struct { long long n; double alpha, beta; } const table[] = {{20, 2.0, 3.0}, {1000, 0.5, 0.5}, {1000000, 3.0, 40.0}, {50, 0.05, 2.0}};
    for (auto const& t : table)
        if (double(t.n) == n && t.alpha == alpha)
            return wr::beta_binomial<Int>({.trials = Int(t.n), .alpha = t.alpha, .beta = t.beta});
    return wr::beta_binomial<Int>({.trials = Int(n), .alpha = alpha, .beta = 1.0});
}

// Hypergéométrique : (N, K) portés par (a, b), n retrouvé dans la table des références.
template <class Int>
wr::hypergeometric<Int> hypergeometric_of(double N, double Ks)
{
    struct { long long N, K, n; } const table[] = {{50, 20, 10}, {1000, 300, 400}, {1000000, 1000, 50000}, {1000000000, 500000000, 1000000}, {30, 29, 15}};
    for (auto const& t : table)
        if (double(t.N) == N && double(t.K) == Ks)
            return wr::hypergeometric<Int>({.population = Int(t.N), .successes = Int(t.K), .draws = Int(t.n)});
    return wr::hypergeometric<Int>({.population = Int(N), .successes = Int(Ks), .draws = 0});
}

// Paramètres des références des lois de compléments (a, b ne portent que deux paramètres).
double frechet_location(double shape, double scale)
{
    struct { double shape, scale, location; } const table[] = {{3.0, 2.0, 0.5}, {0.5, 1.0, 0.0}, {1.5, 3.0, -2.0}, {50.0, 1.0, 0.0}, {2.5, 1e-5, 1.0}};
    for (auto const& t : table)
        if (float(t.shape) == float(shape) && float(t.scale) == float(scale))
            return t.location;
    return 0;
}
template <class T>
wr::triangular<T> triangle(double low, double high)
{
    struct { double low, mode, high; } const table[] = {{0.0, 1.0, 3.0}, {-1.0, -1.0, 2.0}, {0.0, 5.0, 5.0}, {1e3, 1e3 + 0.25, 1e3 + 1.0}};
    for (auto const& t : table)
        if (float(t.low) == float(low) && float(t.high) == float(high))
            return wr::triangular<T>({.lower = T(t.low), .mode = T(t.mode), .upper = T(t.high)});
    return wr::triangular<T>({.lower = T(low), .mode = T(low), .upper = T(high)});
}

int main()
{
    namespace d = law_reference_data;
    // En double : un seul arrondi final, sur un résultat intermédiaire précis à
    // 2^-60 près ou mieux ; log_pdf perd en relatif quand ses termes se compensent.
    // En float : le résultat double arrondi en float.
    constexpr double tight = 0.51, loose = 0.5 + 0x1p-20;
#define LAW_CHECKS(T, bound)                                                                                                  \
    check<T>("exponential pdf " #T, d::exponential_pdf_##T, [](double a, double b, double x) { return wr::pdf(exponential<T>(a, b), T(x)); }, bound); \
    check<T>("exponential log_pdf " #T, d::exponential_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(exponential<T>(a, b), T(x)); }, 1.0); \
    check<T>("exponential cdf " #T, d::exponential_cdf_##T, [](double a, double b, double x) { return wr::cdf(exponential<T>(a, b), T(x)); }, bound); \
    check<T>("exponential ccdf " #T, d::exponential_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(exponential<T>(a, b), T(x)); }, bound); \
    check<T>("exponential quantile " #T, d::exponential_quantile_##T, [](double a, double b, double p) { return wr::quantile(exponential<T>(a, b), p); }, bound); \
    check<T>("normal pdf " #T, d::normal_pdf_##T, [](double a, double b, double x) { return wr::pdf(normal<T>(a, b), T(x)); }, bound); \
    check<T>("normal log_pdf " #T, d::normal_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(normal<T>(a, b), T(x)); }, 1.0); \
    check<T>("normal cdf " #T, d::normal_cdf_##T, [](double a, double b, double x) { return wr::cdf(normal<T>(a, b), T(x)); }, bound); \
    check<T>("normal ccdf " #T, d::normal_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(normal<T>(a, b), T(x)); }, bound); \
    check<T>("normal quantile " #T, d::normal_quantile_##T, [](double a, double b, double p) { return wr::quantile(normal<T>(a, b), p); }, bound);
#define GAMMA_BETA_CHECKS(T, bound)                                                                                           \
    check<T>("gamma pdf " #T, d::gamma_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::gamma<T>({.rate = T(b), .shape = T(a)}), T(x)); }, bound); \
    check<T>("gamma log_pdf " #T, d::gamma_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::gamma<T>({.rate = T(b), .shape = T(a)}), T(x)); }, 1.0); \
    check<T>("gamma (échelle) pdf " #T, d::gamma_scale_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::gamma<T>({.scale = T(b), .shape = T(a)}), T(x)); }, bound); \
    check<T>("gamma (échelle) log_pdf " #T, d::gamma_scale_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::gamma<T>({.scale = T(b), .shape = T(a)}), T(x)); }, 1.0); \
    check<T>("beta pdf " #T, d::beta_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::beta<T>({.alpha = T(a), .beta = T(b)}), T(x)); }, bound); \
    check<T>("beta log_pdf " #T, d::beta_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::beta<T>({.alpha = T(a), .beta = T(b)}), T(x)); }, 1.0);
    GAMMA_BETA_CHECKS(double, tight)
    GAMMA_BETA_CHECKS(float, loose)
#undef GAMMA_BETA_CHECKS
#define PHASE4_CHECKS(T, bound)                                                                                              \
    check<T>("gamma cdf " #T, d::gamma_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::gamma<T>({.rate = T(b), .shape = T(a)}), T(x)); }, bound); \
    check<T>("gamma ccdf " #T, d::gamma_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::gamma<T>({.rate = T(b), .shape = T(a)}), T(x)); }, bound); \
    check<T>("gamma quantile " #T, d::gamma_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::gamma<T>({.rate = T(b), .shape = T(a)}), p); }, bound); \
    check<T>("gamma (éch.) cdf " #T, d::gamma_scale_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::gamma<T>({.scale = T(b), .shape = T(a)}), T(x)); }, bound); \
    check<T>("gamma (éch.) ccdf " #T, d::gamma_scale_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::gamma<T>({.scale = T(b), .shape = T(a)}), T(x)); }, bound); \
    check<T>("gamma (éch.) quantile " #T, d::gamma_scale_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::gamma<T>({.scale = T(b), .shape = T(a)}), p); }, bound); \
    check<T>("gamma entropy " #T, d::gamma_entropy_##T, [](double a, double b, double) { return wr::gamma<T>({.rate = T(b), .shape = T(a)}).entropy(); }, 1.0); \
    check<T>("gamma (éch.) entropy " #T, d::gamma_scale_entropy_##T, [](double a, double b, double) { return wr::gamma<T>({.scale = T(b), .shape = T(a)}).entropy(); }, 1.0); \
    check<T>("beta cdf " #T, d::beta_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::beta<T>({.alpha = T(a), .beta = T(b)}), T(x)); }, bound); \
    check<T>("beta ccdf " #T, d::beta_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::beta<T>({.alpha = T(a), .beta = T(b)}), T(x)); }, bound); \
    check<T>("beta quantile " #T, d::beta_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::beta<T>({.alpha = T(a), .beta = T(b)}), p); }, bound); \
    check<T>("beta entropy " #T, d::beta_entropy_##T, [](double a, double b, double) { return wr::beta<T>({.alpha = T(a), .beta = T(b)}).entropy(); }, 1.0);
    PHASE4_CHECKS(double, tight)
    PHASE4_CHECKS(float, loose)
#undef PHASE4_CHECKS
#define DERIVED_CHECKS(T, bound)                                                                                              \
    check<T>("chi_squared cdf " #T, d::chi_squared_cdf_##T, [](double a, double, double x) { return wr::cdf(wr::chi_squared<T>({.degrees_of_freedom = T(a)}), T(x)); }, bound); \
    check<T>("chi_squared quantile " #T, d::chi_squared_quantile_##T, [](double a, double, double p) { return wr::quantile(wr::chi_squared<T>({.degrees_of_freedom = T(a)}), p); }, bound); \
    check<T>("student_t pdf " #T, d::student_t_pdf_##T, [](double a, double, double x) { return wr::pdf(wr::student_t<T>({.degrees_of_freedom = T(a)}), T(x)); }, bound); \
    check<T>("student_t log_pdf " #T, d::student_t_log_pdf_##T, [](double a, double, double x) { return wr::log_pdf(wr::student_t<T>({.degrees_of_freedom = T(a)}), T(x)); }, 1.0); \
    check<T>("student_t cdf " #T, d::student_t_cdf_##T, [](double a, double, double x) { return wr::cdf(wr::student_t<T>({.degrees_of_freedom = T(a)}), T(x)); }, bound); \
    check<T>("student_t ccdf " #T, d::student_t_ccdf_##T, [](double a, double, double x) { return wr::ccdf(wr::student_t<T>({.degrees_of_freedom = T(a)}), T(x)); }, bound); \
    check<T>("student_t quantile " #T, d::student_t_quantile_##T, [](double a, double, double p) { return wr::quantile(wr::student_t<T>({.degrees_of_freedom = T(a)}), p); }, bound); \
    check<T>("student_t entropy " #T, d::student_t_entropy_##T, [](double a, double, double) { return wr::student_t<T>({.degrees_of_freedom = T(a)}).entropy(); }, 1.0); \
    check<T>("lognormal pdf " #T, d::lognormal_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::lognormal<T>({.log_mean = T(a), .log_stddev = T(b)}), T(x)); }, bound); \
    check<T>("lognormal log_pdf " #T, d::lognormal_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::lognormal<T>({.log_mean = T(a), .log_stddev = T(b)}), T(x)); }, 1.0); \
    check<T>("lognormal cdf " #T, d::lognormal_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::lognormal<T>({.log_mean = T(a), .log_stddev = T(b)}), T(x)); }, bound); \
    check<T>("lognormal ccdf " #T, d::lognormal_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::lognormal<T>({.log_mean = T(a), .log_stddev = T(b)}), T(x)); }, bound); \
    check<T>("lognormal quantile " #T, d::lognormal_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::lognormal<T>({.log_mean = T(a), .log_stddev = T(b)}), p); }, bound); \
    check<T>("lognormal mean " #T, d::lognormal_mean_##T, [](double a, double b, double) { return wr::lognormal<T>({.log_mean = T(a), .log_stddev = T(b)}).mean(); }, bound); \
    check<T>("lognormal variance " #T, d::lognormal_variance_##T, [](double a, double b, double) { return wr::lognormal<T>({.log_mean = T(a), .log_stddev = T(b)}).variance(); }, bound); \
    check<T>("lognormal entropy " #T, d::lognormal_entropy_##T, [](double a, double b, double) { return wr::lognormal<T>({.log_mean = T(a), .log_stddev = T(b)}).entropy(); }, 1.0); \
    check<T>("cauchy pdf " #T, d::cauchy_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::cauchy<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("cauchy log_pdf " #T, d::cauchy_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::cauchy<T>({.location = T(a), .scale = T(b)}), T(x)); }, 1.0); \
    check<T>("cauchy cdf " #T, d::cauchy_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::cauchy<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("cauchy ccdf " #T, d::cauchy_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::cauchy<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("cauchy quantile " #T, d::cauchy_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::cauchy<T>({.location = T(a), .scale = T(b)}), p); }, bound); \
    check<T>("cauchy entropy " #T, d::cauchy_entropy_##T, [](double a, double b, double) { return wr::cauchy<T>({.location = T(a), .scale = T(b)}).entropy(); }, 1.0);
    DERIVED_CHECKS(double, tight)
    DERIVED_CHECKS(float, loose)
#undef DERIVED_CHECKS
#define SECOND_CHECKS(T, bound)                                                                                               \
    check<T>("fisher_f pdf " #T, d::fisher_f_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::fisher_f<T>({.d1 = T(a), .d2 = T(b)}), T(x)); }, bound); \
    check<T>("fisher_f log_pdf " #T, d::fisher_f_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::fisher_f<T>({.d1 = T(a), .d2 = T(b)}), T(x)); }, 1.0); \
    check<T>("fisher_f cdf " #T, d::fisher_f_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::fisher_f<T>({.d1 = T(a), .d2 = T(b)}), T(x)); }, bound); \
    check<T>("fisher_f ccdf " #T, d::fisher_f_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::fisher_f<T>({.d1 = T(a), .d2 = T(b)}), T(x)); }, bound); \
    check<T>("fisher_f quantile " #T, d::fisher_f_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::fisher_f<T>({.d1 = T(a), .d2 = T(b)}), p); }, bound); \
    check<T>("fisher_f entropy " #T, d::fisher_f_entropy_##T, [](double a, double b, double) { return wr::fisher_f<T>({.d1 = T(a), .d2 = T(b)}).entropy(); }, 1.0); \
    check<T>("weibull pdf " #T, d::weibull_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::weibull<T>({.scale = T(b), .shape = T(a)}), T(x)); }, bound); \
    check<T>("weibull log_pdf " #T, d::weibull_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::weibull<T>({.scale = T(b), .shape = T(a)}), T(x)); }, 1.0); \
    check<T>("weibull cdf " #T, d::weibull_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::weibull<T>({.scale = T(b), .shape = T(a)}), T(x)); }, bound); \
    check<T>("weibull ccdf " #T, d::weibull_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::weibull<T>({.scale = T(b), .shape = T(a)}), T(x)); }, bound); \
    check<T>("weibull quantile " #T, d::weibull_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::weibull<T>({.scale = T(b), .shape = T(a)}), p); }, bound); \
    check<T>("weibull mean " #T, d::weibull_mean_##T, [](double a, double b, double) { return wr::weibull<T>({.scale = T(b), .shape = T(a)}).mean(); }, bound); \
    check<T>("weibull variance " #T, d::weibull_variance_##T, [](double a, double b, double) { return wr::weibull<T>({.scale = T(b), .shape = T(a)}).variance(); }, bound); \
    check<T>("weibull entropy " #T, d::weibull_entropy_##T, [](double a, double b, double) { return wr::weibull<T>({.scale = T(b), .shape = T(a)}).entropy(); }, 1.0); \
    check<T>("laplace pdf " #T, d::laplace_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::laplace<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("laplace log_pdf " #T, d::laplace_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::laplace<T>({.location = T(a), .scale = T(b)}), T(x)); }, 1.0); \
    check<T>("laplace cdf " #T, d::laplace_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::laplace<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("laplace ccdf " #T, d::laplace_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::laplace<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("laplace quantile " #T, d::laplace_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::laplace<T>({.location = T(a), .scale = T(b)}), p); }, bound); \
    check<T>("logistic pdf " #T, d::logistic_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::logistic<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("logistic log_pdf " #T, d::logistic_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::logistic<T>({.location = T(a), .scale = T(b)}), T(x)); }, 1.0); \
    check<T>("logistic cdf " #T, d::logistic_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::logistic<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("logistic ccdf " #T, d::logistic_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::logistic<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("logistic quantile " #T, d::logistic_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::logistic<T>({.location = T(a), .scale = T(b)}), p); }, bound); \
    check<T>("pareto pdf " #T, d::pareto_pdf_##T, [](double a, double b, double x) { return wr::pdf(wr::pareto<T>({.scale = T(a), .shape = T(b)}), T(x)); }, bound); \
    check<T>("pareto log_pdf " #T, d::pareto_log_pdf_##T, [](double a, double b, double x) { return wr::log_pdf(wr::pareto<T>({.scale = T(a), .shape = T(b)}), T(x)); }, 1.0); \
    check<T>("pareto cdf " #T, d::pareto_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::pareto<T>({.scale = T(a), .shape = T(b)}), T(x)); }, bound); \
    check<T>("pareto ccdf " #T, d::pareto_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::pareto<T>({.scale = T(a), .shape = T(b)}), T(x)); }, bound); \
    check<T>("pareto quantile " #T, d::pareto_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::pareto<T>({.scale = T(a), .shape = T(b)}), p); }, bound); \
    check<T>("pareto entropy " #T, d::pareto_entropy_##T, [](double a, double b, double) { return wr::pareto<T>({.scale = T(a), .shape = T(b)}).entropy(); }, 1.0);
#define CQ(T, bound, name, label, make)                                                                                                \
    check<T>(label " cquantile " #T, d::name##_cquantile_##T, [](double a, double b, double q) { (void)a; (void)b; return wr::cquantile(make, q); }, bound);
#define SURVIVAL_CHECKS(T, bound)                                                                                            \
    CQ(T, bound, exponential, "exponential", (exponential<T>(a, b)))                                                                  \
    CQ(T, bound, normal, "normal", (normal<T>(a, b)))                                                                                 \
    CQ(T, bound, gamma, "gamma", (wr::gamma<T>({.rate = T(b), .shape = T(a)})))                                                      \
    CQ(T, bound, gamma_scale, "gamma (éch.)", (wr::gamma<T>({.scale = T(b), .shape = T(a)})))                                         \
    CQ(T, bound, beta, "beta", (wr::beta<T>({.alpha = T(a), .beta = T(b)})))                                                         \
    CQ(T, bound, chi_squared, "chi_squared", (wr::chi_squared<T>({.degrees_of_freedom = T(a)})))                                     \
    CQ(T, bound, student_t, "student_t", (wr::student_t<T>({.degrees_of_freedom = T(a)})))                                           \
    CQ(T, bound, lognormal, "lognormal", (wr::lognormal<T>({.log_mean = T(a), .log_stddev = T(b)})))                                            \
    CQ(T, bound, cauchy, "cauchy", (wr::cauchy<T>({.location = T(a), .scale = T(b)})))                                               \
    CQ(T, bound, fisher_f, "fisher_f", (wr::fisher_f<T>({.d1 = T(a), .d2 = T(b)})))                                                  \
    CQ(T, bound, weibull, "weibull", (wr::weibull<T>({.scale = T(b), .shape = T(a)})))                                               \
    CQ(T, bound, laplace, "laplace", (wr::laplace<T>({.location = T(a), .scale = T(b)})))                                            \
    CQ(T, bound, logistic, "logistic", (wr::logistic<T>({.location = T(a), .scale = T(b)})))                                         \
    CQ(T, bound, pareto, "pareto", (wr::pareto<T>({.scale = T(a), .shape = T(b)})))
#define STANDARD_CHECKS(T, bound, name, make)                                                                                \
    check<T>(#name " pdf " #T, d::name##_pdf_##T, [](double a, double b, double x) { (void)b; return wr::pdf(make, T(x)); }, bound);              \
    check<T>(#name " log_pdf " #T, d::name##_log_pdf_##T, [](double a, double b, double x) { (void)b; return wr::log_pdf(make, T(x)); }, 1.0);    \
    check<T>(#name " cdf " #T, d::name##_cdf_##T, [](double a, double b, double x) { (void)b; return wr::cdf(make, T(x)); }, bound);              \
    check<T>(#name " ccdf " #T, d::name##_ccdf_##T, [](double a, double b, double x) { (void)b; return wr::ccdf(make, T(x)); }, bound);            \
    check<T>(#name " quantile " #T, d::name##_quantile_##T, [](double a, double b, double p) { (void)b; return wr::quantile(make, p); }, bound);   \
    check<T>(#name " cquantile " #T, d::name##_cquantile_##T, [](double a, double b, double q) { (void)b; return wr::cquantile(make, q); }, bound);
#define COMPLEMENT_CHECKS(T, bound)                                                                                          \
    STANDARD_CHECKS(T, bound, gumbel, (wr::gumbel<T>({.location = T(a), .scale = T(b)})))                                   \
    STANDARD_CHECKS(T, bound, frechet, (wr::frechet<T>({.location = T(frechet_location(a, b)), .scale = T(b), .shape = T(a)})))   \
    STANDARD_CHECKS(T, bound, rayleigh, (wr::rayleigh<T>({.scale = T(a)})))                                                 \
    STANDARD_CHECKS(T, bound, arcsine, (wr::arcsine<T>({.lower = T(a), .upper = T(b)})))                                        \
    check<T>("triangular pdf " #T, d::triangular_pdf_##T, [](double a, double b, double x) { return wr::pdf(triangle<T>(a, b), T(x)); }, bound);   \
    check<T>("triangular cdf " #T, d::triangular_cdf_##T, [](double a, double b, double x) { return wr::cdf(triangle<T>(a, b), T(x)); }, bound);   \
    check<T>("triangular ccdf " #T, d::triangular_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(triangle<T>(a, b), T(x)); }, bound); \
    check<T>("triangular quantile " #T, d::triangular_quantile_##T, [](double a, double b, double p) { return wr::quantile(triangle<T>(a, b), p); }, bound); \
    check<T>("triangular cquantile " #T, d::triangular_cquantile_##T, [](double a, double b, double q) { return wr::cquantile(triangle<T>(a, b), q); }, bound); \
    check<T>("triangular variance " #T, d::triangular_variance_##T, [](double a, double b, double) { return triangle<T>(a, b).variance(); }, bound);       \
    check<T>("triangular mean " #T, d::triangular_mean_##T, [](double a, double b, double) { return triangle<T>(a, b).mean(); }, bound);
#define INVERSE_GAUSSIAN_CHECKS(T, bound)                                                                                    \
    STANDARD_CHECKS(T, bound, inverse_gaussian, (wr::inverse_gaussian<T>({.mean = T(a), .shape = T(b)})))                   \
    check<T>("inverse_gaussian variance " #T, d::inverse_gaussian_variance_##T, [](double a, double b, double) { return wr::inverse_gaussian<T>({.mean = T(a), .shape = T(b)}).variance(); }, bound); \
    check<T>("inverse_gaussian mode " #T, d::inverse_gaussian_mode_##T, [](double a, double b, double) { return wr::inverse_gaussian<T>({.mean = T(a), .shape = T(b)}).mode(); }, bound);           \
    check<T>("inverse_gaussian skewness " #T, d::inverse_gaussian_skewness_##T, [](double a, double b, double) { return wr::inverse_gaussian<T>({.mean = T(a), .shape = T(b)}).skewness(); }, bound); \
    check<T>("inverse_gaussian excess_kurtosis " #T, d::inverse_gaussian_excess_kurtosis_##T, [](double a, double b, double) { return wr::inverse_gaussian<T>({.mean = T(a), .shape = T(b)}).excess_kurtosis(); }, bound);
#define VON_MISES_CHECKS(T, bound)                                                                                            \
    check<T>("von_mises cdf " #T, d::von_mises_cdf_##T, [](double a, double b, double x) { return wr::cdf(wr::von_mises<T>({.location = T(a), .concentration = T(b)}), T(x)); }, bound);      \
    check<T>("von_mises ccdf " #T, d::von_mises_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(wr::von_mises<T>({.location = T(a), .concentration = T(b)}), T(x)); }, bound);    \
    check<T>("von_mises quantile " #T, d::von_mises_quantile_##T, [](double a, double b, double p) { return wr::quantile(wr::von_mises<T>({.location = T(a), .concentration = T(b)}), p); }, bound); \
    check<T>("von_mises cquantile " #T, d::von_mises_cquantile_##T, [](double a, double b, double q) { return wr::cquantile(wr::von_mises<T>({.location = T(a), .concentration = T(b)}), q); }, bound);
#define SKEW_NORMAL_CHECKS(T, bound)                                                                                              check<T>("skew_normal pdf " #T, d::skew_normal_pdf_##T, [](double a, double b, double x) { return wr::pdf(skew_of<T>(a, b), T(x)); }, bound);             check<T>("skew_normal cdf " #T, d::skew_normal_cdf_##T, [](double a, double b, double x) { return wr::cdf(skew_of<T>(a, b), T(x)); }, bound);             check<T>("skew_normal ccdf " #T, d::skew_normal_ccdf_##T, [](double a, double b, double x) { return wr::ccdf(skew_of<T>(a, b), T(x)); }, bound);           check<T>("skew_normal quantile " #T, d::skew_normal_quantile_##T, [](double a, double b, double p) { return wr::quantile(skew_of<T>(a, b), p); }, bound);     check<T>("skew_normal cquantile " #T, d::skew_normal_cquantile_##T, [](double a, double b, double q) { return wr::cquantile(skew_of<T>(a, b), q); }, bound);     check<T>("skew_normal mean " #T, d::skew_normal_mean_##T, [](double a, double b, double) { return skew_of<T>(a, b).mean(); }, bound);                  check<T>("skew_normal variance " #T, d::skew_normal_variance_##T, [](double a, double b, double) { return skew_of<T>(a, b).variance(); }, bound);      check<T>("skew_normal skewness " #T, d::skew_normal_skewness_##T, [](double a, double b, double) { return skew_of<T>(a, b).skewness(); }, bound);
    SKEW_NORMAL_CHECKS(double, tight)
    SKEW_NORMAL_CHECKS(float, loose)
#undef SKEW_NORMAL_CHECKS
    {
        auto const St = [](double a, double b) { return wr::stable({.alpha = a, .beta = b}); };
        check<double>("stable pdf", d::stable_pdf_double, [&](double a, double b, double x) { return wr::pdf(St(a, b), x); }, tight);
        check<double>("stable cdf", d::stable_cdf_double, [&](double a, double b, double x) { return wr::cdf(St(a, b), x); }, tight);
        check<double>("stable ccdf", d::stable_ccdf_double, [&](double a, double b, double x) { return wr::ccdf(St(a, b), x); }, tight);
    }
    VON_MISES_CHECKS(double, tight)
    VON_MISES_CHECKS(float, loose)
#undef VON_MISES_CHECKS
    INVERSE_GAUSSIAN_CHECKS(double, tight)
    INVERSE_GAUSSIAN_CHECKS(float, loose)
#undef INVERSE_GAUSSIAN_CHECKS
    COMPLEMENT_CHECKS(double, tight)
    COMPLEMENT_CHECKS(float, loose)
    check<double>("gumbel mean", d::gumbel_mean_double, [](double a, double b, double) { return wr::gumbel({.location = a, .scale = b}).mean(); }, tight);
    check<double>("gumbel variance", d::gumbel_variance_double, [](double a, double b, double) { return wr::gumbel({.location = a, .scale = b}).variance(); }, tight);
    check<double>("gumbel median", d::gumbel_median_double, [](double a, double b, double) { return wr::gumbel({.location = a, .scale = b}).median(); }, tight);
    check<double>("gumbel entropy", d::gumbel_entropy_double, [](double a, double b, double) { return wr::gumbel({.location = a, .scale = b}).entropy(); }, 1.0);
    check<double>("frechet mean", d::frechet_mean_double, [](double a, double b, double) { return wr::frechet({.location = frechet_location(a, b), .scale = b, .shape = a}).mean(); }, tight);
    check<double>("frechet variance", d::frechet_variance_double, [](double a, double b, double) { return wr::frechet({.location = frechet_location(a, b), .scale = b, .shape = a}).variance(); }, tight);
    check<double>("frechet median", d::frechet_median_double, [](double a, double b, double) { return wr::frechet({.location = frechet_location(a, b), .scale = b, .shape = a}).median(); }, tight);
    check<double>("frechet mode", d::frechet_mode_double, [](double a, double b, double) { return wr::frechet({.location = frechet_location(a, b), .scale = b, .shape = a}).mode(); }, tight);
    check<double>("frechet entropy", d::frechet_entropy_double, [](double a, double b, double) { return wr::frechet({.location = frechet_location(a, b), .scale = b, .shape = a}).entropy(); }, 1.0);
    check<double>("rayleigh mean", d::rayleigh_mean_double, [](double a, double, double) { return wr::rayleigh({.scale = a}).mean(); }, tight);
    check<double>("rayleigh variance", d::rayleigh_variance_double, [](double a, double, double) { return wr::rayleigh({.scale = a}).variance(); }, tight);
    check<double>("rayleigh median", d::rayleigh_median_double, [](double a, double, double) { return wr::rayleigh({.scale = a}).median(); }, tight);
    check<double>("rayleigh entropy", d::rayleigh_entropy_double, [](double a, double, double) { return wr::rayleigh({.scale = a}).entropy(); }, 1.0);
    check<double>("arcsine variance", d::arcsine_variance_double, [](double a, double b, double) { return wr::arcsine({.lower = a, .upper = b}).variance(); }, tight);
    check<double>("arcsine entropy", d::arcsine_entropy_double, [](double a, double b, double) { return wr::arcsine({.lower = a, .upper = b}).entropy(); }, 1.0);
#undef COMPLEMENT_CHECKS
#undef STANDARD_CHECKS
    SURVIVAL_CHECKS(double, tight)
    {
        // von Mises et von Mises–Fisher (μ = e_p, x = (√(1 − t²), 0, …, 0, t)).
        auto const V = [](double m, double k) { return wr::von_mises({.location = m, .concentration = k}); };
        check<double>("von_mises pdf", d::von_mises_pdf_double, [&](double a, double b, double x) { return wr::pdf(V(a, b), x); }, tight);
        check<double>("von_mises log_pdf", d::von_mises_log_pdf_double, [&](double a, double b, double x) { return wr::log_pdf(V(a, b), x); }, 1.0);
        check<double>("von_mises resultant", d::von_mises_resultant_double, [&](double a, double b, double) { return V(a, b).mean_resultant_length(); }, tight);
        check<double>("von_mises circ. var.", d::von_mises_circular_variance_double, [&](double a, double b, double) { return V(a, b).circular_variance(); }, tight);
        check<double>("von_mises entropy", d::von_mises_entropy_double, [&](double a, double b, double) { return V(a, b).entropy(); }, 1.0);
        // (dimension p >= 3 : vecteurs construits par push_back, sans front()/back() que GCC -O2
        // croit pouvoir viser un vecteur vide, -Warray-bounds)
        auto const F = [](double k, double p) {
            std::vector<double> mu(static_cast<std::size_t>(p) - 1, 0.0);
            mu.push_back(1.0);
            return wr::von_mises_fisher<double>({.mean_direction = mu, .concentration = k});
        };
        auto const point = [](double p, double t) {
            std::vector<double> x{wsl::special_functions::detail::pinned_portable::sqrt(wsl::floating_point_arithmetic::detail::fmul(1 - t, 1 + t))};
            x.resize(static_cast<std::size_t>(p) - 1, 0.0);
            x.push_back(t);
            return x;
        };
        check<double>("vmf pdf", d::von_mises_fisher_pdf_double, [&](double a, double b, double t) { return wr::pdf(F(a, b), point(b, t)); }, tight);
        check<double>("vmf log_pdf", d::von_mises_fisher_log_pdf_double, [&](double a, double b, double t) { return wr::log_pdf(F(a, b), point(b, t)); }, 1.0);
        check<double>("vmf resultant", d::von_mises_fisher_resultant_double, [&](double a, double b, double) { return F(a, b).mean_resultant_length(); }, tight);
        check<double>("vmf entropy", d::von_mises_fisher_entropy_double, [&](double a, double b, double) { return F(a, b).entropy(); }, 1.0);
    }
    {
        // Lois construites : erreur de composition (design/random.md), bornes mesurées.
        auto const Z = [](double nu) { return wr::affine(wr::student_t({.degrees_of_freedom = nu}), {.scale = -2, .shift = 5}); };
        auto const N = [](double a, double b) { return wr::truncate(wr::normal({.mean = 0.0, .stddev = 1.0}), {.lower = a, .upper = b}); };
        wr::mixture const M(wr::weighted(0.3, wr::normal({.mean = -2.0, .stddev = 1.0})), wr::weighted(0.7, wr::gamma({.rate = 1.0, .shape = 3.0})));
        // Bornes : mesures arrondies au-dessus ; le quantile tronqué près d'un zéro du résultat.
        check<double>("affine student pdf", d::affine_student_pdf_double, [&](double a, double, double y) { return wr::pdf(Z(a), y); }, 0.51);
        check<double>("affine student cdf", d::affine_student_cdf_double, [&](double a, double, double y) { return wr::cdf(Z(a), y); }, 0.51);
        check<double>("affine student ccdf", d::affine_student_ccdf_double, [&](double a, double, double y) { return wr::ccdf(Z(a), y); }, 0.51);
        check<double>("affine student quantile", d::affine_student_quantile_double, [&](double a, double, double p) { return wr::quantile(Z(a), p); }, 4.0);
        check<double>("truncated normal cdf", d::truncated_normal_cdf_double, [&](double a, double b, double x) { return wr::cdf(N(a, b), x); }, 3.0);
        check<double>("truncated normal ccdf", d::truncated_normal_ccdf_double, [&](double a, double b, double x) { return wr::ccdf(N(a, b), x); }, 6.0);
        check<double>("truncated normal quantile", d::truncated_normal_quantile_double, [&](double a, double b, double p) { return wr::quantile(N(a, b), p); }, 64.0);
        auto const PB = [](double mu, double p) { return wr::convolution(wr::poisson({.mean = mu}), wr::binomial({.trials = 20, .p = p})); };
        check<double>("sum poisson+binomial pmf", d::sum_poisson_binomial_pmf_double, [&](double a, double b, double k) { return wr::pmf(PB(a, b), static_cast<std::int64_t>(k)); }, 1.5);
        check<double>("sum poisson+binomial cdf", d::sum_poisson_binomial_cdf_double, [&](double a, double b, double k) { return wr::cdf(PB(a, b), static_cast<std::int64_t>(k)); }, 1.0);
        check<double>("sum poisson+binomial ccdf", d::sum_poisson_binomial_ccdf_double, [&](double a, double b, double k) { return wr::ccdf(PB(a, b), static_cast<std::int64_t>(k)); }, 1.0);
        check<double>("mixture pdf", d::mixture_pdf_double, [&](double, double, double x) { return wr::pdf(M, x); }, 1.5);
        check<double>("mixture cdf", d::mixture_cdf_double, [&](double, double, double x) { return wr::cdf(M, x); }, 1.5);
        check<double>("mixture quantile", d::mixture_quantile_double, [&](double, double, double p) { return wr::quantile(M, p); }, 2.0);
    }
    SURVIVAL_CHECKS(float, loose)
#undef SURVIVAL_CHECKS
#undef CQ
    SECOND_CHECKS(double, tight)
    SECOND_CHECKS(float, loose)
#undef SECOND_CHECKS
    LAW_CHECKS(double, tight)
    LAW_CHECKS(float, loose)
#undef LAW_CHECKS

    // Lois discrètes (masses en double).
    using K = long long;
    auto const G = [](double a) { return wr::geometric<K>({.p = a}); };
    auto const P = [](double a) { return wr::poisson<K>({.mean = a}); };
    auto const B = [](double a, double b) { return wr::binomial<K>({.trials = K(a), .p = b}); };
    check<double>("geometric pmf", d::geometric_pmf_double, [&](double a, double, double x) { return wr::pmf(G(a), K(x)); }, tight);
    check<double>("geometric log_pmf", d::geometric_log_pmf_double, [&](double a, double, double x) { return wr::log_pmf(G(a), K(x)); }, 1.0);
    check<double>("geometric cdf", d::geometric_cdf_double, [&](double a, double, double x) { return wr::cdf(G(a), K(x)); }, tight);
    check<double>("geometric ccdf", d::geometric_ccdf_double, [&](double a, double, double x) { return wr::ccdf(G(a), K(x)); }, tight);
    check<double>("poisson pmf", d::poisson_pmf_double, [&](double a, double, double x) { return wr::pmf(P(a), K(x)); }, tight);
    check<double>("poisson log_pmf", d::poisson_log_pmf_double, [&](double a, double, double x) { return wr::log_pmf(P(a), K(x)); }, 1.0);
    check<double>("binomial pmf", d::binomial_pmf_double, [&](double a, double b, double x) { return wr::pmf(B(a, b), K(x)); }, tight);
    check<double>("binomial log_pmf", d::binomial_log_pmf_double, [&](double a, double b, double x) { return wr::log_pmf(B(a, b), K(x)); }, 1.0);

    check<double>("poisson cdf", d::poisson_cdf_double, [](double a, double, double x) { return wr::cdf(wr::poisson<K>({.mean = a}), K(x)); }, tight);
    check<double>("poisson ccdf", d::poisson_ccdf_double, [](double a, double, double x) { return wr::ccdf(wr::poisson<K>({.mean = a}), K(x)); }, tight);
    auto const NB = [](double a, double b) { return wr::negative_binomial<K>({.successes = a, .p = b}); };
    check<double>("negative_binomial pmf", d::negative_binomial_pmf_double, [&](double a, double b, double x) { return wr::pmf(NB(a, b), K(x)); }, tight);
    check<double>("negative_binomial log_pmf", d::negative_binomial_log_pmf_double, [&](double a, double b, double x) { return wr::log_pmf(NB(a, b), K(x)); }, 1.0);
    check<double>("negative_binomial cdf", d::negative_binomial_cdf_double, [&](double a, double b, double x) { return wr::cdf(NB(a, b), K(x)); }, tight);
    check<double>("negative_binomial ccdf", d::negative_binomial_ccdf_double, [&](double a, double b, double x) { return wr::ccdf(NB(a, b), K(x)); }, tight);
    check<double>("beta_binomial pmf", d::beta_binomial_pmf_double, [](double a, double b, double x) { return wr::pmf(beta_binomial_of<K>(a, b), K(x)); }, tight);
    check<double>("beta_binomial log_pmf", d::beta_binomial_log_pmf_double, [](double a, double b, double x) { return wr::log_pmf(beta_binomial_of<K>(a, b), K(x)); }, 1.0);
    check<double>("beta_binomial cdf", d::beta_binomial_cdf_double, [](double a, double b, double x) { return wr::cdf(beta_binomial_of<K>(a, b), K(x)); }, tight);
    check<double>("beta_binomial ccdf", d::beta_binomial_ccdf_double, [](double a, double b, double x) { return wr::ccdf(beta_binomial_of<K>(a, b), K(x)); }, tight);
    check<double>("beta_binomial variance", d::beta_binomial_variance_double, [](double a, double b, double) { return beta_binomial_of<K>(a, b).variance(); }, tight);
    check<double>("logarithmic pmf", d::logarithmic_pmf_double, [](double a, double, double x) { return wr::pmf(wr::logarithmic<K>({.p = a}), K(x)); }, tight);
    check<double>("logarithmic cdf", d::logarithmic_cdf_double, [](double a, double, double x) { return wr::cdf(wr::logarithmic<K>({.p = a}), K(x)); }, tight);
    check<double>("logarithmic ccdf", d::logarithmic_ccdf_double, [](double a, double, double x) { return wr::ccdf(wr::logarithmic<K>({.p = a}), K(x)); }, tight);
    check<double>("logarithmic mean", d::logarithmic_mean_double, [](double a, double, double) { return wr::logarithmic<K>({.p = a}).mean(); }, tight);
    check<double>("logarithmic variance", d::logarithmic_variance_double, [](double a, double, double) { return wr::logarithmic<K>({.p = a}).variance(); }, tight);
    auto const DG = [](double a, double b) { return wr::discrete_gaussian<K>({.location = a, .scale = b}); };
    check<double>("discrete_gaussian pmf", d::discrete_gaussian_pmf_double, [&](double a, double b, double x) { return wr::pmf(DG(a, b), K(x)); }, tight);
    check<double>("discrete_gaussian cdf", d::discrete_gaussian_cdf_double, [&](double a, double b, double x) { return wr::cdf(DG(a, b), K(x)); }, tight);
    check<double>("discrete_gaussian ccdf", d::discrete_gaussian_ccdf_double, [&](double a, double b, double x) { return wr::ccdf(DG(a, b), K(x)); }, tight);
    check<double>("discrete_gaussian mean", d::discrete_gaussian_mean_double, [&](double a, double b, double) { return DG(a, b).mean(); }, tight);
    check<double>("discrete_gaussian variance", d::discrete_gaussian_variance_double, [&](double a, double b, double) { return DG(a, b).variance(); }, tight);
    auto const RI = [](double a, double b) { return wr::rice<double>({.scale = a, .noncentrality = b}); };
    check<double>("rice pdf", d::rice_pdf_double, [&](double a, double b, double x) { return RI(a, b).pdf(x); }, tight);
    check<double>("rice log_pdf", d::rice_log_pdf_double, [&](double a, double b, double x) { return RI(a, b).log_pdf(x); }, 1.0);
    check<double>("rice cdf", d::rice_cdf_double, [&](double a, double b, double x) { return RI(a, b).cdf(x); }, tight);
    check<double>("rice ccdf", d::rice_ccdf_double, [&](double a, double b, double x) { return RI(a, b).ccdf(x); }, tight);
    check<double>("rice quantile", d::rice_quantile_double, [&](double a, double b, double p) { return RI(a, b).quantile(p); }, tight);
    check<double>("rice cquantile", d::rice_cquantile_double, [&](double a, double b, double q) { return RI(a, b).cquantile(q); }, tight);
    check<double>("rice mean", d::rice_mean_double, [&](double a, double b, double) { return RI(a, b).mean(); }, tight);
    check<double>("rice variance", d::rice_variance_double, [&](double a, double b, double) { return RI(a, b).variance(); }, tight);
    auto const ZI = [](double a, double b) { return wr::zipf<K>({.size = K(a), .exponent = b}); };
    check<double>("zipf pmf", d::zipf_pmf_double, [&](double a, double b, double x) { return wr::pmf(ZI(a, b), K(x)); }, tight);
    check<double>("zipf log_pmf", d::zipf_log_pmf_double, [&](double a, double b, double x) { return wr::log_pmf(ZI(a, b), K(x)); }, 1.0);
    check<double>("zipf cdf", d::zipf_cdf_double, [&](double a, double b, double x) { return wr::cdf(ZI(a, b), K(x)); }, tight);
    check<double>("zipf ccdf", d::zipf_ccdf_double, [&](double a, double b, double x) { return wr::ccdf(ZI(a, b), K(x)); }, tight);
    check<double>("zipf mean", d::zipf_mean_double, [&](double a, double b, double) { return ZI(a, b).mean(); }, tight);
    check<double>("zipf variance", d::zipf_variance_double, [&](double a, double b, double) { return ZI(a, b).variance(); }, tight);
    auto const SK = [](double a, double b) { return wr::skellam<K>({.mean1 = a, .mean2 = b}); };
    check<double>("skellam pmf", d::skellam_pmf_double, [&](double a, double b, double x) { return wr::pmf(SK(a, b), K(x)); }, tight);
    check<double>("skellam log_pmf", d::skellam_log_pmf_double, [&](double a, double b, double x) { return wr::log_pmf(SK(a, b), K(x)); }, 1.0);
    check<double>("skellam cdf", d::skellam_cdf_double, [&](double a, double b, double x) { return wr::cdf(SK(a, b), K(x)); }, tight);
    check<double>("skellam ccdf", d::skellam_ccdf_double, [&](double a, double b, double x) { return wr::ccdf(SK(a, b), K(x)); }, tight);
    check<double>("skellam skewness", d::skellam_skewness_double, [&](double a, double b, double) { return SK(a, b).skewness(); }, tight);
    check<double>("skellam excess_kurtosis", d::skellam_excess_kurtosis_double, [&](double a, double b, double) { return SK(a, b).excess_kurtosis(); }, tight);
    check<double>("hypergeometric pmf", d::hypergeometric_pmf_double, [](double a, double b, double x) { return wr::pmf(hypergeometric_of<K>(a, b), K(x)); }, tight);
    check<double>("hypergeometric log_pmf", d::hypergeometric_log_pmf_double, [](double a, double b, double x) { return wr::log_pmf(hypergeometric_of<K>(a, b), K(x)); }, 1.0);
    check<double>("hypergeometric cdf", d::hypergeometric_cdf_double, [](double a, double b, double x) { return wr::cdf(hypergeometric_of<K>(a, b), K(x)); }, tight);
    check<double>("hypergeometric ccdf", d::hypergeometric_ccdf_double, [](double a, double b, double x) { return wr::ccdf(hypergeometric_of<K>(a, b), K(x)); }, tight);
    check<double>("hypergeometric mean", d::hypergeometric_mean_double, [](double a, double b, double) { return hypergeometric_of<K>(a, b).mean(); }, tight);
    check<double>("hypergeometric variance", d::hypergeometric_variance_double, [](double a, double b, double) { return hypergeometric_of<K>(a, b).variance(); }, tight);
    check<double>("binomial cdf", d::binomial_cdf_double, [](double a, double b, double x) { return wr::cdf(wr::binomial<K>({.trials = K(a), .p = b}), K(x)); }, tight);
    check<double>("binomial ccdf", d::binomial_ccdf_double, [](double a, double b, double x) { return wr::ccdf(wr::binomial<K>({.trials = K(a), .p = b}), K(x)); }, tight);
    {
        std::size_t wrong = 0;
        for (auto const& q : d::poisson_quantile)
            wrong += wr::quantile(wr::poisson<K>({.mean = q.a}), q.u) != q.k;
        std::printf("%-28s %5zu points : %zu faux\n", "poisson quantile", std::size(d::poisson_quantile), wrong);
        failed |= wrong != 0;
        wrong = 0;
        for (auto const& q : d::binomial_quantile)
            wrong += wr::quantile(wr::binomial<K>({.trials = K(q.a), .p = q.b}), q.u) != q.k;
        std::printf("%-28s %5zu points : %zu faux\n", "binomial quantile", std::size(d::binomial_quantile), wrong);
        failed |= wrong != 0;
        wrong = 0;
        for (auto const& q : d::negative_binomial_quantile)
            wrong += wr::quantile(wr::negative_binomial<K>({.successes = q.a, .p = q.b}), q.u) != q.k;
        std::printf("%-28s %5zu points : %zu faux\n", "negative_binomial quantile", std::size(d::negative_binomial_quantile), wrong);
        failed |= wrong != 0;
        wrong = 0;
        for (auto const& q : d::hypergeometric_quantile)
            wrong += wr::quantile(hypergeometric_of<K>(q.a, q.b), q.u) != q.k;
        std::printf("%-28s %5zu points : %zu faux\n", "hypergeometric quantile", std::size(d::hypergeometric_quantile), wrong);
        failed |= wrong != 0;
        auto const quantiles = [&](char const* name, auto const& table, auto make) {
            std::size_t w = 0;
            for (auto const& q : table)
                w += wr::quantile(make(q.a, q.b), q.u) != q.k;
            std::printf("%-28s %5zu points : %zu faux\n", name, std::size(table), w);
            failed |= w != 0;
        };
        quantiles("beta_binomial quantile", d::beta_binomial_quantile, [](double a, double b) { return beta_binomial_of<K>(a, b); });
        quantiles("logarithmic quantile", d::logarithmic_quantile, [](double a, double) { return wr::logarithmic<K>({.p = a}); });
        quantiles("zipf quantile", d::zipf_quantile, [](double a, double b) { return wr::zipf<K>({.size = K(a), .exponent = b}); });
        quantiles("skellam quantile", d::skellam_quantile, [](double a, double b) { return wr::skellam<K>({.mean1 = a, .mean2 = b}); });
        quantiles("discrete_gaussian quantile", d::discrete_gaussian_quantile, [](double a, double b) { return wr::discrete_gaussian<K>({.location = a, .scale = b}); });
    }

    // Quantile de la géométrique : exact.
    std::size_t wrong = 0;
    for (auto const& q : d::geometric_quantile)
        wrong += wr::quantile(G(q.p), q.u) != q.k;
    std::printf("%-28s %5zu points : %zu faux\n", "geometric quantile", std::size(d::geometric_quantile), wrong);
    failed |= wrong != 0;
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
