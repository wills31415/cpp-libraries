// Précision des lois de wsl::random : erreur en ulp contre les valeurs de
// référence de mpmath (tools/random/law_references.py, engendrées par CMake).

#include <wsl/random.hpp>

#include "data/random_law_reference.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <iterator>
#include <span>

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
    wsl::exponential<T> exponential(double rate, double) { return wsl::exponential<T>({.rate = T(rate)}); }
    template <class T>
    wsl::normal<T> normal(double mean, double sd) { return wsl::normal<T>({.mean = T(mean), .sd = T(sd)}); }
}

int main()
{
    namespace d = law_reference_data;
    // En double : un seul arrondi final, sur un résultat intermédiaire précis à
    // 2^-60 près ou mieux ; log_pdf perd en relatif quand ses termes se compensent.
    // En float : le résultat double arrondi en float.
    constexpr double tight = 0.51, loose = 0.5 + 0x1p-20;
#define LAW_CHECKS(T, bound)                                                                                                  \
    check<T>("exponential pdf " #T, d::exponential_pdf_##T, [](double a, double b, double x) { return pdf(exponential<T>(a, b), T(x)); }, bound); \
    check<T>("exponential log_pdf " #T, d::exponential_log_pdf_##T, [](double a, double b, double x) { return log_pdf(exponential<T>(a, b), T(x)); }, 1.0); \
    check<T>("exponential cdf " #T, d::exponential_cdf_##T, [](double a, double b, double x) { return cdf(exponential<T>(a, b), T(x)); }, bound); \
    check<T>("exponential ccdf " #T, d::exponential_ccdf_##T, [](double a, double b, double x) { return ccdf(exponential<T>(a, b), T(x)); }, bound); \
    check<T>("exponential quantile " #T, d::exponential_quantile_##T, [](double a, double b, double p) { return quantile(exponential<T>(a, b), p); }, bound); \
    check<T>("normal pdf " #T, d::normal_pdf_##T, [](double a, double b, double x) { return pdf(normal<T>(a, b), T(x)); }, bound); \
    check<T>("normal log_pdf " #T, d::normal_log_pdf_##T, [](double a, double b, double x) { return log_pdf(normal<T>(a, b), T(x)); }, 1.0); \
    check<T>("normal cdf " #T, d::normal_cdf_##T, [](double a, double b, double x) { return cdf(normal<T>(a, b), T(x)); }, bound); \
    check<T>("normal ccdf " #T, d::normal_ccdf_##T, [](double a, double b, double x) { return ccdf(normal<T>(a, b), T(x)); }, bound); \
    check<T>("normal quantile " #T, d::normal_quantile_##T, [](double a, double b, double p) { return quantile(normal<T>(a, b), p); }, bound);
#define GAMMA_BETA_CHECKS(T, bound)                                                                                           \
    check<T>("gamma pdf " #T, d::gamma_pdf_##T, [](double a, double b, double x) { return pdf(wsl::gamma<T>({.shape = T(a), .rate = T(b)}), T(x)); }, bound); \
    check<T>("gamma log_pdf " #T, d::gamma_log_pdf_##T, [](double a, double b, double x) { return log_pdf(wsl::gamma<T>({.shape = T(a), .rate = T(b)}), T(x)); }, 1.0); \
    check<T>("gamma (échelle) pdf " #T, d::gamma_scale_pdf_##T, [](double a, double b, double x) { return pdf(wsl::gamma<T>({.shape = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("gamma (échelle) log_pdf " #T, d::gamma_scale_log_pdf_##T, [](double a, double b, double x) { return log_pdf(wsl::gamma<T>({.shape = T(a), .scale = T(b)}), T(x)); }, 1.0); \
    check<T>("beta pdf " #T, d::beta_pdf_##T, [](double a, double b, double x) { return pdf(wsl::beta<T>({.alpha = T(a), .beta = T(b)}), T(x)); }, bound); \
    check<T>("beta log_pdf " #T, d::beta_log_pdf_##T, [](double a, double b, double x) { return log_pdf(wsl::beta<T>({.alpha = T(a), .beta = T(b)}), T(x)); }, 1.0);
    GAMMA_BETA_CHECKS(double, tight)
    GAMMA_BETA_CHECKS(float, loose)
#undef GAMMA_BETA_CHECKS
#define PHASE4_CHECKS(T, bound)                                                                                              \
    check<T>("gamma cdf " #T, d::gamma_cdf_##T, [](double a, double b, double x) { return cdf(wsl::gamma<T>({.shape = T(a), .rate = T(b)}), T(x)); }, bound); \
    check<T>("gamma ccdf " #T, d::gamma_ccdf_##T, [](double a, double b, double x) { return ccdf(wsl::gamma<T>({.shape = T(a), .rate = T(b)}), T(x)); }, bound); \
    check<T>("gamma quantile " #T, d::gamma_quantile_##T, [](double a, double b, double p) { return quantile(wsl::gamma<T>({.shape = T(a), .rate = T(b)}), p); }, bound); \
    check<T>("gamma (éch.) cdf " #T, d::gamma_scale_cdf_##T, [](double a, double b, double x) { return cdf(wsl::gamma<T>({.shape = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("gamma (éch.) ccdf " #T, d::gamma_scale_ccdf_##T, [](double a, double b, double x) { return ccdf(wsl::gamma<T>({.shape = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("gamma (éch.) quantile " #T, d::gamma_scale_quantile_##T, [](double a, double b, double p) { return quantile(wsl::gamma<T>({.shape = T(a), .scale = T(b)}), p); }, bound); \
    check<T>("gamma entropy " #T, d::gamma_entropy_##T, [](double a, double b, double) { return wsl::gamma<T>({.shape = T(a), .rate = T(b)}).entropy(); }, 1.0); \
    check<T>("gamma (éch.) entropy " #T, d::gamma_scale_entropy_##T, [](double a, double b, double) { return wsl::gamma<T>({.shape = T(a), .scale = T(b)}).entropy(); }, 1.0); \
    check<T>("beta cdf " #T, d::beta_cdf_##T, [](double a, double b, double x) { return cdf(wsl::beta<T>({.alpha = T(a), .beta = T(b)}), T(x)); }, bound); \
    check<T>("beta ccdf " #T, d::beta_ccdf_##T, [](double a, double b, double x) { return ccdf(wsl::beta<T>({.alpha = T(a), .beta = T(b)}), T(x)); }, bound); \
    check<T>("beta quantile " #T, d::beta_quantile_##T, [](double a, double b, double p) { return quantile(wsl::beta<T>({.alpha = T(a), .beta = T(b)}), p); }, bound); \
    check<T>("beta entropy " #T, d::beta_entropy_##T, [](double a, double b, double) { return wsl::beta<T>({.alpha = T(a), .beta = T(b)}).entropy(); }, 1.0);
    PHASE4_CHECKS(double, tight)
    PHASE4_CHECKS(float, loose)
#undef PHASE4_CHECKS
#define DERIVED_CHECKS(T, bound)                                                                                              \
    check<T>("chi_squared cdf " #T, d::chi_squared_cdf_##T, [](double a, double, double x) { return cdf(wsl::chi_squared<T>({.degrees_of_freedom = T(a)}), T(x)); }, bound); \
    check<T>("chi_squared quantile " #T, d::chi_squared_quantile_##T, [](double a, double, double p) { return quantile(wsl::chi_squared<T>({.degrees_of_freedom = T(a)}), p); }, bound); \
    check<T>("student_t pdf " #T, d::student_t_pdf_##T, [](double a, double, double x) { return pdf(wsl::student_t<T>({.degrees_of_freedom = T(a)}), T(x)); }, bound); \
    check<T>("student_t log_pdf " #T, d::student_t_log_pdf_##T, [](double a, double, double x) { return log_pdf(wsl::student_t<T>({.degrees_of_freedom = T(a)}), T(x)); }, 1.0); \
    check<T>("student_t cdf " #T, d::student_t_cdf_##T, [](double a, double, double x) { return cdf(wsl::student_t<T>({.degrees_of_freedom = T(a)}), T(x)); }, bound); \
    check<T>("student_t ccdf " #T, d::student_t_ccdf_##T, [](double a, double, double x) { return ccdf(wsl::student_t<T>({.degrees_of_freedom = T(a)}), T(x)); }, bound); \
    check<T>("student_t quantile " #T, d::student_t_quantile_##T, [](double a, double, double p) { return quantile(wsl::student_t<T>({.degrees_of_freedom = T(a)}), p); }, bound); \
    check<T>("student_t entropy " #T, d::student_t_entropy_##T, [](double a, double, double) { return wsl::student_t<T>({.degrees_of_freedom = T(a)}).entropy(); }, 1.0); \
    check<T>("lognormal pdf " #T, d::lognormal_pdf_##T, [](double a, double b, double x) { return pdf(wsl::lognormal<T>({.mu = T(a), .sigma = T(b)}), T(x)); }, bound); \
    check<T>("lognormal log_pdf " #T, d::lognormal_log_pdf_##T, [](double a, double b, double x) { return log_pdf(wsl::lognormal<T>({.mu = T(a), .sigma = T(b)}), T(x)); }, 1.0); \
    check<T>("lognormal cdf " #T, d::lognormal_cdf_##T, [](double a, double b, double x) { return cdf(wsl::lognormal<T>({.mu = T(a), .sigma = T(b)}), T(x)); }, bound); \
    check<T>("lognormal ccdf " #T, d::lognormal_ccdf_##T, [](double a, double b, double x) { return ccdf(wsl::lognormal<T>({.mu = T(a), .sigma = T(b)}), T(x)); }, bound); \
    check<T>("lognormal quantile " #T, d::lognormal_quantile_##T, [](double a, double b, double p) { return quantile(wsl::lognormal<T>({.mu = T(a), .sigma = T(b)}), p); }, bound); \
    check<T>("lognormal mean " #T, d::lognormal_mean_##T, [](double a, double b, double) { return wsl::lognormal<T>({.mu = T(a), .sigma = T(b)}).mean(); }, bound); \
    check<T>("lognormal variance " #T, d::lognormal_variance_##T, [](double a, double b, double) { return wsl::lognormal<T>({.mu = T(a), .sigma = T(b)}).variance(); }, bound); \
    check<T>("lognormal entropy " #T, d::lognormal_entropy_##T, [](double a, double b, double) { return wsl::lognormal<T>({.mu = T(a), .sigma = T(b)}).entropy(); }, 1.0); \
    check<T>("cauchy pdf " #T, d::cauchy_pdf_##T, [](double a, double b, double x) { return pdf(wsl::cauchy<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("cauchy log_pdf " #T, d::cauchy_log_pdf_##T, [](double a, double b, double x) { return log_pdf(wsl::cauchy<T>({.location = T(a), .scale = T(b)}), T(x)); }, 1.0); \
    check<T>("cauchy cdf " #T, d::cauchy_cdf_##T, [](double a, double b, double x) { return cdf(wsl::cauchy<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("cauchy ccdf " #T, d::cauchy_ccdf_##T, [](double a, double b, double x) { return ccdf(wsl::cauchy<T>({.location = T(a), .scale = T(b)}), T(x)); }, bound); \
    check<T>("cauchy quantile " #T, d::cauchy_quantile_##T, [](double a, double b, double p) { return quantile(wsl::cauchy<T>({.location = T(a), .scale = T(b)}), p); }, bound); \
    check<T>("cauchy entropy " #T, d::cauchy_entropy_##T, [](double a, double b, double) { return wsl::cauchy<T>({.location = T(a), .scale = T(b)}).entropy(); }, 1.0);
    DERIVED_CHECKS(double, tight)
    DERIVED_CHECKS(float, loose)
#undef DERIVED_CHECKS
    LAW_CHECKS(double, tight)
    LAW_CHECKS(float, loose)
#undef LAW_CHECKS

    // Lois discrètes (masses en double).
    using K = long long;
    auto const G = [](double a) { return wsl::geometric<K>({.p = a}); };
    auto const P = [](double a) { return wsl::poisson<K>({.mean = a}); };
    auto const B = [](double a, double b) { return wsl::binomial<K>({.trials = K(a), .p = b}); };
    check<double>("geometric pmf", d::geometric_pmf_double, [&](double a, double, double x) { return pmf(G(a), K(x)); }, tight);
    check<double>("geometric log_pmf", d::geometric_log_pmf_double, [&](double a, double, double x) { return log_pmf(G(a), K(x)); }, 1.0);
    check<double>("geometric cdf", d::geometric_cdf_double, [&](double a, double, double x) { return cdf(G(a), K(x)); }, tight);
    check<double>("geometric ccdf", d::geometric_ccdf_double, [&](double a, double, double x) { return ccdf(G(a), K(x)); }, tight);
    check<double>("poisson pmf", d::poisson_pmf_double, [&](double a, double, double x) { return pmf(P(a), K(x)); }, tight);
    check<double>("poisson log_pmf", d::poisson_log_pmf_double, [&](double a, double, double x) { return log_pmf(P(a), K(x)); }, 1.0);
    check<double>("binomial pmf", d::binomial_pmf_double, [&](double a, double b, double x) { return pmf(B(a, b), K(x)); }, tight);
    check<double>("binomial log_pmf", d::binomial_log_pmf_double, [&](double a, double b, double x) { return log_pmf(B(a, b), K(x)); }, 1.0);

    check<double>("poisson cdf", d::poisson_cdf_double, [](double a, double, double x) { return cdf(wsl::poisson<K>({.mean = a}), K(x)); }, tight);
    check<double>("poisson ccdf", d::poisson_ccdf_double, [](double a, double, double x) { return ccdf(wsl::poisson<K>({.mean = a}), K(x)); }, tight);
    check<double>("binomial cdf", d::binomial_cdf_double, [](double a, double b, double x) { return cdf(wsl::binomial<K>({.trials = K(a), .p = b}), K(x)); }, tight);
    check<double>("binomial ccdf", d::binomial_ccdf_double, [](double a, double b, double x) { return ccdf(wsl::binomial<K>({.trials = K(a), .p = b}), K(x)); }, tight);
    {
        std::size_t wrong = 0;
        for (auto const& q : d::poisson_quantile)
            wrong += quantile(wsl::poisson<K>({.mean = q.a}), q.u) != q.k;
        std::printf("%-28s %5zu points : %zu faux\n", "poisson quantile", std::size(d::poisson_quantile), wrong);
        failed |= wrong != 0;
        wrong = 0;
        for (auto const& q : d::binomial_quantile)
            wrong += quantile(wsl::binomial<K>({.trials = K(q.a), .p = q.b}), q.u) != q.k;
        std::printf("%-28s %5zu points : %zu faux\n", "binomial quantile", std::size(d::binomial_quantile), wrong);
        failed |= wrong != 0;
    }

    // Quantile de la géométrique : exact.
    std::size_t wrong = 0;
    for (auto const& q : d::geometric_quantile)
        wrong += quantile(G(q.p), q.u) != q.k;
    std::printf("%-28s %5zu points : %zu faux\n", "geometric quantile", std::size(d::geometric_quantile), wrong);
    failed |= wrong != 0;
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
