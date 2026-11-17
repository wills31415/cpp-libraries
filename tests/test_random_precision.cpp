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

    // Quantile de la géométrique : exact.
    std::size_t wrong = 0;
    for (auto const& q : d::geometric_quantile)
        wrong += quantile(G(q.p), q.u) != q.k;
    std::printf("%-28s %5zu points : %zu faux\n", "geometric quantile", std::size(d::geometric_quantile), wrong);
    failed |= wrong != 0;
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
