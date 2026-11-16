// Précision et cas limites des fonctions élémentaires de wsl::special.
//
// Précision : erreur en ulp contre les valeurs de référence de mpmath
// (tests/data, engendrées par tools/special/generate.py), avec la position exacte
// de la valeur vraie entre deux flottants : l'erreur est mesurée, pas seulement
// « arrondi correct ou non ».

#include <wsl/special.hpp>

#include "data/special_elementary_reference.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <span>

namespace sp = wsl::special;
using reference_data::reference;

namespace
{
    template <class T>
    T ulp(T x)
    {
        T const a = std::abs(x);
        return std::nextafter(a, std::numeric_limits<T>::infinity()) - a;
    }

    struct report
    {
        double max_error = 0;
        double worst_x = 0;
        std::size_t not_rounded = 0;   // résultats différents de l'arrondi correct
        std::size_t points = 0;
    };

    // Erreur de f(x) en ulp de l'arrondi correct : |f(x) − valeur exacte| / ulp.
    template <class T, class F>
    report measure(std::span<reference<T> const> references, F f)
    {
        report r;
        for (auto const& ref : references)
        {
            T const y = f(ref.x);
            double const error = std::abs(double(y - ref.rounded) / double(ulp(ref.rounded)) - ref.offset);
            if (y != ref.rounded)
                ++r.not_rounded;
            if (error > r.max_error)
            {
                r.max_error = error;
                r.worst_x = double(ref.x);
            }
            ++r.points;
        }
        return r;
    }

    template <class T, class F>
    void check(char const* name, std::span<reference<T> const> references, F f, double bound)
    {
        report const r = measure(references, f);
        std::printf("%-12s %5zu points : erreur max %.4f ulp (x = %a), arrondi correct %.3f %%\n", name,
                    r.points, r.max_error, r.worst_x,
                    100.0 * double(r.points - r.not_rounded) / double(r.points));
        assert(r.max_error < bound);
    }

    constexpr double inf = std::numeric_limits<double>::infinity();
    constexpr float inff = std::numeric_limits<float>::infinity();

    bool is_positive_zero(double x) { return x == 0 && !std::signbit(x); }
}

// ---------------------------------------------------------------------------
// Évaluation à la compilation : mêmes résultats qu'à l'exécution.
// ---------------------------------------------------------------------------
static_assert(sp::exp(0.0) == 1.0);
static_assert(sp::exp(1.0) == 0x1.5bf0a8b145769p+1);    // e, arrondi correct
static_assert(sp::log(1.0) == 0.0);
static_assert(sp::log(2.0) == 0x1.62e42fefa39efp-1);    // ln 2, arrondi correct
static_assert(sp::exp(0.0f) == 1.0f);
static_assert(sp::log(1.0f) == 0.0f);

constexpr double compile_time[] = {sp::exp(-3.7), sp::exp(123.456), sp::exp(-720.5), sp::log(0x1p-1070),
                                   sp::log(3.0), sp::log(0.999999), sp::log(1e300)};

int main()
{
    // Précision. Bornes : celles annoncées dans elementary.hpp.
    check<double>("exp double", reference_data::exp_double, [](double x) { return sp::exp(x); }, 0.51);
    check<double>("log double", reference_data::log_double, [](double x) { return sp::log(x); }, 0.51);
    check<float>("exp float", reference_data::exp_float, [](float x) { return sp::exp(x); }, 0.5 + 0x1p-28);
    check<float>("log float", reference_data::log_float, [](float x) { return sp::log(x); }, 0.5 + 0x1p-28);

    // L'évaluation à la compilation donne les mêmes bits qu'à l'exécution.
    volatile double inputs[] = {-3.7, 123.456, -720.5, 0x1p-1070, 3.0, 0.999999, 1e300};
    for (int i = 0; i < 3; ++i)
        assert(sp::exp(double(inputs[i])) == compile_time[i]);
    for (int i = 3; i < 7; ++i)
        assert(sp::log(double(inputs[i])) == compile_time[i]);

    // exp : cas limites (Annexe F de C).
    assert(sp::exp(0.0) == 1.0 && sp::exp(-0.0) == 1.0);
    assert(sp::exp(inf) == inf);
    assert(is_positive_zero(sp::exp(-inf)));
    assert(std::isnan(sp::exp(std::numeric_limits<double>::quiet_NaN())));
    assert(sp::exp(0x1.62e42fefa39efp+9) == 0x1.fffffffffff2ap+1023);   // plus grand x fini
    assert(sp::exp(std::nextafter(0x1.62e42fefa39efp+9, inf)) == inf);
    assert(sp::exp(1000.0) == inf);
    assert(sp::exp(-745.0) == 0x1p-1074);                 // plus petit sous-normal
    assert(is_positive_zero(sp::exp(-746.0)));
    assert(is_positive_zero(sp::exp(-1000.0)));
    assert(sp::exp(0x1p-60) == 1.0 && sp::exp(-0x1p-60) == 1.0);
    assert(sp::exp(88.8f) == inff && sp::exp(-104.0f) == 0.0f && sp::exp(-inff) == 0.0f);
    assert(std::isnan(sp::exp(std::numeric_limits<float>::quiet_NaN())));

    // log : cas limites.
    assert(is_positive_zero(sp::log(1.0)));
    assert(sp::log(0.0) == -inf && sp::log(-0.0) == -inf);
    assert(std::isnan(sp::log(-1.0)) && std::isnan(sp::log(-inf)));
    assert(sp::log(inf) == inf);
    assert(std::isnan(sp::log(std::numeric_limits<double>::quiet_NaN())));
    assert(sp::log(0x1p-1074) == -0x1.74385446d71c3p+9);   // plus petit sous-normal
    assert(sp::log(std::numeric_limits<double>::max()) == 0x1.62e42fefa39efp+9);
    assert(sp::log(0.0f) == -inff && std::isnan(sp::log(-1.0f)) && sp::log(inff) == inff);

    // Entiers, comme <cmath>.
    assert(sp::exp(0) == 1.0 && sp::log(1) == 0.0);

    // Cohérence : log(exp(x)) ≈ x.
    for (double x = -700; x < 700; x += 0.37)
        assert(std::abs(sp::log(sp::exp(x)) - x) <= 4 * ulp(x) + 1e-15);

    return 0;
}
