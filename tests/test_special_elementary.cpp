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
using reference_data::reference2;

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

    // Fonctions à deux arguments.
    template <class T, class F>
    void check2(char const* name, std::span<reference2<T> const> references, F f, double bound)
    {
        report r;
        for (auto const& ref : references)
        {
            T const v = f(ref.x, ref.y);
            double const error = std::abs(double(v - ref.rounded) / double(ulp(ref.rounded)) - ref.offset);
            if (v != ref.rounded)
                ++r.not_rounded;
            if (error > r.max_error)
            {
                r.max_error = error;
                r.worst_x = double(ref.x);
            }
            ++r.points;
        }
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
static_assert(sp::exp2(10.0) == 1024.0 && sp::log2(1024.0) == 10.0 && sp::log10(1000.0) == 3.0);
static_assert(sp::expm1(0x1p-60) == 0x1p-60 && sp::log1p(0x1p-60) == 0x1p-60);
static_assert(sp::pow(2.0, 10.0) == 1024.0 && sp::pow(-3.0, 3.0) == -27.0);

constexpr double compile_time[] = {sp::exp(-3.7), sp::exp(123.456), sp::exp(-720.5), sp::log(0x1p-1070),
                                   sp::log(3.0), sp::log(0.999999), sp::log(1e300)};

int main()
{
    // Précision. Bornes : celles annoncées dans elementary.hpp.
    check<double>("exp double", reference_data::exp_double, [](double x) { return sp::exp(x); }, 0.51);
    check<double>("log double", reference_data::log_double, [](double x) { return sp::log(x); }, 0.51);
    check<float>("exp float", reference_data::exp_float, [](float x) { return sp::exp(x); }, 0.5 + 0x1p-28);
    check<float>("log float", reference_data::log_float, [](float x) { return sp::log(x); }, 0.5 + 0x1p-28);
    check<double>("exp2 double", reference_data::exp2_double, [](double x) { return sp::exp2(x); }, 0.51);
    check<double>("expm1 double", reference_data::expm1_double, [](double x) { return sp::expm1(x); }, 0.51);
    check<double>("log2 double", reference_data::log2_double, [](double x) { return sp::log2(x); }, 0.51);
    check<double>("log10 double", reference_data::log10_double, [](double x) { return sp::log10(x); }, 0.51);
    check<double>("log1p double", reference_data::log1p_double, [](double x) { return sp::log1p(x); }, 0.51);
    check<float>("exp2 float", reference_data::exp2_float, [](float x) { return sp::exp2(x); }, 0.5 + 0x1p-28);
    check<float>("expm1 float", reference_data::expm1_float, [](float x) { return sp::expm1(x); }, 0.5 + 0x1p-28);
    check<float>("log2 float", reference_data::log2_float, [](float x) { return sp::log2(x); }, 0.5 + 0x1p-28);
    check<float>("log10 float", reference_data::log10_float, [](float x) { return sp::log10(x); }, 0.5 + 0x1p-28);
    check<float>("log1p float", reference_data::log1p_float, [](float x) { return sp::log1p(x); }, 0.5 + 0x1p-28);
    check2<double>("pow double", reference_data::pow_double, [](double x, double y) { return sp::pow(x, y); }, 0.51);
    check2<float>("pow float", reference_data::pow_float, [](float x, float y) { return sp::pow(x, y); }, 0.5 + 0x1p-28);

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

    // exp2 : puissances exactes, bornes.
    for (int n = -1074; n < 1024; ++n)
        assert(sp::exp2(double(n)) == std::ldexp(1.0, n));
    assert(sp::exp2(1024.0) == inf && is_positive_zero(sp::exp2(-1080.0)) && is_positive_zero(sp::exp2(-inf)));
    assert(sp::exp2(inf) == inf && std::isnan(sp::exp2(std::numeric_limits<double>::quiet_NaN())));
    assert(sp::exp2(0.5) == 0x1.6a09e667f3bcdp+0);              // √2, arrondi correct

    // expm1 : signe des zéros, petits arguments, bornes.
    assert(is_positive_zero(sp::expm1(0.0)) && std::signbit(sp::expm1(-0.0)) && sp::expm1(-0.0) == 0);
    assert(sp::expm1(0x1p-60) == 0x1p-60 && sp::expm1(-0x1p-60) == -0x1p-60);
    assert(sp::expm1(0x1p-1074) == 0x1p-1074);
    assert(sp::expm1(-inf) == -1.0 && sp::expm1(-40.0) == -1.0 && sp::expm1(inf) == inf);
    assert(sp::expm1(1000.0) == inf);

    // log2, log10 : puissances exactes.
    for (int n = -1074; n < 1024; ++n)
        assert(sp::log2(std::ldexp(1.0, n)) == n);
    double power_of_ten = 1;
    for (int n = 0; n <= 22; ++n, power_of_ten *= 10)            // 10^n exact en double jusqu'à 10^22
        assert(sp::log10(power_of_ten) == n);
    assert(sp::log2(0.0) == -inf && std::isnan(sp::log2(-1.0)) && sp::log10(inf) == inf);
    assert(is_positive_zero(sp::log2(1.0)) && is_positive_zero(sp::log10(1.0)));

    // log1p : signe des zéros, petits arguments, bornes.
    assert(is_positive_zero(sp::log1p(0.0)) && std::signbit(sp::log1p(-0.0)) && sp::log1p(-0.0) == 0);
    assert(sp::log1p(0x1p-60) == 0x1p-60 && sp::log1p(-0x1p-1074) == -0x1p-1074);
    assert(sp::log1p(-1.0) == -inf && std::isnan(sp::log1p(-1.5)) && sp::log1p(inf) == inf);

    // pow : cas à part de l'Annexe F de C, comparés à <cmath> sur toutes les
    // combinaisons de valeurs remarquables (le signe des zéros compris).
    {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        double const values[] = {0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 3.0, -3.0, 2.5, -2.5, 0x1p-1074,
                                 1e300, -1e300, inf, -inf, nan, 0x1p53, -0x1p53, 0x1p53 + 2, 1e20, -1e20};
        for (double x : values)
            for (double y : values)
            {
                double const ours = sp::pow(x, y), theirs = std::pow(x, y);
                bool const same = (std::isnan(ours) && std::isnan(theirs))
                               || (ours == theirs && std::signbit(ours) == std::signbit(theirs))
                               || (std::isfinite(theirs) && theirs != 0 && std::abs(ours - theirs) <= ulp(theirs));
                if (!same)
                    std::printf("pow(%a, %a) : %a, <cmath> %a\n", x, y, ours, theirs);
                assert(same);
            }
    }
    assert(sp::pow(2.0, 10.0) == 1024.0 && sp::pow(2.0, -1074.0) == 0x1p-1074 && sp::pow(10.0, 22.0) == 1e22);
    assert(sp::pow(-2.0, 3.0) == -8.0 && sp::pow(-2.0, 2.0) == 4.0 && std::isnan(sp::pow(-2.0, 0.5)));
    assert(sp::pow(4.0, 0.5) == 2.0 && sp::pow(0.25, -0.5) == 2.0);
    for (double x = 0.001; x < 1000; x *= 1.37)
        assert(sp::pow(x, 1.0) == x);
    assert(sp::pow(2, 3) == 8.0 && sp::pow(2.0f, 3.0f) == 8.0f && sp::pow(2.0f, 0.5) == sp::pow(2.0, 0.5));

    // Entiers, comme <cmath>.
    assert(sp::exp(0) == 1.0 && sp::log(1) == 0.0);

    // Cohérence : log(exp(x)) ≈ x, log1p(expm1(x)) ≈ x, log2(exp2(x)) ≈ x.
    for (double x = -700; x < 700; x += 0.37)
    {
        assert(std::abs(sp::log(sp::exp(x)) - x) <= 4 * ulp(x) + 1e-15);
        assert(std::abs(sp::log2(sp::exp2(x)) - x) <= 4 * ulp(x) + 1e-15);
    }
    // (Pour x < −1, expm1(x) est proche de −1 et ne garde que peu de bits de
    // e^x : l'aller-retour y est mal conditionné, quelle que soit la précision.)
    for (double x = -1; x < 700; x += 0.173)
        assert(std::abs(sp::log1p(sp::expm1(x)) - x) <= 4 * ulp(x));

    return 0;
}
