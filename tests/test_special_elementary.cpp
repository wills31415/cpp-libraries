// Précision et cas limites des fonctions élémentaires de wsl::special.
//
// Précision : erreur en ulp contre les valeurs de référence de mpmath
// (tests/data, engendrées par tools/special/generate.py), avec la position exacte
// de la valeur vraie entre deux flottants : l'erreur est mesurée, pas seulement
// « arrondi correct ou non ».

#include <wsl/special.hpp>

#include "data/special_elementary_reference.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <iterator>
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
static_assert(sp::sin(0.0) == 0.0 && sp::cos(0.0) == 1.0 && sp::atan(1.0) == 0x1.921fb54442d18p-1);
static_assert(sp::sin(1e300) == sp::sin(1e300) && sp::atan2(1.0, -1.0) == 0x1.2d97c7f3321d2p+1);
static_assert(sp::sqrt(2.0) == 0x1.6a09e667f3bcdp+0 && sp::sqrt(0x1p-1074) == 0x1p-537);
static_assert(sp::cbrt(27.0) == 3.0 && sp::cbrt(-8.0) == -2.0 && sp::hypot(3.0, 4.0) == 5.0);
static_assert(sp::erf(0.5) == 0x1.0a7ef5c18edd2p-1 && sp::erfc(10.0) == sp::erfc(10.0));
static_assert(sp::tgamma(5.0) == 24.0 && sp::lgamma(1.0) == 0.0 && sp::erf_inv(0.0) == 0.0 && sp::erfc_inv(1.0) == 0.0);

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
    check<double>("sin double", reference_data::sin_double, [](double x) { return sp::sin(x); }, 0.51);
    check<double>("cos double", reference_data::cos_double, [](double x) { return sp::cos(x); }, 0.51);
    check<double>("tan double", reference_data::tan_double, [](double x) { return sp::tan(x); }, 0.51);
    check<double>("atan double", reference_data::atan_double, [](double x) { return sp::atan(x); }, 0.51);
    check2<double>("atan2 double", reference_data::atan2_double, [](double y, double x) { return sp::atan2(y, x); }, 0.51);
    check<float>("sin float", reference_data::sin_float, [](float x) { return sp::sin(x); }, 0.5 + 0x1p-28);
    check<float>("cos float", reference_data::cos_float, [](float x) { return sp::cos(x); }, 0.5 + 0x1p-28);
    check<float>("tan float", reference_data::tan_float, [](float x) { return sp::tan(x); }, 0.5 + 0x1p-28);
    check<float>("atan float", reference_data::atan_float, [](float x) { return sp::atan(x); }, 0.5 + 0x1p-28);
    check2<float>("atan2 float", reference_data::atan2_float, [](float y, float x) { return sp::atan2(y, x); }, 0.5 + 0x1p-28);
    check<double>("cbrt double", reference_data::cbrt_double, [](double x) { return sp::cbrt(x); }, 0.51);
    check<float>("cbrt float", reference_data::cbrt_float, [](float x) { return sp::cbrt(x); }, 0.5 + 0x1p-28);
    check2<double>("hypot double", reference_data::hypot_double, [](double x, double y) { return sp::hypot(x, y); }, 0.51);
    check2<float>("hypot float", reference_data::hypot_float, [](float x, float y) { return sp::hypot(x, y); }, 0.5 + 0x1p-28);
    check<double>("erf double", reference_data::erf_double, [](double x) { return sp::erf(x); }, 0.51);
    check<double>("erfc double", reference_data::erfc_double, [](double x) { return sp::erfc(x); }, 0.51);
    check<float>("erf float", reference_data::erf_float, [](float x) { return sp::erf(x); }, 0.5 + 0x1p-28);
    check<float>("erfc float", reference_data::erfc_float, [](float x) { return sp::erfc(x); }, 0.5 + 0x1p-28);
    check<double>("lgamma x>0", reference_data::lgamma_positive_double, [](double x) { return sp::lgamma(x); }, 0.51);
    // lgamma, x < 0 : la réflexion retranche deux termes presque égaux près des
    // zéros de lgamma. Garantie : chaque résultat est à moins de 0,51 ulp, ou, à
    // défaut (près d'un zéro), à moins de 2^-58 en valeur absolue.
    {
        double relative = 0, absolute = 0;
        std::size_t beyond = 0;
        for (auto const& ref : reference_data::lgamma_negative_double)
        {
            double const y = sp::lgamma(ref.x);
            double const error_ulp = std::abs((y - ref.rounded) / ulp(ref.rounded) - ref.offset);
            if (error_ulp < 0.51)
                relative = std::max(relative, error_ulp);
            else
            {
                ++beyond;
                absolute = std::max(absolute, error_ulp * ulp(ref.rounded));
            }
        }
        std::printf("lgamma x<0   %5zu points : %zu au-delà de 0,51 ulp (près des zéros), erreur absolue max %a\n",
                    std::size(reference_data::lgamma_negative_double), beyond, absolute);
        assert(absolute < 0x1p-58);
    }
    check<double>("tgamma double", reference_data::tgamma_double, [](double x) { return sp::tgamma(x); }, 0.51);
    check<double>("digamma x>0", reference_data::digamma_double, [](double x) { return sp::digamma(x); }, 0.51);
    check<float>("digamma float", reference_data::digamma_float, [](float x) { return sp::digamma(x); }, 0.5 + 0x1p-28);
    // digamma, x < 0 : même garantie que lgamma (0,51 ulp, ou 2^-58 en absolu près des zéros).
    {
        double absolute = 0;
        std::size_t beyond = 0;
        for (auto const& ref : reference_data::digamma_negative_double)
        {
            double const y = sp::digamma(ref.x);
            double const error_ulp = std::abs((y - ref.rounded) / ulp(ref.rounded) - ref.offset);
            if (error_ulp >= 0.51)
            {
                ++beyond;
                absolute = std::max(absolute, error_ulp * ulp(ref.rounded));
            }
        }
        std::printf("digamma x<0  %5zu points : %zu au-delà de 0,51 ulp (près des zéros), erreur absolue max %a\n",
                    std::size(reference_data::digamma_negative_double), beyond, absolute);
        assert(absolute < 0x1p-58);
    }
    check<float>("lgamma float", reference_data::lgamma_float, [](float x) { return sp::lgamma(x); }, 0.5 + 0x1p-28);
    check<float>("tgamma float", reference_data::tgamma_float, [](float x) { return sp::tgamma(x); }, 0.5 + 0x1p-28);
    check<double>("erf_inv dbl", reference_data::erf_inv_double, [](double y) { return sp::erf_inv(y); }, 0.51);
    check<double>("erfc_inv dbl", reference_data::erfc_inv_double, [](double z) { return sp::erfc_inv(z); }, 0.51);
    check<float>("erf_inv float", reference_data::erf_inv_float, [](float y) { return sp::erf_inv(y); }, 0.5 + 0x1p-28);
    check<float>("erfc_inv flt", reference_data::erfc_inv_float, [](float z) { return sp::erfc_inv(z); }, 0.5 + 0x1p-28);

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

    // Trigonométrie : cas à part comparés à <cmath>, signe des zéros compris.
    {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        auto const same = [](double ours, double theirs) {
            return (std::isnan(ours) && std::isnan(theirs)) || (ours == theirs && std::signbit(ours) == std::signbit(theirs));
        };
        double const singles[] = {0.0, -0.0, inf, -inf, nan, 0x1p-1074, -0x1p-1074, 0x1p-30, -0x1p-30};
        for (double x : singles)
        {
            assert(same(sp::sin(x), std::sin(x)) && same(sp::cos(x), std::cos(x)));
            assert(same(sp::tan(x), std::tan(x)) && same(sp::atan(x), std::atan(x)));
        }
        assert(sp::atan(inf) == 0x1.921fb54442d18p+0 && sp::atan(-inf) == -0x1.921fb54442d18p+0);
        double const values[] = {0.0, -0.0, 1.0, -1.0, 0x1p-1074, -0x1p-1074, 1e300, -1e300, inf, -inf, nan};
        for (double y : values)
            for (double x : values)
            {
                double const ours = sp::atan2(y, x), theirs = std::atan2(y, x);
                bool const ok = same(ours, theirs) || std::abs(ours - theirs) <= ulp(theirs);
                if (!ok)
                    std::printf("atan2(%a, %a) : %a, <cmath> %a\n", y, x, ours, theirs);
                assert(ok);
            }
    }
    // atan2 aux rapports extrêmes, autour des seuils de mise à l'échelle : jamais
    // NaN, et comme <cmath> (un premier essai rendait NaN pour y ≈ 2^-835,
    // x ≈ 2^797, le découpage de Veltkamp débordant).
    for (double big : {0x1p650, 0x1p699, 0x1p701, 0x1p797, 0x1p799, 0x1p801, 0x1p899, 0x1p901, 0x1p1000})
        for (double small : {0x1p-1074, 0x1p-1000, 0x1p-835, 0x1p-801, 0x1p-799, 0x1p-500})
            for (double sy : {1.0, -1.0})
                for (double sx : {1.0, -1.0})
                {
                    double const a = sp::atan2(sy * small, sx * big), b = sp::atan2(sy * big, sx * small);
                    assert(!std::isnan(a) && !std::isnan(b));
                    assert(std::abs(a - std::atan2(sy * small, sx * big)) <= ulp(std::atan2(sy * small, sx * big)));
                    assert(std::abs(b - std::atan2(sy * big, sx * small)) <= ulp(std::atan2(sy * big, sx * small)));
                }

    // Identités : sin² + cos² ≈ 1, tan = sin/cos, atan(tan x) = x sur (−π/2, π/2),
    // symétries exactes.
    for (double x = -50; x < 50; x += 0.0137)
    {
        double const s = sp::sin(x), c = sp::cos(x);
        assert(std::abs(s * s + c * c - 1) <= 4e-16);
        assert(sp::sin(-x) == -s && sp::cos(-x) == c && sp::tan(-x) == -sp::tan(x));
        if (std::abs(c) > 1e-3)
            assert(std::abs(sp::tan(x) - s / c) <= 4 * ulp(s / c));
    }
    for (double x = -1.5; x < 1.5; x += 0.00731)
        assert(std::abs(sp::atan(sp::tan(x)) - x) <= 2 * ulp(x) + 1e-300);
    assert(sp::sin(0x1.921fb54442d18p+1) == 0x1.1a62633145c07p-53);   // sin(π arrondi) = π − π arrondi
    assert(sp::cos(0x1.921fb54442d18p+0) == 0x1.1a62633145c07p-54);

    // sqrt : identique à la racine matérielle ; cbrt et hypot : cas exacts et à part.
    for (double x = 1e-300; x < 1e300; x *= 1.7)
        assert(sp::sqrt(x) == std::sqrt(x));
    assert(std::signbit(sp::sqrt(-0.0)) && std::isnan(sp::sqrt(-1.0)) && sp::sqrt(inf) == inf);
    for (int n = -100; n <= 100; ++n)
        assert(sp::cbrt(double(n) * n * n) == n);
    assert(std::signbit(sp::cbrt(-0.0)) && sp::cbrt(-inf) == -inf && std::isnan(sp::cbrt(std::nan(""))));
    assert(sp::cbrt(0x1p-1074) == 0x1p-358);                     // exact (1074 = 3 · 358)
    assert(sp::hypot(3.0, 4.0) == 5.0 && sp::hypot(-5.0, 12.0) == 13.0 && sp::hypot(-7.0, 0.0) == 7.0);
    assert(sp::hypot(inf, std::nan("")) == inf && sp::hypot(std::nan(""), -inf) == inf);
    assert(std::isnan(sp::hypot(std::nan(""), 1.0)));
    // Valeurs exactes (mpmath) : les libm ne sont pas toutes correctement arrondies.
    assert(sp::hypot(1e308, 1e308) == 0x1.92c80954c51f5p+1023 && sp::hypot(1.7e308, 1.7e308) == inf);
    assert(sp::hypot(0x1p-1074, 0x1p-1074) == 0x1p-1074);

    // erf, erfc : cas à part, symétries, compléments.
    assert(sp::erf(0.0) == 0.0 && std::signbit(sp::erf(-0.0)) && sp::erf(inf) == 1.0 && sp::erf(-inf) == -1.0);
    assert(sp::erfc(0.0) == 1.0 && sp::erfc(-0.0) == 1.0 && is_positive_zero(sp::erfc(inf)) && sp::erfc(-inf) == 2.0);
    assert(std::isnan(sp::erf(std::nan(""))) && std::isnan(sp::erfc(std::nan(""))));
    assert(sp::erf(0x1p-1074) == 0x1p-1074 && sp::erf(1e-300) == sp::erf(1e-300));
    assert(is_positive_zero(sp::erfc(28.0)) && sp::erfc(27.0) > 0 && sp::erf(6.0) == 1.0 && sp::erfc(-6.0) == 2.0);
    for (double x = -7; x < 7; x += 0.0123)
    {
        assert(sp::erf(-x) == -sp::erf(x));
        assert(std::abs(sp::erf(x) + sp::erfc(x) - 1) <= 3e-16);
        assert(std::abs(sp::erfc(x) + sp::erfc(-x) - 2) <= 6e-16);
    }

    // lgamma, tgamma : cas à part, factorielles exactes, signes.
    assert(is_positive_zero(sp::lgamma(1.0)) && is_positive_zero(sp::lgamma(2.0)));
    assert(sp::lgamma(0.0) == inf && sp::lgamma(-0.0) == inf && sp::lgamma(-3.0) == inf && sp::lgamma(-inf) == inf);
    assert(sp::lgamma(inf) == inf && sp::lgamma(2.6e305) == inf && std::isnan(sp::lgamma(std::nan(""))));
    assert(sp::tgamma(0.0) == inf && sp::tgamma(-0.0) == -inf && std::isnan(sp::tgamma(-3.0)));
    assert(std::isnan(sp::tgamma(-inf)) && sp::tgamma(inf) == inf && sp::tgamma(172.0) == inf);
    double factorial = 1;
    for (int n = 1; n <= 23; ++n)                                    // (n−1)! exact en double jusqu'à 22!
    {
        assert(sp::tgamma(double(n)) == factorial);
        factorial *= n;
    }
    assert(sp::tgamma(0.5) == 0x1.c5bf891b4ef6bp+0 && sp::tgamma(-0.5) == -0x1.c5bf891b4ef6bp+1);
    assert(sp::tgamma(-1.5) == 0x1.2e7fb0bcdf4f2p+1 && sp::lgamma(0.5) == 0x1.250d048e7a1bdp-1);
    for (double x = -20.75; x < 0; x += 0.5)                         // signe : négatif sur (−1, 0), (−3, −2)…
    {
        bool const negative = (static_cast<long long>(std::floor(x)) & 1) != 0;
        assert((sp::tgamma(x) < 0) == negative);
    }
    assert(sp::tgamma(1e-300) == 1.0 / 1e-300 && sp::tgamma(-1e-300) == -1.0 / 1e-300 && sp::lgamma(1e-300) == -sp::log(1e-300));

    // erf_inv, erfc_inv : cas à part, allers-retours.
    assert(sp::erf_inv(0.0) == 0.0 && std::signbit(sp::erf_inv(-0.0)));
    assert(sp::erf_inv(1.0) == inf && sp::erf_inv(-1.0) == -inf && std::isnan(sp::erf_inv(1.5)));
    assert(std::isnan(sp::erf_inv(std::nan(""))) && std::isnan(sp::erf_inv(-inf)));
    assert(sp::erfc_inv(0.0) == inf && sp::erfc_inv(2.0) == -inf && is_positive_zero(sp::erfc_inv(1.0)));
    assert(std::isnan(sp::erfc_inv(-0.5)) && std::isnan(sp::erfc_inv(2.5)) && std::isnan(sp::erfc_inv(std::nan(""))));
    assert(sp::erf_inv(0x1p-1074) == 0x1p-1074 && sp::erfc_inv(0x1p-1074) > 27 && sp::erfc_inv(0x1p-1074) < 27.3);
    for (double y = -0.999; y < 1; y += 0.00137)
        assert(std::abs(sp::erf(sp::erf_inv(y)) - y) <= 2 * ulp(y) + 1e-17);
    for (double z = 1e-300; z < 0.5; z *= 1.9)
    {
        // Conditionnement : une erreur relative ε sur x en devient une de 2x² ε sur z.
        double const x = sp::erfc_inv(z);
        assert(std::abs(sp::erfc(x) / z - 1) <= 4 * 0x1p-53 * (2 * x * x + 1));
    }

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
