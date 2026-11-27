// Étiquettes : precision long
// Précision des fonctions élémentaires dans les deux modes (design/interfaces.md, décisions 24
// à 27, 38) : déterministe (défaut ; le compilateur peut fusionner des produits) et portable
// (mêmes bits partout). Une fusion supprime un arrondi : les bornes, établies opération par
// opération, doivent tenir dans les deux modes. Les tests compilés avec -ffp-contract=fast
// (configurations FMA) l'éprouvent vraiment.

#include <wsl/elementary_functions.hpp>

#include "precision_common.hpp"

namespace ef = wsl::elementary_functions;

namespace
{
    template <class Mode>
    void check_mode(char const* mode)
    {
        std::printf("--- mode %s\n", mode);
        check<double>("exp double", reference_data::exp_double, [](double x) { return ef::exp(x, Mode{}); }, 0.51);
        check<float>("exp float", reference_data::exp_float, [](float x) { return ef::exp(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("exp2 double", reference_data::exp2_double, [](double x) { return ef::exp2(x, Mode{}); }, 0.51);
        check<float>("exp2 float", reference_data::exp2_float, [](float x) { return ef::exp2(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("expm1 double", reference_data::expm1_double, [](double x) { return ef::expm1(x, Mode{}); }, 0.51);
        check<float>("expm1 float", reference_data::expm1_float, [](float x) { return ef::expm1(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("log double", reference_data::log_double, [](double x) { return ef::log(x, Mode{}); }, 0.51);
        check<float>("log float", reference_data::log_float, [](float x) { return ef::log(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("log2 double", reference_data::log2_double, [](double x) { return ef::log2(x, Mode{}); }, 0.51);
        check<float>("log2 float", reference_data::log2_float, [](float x) { return ef::log2(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("log10 double", reference_data::log10_double, [](double x) { return ef::log10(x, Mode{}); }, 0.51);
        check<float>("log10 float", reference_data::log10_float, [](float x) { return ef::log10(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("log1p double", reference_data::log1p_double, [](double x) { return ef::log1p(x, Mode{}); }, 0.51);
        check<float>("log1p float", reference_data::log1p_float, [](float x) { return ef::log1p(x, Mode{}); }, 0.5 + 0x1p-28);
        check2<double>("pow double", reference_data::pow_double, [](double x, double y) { return ef::pow(x, y, Mode{}); }, 0.51);
        check2<float>("pow float", reference_data::pow_float, [](float x, float y) { return ef::pow(x, y, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("sin double", reference_data::sin_double, [](double x) { return ef::sin(x, Mode{}); }, 0.51);
        check<float>("sin float", reference_data::sin_float, [](float x) { return ef::sin(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("cos double", reference_data::cos_double, [](double x) { return ef::cos(x, Mode{}); }, 0.51);
        check<float>("cos float", reference_data::cos_float, [](float x) { return ef::cos(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("tan double", reference_data::tan_double, [](double x) { return ef::tan(x, Mode{}); }, 0.51);
        check<float>("tan float", reference_data::tan_float, [](float x) { return ef::tan(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("atan double", reference_data::atan_double, [](double x) { return ef::atan(x, Mode{}); }, 0.51);
        check<float>("atan float", reference_data::atan_float, [](float x) { return ef::atan(x, Mode{}); }, 0.5 + 0x1p-28);
        check<double>("cbrt double", reference_data::cbrt_double, [](double x) { return ef::cbrt(x, Mode{}); }, 0.51);
        check<float>("cbrt float", reference_data::cbrt_float, [](float x) { return ef::cbrt(x, Mode{}); }, 0.5 + 0x1p-28);
        check2<double>("atan2 double", reference_data::atan2_double, [](double y, double x) { return ef::atan2(y, x, Mode{}); }, 0.51);
        check2<float>("atan2 float", reference_data::atan2_float, [](float y, float x) { return ef::atan2(y, x, Mode{}); }, 0.5 + 0x1p-28);
        check2<double>("hypot double", reference_data::hypot_double, [](double x, double y) { return ef::hypot(x, y, Mode{}); }, 0.51);
        check2<float>("hypot float", reference_data::hypot_float, [](float x, float y) { return ef::hypot(x, y, Mode{}); }, 0.5 + 0x1p-28);
    }
}

// Les étiquettes : défaut, explicites, refus de wsl::exact (pas l'arrondi correct).
static_assert(ef::exp(0.0) == 1.0 && ef::exp(0.0, wsl::portable) == 1.0 && ef::exp(0.0, wsl::deterministic) == 1.0);
static_assert(ef::exp(1.0, wsl::portable) == 0x1.5bf0a8b145769p+1);
static_assert(ef::exp(1.0, wsl::portable, wsl::approximate) == ef::exp(1.0, wsl::approximate, wsl::portable));
static_assert(ef::exp2(10.0) == 1024.0 && ef::expm1(0x1p-60, wsl::portable) == 0x1p-60);
template <class... M>
constexpr bool exp_accepts = requires { ef::exp(1.0, M{}...); };
static_assert(exp_accepts<> && exp_accepts<wsl::portable_t> && exp_accepts<wsl::deterministic_t, wsl::approximate_t>);
static_assert(!exp_accepts<wsl::portable_t, wsl::deterministic_t>);   // deux étiquettes d'un même axe
static_assert(!exp_accepts<int>);                                     // pas une étiquette
template <class... M>
constexpr bool sqrt_accepts = requires { ef::sqrt(2.0, M{}...); };
static_assert(sqrt_accepts<wsl::exact_t> && sqrt_accepts<wsl::exact_t, wsl::portable_t>);   // sqrt l'est
static_assert(ef::sqrt(2.0, wsl::exact) == 0x1.6a09e667f3bcdp+0);

int main()
{
    check_mode<wsl::deterministic_t>("déterministe");
    check_mode<wsl::portable_t>("portable");
    std::puts("elementary_functions.precision : OK");
}
