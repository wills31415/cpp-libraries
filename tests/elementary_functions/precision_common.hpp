#pragma once

// Mesure de précision commune aux tests de précision (elementary_functions, special_functions) :
// erreur en ulp contre les valeurs de référence de mpmath, avec la position exacte de la valeur
// vraie entre deux flottants.

#include "data/special_elementary_reference.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <span>

using reference_data::reference;
using reference_data::reference2;
using reference_data::reference3;

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

    template <class T, class F>
    void check3(char const* name, std::span<reference3<T> const> references, F f, double bound)
    {
        report r;
        for (auto const& ref : references)
        {
            T const v = f(ref.a, ref.b, ref.x);
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

    [[maybe_unused]] bool is_positive_zero(double x) { return x == 0 && !std::signbit(x); }
}
