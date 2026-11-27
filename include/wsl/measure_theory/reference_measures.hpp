#pragma once

// Mesures de référence (design/interfaces.md, décision 23) : la densité d'une loi est sa
// dérivée de Radon-Nikodym par rapport à l'une d'elles. Pour l'instant de simples étiquettes ;
// la théorie de la mesure (intégrale, mesures à densité, produits, images : décision 39) reste
// à concevoir.

#include <wsl/modes.hpp>

namespace wsl::inline WSL_DEFAULTS::measure_theory
{
    // Mesure de Lebesgue (longueur, aire, volume).
    struct lebesgue
    {
    };
    // Mesure de comptage.
    struct counting
    {
    };
    // Mesure de surface (de Hausdorff) d'une sphère.
    struct surface
    {
    };
}
