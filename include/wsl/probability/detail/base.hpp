#pragma once

// Socle commun des en-têtes de probability : les couches du dessous (moteurs, primitives et
// tirages exacts de random, calcul flottant portable, mesures de référence), et l'accès des
// outils internes de probability à ceux de random. Rien ici n'est une interface.

#include <wsl/measure_theory/reference_measures.hpp>
#include <wsl/modes.hpp>
#include <wsl/random/detail/portable_math.hpp>
#include <wsl/random/detail/rational.hpp>
#include <wsl/random/engines.hpp>
#include <wsl/random/uniform.hpp>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
}

// Les outils internes de probability voient ceux de random, et les lois (types et concepts
// qu'ils manipulent).
namespace wsl::inline WSL_DEFAULTS::probability::detail
{
    using namespace random::detail;
    using namespace random::engines;
    using namespace laws;
}
