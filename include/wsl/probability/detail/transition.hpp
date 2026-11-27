#pragma once

// Transition (refonte, étape 5a) : random est découpée en wsl::random (moteurs, primitives,
// tirages exacts), wsl::probability (lois, variables, algèbre des lois),
// wsl::probability::stochastic_processes et wsl::measure_theory (design/interfaces.md,
// décisions 37, 39, 40).
// Le code déplacé garde ses noms non qualifiés : ces déclarations les rendent visibles
// jusqu'à la passe 5b, qui les appellera qualifiés. Rien ici n'est une interface.

#include <wsl/measure_theory/reference_measures.hpp>
#include <wsl/modes.hpp>
#include <wsl/random/detail/rational.hpp>
#include <wsl/random/engines.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

namespace wsl::inline WSL_DEFAULTS::random::detail
{
}

namespace wsl::inline WSL_DEFAULTS::probability
{
    using random::bernoulli_trial;
    using random::ratio;
    using random::random_bits;
    using random::shuffle;
    using random::uniform_below;
    using random::uniform_int;
    using random::uniform_real;
    using random::uniform_real_fast;

    using measure_theory::counting;
    using measure_theory::lebesgue;
    using measure_theory::surface;

    namespace laws
    {
    }
    using namespace laws;

    namespace detail
    {
        using namespace random::detail;
        using namespace random::engines;
        using namespace special::detail;
    }
    namespace stochastic_processes
    {
        using namespace laws;
    }
}
