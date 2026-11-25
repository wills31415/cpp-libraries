// wsl::mdarray, mode reproductible par défaut (WSL_MDARRAY_REPRODUCIBLE, I9) : sans
// marque, les instructions sont reproductibles, même empreinte que wsl::reproducible ;
// wsl::fast(…) reste possible. (Un exécutable à part : le mode par défaut doit être le
// même dans tout un programme.)

#define WSL_MDARRAY_REPRODUCIBLE 1
#include "test_mdarray_modes_common.hpp"

#include <cinttypes>
#include <cstdio>

int main()
{
    std::uint64_t const unmarked = modes::compute(modes::unmarked);
    std::uint64_t const explicit_mode = modes::compute(modes::reproducible);
    (void)modes::compute(modes::fast);
    if (unmarked != explicit_mode || unmarked != 0xf33b15d552bca1bcULL)
    {
        std::printf("empreintes différentes : 0x%016" PRIx64 ", 0x%016" PRIx64 "\n", unmarked, explicit_mode);
        return 1;
    }
    std::puts("test_mdarray_modes_default : OK");
}
