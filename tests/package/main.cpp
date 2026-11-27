// Les deux écritures des inclusions (design/interfaces.md, décision 30) : sans extension,
// la forme documentée, et avec .hpp.
#include <wsl/mdarray>                 // parapluie : fichier et dossier du même nom
#include <wsl/permutation.hpp>
#include <wsl/random>
#include <wsl/random/engines>          // notion sous un parapluie
#include <wsl/special_functions/gamma>

#include <cstdlib>

int main()
{
    static_assert(__cplusplus > 202002L, "wsl::mdarray doit imposer C++23");
    wsl::random::engines::pcg64 g(42);
    bool const ok = wsl::special_functions::lgamma(1.0) == 0.0 && g() != g();
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
