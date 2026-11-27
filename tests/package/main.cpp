// Les deux écritures des inclusions (docs/interfaces.md, décision 30) : sans extension,
// la forme documentée, et avec .hpp.
#include <wsl/array2d>
#include <wsl/graded_io.hpp>
#include <wsl/random>                  // parapluie : fichier et dossier du même nom
#include <wsl/random/engines>          // notion sous un parapluie
#include <wsl/special/gamma>
#include <wsl/mdarray>

#include <cstdlib>

int main()
{
    static_assert(__cplusplus > 202002L, "wsl::array2d doit imposer C++23");
    auto m = wsl::array2d<int, 2, 2>::from_rows({{1, 2}, {3, 4}});
    wsl::random::pcg64 g(42);
    bool const ok = m.transposed()[0, 1] == 3 && wsl::to_string(wsl::jet<int, 3>{1, 2}) == "{1, 2}"
                 && wsl::special::lgamma(1.0) == 0.0 && g() != g();
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
