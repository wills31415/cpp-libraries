#include <wsl/array2d.hpp>
#include <wsl/graded_io.hpp>

#include <cstdlib>

int main()
{
    static_assert(__cplusplus > 202002L, "wsl::array2d doit imposer C++23");
    auto m = wsl::array2d<int, 2, 2>::from_rows({{1, 2}, {3, 4}});
    bool const ok = m.transposed()[0, 1] == 3 && wsl::to_string(wsl::jet<int, 3>{1, 2}) == "{1, 2}";
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
