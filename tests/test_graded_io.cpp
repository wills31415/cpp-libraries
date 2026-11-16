#include <wsl/graded_io.hpp>

#include <cassert>
#include <sstream>

using namespace wsl;
int main()
{
    // forme compacte : les equivalents s'affichent identiquement
    assert(std::format("{}", series<int, 2>{1, 2}) == "{1, 2}");
    assert(std::format("{}", series<int, 9>{1, 2, 0, 0}) == "{1, 2}");
    assert(std::format("{}", jet<int, 4>{1, 2}) == "{1, 2}");

    // les deux arites s'affichent pareil : la distinction est nominale
    assert(std::format("{}", jet<int, 5>{0, 3, 7}) ==
           std::format("{}", series<int, 5>{0, 3, 7}));
    assert(std::format("{}", index_jet<3>{4, 5}) == "{4, 5}");

    // les zeros de tete restent
    assert(std::format("{}", jet<int, 4>{0, 0, 3}) == "{0, 0, 3}");

    // neutre
    assert(std::format("{}", jet<int, 4>{}) == "{}");
    assert(std::format("{}", jet<int, 0>{}) == "{}");

    // signes et types
    assert(std::format("{}", jet<int, 3>{-1, 0, -7}) == "{-1, 0, -7}");
    assert(std::format("{}", jet<unsigned, 2>{4294967295u, 1}) == "{4294967295, 1}");

    // forme complete : taille pour le fixe, plafond pour le flexible
    assert(std::format("{:#}", jet<int, 4>{1, 2}) == "{1, 2, 0, 0}");
    assert(std::format("{:#}", series<int, 4>{1, 2}) == "{1, 2, 0, 0}");
    assert(std::format("{:#}", jet<int, 4>{}) == "{0, 0, 0, 0}");
    assert(std::format("{:#}", series<int, 2>{1, 2}) == "{1, 2}");

    // flux, trouve par ADL
    std::ostringstream stream;
    stream << jet<int, 4>{1, 2, 0, 0} << ' ' << series<int, 3>{0, 5};
    assert(stream.str() == "{1, 2} {0, 5}");

    assert(to_string(jet<int, 3>{7}) == "{7}");

    // imbrication dans un format plus large
    assert(std::format("de {} a {}", jet<int, 2>{}, jet<int, 2>{3, 4}) ==
           "de {} a {3, 4}");

    return 0;
}
