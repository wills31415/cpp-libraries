#include <iostream>
#include <string>

#include <wsl/dlist.hpp>

int main()
{
    std::cout << "== Constructions ==\n";

    wsl::dlist<std::string> empty;
    std::cout << empty << "\n\n";

    wsl::dlist<std::string> single{"une note"};
    std::cout << single << "\n\n";

    wsl::dlist<std::string> tail_only(wsl::dlist<std::string>::with_tail, 10);
    std::cout << tail_only << "\n\n";

    std::cout << "== Insertion de notes ==\n";
    wsl::dlist<std::string> document(
        wsl::dlist<std::string>::with_tail, 12);
    document.emplace(2, "premiere note");
    document.emplace(7, "derniere note");
    std::cout << document << "\n\n";

    std::cout << "== Insertion d'un fragment ==\n";
    wsl::dlist<std::string> fragment(
        wsl::dlist<std::string>::with_tail, 4);
    fragment.emplace(1, "note inseree");
    fragment.emplace(2, "seconde note");
    document.insert(5, fragment);
    std::cout << document << "\n\n";

    std::cout << "== Extraction d'un fragment ==\n";
    auto extracted = document.extract(4, 4);
    std::cout << "document:\n"
              << document << "\n\n";
    std::cout << "extrait:\n"
              << extracted << '\n';

    return 0;
}
