

#include <cassert>
#include <compare>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>

#include <wsl/dlist.hpp>

int main()
{
    wsl::dlist<std::string> empty;
    assert(empty.empty());
    assert(empty.extent() == 0);

    wsl::dlist<std::string> single{"single"};
    assert(single.size() == 1);
    assert(*single.begin() == "single");
    assert(single.extent() == 0);

    wsl::dlist<std::string> tail_only(wsl::dlist<std::string>::with_tail, 10);
    assert(tail_only.empty());
    assert(tail_only.extent() == 10);

    wsl::dlist<std::string> values(wsl::dlist<std::string>::with_tail, 8);
    values.emplace(2, "left");
    values.emplace(5, "right");

    wsl::dlist<std::string> inserted(wsl::dlist<std::string>::with_tail, 4);
    inserted.emplace(1, "middle");
    inserted.emplace(3, "same fragment");

    auto inserted_node = values.insert(3, inserted);
    assert(*inserted_node == "middle");
    assert(values.extent() == 12);

    auto extracted = values.extract(3, 4);
    assert(extracted.size() == 2);
    assert(values.size() == 2);
    assert(values.extent() == 8);

    wsl::dlist<std::string> suffix(wsl::dlist<std::string>::with_tail, 2);
    suffix.emplace(1, "end");
    auto suffix_node = values.insert(5, suffix);
    assert(*suffix_node == "end");
    assert(values.extent() == 10);

    wsl::dlist<std::unique_ptr<int>> movable(
        wsl::dlist<std::unique_ptr<int>>::with_tail, 1);
    movable.emplace(0, std::make_unique<int>(42));
    auto movable_part = movable.extract(0, 1);
    assert(**movable_part.begin() == 42);

    std::ostringstream output;
    output << values;
    assert(output.str() ==
           "dlist{size=3, extent=10}\n"
           "\xE2\x95\x9E\xE2\x95\x90 head  2\n"
           "\xE2\x94\x9C\xE2\x94\x80 note  \"left\"\n"
           "\xE2\x95\x9F\xE2\x94\x80 gap   4\n"
           "\xE2\x94\x9C\xE2\x94\x80 note  \"end\"\n"
           "\xE2\x95\x9F\xE2\x94\x80 gap   1\n"
           "\xE2\x94\x9C\xE2\x94\x80 note  \"right\"\n"
           "\xE2\x95\x98\xE2\x95\x90 tail  3");

    std::ostringstream empty_output;
    empty_output << empty;
    assert(empty_output.str() == "dlist{size=0, extent=0}");

    return 0;
}