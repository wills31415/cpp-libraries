

#include <algorithm>
#include <cassert>
#include <compare>
#include <iterator>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>
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

    // Positions publiques, effacement par itérateur (un nœud, une plage) : les autres
    // positions et l'étendue sont inchangées.
    {
        wsl::dlist<int> d(wsl::dlist<int>::with_tail, 20);                  // étendue 20, sans valeur
        for (int k : {3, 7, 8, 15, 20})
            d.insert(static_cast<std::size_t>(k), k * 10);
        auto positions = [&] {
            std::vector<std::size_t> r;
            for (auto it = d.cbegin(); it != d.cend(); ++it)
                r.push_back(d.position(it));
            return r;
        };
        assert((positions() == std::vector<std::size_t>{3, 7, 8, 15, 20}) && d.extent() == 20);
        auto next = d.erase(std::next(d.cbegin()));                          // retire 7
        assert(*next == 80 && (positions() == std::vector<std::size_t>{3, 8, 15, 20}) && d.extent() == 20);
        next = d.erase(std::next(d.cbegin()), std::prev(d.cend()));         // retire 8 et 15
        assert(*next == 200 && (positions() == std::vector<std::size_t>{3, 20}) && d.extent() == 20);
        next = d.erase(d.cbegin(), d.cend());
        assert(next == d.end() && d.empty() && d.extent() == 20);
    }
    // extract sur une longue liste, contre un modèle (positions absolues) : la plage est
    // retirée, la suite recule de sa longueur.
    {
        wsl::dlist<int> d(wsl::dlist<int>::with_tail, 7000);
        std::vector<std::pair<std::size_t, int>> model;
        for (int k = 0; k < 2000; ++k)
        {
            std::size_t const pos = static_cast<std::size_t>(3 * k + (k % 5));
            d.insert(pos, k);
            model.emplace_back(pos, k);
        }
        std::ranges::sort(model);                                           // la liste range par position
        std::size_t const from = 1000, length = 2500;
        auto piece = d.extract(from, length);
        std::vector<std::pair<std::size_t, int>> kept, cut;
        for (auto [pos, v] : model)
            if (pos < from)
                kept.emplace_back(pos, v);
            else if (pos < from + length)
                cut.emplace_back(pos - from, v);
            else
                kept.emplace_back(pos - length, v);
        auto as_pairs = [](wsl::dlist<int> const& l) {
            std::vector<std::pair<std::size_t, int>> r;
            for (auto it = l.cbegin(); it != l.cend(); ++it)
                r.emplace_back(l.position(it), *it);
            return r;
        };
        assert(as_pairs(d) == kept && as_pairs(piece) == cut && piece.extent() == length);
    }

    return 0;
}