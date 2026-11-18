// wsl::permutation : construction, composition, inverse, cycles, signature,
// inversions, rang lexicographique.

#include <wsl/permutation.hpp>

#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace
{
    template <class F>
    bool throws(F f)
    {
        try { f(); } catch (std::invalid_argument const&) { return true; } catch (std::out_of_range const&) { return true; }
        return false;
    }

    // Inversions par la définition.
    std::uint64_t naive_inversions(wsl::permutation const& p)
    {
        std::uint64_t c = 0;
        for (std::size_t i = 0; i < p.size(); ++i)
            for (std::size_t j = i + 1; j < p.size(); ++j)
                c += p(i) > p(j);
        return c;
    }
}

static_assert(wsl::permutation{2, 0, 1}.inverse() == wsl::permutation{1, 2, 0});
static_assert(wsl::permutation{1, 0, 3, 2}.sign() == 1 && wsl::permutation{1, 2, 0}.cycle_count() == 1);
static_assert(wsl::permutation{3, 1, 2, 0}.inversions() == 5 && wsl::permutation(5).is_identity());
static_assert((wsl::permutation{1, 2, 0} * wsl::permutation{1, 0, 2}) == wsl::permutation{2, 1, 0});   // s(t(i))

int main()
{
    wsl::permutation const s{2, 0, 1}, t = wsl::permutation::from_cycles(3, {{0, 1}});
    assert(t == (wsl::permutation{1, 0, 2}) && (s * t)(0) == s(t(0)) && (s * s.inverse()).is_identity());
    assert(s.power(3).is_identity() && s.power(-1) == s.inverse() && s.power(2) == s * s && s.power(0).is_identity());
    auto const c = wsl::permutation{1, 0, 3, 4, 2, 5}.cycles();
    assert(c.size() == 3 && c[0] == (std::vector<std::size_t>{0, 1}) && c[1] == (std::vector<std::size_t>{2, 3, 4}) && c[2] == (std::vector<std::size_t>{5}));
    wsl::permutation const u{1, 0, 3, 4, 2, 5};
    assert((u.cycle_type() == std::vector<std::size_t>{3, 2, 1}) && u.fixed_points() == 1 && u.sign() == -1);
    // Inversions (tri fusion) et rang lexicographique, sur toutes les permutations de 6 éléments.
    wsl::permutation p(6);
    std::uint64_t rank = 0;
    do
    {
        assert(p.inversions() == naive_inversions(p) && p.lexicographic_rank() == rank);
        assert(wsl::permutation::from_lexicographic_rank(6, rank) == p);
        assert(p.sign() == (naive_inversions(p) % 2 == 0 ? 1 : -1));
        ++rank;
    } while (p.next());
    assert(rank == 720 && p.is_identity());
    // Grande permutation : inversions de l'inverse = inversions ; signature multiplicative.
    std::vector<std::size_t> images(1000);
    for (std::size_t i = 0; i < 1000; ++i)
        images[i] = (i * 389 + 17) % 1000;
    wsl::permutation const q(images);
    assert(q.inverse().inversions() == q.inversions());
    assert((q * q.inverse()).is_identity() && (q * q).sign() == 1);
    std::vector<int> const v{10, 20, 30};
    assert((wsl::permute(s, v) == std::vector<int>{20, 30, 10}));        // r[s(i)] = v[i]
    assert(throws([] { wsl::permutation x{0, 0, 1}; }) && throws([] { wsl::permutation x{0, 3}; }));
    assert(throws([] { (void)(wsl::permutation(3) * wsl::permutation(4)); }) && throws([] { (void)wsl::permutation::from_cycles(3, {{0, 1}, {1, 2}}); }));
    assert(throws([] { (void)wsl::permutation(21).lexicographic_rank(); }) && throws([] { (void)wsl::permutation(3).at(3); }));
    std::printf("test_permutation : ok\n");
}
