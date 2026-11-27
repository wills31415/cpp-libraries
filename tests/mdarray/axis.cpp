// wsl::mdarray, axes (design/mdarray.md, L1) : identité par le type, étiquettes,
// groupes, axes prédéfinis, absence de copie. Les refus qui sont des erreurs de
// compilation sont dans tests/mdarray/compile_fail/axis.cpp.

#include <wsl/mdarray.hpp>

#include <cstdio>
#include <string_view>
#include <type_traits>

namespace
{
    template <class A, class B>
    constexpr bool same(A const&, B const&)
    {
        return std::is_same_v<A, B>;
    }

    // Une fonction constexpr, un modèle : les axes y sont distincts (évaluation constante).
    template <class T>
    constexpr bool distinct_in_template(T)
    {
        auto i = wsl::axis();
        auto j = wsl::axis();
        auto [p, q] = wsl::axis<2>();
        auto [mu, nu] = wsl::axis<"mu", "nu">();
        auto mu2 = wsl::axis<"mu">();
        return !same(i, j) && !same(p, q) && !same(i, p) && !same(mu, nu) && !same(mu, mu2) && mu.name() == "mu"
               && name_of(nu) == "nu" && name_of(i).empty();
    }
    static_assert(distinct_in_template(1) && distinct_in_template(1.0));

    // Deux appels identiques donnent des groupes distincts.
    constexpr bool groups()
    {
        auto [a, b] = wsl::axis<2>();
        auto [c, d] = wsl::axis<2>();
        auto [x] = wsl::axis<1>();
        return !same(a, c) && !same(b, d) && !same(a, x);
    }
    static_assert(groups());

    // Axes prédéfinis, étiquetés.
    using namespace wsl::axes;
    static_assert(!same(i, j) && i.name() == "i" && name_of(s) == "s");

    // Un axe ne se copie pas (deux noms ne désignent jamais le même axe).
    using A = decltype(wsl::axis());
    static_assert(!std::is_copy_constructible_v<A> && !std::is_copy_assignable_v<A>);
    static_assert(!std::is_move_constructible_v<A>);
    static_assert(std::is_empty_v<A>);
}

int main()
{
    auto [r1, r2] = wsl::axis<"r1", "r2">();
    std::printf("mdarray.axis : OK (%.*s, %.*s)\n", int(r1.name().size()), r1.name().data(), int(r2.name().size()),
                r2.name().data());
}
