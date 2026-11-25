#pragma once

// Listes de types, pour l'algèbre des axes à la compilation (docs/mdarray.md, I5).

#include <cstddef>
#include <type_traits>

namespace wsl::detail::mdarray
{
    template <class... T>
    struct tlist
    {
        static constexpr std::size_t size = sizeof...(T);
    };

    template <class... L>
    struct tconcat;
    template <>
    struct tconcat<>
    {
        using type = tlist<>;
    };
    template <class... A>
    struct tconcat<tlist<A...>>
    {
        using type = tlist<A...>;
    };
    template <class... A, class... B, class... Rest>
    struct tconcat<tlist<A...>, tlist<B...>, Rest...>
    {
        using type = typename tconcat<tlist<A..., B...>, Rest...>::type;
    };
    template <class... L>
    using tconcat_t = typename tconcat<L...>::type;

    template <class X, class L>
    inline constexpr std::size_t tcount = 0;
    template <class X, class... T>
    inline constexpr std::size_t tcount<X, tlist<T...>> = (std::size_t{std::is_same_v<X, T>} + ... + 0);

    template <class X, class L>
    inline constexpr bool tcontains = tcount<X, L> > 0;

    // Sans répétition, dans l'ordre de première apparition.
    template <class Done, class Todo>
    struct tunique_impl;
    template <class... D>
    struct tunique_impl<tlist<D...>, tlist<>>
    {
        using type = tlist<D...>;
    };
    template <class... D, class X, class... Rest>
    struct tunique_impl<tlist<D...>, tlist<X, Rest...>>
    {
        using type = typename tunique_impl<std::conditional_t<tcontains<X, tlist<D...>>, tlist<D...>, tlist<D..., X>>,
                                           tlist<Rest...>>::type;
    };
    template <class L>
    using tunique_t = typename tunique_impl<tlist<>, L>::type;

    // Éléments de L absents de S (répétitions gardées).
    template <class L, class S>
    struct tminus;
    template <class S>
    struct tminus<tlist<>, S>
    {
        using type = tlist<>;
    };
    template <class X, class... Rest, class S>
    struct tminus<tlist<X, Rest...>, S>
    {
        using rest = typename tminus<tlist<Rest...>, S>::type;
        using type = std::conditional_t<tcontains<X, S>, rest, tconcat_t<tlist<X>, rest>>;
    };
    template <class L, class S>
    using tminus_t = typename tminus<L, S>::type;

    // Position de X dans une liste sans répétition.
    template <class X, class L>
    struct tindex;
    template <class X, class... T>
    struct tindex<X, tlist<X, T...>>
    {
        static constexpr std::size_t value = 0;
    };
    template <class X, class Y, class... T>
    struct tindex<X, tlist<Y, T...>>
    {
        static constexpr std::size_t value = 1 + tindex<X, tlist<T...>>::value;
    };
    template <class X, class L>
    inline constexpr std::size_t tindex_v = tindex<X, L>::value;

    // Vrai si f<X>() pour tout X de L.
    template <class L, template <class> class P>
    inline constexpr bool tall = false;
    template <class... T, template <class> class P>
    inline constexpr bool tall<tlist<T...>, P> = (P<T>::value && ...);
}
