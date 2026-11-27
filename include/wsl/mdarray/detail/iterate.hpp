#pragma once

#include <wsl/modes.hpp>

// Parcours interne de toutes les positions d'un ensemble d'indices, en ordre
// lexicographique (ligne majeure), dernier indice dans la boucle intérieure.

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace wsl::inline WSL_DEFAULTS::detail::mdarray
{
    // f(i₀, i₁, …) pour chaque position de e ; s'arrête et rend faux dès que f
    // rend faux (f peut aussi rendre void).
    template <class Extents, class F>
    constexpr bool for_each_index(Extents const& e, F&& f)
    {
        using I = typename Extents::index_type;
        constexpr std::size_t R = Extents::rank();
        if constexpr (R == 0)
        {
            if constexpr (std::is_void_v<decltype(f())>)
            {
                f();
                return true;
            }
            else
                return static_cast<bool>(f());
        }
        else
        {
            for (std::size_t r = 0; r < R; ++r)
                if (e.extent(r) == 0)
                    return true;
            std::array<I, R> i{};
            while (true)
            {
                bool const keep = [&]<std::size_t... K>(std::index_sequence<K...>) {
                    if constexpr (std::is_void_v<decltype(f(i[K]...))>)
                    {
                        f(i[K]...);
                        return true;
                    }
                    else
                        return static_cast<bool>(f(i[K]...));
                }(std::make_index_sequence<R>{});
                if (!keep)
                    return false;
                std::size_t r = R;
                while (r-- > 0)
                {
                    if (++i[r] < e.extent(r))
                        break;
                    i[r] = 0;
                    if (r == 0)
                        return true;
                }
            }
        }
    }
}
