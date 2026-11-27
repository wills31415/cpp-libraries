#pragma once

#include <wsl/modes.hpp>

// Axes (design/mdarray.md, L1) : un nom pour toutes les positions d'une dimension.
//
//   auto i = wsl::axis();                       // un axe
//   auto mu = wsl::axis<"mu">();                // étiqueté (messages, affichage)
//   auto [mu, nu] = wsl::axis<"mu", "nu">();    // plusieurs, étiquetés
//   auto [i, j] = wsl::axis<2>();               // plusieurs
//   using namespace wsl::axes;                  // i, j, k, l, m, n, p, q, r, s prédéfinis
//
// L'identité d'un axe est son type : chaque appel de wsl::axis crée des axes
// distincts (lambda en argument par défaut d'un modèle de fonction). Un axe ne se
// copie pas (deux noms ne désignent jamais le même axe) ; `auto i = axis(), j =
// axis();` est refusé par la règle de `auto` (types différents). Au niveau d'un
// espace de noms dans un en-tête, donner une étiquette de type explicite
// (`inline constexpr wsl::mdaxis<tag> i;`) : le lambda n'y est fiable sur aucun
// compilateur.

#include <algorithm>
#include <cstddef>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

namespace wsl::inline WSL_DEFAULTS
{
    namespace detail::mdarray
    {
        // Chaîne en paramètre de modèle (étiquette d'un axe).
        template <std::size_t N>
        struct fixed_string
        {
            char value[N]{};

            constexpr fixed_string(char const (&s)[N]) noexcept { std::copy_n(s, N, value); }
            constexpr std::string_view view() const noexcept { return {value, N - 1}; }
        };

        template <fixed_string Label, class Id>
        struct labelled
        {
            static constexpr std::string_view name = Label.view();
        };

        template <class Id, std::size_t K>
        struct numbered
        {
        };

        template <fixed_string... L>
        constexpr bool distinct_labels() noexcept
        {
            std::string_view const v[] = {std::string_view{}, L.view()...};
            for (std::size_t a = 1; a < sizeof...(L) + 1; ++a)
                for (std::size_t b = a + 1; b < sizeof...(L) + 1; ++b)
                    if (v[a] == v[b])
                        return false;
            return true;
        }
    }

    // Un axe : un objet sans valeur, dont l'identité est le type.
    template <class Id>
    class mdaxis
    {
    public:
        using id_type = Id;

        constexpr mdaxis() noexcept = default;
        mdaxis(mdaxis const&) = delete;
        mdaxis& operator=(mdaxis const&) = delete;

        // Étiquette, vide si l'axe n'en a pas.
        static constexpr std::string_view name() noexcept
        {
            if constexpr (requires { Id::name; })
                return Id::name;
            else
                return {};
        }
    };

    namespace detail::mdarray
    {
        template <class T>
        inline constexpr bool is_axis = false;
        template <class Id>
        inline constexpr bool is_axis<wsl::mdaxis<Id>> = true;

        template <class T>
        concept axis_like = is_axis<std::remove_cvref_t<T>>;

        template <class Id, std::size_t... K>
        constexpr auto numbered_axes(std::index_sequence<K...>) noexcept
        {
            return std::tuple<wsl::mdaxis<numbered<Id, K>>...>{};
        }
    }

    // wsl::axis(), wsl::axis<"mu">(), wsl::axis<"mu", "nu">() : un axe distinct à
    // chaque appel ; plusieurs étiquettes : un groupe, à décomposer.
    template <detail::mdarray::fixed_string... Labels, class Id = decltype([] {})>
    constexpr auto axis() noexcept
    {
        static_assert(detail::mdarray::distinct_labels<Labels...>(), "wsl::axis : étiquettes en double");
        if constexpr (sizeof...(Labels) == 0)
            return mdaxis<Id>{};
        else if constexpr (sizeof...(Labels) == 1)
            return mdaxis<detail::mdarray::labelled<Labels..., Id>>{};
        else
            return std::tuple<mdaxis<detail::mdarray::labelled<Labels, Id>>...>{};
    }

    // wsl::axis<N>() : un groupe de N axes distincts, à décomposer.
    template <std::size_t N, class Id = decltype([] {})>
    constexpr auto axis() noexcept
    {
        static_assert(N > 0, "wsl::axis<N> : au moins un axe");
        return detail::mdarray::numbered_axes<Id>(std::make_index_sequence<N>{});
    }

    // Étiquette d'un axe (vide sans étiquette).
    template <class Id>
    constexpr std::string_view name_of(mdaxis<Id> const&) noexcept
    {
        return mdaxis<Id>::name();
    }

    // Axes prédéfinis, étiquetés : using namespace wsl::axes;
    namespace axes
    {
        namespace tags
        {
#define WSL_MDARRAY_AXIS_TAG(x)                              \
    struct x                                                 \
    {                                                        \
        static constexpr std::string_view name = #x;         \
    };
            WSL_MDARRAY_AXIS_TAG(i)
            WSL_MDARRAY_AXIS_TAG(j)
            WSL_MDARRAY_AXIS_TAG(k)
            WSL_MDARRAY_AXIS_TAG(l)
            WSL_MDARRAY_AXIS_TAG(m)
            WSL_MDARRAY_AXIS_TAG(n)
            WSL_MDARRAY_AXIS_TAG(p)
            WSL_MDARRAY_AXIS_TAG(q)
            WSL_MDARRAY_AXIS_TAG(r)
            WSL_MDARRAY_AXIS_TAG(s)
#undef WSL_MDARRAY_AXIS_TAG
        }

        inline constexpr mdaxis<tags::i> i;
        inline constexpr mdaxis<tags::j> j;
        inline constexpr mdaxis<tags::k> k;
        inline constexpr mdaxis<tags::l> l;
        inline constexpr mdaxis<tags::m> m;
        inline constexpr mdaxis<tags::n> n;
        inline constexpr mdaxis<tags::p> p;
        inline constexpr mdaxis<tags::q> q;
        inline constexpr mdaxis<tags::r> r;
        inline constexpr mdaxis<tags::s> s;
    }
}
