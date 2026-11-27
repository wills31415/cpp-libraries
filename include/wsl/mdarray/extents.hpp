#pragma once

#include <wsl/modes.hpp>

// Étendues d'un tableau multidimensionnel : l'ensemble d'indices
// E = [0, n₁) × ⋯ × [0, n_N) (design/mdarray.md, A5, I2).
//
//   wsl::extents<std::size_t, 3, std::dynamic_extent> e(4);   // 3 × 4, la seconde dynamique
//   wsl::dims<2> f(3, 4);                                     // toutes dynamiques (C++26)
//   e.contains(2, 3);  e.contains(p);  e.contains(p, d);      // ajouts
//
// Noms et sémantique de std::extents (C++23), ajouts seulement : contains, et
// wrap, clamp (positions ramenées dans E, D22, D23).

#include <wsl/mdarray/detail/precondition.hpp>
#include <wsl/mdarray/integer.hpp>
#include <wsl/mdarray/position.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>

namespace wsl::inline WSL_DEFAULTS
{
    using std::dynamic_extent;

    template <class IndexType, std::size_t... Extents>
    class extents;

    namespace detail::mdarray
    {
        template <class T>
        inline constexpr bool is_extents = false;
        template <class I, std::size_t... E>
        inline constexpr bool is_extents<wsl::extents<I, E...>> = true;

        template <class T>
        concept extents_like = is_extents<std::remove_cvref_t<T>>;

        // Nombre d'étendues dynamiques parmi les r premières.
        template <std::size_t... E>
        constexpr std::size_t dynamic_before(std::size_t r) noexcept
        {
            constexpr std::array<std::size_t, sizeof...(E)> values{E...};
            std::size_t count = 0;
            for (std::size_t k = 0; k < r; ++k)
                count += values[k] == dynamic_extent;
            return count;
        }

        // Étendues dynamiques ; vide (sans octet) quand il n'y en a aucune, pour que
        // des étendues fixes n'occupent aucune place.
        template <class I, std::size_t N>
        struct dynamic_storage
        {
            std::array<I, N> values{};
            constexpr I operator[](std::size_t k) const noexcept { return values[k]; }
            constexpr void set(std::size_t k, I v) noexcept { values[k] = v; }
        };
        template <class I>
        struct dynamic_storage<I, 0>
        {
            constexpr I operator[](std::size_t) const noexcept { return I{}; }
            constexpr void set(std::size_t, I) noexcept {}
        };

        template <class To, class From>
        constexpr To index_cast(From value, char const* message) noexcept
        {
            precondition(std::cmp_greater_equal(value, 0) && std::in_range<To>(value), message);
            return static_cast<To>(value);
        }
    }

    template <class IndexType, std::size_t... Extents>
    class extents
    {
        static_assert(detail::mdarray::index_integer<IndexType>,
                      "wsl::extents : IndexType doit être un type entier (ni bool ni caractère)");
        static_assert(((Extents == dynamic_extent || std::in_range<IndexType>(Extents)) && ...),
                      "wsl::extents : étendue fixe non représentable dans IndexType");

    public:
        using index_type = IndexType;
        using size_type = std::make_unsigned_t<IndexType>;
        using rank_type = std::size_t;

        static constexpr rank_type rank() noexcept { return sizeof...(Extents); }
        static constexpr rank_type rank_dynamic() noexcept { return ((Extents == dynamic_extent) + ... + 0); }

        static constexpr std::size_t static_extent(rank_type r) noexcept
        {
            constexpr std::array<std::size_t, sizeof...(Extents)> values{Extents...};
            detail::mdarray::precondition(r < rank(), "extents::static_extent : rang hors limites");
            return values[r];
        }

        constexpr index_type extent(rank_type r) const noexcept
        {
            detail::mdarray::precondition(r < rank(), "extents::extent : rang hors limites");
            if (static_extent(r) == dynamic_extent)
                return dynamic_[dynamic_index(r)];
            return static_cast<index_type>(static_extent(r));
        }

        constexpr extents() noexcept = default;

        // Depuis d'autres étendues de même rang, compatibles.
        template <class OtherIndexType, std::size_t... OtherExtents>
            requires(sizeof...(OtherExtents) == rank()
                     && ((OtherExtents == dynamic_extent || Extents == dynamic_extent || OtherExtents == Extents) && ...))
        constexpr explicit((((Extents != dynamic_extent) && (OtherExtents == dynamic_extent)) || ...)
                           || std::cmp_less(std::numeric_limits<index_type>::max(), std::numeric_limits<OtherIndexType>::max()))
            extents(extents<OtherIndexType, OtherExtents...> const& other) noexcept
        {
            for (rank_type r = 0; r < rank(); ++r)
            {
                auto const v = other.extent(r);
                if (static_extent(r) == dynamic_extent)
                    dynamic_.set(dynamic_index(r), detail::mdarray::index_cast<index_type>(v, "extents : étendue non représentable"));
                else
                    detail::mdarray::precondition(std::cmp_equal(v, static_extent(r)), "extents : étendue fixe différente");
            }
        }

        // Les étendues dynamiques seules, ou toutes.
        template <class... OtherIndexTypes>
            requires((std::is_convertible_v<OtherIndexTypes, index_type> && ...)
                     && (std::is_nothrow_constructible_v<index_type, OtherIndexTypes> && ...)
                     && (sizeof...(OtherIndexTypes) == rank_dynamic() || sizeof...(OtherIndexTypes) == rank()))
        constexpr explicit extents(OtherIndexTypes... values) noexcept
        {
            std::array<index_type, sizeof...(OtherIndexTypes)> const v{
                detail::mdarray::index_cast<index_type>(values, "extents : étendue négative ou non représentable")...};
            assign(std::span<index_type const, sizeof...(OtherIndexTypes)>(v));
        }

        template <class OtherIndexType, std::size_t N>
            requires(std::is_convertible_v<OtherIndexType const&, index_type>
                     && std::is_nothrow_constructible_v<index_type, OtherIndexType const&>
                     && (N == rank_dynamic() || N == rank()))
        constexpr explicit(N != rank_dynamic()) extents(std::span<OtherIndexType, N> values) noexcept
        {
            std::array<index_type, N> v{};
            for (std::size_t k = 0; k < N; ++k)
                v[k] = detail::mdarray::index_cast<index_type>(values[k], "extents : étendue négative ou non représentable");
            assign(std::span<index_type const, N>(v));
        }

        template <class OtherIndexType, std::size_t N>
            requires(std::is_convertible_v<OtherIndexType const&, index_type>
                     && std::is_nothrow_constructible_v<index_type, OtherIndexType const&>
                     && (N == rank_dynamic() || N == rank()))
        constexpr explicit(N != rank_dynamic()) extents(std::array<OtherIndexType, N> const& values) noexcept
            : extents(std::span<OtherIndexType const, N>(values))
        {
        }

        template <class OtherIndexType, std::size_t... OtherExtents>
        friend constexpr bool operator==(extents const& a, extents<OtherIndexType, OtherExtents...> const& b) noexcept
        {
            if constexpr (rank() != sizeof...(OtherExtents))
                return false;
            else
            {
                for (rank_type r = 0; r < rank(); ++r)
                    if (!std::cmp_equal(a.extent(r), b.extent(r)))
                        return false;
                return true;
            }
        }

        // --- Ajouts (D19) : appartenance à E, toujours vérifiée ----------------

        template <detail::mdarray::index_integer... I>
            requires(sizeof...(I) == rank())
        constexpr bool contains(I... indices) const noexcept
        {
            [[maybe_unused]] rank_type r = 0;
            return ((std::cmp_greater_equal(indices, 0) && std::cmp_less(indices, extent(r++))) && ...);
        }

        template <class I>
        constexpr bool contains(mdposition<sizeof...(Extents), I> const& p) const noexcept
        {
            for (rank_type r = 0; r < rank(); ++r)
                if (!std::cmp_less(p[r], extent(r)))
                    return false;
            return true;
        }

        template <detail::mdarray::index_integer J, std::size_t N>
            requires(N == sizeof...(Extents))
        constexpr bool contains(std::array<J, N> const& p) const noexcept
        {
            for (rank_type r = 0; r < rank(); ++r)
                if (!(std::cmp_greater_equal(p[r], 0) && std::cmp_less(p[r], extent(r))))
                    return false;
            return true;
        }

        // p + d ∈ E, sans former p + d (qui peut être négatif).
        template <class I, class D>
        constexpr bool contains(mdposition<sizeof...(Extents), I> const& p,
                                mddisplacement<sizeof...(Extents), D> const& d) const noexcept
        {
            for (rank_type r = 0; r < rank(); ++r)
            {
                auto const q = p[r];
                auto const s = d[r];
                if (s < 0)
                {
                    if (std::cmp_less(q, detail::mdarray::magnitude(s)))
                        return false;
                    if (!std::cmp_less(q - static_cast<I>(detail::mdarray::magnitude(s)), extent(r)))
                        return false;
                }
                else
                {
                    // q + s < n  ⇔  s < n − q (n − q > 0 quand q < n).
                    if (!std::cmp_less(q, extent(r)) || !std::cmp_less(s, extent(r) - static_cast<index_type>(q)))
                        return false;
                }
            }
            return true;
        }

    private:
        static constexpr rank_type dynamic_index(rank_type r) noexcept
        {
            return detail::mdarray::dynamic_before<Extents...>(r);
        }

        // N == rank_dynamic() : les étendues dynamiques ; N == rank() : toutes, les
        // fixes devant coïncider (précondition).
        template <std::size_t N>
        constexpr void assign(std::span<index_type const, N> v) noexcept
        {
            if constexpr (N == rank_dynamic())
            {
                for (std::size_t k = 0; k < N; ++k)
                    dynamic_.set(k, v[k]);
            }
            else
            {
                for (rank_type r = 0; r < rank(); ++r)
                {
                    if (static_extent(r) == dynamic_extent)
                        dynamic_.set(dynamic_index(r), v[r]);
                    else
                        detail::mdarray::precondition(std::cmp_equal(v[r], static_extent(r)), "extents : étendue fixe différente");
                }
            }
        }

        [[no_unique_address]] detail::mdarray::dynamic_storage<index_type, rank_dynamic()> dynamic_{};
    };

    namespace detail::mdarray
    {
        template <class IndexType, class Seq>
        struct dextents_of;
        template <class IndexType, std::size_t... K>
        struct dextents_of<IndexType, std::index_sequence<K...>>
        {
            using type = wsl::extents<IndexType, ((void)K, dynamic_extent)...>;
        };
    }

    template <class IndexType, std::size_t Rank>
    using dextents = typename detail::mdarray::dextents_of<IndexType, std::make_index_sequence<Rank>>::type;

    // C++26 (P2389).
    template <std::size_t Rank, class IndexType = std::size_t>
    using dims = dextents<IndexType, Rank>;

    template <class... Integrals>
        requires(std::is_convertible_v<Integrals, std::size_t> && ...)
    explicit extents(Integrals...) -> extents<std::size_t, ((void)sizeof(Integrals), dynamic_extent)...>;

    // --- Ajouts (D22, D23) : positions ramenées dans E -------------------------

    // Bord périodique : p + d modulo E (tore), précondition : étendues non nulles.
    template <class I, class D, class IndexType, std::size_t... Extents>
    constexpr mdposition<sizeof...(Extents), I> wrap(mdposition<sizeof...(Extents), I> const& p,
                                                     mddisplacement<sizeof...(Extents), D> const& d,
                                                     extents<IndexType, Extents...> const& e) noexcept
    {
        std::array<I, sizeof...(Extents)> r{};
        for (std::size_t k = 0; k < sizeof...(Extents); ++k)
        {
            detail::mdarray::precondition(e.extent(k) > 0, "wrap : étendue nulle");
            // (p + d) mod n, calculé sans débordement : p mod n + d mod n, dans [0, 2n).
            using U = std::make_unsigned_t<std::common_type_t<std::make_unsigned_t<I>, std::make_unsigned_t<IndexType>>>;
            U const n = static_cast<U>(e.extent(k));
            U const a = static_cast<U>(p[k]) % n;
            U const b = d[k] >= 0 ? static_cast<U>(detail::mdarray::magnitude(d[k])) % n
                                  : static_cast<U>((n - static_cast<U>(detail::mdarray::magnitude(d[k])) % n) % n);
            U const s = a >= n - b ? a - (n - b) : a + b;
            r[k] = static_cast<I>(s);
        }
        return mdposition<sizeof...(Extents), I>(r);
    }

    template <class I, class IndexType, std::size_t... Extents>
    constexpr mdposition<sizeof...(Extents), I> wrap(mdposition<sizeof...(Extents), I> const& p,
                                                     extents<IndexType, Extents...> const& e) noexcept
    {
        return wrap(p, mddisplacement<sizeof...(Extents)>{}, e);
    }

    // Bord recopié : p + d ramené au bord le plus proche, précondition : étendues non nulles.
    template <class I, class D, class IndexType, std::size_t... Extents>
    constexpr mdposition<sizeof...(Extents), I> clamp(mdposition<sizeof...(Extents), I> const& p,
                                                      mddisplacement<sizeof...(Extents), D> const& d,
                                                      extents<IndexType, Extents...> const& e) noexcept
    {
        std::array<I, sizeof...(Extents)> r{};
        for (std::size_t k = 0; k < sizeof...(Extents); ++k)
        {
            detail::mdarray::precondition(e.extent(k) > 0, "clamp : étendue nulle");
            auto const last = e.extent(k) - 1;
            auto const q = p[k];
            auto const s = d[k];
            if (s < 0)
            {
                auto const m = detail::mdarray::magnitude(s);
                r[k] = std::cmp_less_equal(q, m) ? I{0}
                                                 : static_cast<I>(std::cmp_less(q - static_cast<I>(m), last) ? q - static_cast<I>(m)
                                                                                                             : static_cast<I>(last));
            }
            else
            {
                // q + s ≥ last  ⇔  s ≥ last − q (si q ≤ last).
                bool const beyond = std::cmp_greater_equal(q, last)
                                    || std::cmp_greater_equal(s, last - static_cast<IndexType>(q));
                r[k] = beyond ? static_cast<I>(last) : static_cast<I>(q + static_cast<I>(s));
            }
        }
        return mdposition<sizeof...(Extents), I>(r);
    }

    template <class I, class IndexType, std::size_t... Extents>
    constexpr mdposition<sizeof...(Extents), I> clamp(mdposition<sizeof...(Extents), I> const& p,
                                                      extents<IndexType, Extents...> const& e) noexcept
    {
        return clamp(p, mddisplacement<sizeof...(Extents)>{}, e);
    }
}
