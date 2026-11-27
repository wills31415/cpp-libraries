#pragma once

// Vue multidimensionnelle (design/mdarray.md, A5) : noms et sémantique de
// std::mdspan (C++23), ajouts seulement :
//   m[p] avec p un wsl::mdposition (L4) ; m.at(…) toujours vérifié,
//   std::out_of_range (D20, comme std::mdspan::at de C++26) ; m.contains(…) (D19) ;
//   conversions avec std::mdspan là où il existe (__cpp_lib_mdspan).

#include <wsl/mdarray/axis.hpp>
#include <wsl/mdarray/detail/expression_base.hpp>
#include <wsl/mdarray/detail/precondition.hpp>
#include <wsl/mdarray/extents.hpp>
#include <wsl/mdarray/layout.hpp>
#include <wsl/mdarray/position.hpp>

#include <array>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>

#if __has_include(<mdspan>)
#    include <mdspan>
#endif

namespace wsl
{
    template <class ElementType>
    struct default_accessor
    {
        using offset_policy = default_accessor;
        using element_type = ElementType;
        using reference = ElementType&;
        using data_handle_type = ElementType*;

        constexpr default_accessor() noexcept = default;

        template <class OtherElementType>
            requires std::is_convertible_v<OtherElementType (*)[], element_type (*)[]>
        constexpr default_accessor(default_accessor<OtherElementType>) noexcept
        {
        }

        constexpr reference access(data_handle_type p, std::size_t i) const noexcept { return p[i]; }
        constexpr data_handle_type offset(data_handle_type p, std::size_t i) const noexcept { return p + i; }
    };

#if defined(__cpp_lib_mdspan)
    namespace detail::mdarray
    {
        template <class Layout>
        struct std_layout;
        template <>
        struct std_layout<wsl::layout_right>
        {
            using type = std::layout_right;
        };
        template <>
        struct std_layout<wsl::layout_left>
        {
            using type = std::layout_left;
        };
        template <>
        struct std_layout<wsl::layout_stride>
        {
            using type = std::layout_stride;
        };

        template <class Layout>
        struct wsl_layout;
        template <>
        struct wsl_layout<std::layout_right>
        {
            using type = wsl::layout_right;
        };
        template <>
        struct wsl_layout<std::layout_left>
        {
            using type = wsl::layout_left;
        };
        template <>
        struct wsl_layout<std::layout_stride>
        {
            using type = wsl::layout_stride;
        };

        template <class E>
        struct std_extents_of;
        template <class I, std::size_t... E>
        struct std_extents_of<wsl::extents<I, E...>>
        {
            using type = std::extents<I, E...>;
        };

        template <class E>
        struct wsl_extents_of;
        template <class I, std::size_t... E>
        struct wsl_extents_of<std::extents<I, E...>>
        {
            using type = wsl::extents<I, E...>;
        };

        template <class To, class From>
        constexpr To convert_extents(From const& e) noexcept
        {
            std::array<typename To::index_type, From::rank()> v{};
            for (std::size_t r = 0; r < From::rank(); ++r)
                v[r] = e.extent(r);
            return To(v);
        }
    }
#endif

    template <class ElementType, class Extents, class LayoutPolicy = layout_right,
              class AccessorPolicy = default_accessor<ElementType>>
    class mdspan
    {
        static_assert(detail::mdarray::is_extents<Extents>, "wsl::mdspan : Extents doit être un wsl::extents");

    public:
        using extents_type = Extents;
        using layout_type = LayoutPolicy;
        using accessor_type = AccessorPolicy;
        using mapping_type = typename layout_type::template mapping<extents_type>;
        using element_type = ElementType;
        using value_type = std::remove_cv_t<element_type>;
        using index_type = typename extents_type::index_type;
        using size_type = typename extents_type::size_type;
        using rank_type = typename extents_type::rank_type;
        using data_handle_type = typename accessor_type::data_handle_type;
        using reference = typename accessor_type::reference;

        static constexpr rank_type rank() noexcept { return extents_type::rank(); }
        static constexpr rank_type rank_dynamic() noexcept { return extents_type::rank_dynamic(); }
        static constexpr std::size_t static_extent(rank_type r) noexcept { return extents_type::static_extent(r); }
        constexpr index_type extent(rank_type r) const noexcept { return extents().extent(r); }

        constexpr mdspan()
            requires(rank_dynamic() > 0 && std::is_default_constructible_v<data_handle_type>
                     && std::is_default_constructible_v<mapping_type> && std::is_default_constructible_v<accessor_type>)
        = default;
        constexpr mdspan(mdspan const&) = default;
        constexpr mdspan(mdspan&&) = default;

        template <class... OtherIndexTypes>
            requires((std::is_convertible_v<OtherIndexTypes, index_type> && ...)
                     && (std::is_nothrow_constructible_v<index_type, OtherIndexTypes> && ...)
                     && (sizeof...(OtherIndexTypes) == rank() || sizeof...(OtherIndexTypes) == rank_dynamic())
                     && std::is_constructible_v<mapping_type, extents_type> && std::is_default_constructible_v<accessor_type>)
        constexpr explicit mdspan(data_handle_type p, OtherIndexTypes... exts)
            : ptr_(std::move(p)), map_(extents_type(static_cast<index_type>(std::move(exts))...))
        {
        }

        template <class OtherIndexType, std::size_t N>
            requires(std::is_convertible_v<OtherIndexType const&, index_type>
                     && std::is_nothrow_constructible_v<index_type, OtherIndexType const&>
                     && (N == rank() || N == rank_dynamic()) && std::is_constructible_v<mapping_type, extents_type>
                     && std::is_default_constructible_v<accessor_type>)
        constexpr explicit(N != rank_dynamic()) mdspan(data_handle_type p, std::span<OtherIndexType, N> exts)
            : ptr_(std::move(p)), map_(extents_type(exts))
        {
        }

        template <class OtherIndexType, std::size_t N>
            requires(std::is_convertible_v<OtherIndexType const&, index_type>
                     && std::is_nothrow_constructible_v<index_type, OtherIndexType const&>
                     && (N == rank() || N == rank_dynamic()) && std::is_constructible_v<mapping_type, extents_type>
                     && std::is_default_constructible_v<accessor_type>)
        constexpr explicit(N != rank_dynamic()) mdspan(data_handle_type p, std::array<OtherIndexType, N> const& exts)
            : ptr_(std::move(p)), map_(extents_type(exts))
        {
        }

        constexpr mdspan(data_handle_type p, extents_type const& ext)
            requires(std::is_constructible_v<mapping_type, extents_type const&> && std::is_default_constructible_v<accessor_type>)
            : ptr_(std::move(p)), map_(ext)
        {
        }

        constexpr mdspan(data_handle_type p, mapping_type const& m)
            requires std::is_default_constructible_v<accessor_type>
            : ptr_(std::move(p)), map_(m)
        {
        }

        constexpr mdspan(data_handle_type p, mapping_type const& m, accessor_type const& a)
            : ptr_(std::move(p)), map_(m), acc_(a)
        {
        }

        template <class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherAccessor>
            requires(std::is_constructible_v<mapping_type, typename OtherLayoutPolicy::template mapping<OtherExtents> const&>
                     && std::is_constructible_v<accessor_type, OtherAccessor const&>)
        constexpr explicit(!std::is_convertible_v<typename OtherLayoutPolicy::template mapping<OtherExtents> const&, mapping_type>
                           || !std::is_convertible_v<OtherAccessor const&, accessor_type>)
            mdspan(mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherAccessor> const& other)
            : ptr_(other.data_handle()), map_(other.mapping()), acc_(other.accessor())
        {
            static_assert(std::is_constructible_v<data_handle_type, typename OtherAccessor::data_handle_type const&>);
            static_assert(std::is_constructible_v<extents_type, OtherExtents>);
        }

        constexpr mdspan& operator=(mdspan const&) = default;
        constexpr mdspan& operator=(mdspan&&) = default;

        // --- Accès (C++23) -----------------------------------------------------

        template <class... OtherIndexTypes>
            requires((std::is_convertible_v<OtherIndexTypes, index_type> && ...)
                     && (std::is_nothrow_constructible_v<index_type, OtherIndexTypes> && ...) && sizeof...(OtherIndexTypes) == rank())
        constexpr reference operator[](OtherIndexTypes... indices) const
        {
            return acc_.access(ptr_, static_cast<std::size_t>(map_(static_cast<index_type>(std::move(indices))...)));
        }

        template <class OtherIndexType>
            requires(std::is_convertible_v<OtherIndexType const&, index_type>
                     && std::is_nothrow_constructible_v<index_type, OtherIndexType const&>)
        constexpr reference operator[](std::span<OtherIndexType, rank()> indices) const
        {
            return [&]<std::size_t... K>(std::index_sequence<K...>) -> reference {
                return (*this)[static_cast<index_type>(std::as_const(indices[K]))...];
            }(std::make_index_sequence<rank()>{});
        }

        template <class OtherIndexType>
            requires(std::is_convertible_v<OtherIndexType const&, index_type>
                     && std::is_nothrow_constructible_v<index_type, OtherIndexType const&>)
        constexpr reference operator[](std::array<OtherIndexType, rank()> const& indices) const
        {
            return (*this)[std::span<OtherIndexType const, rank()>(indices)];
        }

        // --- Ajouts : axes (L2, L3), position (L4), at (D20), contains (D19) -------

        // m[i, j], m[1, j] : une expression indicée (monde des axes) ; une cible si la
        // vue est modifiable (m[i, j] = …).
        template <class... A>
            requires(sizeof...(A) == rank() && (detail::mdarray::axis_like<A> || ...)
                     && ((detail::mdarray::axis_like<A> || detail::mdarray::index_integer<std::remove_cvref_t<A>>) && ...))
        constexpr auto operator[](A const&... a) const
        {
            return detail::mdarray::make_indexed(*this, a...);
        }

        // u[i + 1] : refusé, i + 1 est une valeur, pas un axe (V2).
        template <class... A>
            requires(detail::mdarray::indexed<A> || ...)
        void operator[](A const&...) const
        {
            static_assert(sizeof...(A) == 0, "wsl::mdarray : i + 1 est une valeur, pas un axe (V2) ; pour décaler, "
                                             "utiliser une vue (roll, shift, slice : lot 2)");
        }

        template <class I>
        constexpr reference operator[](mdposition<rank(), I> const& p) const
        {
            return (*this)[p.to_array()];
        }

        template <detail::mdarray::index_integer... OtherIndexTypes>
            requires(sizeof...(OtherIndexTypes) == rank())
        constexpr reference at(OtherIndexTypes... indices) const
        {
            if (!extents().contains(indices...))
                throw std::out_of_range("wsl::mdspan::at : position hors du tableau");
            return (*this)[static_cast<index_type>(indices)...];
        }

        template <detail::mdarray::index_integer OtherIndexType>
        constexpr reference at(std::span<OtherIndexType, rank()> indices) const
        {
            return [&]<std::size_t... K>(std::index_sequence<K...>) -> reference {
                return at(indices[K]...);
            }(std::make_index_sequence<rank()>{});
        }

        template <detail::mdarray::index_integer OtherIndexType>
        constexpr reference at(std::array<OtherIndexType, rank()> const& indices) const
        {
            return at(std::span<OtherIndexType const, rank()>(indices));
        }

        template <class I>
        constexpr reference at(mdposition<rank(), I> const& p) const
        {
            return at(p.to_array());
        }

        template <class... A>
        constexpr bool contains(A const&... a) const noexcept
            requires requires(extents_type const& e) { e.contains(a...); }
        {
            return extents().contains(a...);
        }

        // --- Observateurs (C++23) ----------------------------------------------

        constexpr size_type size() const noexcept
        {
            size_type s{1};
            for (rank_type r = 0; r < rank(); ++r)
            {
                detail::mdarray::precondition(!__builtin_mul_overflow(s, static_cast<size_type>(extent(r)), &s),
                                              "mdspan::size : taille non représentable");
            }
            return s;
        }

        [[nodiscard]] constexpr bool empty() const noexcept
        {
            for (rank_type r = 0; r < rank(); ++r)
                if (extent(r) == 0)
                    return true;
            return false;
        }

        friend constexpr void swap(mdspan& x, mdspan& y) noexcept
        {
            std::swap(x.ptr_, y.ptr_);
            std::swap(x.map_, y.map_);
            std::swap(x.acc_, y.acc_);
        }

        constexpr extents_type const& extents() const noexcept { return map_.extents(); }
        constexpr data_handle_type const& data_handle() const noexcept { return ptr_; }
        constexpr mapping_type const& mapping() const noexcept { return map_; }
        constexpr accessor_type const& accessor() const noexcept { return acc_; }

        static constexpr bool is_always_unique() { return mapping_type::is_always_unique(); }
        static constexpr bool is_always_exhaustive() { return mapping_type::is_always_exhaustive(); }
        static constexpr bool is_always_strided() { return mapping_type::is_always_strided(); }

        constexpr bool is_unique() const { return map_.is_unique(); }
        constexpr bool is_exhaustive() const { return map_.is_exhaustive(); }
        constexpr bool is_strided() const { return map_.is_strided(); }
        constexpr index_type stride(rank_type r) const { return map_.stride(r); }

#if defined(__cpp_lib_mdspan)
        // --- Conversions avec std::mdspan (A5) -----------------------------------

        template <class StdLayout>
            requires(std::is_same_v<accessor_type, default_accessor<element_type>>
                     && std::is_same_v<layout_type, typename detail::mdarray::wsl_layout<StdLayout>::type>)
        constexpr mdspan(std::mdspan<element_type, typename detail::mdarray::std_extents_of<extents_type>::type, StdLayout> const& s)
            : ptr_(s.data_handle()), map_(from_std_mapping(s.mapping()))
        {
        }

        // Vers la disposition correspondante de la norme, ou vers std::layout_stride.
        template <class StdLayout>
            requires(std::is_same_v<accessor_type, default_accessor<element_type>>
                     && (std::is_same_v<StdLayout, typename detail::mdarray::std_layout<layout_type>::type>
                         || std::is_same_v<StdLayout, std::layout_stride>))
        constexpr operator std::mdspan<element_type, typename detail::mdarray::std_extents_of<extents_type>::type, StdLayout>() const
        {
            using std_extents = typename detail::mdarray::std_extents_of<extents_type>::type;
            auto const e = detail::mdarray::convert_extents<std_extents>(extents());
            if constexpr (std::is_same_v<StdLayout, std::layout_stride>)
            {
                std::array<index_type, rank()> s{};
                for (rank_type r = 0; r < rank(); ++r)
                    s[r] = map_.stride(r);
                return std::mdspan<element_type, std_extents, StdLayout>(ptr_, typename StdLayout::template mapping<std_extents>(e, s));
            }
            else
                return std::mdspan<element_type, std_extents, StdLayout>(ptr_, e);
        }

    private:
        template <class StdMapping>
        static constexpr mapping_type from_std_mapping(StdMapping const& m)
        {
            auto const e = detail::mdarray::convert_extents<extents_type>(m.extents());
            if constexpr (std::is_same_v<layout_type, layout_stride>)
            {
                std::array<index_type, rank()> s{};
                for (rank_type r = 0; r < rank(); ++r)
                    s[r] = m.stride(r);
                return mapping_type(e, s);
            }
            else
                return mapping_type(e);
        }

    public:
#endif

    private:
        data_handle_type ptr_{};
        [[no_unique_address]] mapping_type map_{};
        [[no_unique_address]] accessor_type acc_{};
    };

    // --- Guides de déduction (C++23) -------------------------------------------

    template <class CArray>
        requires(std::is_array_v<CArray> && std::rank_v<CArray> == 1)
    mdspan(CArray&) -> mdspan<std::remove_all_extents_t<CArray>, extents<std::size_t, std::extent_v<CArray, 0>>>;

    template <class Pointer>
        requires(std::is_pointer_v<std::remove_reference_t<Pointer>>)
    mdspan(Pointer&&) -> mdspan<std::remove_pointer_t<std::remove_reference_t<Pointer>>, extents<std::size_t>>;

    template <class ElementType, class... Integrals>
        requires((std::is_convertible_v<Integrals, std::size_t> && ...) && sizeof...(Integrals) > 0)
    explicit mdspan(ElementType*, Integrals...) -> mdspan<ElementType, dextents<std::size_t, sizeof...(Integrals)>>;

    template <class ElementType, class OtherIndexType, std::size_t N>
    mdspan(ElementType*, std::span<OtherIndexType, N>) -> mdspan<ElementType, dextents<std::size_t, N>>;

    template <class ElementType, class OtherIndexType, std::size_t N>
    mdspan(ElementType*, std::array<OtherIndexType, N> const&) -> mdspan<ElementType, dextents<std::size_t, N>>;

    template <class ElementType, class IndexType, std::size_t... ExtentsPack>
    mdspan(ElementType*, extents<IndexType, ExtentsPack...> const&) -> mdspan<ElementType, extents<IndexType, ExtentsPack...>>;

    template <class ElementType, class MappingType>
    mdspan(ElementType*, MappingType const&)
        -> mdspan<ElementType, typename MappingType::extents_type, typename MappingType::layout_type>;

    template <class MappingType, class AccessorType>
    mdspan(typename AccessorType::data_handle_type const&, MappingType const&, AccessorType const&)
        -> mdspan<typename AccessorType::element_type, typename MappingType::extents_type, typename MappingType::layout_type,
                  AccessorType>;

#if defined(__cpp_lib_mdspan)
    template <class T, class I, std::size_t... E, class StdLayout>
    mdspan(std::mdspan<T, std::extents<I, E...>, StdLayout> const&)
        -> mdspan<T, extents<I, E...>, typename detail::mdarray::wsl_layout<StdLayout>::type>;
#endif
}
