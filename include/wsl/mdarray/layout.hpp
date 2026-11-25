#pragma once

// Dispositions : correspondance entre une position et une place en mémoire
// (docs/mdarray.md, A5). Noms et sémantique de <mdspan> (C++23) : layout_right
// (ligne majeure, le dernier indice varie le plus vite), layout_left (colonne
// majeure), layout_stride (pas quelconques).

#include <wsl/mdarray/detail/precondition.hpp>
#include <wsl/mdarray/extents.hpp>

#include <array>
#include <cstddef>
#include <span>
#include <type_traits>
#include <utility>

namespace wsl
{
    struct layout_right
    {
        template <class Extents>
        class mapping;
    };

    struct layout_left
    {
        template <class Extents>
        class mapping;
    };

    struct layout_stride
    {
        template <class Extents>
        class mapping;
    };

    namespace detail::mdarray
    {
        // Produit des étendues [begin, end), précondition : représentable.
        template <class Extents>
        constexpr typename Extents::index_type extents_product(Extents const& e, std::size_t begin, std::size_t end) noexcept
        {
            using I = typename Extents::index_type;
            I p{1};
            for (std::size_t r = begin; r < end; ++r)
            {
                I next{};
                precondition(!__builtin_mul_overflow(p, e.extent(r), &next),
                             "disposition : taille de l'ensemble d'indices non représentable");
                p = next;
            }
            return p;
        }

        // Toutes les étendues sont non nulles et leur produit, représentable.
        template <class Extents>
        constexpr bool size_representable(Extents const& e) noexcept
        {
            using I = typename Extents::index_type;
            I p{1};
            for (std::size_t r = 0; r < Extents::rank(); ++r)
            {
                if (e.extent(r) == 0)
                    return true;
                if (__builtin_mul_overflow(p, e.extent(r), &p))
                    return false;
            }
            return true;
        }

        template <class Mapping, class... Indices>
        constexpr void check_indices(typename Mapping::extents_type const& e, Indices... i) noexcept
        {
            using I = typename Mapping::index_type;
            precondition(e.contains(static_cast<I>(i)...), "disposition : position hors de l'ensemble d'indices");
        }

        template <class T, template <class> class M>
        inline constexpr bool is_mapping_of = false;
        template <class E, template <class> class M>
        inline constexpr bool is_mapping_of<M<E>, M> = true;

        template <class M>
        inline constexpr bool is_wsl_mapping =
            std::is_same_v<typename M::layout_type, wsl::layout_right> || std::is_same_v<typename M::layout_type, wsl::layout_left>
            || std::is_same_v<typename M::layout_type, wsl::layout_stride>;

        template <class M>
        concept strided_wsl_mapping = requires { typename M::layout_type; } && is_wsl_mapping<M> && M::is_always_strided();
    }

    // -----------------------------------------------------------------------
    // layout_right
    // -----------------------------------------------------------------------

    template <class Extents>
    class layout_right::mapping
    {
        static_assert(detail::mdarray::is_extents<Extents>, "layout_right::mapping : Extents doit être un wsl::extents");

    public:
        using extents_type = Extents;
        using index_type = typename extents_type::index_type;
        using size_type = typename extents_type::size_type;
        using rank_type = typename extents_type::rank_type;
        using layout_type = layout_right;

        constexpr mapping() noexcept = default;
        constexpr mapping(mapping const&) noexcept = default;

        constexpr mapping(extents_type const& e) noexcept : extents_(e)
        {
            detail::mdarray::precondition(detail::mdarray::size_representable(e),
                                          "layout_right : taille de l'ensemble d'indices non représentable");
        }

        template <class OtherExtents>
            requires std::is_constructible_v<extents_type, OtherExtents>
        constexpr explicit(!std::is_convertible_v<OtherExtents, extents_type>) mapping(mapping<OtherExtents> const& other) noexcept
            : mapping(extents_type(other.extents()))
        {
        }

        template <class OtherExtents>
            requires(extents_type::rank() <= 1 && std::is_constructible_v<extents_type, OtherExtents>)
        constexpr explicit(!std::is_convertible_v<OtherExtents, extents_type>)
            mapping(layout_left::mapping<OtherExtents> const& other) noexcept
            : mapping(extents_type(other.extents()))
        {
        }

        template <class OtherExtents>
            requires std::is_constructible_v<extents_type, OtherExtents>
        constexpr explicit(extents_type::rank() > 0) mapping(layout_stride::mapping<OtherExtents> const& other) noexcept
            : mapping(extents_type(other.extents()))
        {
            for (rank_type r = 0; r < extents_type::rank(); ++r)
                detail::mdarray::precondition(std::cmp_equal(other.stride(r), stride(r)), "layout_right : pas incompatibles");
        }

        constexpr mapping& operator=(mapping const&) noexcept = default;

        constexpr extents_type const& extents() const noexcept { return extents_; }

        constexpr index_type required_span_size() const noexcept
        {
            return detail::mdarray::extents_product(extents_, 0, extents_type::rank());
        }

        template <class... Indices>
            requires(sizeof...(Indices) == extents_type::rank() && (std::is_convertible_v<Indices, index_type> && ...)
                     && (std::is_nothrow_constructible_v<index_type, Indices> && ...))
        constexpr index_type operator()(Indices... indices) const noexcept
        {
            detail::mdarray::check_indices<mapping>(extents_, indices...);
            index_type offset{0};
            [[maybe_unused]] rank_type r = 0;
            ((offset = static_cast<index_type>(offset * extents_.extent(r++) + static_cast<index_type>(indices))), ...);
            return offset;
        }

        static constexpr bool is_always_unique() noexcept { return true; }
        static constexpr bool is_always_exhaustive() noexcept { return true; }
        static constexpr bool is_always_strided() noexcept { return true; }
        static constexpr bool is_unique() noexcept { return true; }
        static constexpr bool is_exhaustive() noexcept { return true; }
        static constexpr bool is_strided() noexcept { return true; }

        constexpr index_type stride(rank_type r) const noexcept
            requires(extents_type::rank() > 0)
        {
            detail::mdarray::precondition(r < extents_type::rank(), "layout_right::stride : rang hors limites");
            return detail::mdarray::extents_product(extents_, r + 1, extents_type::rank());
        }

        template <class OtherExtents>
            requires(OtherExtents::rank() == extents_type::rank())
        friend constexpr bool operator==(mapping const& a, mapping<OtherExtents> const& b) noexcept
        {
            return a.extents() == b.extents();
        }

    private:
        [[no_unique_address]] extents_type extents_{};
    };

    // -----------------------------------------------------------------------
    // layout_left
    // -----------------------------------------------------------------------

    template <class Extents>
    class layout_left::mapping
    {
        static_assert(detail::mdarray::is_extents<Extents>, "layout_left::mapping : Extents doit être un wsl::extents");

    public:
        using extents_type = Extents;
        using index_type = typename extents_type::index_type;
        using size_type = typename extents_type::size_type;
        using rank_type = typename extents_type::rank_type;
        using layout_type = layout_left;

        constexpr mapping() noexcept = default;
        constexpr mapping(mapping const&) noexcept = default;

        constexpr mapping(extents_type const& e) noexcept : extents_(e)
        {
            detail::mdarray::precondition(detail::mdarray::size_representable(e),
                                          "layout_left : taille de l'ensemble d'indices non représentable");
        }

        template <class OtherExtents>
            requires std::is_constructible_v<extents_type, OtherExtents>
        constexpr explicit(!std::is_convertible_v<OtherExtents, extents_type>) mapping(mapping<OtherExtents> const& other) noexcept
            : mapping(extents_type(other.extents()))
        {
        }

        template <class OtherExtents>
            requires(extents_type::rank() <= 1 && std::is_constructible_v<extents_type, OtherExtents>)
        constexpr explicit(!std::is_convertible_v<OtherExtents, extents_type>)
            mapping(layout_right::mapping<OtherExtents> const& other) noexcept
            : mapping(extents_type(other.extents()))
        {
        }

        template <class OtherExtents>
            requires std::is_constructible_v<extents_type, OtherExtents>
        constexpr explicit(extents_type::rank() > 0) mapping(layout_stride::mapping<OtherExtents> const& other) noexcept
            : mapping(extents_type(other.extents()))
        {
            for (rank_type r = 0; r < extents_type::rank(); ++r)
                detail::mdarray::precondition(std::cmp_equal(other.stride(r), stride(r)), "layout_left : pas incompatibles");
        }

        constexpr mapping& operator=(mapping const&) noexcept = default;

        constexpr extents_type const& extents() const noexcept { return extents_; }

        constexpr index_type required_span_size() const noexcept
        {
            return detail::mdarray::extents_product(extents_, 0, extents_type::rank());
        }

        template <class... Indices>
            requires(sizeof...(Indices) == extents_type::rank() && (std::is_convertible_v<Indices, index_type> && ...)
                     && (std::is_nothrow_constructible_v<index_type, Indices> && ...))
        constexpr index_type operator()(Indices... indices) const noexcept
        {
            detail::mdarray::check_indices<mapping>(extents_, indices...);
            std::array<index_type, sizeof...(Indices)> const i{static_cast<index_type>(indices)...};
            index_type offset{0};
            for (rank_type r = extents_type::rank(); r-- > 0;)
                offset = static_cast<index_type>(offset * extents_.extent(r) + i[r]);
            return offset;
        }

        static constexpr bool is_always_unique() noexcept { return true; }
        static constexpr bool is_always_exhaustive() noexcept { return true; }
        static constexpr bool is_always_strided() noexcept { return true; }
        static constexpr bool is_unique() noexcept { return true; }
        static constexpr bool is_exhaustive() noexcept { return true; }
        static constexpr bool is_strided() noexcept { return true; }

        constexpr index_type stride(rank_type r) const noexcept
            requires(extents_type::rank() > 0)
        {
            detail::mdarray::precondition(r < extents_type::rank(), "layout_left::stride : rang hors limites");
            return detail::mdarray::extents_product(extents_, 0, r);
        }

        template <class OtherExtents>
            requires(OtherExtents::rank() == extents_type::rank())
        friend constexpr bool operator==(mapping const& a, mapping<OtherExtents> const& b) noexcept
        {
            return a.extents() == b.extents();
        }

    private:
        [[no_unique_address]] extents_type extents_{};
    };

    // -----------------------------------------------------------------------
    // layout_stride
    // -----------------------------------------------------------------------

    template <class Extents>
    class layout_stride::mapping
    {
        static_assert(detail::mdarray::is_extents<Extents>, "layout_stride::mapping : Extents doit être un wsl::extents");

    public:
        using extents_type = Extents;
        using index_type = typename extents_type::index_type;
        using size_type = typename extents_type::size_type;
        using rank_type = typename extents_type::rank_type;
        using layout_type = layout_stride;

        // Les pas de layout_right pour extents_type() (comme la norme).
        constexpr mapping() noexcept : mapping(layout_right::mapping<extents_type>()) {}
        constexpr mapping(mapping const&) noexcept = default;

        template <class OtherIndexType>
            requires(std::is_convertible_v<OtherIndexType const&, index_type>
                     && std::is_nothrow_constructible_v<index_type, OtherIndexType const&>)
        constexpr mapping(extents_type const& e, std::span<OtherIndexType, extents_type::rank()> s) noexcept : extents_(e)
        {
            for (rank_type r = 0; r < extents_type::rank(); ++r)
            {
                detail::mdarray::precondition(std::cmp_greater(s[r], 0), "layout_stride : pas non positif");
                strides_[r] = detail::mdarray::index_cast<index_type>(s[r], "layout_stride : pas non représentable");
            }
            detail::mdarray::precondition(unique_strides(), "layout_stride : correspondance non injective");
        }

        template <class OtherIndexType>
            requires(std::is_convertible_v<OtherIndexType const&, index_type>
                     && std::is_nothrow_constructible_v<index_type, OtherIndexType const&>)
        constexpr mapping(extents_type const& e, std::array<OtherIndexType, extents_type::rank()> const& s) noexcept
            : mapping(e, std::span<OtherIndexType const, extents_type::rank()>(s))
        {
        }

        // Depuis toute disposition à pas de wsl (layout_right, layout_left, layout_stride).
        template <class StridedMapping>
            requires(detail::mdarray::strided_wsl_mapping<StridedMapping>
                     && std::is_constructible_v<extents_type, typename StridedMapping::extents_type>
                     && StridedMapping::extents_type::rank() == extents_type::rank())
        constexpr explicit(!std::is_convertible_v<typename StridedMapping::extents_type, extents_type>)
            mapping(StridedMapping const& other) noexcept
            : extents_(other.extents())
        {
            if constexpr (extents_type::rank() > 0)
                for (rank_type r = 0; r < extents_type::rank(); ++r)
                    strides_[r] = detail::mdarray::index_cast<index_type>(other.stride(r), "layout_stride : pas non représentable");
        }

        constexpr mapping& operator=(mapping const&) noexcept = default;

        constexpr extents_type const& extents() const noexcept { return extents_; }
        constexpr std::array<index_type, extents_type::rank()> strides() const noexcept { return strides_; }

        constexpr index_type required_span_size() const noexcept
        {
            if constexpr (extents_type::rank() == 0)
                return 1;
            else
            {
                for (rank_type r = 0; r < extents_type::rank(); ++r)
                    if (extents_.extent(r) == 0)
                        return 0;
                index_type size{1};
                for (rank_type r = 0; r < extents_type::rank(); ++r)
                    size = static_cast<index_type>(size + (extents_.extent(r) - 1) * strides_[r]);
                return size;
            }
        }

        template <class... Indices>
            requires(sizeof...(Indices) == extents_type::rank() && (std::is_convertible_v<Indices, index_type> && ...)
                     && (std::is_nothrow_constructible_v<index_type, Indices> && ...))
        constexpr index_type operator()(Indices... indices) const noexcept
        {
            detail::mdarray::check_indices<mapping>(extents_, indices...);
            index_type offset{0};
            [[maybe_unused]] rank_type r = 0;
            ((offset = static_cast<index_type>(offset + static_cast<index_type>(indices) * strides_[r++])), ...);
            return offset;
        }

        static constexpr bool is_always_unique() noexcept { return true; }
        static constexpr bool is_always_exhaustive() noexcept { return false; }
        static constexpr bool is_always_strided() noexcept { return true; }
        static constexpr bool is_unique() noexcept { return true; }
        static constexpr bool is_strided() noexcept { return true; }

        // Exhaustive : toute place de [0, required_span_size()) est atteinte
        // ([mdspan.layout.reqmts]). La correspondance étant injective (précondition),
        // c'est exactement required_span_size() == taille de l'ensemble d'indices ;
        // vrai pour un ensemble vide (0 == 0, [mdspan.layout.stride.obs] 6.1). La
        // formulation par permutation (6.2) dit faux pour une longueur 1 de pas
        // quelconque, où rien ne manque pourtant : on suit la propriété, comme libc++.
        constexpr bool is_exhaustive() const noexcept
        {
            if constexpr (extents_type::rank() == 0)
                return true;
            else
            {
                index_type size{1};
                for (rank_type r = 0; r < extents_type::rank(); ++r)
                    size = static_cast<index_type>(size * extents_.extent(r));
                return required_span_size() == size;
            }
        }

        constexpr index_type stride(rank_type r) const noexcept
        {
            detail::mdarray::precondition(r < extents_type::rank(), "layout_stride::stride : rang hors limites");
            return strides_[r];
        }

        template <class OtherMapping>
            requires(detail::mdarray::strided_wsl_mapping<OtherMapping>
                     && OtherMapping::extents_type::rank() == extents_type::rank())
        friend constexpr bool operator==(mapping const& a, OtherMapping const& b) noexcept
        {
            if (!(a.extents() == b.extents()))
                return false;
            if constexpr (extents_type::rank() > 0)
                for (rank_type r = 0; r < extents_type::rank(); ++r)
                    if (!std::cmp_equal(a.stride(r), b.stride(r)))
                        return false;
            return true;
        }

    private:
        constexpr bool less(rank_type a, rank_type b) const noexcept
        {
            return strides_[a] < strides_[b] || (strides_[a] == strides_[b] && extents_.extent(a) < extents_.extent(b));
        }

        // Précondition de la norme : il existe un ordre P des rangs où
        // s[pᵢ] ≥ s[pᵢ₋₁] · e(pᵢ₋₁) (pas emboîtés, donc correspondance injective).
        constexpr bool unique_strides() const noexcept
        {
            if constexpr (extents_type::rank() <= 1)
                return true;
            else
            {
                std::array<rank_type, extents_type::rank()> order{};
                for (rank_type r = 0; r < extents_type::rank(); ++r)
                    order[r] = r;
                for (rank_type a = 1; a < extents_type::rank(); ++a)
                    for (rank_type b = a; b > 0 && less(order[b], order[b - 1]); --b)
                        std::swap(order[b], order[b - 1]);
                for (rank_type k = 1; k < extents_type::rank(); ++k)
                {
                    index_type previous{};
                    if (__builtin_mul_overflow(strides_[order[k - 1]], extents_.extent(order[k - 1]), &previous))
                        return false;
                    if (strides_[order[k]] < previous)
                        return false;
                }
                return true;
            }
        }

        [[no_unique_address]] extents_type extents_{};
        std::array<index_type, extents_type::rank()> strides_{};
    };
}
