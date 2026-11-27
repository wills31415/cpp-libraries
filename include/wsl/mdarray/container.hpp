#pragma once

#include <wsl/modes.hpp>

// Tableau multidimensionnel possédant ses éléments (design/mdarray.md, I13) : la
// sémantique de std::mdarray proposée (P1684), valeur (copie profonde).
//
//   wsl::mdarray<double, wsl::dims<2>> a(3, 4);                       // 3 × 4, zéros
//   wsl::mdarray<double, wsl::extents<std::size_t, 3, 3>> g;          // 3 × 3, sans allocation
//   wsl::mdarray<double, wsl::dims<2>> m{{1, 2, 3}, {4, 5, 6}};       // littéral (rang ≥ 2)
//   m[1, 2];  m.at(1, 2);  m.to_mdspan();  equal(a, b);
//
// Conteneur par défaut : std::array si toutes les étendues sont fixes et la
// disposition toujours exhaustive (aucune allocation, constexpr complet), sinon
// std::vector (écart assumé avec P1684). Le littéral imbriqué n'existe qu'en
// rang ≥ 2 : en rang 1, {3} resterait ambigu (une étendue ou un élément ?).

#include <wsl/mdarray/detail/expression_base.hpp>
#include <wsl/mdarray/detail/iterate.hpp>
#include <wsl/mdarray/detail/loop.hpp>
#include <wsl/mdarray/detail/precondition.hpp>
#include <wsl/mdarray/extents.hpp>
#include <wsl/mdarray/layout.hpp>
#include <wsl/mdarray/mdspan.hpp>
#include <wsl/mdarray/position.hpp>

#include <algorithm>
#include <functional>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

// WSL_MDARRAY_INDEPENDENT_LOOP demande la vectorisation : Clang avertit quand elle est
// impossible (mode reproductible, fonctions), ce qui est attendu (detail/loop.hpp).
#if defined(__clang__)
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wpass-failed"
#endif

namespace wsl::inline WSL_DEFAULTS
{
    namespace detail::mdarray
    {
        template <class Extents>
        constexpr std::size_t static_size() noexcept
        {
            std::size_t n = 1;
            for (std::size_t r = 0; r < Extents::rank(); ++r)
                n *= Extents::static_extent(r);
            return n;
        }

        // Tampon contigu sur le tas, pour bool (std::vector<bool> ne stocke pas de
        // vrais bool et n'offre aucun pointeur sur ses éléments).
        template <class T>
        class heap_buffer
        {
        public:
            using value_type = T;
            using iterator = T*;
            using const_iterator = T const*;

            constexpr heap_buffer() = default;
            constexpr explicit heap_buffer(std::size_t n) : data_(n ? new T[n]() : nullptr), size_(n) {}
            constexpr heap_buffer(std::size_t n, T const& v) : heap_buffer(n) { std::ranges::fill(*this, v); }
            constexpr heap_buffer(heap_buffer const& o) : heap_buffer(o.size_) { std::ranges::copy(o, begin()); }
            constexpr heap_buffer(heap_buffer&& o) noexcept : data_(std::exchange(o.data_, nullptr)), size_(std::exchange(o.size_, 0)) {}
            constexpr heap_buffer& operator=(heap_buffer o) noexcept
            {
                std::swap(data_, o.data_);
                std::swap(size_, o.size_);
                return *this;
            }
            constexpr ~heap_buffer() { delete[] data_; }

            constexpr T* begin() noexcept { return data_; }
            constexpr T* end() noexcept { return data_ + size_; }
            constexpr T const* begin() const noexcept { return data_; }
            constexpr T const* end() const noexcept { return data_ + size_; }
            constexpr std::size_t size() const noexcept { return size_; }
            constexpr T& operator[](std::size_t k) noexcept { return data_[k]; }
            constexpr T const& operator[](std::size_t k) const noexcept { return data_[k]; }

        private:
            T* data_ = nullptr;
            std::size_t size_ = 0;
        };

        template <class T, class Extents, class Layout>
        struct default_container
        {
            using type = std::conditional_t<std::is_same_v<T, bool>, heap_buffer<bool>, std::vector<T>>;
        };

        template <class T, class Extents, class Layout>
            requires(Extents::rank_dynamic() == 0 && Layout::template mapping<Extents>::is_always_exhaustive())
        struct default_container<T, Extents, Layout>
        {
            using type = std::array<T, static_size<Extents>()>;
        };

        // std::initializer_list imbriquée R fois.
        template <class T, std::size_t R>
        struct nested_list
        {
            using type = std::initializer_list<typename nested_list<T, R - 1>::type>;
        };
        template <class T>
        struct nested_list<T, 1>
        {
            using type = std::initializer_list<T>;
        };
        struct no_literal
        {
        };
        template <class T>
        struct nested_list<T, 0>
        {
            using type = no_literal;
        };

        // Forme d'un littéral imbriqué (précondition : rectangulaire).
        template <std::size_t R, class List, class I>
        constexpr void list_shape(List const& list, std::array<I, R>& shape, std::size_t depth)
        {
            I const n = static_cast<I>(list.size());
            if (depth == 0 || shape[depth] == static_cast<I>(-1))
                shape[depth] = n;
            precondition(shape[depth] == n, "mdarray : littéral non rectangulaire");
            if constexpr (requires { list.begin()->begin(); })
                for (auto const& sub : list)
                    list_shape<R>(sub, shape, depth + 1);
        }

        template <class T>
        inline constexpr bool is_mdarray = false;

        template <class T>
        inline constexpr bool is_mdspan = false;
        template <class T, class E, class L, class A>
        inline constexpr bool is_mdspan<wsl::mdspan<T, E, L, A>> = true;
    }

    template <class ElementType, class Extents, class LayoutPolicy = layout_right,
              class Container = typename detail::mdarray::default_container<ElementType, Extents, LayoutPolicy>::type>
    class mdarray
    {
        static_assert(detail::mdarray::is_extents<Extents>, "wsl::mdarray : Extents doit être un wsl::extents");
        static_assert(std::is_same_v<typename Container::value_type, ElementType>,
                      "wsl::mdarray : le conteneur doit contenir des ElementType");

    public:
        using extents_type = Extents;
        using layout_type = LayoutPolicy;
        using container_type = Container;
        using mapping_type = typename layout_type::template mapping<extents_type>;
        using element_type = ElementType;
        using value_type = ElementType;
        using index_type = typename extents_type::index_type;
        using size_type = typename extents_type::size_type;
        using rank_type = typename extents_type::rank_type;
        using mdspan_type = wsl::mdspan<element_type, extents_type, layout_type>;
        using const_mdspan_type = wsl::mdspan<element_type const, extents_type, layout_type>;
        using pointer = element_type*;
        using const_pointer = element_type const*;
        using reference = element_type&;
        using const_reference = element_type const&;

        static constexpr rank_type rank() noexcept { return extents_type::rank(); }
        static constexpr rank_type rank_dynamic() noexcept { return extents_type::rank_dynamic(); }
        static constexpr std::size_t static_extent(rank_type r) noexcept { return extents_type::static_extent(r); }
        constexpr index_type extent(rank_type r) const noexcept { return extents().extent(r); }

        // --- Construction (P1684), éléments initialisés par valeur ------------------

        constexpr mdarray()
            requires(std::is_default_constructible_v<mapping_type>)
            : map_(), ctr_(make_container(map_))
        {
        }

        constexpr mdarray(mdarray const&) = default;
        constexpr mdarray(mdarray&&) = default;
        constexpr mdarray& operator=(mdarray const&) = default;
        constexpr mdarray& operator=(mdarray&&) = default;

        template <class... OtherIndexTypes>
            requires((std::is_convertible_v<OtherIndexTypes, index_type> && ...)
                     && (std::is_nothrow_constructible_v<index_type, OtherIndexTypes> && ...)
                     && (sizeof...(OtherIndexTypes) == rank() || sizeof...(OtherIndexTypes) == rank_dynamic())
                     && sizeof...(OtherIndexTypes) > 0 && std::is_constructible_v<mapping_type, extents_type>)
        constexpr explicit mdarray(OtherIndexTypes... exts)
            : map_(extents_type(static_cast<index_type>(exts)...)), ctr_(make_container(map_))
        {
        }

        constexpr explicit mdarray(extents_type const& e)
            requires std::is_constructible_v<mapping_type, extents_type const&>
            : map_(e), ctr_(make_container(map_))
        {
        }
        constexpr explicit mdarray(mapping_type const& m) : map_(m), ctr_(make_container(map_)) {}

        constexpr mdarray(extents_type const& e, value_type const& v)
            requires std::is_constructible_v<mapping_type, extents_type const&>
            : map_(e), ctr_(make_container(map_, v))
        {
        }
        constexpr mdarray(mapping_type const& m, value_type const& v) : map_(m), ctr_(make_container(map_, v)) {}

        constexpr mdarray(extents_type const& e, container_type const& c)
            requires std::is_constructible_v<mapping_type, extents_type const&>
            : map_(e), ctr_(c)
        {
            check_size();
        }
        constexpr mdarray(mapping_type const& m, container_type const& c) : map_(m), ctr_(c) { check_size(); }
        constexpr mdarray(extents_type const& e, container_type&& c)
            requires std::is_constructible_v<mapping_type, extents_type const&>
            : map_(e), ctr_(std::move(c))
        {
            check_size();
        }
        constexpr mdarray(mapping_type const& m, container_type&& c) : map_(m), ctr_(std::move(c)) { check_size(); }

        // Littéral imbriqué, en rang ≥ 2 : {{1, 2, 3}, {4, 5, 6}} ; les étendues
        // dynamiques sont prises au littéral, les fixes doivent coïncider (précondition).
        constexpr mdarray(typename detail::mdarray::nested_list<value_type, rank()>::type values)
            requires(rank() >= 2 && std::is_constructible_v<mapping_type, extents_type>)
            : map_(literal_extents(values)), ctr_(make_container(map_))
        {
            std::array<index_type, rank()> i{};
            fill_literal<0>(values, i);
        }

        // Copie des éléments d'une vue (P1684).
        template <class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherAccessor>
            requires(std::is_constructible_v<extents_type, OtherExtents>
                     && std::is_constructible_v<value_type, typename OtherAccessor::reference>)
        constexpr explicit(!std::is_convertible_v<OtherExtents, extents_type>)
            mdarray(wsl::mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherAccessor> const& other)
            : map_(extents_type(other.extents())), ctr_(make_container(map_))
        {
            copy_from(other);
        }

        // Autre type d'éléments, d'étendues ou de conteneur.
        template <class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherContainer>
            requires(!std::is_same_v<mdarray, wsl::mdarray<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherContainer>>
                     && std::is_constructible_v<extents_type, OtherExtents> && std::is_constructible_v<value_type, OtherElementType const&>)
        constexpr explicit(!std::is_convertible_v<OtherExtents, extents_type> || !std::is_convertible_v<OtherElementType const&, value_type>)
            mdarray(wsl::mdarray<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherContainer> const& other)
            : map_(extents_type(other.extents())), ctr_(make_container(map_))
        {
            copy_from(other.to_mdspan());
        }

        // Depuis une expression de tableaux entiers : mdarray c = a * b + 1; (V8)
        template <detail::mdarray::expression E>
            requires std::is_constructible_v<extents_type, typename std::remove_cvref_t<E>::extents_type>
        constexpr mdarray(E&& e) : map_(extents_type(e.extents())), ctr_(make_container(map_))
        {
            static_assert(!std::is_lvalue_reference_v<E>, "wsl::mdarray : une expression ne se garde pas (L2) ; std::move(e) si c'est voulu");
            assign_elements(detail::mdarray::default_mode(std::move(e)));
        }

        // y = expr : formes égales (précondition), pas de redimensionnement (I13).
        template <detail::mdarray::expression E>
        constexpr mdarray& operator=(E&& e)
        {
            static_assert(!std::is_lvalue_reference_v<E>, "wsl::mdarray : une expression ne se garde pas (L2) ; std::move(e) si c'est voulu");
            check_shape(e.extents());
            assign_elements(detail::mdarray::default_mode(std::move(e)));
            return *this;
        }

        // --- Accès : ceux de mdspan ----------------------------------------------

        template <class... I>
            requires requires(mdspan_type m, I const&... i) { m[i...]; }
        constexpr decltype(auto) operator[](I const&... indices)
        {
            return to_mdspan()[indices...];
        }

        template <class... I>
            requires requires(const_mdspan_type m, I const&... i) { m[i...]; }
        constexpr decltype(auto) operator[](I const&... indices) const
        {
            return to_mdspan()[indices...];
        }

        template <class... I>
            requires requires(mdspan_type m, I const&... i) { m.at(i...); }
        constexpr reference at(I const&... indices)
        {
            return to_mdspan().at(indices...);
        }

        template <class... I>
            requires requires(const_mdspan_type m, I const&... i) { m.at(i...); }
        constexpr const_reference at(I const&... indices) const
        {
            return to_mdspan().at(indices...);
        }

        template <class... A>
        constexpr bool contains(A const&... a) const noexcept
            requires requires(extents_type const& e) { e.contains(a...); }
        {
            return extents().contains(a...);
        }

        // --- Vue ---------------------------------------------------------------

        constexpr mdspan_type to_mdspan() { return mdspan_type(data(), map_); }
        constexpr const_mdspan_type to_mdspan() const { return const_mdspan_type(data(), map_); }

        template <class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherAccessor>
            requires std::is_assignable_v<wsl::mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherAccessor>&, mdspan_type>
        constexpr operator wsl::mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherAccessor>()
        {
            return to_mdspan();
        }

        template <class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherAccessor>
            requires std::is_assignable_v<wsl::mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherAccessor>&,
                                          const_mdspan_type>
        constexpr operator wsl::mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherAccessor>() const
        {
            return to_mdspan();
        }

        // --- Observateurs --------------------------------------------------------

        constexpr extents_type const& extents() const noexcept { return map_.extents(); }
        constexpr mapping_type const& mapping() const noexcept { return map_; }
        constexpr size_type size() const noexcept { return static_cast<size_type>(to_mdspan().size()); }
        [[nodiscard]] constexpr bool empty() const noexcept { return to_mdspan().empty(); }

        constexpr pointer data() noexcept { return std::to_address(ctr_.begin()); }
        constexpr const_pointer data() const noexcept { return std::to_address(ctr_.begin()); }
        constexpr container_type& container() noexcept { return ctr_; }
        constexpr container_type const& container() const noexcept { return ctr_; }
        constexpr container_type&& extract_container() && noexcept { return std::move(ctr_); }

        static constexpr bool is_always_unique() { return mapping_type::is_always_unique(); }
        static constexpr bool is_always_exhaustive() { return mapping_type::is_always_exhaustive(); }
        static constexpr bool is_always_strided() { return mapping_type::is_always_strided(); }
        constexpr bool is_unique() const { return map_.is_unique(); }
        constexpr bool is_exhaustive() const { return map_.is_exhaustive(); }
        constexpr bool is_strided() const { return map_.is_strided(); }
        constexpr index_type stride(rank_type r) const { return map_.stride(r); }

        // Parcours à plat du stockage (ajout, comme array2d).
        constexpr auto begin() noexcept { return ctr_.begin(); }
        constexpr auto end() noexcept { return ctr_.end(); }
        constexpr auto begin() const noexcept { return ctr_.begin(); }
        constexpr auto end() const noexcept { return ctr_.end(); }

        constexpr void fill(value_type const& v) { std::ranges::fill(ctr_, v); }

        friend constexpr void swap(mdarray& x, mdarray& y) noexcept
        {
            using std::swap;
            swap(x.map_, y.map_);
            swap(x.ctr_, y.ctr_);
        }

    private:
        static constexpr container_type make_container(mapping_type const& m)
        {
            if constexpr (requires { container_type(std::size_t{}); })
                return container_type(static_cast<std::size_t>(m.required_span_size()));
            else
            {
                detail::mdarray::precondition(std::tuple_size_v<container_type> >= static_cast<std::size_t>(m.required_span_size()),
                                              "mdarray : conteneur trop petit");
                return container_type{};
            }
        }

        static constexpr container_type make_container(mapping_type const& m, value_type const& v)
        {
            if constexpr (requires { container_type(std::size_t{}, v); })
                return container_type(static_cast<std::size_t>(m.required_span_size()), v);
            else
            {
                container_type c = make_container(m);
                std::ranges::fill(c, v);
                return c;
            }
        }

        constexpr void check_size() const noexcept
        {
            detail::mdarray::precondition(std::ranges::size(ctr_) >= static_cast<std::size_t>(map_.required_span_size()),
                                          "mdarray : conteneur trop petit pour la disposition");
        }

        template <class E>
        constexpr void check_shape(E const& e) const
        {
            static_assert(E::rank() == rank(), "wsl::mdarray : expression de rang différent");
            for (rank_type r = 0; r < rank(); ++r)
                detail::mdarray::precondition(std::cmp_equal(e.extent(r), extent(r)),
                                              "wsl::mdarray : affectation d'une expression d'une autre forme (I13)");
        }

        // Une opérande lit-elle notre stockage à d'autres positions que la nôtre ? (I12)
        template <class E>
        constexpr bool aliased(E const& e) const
        {
            if consteval
            {
                return true; // pas d'ordre entre pointeurs en évaluation constante : par prudence
            }
            else
            {
                bool result = false;
                void const* lo = data();
                void const* hi = data() + std::ranges::size(ctr_);
                e.for_each_storage([&](auto const& view) {
                    void const* p = view.data_handle();
                    void const* q = view.data_handle() + view.mapping().required_span_size();
                    std::less<void const*> less;
                    bool const overlap = less(p, hi) && less(lo, q);
                    if (!overlap)
                        return;
                    if constexpr (std::is_same_v<typename std::remove_cvref_t<decltype(view)>::mapping_type, mapping_type>)
                        if (view.data_handle() == data() && view.mapping() == map_)
                            return; // mêmes positions : sûr
                    result = true;
                });
                return result;
            }
        }

        template <class E>
        constexpr void assign_elements(E const& e)
        {
            if (aliased(e))
            {
                container_type temporary = make_container(map_);
                mdspan_type t(std::to_address(temporary.begin()), map_);
                detail::mdarray::for_each_index(extents(), [&](auto... i) { t[i...] = static_cast<value_type>(e(i...)); });
                ctr_ = std::move(temporary);
            }
            else if constexpr (std::is_same_v<layout_type, layout_right> && requires { E::flat_ok; })
            {
                if constexpr (E::flat_ok)
                {
                    // Ligne majeure de part et d'autre : une boucle sur le stockage.
                    pointer const d = data();
                    std::size_t const n = static_cast<std::size_t>(size());
                    std::size_t k = 0;
                    if constexpr (requires { E::lane_barrier; })
                        if constexpr (E::lane_barrier)
                        {
                            // Produit reproductible : par paquets de 8, une barrière par
                            // registre vectoriel (mêmes bits que la boucle simple).
                            constexpr std::size_t W = 8;
                            for (; k + W <= n; k += W)
                            {
                                auto const p = e.template flat_lanes<W>(k);
                                for (std::size_t l = 0; l < W; ++l)
                                    d[k + l] = static_cast<value_type>(p.v[l]);
                            }
                        }
                    WSL_MDARRAY_INDEPENDENT_LOOP
                    for (; k < n; ++k)
                        d[k] = static_cast<value_type>(e.flat(k));
                }
                else
                {
                    auto self = to_mdspan();
                    detail::mdarray::for_each_index(extents(), [&](auto... i) { self[i...] = static_cast<value_type>(e(i...)); });
                }
            }
            else
            {
                auto self = to_mdspan();
                detail::mdarray::for_each_index(extents(), [&](auto... i) { self[i...] = static_cast<value_type>(e(i...)); });
            }
        }

        template <class View>
        constexpr void copy_from(View const& other)
        {
            auto self = to_mdspan();
            detail::mdarray::for_each_index(extents(), [&](auto... i) { self[i...] = static_cast<value_type>(other[i...]); });
        }

        static constexpr extents_type literal_extents(typename detail::mdarray::nested_list<value_type, rank()>::type values)
        {
            std::array<index_type, rank()> shape{};
            for (auto& s : shape)
                s = static_cast<index_type>(-1);
            detail::mdarray::list_shape<rank()>(values, shape, 0);
            for (auto& s : shape)
                if (s == static_cast<index_type>(-1))
                    s = 0;
            return extents_type(shape);
        }

        template <std::size_t Depth, class List>
        constexpr void fill_literal(List const& list, std::array<index_type, rank()>& i)
        {
            i[Depth] = 0;
            for (auto const& item : list)
            {
                if constexpr (Depth + 1 == rank())
                    to_mdspan()[i] = item;
                else
                    fill_literal<Depth + 1>(item, i);
                ++i[Depth];
            }
        }

        [[no_unique_address]] mapping_type map_;
        container_type ctr_;
    };

    namespace detail::mdarray
    {
        template <class T, class E, class L, class C>
        inline constexpr bool is_mdarray<wsl::mdarray<T, E, L, C>> = true;

        template <class T>
        concept array_like = is_mdarray<std::remove_cvref_t<T>> || is_mdspan<std::remove_cvref_t<T>>;

        template <class T>
        constexpr auto as_view(T const& t)
        {
            if constexpr (is_mdarray<T>)
                return t.to_mdspan();
            else
                return t;
        }
    }

    template <class T, class E, class L, class A>
    mdarray(wsl::mdspan<T, E, L, A> const&) -> mdarray<std::remove_cv_t<T>, E, L>;

    // Égalité de deux tableaux ou vues : mêmes étendues, mêmes éléments aux mêmes
    // positions (I13 ; `==` sera élément par élément, comme std::valarray).
    template <detail::mdarray::array_like A, detail::mdarray::array_like B>
        requires(std::remove_cvref_t<A>::rank() == std::remove_cvref_t<B>::rank())
    constexpr bool equal(A const& a, B const& b)
    {
        auto const x = detail::mdarray::as_view(a);
        auto const y = detail::mdarray::as_view(b);
        if (!(x.extents() == y.extents()))
            return false;
        return detail::mdarray::for_each_index(x.extents(), [&](auto... i) { return x[i...] == y[i...]; });
    }
}

#if defined(__clang__)
#    pragma clang diagnostic pop
#endif
