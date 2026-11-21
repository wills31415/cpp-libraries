#pragma once

// Tableau 2D de taille fixe R × C (agrégat, comme std::array), disposition par lignes
// (row_major, défaut) ou par colonnes (column_major). Accès m[i, j] (C++23), dimensions
// extent(0) = R et extent(1) = C comme std::mdspan ; vues de lignes et de colonnes
// (std::span si contiguës, strided_view sinon), to_mdspan().
//
//   wsl::array2d<int, 2, 3> m = wsl::array2d<int, 2, 3>::from_rows({{1, 2, 3}, {4, 5, 6}});
//   m[1, 2];                      // 6
//   for (int x : m.column(1)) …   // 2, 5
//   auto t = m.transposed();      // 3 × 2

#include <algorithm>
#include <compare>          // IWYU pragma: keep
#include <concepts>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <version>
#if __has_include(<mdspan>)
#include <mdspan>
#endif

namespace wsl
{

    struct row_major
    {
        template <std::size_t R, std::size_t C>
        static constexpr std::size_t index(std::size_t i, std::size_t j) noexcept {return i * C + j;}

#ifdef __cpp_lib_mdspan
        using mdspan_layout = std::layout_right;
#endif
    };

    struct column_major
    {
        template <std::size_t R, std::size_t C>
        static constexpr std::size_t index(std::size_t i, std::size_t j) noexcept {return j * R + i;}

#ifdef __cpp_lib_mdspan
        using mdspan_layout = std::layout_left;
#endif
    };


    template <typename L, std::size_t R, std::size_t C>
    concept layout2d = requires(std::size_t i, std::size_t j)
    {
        {L::template index<R, C>(i, j)} noexcept -> std::same_as<std::size_t>;
    };


    template <typename T>
    class strided_view : public std::ranges::view_interface<strided_view<T>>
    {
    public:
        class iterator
        {
        public:
            using iterator_concept  = std::random_access_iterator_tag;
            using iterator_category = std::random_access_iterator_tag;
            using value_type        = std::remove_cv_t<T>;
            using difference_type   = std::ptrdiff_t;
            using reference         = T&;

            iterator() = default;

            constexpr iterator(T* base, difference_type index, difference_type stride) noexcept
            : base_(base), index_(index), stride_(stride)
            {}

            constexpr T& operator * () const noexcept {return base_[index_ * stride_];}
            constexpr T* operator -> () const noexcept {return base_ + index_ * stride_;}
            constexpr T& operator [] (difference_type n) const noexcept {return base_[(index_ + n) * stride_];}

            constexpr iterator& operator ++ () noexcept                     {++index_; return *this;}
            constexpr iterator  operator ++ (int) noexcept                  {auto tmp = *this; ++index_; return tmp;}
            constexpr iterator& operator -- () noexcept                     {--index_; return *this;}
            constexpr iterator  operator -- (int) noexcept                  {auto tmp = *this; --index_; return tmp;}
            constexpr iterator& operator += (difference_type n) noexcept    {index_ += n; return *this;}
            constexpr iterator& operator -= (difference_type n) noexcept    {index_ -= n; return *this;}

            friend constexpr iterator           operator + (iterator it, difference_type n) noexcept {return it += n;}
            friend constexpr iterator           operator + (difference_type n, iterator it) noexcept {return it += n;}
            friend constexpr iterator           operator - (iterator it, difference_type n) noexcept {return it -= n;}
            friend constexpr difference_type    operator -   (const iterator& a, const iterator& b) noexcept {return a.index_ - b.index_;}
            friend constexpr bool               operator ==  (const iterator& a, const iterator& b) noexcept {return a.index_ == b.index_;}
            friend constexpr auto               operator <=> (const iterator& a, const iterator& b) noexcept {return a.index_ <=> b.index_;}

        private:
            T* base_ = nullptr;
            difference_type index_ = 0;
            difference_type stride_ = 1;
        };

        strided_view() = default;

        constexpr strided_view(T* base, std::size_t count, std::ptrdiff_t stride) noexcept
        : base_(base), count_(count), stride_(stride)
        {}

        constexpr iterator begin() const noexcept {return iterator(base_, 0, stride_);}
        constexpr iterator end() const noexcept {return iterator(base_, static_cast<std::ptrdiff_t>(count_), stride_);}
        constexpr std::size_t size() const noexcept {return count_;}
        constexpr std::ptrdiff_t stride() const noexcept {return stride_;}

    private:
        T* base_ = nullptr;
        std::size_t count_ = 0;
        std::ptrdiff_t stride_ = 1;
    };


    template <typename T, std::size_t R, std::size_t C, typename Layout = row_major>
    requires layout2d<Layout, R, C>
    struct array2d
    {
        static_assert(R > 0 && C > 0, "wsl::array2d : dimensions nulles non prises en charge");

        using value_type                = T;
        using layout_type               = Layout;
        using size_type                 = std::size_t;
        using difference_type           = std::ptrdiff_t;
        using reference                 = T&;
        using const_reference           = T const&;
        using pointer                   = T*;
        using const_pointer             = T const*;
        using iterator                  = T*;
        using const_iterator            = T const*;
        using reverse_iterator          = std::reverse_iterator<iterator>;
        using const_reverse_iterator    = std::reverse_iterator<const_iterator>;

        T elems_[R * C];

        static constexpr size_type rank() noexcept              {return 2;}
        static constexpr size_type extent(size_type i) noexcept {return i == 0 ? R : C;}
        static constexpr size_type size() noexcept              {return R * C;}
        static constexpr size_type max_size() noexcept          {return R * C;}
        [[nodiscard]] static constexpr bool empty() noexcept    {return false;}

        constexpr reference         operator [] (size_type i, size_type j) noexcept {return elems_[Layout::template index<R, C>(i, j)];}
        constexpr const_reference   operator [] (size_type i, size_type j) const noexcept {return elems_[Layout::template index<R, C>(i, j)];}
        constexpr reference         at(size_type i, size_type j) {check_bounds(i, j); return (*this)[i, j];}
        constexpr const_reference   at(size_type i, size_type j) const {check_bounds(i, j); return (*this)[i, j];}
        constexpr T*                data() noexcept {return elems_;}
        constexpr T const*          data() const noexcept {return elems_;}

        constexpr iterator begin() noexcept {return elems_;}
        constexpr const_iterator begin() const noexcept {return elems_;}
        constexpr iterator end() noexcept {return elems_ + R * C;}
        constexpr const_iterator end() const noexcept {return elems_ + R * C;}
        constexpr const_iterator cbegin() const noexcept {return begin();}
        constexpr const_iterator cend() const noexcept {return end();}
        constexpr reverse_iterator rbegin() noexcept {return reverse_iterator(end());}
        constexpr const_reverse_iterator rbegin() const noexcept {return const_reverse_iterator(end());}
        constexpr reverse_iterator rend() noexcept {return reverse_iterator(begin());}
        constexpr const_reverse_iterator rend() const noexcept {return const_reverse_iterator(begin());}
        constexpr const_reverse_iterator crbegin() const noexcept {return rbegin();}
        constexpr const_reverse_iterator crend() const noexcept {return rend();}

        static constexpr difference_type row_stride = static_cast<difference_type>(Layout::template index<R, C>(1, 0)) - static_cast<difference_type>(Layout::template index<R, C>(0, 0));
        static constexpr difference_type column_stride = static_cast<difference_type>(Layout::template index<R, C>(0, 1)) - static_cast<difference_type>(Layout::template index<R, C>(0, 0));

        template <typename Self>
        constexpr auto row(this Self& self, size_type i) noexcept
        {
            using E = std::remove_reference_t<decltype(self.elems_[0])>;
            E* first = self.elems_ + Layout::template index<R, C>(i, 0);
            if constexpr (column_stride == 1) return std::span<E, C>(first, C);
            else return strided_view<E>(first, C, column_stride);
        }

        template <typename Self>
        constexpr auto column(this Self& self, size_type j) noexcept
        {
            using E = std::remove_reference_t<decltype(self.elems_[0])>;
            E* first = self.elems_ + Layout::template index<R, C>(0, j);
            if constexpr (row_stride == 1) return std::span<E, R>(first, R);
            else return strided_view<E>(first, R, row_stride);
        }

        template <typename Self>
        constexpr auto rows(this Self& self) noexcept
        {
            return std::views::iota(size_type{0}, R) | std::views::transform([&self](size_type i){return self.row(i);});
        }

        template <typename Self>
        constexpr auto columns(this Self& self) noexcept
        {
            return std::views::iota(size_type{0}, C) | std::views::transform([&self](size_type j){return self.column(j);});
        }

#ifdef __cpp_lib_mdspan
        constexpr auto to_mdspan() noexcept requires requires {typename Layout::mdspan_layout;}
        {
            return std::mdspan<T, std::extents<size_type, R, C>, typename Layout::mdspan_layout>(elems_);
        }
        constexpr auto to_mdspan() const noexcept requires requires {typename Layout::mdspan_layout;}
        {
            return std::mdspan<T const, std::extents<size_type, R, C>, typename Layout::mdspan_layout>(elems_);
        }
#endif

        static constexpr array2d from_rows(T const (&values)[R][C])
        requires std::default_initializable<T> && std::copyable<T>
        {
            array2d result{};
            for (size_type i = 0; i < R; ++i)
                for (size_type j = 0; j < C; ++j)
                    result[i, j] = values[i][j];
            return result;
        }

        constexpr void fill(T const& value) {std::ranges::fill(elems_, value);}

        constexpr void swap(array2d& other) noexcept(std::is_nothrow_swappable_v<T>)
        {
            std::ranges::swap_ranges(elems_, other.elems_);
        }
        friend constexpr void swap(array2d& a, array2d& b) noexcept(std::is_nothrow_swappable_v<T>)
        {
            a.swap(b);
        }

        constexpr array2d<T, C, R, Layout> transposed() const
        requires std::default_initializable<T> && std::copyable<T>
        {
            array2d<T, C, R, Layout> result{};
            for (size_type i = 0; i < R; ++i)
                for (size_type j = 0; j < C; ++j)
                    result[j, i] = (*this)[i, j];
            return result;
        }

        friend constexpr bool operator ==  (array2d const&, array2d const&) = default;

        // Pas de relation d'ordre canonique entre deux tableaux 2D : l'opérateur est volontairement supprimé.
        friend constexpr auto operator <=> (array2d const&, array2d const&) = delete;

    private:
        static constexpr void check_bounds(size_type i, size_type j)
        {
            if (i >= R || j >= C) throw std::out_of_range("wsl::array2d::at : indice hors bornes");
        }
    };

}   // namespace wsl


template <typename T>
inline constexpr bool std::ranges::enable_borrowed_range<wsl::strided_view<T>> = true;
