#pragma once

// Positions et déplacements dans un tableau multidimensionnel (docs/mdarray.md, §4).
//
//   wsl::mdposition p{2, 1};                     // mdposition<2> en std::size_t
//   wsl::mddisplacement d{0, -1};                // mddisplacement<2> en std::ptrdiff_t
//   auto q = p + d;                              // position (2, 0)
//   auto e = q - p;                              // déplacement (0, -1)
//   p + wsl::mddisplacement<2>::unit(0);         // voisin dans la direction 0
//
// Modèle affine, celui de std::chrono (time_point, duration) : position − position
// = déplacement, position ± déplacement = position. On n'ajoute, ne multiplie ni ne
// divise des positions. Les déplacements forment un ℤ-module : d ± e, −d, et `*`,
// `/`, `%` par un entier (troncature vers zéro, comme les entiers de C++).
// Ordre : lexicographique, total (comme std::array) ; ordre composante par
// composante : all_less, all_less_equal. Préconditions (P3) : résultat
// représentable, composantes d'une position positives.

#include <wsl/mdarray/detail/precondition.hpp>
#include <wsl/mdarray/integer.hpp>

#include <array>
#include <compare>
#include <cstddef>
#include <format>
#include <functional>
#include <limits>
#include <ostream>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>

namespace wsl
{
    template <std::size_t N, class IndexType = std::size_t>
    class mdposition;

    template <std::size_t N, class DiffType = std::ptrdiff_t>
    class mddisplacement;

    namespace detail::mdarray
    {
        template <class T>
        inline constexpr bool is_mdposition = false;
        template <std::size_t N, class I>
        inline constexpr bool is_mdposition<wsl::mdposition<N, I>> = true;

        template <class T>
        inline constexpr bool is_mddisplacement = false;
        template <std::size_t N, class D>
        inline constexpr bool is_mddisplacement<wsl::mddisplacement<N, D>> = true;

        template <class T>
        concept position_like = is_mdposition<std::remove_cvref_t<T>>;
        template <class T>
        concept displacement_like = is_mddisplacement<std::remove_cvref_t<T>>;

        // Toute valeur positive de From est représentable dans To.
        template <class From, class To>
        inline constexpr bool lossless_positive =
            std::cmp_less_equal(std::numeric_limits<From>::max(), std::numeric_limits<To>::max());

        // Toute valeur de From est représentable dans To.
        template <class From, class To>
        inline constexpr bool lossless =
            std::cmp_less_equal(std::numeric_limits<From>::max(), std::numeric_limits<To>::max())
            && std::cmp_greater_equal(std::numeric_limits<From>::min(), std::numeric_limits<To>::min());

        // Valeur absolue de d, sans débordement (dans le type non signé).
        template <class D>
        constexpr std::make_unsigned_t<D> magnitude(D d) noexcept
        {
            using U = std::make_unsigned_t<D>;
            return d < 0 ? static_cast<U>(U{0} - static_cast<U>(d)) : static_cast<U>(d);
        }

        // Conversion en To, précondition : représentable.
        template <class To, class From>
        constexpr To narrow(From value, char const* message) noexcept
        {
            precondition(std::in_range<To>(value), message);
            return static_cast<To>(value);
        }

        // p + d, précondition : résultat dans [0, max(I)].
        template <class I, class D>
        constexpr I shift(I p, D d) noexcept
        {
            if (d >= 0)
            {
                precondition(std::cmp_less_equal(d, std::numeric_limits<I>::max() - p),
                             "position + déplacement : dépassement du type des composantes");
                return static_cast<I>(p + static_cast<I>(d));
            }
            auto const m = magnitude(d);
            precondition(std::cmp_greater_equal(p, m), "position + déplacement : composante négative");
            return static_cast<I>(p - static_cast<I>(m));
        }

        // p − d (sans former −d, qui peut déborder).
        template <class I, class D>
        constexpr I unshift(I p, D d) noexcept
        {
            if (d <= 0)
            {
                auto const m = magnitude(d);
                precondition(std::cmp_less_equal(m, std::numeric_limits<I>::max() - p),
                             "position − déplacement : dépassement du type des composantes");
                return static_cast<I>(p + static_cast<I>(m));
            }
            precondition(std::cmp_greater_equal(p, d), "position − déplacement : composante négative");
            return static_cast<I>(p - static_cast<I>(d));
        }

        // a − b pour deux positions, dans le type signé S.
        template <class S, class A, class B>
        constexpr S difference(A a, B b) noexcept
        {
            if (std::cmp_greater_equal(a, b))
            {
                using C = std::make_unsigned_t<std::common_type_t<A, B>>;
                C const m = static_cast<C>(static_cast<C>(a) - static_cast<C>(b));
                precondition(std::cmp_less_equal(m, std::numeric_limits<S>::max()),
                             "position − position : différence non représentable");
                return static_cast<S>(m);
            }
            using C = std::make_unsigned_t<std::common_type_t<A, B>>;
            C const m = static_cast<C>(static_cast<C>(b) - static_cast<C>(a));
            precondition(std::cmp_less_equal(m, magnitude(std::numeric_limits<S>::min())),
                         "position − position : différence non représentable");
            // −m, écrit pour que m == |min(S)| ne déborde pas.
            return static_cast<S>(-static_cast<S>(m - 1) - 1);
        }

        template <class D>
        constexpr D checked_add(D a, D b, char const* message) noexcept
        {
            D r{};
            precondition(!__builtin_add_overflow(a, b, &r), message);
            return r;
        }

        template <class D>
        constexpr D checked_sub(D a, D b, char const* message) noexcept
        {
            D r{};
            precondition(!__builtin_sub_overflow(a, b, &r), message);
            return r;
        }

        template <class D>
        constexpr D checked_mul(D a, D b, char const* message) noexcept
        {
            D r{};
            precondition(!__builtin_mul_overflow(a, b, &r), message);
            return r;
        }

        // Comparaison lexicographique de deux suites de valeurs entières
        // (valeurs mathématiques : jamais le piège -1 < 0u).
        template <std::size_t N, class A, class B>
        constexpr std::strong_ordering lexicographic(std::array<A, N> const& a, std::array<B, N> const& b) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
            {
                if (std::cmp_less(a[k], b[k]))
                    return std::strong_ordering::less;
                if (std::cmp_greater(a[k], b[k]))
                    return std::strong_ordering::greater;
            }
            return std::strong_ordering::equal;
        }

        inline constexpr std::size_t hash_mix(std::size_t seed, std::size_t value) noexcept
        {
            return seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
        }
    }

    // -----------------------------------------------------------------------
    // mddisplacement : un vecteur de ℤᴺ (déplacement entre deux positions)
    // -----------------------------------------------------------------------

    template <std::size_t N, class DiffType>
    class mddisplacement
    {
        static_assert(detail::mdarray::index_integer<DiffType> && std::is_signed_v<DiffType>,
                      "wsl::mddisplacement : DiffType doit être un type entier signé");

    public:
        using difference_type = DiffType;
        using value_type = DiffType;

        static constexpr std::size_t rank() noexcept { return N; }

        // Déplacement nul.
        constexpr mddisplacement() noexcept = default;

        // Exactement N composantes, de tout type entier (précondition : représentables).
        template <detail::mdarray::index_integer... V>
            requires(sizeof...(V) == N && N > 0)
        constexpr explicit(N == 1) mddisplacement(V... values) noexcept
            : values_{detail::mdarray::narrow<DiffType>(values, "mddisplacement : composante non représentable")...}
        {
        }

        template <detail::mdarray::index_integer V>
        constexpr mddisplacement(std::array<V, N> const& values) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = detail::mdarray::narrow<DiffType>(values[k], "mddisplacement : composante non représentable");
        }

        template <detail::mdarray::index_integer V>
        constexpr mddisplacement(std::span<V, N> values) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = detail::mdarray::narrow<DiffType>(values[k], "mddisplacement : composante non représentable");
        }

        template <detail::mdarray::index_integer V>
        constexpr explicit mddisplacement(std::span<V> values) noexcept
        {
            detail::mdarray::precondition(values.size() == N, "mddisplacement : nombre de composantes différent du rang");
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = detail::mdarray::narrow<DiffType>(values[k], "mddisplacement : composante non représentable");
        }

        // Autre type de composantes : implicite sans perte, explicite sinon (D28).
        template <class E>
        constexpr explicit(!detail::mdarray::lossless<E, DiffType>) mddisplacement(mddisplacement<N, E> const& other) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = detail::mdarray::narrow<DiffType>(other[k], "mddisplacement : composante non représentable");
        }

        // Toutes les composantes à value.
        template <detail::mdarray::index_integer V>
        static constexpr mddisplacement filled(V value) noexcept
        {
            mddisplacement d;
            DiffType const v = detail::mdarray::narrow<DiffType>(value, "mddisplacement::filled : valeur non représentable");
            for (auto& x : d.values_)
                x = v;
            return d;
        }

        // Vecteur de base e_k de ℤᴺ.
        static constexpr mddisplacement unit(std::size_t k) noexcept
        {
            detail::mdarray::precondition(k < N, "mddisplacement::unit : direction hors du rang");
            mddisplacement d;
            d.values_[k] = 1;
            return d;
        }

        constexpr DiffType operator[](std::size_t k) const noexcept
        {
            detail::mdarray::precondition(k < N, "mddisplacement : composante hors du rang");
            return values_[k];
        }

        constexpr std::array<DiffType, N> to_array() const noexcept { return values_; }

        template <std::size_t K>
        constexpr mddisplacement<K, DiffType> first() const noexcept
        {
            static_assert(K <= N, "mddisplacement::first : plus de composantes que le rang");
            std::array<DiffType, K> r{};
            for (std::size_t k = 0; k < K; ++k)
                r[k] = values_[k];
            return mddisplacement<K, DiffType>(r);
        }

        template <std::size_t K>
        constexpr mddisplacement<K, DiffType> last() const noexcept
        {
            static_assert(K <= N, "mddisplacement::last : plus de composantes que le rang");
            std::array<DiffType, K> r{};
            for (std::size_t k = 0; k < K; ++k)
                r[k] = values_[N - K + k];
            return mddisplacement<K, DiffType>(r);
        }

        template <std::size_t K>
        friend constexpr DiffType get(mddisplacement const& d) noexcept
        {
            static_assert(K < N, "get : composante hors du rang");
            return d.values_[K];
        }

        // ℤ-module.
        constexpr mddisplacement operator-() const noexcept
        {
            mddisplacement r;
            for (std::size_t k = 0; k < N; ++k)
                r.values_[k] = detail::mdarray::checked_sub(DiffType{0}, values_[k], "−déplacement : dépassement");
            return r;
        }

        template <class E>
        constexpr mddisplacement& operator+=(mddisplacement<N, E> const& e) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = detail::mdarray::checked_add(
                    values_[k], detail::mdarray::narrow<DiffType>(e[k], "déplacement += : composante non représentable"),
                    "déplacement += : dépassement");
            return *this;
        }

        template <class E>
        constexpr mddisplacement& operator-=(mddisplacement<N, E> const& e) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = detail::mdarray::checked_sub(
                    values_[k], detail::mdarray::narrow<DiffType>(e[k], "déplacement −= : composante non représentable"),
                    "déplacement −= : dépassement");
            return *this;
        }

        template <detail::mdarray::index_integer V>
        constexpr mddisplacement& operator*=(V factor) noexcept
        {
            DiffType const f = detail::mdarray::narrow<DiffType>(factor, "déplacement *= : facteur non représentable");
            for (auto& x : values_)
                x = detail::mdarray::checked_mul(x, f, "déplacement *= : dépassement");
            return *this;
        }

        template <detail::mdarray::index_integer V>
        constexpr mddisplacement& operator/=(V divisor) noexcept
        {
            DiffType const q = detail::mdarray::narrow<DiffType>(divisor, "déplacement /= : diviseur non représentable");
            detail::mdarray::precondition(q != 0, "déplacement /= : division par zéro");
            for (auto& x : values_)
            {
                detail::mdarray::precondition(!(q == -1 && x == std::numeric_limits<DiffType>::min()),
                                              "déplacement /= : dépassement");
                x = static_cast<DiffType>(x / q);
            }
            return *this;
        }

        template <detail::mdarray::index_integer V>
        constexpr mddisplacement& operator%=(V divisor) noexcept
        {
            DiffType const q = detail::mdarray::narrow<DiffType>(divisor, "déplacement %= : diviseur non représentable");
            detail::mdarray::precondition(q != 0, "déplacement %= : division par zéro");
            for (auto& x : values_)
                x = q == -1 ? DiffType{0} : static_cast<DiffType>(x % q);
            return *this;
        }

        template <class E>
        friend constexpr auto operator+(mddisplacement const& d, mddisplacement<N, E> const& e) noexcept
        {
            mddisplacement<N, std::common_type_t<DiffType, E>> r(d);
            r += e;
            return r;
        }

        template <class E>
        friend constexpr auto operator-(mddisplacement const& d, mddisplacement<N, E> const& e) noexcept
        {
            mddisplacement<N, std::common_type_t<DiffType, E>> r(d);
            r -= e;
            return r;
        }

        template <detail::mdarray::index_integer V>
        friend constexpr mddisplacement operator*(mddisplacement d, V factor) noexcept { return d *= factor; }
        template <detail::mdarray::index_integer V>
        friend constexpr mddisplacement operator*(V factor, mddisplacement d) noexcept { return d *= factor; }
        template <detail::mdarray::index_integer V>
        friend constexpr mddisplacement operator/(mddisplacement d, V divisor) noexcept { return d /= divisor; }
        template <detail::mdarray::index_integer V>
        friend constexpr mddisplacement operator%(mddisplacement d, V divisor) noexcept { return d %= divisor; }

        template <class E>
        friend constexpr bool operator==(mddisplacement const& d, mddisplacement<N, E> const& e) noexcept
        {
            return detail::mdarray::lexicographic(d.values_, e.to_array()) == std::strong_ordering::equal;
        }

        template <class E>
        friend constexpr std::strong_ordering operator<=>(mddisplacement const& d, mddisplacement<N, E> const& e) noexcept
        {
            return detail::mdarray::lexicographic(d.values_, e.to_array());
        }

    private:
        std::array<DiffType, N> values_{};
    };

    template <detail::mdarray::index_integer... V>
    mddisplacement(V...) -> mddisplacement<sizeof...(V), std::ptrdiff_t>;

    template <detail::mdarray::index_integer V, std::size_t N>
    mddisplacement(std::array<V, N>) -> mddisplacement<N, std::ptrdiff_t>;

    // -----------------------------------------------------------------------
    // mdposition : une position dans un tableau (un point de ℕᴺ ⊂ ℤᴺ)
    // -----------------------------------------------------------------------

    template <std::size_t N, class IndexType>
    class mdposition
    {
        static_assert(detail::mdarray::index_integer<IndexType>,
                      "wsl::mdposition : IndexType doit être un type entier (ni bool ni caractère)");

    public:
        using index_type = IndexType;
        using value_type = IndexType;
        using displacement_type = mddisplacement<N, std::make_signed_t<IndexType>>;

        static constexpr std::size_t rank() noexcept { return N; }

        // L'origine.
        constexpr mdposition() noexcept = default;

        // Exactement N composantes, de tout type entier, comme std::extents
        // (précondition : positives et représentables).
        template <detail::mdarray::index_integer... V>
            requires(sizeof...(V) == N && N > 0)
        constexpr explicit(N == 1) mdposition(V... values) noexcept
            : values_{component(values)...}
        {
        }

        template <detail::mdarray::index_integer V>
        constexpr mdposition(std::array<V, N> const& values) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = component(values[k]);
        }

        template <detail::mdarray::index_integer V>
        constexpr mdposition(std::span<V, N> values) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = component(values[k]);
        }

        template <detail::mdarray::index_integer V>
        constexpr explicit mdposition(std::span<V> values) noexcept
        {
            detail::mdarray::precondition(values.size() == N, "mdposition : nombre de composantes différent du rang");
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = component(values[k]);
        }

        // Autre type de composantes : implicite sans perte, explicite sinon (D28).
        template <class J>
        constexpr explicit(!detail::mdarray::lossless_positive<J, IndexType>) mdposition(mdposition<N, J> const& other) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = component(other[k]);
        }

        // Toutes les composantes à value.
        template <detail::mdarray::index_integer V>
        static constexpr mdposition filled(V value) noexcept
        {
            mdposition p;
            IndexType const v = component(value);
            for (auto& x : p.values_)
                x = v;
            return p;
        }

        constexpr IndexType operator[](std::size_t k) const noexcept
        {
            detail::mdarray::precondition(k < N, "mdposition : composante hors du rang");
            return values_[k];
        }

        constexpr std::array<IndexType, N> to_array() const noexcept { return values_; }

        template <std::size_t K>
        constexpr mdposition<K, IndexType> first() const noexcept
        {
            static_assert(K <= N, "mdposition::first : plus de composantes que le rang");
            std::array<IndexType, K> r{};
            for (std::size_t k = 0; k < K; ++k)
                r[k] = values_[k];
            return mdposition<K, IndexType>(r);
        }

        template <std::size_t K>
        constexpr mdposition<K, IndexType> last() const noexcept
        {
            static_assert(K <= N, "mdposition::last : plus de composantes que le rang");
            std::array<IndexType, K> r{};
            for (std::size_t k = 0; k < K; ++k)
                r[k] = values_[N - K + k];
            return mdposition<K, IndexType>(r);
        }

        template <std::size_t K>
        friend constexpr IndexType get(mdposition const& p) noexcept
        {
            static_assert(K < N, "get : composante hors du rang");
            return p.values_[K];
        }

        // Action des déplacements.
        template <class D>
        constexpr mdposition& operator+=(mddisplacement<N, D> const& d) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = detail::mdarray::shift(values_[k], d[k]);
            return *this;
        }

        template <class D>
        constexpr mdposition& operator-=(mddisplacement<N, D> const& d) noexcept
        {
            for (std::size_t k = 0; k < N; ++k)
                values_[k] = detail::mdarray::unshift(values_[k], d[k]);
            return *this;
        }

        template <class D>
        friend constexpr mdposition operator+(mdposition p, mddisplacement<N, D> const& d) noexcept { return p += d; }
        template <class D>
        friend constexpr mdposition operator+(mddisplacement<N, D> const& d, mdposition p) noexcept { return p += d; }
        template <class D>
        friend constexpr mdposition operator-(mdposition p, mddisplacement<N, D> const& d) noexcept { return p -= d; }

        // position − position = déplacement, dans le type signé commun.
        template <class J>
        friend constexpr auto operator-(mdposition const& p, mdposition<N, J> const& q) noexcept
        {
            using S = std::make_signed_t<std::common_type_t<IndexType, J>>;
            std::array<S, N> r{};
            for (std::size_t k = 0; k < N; ++k)
                r[k] = detail::mdarray::difference<S>(p.values_[k], q[k]);
            return mddisplacement<N, S>(r);
        }

        template <class J>
        friend constexpr bool operator==(mdposition const& p, mdposition<N, J> const& q) noexcept
        {
            return detail::mdarray::lexicographic(p.values_, q.to_array()) == std::strong_ordering::equal;
        }

        template <class J>
        friend constexpr std::strong_ordering operator<=>(mdposition const& p, mdposition<N, J> const& q) noexcept
        {
            return detail::mdarray::lexicographic(p.values_, q.to_array());
        }

    private:
        template <class V>
        static constexpr IndexType component(V value) noexcept
        {
            detail::mdarray::precondition(std::cmp_greater_equal(value, 0), "mdposition : composante négative");
            return detail::mdarray::narrow<IndexType>(value, "mdposition : composante non représentable");
        }

        std::array<IndexType, N> values_{};
    };

    template <detail::mdarray::index_integer... V>
    mdposition(V...) -> mdposition<sizeof...(V), std::size_t>;

    template <detail::mdarray::index_integer V, std::size_t N>
    mdposition(std::array<V, N>) -> mdposition<N, std::size_t>;

    // -----------------------------------------------------------------------
    // Fonctions (trouvées par ADL)
    // -----------------------------------------------------------------------

    // Ordre composante par composante (ordre partiel, D18).
    template <class A, class B>
        requires((detail::mdarray::position_like<A> && detail::mdarray::position_like<B>)
                 || (detail::mdarray::displacement_like<A> && detail::mdarray::displacement_like<B>))
                && (std::remove_cvref_t<A>::rank() == std::remove_cvref_t<B>::rank())
    constexpr bool all_less(A const& q, B const& p) noexcept
    {
        for (std::size_t k = 0; k < A::rank(); ++k)
            if (!std::cmp_less(q[k], p[k]))
                return false;
        return true;
    }

    template <class A, class B>
        requires((detail::mdarray::position_like<A> && detail::mdarray::position_like<B>)
                 || (detail::mdarray::displacement_like<A> && detail::mdarray::displacement_like<B>))
                && (std::remove_cvref_t<A>::rank() == std::remove_cvref_t<B>::rank())
    constexpr bool all_less_equal(A const& q, B const& p) noexcept
    {
        for (std::size_t k = 0; k < A::rank(); ++k)
            if (!std::cmp_less_equal(q[k], p[k]))
                return false;
        return true;
    }

    // Assemblage (D24) : les composantes de chacun, dans l'ordre.
    template <std::size_t... N, class... I>
    constexpr auto concat(mdposition<N, I> const&... p) noexcept
    {
        using T = std::common_type_t<I...>;
        std::array<T, (N + ...)> r{};
        std::size_t k = 0;
        ((void)[&] {
            for (std::size_t j = 0; j < N; ++j)
                r[k++] = static_cast<T>(p[j]);
        }(), ...);
        return mdposition<(N + ...), T>(r);
    }

    template <std::size_t... N, class... D>
    constexpr auto concat(mddisplacement<N, D> const&... d) noexcept
    {
        using T = std::common_type_t<D...>;
        std::array<T, (N + ...)> r{};
        std::size_t k = 0;
        ((void)[&] {
            for (std::size_t j = 0; j < N; ++j)
                r[k++] = static_cast<T>(d[j]);
        }(), ...);
        return mddisplacement<(N + ...), T>(r);
    }

    // Normes d'un déplacement (D26), dans DiffType (précondition : représentables).
    template <std::size_t N, class D>
    constexpr D norm1(mddisplacement<N, D> const& d) noexcept
    {
        D s{0};
        for (std::size_t k = 0; k < N; ++k)
            s = detail::mdarray::checked_add(s, detail::mdarray::narrow<D>(detail::mdarray::magnitude(d[k]), "norm1 : dépassement"),
                                             "norm1 : dépassement");
        return s;
    }

    template <std::size_t N, class D>
    constexpr D norm_inf(mddisplacement<N, D> const& d) noexcept
    {
        D m{0};
        for (std::size_t k = 0; k < N; ++k)
        {
            D const a = detail::mdarray::narrow<D>(detail::mdarray::magnitude(d[k]), "norm_inf : dépassement");
            m = a > m ? a : m;
        }
        return m;
    }

    template <std::size_t N, class D>
    constexpr D norm2_squared(mddisplacement<N, D> const& d) noexcept
    {
        D s{0};
        for (std::size_t k = 0; k < N; ++k)
            s = detail::mdarray::checked_add(s, detail::mdarray::checked_mul(d[k], d[k], "norm2_squared : dépassement"),
                                             "norm2_squared : dépassement");
        return s;
    }

    // Division euclidienne composante par composante (D15), k > 0.
    template <std::size_t N, class D, detail::mdarray::index_integer V>
    constexpr mddisplacement<N, D> div_euclid(mddisplacement<N, D> const& d, V k) noexcept
    {
        D const q = detail::mdarray::narrow<D>(k, "div_euclid : diviseur non représentable");
        std::array<D, N> r{};
        for (std::size_t j = 0; j < N; ++j)
            r[j] = div_euclid(d[j], q);
        return mddisplacement<N, D>(r);
    }

    template <std::size_t N, class D, detail::mdarray::index_integer V>
    constexpr mddisplacement<N, D> rem_euclid(mddisplacement<N, D> const& d, V k) noexcept
    {
        D const q = detail::mdarray::narrow<D>(k, "rem_euclid : diviseur non représentable");
        std::array<D, N> r{};
        for (std::size_t j = 0; j < N; ++j)
            r[j] = rem_euclid(d[j], q);
        return mddisplacement<N, D>(r);
    }

    // Affichage « (2, 1) » par std::format (D27).
    template <std::size_t N, class I>
    std::ostream& operator<<(std::ostream& out, mdposition<N, I> const& p);
    template <std::size_t N, class D>
    std::ostream& operator<<(std::ostream& out, mddisplacement<N, D> const& d);
}

// Protocole des tuples : auto [i, j] = p;
template <std::size_t N, class I>
struct std::tuple_size<wsl::mdposition<N, I>> : std::integral_constant<std::size_t, N>
{
};
template <std::size_t K, std::size_t N, class I>
struct std::tuple_element<K, wsl::mdposition<N, I>>
{
    using type = I;
};
template <std::size_t N, class D>
struct std::tuple_size<wsl::mddisplacement<N, D>> : std::integral_constant<std::size_t, N>
{
};
template <std::size_t K, std::size_t N, class D>
struct std::tuple_element<K, wsl::mddisplacement<N, D>>
{
    using type = D;
};

// Hachage (D27) : positions et déplacements dans std::unordered_map, std::unordered_set.
template <std::size_t N, class I>
struct std::hash<wsl::mdposition<N, I>>
{
    constexpr std::size_t operator()(wsl::mdposition<N, I> const& p) const noexcept
    {
        std::size_t seed = N;
        for (std::size_t k = 0; k < N; ++k)
            seed = wsl::detail::mdarray::hash_mix(seed, static_cast<std::size_t>(p[k]));
        return seed;
    }
};
template <std::size_t N, class D>
struct std::hash<wsl::mddisplacement<N, D>>
{
    constexpr std::size_t operator()(wsl::mddisplacement<N, D> const& d) const noexcept
    {
        std::size_t seed = N;
        for (std::size_t k = 0; k < N; ++k)
            seed = wsl::detail::mdarray::hash_mix(seed, static_cast<std::size_t>(d[k]));
        return seed;
    }
};

// Formatage (D27) : « (2, 1) », les options s'appliquent à chaque composante.
namespace wsl::detail::mdarray
{
    template <class Tuple, class T>
    struct tuple_formatter
    {
        std::formatter<T, char> element;

        constexpr auto parse(std::format_parse_context& context) { return element.parse(context); }

        template <class Context>
        auto format(Tuple const& t, Context& context) const
        {
            auto out = context.out();
            *out++ = '(';
            for (std::size_t k = 0; k < Tuple::rank(); ++k)
            {
                if (k > 0)
                {
                    *out++ = ',';
                    *out++ = ' ';
                }
                context.advance_to(out);
                out = element.format(t[k], context);
            }
            *out++ = ')';
            return out;
        }
    };
}

template <std::size_t N, class I>
struct std::formatter<wsl::mdposition<N, I>, char> : wsl::detail::mdarray::tuple_formatter<wsl::mdposition<N, I>, I>
{
};
template <std::size_t N, class D>
struct std::formatter<wsl::mddisplacement<N, D>, char> : wsl::detail::mdarray::tuple_formatter<wsl::mddisplacement<N, D>, D>
{
};

template <std::size_t N, class I>
std::ostream& wsl::operator<<(std::ostream& out, wsl::mdposition<N, I> const& p)
{
    return out << std::format("{}", p);
}
template <std::size_t N, class D>
std::ostream& wsl::operator<<(std::ostream& out, wsl::mddisplacement<N, D> const& d)
{
    return out << std::format("{}", d);
}
