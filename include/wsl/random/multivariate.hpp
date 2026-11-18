#pragma once

// Socle des lois multivariées (voir docs/random.md §12) : dimension fixe ou dynamique,
// sur le modèle de std::span, vecteurs et matrices carrées.
//
//   wsl::dirichlet<double, 3>                // valeurs std::array<double, 3>
//   wsl::dirichlet<double, wsl::dynamic>     // valeurs std::vector<double>, dimension donnée par les paramètres
//   wsl::square_matrix<double, 2> S{{1, 0.5}, {0.5, 2}};   // ligne par ligne
//
// Calcul (detail) : produits scalaires et sommes en double-double, produits par fmul ;
// factorisation de Cholesky (Σ = L Lᵀ, L triangulaire inférieure) et résolution
// triangulaire, bit-exactes.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/special.hpp>

#include <array>
#include <cstddef>
#include <initializer_list>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace wsl
{
    inline constexpr std::size_t dynamic = std::dynamic_extent;

    // Vecteur de dimension N (std::array) ou dynamique (std::vector).
    template <class T, std::size_t N>
    using vector_of = std::conditional_t<N == dynamic, std::vector<T>, std::array<T, N>>;

    namespace detail
    {
        template <class T, std::size_t N>
        constexpr vector_of<T, N> make_vector(std::size_t n, T value = T{})
        {
            if constexpr (N == dynamic)
                return std::vector<T>(n, value);
            else
            {
                std::array<T, N> v{};
                v.fill(value);
                return v;
            }
        }
    }

    // Matrice carrée n × n, ligne par ligne ; n = N, ou donné à la construction.
    template <class T, std::size_t N = dynamic>
    class square_matrix
    {
    public:
        using value_type = T;

        constexpr square_matrix() requires(N != dynamic) = default;
        constexpr explicit square_matrix(std::size_t n) requires(N == dynamic) : n_(n), data_(n * n) {}

        // {{a₀₀, a₀₁, …}, {a₁₀, …}, …}
        constexpr square_matrix(std::initializer_list<std::initializer_list<T>> rows) : n_(rows.size())
        {
            if constexpr (N != dynamic)
                if (rows.size() != N)
                    throw std::domain_error("wsl::square_matrix : nombre de lignes différent de la dimension");
            if constexpr (N == dynamic)
                data_.resize(n_ * n_);
            std::size_t i = 0;
            for (auto const& row : rows)
            {
                if (row.size() != n_)
                    throw std::domain_error("wsl::square_matrix : la matrice doit être carrée");
                std::size_t j = 0;
                for (T const& v : row)
                    data_[i * n_ + j++] = v;
                ++i;
            }
        }

        constexpr std::size_t size() const noexcept { return N == dynamic ? n_ : N; }
        constexpr T& operator()(std::size_t i, std::size_t j) noexcept { return data_[i * size() + j]; }
        constexpr T const& operator()(std::size_t i, std::size_t j) const noexcept { return data_[i * size() + j]; }

        friend constexpr bool operator==(square_matrix const& a, square_matrix const& b) noexcept { return a.size() == b.size() && a.data_ == b.data_; }

    private:
        std::size_t n_ = N == dynamic ? 0 : N;
        std::conditional_t<N == dynamic, std::vector<T>, std::array<T, N == dynamic ? 1 : N * N>> data_{};
    };

    namespace detail
    {
        // Σᵢ aᵢ bᵢ en double-double.
        template <class F>
        constexpr sd::double_double dot(std::size_t n, F term) noexcept
        {
            sd::double_double s{0.0, 0.0};
            for (std::size_t i = 0; i < n; ++i)
            {
                auto const [a, b] = term(i);
                s = sum(s, product(a, b));
            }
            return s;
        }

        // Σ = L Lᵀ, L triangulaire inférieure (diagonale > 0) ; faux si Σ n'est pas
        // symétrique définie positive. Sommes en double-double, un arrondi par coefficient.
        template <class T, std::size_t N>
        constexpr bool cholesky(square_matrix<T, N> const& s, square_matrix<double, N>& l)
        {
            std::size_t const n = s.size();
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    if (!(double(s(i, j)) == double(s(j, i))) || !is_finite(double(s(i, j))))
                        return false;
            for (std::size_t j = 0; j < n; ++j)
            {
                auto d = dot(j, [&](std::size_t k) { return std::pair{l(j, k), l(j, k)}; });
                auto const r = sum({double(s(j, j)), 0.0}, {-d.hi, -d.lo});
                double const v = r.hi + r.lo;
                if (!(v > 0))
                    return false;
                double const ljj = wsl::special::sqrt(v);
                l(j, j) = ljj;
                for (std::size_t i = j + 1; i < n; ++i)
                {
                    auto const e = dot(j, [&](std::size_t k) { return std::pair{l(i, k), l(j, k)}; });
                    auto const t = sum({double(s(i, j)), 0.0}, {-e.hi, -e.lo});
                    l(i, j) = (t.hi + t.lo) / ljj;
                }
                for (std::size_t i = 0; i < j; ++i)
                    l(i, j) = 0;
            }
            return true;
        }

        // z tel que L z = b (L triangulaire inférieure).
        template <std::size_t N, class V>
        constexpr V forward_substitution(square_matrix<double, N> const& l, V b)
        {
            std::size_t const n = l.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                auto const e = dot(i, [&](std::size_t k) { return std::pair{l(i, k), double(b[k])}; });
                auto const t = sum({double(b[i]), 0.0}, {-e.hi, -e.lo});
                b[i] = static_cast<typename V::value_type>((t.hi + t.lo) / l(i, i));
            }
            return b;
        }
    }
}
