#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Socle des lois multivariées (voir design/random.md §12) : dimension fixe ou dynamique,
// sur le modèle de std::span, vecteurs et matrices carrées.
//
//   wsl::probability::laws::dirichlet<double, 3>                // valeurs std::array<double, 3>
//   wsl::probability::laws::dirichlet<double, wsl::probability::dynamic>     // valeurs std::vector<double>, dimension donnée par les paramètres
//   wsl::probability::square_matrix<double, 2> S{{1, 0.5}, {0.5, 2}};   // ligne par ligne
//
// Calcul (detail) : produits scalaires et sommes en double-double, produits par fmul ;
// factorisation de Cholesky (Σ = L Lᵀ, L triangulaire inférieure) et résolution
// triangulaire, bit-exactes.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/loader.hpp>
#include <wsl/random/detail/portable_math.hpp>

#include <array>
#include <cstddef>
#include <initializer_list>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace wsl::inline WSL_DEFAULTS::probability
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
                    throw std::domain_error("wsl::probability::square_matrix : nombre de lignes différent de la dimension");
            if constexpr (N == dynamic)
                data_.resize(n_ * n_);
            std::size_t i = 0;
            for (auto const& row : rows)
            {
                if (row.size() != n_)
                    throw std::domain_error("wsl::probability::square_matrix : la matrice doit être carrée");
                std::size_t j = 0;
                for (T const& v : row)
                    data_[i * n_ + j++] = v;
                ++i;
            }
        }

        // Dimensions comme std::mdspan : extent(0) = extent(1) = n ; accès m[i, j] (par lignes).
        static constexpr std::size_t rank() noexcept { return 2; }
        constexpr std::size_t extent(std::size_t) const noexcept { return N == dynamic ? n_ : N; }
        constexpr T& operator[](std::size_t i, std::size_t j) noexcept { return data_[i * extent(0) + j]; }
        constexpr T const& operator[](std::size_t i, std::size_t j) const noexcept { return data_[i * extent(0) + j]; }

        friend constexpr bool operator==(square_matrix const& a, square_matrix const& b) noexcept { return a.extent(0) == b.extent(0) && a.data_ == b.data_; }

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
            std::size_t const n = s.extent(0);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    if (!(double(s[i, j]) == double(s[j, i])) || !is_finite(double(s[i, j])))
                        return false;
            for (std::size_t j = 0; j < n; ++j)
            {
                auto d = dot(j, [&](std::size_t k) { return std::pair{l[j, k], l[j, k]}; });
                auto const r = sum({double(s[j, j]), 0.0}, {-d.hi, -d.lo});
                double const v = r.hi + r.lo;
                if (!(v > 0))
                    return false;
                double const ljj = wsl::random::detail::pm::sqrt(v);
                l[j, j] = ljj;
                for (std::size_t i = j + 1; i < n; ++i)
                {
                    auto const e = dot(j, [&](std::size_t k) { return std::pair{l[i, k], l[j, k]}; });
                    auto const t = sum({double(s[i, j]), 0.0}, {-e.hi, -e.lo});
                    l[i, j] = (t.hi + t.lo) / ljj;
                }
                for (std::size_t i = 0; i < j; ++i)
                    l[i, j] = 0;
            }
            return true;
        }

        template <class T, std::size_t N>
        constexpr square_matrix<T, N> make_matrix(std::size_t n)
        {
            if constexpr (N == dynamic)
                return square_matrix<T, N>(n);
            else
                return square_matrix<T, N>();
        }

        // Σ log Lᵢᵢ (½ log det de L Lᵀ), en double-double.
        template <std::size_t N>
        constexpr sd::double_double half_log_determinant(square_matrix<double, N> const& l) noexcept
        {
            sd::double_double s{0.0, 0.0};
            for (std::size_t i = 0; i < l.extent(0); ++i)
                s = sum(s, log_precise(l[i, i]));
            return s;
        }

        // A B, A et B quelconques ; chaque coefficient arrondi une fois.
        template <std::size_t N>
        constexpr square_matrix<double, N> multiply(square_matrix<double, N> const& a, square_matrix<double, N> const& b)
        {
            std::size_t const n = a.extent(0);
            auto c = make_matrix<double, N>(n);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                {
                    auto const t = dot(n, [&](std::size_t k) { return std::pair{a[i, k], b[k, j]}; });
                    c[i, j] = t.hi + t.lo;
                }
            return c;
        }

        // B Bᵀ (symétrique exactement : (i, j) calculé une fois).
        template <class T, std::size_t N>
        constexpr square_matrix<T, N> gram(square_matrix<double, N> const& b)
        {
            std::size_t const n = b.extent(0);
            auto c = make_matrix<T, N>(n);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j <= i; ++j)
                {
                    auto const t = dot(n, [&](std::size_t k) { return std::pair{b[i, k], b[j, k]}; });
                    c[i, j] = c[j, i] = static_cast<T>(t.hi + t.lo);
                }
            return c;
        }

        // A⁻¹, A triangulaire inférieure inversible (colonne par colonne).
        template <std::size_t N>
        constexpr square_matrix<double, N> lower_inverse(square_matrix<double, N> const& a)
        {
            std::size_t const n = a.extent(0);
            auto r = make_matrix<double, N>(n);
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t i = j; i < n; ++i)
                {
                    auto const e = dot(i - j, [&](std::size_t k) { return std::pair{a[i, j + k], r[j + k, j]}; });
                    auto const t = sum({i == j ? 1.0 : 0.0, 0.0}, {-e.hi, -e.lo});
                    r[i, j] = (t.hi + t.lo) / a[i, i];
                }
            return r;
        }

        // ‖L⁻¹ M‖²_F, L triangulaire inférieure : tr(Mᵀ (L Lᵀ)⁻¹ M), en double-double.
        template <std::size_t N>
        constexpr sd::double_double trace_solve(square_matrix<double, N> const& l, square_matrix<double, N> const& m)
        {
            std::size_t const n = l.extent(0);
            sd::double_double s{0.0, 0.0};
            auto column = make_vector<double, N>(n);
            for (std::size_t j = 0; j < n; ++j)
            {
                for (std::size_t i = 0; i < n; ++i)
                {
                    auto const e = dot(i, [&](std::size_t k) { return std::pair{l[i, k], column[k]}; });
                    auto const t = sum({m[i, j], 0.0}, {-e.hi, -e.lo});
                    column[i] = (t.hi + t.lo) / l[i, i];
                }
                s = sum(s, dot(n, [&](std::size_t k) { return std::pair{column[k], column[k]}; }));
            }
            return s;
        }

        // log Γ_p(a) = p(p − 1)/4 log π + Σⱼ lgamma(a − j/2), a > (p − 1)/2, en double-double.
        constexpr sd::double_double log_multivariate_gamma(std::size_t p, double a) noexcept
        {
            double const c = sd::fmul(double(p) * double(p - 1), 0.25);    // exact
            sd::double_double s = mul({c, 0.0}, {sd::log_pi_hi, sd::log_pi_lo});
            for (std::size_t j = 0; j < p; ++j)
                s = sum(s, log_gamma(a - sd::fmul(double(j), 0.5)));
            return s;
        }
        // ψ_p(a) = Σⱼ ψ(a − j/2).
        constexpr double multivariate_digamma(std::size_t p, double a) noexcept
        {
            sd::double_double s{0.0, 0.0};
            for (std::size_t j = 0; j < p; ++j)
                s = sum(s, {wsl::random::detail::pm::digamma(a - sd::fmul(double(j), 0.5)), 0.0});
            return s.hi + s.lo;
        }

        // z tel que L z = b (L triangulaire inférieure).
        template <std::size_t N, class V>
        constexpr V forward_substitution(square_matrix<double, N> const& l, V b)
        {
            std::size_t const n = l.extent(0);
            for (std::size_t i = 0; i < n; ++i)
            {
                auto const e = dot(i, [&](std::size_t k) { return std::pair{l[i, k], double(b[k])}; });
                auto const t = sum({double(b[i]), 0.0}, {-e.hi, -e.lo});
                b[i] = static_cast<typename V::value_type>((t.hi + t.lo) / l[i, i]);
            }
            return b;
        }
    }
}
