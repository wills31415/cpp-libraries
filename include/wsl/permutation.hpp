#pragma once

// Permutations de {0, …, n − 1}, constexpr.
//
//   wsl::permutation s{2, 0, 1};              // s(0) = 2, s(1) = 0, s(2) = 1
//   wsl::permutation t = wsl::permutation::from_cycles(3, {{0, 1}});
//   (s * t)(0);  s.inverse();  s.cycles();  s.sign();  s.inversions();
//
// Convention : (s * t)(i) = s(t(i)) (composition des applications, t d'abord).
// permute(s, v) déplace l'élément d'indice i en position s(i) : r[s(i)] = v[i].
// Les constructeurs vérifient que les images forment une permutation
// (std::invalid_argument sinon) ; composer des permutations de tailles différentes
// lève std::invalid_argument.

#include <algorithm>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wsl
{
    class permutation
    {
    public:
        using size_type = std::size_t;
        using const_iterator = std::vector<size_type>::const_iterator;

        constexpr permutation() = default;
        // Identité de taille n.
        constexpr explicit permutation(size_type n) : images_(n)
        {
            for (size_type i = 0; i < n; ++i)
                images_[i] = i;
        }
        // Images s(0), …, s(n − 1).
        constexpr explicit permutation(std::vector<size_type> images) : images_(std::move(images)) { check(); }
        constexpr permutation(std::initializer_list<size_type> images) : images_(images) { check(); }

        // Produit de cycles disjoints (les éléments absents sont fixes).
        static constexpr permutation from_cycles(size_type n, std::vector<std::vector<size_type>> const& cycles)
        {
            permutation p(n);
            std::vector<bool> seen(n, false);
            for (auto const& c : cycles)
                for (size_type k = 0; k < c.size(); ++k)
                {
                    size_type const a = c[k], b = c[(k + 1) % c.size()];
                    if (a >= n || b >= n || seen[a])
                        throw std::invalid_argument("wsl::permutation : cycles non disjoints ou hors de {0, …, n − 1}");
                    seen[a] = true;
                    p.images_[a] = b;
                }
            return p;
        }

        constexpr size_type size() const noexcept { return images_.size(); }
        constexpr bool empty() const noexcept { return images_.empty(); }
        constexpr size_type operator()(size_type i) const noexcept { return images_[i]; }
        constexpr size_type operator[](size_type i) const noexcept { return images_[i]; }
        constexpr size_type at(size_type i) const
        {
            if (i >= size())
                throw std::out_of_range("wsl::permutation::at");
            return images_[i];
        }
        constexpr std::span<size_type const> images() const noexcept { return images_; }
        constexpr const_iterator begin() const noexcept { return images_.begin(); }
        constexpr const_iterator end() const noexcept { return images_.end(); }

        // (s * t)(i) = s(t(i))
        friend constexpr permutation operator*(permutation const& s, permutation const& t)
        {
            if (s.size() != t.size())
                throw std::invalid_argument("wsl::permutation : composition de tailles différentes");
            permutation r;
            r.images_.resize(s.size());
            for (size_type i = 0; i < s.size(); ++i)
                r.images_[i] = s.images_[t.images_[i]];
            return r;
        }
        constexpr permutation& operator*=(permutation const& t) { return *this = *this * t; }

        constexpr permutation inverse() const
        {
            permutation r;
            r.images_.resize(size());
            for (size_type i = 0; i < size(); ++i)
                r.images_[images_[i]] = i;
            return r;
        }

        // sᵏ, k entier relatif (exponentiation rapide).
        constexpr permutation power(long long k) const
        {
            permutation base = k < 0 ? inverse() : *this, r(size());
            unsigned long long e = k < 0 ? 0ULL - static_cast<unsigned long long>(k) : static_cast<unsigned long long>(k);
            while (e != 0)
            {
                if (e & 1)
                    r = r * base;
                base = base * base;
                e >>= 1;
            }
            return r;
        }

        constexpr bool is_identity() const noexcept
        {
            for (size_type i = 0; i < size(); ++i)
                if (images_[i] != i)
                    return false;
            return true;
        }
        constexpr size_type fixed_points() const noexcept
        {
            size_type c = 0;
            for (size_type i = 0; i < size(); ++i)
                c += images_[i] == i;
            return c;
        }

        // Cycles, chacun commençant par son plus petit élément, rangés par celui-ci
        // (points fixes compris).
        constexpr std::vector<std::vector<size_type>> cycles() const
        {
            std::vector<std::vector<size_type>> r;
            std::vector<bool> seen(size(), false);
            for (size_type i = 0; i < size(); ++i)
            {
                if (seen[i])
                    continue;
                std::vector<size_type> c;
                for (size_type j = i; !seen[j]; j = images_[j])
                {
                    seen[j] = true;
                    c.push_back(j);
                }
                r.push_back(std::move(c));
            }
            return r;
        }
        constexpr size_type cycle_count() const
        {
            size_type c = 0;
            std::vector<bool> seen(size(), false);
            for (size_type i = 0; i < size(); ++i)
                if (!seen[i])
                {
                    ++c;
                    for (size_type j = i; !seen[j]; j = images_[j])
                        seen[j] = true;
                }
            return c;
        }
        // Longueurs des cycles, décroissantes (une partition de n).
        constexpr std::vector<size_type> cycle_type() const
        {
            std::vector<size_type> r;
            for (auto const& c : cycles())
                r.push_back(c.size());
            std::ranges::sort(r, std::greater<>{});
            return r;
        }
        // (−1)^(n − nombre de cycles)
        constexpr int sign() const { return (size() - cycle_count()) % 2 == 0 ? 1 : -1; }

        // Nombre de paires i < j avec s(i) > s(j), en O(n log n) (tri fusion).
        constexpr std::uint64_t inversions() const
        {
            std::vector<size_type> a = images_, buffer(size());
            std::uint64_t count = 0;
            for (size_type width = 1; width < size(); width *= 2)
            {
                for (size_type lo = 0; lo < size(); lo += 2 * width)
                {
                    size_type const mid = std::min(lo + width, size()), hi = std::min(lo + 2 * width, size());
                    size_type i = lo, j = mid, k = lo;
                    while (i < mid && j < hi)
                    {
                        if (a[j] < a[i])
                        {
                            count += mid - i;
                            buffer[k++] = a[j++];
                        }
                        else
                            buffer[k++] = a[i++];
                    }
                    while (i < mid)
                        buffer[k++] = a[i++];
                    while (j < hi)
                        buffer[k++] = a[j++];
                }
                std::swap(a, buffer);
            }
            return count;
        }

        // Rang dans l'ordre lexicographique des images (0 pour l'identité), n <= 20.
        constexpr std::uint64_t lexicographic_rank() const
        {
            if (size() > 20)
                throw std::out_of_range("wsl::permutation::lexicographic_rank : n > 20");
            std::uint64_t rank = 0;
            for (size_type i = 0; i < size(); ++i)
            {
                std::uint64_t smaller = 0;
                for (size_type j = i + 1; j < size(); ++j)
                    smaller += images_[j] < images_[i];
                rank = rank * (size() - i) + smaller;
            }
            return rank;
        }
        static constexpr permutation from_lexicographic_rank(size_type n, std::uint64_t rank)
        {
            if (n > 20)
                throw std::out_of_range("wsl::permutation::from_lexicographic_rank : n > 20");
            std::uint64_t factorial = 1;
            for (size_type i = 2; i <= n; ++i)
                factorial *= i;
            if (n > 0 && rank >= factorial)
                throw std::out_of_range("wsl::permutation::from_lexicographic_rank : rang >= n!");
            std::vector<size_type> free(n), images;
            for (size_type i = 0; i < n; ++i)
                free[i] = i;
            for (size_type i = n; i > 0; --i)
            {
                factorial /= i;
                size_type const k = static_cast<size_type>(rank / factorial);
                rank %= factorial;
                images.push_back(free[k]);
                free.erase(free.begin() + static_cast<std::ptrdiff_t>(k));
            }
            permutation p;
            p.images_ = std::move(images);
            return p;
        }

        // Permutation suivante dans l'ordre lexicographique (faux après la dernière, comme
        // std::next_permutation, qui revient alors à l'identité).
        constexpr bool next() { return std::next_permutation(images_.begin(), images_.end()); }

        friend constexpr bool operator==(permutation const&, permutation const&) = default;
        friend constexpr auto operator<=>(permutation const& a, permutation const& b) { return a.images_ <=> b.images_; }

    private:
        constexpr void check() const
        {
            std::vector<bool> seen(size(), false);
            for (size_type v : images_)
            {
                if (v >= size() || seen[v])
                    throw std::invalid_argument("wsl::permutation : les images doivent être 0, …, n − 1, chacune une fois");
                seen[v] = true;
            }
        }

        std::vector<size_type> images_;
    };

    // r[s(i)] = v[i] : l'élément d'indice i va en position s(i).
    template <class T>
    constexpr std::vector<T> permute(permutation const& s, std::span<T const> v)
    {
        if (v.size() != s.size())
            throw std::invalid_argument("wsl::permute : taille différente de celle de la permutation");
        std::vector<T> r(v.size());
        for (std::size_t i = 0; i < v.size(); ++i)
            r[s(i)] = v[i];
        return r;
    }
    template <class T>
    constexpr std::vector<T> permute(permutation const& s, std::vector<T> const& v)
    {
        return permute(s, std::span<T const>(v));
    }
}
