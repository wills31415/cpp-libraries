#pragma once

// Matrices aléatoires : groupes orthogonal et unitaire sous la mesure de Haar, ensembles
// gaussiens orthogonal (GOE) et unitaire (GUE).
//
//   wsl::random::haar_orthogonal<double, 3> O;   wsl::random::haar_unitary<double> U({.dimension = 4});
//   wsl::random::goe<double, 5> H;               wsl::random::gue<double, 5> K;
//   auto Q = O(g);  H.log_pdf(H(g));
//
// Haar : Q de la factorisation QR d'une matrice gaussienne (réelle, ou complexe de
// parties réelle et imaginaire N(0, ½)) à diagonale de R positive, par Gram–Schmidt
// réorthogonalisé (deux projections par colonne, produits scalaires en double-double)
// — loi exactement de Haar (Mezzadri 2007), orthogonalité à quelques ulp.
// GOE, GUE : densité ∝ e^{−tr H²/2} sur les matrices symétriques (hermitiennes) :
// diagonale N(0, 1), hors diagonale N(0, ½) (parties réelle et imaginaire N(0, ½) pour
// le GUE) ; log-densité −½ Σ hᵢᵢ² − Σᵢ<ⱼ |hᵢⱼ|² − (n/2) log 2π − (n(n − 1)/4) log π
// (n(n − 1)/2 pour le GUE) ; −∞ si H n'est pas symétrique (hermitienne).
//
// Valeurs : wsl::random::square_matrix<T, N> ou wsl::random::square_matrix<std::complex<T>, N>. Les
// produits complexes sont calculés ici (fmul), jamais par l'opérateur * de
// std::complex, qui se contracterait en fma selon les options du projet.
// Tirage (approché, bit-exact) : normales exactes, colonne par colonne (Haar), ou
// triangle supérieur ligne par ligne (GOE, GUE), parties réelles puis imaginaires.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/special.hpp>

#include <complex>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    namespace detail
    {
        // Gram–Schmidt réorthogonalisé, colonnes réelles (Re) ou complexes (Re, Im).
        // a : n colonnes de n coefficients (re[j][i], im[j][i]).
        template <bool Complex>
        constexpr void orthonormalize(std::size_t n, std::vector<std::vector<double>>& re, std::vector<std::vector<double>>& im)
        {
            for (std::size_t j = 0; j < n; ++j)
            {
                for (int pass = 0; pass < 2; ++pass)
                    for (std::size_t k = 0; k < j; ++k)
                    {
                        // c = ⟨q_k, a_j⟩ = Σ conj(q) a
                        auto const cr = dot(n, [&](std::size_t i) { return std::pair{re[k][i], re[j][i]}; });
                        double c_re = cr.hi + cr.lo, c_im = 0;
                        if constexpr (Complex)
                        {
                            auto const ci = dot(n, [&](std::size_t i) { return std::pair{im[k][i], im[j][i]}; });
                            auto const r2 = sum(cr, ci);
                            auto const t1 = dot(n, [&](std::size_t i) { return std::pair{re[k][i], im[j][i]}; });
                            auto const t2 = dot(n, [&](std::size_t i) { return std::pair{im[k][i], re[j][i]}; });
                            auto const r3 = sum(t1, {-t2.hi, -t2.lo});
                            c_re = r2.hi + r2.lo;
                            c_im = r3.hi + r3.lo;
                        }
                        for (std::size_t i = 0; i < n; ++i)                // a −= c q
                        {
                            if constexpr (Complex)
                            {
                                re[j][i] = (re[j][i] - sd::fmul(c_re, re[k][i])) + sd::fmul(c_im, im[k][i]);
                                im[j][i] = (im[j][i] - sd::fmul(c_re, im[k][i])) - sd::fmul(c_im, re[k][i]);
                            }
                            else
                                re[j][i] = re[j][i] - sd::fmul(c_re, re[k][i]);
                        }
                    }
                auto r2 = dot(n, [&](std::size_t i) { return std::pair{re[j][i], re[j][i]}; });
                if constexpr (Complex)
                    r2 = sum(r2, dot(n, [&](std::size_t i) { return std::pair{im[j][i], im[j][i]}; }));
                scaled const r = sqrt(normalized(r2.hi, r2.lo));
                for (std::size_t i = 0; i < n; ++i)
                {
                    auto const divide = [&](double v) {
                        double const a = magnitude(v);
                        double const q = a == 0 ? 0.0 : rounded(normalized(a, 0.0) / r);
                        return v < 0 ? -q : q;
                    };
                    re[j][i] = divide(re[j][i]);
                    if constexpr (Complex)
                        im[j][i] = divide(im[j][i]);
                }
            }
        }

        template <class T, std::size_t N>
        constexpr interval<square_matrix<T, N>> box(std::size_t n, T lo, T up)
        {
            auto a = make_matrix<T, N>(n), b = make_matrix<T, N>(n);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                {
                    a(i, j) = lo;
                    b(i, j) = up;
                }
            return {a, b};
        }

        // Dimension : N, ou donnée (et vérifiée).
        template <std::size_t N>
        constexpr std::size_t checked(std::size_t n, char const* name)
        {
            if (n < 1 || (N != dynamic && n != N))
                throw std::domain_error(name);
            return n;
        }
    }

    struct matrix_dimension
    {
        std::size_t dimension;
    };

    template <law_real T = double, std::size_t N = dynamic>
    class haar_orthogonal : public law_interface
    {
    public:
        using value_type = square_matrix<T, N>;
        using parameters = matrix_dimension;
        static constexpr exactness sampling = exactness::approximate;

        constexpr haar_orthogonal() requires(N != dynamic) : n_(N) {}
        constexpr haar_orthogonal(parameters p) : n_(detail::checked<N>(p.dimension, "wsl::random::haar_orthogonal : dimension invalide")) {}

        constexpr std::size_t dimension() const noexcept { return n_; }
        constexpr interval<value_type> support() const { return detail::box<T, N>(n_, T(-1), T(1)); }
        constexpr value_type mean() const { return detail::box<T, N>(n_, T(0), T(0)).lower; }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& O = *this;
            std::size_t const n = O.n_;
            normal<double> const standard({.mean = 0.0, .stddev = 1.0});
            std::vector<std::vector<double>> re(n, std::vector<double>(n)), im;
            for (auto& column : re)
                for (auto& v : column)
                    v = wsl::random::sample(standard, g);
            detail::orthonormalize<false>(n, re, im);
            auto q = detail::make_matrix<T, N>(n);
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t i = 0; i < n; ++i)
                    q(i, j) = static_cast<T>(re[j][i]);
            return q;
        }

        friend constexpr bool operator==(haar_orthogonal const&, haar_orthogonal const&) = default;

    private:
        std::size_t n_;
    };

    template <law_real T = double, std::size_t N = dynamic>
    class haar_unitary : public law_interface
    {
    public:
        using value_type = square_matrix<std::complex<T>, N>;
        using parameters = matrix_dimension;
        static constexpr exactness sampling = exactness::approximate;

        constexpr haar_unitary() requires(N != dynamic) : n_(N) {}
        constexpr haar_unitary(parameters p) : n_(detail::checked<N>(p.dimension, "wsl::random::haar_unitary : dimension invalide")) {}

        constexpr std::size_t dimension() const noexcept { return n_; }
        constexpr interval<value_type> support() const { return detail::box<std::complex<T>, N>(n_, {T(-1), T(-1)}, {T(1), T(1)}); }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& U = *this;
            std::size_t const n = U.n_;
            normal<double> const half({.mean = 0.0, .stddev = detail::inv_sqrt2_hi});
            std::vector<std::vector<double>> re(n, std::vector<double>(n)), im(n, std::vector<double>(n));
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t i = 0; i < n; ++i)
                {
                    re[j][i] = wsl::random::sample(half, g);
                    im[j][i] = wsl::random::sample(half, g);
                }
            detail::orthonormalize<true>(n, re, im);
            auto q = detail::make_matrix<std::complex<T>, N>(n);
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t i = 0; i < n; ++i)
                    q(i, j) = {static_cast<T>(re[j][i]), static_cast<T>(im[j][i])};
            return q;
        }

        friend constexpr bool operator==(haar_unitary const&, haar_unitary const&) = default;

    private:
        std::size_t n_;
    };

    // GOE, GUE : log-densité commune, Σ des carrés en double-double.
    template <law_real T = double, std::size_t N = dynamic>
    class goe : public law_interface
    {
    public:
        using value_type = square_matrix<T, N>;
        using parameters = matrix_dimension;
        static constexpr exactness sampling = exactness::approximate;

        constexpr goe() requires(N != dynamic) : n_(N) {}
        constexpr goe(parameters p) : n_(detail::checked<N>(p.dimension, "wsl::random::goe : dimension invalide")) {}

        constexpr std::size_t dimension() const noexcept { return n_; }
        constexpr interval<value_type> support() const { return detail::box<T, N>(n_, -std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()); }
        constexpr value_type mean() const { return detail::box<T, N>(n_, T(0), T(0)).lower; }

        constexpr double log_pdf(value_type const& h) const noexcept
        {
            auto const& H = *this;
            auto const l = H.log_density_dd(h);
            return l.hi + l.lo;
        }
        constexpr double pdf(value_type const& h) const noexcept { auto const& H = *this; return detail::rounded_exp(H.log_density_dd(h)); }

        constexpr detail::sd::double_double log_density_dd(value_type const& h) const noexcept
        {
            namespace fp = detail::sd;
            std::size_t const n = n_;
            constexpr fp::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            if (h.size() != n)
                return none;
            fp::double_double s{0.0, 0.0};
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = i; j < n; ++j)
                {
                    double const v = h(i, j);
                    if (!(v == double(h(j, i))) || !detail::is_finite(v))
                        return none;
                    auto const q = detail::product(v, v);
                    s = detail::sum(s, i == j ? fp::double_double{fp::fmul(q.hi, 0.5), fp::fmul(q.lo, 0.5)} : q);
                }
            return normalized_log(n, s, 0.25);
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& H = *this;
            std::size_t const n = H.n_;
            normal<double> const one({.mean = 0.0, .stddev = 1.0}), half({.mean = 0.0, .stddev = detail::inv_sqrt2_hi});
            auto h = detail::make_matrix<T, N>(n);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = i; j < n; ++j)
                    h(i, j) = h(j, i) = static_cast<T>(wsl::random::sample(i == j ? one : half, g));
            return h;
        }

        friend constexpr bool operator==(goe const&, goe const&) = default;

        // −s − (n/2) log 2π − c n(n − 1) log π
        static constexpr detail::sd::double_double normalized_log(std::size_t n, detail::sd::double_double s, double c) noexcept
        {
            namespace fp = detail::sd;
            auto l = detail::sum({-s.hi, -s.lo}, detail::mul({-double(n), 0.0}, {detail::half_log_2pi_hi, detail::half_log_2pi_lo}));
            auto const k = detail::mul({-fp::fmul(c, double(n) * double(n - 1)), 0.0}, {fp::log_pi_hi, fp::log_pi_lo});
            return detail::sum(l, k);
        }
        // e^{l}, un arrondi.

    private:
        std::size_t n_;
    };

    template <law_real T = double, std::size_t N = dynamic>
    class gue : public law_interface
    {
    public:
        using value_type = square_matrix<std::complex<T>, N>;
        using parameters = matrix_dimension;
        static constexpr exactness sampling = exactness::approximate;

        constexpr gue() requires(N != dynamic) : n_(N) {}
        constexpr gue(parameters p) : n_(detail::checked<N>(p.dimension, "wsl::random::gue : dimension invalide")) {}

        constexpr std::size_t dimension() const noexcept { return n_; }
        constexpr interval<value_type> support() const
        {
            constexpr T inf = std::numeric_limits<T>::infinity();
            return detail::box<std::complex<T>, N>(n_, {-inf, -inf}, {inf, inf});
        }

        constexpr double log_pdf(value_type const& h) const noexcept
        {
            auto const& H = *this;
            auto const l = H.log_density_dd(h);
            return l.hi + l.lo;
        }
        constexpr double pdf(value_type const& h) const noexcept { auto const& H = *this; return detail::rounded_exp(H.log_density_dd(h)); }

        constexpr detail::sd::double_double log_density_dd(value_type const& h) const noexcept
        {
            namespace fp = detail::sd;
            std::size_t const n = n_;
            constexpr fp::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            if (h.size() != n)
                return none;
            fp::double_double s{0.0, 0.0};
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = i; j < n; ++j)
                {
                    double const a = h(i, j).real(), b = h(i, j).imag();
                    if (!(a == double(h(j, i).real())) || !(b == -double(h(j, i).imag())) || !detail::is_finite(a) || !detail::is_finite(b) || (i == j && b != 0))
                        return none;
                    auto const q = detail::sum(detail::product(a, a), detail::product(b, b));
                    s = detail::sum(s, i == j ? fp::double_double{fp::fmul(q.hi, 0.5), fp::fmul(q.lo, 0.5)} : q);
                }
            return goe<T, N>::normalized_log(n, s, 0.5);
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& H = *this;
            std::size_t const n = H.n_;
            normal<double> const one({.mean = 0.0, .stddev = 1.0}), half({.mean = 0.0, .stddev = detail::inv_sqrt2_hi});
            auto h = detail::make_matrix<std::complex<T>, N>(n);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = i; j < n; ++j)
                {
                    if (i == j)
                        h(i, i) = {static_cast<T>(wsl::random::sample(one, g)), T(0)};
                    else
                    {
                        T const a = static_cast<T>(wsl::random::sample(half, g));
                        T const b = static_cast<T>(wsl::random::sample(half, g));
                        h(i, j) = {a, b};
                        h(j, i) = {a, -b};
                    }
                }
            return h;
        }

        friend constexpr bool operator==(gue const&, gue const&) = default;

    private:
        std::size_t n_;
    };
}
