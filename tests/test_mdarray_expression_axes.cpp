// wsl::mdarray, monde des axes (docs/mdarray.md, V2 à V4, V9, L1 à L14, I6, I10 à I12) :
// les exemples du document, chacun comparé à une boucle écrite à la main. Les refus
// sont dans tests/compile_fail/mdarray/axes.cpp.

#include <wsl/mdarray.hpp>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <type_traits>

namespace
{
    using D1 = wsl::dims<1>;
    using D2 = wsl::dims<2>;
    using D3 = wsl::dims<3>;
    using V = wsl::mdarray<double, D1>;
    using M = wsl::mdarray<double, D2>;
    using T3 = wsl::mdarray<double, D3>;

    bool close(double a, double b) { return std::fabs(a - b) <= 1e-12 * (1 + std::fabs(a) + std::fabs(b)); }

    // Remplissage déterministe : m(i, j) = f(i, j).
    template <class A, class F>
    void fill2(A& a, F f)
    {
        for (std::size_t i = 0; i < a.extent(0); ++i)
            for (std::size_t j = 0; j < a.extent(1); ++j)
                a[i, j] = f(double(i), double(j));
    }

    // Calculs à la compilation : le monde des axes est constexpr.
    constexpr bool compile_time()
    {
        auto [i, j, k] = wsl::axis<3>();
        wsl::mdarray<int, wsl::extents<std::size_t, 2, 2>> a{{1, 2}, {3, 4}}, b{{5, 6}, {7, 8}}, c;
        c[i, k] = a[i, j] * b[j, k];
        int const trace = a[i, i];
        int const frobenius = a[i, j] * a[i, j];
        return c[0, 0] == 19 && c[1, 1] == 50 && trace == 5 && frobenius == 30;
    }
    static_assert(compile_time());

    void linear_algebra()
    {
        auto [i, j, k, l] = wsl::axis<4>();
        M A(3, 4), B(4, 2), C(3, 2);
        fill2(A, [](double r, double s) { return r + 2 * s + 1; });
        fill2(B, [](double r, double s) { return r * s - 1; });
        V x(4), y(3), b(3);
        for (std::size_t r = 0; r < 4; ++r)
            x[r] = 0.5 * double(r) - 1;
        for (std::size_t r = 0; r < 3; ++r)
            b[r] = double(r) * double(r);

        // Produit matriciel et application affine (V3).
        C[i, k] = A[i, j] * B[j, k];
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t s = 0; s < 2; ++s)
            {
                double e = 0;
                for (std::size_t t = 0; t < 4; ++t)
                    e += A[r, t] * B[t, s];
                assert((close(C[r, s], e)));
            }
        y[i] = A[i, j] * x[j] + b[i];
        for (std::size_t r = 0; r < 3; ++r)
        {
            double e = b[r];
            for (std::size_t t = 0; t < 4; ++t)
                e += A[r, t] * x[t];
            assert((close(y[r], e)));
        }

        // Scalaires (V9) : norme, forme bilinéaire, wsl::value.
        double const n = sqrt(x[i] * x[i]);
        double e = 0;
        for (std::size_t t = 0; t < 4; ++t)
            e += x[t] * x[t];
        assert((close(n, std::sqrt(e))));
        M Q(4, 4);
        fill2(Q, [](double r, double s) { return r == s ? 2.0 : 0.25; });
        double const q = x[k] * Q[k, l] * x[l];
        double eq = 0;
        for (std::size_t r = 0; r < 4; ++r)
            for (std::size_t s = 0; s < 4; ++s)
                eq += x[r] * Q[r, s] * x[s];
        assert((close(q, eq)));
        auto const v = wsl::value(x[i] * x[i]);
        static_assert(std::is_same_v<std::remove_const_t<decltype(v)>, double>);
        assert((close(v, e)));
        if (x[i] * x[i] > 1) // comparaison d'un scalaire
            assert((e > 1));

        // Diagonale (lecture, écriture), lignes multipliées, transposée, produit tensoriel (V4).
        M S(4, 4);
        S[i, j] = 0;
        S[i, i] = 1;
        assert((S[2, 2] == 1 && S[2, 1] == 0));
        V d(4);
        d[i] = Q[i, i];
        assert((d[3] == 2));
        M R(3, 4);
        R[i, j] = A[i, j] * b[i];
        assert((close(R[2, 3], A[2, 3] * b[2])));
        M At(4, 3);
        At[i, j] = A[j, i];
        assert((At[3, 1] == A[1, 3]));
        M O(3, 4);
        O[i, j] = b[i] * x[j];
        assert((close(O[2, 1], b[2] * x[1])));

        // Diffusion le long des axes de gauche, remplissage d'une colonne (V4).
        M W(3, 4);
        W[i, j] = b[i];
        W[i, 0] = 0;
        assert((W[2, 3] == 4 && W[2, 0] == 0));

        // Affectations composées : vrais raccourcis, i à gauche jamais sommé (V9).
        M Y(3, 4);
        Y[i, j] = A[i, j];
        Y[i, j] *= b[i];
        assert((close(Y[2, 3], R[2, 3])));
        Y[i, j] += 1;
        assert((close(Y[0, 0], A[0, 0] * b[0] + 1)));
    }

    void reductions()
    {
        auto [i, j] = wsl::axis<2>();
        M a(3, 4);
        fill2(a, [](double r, double s) { return (r + 1) * (s - 1.5); });
        V r(3), m(3);
        wsl::mdarray<std::ptrdiff_t, D1> idx(3);
        r[i] = sum(j)(a[i, j]); // sommes des lignes (L8)
        m[i] = max(j)(a[i, j]);
        idx[i] = argmax(j)(a[i, j]);
        for (std::size_t p = 0; p < 3; ++p)
        {
            double s = 0, mx = a[p, 0];
            for (std::size_t q = 0; q < 4; ++q)
            {
                s += a[p, q];
                mx = std::fmax(mx, a[p, q]);
            }
            assert((close(r[p], s) && m[p] == mx && idx[p] == 3));
        }
        // Entièrement liée : calculée tout de suite, auto permis (V9).
        auto const total = sum(i, j)(a[i, j]);
        static_assert(std::is_same_v<std::remove_const_t<decltype(total)>, double>);
        assert((close(total, sum(a))));
        auto const positives = count(i, j)(a[i, j] > 0);
        assert((positives == 6));
        assert((any(i, j)(a[i, j] > 3) && !all(i, j)(a[i, j] > 0)));
        // Softmax (une ligne).
        V x(4), y(4);
        for (std::size_t p = 0; p < 4; ++p)
            x[p] = double(p) / 4;
        y[i] = exp(x[i]) / sum(j)(exp(x[j]));
        double z = 0;
        for (std::size_t p = 0; p < 4; ++p)
            z += std::exp(x[p]);
        assert((close(y[2], std::exp(x[2]) / z)));
        // max(i, j) : élément par élément (des positions).
        wsl::mdarray<std::ptrdiff_t, D2> mx(3, 3);
        mx[i, j] = max(i, j);
        assert((mx[1, 2] == 2 && mx[2, 0] == 2));
    }

    void values_and_functions()
    {
        auto [i, j] = wsl::axis<2>();
        // Axes comme valeurs (L12), Kronecker (L13).
        M h(3, 3);
        h[i, j] = 1.0 / (i + j + 1);
        assert((close(h[1, 2], 1.0 / 4)));
        M I3(3, 3);
        I3[i, j] = wsl::kronecker(i, j);
        assert((I3[1, 1] == 1 && I3[1, 2] == 0));
        V p(4);
        for (std::size_t r = 0; r < 4; ++r)
            p[r] = 0.25;
        double const mean = i * p[i];
        assert((close(mean, 1.5)));
        M T(3, 3);
        T[i, j] = exp(-abs(i - j));
        assert((close(T[0, 2], std::exp(-2.0))));
        M L(3, 3);
        L[i, j] = where(i >= j, h[i, j], 0);
        assert((L[0, 1] == 0 && close(L[2, 1], h[2, 1])));

        // Matrice des distances, noyau gaussien, fonction de l'utilisateur (I11).
        V x(3), y(3);
        for (std::size_t r = 0; r < 3; ++r)
        {
            x[r] = double(r);
            y[r] = 2 * double(r);
        }
        M D(3, 3), K(3, 3), F2(3, 3);
        D[i, j] = hypot(x[i] - x[j], y[i] - y[j]);
        assert((close(D[0, 2], std::hypot(2.0, 4.0))));
        K[i, j] = exp(-(x[i] - x[j]) * (x[i] - x[j]));
        assert((close(K[0, 1], std::exp(-1.0))));
        auto F = wsl::vectorize([](double u, double v) { return u * 10 + v; });
        F2[i, j] = F(x[i], x[j]);
        assert((F2[2, 1] == 21));
        // Perte totale : huber(p[i], t[i]) somme i (convention C).
        auto huber = wsl::vectorize([](double u, double v) { double d = u - v; return std::fabs(d) < 1 ? d * d / 2 : std::fabs(d) - 0.5; });
        double const loss = huber(x[i], y[i]);
        assert((close(loss, 0 + 0.5 + 1.5)));
        // Expression de tableaux entiers indicée.
        V s(3);
        s[i] = (x + y)[i] * 2;
        assert((s[2] == 12));
        // Comparaisons par les axes : élément par élément (I10).
        wsl::mdarray<bool, D1> mask(3);
        mask[i] = x[i] < 1.5;
        assert((mask[1] && !mask[2]));
        double const masked = (x[i] < 1.5) * y[i];
        assert((masked == 2));
    }

    void fields_and_mechanics()
    {
        auto [p, i, j, k, l] = wsl::axis<5>();
        // En chaque point p : gₚᵢⱼ vᵢ vⱼ (V3 : p à gauche, jamais sommé).
        T3 g(4, 3, 3);
        M v(4, 3);
        for (std::size_t a = 0; a < 4; ++a)
            for (std::size_t b = 0; b < 3; ++b)
            {
                v[a, b] = double(a) + double(b) / 10;
                for (std::size_t c = 0; c < 3; ++c)
                    g[a, b, c] = b == c ? 1.0 + double(a) : 0.1;
            }
        V n2(4);
        n2[p] = g[p, i, j] * v[p, i] * v[p, j];
        for (std::size_t a = 0; a < 4; ++a)
        {
            double e = 0;
            for (std::size_t b = 0; b < 3; ++b)
                for (std::size_t c = 0; c < 3; ++c)
                    e += g[a, b, c] * v[a, b] * v[a, c];
            assert((close(n2[a], e)));
        }
        // Hooke isotrope : σᵢⱼ = λ εₖₖ δᵢⱼ + 2μ εᵢⱼ.
        M eps{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}, sigma(3, 3);
        double const lambda = 2, mu = 3;
        sigma[i, j] = lambda * eps[k, k] * wsl::kronecker(i, j) + 2 * mu * eps[i, j];
        assert((sigma[0, 0] == 2 * 11 + 6 * 1 && sigma[0, 1] == 6 * 2));
        // Matrice de masse : Mᵢⱼ = Σ_q w_q φᵢ(q) φⱼ(q).
        V w{3};
        w[0] = 0.5;
        w[1] = 1;
        w[2] = 0.5;
        M phi{{1, 0.5, 0}, {0, 0.5, 1}}, mass(2, 2);
        mass[i, j] = w[l] * phi[i, l] * phi[j, l];
        assert((close(mass[0, 0], 0.5 + 0.25) && close(mass[0, 1], 0.25)));
    }

    void aliasing()
    {
        auto [i, j] = wsl::axis<2>();
        // y = A y sur place (I12) : un temporaire.
        M A{{1, 2}, {3, 4}};
        V y(2);
        y[0] = 1;
        y[1] = 1;
        y[i] = A[i, j] * y[j];
        assert((y[0] == 3 && y[1] == 7));
        // Transposée sur place.
        M B{{1, 2}, {3, 4}};
        B[i, j] = B[j, i];
        assert((B[0, 1] == 3 && B[1, 0] == 2));
        // Mêmes positions : sur place, sans temporaire nécessaire.
        B[i, j] = B[i, j] * 2;
        assert((B[1, 1] == 8));
    }
}

int main()
{
    linear_algebra();
    reductions();
    values_and_functions();
    fields_and_mechanics();
    aliasing();
    std::puts("test_mdarray_expression_axes : OK");
}
