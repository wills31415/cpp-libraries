// wsl::mdarray, monde des axes (docs/mdarray.md, V2 à V4, V9, L1 à L14, I6, I10 à I12) :
// les exemples du document, chacun comparé à une boucle écrite à la main. Les refus
// sont dans tests/compile_fail/mdarray/axes.cpp.

#include <wsl/mdarray.hpp>

#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <type_traits>
#include <vector>

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

    // Accès à pas fixe : mêmes bits que l'évaluation générale, qu'on force en lisant
    // (A * 1.0)[i, j] (une expression, hors du chemin à pas fixe ; × 1 est exact). Mode
    // reproductible : sans fma, les deux formes font exactement les mêmes opérations.
    template <class X, class Y>
    bool same_bits(X const& x, Y const& y)
    {
        for (std::size_t q = 0; q < x.size(); ++q)
            if (std::bit_cast<std::uint64_t>(x.data()[q]) != std::bit_cast<std::uint64_t>(y.data()[q]))
                return false;
        return true;
    }

    void strided_access()
    {
        auto [i, j, k, p] = wsl::axis<4>();
        constexpr std::size_t n = 203; // au-delà de 128 termes : somme partagée en deux moitiés
        M A(n, n), B(n, 5);
        fill2(A, [](double r, double c) { return 1 / (1 + r + 2 * c); });
        fill2(B, [](double r, double c) { return std::sqrt(1 + r * c); });
        V x(n), y1(n), y2(n);
        for (std::size_t q = 0; q < n; ++q)
            x[q] = 1 / (3.0 + q);

        // Produit matrice × vecteur, produit matriciel.
        y1[i] = wsl::reproducible(A[i, j] * x[j]);
        y2[i] = wsl::reproducible((A * 1.0)[i, j] * x[j]);
        assert(same_bits(y1, y2));
        M C1(n, 5), C2(n, 5);
        C1[i, k] = wsl::reproducible(A[i, j] * B[j, k]);
        C2[i, k] = wsl::reproducible(A[i, j] * (B * 1.0)[j, k]);
        assert(same_bits(C1, C2));

        // Voies : 21 colonnes (deux paquets de 8 et un reste), valeurs d'axes, position fixe,
        // scalaire, cible relue.
        M B21(n, 21), E1(n, 21), E2(n, 21);
        fill2(B21, [](double r, double c) { return 1 / (2 + r + 0.5 * c); });
        E1[i, k] = wsl::reproducible((A[i, j] + k) * B21[j, k] * 2.5 + i + A[3, j] * B21[j, k]);
        E2[i, k] = wsl::reproducible((A[i, j] + k) * (B21 * 1.0)[j, k] * 2.5 + i + A[3, j] * B21[j, k]);
        assert(same_bits(E1, E2));
        E1[i, k] = wsl::reproducible(A[i, j] * B21[j, k] + E1[i, k]);
        E2[i, k] = wsl::reproducible(A[i, j] * (B21 * 1.0)[j, k] + E2[i, k]);
        assert(same_bits(E1, E2));
        M Q1(n, 21), Q2(n, 21);
        Q1[i, k] = wsl::reproducible(wsl::kronecker(i, k) * B21[j, k] * x[j]);
        Q2[i, k] = wsl::reproducible(wsl::kronecker(i, k) * (B21 * 1.0)[j, k] * x[j]);
        assert(same_bits(Q1, Q2));

        // Découpage du cache le long de l'axe sommé (plus de 128 termes : tranches de la
        // somme par paires, 37 lignes = un bloc de 32 et un reste) ; rang 3 aussi.
        M A300(37, 300), B300(300, 21), G1(37, 21), G2(37, 21);
        fill2(A300, [](double r, double c) { return 1 / (1 + r + 0.25 * c); });
        fill2(B300, [](double r, double c) { return std::sqrt(2 + r + c); });
        G1[i, k] = wsl::reproducible(A300[i, j] * B300[j, k]);
        G2[i, k] = wsl::reproducible(A300[i, j] * (B300 * 1.0)[j, k]);
        assert(same_bits(G1, G2));
        V w300(300);
        for (std::size_t q = 0; q < 300; ++q)
            w300[q] = 1 / (7.0 + q);
        G1[i, k] = wsl::reproducible(A300[i, j] * w300[j] * B300[j, k]); // w[j] recopié aussi
        G2[i, k] = wsl::reproducible(A300[i, j] * w300[j] * (B300 * 1.0)[j, k]);
        assert(same_bits(G1, G2));
        T3 X(3, 5, 300), H1(3, 5, 21), H2(3, 5, 21);
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 5; ++b)
                for (std::size_t c = 0; c < 300; ++c)
                    X[a, b, c] = 1 / (1.0 + a + 2 * b + 0.5 * c);
        H1[p, i, k] = wsl::reproducible(X[p, i, j] * B300[j, k]);
        H2[p, i, k] = wsl::reproducible(X[p, i, j] * (B300 * 1.0)[j, k]);
        assert(same_bits(H1, H2));

        // Somme le long de deux axes (indice aplati, 12 × 25 = 300 termes) : pas cohérents,
        // somme aplatie, t[i, j, k] recopié ; uj[j, i], contigu le long de i (pas 1 fixé) :
        // ni aplatie ni recopiée.
        auto [l] = wsl::axis<1>();
        T3 gij(9, 12, 25), t(12, 25, 21);
        M uj(25, 12), c1(9, 21), c2(9, 21);
        for (std::size_t a = 0; a < 12; ++a)
            for (std::size_t b = 0; b < 25; ++b)
            {
                uj[b, a] = 1 / (3.0 + a + b);
                for (std::size_t c = 0; c < 9; ++c)
                    gij[c, a, b] = 1 / (1.0 + c + 0.5 * a + 0.25 * b);
                for (std::size_t c = 0; c < 21; ++c)
                    t[a, b, c] = std::sqrt(1.0 + a + 2 * b + 3 * c);
            }
        c1[p, l] = wsl::reproducible(gij[p, i, j] * t[i, j, l]);
        c2[p, l] = wsl::reproducible(gij[p, i, j] * (t * 1.0)[i, j, l]);
        assert(same_bits(c1, c2));
        c1[p, l] = wsl::reproducible(gij[p, i, j] * uj[j, i] * t[i, j, l]);
        c2[p, l] = wsl::reproducible(gij[p, i, j] * uj[j, i] * (t * 1.0)[i, j, l]);
        assert(same_bits(c1, c2));
        // Valeur d'axe dans le corps : somme non aplatie.
        c1[p, l] = wsl::reproducible(gij[p, i, j] * t[i, j, l] * (j + 1));
        c2[p, l] = wsl::reproducible(gij[p, i, j] * (t * 1.0)[i, j, l] * (j + 1));
        assert(same_bits(c1, c2));

        // Dispositions layout_left et layout_stride (une colonne sur deux).
        std::vector<double> buffer(2 * n * 5);
        for (std::size_t q = 0; q < buffer.size(); ++q)
            buffer[q] = 1 / (1.0 + q);
        wsl::mdspan<double, D2, wsl::layout_left> L(buffer.data(), n, 5);
        wsl::mdspan<double, D2, wsl::layout_stride> S(buffer.data(), wsl::layout_stride::mapping<D2>(D2(n, 5), std::array{2, 2 * int(n)}));
        V u(5), v1(5), v2(5);
        u[i] = wsl::reproducible(L[j, i] * x[j]);
        v1[i] = wsl::reproducible(S[j, i] * x[j]);
        for (std::size_t c = 0; c < 5; ++c)
        {
            double sl = 0, ss = 0;
            for (std::size_t r = 0; r < n; ++r)
            {
                sl += buffer[r + n * c] * x[r];
                ss += buffer[2 * r + 2 * n * c] * x[r];
            }
            assert((close(u[c], sl) && close(v1[c], ss)));
        }
        v2[i] = wsl::reproducible((S * 1.0)[j, i] * x[j]);
        assert(same_bits(v1, v2));

        // Position fixe, diagonale, deux axes sommés, valeurs d'axes, δ.
        double const row = wsl::reproducible(A[2, j] * x[j]);
        assert(row == wsl::reproducible((A * 1.0)[2, j] * x[j]));
        V d(n);
        d[i] = A[i, i];
        assert((d[7] == A[7, 7]));
        double const frobenius = wsl::reproducible(C1[i, k] * B[i, k] + 0 * B[i, k] * B[i, k]);
        double const reference = wsl::reproducible((C1 * 1.0)[i, k] * B[i, k] + 0 * B[i, k] * B[i, k]);
        assert(frobenius == reference);
        double const moment = x[i] * i;
        double hand = 0;
        for (std::size_t q = 0; q < n; ++q)
            hand += x[q] * double(q);
        assert(close(moment, hand));
        y1[i] = wsl::kronecker(i, j) * x[j];
        assert(same_bits(y1, x));

        // Sommes imbriquées : n2[p] = g[p, i, j] v[p, i] v[p, j], longueurs statiques et dynamiques.
        constexpr std::size_t P = 11;
        wsl::mdarray<double, wsl::extents<std::size_t, wsl::dynamic_extent, 3, 3>> g(P);
        wsl::mdarray<double, wsl::extents<std::size_t, wsl::dynamic_extent, 3>> w(P);
        T3 gd(P, 3, 3);
        M wd(P, 3);
        for (std::size_t q = 0; q < P; ++q)
            for (std::size_t r = 0; r < 3; ++r)
            {
                w[q, r] = wd[q, r] = 1 / (1.0 + q + r);
                for (std::size_t t = 0; t < 3; ++t)
                    g[q, r, t] = gd[q, r, t] = r == t ? 1 + 0.1 * q : 0.1 * (r + t);
            }
        V n1(P), n2(P), n3(P);
        n1[p] = wsl::reproducible(g[p, i, j] * w[p, i] * w[p, j]);
        n2[p] = wsl::reproducible(gd[p, i, j] * wd[p, i] * wd[p, j]);
        n3[p] = wsl::reproducible((gd * 1.0)[p, i, j] * wd[p, i] * wd[p, j]);
        assert((same_bits(n1, n2) && same_bits(n2, n3)));

        // Réductions explicites, cible relue, longueur nulle.
        V m1(n), s1(n), s2(n);
        m1[i] = wsl::max(j)(A[i, j]);
        assert((m1[4] == A[4, 0]));
        s1[i] = wsl::reproducible(wsl::sum(j)(A[i, j] * x[j]));
        s2[i] = wsl::reproducible(A[i, j] * x[j]);
        assert(same_bits(s1, s2));
        y1[i] = 2 * x[i] + y1[i];
        assert((y1[5] == 3 * x[5]));
        V z(0);
        M A0(3, 0);
        V y3(D1(3), 1.0);
        double const empty = z[i] * z[i];
        y3[i] = A0[i, j] * z[j];
        assert((empty == 0 && y3[0] == 0 && y3[2] == 0));
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
    strided_access();
    aliasing();
    std::puts("test_mdarray_expression_axes : OK");
}
