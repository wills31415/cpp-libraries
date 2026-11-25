// Banc d'essai de wsl::mdarray : chaque noyau (-DKERNEL=…) contre la même boucle
// écrite à la main. Un programme par noyau (tools/mdarray/benchmark.sh).
#include <wsl/mdarray.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <vector>

using V = wsl::mdarray<double, wsl::dims<1>>;
using M = wsl::mdarray<double, wsl::dims<2>>;
// Champ : composantes de longueur 3 connues à la compilation, comme dans la boucle à la main.
using F2 = wsl::mdarray<double, wsl::extents<std::size_t, wsl::dynamic_extent, 3>>;
using F3 = wsl::mdarray<double, wsl::extents<std::size_t, wsl::dynamic_extent, 3, 3>>;

template <class F>
double best_ns(F f, int repeat)
{
    double best = 1e300;
    for (int r = 0; r < repeat; ++r)
    {
        auto const t0 = std::chrono::steady_clock::now();
        f();
        auto const t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0).count());
    }
    return best;
}

volatile double sink;

int main()
{
    auto [i, j, k, p] = wsl::axis<4>();
    (void)i; (void)j; (void)k; (void)p;
#if KERNEL == 1 // y = a x + y, tableaux entiers
    constexpr std::size_t n = 1 << 20;
    V x(n), y(n);
    for (std::size_t q = 0; q < n; ++q) { x[q] = 1.0 / (q + 1); y[q] = 0.5; }
    double a = 1.0000001;
    double lib = best_ns([&] { y = a * x + y; }, 20);
    double hand = best_ns([&] { double* yd = y.data(); double const* xd = x.data(); for (std::size_t q = 0; q < n; ++q) yd[q] = a * xd[q] + yd[q]; }, 20);
    std::printf("axpy (tableaux)   %8.3f ns/élément  à la main %8.3f  rapport %5.1f\n", lib / n, hand / n, lib / hand);
#elif KERNEL == 2 // y[i] = a x[i] + y[i], axes
    constexpr std::size_t n = 1 << 20;
    V x(n), y(n);
    for (std::size_t q = 0; q < n; ++q) { x[q] = 1.0 / (q + 1); y[q] = 0.5; }
    double a = 1.0000001;
    double lib = best_ns([&] { y[i] = a * x[i] + y[i]; }, 20);
    double hand = best_ns([&] { double* yd = y.data(); double const* xd = x.data(); for (std::size_t q = 0; q < n; ++q) yd[q] = a * xd[q] + yd[q]; }, 20);
    std::printf("axpy (axes)       %8.3f ns/élément  à la main %8.3f  rapport %5.1f\n", lib / n, hand / n, lib / hand);
#elif KERNEL == 3 // s = sum(x)
    constexpr std::size_t n = 1 << 20;
    V x(n);
    for (std::size_t q = 0; q < n; ++q) x[q] = 1.0 / (q + 1);
    double lib = best_ns([&] { sink = sum(x); }, 20);
    double hand = best_ns([&] { double s = 0; double const* xd = x.data(); for (std::size_t q = 0; q < n; ++q) s += xd[q]; sink = s; }, 20);
    std::printf("somme (tableau)   %8.3f ns/élément  à la main %8.3f  rapport %5.1f\n", lib / n, hand / n, lib / hand);
#elif KERNEL == 4 // s = x[i] x[i]
    constexpr std::size_t n = 1 << 20;
    V x(n);
    for (std::size_t q = 0; q < n; ++q) x[q] = 1.0 / (q + 1);
    double lib = best_ns([&] { double s = x[i] * x[i]; sink = s; }, 20);
    double hand = best_ns([&] { double s = 0; double const* xd = x.data(); for (std::size_t q = 0; q < n; ++q) s += xd[q] * xd[q]; sink = s; }, 20);
    std::printf("produit scalaire  %8.3f ns/élément  à la main %8.3f  rapport %5.1f\n", lib / n, hand / n, lib / hand);
#elif KERNEL == 5 // y[i] = A[i, j] x[j]
    constexpr std::size_t n = 1024;
    M A(n, n);
    V x(n), y(n);
    for (std::size_t q = 0; q < n; ++q) { x[q] = 1.0 / (q + 1); for (std::size_t r = 0; r < n; ++r) A[q, r] = double(q + r) / n; }
    double lib = best_ns([&] { y[i] = A[i, j] * x[j]; }, 10);
    double hand = best_ns([&] { for (std::size_t q = 0; q < n; ++q) { double s = 0; for (std::size_t r = 0; r < n; ++r) s += A.data()[q * n + r] * x.data()[r]; y.data()[q] = s; } }, 10);
    std::printf("matrice × vecteur %8.3f ns/terme    à la main %8.3f  rapport %5.1f\n", lib / (n * n), hand / (n * n), lib / hand);
#elif KERNEL == 6 // C[i, k] = A[i, j] B[j, k]
    constexpr std::size_t n = 192;
    M A(n, n), B(n, n), C(n, n);
    for (std::size_t q = 0; q < n; ++q) for (std::size_t r = 0; r < n; ++r) { A[q, r] = double(q + r) / n; B[q, r] = double(q) - double(r); }
    double lib = best_ns([&] { C[i, k] = A[i, j] * B[j, k]; }, 5);
    double hand = best_ns([&] { for (std::size_t q = 0; q < n; ++q) for (std::size_t r = 0; r < n; ++r) { double s = 0; for (std::size_t t = 0; t < n; ++t) s += A.data()[q * n + t] * B.data()[t * n + r]; C.data()[q * n + r] = s; } }, 5);
    std::printf("produit matriciel %8.3f ns/terme    à la main %8.3f  rapport %5.1f\n", lib / (n * n * n), hand / (n * n * n), lib / hand);
#elif KERNEL == 7 // n2[p] = g[p, i, j] v[p, i] v[p, j]
    constexpr std::size_t P = 1 << 16;
    F3 g(P);
    F2 v(P);
    V n2(P);
    for (std::size_t q = 0; q < P; ++q) for (std::size_t r = 0; r < 3; ++r) { v[q, r] = 1.0 / (q + r + 1); for (std::size_t t = 0; t < 3; ++t) g[q, r, t] = r == t ? 1.0 : 0.1; }
    double lib = best_ns([&] { n2[p] = g[p, i, j] * v[p, i] * v[p, j]; }, 10);
    double hand = best_ns([&] { for (std::size_t q = 0; q < P; ++q) { double s = 0; for (std::size_t r = 0; r < 3; ++r) for (std::size_t t = 0; t < 3; ++t) s += g.data()[q * 9 + r * 3 + t] * v.data()[q * 3 + r] * v.data()[q * 3 + t]; n2.data()[q] = s; } }, 10);
    std::printf("champ g v v       %8.3f ns/point    à la main %8.3f  rapport %5.1f\n", lib / P, hand / P, lib / hand);
#elif KERNEL == 8 // C = A ∘ B + 1, 2D tableaux entiers
    constexpr std::size_t n = 1024;
    M A(n, n), B(n, n), C(n, n);
    for (std::size_t q = 0; q < n; ++q) for (std::size_t r = 0; r < n; ++r) { A[q, r] = double(q + r); B[q, r] = 0.5; }
    double lib = best_ns([&] { C = A * B + 1; }, 10);
    double hand = best_ns([&] { for (std::size_t q = 0; q < n * n; ++q) C.data()[q] = A.data()[q] * B.data()[q] + 1; }, 10);
    std::printf("A∘B + 1 (2D)      %8.3f ns/élément  à la main %8.3f  rapport %5.1f\n", lib / (n * n), hand / (n * n), lib / hand);
#endif
    return 0;
}
