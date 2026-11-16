// Banc d'essai informatif : wsl::special contre <cmath>, en nanosecondes par appel.
// À lancer en Release : ./build/release/examples/special_benchmark

#include <wsl/special.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

namespace sp = wsl::special;

namespace
{
    // Minimum sur plusieurs mesures : sur une machine partagée (fréquence variable,
    // autres processus), le minimum est la mesure la plus fidèle du coût propre.
    template <class T, class F>
    double nanoseconds_per_call(std::vector<T> const& inputs, F f, T& sink)
    {
        double best = 1e300;
        T sum = 0;
        for (int repetition = 0; repetition < 9; ++repetition)
        {
            auto const start = std::chrono::steady_clock::now();
            for (int round = 0; round < 10; ++round)
                for (T x : inputs)
                    sum += f(x);
            auto const stop = std::chrono::steady_clock::now();
            double const ns = std::chrono::duration<double, std::nano>(stop - start).count() / (10.0 * double(inputs.size()));
            best = ns < best ? ns : best;
        }
        sink += sum;
        return best;
    }

    template <class T>
    std::vector<T> inputs(T low, T high)
    {
        std::vector<T> v(1 << 16);
        for (std::size_t i = 0; i < v.size(); ++i)
            v[i] = low + (high - low) * T(i) / T(v.size());
        return v;
    }
}

int main()
{
    double sink = 0;
    float sinkf = 0;
    auto const xd = inputs<double>(-700, 700), ld = inputs<double>(1e-6, 1e6);
    auto const xf = inputs<float>(-80, 80), lf = inputs<float>(1e-6f, 1e6f);

    std::printf("%-12s %10s %10s\n", "", "wsl", "<cmath>");
    std::printf("%-12s %8.2f ns %8.2f ns\n", "exp double",
                nanoseconds_per_call(xd, [](double x) { return sp::exp(x); }, sink),
                nanoseconds_per_call(xd, [](double x) { return std::exp(x); }, sink));
    std::printf("%-12s %8.2f ns %8.2f ns\n", "log double",
                nanoseconds_per_call(ld, [](double x) { return sp::log(x); }, sink),
                nanoseconds_per_call(ld, [](double x) { return std::log(x); }, sink));
    std::printf("%-12s %8.2f ns %8.2f ns\n", "exp float",
                nanoseconds_per_call(xf, [](float x) { return sp::exp(x); }, sinkf),
                nanoseconds_per_call(xf, [](float x) { return std::exp(x); }, sinkf));
    std::printf("%-12s %8.2f ns %8.2f ns\n", "log float",
                nanoseconds_per_call(lf, [](float x) { return sp::log(x); }, sinkf),
                nanoseconds_per_call(lf, [](float x) { return std::log(x); }, sinkf));
    std::printf("(somme de contrôle : %g)\n", sink + double(sinkf));
}
