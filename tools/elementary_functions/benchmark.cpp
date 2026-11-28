// Mesure d'une seule fonction, compilée par tools/elementary_functions/benchmark.sh avec
// FUNCTION, TYPE, LOW et HIGH définis. Un programme par fonction : dans un même
// programme, les boucles de mesure de plusieurs fonctions se sont montrées
// capables de fausser les chiffres d'un facteur 3 ou plus (constaté sur Intel).

#include <wsl/elementary_functions.hpp>
#include <wsl/special_functions.hpp>
#include <wsl/special_functions/detail/pinned_portable.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

int main()
{
    std::vector<TYPE> inputs(1 << 16);
    for (std::size_t i = 0; i < inputs.size(); ++i)
        inputs[i] = TYPE(LOW) + (TYPE(HIGH) - TYPE(LOW)) * TYPE(i) / TYPE(inputs.size());

    double best = 1e300;
    TYPE sum = 0;
    for (int repetition = 0; repetition < 9; ++repetition)
    {
        auto const start = std::chrono::steady_clock::now();
        for (int round = 0; round < 10; ++round)
            for (TYPE x : inputs)
                sum += FUNCTION(x);
        auto const stop = std::chrono::steady_clock::now();
        double const ns = std::chrono::duration<double, std::nano>(stop - start).count() / (10.0 * double(inputs.size()));
        best = ns < best ? ns : best;
    }
    std::printf("%.2f %g\n", best, double(sum));
}
