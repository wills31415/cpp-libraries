#include <wsl/progress.hpp>
#include <cassert>
#include <chrono>
#include <cmath>

using namespace wsl;
using namespace std::chrono_literals;

// Le suivi des tâches longues (util/progress.hpp) : la vitesse récente, sur une fenêtre glissante.

namespace
{
    bool near(double a, double b) {return std::abs(a - b) < 1e-9;}
}

void readable()
{
    assert(readable_duration(42) == "42 s" && readable_duration(185) == "3 min 05 s" && readable_duration(4020) == "1 h 07 min");
    assert(readable_rate(12.3) == "12/s" && readable_rate(0.583) == "0.58/s" && readable_rate(0.05) == "3.0/min");
}

void sliding_window()
{
    using clock = recent_rate::clock;
    clock::time_point const t0{};
    recent_rate r(10min);
    r.reset(t0);
    assert(r.per_second(t0, 0) == 0);

    // Avant 10 minutes : la moyenne depuis le début. 1 unité par seconde pendant 5 minutes.
    for(int s = 1; s <= 300; ++s) r.add(t0 + std::chrono::seconds(s), std::size_t(s));
    assert(near(r.per_second(t0 + 300s, 300), 1));

    // Puis 3 fois plus lent pendant 20 minutes : la vitesse passe à 1/3 dès que la fenêtre ne voit plus que la
    // partie lente, alors que la moyenne depuis le début serait encore de 0,5.
    std::size_t done = 300;
    for(int s = 301; s <= 1500; ++s)
    {
        if(s % 3 == 0) ++done;
        r.add(t0 + std::chrono::seconds(s), done);
    }
    double const recent = r.per_second(t0 + 1500s, done);
    assert(std::abs(recent - 1.0 / 3) < 0.01);
    assert(std::abs(double(done) / 1500 - 0.47) < 0.01);

    // Plusieurs appels dans la même seconde ne gardent que le dernier avancement.
    r.add(t0 + 1500s + 200ms, done + 5);
    assert(r.per_second(t0 + 1500s + 200ms, done + 5) > recent);
}

int main()
{
    readable();
    sliding_window();
}
