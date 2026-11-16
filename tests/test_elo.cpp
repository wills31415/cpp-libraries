#include <wsl/elo.hpp>
#include <cassert>
#include <cmath>
#include <vector>

using namespace wsl;

namespace
{
    bool near(double a, double b, double tolerance) {return std::abs(a - b) <= tolerance;}
}

void scores()
{
    assert(elo_from_score(0.5) == 0);
    assert(near(elo_from_score(0.64), 100, 1));                 // +100 Elo : environ 64 %
    assert(near(elo_from_score(0.91), 402, 2));                 // +400 Elo : environ 91 %
    assert(near(score_from_elo(elo_from_score(0.7)), 0.7, 1e-9));
    assert(near(elo_from_score(0.3), -elo_from_score(0.7), 1e-9));
    assert(std::isfinite(elo_from_score(1)) && std::isfinite(elo_from_score(0)));
}

void estimates()
{
    // 60 % sur 100 parties : environ +70 Elo, intervalle d'environ ±70.
    tally const t{60, 0, 40};
    elo_estimate const e = estimate_elo(t);
    assert(near(e.elo, 70.4, 0.5) && e.low > 0 && e.high < 150 && e.low < e.elo && e.elo < e.high);
    // Les nuls réduisent l'incertitude.
    elo_estimate const d = estimate_elo(tally{20, 80, 0});
    assert(near(d.elo, 70.4, 0.5) && d.high - d.low < e.high - e.low);
    // Le bilan inverse donne l'écart inverse.
    assert(near(estimate_elo(t.reversed()).elo, -e.elo, 1e-9));
    // Tout gagner : écart infini, et seule la borne basse compte. Sur 200 parties, au moins +692 Elo.
    elo_estimate const all = estimate_elo(tally{200, 0, 0});
    assert(std::isinf(all.elo) && all.elo > 0 && near(all.low, 692, 1));
    assert(std::isinf(estimate_elo(tally{0, 0, 200}).elo) && near(estimate_elo(tally{0, 0, 200}).high, -all.low, 1e-9));
}

void ratings()
{
    // Deux agents : Bradley-Terry retrouve l'écart du score, à la partie nulle fictive près.
    std::vector<std::vector<tally>> two(2, std::vector<tally>(2));
    two[0][1] = {36, 0, 64};
    two[1][0] = two[0][1].reversed();
    auto const elo = bradley_terry(two);
    assert(elo[0] == 0 && near(elo[1], elo_from_score(64.5 / 101), 1));

    // Trois agents, A > B > C.
    std::vector<std::vector<tally>> three(3, std::vector<tally>(3));
    auto const set = [&](std::size_t i, std::size_t j, tally t) {three[i][j] = t; three[j][i] = t.reversed();};
    set(0, 1, {75, 0, 25});
    set(1, 2, {75, 0, 25});
    set(0, 2, {90, 0, 10});
    auto const r = rate(three);
    assert(r[0].elo == 0 && r[0].elo > r[1].elo && r[1].elo > r[2].elo);
    for(rating const& x : r) assert(x.low <= x.elo + 1e-9 && x.elo <= x.high + 1e-9);

    // Un agent qui gagne tout garde une note finie.
    set(0, 1, {100, 0, 0});
    assert(std::isfinite(bradley_terry(three)[1]));
}

void convergence()
{
    // La référence (agent 0) très en dessous des autres : l'algorithme MM doit converger vite. Avec 200
    // itérations, les notes vérifient déjà les équations de Bradley-Terry : les points de chacun (partie nulle
    // fictive comprise) valent ceux qu'attend le modèle. Une normalisation qui ne divisait que la référence
    // demandait plus de 2 000 itérations.
    std::size_t const n = 5;
    std::vector<std::vector<tally>> results(n, std::vector<tally>(n));
    auto const set = [&](std::size_t i, std::size_t j, tally t) {results[i][j] = t; results[j][i] = t.reversed();};
    for(std::size_t j = 1; j < n; ++j) set(0, j, {0, 1, 49});
    set(1, 2, {22, 3, 25});
    set(1, 3, {20, 2, 28});
    set(1, 4, {18, 4, 28});
    set(2, 3, {24, 2, 24});
    set(2, 4, {21, 3, 26});
    set(3, 4, {23, 4, 23});

    auto const elo = bradley_terry(results, 200);
    assert(elo[0] == 0 && elo[1] > 500);
    for(std::size_t i = 0; i < n; ++i)
    {
        double points = 0, expected = 0;
        for(std::size_t j = 0; j < n; ++j)
        {
            if(i == j) continue;
            points += results[i][j].points() + 0.5;
            expected += double(results[i][j].games() + 1) * score_from_elo(elo[i] - elo[j]);
        }
        assert(near(points, expected, 1e-6 * points));
    }
}

int main()
{
    scores();
    estimates();
    ratings();
    convergence();
}
