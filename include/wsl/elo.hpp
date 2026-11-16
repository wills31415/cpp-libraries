#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <random>
#include <string>
#include <vector>


// Classements Elo, pour mesurer la force des agents.
//
// Une partie rapporte 1 point au vainqueur, 0,5 à chacun en cas de nul. Un écart d'Elo D correspond à un
// score moyen attendu s = 1 / (1 + 10^(-D/400)) : +100 Elo, environ 64 % des points ; +400, environ 91 %.

namespace wsl
{

    // Bilan d'un agent face à un autre.
    struct tally
    {
        std::size_t wins = 0, draws = 0, losses = 0;

        constexpr std::size_t games() const noexcept {return wins + draws + losses;}
        constexpr double points() const noexcept {return double(wins) + 0.5 * double(draws);}
        constexpr double score() const noexcept {return games() ? points() / double(games()) : 0.5;}

        constexpr tally reversed() const noexcept {return {losses, draws, wins};}
        friend constexpr bool operator == (tally const&, tally const&) = default;
        constexpr tally& operator += (tally const& other) noexcept
        {
            wins += other.wins;
            draws += other.draws;
            losses += other.losses;
            return *this;
        }
    };

    // Écart d'Elo qui correspond à un score moyen. Borné pour rester fini à 0 % et 100 %.
    inline double elo_from_score(double score) noexcept
    {
        double const s = std::clamp(score, 1e-4, 1 - 1e-4);
        return -400 * std::log10(1 / s - 1);
    }

    inline double score_from_elo(double elo) noexcept {return 1 / (1 + std::pow(10.0, -elo / 400));}

    // Un écart Elo arrondi, toujours signé : « +0 » plutôt que « -0 ».
    inline std::string signed_elo(double elo) {return std::abs(elo) < 0.5 ? "+0" : std::format("{:+.0f}", elo);}

    // Écart d'Elo entre deux agents, avec son intervalle de confiance à 95 %, d'après leurs parties.
    struct elo_estimate
    {
        double elo, low, high;
    };

    // Si un agent gagne toutes les parties, l'écart n'est pas mesurable : il est infini, et seule sa borne
    // basse a un sens, celle du score p tel que gagner n parties sur n arrive avec 2,5 % de chances (p^n).
    inline elo_estimate estimate_elo(tally const& t) noexcept
    {
        constexpr double infinity = std::numeric_limits<double>::infinity();
        double const n = double(t.games());
        if(n == 0) return {0, -infinity, infinity};
        if(t.wins == t.games()) return {infinity, elo_from_score(std::pow(0.025, 1 / n)), infinity};
        if(t.losses == t.games()) return {-infinity, -infinity, -elo_from_score(std::pow(0.025, 1 / n))};
        double const s = t.score();
        // Variance du score d'une partie, puis erreur type de la moyenne.
        double const variance = (double(t.wins) * (1 - s) * (1 - s) + double(t.draws) * (0.5 - s) * (0.5 - s)
                                 + double(t.losses) * s * s) / n;
        double const error = std::sqrt(variance / n);
        return {elo_from_score(s), elo_from_score(s - 1.96 * error), elo_from_score(s + 1.96 * error)};
    }

    // Classement de plusieurs agents d'après leurs bilans deux à deux (`results[i][j]` : bilan de i face à
    // j), par le modèle de Bradley-Terry, qui est le modèle de l'Elo. Calculé par l'algorithme MM de Hunter
    // (2004), répété jusqu'à convergence. Chaque paire reçoit en plus une partie nulle fictive, pour qu'un
    // agent qui gagne tout garde une note finie. Le premier agent sert de référence : sa note vaut 0.
    // `results[j][i]` doit être le bilan inverse de `results[i][j]`.
    inline std::vector<double> bradley_terry(std::vector<std::vector<tally>> const& results, int max_iterations = 100000)
    {
        std::size_t const n = results.size();
        for(std::size_t i = 0; i < n; ++i)
        {
            assert(results[i].size() == n);
            for(std::size_t j = 0; j < n; ++j) assert(results[i][j] == results[j][i].reversed());
        }
        std::vector<double> strength(n, 1.0);
        auto const points = [&](std::size_t i, std::size_t j) {return results[i][j].points() + (i != j ? 0.5 : 0.0);};
        auto const games = [&](std::size_t i, std::size_t j) {return double(results[i][j].games()) + (i != j ? 1.0 : 0.0);};

        for(int it = 0; it < max_iterations; ++it)
        {
            std::vector<double> next(n);
            for(std::size_t i = 0; i < n; ++i)
            {
                double won = 0, expected = 0;
                for(std::size_t j = 0; j < n; ++j)
                {
                    if(i == j) continue;
                    won += points(i, j);
                    expected += games(i, j) / (strength[i] + strength[j]);
                }
                next[i] = expected > 0 ? won / expected : strength[i];
            }
            // Le modèle ne dépend que des rapports de forces : on ramène la référence à 1, en divisant toutes
            // les forces par la même valeur.
            double const reference = next[0];
            double change = 0;
            for(std::size_t i = 0; i < n; ++i)
            {
                next[i] /= reference;
                change = std::max(change, std::abs(std::log(next[i] / strength[i])));
            }
            strength = next;
            if(change < 1e-10) break;
        }

        std::vector<double> elo(n);
        for(std::size_t i = 0; i < n; ++i) elo[i] = 400 * std::log10(strength[i]);
        return elo;
    }

    // Intervalles à 95 % du classement, par rééchantillonnage des parties de chaque paire. Pour tirer les
    // parties rééchantillonnées, chaque issue (victoire, nul, défaite) reçoit un tiers de partie fictive : un
    // bilan parfait garde ainsi un peu d'incertitude, au lieu d'être reproduit à l'identique.
    //
    // Comme les notes, les intervalles sont relatifs au premier agent : celui d'un agent est l'intervalle de
    // son écart avec la référence (qui a donc [0, 0]), et comprend l'incertitude sur la référence elle-même.
    struct rating
    {
        double elo, low, high;
    };

    inline std::vector<rating> rate(std::vector<std::vector<tally>> const& results, unsigned seed = 1, int samples = 200)
    {
        assert(samples >= 40);      // au moins un tirage de chaque côté des quantiles à 2,5 % et 97,5 %
        std::size_t const n = results.size();
        std::vector<double> const central = bradley_terry(results);
        std::vector<std::vector<double>> resampled_elos(n);      // les notes de chaque agent, tirage par tirage
        std::mt19937 gen(seed);

        for(int sample = 0; sample < samples; ++sample)
        {
            std::vector<std::vector<tally>> resampled(n, std::vector<tally>(n));
            for(std::size_t i = 0; i < n; ++i)
                for(std::size_t j = i + 1; j < n; ++j)
                {
                    tally const& t = results[i][j];
                    std::size_t const g = t.games();
                    if(g == 0) continue;
                    std::discrete_distribution<int> outcome({double(t.wins) + 1.0 / 3, double(t.draws) + 1.0 / 3, double(t.losses) + 1.0 / 3});
                    tally r;
                    for(std::size_t k = 0; k < g; ++k)
                    {
                        int const o = outcome(gen);
                        if(o == 0) ++r.wins;
                        else if(o == 1) ++r.draws;
                        else ++r.losses;
                    }
                    resampled[i][j] = r;
                    resampled[j][i] = r.reversed();
                }
            std::vector<double> const elo = bradley_terry(resampled);
            for(std::size_t i = 0; i < n; ++i) resampled_elos[i].push_back(elo[i]);
        }

        std::vector<rating> out(n);
        for(std::size_t i = 0; i < n; ++i)
        {
            std::ranges::sort(resampled_elos[i]);
            std::size_t const low = std::size_t(0.025 * double(samples)), high = std::size_t(0.975 * double(samples)) - 1;
            out[i] = {central[i], std::min(central[i], resampled_elos[i][low]), std::max(central[i], resampled_elos[i][high])};
        }
        return out;
    }

}   // namespace wsl
