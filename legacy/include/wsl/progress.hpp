#pragma once

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <deque>
#include <format>
#include <limits>
#include <print>
#include <string>
#include <utility>

#include <unistd.h>


// Suivi d'une tâche longue, sur la sortie d'erreur : part faite, temps écoulé, vitesse récente, temps restant
// estimé.
//
// Le temps restant suit la vitesse des 10 dernières minutes, et non la moyenne depuis le début : un
// ralentissement passager (un autre programme qui prend la carte graphique, par exemple) ne fausse plus
// l'estimation pendant des heures. Avant 10 minutes, c'est la moyenne depuis le début.
//
// Sur un terminal, une ligne réécrite sur place (au plus quatre fois par seconde). Sinon (sortie redirigée
// vers un fichier, tâche lancée en arrière-plan), une ligne à chaque tranche de 5 % ou toutes les 30 s, pour
// qu'un journal reste lisible ; rien pendant les 2 premières secondes, pour ne pas encombrer de tâche courte. L'appelant protège les appels s'ils viennent de plusieurs fils.

namespace wsl
{

    // Une durée lisible : « 42 s », « 3 min 05 s », « 1 h 07 min ».
    inline std::string readable_duration(double seconds)
    {
        auto const s = static_cast<long long>(seconds + 0.5);
        if (s < 60) return std::format("{} s", s);
        if (s < 3600) return std::format("{} min {:02} s", s / 60, s % 60);
        return std::format("{} h {:02} min", s / 3600, s % 3600 / 60);
    }

    // Une vitesse lisible, en unités par seconde ou par minute : « 12/s », « 0.58/s », « 3.2/min ».
    inline std::string readable_rate(double per_second)
    {
        if (per_second >= 10) return std::format("{:.0f}/s", per_second);
        if (per_second >= 0.1) return std::format("{:.2f}/s", per_second);
        return std::format("{:.1f}/min", 60 * per_second);
    }

    // La vitesse d'avancement sur une fenêtre de temps glissante : l'avancement retenu au plus une fois par
    // seconde, et oublié au-delà de la fenêtre.
    class recent_rate
    {
    public:
        using clock = std::chrono::steady_clock;

        explicit recent_rate(clock::duration window = std::chrono::minutes(10)) noexcept : window_(window) {}

        void reset(clock::time_point start)
        {
            samples_.assign(1, {start, 0});
        }

        // `done` unités faites à l'instant `now`.
        void add(clock::time_point now, std::size_t done)
        {
            if (samples_.empty() || now - samples_.back().time >= std::chrono::seconds(1)) samples_.push_back({now, done});
            else samples_.back().done = done;
            // Le plus ancien point gardé est le dernier d'avant la fenêtre : elle est ainsi toujours couverte.
            while (samples_.size() > 2 && now - samples_[1].time >= window_) samples_.pop_front();
        }

        // Unités par seconde sur la fenêtre (depuis le début si elle est plus courte), jusqu'à `now` ; 0 sans
        // avancement mesurable.
        double per_second(clock::time_point now, std::size_t done) const
        {
            if (samples_.empty()) return 0;
            double const seconds = std::chrono::duration<double>(now - samples_.front().time).count();
            return seconds > 0 && done > samples_.front().done ? double(done - samples_.front().done) / seconds : 0;
        }

    private:
        struct sample
        {
            clock::time_point time;
            std::size_t done;
        };

        clock::duration window_;
        std::deque<sample> samples_;
    };

    class progress_meter
    {
    public:
        using clock = std::chrono::steady_clock;

        explicit progress_meter(std::string label) : label_(std::move(label)), terminal_(isatty(fileno(stderr)) != 0)
        {
            rate_.reset(start_);
        }

        // Reprend le compte du temps (pour une tâche qui commence après la construction).
        void restart()
        {
            start_ = last_ = clock::now();
            last_tranche_ = 0;
            printed_ = false;
            rate_.reset(start_);
        }

        // `done` unités faites sur `total`.
        void operator()(std::size_t done, std::size_t total)
        {
            auto const now = clock::now();
            rate_.add(now, done);
            // tranches de 5 % : 20 done/total, sans débordement (division préalable pour les totaux > 2^59)
            std::size_t const tranche = !total || done >= total ? 20
                                      : total <= std::numeric_limits<std::size_t>::max() / 20 ? done * 20 / total
                                      : done / (total / 20);
            bool const end = done >= total;
            if (terminal_)
            {
                if (!end && now - last_ < std::chrono::milliseconds(250)) return;
                std::print(stderr, "\r\x1b[K{}", line(done, total, now));
            }
            else
            {
                if (end ? !printed_ : now - start_ < std::chrono::seconds(2)) return;
                if (!end && tranche == last_tranche_ && now - last_ < std::chrono::seconds(30)) return;
                std::println(stderr, "{}", line(done, total, now));
                printed_ = true;
            }
            last_ = now;
            last_tranche_ = tranche;
            if (end) clear();
        }

        // Efface la ligne sur place, avant d'écrire autre chose sur le terminal.
        void clear() const
        {
            if (terminal_) std::print(stderr, "\r\x1b[K");
        }

    private:
        std::string line(std::size_t done, std::size_t total, clock::time_point now) const
        {
            double const elapsed = std::chrono::duration<double>(now - start_).count();
            double const part = total ? double(done) / double(total) : 1.0;
            std::string text = std::format("{} : {} sur {} ({:.0f} %), {}", label_, done, total, 100 * part, readable_duration(elapsed));
            if (done > 0 && done < total)
            {
                double const rate = rate_.per_second(now, done);
                if (rate > 0) text += std::format(", {}, encore environ {}", readable_rate(rate), readable_duration(double(total - done) / rate));
            }
            return text;
        }

        std::string label_;
        bool terminal_;
        clock::time_point start_ = clock::now();
        clock::time_point last_ = start_;
        std::size_t last_tranche_ = 0;
        bool printed_ = false;
        recent_rate rate_;
    };

}   // namespace wsl
