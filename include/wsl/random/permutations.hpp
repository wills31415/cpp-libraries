#pragma once

// Lois sur les permutations de {0, …, n − 1} (valeurs wsl::permutation).
//
//   wsl::uniform_permutation U({.size = 10});
//   wsl::ewens E({.size = 10, .theta = 2.5});      // P(σ) = θ^{cycles(σ)} / θ(θ + 1)…(θ + n − 1)
//   wsl::mallows M({.size = 10, .q = 0.7});        // P(σ) = q^{inversions(σ)} / Π [i]_q
//   sample(U, g); pmf(E, sigma); log_pmf(M, sigma); E.mean_cycles(); M.mean_inversions();
//
// Tirages exacts (la permutation tirée a exactement la probabilité annoncée) :
// - uniforme : Fisher–Yates (n − 1 entiers uniformes) ;
// - Ewens : restaurant chinois — l'élément i ouvre un cycle avec probabilité θ/(θ + i),
//   sinon s'insère après un élément uniforme des i premiers ; θ/(θ + i) décidé exactement
//   (detail/choice.hpp) ;
// - Mallows : code d'insertion V_i ∈ {0, …, i − 1}, P(V_i = v) ∝ q^v, par épreuves de
//   Bernoulli exactes (géométrique tronquée si i(1 − q) >= 1, sinon uniforme acceptée
//   avec probabilité q^v ; q > 1 par symétrie, de rapport 1/q) ; coût O(min(i, 1/|1 − q|))
//   épreuves par position.
// Masses en logarithme double-double, un arrondi (0 si elles sous-débordent).

#include <wsl/permutation.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/choice.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace wsl
{
    namespace detail
    {
        constexpr double exp_of_dd(sd::double_double l) noexcept
        {
            if (!(l.hi > -1500))
                return l.hi == l.hi ? 0.0 : l.hi;
            if (l.hi > 1400)
                return std::numeric_limits<double>::infinity();
            return rounded(exp_scaled(l.hi, l.lo));
        }

        // log Γ(y), y > 0, en double-double (y < 1 : lgamma(1 + y) − log y).
        constexpr sd::double_double log_gamma_any(double y) noexcept
        {
            if (y >= 1)
                return log_gamma(y);
            auto const [s_hi, s_lo] = sd::two_sum(1.0, y);
            auto const l = y <= 0.25 ? sd::lgamma1p_small(y) : log_gamma(s_hi, s_lo);
            auto const ly = log_precise(y);
            return sum(l, {-ly.hi, -ly.lo});
        }

        // log n!
        constexpr sd::double_double log_factorial(std::uint64_t n) noexcept
        {
            return n < 2 ? sd::double_double{0.0, 0.0} : log_gamma(double(n) + 1);
        }

        // 1 − e^x, x < 0 en double-double.
        constexpr sd::double_double one_minus_exp(sd::double_double x) noexcept
        {
            if (x.hi > -40)
            {
                auto const m = expm1_double_double(x.hi, x.lo);
                return {-m.hi, -m.lo};
            }
            if (x.hi < -1500)
                return {1.0, 0.0};
            auto const e = unscaled(exp_scaled(x.hi, x.lo));
            return sum({1.0, 0.0}, {-e.hi, -e.lo});
        }
    }

    struct permutation_size
    {
        std::size_t size;
    };

    class uniform_permutation
    {
    public:
        using value_type = permutation;
        using parameters = permutation_size;
        static constexpr exactness sampling = exactness::exact;

        constexpr explicit uniform_permutation(parameters p) : n_(p.size) {}

        constexpr std::size_t size() const noexcept { return n_; }
        constexpr double entropy() const noexcept
        {
            auto const l = detail::log_factorial(n_);
            return l.hi + l.lo;
        }

        friend constexpr double log_pmf(uniform_permutation const& U, permutation const& s) noexcept
        {
            return s.size() == U.n_ ? -U.entropy() : -std::numeric_limits<double>::infinity();
        }
        friend constexpr double pmf(uniform_permutation const& U, permutation const& s) noexcept
        {
            if (s.size() != U.n_)
                return 0;
            auto const l = detail::log_factorial(U.n_);
            return detail::exp_of_dd({-l.hi, -l.lo});
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr permutation sample(uniform_permutation const& U, G& g)
        {
            std::vector<std::size_t> images(U.n_);
            for (std::size_t i = 0; i < U.n_; ++i)
                images[i] = i;
            shuffle(images, g);
            return permutation(std::move(images));
        }

        friend constexpr bool operator==(uniform_permutation const&, uniform_permutation const&) = default;

    private:
        std::size_t n_;
    };

    class ewens
    {
    public:
        using value_type = permutation;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t size;
            double theta;
        };

        constexpr explicit ewens(parameters p) : n_(p.size), theta_(p.theta)
        {
            if (!(theta_ > 0) || !detail::is_finite(theta_))
                throw std::domain_error("wsl::ewens : il faut theta > 0, fini");
            // log θ(θ + 1)…(θ + n − 1) = log Γ(θ + n) − log Γ(θ)
            if (n_ > 0)
            {
                auto const [s_hi, s_lo] = detail::sd::two_sum(theta_, double(n_));
                auto const a = detail::log_gamma(s_hi, s_lo), b = detail::log_gamma_any(theta_);
                log_rising_ = detail::sum(a, {-b.hi, -b.lo});
            }
            log_theta_ = detail::log_precise(theta_);
        }

        constexpr std::size_t size() const noexcept { return n_; }
        constexpr double theta() const noexcept { return theta_; }
        // E[nombre de cycles] = Σ_{i<n} θ/(θ + i)
        constexpr double mean_cycles() const noexcept
        {
            detail::sd::double_double s{0.0, 0.0};
            for (std::size_t i = 0; i < n_; ++i)
            {
                auto const [d_hi, d_lo] = detail::sd::two_sum(theta_, double(i));
                s = detail::sum(s, detail::sd::quotient_dd({theta_, 0.0}, {d_hi, d_lo}));
            }
            return s.hi + s.lo;
        }

        // k log θ − log θ^{(n)}
        friend constexpr double log_pmf(ewens const& E, permutation const& s) noexcept
        {
            auto const l = E.log_mass(s);
            return l.hi + l.lo;
        }
        friend constexpr double pmf(ewens const& E, permutation const& s) noexcept { return detail::exp_of_dd(E.log_mass(s)); }

        template <std::uniform_random_bit_generator G>
        friend constexpr permutation sample(ewens const& E, G& g)
        {
            std::vector<std::size_t> images(E.n_);
            for (std::size_t i = 0; i < E.n_; ++i)
            {
                if (i == 0 || detail::choose(g, E.theta_, double(i)))       // nouveau cycle : θ/(θ + i)
                    images[i] = i;
                else                                                        // après j : i → σ(j), j → i
                {
                    auto const j = static_cast<std::size_t>(uniform_below(g, i));
                    images[i] = images[j];
                    images[j] = i;
                }
            }
            return permutation(std::move(images));
        }

        friend constexpr bool operator==(ewens const& a, ewens const& b) noexcept { return a.n_ == b.n_ && a.theta_ == b.theta_; }

    private:
        constexpr detail::sd::double_double log_mass(permutation const& s) const noexcept
        {
            if (s.size() != n_)
                return {-std::numeric_limits<double>::infinity(), 0.0};
            auto const l = detail::mul({double(s.cycle_count()), 0.0}, log_theta_);
            return detail::sum(l, {-log_rising_.hi, -log_rising_.lo});
        }

        std::size_t n_;
        double theta_;
        detail::sd::double_double log_rising_{}, log_theta_{};
    };

    class mallows
    {
    public:
        using value_type = permutation;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t size;
            double q;
        };

        constexpr explicit mallows(parameters p) : n_(p.size), q_(p.q)
        {
            if (!(q_ > 0) || !detail::is_finite(q_))
                throw std::domain_error("wsl::mallows : il faut q > 0, fini");
            log_q_ = detail::log_precise(q_);
            // log Z = Σ_{i=1}^{n} log [i]_q
            for (std::size_t i = 1; i <= n_; ++i)
                log_z_ = detail::sum(log_z_, log_q_integer(i));
        }

        constexpr std::size_t size() const noexcept { return n_; }
        constexpr double q() const noexcept { return q_; }
        // E[inversions] = Σ_i (Σ_v v q^v)/[i]_q
        constexpr double mean_inversions() const noexcept
        {
            detail::sd::double_double s{0.0, 0.0};
            for (std::size_t i = 2; i <= n_; ++i)
                s = detail::sum(s, mean_code(i));
            return s.hi + s.lo;
        }

        friend constexpr double log_pmf(mallows const& M, permutation const& s) noexcept
        {
            auto const l = M.log_mass(s);
            return l.hi + l.lo;
        }
        friend constexpr double pmf(mallows const& M, permutation const& s) noexcept { return detail::exp_of_dd(M.log_mass(s)); }

        template <std::uniform_random_bit_generator G>
        friend constexpr permutation sample(mallows const& M, G& g)
        {
            // La valeur i s'insère en position i − V_i (V_i inversions avec les plus petites).
            std::vector<std::size_t> line;
            line.reserve(M.n_);
            for (std::size_t i = 0; i < M.n_; ++i)
            {
                std::size_t const v = M.code(g, i + 1);
                line.insert(line.begin() + static_cast<std::ptrdiff_t>(i - v), i);
            }
            return permutation(std::move(line));
        }

        friend constexpr bool operator==(mallows const& a, mallows const& b) noexcept { return a.n_ == b.n_ && a.q_ == b.q_; }

    private:
        // V ∈ {0, …, m − 1}, P(V = v) ∝ q^v, exactement.
        template <class G>
        constexpr std::size_t code(G& g, std::size_t m) const
        {
            if (m == 1)
                return 0;
            if (q_ == 1)
                return static_cast<std::size_t>(uniform_below(g, m));
            bool const below = q_ < 1;                                      // rapport r = q, ou 1/q (V renversé)
            auto const success = [&] { return below ? bernoulli_trial(g, q_) : detail::fraction(g, 1.0, q_); };
            double const gap = below ? 1 - q_ : 1 - 1 / q_;                 // 1 − r (seulement pour choisir la méthode)
            std::size_t w = 0;
            if (double(m) * gap >= 1)                                       // géométrique tronquée
                for (;;)
                {
                    w = 0;
                    while (w < m && success())
                        ++w;
                    if (w < m)
                        break;
                }
            else                                                            // uniforme acceptée avec probabilité r^w
                for (;;)
                {
                    w = static_cast<std::size_t>(uniform_below(g, m));
                    std::size_t k = 0;
                    while (k < w && success())
                        ++k;
                    if (k == w)
                        break;
                }
            return below ? w : m - 1 - w;
        }

        // log [i]_q = log((1 − q^i)/(1 − q)) ; q > 1 : (i − 1) log q + log [i]_{1/q}.
        constexpr detail::sd::double_double log_q_integer(std::size_t i) const noexcept
        {
            if (i == 1)
                return {0.0, 0.0};
            if (q_ == 1)
                return detail::log_precise(double(i));
            detail::sd::double_double const l = q_ < 1 ? log_q_ : detail::sd::double_double{-log_q_.hi, -log_q_.lo};   // < 0
            auto const a = detail::one_minus_exp(detail::mul({double(i), 0.0}, l)), b = detail::one_minus_exp(l);
            auto r = detail::sum(detail::log(detail::normalized(a.hi, a.lo)), [&] { auto const lb = detail::log(detail::normalized(b.hi, b.lo)); return detail::sd::double_double{-lb.hi, -lb.lo}; }());
            if (q_ > 1)
                r = detail::sum(r, detail::mul({double(i - 1), 0.0}, log_q_));
            return r;
        }

        // E[V_i] = Σ_{v<i} v q^v / [i]_q
        constexpr detail::sd::double_double mean_code(std::size_t i) const noexcept
        {
            if (q_ == 1)
                return {detail::sd::fmul(double(i - 1), 0.5), 0.0};
            // r/(1 − r) − i rⁱ/(1 − rⁱ), r = q < 1 ; q > 1 : (i − 1) − celle de r = 1/q (en double)
            bool const below = q_ < 1;
            double const r = below ? q_ : 1 / q_;
            double const ri = wsl::special::exp(detail::sd::fmul(double(i), below ? log_q_.hi : -log_q_.hi));
            double const v = r / (1 - r) - detail::sd::fmul(double(i), ri) / (1 - ri);
            return {below ? v : double(i - 1) - v, 0.0};
        }

        constexpr detail::sd::double_double log_mass(permutation const& s) const noexcept
        {
            if (s.size() != n_)
                return {-std::numeric_limits<double>::infinity(), 0.0};
            auto const l = detail::mul({double(s.inversions()), 0.0}, log_q_);
            return detail::sum(l, {-log_z_.hi, -log_z_.lo});
        }

        std::size_t n_;
        double q_;
        detail::sd::double_double log_q_{}, log_z_{0.0, 0.0};
    };
}
