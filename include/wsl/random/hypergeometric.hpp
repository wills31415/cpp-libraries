#pragma once

// Loi hypergéométrique : nombre de succès en n tirages sans remise dans une population de
// N éléments dont K succès : P(X = k) = C(K, k) C(N − K, n − k)/C(N, n),
// max(0, n − (N − K)) <= k <= min(n, K).
//
//   wsl::hypergeometric X({.population = 50, .successes = 20, .draws = 10});
//   pmf(X, 4); cdf(X, 4); quantile(X, 0.5); sample(X, g);
//
// Masse (Loader, comme dhyper de R) : avec π = n/N en double-double, le produit de deux
// masses binomiales b(k; K, π) b(n − k; N − K, π) divisé par b(n; N, π), en logarithme
// (δ, bd0), un arrondi. Répartition : la queue la plus petite (sous ou au-dessus du mode),
// sommée depuis k par le rapport de masses consécutives, en double-double (O(σ) termes) ;
// l'autre par complément. Quantile exact (detail/discrete.hpp).
//
// Tirage (approché, bit-exact) : inversion depuis le mode (comparaisons exactes de
// l'uniforme aux sommes partielles en double-double), O(σ) en moyenne.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/binomial.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/discrete.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>

namespace wsl
{
    template <law_integer I = std::int64_t>
    class hypergeometric
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            I population, successes, draws;
        };

        constexpr hypergeometric(parameters q) : N_(q.population), K_(q.successes), n_(q.draws)
        {
            if (N_ < 0 || K_ < 0 || n_ < 0 || K_ > N_ || n_ > N_ || static_cast<std::uint64_t>(N_) > (std::uint64_t{1} << 53))
                throw std::domain_error("wsl::hypergeometric : il faut 0 <= successes, draws <= population <= 2^53");
        }

        constexpr I population() const noexcept { return N_; }
        constexpr I successes() const noexcept { return K_; }
        constexpr I draws() const noexcept { return n_; }
        constexpr I low() const noexcept { return n_ - (N_ - K_) > 0 ? n_ - (N_ - K_) : I(0); }
        constexpr I high() const noexcept { return n_ < K_ ? n_ : K_; }
        constexpr interval<I> support() const noexcept { return {low(), high()}; }
        // n K/N
        constexpr double mean() const noexcept
        {
            if (n_ == 0 || K_ == 0)
                return 0;
            auto const q = detail::sd::quotient_dd(exact_product(n_, K_), {double(N_), 0.0});
            return q.hi + q.lo;
        }
        // n K (N − K)(N − n)/(N² (N − 1))
        constexpr double variance() const noexcept
        {
            if (n_ == 0 || K_ == 0 || K_ == N_ || n_ == N_)
                return 0;
            auto const a = exact_product(n_, K_), b = exact_product(N_ - K_, N_ - n_);
            detail::scaled const num = detail::normalized(a.hi, a.lo) * detail::normalized(b.hi, b.lo);
            return detail::rounded(num / double(N_) / double(N_) / double(N_ - 1));
        }
        // ⌊(n + 1)(K + 1)/(N + 2)⌋, exact
        constexpr I mode() const noexcept
        {
            auto const num = exact_product(n_ + 1, K_ + 1);
            double const d = double(N_ + 2);
            I f = static_cast<I>((num.hi + num.lo) / d);
            auto const le = [&](I c) {                                      // c d <= num, exact
                auto const t = exact_product(c, N_ + 2);
                return t.hi < num.hi || (t.hi == num.hi && t.lo <= num.lo);
            };
            while (f > 0 && !le(f))
                --f;
            while (le(f + 1))
                ++f;
            return f;
        }

        friend constexpr double pmf(hypergeometric const& X, I k) noexcept
        {
            if (k < X.low() || k > X.high())
                return 0;
            if (X.low() == X.high())
                return 1;
            auto const l = X.log_mass(k);
            if (!(l.hi > -1500))
                return 0;
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }
        friend constexpr double log_pmf(hypergeometric const& X, I k) noexcept
        {
            if (k < X.low() || k > X.high())
                return -std::numeric_limits<double>::infinity();
            if (X.low() == X.high())
                return 0;
            auto const l = X.log_mass(k);
            return l.hi + l.lo;
        }
        friend constexpr double cdf(hypergeometric const& X, I k) noexcept
        {
            if (k < X.low())
                return 0;
            if (k >= X.high())
                return 1;
            auto const [v, is_f] = X.distribution(k);
            return is_f ? detail::rounded(v) : detail::sd::one_minus(v);
        }
        friend constexpr double ccdf(hypergeometric const& X, I k) noexcept
        {
            if (k < X.low())
                return 1;
            if (k >= X.high())
                return 0;
            auto const [v, is_f] = X.distribution(k);
            return is_f ? detail::sd::one_minus(v) : detail::rounded(v);
        }

        // Plus petit k tel que P(X <= k) >= u.
        friend constexpr I quantile(hypergeometric const& X, double u)
        {
            if (!(u >= 0 && u <= 1))
                throw std::domain_error("wsl::quantile : p doit être dans [0, 1]");
            if (u == 0 || X.low() == X.high())
                return X.low();
            if (u == 1)
                return X.high();
            auto const ge = [&](I k) {
                if (k >= X.high())
                    return true;
                auto const [v, is_f] = X.distribution(k);
                return detail::at_least(v, is_f, u);
            };
            return detail::smallest_true(X.mode(), X.low(), X.high(), ge);
        }

        // Inversion depuis le mode : u comparé exactement aux sommes partielles de masses
        // (double-double), en s'éloignant du mode de part et d'autre.
        template <std::uniform_random_bit_generator G>
        friend constexpr I sample(hypergeometric const& X, G& g)
        {
            if (X.low() == X.high())
                return X.low();
            double const u = uniform_real<double>(g);
            I const m = X.mode() < X.low() ? X.low() : X.mode() > X.high() ? X.high() : X.mode();
            // F(m − 1) et masses vers le bas, puis vers le haut
            auto const [v, is_f] = X.distribution(m > X.low() ? m - 1 : m);
            detail::sd::double_double f_below{0.0, 0.0};                   // F(m − 1)
            if (m > X.low())
            {
                auto const w = detail::unscaled(v);
                f_below = is_f ? w : detail::sum({1.0, 0.0}, {-w.hi, -w.lo});
            }
            auto const lm = X.log_mass(m);
            detail::sd::double_double const pm = detail::unscaled(detail::exp_scaled(lm.hi, lm.lo));
            auto const below = [u](detail::sd::double_double f) { return f.hi > u || (f.hi == u && f.lo > 0); };   // u < f
            if (below(f_below))                                             // u < F(m − 1) : vers le bas
            {
                detail::sd::double_double f = f_below, mass = pm;
                for (I k = m - 1; k > X.low(); --k)
                {
                    mass = detail::mul(mass, X.ratio_down(k + 1));          // P(k)
                    f = detail::sum(f, {-mass.hi, -mass.lo});               // F(k − 1)
                    if (!below(f))
                        return k;
                }
                return X.low();
            }
            detail::sd::double_double f = detail::sum(f_below, pm), mass = pm;   // F(m)
            for (I k = m; k < X.high(); ++k)
            {
                if (below(f))
                    return k;
                mass = detail::mul(mass, X.ratio_up(k));                    // P(k + 1)
                f = detail::sum(f, mass);
            }
            return X.high();
        }

        friend constexpr bool operator==(hypergeometric const&, hypergeometric const&) = default;

    private:
        // log b(x; n, π) (masse binomiale à π = n/N), m = nπ, n − m en double-double.
        static constexpr detail::sd::double_double binomial_log(double x, double n, detail::sd::double_double m, detail::sd::double_double rest) noexcept
        {
            namespace sd = detail::sd;
            if (n == 0)
                return {0.0, 0.0};
            if (x == 0)                                                     // n log((n − m)/n)
                return detail::mul({n, 0.0}, detail::log(detail::normalized(rest.hi, rest.lo) / n));
            if (x == n)
                return detail::mul({n, 0.0}, detail::log(detail::normalized(m.hi, m.lo) / n));
            auto e = detail::stirling_error(n);
            e = detail::sum(e, detail::negated(detail::stirling_error(x)));
            e = detail::sum(e, detail::negated(detail::stirling_error(n - x)));
            e = detail::sum(e, detail::negated(detail::deviance(x, m.hi, m.lo)));
            e = detail::sum(e, detail::negated(detail::deviance(n - x, rest.hi, rest.lo)));
            detail::scaled const r = detail::normalized(n, 0.0) / (detail::normalized(x, 0.0) * detail::normalized(n - x, 0.0));
            auto const l = detail::log(r);                                  // ½ log(n/(2π x (n − x)))
            e = detail::sum(e, {sd::fmul(l.hi, 0.5), sd::fmul(l.lo, 0.5)});
            return detail::sum(e, {-detail::half_log_2pi_hi, -detail::half_log_2pi_lo});
        }

        // log P(X = k) = log b(k; K, π) + log b(n − k; N − K, π) − log b(n; N, π), π = n/N.
        constexpr detail::sd::double_double log_mass(I k) const noexcept
        {
            double const N = double(N_), K = double(K_), n = double(n_), x = double(k);
            auto const frac = [&](I c) {                                    // c π et c (1 − π), c entier
                return std::pair{detail::sd::quotient_dd(exact_product(c, n_), {N, 0.0}), detail::sd::quotient_dd(exact_product(c, N_ - n_), {N, 0.0})};
            };
            auto const [k1, k2] = frac(K_);
            auto const [j1, j2] = frac(N_ - K_);
            auto e = binomial_log(x, K, k1, k2);
            e = detail::sum(e, binomial_log(n - x, N - K, j1, j2));
            auto const total = binomial_log(n, N, {n, 0.0}, {N - n, 0.0});
            return detail::sum(e, detail::negated(total));
        }

        // P(k − 1)/P(k) = k (N − K − n + k)/((K − k + 1)(n − k + 1))
        constexpr detail::sd::double_double ratio_down(I k) const noexcept
        {
            return detail::sd::quotient_dd(exact_product(k, N_ - K_ - n_ + k), exact_product(K_ - k + 1, n_ - k + 1));
        }
        // P(k + 1)/P(k) = (K − k)(n − k)/((k + 1)(N − K − n + k + 1))
        constexpr detail::sd::double_double ratio_up(I k) const noexcept
        {
            return detail::sd::quotient_dd(exact_product(K_ - k, n_ - k), exact_product(k + 1, N_ - K_ - n_ + k + 1));
        }

        // a b exact en double-double (a, b <= 2^53).
        static constexpr detail::sd::double_double exact_product(I a, I b) noexcept
        {
            auto const [h, l] = detail::sd::two_product(double(a), double(b));
            return {h, l};
        }

        // La plus petite part : F(k) (is_f) si k est sous le mode, sinon 1 − F(k) ; low <= k < high.
        struct part
        {
            detail::scaled value;
            bool is_f;
        };
        constexpr part distribution(I k) const noexcept
        {
            bool const lower = k < mode();
            I j = lower ? k : k + 1;
            auto const l = log_mass(j);
            if (!(l.hi > -1500))
                return {{0.0, 0.0, 0}, lower};
            detail::sd::double_double total{1.0, 0.0}, term{1.0, 0.0};      // Σ P(·)/P(j)
            for (;;)
            {
                if (lower ? j <= low() : j >= high())
                    break;
                term = detail::mul(term, lower ? ratio_down(j) : ratio_up(j));
                j = lower ? j - 1 : j + 1;
                total = detail::sum(total, term);
                if (detail::magnitude(term.hi) < detail::sd::fmul(total.hi, 0x1p-112))
                    break;
            }
            return {detail::exp_scaled(l.hi, l.lo) * detail::normalized(total.hi, total.lo), lower};
        }

        I N_, K_, n_;
    };
}
