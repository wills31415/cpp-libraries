#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Loi bêta-binomiale : binomiale de n épreuves dont le paramètre suit une loi bêta(α, β),
// P(X = k) = C(n, k) B(k + α, n − k + β)/B(α, β), 0 <= k <= n.
//
//   wsl::probability::beta_binomial X({.trials = 20, .alpha = 2, .beta = 3});
//   X.pmf(6); X.cdf(6); X.quantile(0.9); X.variance(); X(g);
//
// Masse : somme des huit log Γ (double-double, wsl::special), les termes −y de Stirling
// s'annulant exactement ; erreur absolue du logarithme ~ 2^-106 n log n (0,5 ulp mesuré
// jusqu'à n = 10^6). Répartition : la queue du côté de k (sous ou au-dessus de la moyenne),
// sommée depuis k par le rapport de masses consécutives, en double-double — arrêt anticipé
// quand α, β >= 1 (masse log-concave), sinon jusqu'au bord : O(n) ; l'autre par complément.
// Quantile exact (detail/discrete.hpp).
//
// Tirage (approché, bit-exact) : P ~ bêta(α, β), puis X ~ binomiale(n, P).

#include <wsl/probability/laws/basic.hpp>
#include <wsl/probability/laws/beta.hpp>
#include <wsl/probability/laws/binomial.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/discrete.hpp>
#include <wsl/probability/detail/loader.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/laws/permutations.hpp>
#include <wsl/special.hpp>

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_integer I = std::int64_t>
    class beta_binomial : public law_interface
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            I trials;
            double alpha, beta;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.trials = trials(), .alpha = alpha(), .beta = beta()}; }

        constexpr beta_binomial(parameters q) : n_(q.trials), a_(q.alpha), b_(q.beta)
        {
            if (n_ < 0 || static_cast<std::uint64_t>(n_) > (std::uint64_t{1} << 53) || !(a_ > 0) || !(b_ > 0) || !detail::is_finite(a_) || !detail::is_finite(b_))
                throw std::domain_error("wsl::probability::beta_binomial : il faut 0 <= trials <= 2^53, alpha > 0 et beta > 0 finis");
        }

        constexpr I trials() const noexcept { return n_; }
        constexpr double alpha() const noexcept { return a_; }
        constexpr double beta() const noexcept { return b_; }
        constexpr interval<I> support() const noexcept { return {I(0), n_}; }
        // n α/(α + β)
        constexpr double mean() const noexcept
        {
            if (n_ == 0)
                return 0;
            auto const [s, s_error] = detail::sd::two_sum(a_, b_);
            return detail::rounded(detail::normalized(double(n_), 0.0) * detail::normalized(a_, 0.0) / detail::normalized(s, s_error));
        }
        // n α β (α + β + n)/((α + β)² (α + β + 1))
        constexpr double variance() const noexcept
        {
            if (n_ == 0)
                return 0;
            auto const s = detail::sum({a_, 0.0}, {b_, 0.0});
            auto const s1 = detail::sum(s, {1.0, 0.0}), sn = detail::sum(s, {double(n_), 0.0});
            detail::scaled const num = detail::normalized(double(n_), 0.0) * detail::normalized(a_, 0.0) * detail::normalized(b_, 0.0) * detail::normalized(sn.hi, sn.lo);
            detail::scaled const ss = detail::normalized(s.hi, s.lo);
            return detail::rounded(num / (ss * ss * detail::normalized(s1.hi, s1.lo)));
        }

        constexpr double pmf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 0 || k > X.n_)
                return 0;
            auto const l = X.log_mass(k);
            if (!(l.hi > -1500))
                return 0;
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }
        constexpr double log_pmf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 0 || k > X.n_)
                return -std::numeric_limits<double>::infinity();
            auto const l = X.log_mass(k);
            return l.hi + l.lo;
        }
        constexpr double cdf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 0)
                return 0;
            if (k >= X.n_)
                return 1;
            auto const [v, is_f] = X.distribution(k);
            return is_f ? detail::rounded(v) : detail::sd::one_minus(v);
        }
        constexpr double ccdf(I k) const noexcept
        {
            auto const& X = *this;
            if (k < 0)
                return 1;
            if (k >= X.n_)
                return 0;
            auto const [v, is_f] = X.distribution(k);
            return is_f ? detail::sd::one_minus(v) : detail::rounded(v);
        }

        // Plus petit k tel que P(X <= k) >= u.
        constexpr I quantile(double u) const
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            if (!(u >= 0 && u <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            if (u == 0 || X.n_ == 0)
                return 0;
            if (u == 1)
                return X.n_;
            double guess = X.mean() + sd::fmul(detail::normal_quantile_guess(u), wsl::special::sqrt(X.variance()));
            guess = guess < 0 ? 0 : guess > double(X.n_) ? double(X.n_) : guess;
            auto const ge = [&](I k) {
                if (k >= X.n_)
                    return true;
                auto const [v, is_f] = X.distribution(k);
                return detail::at_least(v, is_f, u);
            };
            return detail::smallest_true(static_cast<I>(guess), I(0), X.n_, ge);
        }

        template <std::uniform_random_bit_generator G>
        constexpr I sample(G& g) const
        {
            auto const& X = *this;
            double const p = detail::native_sample(wsl::probability::beta<double>({.alpha = X.a_, .beta = X.b_}), g);
            return detail::native_sample(binomial<I>({.trials = X.n_, .p = p}), g);
        }

        friend constexpr bool operator==(beta_binomial const& a, beta_binomial const& b) noexcept { return a.n_ == b.n_ && a.a_ == b.a_ && a.b_ == b.b_; }

    private:
        static constexpr detail::sd::double_double lg(double y) noexcept { return detail::log_gamma_any(y); }
        static constexpr detail::sd::double_double lg(detail::sd::double_double y) noexcept
        {
            return y.lo == 0 ? detail::log_gamma_any(y.hi) : detail::log_gamma(y.hi, y.lo);
        }

        // log Γ(n+1) + log Γ(k+α) + log Γ(n−k+β) + log Γ(α+β) − log Γ(k+1) − log Γ(n−k+1) − log Γ(n+α+β) − log Γ(α) − log Γ(β)
        constexpr detail::sd::double_double log_mass(I k) const noexcept
        {
            double const n = double(n_), x = double(k), y = double(n_ - k);
            auto const xa = detail::sum({x, 0.0}, {a_, 0.0}), yb = detail::sum({y, 0.0}, {b_, 0.0});
            auto const ab = detail::sum({a_, 0.0}, {b_, 0.0}), nab = detail::sum(ab, {n, 0.0});
            auto e = detail::sum(lg(n + 1), lg(xa));
            e = detail::sum(e, lg(yb));
            e = detail::sum(e, lg(ab));
            e = detail::sum(e, detail::negated(lg(x + 1)));
            e = detail::sum(e, detail::negated(lg(y + 1)));
            e = detail::sum(e, detail::negated(lg(nab)));
            e = detail::sum(e, detail::negated(lg(a_)));
            return detail::sum(e, detail::negated(lg(b_)));
        }

        // P(j − 1)/P(j) = j (n − j + β)/((n − j + 1)(j − 1 + α))
        constexpr detail::sd::double_double ratio_down(I j) const noexcept
        {
            auto const num = detail::mul({double(j), 0.0}, detail::sum({double(n_ - j), 0.0}, {b_, 0.0}));
            auto const den = detail::mul({double(n_ - j + 1), 0.0}, detail::sum({double(j - 1), 0.0}, {a_, 0.0}));
            return detail::sd::quotient_dd(num, den);
        }
        // P(j + 1)/P(j) = (n − j)(j + α)/((j + 1)(n − j − 1 + β))
        constexpr detail::sd::double_double ratio_up(I j) const noexcept
        {
            auto const num = detail::mul({double(n_ - j), 0.0}, detail::sum({double(j), 0.0}, {a_, 0.0}));
            auto const den = detail::mul({double(j + 1), 0.0}, detail::sum({double(n_ - j - 1), 0.0}, {b_, 0.0}));
            return detail::sd::quotient_dd(num, den);
        }

        struct part
        {
            detail::scaled value;
            bool is_f;
        };
        // F(k) si k est sous la moyenne, sinon 1 − F(k) ; 0 <= k < n.
        constexpr part distribution(I k) const noexcept
        {
            bool const lower = double(k) < mean();
            bool const concave = a_ >= 1 && b_ >= 1;
            I j = lower ? k : k + 1;
            auto const l = log_mass(j);
            detail::sd::double_double total{1.0, 0.0}, term{1.0, 0.0};      // Σ P(·)/P(j)
            while (lower ? j > 0 : j < n_)
            {
                auto const r = lower ? ratio_down(j) : ratio_up(j);
                term = detail::mul(term, r);
                j = lower ? j - 1 : j + 1;
                total = detail::sum(total, term);
                if (concave && r.hi < 1 && detail::magnitude(term.hi) < detail::sd::fmul(total.hi, 0x1p-112))
                    break;
            }
            if (!(l.hi > -1500))
                return {{0.0, 0.0, 0}, lower};
            return {detail::exp_scaled(l.hi, l.lo) * detail::normalized(total.hi, total.lo), lower};
        }

        I n_;
        double a_, b_;
    };
}
