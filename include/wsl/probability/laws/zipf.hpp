#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Loi de Zipf (finie) sur {1, …, N} : P(X = k) = k^{−s}/H_{N,s}, H_{N,s} = Σ_{j=1}^{N} j^{−s},
// s >= 0 (s = 0 : uniforme), N <= 2^53.
//
//   wsl::probability::laws::zipf X({.size = 1000, .exponent = 1.2});
//   X.pmf(3); X.cdf(10); X.quantile(0.9); X.mean(); X(g);
//
// Sommes S_p(a, b) = Σ_{j=a}^{b} j^p (p = −s, 1 − s, 2 − s) en double-double : directes pour
// j < a₀ = max(64, 8|p|), puis Euler–Maclaurin à 12 termes de Bernoulli, tout relativement à
// a₀^p (pas de débordement) ; l'intégrale a^{p+1}·expm1((p + 1) log(b/a))/(p + 1) reste
// stable près de s = 1. Répartition H_{k,s}/H_{N,s}, survie S_{−s}(k + 1, N)/H_{N,s} (la queue
// directement, sans différence) ; quantiles exacts (detail/discrete.hpp). Moyenne
// 1 + A/H et variance (B H − A²)/H², A = Σ (k − 1) k^{−s}, B = Σ (k − 1)² k^{−s} (pas de
// compensation H_{s−1} − H_s pour s grand).
//
// Tirage (approché, bit-exact) : rejet par inversion de Hörmann et Derflinger (1996) ;
// uniforme pour s = 0.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/discrete.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/random/detail/portable_math.hpp>

#include <array>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
}

// Outils internes : dans wsl::probability::detail (un seul « detail » visible depuis probability).
namespace wsl::inline WSL_DEFAULTS::probability::detail::zipf_sums
{
        namespace sd = random::detail::sd;
        using dd = sd::double_double;

        // B_{2m}/(2m)!, m = 1 … 12, en double-double (numérateur, dénominateur de B_{2m}).
        constexpr std::array<dd, 12> bernoulli_ratios() noexcept
        {
            constexpr double num[12] = {1, -1, 1, -1, 5, -691, 7, -3617, 43867, -174611, 854513, -236364091};
            constexpr double den[12] = {6, 30, 42, 30, 66, 2730, 6, 510, 798, 330, 138, 2730};
            std::array<dd, 12> out{};
            dd factorial{1.0, 0.0};
            for (int m = 1; m <= 12; ++m)
            {
                factorial = mul(mul(factorial, {double(2 * m - 1), 0.0}), {double(2 * m), 0.0});   // (2m)!, exact
                out[m - 1] = sd::quotient_dd(sd::quotient_dd({num[m - 1], 0.0}, {den[m - 1], 0.0}), factorial);
            }
            return out;
        }
        inline constexpr std::array<dd, 12> bernoulli = bernoulli_ratios();

        // x^p à l'échelle (x >= 1 entier exact), 0 si sous 2^-2000.
        constexpr scaled power(double x, double p) noexcept
        {
            auto const l = mul(log_precise(x), {p, 0.0});
            if (!(l.hi > -1500))
                return {0.0, 0.0, 0};
            return exp_scaled(l.hi, l.lo);
        }

        // (b/a)^p à l'échelle.
        constexpr scaled power_ratio(double b, double a, double p) noexcept
        {
            auto const l = mul(sum(log_precise(b), negated(log_precise(a))), {p, 0.0});
            if (!(l.hi > -1500))
                return {0.0, 0.0, 0};
            return exp_scaled(l.hi, l.lo);
        }

        // Σ_{j=a}^{b} j^p par Euler–Maclaurin, a assez grand (a >= 8|p|, a >= 64), relativement à
        // a^p : la valeur est a^p · r, r rendu en double-double (r >= 1/2).
        constexpr dd euler_maclaurin_ratio(double a, double b, double p) noexcept
        {
            dd const L = sum(log_precise(b), negated(log_precise(a)));     // log(b/a)
            dd const q = sum({p, 0.0}, {1.0, 0.0});                          // p + 1
            dd const t = mul(q, L);
            // ∫_a^b x^p dx / a^p = a · (e^t − 1)/(p + 1) (= a · L si p = −1)
            dd integral{};
            if (magnitude(t.hi) < 0x1p-12)
            {
                // (e^t − 1)/t = 1 + t/2 + t²/6 + t³/24 + t⁴/120 + t⁵/720
                dd e{1.0 / 720, 0.0};
                for (double c : {120.0, 24.0, 6.0, 2.0, 1.0})
                    e = sum(mul(e, t), sd::quotient_dd({1.0, 0.0}, {c, 0.0}));
                integral = mul(mul({a, 0.0}, L), e);
            }
            else
            {
                dd em1{};
                if (t.hi > 700)
                    em1 = unscaled(exp_scaled(t.hi, t.lo));                 // −1 négligeable
                else
                {
                    auto const e = t.hi < -745 ? dd{0.0, 0.0} : unscaled(exp_scaled(t.hi, t.lo));
                    em1 = sum(e, {-1.0, 0.0});
                }
                integral = mul({a, 0.0}, sd::quotient_dd(em1, q));
            }
            // (f(a) + f(b))/2 / a^p = (1 + (b/a)^p)/2
            auto const ratio_b = unscaled(power_ratio(b, a, p));
            dd r = sum(integral, {0.5, 0.0});
            r = sum(r, {sd::fmul(ratio_b.hi, 0.5), sd::fmul(ratio_b.lo, 0.5)});
            // Σ_m B_{2m}/(2m)! (f^{(2m−1)}(b) − f^{(2m−1)}(a)) / a^p,
            // f^{(n)}(x) = p(p − 1)…(p − n + 1) x^{p−n} : (p)_n a^{−n} ((b/a)^{p−n} − 1)
            dd falling{p, 0.0};                                             // (p)_1
            dd const inv_a = sd::quotient_dd({1.0, 0.0}, {a, 0.0});
            dd scale = inv_a;                                               // a^{−n}, n = 1
            for (int m = 1; m <= 12; ++m)
            {
                int const n = 2 * m - 1;
                auto const rb = unscaled(power_ratio(b, a, p - n));         // (b/a)^{p−n}
                dd const term = mul(mul(bernoulli[m - 1], mul(falling, scale)), sum(rb, {-1.0, 0.0}));
                r = sum(r, term);
                if (magnitude(term.hi) < sd::fmul(magnitude(r.hi), 0x1p-118))
                    break;
                // (p)_{n+2} = (p)_n (p − n)(p − n − 1), a^{−n−2}
                falling = mul(mul(falling, {p - n, 0.0}), {p - n - 1, 0.0});
                scale = mul(mul(scale, inv_a), inv_a);
            }
            return r;
        }

        // Σ_{j=a}^{b} j^p, à l'échelle (a, b entiers exacts, 1 <= a).
        constexpr scaled power_sum(double a, double b, double p) noexcept
        {
            if (b < a)
                return {0.0, 0.0, 0};
            double const a0 = magnitude(p) * 8 > 64 ? double(std::int64_t(magnitude(p) * 8) + 1) : 64.0;   // ((|p| + 2m)/(2πa))² < 4·10⁻⁴
            scaled s{0.0, 0.0, 0};
            double j = a;
            for (; j <= b && j < a0; ++j)
                s = add_scaled(s, power(j, p));
            if (j > b)
                return s;
            auto const base = power(j, p);
            if (base.hi == 0)
                return s;
            auto const r = euler_maclaurin_ratio(j, b, p);
            return add_scaled(s, base * normalized(r.hi, r.lo));
        }
}

namespace wsl::inline WSL_DEFAULTS::probability::laws
{

    template <law_integer I = std::int64_t>
    class zipf : public law_interface
    {
    public:
        using value_type = I;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            I size;
            double exponent;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.size = n_, .exponent = s_}; }

        constexpr zipf(parameters p) : n_(p.size), s_(p.exponent)
        {
            if (!(n_ >= 1) || double(n_) > 0x1p53 || !(s_ >= 0) || !detail::is_finite(s_))
                throw std::domain_error("wsl::probability::laws::zipf : il faut 1 <= size <= 2^53 et exponent >= 0, fini");
            h_ = detail::zipf_sums::power_sum(1.0, double(n_), -s_);
            if (s_ > 0)
                setup_sampler();
        }

        constexpr I size() const noexcept { return n_; }
        constexpr double exponent() const noexcept { return s_; }
        constexpr interval<I> support() const noexcept { return {I(1), n_}; }
        constexpr I mode() const noexcept { return 1; }
        constexpr I median() const noexcept { return quantile(0.5); }

        // 1 + A/H, A = Σ (k − 1) k^{−s}
        constexpr double mean() const noexcept { return detail::rounded(detail::add_scaled(detail::normalized(1.0, 0.0), moment_a() / h_)); }
        // (B H − A²)/H²
        constexpr double variance() const noexcept
        {
            auto const a = moment_a(), b = moment_b();
            auto const ar = detail::unscaled(a / h_), br = detail::unscaled(b / h_);
            auto const v = detail::sum(br, detail::negated(detail::mul(ar, ar)));
            return v.hi > 0 ? v.hi + v.lo : 0.0;
        }

        constexpr double pmf(I k) const noexcept
        {
            if (k < 1 || k > n_)
                return 0;
            return detail::rounded(detail::zipf_sums::power(double(k), -s_) / h_);
        }
        constexpr double log_pmf(I k) const noexcept
        {
            if (k < 1 || k > n_)
                return -std::numeric_limits<double>::infinity();
            auto const l = detail::sum(detail::mul(detail::log_precise(double(k)), {-s_, 0.0}), detail::negated(detail::log(h_)));
            return l.hi + l.lo;
        }
        constexpr double cdf(I k) const noexcept
        {
            auto const [v, is_f] = distribution(k);
            return is_f ? detail::rounded(v) : detail::sd::one_minus(v);
        }
        constexpr double ccdf(I k) const noexcept
        {
            auto const [v, is_f] = distribution(k);
            return is_f ? detail::sd::one_minus(v) : detail::rounded(v);
        }

        constexpr I quantile(double u) const
        {
            auto const& X = *this;
            if (!(u >= 0 && u <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            if (u == 0)
                return 1;
            if (u == 1)
                return X.n_;
            auto const ge = [&](I k) {
                auto const [v, is_f] = X.distribution(k);
                return detail::at_least(v, is_f, u);
            };
            return detail::smallest_true(I(1), I(1), X.n_, ge);
        }

        template <std::uniform_random_bit_generator G>
        constexpr I sample(G& g) const
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            if (X.s_ == 0)
                return random::uniform_int(I(1), X.n_, g);
            for (;;)
            {
                double const u = X.h_integral_n_ + sd::fmul(random::uniform_real<double>(g, wsl::exact), X.h_integral_x1_ - X.h_integral_n_);
                double const x = X.h_integral_inverse(u);
                double k = x + 0.5;
                k = k < 1 ? 1 : k > double(X.n_) ? double(X.n_) : double(std::int64_t(k));
                if (k - x <= X.threshold_ || u >= X.h_integral(k + 0.5) - X.h(k))
                    return static_cast<I>(k);
            }
        }

        friend constexpr bool operator==(zipf const& a, zipf const& b) noexcept { return a.n_ == b.n_ && a.s_ == b.s_; }

    private:
        // F(k) (is_f) si la queue basse est la plus petite, sinon 1 − F(k), sommée directement.
        struct part
        {
            detail::scaled value;
            bool is_f;
        };
        constexpr part distribution(I k) const noexcept
        {
            if (k < 1)
                return {{0.0, 0.0, 0}, true};
            if (k >= n_)
                return {{0.0, 0.0, 0}, false};
            auto const lower = detail::zipf_sums::power_sum(1.0, double(k), -s_);
            auto const r = lower / h_;
            if (detail::compare(r, {0.5, 0.0}) <= 0)
                return {r, true};
            return {detail::zipf_sums::power_sum(double(k) + 1, double(n_), -s_) / h_, false};
        }

        // A = Σ_{k=2}^{N} (k − 1) k^{−s}, B = Σ_{k=2}^{N} (k − 1)² k^{−s}
        constexpr detail::scaled moment_a() const noexcept { return moment(1); }
        constexpr detail::scaled moment_b() const noexcept { return moment(2); }
        constexpr detail::scaled moment(int order) const noexcept
        {
            namespace z = detail::zipf_sums;
            double const n = double(n_);
            if (n < 2)
                return {0.0, 0.0, 0};
            constexpr double a0 = 64;                                       // queue : power_sum choisit son seuil
            detail::scaled direct{0.0, 0.0, 0};
            double j = 2;
            for (; j <= n && j < a0; ++j)                                   // termes positifs
            {
                double const c = order == 1 ? j - 1 : detail::sd::fmul(j - 1, j - 1);
                direct = detail::add_scaled(direct, z::power(j, -s_) * detail::normalized(c, 0.0));
            }
            if (j > n)
                return direct;
            // (k − 1) k^{−s} = k^{1−s} − k^{−s} ; (k − 1)² k^{−s} = k^{2−s} − 2 k^{1−s} + k^{−s}
            auto const p0 = z::power_sum(j, n, -s_), p1 = z::power_sum(j, n, 1 - s_);
            detail::scaled tail{};
            if (order == 1)
                tail = detail::scaled_difference(p1, p0);
            else
            {
                auto const p2 = z::power_sum(j, n, 2 - s_);
                tail = detail::scaled_difference(detail::add_scaled(p2, p0), p1 * detail::normalized(2.0, 0.0));
            }
            return detail::add_scaled(direct, tail);
        }

        // Rejet par inversion (Hörmann et Derflinger) : H(x) = ∫ x^{−s}, h(x) = x^{−s}.
        static constexpr double helper1(double x) noexcept                 // log1p(x)/x
        {
            namespace sd = detail::sd;
            return detail::magnitude(x) > 1e-8 ? wsl::random::detail::pm::log1p(x) / x : 1 - sd::fmul(x, 0.5 - sd::fmul(x, 1.0 / 3 - sd::fmul(0.25, x)));
        }
        static constexpr double helper2(double x) noexcept                 // expm1(x)/x
        {
            namespace sd = detail::sd;
            return detail::magnitude(x) > 1e-8 ? wsl::random::detail::pm::expm1(x) / x : 1 + sd::fmul(sd::fmul(x, 0.5), 1 + sd::fmul(x / 3, 1 + sd::fmul(0.25, x)));
        }
        constexpr double h_integral(double x) const noexcept
        {
            double const l = wsl::random::detail::pm::log(x);
            return detail::sd::fmul(helper2(detail::sd::fmul(1 - s_, l)), l);
        }
        constexpr double h(double x) const noexcept { return wsl::random::detail::pm::exp(detail::sd::fmul(-s_, wsl::random::detail::pm::log(x))); }
        constexpr double h_integral_inverse(double x) const noexcept
        {
            double t = detail::sd::fmul(x, 1 - s_);
            if (t < -1)
                t = -1;
            return wsl::random::detail::pm::exp(detail::sd::fmul(helper1(t), x));
        }
        constexpr void setup_sampler() noexcept
        {
            h_integral_x1_ = h_integral(1.5) - 1;
            h_integral_n_ = h_integral(double(n_) + 0.5);
            threshold_ = 2 - h_integral_inverse(h_integral(2.5) - h(2));
        }

        I n_;
        double s_;
        detail::scaled h_{1.0, 0.0, 0};
        double h_integral_x1_ = 0, h_integral_n_ = 0, threshold_ = 0;
    };
}
