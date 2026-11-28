#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Loi de Rice : module R = |(ν + σ Z₁, σ Z₂)| d'une gaussienne bidimensionnelle décentrée,
// ν >= 0, σ > 0 (ν = 0 : Rayleigh). En unités réduites y = x/σ, c = ν/σ :
// densité g(y) = y e^{−(y² + c²)/2} I₀(c y) ; répartition 1 − Q₁(c, y) (Marcum).
//
//   wsl::probability::rice R({.scale = 1.5, .noncentrality = 2.0});
//   R.pdf(2.0); R.cdf(1.0); R.quantile(0.9); R.mean(); R(g);
//
// Densité en logarithme double-double (log I₀ de wsl::special, corrigé au premier ordre de
// l'arrondi de z = c y par I₁/I₀). Répartition par les séries de Marcum, termes positifs :
// S = e^{−(y² + c²)/2} Σ_{k>=0} (c/y)^k I_k(c y), F = e^{−(y² + c²)/2} Σ_{k>=1} (y/c)^k I_k(c y),
// la plus petite des deux (y face à √(c² + 2 log 2)), l'autre par complément ; quantiles
// par le solveur commun (detail/root.hpp). Moyenne σ √(π/2) L_{1/2}(−c²/2), variance
// σ² (2 + c² − m²) ; pour c >= 12, développement asymptotique m = c (1 + Σ c_k (2/c²)^k),
// c_k = ((−1/2)_k)²/k!, d'où δ = m − c et Var/σ² = 2 − 2cδ − δ² sans compensation.
// ν = 0 : fonctions de rayleigh (formes closes). Tirage (approché, bit-exact) :
// hypot(ν + σ Z₁, σ Z₂), Z exactes.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/root.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/laws/normal.hpp>
#include <wsl/probability/laws/rayleigh.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_real T = double>
    class rice : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T scale = 1;
            T noncentrality;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.scale = sigma_, .noncentrality = nu_}; }

        constexpr rice(parameters p) : sigma_(p.scale), nu_(p.noncentrality)
        {
            if (!(sigma_ > 0) || !detail::is_finite(sigma_) || !(nu_ >= 0) || !detail::is_finite(nu_))
                throw std::domain_error("wsl::probability::rice : il faut scale > 0 et noncentrality >= 0, finis");
            c_ = nu_ == 0 ? sd_t{0.0, 0.0} : detail::unscaled(detail::normalized(double(nu_), 0.0) / double(sigma_));
            log_sigma_ = detail::log_precise(double(sigma_));
            if (nu_ == 0)
                rayleigh_ = rayleigh<T>({.scale = sigma_});
        }

        constexpr T scale() const noexcept { return sigma_; }
        constexpr T noncentrality() const noexcept { return nu_; }
        constexpr interval<T> support() const noexcept { return {T(0), std::numeric_limits<T>::infinity()}; }

        constexpr T mean() const noexcept
        {
            if (nu_ == 0)
                return rayleigh_.mean();
            auto const m = reduced_mean();
            return static_cast<T>(detail::rounded(detail::normalized(m.hi, m.lo) * detail::normalized(double(sigma_), 0.0)));
        }
        constexpr T variance() const noexcept
        {
            if (nu_ == 0)
                return rayleigh_.variance();
            auto const v = reduced_variance();
            return static_cast<T>(detail::rounded(detail::normalized(v.hi, v.lo) * detail::normalized(double(sigma_), 0.0) * detail::normalized(double(sigma_), 0.0)));
        }
        constexpr T median() const noexcept { return quantile(0.5); }

        constexpr T pdf(T x) const noexcept
        {
            if (nu_ == 0)
                return rayleigh_.pdf(x);
            double const v = x;
            if (!(v > 0))
                return v == v ? T(0) : x;
            if (!detail::is_finite(v))
                return 0;
            auto const l = detail::sum(log_reduced_density(reduced(v)), detail::negated(log_sigma_));
            return l.hi > -1500 ? static_cast<T>(detail::rounded(detail::exp_scaled(l.hi, l.lo))) : T(0);
        }
        constexpr T log_pdf(T x) const noexcept
        {
            if (nu_ == 0)
                return rayleigh_.log_pdf(x);
            double const v = x;
            if (!(v > 0))
                return v == v ? -std::numeric_limits<T>::infinity() : x;
            if (!detail::is_finite(v))
                return -std::numeric_limits<T>::infinity();
            auto const l = detail::sum(log_reduced_density(reduced(v)), detail::negated(log_sigma_));
            return static_cast<T>(l.hi + l.lo);
        }
        constexpr T cdf(T x) const noexcept { return nu_ == 0 ? rayleigh_.cdf(x) : tail(x, false); }
        constexpr T ccdf(T x) const noexcept { return nu_ == 0 ? rayleigh_.ccdf(x) : tail(x, true); }
        constexpr T quantile(double p) const { return nu_ == 0 ? rayleigh_.quantile(p) : inverse(p, false); }
        constexpr T cquantile(double q) const { return nu_ == 0 ? rayleigh_.cquantile(q) : inverse(q, true); }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            double const a = detail::native_sample(normal<double>({.mean = double(X.nu_), .stddev = double(X.sigma_)}), g);
            double const b = detail::native_sample(normal<double>({.mean = 0.0, .stddev = double(X.sigma_)}), g);
            return static_cast<T>(wsl::special::hypot(a, b));
        }
        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g, approximate_t) const
        {
            auto const& X = *this;
            double const a = normal<double>({.mean = double(X.nu_), .stddev = double(X.sigma_)})(g, wsl::approximate);
            double const b = normal<double>({.mean = 0.0, .stddev = double(X.sigma_)})(g, wsl::approximate);
            return static_cast<T>(wsl::special::hypot(a, b));
        }

        friend constexpr bool operator==(rice const& a, rice const& b) noexcept { return a.sigma_ == b.sigma_ && a.nu_ == b.nu_; }

    private:
        using sd_t = detail::sd::double_double;

        // y = x/σ en double-double (x > 0 fini).
        constexpr sd_t reduced(double x) const noexcept { return detail::unscaled(detail::normalized(x, 0.0) / double(sigma_)); }

        // log I_n(z), z en double-double : log I_n(z_hi) + z_lo (I_{n+1}/I_n + n/z_hi) ; mémoire des
        // trois derniers (les séries avancent d'un pas).
        struct bessel_cache
        {
            double n[3] = {-1, -1, -1};
            sd_t v[3]{};
            int next = 0;
        };
        static constexpr sd_t log_bessel_hi(double n, double z, bessel_cache& c) noexcept
        {
            for (int i = 0; i < 3; ++i)
                if (c.n[i] == n)
                    return c.v[i];
            auto const v = wsl::special::detail::log_cyl_bessel_i(n, z);
            c.n[c.next] = n;
            c.v[c.next] = v;
            c.next = (c.next + 1) % 3;
            return v;
        }
        static constexpr sd_t log_bessel(double n, sd_t z, bessel_cache& c) noexcept
        {
            namespace sd = detail::sd;
            auto const l = log_bessel_hi(n, z.hi, c);
            if (z.lo == 0)
                return l;
            auto const up = log_bessel_hi(n + 1, z.hi, c);
            double const ratio = wsl::special::exp((up.hi - l.hi) + (up.lo - l.lo));   // I_{n+1}/I_n
            return detail::sum(l, {sd::fmul(z.lo, ratio + n / z.hi), 0.0});
        }

        // −(y² + c²)/2, en double-double.
        constexpr sd_t minus_half_square(sd_t y) const noexcept
        {
            namespace sd = detail::sd;
            auto const s = detail::sum(detail::mul(y, y), detail::mul(c_, c_));
            return {sd::fmul(-s.hi, 0.5), sd::fmul(-s.lo, 0.5)};
        }

        // log g(y) = log y − (y² + c²)/2 + log I₀(c y).
        constexpr sd_t log_reduced_density(sd_t y) const noexcept
        {
            auto l = detail::sum(detail::log(detail::normalized(y.hi, y.lo)), minus_half_square(y));
            if (c_.hi != 0)
            {
                bessel_cache cache;
                l = detail::sum(l, log_bessel(0, detail::mul(c_, y), cache));
            }
            return l;
        }

        // F(y) (upper = false) ou S(y), et la plus petite des deux sommée directement (is_lower :
        // la valeur rendue est F).
        struct part
        {
            detail::scaled value;
            bool is_lower;
        };
        constexpr part parts(sd_t y) const noexcept                         // c > 0
        {
            namespace sd = detail::sd;
            auto const e0 = minus_half_square(y);
            double const median_guess = wsl::special::sqrt(sd::fmul(c_.hi, c_.hi) + 1.3862943611198906);
            bool const lower = y.hi < median_guess;
            auto const z = detail::mul(c_, y);
            // log r : r = c/y (S) ou y/c (F)
            auto const lr = detail::sum(detail::log(detail::normalized(c_.hi, c_.lo)), detail::negated(detail::log(detail::normalized(y.hi, y.lo))));
            sd_t const log_ratio = lower ? detail::negated(lr) : lr;
            bessel_cache cache;
            double k = lower ? 1 : 0;
            auto const term_log = [&](double n) { return detail::sum(detail::sum(e0, detail::mul({n, 0.0}, log_ratio)), log_bessel(n, z, cache)); };
            auto const base = term_log(k);
            sd_t total{0.0, 0.0};
            double previous = std::numeric_limits<double>::infinity();
            for (int step = 0; step < 1000000; ++step, ++k)
            {
                auto const e = detail::sum(term_log(k), detail::negated(base));
                if (!(e.hi > -1500))
                {
                    if (k > z.hi)
                        break;
                    continue;
                }
                auto const w = detail::unscaled(detail::exp_scaled(e.hi, e.lo));
                total = detail::sum(total, w);
                if (step > 0 && w.hi < previous && w.hi < sd::fmul(total.hi, 0x1p-112))
                    break;
                previous = w.hi;
            }
            if (!(base.hi > -1500) || total.hi == 0)
                return {{0.0, 0.0, 0}, lower};
            return {detail::exp_scaled(base.hi, base.lo) * detail::normalized(total.hi, total.lo), lower};
        }

        constexpr T tail(T x, bool upper) const noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0))
                return upper ? T(1) : T(0);
            if (!detail::is_finite(v))
                return upper ? T(0) : T(1);
            auto const [value, is_lower] = parts(reduced(v));
            return static_cast<T>(is_lower == !upper ? detail::rounded(value) : detail::sd::one_minus(value));
        }

        // y tel que F(y) = t (upper = false) ou S(y) = t, t <= ½.
        constexpr sd_t solve(double t, bool upper) const
        {
            namespace sd = detail::sd;
            auto const evaluate = [&](sd_t y) {
                detail::part_value r{{0.0, 0.0, 0}, {0.0, 0.0, 0}, y};
                if (!(y.hi > 0))
                {
                    r.value = upper ? detail::normalized(1.0, 0.0) : detail::scaled{0.0, 0.0, 0};
                    return r;
                }
                auto const [value, is_lower] = parts(y);
                bool const want_lower = !upper;
                if (is_lower == want_lower)
                    r.value = value;
                else
                {
                    auto const one = detail::unscaled(value);
                    auto const c = detail::sum({1.0, 0.0}, detail::negated(one));
                    r.value = c.hi > 0 ? detail::normalized(c.hi, c.lo) : detail::scaled{0.0, 0.0, 0};
                }
                auto const l = log_reduced_density(y);
                if (l.hi > -1500)
                    r.slope = detail::exp_scaled(l.hi, l.lo);
                return r;
            };
            // départ : la médiane approchée ; petite probabilité de répartition : F ≈ e^{−c²/2} y²/2
            // (c y petit), y₀ = √(2t) e^{c²/4}, calculé en logarithme (t peut être sous-normal)
            double start = wsl::special::sqrt(sd::fmul(c_.hi, c_.hi) + 1.3862943611198906);
            if (!upper)
            {
                double const log_y0 = sd::fmul(0.5, wsl::special::log(t) + 0.6931471805599453) + sd::fmul(sd::fmul(c_.hi, c_.hi), 0.25);
                if (log_y0 < wsl::special::log(start))
                    start = wsl::special::exp(log_y0);
            }
            return detail::solve_part(evaluate, !upper, t, {.start = start, .lo = 0.0, .max_step = 4},
                                      [](double y, int side) { return side > 0 ? y * 0.5 : detail::sd::fmul(y, 2.0) + 1; }, detail::arithmetic_midpoint, 400);
        }

        constexpr T inverse(double t, bool upper) const
        {
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::probability::quantile : probabilité hors de [0, 1]");
            if (t == 0)
                return upper ? std::numeric_limits<T>::infinity() : T(0);
            if (t == 1)
                return upper ? T(0) : std::numeric_limits<T>::infinity();
            auto const y = t <= 0.5 ? solve(t, upper) : solve(1.0 - t, !upper);
            if (!(y.hi > 0))
                return 0;
            return static_cast<T>(detail::rounded(detail::normalized(y.hi, y.lo) * detail::normalized(double(sigma_), 0.0)));
        }

        // Moyenne réduite m(c) = √(π/2) L_{1/2}(−c²/2), en double-double.
        constexpr sd_t reduced_mean() const noexcept
        {
            namespace sd = detail::sd;
            constexpr sd_t root_half_pi{0x1.40d931ff62706p+0, -0x1.a6a0d6f814637p-54};   // √(π/2) (mpmath)
            if (c_.hi == 0)
                return root_half_pi;
            if (c_.hi >= 12)
                return detail::sum(c_, asymptotic_delta());
            // e^{−w}[(1 + 2w) I₀(w) + 2w I₁(w)], w = c²/4
            auto const c2 = detail::mul(c_, c_);
            sd_t const w{sd::fmul(c2.hi, 0.25), sd::fmul(c2.lo, 0.25)};
            bessel_cache cache;
            auto const l0 = detail::sum(log_bessel(0, w, cache), detail::negated(w));
            auto const l1 = detail::sum(log_bessel(1, w, cache), detail::negated(w));
            auto const i0 = detail::unscaled(detail::exp_scaled(l0.hi, l0.lo)), i1 = detail::unscaled(detail::exp_scaled(l1.hi, l1.lo));
            sd_t const two_w{sd::fmul(w.hi, 2.0), sd::fmul(w.lo, 2.0)};
            auto const inner = detail::sum(detail::mul(detail::sum({1.0, 0.0}, two_w), i0), detail::mul(two_w, i1));
            return detail::mul(root_half_pi, inner);
        }
        // δ = m − c = c Σ_{k>=1} c_k (2/c²)^k, c_{k+1} = c_k (k − ½)²/(k + 1) ; c >= 12.
        constexpr sd_t asymptotic_delta() const noexcept
        {
            auto const c2 = detail::mul(c_, c_);
            sd_t const x = detail::sd::quotient_dd({2.0, 0.0}, c2);         // 2/c²
            sd_t coefficient{0.25, 0.0};                                    // c₁
            sd_t power = x;
            sd_t s{0.0, 0.0};
            for (int k = 1; k < 400; ++k)
            {
                auto const term = detail::mul(coefficient, power);
                s = detail::sum(s, term);
                if (term.hi < detail::sd::fmul(s.hi, 0x1p-112))
                    break;
                double const h = double(k) - 0.5;
                coefficient = detail::sd::quotient_dd(detail::mul(coefficient, {detail::sd::fmul(h, h), 0.0}), {double(k + 1), 0.0});
                power = detail::mul(power, x);
                if (detail::mul(coefficient, power).hi > term.hi)           // série asymptotique : arrêt optimal
                    break;
            }
            return detail::mul(c_, s);
        }
        // Variance réduite : 2 + c² − m², ou 2 − 2cδ − δ² (c >= 12).
        constexpr sd_t reduced_variance() const noexcept
        {
            namespace sd = detail::sd;
            if (c_.hi >= 12)
            {
                auto const d = asymptotic_delta();
                auto const two_cd = detail::mul(c_, d);
                return detail::sum(detail::sum({2.0, 0.0}, sd_t{sd::fmul(-two_cd.hi, 2.0), sd::fmul(-two_cd.lo, 2.0)}), detail::negated(detail::mul(d, d)));
            }
            auto const m = reduced_mean();
            return detail::sum(detail::sum({2.0, 0.0}, detail::mul(c_, c_)), detail::negated(detail::mul(m, m)));
        }

        T sigma_, nu_;
        rayleigh<T> rayleigh_{};
        sd_t c_{0.0, 0.0}, log_sigma_{0.0, 0.0};
    };
}
