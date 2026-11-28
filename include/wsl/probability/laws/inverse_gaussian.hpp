#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Loi inverse-gaussienne (de Wald) de moyenne μ et de forme λ :
// densité √(λ/(2πx³)) e^{−λ(x − μ)²/(2μ²x)}, x > 0.
//
//   wsl::probability::laws::inverse_gaussian W({.mean = 1, .shape = 3});
//   W.cdf(2.0); W.ccdf(50.0); W.quantile(0.99); W.variance(); W(g);
//
// Avec c = √(λ/x), a = c (x − μ)/μ, δ = 2c et R(s) = e^{s²/2} Φ(−s) (rapport de Mills) :
//     P(X > x) = e^{−a²/2} (R(a) − R(a + δ)),   P(X <= x) = e^{−a²/2} (R(|a|) + R(a + δ)) si a < 0.
// La survie se compense quand δ est petit devant 1 + a⁺ (queue droite, x >> μ, ou λ/x
// petit) : la différence est alors la série −Σ_{k>=1} r_k(a) δᵏ des coefficients de Taylor
// de R en a, (k + 1) r_{k+1} = a r_k + r_{k−1}, r_1 = a r_0 − 1/√(2π) — récurrence avancée
// pour a < 4,5, rétrograde (Miller, 90 termes) au-delà. Sinon, différence directe (perte
// d'un facteur 6 au plus). R lui-même : séries en 0 (½ e^{s²/2} − Σ s^{2n+1}/((2n+1)!! √(2π)))
// sous 4,5, Miller jusqu'à 2048, développement asymptotique au-delà ; à 2⁻⁸⁵ près.
// Répartition : la part la plus petite est calculée, l'autre en est le complément.
// Quantiles : Newton en log x sur log F (ou log S), encadré, puis un pas de Newton à
// résidu précis depuis le flottant trouvé ; un arrondi.
//
// Tirage (approché, bit-exact) : Michael, Schucany et Haas (1976), une normale exacte et
// une épreuve de Bernoulli exacte.
//
// Pas d'entropie : elle demande E₁(2λ/μ), absente de wsl::special.

#include <wsl/probability/laws/basic.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/discrete.hpp>
#include <wsl/probability/detail/root.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/laws/normal.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/random/detail/portable_math.hpp>

#include <bit>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
}

// Outils internes : dans wsl::probability::detail (un seul « detail » visible depuis probability).
namespace wsl::inline WSL_DEFAULTS::probability::detail
{
        inline constexpr int mills_order = 90;                              // départ de Miller

        // R(s) = e^{s²/2} Φ(−s) pour |s| < 4,5 : ½ Σ (s²/2)ⁿ/n! − (1/√(2π)) Σ s^{2n+1}/(2n+1)!!,
        // termes positifs en double-double (perte de 2²² au plus, à 4,5).
        constexpr sd::double_double mills_near_zero(sd::double_double s) noexcept
        {
            auto const s2 = mul(s, s);
            auto const h = times_double(s2, 0.5);
            sd::double_double e{1.0, 0.0}, et{1.0, 0.0}, o = s, ot = s;
            for (int n = 1; n < 200; ++n)
            {
                et = sd::quotient_dd(mul(et, h), {static_cast<double>(n), 0.0});
                ot = sd::quotient_dd(mul(ot, s2), {static_cast<double>(2 * n + 1), 0.0});
                e = sum(e, et);
                o = sum(o, ot);
                if (magnitude(et.hi) < 0x1p-114 && magnitude(ot.hi) < 0x1p-114)
                    break;
            }
            auto const odd = mul(o, {inv_sqrt_2pi_hi, inv_sqrt_2pi_lo});
            return sum(times_double(e, 0.5), negated(odd));
        }

        // Coefficients de Taylor r_0 … r_{n−1} de R en a, 4,5 <= a < 2048 : récurrence
        // rétrograde depuis r_{N+1} = 0, r_N = 1 (solution minimale), renormalisée par 2^-300
        // quand elle grandit, puis normalisée par a r_0 − r_1 = 1/√(2π) (termes positifs).
        constexpr void mills_miller(sd::double_double a, sd::double_double* r, int n) noexcept
        {
            sd::double_double y[mills_order + 2]{};
            y[mills_order] = {1.0, 0.0};
            for (int k = mills_order; k >= 1; --k)
            {
                y[k - 1] = sum(times_double(y[k + 1], static_cast<double>(k + 1)), negated(mul(a, y[k])));
                if (magnitude(y[k - 1].hi) > 0x1p300)                       // reste sous 2^995 (produits de Dekker)
                    for (int j = k - 1; j <= mills_order; ++j)
                        y[j] = {sd::scale_by_power_of_two(y[j].hi, -300), sd::scale_by_power_of_two(y[j].lo, -300)};
            }
            auto const norm = sd::quotient_dd({inv_sqrt_2pi_hi, inv_sqrt_2pi_lo}, sum(mul(a, y[0]), negated(y[1])));
            for (int k = 0; k < n; ++k)
                r[k] = mul(y[k], norm);
        }

        // R(s), s >= 2048 à l'échelle : (1/(√(2π) s)) (1 − v + 3v² − 15v³ + 105v⁴ − 945v⁵), v = 1/s².
        constexpr scaled mills_asymptotic(scaled const& s) noexcept
        {
            double v = 0;
            if (s.e < 500)
            {
                auto const [s_hi, s_lo] = unscaled(s);
                v = 1.0 / sd::fmul(s_hi, s_hi);
            }
            double series = 0;
            for (int k = 5; k >= 1; --k)                                    // −v (1 − 3v (1 − 5v …))
                series = sd::fmul(-static_cast<double>(2 * k - 1), sd::fmul(v, 1.0 + series));
            auto const [t_hi, t_lo] = sd::fast_two_sum(1.0, series);
            return normalized(inv_sqrt_2pi_hi, inv_sqrt_2pi_lo) * normalized(t_hi, t_lo) / s;
        }

        // R(s), s >= −1/8, à l'échelle.
        constexpr scaled mills(sd::double_double s) noexcept
        {
            if (s.hi < 4.5)
            {
                auto const r = mills_near_zero(s);
                return normalized(r.hi, r.lo);
            }
            if (s.hi < 2048)
            {
                sd::double_double r[1]{};
                mills_miller(s, r, 1);
                return normalized(r[0].hi, r[0].lo);
            }
            return mills_asymptotic(normalized(s.hi, s.lo));
        }

        // R(a) − R(a + δ) = −Σ_{k>=1} r_k δᵏ, −1/8 < a < 64, 0 < δ < (1 + a⁺)/4 à l'échelle :
        // termes décroissants d'un facteur 4 au moins, sans compensation.
        constexpr scaled mills_difference(sd::double_double a, scaled const& delta) noexcept
        {
            constexpr int count = mills_order - 5;
            sd::double_double r[count]{};
            if (a.hi < 4.5)
            {
                r[0] = mills_near_zero(a);
                r[1] = sum(mul(a, r[0]), {-inv_sqrt_2pi_hi, -inv_sqrt_2pi_lo});
                for (int k = 1; k + 1 < count; ++k)
                    r[k + 1] = sd::quotient_dd(sum(mul(a, r[k]), r[k - 1]), {static_cast<double>(k + 1), 0.0});
            }
            else
                mills_miller(a, r, count);
            // −Σ r_k δ^{k−1}, puis le facteur δ à l'échelle (δ minuscule : r_1 seul)
            sd::double_double t = negated(r[1]);
            if (delta.e > -500)
            {
                auto const d = unscaled(delta);
                sd::double_double power{1.0, 0.0};
                for (int k = 2; k < count; ++k)
                {
                    power = mul(power, d);
                    auto const term = negated(mul(r[k], power));
                    t = sum(t, term);
                    if (magnitude(term.hi) < sd::fmul(magnitude(t.hi), 0x1p-82))
                        break;
                }
            }
            return normalized(t.hi, t.lo) * delta;
        }

        // Fonctions de l'inverse-gaussienne en x > 0 fini : répartition, survie (à
        // l'échelle, 0 : hi = 0) et éléments de la densité.
        struct inverse_gaussian_parts
        {
            scaled lower, upper;
            scaled c_over_x;                                                // √(λ/x³)
            sd::double_double half_a2;                                      // a²/2 (+∞ si hors de portée)
        };

        constexpr inverse_gaussian_parts inverse_gaussian_at(double x, double mu, double lambda) noexcept
        {
            constexpr double infinity = std::numeric_limits<double>::infinity();
            inverse_gaussian_parts out{};
            scaled const c = sqrt(normalized(lambda, 0.0) / x);
            out.c_over_x = c / x;
            auto const [d_hi, d_lo] = sd::two_sum(x, -mu);
            bool const negative = d_hi < 0;
            scaled a_mag{0.0, 0.0, 0};
            if (d_hi != 0)
                a_mag = c * (normalized(negative ? -d_hi : d_hi, negative ? -d_lo : d_lo) / mu);
            if (a_mag.hi != 0 && a_mag.e >= 6)                              // |a| >= 64 : e^{−a²/2} < 2^-2900
            {
                scaled const a2 = a_mag * a_mag;
                out.half_a2 = a2.e > 1000 ? sd::double_double{infinity, 0.0} : times_double(unscaled(a2), 0.5);
                if (negative)
                    out.upper = {1.0, 0.0, 0};
                else
                    out.lower = {1.0, 0.0, 0};
                return out;
            }
            sd::double_double a{0.0, 0.0};
            if (a_mag.hi != 0)
            {
                auto const u = unscaled(a_mag);
                a = negative ? negated(u) : u;
            }
            out.half_a2 = times_double(mul(a, a), 0.5);
            scaled const e = exp_scaled(-out.half_a2.hi, -out.half_a2.lo);
            scaled delta = c;
            delta.e += 1;
            // régime de la série : δ < (1 + a⁺)/4
            sd::double_double const bound = times_double(sum({1.0, 0.0}, negative ? sd::double_double{0.0, 0.0} : a), 0.25);
            bool const by_series = compare(delta, bound) < 0;
            auto const r_b = [&] {
                if (delta.e >= 1000)
                    return mills_asymptotic(delta);
                return mills(sum(a, unscaled(delta)));
            };
            if (!negative)
            {
                scaled const d = by_series ? mills_difference(a, delta) : scaled_difference(mills(a), r_b());
                out.upper = d.hi == 0 ? d : e * d;
                out.lower = scaled_complement(out.upper);
            }
            else
            {
                out.lower = e * add_scaled(mills(negated(a)), r_b());
                out.upper = by_series ? e * mills_difference(a, delta) : scaled_complement(out.lower);
            }
            return out;
        }
}

namespace wsl::inline WSL_DEFAULTS::probability::laws
{

    template <law_real T = double>
    class inverse_gaussian : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T mean, shape;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.mean = mean(), .shape = shape()}; }

        constexpr inverse_gaussian(parameters p) : mu_(p.mean), lambda_(p.shape)
        {
            if (!(mu_ > 0) || !detail::is_finite(mu_) || !(lambda_ > 0) || !detail::is_finite(lambda_))
                throw std::domain_error("wsl::probability::laws::inverse_gaussian : il faut mean > 0 et shape > 0, finis");
        }

        constexpr T mean() const noexcept { return mu_; }
        constexpr T shape() const noexcept { return lambda_; }
        constexpr interval<T> support() const noexcept { return {T(0), std::numeric_limits<T>::infinity()}; }
        // μ³/λ
        constexpr T variance() const noexcept
        {
            detail::scaled const m = detail::normalized(double(mu_), 0.0);
            return static_cast<T>(detail::rounded(m * m * m / double(lambda_)));
        }
        // 3 √(μ/λ)
        constexpr T skewness() const noexcept
        {
            return static_cast<T>(detail::rounded(detail::sqrt(detail::normalized(double(mu_), 0.0) / double(lambda_)) * detail::normalized(3.0, 0.0)));
        }
        // 15 μ/λ
        constexpr T excess_kurtosis() const noexcept
        {
            return static_cast<T>(detail::rounded(detail::normalized(double(mu_), 0.0) / double(lambda_) * detail::normalized(15.0, 0.0)));
        }
        // μ (√(1 + t²) − t) = μ/(√(1 + t²) + t), t = 3μ/(2λ)
        constexpr T mode() const noexcept
        {
            detail::scaled t = detail::normalized(double(mu_), 0.0) / double(lambda_) * detail::normalized(3.0, 0.0);
            t.e -= 1;
            detail::scaled w = t;
            if (t.e > 400)
                w.e += 1;
            else
            {
                auto const td = detail::unscaled(t);
                auto const one_plus = detail::sum({1.0, 0.0}, detail::mul(td, td));
                w = detail::add_scaled(detail::sqrt(detail::normalized(one_plus.hi, one_plus.lo)), t);
            }
            return static_cast<T>(detail::rounded(detail::normalized(double(mu_), 0.0) / w));
        }
        constexpr T median() const noexcept { return wsl::probability::quantile(*this, 0.5); }

        constexpr T pdf(T x) const noexcept
        {
            auto const& W = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0) || v == std::numeric_limits<double>::infinity())
                return 0;
            auto const p = detail::inverse_gaussian_at(v, double(W.mu_), double(W.lambda_));
            if (!(p.half_a2.hi < 1500))
                return 0;
            detail::scaled const e = detail::exp_scaled(-p.half_a2.hi, -p.half_a2.lo);
            return static_cast<T>(detail::rounded(p.c_over_x * e * detail::normalized(detail::inv_sqrt_2pi_hi, detail::inv_sqrt_2pi_lo)));
        }
        // log √(λ/x³) − ½ log 2π − a²/2
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& W = *this;
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0) || v == std::numeric_limits<double>::infinity())
                return -std::numeric_limits<T>::infinity();
            auto const p = detail::inverse_gaussian_at(v, double(W.mu_), double(W.lambda_));
            auto l = detail::sum(detail::log(p.c_over_x), {-detail::half_log_2pi_hi, -detail::half_log_2pi_lo});
            l = detail::sum(l, detail::negated(p.half_a2));
            return static_cast<T>(l.hi + l.lo);
        }
        constexpr T cdf(T x) const noexcept { auto const& W = *this; return W.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& W = *this; return W.distribution(x, true); }
        constexpr T quantile(double p) const { auto const& W = *this; return W.inverse(p, false); }
        constexpr T cquantile(double q) const { auto const& W = *this; return W.inverse(q, true); }

        // Michael, Schucany et Haas : y = N², x₁ = μ²/x₂, x₂ = μ (1 + (μy + √(μy (4λ + μy)))/(2λ)) ;
        // x₁ avec probabilité μ/(μ + x₁), x₂ sinon.
        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& W = *this;
            namespace fp = detail::sd;
            double const mu = W.mu_, lambda = W.lambda_;
            double const z = detail::native_sample(normal<double>({.mean = 0.0, .stddev = 1.0}), g);
            double const w = fp::fmul(mu, fp::fmul(z, z));
            double const s = wsl::random::detail::pm::sqrt(fp::fmul(w, fp::fmul(4.0, lambda) + w));
            double const x2 = mu + fp::fmul(mu, (w + s) / fp::fmul(2.0, lambda));
            double const x1 = fp::fmul(mu, mu / x2);
            return static_cast<T>(random::bernoulli_trial(mu / (mu + x1), g) ? x1 : x2);
        }

        friend constexpr bool operator==(inverse_gaussian const&, inverse_gaussian const&) = default;

    private:
        constexpr T distribution(T x, bool upper) const noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (!(v > 0))
                return upper ? T(1) : T(0);
            if (v == std::numeric_limits<double>::infinity())
                return upper ? T(0) : T(1);
            auto const p = detail::inverse_gaussian_at(v, double(mu_), double(lambda_));
            return static_cast<T>(detail::rounded(upper ? p.upper : p.lower));
        }

        // Plus petit côté : F(x) = t (upper = false) ou S(x) = t, t <= ½ ; solveur commun
        // (detail/root.hpp) dans la variable u = log x (dV/du = x f), pas d'au plus 16.
        constexpr double solve(double t, bool upper) const
        {
            double const mu = mu_, lambda = lambda_;
            auto const evaluate = [&](detail::sd::double_double u) {
                double x = detail::rounded(detail::exp_scaled(u.hi, u.lo));
                if (!(x > 0))
                    x = std::numeric_limits<double>::denorm_min();
                if (!(x < std::numeric_limits<double>::infinity()))
                    x = std::numeric_limits<double>::max();
                auto const p = detail::inverse_gaussian_at(x, mu, lambda);
                detail::part_value r{upper ? p.upper : p.lower, {0.0, 0.0, 0}, detail::log_precise(x)};
                if (r.value.hi != 0 && p.half_a2.hi < 1500)
                    r.slope = p.c_over_x * detail::exp_scaled(-p.half_a2.hi, -p.half_a2.lo)
                            * detail::normalized(detail::inv_sqrt_2pi_hi, detail::inv_sqrt_2pi_lo) * detail::normalized(x, 0.0);
                return r;
            };
            auto const outward = [](double u, int side) {
                double const next = side > 0 ? u - 16 : u + 16;
                return next < -745.0 ? -745.0 : next > 709.7 ? 709.7 : next;
            };
            auto const u = detail::solve_part(evaluate, !upper, t, {.start = wsl::random::detail::pm::log(mu), .max_step = 16}, outward, detail::arithmetic_midpoint);
            return nearest_double(detail::rounded(detail::exp_scaled(u.hi, u.lo)), t, upper);
        }

        // Finition : V (F ou S) n'est évaluée qu'en des doubles ; quand elle varie de plus d'un
        // ulp par ulp de x (x f/V grand : queues des lois concentrées), le solveur peut s'arrêter
        // à quelques ulp. Encadrer la racine entre deux doubles voisins x_a <= x* < x_b, puis
        // garder le plus proche : V linéaire à l'échelle d'un ulp, x_a si t est plus près de V(x_a)
        // que du milieu. Au plus 8 pas de chaque côté, sinon le résultat du solveur.
        constexpr double nearest_double(double x, double t, bool upper) const noexcept
        {
            if (!(x > std::numeric_limits<double>::min() && x < std::numeric_limits<double>::max()))
                return x;
            double const mu = mu_, lambda = lambda_;
            auto const value = [&](double y) {
                auto const p = detail::inverse_gaussian_at(y, mu, lambda);
                return upper ? p.upper : p.lower;
            };
            auto const step = [](double y, int d) { return std::bit_cast<double>(std::bit_cast<std::uint64_t>(y) + static_cast<std::uint64_t>(std::int64_t{d})); };
            // y à gauche de la racine (ou dessus) : F(y) <= t, S(y) >= t
            auto const left = [&](detail::scaled const& v) { int const c = detail::compare(v, {t, 0.0}); return upper ? c >= 0 : c <= 0; };
            double a = x;
            detail::scaled va = value(a);
            int moves = 0;
            for (; moves < 8 && !left(va); ++moves)
                va = value(a = step(a, -1));
            detail::scaled vb = value(step(a, 1));
            for (; moves < 16 && left(vb); ++moves)
            {
                a = step(a, 1);
                va = vb;
                vb = value(step(a, 1));
            }
            if (moves >= 16 || !left(va))
                return x;
            detail::scaled middle = detail::add_scaled(va, vb);
            middle.e -= 1;
            int const c = detail::compare(middle, {t, 0.0});                // signe de milieu − t
            bool const keep_a = upper ? c < 0 : c > 0;                      // t entre V(x_a) et le milieu
            return keep_a ? a : step(a, 1);
        }

        constexpr T inverse(double t, bool upper) const
        {
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::probability::quantile : probabilité hors de [0, 1]");
            if (t == 0)
                return upper ? std::numeric_limits<T>::infinity() : T(0);
            if (t == 1)
                return upper ? T(0) : std::numeric_limits<T>::infinity();
            if (t <= 0.5)
                return static_cast<T>(solve(t, upper));
            return static_cast<T>(solve(1.0 - t, !upper));                  // 1 − t exact
        }

        T mu_, lambda_;
    };
}
