#pragma once

// Loi normale asymétrique (Azzalini) : densité (2/ω) φ(z) Φ(αz), z = (x − ξ)/ω.
//
//   wsl::random::skew_normal X({.location = 0, .scale = 1, .shape = 4});
//   X.pdf(0.5); X.cdf(-1); X.quantile(0.99); X.mean(); X(g);
//
// Densité : φ(z) et Φ(αz) à l'échelle (z en double-double, Φ par erfc de wsl::special), un
// arrondi. Répartition : Φ(z) − 2T(z, α) (T d'Owen) se compense dans la queue (α = 1 :
// F = Φ², calculée Φ − Φ(1 − Φ)) ; on intègre donc la plus petite queue de la densité, du
// côté de z par rapport à la moyenne, par morceaux et Gauss–Legendre à 12 nœuds (ceux de
// von Mises) : largeur <= 1, 4/|(log f)′|, 1/√|(log f)″| (erreur de la règle sous 2^-78), et
// max(1/|α|, |t|/2), car Φ(αt) bascule sur une échelle 1/|α| autour de 0 ; la densité est log-concave, la
// queue restante est bornée par f/|(log f)′| au-delà du mode : arrêt quand elle devient
// négligeable. L'autre part par complément. Quantiles : Newton sur le logarithme de la
// part visée, encadré, puis deux pas à résidu précis.
// Moments : δ = α/√(1 + α²), moyenne ξ + ωδ√(2/π), variance ω²(1 − 2δ²/π), asymétrie et
// aplatissement en forme close.
//
// Tirage (approché, bit-exact) : U, V normales exactes ; X = ξ + ωU si V < αU, ξ − ωU sinon.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/trigonometry.hpp>
#include <wsl/random/detail/quadrature.hpp>
#include <wsl/random/detail/root.hpp>
#include <wsl/random/detail/gaussian.hpp>
#include <wsl/random/detail/discrete.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/special.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    namespace detail
    {

        // f(t) = 2 φ(t) Φ(αt), t en double-double, à l'échelle.
        constexpr scaled skew_normal_density(double alpha, sd::double_double t) noexcept
        {
            auto const t2 = mul(t, t);
            double const e_hi = sd::fmul(t2.hi, -0.5), e_lo = sd::fmul(t2.lo, -0.5);
            if (!(e_hi > -1500))
                return {0.0, 0.0, 0};
            scaled const phi = normal_lower_scaled(mul({alpha, 0.0}, t));
            if (phi.hi == 0)
                return phi;
            scaled g = exp_scaled(e_hi, e_lo) * normalized(inv_sqrt_2pi_hi, inv_sqrt_2pi_lo) * phi;
            g.e += 1;
            return g;
        }

        // (log f)′(t) = −t + α φ(αt)/Φ(αt), (log f)″ (valeurs approchées, pour les pas).
        struct log_slope
        {
            double first, second;
        };
        constexpr log_slope skew_normal_slope(double alpha, double t) noexcept
        {
            double const u = sd::fmul(alpha, t);
            double r = 0;                                                   // φ(u)/Φ(u)
            if (u < -30)
                r = -u - 1 / u;
            else
            {
                double const big_phi = 0.5 * wsl::special::erfc(-sd::fmul(u, inv_sqrt2_hi));
                r = big_phi > 0 ? sd::fmul(inv_sqrt_2pi_hi, wsl::special::exp(sd::fmul(-0.5, sd::fmul(u, u)))) / big_phi : -u;
            }
            double const first = sd::fmul(alpha, r) - t;
            double const second = -1 - sd::fmul(sd::fmul(alpha, alpha), sd::fmul(r, u + r));
            return {first, second};
        }

        // ∫ f sur la queue issue de z (vers −∞ si lower, vers +∞ sinon), à l'échelle.
        constexpr scaled skew_normal_tail(double alpha, sd::double_double z, bool lower) noexcept
        {
            namespace fp = sd;
            auto const f = [alpha](sd::double_double t) { return skew_normal_density(alpha, t); };
            auto const width = [alpha](sd::double_double t0) {
                auto const [first, second] = skew_normal_slope(alpha, t0.hi);
                double w = 1.0;
                if (fp::fmul(magnitude(first), w) > 4.0)
                    w = 4.0 / magnitude(first);
                if (fp::fmul(magnitude(second), fp::fmul(w, w)) > 1.0)
                    w = 1.0 / wsl::special::sqrt(magnitude(second));
                if (alpha != 0)                                             // Φ(αt) bascule sur 1/|α| autour de 0
                {
                    double const near = 1.0 / magnitude(alpha), half_distance = fp::fmul(magnitude(t0.hi), 0.5);
                    double const limit = near > half_distance ? near : half_distance;
                    if (w > limit)
                        w = limit;
                }
                return w;
            };
            auto const stop = [&](sd::double_double t0, scaled const& total, sd::double_double) {
                if (total.hi == 0)
                    return false;
                scaled const level = f(t0);
                if (level.hi == 0)
                    return true;
                double const first = skew_normal_slope(alpha, t0.hi).first;
                double const outward = lower ? first : -first;              // > 0 : f décroît vers l'extérieur
                return outward > 0 && (level / outward).e < total.e - 80;   // reste <= f(t₀)/|(log f)′(t₀)|
            };
            auto const piece = [&](sd::double_double lo, sd::double_double hi, scaled const&) { return legendre12(f, lo, hi); };
            sd::double_double const end{lower ? -60.0 : 60.0, 0.0};        // e^{−t²/2} < 2^-2500 au-delà
            return march_integral(z, end, !lower, width, piece, stop);
        }
    }

    template <law_real T = double>
    class skew_normal : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            T location = 0, scale = 1, shape = 0;
        };

        constexpr skew_normal() : skew_normal(parameters{}) {}

        constexpr skew_normal(parameters p) : xi_(p.location), omega_(p.scale), alpha_(p.shape)
        {
            if (!detail::is_finite(xi_) || !(omega_ > 0) || !detail::is_finite(omega_) || !detail::is_finite(alpha_))
                throw std::domain_error("wsl::random::skew_normal : il faut location et shape finis, scale > 0 fini");
            // δ = α/√(1 + α²) et δ√(2/π), en double-double
            namespace sd = detail::sd;
            double const a = alpha_;
            if (a != 0)
            {
                auto const a2 = detail::mul({a, 0.0}, {a, 0.0});
                auto const one_a2 = detail::sum({1.0, 0.0}, a2);
                detail::scaled const root = detail::sqrt(detail::normalized(one_a2.hi, one_a2.lo));
                auto const d = detail::unscaled(detail::normalized(detail::magnitude(a), 0.0) / root);
                delta_ = a < 0 ? detail::negated(d) : d;
            }
            // √(2/π) = 2 · (1/√(2π))
            auto const m = detail::mul(delta_, {sd::fmul(detail::inv_sqrt_2pi_hi, 2.0), sd::fmul(detail::inv_sqrt_2pi_lo, 2.0)});
            mu_z_ = m;
            // 1 − 2δ²/π = 1 − μ_z²
            variance_z_ = detail::sum({1.0, 0.0}, detail::negated(detail::mul(m, m)));
        }

        constexpr T location() const noexcept { return xi_; }
        constexpr T scale() const noexcept { return omega_; }
        constexpr T shape() const noexcept { return alpha_; }
        constexpr interval<T> support() const noexcept { return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()}; }
        constexpr T mean() const noexcept
        {
            if (mu_z_.hi == 0)
                return xi_;
            return static_cast<T>(detail::affine(double(xi_), double(omega_), mu_z_.hi, mu_z_.lo));
        }
        constexpr T variance() const noexcept
        {
            detail::scaled const w = detail::normalized(double(omega_), 0.0);
            return static_cast<T>(detail::rounded(w * w * detail::normalized(variance_z_.hi, variance_z_.lo)));
        }
        // ((4 − π)/2) μ_z³/σ_z³
        constexpr T skewness() const noexcept
        {
            if (mu_z_.hi == 0)
                return 0;
            detail::scaled const m = detail::normalized(detail::magnitude(mu_z_.hi), mu_z_.hi < 0 ? -mu_z_.lo : mu_z_.lo);
            detail::scaled const v = detail::normalized(variance_z_.hi, variance_z_.lo);
            double const r = detail::rounded(detail::normalized(detail::four_minus_pi_half_hi, detail::four_minus_pi_half_lo) * m * m * m / (v * detail::sqrt(v)));
            return static_cast<T>(mu_z_.hi < 0 ? -r : r);
        }
        // 2(π − 3) μ_z⁴/σ_z⁴
        constexpr T excess_kurtosis() const noexcept
        {
            if (mu_z_.hi == 0)
                return 0;
            detail::scaled const m = detail::normalized(detail::magnitude(mu_z_.hi), mu_z_.hi < 0 ? -mu_z_.lo : mu_z_.lo);
            detail::scaled const v = detail::normalized(variance_z_.hi, variance_z_.lo);
            auto const c = detail::sum({detail::pi_hi, detail::pi_lo}, {-3.0, 0.0});    // π − 3
            detail::scaled k = detail::normalized(c.hi, c.lo) * m * m * m * m / (v * v);
            k.e += 1;
            return static_cast<T>(detail::rounded(k));
        }

        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x;
            if (!(v == v))
                return x;
            auto const z = detail::standardize(v, double(X.xi_), double(X.omega_));
            detail::scaled const f = detail::skew_normal_density(double(X.alpha_), z);
            if (f.hi == 0)
                return 0;
            return static_cast<T>(detail::rounded(f / double(X.omega_)));
        }
        // log 2 − z²/2 − ½ log 2π + log Φ(αz) − log ω
        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            namespace sd = detail::sd;
            double const v = x;
            if (!(v == v))
                return x;
            auto const z = detail::standardize(v, double(X.xi_), double(X.omega_));
            auto const z2 = detail::mul(z, z);
            detail::scaled const phi = detail::normal_lower_scaled(detail::mul({double(X.alpha_), 0.0}, z));
            if (phi.hi == 0 || !detail::is_finite(z2.hi))
                return -std::numeric_limits<T>::infinity();
            auto l = detail::sum({sd::ln2_hi, sd::ln2_lo}, {sd::fmul(z2.hi, -0.5), sd::fmul(z2.lo, -0.5)});
            l = detail::sum(l, {-detail::half_log_2pi_hi, -detail::half_log_2pi_lo});
            l = detail::sum(l, detail::log(phi));
            auto const lw = detail::log_precise(double(X.omega_));
            l = detail::sum(l, detail::negated(lw));
            return static_cast<T>(l.hi + l.lo);
        }
        constexpr T cdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, false); }
        constexpr T ccdf(T x) const noexcept { auto const& X = *this; return X.distribution(x, true); }
        constexpr T quantile(double p) const { auto const& X = *this; return X.inverse(p, false); }
        constexpr T cquantile(double q) const { auto const& X = *this; return X.inverse(q, true); }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            normal<double> const N({.mean = 0.0, .stddev = 1.0});
            double const u = wsl::random::sample(N, g), v = wsl::random::sample(N, g);
            double const z = v < detail::sd::fmul(double(X.alpha_), u) ? u : -u;
            return static_cast<T>(double(X.xi_) + detail::sd::fmul(double(X.omega_), z));
        }

        friend constexpr bool operator==(skew_normal const& a, skew_normal const& b) noexcept { return a.xi_ == b.xi_ && a.omega_ == b.omega_ && a.alpha_ == b.alpha_; }

    private:
        // (F(z), S(z)) à l'échelle : la queue du côté de z par rapport à la moyenne, l'autre par complément.
        struct parts
        {
            detail::scaled lower, upper;
        };
        constexpr parts both(detail::sd::double_double z) const noexcept
        {
            bool const lower = z.hi < mu_z_.hi || (z.hi == mu_z_.hi && z.lo < mu_z_.lo);
            detail::scaled const tail = detail::skew_normal_tail(double(alpha_), z, lower);
            detail::scaled other{1.0, 0.0, 0};
            if (tail.hi != 0)
            {
                auto const [t_hi, t_lo] = detail::unscaled(tail);
                auto const [s, s_error] = detail::sd::two_sum(1.0, -t_hi);
                auto const [h, l] = detail::sd::fast_two_sum(s, s_error - t_lo);
                other = h > 0 ? detail::normalized(h, l) : detail::scaled{0.0, 0.0, 0};
            }
            return lower ? parts{tail, other} : parts{other, tail};
        }

        constexpr T distribution(T x, bool upper) const noexcept
        {
            double const v = x;
            if (!(v == v))
                return x;
            if (v == std::numeric_limits<double>::infinity() || v == -std::numeric_limits<double>::infinity())
                return (v > 0) != upper ? T(1) : T(0);
            auto const z = detail::standardize(v, double(xi_), double(omega_));
            if (!detail::is_finite(z.hi))
                return (z.hi > 0) != upper ? T(1) : T(0);
            auto const p = both(z);
            return static_cast<T>(detail::rounded(upper ? p.upper : p.lower));
        }

        // z tel que F(z) = t (upper = false) ou S(z) = t, t <= ½ : solveur commun
        // (detail/root.hpp), départ μ_z + σ_z Φ⁻¹, pas d'au plus 8, ±4 hors encadrement.
        constexpr detail::sd::double_double solve(double t, bool upper) const
        {
            namespace sd = detail::sd;
            double const a = alpha_;
            auto const evaluate = [&](sd::double_double z) {
                auto const p = both(z);
                detail::part_value r{upper ? p.upper : p.lower, {0.0, 0.0, 0}, z};
                if (r.value.hi != 0)
                    r.slope = detail::skew_normal_density(a, z);
                return r;
            };
            double const sigma = wsl::special::sqrt(variance_z_.hi);
            double const q = detail::standard_quantile(t).hi;
            double const start = mu_z_.hi + sd::fmul(sigma, upper ? -q : q);
            return detail::solve_part(evaluate, !upper, t, {.start = start, .max_step = 8},
                                      [](double z, int side) { return side > 0 ? z - 4 : z + 4; }, detail::arithmetic_midpoint, 300);
        }

        constexpr T inverse(double t, bool upper) const
        {
            if (!(t >= 0 && t <= 1))
                throw std::domain_error("wsl::random::quantile : probabilité hors de [0, 1]");
            if (t == 0)
                return upper ? std::numeric_limits<T>::infinity() : -std::numeric_limits<T>::infinity();
            if (t == 1)
                return upper ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity();
            auto const z = t <= 0.5 ? solve(t, upper) : solve(1.0 - t, !upper);
            if (z.hi == 0)
                return xi_;
            return static_cast<T>(detail::affine(double(xi_), double(omega_), z.hi, z.lo));
        }

        T xi_, omega_, alpha_;
        detail::sd::double_double delta_{0.0, 0.0}, mu_z_{0.0, 0.0}, variance_z_{1.0, 0.0};
    };
}
