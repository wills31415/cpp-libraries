#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Loi bêta de paramètres α > 0 et β > 0 : densité x^{α−1} (1 − x)^{β−1} / B(α, β)
// sur [0, 1].
//
//   wsl::probability::laws::beta X({.alpha = 2, .beta = 5});
//   X.mean(); X.a(); X.b(); X.pdf(0.3); X.log_pdf(0.3); X(g);
//
// (Une fonction membre ne peut pas s'appeler comme sa classe : les paramètres se
// lisent par a() et b().)
//
// Densité : comme dbeta de R. Si α > 2 et β > 2, par la binomiale généralisée,
//   (n + 1) e^{δ(n) − δ(k) − δ(n−k) − bd0(k, nx) − bd0(n−k, n(1−x))} √(n / (2π k (n−k))),
// k = α − 1, n = α + β − 2 (décomposition de Loader, detail/loader.hpp) ; sinon
// directement, (α − 1) log x + (β − 1) log1p(−x) − log B(α, β), log B en
// double-double sans compensation (detail::log_beta). Un seul arrondi.
// Répartition, survie, quantile : ibeta, ibetac, ibeta_inv de wsl::special.
// Entropie réécrite avec l'erreur de Stirling δ (n = α + β) :
//     ½ log(αβ/n³) + ½ log 2π + ½(1 − 1/α − 1/β + 2/n) + δ(α) + δ(β) − δ(n)
//     − (α − 1)δ'(α) − (β − 1)δ'(β) + (n − 2)δ'(n),
// sans les compensations de la forme usuelle (log B − (α − 1)ψ(α) − …).
//
// Tirage (approché, bit-exact) : X = G₁/(G₁ + G₂), G₁, G₂ de lois gamma(α) et
// gamma(β) (Marsaglia–Tsang), en logarithmes si α < 1 ou β < 1.

#include <wsl/probability/laws/basic.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/probability/detail/loader.hpp>
#include <wsl/probability/laws/gamma.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/random/detail/portable_math.hpp>

#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
}

// Outils internes : dans wsl::probability::detail (un seul « detail » visible depuis probability).
namespace wsl::inline WSL_DEFAULTS::probability::detail
{
        // Entropie de la loi bêta(a, b), en double-double (voir beta::entropy).
        constexpr sd::double_double beta_entropy(double a, double b) noexcept
        {
            namespace fp = sd;
            auto const [n_hi, n_lo] = fp::two_sum(a, b);
            fp::double_double const n{n_hi, n_lo};
            auto const half = [](fp::double_double v) { return fp::double_double{fp::fmul(v.hi, 0.5), fp::fmul(v.lo, 0.5)}; };
            auto const negated = [](fp::double_double v) { return fp::double_double{-v.hi, -v.lo}; };
            auto const reciprocal = [](fp::double_double v) {
                double const q = 1.0 / v.hi;
                auto const [m, m_error] = fp::two_product(q, v.hi);
                return fp::double_double{q, (((1.0 - m) - m_error) - fp::fmul(q, v.lo)) / v.hi};
            };
            auto const ln = detail::log(detail::normalized(n_hi, n_lo));
            fp::double_double h = half(detail::log_precise(a));
            h = detail::sum(h, half(detail::log_precise(b)));
            h = detail::sum(h, detail::mul(ln, {-1.5, 0.0}));                  // −(3/2) log n : 1,5 log n n'est pas exact
            h = detail::sum(h, {fp::half_log_two_pi_hi, fp::half_log_two_pi_lo});
            fp::double_double c{1.0, 0.0};                                  // 1 − 1/a − 1/b + 2/n
            c = detail::sum(c, negated(reciprocal({a, 0.0})));
            c = detail::sum(c, negated(reciprocal({b, 0.0})));
            fp::double_double const two_over_n = reciprocal(n);
            c = detail::sum(c, {fp::fmul(two_over_n.hi, 2.0), fp::fmul(two_over_n.lo, 2.0)});
            h = detail::sum(h, half(c));
            h = detail::sum(h, detail::stirling_error(a));
            h = detail::sum(h, detail::stirling_error(b));
            h = detail::sum(h, negated(detail::stirling_error(n_hi, n_lo)));
            auto const [am_hi, am_lo] = fp::two_sum(a, -1.0);
            auto const [bm_hi, bm_lo] = fp::two_sum(b, -1.0);
            h = detail::sum(h, negated(detail::mul({am_hi, am_lo}, detail::stirling_error_derivative(a))));
            h = detail::sum(h, negated(detail::mul({bm_hi, bm_lo}, detail::stirling_error_derivative(b))));
            h = detail::sum(h, detail::mul(detail::sum(n, {-2.0, 0.0}), detail::stirling_error_derivative(n_hi)));
            return h;
        }
}

namespace wsl::inline WSL_DEFAULTS::probability::laws
{

    template <law_real T = double>
    class beta : public law_interface
    {
    public:
        using value_type = T;
        static constexpr exactness sampling = exactness::approximate;

        struct shapes
        {
            T alpha, beta;
        };
        using param_type = shapes;
        constexpr param_type param() const { return {.alpha = a_, .beta = b_}; }   // pas d'accesseur beta() : c'est le nom de la classe

        constexpr beta(shapes p) : a_(p.alpha), b_(p.beta)
        {
            namespace fp = detail::sd;
            if (!(a_ > 0) || !detail::is_finite(a_) || !(b_ > 0) || !detail::is_finite(b_))
                throw std::domain_error("wsl::probability::laws::beta : il faut alpha > 0 et beta > 0, finis");
            double const a = a_ >= 1 ? double(a_) : double(a_) + 1, b = b_ >= 1 ? double(b_) : double(b_) + 1;
            d_a_ = a - 1.0 / 3;
            c_a_ = 1 / wsl::random::detail::pm::sqrt(fp::fmul(9.0, d_a_));
            d_b_ = b - 1.0 / 3;
            c_b_ = 1 / wsl::random::detail::pm::sqrt(fp::fmul(9.0, d_b_));
        }

        constexpr T a() const noexcept { return a_; }
        constexpr T b() const noexcept { return b_; }
        constexpr interval<T> support() const noexcept { return {T(0), T(1)}; }

        // α/(α + β).
        constexpr T mean() const noexcept { return static_cast<T>(detail::rounded(detail::normalized(a_, 0.0) / total())); }
        // αβ / ((α + β)² (α + β + 1)).
        constexpr T variance() const noexcept
        {
            detail::scaled const s = total();
            auto const [s_hi, s_lo] = detail::sd::two_sum(double(a_), double(b_));
            auto const [u_hi, u_lo] = detail::sum({s_hi, s_lo}, {1.0, 0.0});
            return static_cast<T>(detail::rounded(detail::normalized(a_, 0.0) * detail::normalized(b_, 0.0) / (s * s * detail::normalized(u_hi, u_lo))));
        }
        // (α − 1)/(α + β − 2) si α, β > 1 ; 0 ou 1 au bord ; NaN si le mode n'est pas unique.
        constexpr T mode() const noexcept
        {
            double const a = a_, b = b_;
            if (a > 1 && b > 1)
            {
                auto const [p1, e1] = detail::sd::two_sum(a, -1.0);
                auto const [p2, e2] = detail::sd::two_sum(b, -1.0);
                auto const [n_hi, n_lo] = detail::sum({p1, e1}, {p2, e2});
                return static_cast<T>(detail::rounded(detail::normalized(p1, e1) / detail::normalized(n_hi, n_lo)));
            }
            if (a <= 1 && b > 1)
                return 0;
            if (a > 1 && b <= 1)
                return 1;
            return std::numeric_limits<T>::quiet_NaN();
        }
        constexpr T skewness() const noexcept
        {
            double const a = a_, b = b_, s = a + b;
            return static_cast<T>(2 * (b - a) * wsl::random::detail::pm::sqrt(s + 1) / ((s + 2) * wsl::random::detail::pm::sqrt(a * b)));
        }
        constexpr T excess_kurtosis() const noexcept
        {
            double const a = a_, b = b_, s = a + b;
            double const d = a - b;
            namespace fp = detail::sd;
            double const numerator = fp::fmul(fp::fmul(d, d), s + 1) - fp::fmul(fp::fmul(a, b), s + 2);
            return static_cast<T>(6 * numerator / (a * b * (s + 2) * (s + 3)));
        }

        constexpr T entropy() const noexcept
        {
            auto const [hi, lo] = detail::beta_entropy(double(a_), double(b_));
            return static_cast<T>(hi + lo);
        }

        constexpr T cdf(T x) const noexcept
        {
            auto const& X = *this;
            if (!(x == x))
                return x;
            return static_cast<T>(x <= 0 ? 0.0 : x >= 1 ? 1.0 : wsl::random::detail::pm::ibeta(double(X.a_), double(X.b_), double(x)));
        }
        constexpr T ccdf(T x) const noexcept
        {
            auto const& X = *this;
            if (!(x == x))
                return x;
            return static_cast<T>(x <= 0 ? 1.0 : x >= 1 ? 0.0 : wsl::random::detail::pm::ibetac(double(X.a_), double(X.b_), double(x)));
        }
        constexpr T quantile(double p) const
        {
            auto const& X = *this;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            return static_cast<T>(wsl::random::detail::pm::ibeta_inv(double(X.a_), double(X.b_), p));
        }

        // Quantile de survie : x tel que P(X > x) = q (q exact, sans passer par 1 − q).
        constexpr T cquantile(double q) const
        {
            auto const& X = *this;
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::probability::cquantile : q doit être dans [0, 1]");
            if (q >= 0.5)
                return wsl::probability::quantile(X, 1.0 - q);                                // 1 − q exact
            return static_cast<T>(wsl::random::detail::pm::ibetac_inv(double(X.a_), double(X.b_), q));
        }

        constexpr T pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x, a = X.a_, b = X.b_;
            if (!(v >= 0 && v <= 1))
                return v == v ? T(0) : x;
            if (v == 0 || v == 1)
                return X.edge(v == 0 ? a : b, v == 0 ? b : a, false);
            auto const [e_hi, e_lo, factor] = X.density_exponent(v);
            if (e_hi < -1500)
                return 0;
            if (e_hi > 1400)
                return std::numeric_limits<T>::infinity();
            return static_cast<T>(detail::rounded(detail::exp_scaled(e_hi, e_lo) * factor));
        }

        constexpr T log_pdf(T x) const noexcept
        {
            auto const& X = *this;
            double const v = x, a = X.a_, b = X.b_;
            if (!(v >= 0 && v <= 1))
                return v == v ? -std::numeric_limits<T>::infinity() : x;
            if (v == 0 || v == 1)
                return X.edge(v == 0 ? a : b, v == 0 ? b : a, true);
            auto const [e_hi, e_lo, factor] = X.density_exponent(v);
            auto const e = detail::sum({e_hi, e_lo}, detail::log(factor));
            return static_cast<T>(e.hi + e.lo);
        }

        template <std::uniform_random_bit_generator G>
        constexpr T sample(G& g) const
        {
            auto const& X = *this;
            double const a = X.a_, b = X.b_;
            double const ga = detail::marsaglia_tsang(X.d_a_, X.c_a_, g);
            double const gb = detail::marsaglia_tsang(X.d_b_, X.c_b_, g);
            if (a >= 1 && b >= 1)
                return static_cast<T>(ga / (ga + gb));
            // log G = log G(s + 1) + log(U)/s pour s < 1, puis 1/(1 + e^{log G₂ − log G₁}).
            double la = wsl::random::detail::pm::log(ga), lb = wsl::random::detail::pm::log(gb);
            if (a < 1)
                la = la + wsl::random::detail::pm::log(detail::uniform_positive(g)) / a;
            if (b < 1)
                lb = lb + wsl::random::detail::pm::log(detail::uniform_positive(g)) / b;
            return static_cast<T>(1 / (1 + wsl::random::detail::pm::exp(lb - la)));
        }

        friend constexpr bool operator==(beta const& p, beta const& q) noexcept { return p.a_ == q.a_ && p.b_ == q.b_; }

    private:
        struct exponent_and_factor
        {
            double hi, lo;
            detail::scaled factor;
        };

        // pdf = e^{hi + lo} · factor.
        constexpr exponent_and_factor density_exponent(double x) const noexcept
        {
            namespace fp = detail::sd;
            double const a = a_, b = b_;
            if (a > 2 && b > 2)                                             // binomiale généralisée
            {
                double const k = a - 1, m = b - 1;                          // exacts sous 2^53
                auto const [n_hi, n_lo] = fp::two_sum(k, m);
                auto const [q_hi, q_lo] = fp::two_sum(1.0, -x);
                auto const np = detail::mul({n_hi, n_lo}, {x, 0.0});
                auto const nq = detail::mul({n_hi, n_lo}, {q_hi, q_lo});
                auto e = detail::stirling_error(n_hi);
                auto const dk = detail::stirling_error(k), dm = detail::stirling_error(m);
                auto const bk = detail::deviance(k, np.hi, np.lo), bm = detail::deviance(m, nq.hi, nq.lo);
                e = detail::sum(e, {-dk.hi, -dk.lo});
                e = detail::sum(e, {-dm.hi, -dm.lo});
                e = detail::sum(e, {-bk.hi, -bk.lo});
                e = detail::sum(e, {-bm.hi, -bm.lo});
                detail::scaled const n = detail::normalized(n_hi, n_lo);
                auto const [u_hi, u_lo] = detail::sum({n_hi, n_lo}, {1.0, 0.0});
                detail::scaled const factor = detail::normalized(u_hi, u_lo)
                                            * detail::sqrt(n / (detail::normalized(k, 0.0) * detail::normalized(m, 0.0)))
                                            * detail::normalized(detail::inv_sqrt_2pi_hi, detail::inv_sqrt_2pi_lo);
                return {e.hi, e.lo, factor};
            }
            // (α − 1) log x + (β − 1) log1p(−x) − log B(α, β)
            auto const [am_hi, am_lo] = fp::two_sum(a, -1.0);
            auto const [bm_hi, bm_lo] = fp::two_sum(b, -1.0);
            auto e = detail::log_beta(a, b);
            e = {-e.hi, -e.lo};
            if (am_hi != 0)
                e = detail::sum(e, detail::mul({am_hi, am_lo}, detail::log_precise(x)));
            if (bm_hi != 0)
                e = detail::sum(e, detail::mul({bm_hi, bm_lo}, detail::log1p_dd(-x)));
            return {e.hi, e.lo, detail::normalized(1.0, 0.0)};
        }

        // Densité en 0 (s = α, t = β) ou en 1 (s = β, t = α).
        constexpr T edge(double s, double t, bool logarithm) const noexcept
        {
            if (s < 1)
                return std::numeric_limits<T>::infinity();
            if (s > 1)
                return logarithm ? -std::numeric_limits<T>::infinity() : T(0);
            if (logarithm)                                                  // densité t
            {
                auto const [hi, lo] = detail::log_precise(t);
                return static_cast<T>(hi + lo);
            }
            return static_cast<T>(t);
        }

        // α + β, à l'échelle.
        constexpr detail::scaled total() const noexcept
        {
            auto const [s, s_error] = detail::sd::two_sum(double(a_), double(b_));
            return detail::normalized(s, s_error);
        }

        T a_, b_;
        double d_a_ = 0, c_a_ = 0, d_b_ = 0, c_b_ = 0;                      // Marsaglia–Tsang
    };
}
