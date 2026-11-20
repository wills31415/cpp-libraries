#pragma once

// Processus en temps continu (docs/random.md §14), exacts point par point (option B) :
// chaque valeur rendue est l'arrondi correct d'une trajectoire tirée exactement selon la
// loi dont les paramètres irrationnels (√Δt, e^{−θt}…) sont arrondis une fois.
//
//   wsl::random::brownian_motion B({.drift = 0, .volatility = 1});
//   B.at(2.0);                                             // normal(0, √2)
//   auto x = wsl::random::sample(B.on_grid({0.5, 1.0, 2.0}), g);        // B aux instants donnés
//   for (auto [t, v] : wsl::random::trajectory(B, g, {.step = 0.01}) | std::views::take(100)) …
//   wsl::random::brownian_bridge D({.from = {0, 0}, .to = {1, 0}});  wsl::random::sample(D.on_grid(ts), g);
//   wsl::random::ornstein_uhlenbeck O({.theta = 1, .mean = 0, .volatility = 1, .start = 2});
//   wsl::random::geometric_brownian_motion S({.drift = 0.05, .volatility = 0.2, .start = 100});
//   wsl::random::poisson_process N({.rate = 2});    N.at(3.0);   wsl::random::sample(N.on(10.0), g);
//   for (double t : wsl::random::trajectory(N, g) | std::views::take(10)) …   // instants des sauts
//   wsl::random::inhomogeneous_poisson_process M(f, /* borne */ 5.0);
//   wsl::random::compound_poisson C(2.0, wsl::random::exponential({.rate = 1.0}));     // (instant, valeur)
//   wsl::random::continuous_markov_chain Q({.generator = {{-1, 1}, {2, -2}}, .start = 0});
//
// Construction (detail/lazy_sum.hpp) :
// - brownien : X(t) = x₀ + μ t + σ Σ √Δt_j Z_j, Z_j normales exactes ; μ t exact
//   (produit de doubles), √Δt_j arrondi une fois (Δt_j exact) ;
// - pont : a + r (b − a) + σ (W(t) − r W(T)), r = (t − t₀)/(T − t₀) arrondi une fois ;
// - Ornstein–Uhlenbeck, par changement de temps : X(t) = μ + a(t)(x₀ − μ) + r(t) W(u(t)),
//   a = e^{−θt}, r = σ a/√(2θ), u = e^{2θt} − 1 (arrondis une fois par instant ; θt <= 350) ;
// - brownien géométrique : s₀ e^{X(t)}, X brownien de dérive μ − σ²/2 exact, l'exponentielle
//   arrondie par wsl::special::exp (approché : 0,5 ulp, pas d'arrondi correct garanti) ;
// - Poisson : Tₖ = Σ Eⱼ/λ, exponentielles exactes ; inhomogène : amincissement, candidat
//   accepté avec probabilité f(t)/Λ exactement, f évaluée à l'instant arrondi (approché) ;
//   composé : instants exacts, valeurs exactes si la loi des sauts se tire paresseusement
//   (normale, exponentielle, Laplace, uniforme), sinon somme exacte des sauts arrondis ;
// - chaîne de Markov en temps continu : séjours Eⱼ/q_x exacts (q_x = −Q_xx), sauts
//   catégoriels exacts ; temps +∞ après absorption.
// Coût : une valeur coûte O(1) en moyenne (voie rapide) ; la mémoire croît avec la
// longueur de la trajectoire (les accroissements restent raffinables).

#include <wsl/random/affine.hpp>
#include <wsl/random/categorical.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/detail/lazy_sum.hpp>
#include <wsl/random/exponential.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/lognormal.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/random/poisson.hpp>
#include <wsl/random/process.hpp>
#include <wsl/special.hpp>

#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wsl::random
{
    struct grid_step
    {
        double step;
    };

    namespace detail
    {
        // √x arrondi une fois, x dyadique exact > 0.
        constexpr double sqrt_of(dyadic_value const& x)
        {
            auto const d = to_dd(x);
            return rounded(sqrt(normalized(d.hi, d.lo)));
        }

        // ±Z, Z normale centrée réduite exacte (demi-normale paresseuse, puis le signe).
        template <class G>
        constexpr lazy_draw standard_normal(G& g)
        {
            return lazy_sample(normal<double>({.mean = 0.0, .sd = 1.0}), g);
        }

        // Brownien standard W aux instants croissants donnés (W(0) = 0) : accroissements
        // √Δ Z dans une somme paresseuse.
        struct brownian_sum
        {
            lazy_sum sum{};
            dyadic_value last{};                                            // dernier instant (exact)

            template <class G>
            constexpr void advance_to(dyadic_value const& t, G& g)
            {
                dyadic_value const delta = t - last;
                if (delta.is_zero())
                    return;
                if (delta.m.negative)
                    throw std::domain_error("wsl : instants décroissants");
                lazy_draw z = standard_normal(g);
                dyadic_value const c = z.coefficient * to_dyadic(sqrt_of(delta));
                sum.add(c, 0, std::move(z.v), g);
                last = t;
            }
        };

        constexpr void check_times(std::vector<double> const& times)
        {
            for (std::size_t i = 0; i < times.size(); ++i)
                if (!(times[i] >= 0) || !is_finite(times[i]) || (i > 0 && times[i] < times[i - 1]))
                    throw std::domain_error("wsl : instants >= 0, finis, croissants");
        }
    }

    // -----------------------------------------------------------------------
    // Mouvement brownien X(t) = x₀ + μ t + σ W(t).
    // -----------------------------------------------------------------------
    class brownian_motion
    {
    public:
        using state_type = std::pair<double, double>;                       // (instant, valeur)

        struct parameters
        {
            double drift = 0;
            double volatility = 1;
            double start = 0;
        };

        constexpr explicit brownian_motion(parameters p) : mu_(p.drift), sigma_(p.volatility), x0_(p.start)
        {
            if (!detail::is_finite(mu_) || !(sigma_ > 0) || !detail::is_finite(sigma_) || !detail::is_finite(x0_))
                throw std::domain_error("wsl::random::brownian_motion : dérive et départ finis, volatilité > 0 finie");
        }

        constexpr double drift() const noexcept { return mu_; }
        constexpr double volatility() const noexcept { return sigma_; }
        constexpr double start() const noexcept { return x0_; }
        // X(t) ~ N(x₀ + μ t, σ √t)
        constexpr normal<> at(double t) const
        {
            if (!(t > 0))
                throw std::domain_error("wsl::random::brownian_motion::at : il faut t > 0");
            double const m = detail::linear(mu_, t, x0_);
            return normal<>({.mean = m, .sd = detail::sd::fmul(sigma_, wsl::special::sqrt(t))});
        }

        // Valeurs aux instants donnés (>= 0, croissants).
        class grid_law : public law_interface
        {
        public:
            using value_type = std::vector<double>;
            static constexpr exactness sampling = exactness::exact;

            constexpr grid_law(parameters b, std::vector<double> times) : b_(b), times_(std::move(times)) { detail::check_times(times_); }

            template <std::uniform_random_bit_generator G>
            constexpr value_type sample(G& g) const
            {
                auto const& L = *this;
                detail::brownian_sum w;
                value_type x;
                for (double t : L.times_)
                    x.push_back(value_at(L.b_, w, detail::to_dyadic(t), g));
                return x;
            }

        private:
            parameters b_;
            std::vector<double> times_;
        };
        constexpr grid_law on_grid(std::vector<double> times) const { return grid_law(parameters{mu_, sigma_, x0_}, std::move(times)); }

        // Marcheur à pas régulier h : (k h, X(k h)), k = 0, 1, …
        struct walker_type
        {
            parameters b;
            detail::dyadic_value h;
            double step;
            detail::brownian_sum w{};
            std::uint64_t k = 0;
            template <class G>
            constexpr state_type first(G& g) { return at_k(g); }
            template <class G>
            constexpr state_type next(G& g)
            {
                ++k;
                return at_k(g);
            }
            template <class G>
            constexpr state_type at_k(G& g)
            {
                detail::dyadic_value const t = h * detail::to_dyadic(double(k));
                return {detail::round_rational<double>(t, detail::big(1)), value_at(b, w, t, g)};
            }
        };
        constexpr walker_type walker(grid_step s) const
        {
            if (!(s.step > 0) || !detail::is_finite(s.step))
                throw std::domain_error("wsl::random::brownian_motion : pas > 0 fini");
            return {parameters{mu_, sigma_, x0_}, detail::to_dyadic(s.step), s.step};
        }

        // Variante rapide (wsl::random::fast) : accroissements par la ziggourat, somme double-double.
        struct fast_walker_type
        {
            double mu, sigma, h, c;
            detail::sd::double_double x;
            std::uint64_t k = 0;
            template <class G>
            constexpr state_type first(G&) { return {0.0, x.hi + x.lo}; }
            template <class G>
            constexpr state_type next(G& g)
            {
                ++k;
                double const z = detail::ziggurat_normal(g);
                x = detail::sum(x, detail::sum(detail::product(mu, h), detail::product(c, z)));
                return {detail::sd::fmul(double(k), h), x.hi + x.lo};
            }
        };
        constexpr fast_walker_type fast_walker(grid_step s) const
        {
            return {mu_, sigma_, s.step, detail::sd::fmul(sigma_, wsl::special::sqrt(s.step)), {x0_, 0.0}};
        }

        // x₀ + μ t + σ W(t), arrondi une fois.
        template <class G>
        static constexpr double value_at(parameters const& b, detail::brownian_sum& w, detail::dyadic_value const& t, G& g)
        {
            w.advance_to(t, g);
            detail::dyadic_value const offset = detail::to_dyadic(b.start) + detail::to_dyadic(b.drift) * t;
            return w.sum.template rounded<double>(g, offset, {}, 0, detail::to_dyadic(b.volatility));
        }

    private:
        double mu_, sigma_, x0_;
    };

    template <std::uniform_random_bit_generator G>
    constexpr auto trajectory(brownian_motion const& b, G& g, grid_step s)
    {
        return process_view<brownian_motion::walker_type, G>(b.walker(s), g);
    }
    template <std::uniform_random_bit_generator G>
    constexpr auto trajectory(brownian_motion const& b, G& g, grid_step s, fast_t)
    {
        return process_view<brownian_motion::fast_walker_type, G>(b.fast_walker(s), g);
    }

    // -----------------------------------------------------------------------
    // Pont brownien de (t₀, a) à (T, b), volatilité σ.
    // -----------------------------------------------------------------------
    class brownian_bridge
    {
    public:
        struct point
        {
            double time, value;
        };
        struct parameters
        {
            point from, to;
            double volatility = 1;
        };

        constexpr explicit brownian_bridge(parameters p) : from_(p.from), to_(p.to), sigma_(p.volatility)
        {
            if (!(to_.time > from_.time) || !detail::is_finite(from_.time) || !detail::is_finite(to_.time) || !detail::is_finite(from_.value)
                || !detail::is_finite(to_.value) || !(sigma_ > 0) || !detail::is_finite(sigma_))
                throw std::domain_error("wsl::random::brownian_bridge : il faut t₀ < T, valeurs finies, volatilité > 0");
        }

        // X(t) ~ N(a + r (b − a), σ √((t − t₀)(T − t)/(T − t₀)))
        constexpr normal<> at(double t) const
        {
            if (!(t > from_.time && t < to_.time))
                throw std::domain_error("wsl::random::brownian_bridge::at : il faut t₀ < t < T");
            double const r = (t - from_.time) / (to_.time - from_.time);
            double const v = detail::sd::fmul(detail::sd::fmul(t - from_.time, to_.time - t) / (to_.time - from_.time), detail::sd::fmul(sigma_, sigma_));
            return normal<>({.mean = detail::linear(r, to_.value - from_.value, from_.value), .sd = wsl::special::sqrt(v)});
        }

        class grid_law : public law_interface
        {
        public:
            using value_type = std::vector<double>;
            static constexpr exactness sampling = exactness::exact;

            constexpr grid_law(parameters b, std::vector<double> times) : b_(b), times_(std::move(times))
            {
                for (std::size_t i = 0; i < times_.size(); ++i)
                    if (!(times_[i] >= b_.from.time && times_[i] <= b_.to.time) || (i > 0 && times_[i] < times_[i - 1]))
                        throw std::domain_error("wsl::random::brownian_bridge : instants croissants dans [t₀, T]");
            }

            // W(tᵢ − t₀) pour tous les instants, puis W(T − t₀) ; X(tᵢ) = a + rᵢ (b − a) + σ (Wᵢ − rᵢ W_T).
            template <std::uniform_random_bit_generator G>
            constexpr value_type sample(G& g) const
            {
                auto const& L = *this;
                namespace d = detail;
                d::brownian_sum w;
                d::dyadic_value const t0 = d::to_dyadic(L.b_.from.time), span = d::to_dyadic(L.b_.to.time) - t0;
                std::vector<std::size_t> prefix;
                for (double t : L.times_)
                {
                    w.advance_to(d::to_dyadic(t) - t0, g);
                    prefix.push_back(w.sum.size());
                }
                w.advance_to(span, g);
                d::dyadic_value const a = d::to_dyadic(L.b_.from.value), ba = d::to_dyadic(L.b_.to.value) - a, s = d::to_dyadic(L.b_.volatility);
                value_type x;
                for (std::size_t i = 0; i < L.times_.size(); ++i)
                {
                    double const r = L.times_[i] == L.b_.to.time ? 1.0 : (L.times_[i] - L.b_.from.time) / (L.b_.to.time - L.b_.from.time);
                    d::dyadic_value const rd = d::to_dyadic(r);
                    x.push_back(w.sum.template rounded<double>(g, a + rd * ba, s, prefix[i], -(s * rd)));
                }
                return x;
            }

        private:
            parameters b_;
            std::vector<double> times_;
        };
        constexpr grid_law on_grid(std::vector<double> times) const { return grid_law(parameters{from_, to_, sigma_}, std::move(times)); }

    private:
        point from_, to_;
        double sigma_;
    };

    // -----------------------------------------------------------------------
    // Ornstein–Uhlenbeck dX = θ(μ − X) dt + σ dW, X(0) = x₀.
    // -----------------------------------------------------------------------
    class ornstein_uhlenbeck
    {
    public:
        using state_type = std::pair<double, double>;

        struct parameters
        {
            double theta, mean = 0, volatility = 1, start = 0;
        };

        constexpr explicit ornstein_uhlenbeck(parameters p) : theta_(p.theta), mu_(p.mean), sigma_(p.volatility), x0_(p.start)
        {
            if (!(theta_ > 0) || !detail::is_finite(theta_) || !detail::is_finite(mu_) || !(sigma_ > 0) || !detail::is_finite(sigma_) || !detail::is_finite(x0_))
                throw std::domain_error("wsl::random::ornstein_uhlenbeck : θ > 0, σ > 0, moyenne et départ finis");
        }

        // X(t) ~ N(μ + e^{−θt}(x₀ − μ), σ √(−expm1(−2θt)/(2θ)))
        constexpr normal<> at(double t) const
        {
            namespace sp = wsl::special;
            if (!(t > 0))
                throw std::domain_error("wsl::random::ornstein_uhlenbeck::at : il faut t > 0");
            double const a = sp::exp(-detail::sd::fmul(theta_, t));
            double const v = -sp::expm1(-detail::sd::fmul(detail::sd::fmul(2.0, theta_), t)) / detail::sd::fmul(2.0, theta_);
            return normal<>({.mean = detail::linear(a, x0_ - mu_, mu_), .sd = detail::sd::fmul(sigma_, sp::sqrt(v))});
        }
        constexpr normal<> stationary() const { return normal<>({.mean = mu_, .sd = sigma_ / wsl::special::sqrt(detail::sd::fmul(2.0, theta_))}); }

        class grid_law : public law_interface
        {
        public:
            using value_type = std::vector<double>;
            static constexpr exactness sampling = exactness::exact;

            constexpr grid_law(parameters o, std::vector<double> times) : o_(o), times_(std::move(times)) { detail::check_times(times_); }

            template <std::uniform_random_bit_generator G>
            constexpr value_type sample(G& g) const
            {
                auto const& L = *this;
                detail::brownian_sum w;
                value_type x;
                for (double t : L.times_)
                    x.push_back(value_at(L.o_, w, t, g));
                return x;
            }

        private:
            parameters o_;
            std::vector<double> times_;
        };
        constexpr grid_law on_grid(std::vector<double> times) const { return grid_law(parameters{theta_, mu_, sigma_, x0_}, std::move(times)); }

        struct walker_type
        {
            parameters o;
            double step;
            detail::brownian_sum w{};
            std::uint64_t k = 0;
            template <class G>
            constexpr state_type first(G& g) { return at_k(g); }
            template <class G>
            constexpr state_type next(G& g)
            {
                ++k;
                return at_k(g);
            }
            template <class G>
            constexpr state_type at_k(G& g)
            {
                double const t = detail::sd::fmul(double(k), step);
                return {t, value_at(o, w, t, g)};
            }
        };
        constexpr walker_type walker(grid_step s) const
        {
            if (!(s.step > 0) || !detail::is_finite(s.step))
                throw std::domain_error("wsl::random::ornstein_uhlenbeck : pas > 0 fini");
            return {parameters{theta_, mu_, sigma_, x0_}, s.step};
        }

        // μ + a (x₀ − μ) + r W(u), a = e^{−θt}, r = σ a/√(2θ), u = e^{2θt} − 1.
        template <class G>
        static constexpr double value_at(parameters const& o, detail::brownian_sum& w, double t, G& g)
        {
            double const theta_ = o.theta, mu_ = o.mean, sigma_ = o.volatility, x0_ = o.start;
            namespace sp = wsl::special;
            namespace d = detail;
            double const theta_t = d::sd::fmul(theta_, t);
            if (!(theta_t <= 350))
                throw std::domain_error("wsl::random::ornstein_uhlenbeck : θ t > 350 (changement de temps hors des doubles)");
            if (t == 0)
                return x0_;
            double const u = sp::expm1(d::sd::fmul(2.0, theta_t)), a = sp::exp(-theta_t);
            double const r = d::sd::fmul(sigma_ / sp::sqrt(d::sd::fmul(2.0, theta_)), a);
            w.advance_to(d::to_dyadic(u), g);
            d::dyadic_value const offset = d::to_dyadic(mu_) + d::to_dyadic(a) * (d::to_dyadic(x0_) - d::to_dyadic(mu_));
            return w.sum.template rounded<double>(g, offset, {}, 0, d::to_dyadic(r));
        }

    private:
        double theta_, mu_, sigma_, x0_;
    };

    template <std::uniform_random_bit_generator G>
    constexpr auto trajectory(ornstein_uhlenbeck const& o, G& g, grid_step s)
    {
        return process_view<ornstein_uhlenbeck::walker_type, G>(o.walker(s), g);
    }

    // -----------------------------------------------------------------------
    // Brownien géométrique S(t) = s₀ e^{(μ − σ²/2) t + σ W(t)}.
    // -----------------------------------------------------------------------
    class geometric_brownian_motion
    {
    public:
        using state_type = std::pair<double, double>;

        struct parameters
        {
            double drift = 0, volatility = 1, start = 1;
        };

        constexpr explicit geometric_brownian_motion(parameters p)
            : log_(brownian_motion::parameters{.drift = p.drift - detail::sd::fmul(0.5, detail::sd::fmul(p.volatility, p.volatility)), .volatility = p.volatility, .start = 0}), s0_(p.start)
        {
            if (!(s0_ > 0) || !detail::is_finite(s0_))
                throw std::domain_error("wsl::random::geometric_brownian_motion : départ > 0 fini");
        }

        constexpr brownian_motion const& logarithm() const noexcept { return log_; }   // log(S/s₀), exact
        // S(t) ~ lognormale(log s₀ + (μ − σ²/2) t, σ √t)
        constexpr lognormal<> at(double t) const
        {
            auto const n = log_.at(t);
            return lognormal<>({.mu = wsl::special::log(s0_) + n.mean(), .sigma = n.sd()});
        }

        class grid_law : public law_interface
        {
        public:
            using value_type = std::vector<double>;
            static constexpr exactness sampling = exactness::approximate;

            constexpr grid_law(brownian_motion log, double s0, std::vector<double> times) : log_(std::move(log)), s0_(s0), times_(std::move(times)) {}

            template <std::uniform_random_bit_generator G>
            constexpr value_type sample(G& g) const
            {
                auto const& L = *this;
                auto x = wsl::random::sample(L.log_.on_grid(L.times_), g);
                for (auto& v : x)
                    v = detail::sd::fmul(L.s0_, wsl::special::exp(v));
                return x;
            }

        private:
            brownian_motion log_;
            double s0_;
            std::vector<double> times_;
        };
        constexpr grid_law on_grid(std::vector<double> times) const { return grid_law(log_, s0_, std::move(times)); }

    private:
        brownian_motion log_;
        double s0_;
    };

    // -----------------------------------------------------------------------
    // Processus de Poisson homogène : instants des sauts T₁ < T₂ < …
    // -----------------------------------------------------------------------
    class poisson_process
    {
    public:
        using state_type = double;

        struct parameters
        {
            double rate;
        };

        constexpr explicit poisson_process(parameters p) : lambda_(p.rate)
        {
            if (!(lambda_ > 0) || !detail::is_finite(lambda_))
                throw std::domain_error("wsl::random::poisson_process : taux > 0 fini");
        }

        constexpr double rate() const noexcept { return lambda_; }
        constexpr poisson<> at(double t) const { return poisson<>({.mean = detail::sd::fmul(lambda_, t)}); }   // N(t)

        struct walker_type
        {
            double lambda;
            detail::lazy_sum times{std::vector<double>{lambda}};
            template <class G>
            constexpr double first(G& g) { return next(g); }
            template <class G>
            constexpr double next(G& g)
            {
                times.add(detail::lazy_sample(exponential<double>({.rate = lambda}), g), g);
                return times.rounded(g);
            }
            // T_last <= t exactement ?
            template <class G>
            constexpr bool not_after(G& g, double t) { return times.compare(g, {}, {}, 0, detail::to_dyadic(1.0), detail::to_dyadic(t)) <= 0; }
        };
        constexpr walker_type walker() const { return {lambda_}; }

        struct fast_walker_type
        {
            double lambda;
            detail::sd::double_double t{0.0, 0.0};
            template <class G>
            constexpr double first(G& g) { return next(g); }
            template <class G>
            constexpr double next(G& g)
            {
                t = detail::sum(t, {wsl::random::sample(exponential<double>({.rate = lambda}), g, fast), 0.0});
                return t.hi + t.lo;
            }
        };

        // Instants dans [0, T] : loi d'un ensemble fini.
        class window_law : public law_interface
        {
        public:
            using value_type = std::vector<double>;
            static constexpr exactness sampling = exactness::exact;

            constexpr window_law(double rate, double horizon) : lambda_(rate), horizon_(horizon)
            {
                if (!(horizon_ >= 0) || !detail::is_finite(horizon_))
                    throw std::domain_error("wsl::random::poisson_process::on : horizon >= 0 fini");
            }

            template <std::uniform_random_bit_generator G>
            constexpr value_type sample(G& g) const
            {
                auto const& W = *this;
                walker_type w{W.lambda_};
                value_type times;
                for (;;)
                {
                    double const t = w.next(g);
                    if (!w.not_after(g, W.horizon_))
                        return times;
                    times.push_back(t);
                }
            }

        private:
            double lambda_, horizon_;
        };
        constexpr window_law on(double horizon) const { return {lambda_, horizon}; }

    private:
        double lambda_;
    };

    template <std::uniform_random_bit_generator G>
    constexpr auto trajectory(poisson_process const& p, G& g, fast_t)
    {
        return process_view<poisson_process::fast_walker_type, G>({p.rate()}, g);
    }

    // Poisson d'intensité f(t) <= Λ, par amincissement (candidats exacts de taux Λ,
    // acceptés avec probabilité f(t)/Λ exactement, f évaluée à l'instant arrondi).
    template <class F>
        requires std::invocable<F const&, double>
    class inhomogeneous_poisson_process : public law_interface
    {
    public:
        using state_type = double;
        static constexpr exactness sampling = exactness::approximate;

        constexpr inhomogeneous_poisson_process(F intensity, double bound) : f_(std::move(intensity)), bound_(bound)
        {
            if (!(bound_ > 0) || !detail::is_finite(bound_))
                throw std::domain_error("wsl::random::inhomogeneous_poisson_process : borne > 0 finie");
        }

        struct walker_type
        {
            inhomogeneous_poisson_process const* p;
            poisson_process::walker_type candidates;
            template <class G>
            constexpr double first(G& g) { return next(g); }
            template <class G>
            constexpr double next(G& g)
            {
                for (;;)
                {
                    double const t = candidates.next(g);
                    double const f = p->f_(t);
                    if (!(f >= 0 && f <= p->bound_))
                        throw std::domain_error("wsl::random::inhomogeneous_poisson_process : il faut 0 <= f(t) <= borne");
                    if (detail::fraction(g, f, p->bound_))
                        return t;
                }
            }
        };
        constexpr walker_type walker() const { return {this, poisson_process::walker_type{bound_}}; }

        class window_law : public law_interface
        {
        public:
            using value_type = std::vector<double>;
            static constexpr exactness sampling = exactness::approximate;

            constexpr window_law(inhomogeneous_poisson_process const& p, double horizon) : p_(p), horizon_(horizon) {}

            template <std::uniform_random_bit_generator G>
            constexpr value_type sample(G& g) const
            {
                auto const& W = *this;
                poisson_process::walker_type c{W.p_.bound_};
                value_type times;
                for (;;)
                {
                    double const t = c.next(g);
                    if (!c.not_after(g, W.horizon_))
                        return times;
                    double const f = W.p_.f_(t);
                    if (!(f >= 0 && f <= W.p_.bound_))
                        throw std::domain_error("wsl::random::inhomogeneous_poisson_process : il faut 0 <= f(t) <= borne");
                    if (detail::fraction(g, f, W.p_.bound_))
                        times.push_back(t);
                }
            }

        private:
            inhomogeneous_poisson_process p_;
            double horizon_;
        };
        constexpr window_law on(double horizon) const { return window_law(*this, horizon); }

    private:
        F f_;
        double bound_;
    };

    // -----------------------------------------------------------------------
    // Poisson composé : (instant, somme des sauts) à chaque saut, depuis (0, 0).
    // -----------------------------------------------------------------------
    template <probability_law L>
        requires std::floating_point<typename L::value_type>
    class compound_poisson : public law_interface
    {
    public:
        using state_type = std::pair<double, double>;
        static constexpr exactness sampling = detail::lazy_samplable<L, std::mt19937_64> ? exactness::exact : exactness::approximate;

        constexpr compound_poisson(double rate, L jump) : lambda_(rate), jump_(std::move(jump))
        {
            if (!(lambda_ > 0) || !detail::is_finite(lambda_))
                throw std::domain_error("wsl::random::compound_poisson : taux > 0 fini");
        }

        constexpr double rate() const noexcept { return lambda_; }
        constexpr L const& jump() const noexcept { return jump_; }
        // E[S(t)] = λ t E[J], Var S(t) = λ t E[J²]
        constexpr double mean(double t) const requires has_mean<L> { return detail::sd::fmul(detail::sd::fmul(lambda_, t), double(jump_.mean())); }
        constexpr double variance(double t) const requires has_mean<L> && has_variance<L>
        {
            double const m = jump_.mean();
            return detail::sd::fmul(detail::sd::fmul(lambda_, t), double(jump_.variance()) + detail::sd::fmul(m, m));
        }

        struct walker_type
        {
            compound_poisson const* p;
            poisson_process::walker_type times;
            detail::lazy_sum values{};
            template <class G>
            constexpr state_type first(G&) { return {0.0, 0.0}; }
            template <class G>
            constexpr state_type next(G& g)
            {
                double const t = times.next(g);
                if constexpr (detail::lazy_samplable<L, G>)
                    values.add(detail::lazy_sample(p->jump_, g), g);
                else
                    values.add_constant(detail::to_dyadic(double(wsl::random::sample(p->jump_, g))));
                return {t, values.rounded(g)};
            }
        };
        constexpr walker_type walker() const { return {this, poisson_process::walker_type{lambda_}}; }

    private:
        double lambda_;
        L jump_;
    };

    // -----------------------------------------------------------------------
    // Chaîne de Markov en temps continu, de générateur Q (lignes de somme nulle).
    // -----------------------------------------------------------------------
    class continuous_markov_chain
    {
    public:
        using state_type = std::pair<double, std::size_t>;                 // (instant du saut, état)

        struct parameters
        {
            square_matrix<double> generator;
            std::size_t start = 0;
        };

        constexpr explicit continuous_markov_chain(parameters p) : q_(std::move(p.generator)), start_(p.start)
        {
            std::size_t const k = q_.size();
            if (k < 1 || start_ >= k)
                throw std::domain_error("wsl::random::continuous_markov_chain : au moins un état, départ valide");
            for (std::size_t i = 0; i < k; ++i)
            {
                std::vector<double> row(k, 0.0);
                detail::sd::double_double out{0.0, 0.0};
                for (std::size_t j = 0; j < k; ++j)
                {
                    if (j == i)
                        continue;
                    if (!(q_(i, j) >= 0) || !detail::is_finite(q_(i, j)))
                        throw std::domain_error("wsl::random::continuous_markov_chain : taux hors diagonale >= 0, finis");
                    row[j] = q_(i, j);
                    out = detail::sum(out, {q_(i, j), 0.0});
                }
                auto const d = detail::sum(out, {q_(i, i), 0.0});
                if (!(detail::magnitude(d.hi) <= detail::sd::fmul(out.hi + 0x1p-1000, 0x1p-40)))
                    throw std::domain_error("wsl::random::continuous_markov_chain : chaque ligne du générateur doit sommer à 0");
                rates_.push_back(-q_(i, i));
                if (out.hi > 0)
                    jumps_.push_back(categorical<>({.p = row}));
                else
                    jumps_.push_back(categorical<>({.p = std::vector<double>(k, 1.0)}));   // absorbant : jamais tiré
            }
        }

        constexpr std::size_t size() const noexcept { return q_.size(); }
        constexpr square_matrix<double> const& generator() const noexcept { return q_; }

        // Loi de X(t) : π₀ e^{Qt}, par mise à l'échelle et carrés (Taylor de degré 24).
        constexpr categorical<> at(double t) const
        {
            std::size_t const k = size();
            double norm = 0;
            for (std::size_t i = 0; i < k; ++i)
                norm = std::max(norm, detail::sd::fmul(detail::magnitude(q_(i, i)), 2.0));
            int squarings = 0;
            double scale = t;
            while (detail::sd::fmul(norm, scale) > 0.5)
            {
                scale = detail::sd::fmul(scale, 0.5);
                ++squarings;
            }
            square_matrix<double> a(k), term(k), e(k);
            for (std::size_t i = 0; i < k; ++i)
                for (std::size_t j = 0; j < k; ++j)
                {
                    a(i, j) = detail::sd::fmul(q_(i, j), scale);
                    term(i, j) = e(i, j) = i == j ? 1.0 : 0.0;
                }
            for (int n = 1; n <= 24; ++n)
            {
                term = detail::multiply(term, a);
                for (std::size_t i = 0; i < k; ++i)
                    for (std::size_t j = 0; j < k; ++j)
                    {
                        term(i, j) = term(i, j) / n;
                        e(i, j) = e(i, j) + term(i, j);
                    }
            }
            for (int s = 0; s < squarings; ++s)
                e = detail::multiply(e, e);
            std::vector<double> pi(k);
            for (std::size_t j = 0; j < k; ++j)
                pi[j] = e(start_, j) < 0 ? 0.0 : e(start_, j);
            return categorical<>({.p = pi});
        }

        struct walker_type
        {
            continuous_markov_chain const* c;
            detail::lazy_sum times;
            std::size_t x;
            bool absorbed = false;
            template <class G>
            constexpr state_type first(G&) { return {0.0, x}; }
            template <class G>
            constexpr state_type next(G& g)
            {
                double const q = c->rates_[x];
                if (absorbed || !(q > 0))
                {
                    absorbed = true;
                    return {std::numeric_limits<double>::infinity(), x};
                }
                times.add(detail::lazy_sample(exponential<double>({.rate = q}), g), g);
                x = static_cast<std::size_t>(wsl::random::sample(c->jumps_[x], g));
                return {times.rounded(g), x};
            }
        };
        constexpr walker_type walker() const { return {this, detail::lazy_sum{}, start_}; }

    private:
        square_matrix<double> q_;
        std::size_t start_;
        std::vector<double> rates_;
        std::vector<categorical<>> jumps_;
    };
}
