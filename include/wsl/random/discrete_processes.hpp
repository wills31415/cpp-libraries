#pragma once

// Processus en temps discret (voir design/random.md §14) : marches aléatoires à pas
// entiers, chaînes de Markov, urnes de Pólya, restaurant chinois (Pitman–Yor).
//
//   wsl::random::random_walk W(wsl::random::rademacher({.p = 0.5}));                       // marche simple
//   wsl::random::markov_chain C({.transition = {{0.9, 0.1}, {0.5, 0.5}}, .start = 0});
//   wsl::random::polya_urn U({.balls = {1, 2}, .added = 1});
//   wsl::random::chinese_restaurant R({.theta = 1.5});
//   for (auto x : wsl::random::trajectory(C, g) | std::views::take(100)) …
//   C.at(5); C.stationary(); wsl::random::sample(C.steps(100), g); wsl::random::pmf(C.steps(3), {0, 0, 1, 1});
//   U.at(10);  R.at(100);                                                  // lois marginales
//
// Tous les tirages sont exacts : pas tirés exactement (sommes entières), transitions
// catégorielles exactes (alias), urnes et restaurant par poids entiers exacts (les
// poids flottants sont des entiers à une même puissance de 2 près ; grands entiers).
// Lois des trajectoires finies (steps(n)) : masse produit, en logarithme double-double.
//
// Lois marginales : chaîne de Markov, catégorielle π₀ Pⁿ (produits en double-double,
// carrés successifs) ; urne de Pólya, Dirichlet-multinomiale (bêta-binomiale pour deux
// couleurs) ; restaurant chinois, partition d'Ewens–Pitman (fonction de probabilité de
// partition échangeable, EPPF).

#include <wsl/random/categorical.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/big.hpp>
#include <wsl/random/detail/choice.hpp>
#include <wsl/random/detail/lazy_sum.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/multinomial.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/permutations.hpp>
#include <wsl/random/process.hpp>
#include <wsl/random/convolution.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wsl::random
{
    // ±1 : P(+1) = p, P(−1) = 1 − p (pas de la marche simple).
    class rademacher : public law_interface
    {
    public:
        using value_type = std::int64_t;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            double p = 0.5;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.p = p()}; }

        constexpr rademacher() noexcept : p_(0.5) {}
        constexpr explicit rademacher(parameters q) : p_(q.p)
        {
            if (!(p_ >= 0 && p_ <= 1))
                throw std::domain_error("wsl::random::rademacher : il faut p dans [0, 1]");
        }

        constexpr double p() const noexcept { return p_; }
        constexpr interval<std::int64_t> support() const noexcept { return {-1, 1}; }
        constexpr double mean() const noexcept { return detail::sd::fmul(2.0, p_) - 1; }
        constexpr double variance() const noexcept { return detail::sd::fmul(4.0, detail::sd::fmul(p_, 1 - p_)); }

        constexpr double pmf(std::int64_t k) const noexcept { auto const& X = *this; return k == 1 ? X.p_ : k == -1 ? 1 - X.p_ : 0.0; }
        constexpr double cdf(std::int64_t k) const noexcept { auto const& X = *this; return k < -1 ? 0.0 : k < 1 ? 1 - X.p_ : 1.0; }
        constexpr double ccdf(std::int64_t k) const noexcept { auto const& X = *this; return k < -1 ? 1.0 : k < 1 ? X.p_ : 0.0; }

        template <std::uniform_random_bit_generator G>
        constexpr std::int64_t sample(G& g) const
        {
            auto const& X = *this;
            return bernoulli_trial(X.p_, g) ? 1 : -1;
        }

        friend constexpr bool operator==(rademacher const&, rademacher const&) = default;

    private:
        double p_;
    };

    namespace detail
    {
        // Loi d'une trajectoire finie : les n + 1 premiers états du marcheur de P.
        template <class P>
        class finite_path : public law_interface
        {
        public:
            using value_type = std::vector<typename P::state_type>;
            static constexpr exactness sampling = exactness::exact;

            constexpr finite_path(P const& process, std::size_t n) : process_(process), n_(n) {}

            constexpr std::size_t steps() const noexcept { return n_; }

            constexpr double log_pmf(value_type const& path) const requires requires(P const& p, value_type const& v) { p.log_path_mass(v); }
            {
                auto const& F = *this;
                if (path.size() != F.n_ + 1)
                    return -std::numeric_limits<double>::infinity();
                auto const l = F.process_.log_path_mass(path);
                return l.hi + l.lo;
            }
            constexpr double pmf(value_type const& path) const requires requires(P const& p, value_type const& v) { p.log_path_mass(v); }
            {
                auto const& F = *this;
                if (path.size() != F.n_ + 1)
                    return 0;
                return exp_of_dd(F.process_.log_path_mass(path));
            }

            template <std::uniform_random_bit_generator G>
            constexpr value_type sample(G& g) const
            {
                auto const& F = *this;
                auto w = F.process_.walker();
                value_type path;
                path.reserve(F.n_ + 1);
                path.push_back(w.first(g));
                for (std::size_t k = 0; k < F.n_; ++k)
                    path.push_back(w.next(g));
                return path;
            }

        private:
            P process_;
            std::size_t n_;
        };
    }

    // -----------------------------------------------------------------------
    // Marche aléatoire à pas entiers : X₀ = start, X_{k+1} = X_k + pas.
    // -----------------------------------------------------------------------
    template <probability_law L>
        requires std::integral<typename L::value_type> || std::floating_point<typename L::value_type>
    class random_walk
    {
    public:
        using state_type = typename L::value_type;
        using step_law = L;

        struct origin
        {
            state_type start = 0;
        };

        constexpr explicit random_walk(L step) : step_(std::move(step)), start_(0) {}
        constexpr random_walk(L step, origin o) : step_(std::move(step)), start_(o.start) {}

        constexpr L const& step() const noexcept { return step_; }
        constexpr state_type start() const noexcept { return start_; }
        // E[X_n] = start + n E[pas], Var X_n = n Var pas.
        constexpr double mean(std::size_t n) const noexcept requires has_mean<L> { return double(start_) + detail::sd::fmul(double(n), double(step_.mean())); }
        constexpr double variance(std::size_t n) const noexcept requires has_variance<L> { return detail::sd::fmul(double(n), double(step_.variance())); }

        // Pas entiers : somme entière (saturée). Pas réels : somme paresseuse exacte, arrondie
        // une fois (pas tirés paresseusement si la loi le permet, sinon pas arrondis sommés
        // exactement).
        struct walker_type
        {
            L step;
            state_type x;
            [[no_unique_address]] std::conditional_t<std::floating_point<state_type>, detail::lazy_sum, detail::sd::double_double> sum{};
            template <class G>
            constexpr state_type first(G&) { return x; }
            template <class G>
            constexpr state_type next(G& g)
            {
                if constexpr (std::integral<state_type>)
                    return x = detail::saturated_sum(x, static_cast<state_type>(wsl::random::sample(step, g)));
                else
                {
                    if constexpr (detail::lazy_samplable<L, G>)
                        sum.add(detail::lazy_sample(step, g), g);
                    else
                        sum.add_constant(detail::to_dyadic(double(wsl::random::sample(step, g))));
                    return static_cast<state_type>(sum.template rounded<double>(g, detail::to_dyadic(double(x)), {}, 0, detail::to_dyadic(1.0)));
                }
            }
        };
        constexpr walker_type walker() const { return {step_, start_}; }
        constexpr detail::finite_path<random_walk> steps(std::size_t n) const { return {*this, n}; }

        // Σ log P(pas = x_{k+1} − x_k), x₀ = start.
        constexpr detail::sd::double_double log_path_mass(std::vector<state_type> const& path) const requires has_pmf<L>
        {
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            if (path.empty() || path[0] != start_)
                return none;
            detail::sd::double_double l{0.0, 0.0};
            for (std::size_t k = 1; k < path.size(); ++k)
            {
                double const p = wsl::random::pmf(step_, static_cast<state_type>(path[k] - path[k - 1]));
                if (!(p > 0))
                    return none;
                l = detail::sum(l, detail::log_precise(p));
            }
            return l;
        }

    private:
        L step_;
        state_type start_;
    };

    // -----------------------------------------------------------------------
    // Chaîne de Markov homogène sur {0, …, K − 1}.
    // -----------------------------------------------------------------------
    class markov_chain
    {
    public:
        using state_type = std::size_t;

        struct from_state
        {
            square_matrix<double> transition;
            std::size_t start = 0;
        };
        struct from_law
        {
            square_matrix<double> transition;
            std::vector<double> initial;
        };

        constexpr explicit markov_chain(from_state p) : p_(std::move(p.transition)), initial_(point_mass(p.start, p_.extent(0)))
        {
            if (p.start >= p_.extent(0))
                throw std::domain_error("wsl::random::markov_chain : état initial hors de {0, …, K − 1}");
            build();
        }
        constexpr explicit markov_chain(from_law p) : p_(std::move(p.transition)), initial_(std::move(p.initial))
        {
            if (initial_.size() != p_.extent(0))
                throw std::domain_error("wsl::random::markov_chain : loi initiale de dimension différente");
            build();
        }

        constexpr std::size_t size() const noexcept { return p_.extent(0); }
        constexpr square_matrix<double> const& transition() const noexcept { return p_; }
        constexpr double probability(std::size_t i, std::size_t j) const noexcept { return wsl::random::pmf(rows_[i], j); }
        constexpr bool is_absorbing(std::size_t i) const noexcept { return probability(i, i) == 1; }

        // Loi de X_n : π₀ Pⁿ (carrés successifs de P, produits en double-double).
        constexpr categorical<> at(std::size_t n) const
        {
            std::vector<double> pi = normalized_initial();
            square_matrix<double> power = probabilities();
            for (std::size_t e = n; e > 0; e >>= 1)
            {
                if (e & 1)
                    pi = times(pi, power);
                if (e > 1)
                    power = detail::multiply(power, power);
            }
            return categorical<>({.p = pi});
        }
        // Loi stationnaire π = π P (système (Pᵀ − I) π = 0, Σ π = 1, pivot partiel) ;
        // std::domain_error si elle n'est pas unique.
        constexpr categorical<> stationary() const
        {
            std::size_t const k = size();
            auto const q = probabilities();
            std::vector<std::vector<double>> a(k, std::vector<double>(k + 1, 0.0));
            for (std::size_t i = 0; i < k; ++i)                             // ligne i : Σ_j (P_ji − δ_ij) π_j = 0
                for (std::size_t j = 0; j < k; ++j)
                    a[i][j] = q[j, i] - (i == j ? 1.0 : 0.0);
            for (std::size_t j = 0; j <= k; ++j)                            // dernière ligne : Σ π = 1
                a[k - 1][j] = 1.0;
            for (std::size_t c = 0; c < k; ++c)
            {
                std::size_t pivot = c;
                for (std::size_t r = c + 1; r < k; ++r)
                    if (detail::magnitude(a[r][c]) > detail::magnitude(a[pivot][c]))
                        pivot = r;
                if (!(detail::magnitude(a[pivot][c]) > 0x1p-40))
                    throw std::domain_error("wsl::random::markov_chain::stationary : loi stationnaire non unique");
                std::swap(a[c], a[pivot]);
                for (std::size_t r = 0; r < k; ++r)
                    if (r != c && a[r][c] != 0)
                    {
                        double const f = a[r][c] / a[c][c];
                        for (std::size_t j = c; j <= k; ++j)
                            a[r][j] = a[r][j] - detail::sd::fmul(f, a[c][j]);
                    }
            }
            std::vector<double> pi(k);
            for (std::size_t i = 0; i < k; ++i)
            {
                double const v = a[i][k] / a[i][i];
                pi[i] = v < 0 ? 0.0 : v;                                    // arrondis : −0 ou −ε
            }
            return categorical<>({.p = pi});
        }

        struct walker_type
        {
            markov_chain const* chain;
            std::size_t x = 0;
            template <class G>
            constexpr std::size_t first(G& g) { return x = static_cast<std::size_t>(wsl::random::sample(chain->initial_law_, g)); }
            template <class G>
            constexpr std::size_t next(G& g) { return x = static_cast<std::size_t>(wsl::random::sample(chain->rows_[x], g)); }
        };
        constexpr walker_type walker() const { return {this}; }
        constexpr detail::finite_path<markov_chain> steps(std::size_t n) const { return {*this, n}; }

        // log π₀(x₀) + Σ log P(x_k, x_{k+1})
        constexpr detail::sd::double_double log_path_mass(std::vector<std::size_t> const& path) const
        {
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            if (path.empty() || path[0] >= size())
                return none;
            double const p0 = wsl::random::pmf(initial_law_, path[0]);
            if (!(p0 > 0))
                return none;
            auto l = detail::log_precise(p0);
            for (std::size_t k = 1; k < path.size(); ++k)
            {
                if (path[k] >= size())
                    return none;
                double const p = probability(path[k - 1], path[k]);
                if (!(p > 0))
                    return none;
                l = detail::sum(l, detail::log_precise(p));
            }
            return l;
        }

    private:
        static constexpr std::vector<double> point_mass(std::size_t i, std::size_t k)
        {
            std::vector<double> v(k, 0.0);
            if (i < k)
                v[i] = 1;
            return v;
        }

        constexpr void build()
        {
            std::size_t const k = p_.extent(0);
            if (k < 1)
                throw std::domain_error("wsl::random::markov_chain : au moins un état");
            for (std::size_t i = 0; i < k; ++i)
            {
                std::vector<double> row(k);
                detail::sd::double_double s{0.0, 0.0};
                for (std::size_t j = 0; j < k; ++j)
                {
                    row[j] = p_[i, j];
                    if (!(row[j] >= 0) || !detail::is_finite(row[j]))
                        throw std::domain_error("wsl::random::markov_chain : probabilités de transition >= 0, finies");
                    s = detail::sum(s, {row[j], 0.0});
                }
                auto const d = detail::sum(s, {-1.0, 0.0});
                if (!(detail::magnitude(d.hi) <= 0x1p-40))
                    throw std::domain_error("wsl::random::markov_chain : chaque ligne de transition doit sommer à 1");
                rows_.push_back(categorical<>({.p = row}));
            }
            initial_law_ = categorical<>({.p = initial_});
        }

        // P normalisée ligne par ligne (probabilités exactes arrondies).
        constexpr square_matrix<double> probabilities() const
        {
            square_matrix<double> q(size());
            for (std::size_t i = 0; i < size(); ++i)
                for (std::size_t j = 0; j < size(); ++j)
                    q[i, j] = probability(i, j);
            return q;
        }
        constexpr std::vector<double> normalized_initial() const
        {
            std::vector<double> v(size());
            for (std::size_t i = 0; i < size(); ++i)
                v[i] = wsl::random::pmf(initial_law_, i);
            return v;
        }
        static constexpr std::vector<double> times(std::vector<double> const& v, square_matrix<double> const& m)
        {
            std::vector<double> r(v.size());
            for (std::size_t j = 0; j < v.size(); ++j)
            {
                auto const t = detail::dot(v.size(), [&](std::size_t i) { return std::pair{v[i], m[i, j]}; });
                r[j] = t.hi + t.lo;
            }
            return r;
        }

        square_matrix<double> p_;
        std::vector<double> initial_;
        std::vector<categorical<>> rows_;
        categorical<> initial_law_{{.p = {1.0}}};
    };

    // -----------------------------------------------------------------------
    // Loi de Dirichlet-multinomiale (n tirages d'une urne de Pólya) ; bêta-binomiale
    // pour K = 2.
    // -----------------------------------------------------------------------
    class dirichlet_multinomial : public law_interface
    {
    public:
        using value_type = std::vector<std::int64_t>;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::int64_t trials;
            std::vector<double> alpha;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.trials = trials(), .alpha = alpha()}; }

        constexpr explicit dirichlet_multinomial(parameters q) : n_(q.trials), alpha_(std::move(q.alpha))
        {
            if (n_ < 0 || alpha_.size() < 2)
                throw std::domain_error("wsl::random::dirichlet_multinomial : il faut trials >= 0 et au moins deux catégories");
            for (double const a : alpha_)
            {
                if (!(a > 0) || !detail::is_finite(a))
                    throw std::domain_error("wsl::random::dirichlet_multinomial : il faut alpha > 0, fini");
                total_ = detail::sum(total_, {a, 0.0});
            }
        }

        constexpr std::int64_t trials() const noexcept { return n_; }
        constexpr std::vector<double> const& alpha() const noexcept { return alpha_; }
        constexpr std::size_t dimension() const noexcept { return alpha_.size(); }
        // n αᵢ/A
        constexpr std::vector<double> mean() const noexcept
        {
            std::vector<double> m(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
            {
                auto const q = detail::sd::quotient_dd(detail::product(double(n_), alpha_[i]), total_);
                m[i] = q.hi + q.lo;
            }
            return m;
        }
        // n pᵢ (1 − pᵢ) (n + A)/(1 + A)
        constexpr std::vector<double> variance() const noexcept
        {
            namespace fp = detail::sd;
            double const a = total_.hi + total_.lo, n = double(n_);
            std::vector<double> v(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
            {
                double const p = alpha_[i] / a;
                v[i] = fp::fmul(fp::fmul(fp::fmul(n, p), 1 - p), (n + a) / (1 + a));
            }
            return v;
        }

        // log n! − Σ log kᵢ! + lgamma(A) − lgamma(n + A) + Σ (lgamma(kᵢ + αᵢ) − lgamma(αᵢ))
        constexpr double log_pmf(value_type const& k) const noexcept
        {
            auto const& D = *this;
            auto const l = D.log_mass(k);
            return l.hi + l.lo;
        }
        constexpr double pmf(value_type const& k) const noexcept { auto const& D = *this; return detail::exp_of_dd(D.log_mass(k)); }

        // Urne de Pólya : poids αᵢ + kᵢ exacts (entiers à une puissance de 2 près).
        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& D = *this;
            value_type k(D.dimension(), 0);
            int e = 0;
            for (double const a : D.alpha_)
            {
                int const ea = detail::dyadic_of(a).e;
                e = ea < e ? ea : e;
            }
            detail::growing_weights w;
            for (double const a : D.alpha_)
                w.push(detail::integer_at(a, e));
            detail::big const one = detail::big(1).shifted_left(-e);
            for (std::int64_t t = 0; t < D.n_; ++t)
            {
                std::size_t const i = w.pick(g);
                ++k[i];
                w.add(i, one);
            }
            return k;
        }

        friend constexpr bool operator==(dirichlet_multinomial const& a, dirichlet_multinomial const& b) noexcept { return a.n_ == b.n_ && a.alpha_ == b.alpha_; }

    private:
        constexpr detail::sd::double_double log_mass(value_type const& k) const noexcept
        {
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            if (k.size() != dimension())
                return none;
            std::int64_t s = 0;
            for (std::int64_t v : k)
            {
                if (v < 0 || v > n_ - s)
                    return none;
                s += v;
            }
            if (s != n_)
                return none;
            auto l = detail::log_factorial(std::uint64_t(n_));
            l = detail::sum(l, detail::log_gamma_any(total_.hi + total_.lo));
            auto const [t_hi, t_lo] = detail::sd::two_sum(total_.hi + total_.lo, double(n_));
            auto const ln = detail::log_gamma(t_hi, t_lo);
            l = detail::sum(l, {-ln.hi, -ln.lo});
            for (std::size_t i = 0; i < dimension(); ++i)
            {
                if (k[i] == 0)
                    continue;
                auto const f = detail::log_factorial(std::uint64_t(k[i]));
                l = detail::sum(l, {-f.hi, -f.lo});
                auto const [a_hi, a_lo] = detail::sd::two_sum(alpha_[i], double(k[i]));
                l = detail::sum(l, detail::log_gamma(a_hi, a_lo));
                auto const la = detail::log_gamma_any(alpha_[i]);
                l = detail::sum(l, {-la.hi, -la.lo});
            }
            return l;
        }

        std::int64_t n_;
        std::vector<double> alpha_;
        detail::sd::double_double total_{0.0, 0.0};
    };

    // -----------------------------------------------------------------------
    // Urne de Pólya : à chaque tirage, une boule de couleur i (probabilité nᵢ/Σ n), puis
    // `added` boules de cette couleur ajoutées. État : la couleur tirée.
    // -----------------------------------------------------------------------
    class polya_urn
    {
    public:
        using state_type = std::size_t;

        struct parameters
        {
            std::vector<double> balls;
            double added = 1;
        };

        constexpr explicit polya_urn(parameters q) : balls_(std::move(q.balls)), added_(q.added)
        {
            if (balls_.size() < 2 || !(added_ >= 0) || !detail::is_finite(added_))
                throw std::domain_error("wsl::random::polya_urn : au moins deux couleurs, added >= 0 fini");
            double total = 0;
            for (double const b : balls_)
            {
                if (!(b >= 0) || !detail::is_finite(b))
                    throw std::domain_error("wsl::random::polya_urn : nombres de boules >= 0, finis");
                total += b;
            }
            if (!(total > 0))
                throw std::domain_error("wsl::random::polya_urn : urne vide");
            exponent_ = added_ > 0 ? detail::dyadic_of(added_).e : 0;
            for (double const b : balls_)
                if (b > 0)
                    exponent_ = std::min(exponent_, detail::dyadic_of(b).e);
        }

        constexpr std::vector<double> const& balls() const noexcept { return balls_; }
        constexpr double added() const noexcept { return added_; }

        // Couleurs comptées après n tirages : Dirichlet-multinomiale(n, balls/added)
        // (added > 0 ; les paramètres sont les rapports, arrondis une fois).
        constexpr dirichlet_multinomial at(std::size_t n) const
        {
            if (!(added_ > 0))
                throw std::domain_error("wsl::random::polya_urn::at : added = 0 (tirages indépendants : multinomiale)");
            std::vector<double> alpha;
            for (double const b : balls_)
                alpha.push_back(b / added_);
            return dirichlet_multinomial({.trials = static_cast<std::int64_t>(n), .alpha = alpha});
        }

        struct walker_type
        {
            detail::growing_weights w;
            detail::big step;
            template <class G>
            constexpr std::size_t first(G& g) { return draw(g); }
            template <class G>
            constexpr std::size_t next(G& g) { return draw(g); }
            template <class G>
            constexpr std::size_t draw(G& g)
            {
                std::size_t const i = w.pick(g);
                w.add(i, step);
                return i;
            }
        };
        constexpr walker_type walker() const
        {
            walker_type r;
            for (double const b : balls_)
                r.w.push(detail::integer_at(b, exponent_));
            r.step = detail::integer_at(added_, exponent_);
            return r;
        }
        constexpr detail::finite_path<polya_urn> steps(std::size_t n) const { return {*this, n}; }

        // Π_k (b_{c_k} + a m_k)/(B + a k), m_k le nombre de c_k déjà tirés.
        constexpr detail::sd::double_double log_path_mass(std::vector<std::size_t> const& colors) const
        {
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            std::vector<double> seen(balls_.size(), 0.0);
            detail::sd::double_double total{0.0, 0.0};
            for (double const b : balls_)
                total = detail::sum(total, {b, 0.0});
            detail::sd::double_double l{0.0, 0.0};
            for (std::size_t k = 0; k < colors.size(); ++k)
            {
                std::size_t const c = colors[k];
                if (c >= balls_.size())
                    return none;
                auto const num = detail::sum({balls_[c], 0.0}, detail::product(added_, seen[c]));
                auto const den = detail::sum(total, detail::product(added_, double(k)));
                if (!(num.hi > 0))
                    return none;
                l = detail::sum(l, detail::log(detail::normalized(num.hi, num.lo)));
                auto const ld = detail::log(detail::normalized(den.hi, den.lo));
                l = detail::sum(l, {-ld.hi, -ld.lo});
                seen[c] += 1;
            }
            return l;
        }

    private:
        std::vector<double> balls_;
        double added_;
        int exponent_ = 0;
    };

    // -----------------------------------------------------------------------
    // Partition d'Ewens–Pitman de {0, …, n − 1} (mot à croissance restreinte), et
    // restaurant chinois (Pitman–Yor) : le client n s'assoit à la table j avec
    // probabilité (n_j − d)/(n + θ), à une nouvelle table avec (θ + K d)/(n + θ).
    // -----------------------------------------------------------------------
    namespace detail
    {
        // Tirage exact du client suivant : d = 0, un client précédent uniforme (O(1)) ;
        // sinon poids entiers (grands entiers), O(K).
        struct restaurant_state
        {
            double theta, d;
            std::vector<std::size_t> tables;                                // effectifs
            std::vector<std::size_t> seats;                                 // table de chaque client (d = 0)
            std::size_t customers = 0;
            // d > 0 : poids entiers, à l'échelle 2^{−e}, des tables (n_j − d) puis de la
            // nouvelle table (θ + K d), tenus à jour.
            growing_weights w{};
            big unit{}, dd{}, th{};

            template <class G>
            constexpr std::size_t seat(G& g)
            {
                std::size_t table = 0;
                if (customers == 0)
                    table = 0;
                else if (d == 0)                                            // nouvelle table : θ/(θ + n) ; sinon un client uniforme
                    table = choose(g, theta, double(customers)) ? tables.size() : seats[static_cast<std::size_t>(uniform_below(customers, g))];
                else
                {
                    std::size_t const c = w.pick(g);                        // dernier indice : nouvelle table
                    table = c + 1 == w.size() ? tables.size() : c;
                }
                if (d > 0)
                    record(table);
                if (table == tables.size())
                    tables.push_back(0);
                ++tables[table];
                ++customers;
                if (d == 0)
                    seats.push_back(table);
                return table;
            }

        private:
            // Avant la mise à jour des effectifs : un client de plus à `table`.
            constexpr void record(std::size_t table)
            {
                if (customers == 0)
                {
                    int e = std::min(dyadic_of(d).e, 0);
                    if (theta != 0)
                        e = std::min(e, dyadic_of(magnitude(theta)).e);
                    unit = big(1).shifted_left(-e);
                    dd = integer_at(d, e);
                    th = integer_at(magnitude(theta), e);
                }
                if (table < tables.size())
                {
                    w.add(table, unit);                                     // n_j − d gagne 1
                    return;
                }
                growing_weights r;                                          // une table s'ouvre : K + 1 tables
                for (std::size_t t : tables)
                    r.push(big(t) * unit - dd);
                r.push(unit - dd);
                big const kd = big(tables.size() + 1) * dd;
                r.push(theta >= 0 ? kd + th : kd - th);                     // θ + (K + 1) d > 0
                w = r;
            }
        };
    }

    class ewens_partition : public law_interface
    {
    public:
        using value_type = std::vector<std::size_t>;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t size;
            double theta;
            double discount = 0;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.size = size(), .theta = theta(), .discount = discount()}; }

        constexpr explicit ewens_partition(parameters q) : n_(q.size), theta_(q.theta), d_(q.discount)
        {
            if (!(d_ >= 0 && d_ < 1) || !(theta_ > -d_) || !detail::is_finite(theta_))
                throw std::domain_error("wsl::random::ewens_partition : il faut 0 <= discount < 1 et theta > −discount");
        }

        constexpr std::size_t size() const noexcept { return n_; }
        constexpr double theta() const noexcept { return theta_; }
        constexpr double discount() const noexcept { return d_; }
        // E[nombre de blocs] : Σ_{i<n} θ/(θ + i) (d = 0) ; (θ/d)((θ + d)_n/(θ)_n − 1) sinon.
        constexpr double mean_blocks() const noexcept
        {
            if (d_ == 0)
                return ewens({.size = n_, .theta = theta_}).mean_cycles();
            // Σ_{i<n} P(le client i ouvre une table) par la récurrence E[K_{i+1}] = E[K_i] + (θ + d E[K_i])/(θ + i)
            double k = 0;
            for (std::size_t i = 0; i < n_; ++i)
                k = k + (theta_ + detail::sd::fmul(d_, k)) / (theta_ + double(i));
            return k;
        }

        // EPPF : Π_{i=1}^{K−1} (θ + i d) / (θ + 1)_{n−1} · Π_j (1 − d)_{n_j − 1}
        constexpr double log_pmf(value_type const& a) const noexcept
        {
            auto const& P = *this;
            auto const l = P.log_mass(a);
            return l.hi + l.lo;
        }
        constexpr double pmf(value_type const& a) const noexcept { auto const& P = *this; return detail::exp_of_dd(P.log_mass(a)); }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& P = *this;
            detail::restaurant_state r{P.theta_, P.d_, {}, {}, 0, {}, {}, {}, {}};
            value_type a;
            a.reserve(P.n_);
            for (std::size_t i = 0; i < P.n_; ++i)
                a.push_back(r.seat(g));
            return a;
        }

        friend constexpr bool operator==(ewens_partition const&, ewens_partition const&) = default;

    private:
        constexpr detail::sd::double_double log_mass(value_type const& a) const noexcept
        {
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            if (a.size() != n_)
                return none;
            std::vector<std::size_t> sizes;
            for (std::size_t v : a)
            {
                if (v > sizes.size())
                    return none;
                if (v == sizes.size())
                    sizes.push_back(0);
                ++sizes[v];
            }
            detail::sd::double_double l{0.0, 0.0};
            auto const add_log = [&](double x_hi, double x_lo, bool minus) {
                auto const t = detail::log(detail::normalized(x_hi, x_lo));
                l = detail::sum(l, minus ? detail::sd::double_double{-t.hi, -t.lo} : t);
            };
            for (std::size_t i = 1; i < sizes.size(); ++i)                  // θ + i d
            {
                auto const [h, e] = detail::sd::two_sum(theta_, detail::sd::fmul(double(i), d_));
                add_log(h, e, false);
            }
            for (std::size_t i = 1; i < n_; ++i)                            // (θ + 1)_{n−1}
            {
                auto const [h, e] = detail::sd::two_sum(theta_, double(i));
                add_log(h, e, true);
            }
            for (std::size_t s : sizes)                                     // (1 − d)_{s − 1}
                for (std::size_t m = 1; m < s; ++m)
                {
                    auto const [h, e] = detail::sd::two_sum(double(m), -d_);
                    add_log(h, e, false);
                }
            return l;
        }

        std::size_t n_;
        double theta_, d_;
    };

    class chinese_restaurant
    {
    public:
        using state_type = std::size_t;                                     // table du client

        struct parameters
        {
            double theta;
            double discount = 0;
        };

        constexpr explicit chinese_restaurant(parameters q) : theta_(q.theta), d_(q.discount)
        {
            if (!(d_ >= 0 && d_ < 1) || !(theta_ > -d_) || !detail::is_finite(theta_))
                throw std::domain_error("wsl::random::chinese_restaurant : il faut 0 <= discount < 1 et theta > −discount");
        }

        constexpr double theta() const noexcept { return theta_; }
        constexpr double discount() const noexcept { return d_; }
        // Tables des n premiers clients : la partition d'Ewens–Pitman.
        constexpr ewens_partition at(std::size_t n) const { return ewens_partition({.size = n, .theta = theta_, .discount = d_}); }
        constexpr ewens_partition steps(std::size_t n) const { return at(n + 1); }   // n + 1 clients

        struct walker_type
        {
            detail::restaurant_state r;
            template <class G>
            constexpr std::size_t first(G& g) { return r.seat(g); }
            template <class G>
            constexpr std::size_t next(G& g) { return r.seat(g); }
        };
        constexpr walker_type walker() const { return {detail::restaurant_state{theta_, d_, {}, {}, 0, {}, {}, {}, {}}}; }

    private:
        double theta_, d_;
    };
}
