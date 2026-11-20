#pragma once

// Graphes et arbres aléatoires, sur les sommets {0, …, n − 1}.
//
//   wsl::random::erdos_renyi G({.vertices = 100, .p = 0.05});           // G(n, p)
//   wsl::random::uniform_graph H({.vertices = 100, .edges = 250});       // G(n, m)
//   wsl::random::barabasi_albert B({.vertices = 1000, .edges = 3});      // attachement préférentiel
//   wsl::random::random_regular_graph R({.vertices = 20, .degree = 3});  // uniforme parmi les d-réguliers
//   wsl::random::uniform_tree T({.vertices = 50});                       // arbres de Cayley : n^{n−2}
//   wsl::random::uniform_binary_tree A({.vertices = 30});                   // arbres binaires : Catalan
//
// Valeurs : liste d'arêtes (u, v), u < v, triée (wsl::random::edge_list) ; arbre binaire : enfants
// {gauche, droit} de chaque nœud, numérotés en préordre (racine 0), wsl::random::binary_tree::none
// pour l'absence d'enfant (wsl::random::binary_tree).
//
// Tirages exacts :
// - G(n, p) : une épreuve de Bernoulli exacte par paire, dans l'ordre lexicographique
//   (O(n²)) ; masse p^m (1 − p)^{N − m}, N = n(n − 1)/2 ;
// - G(n, m) : m paires par l'algorithme de Floyd ; 1/C(N, m) ;
// - Barabási–Albert : K_{m+1} au départ, puis chaque sommet s'attache à m sommets
//   distincts tirés proportionnellement au degré (extrémité uniforme d'une arête,
//   doublons rejetés) ; pas de masse calculée ;
// - d-régulier : modèle de configuration (couplage uniforme des n d demi-arêtes) rejeté
//   tant qu'il a une boucle ou une arête double : uniforme sur les graphes simples
//   (essais e^{(d² − 1)/4} en moyenne : petits d) ;
// - Cayley : suite de Prüfer uniforme, décodée en O(n) ; n^{−(n−2)} ;
// - binaire : algorithme de Rémy (1985) sur les arbres binaires complets à n nœuds
//   internes, puis feuilles retirées ; 1/C_n, C_n = C(2n, n)/(n + 1).

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/permutations.hpp>
#include <wsl/random/subsets.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wsl::random
{
    using edge_list = std::vector<std::pair<std::size_t, std::size_t>>;

    struct binary_tree
    {
        static constexpr std::size_t none = static_cast<std::size_t>(-1);
        std::vector<std::array<std::size_t, 2>> children;                   // {gauche, droit}, préordre

        constexpr std::size_t size() const noexcept { return children.size(); }
        friend constexpr bool operator==(binary_tree const&, binary_tree const&) = default;
        friend constexpr auto operator<=>(binary_tree const&, binary_tree const&) = default;
    };

    namespace detail
    {
        // u < v, croissantes, sur n sommets.
        constexpr bool is_edge_list(edge_list const& e, std::size_t n) noexcept
        {
            for (std::size_t i = 0; i < e.size(); ++i)
                if (e[i].first >= e[i].second || e[i].second >= n || (i > 0 && !(e[i - 1] < e[i])))
                    return false;
            return true;
        }

        // k-ième paire (u, v), u < v, dans l'ordre lexicographique.
        constexpr std::pair<std::size_t, std::size_t> pair_of_index(std::uint64_t k, std::size_t n) noexcept
        {
            std::size_t u = 0;
            std::uint64_t row = n - 1;
            while (k >= row)
            {
                k -= row;
                ++u;
                --row;
            }
            return {u, u + 1 + static_cast<std::size_t>(k)};
        }

        constexpr std::uint64_t pair_count(std::size_t n) noexcept { return std::uint64_t(n) * (n - 1) / 2; }

        constexpr void normalize(edge_list& e)
        {
            for (auto& [u, v] : e)
                if (u > v)
                    std::swap(u, v);
            std::ranges::sort(e);
        }

        // Connexe et sans cycle (n − 1 arêtes), par union-find.
        constexpr bool is_tree(edge_list const& e, std::size_t n)
        {
            if (n == 0 || e.size() != n - 1 || !is_edge_list(e, n))
                return false;
            std::vector<std::size_t> parent(n);
            for (std::size_t i = 0; i < n; ++i)
                parent[i] = i;
            auto const find = [&](std::size_t x) {
                while (parent[x] != x)
                    x = parent[x] = parent[parent[x]];
                return x;
            };
            for (auto const& [u, v] : e)
            {
                std::size_t const a = find(u), b = find(v);
                if (a == b)
                    return false;
                parent[a] = b;
            }
            return true;
        }
    }

    class erdos_renyi : public law_interface
    {
    public:
        using value_type = edge_list;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t vertices;
            double p;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.vertices = vertices(), .p = p()}; }

        constexpr explicit erdos_renyi(parameters q) : n_(q.vertices), p_(q.p)
        {
            if (!(p_ >= 0 && p_ <= 1))
                throw std::domain_error("wsl::random::erdos_renyi : il faut p dans [0, 1]");
        }

        constexpr std::size_t vertices() const noexcept { return n_; }
        constexpr double p() const noexcept { return p_; }
        constexpr double mean_edges() const noexcept { return detail::sd::fmul(double(detail::pair_count(n_)), p_); }

        constexpr double log_pmf(edge_list const& e) const noexcept
        {
            auto const& G = *this;
            auto const l = G.log_mass(e);
            return l.hi + l.lo;
        }
        constexpr double pmf(edge_list const& e) const noexcept { auto const& G = *this; return detail::exp_of_dd(G.log_mass(e)); }

        template <std::uniform_random_bit_generator G>
        constexpr edge_list sample(G& g) const
        {
            auto const& R = *this;
            edge_list e;
            for (std::size_t u = 0; u < R.n_; ++u)
                for (std::size_t v = u + 1; v < R.n_; ++v)
                    if (bernoulli_trial(R.p_, g))
                        e.emplace_back(u, v);
            return e;
        }

        friend constexpr bool operator==(erdos_renyi const&, erdos_renyi const&) = default;

    private:
        constexpr detail::sd::double_double log_mass(edge_list const& e) const noexcept
        {
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            if (!detail::is_edge_list(e, n_))
                return none;
            double const m = double(e.size()), rest = double(detail::pair_count(n_)) - m;
            detail::sd::double_double l{0.0, 0.0};
            if (m > 0)
            {
                if (p_ == 0)
                    return none;
                l = detail::mul({m, 0.0}, detail::log_precise(p_));
            }
            if (rest > 0)
            {
                if (p_ == 1)
                    return none;
                l = detail::sum(l, detail::mul({rest, 0.0}, detail::log1p_dd(-p_)));
            }
            return l;
        }

        std::size_t n_;
        double p_;
    };

    class uniform_graph : public law_interface
    {
    public:
        using value_type = edge_list;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t vertices, edges;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.vertices = n_, .edges = m_}; }

        constexpr explicit uniform_graph(parameters q) : n_(q.vertices), m_(q.edges)
        {
            if (m_ > detail::pair_count(n_))
                throw std::domain_error("wsl::random::uniform_graph : plus d'arêtes que de paires");
            log_count_ = detail::log_binomial(detail::pair_count(n_), m_);
        }

        constexpr std::size_t vertices() const noexcept { return n_; }
        constexpr std::size_t edges() const noexcept { return m_; }
        constexpr double entropy() const noexcept { return log_count_.hi + log_count_.lo; }

        constexpr double log_pmf(edge_list const& e) const noexcept
        {
            auto const& G = *this;
            return e.size() == G.m_ && detail::is_edge_list(e, G.n_) ? -G.entropy() : -std::numeric_limits<double>::infinity();
        }
        constexpr double pmf(edge_list const& e) const noexcept
        {
            auto const& G = *this;
            return e.size() == G.m_ && detail::is_edge_list(e, G.n_) ? detail::exp_of_dd({-G.log_count_.hi, -G.log_count_.lo}) : 0.0;
        }

        // Floyd sur les indices des paires.
        template <std::uniform_random_bit_generator G>
        constexpr edge_list sample(G& g) const
        {
            auto const& R = *this;
            std::uint64_t const total = detail::pair_count(R.n_);
            std::vector<std::uint64_t> chosen;
            chosen.reserve(R.m_);
            for (std::uint64_t j = total - R.m_; j < total; ++j)
            {
                std::uint64_t const t = uniform_below(j + 1, g);
                auto const at = std::ranges::lower_bound(chosen, t);
                if (at != chosen.end() && *at == t)
                    chosen.insert(std::ranges::lower_bound(chosen, j), j);
                else
                    chosen.insert(at, t);
            }
            edge_list e;
            e.reserve(R.m_);
            for (std::uint64_t k : chosen)
                e.push_back(detail::pair_of_index(k, R.n_));
            return e;
        }

        friend constexpr bool operator==(uniform_graph const& a, uniform_graph const& b) noexcept { return a.n_ == b.n_ && a.m_ == b.m_; }

    private:
        std::size_t n_, m_;
        detail::sd::double_double log_count_{};
    };

    class barabasi_albert : public law_interface
    {
    public:
        using value_type = edge_list;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t vertices, edges;                                    // m arêtes par nouveau sommet
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.vertices = n_, .edges = m_}; }

        constexpr explicit barabasi_albert(parameters q) : n_(q.vertices), m_(q.edges)
        {
            if (m_ < 1 || n_ < m_ + 1)
                throw std::domain_error("wsl::random::barabasi_albert : il faut edges >= 1 et vertices >= edges + 1");
        }

        constexpr std::size_t vertices() const noexcept { return n_; }
        // C(m + 1, 2) + (n − m − 1) m
        constexpr std::size_t edge_count() const noexcept { return (m_ + 1) * m_ / 2 + (n_ - m_ - 1) * m_; }

        template <std::uniform_random_bit_generator G>
        constexpr edge_list sample(G& g) const
        {
            auto const& B = *this;
            edge_list e;
            std::vector<std::size_t> ends;                                  // chaque sommet, autant de fois que son degré
            for (std::size_t u = 0; u <= B.m_; ++u)
                for (std::size_t v = u + 1; v <= B.m_; ++v)
                {
                    e.emplace_back(u, v);
                    ends.push_back(u);
                    ends.push_back(v);
                }
            std::vector<std::size_t> targets;
            for (std::size_t v = B.m_ + 1; v < B.n_; ++v)
            {
                targets.clear();
                std::size_t const available = ends.size();                 // degrés avant l'arrivée de v
                while (targets.size() < B.m_)
                {
                    std::size_t const t = ends[static_cast<std::size_t>(uniform_below(available, g))];
                    if (std::ranges::find(targets, t) == targets.end())
                        targets.push_back(t);
                }
                for (std::size_t t : targets)
                {
                    e.emplace_back(t, v);
                    ends.push_back(t);
                    ends.push_back(v);
                }
            }
            detail::normalize(e);
            return e;
        }

        friend constexpr bool operator==(barabasi_albert const&, barabasi_albert const&) = default;

    private:
        std::size_t n_, m_;
    };

    class random_regular_graph : public law_interface
    {
    public:
        using value_type = edge_list;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t vertices, degree;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.vertices = vertices(), .degree = degree()}; }

        constexpr explicit random_regular_graph(parameters q) : n_(q.vertices), d_(q.degree)
        {
            if (d_ >= n_ || (n_ * d_) % 2 != 0)
                throw std::domain_error("wsl::random::random_regular_graph : il faut degree < vertices et vertices · degree pair");
        }

        constexpr std::size_t vertices() const noexcept { return n_; }
        constexpr std::size_t degree() const noexcept { return d_; }

        template <std::uniform_random_bit_generator G>
        constexpr edge_list sample(G& g) const
        {
            auto const& R = *this;
            std::size_t const points = R.n_ * R.d_;
            std::vector<std::size_t> owner(points);
            for (;;)
            {
                for (std::size_t i = 0; i < points; ++i)
                    owner[i] = i / R.d_;
                shuffle(owner, g);                                          // couplage uniforme : paires consécutives
                edge_list e;
                e.reserve(points / 2);
                bool simple = true;
                for (std::size_t i = 0; i < points && simple; i += 2)
                {
                    if (owner[i] == owner[i + 1])
                        simple = false;
                    else
                        e.emplace_back(owner[i], owner[i + 1]);
                }
                if (!simple)
                    continue;
                detail::normalize(e);
                if (std::ranges::adjacent_find(e) == e.end())
                    return e;
            }
        }

        friend constexpr bool operator==(random_regular_graph const&, random_regular_graph const&) = default;

    private:
        std::size_t n_, d_;
    };

    class uniform_tree : public law_interface
    {
    public:
        using value_type = edge_list;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t vertices;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.vertices = vertices()}; }

        constexpr explicit uniform_tree(parameters q) : n_(q.vertices)
        {
            if (n_ < 1)
                throw std::domain_error("wsl::random::uniform_tree : il faut au moins un sommet");
        }

        constexpr std::size_t vertices() const noexcept { return n_; }
        // (n − 2) log n
        constexpr double entropy() const noexcept
        {
            if (n_ <= 2)
                return 0;
            auto const l = detail::mul({double(n_ - 2), 0.0}, detail::log_precise(double(n_)));
            return l.hi + l.lo;
        }

        constexpr double log_pmf(edge_list const& e) const { auto const& T = *this; return detail::is_tree(e, T.n_) ? -T.entropy() : -std::numeric_limits<double>::infinity(); }
        constexpr double pmf(edge_list const& e) const
        {
            auto const& T = *this;
            if (!detail::is_tree(e, T.n_))
                return 0;
            if (T.n_ <= 2)
                return 1;
            return detail::exp_of_dd(detail::mul({-double(T.n_ - 2), 0.0}, detail::log_precise(double(T.n_))));
        }

        // Suite de Prüfer uniforme (n − 2 lettres), décodage linéaire.
        template <std::uniform_random_bit_generator G>
        constexpr edge_list sample(G& g) const
        {
            auto const& T = *this;
            std::size_t const n = T.n_;
            edge_list e;
            if (n < 2)
                return e;
            std::vector<std::size_t> code(n - 2), degree(n, 1);
            for (auto& c : code)
            {
                c = static_cast<std::size_t>(uniform_below(n, g));
                ++degree[c];
            }
            std::size_t pointer = 0;
            while (degree[pointer] != 1)
                ++pointer;
            std::size_t leaf = pointer;
            for (std::size_t a : code)
            {
                e.emplace_back(leaf, a);
                if (--degree[a] == 1 && a < pointer)
                    leaf = a;
                else
                {
                    ++pointer;
                    while (degree[pointer] != 1)
                        ++pointer;
                    leaf = pointer;
                }
            }
            e.emplace_back(leaf, n - 1);
            detail::normalize(e);
            return e;
        }

        friend constexpr bool operator==(uniform_tree const&, uniform_tree const&) = default;

    private:
        std::size_t n_;
    };

    class uniform_binary_tree : public law_interface
    {
    public:
        using value_type = binary_tree;
        static constexpr exactness sampling = exactness::exact;

        struct parameters
        {
            std::size_t vertices;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.vertices = vertices()}; }

        constexpr explicit uniform_binary_tree(parameters q) : n_(q.vertices)
        {
            // log C_n = log C(2n, n) − log(n + 1)
            auto const l = detail::log_binomial(2 * std::uint64_t(n_), n_);
            auto const ln = detail::log_precise(double(n_) + 1);
            log_count_ = detail::sum(l, {-ln.hi, -ln.lo});
        }

        constexpr std::size_t vertices() const noexcept { return n_; }
        constexpr double entropy() const noexcept { return log_count_.hi + log_count_.lo; }

        constexpr double log_pmf(binary_tree const& t) const { auto const& T = *this; return T.valid(t) ? -T.entropy() : -std::numeric_limits<double>::infinity(); }
        constexpr double pmf(binary_tree const& t) const
        {
            auto const& T = *this;
            return T.valid(t) ? detail::exp_of_dd({-T.log_count_.hi, -T.log_count_.lo}) : 0.0;
        }

        // Rémy : arbre complet, chaque étape choisit un nœud uniforme parmi 2i − 1 et un côté.
        template <std::uniform_random_bit_generator G>
        constexpr binary_tree sample(G& g) const
        {
            auto const& T = *this;
            std::size_t const n = T.n_;
            constexpr std::size_t none = binary_tree::none;
            std::vector<std::array<std::size_t, 2>> child(2 * n + 1, {none, none});
            std::vector<std::size_t> parent(2 * n + 1, none);
            std::size_t root = 0, count = 1;                                // le nœud 0 : une feuille
            for (std::size_t i = 1; i <= n; ++i)
            {
                std::size_t const x = static_cast<std::size_t>(uniform_below(count, g));
                std::size_t const side = static_cast<std::size_t>(uniform_below(2, g));
                std::size_t const y = count, z = count + 1;                 // y interne, z feuille
                count += 2;
                std::size_t const p = parent[x];
                if (p == none)
                    root = y;
                else
                    child[p][child[p][0] == x ? 0 : 1] = y;
                parent[y] = p;
                child[y][side] = z;
                child[y][1 - side] = x;
                parent[x] = parent[z] = y;
            }
            // Nœuds internes en préordre ; les feuilles disparaissent.
            binary_tree t;
            if (n == 0)
                return t;
            t.children.reserve(n);
            std::vector<std::size_t> stack{root}, label(2 * n + 1, none);
            std::vector<std::size_t> order;
            while (!stack.empty())
            {
                std::size_t const v = stack.back();
                stack.pop_back();
                if (child[v][0] == none)
                    continue;                                               // feuille
                label[v] = order.size();
                order.push_back(v);
                stack.push_back(child[v][1]);
                stack.push_back(child[v][0]);
            }
            for (std::size_t v : order)
            {
                std::array<std::size_t, 2> c{none, none};
                for (int s = 0; s < 2; ++s)
                    if (child[child[v][s]][0] != none)
                        c[s] = label[child[v][s]];
                t.children.push_back(c);
            }
            return t;
        }

        friend constexpr bool operator==(uniform_binary_tree const& a, uniform_binary_tree const& b) noexcept { return a.n_ == b.n_; }

    private:
        // n nœuds numérotés en préordre depuis 0.
        constexpr bool valid(binary_tree const& t) const
        {
            if (t.size() != n_)
                return false;
            if (n_ == 0)
                return true;
            std::size_t next = 1;
            std::vector<std::size_t> stack{0};
            while (!stack.empty())
            {
                std::size_t const v = stack.back();
                stack.pop_back();
                auto const [l, r] = t.children[v];                          // préordre : le prochain dépilé porte le numéro suivant
                if (r != binary_tree::none)
                    stack.push_back(r);
                if (l != binary_tree::none)
                    stack.push_back(l);
                if (!stack.empty())
                {
                    if (stack.back() != next || next >= n_)
                        return false;
                    ++next;
                }
            }
            return next == n_;
        }

        std::size_t n_;
        detail::sd::double_double log_count_{};
    };
}
