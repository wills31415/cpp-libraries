#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Loi multinomiale : n épreuves indépendantes à K issues de probabilités p₁, …, p_K
// (K >= 2, Σ pᵢ = 1 à 2^-40 près) ; Nᵢ le nombre d'issues i.
//
//   wsl::probability::multinomial<std::int64_t, 3> M({.trials = 10, .p = {0.2, 0.3, 0.5}});
//   M.pmf({2, 3, 5}); M.mean(); M.covariance(); M.marginal(2); M(g);
//
// Masse n!/Π kᵢ! Π pᵢ^{kᵢ} en logarithme double-double (lgamma, log pᵢ), un arrondi ;
// moyenne n pᵢ, covariance n pᵢ (δᵢⱼ − pⱼ) ; marginale i : binomiale(n, pᵢ).
//
// Tirage (approché comme la binomiale, bit-exact) : N₁ ~ B(n, p₁), puis
// Nᵢ ~ B(n − N₁ − … − Nᵢ₋₁, pᵢ / (pᵢ + … + p_K)), les restes en double-double ; K − 1
// tirages binomiaux, dans l'ordre.

#include <wsl/probability/laws/binomial.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/detail/loader.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/multivariate.hpp>
#include <wsl/special.hpp>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace wsl::inline WSL_DEFAULTS::probability::laws
{
    template <law_integer I = std::int64_t, std::size_t N = dynamic>
    class multinomial : public law_interface
    {
    public:
        using value_type = vector_of<I, N>;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            I trials;
            vector_of<double, N> p;
        };
        using param_type = parameters;
        constexpr param_type param() const { return {.trials = trials(), .p = p()}; }

        constexpr multinomial(parameters q) : n_(q.trials), p_(q.p)
        {
            if (n_ < 0)
                throw std::domain_error("wsl::probability::multinomial : il faut trials >= 0");
            if (p_.size() < 2)
                throw std::domain_error("wsl::probability::multinomial : il faut au moins deux issues");
            detail::sd::double_double s{0.0, 0.0};
            for (double const v : p_)
            {
                if (!(v >= 0) || !(v <= 1))
                    throw std::domain_error("wsl::probability::multinomial : il faut p dans [0, 1]");
                s = detail::sum(s, {v, 0.0});
            }
            auto const d = detail::sum(s, {-1.0, 0.0});
            if (!((d.hi < 0 ? -d.hi : d.hi) <= 0x1p-40))
                throw std::domain_error("wsl::probability::multinomial : la somme des p doit valoir 1");
            // Probabilités conditionnelles pᵢ / (pᵢ + … + p_K), restes en double-double.
            std::size_t const k = p_.size();
            conditional_.assign(k, 1.0);
            detail::sd::double_double rest{0.0, 0.0};
            for (std::size_t i = k; i-- > 0;)
            {
                rest = detail::sum(rest, {p_[i], 0.0});
                double const r = rest.hi + rest.lo;
                double const c = r > 0 ? p_[i] / r : 0.0;
                conditional_[i] = c > 1 ? 1.0 : c;
            }
        }

        constexpr I trials() const noexcept { return n_; }
        constexpr vector_of<double, N> const& p() const noexcept { return p_; }
        constexpr std::size_t dimension() const noexcept { return p_.size(); }
        constexpr interval<value_type> support() const noexcept
        {
            return {detail::make_vector<I, N>(dimension(), I(0)), detail::make_vector<I, N>(dimension(), n_)};
        }

        constexpr vector_of<double, N> mean() const noexcept
        {
            auto m = detail::make_vector<double, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                m[i] = detail::sd::fmul(double(n_), p_[i]);
            return m;
        }
        constexpr vector_of<double, N> variance() const noexcept
        {
            auto v = detail::make_vector<double, N>(dimension());
            for (std::size_t i = 0; i < dimension(); ++i)
                v[i] = covariance_entry(i, i);
            return v;
        }
        constexpr square_matrix<double, N> covariance() const
        {
            square_matrix<double, N> c = [&] {
                if constexpr (N == dynamic)
                    return square_matrix<double, N>(dimension());
                else
                    return square_matrix<double, N>();
            }();
            for (std::size_t i = 0; i < dimension(); ++i)
                for (std::size_t j = 0; j < dimension(); ++j)
                    c[i, j] = covariance_entry(i, j);
            return c;
        }
        constexpr binomial<I> marginal(std::size_t i) const { return binomial<I>({.trials = n_, .p = p_[i]}); }

        constexpr double log_pmf(value_type const& k) const noexcept
        {
            auto const& M = *this;
            auto const l = M.log_mass(k);
            return l.hi + l.lo;
        }
        constexpr double pmf(value_type const& k) const noexcept
        {
            auto const& M = *this;
            auto const l = M.log_mass(k);
            if (!(l.hi > -1500))
                return 0;
            return detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }

        template <std::uniform_random_bit_generator G>
        constexpr value_type sample(G& g) const
        {
            auto const& M = *this;
            value_type k = detail::make_vector<I, N>(M.dimension(), I(0));
            I rest = M.n_;
            std::size_t const last = M.dimension() - 1;
            for (std::size_t i = 0; i < last; ++i)
            {
                if (rest == 0)
                    break;
                double const c = M.conditional_[i];
                I const v = c >= 1 ? rest : c <= 0 ? I(0) : detail::native_sample(binomial<I>({.trials = rest, .p = c}), g);
                k[i] = v;
                rest = static_cast<I>(rest - v);
            }
            k[last] = static_cast<I>(k[last] + rest);
            return k;
        }

        friend constexpr bool operator==(multinomial const& a, multinomial const& b) noexcept { return a.n_ == b.n_ && a.p_ == b.p_; }

    private:
        // n pᵢ (δᵢⱼ − pⱼ)
        constexpr double covariance_entry(std::size_t i, std::size_t j) const noexcept
        {
            double const pi = p_[i], pj = p_[j];
            if (pi == 0 || pj == 0)
                return 0;
            if (n_ == 0)
                return 0;
            detail::scaled const npi = detail::normalized(double(n_), 0.0) * detail::normalized(pi, 0.0);
            if (i == j)
            {
                auto const [q, q_error] = detail::sd::two_sum(1.0, -pi);
                if (!(q > 0))
                    return 0;
                return detail::rounded(npi * detail::normalized(q, q_error));
            }
            return -detail::rounded(npi * detail::normalized(pj, 0.0));
        }

        // log n! − Σ log kᵢ! + Σ kᵢ log pᵢ ; −∞ hors du support.
        constexpr detail::sd::double_double log_mass(value_type const& k) const noexcept
        {
            constexpr detail::sd::double_double none{-std::numeric_limits<double>::infinity(), 0.0};
            if (k.size() != dimension())
                return none;
            I total = 0;
            detail::sd::double_double l = detail::log_gamma(double(n_) + 1);
            for (std::size_t i = 0; i < dimension(); ++i)
            {
                if (k[i] < 0 || k[i] > n_ - total)
                    return none;
                total = static_cast<I>(total + k[i]);
                if (k[i] == 0)
                    continue;
                if (p_[i] == 0)
                    return none;
                auto const lg = detail::log_gamma(double(k[i]) + 1);
                l = detail::sum(l, {-lg.hi, -lg.lo});
                l = detail::sum(l, detail::mul({double(k[i]), 0.0}, detail::log_precise(p_[i])));
            }
            if (total != n_)
                return none;
            return l;
        }

        I n_;
        vector_of<double, N> p_;
        std::vector<double> conditional_;
    };
}
