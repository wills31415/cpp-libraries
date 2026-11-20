#pragma once

// Lois uniformes sur la sphère unité Sⁿ⁻¹ = {x ∈ Rⁿ : |x| = 1} et sur la boule unité.
//
//   wsl::random::uniform_sphere<double, 3> S;                                  // dimension fixe
//   wsl::random::uniform_ball<double, wsl::random::dynamic> B({.dimension = 10});      // dimension dynamique
//   sample(S, g); pdf(B, x); S.surface_area(); B.volume();
//
// Sphère : densité par rapport à la mesure de surface, 1/|Sⁿ⁻¹| = Γ(n/2)/(2π^{n/2}),
// sur les x de norme 1 à n 2^-50 près (0 ailleurs) ; boule : 1/volume, volume
// π^{n/2}/Γ(n/2 + 1). Moyenne nulle, covariance I/n (sphère), I/(n + 2) (boule).
//
// Tirage (approché, bit-exact) : sphère, z/|z|, z de composantes normales centrées
// réduites tirées exactement, |z| en double-double, chaque composante arrondie une
// fois ; boule, la même direction multipliée par U^{1/n}, U uniforme sur (0, 1).

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/detail/loader.hpp>
#include <wsl/random/law.hpp>
#include <wsl/random/multivariate.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/special.hpp>

#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>

namespace wsl::random
{
    namespace detail
    {
        // log |Sⁿ⁻¹| = log 2 + (n/2) log π − lgamma(n/2), en double-double.
        constexpr sd::double_double log_sphere_area(std::size_t n) noexcept
        {
            double const h = sd::fmul(double(n), 0.5);
            auto l = sum({sd::ln2_hi, sd::ln2_lo}, mul({h, 0.0}, {sd::log_pi_hi, sd::log_pi_lo}));
            auto const lg = log_gamma(h);
            return sum(l, {-lg.hi, -lg.lo});
        }

        // Direction uniforme : z/|z|.
        template <class T, std::size_t N, class G>
        constexpr vector_of<T, N> uniform_direction(std::size_t n, G& g)
        {
            normal<double> const standard({.mean = 0.0, .sd = 1.0});
            auto z = make_vector<double, N>(n);
            for (;;)
            {
                for (auto& v : z)
                    v = sample(standard, g);
                auto const r2 = dot(n, [&](std::size_t k) { return std::pair{z[k], z[k]}; });
                if (r2.hi > 0)                                              // (probabilité nulle d'échouer)
                {
                    scaled const r = sqrt(normalized(r2.hi, r2.lo));
                    auto x = make_vector<T, N>(n);
                    for (std::size_t i = 0; i < n; ++i)
                    {
                        double const a = z[i] < 0 ? -z[i] : z[i];
                        double const q = a == 0 ? 0.0 : rounded(normalized(a, 0.0) / r);
                        x[i] = static_cast<T>(z[i] < 0 ? -q : q);
                    }
                    return x;
                }
            }
        }

        // |x|² − 1 en double-double.
        template <class V>
        constexpr double squared_norm_minus_one(V const& x) noexcept
        {
            auto const s = dot(x.size(), [&](std::size_t k) { return std::pair{double(x[k]), double(x[k])}; });
            auto const d = sum(s, {-1.0, 0.0});
            return d.hi + d.lo;
        }
    }

    template <law_real T = double, std::size_t N = dynamic>
    class uniform_sphere
    {
    public:
        using value_type = vector_of<T, N>;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            std::size_t dimension;
        };

        constexpr uniform_sphere() requires(N != dynamic) : n_(N) { check(); }
        constexpr uniform_sphere(parameters p) : n_(p.dimension)
        {
            if (N != dynamic && n_ != N)
                throw std::domain_error("wsl::random::uniform_sphere : dimension différente de N");
            check();
        }

        constexpr std::size_t dimension() const noexcept { return n_; }
        constexpr interval<value_type> support() const noexcept
        {
            return {detail::make_vector<T, N>(n_, T(-1)), detail::make_vector<T, N>(n_, T(1))};
        }
        constexpr value_type mean() const noexcept { return detail::make_vector<T, N>(n_, T(0)); }
        constexpr value_type variance() const noexcept { return detail::make_vector<T, N>(n_, static_cast<T>(1.0 / double(n_))); }
        constexpr double surface_area() const noexcept
        {
            auto const l = detail::log_sphere_area(n_);
            return l.hi > 1400 ? std::numeric_limits<double>::infinity() : detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }
        constexpr double entropy() const noexcept
        {
            auto const l = detail::log_sphere_area(n_);
            return l.hi + l.lo;
        }

        friend constexpr double pdf(uniform_sphere const& S, value_type const& x) noexcept
        {
            if (!S.on_sphere(x))
                return 0;
            auto const l = detail::log_sphere_area(S.n_);
            return detail::rounded(detail::exp_scaled(-l.hi, -l.lo));
        }
        friend constexpr double log_pdf(uniform_sphere const& S, value_type const& x) noexcept
        {
            if (!S.on_sphere(x))
                return -std::numeric_limits<double>::infinity();
            return -S.entropy();
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr value_type sample(uniform_sphere const& S, G& g)
        {
            return detail::uniform_direction<T, N>(S.n_, g);
        }

        friend constexpr bool operator==(uniform_sphere const&, uniform_sphere const&) = default;

    private:
        constexpr void check() const
        {
            if (n_ < 1)
                throw std::domain_error("wsl::random::uniform_sphere : il faut une dimension >= 1");
        }
        constexpr bool on_sphere(value_type const& x) const noexcept
        {
            if (x.size() != n_)
                return false;
            double const d = detail::squared_norm_minus_one(x);
            return (d < 0 ? -d : d) <= double(n_) * 0x1p-50;
        }

        std::size_t n_;
    };

    template <law_real T = double, std::size_t N = dynamic>
    class uniform_ball
    {
    public:
        using value_type = vector_of<T, N>;
        static constexpr exactness sampling = exactness::approximate;

        struct parameters
        {
            std::size_t dimension;
        };

        constexpr uniform_ball() requires(N != dynamic) : n_(N) { check(); }
        constexpr uniform_ball(parameters p) : n_(p.dimension)
        {
            if (N != dynamic && n_ != N)
                throw std::domain_error("wsl::random::uniform_ball : dimension différente de N");
            check();
        }

        constexpr std::size_t dimension() const noexcept { return n_; }
        constexpr interval<value_type> support() const noexcept
        {
            return {detail::make_vector<T, N>(n_, T(-1)), detail::make_vector<T, N>(n_, T(1))};
        }
        constexpr value_type mean() const noexcept { return detail::make_vector<T, N>(n_, T(0)); }
        constexpr value_type variance() const noexcept { return detail::make_vector<T, N>(n_, static_cast<T>(1.0 / (double(n_) + 2))); }
        // π^{n/2}/Γ(n/2 + 1) = |Sⁿ⁻¹|/n
        constexpr double volume() const noexcept
        {
            auto const l = log_volume();
            return l.hi > 1400 ? std::numeric_limits<double>::infinity() : l.hi < -1500 ? 0.0 : detail::rounded(detail::exp_scaled(l.hi, l.lo));
        }
        constexpr double entropy() const noexcept
        {
            auto const l = log_volume();
            return l.hi + l.lo;
        }

        friend constexpr double pdf(uniform_ball const& B, value_type const& x) noexcept
        {
            if (!B.inside(x))
                return 0;
            auto const l = B.log_volume();
            return l.hi < -1400 ? std::numeric_limits<double>::infinity() : detail::rounded(detail::exp_scaled(-l.hi, -l.lo));
        }
        friend constexpr double log_pdf(uniform_ball const& B, value_type const& x) noexcept
        {
            return B.inside(x) ? -B.entropy() : -std::numeric_limits<double>::infinity();
        }

        template <std::uniform_random_bit_generator G>
        friend constexpr value_type sample(uniform_ball const& B, G& g)
        {
            auto x = detail::uniform_direction<T, N>(B.n_, g);
            double const u = detail::uniform_positive(g);
            double const r = B.n_ == 1 ? u : wsl::special::exp(wsl::special::log(u) / double(B.n_));
            for (auto& v : x)
                v = static_cast<T>(detail::sd::fmul(double(v), r));
            return x;
        }

        friend constexpr bool operator==(uniform_ball const&, uniform_ball const&) = default;

    private:
        constexpr void check() const
        {
            if (n_ < 1)
                throw std::domain_error("wsl::random::uniform_ball : il faut une dimension >= 1");
        }
        constexpr detail::sd::double_double log_volume() const noexcept
        {
            auto const a = detail::log_sphere_area(n_);
            auto const ln = detail::log_precise(double(n_));
            return detail::sum(a, {-ln.hi, -ln.lo});
        }
        constexpr bool inside(value_type const& x) const noexcept
        {
            return x.size() == n_ && detail::squared_norm_minus_one(x) <= 0;
        }

        std::size_t n_;
    };
}
