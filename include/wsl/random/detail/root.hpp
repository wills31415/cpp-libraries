#pragma once

// Solveur commun des quantiles sans forme close : u tel que V(u) = t, V la plus petite part
// d'une loi (répartition ou survie, t <= ½ en pratique), monotone en u.
//
// Newton sur log V : u ← u − log(V/t) · V/V′, encadré (bissection dès que le pas sort de
// l'encadrement, ou dépasse max_step) ; tant qu'un côté de l'encadrement manque, la loi
// fournit le pas vers l'extérieur. Convergence quand |log(V/t)| < 2^-36 ou que le pas est
// plus fin que les doubles : deux pas à résidu précis, en double-double, depuis le point
// effectivement évalué (resultat en double-double).

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/discrete.hpp>

#include <limits>

namespace wsl::random::detail
{
    // Ce que la loi rend en un point : V(u), |dV/du| (0 : pas de pas de Newton), et le point
    // effectivement évalué (u arrondi par la loi, par exemple x = e^u en double).
    struct part_value
    {
        scaled value, slope;
        sd::double_double at;
    };

    struct root_bounds
    {
        double start;
        double lo = -std::numeric_limits<double>::infinity(), hi = std::numeric_limits<double>::infinity();
        double max_step = std::numeric_limits<double>::infinity();
    };

    // evaluate(u) -> part_value ; outward(u, side) -> u suivant (side > 0 : u trop grand) ;
    // midpoint(lo, hi) -> point de bissection.
    template <class Evaluate, class Outward, class Midpoint>
    constexpr sd::double_double solve_part(Evaluate const& evaluate, bool increasing, double t, root_bounds b,
                                           Outward const& outward, Midpoint const& midpoint, int iterations = 400)
    {
        auto const log_t = log_precise(t);
        struct step
        {
            int side;                                                       // > 0 : u trop grand
            sd::double_double du;
            double residual;                                                // |log(V/t)|
            bool finite;
            sd::double_double at;
        };
        auto const newton = [&](sd::double_double u) {
            part_value const p = evaluate(u);
            int const c = p.value.hi == 0 ? -1 : compare(p.value, {t, 0.0});
            step s{increasing ? c : -c, {0.0, 0.0}, 0.0, false, p.at};
            if (p.value.hi == 0 || p.slope.hi == 0)
                return s;
            auto const h = sum(log(p.value), negated(log_t));               // log(V/t)
            auto const r = unscaled(p.value / p.slope);                     // V/|V′|
            auto const d = mul(h, r);
            s.du = increasing ? negated(d) : d;                             // u − h V/V′
            s.residual = magnitude(h.hi);
            s.finite = is_finite(s.du.hi);
            return s;
        };
        double lo = b.lo, hi = b.hi, u = b.start;
        constexpr double infinity = std::numeric_limits<double>::infinity();
        for (int iteration = 0; iteration < iterations; ++iteration)
        {
            step const s = newton({u, 0.0});
            if (s.side == 0)
                return s.at;
            if (s.side > 0)
                hi = u < hi ? u : hi;
            else
                lo = u > lo ? u : lo;
            double next = u + s.du.hi;
            if (s.finite && (0x1p-36 > s.residual || next == u))            // deux pas à résidu précis
            {
                sd::double_double r = sum(s.at, s.du);
                step const s2 = newton(r);
                if (s2.finite)
                    r = sum(s2.at, s2.du);
                return r;
            }
            bool const inside = s.finite && next > lo && next < hi && !(s.du.hi > b.max_step) && !(-b.max_step > s.du.hi);
            if (!inside)
                next = lo > -infinity && hi < infinity ? midpoint(lo, hi) : outward(u, s.side);
            if (next == u)
                return s.at;
            u = next;
        }
        return {u, 0.0};
    }

    // Bissection arithmétique (défaut).
    inline constexpr auto arithmetic_midpoint = [](double lo, double hi) { return sd::fmul(lo, 0.5) + sd::fmul(hi, 0.5); };
}
