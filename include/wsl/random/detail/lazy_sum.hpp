#pragma once

// Sommes exactes de réels paresseux, arrondies correctement (processus en temps
// continu, marches à pas réels ; docs/random.md §14, option B).
//
// Une somme paresseuse tient des termes c_j V_j / q_j : c_j dyadique exact (grand entier
// fois puissance de 2), V_j >= 0 réel paresseux (k + u, u aux bits tirés à la demande),
// q_j > 0 un diviseur flottant pris dans une petite liste (taux d'exponentielles), plus une
// constante dyadique exacte. On arrondit
//
//   offset + α Σ_{j<k} c_j V_j / q_j + β (constante + Σ_j c_j V_j / q_j)
//
// (α, β, offset dyadiques exacts) au plus proche : c'est l'arrondi correct de la valeur
// tirée exactement. Voie rapide : milieux des intervalles en double-double, avec une borne
// d'erreur rigoureuse ; si toute la boule tombe dans la cellule d'arrondi d'un même double,
// c'est lui. Sinon, si la largeur des intervalles en est la cause, tous les termes gagnent
// un mot (64 bits) et l'on recommence (coût amorti : la largeur cumulée baisse de 2^64) ;
// sinon (valeur presque à mi-chemin de deux doubles), voie exacte en grands entiers sur les
// bornes, et raffinement jusqu'à décision. Les deux voies donnent le même résultat et le
// même nombre de mots tirés : la voie rapide ne rend que ce qui est certain.

#include <wsl/random/basic_laws.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/big.hpp>
#include <wsl/random/detail/lazy.hpp>
#include <wsl/random/detail/lazy_root.hpp>
#include <wsl/random/exponential.hpp>
#include <wsl/random/laplace.hpp>
#include <wsl/random/normal.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wsl::random::detail
{
    // ±m 2^e exactement.
    struct dyadic_value
    {
        signed_big m;
        int e = 0;

        constexpr bool is_zero() const noexcept { return m.magnitude.is_zero(); }
    };

    constexpr dyadic_value to_dyadic(double x)
    {
        auto const d = decompose(x);
        if (d.mantissa == 0)
            return {};
        return {{d.negative, big(d.mantissa)}, d.exponent};
    }
    constexpr dyadic_value operator*(dyadic_value const& a, dyadic_value const& b)
    {
        if (a.is_zero() || b.is_zero())
            return {};
        return {a.m * b.m, a.e + b.e};
    }
    constexpr dyadic_value operator-(dyadic_value a)
    {
        if (!a.is_zero())
            a.m.negative = !a.m.negative;
        return a;
    }
    constexpr dyadic_value operator+(dyadic_value const& a, dyadic_value const& b)
    {
        if (a.is_zero())
            return b;
        if (b.is_zero())
            return a;
        int const e = std::min(a.e, b.e);
        signed_big const x{a.m.negative, a.m.magnitude.shifted_left(a.e - e)}, y{b.m.negative, b.m.magnitude.shifted_left(b.e - e)};
        return {x + y, e};
    }
    constexpr dyadic_value operator-(dyadic_value const& a, dyadic_value const& b) { return a + (-b); }

    // Approximation double-double (106 bits de tête) ; erreur relative < 2^-104.
    constexpr sd::double_double to_dd(dyadic_value const& a) noexcept
    {
        if (a.is_zero())
            return {0.0, 0.0};
        int const length = a.m.magnitude.bit_length();
        int const shift = length > 106 ? length - 106 : 0;
        big const top = a.m.magnitude.shifted_right(shift);
        std::uint64_t const high = top.shifted_right(top.bit_length() > 53 ? top.bit_length() - 53 : 0).bits_at(0);
        int const low_bits = top.bit_length() > 53 ? top.bit_length() - 53 : 0;
        std::uint64_t const low = low_bits == 0 ? 0 : top.bits_at(0) & ((low_bits >= 64 ? ~std::uint64_t{0} : (std::uint64_t{1} << low_bits) - 1));
        double const h = sd::scale_by_power_of_two(static_cast<double>(high), a.e + shift + low_bits);
        double const l = sd::scale_by_power_of_two(static_cast<double>(low), a.e + shift);
        auto const [s, t] = sd::fast_two_sum(h, l);
        return a.m.negative ? sd::double_double{-s, -t} : sd::double_double{s, t};
    }

    // ±m 2^e / d au plus proche (pair) en T, d > 0.
    template <class T>
    constexpr T round_rational(dyadic_value const& n, big const& d)
    {
        if (n.is_zero())
            return T(0);
        if (d == big(1))
            return round_to<T>(n.m.negative, n.m.magnitude, false, n.e);
        int const extra = std::max(0, 120 + d.bit_length() - n.m.magnitude.bit_length());
        auto const [q, r] = divide(n.m.magnitude.shifted_left(extra), d);
        return round_to<T>(n.m.negative, q, !r.is_zero(), n.e - extra);
    }

    // Signe de ±m 2^e / d − t (d > 0).
    constexpr int sign_versus(dyadic_value const& n, big const& d, dyadic_value const& t)
    {
        dyadic_value const td = t * dyadic_value{{false, d}, 0};
        dyadic_value const diff = n - td;
        return diff.is_zero() ? 0 : diff.m.negative ? -1 : 1;
    }

    // Réel paresseux tiré exactement, sous la forme offset + c V / q, ou offset + c √v V si
    // root (normale de σ = √v irrationnel : v exact, à l'échelle).
    struct lazy_draw
    {
        dyadic_value offset;
        dyadic_value coefficient;
        double divisor = 1;
        lazy_real v;
        bool root = false;
        scaled variance{1.0, 0.0, 0};
    };

    template <class T, class G>
    constexpr lazy_draw lazy_sample(normal<T> const& X, G& g)
    {
        lazy_real v = exact_half_normal(g);
        bool const negative = (random_bits(g) >> 63) != 0;                  // même ordre que sample(normal)
        bool const root = !X.has_double_stddev();                           // σ = √v irrationnel : terme ±√v V
        dyadic_value c = to_dyadic(root ? 1.0 : double(X.stddev()));
        if (negative)
            c = -c;
        return {to_dyadic(double(X.mean())), c, 1.0, std::move(v), root, X.exact_variance()};
    }
    template <class T, class G>
    constexpr lazy_draw lazy_sample(exponential<T> const& X, G& g)
    {
        if (double const r = X.exact_rate().as_double(); r != 0)
            return {{}, to_dyadic(1.0), r, exact_exponential(g), false, {1.0, 0.0, 0}};
        auto const f = X.exact_rate().sampling_factors();                    // V/λ = V σ / d, exact
        return {{}, to_dyadic(f.sigma), f.divisor, exact_exponential(g), false, {1.0, 0.0, 0}};
    }
    template <class T, class G>
    constexpr lazy_draw lazy_sample(laplace<T> const& X, G& g)
    {
        lazy_real v = exact_exponential(g);
        bool const negative = (random_bits(g) >> 63) != 0;
        dyadic_value c = to_dyadic(double(X.scale()));
        if (negative)
            c = -c;
        return {to_dyadic(double(X.location())), c, 1.0, std::move(v), false, {1.0, 0.0, 0}};
    }
    template <class T, class G>
    constexpr lazy_draw lazy_sample(uniform<T> const& X, G&)
    {
        return {to_dyadic(double(X.lower())), to_dyadic(double(X.upper())) - to_dyadic(double(X.lower())), 1.0, lazy_real{}, false, {1.0, 0.0, 0}};
    }

    template <class L, class G>
    concept lazy_samplable = requires(L const& law, G& g) {
        { lazy_sample(law, g) } -> std::same_as<lazy_draw>;
    };

    class lazy_sum
    {
    public:
        constexpr lazy_sum() : lazy_sum(std::vector<double>{1.0}) {}
        constexpr explicit lazy_sum(std::vector<double> divisors) : divisors_(std::move(divisors))
        {
            for (double const q : divisors_)
            {
                auto const d = decompose(q);
                divisor_mantissa_.push_back(d.mantissa);
                divisor_exponent_.push_back(d.exponent);
            }
        }

        constexpr std::size_t size() const noexcept { return atoms_.size(); }

        // + c V / q_index (ou c √v V si root) ; V raffiné au niveau courant.
        template <class G>
        constexpr void add(dyadic_value c, std::size_t divisor_index, lazy_real v, G& g, bool root = false)
        {
            while (v.fraction.size() < level_)
                v.fraction.refine(g);
            atoms_.push_back({std::move(c), divisor_index, std::move(v), {}, 0, 0, root});
            refresh(atoms_.size() - 1);
            append_prefix();
        }
        template <class G>
        constexpr void add(lazy_draw d, G& g)
        {
            std::size_t index = 0;
            while (index < divisors_.size() && divisors_[index] != d.divisor)
                ++index;
            if (index == divisors_.size())
            {
                divisors_.push_back(d.divisor);
                auto const q = decompose(d.divisor);
                divisor_mantissa_.push_back(q.mantissa);
                divisor_exponent_.push_back(q.exponent);
            }
            if (d.root)
                set_root(d.variance);
            add_constant(d.offset);
            add(std::move(d.coefficient), index, std::move(d.v), g, d.root);
        }
        // Constante exacte (pondérée par β).
        constexpr void add_constant(dyadic_value const& d)
        {
            if (d.is_zero())
                return;
            constant_ = constant_ + d;
            constant_dd_ = to_dd(constant_);
        }

        // Arrondi correct de offset + α Σ_{j<k} c_j V_j/q_j + β (constante + Σ_j c_j V_j/q_j).
        template <class T, class G>
        constexpr T rounded(G& g, dyadic_value const& offset, dyadic_value const& alpha, std::size_t k, dyadic_value const& beta)
        {
            for (;;)
            {
                T quick{};
                auto const verdict = level_ > 40 ? quick_verdict::undecided : quickly(offset, alpha, k, beta, quick);   // au-delà : voie exacte
                if (verdict == quick_verdict::decided)
                    return quick;
                if (verdict == quick_verdict::too_wide)
                {
                    level_up(g);
                    continue;
                }
                auto const [low, width, root_low, root_width, denominator] = exact_bounds(offset, alpha, k, beta);
                T const a = round_end<T>(low, root_low, denominator), c = round_end<T>(low + width, root_low + root_width, denominator);
                using bits_type = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;
                if (std::bit_cast<bits_type>(a) == std::bit_cast<bits_type>(c))
                    return a;
                level_up(g);
            }
        }
        template <class G>
        constexpr double rounded(G& g) { return rounded<double>(g, {}, {}, 0, to_dyadic(1.0)); }

        // Signe de (valeur − t), la valeur comme dans rounded.
        template <class G>
        constexpr int compare(G& g, dyadic_value const& offset, dyadic_value const& alpha, std::size_t k, dyadic_value const& beta, dyadic_value const& t)
        {
            for (;;)
            {
                auto const [y, bound, refinable] = approximate(offset, alpha, k, beta);
                auto const td = to_dd(t);
                bool const usable = level_ <= 40 && bound < std::numeric_limits<double>::infinity() && is_finite(y.hi);
                auto const d = sum(y, {-td.hi, -td.lo});
                double const gap = d.hi + d.lo, margin = bound + sd::fmul(magnitude(td.hi) + magnitude(y.hi), 0x1p-100);
                if (usable && gap > margin)
                    return 1;
                if (usable && gap < -margin)
                    return -1;
                if (usable && refinable > sd::fmul(magnitude(td.hi) + 0x1p-1000, 0x1p-40) && refinable > sd::fmul(bound, 0.5))
                {
                    level_up(g);
                    continue;
                }
                auto const [low, width, root_low, root_width, denominator] = exact_bounds(offset, alpha, k, beta);
                int const a = sign_end(low, root_low, denominator, t), b = sign_end(low + width, root_low + root_width, denominator, t);
                if (a == b && a != 0)
                    return a;
                level_up(g);
            }
        }

        // Contrôle (tests) : après extra mots de plus sur chaque terme, les deux bornes exactes
        // s'arrondissent-elles à r ?
        template <class T, class G>
        constexpr bool confirms(G& g, T r, int extra, dyadic_value const& offset, dyadic_value const& alpha, std::size_t k, dyadic_value const& beta)
        {
            for (int i = 0; i < extra; ++i)
                level_up(g);
            auto const [low, width, root_low, root_width, denominator] = exact_bounds(offset, alpha, k, beta);
            using bits_type = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;
            return std::bit_cast<bits_type>(round_end<T>(low, root_low, denominator)) == std::bit_cast<bits_type>(r)
                && std::bit_cast<bits_type>(round_end<T>(low + width, root_low + root_width, denominator)) == std::bit_cast<bits_type>(r);
        }

    private:
        struct atom
        {
            dyadic_value c;
            std::size_t divisor;
            lazy_real v;
            sd::double_double mid;                                          // ≈ c mid(V) / q
            double width;                                                   // >= |c| 2^-b / q (réductible)
            double error;                                                   // >= |mid − c mid(V)/q| (fixe)
            bool root = false;                                              // terme c √v V
        };

        // √v commun des termes racine (une seule variance par somme : pas de racines mêlées).
        constexpr void set_root(scaled const& v)
        {
            if (has_root_)
            {
                if (!(v.hi == root_variance_.hi && v.lo == root_variance_.lo && v.e == root_variance_.e))
                    throw std::domain_error("wsl::random : somme exacte de normales de variances différentes, à racines irrationnelles");
                return;
            }
            has_root_ = true;
            root_variance_ = v;
            root_square_ = exact_square_of(v);
            auto const r = sqrt(v);
            sigma_dd_ = unscaled(r);
            if (!is_finite(sigma_dd_.hi) || sigma_dd_.hi == 0)
                sigma_dd_ = {std::numeric_limits<double>::infinity(), 0.0};   // hors de portée : voie exacte
        }

        // Arrondi correct de (A + √v B) / d, et signe de (A + √v B) / d − t (B nul : rationnel).
        template <class T>
        constexpr T round_end(dyadic_value const& a, dyadic_value const& b, big const& d) const
        {
            if (b.is_zero())
                return round_rational<T>(a, d);
            exact_dyadic const ea{a.m, a.e}, eb{b.m, b.e};
            auto const side = [&](exact_dyadic const& x) {
                exact_dyadic const xd{{x.m.negative, x.m.magnitude * d}, x.e};   // x · d
                return sign_plus_root(exact_sum(ea, exact_negated(xd)), root_square_, eb);
            };
            auto const ad = to_dd(a), bd = to_dd(b), dd = to_dd(dyadic_value{{false, d}, 0});
            double const candidate = (ad.hi + sd::fmul(sigma_dd_.hi, bd.hi)) / dd.hi;
            return round_by_side<T>(side, static_cast<T>(is_finite(candidate) ? candidate : 0.0));
        }
        constexpr int sign_end(dyadic_value const& a, dyadic_value const& b, big const& d, dyadic_value const& t) const
        {
            if (b.is_zero())
                return sign_versus(a, d, t);
            dyadic_value const td = t * dyadic_value{{false, d}, 0};
            dyadic_value const p = a - td;
            return sign_plus_root(exact_dyadic{p.m, p.e}, root_square_, exact_dyadic{b.m, b.e});
        }

        enum class quick_verdict
        {
            decided,
            too_wide,
            undecided
        };

        // Milieu et largeur du terme j (bornes rigoureuses : largeur arrondie vers le haut).
        constexpr void refresh(std::size_t j)
        {
            atom& a = atoms_[j];
            int const b = a.v.fraction.bits();
            // milieu de V : k + (A + ½) 2^-b, sur les deux premiers mots (reste dans l'erreur)
            double const w0 = a.v.fraction.size() > 0 ? double(a.v.fraction.word(0) >> 11) : 0.0;
            double const w0_low = a.v.fraction.size() > 0 ? double(a.v.fraction.word(0) & 0x7ff) : 0.0;
            double const w1 = a.v.fraction.size() > 1 ? double(a.v.fraction.word(1) >> 11) : 0.0;
            sd::double_double mv{double(a.v.integer), 0.0};
            mv = sum(mv, {sd::scale_by_power_of_two(w0, -53), 0.0});
            mv = sum(mv, {sd::scale_by_power_of_two(w0_low, -64), 0.0});
            mv = sum(mv, {sd::scale_by_power_of_two(w1, -117), 0.0});
            mv = sum(mv, {sd::scale_by_power_of_two(1.0, -b - 1), 0.0});    // + ½ ulp de l'intervalle
            auto c = to_dd(a.c);
            if (a.root)
                c = mul(c, sigma_dd_);                                      // √v à 2^-100 près : erreur ci-dessous
            auto m = mul(c, mv);
            double const q = divisors_[a.divisor];
            if (q != 1)
                m = sd::quotient_dd(m, {q, 0.0});
            a.mid = m;
            double const cw = sd::fmul(sd::fmul(magnitude(c.hi), 1 + 0x1p-50) / q, 1 + 0x1p-50);
            // largeur de l'intervalle (réductible) ; erreur du milieu (bits au-delà du deuxième mot
            // ignorés, arrondis du double-double, partie entière grande), qui ne diminue pas
            a.width = sd::scale_by_power_of_two(cw, -b);
            a.error = sd::fmul(magnitude(m.hi), a.root ? 0x1p-96 : 0x1p-100) + (a.v.fraction.size() > 2 ? sd::scale_by_power_of_two(cw, -116) : 0.0);
            if (a.v.integer >= (std::uint64_t{1} << 50))
                a.error = a.error + sd::fmul(magnitude(m.hi), 0x1p-90);
            if (!is_finite(a.width) || !is_finite(a.error) || !is_finite(m.hi) || (c.hi == 0 && !a.c.is_zero()))
                a.error = std::numeric_limits<double>::infinity();             // hors de portée des doubles : voie exacte
        }

        constexpr void append_prefix()
        {
            std::size_t const j = atoms_.size() - 1;
            auto const& a = atoms_[j];
            prefix_mid_.push_back(j == 0 ? a.mid : sum(prefix_mid_.back(), a.mid));
            double const previous_error = j == 0 ? 0.0 : prefix_error_.back();
            double const previous_width = j == 0 ? 0.0 : prefix_width_.back();
            prefix_error_.push_back(sd::fmul(previous_error + a.error + sd::fmul(magnitude(prefix_mid_.back().hi) + magnitude(a.mid.hi), 0x1p-102), 1 + 0x1p-50));
            prefix_width_.push_back(sd::fmul(previous_width + a.width, 1 + 0x1p-50));
        }

        constexpr void rebuild_prefix()
        {
            prefix_mid_.clear();
            prefix_error_.clear();
            prefix_width_.clear();
            for (std::size_t j = 0; j < atoms_.size(); ++j)
            {
                refresh(j);
                auto const& a = atoms_[j];
                prefix_mid_.push_back(j == 0 ? a.mid : sum(prefix_mid_.back(), a.mid));
                double const previous_error = j == 0 ? 0.0 : prefix_error_.back();
                double const previous_width = j == 0 ? 0.0 : prefix_width_.back();
                prefix_error_.push_back(sd::fmul(previous_error + a.error + sd::fmul(magnitude(prefix_mid_.back().hi) + magnitude(a.mid.hi), 0x1p-102), 1 + 0x1p-50));
                prefix_width_.push_back(sd::fmul(previous_width + a.width, 1 + 0x1p-50));
            }
        }

        // Tous les termes gagnent un mot.
        template <class G>
        constexpr void level_up(G& g)
        {
            ++level_;
            for (auto& a : atoms_)
                while (a.v.fraction.size() < level_)
                    a.v.fraction.refine(g);
            rebuild_prefix();
        }

        // Valeur approchée et borne rigoureuse de l'écart à la vraie valeur.
        struct approximation
        {
            sd::double_double y;
            double bound;                                                   // écart total à la vraie valeur
            double refinable;                                               // part due aux largeurs (raffiner la réduit)
        };
        // Termes j < k de poids α + β, les autres de poids β : y = offset + (α + β) P_k
        // + β (P_n − P_k + constante), largeurs pondérées de même (les poids qui s'annulent
        // n'apportent aucune largeur).
        constexpr approximation approximate(dyadic_value const& offset, dyadic_value const& alpha, std::size_t k, dyadic_value const& beta) const
        {
            auto const o = to_dd(offset), a1 = to_dd(alpha + beta), be = to_dd(beta);
            std::size_t const n = atoms_.size();
            sd::double_double const pk = k == 0 ? sd::double_double{0.0, 0.0} : prefix_mid_[k - 1];
            sd::double_double const pn = n == 0 ? sd::double_double{0.0, 0.0} : prefix_mid_[n - 1];
            sd::double_double const rest = sum(sum(pn, {-pk.hi, -pk.lo}), constant_dd_);
            double const wk = k == 0 ? 0.0 : prefix_width_[k - 1], wn = n == 0 ? 0.0 : prefix_width_[n - 1];
            double const ek = k == 0 ? 0.0 : prefix_error_[k - 1], en = n == 0 ? 0.0 : prefix_error_[n - 1];
            // largeur des termes j >= k : W_n − W_k, W_k calculé majoré d'au plus 2^-20 relativement
            double const w_rest = k >= n ? 0.0 : wn - wk / (1 + 0x1p-20) > 0 ? sd::fmul(wn - wk / (1 + 0x1p-20), 1 + 0x1p-50) : 0.0;
            sd::double_double y = o;
            if (!(alpha + beta).is_zero())
                y = sum(y, mul(a1, pk));
            if (!beta.is_zero())
                y = sum(y, mul(be, rest));
            double const refinable = sd::fmul(sd::fmul(magnitude(a1.hi), sd::fmul(wk, 0.5)) + sd::fmul(magnitude(be.hi), sd::fmul(w_rest, 0.5)), 1 + 0x1p-40);
            double const magnitude_sum = magnitude(o.hi) + sd::fmul(magnitude(a1.hi), magnitude(pk.hi))
                                       + sd::fmul(magnitude(be.hi), magnitude(pn.hi) + magnitude(pk.hi) + magnitude(constant_dd_.hi));
            double const fixed = sd::fmul(sd::fmul(magnitude(a1.hi), ek) + sd::fmul(magnitude(be.hi), en + ek), 1 + 0x1p-40) + sd::fmul(magnitude_sum, 0x1p-97);
            return {y, refinable + fixed, refinable};
        }

        template <class T>
        constexpr quick_verdict quickly(dyadic_value const& offset, dyadic_value const& alpha, std::size_t k, dyadic_value const& beta, T& result) const
        {
            using F = format<T>;
            auto const [y, bound, refinable] = approximate(offset, alpha, k, beta);
            if (!(bound < std::numeric_limits<double>::infinity()) || !is_finite(y.hi))
                return quick_verdict::undecided;                            // voie exacte
            auto const r = static_cast<T>(y.hi);
            double const a = magnitude(double(r));
            if (!(a >= double(std::numeric_limits<T>::min()) * 2 && a < double(std::numeric_limits<T>::max()) / 2))
                return a == 0 && refinable > 0x1p-1000 ? quick_verdict::too_wide : quick_verdict::undecided;
            int const e = sd::exponent_of(a);
            double const up = sd::power_of_two(e - F::digits + 1);
            double const down = a == sd::power_of_two(e) ? sd::fmul(up, 0.5) : up;
            double deviation = (y.hi - double(r)) + y.lo;
            if (r < 0)
                deviation = -deviation;
            double const margin = 1 - 0x1p-30;
            if (deviation + bound < sd::fmul(sd::fmul(up, 0.5), margin) && bound - deviation < sd::fmul(sd::fmul(down, 0.5), margin))
            {
                result = r;
                return quick_verdict::decided;
            }
            // Raffiner seulement si les largeurs dominent (les erreurs d'arrondi du double-double ne
            // diminuent pas : près de 0, c'est la voie exacte qui tranche).
            return refinable > sd::fmul(down, 0x1p-12) && refinable > sd::fmul(bound, 0.5) ? quick_verdict::too_wide : quick_verdict::undecided;
        }

        // Bornes exactes : valeur ∈ [low, low + width] / denominator.
        struct exact
        {
            dyadic_value low, width;                                        // partie dyadique
            dyadic_value root_low, root_width;                              // coefficient de √v
            big denominator;
        };
        constexpr exact exact_bounds(dyadic_value const& offset, dyadic_value const& alpha, std::size_t k, dyadic_value const& beta) const
        {
            // dénominateur commun : produit des mantisses des diviseurs employés
            std::vector<bool> used(divisors_.size(), false);
            for (auto const& a : atoms_)
                used[a.divisor] = true;
            big denominator(1);
            for (std::size_t i = 0; i < divisors_.size(); ++i)
                if (used[i])
                    denominator = denominator * big(divisor_mantissa_[i]);
            dyadic_value const den{{false, denominator}, 0};
            dyadic_value low = offset * den + beta * constant_ * den, width, root_low, root_width;
            for (std::size_t j = 0; j < atoms_.size(); ++j)
            {
                auto const& a = atoms_[j];
                dyadic_value weight = beta;
                if (j < k)
                    weight = weight + alpha;
                if (weight.is_zero())
                    continue;
                dyadic_value const gamma = weight * a.c;
                int const b = a.v.fraction.bits();
                big const lo = big(a.v.integer).shifted_left(b) + a.v.fraction.numerator();
                // γ V / q = γ (lo + θ) 2^{−b} 2^{−e_q} (D/m_q) / D, θ ∈ [0, 1]
                big const share = divide(denominator, big(divisor_mantissa_[a.divisor])).quotient;
                dyadic_value const unit = gamma * dyadic_value{{false, share}, -b - divisor_exponent_[a.divisor]};
                dyadic_value const bottom = unit * dyadic_value{{false, gamma.m.negative ? lo + big(1) : lo}, 0};   // borne basse de γ V
                dyadic_value magnitude_unit = unit;
                magnitude_unit.m.negative = false;
                if (a.root)
                {
                    root_low = root_low + bottom;
                    root_width = root_width + magnitude_unit;
                }
                else
                {
                    low = low + bottom;
                    width = width + magnitude_unit;
                }
            }
            return {low, width, root_low, root_width, denominator};
        }

        std::vector<double> divisors_;
        std::vector<std::uint64_t> divisor_mantissa_;
        std::vector<int> divisor_exponent_;
        std::vector<atom> atoms_;
        std::vector<sd::double_double> prefix_mid_;
        std::vector<double> prefix_error_, prefix_width_;
        dyadic_value constant_;
        sd::double_double constant_dd_{0.0, 0.0};
        int level_ = 1;
        bool has_root_ = false;
        scaled root_variance_{1.0, 0.0, 0};
        exact_square root_square_{};
        sd::double_double sigma_dd_{1.0, 0.0};
    };
}
