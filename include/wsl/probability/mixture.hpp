#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/transition.hpp>   // transition (étape 5a)

// Algèbre des lois : mélange fini de lois (voir design/random.md §4.6).
//
//   wsl::probability::mixture M(wsl::probability::weighted(0.3, wsl::probability::normal({.mean = -2, .stddev = 1})),
//                  wsl::probability::weighted(0.7, wsl::probability::gamma({.rate = 1, .shape = 3})));
//   M.pdf(0.5); M.cdf(1.0); M.quantile(0.9); M.mean(); M.variance(); M(g);
//
// Composantes de types quelconques, de même value_type (réel ou entier) ; poids > 0
// quelconques, normalisés par leur somme W. Densité (masse), répartition et survie
// Σ wᵢ fᵢ / W en double-double, un arrondi final sur les valeurs des composantes ;
// moyenne Σ wᵢ mᵢ/W, variance Σ wᵢ (vᵢ + (mᵢ − m)²)/W (sans compensation),
// asymétrie et kurtosis par les moments centrés ; quantile et cquantile par
// dichotomie sur les flottants (les entiers pour une loi discrète) entre les quantiles
// extrêmes des composantes, comparaisons en double-double sur la répartition (ou la
// survie, du côté où elle est petite).
//
// Tirage exact si celui de chaque composante l'est : la composante est choisie avec
// probabilité wᵢ/W exactement (les poids, flottants, sont des entiers à un facteur 2^e
// commun près : un entier uniforme sous leur somme, sur 64 bits ou en grands entiers),
// puis tirée.

#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/random/detail/big.hpp>
#include <wsl/random/detail/choice.hpp>
#include <wsl/probability/detail/discrete.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/random/uniform.hpp>
#include <wsl/special.hpp>

#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace wsl::inline WSL_DEFAULTS::probability
{
    // Composante pondérée d'un mélange.
    template <probability_law L>
    struct weighted
    {
        double weight;
        L law;

        constexpr weighted(double w, L const& l) : weight(w), law(l) {}
    };

    namespace detail
    {
        // Clé entière croissante avec le flottant (−0 juste sous +0).
        template <std::floating_point T>
        constexpr auto ordered_key(T x) noexcept
        {
            using U = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;
            constexpr U sign = U(1) << (8 * sizeof(T) - 1);
            U const b = std::bit_cast<U>(x);
            return (b & sign) != 0 ? U(~b) : U(b | sign);
        }
        template <std::floating_point T, class U>
        constexpr T from_ordered_key(U k) noexcept
        {
            constexpr U sign = U(1) << (8 * sizeof(T) - 1);
            return std::bit_cast<T>((k & sign) != 0 ? U(k & ~sign) : U(~k));
        }

        // Plus petit flottant x de [low, high] tel que pred(x) (pred croissant ; high
        // rendu si pred n'est jamais vrai) : dichotomie sur les clés, au plus 64 pas.
        template <std::floating_point T, class P>
        constexpr T smallest_real(T low, T high, P pred)
        {
            if (pred(low))
                return low;
            auto a = ordered_key(low), b = ordered_key(high);                 // pred faux en a
            while (b - a > 1)
            {
                auto const m = a + (b - a) / 2;
                if (pred(from_ordered_key<T>(m)))
                    b = m;
                else
                    a = m;
            }
            return from_ordered_key<T>(b);
        }


    }

    template <probability_law... L>
        requires(sizeof...(L) >= 1)
    class mixture : public law_interface
    {
        using first = std::tuple_element_t<0, std::tuple<L...>>;
        static constexpr std::size_t n = sizeof...(L);

    public:
        using value_type = typename first::value_type;
        static_assert((std::is_same_v<typename L::value_type, value_type> && ...), "wsl::probability::mixture : composantes de même value_type");
        static constexpr bool continuous = std::is_floating_point_v<value_type>;
        using moment_type = std::conditional_t<continuous, value_type, double>;   // comme les lois de base
        static constexpr exactness sampling = ((L::sampling == exactness::exact) && ...) ? exactness::exact : exactness::approximate;

        constexpr explicit mixture(weighted<L> const&... c) : laws_(c.law...), w_{c.weight...}
        {
            for (double const w : w_)
                if (!(w > 0) || !detail::is_finite(w))
                    throw std::domain_error("wsl::probability::mixture : les poids doivent être > 0, finis");
            total_weight_ = {0.0, 0.0};
            for (double const w : w_)
                total_weight_ = detail::sum(total_weight_, {w, 0.0});
            integer_weights();
        }

        static constexpr std::size_t size() noexcept { return n; }
        constexpr double weight(std::size_t i) const noexcept { return w_[i]; }   // poids donné (non normalisé)
        template <std::size_t I>
        constexpr auto const& component() const noexcept { return std::get<I>(laws_); }

        constexpr interval<value_type> support() const noexcept
        {
            interval<value_type> s = std::get<0>(laws_).support();
            apply_all([&](auto const& c) {
                auto const t = c.support();
                s.lower = t.lower < s.lower ? t.lower : s.lower;
                s.upper = t.upper > s.upper ? t.upper : s.upper;
            });
            return s;
        }

        constexpr moment_type mean() const noexcept requires(has_mean<L> && ...)
        {
            auto const s = average([](auto const& c) { return double(c.mean()); });
            return static_cast<moment_type>(s.hi + s.lo);
        }
        // Σ wᵢ (vᵢ + (mᵢ − m)²) / W
        constexpr moment_type variance() const noexcept requires(has_mean<L> && ...) && (has_variance<L> && ...)
        {
            return static_cast<moment_type>(variance_double());
        }
        constexpr moment_type skewness() const noexcept requires(has_mean<L> && ...) && (has_variance<L> && ...) && (has_skewness<L> && ...)
        {
            namespace fp = detail::sd;
            double const m = mean_double(), v = variance_double();
            auto const s = average([m](auto const& c) {                    // E[(X − m)³ | i]
                double const d = double(c.mean()) - m, vi = c.variance();
                double const si = fp::fmul(double(c.skewness()), fp::fmul(vi, wsl::special::sqrt(vi)));
                return (si + fp::fmul(fp::fmul(3.0, vi), d)) + fp::fmul(fp::fmul(d, d), d);
            });
            return static_cast<moment_type>((s.hi + s.lo) / (v * wsl::special::sqrt(v)));
        }
        constexpr moment_type excess_kurtosis() const noexcept
            requires(has_mean<L> && ...) && (has_variance<L> && ...) && (has_skewness<L> && ...) && (has_excess_kurtosis<L> && ...)
        {
            namespace fp = detail::sd;
            double const m = mean_double(), v = variance_double();
            auto const s = average([m](auto const& c) {                    // E[(X − m)⁴ | i]
                double const d = double(c.mean()) - m, vi = c.variance();
                double const k = fp::fmul(double(c.excess_kurtosis()) + 3.0, fp::fmul(vi, vi));
                double const t3 = fp::fmul(fp::fmul(4.0, d), fp::fmul(double(c.skewness()), fp::fmul(vi, wsl::special::sqrt(vi))));
                double const d2 = fp::fmul(d, d);
                return ((k + t3) + fp::fmul(fp::fmul(6.0, vi), d2)) + fp::fmul(d2, d2);
            });
            return static_cast<moment_type>((s.hi + s.lo) / fp::fmul(v, v) - 3.0);
        }
        constexpr value_type median() const requires(has_quantile<L> && ...) && (has_cdf<L> && ...) && (has_ccdf<L> && ...)
        {
            return wsl::probability::quantile(*this, 0.5);
        }

        constexpr value_type pdf(value_type x) const noexcept requires continuous && (has_pdf<L> && ...)
        {
            auto const& M = *this;
            auto const s = M.average([x](auto const& c) { return double(wsl::probability::pdf(c, x)); });
            return static_cast<value_type>(s.hi + s.lo);
        }
        constexpr double pmf(value_type k) const noexcept requires(!continuous) && (has_pmf<L> && ...)
        {
            auto const& M = *this;
            auto const s = M.average([k](auto const& c) { return double(wsl::probability::pmf(c, k)); });
            return s.hi + s.lo;
        }
        // log Σ wᵢ e^{lᵢ} − log W, centré sur le plus grand terme.
        constexpr value_type log_pdf(value_type x) const noexcept
            requires continuous && (requires(L const& c, value_type v) { wsl::probability::log_pdf(c, v); } && ...)
        {
            auto const& M = *this;
            std::array<double, n> l{};
            std::size_t i = 0;
            M.apply_all([&](auto const& c) {
                double const li = wsl::probability::log_pdf(c, x);
                if (!(li > -std::numeric_limits<double>::infinity()))       // densité nulle (ou NaN) : à part
                {
                    l[i++] = li;
                    return;
                }
                auto const t = detail::sum({li, 0.0}, detail::log_precise(M.w_[i]));
                l[i++] = t.hi + t.lo;
            });
            double top = -std::numeric_limits<double>::infinity();
            for (double const v : l)
                top = v > top ? v : top;
            if (!(top > -std::numeric_limits<double>::infinity()))
                return static_cast<value_type>(top);
            double s = 0;
            for (double const v : l)
            {
                if (!(v == v))
                    return static_cast<value_type>(v);
                if (v > -std::numeric_limits<double>::infinity())
                    s = s + wsl::special::exp(v - top);
            }
            auto const lw = detail::log(detail::normalized(M.total_weight_.hi, M.total_weight_.lo));
            auto const r = detail::sum({top, 0.0}, detail::sum(detail::log_precise(s), {-lw.hi, -lw.lo}));
            return static_cast<value_type>(r.hi + r.lo);
        }

        constexpr std::conditional_t<continuous, value_type, double> cdf(value_type x) const noexcept requires(has_cdf<L> && ...)
        {
            auto const& M = *this;
            auto const s = M.cdf_dd(x);
            return static_cast<std::conditional_t<continuous, value_type, double>>(clamp01(s.hi + s.lo));
        }
        constexpr std::conditional_t<continuous, value_type, double> ccdf(value_type x) const noexcept requires(has_ccdf<L> && ...)
        {
            auto const& M = *this;
            auto const s = M.ccdf_dd(x);
            return static_cast<std::conditional_t<continuous, value_type, double>>(clamp01(s.hi + s.lo));
        }

        // Plus petit x tel que F(x) >= p, entre les quantiles extrêmes des composantes.
        constexpr value_type quantile(double p) const requires(has_quantile<L> && ...) && (has_cdf<L> && ...) && (has_ccdf<L> && ...)
        {
            auto const& M = *this;
            if (!(p >= 0 && p <= 1))
                throw std::domain_error("wsl::probability::quantile : p doit être dans [0, 1]");
            auto const [low, high] = M.bracket([p](auto const& c) { return static_cast<value_type>(wsl::probability::quantile(c, p)); });
            if (p == 0 || p == 1 || low == high)
                return p == 0 ? low : high;
            double const q = 1.0 - p;                                       // utilisé pour p > ½ : exact
            auto const ge = [&](value_type x) { return p <= 0.5 ? detail::at_least(M.cdf_dd(x), p) : detail::at_most(M.ccdf_dd(x), q); };
            if constexpr (continuous)
                return detail::smallest_real(low, high, ge);
            else
                return detail::smallest_true(low, low, high, ge);
        }
        // Plus petit x tel que P(X > x) <= q.
        constexpr value_type cquantile(double q) const requires(has_cquantile<L> && ...) && (has_cdf<L> && ...) && (has_ccdf<L> && ...)
        {
            auto const& M = *this;
            if (!(q >= 0 && q <= 1))
                throw std::domain_error("wsl::probability::cquantile : q doit être dans [0, 1]");
            auto const [low, high] = M.bracket([q](auto const& c) { return static_cast<value_type>(wsl::probability::cquantile(c, q)); });
            if (q == 0 || q == 1 || low == high)
                return q == 1 ? low : high;
            double const p = 1.0 - q;                                       // utilisé pour q > ½ : exact
            auto const ok = [&](value_type x) { return q <= 0.5 ? detail::at_most(M.ccdf_dd(x), q) : detail::at_least(M.cdf_dd(x), p); };
            if constexpr (continuous)
                return detail::smallest_real(low, high, ok);
            else
                return detail::smallest_true(low, low, high, ok);
        }

        template <std::uniform_random_bit_generator G>
            requires(samplable<L, G> && ...)
        constexpr value_type sample(G& g) const
        {
            auto const& M = *this;
            return M.visit(M.choose(g), [&g](auto const& c) { return static_cast<value_type>(detail::native_sample(c, g)); });
        }
        template <std::uniform_random_bit_generator G>
            requires(samplable<L, G> && ...)
        constexpr value_type sample(G& g, approximate_t) const
        {
            auto const& M = *this;
            return M.visit(M.choose(g), [&g](auto const& c) { return static_cast<value_type>(wsl::probability::sample(c, g, wsl::approximate)); });
        }

        friend constexpr bool operator==(mixture const& a, mixture const& b) noexcept { return a.laws_ == b.laws_ && a.w_ == b.w_; }

    private:
        constexpr double mean_double() const noexcept
        {
            auto const s = average([](auto const& c) { return double(c.mean()); });
            return s.hi + s.lo;
        }
        constexpr double variance_double() const noexcept
        {
            double const m = mean_double();
            auto const s = average([m](auto const& c) {
                double const d = double(c.mean()) - m;
                return double(c.variance()) + detail::sd::fmul(d, d);
            });
            return s.hi + s.lo;
        }

        static constexpr double clamp01(double v) noexcept { return v < 0 ? 0.0 : v > 1 ? 1.0 : v; }

        template <class F>
        constexpr void apply_all(F&& f) const
        {
            std::apply([&](auto const&... c) { (f(c), ...); }, laws_);
        }

        template <class F>
        constexpr value_type visit(std::size_t i, F&& f) const
        {
            value_type r{};
            [&]<std::size_t... I>(std::index_sequence<I...>) {
                ((i == I ? (void)(r = f(std::get<I>(laws_))) : void()), ...);
            }(std::index_sequence_for<L...>{});
            return r;
        }

        // Σ wᵢ f(composante i) / W, en double-double.
        template <class F>
        constexpr detail::sd::double_double average(F&& f) const
        {
            detail::sd::double_double s{0.0, 0.0};
            std::size_t i = 0;
            apply_all([&](auto const& c) {
                double const v = f(c);
                s = detail::sum(s, detail::product(w_[i++], v));
            });
            if (!detail::is_finite(s.hi))
                return {s.hi, 0.0};
            return detail::sd::quotient_dd(s, total_weight_);
        }

        constexpr detail::sd::double_double cdf_dd(value_type x) const { return average([x](auto const& c) { return double(wsl::probability::cdf(c, x)); }); }
        constexpr detail::sd::double_double ccdf_dd(value_type x) const { return average([x](auto const& c) { return double(wsl::probability::ccdf(c, x)); }); }

        // [min, max] des quantiles des composantes.
        template <class Q>
        constexpr interval<value_type> bracket(Q q) const
        {
            interval<value_type> b{};
            bool started = false;
            apply_all([&](auto const& c) {
                value_type const v = q(c);
                if (!started || v < b.lower)
                    b.lower = v;
                if (!started || v > b.upper)
                    b.upper = v;
                started = true;
            });
            return b;
        }

        // Poids wᵢ = mᵢ 2^{eᵢ}, mᵢ entiers : Nᵢ = mᵢ 2^{eᵢ − min e}, cumuls sur 64 bits s'ils
        // y tiennent, en grands entiers sinon.
        constexpr void integer_weights()
        {
            std::array<std::uint64_t, n> m{};
            std::array<int, n> e{};
            int low = std::numeric_limits<int>::max(), top = std::numeric_limits<int>::min();
            for (std::size_t i = 0; i < n; ++i)
            {
                std::uint64_t const b = std::bit_cast<std::uint64_t>(w_[i]);
                int const biased = static_cast<int>(b >> 52);
                std::uint64_t const f = b & 0x000fffffffffffff;
                m[i] = biased == 0 ? f : f | 0x0010000000000000;
                e[i] = biased == 0 ? -1074 : biased - 1075;
                int const z = std::countr_zero(m[i]);                       // mantisse impaire : moins de bits
                m[i] >>= z;
                e[i] += z;
                low = e[i] < low ? e[i] : low;
            }
            for (std::size_t i = 0; i < n; ++i)
            {
                int const t = static_cast<int>(std::bit_width(m[i])) + e[i];
                top = t > top ? t : top;
            }
            if (top - low + static_cast<int>(std::bit_width(n)) <= 63)
            {
                std::uint64_t c = 0;
                for (std::size_t i = 0; i < n; ++i)
                {
                    c += m[i] << (e[i] - low);
                    cumulative_[i] = c;
                }
            }
            else
            {
                detail::big c;
                for (std::size_t i = 0; i < n; ++i)
                {
                    c = c + detail::big(m[i]).shifted_left(e[i] - low);
                    cumulative_big_.push_back(c);
                }
            }
        }

        // Indice de composante, de probabilité wᵢ/W exactement.
        template <class G>
        constexpr std::size_t choose(G& g) const
        {
            if constexpr (n == 1)
                return 0;
            else if (cumulative_big_.empty())
            {
                std::uint64_t const r = uniform_below(cumulative_[n - 1], g);
                for (std::size_t i = 0; i + 1 < n; ++i)
                    if (r < cumulative_[i])
                        return i;
                return n - 1;
            }
            else
            {
                detail::big const r = detail::uniform_below_big(g, cumulative_big_.back());
                for (std::size_t i = 0; i + 1 < n; ++i)
                    if (r < cumulative_big_[i])
                        return i;
                return n - 1;
            }
        }

        std::tuple<L...> laws_;
        std::array<double, n> w_;
        detail::sd::double_double total_weight_{};
        std::array<std::uint64_t, n> cumulative_{};
        std::vector<detail::big> cumulative_big_;
    };

    template <class... L>
    mixture(weighted<L>...) -> mixture<L...>;
}
