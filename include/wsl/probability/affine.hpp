#pragma once

#include <wsl/modes.hpp>
#include <wsl/probability/detail/base.hpp>

// Algèbre des lois : image affine Y = aX + b d'une loi réelle (voir design/random.md §4.6).
//
//   wsl::probability::laws::normal N({.mean = 1, .stddev = 2});
//   auto Y = wsl::probability::affine(N, {.scale = 3, .shift = -1});   // wsl::probability::laws::normal<> : N(2, 6), paramètres calculés
//   wsl::probability::laws::student_t T({.degrees_of_freedom = 4});
//   auto Z = wsl::probability::affine(T, {.scale = -2, .shift = 5});   // wsl::probability::affine_image<student_t<>>
//   Z.cdf(1.0); Z.quantile(0.99); Z.mean(); Z(g);
//
// Identités : une famille fermée par image affine le déclare par
//   affine_closure(X, a, b) -> L               (normale, Cauchy, Laplace, logistique, uniforme),
// une famille fermée par changement d'échelle positif par
//   scale_closure(X, c) -> std::optional<L>    (exponentielle, gamma, Weibull, Pareto, log-normale) ;
// fonctions libres trouvées par ADL, que toute loi peut fournir. Avec affine_closure,
// a·X + b est une loi de la famille, dont les paramètres sont arrondis une fois (à
// chaque opération : 3 * N - 1 arrondit 3μ puis 3μ − 1 ; wsl::probability::affine(N, {3, -1}), une
// fois). Avec scale_closure, aX (a > 0, b = 0) reste de type affine_image<L> mais
// toutes ses fonctions sont celles de la loi de la famille : même précision.
//
// Sinon (affine_image<L>), composition : pdf_Y(y) = pdf_X(z)/|a|, z = (y − b)/a
// arrondi une fois ; cdf et ccdf échangées si a < 0 ; quantile a·Q_X(p) + b (ou a·cQ_X(p)
// + b), un arrondi ; moments transformés (a·m + b, a²v, ±asymétrie…). Bit-exact ;
// l'erreur ajoutée est celle de l'arrondi de z, amplifiée par le conditionnement
// |z f(z)/F(z)| de la répartition (mesurée : design/random.md).
//
// Tirage : a·x + b, x tiré selon X, arrondi une fois (approché).

#include <wsl/probability/laws/basic.hpp>
#include <wsl/probability/laws/cauchy.hpp>
#include <wsl/random/detail/arithmetic.hpp>
#include <wsl/probability/laws/exponential.hpp>
#include <wsl/probability/laws/gamma.hpp>
#include <wsl/probability/laws/laplace.hpp>
#include <wsl/probability/law.hpp>
#include <wsl/probability/laws/logistic.hpp>
#include <wsl/probability/laws/lognormal.hpp>
#include <wsl/probability/laws/normal.hpp>
#include <wsl/probability/laws/pareto.hpp>
#include <wsl/probability/laws/weibull.hpp>

#include <concepts>
#include <optional>
#include <random>
#include <stdexcept>
#include <type_traits>

namespace wsl::inline WSL_DEFAULTS::probability
{
    namespace detail
    {
        // a x + b, arrondi une fois (x, a quelconques ; infinis et NaN propagés).
        constexpr double linear(double a, double x, double b) noexcept
        {
            if (b == 0)
                return sd::fmul(a, x);
            if (a == 0 || x == 0)
                return b;
            if (!is_finite(a) || !is_finite(x))
                return sd::fmul(a, x) + b;
            bool const negative = (a < 0) != (x < 0);
            return affine_scaled(b, negative, normalized(magnitude(a), 0.0) * normalized(magnitude(x), 0.0));
        }

        // (y − b)/a, arrondi une fois.
        constexpr double unshift(double y, double b, double a) noexcept
        {
            if (!is_finite(y))
                return y / a;
            auto const [d, e] = sd::two_sum(y, -b);
            if (!is_finite(d))                                              // y − b déborde
                return y / a - b / a;
            if (d == 0)
                return 0.0;
            bool const negative = (d < 0) != (a < 0);
            double const r = rounded(normalized(magnitude(d), d < 0 ? -e : e) / magnitude(a));
            return negative ? -r : r;
        }
    }

    // Lois à valeurs réelles.
    template <class L>
    concept real_law = has_support<L> && std::floating_point<typename L::value_type>;

    template <class L>
    concept affine_closed = real_law<L> && requires(L const& X) {
        { affine_closure(X, 1.0, 0.0) } -> std::same_as<L>;
    };

    template <class L>
    concept scale_closed = real_law<L> && requires(L const& X) {
        { scale_closure(X, 1.0) } -> std::same_as<std::optional<L>>;
    };

    // Paramètres de l'image affine x ↦ scale·x + shift.
    struct affine_map
    {
        double scale = 1;
        double shift = 0;
    };

    // -----------------------------------------------------------------------
    // Identités : familles fermées.
    // -----------------------------------------------------------------------
    template <laws::law_real T>
    constexpr laws::normal<T> affine_closure(laws::normal<T> const& X, double a, double b)
    {
        // variance a² v (double-double : un arrondi à 2^-106 au plus), moyenne a μ + b arrondie une fois
        detail::scaled const a2 = detail::normalized(detail::magnitude(a), 0.0) * detail::normalized(detail::magnitude(a), 0.0);
        return laws::normal<T>(static_cast<T>(detail::linear(a, double(X.mean()), b)), a2 * X.exact_variance());
    }
    template <laws::law_real T>
    constexpr laws::cauchy<T> affine_closure(laws::cauchy<T> const& X, double a, double b)
    {
        return laws::cauchy<T>({.location = static_cast<T>(detail::linear(a, double(X.location()), b)), .scale = static_cast<T>(detail::sd::fmul(detail::magnitude(a), double(X.scale())))});
    }
    template <laws::law_real T>
    constexpr laws::laplace<T> affine_closure(laws::laplace<T> const& X, double a, double b)
    {
        return laws::laplace<T>({.location = static_cast<T>(detail::linear(a, double(X.location()), b)), .scale = static_cast<T>(detail::sd::fmul(detail::magnitude(a), double(X.scale())))});
    }
    template <laws::law_real T>
    constexpr laws::logistic<T> affine_closure(laws::logistic<T> const& X, double a, double b)
    {
        return laws::logistic<T>({.location = static_cast<T>(detail::linear(a, double(X.location()), b)), .scale = static_cast<T>(detail::sd::fmul(detail::magnitude(a), double(X.scale())))});
    }
    template <laws::law_real T>
    constexpr laws::uniform<T> affine_closure(laws::uniform<T> const& X, double a, double b)
    {
        auto const u = static_cast<T>(detail::linear(a, double(X.lower()), b)), v = static_cast<T>(detail::linear(a, double(X.upper()), b));
        return a > 0 ? laws::uniform<T>({.lower = u, .upper = v}) : laws::uniform<T>({.lower = v, .upper = u});
    }

    template <laws::law_real T>
    constexpr std::optional<laws::exponential<T>> scale_closure(laws::exponential<T> const& X, double c) noexcept
    {
        auto const l = X.exact_rate().divided(c);                           // λ/c exact
        if (!l)
            return std::nullopt;
        auto const [r, s] = detail::parameters_of(*l);
        return laws::exponential<T>({.rate = r, .scale = s});
    }
    template <laws::law_real T>
    constexpr std::optional<laws::weibull<T>> scale_closure(laws::weibull<T> const& X, double c) noexcept
    {
        auto const l = static_cast<T>(detail::sd::fmul(double(X.scale()), c));
        if (!(l > 0) || !detail::is_finite(l))
            return std::nullopt;
        return laws::weibull<T>({.scale = l, .shape = X.shape()});
    }
    template <laws::law_real T>
    constexpr std::optional<laws::pareto<T>> scale_closure(laws::pareto<T> const& X, double c) noexcept
    {
        auto const x = static_cast<T>(detail::sd::fmul(double(X.scale()), c));
        if (!(x > 0) || !detail::is_finite(x))
            return std::nullopt;
        return laws::pareto<T>({.scale = x, .shape = X.shape()});
    }
    // c e^Y = e^{Y + log c} : μ + log c, arrondi une fois.
    template <laws::law_real T>
    constexpr std::optional<laws::lognormal<T>> scale_closure(laws::lognormal<T> const& X, double c) noexcept
    {
        auto const m = detail::sum({double(X.log_mean()), 0.0}, detail::log_precise(c));
        auto const mu = static_cast<T>(m.hi + m.lo);
        if (!detail::is_finite(mu))
            return std::nullopt;
        return laws::lognormal<T>({.log_mean = mu, .log_stddev = X.log_stddev()});
    }

    // -----------------------------------------------------------------------
    // Image affine générique.
    // -----------------------------------------------------------------------
    template <real_law L>
    class affine_image : public law_interface
    {
    public:
        using value_type = typename L::value_type;
        using law_type = L;
        static constexpr exactness sampling = exactness::approximate;

        constexpr affine_image(L const& law, affine_map m) : law_(law), a_(m.scale), b_(m.shift)
        {
            if (!(a_ != 0) || !detail::is_finite(a_) || !detail::is_finite(b_))
                throw std::domain_error("wsl::probability::affine : il faut scale non nul et shift fini");
            if constexpr (scale_closed<L>)
                if (b_ == 0 && a_ > 0)
                    family_ = scale_closure(law_, a_);
        }

        constexpr L const& law() const noexcept { return law_; }
        constexpr double scale() const noexcept { return a_; }
        constexpr double shift() const noexcept { return b_; }

        constexpr interval<value_type> support() const noexcept
        {
            auto const [lo, up] = law_.support();
            auto const u = image(lo), v = image(up);
            return a_ > 0 ? interval<value_type>{u, v} : interval<value_type>{v, u};
        }

        constexpr value_type mean() const noexcept requires has_mean<L>
        {
            if constexpr (scale_closed<L>)
                if (family_)
                    return family_->mean();
            return image(law_.mean());
        }
        constexpr value_type median() const noexcept requires has_median<L>
        {
            if constexpr (scale_closed<L>)
                if (family_)
                    return family_->median();
            return image(law_.median());
        }
        constexpr value_type mode() const noexcept requires has_mode<L>
        {
            if constexpr (scale_closed<L>)
                if (family_)
                    return family_->mode();
            return image(law_.mode());
        }
        // a² v, arrondi une fois.
        constexpr value_type variance() const noexcept requires has_variance<L>
        {
            if constexpr (scale_closed<L>)
                if (family_)
                    return family_->variance();
            double const v = law_.variance();
            if (!(v > 0) || !detail::is_finite(v))
                return static_cast<value_type>(v);
            detail::scaled const s = detail::normalized(detail::magnitude(a_), 0.0);
            return static_cast<value_type>(detail::rounded(s * s * detail::normalized(v, 0.0)));
        }
        constexpr value_type skewness() const noexcept requires has_skewness<L>
        {
            auto const s = law_.skewness();
            return a_ > 0 || s == 0 ? s : -s;
        }
        constexpr value_type excess_kurtosis() const noexcept requires has_excess_kurtosis<L> { return law_.excess_kurtosis(); }
        // H + log|a|.
        constexpr value_type entropy() const noexcept requires has_entropy<L>
        {
            if constexpr (scale_closed<L>)
                if (family_)
                    return family_->entropy();
            auto const h = detail::sum({double(law_.entropy()), 0.0}, detail::log_precise(detail::magnitude(a_)));
            return static_cast<value_type>(h.hi + h.lo);
        }

        constexpr value_type pdf(value_type y) const noexcept requires has_pdf<L>
        {
            auto const& Y = *this;
            if constexpr (scale_closed<L>)
                if (Y.family_)
                    return wsl::probability::pdf(*Y.family_, y);
            double const f = wsl::probability::pdf(Y.law_, Y.preimage(y));
            return static_cast<value_type>(f / detail::magnitude(Y.a_));
        }
        constexpr value_type log_pdf(value_type y) const noexcept requires requires(L const& X, value_type x) { wsl::probability::log_pdf(X, x); }
        {
            auto const& Y = *this;
            if constexpr (scale_closed<L>)
                if (Y.family_)
                    return wsl::probability::log_pdf(*Y.family_, y);
            auto const la = detail::log_precise(detail::magnitude(Y.a_));
            auto const l = detail::sum({double(wsl::probability::log_pdf(Y.law_, Y.preimage(y))), 0.0}, {-la.hi, -la.lo});
            return static_cast<value_type>(l.hi + l.lo);
        }

        // P(aX + b <= y) = P(X <= z) si a > 0, P(X >= z) sinon (loi continue).
        constexpr value_type cdf(value_type y) const noexcept requires has_cdf<L> && has_ccdf<L>
        {
            auto const& Y = *this;
            if constexpr (scale_closed<L>)
                if (Y.family_)
                    return wsl::probability::cdf(*Y.family_, y);
            auto const z = Y.preimage(y);
            return Y.a_ > 0 ? wsl::probability::cdf(Y.law_, z) : wsl::probability::ccdf(Y.law_, z);
        }
        constexpr value_type ccdf(value_type y) const noexcept requires has_cdf<L> && has_ccdf<L>
        {
            auto const& Y = *this;
            if constexpr (scale_closed<L>)
                if (Y.family_)
                    return wsl::probability::ccdf(*Y.family_, y);
            auto const z = Y.preimage(y);
            return Y.a_ > 0 ? wsl::probability::ccdf(Y.law_, z) : wsl::probability::cdf(Y.law_, z);
        }

        constexpr value_type quantile(double p) const requires has_quantile<L> && has_cquantile<L>
        {
            auto const& Y = *this;
            if constexpr (scale_closed<L>)
                if (Y.family_)
                    return wsl::probability::quantile(*Y.family_, p);
            return Y.image(Y.a_ > 0 ? wsl::probability::quantile(Y.law_, p) : wsl::probability::cquantile(Y.law_, p));
        }
        constexpr value_type cquantile(double q) const requires has_quantile<L> && has_cquantile<L>
        {
            auto const& Y = *this;
            if constexpr (scale_closed<L>)
                if (Y.family_)
                    return wsl::probability::cquantile(*Y.family_, q);
            return Y.image(Y.a_ > 0 ? wsl::probability::cquantile(Y.law_, q) : wsl::probability::quantile(Y.law_, q));
        }

        template <std::uniform_random_bit_generator G>
            requires samplable<L, G>
        constexpr value_type sample(G& g) const
        {
            auto const& Y = *this;
            return Y.image(detail::native_sample(Y.law_, g));
        }
        template <std::uniform_random_bit_generator G>
            requires samplable<L, G>
        constexpr value_type sample(G& g, approximate_t) const
        {
            auto const& Y = *this;
            return Y.image(wsl::probability::sample(Y.law_, g, wsl::approximate));
        }

        friend constexpr bool operator==(affine_image const& p, affine_image const& q) noexcept
        {
            return p.law_ == q.law_ && p.a_ == q.a_ && p.b_ == q.b_;
        }

    private:
        constexpr value_type image(value_type x) const noexcept { return static_cast<value_type>(detail::linear(a_, double(x), b_)); }
        constexpr value_type preimage(value_type y) const noexcept { return static_cast<value_type>(detail::unshift(double(y), b_, a_)); }

        L law_;
        double a_, b_;
        std::optional<L> family_;                                           // aX dans la famille (scale_closure)
    };

    // a·X + b : loi de la famille si elle est fermée, image affine sinon ; une image
    // affine se compose (a'(aX + b) + b' = (a'a)X + (a'b + b'), coefficients arrondis une fois).
    template <real_law L>
    constexpr auto affine(L const& X, affine_map m)
    {
        if constexpr (affine_closed<L>)
        {
            if (!(m.scale != 0) || !detail::is_finite(m.scale) || !detail::is_finite(m.shift))
                throw std::domain_error("wsl::probability::affine : il faut scale non nul et shift fini");
            return affine_closure(X, m.scale, m.shift);
        }
        else
            return affine_image<L>(X, m);
    }
    template <real_law L>
    constexpr affine_image<L> affine(affine_image<L> const& Y, affine_map m)
    {
        if (!(m.scale != 0) || !detail::is_finite(m.scale) || !detail::is_finite(m.shift))
            throw std::domain_error("wsl::probability::affine : il faut scale non nul et shift fini");
        return affine_image<L>(Y.law(), {.scale = detail::sd::fmul(m.scale, Y.scale()), .shift = detail::linear(m.scale, Y.shift(), m.shift)});
    }
}
