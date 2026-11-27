"""Exponentielle, normale, gamma, bêta, géométrique, Poisson, binomiale ; constantes de extended.hpp."""

from common import *                                         # noqa: F401,F403


def laws(kind):
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    out = []
    for rate in EXPONENTIAL_RATES:
        if kind == "float" and not (2.0 ** -120 < rate < 2.0 ** 120):
            continue
        r = mpf(to_type(rate))
        xs = exponential_points(rate, 250 if kind == "double" else 100, to_type)
        out.append(("exponential_pdf", rate, 0.0, xs, lambda x, r=r: r * mp.exp(-r * x)))
        out.append(("exponential_log_pdf", rate, 0.0, xs, lambda x, r=r: mp.log(r) - r * x))
        out.append(("exponential_cdf", rate, 0.0, xs, lambda x, r=r: -mp.expm1(-r * x)))
        out.append(("exponential_ccdf", rate, 0.0, xs, lambda x, r=r: mp.exp(-r * x)))
        out.append(("exponential_quantile", rate, 0.0, probabilities(100 if kind == "double" else 40),
                    lambda p, r=r: -mp.log1p(-p) / r))
    for mu, sigma in NORMAL_PARAMETERS:
        if kind == "float" and not (2.0 ** -120 < sigma < 2.0 ** 120 and abs(mu) < 2.0 ** 120):
            continue
        m, s = mpf(to_type(mu)), mpf(to_type(sigma))
        xs = normal_points(mu, sigma, 250 if kind == "double" else 100, to_type)
        z = lambda x, m=m, s=s: (x - m) / s
        out.append(("normal_pdf", mu, sigma, xs, lambda x, z=z, s=s: mp.exp(-z(x) ** 2 / 2) / (s * mp.sqrt(2 * mp.pi))))
        out.append(("normal_log_pdf", mu, sigma, xs, lambda x, z=z, s=s: -z(x) ** 2 / 2 - mp.log(s) - mp.log(2 * mp.pi) / 2))
        out.append(("normal_cdf", mu, sigma, xs, lambda x, z=z: mp.erfc(-z(x) / mp.sqrt(2)) / 2))
        out.append(("normal_ccdf", mu, sigma, xs, lambda x, z=z: mp.erfc(z(x) / mp.sqrt(2)) / 2))
        out.append(("normal_quantile", mu, sigma, probabilities(100 if kind == "double" else 40),
                    lambda p, m=m, s=s: m + s * normal_quantile(p)))
    for shape, rate, by_scale in GAMMA:
        if kind == "float" and not all(2.0 ** -120 < v < 2.0 ** 120 for v in (shape, rate)):
            continue
        a, r = mpf(to_type(shape)), mpf(to_type(rate))
        lam = 1 / r if by_scale else r
        xs = []
        for _ in range(250 if kind == "double" else 100):
            t = shape * log_uniform(2.0 ** -40, 40.0) if rng.random() < 0.8 else log_uniform(2.0 ** -60, 2.0)
            xs.append(to_type(float(t / lam)))
        xs = [x for x in xs if 0 < x < math.inf]
        logpdf = lambda x, a=a, lam=lam: a * mp.log(lam) + (a - 1) * mp.log(x) - lam * x - mp.loggamma(a)
        name = "gamma_scale" if by_scale else "gamma"
        out.append((f"{name}_pdf", shape, rate, xs, lambda x, f=logpdf: mp.exp(f(x))))
        out.append((f"{name}_log_pdf", shape, rate, xs, logpdf))
    for alpha, beta in BETA:
        if kind == "float" and not all(2.0 ** -120 < v < 2.0 ** 120 for v in (alpha, beta)):
            continue
        a, b = mpf(to_type(alpha)), mpf(to_type(beta))
        xs = []
        for _ in range(250 if kind == "double" else 100):
            u = rng.random()
            x = rng.random() if u < 0.4 else log_uniform(2.0 ** -60, 0.5) if u < 0.7 else 1 - log_uniform(2.0 ** -53, 0.5)
            xs.append(to_type(x))
        m = alpha / (alpha + beta)
        xs += [to_type(m * f) for f in (0.5, 0.9, 1.0, 1.1, 2.0)]
        xs = [x for x in xs if 0 < x < 1]
        lb = mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)
        logpdf = lambda x, a=a, b=b, lb=lb: (a - 1) * mp.log(x) + (b - 1) * mp.log1p(-x) - lb
        out.append(("beta_pdf", alpha, beta, xs, lambda x, f=logpdf: mp.exp(f(x))))
        out.append(("beta_log_pdf", alpha, beta, xs, logpdf))
    return out


GAMMA = [(2.5, 2.0, False), (0.3, 1.0, False), (1.0, 3.0, False), (1e6, 1.0, False), (0.999, 1e3, False),
         (1e-3, 5.0, True), (50.0, 1e-200, True), (7.25, 1e100, True), (3.0, 0.5, True), (1e-300, 1.0, False)]
BETA = [(2.0, 5.0), (0.5, 0.5), (1.5, 1e8), (30.0, 70.0), (1e-3, 2.0), (2.5, 1e-2), (1e5, 3e5), (1.0, 1.0),
        (1.0, 3.0), (0.2, 50.0), (3.0, 1e12), (2.0, 2.0), (12.0, 0.7)]


GEOMETRIC_P = [0.2, 0.5, 1e-3, 1e-12, 0.999, 2.0 ** -60, 0.75]
POISSON_MEANS = [0.5, 3.5, 25.0, 1e3, 1e6, 1e-5, 1e12, 9.99, 10.0]
BINOMIAL = [(20, 0.3), (1000, 0.5), (10 ** 6, 1e-3), (50, 0.999), (10 ** 9, 0.37), (3, 1e-10), (7, 0.5)]


def integer_points(center, spread, low, high, count):
    ks = {low, high, int(center)}
    for _ in range(count):
        k = int(round(center + spread * rng.uniform(-40, 40) * rng.random() ** 3))
        ks.add(min(max(k, low), high))
    return sorted(ks)


def discrete_laws():
    out = []
    for p in GEOMETRIC_P:
        pp = mpf(p)
        mean = (1 - p) / p
        ks = sorted({min(int(mean * log_uniform(2.0 ** -20, 600.0)), 2 ** 62) for _ in range(250)} | {0, 1})
        lq = mp.log1p(-pp)
        out.append(("geometric_pmf", p, 0.0, ks, lambda k, pp=pp: pp * (1 - pp) ** k))
        out.append(("geometric_log_pmf", p, 0.0, ks, lambda k, pp=pp, lq=lq: mp.log(pp) + k * lq))
        out.append(("geometric_cdf", p, 0.0, ks, lambda k, lq=lq: -mp.expm1((k + 1) * lq)))
        out.append(("geometric_ccdf", p, 0.0, ks, lambda k, lq=lq: mp.exp((k + 1) * lq)))
    for mu in POISSON_MEANS:
        m = mpf(mu)
        ks = integer_points(mu, math.sqrt(mu) + 1, 0, 2 ** 62, 250)
        ks = [k for k in ks if k >= 0]
        logpmf = lambda k, m=m: -m + k * mp.log(m) - mp.loggamma(k + 1)
        out.append(("poisson_pmf", mu, 0.0, ks, lambda k, f=logpmf: mp.exp(f(k))))
        out.append(("poisson_log_pmf", mu, 0.0, ks, logpmf))
    for n, p in BINOMIAL:
        pp = mpf(p)
        sd = math.sqrt(n * p * (1 - p)) + 1
        ks = list(range(n + 1)) if n <= 60 else integer_points(n * p, sd, 0, n, 250)
        logpmf = lambda k, n=n, pp=pp: (mp.loggamma(n + 1) - mp.loggamma(k + 1) - mp.loggamma(n - k + 1)
                                        + k * mp.log(pp) + (n - k) * mp.log1p(-pp))
        out.append(("binomial_pmf", float(n), p, ks, lambda k, f=logpmf: mp.exp(f(k))))
        out.append(("binomial_log_pmf", float(n), p, ks, logpmf))
    return out


def geometric_quantiles():
    rows = []
    for p in GEOMETRIC_P:
        lq = mp.log1p(-mpf(p))
        for u in probabilities(60):
            t = mp.log1p(-mpf(u)) / lq
            k = int(mp.ceil(t)) - 1
            if abs(t - mp.nint(t)) < mpf(2) ** -90 * t or k >= 2 ** 62:
                continue
            rows.append(f"        {{{g.hexd(p)}, {g.hexd(u)}, {k}}},\n")
    return rows


# Phase 4 : répartitions, survies, quantiles (gamma et bêta incomplètes), entropies.
P4_POISSON = [0.5, 3.5, 25.0, 1e3, 1e6, 1e-5, 9.99]
P4_BINOMIAL = [(20, 0.3), (1000, 0.5), (10 ** 6, 1e-3), (50, 0.999), (3, 1e-10), (7, 0.5)]


def phase4_laws(kind):
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    out = []
    count = 120 if kind == "double" else 50
    for shape, rate, by_scale in GAMMA:
        if kind == "float" and not all(2.0 ** -120 < v < 2.0 ** 120 for v in (shape, rate)):
            continue
        a, r = to_type(shape), to_type(rate)
        lam = 1 / mpf(r) if by_scale else mpf(r)
        xs = [to_type(float(shape * log_uniform(2.0 ** -20, 8.0) / lam)) for _ in range(count)]
        xs = [x for x in xs if 0 < x < math.inf]
        t = lambda x, lam=lam, r=r, by_scale=by_scale: (mpf(x) / mpf(r)) if by_scale else mpf(x) * lam
        name = "gamma_scale" if by_scale else "gamma"
        out.append((f"{name}_cdf", shape, rate, xs, lambda x, a=a, t=t: g.igamma_reference(a, t(x))[0]))
        out.append((f"{name}_ccdf", shape, rate, xs, lambda x, a=a, t=t: g.igamma_reference(a, t(x))[1]))
        ps = probabilities(count // 3)
        if shape <= 1e4:                                                # coût de la référence (inverse à 320 bits)
            out.append((f"{name}_quantile", shape, rate, ps,
                        lambda u, a=a, r=r, by_scale=by_scale: g.igamma_inverse_exact(a, u, False) * (mpf(r) if by_scale else 1 / mpf(r))))
    for alpha, beta in BETA:
        if kind == "float" and not all(2.0 ** -120 < v < 2.0 ** 120 for v in (alpha, beta)):
            continue
        a, b = to_type(alpha), to_type(beta)
        xs = []
        for _ in range(count):
            u = rng.random()
            xs.append(to_type(rng.random() if u < 0.4 else log_uniform(2.0 ** -40, 0.5) if u < 0.7 else 1 - log_uniform(2.0 ** -40, 0.5)))
        xs = [x for x in xs if 0 < x < 1]
        out.append(("beta_cdf", alpha, beta, xs, lambda x, a=a, b=b: g.ibeta_reference(a, b, x)[0]))
        out.append(("beta_ccdf", alpha, beta, xs, lambda x, a=a, b=b: g.ibeta_reference(a, b, x)[1]))
        if max(alpha, beta) <= 2.0 ** 12:
            out.append(("beta_quantile", alpha, beta, probabilities(count // 3),
                        lambda u, a=a, b=b: g.ibeta_inverse_exact(a, b, u, False)))
    # Entropies : une valeur par loi (x sans objet).
    for shape, rate, by_scale in GAMMA:
        if kind == "float" and not all(2.0 ** -120 < v < 2.0 ** 120 for v in (shape, rate)):
            continue
        a, r = mpf(to_type(shape)), mpf(to_type(rate))
        lam = 1 / r if by_scale else r
        name = "gamma_scale" if by_scale else "gamma"
        out.append((f"{name}_entropy", shape, rate, [0.0],
                    lambda _, a=a, lam=lam: a - mp.log(lam) + mp.loggamma(a) + (1 - a) * mp.digamma(a)))
    for alpha, beta in BETA:
        if kind == "float" and not all(2.0 ** -120 < v < 2.0 ** 120 for v in (alpha, beta)):
            continue
        a, b = mpf(to_type(alpha)), mpf(to_type(beta))
        out.append(("beta_entropy", alpha, beta, [0.0],
                    lambda _, a=a, b=b: (mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)) - (a - 1) * mp.digamma(a)
                    - (b - 1) * mp.digamma(b) + (a + b - 2) * mp.digamma(a + b)))
    if kind == "double":
        for mu in P4_POISSON:
            ks = [k for k in integer_points(mu, math.sqrt(mu) + 1, 0, 2 ** 62, 80) if k >= 0]
            out.append(("poisson_cdf", mu, 0.0, ks, lambda k, mu=mu: g.igamma_reference(float(k) + 1, mu)[1]))
            out.append(("poisson_ccdf", mu, 0.0, ks, lambda k, mu=mu: g.igamma_reference(float(k) + 1, mu)[0]))
        for n, pp in P4_BINOMIAL:
            sd = math.sqrt(n * pp * (1 - pp)) + 1
            ks = [k for k in (list(range(n)) if n <= 60 else integer_points(n * pp, sd, 0, n - 1, 80)) if k < n]
            out.append(("binomial_cdf", float(n), pp, ks, lambda k, n=n, pp=pp: g.ibeta_reference(float(k) + 1, float(n - k), pp)[1]))
            out.append(("binomial_ccdf", float(n), pp, ks, lambda k, n=n, pp=pp: g.ibeta_reference(float(k) + 1, float(n - k), pp)[0]))
    return out


def entries(kind):
    return laws(kind) + (discrete_laws() if kind == "double" else []) + phase4_laws(kind)


def tables():
    poisson_cdf = lambda mu, _, k: g.igamma_reference(float(k) + 1, mu)[1]
    poisson_start = lambda mu, _, u: mu + math.sqrt(mu) * normal_start(u)
    binomial_cdf = lambda n, pp, k: mpf(1) if k >= n else g.ibeta_reference(float(k) + 1, float(n - k), pp)[1]
    binomial_start = lambda n, pp, u: min(n, max(0.0, n * pp + math.sqrt(n * pp * (1 - pp)) * normal_start(u)))
    return [("quantile_point", "geometric_quantile", geometric_quantiles()),
            ("quantile_point2", "poisson_quantile", discrete_quantiles("poisson_quantile", [(mu, 0.0) for mu in P4_POISSON], poisson_cdf, poisson_start)),
            ("quantile_point2", "binomial_quantile", discrete_quantiles("binomial_quantile", [(float(n), pp) for n, pp in P4_BINOMIAL], binomial_cdf, binomial_start))]


def main():
    run('base', entries, tables, check=check_constants)
