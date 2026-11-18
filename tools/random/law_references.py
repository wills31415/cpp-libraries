#!/usr/bin/env python3
"""Valeurs de référence des lois de wsl::random (densités, répartitions, quantiles).

    uv run --no-project --with mpmath python tools/random/law_references.py FICHIER

Comme pour wsl::special : calculées avec mpmath à 256 bits, arrondies au plus proche en double ou en float, avec la position
exacte de la valeur entre deux flottants (en ulp). Non versionnées : CMake les
engendre dans le dossier de build. Déterministe (graine fixe).

Vérifie aussi les constantes en double-double de include/wsl/special/detail/extended.hpp.
"""

import math
import random
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/special"))
import generate as g                                         # noqa: E402  (arrondis, inverse d'erfc)
from mpmath import mp, mpf                                   # noqa: E402

mp.prec = 256
rng = random.Random(20261117)

HEADER = "// ENGENDRÉ par tools/random/law_references.py — ne pas modifier à la main.\n"


def check_constants():
    text = (ROOT / "include/wsl/special/detail/extended.hpp").read_text()
    exact = {
        "inv_sqrt2": 1 / mp.sqrt(2), "sqrt2": mp.sqrt(2), "inv_sqrt_2pi": 1 / mp.sqrt(2 * mp.pi),
        "half_log_2pi": mp.log(2 * mp.pi) / 2, "normal_entropy": (1 + mp.log(2 * mp.pi)) / 2,
        "half_sqrt_pi": mp.sqrt(mp.pi) / 2, "inv_sqrt_pi": 1 / mp.sqrt(mp.pi), "one_twelfth": mpf(1) / 12, "inv_pi": 1 / mp.pi, "log_4pi": mp.log(4 * mp.pi), "pi_squared_over_3": mp.pi ** 2 / 3,
    }
    for name, v in exact.items():
        hi = float(v)
        m = re.search(rf"\b{name}(?:_hi)? = ([-0-9a-fx.p+]+)", text)
        assert m and float.fromhex(m.group(1)) == hi, f"constante {name} fausse"
        m = re.search(rf"\b{name}_lo = ([-0-9a-fx.p+]+)", text)
        if m:
            assert float.fromhex(m.group(1)) == float(v - mpf(hi)), f"constante {name}_lo fausse"


def log_uniform(lo, hi):
    return 2.0 ** rng.uniform(math.log2(lo), math.log2(hi))


# ---------------------------------------------------------------------------
# Lois : paramètres, points, valeurs exactes
# ---------------------------------------------------------------------------

EXPONENTIAL_RATES = [1.0, 2.5, 1e-3, 1e5, 2.0 ** -1000, 2.0 ** 1000, 0.1]
NORMAL_PARAMETERS = [(0.0, 1.0), (3.0, 0.5), (-1e3, 7.0), (0.0, 1e-300), (1e10, 3e5), (0.25, 1e200), (-2.0, 0.3)]


def exponential_points(rate, count, to_type):
    xs = []
    for _ in range(count):
        t = log_uniform(2.0 ** -70, 1500.0)                  # λx
        xs.append(to_type(t / rate))
    xs += [to_type(1 / rate), to_type(math.log(2) / rate)]
    return [x for x in xs if x > 0 and math.isfinite(x)]


def normal_points(mu, sigma, count, to_type):
    xs = []
    for i in range(count):
        z = rng.uniform(-40, 40) if i % 2 else rng.choice([-1, 1]) * log_uniform(2.0 ** -60, 40.0)
        xs.append(to_type(mu + sigma * z))
    return [x for x in xs if math.isfinite(x)]


def probabilities(count):
    ps = [rng.random() for _ in range(count)]
    ps += [2.0 ** -rng.uniform(1, 1074) for _ in range(count)]
    ps += [1 - 2.0 ** -rng.uniform(1, 53) for _ in range(count // 4)]
    ps += [0.5 + rng.choice([-1, 1]) * 2.0 ** -rng.uniform(2, 54) for _ in range(count // 4)]
    return sorted({p for p in ps if 0 < p < 1})


def normal_quantile(p):
    return -mp.sqrt(2) * g.erfcinv_exact(2 * mpf(p))       # 2p exact : mpmath garde l'exposant à part


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


CHI2 = [1.0, 3.0, 10.0, 0.5, 100.0, 2500.0]
STUDENT = [1.0, 2.0, 5.0, 0.5, 30.0, 1e3, 3.7, 1e4]
LOGNORMAL = [(0.0, 1.0), (0.3, 0.5), (-5.0, 2.0), (10.0, 0.01), (0.0, 1e-3)]
CAUCHY = [(0.0, 1.0), (1.0, 2.0), (-1e5, 3.0), (0.0, 1e-10), (1e10, 1e3)]


def student_tail(nu, t):
    """P(T > |t|) = I_x(ν/2, ½)/2, x = ν/(ν + t²)."""
    with mp.workprec(450):
        nu, t = mpf(nu), mpf(t)
        x = nu / (nu + t * t)
        return g.ibeta_reference(nu / 2, mpf(1) / 2, x)[0] / 2


def derived_laws(kind):
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 120 if kind == "double" else 50
    out = []
    for k in CHI2:
        kk = mpf(to_type(k))
        xs = [to_type(float(k * log_uniform(2.0 ** -20, 8.0))) for _ in range(count)]
        out.append(("chi_squared_cdf", k, 0.0, xs, lambda x, kk=kk: g.igamma_reference(kk / 2, mpf(x) / 2)[0]))
        if k <= 100:                                                    # coût des références
            out.append(("chi_squared_quantile", k, 0.0, probabilities(count // 3), lambda u, kk=kk: 2 * g.igamma_inverse_exact(kk / 2, u, False)))
    for nu in STUDENT:
        n = mpf(to_type(nu))
        xs = [to_type(r) for r in [rng.choice([-1, 1]) * log_uniform(2.0 ** -30, 2.0 ** 30) for _ in range(count)]]
        c = mp.loggamma((n + 1) / 2) - mp.loggamma(n / 2) - mp.log(n * mp.pi) / 2
        logpdf = lambda x, n=n, c=c: c - (n + 1) / 2 * mp.log1p(mpf(x) ** 2 / n)
        out.append(("student_t_pdf", nu, 0.0, xs, lambda x, f=logpdf: mp.exp(f(x))))
        out.append(("student_t_log_pdf", nu, 0.0, xs, logpdf))
        out.append(("student_t_cdf", nu, 0.0, xs, lambda x, n=n: student_tail(n, x) if x < 0 else 1 - student_tail(n, x)))
        out.append(("student_t_ccdf", nu, 0.0, xs, lambda x, n=n: student_tail(n, x) if x > 0 else 1 - student_tail(n, x)))
        def quantile(u, n=n):
            if u == 0.5:
                return mpf(0)
            tail = min(mpf(u), 1 - mpf(u))
            with mp.workprec(450):
                x = g.ibeta_inverse_exact(n / 2, mpf(1) / 2, 2 * tail, False)
                t = mp.sqrt(n * (1 - x) / x)
            return -t if u < 0.5 else t
        if nu <= 100:                                                   # coût des références (inverse de la bêta)
            out.append(("student_t_quantile", nu, 0.0, probabilities(count // 3), quantile))
        h = n / 2
        out.append(("student_t_entropy", nu, 0.0, [0.0], lambda _, n=n, h=h: (h + mpf(1) / 2) * (mp.digamma(h + mpf(1) / 2) - mp.digamma(h))
                    + mp.log(mp.sqrt(n)) + (mp.loggamma(h) + mp.loggamma(mpf(1) / 2) - mp.loggamma(h + mpf(1) / 2))))
    for mu, sigma in LOGNORMAL:
        m, sg = mpf(to_type(mu)), mpf(to_type(sigma))
        xs = [to_type(float(mp.exp(m + sg * rng.uniform(-30, 30)))) for _ in range(count)]
        xs = [x for x in xs if 0 < x < math.inf]
        z = lambda x, m=m, sg=sg: (mp.log(mpf(x)) - m) / sg
        logpdf = lambda x, z=z, sg=sg: -z(x) ** 2 / 2 - mp.log(sg) - mp.log(mpf(x)) - mp.log(2 * mp.pi) / 2
        out.append(("lognormal_pdf", mu, sigma, xs, lambda x, f=logpdf: mp.exp(f(x))))
        out.append(("lognormal_log_pdf", mu, sigma, xs, logpdf))
        out.append(("lognormal_cdf", mu, sigma, xs, lambda x, z=z: mp.erfc(-z(x) / mp.sqrt(2)) / 2))
        out.append(("lognormal_ccdf", mu, sigma, xs, lambda x, z=z: mp.erfc(z(x) / mp.sqrt(2)) / 2))
        out.append(("lognormal_quantile", mu, sigma, probabilities(count // 3), lambda u, m=m, sg=sg: mp.exp(m + sg * normal_quantile(u))))
        out.append(("lognormal_mean", mu, sigma, [0.0], lambda _, m=m, sg=sg: mp.exp(m + sg ** 2 / 2)))
        out.append(("lognormal_variance", mu, sigma, [0.0], lambda _, m=m, sg=sg: mp.expm1(sg ** 2) * mp.exp(2 * m + sg ** 2)))
        out.append(("lognormal_entropy", mu, sigma, [0.0], lambda _, m=m, sg=sg: m + mp.log(sg) + (1 + mp.log(2 * mp.pi)) / 2))
    for x0, gm in CAUCHY:
        a, b = mpf(to_type(x0)), mpf(to_type(gm))
        xs = [to_type(float(a + b * rng.choice([-1, 1]) * log_uniform(2.0 ** -30, 2.0 ** 60))) for _ in range(count)]
        z = lambda x, a=a, b=b: (mpf(x) - a) / b
        out.append(("cauchy_pdf", x0, gm, xs, lambda x, z=z, b=b: 1 / (mp.pi * b * (1 + z(x) ** 2))))
        out.append(("cauchy_log_pdf", x0, gm, xs, lambda x, z=z, b=b: -mp.log(mp.pi * b * (1 + z(x) ** 2))))
        out.append(("cauchy_cdf", x0, gm, xs, lambda x, z=z: mp.atan2(1, -z(x)) / mp.pi))
        out.append(("cauchy_ccdf", x0, gm, xs, lambda x, z=z: mp.atan2(1, z(x)) / mp.pi))
        # −γ cot(πu) (u < ½), γ cot(π(1 − u)) : tan(π(u − ½)) perdrait u minuscule dans u − ½.
        out.append(("cauchy_quantile", x0, gm, probabilities(count // 3),
                    lambda u, a=a, b=b: a - b / mp.tan(mp.pi * mpf(u)) if u < 0.5 else a + b / mp.tan(mp.pi * (1 - mpf(u)))))
        out.append(("cauchy_entropy", x0, gm, [0.0], lambda _, b=b: mp.log(4 * mp.pi * b)))
    return out


FISHER = [(5.0, 10.0), (1.0, 1.0), (2.0, 7.0), (30.0, 3.0), (0.5, 40.0), (100.0, 200.0)]
WEIBULL = [(1.5, 2.0), (0.5, 1.0), (1.0, 3.0), (5.0, 0.1), (50.0, 1.0), (0.1, 1e3), (300.0, 2.0)]
LAPLACE = [(0.0, 1.0), (1.0, 2.0), (-1e5, 3.0), (0.0, 1e-10)]
LOGISTIC = [(0.0, 1.0), (1.0, 2.0), (-1e5, 3.0), (0.0, 1e-10)]
PARETO = [(1.0, 2.5), (3.0, 0.5), (1e-5, 10.0), (2.0, 1.0), (1.0, 300.0)]


def second_batch(kind):
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 100 if kind == "double" else 40
    out = []
    for d1, d2 in FISHER:
        a1, a2 = mpf(to_type(d1)), mpf(to_type(d2))
        xs = [to_type(float(log_uniform(2.0 ** -30, 2.0 ** 30))) for _ in range(count)]
        u = lambda x, a1=a1, a2=a2: (a1 * mpf(x), a2)
        logpdf = lambda x, a1=a1, a2=a2: ((a1 / 2) * mp.log(a1 * mpf(x) / (a1 * mpf(x) + a2)) + (a2 / 2) * mp.log(a2 / (a1 * mpf(x) + a2))
                                          - (mp.loggamma(a1 / 2) + mp.loggamma(a2 / 2) - mp.loggamma((a1 + a2) / 2)) - mp.log(mpf(x)))
        out.append(("fisher_f_pdf", d1, d2, xs, lambda x, f=logpdf: mp.exp(f(x))))
        out.append(("fisher_f_log_pdf", d1, d2, xs, logpdf))
        def ratio(x, a1=a1, a2=a2):
            with mp.workprec(450):
                return a1 * mpf(x) / (a1 * mpf(x) + a2)
        out.append(("fisher_f_cdf", d1, d2, xs, lambda x, a1=a1, a2=a2, ratio=ratio: g.ibeta_reference(a1 / 2, a2 / 2, ratio(x))[0]))
        out.append(("fisher_f_ccdf", d1, d2, xs, lambda x, a1=a1, a2=a2, ratio=ratio: g.ibeta_reference(a1 / 2, a2 / 2, ratio(x))[1]))
        def quantile(p, a1=a1, a2=a2):
            with mp.workprec(450):
                u_ = g.ibeta_inverse_exact(a1 / 2, a2 / 2, p, False)
                return a2 * u_ / (a1 * (1 - u_))
        out.append(("fisher_f_quantile", d1, d2, probabilities(count // 4), quantile))
        a, b = a1 / 2, a2 / 2
        out.append(("fisher_f_entropy", d1, d2, [0.0], lambda _, a=a, b=b: mp.log(b / a) + (mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b))
                    + (1 - a) * mp.digamma(a) - (1 + b) * mp.digamma(b) + (a + b) * mp.digamma(a + b)))
    for k, lam in WEIBULL:
        kk, ll = mpf(to_type(k)), mpf(to_type(lam))
        xs = [to_type(float(ll * log_uniform(2.0 ** -30, 2.0 ** 3) ** (1 / kk))) for _ in range(count)]
        xs = [x for x in xs if 0 < x < math.inf]
        t = lambda x, kk=kk, ll=ll: (mpf(x) / ll) ** kk
        out.append(("weibull_pdf", k, lam, xs, lambda x, t=t, kk=kk: kk * t(x) * mp.exp(-t(x)) / mpf(x)))
        out.append(("weibull_log_pdf", k, lam, xs, lambda x, t=t, kk=kk: mp.log(kk) + mp.log(t(x)) - t(x) - mp.log(mpf(x))))
        out.append(("weibull_cdf", k, lam, xs, lambda x, t=t: -mp.expm1(-t(x))))
        out.append(("weibull_ccdf", k, lam, xs, lambda x, t=t: mp.exp(-t(x))))
        out.append(("weibull_quantile", k, lam, probabilities(count // 3), lambda p, kk=kk, ll=ll: ll * (-mp.log1p(-mpf(p))) ** (1 / kk)))
        out.append(("weibull_mean", k, lam, [0.0], lambda _, kk=kk, ll=ll: ll * mp.gamma(1 + 1 / kk)))
        out.append(("weibull_variance", k, lam, [0.0], lambda _, kk=kk, ll=ll: ll ** 2 * (mp.gamma(1 + 2 / kk) - mp.gamma(1 + 1 / kk) ** 2)))
        out.append(("weibull_entropy", k, lam, [0.0], lambda _, kk=kk, ll=ll: mp.euler * (1 - 1 / kk) + mp.log(ll / kk) + 1))
    for name, params, tail in (("laplace", LAPLACE, lambda z: mp.exp(-abs(z)) / 2), ("logistic", LOGISTIC, lambda z: mp.exp(-abs(z)) / (1 + mp.exp(-abs(z))))):
        for mu, sc in params:
            m, s_ = mpf(to_type(mu)), mpf(to_type(sc))
            xs = [to_type(float(m + s_ * rng.choice([-1, 1]) * log_uniform(2.0 ** -30, 700.0))) for _ in range(count)]
            z = lambda x, m=m, s_=s_: (mpf(x) - m) / s_
            if name == "laplace":
                pdf = lambda x, z=z, s_=s_: mp.exp(-abs(z(x))) / (2 * s_)
                q = lambda p, m=m, s_=s_: m + s_ * mp.log(2 * mpf(p)) if p < 0.5 else m - s_ * mp.log(2 * (1 - mpf(p)))
            else:
                pdf = lambda x, z=z, s_=s_: mp.exp(-abs(z(x))) / (s_ * (1 + mp.exp(-abs(z(x)))) ** 2)
                q = lambda p, m=m, s_=s_: m + s_ * (mp.log(mpf(p)) - mp.log1p(-mpf(p)))
            out.append((f"{name}_pdf", mu, sc, xs, pdf))
            out.append((f"{name}_log_pdf", mu, sc, xs, lambda x, pdf=pdf: mp.log(pdf(x))))
            out.append((f"{name}_cdf", mu, sc, xs, lambda x, z=z, tail=tail: tail(z(x)) if z(x) < 0 else 1 - tail(z(x))))
            out.append((f"{name}_ccdf", mu, sc, xs, lambda x, z=z, tail=tail: tail(z(x)) if z(x) > 0 else 1 - tail(z(x))))
            out.append((f"{name}_quantile", mu, sc, [p for p in probabilities(count // 3) if p != 0.5], q))
    for xm, al in PARETO:
        x0, a = mpf(to_type(xm)), mpf(to_type(al))
        xs = [to_type(float(x0 * log_uniform(1.0 + 2.0 ** -40, 2.0 ** 40))) for _ in range(count)]
        xs = [x for x in xs if x > float(x0) and x < math.inf]
        out.append(("pareto_pdf", xm, al, xs, lambda x, x0=x0, a=a: a / mpf(x) * (x0 / mpf(x)) ** a))
        out.append(("pareto_log_pdf", xm, al, xs, lambda x, x0=x0, a=a: mp.log(a / mpf(x)) + a * mp.log(x0 / mpf(x))))
        out.append(("pareto_cdf", xm, al, xs, lambda x, x0=x0, a=a: -mp.expm1(a * mp.log(x0 / mpf(x)))))
        out.append(("pareto_ccdf", xm, al, xs, lambda x, x0=x0, a=a: (x0 / mpf(x)) ** a))
        out.append(("pareto_quantile", xm, al, probabilities(count // 3), lambda p, x0=x0, a=a: x0 * mp.exp(-mp.log1p(-mpf(p)) / a)))
        out.append(("pareto_entropy", xm, al, [0.0], lambda _, x0=x0, a=a: mp.log(x0 / a) + 1 + 1 / a))
    return out


def discrete_quantiles(name, cases, cdf, start):
    """Lignes {a, b, u, k} : k = plus petit entier tel que F(k) >= u (F exacte), sauf ambiguïté."""
    tasks = [(a, b, u) for a, b in cases for u in probabilities(40)]       # tirages dans l'ordre, avant le calcul
    def solve(a, b, u):
        k = max(0, int(start(a, b, u)))
        while k > 0 and cdf(a, b, k - 1) >= u:
            k -= 1
        while cdf(a, b, k) < u:
            k += 1
        near = lambda j: j >= 0 and abs(cdf(a, b, j) - u) < mpf(2) ** -90 * u
        return mpf(-1) if near(k) or near(k - 1) or k >= 2 ** 62 else mpf(k)
    rows = []
    for (a, b, u), k in zip(tasks, g.parallel_map(solve, tasks)):
        if k >= 0:
            rows.append(f"        {{{g.hexd(float(a))}, {g.hexd(float(b))}, {g.hexd(u)}, {int(k)}}},\n")
    return rows


def reference_file():
    out = [HEADER, "#pragma once\n\n", "// Valeurs de référence des lois de wsl::random (mpmath).\n\n",
           "namespace law_reference_data\n{\n",
           "    template <class T>\n    struct point\n    {\n        double a, b, x;\n        T rounded;\n        double offset;\n    };\n\n"]
    total = 0
    for kind, rounder, ulp, fmt in (("double", g.round_double, g.ulp_double, g.hexd),
                                    ("float", g.round_float, g.ulp_float, g.hexf)):
        groups = {}
        for name, a, b, xs, f in laws(kind) + (discrete_laws() if kind == "double" else []) + phase4_laws(kind) + derived_laws(kind) + second_batch(kind):
            rows = groups.setdefault(name, [])
            for x, exact in zip(xs, g.parallel_map(lambda v, f=f: f(mpf(v)), [(x,) for x in xs])):
                cr = rounder(exact)
                if math.isinf(cr) or cr == 0:
                    continue
                offset = float((exact - mpf(cr)) / mpf(ulp(cr)))
                rows.append(f"        {{{g.hexd(float(a))}, {g.hexd(float(b))}, {g.hexd(float(x))}, {fmt(cr)}, {offset!r}}},\n")
        for name, rows in groups.items():
            out.append(f"    inline constexpr point<{kind}> {name}_{kind}[] = {{\n")
            out += rows
            out.append("    };\n\n")
            total += len(rows)
    out.append("    struct quantile_point\n    {\n        double p, u;\n        long long k;\n    };\n\n")
    out.append("    // Quantiles exacts de la loi géométrique.\n    inline constexpr quantile_point geometric_quantile[] = {\n")
    rows = geometric_quantiles()
    out += rows
    out.append("    };\n\n")
    total += len(rows)
    out.append("    struct quantile_point2\n    {\n        double a, b, u;\n        long long k;\n    };\n\n")
    import math as _m
    poisson_cdf = lambda mu, _, k: g.igamma_reference(float(k) + 1, mu)[1]
    poisson_start = lambda mu, _, u: mu + _m.sqrt(mu) * float(-mp.sqrt(2) * g.erfcinv_exact(2 * mpf(u)))
    binomial_cdf = lambda n, pp, k: mpf(1) if k >= n else g.ibeta_reference(float(k) + 1, float(n - k), pp)[1]
    binomial_start = lambda n, pp, u: min(n, max(0.0, n * pp + _m.sqrt(n * pp * (1 - pp)) * float(-mp.sqrt(2) * g.erfcinv_exact(2 * mpf(u)))))
    for name, cases, cdf, start in (("poisson_quantile", [(mu, 0.0) for mu in P4_POISSON], poisson_cdf, poisson_start),
                                    ("binomial_quantile", [(float(n), pp) for n, pp in P4_BINOMIAL], binomial_cdf, binomial_start)):
        rows = discrete_quantiles(name, cases, cdf, start)
        out.append(f"    // Quantiles exacts.\n    inline constexpr quantile_point2 {name}[] = {{\n")
        out += rows
        out.append("    };\n\n")
        total += len(rows)
    out.append("}\n")
    return "".join(out), total


def main():
    if len(sys.argv) != 2:
        sys.exit("usage : law_references.py FICHIER")
    check_constants()
    text, total = reference_file()
    output = Path(sys.argv[1])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text)
    print(f"{total} points de référence des lois : {output}")


if __name__ == "__main__":
    main()
