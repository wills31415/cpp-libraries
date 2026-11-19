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
        "half_sqrt_pi": mp.sqrt(mp.pi) / 2, "inv_sqrt_pi": 1 / mp.sqrt(mp.pi), "one_twelfth": mpf(1) / 12, "inv_pi": 1 / mp.pi, "log_4pi": mp.log(4 * mp.pi), "pi_squared_over_3": mp.pi ** 2 / 3, "sqrt_half_pi": mp.sqrt(mp.pi / 2), "four_minus_pi_half": (4 - mp.pi) / 2, "sqrt_two_ln2": mp.sqrt(2 * mp.log(2)), "gumbel_skewness": 12 * mp.sqrt(6) * mp.zeta(3) / mp.pi ** 3, "rayleigh_skewness": 2 * mp.sqrt(mp.pi) * (mp.pi - 3) / (4 - mp.pi) ** mpf(1.5), "rayleigh_kurtosis": -(6 * mp.pi ** 2 - 24 * mp.pi + 16) / (4 - mp.pi) ** 2,
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
CAUCHY = [(0.0, 1.0), (1.0, 2.0), (-1e5, 3.0), (0.0, 1e-10), (1e10, 1e3), (0.0, 2.0 ** -900)]


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
        if kind == "float" and not 2.0 ** -120 < gm < 2.0 ** 120:
            continue
        a, b = mpf(to_type(x0)), mpf(to_type(gm))
        top = 2.0 ** 1000 if kind == "double" else 2.0 ** 90                   # |z| > 2^600 compris (double)
        xs = [to_type(float(a + b * rng.choice([-1, 1]) * log_uniform(2.0 ** -30, top))) for _ in range(count)]
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
            far = (lambda: log_uniform(2.0 ** 600, 2.0 ** 900)) if kind == "double" else (lambda: log_uniform(2.0 ** 60, 2.0 ** 100))
            xs = [to_type(float(m + s_ * rng.choice([-1, 1]) * (log_uniform(2.0 ** -30, 700.0) if rng.random() < 0.9 else far()))) for _ in range(count)]
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


def survival_quantiles(kind):
    """Quantiles de survie : x tel que P(X > x) = q, écrits sur la queue (q exact, jamais 1 − q
    arrondi)."""
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 40 if kind == "double" else 16
    fits = lambda *v: kind == "double" or all(2.0 ** -120 < abs(x) < 2.0 ** 120 or x == 0 for x in v)
    out = []
    qs = lambda: probabilities(count)
    for rate in EXPONENTIAL_RATES:
        if fits(rate):
            r = mpf(to_type(rate))
            out.append(("exponential_cquantile", rate, 0.0, qs(), lambda q, r=r: -mp.log(mpf(q)) / r))
    for mu, sigma in NORMAL_PARAMETERS:
        if fits(mu, sigma):
            m, sg = mpf(to_type(mu)), mpf(to_type(sigma))
            out.append(("normal_cquantile", mu, sigma, qs(), lambda q, m=m, sg=sg: m - sg * normal_quantile(q)))
    for shape, rate, by_scale in GAMMA:
        if fits(shape, rate) and 1e-100 < shape <= 1e4:              # forme 1e-300 : solution hors de toute échelle
            a, r = to_type(shape), to_type(rate)
            out.append(("gamma_scale_cquantile" if by_scale else "gamma_cquantile", shape, rate, probabilities(count // 2),
                        lambda q, a=a, r=r, by_scale=by_scale: g.igamma_inverse_exact(a, q, True) * (mpf(r) if by_scale else 1 / mpf(r))))
    for alpha, beta in BETA:
        if fits(alpha, beta) and max(alpha, beta) <= 2.0 ** 12:
            a, b = to_type(alpha), to_type(beta)
            out.append(("beta_cquantile", alpha, beta, probabilities(count // 2), lambda q, a=a, b=b: g.ibeta_inverse_exact(a, b, q, True)))
    for k in CHI2:
        if k <= 100:
            kk = mpf(to_type(k))
            out.append(("chi_squared_cquantile", k, 0.0, probabilities(count // 2), lambda q, kk=kk: 2 * g.igamma_inverse_exact(kk / 2, q, True)))
    for nu in STUDENT:
        if nu <= 100:
            n = mpf(to_type(nu))
            def survival(q, n=n):
                if q == 0.5:
                    return mpf(0)
                tail = min(mpf(q), 1 - mpf(q))
                with mp.workprec(450):
                    x = g.ibeta_inverse_exact(n / 2, mpf(1) / 2, 2 * tail, False)
                    t = mp.sqrt(n * (1 - x) / x)
                return t if q < 0.5 else -t
            out.append(("student_t_cquantile", nu, 0.0, probabilities(count // 2), survival))
    for mu, sigma in LOGNORMAL:
        m, sg = mpf(to_type(mu)), mpf(to_type(sigma))
        out.append(("lognormal_cquantile", mu, sigma, qs(), lambda q, m=m, sg=sg: mp.exp(m - sg * normal_quantile(q))))
    for x0, gm in CAUCHY:
        if fits(x0, gm):
            a, b = mpf(to_type(x0)), mpf(to_type(gm))
            out.append(("cauchy_cquantile", x0, gm, qs(),
                        lambda q, a=a, b=b: a + b / mp.tan(mp.pi * mpf(q)) if q < 0.5 else a - b / mp.tan(mp.pi * (1 - mpf(q)))))
    for d1, d2 in FISHER:
        a1, a2 = mpf(to_type(d1)), mpf(to_type(d2))
        def survival(q, a1=a1, a2=a2):
            with mp.workprec(450):
                w = g.ibeta_inverse_exact(a2 / 2, a1 / 2, q, False)     # 1 − u, I_w(b, a) = q (1 − u minuscule)
                return a2 * (1 - w) / (a1 * w)
        out.append(("fisher_f_cquantile", d1, d2, probabilities(count // 2), survival))
    for k, lam in WEIBULL:
        kk, ll = mpf(to_type(k)), mpf(to_type(lam))
        out.append(("weibull_cquantile", k, lam, qs(), lambda q, kk=kk, ll=ll: ll * (-mp.log(mpf(q))) ** (1 / kk)))
    for mu, sc in LAPLACE:
        m, b = mpf(to_type(mu)), mpf(to_type(sc))
        out.append(("laplace_cquantile", mu, sc, [q for q in qs() if q != 0.5],
                    lambda q, m=m, b=b: m - b * mp.log(2 * mpf(q)) if q < 0.5 else m + b * mp.log(2 * (1 - mpf(q)))))
    for mu, sc in LOGISTIC:
        m, sl = mpf(to_type(mu)), mpf(to_type(sc))
        out.append(("logistic_cquantile", mu, sc, [q for q in qs() if q != 0.5], lambda q, m=m, sl=sl: m - sl * (mp.log(mpf(q)) - mp.log1p(-mpf(q)))))
    for xm, al in PARETO:
        x0, a = mpf(to_type(xm)), mpf(to_type(al))
        out.append(("pareto_cquantile", xm, al, qs(), lambda q, x0=x0, a=a: x0 * mp.exp(-mp.log(mpf(q)) / a)))
    return out


TRUNCATED = [(-1.0, 2.0), (10.0, 1e300), (-1e300, -5.0), (3.0, 4.0), (-0.5, 0.75)]   # zéro du quantile loin de p = ½


def algebra_laws(kind):
    """Lois construites (erreur de composition mesurée) : Z = −2T + 5, T de Student ; N(0, 1)
    tronquée à [a, b] ; mélange 0,3 N(−2, 1) + 0,7 gamma(3, 1). En double seulement."""
    if kind != "double":
        return []
    out = []
    count = 60
    for nu in [4.0, 1.5, 30.0]:
        n = mpf(nu)
        ys = [float(5 + rng.choice([-1, 1]) * log_uniform(2.0 ** -20, 2.0 ** 40)) for _ in range(count)]
        z = lambda y: (5 - mpf(y)) / 2                                  # cdf_Z(y) = P(T >= z)
        upper = lambda t, n=n: student_tail(n, t) if t > 0 else 1 - student_tail(n, t)
        c = mp.loggamma((n + 1) / 2) - mp.loggamma(n / 2) - mp.log(n * mp.pi) / 2
        out.append(("affine_student_pdf", nu, 0.0, ys, lambda y, n=n, c=c: mp.exp(c - (n + 1) / 2 * mp.log1p(z(y) ** 2 / n)) / 2))
        out.append(("affine_student_cdf", nu, 0.0, ys, lambda y, upper=upper: upper(z(y))))
        out.append(("affine_student_ccdf", nu, 0.0, ys, lambda y, n=n: 1 - student_tail(n, z(y)) if z(y) > 0 else student_tail(n, z(y))))
        def quantile(u, n=n):                                           # 5 − 2 Q_T(1 − u) = 5 + 2 Q_T(u)
            if u == 0.5:
                return mpf(5)
            tail = min(mpf(u), 1 - mpf(u))
            with mp.workprec(450):
                x = g.ibeta_inverse_exact(n / 2, mpf(1) / 2, 2 * tail, False)
                t = mp.sqrt(n * (1 - x) / x)
            return 5 - 2 * t if u < 0.5 else 5 + 2 * t
        out.append(("affine_student_quantile", nu, 0.0, probabilities(count // 3), quantile))
    S = lambda x: mpf(0) if x > 60 else mpf(1) if x < -60 else mp.erfc(mpf(x) / mp.sqrt(2)) / 2   # S(60) < 2^-2600
    for a, b in TRUNCATED:
        # Mesure hors des bornes (conditionnement |x|/(x − a)) et là où la survie de X reste
        # normale (|x| < 37) : ailleurs la composition perd ses bits (docs/random.md).
        lo_, hi_ = max(a, -37.0), min(b, 37.0)
        w = (hi_ - lo_) / 8
        xs = [float(rng.uniform(lo_ + w, hi_ - w)) for _ in range(count)]
        def cdf_t(x, a=a, b=b):
            with mp.workprec(1200):                                     # 1 − S(x) jusqu'à 10⁻³⁰⁰
                return (S(a) - S(x)) / (S(a) - S(b))
        def ccdf_t(x, a=a, b=b):
            with mp.workprec(1200):
                return (S(x) - S(b)) / (S(a) - S(b))
        def quantile_t(p, a=a, b=b):
            with mp.workprec(1200):
                t = S(a) - mpf(p) * (S(a) - S(b))                       # S(x) = t
                if t <= mpf(1) / 2:
                    return mp.sqrt(2) * g.erfcinv_exact(2 * t)
                u = S(-a) + mpf(p) * (S(a) - S(b))                      # Φ(x) = 1 − t, sans compensation
                return -mp.sqrt(2) * g.erfcinv_exact(2 * u)
        out.append(("truncated_normal_cdf", a, b, xs, cdf_t))
        out.append(("truncated_normal_ccdf", a, b, xs, ccdf_t))
        out.append(("truncated_normal_quantile", a, b, [p for p in probabilities(count // 3) if 2.0 ** -900 < p < 1 - 2.0 ** -40], quantile_t))
    def mix_cdf(x):
        x = mpf(x)
        with mp.workprec(400):
            g3 = 0 if x <= 0 else 1 - mp.exp(-x) * (1 + x + x ** 2 / 2)
            return mpf(3) / 10 * mp.erfc(-(x + 2) / mp.sqrt(2)) / 2 + mpf(7) / 10 * g3
    def mix_pdf(x):
        x = mpf(x)
        return mpf(3) / 10 * mp.exp(-(x + 2) ** 2 / 2) / mp.sqrt(2 * mp.pi) + (mpf(7) / 10 * x ** 2 * mp.exp(-x) / 2 if x > 0 else 0)
    def mix_quantile(p):
        lo_, hi_ = mpf(-40), mpf(800)
        with mp.workprec(300):
            for _ in range(300):
                m = (lo_ + hi_) / 2
                if mix_cdf(m) < p:
                    lo_ = m
                else:
                    hi_ = m
            return (lo_ + hi_) / 2
    # Poisson(μ) + binomiale(20, p) : convolution exacte (k entier, porté par x).
    for mu, pp in [(3.5, 0.3), (50.0, 0.9), (0.01, 0.5), (400.0, 0.05)]:
        m_, q_ = mpf(mu), mpf(pp)
        b_pmf = [mp.binomial(20, j) * q_ ** j * (1 - q_) ** (20 - j) for j in range(21)]
        top = int(mu + 20 + 12 * math.sqrt(mu + 5)) + 30
        ks = sorted({rng.randrange(0, top) for _ in range(count)})
        p_pmf = lambda i, m_=m_: mpf(0) if i < 0 else mp.exp(i * mp.log(m_) - m_ - mp.loggamma(i + 1)) if i > 0 else mp.exp(-m_)
        out.append(("sum_poisson_binomial_pmf", mu, pp, [float(k) for k in ks],
                    lambda k, b_pmf=b_pmf, p_pmf=p_pmf: mp.fsum(b_pmf[j] * p_pmf(int(k) - j) for j in range(21))))
        out.append(("sum_poisson_binomial_cdf", mu, pp, [float(k) for k in ks],
                    lambda k, b_pmf=b_pmf, mu=mu: mp.fsum(b_pmf[j] * g.igamma_reference(float(int(k) - j) + 1, mu)[1] for j in range(21) if int(k) - j >= 0)))
        out.append(("sum_poisson_binomial_ccdf", mu, pp, [float(k) for k in ks],
                    lambda k, b_pmf=b_pmf, mu=mu: mp.fsum(b_pmf[j] * (g.igamma_reference(float(int(k) - j) + 1, mu)[0] if int(k) - j >= 0 else mpf(1)) for j in range(21))))
    xs = [float(rng.uniform(-8, 20)) for _ in range(count)]
    out.append(("mixture_pdf", 0.0, 0.0, xs, mix_pdf))
    out.append(("mixture_cdf", 0.0, 0.0, xs, mix_cdf))
    # (loin de p = F(0) = 0,3 Φ(2), où le quantile s'annule : erreur relative sans borne)
    out.append(("mixture_quantile", 0.0, 0.0, [p for p in probabilities(count // 4) if p > 1e-300 and abs(p - 0.29317) > 0.05], mix_quantile))
    return out


VON_MISES = [(0.0, 0.01), (0.5, 0.5), (-2.0, 4.0), (1.0, 30.0), (0.0, 1e3), (3.0, 1e5)]
VON_MISES_FISHER = [(0.5, 3), (10.0, 3), (1e3, 3), (3.5, 5), (50.0, 10), (1e4, 4), (200.0, 120)]


def circular_laws(kind):
    """von Mises (a = μ, b = κ) ; von Mises–Fisher (a = κ, b = p, μ = e_p, x = (√(1 − t²), 0, …, 0, t))."""
    if kind != "double":
        return []
    out = []
    count = 120
    bessel = lambda nu, k: mp.besseli(mpf(nu), mpf(k), maxterms=10 ** 7)
    for mu, kappa in VON_MISES:
        m, k = mpf(mu), mpf(kappa)
        xs = [float(mu + rng.uniform(-math.pi, math.pi) * (2.0 ** -rng.uniform(0, 12) if rng.random() < 0.4 else 1)) for _ in range(count)]
        with mp.workprec(400):
            i0, i1 = bessel(0, k), bessel(1, k)
            lnorm = mp.log(2 * mp.pi) + mp.log(i0)
        out.append(("von_mises_pdf", mu, kappa, xs, lambda x, m=m, k=k, lnorm=lnorm: mp.exp(k * mp.cos(mpf(x) - m) - lnorm)))
        out.append(("von_mises_log_pdf", mu, kappa, xs, lambda x, m=m, k=k, lnorm=lnorm: k * mp.cos(mpf(x) - m) - lnorm))
        out.append(("von_mises_resultant", mu, kappa, [0.0], lambda _, i0=i0, i1=i1: i1 / i0))
        out.append(("von_mises_circular_variance", mu, kappa, [0.0], lambda _, i0=i0, i1=i1: 1 - i1 / i0))
        out.append(("von_mises_entropy", mu, kappa, [0.0], lambda _, k=k, i0=i0, i1=i1, lnorm=lnorm: lnorm - k * i1 / i0))
    for kappa, dim in VON_MISES_FISHER:
        k, d = mpf(kappa), dim
        ts = sorted({1 - 2.0 ** -rng.uniform(1, 40) if rng.random() < 0.5 else rng.uniform(-1, 1) for _ in range(count)})
        with mp.workprec(400):
            nu = mpf(d) / 2 - 1
            lc = nu * mp.log(k) - (mpf(d) / 2) * mp.log(2 * mp.pi) - mp.log(bessel(nu, k))
            a = bessel(nu + 1, k) / bessel(nu, k)
        out.append(("von_mises_fisher_pdf", kappa, float(d), ts, lambda t, k=k, lc=lc: mp.exp(lc + k * mpf(t))))
        out.append(("von_mises_fisher_log_pdf", kappa, float(d), ts, lambda t, k=k, lc=lc: lc + k * mpf(t)))
        out.append(("von_mises_fisher_resultant", kappa, float(d), [0.0], lambda _, a=a: a))
        out.append(("von_mises_fisher_entropy", kappa, float(d), [0.0], lambda _, k=k, lc=lc, a=a: -lc - k * a))
    return out


GUMBEL = [(0.0, 1.0), (1.0, 2.0), (-1e5, 3.0), (0.0, 1e-10)]
FRECHET = [(3.0, 2.0, 0.5), (0.5, 1.0, 0.0), (1.5, 3.0, -2.0), (50.0, 1.0, 0.0), (2.5, 1e-5, 1.0)]
RAYLEIGH = [1.0, 2.0, 1e-10, 3e5]
TRIANGULAR = [(0.0, 1.0, 3.0), (-1.0, -1.0, 2.0), (0.0, 5.0, 5.0), (1e3, 1e3 + 0.25, 1e3 + 1.0)]
ARCSINE = [(0.0, 1.0), (-2.0, 3.0), (1e3, 1e3 + 1e-3)]


def complement_laws(kind):
    """Gumbel, Fréchet, Rayleigh, triangulaire, arcsinus (Erlang : la gamma, déjà mesurée)."""
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 100 if kind == "double" else 40
    out = []
    qs = lambda: probabilities(count // 3)

    def standard(name, a, b, xs, cdf, ccdf, pdf, logpdf, quantile, cquantile, constants):
        out.append((f"{name}_pdf", a, b, xs, pdf))
        out.append((f"{name}_log_pdf", a, b, xs, logpdf))
        out.append((f"{name}_cdf", a, b, xs, cdf))
        out.append((f"{name}_ccdf", a, b, xs, ccdf))
        out.append((f"{name}_quantile", a, b, qs(), quantile))
        out.append((f"{name}_cquantile", a, b, qs(), cquantile))
        for c, f in constants:
            out.append((f"{name}_{c}", a, b, [0.0], lambda _, f=f: f))

    for mu, beta in GUMBEL:
        m, b = mpf(to_type(mu)), mpf(to_type(beta))
        xs = [to_type(float(m + b * rng.choice([-1, 1]) * log_uniform(2.0 ** -20, 40.0))) for _ in range(count)]
        z = lambda x, m=m, b=b: (mpf(x) - m) / b
        standard("gumbel", mu, beta, xs,
                 lambda x, z=z: mp.exp(-mp.exp(-z(x))), lambda x, z=z: -mp.expm1(-mp.exp(-z(x))),
                 lambda x, z=z, b=b: mp.exp(-z(x) - mp.exp(-z(x))) / b, lambda x, z=z, b=b: -z(x) - mp.exp(-z(x)) - mp.log(b),
                 lambda p, m=m, b=b: m - b * mp.log(-mp.log(mpf(p))), lambda q, m=m, b=b: m - b * mp.log(-mp.log1p(-mpf(q))),
                 [("mean", m + b * mp.euler), ("variance", b ** 2 * mp.pi ** 2 / 6), ("median", m - b * mp.log(mp.log(2))), ("entropy", mp.log(b) + mp.euler + 1)])
    for al, sc, loc in FRECHET:
        a, sv, lv = mpf(to_type(al)), mpf(to_type(sc)), mpf(to_type(loc))
        xs = [to_type(float(lv + sv * log_uniform(2.0 ** -1, 2.0 ** 20) ** (1 / a) * (log_uniform(0.2, 3.0) if rng.random() < 0.5 else 1))) for _ in range(count)]
        xs = [x for x in xs if x > float(lv) and x < math.inf]
        zz = lambda x, sv=sv, lv=lv: (mpf(x) - lv) / sv
        const = [("median", lv + sv * mp.log(2) ** (-1 / a)), ("mode", lv + sv * (a / (1 + a)) ** (1 / a)), ("entropy", 1 + mp.euler / a + mp.euler + mp.log(sv / a))]
        if al > 1:
            const.append(("mean", lv + sv * mp.gamma(1 - 1 / a)))
        if al > 2:
            const.append(("variance", sv ** 2 * (mp.gamma(1 - 2 / a) - mp.gamma(1 - 1 / a) ** 2)))
        standard("frechet", al, sc, xs,
                 lambda x, zz=zz, a=a: mp.exp(-zz(x) ** -a), lambda x, zz=zz, a=a: -mp.expm1(-zz(x) ** -a),
                 lambda x, zz=zz, a=a, sv=sv: a / sv * zz(x) ** (-1 - a) * mp.exp(-zz(x) ** -a),
                 lambda x, zz=zz, a=a, sv=sv: mp.log(a / sv) - (1 + a) * mp.log(zz(x)) - zz(x) ** -a,
                 lambda p, a=a, sv=sv, lv=lv: lv + sv * (-mp.log(mpf(p))) ** (-1 / a), lambda q, a=a, sv=sv, lv=lv: lv + sv * (-mp.log1p(-mpf(q))) ** (-1 / a), const)
    for sig in RAYLEIGH:
        sv = mpf(to_type(sig))
        if kind == "float" and not 2.0 ** -100 < sig < 2.0 ** 100:
            continue
        xs = [to_type(float(sv * log_uniform(2.0 ** -30, 30.0))) for _ in range(count)]
        h = lambda x, sv=sv: mpf(x) ** 2 / (2 * sv ** 2)
        standard("rayleigh", sig, 0.0, xs, lambda x, h=h: -mp.expm1(-h(x)), lambda x, h=h: mp.exp(-h(x)),
                 lambda x, h=h, sv=sv: mpf(x) / sv ** 2 * mp.exp(-h(x)), lambda x, h=h, sv=sv: mp.log(mpf(x) / sv ** 2) - h(x),
                 lambda p, sv=sv: sv * mp.sqrt(-2 * mp.log1p(-mpf(p))), lambda q, sv=sv: sv * mp.sqrt(-2 * mp.log(mpf(q))),
                 [("mean", sv * mp.sqrt(mp.pi / 2)), ("variance", (4 - mp.pi) / 2 * sv ** 2), ("median", sv * mp.sqrt(2 * mp.log(2))), ("entropy", 1 + mp.log(sv / mp.sqrt(2)) + mp.euler / 2)])
    for lo_, mo_, hi_ in TRIANGULAR:
        a, c, b = mpf(to_type(lo_)), mpf(to_type(mo_)), mpf(to_type(hi_))
        xs = [to_type(float(a + (b - a) * mpf(rng.random()))) for _ in range(count)] + [to_type(float(a + (c - a) * mpf(2.0 ** -rng.uniform(1, 40)))) for _ in range(count // 4)]
        xs = [x for x in xs if float(a) < x < float(b)]
        def cdf_t(x, a=a, b=b, c=c):
            x = mpf(x)
            return (x - a) ** 2 / ((b - a) * (c - a)) if x <= c else 1 - (b - x) ** 2 / ((b - a) * (b - c))
        def pdf_t(x, a=a, b=b, c=c):
            x = mpf(x)
            return 2 * (x - a) / ((b - a) * (c - a)) if x < c else 2 * (b - x) / ((b - a) * (b - c))
        def q_t(p, a=a, b=b, c=c):
            p = mpf(p)
            return a + mp.sqrt(p * (b - a) * (c - a)) if p * (b - a) <= c - a else b - mp.sqrt((1 - p) * (b - a) * (b - c))
        out.append(("triangular_pdf", lo_, hi_, [x for x in xs if float(x) != float(c)], pdf_t))
        out.append(("triangular_cdf", lo_, hi_, xs, cdf_t))
        out.append(("triangular_ccdf", lo_, hi_, xs, lambda x, a=a, b=b, c=c, cdf_t=cdf_t: (b - mpf(x)) ** 2 / ((b - a) * (b - c)) if mpf(x) > c else 1 - (mpf(x) - a) ** 2 / ((b - a) * (c - a))))
        out.append(("triangular_quantile", lo_, hi_, qs(), q_t))
        out.append(("triangular_cquantile", lo_, hi_, qs(), lambda q, q_t=q_t: q_t(1 - mpf(q))))
        out.append(("triangular_variance", lo_, hi_, [0.0], lambda _, a=a, b=b, c=c: ((a - b) ** 2 + (a - c) ** 2 + (b - c) ** 2) / 36))
        out.append(("triangular_mean", lo_, hi_, [0.0], lambda _, a=a, b=b, c=c: (a + b + c) / 3))
    for lo_, hi_ in ARCSINE:
        a, b = mpf(to_type(lo_)), mpf(to_type(hi_))
        xs = [to_type(float(a + (b - a) * mpf(rng.random()))) for _ in range(count)] + [to_type(float(a + (b - a) * mpf(2.0 ** -rng.uniform(1, 40)))) for _ in range(count // 4)] \
             + [to_type(float(b - (b - a) * mpf(2.0 ** -rng.uniform(1, 40)))) for _ in range(count // 4)]
        xs = [x for x in xs if float(a) < x < float(b)]
        u = lambda x, a=a, b=b: (mpf(x) - a) / (b - a)
        def cdf_a(x, u=u):
            v = u(x)
            return 2 / mp.pi * mp.atan(mp.sqrt(v / (1 - v)))
        def ccdf_a(x, u=u):
            v = u(x)
            return 2 / mp.pi * mp.atan(mp.sqrt((1 - v) / v))
        def q_a(p, a=a, b=b):
            p = mpf(p)
            return a + (b - a) * mp.sin(mp.pi * p / 2) ** 2 if p <= mpf(1) / 2 else b - (b - a) * mp.sin(mp.pi * (1 - p) / 2) ** 2
        def cq_a(q, a=a, b=b):
            q = mpf(q)
            return b - (b - a) * mp.sin(mp.pi * q / 2) ** 2 if q <= mpf(1) / 2 else a + (b - a) * mp.sin(mp.pi * (1 - q) / 2) ** 2
        standard("arcsine", lo_, hi_, xs, cdf_a, ccdf_a,
                 lambda x, a=a, b=b: 1 / (mp.pi * mp.sqrt((mpf(x) - a) * (b - mpf(x)))), lambda x, a=a, b=b: -mp.log(mp.pi) - mp.log((mpf(x) - a) * (b - mpf(x))) / 2,
                 q_a, cq_a, [("variance", (b - a) ** 2 / 8), ("entropy", mp.log(mp.pi * (b - a) / 4))])
    return out


INVERSE_GAUSSIAN = [(1.0, 1.0), (1.0, 0.01), (2.0, 50.0), (1e-3, 1e3), (1e5, 3.0), (1.0, 1e-6)]


def inverse_gaussian_parts(mu, lam, x, prec=320):
    """(densité, répartition, survie) de l'inverse-gaussienne (survie compensée d'au plus 2^35 ici)."""
    with mp.workprec(prec):
        x = mpf(x)
        c = mp.sqrt(lam / x)
        a, b = c * (x / mu - 1), c * (x / mu + 1)
        tail = mp.exp(2 * lam / mu) * mp.erfc(b / mp.sqrt(2)) / 2
        lower = mp.erfc(-a / mp.sqrt(2)) / 2 + tail
        upper = mp.erfc(a / mp.sqrt(2)) / 2 - tail
        density = mp.sqrt(lam / (2 * mp.pi * x ** 3)) * mp.exp(-lam * (x - mu) ** 2 / (2 * mu ** 2 * x))
        return +density, +lower, +upper


def inverse_gaussian_quantile(mu, lam, t, upper):
    """x tel que F(x) = t (ou S(x) = t) : bissection en log x à 160 bits, puis Newton à 320 bits."""
    t = mpf(t)
    value = lambda x: inverse_gaussian_parts(mu, lam, x, 160)[2 if upper else 1]
    lo, hi = mpf(2) ** -1100, mpf(2) ** 1100
    for _ in range(50):
        mid = mp.sqrt(lo * hi)
        if (value(mid) > t) != upper:
            hi = mid
        else:
            lo = mid
    x = mp.sqrt(lo * hi)
    with mp.workprec(320):
        for _ in range(10):
            f, F, S = inverse_gaussian_parts(mu, lam, x, 320)
            step = ((S if upper else F) - t) / f
            x = x + step if upper else x - step
            if abs(step) < x * mpf(2) ** -200:
                break
        return +x


def inverse_gaussian_laws(kind):
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 80 if kind == "double" else 30
    out = []
    for mu_, lam_ in INVERSE_GAUSSIAN:
        mu, lam = mpf(to_type(mu_)), mpf(to_type(lam_))
        xs = [to_type(float(mu * log_uniform(2.0 ** -12, 2.0 ** 30))) for _ in range(count)]
        xs += [to_type(float(mu * (1 + rng.choice([-1, 1]) * log_uniform(2.0 ** -30, 0.5)))) for _ in range(count // 2)]
        xs = sorted({x for x in xs if 0 < x < math.inf})
        out.append(("inverse_gaussian_pdf", mu_, lam_, xs, lambda x, mu=mu, lam=lam: inverse_gaussian_parts(mu, lam, x)[0]))
        out.append(("inverse_gaussian_log_pdf", mu_, lam_, xs, lambda x, mu=mu, lam=lam: mp.log(inverse_gaussian_parts(mu, lam, x)[0])))
        out.append(("inverse_gaussian_cdf", mu_, lam_, xs, lambda x, mu=mu, lam=lam: inverse_gaussian_parts(mu, lam, x)[1]))
        out.append(("inverse_gaussian_ccdf", mu_, lam_, xs, lambda x, mu=mu, lam=lam: inverse_gaussian_parts(mu, lam, x)[2]))
        qs = [q for q in probabilities(count // 5) if q > 1e-300 or kind == "double"]
        def quantile(p, mu=mu, lam=lam):
            p = mpf(p)
            return inverse_gaussian_quantile(mu, lam, p, False) if p <= 0.5 else inverse_gaussian_quantile(mu, lam, 1 - p, True)
        def cquantile(q, mu=mu, lam=lam):
            q = mpf(q)
            return inverse_gaussian_quantile(mu, lam, q, True) if q <= 0.5 else inverse_gaussian_quantile(mu, lam, 1 - q, False)
        out.append(("inverse_gaussian_quantile", mu_, lam_, qs, quantile))
        out.append(("inverse_gaussian_cquantile", mu_, lam_, qs, cquantile))
        t = 3 * mu / (2 * lam)
        for c, v in [("variance", mu ** 3 / lam), ("mode", mu * (mp.sqrt(1 + t * t) - t)), ("skewness", 3 * mp.sqrt(mu / lam)), ("excess_kurtosis", 15 * mu / lam)]:
            out.append((f"inverse_gaussian_{c}", mu_, lam_, [0.0], lambda _, v=v: v))
    return out


def von_mises_norm(k):
    return 2 * mp.pi * mp.besseli(0, k, maxterms=10 ** 7) * mp.exp(-k) if k > 0 else 2 * mp.pi


def von_mises_tail(k, phi, prec):
    """(G(φ), g(φ)), G(φ) = ∫_φ^π e^{κ(cos t − 1)} dt / (2π e^{−κ} I₀(κ)). En w = 2κ (sin²(t/2) − sin²(φ/2)),
    dt = dw/(κ sin t), sin t = 2√(x(1 − x)), x = sin²(t/2) : sur [0, w(π)/2], coupures géométriques
    depuis l'échelle 2κ sin²(φ/2) de la quasi-singularité en w = −2κ sin²(φ/2) ; sur [w(π)/2, w(π)],
    en u = √(w(π) − w), qui lève la singularité en t = π (intégrande e^{−w} √2/√(κ x))."""
    with mp.workprec(prec):
        k, phi = mpf(k), mpf(phi)
        norm = von_mises_norm(k)
        s0, c0 = mp.sin(phi / 2) ** 2, mp.cos(phi / 2) ** 2
        density = mp.exp(-2 * k * s0) / norm
        if k == 0:
            return (mp.pi - phi) / norm, density
        top = 2 * k * c0                                        # w(π)
        cap = mpf(3 * prec)                                     # e^{−w} négligeable au-delà
        half = top / 2
        lower_end = min(half, cap)
        f = lambda w: mp.exp(-w) / (2 * k * mp.sqrt((s0 + w / (2 * k)) * (c0 - w / (2 * k))))
        pts = [mpf(0)]
        b = max(2 * k * s0, mpf(2) ** -prec)
        while b < lower_end:
            if b > pts[-1]:
                pts.append(b)
            b *= 2
        pts.append(lower_end)
        total = mp.quad(f, pts)
        if half < cap:
            g = lambda u: mp.exp(-(top - u * u)) * mp.sqrt(2) / mp.sqrt(k * (1 - u * u / (2 * k)))
            total += mp.quad(g, [0, mp.sqrt(half)])
        return mp.exp(-2 * k * s0) * total / norm, density


def von_mises_head(k, phi, prec):
    """D(φ) = ∫_0^φ e^{κ(cos t − 1)} dt / (2π e^{−κ} I₀(κ)) (φ petit : D = ½ − G sans compensation)."""
    with mp.workprec(prec):
        k, phi = mpf(k), mpf(phi)
        n = int(min(64, 2 + phi * mp.sqrt(k + 1)))
        pts = [phi * j / n for j in range(n + 1)]
        return mp.quad(lambda t: mp.exp(-2 * k * mp.sin(t / 2) ** 2), pts) / von_mises_norm(k)


def von_mises_solve(k, t, low_prec=53, prec=120):
    """φ tel que G(φ) = t (t <= 1/2) : par G pour t <= 1/4, par D(φ) = 1/2 − t au-delà."""
    t = mpf(t)
    head = t > mpf(1) / 4
    r = mpf(1) / 2 - t if head else t
    if r == 0:
        return mpf(0)
    value = (lambda phi, pr: von_mises_head(k, phi, pr)) if head else (lambda phi, pr: von_mises_tail(k, phi, pr)[0])
    lo, hi = mpf(0), +mp.pi
    for _ in range(25):
        mid = (lo + hi) / 2
        if (value(mid, low_prec) > r) != head:
            lo = mid
        else:
            hi = mid
    phi = (lo + hi) / 2
    with mp.workprec(prec):
        for _ in range(10):
            g = von_mises_tail(k, phi, prec)[1]
            step = (value(phi, prec) - r) / g
            phi = phi - step if head else phi + step
            if abs(step) < mpf(2) ** -100 * (phi + mpf(2) ** -1000):
                break
        return +phi


def von_mises_distribution(kind):
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 50 if kind == "double" else 20
    out = []
    for mu_, kappa_ in VON_MISES:
        m, k = mpf(to_type(mu_)), mpf(to_type(kappa_))
        thetas = [rng.uniform(-math.pi, math.pi) for _ in range(count)]
        thetas += [rng.choice([-1, 1]) * math.pi * (1 - 2.0 ** -rng.uniform(1, 45)) for _ in range(count // 4)]
        thetas += [rng.choice([-1, 1]) * 2.0 ** -rng.uniform(1, 40) for _ in range(count // 4)]
        xs = sorted({to_type(float(m) + th) for th in thetas})
        xs = [x for x in xs if abs(mpf(x) - m) < mp.pi]
        def cdf(x, m=m, k=k, upper=False):
            th = mpf(x) - m
            G = von_mises_tail(k, abs(th), 120)[0]
            return G if (th < 0) != upper else 1 - G
        out.append(("von_mises_cdf", mu_, kappa_, xs, cdf))
        out.append(("von_mises_ccdf", mu_, kappa_, xs, lambda x, cdf=cdf: cdf(x, upper=True)))
        qs = probabilities(count // 8)
        def quantile(p, m=m, k=k, upper=False):
            p = mpf(p)
            low = not upper
            if p > 0.5:
                p, low = 1 - p, not low
            phi = von_mises_solve(k, p)
            return m - phi if low else m + phi
        out.append(("von_mises_quantile", mu_, kappa_, qs, quantile))
        out.append(("von_mises_cquantile", mu_, kappa_, qs, lambda q, quantile=quantile: quantile(q, upper=True)))
    return out


NEGATIVE_BINOMIAL = [(3.5, 0.4), (0.3, 0.05), (1e3, 0.999), (50.0, 0.5), (1e-3, 0.9), (2.0, 1e-4)]
HYPERGEOMETRIC = [(50, 20, 10), (1000, 300, 400), (10 ** 6, 10 ** 3, 5 * 10 ** 4), (10 ** 9, 5 * 10 ** 8, 10 ** 6), (30, 29, 15)]


def negative_binomial_log_pmf(r, pp, k):
    with mp.workprec(300):
        r, pp, k = mpf(r), mpf(pp), mpf(k)
        return mp.loggamma(k + r) - mp.loggamma(k + 1) - mp.loggamma(r) + r * mp.log(pp) + k * mp.log1p(-pp)


def hypergeometric_log_pmf(N, K, n, k):
    lc = lambda a, b: mp.loggamma(a + 1) - mp.loggamma(b + 1) - mp.loggamma(a - b + 1)
    with mp.workprec(200):
        return lc(mpf(K), mpf(k)) + lc(mpf(N - K), mpf(n - k)) - lc(mpf(N), mpf(n))


def hypergeometric_tail(N, K, n, k, upper):
    """Σ_{j<=k} P(j) (ou Σ_{j>k}) par le rapport de masses consécutives, à 320 bits."""
    lo, hi = max(0, n - (N - K)), min(n, K)
    with mp.workprec(200):
        j = k + 1 if upper else k
        if j < lo or j > hi:
            return mpf(0)
        base = mp.exp(hypergeometric_log_pmf(N, K, n, j))
        total, term = mpf(1), mpf(1)
        while (j < hi) if upper else (j > lo):
            if upper:
                term *= mpf(K - j) * (n - j) / (mpf(j + 1) * (N - K - n + j + 1))
                j += 1
            else:
                term *= mpf(j) * (N - K - n + j) / (mpf(K - j + 1) * (n - j + 1))
                j -= 1
            total += term
            if term < total * mpf(2) ** -130:
                break
        return base * total


HYPERGEOMETRIC_DRAWS = {(N, K): n for N, K, n in HYPERGEOMETRIC}

BETA_BINOMIAL = [(20, 2.0, 3.0), (1000, 0.5, 0.5), (10 ** 6, 3.0, 40.0), (50, 0.05, 2.0)]
BETA_BINOMIAL_BETA = {(n, a): b for n, a, b in BETA_BINOMIAL}
LOGARITHMIC = [0.9, 0.01, 0.5, 0.999, 1e-6]
DISCRETE_GAUSSIAN = [(0.3, 2.5), (0.3, 0.4), (-7.25, 30.0), (1e6 + 0.5, 1.0), (0.0, 0.05), (2.0, 1e3)]


def beta_binomial_log_pmf(n, a, b, k):
    with mp.workprec(200):
        a, b = mpf(a), mpf(b)
        return mp.log(mp.binomial(n, k)) + mp.log(mp.beta(k + a, n - k + b)) - mp.log(mp.beta(a, b))


_PREFIX_SUMS = {}


def prefix_sum(key, first, ratio, k):
    """Σ_{j<=k} t_j, t_0 = first(), t_{j+1} = t_j ratio(j), à 200 bits ; sommes cumulées gardées par
    processus et prolongées à la demande."""
    with mp.workprec(200):
        if key not in _PREFIX_SUMS:
            t = first()
            _PREFIX_SUMS[key] = ([t], [t])
        terms, sums = _PREFIX_SUMS[key]
        while len(sums) <= k:
            j = len(terms) - 1
            terms.append(terms[-1] * ratio(j))
            sums.append(sums[-1] + terms[-1])
        return sums[k]


def beta_binomial_cdf(n, a, b, k):
    """Σ_{j<=k} P(j) par le rapport de masses, à 200 bits."""
    if k >= n:
        return mpf(1)
    a_, b_ = mpf(a), mpf(b)
    return prefix_sum(("beta_binomial", n, a, b), lambda: mp.exp(beta_binomial_log_pmf(n, a, b, 0)),
                      lambda j: mpf(n - j) * (j + a_) / ((j + 1) * (n - j - 1 + b_)), k)


def logarithmic_cdf(pp, k):
    """Σ_{j=1}^{k} p^j/j / L, à 200 bits (terme d'indice i : p^{i+1}/((i + 1) L))."""
    if k < 1:
        return mpf(0)
    with mp.workprec(200):
        p_ = mpf(pp)
        L = -mp.log1p(-p_)
    return prefix_sum(("logarithmic", pp), lambda: p_ / L, lambda i: p_ * (i + 1) / (i + 2), k - 1)


_DISCRETE_GAUSSIAN_TABLES = {}


def discrete_gaussian_table(mu, sigma):
    """(c, span, termes f(c − span … c + span), sommes cumulées, Z) à 200 bits, une fois par processus ;
    au-delà de c ± (40σ + 40), les termes sont sous 2^-1100 (les quantiles de probabilités
    sous-normales tombent jusque vers 39σ)."""
    key = (mu, sigma)
    if key not in _DISCRETE_GAUSSIAN_TABLES:
        with mp.workprec(200):
            m, s_ = mpf(mu), mpf(sigma)
            c, span = int(mp.nint(m)), int(40 * sigma) + 40
            terms = [mp.exp(-(j - m) ** 2 / (2 * s_ ** 2)) for j in range(c - span, c + span + 1)]
            prefix, acc = [], mpf(0)
            for t in terms:
                acc += t
                prefix.append(acc)
            _DISCRETE_GAUSSIAN_TABLES[key] = (c, span, terms, prefix, acc)
    return _DISCRETE_GAUSSIAN_TABLES[key]


def discrete_gaussian_parts(mu, sigma, k):
    """(P(k), F(k), 1 − F(k)) par sommes directes à 200 bits."""
    c, span, terms, prefix, z = discrete_gaussian_table(mu, sigma)
    with mp.workprec(200):
        i = k - (c - span)
        if i < 0:
            return mpf(0), mpf(0), mpf(1)
        if i >= len(terms):
            return mpf(0), mpf(1), mpf(0)
        lower = prefix[i]
        upper = z - lower if lower > z / 2 else mp.fsum(terms[i + 1:])
        return terms[i] / z, lower / z, upper / z


def hypergeometric_cdf_reference(N, K, k):
    n = HYPERGEOMETRIC_DRAWS[(N, K)]
    if k >= min(n, K):
        return mpf(1)
    return hypergeometric_tail(N, K, n, k, False)


def hypergeometric_start(N, K, u):
    n = HYPERGEOMETRIC_DRAWS[(N, K)]
    mean, sd = n * K / N, math.sqrt(max(n * K * (N - K) * (N - n) / (N * N * (N - 1)), 0.25))
    return min(min(n, K), max(max(0, n - (N - K)), mean + sd * float(-mp.sqrt(2) * g.erfcinv_exact(2 * mpf(u)))))


def new_discrete_laws():
    out = []
    for n, a, b in BETA_BINOMIAL:
        mean = n * a / (a + b); sd = math.sqrt(n * a * b * (a + b + n) / ((a + b) ** 2 * (a + b + 1)))
        ks = sorted({min(n, max(0, int(mean + sd * rng.uniform(-5, 5)))) for _ in range(40)} | {0, n})
        out.append(("beta_binomial_pmf", float(n), a, [float(k) for k in ks], lambda k, n=n, a=a, b=b: mp.exp(beta_binomial_log_pmf(n, a, b, int(k)))))
        out.append(("beta_binomial_log_pmf", float(n), a, [float(k) for k in ks], lambda k, n=n, a=a, b=b: beta_binomial_log_pmf(n, a, b, int(k))))
        if n <= 1000:
            inner = [float(k) for k in ks if k < n]
            out.append(("beta_binomial_cdf", float(n), a, inner, lambda k, n=n, a=a, b=b: beta_binomial_cdf(n, a, b, int(k))))
            out.append(("beta_binomial_ccdf", float(n), a, inner, lambda k, n=n, a=a, b=b: 1 - beta_binomial_cdf(n, a, b, int(k))))
        out.append(("beta_binomial_variance", float(n), a, [0.0], lambda _, n=n, a=a, b=b: mpf(n) * a * b * (a + b + n) / ((mpf(a) + b) ** 2 * (mpf(a) + b + 1))))
    for pp in LOGARITHMIC:
        top = int(min(3000, 40 / -math.log(pp)))
        ks = sorted({max(1, int(2.0 ** rng.uniform(0, math.log2(top + 1)))) for _ in range(40)} | {1, 2})
        L = -mp.log1p(-mpf(pp))
        out.append(("logarithmic_pmf", pp, 0.0, [float(k) for k in ks], lambda k, pp=pp, L=L: mpf(pp) ** int(k) / (int(k) * L)))
        out.append(("logarithmic_cdf", pp, 0.0, [float(k) for k in ks], lambda k, pp=pp: logarithmic_cdf(pp, int(k))))
        out.append(("logarithmic_ccdf", pp, 0.0, [float(k) for k in ks], lambda k, pp=pp: 1 - logarithmic_cdf(pp, int(k))))
        out.append(("logarithmic_mean", pp, 0.0, [0.0], lambda _, pp=pp, L=L: mpf(pp) / ((1 - mpf(pp)) * L)))
        out.append(("logarithmic_variance", pp, 0.0, [0.0], lambda _, pp=pp, L=L: mpf(pp) * (L - pp) / ((1 - mpf(pp)) ** 2 * L ** 2)))
    for mu, sigma in DISCRETE_GAUSSIAN:
        c = round(mu)
        ks = sorted({int(c + sigma * rng.uniform(-9, 9)) for _ in range(40)} | {c, c - 1, c + 1})
        for name, idx in (("pmf", 0), ("cdf", 1), ("ccdf", 2)):
            out.append((f"discrete_gaussian_{name}", mu, sigma, [float(k) for k in ks], lambda k, mu=mu, sigma=sigma, idx=idx: discrete_gaussian_parts(mu, sigma, int(k))[idx]))
        def moments(mu=mu, sigma=sigma, which=0):
            c, span, terms, prefix, z = discrete_gaussian_table(mu, sigma)
            with mp.workprec(200):
                js = range(c - span, c + span + 1)
                m1 = mp.fsum(j * t for j, t in zip(js, terms)) / z
                return m1 if which == 0 else mp.fsum((j - m1) ** 2 * t for j, t in zip(js, terms)) / z
        out.append(("discrete_gaussian_mean", mu, sigma, [0.0], lambda _, moments=moments: moments(which=0)))
        out.append(("discrete_gaussian_variance", mu, sigma, [0.0], lambda _, moments=moments: moments(which=1)))
    for r, pp in NEGATIVE_BINOMIAL:
        q = 1 - pp
        mean, sd = r * q / pp, math.sqrt(r * q) / pp
        ks = sorted({max(0, int(mean + sd * rng.uniform(-6, 12))) for _ in range(60)} | {0, 1, int(mean)} | {int(mean + sd * 2.0 ** rng.uniform(3, 6)) for _ in range(8)})
        out.append(("negative_binomial_pmf", r, pp, [float(k) for k in ks], lambda k, r=r, pp=pp: mp.exp(negative_binomial_log_pmf(r, pp, k))))
        out.append(("negative_binomial_log_pmf", r, pp, [float(k) for k in ks if k != 0 or r != 1], lambda k, r=r, pp=pp: negative_binomial_log_pmf(r, pp, k)))
        out.append(("negative_binomial_cdf", r, pp, [float(k) for k in ks], lambda k, r=r, pp=pp: g.ibeta_reference(r, float(k) + 1, pp)[0]))
        out.append(("negative_binomial_ccdf", r, pp, [float(k) for k in ks], lambda k, r=r, pp=pp: g.ibeta_reference(r, float(k) + 1, pp)[1]))
    for N, K, n in HYPERGEOMETRIC:
        lo, hi = max(0, n - (N - K)), min(n, K)
        mean = n * K / N
        sd = math.sqrt(max(n * K * (N - K) * (N - n) / (N * N * (N - 1)), 0.25))
        ks = sorted({min(hi, max(lo, int(mean + sd * rng.uniform(-9, 9)))) for _ in range(50)} | {lo, hi})
        out.append(("hypergeometric_pmf", float(N), float(K), [float(k) for k in ks], lambda k, N=N, K=K, n=n: mp.exp(hypergeometric_log_pmf(N, K, n, int(k)))))
        out.append(("hypergeometric_log_pmf", float(N), float(K), [float(k) for k in ks], lambda k, N=N, K=K, n=n: hypergeometric_log_pmf(N, K, n, int(k))))
        inner = [float(k) for k in ks if lo <= k < hi]
        out.append(("hypergeometric_cdf", float(N), float(K), inner, lambda k, N=N, K=K, n=n: hypergeometric_tail(N, K, n, int(k), False)))
        out.append(("hypergeometric_ccdf", float(N), float(K), inner, lambda k, N=N, K=K, n=n: hypergeometric_tail(N, K, n, int(k), True)))
        out.append(("hypergeometric_mean", float(N), float(K), [0.0], lambda _, N=N, K=K, n=n: mpf(n) * K / N))
        out.append(("hypergeometric_variance", float(N), float(K), [0.0], lambda _, N=N, K=K, n=n: mpf(n) * K * (N - K) * (N - n) / (mpf(N) ** 2 * (N - 1))))
    return out


SKEW_NORMAL = [(4.0, 0.5, 2.0), (1.0, 0.0, 1.0), (-3.0, -1.0, 0.5), (50.0, 0.0, 1.0), (0.3, 2.0, 3.0)]
SKEW_NORMAL_SCALE = {(a, xi): om for a, xi, om in SKEW_NORMAL}


def skew_normal_small_tail(a, z, prec=110):
    """(queue, lower) : ∫ f, f = 2φ(t)Φ(αt), sur la queue issue de z du côté opposé à la moyenne
    (la plus petite) ; morceaux de largeur min(1, 4/|(log f)′|, 1/√|(log f)″|), Gauss–Legendre
    de mpmath, arrêt quand f/|(log f)′| devient négligeable au-delà du mode."""
    with mp.workprec(prec):
        a, z = mpf(a), mpf(z)
        delta = a / mp.sqrt(1 + a * a)
        lower = z < delta * mp.sqrt(2 / mp.pi)
        f = lambda t: 2 * mp.npdf(t) * mp.ncdf(a * t)
        def slopes(t):
            u = float(a * t)
            big = float(mp.ncdf(a * t))
            r = float(mp.npdf(a * t)) / big if big > 1e-300 else -u - 1 / u
            return float(a) * r - float(t), -1 - float(a) ** 2 * r * (u + r)
        scale = f(z)                                                    # mp.quad juge sa convergence en absolu
        if scale == 0:
            return mpf(0), lower
        g = lambda t: f(t) / scale
        t, total = z, mpf(0)
        for _ in range(100000):
            first, second = slopes(t)
            outward = first if lower else -first
            ft = f(t)
            if total > 0 and outward > 0 and ft / scale / outward < total * mpf(2) ** -100:
                break
            if abs(t) > 80:
                break
            w = 1.0
            if abs(first) * w > 4:
                w = 4 / abs(first)
            if abs(second) * w * w > 1:
                w = 1 / math.sqrt(abs(second))
            if a != 0:                                                      # Φ(αt) bascule sur 1/|α| autour de 0
                w = min(w, max(1 / abs(float(a)), abs(float(t)) / 2))
            nxt = t - w if lower else t + w
            total += abs(mp.quad(g, [t, nxt], method="gauss-legendre"))
            t = nxt
        return total * scale, lower


def skew_normal_parts(a, xi, om, x):
    """(densité, F, S) en x."""
    with mp.workprec(110):
        z = (mpf(x) - xi) / om
        tail, lower = skew_normal_small_tail(a, z)
        other = 1 - tail
        return 2 * mp.npdf(z) * mp.ncdf(mpf(a) * z) / om, (tail if lower else other), (other if lower else tail)


def skew_normal_quantile(a, xi, om, t, upper):
    """z tel que F (ou S) = t, t <= ½ : Newton sur log V, départ normal, encadré."""
    with mp.workprec(110):
        t = mpf(t)
        delta = mpf(a) / mp.sqrt(1 + mpf(a) ** 2)
        m, sd = delta * mp.sqrt(2 / mp.pi), mp.sqrt(1 - 2 * delta ** 2 / mp.pi)
        q = -mp.sqrt(2) * mp.erfinv(1 - 2 * t)
        z = m - sd * q if upper else m + sd * q
        lo, hi = -mp.inf, mp.inf
        for _ in range(200):
            f, F, S = skew_normal_parts(a, 0, 1, z)
            V = S if upper else F
            if (V > t) != upper:
                hi = min(hi, z)
            else:
                lo = max(lo, z)
            if V == 0 or f == 0:
                step = mpf(-4) if (V > t) != upper else mpf(4)
            else:
                step = mp.log(V / t) * V / f * (1 if upper else -1)
            nxt = z + step
            if not (lo < nxt < hi) or abs(step) > 4:
                nxt = (lo + hi) / 2 if lo > -mp.inf and hi < mp.inf else z + (4 if step > 0 else -4)
            if abs(nxt - z) < mpf(2) ** -90 * (1 + abs(z)):
                z = nxt
                break
            z = nxt
        return xi + om * z


def skew_normal_laws(kind):
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 30 if kind == "double" else 12
    out = []
    for a, xi, om in SKEW_NORMAL:
        A, X, W = mpf(to_type(a)), mpf(to_type(xi)), mpf(to_type(om))
        delta = A / mp.sqrt(1 + A * A)
        m = float(X + W * delta * mp.sqrt(2 / mp.pi))
        zs = [rng.uniform(-4, 4) for _ in range(count)] + [rng.choice([-1, 1]) * 2.0 ** rng.uniform(2, 5.2) for _ in range(count // 2)]
        xs = sorted({to_type(m + float(W) * z) for z in zs})
        parts = lambda x, A=A, X=X, W=W: skew_normal_parts(A, X, W, x)
        out.append(("skew_normal_pdf", a, xi, xs, lambda x, A=A, X=X, W=W: 2 * mp.npdf((mpf(x) - X) / W) * mp.ncdf(A * (mpf(x) - X) / W) / W))
        out.append(("skew_normal_cdf", a, xi, xs, lambda x, parts=parts: parts(x)[1]))
        out.append(("skew_normal_ccdf", a, xi, xs, lambda x, parts=parts: parts(x)[2]))
        qs = [q for q in probabilities(3 if kind == "double" else 1)]
        def quantile(p, A=A, X=X, W=W, upper=False):
            p = mpf(p)
            return skew_normal_quantile(A, X, W, p, upper) if p <= 0.5 else skew_normal_quantile(A, X, W, 1 - p, not upper)
        out.append(("skew_normal_quantile", a, xi, qs, quantile))
        out.append(("skew_normal_cquantile", a, xi, qs, lambda q, quantile=quantile: quantile(q, upper=True)))
        out.append(("skew_normal_mean", a, xi, [0.0], lambda _, X=X, W=W, delta=delta: X + W * delta * mp.sqrt(2 / mp.pi)))
        out.append(("skew_normal_variance", a, xi, [0.0], lambda _, W=W, delta=delta: W ** 2 * (1 - 2 * delta ** 2 / mp.pi)))
        out.append(("skew_normal_skewness", a, xi, [0.0], lambda _, delta=delta: (4 - mp.pi) / 2 * (delta * mp.sqrt(2 / mp.pi)) ** 3 / (1 - 2 * delta ** 2 / mp.pi) ** mpf(1.5)))
    return out


def discrete_quantiles(name, cases, cdf, start, low=0):
    """Lignes {a, b, u, k} : k = plus petit entier >= low tel que F(k) >= u (F exacte, F(low − 1) = 0),
    sauf ambiguïté (F(k) ou F(k − 1) à 2^-90 près de u)."""
    tasks = [(a, b, u) for a, b in cases for u in probabilities(40)]       # tirages dans l'ordre, avant le calcul
    ambiguous = mpf(-2) ** 70
    def solve(a, b, u):
        k = max(low, int(start(a, b, u)))                                  # recherche exponentielle, puis dichotomie
        if cdf(a, b, k) >= u:
            hi, step = k, 1
            lo = hi - step
            while lo >= low and cdf(a, b, lo) >= u:
                hi, step = lo, 2 * step
                lo = hi - step
            lo = max(lo, low - 1)                                           # F(lo) < u
        else:
            lo, step = k, 1
            hi = lo + step
            while cdf(a, b, hi) < u:
                lo, step = hi, 2 * step
                hi = lo + step
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if cdf(a, b, mid) >= u:
                hi = mid
            else:
                lo = mid
        k = hi
        near = lambda j: j >= low and abs(cdf(a, b, j) - u) < mpf(2) ** -90 * u
        return ambiguous if near(k) or near(k - 1) or k >= 2 ** 62 else mpf(k)
    rows = []
    for (a, b, u), k in zip(tasks, g.parallel_map(solve, tasks)):
        if k != ambiguous:
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
        entries = laws(kind) + (discrete_laws() if kind == "double" else []) + phase4_laws(kind) + derived_laws(kind) + second_batch(kind) + survival_quantiles(kind) + algebra_laws(kind) + circular_laws(kind) + complement_laws(kind) + inverse_gaussian_laws(kind) + von_mises_distribution(kind) + skew_normal_laws(kind) + (new_discrete_laws() if kind == "double" else [])
        # tous les points de toutes les listes dans un seul calcul parallèle (les listes courtes
        # mais coûteuses, quantiles et constantes, tournaient sinon en série)
        tasks = [(i, x) for i, (_, _, _, xs, _) in enumerate(entries) for x in xs]
        values = iter(g.parallel_map(lambda i, v: entries[i][4](mpf(v)), tasks))
        for name, a, b, xs, f in entries:
            rows = groups.setdefault(name, [])
            for x, exact in zip(xs, [next(values) for _ in xs]):
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
    for name, cases, cdf, start, *low in (("poisson_quantile", [(mu, 0.0) for mu in P4_POISSON], poisson_cdf, poisson_start),
                                    ("binomial_quantile", [(float(n), pp) for n, pp in P4_BINOMIAL], binomial_cdf, binomial_start),
                                    ("negative_binomial_quantile", NEGATIVE_BINOMIAL, lambda r, pp, k: g.ibeta_reference(r, float(k) + 1, pp)[0],
                                     lambda r, pp, u: max(0.0, r * (1 - pp) / pp + _m.sqrt(r * (1 - pp)) / pp * float(-mp.sqrt(2) * g.erfcinv_exact(2 * mpf(u))))),
                                    ("beta_binomial_quantile", [(float(n), a) for n, a, _ in BETA_BINOMIAL if n <= 1000],
                                     lambda n, a, k: beta_binomial_cdf(int(n), a, BETA_BINOMIAL_BETA[(int(n), a)], int(k)), lambda n, a, u: n * u),
                                    ("logarithmic_quantile", [(pp, 0.0) for pp in LOGARITHMIC], lambda pp, _, k: mpf(0) if k < 1 else logarithmic_cdf(pp, int(k)), lambda pp, _, u: 1.0),
                                    ("discrete_gaussian_quantile", [(mu, sigma) for mu, sigma in DISCRETE_GAUSSIAN if sigma <= 30],
                                     lambda mu, sigma, k: discrete_gaussian_parts(mu, sigma, int(k))[1], lambda mu, sigma, u: mu + sigma * float(-mp.sqrt(2) * g.erfcinv_exact(2 * mpf(u))), -10 ** 6),
                                    ("hypergeometric_quantile", [(float(N), float(K)) for N, K, _ in HYPERGEOMETRIC],
                                     lambda N, K, k: hypergeometric_cdf_reference(int(N), int(K), int(k)),
                                     lambda N, K, u: hypergeometric_start(int(N), int(K), u))):
        rows = discrete_quantiles(name, cases, cdf, start, *low)
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
