#!/usr/bin/env python3
"""Valeurs de référence des lois de wsl::random (densités, répartitions, quantiles).

    uv run --no-project --with mpmath python tools/random/law_references.py FICHIER

Comme pour wsl::special : calculées avec mpmath à 256 bits, arrondies au plus proche en double ou en float, avec la position
exacte de la valeur entre deux flottants (en ulp). Non versionnées : CMake les
engendre dans le dossier de build. Déterministe (graine fixe).

Vérifie aussi les constantes en double-double de include/wsl/random/detail/arithmetic.hpp.
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
    text = (ROOT / "include/wsl/random/detail/arithmetic.hpp").read_text()
    exact = {
        "inv_sqrt2": 1 / mp.sqrt(2), "sqrt2": mp.sqrt(2), "inv_sqrt_2pi": 1 / mp.sqrt(2 * mp.pi),
        "half_log_2pi": mp.log(2 * mp.pi) / 2, "normal_entropy": (1 + mp.log(2 * mp.pi)) / 2,
        "half_sqrt_pi": mp.sqrt(mp.pi) / 2, "inv_sqrt_pi": 1 / mp.sqrt(mp.pi), "one_twelfth": mpf(1) / 12,
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


def reference_file():
    out = [HEADER, "#pragma once\n\n", "// Valeurs de référence des lois de wsl::random (mpmath).\n\n",
           "namespace law_reference_data\n{\n",
           "    template <class T>\n    struct point\n    {\n        double a, b, x;\n        T rounded;\n        double offset;\n    };\n\n"]
    total = 0
    for kind, rounder, ulp, fmt in (("double", g.round_double, g.ulp_double, g.hexd),
                                    ("float", g.round_float, g.ulp_float, g.hexf)):
        groups = {}
        for name, a, b, xs, f in laws(kind) + (discrete_laws() if kind == "double" else []):
            rows = groups.setdefault(name, [])
            for x in xs:
                exact = f(mpf(x))
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
    out.append("    };\n")
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
