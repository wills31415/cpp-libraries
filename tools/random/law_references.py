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
        "half_sqrt_pi": mp.sqrt(mp.pi) / 2, "inv_sqrt_pi": 1 / mp.sqrt(mp.pi),
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
    return out


def reference_file():
    out = [HEADER, "#pragma once\n\n", "// Valeurs de référence des lois de wsl::random (mpmath).\n\n",
           "namespace law_reference_data\n{\n",
           "    template <class T>\n    struct point\n    {\n        double a, b, x;\n        T rounded;\n        double offset;\n    };\n\n"]
    total = 0
    for kind, rounder, ulp, fmt in (("double", g.round_double, g.ulp_double, g.hexd),
                                    ("float", g.round_float, g.ulp_float, g.hexf)):
        groups = {}
        for name, a, b, xs, f in laws(kind):
            rows = groups.setdefault(name, [])
            for x in xs:
                exact = f(mpf(x))
                cr = rounder(exact)
                if math.isinf(cr) or cr == 0:
                    continue
                offset = float((exact - mpf(cr)) / mpf(ulp(cr)))
                rows.append(f"        {{{g.hexd(a)}, {g.hexd(b)}, {g.hexd(x)}, {fmt(cr)}, {offset!r}}},\n")
        for name, rows in groups.items():
            out.append(f"    inline constexpr point<{kind}> {name}_{kind}[] = {{\n")
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
