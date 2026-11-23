"""Outils communs des références des lois de wsl::random (une famille de lois par module).

Chaque module engendre un en-tête data/random_law_reference/<famille>.hpp, que CMake ne
recalcule que si le module (ou ce fichier, ou generate.py) change. Valeurs calculées avec
mpmath à 256 bits, arrondies au plus proche en double ou en float, avec la position exacte
de la valeur entre deux flottants (en ulp). Déterministe : graine fixe par famille.
"""

import math
import random
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/special"))
import generate as g                                         # noqa: E402  (arrondis, inverse d'erfc)
from mpmath import mp, mpf                                   # noqa: E402

mp.prec = 256
rng = random.Random(20261117)                                # réensemencé par famille dans run()

HEADER = "// ENGENDRÉ par tools/random/references/{}.py — ne pas modifier à la main.\n"


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


def run(family, entries, tables=lambda: [], check=None):
    """Écrit l'en-tête de la famille : points {a, b, x, arrondi, position} par loi et par type
    (entries(kind), kind = "double" ou "float"), puis les tables de quantiles exacts
    (tables() : (type de ligne, nom, lignes))."""
    if len(sys.argv) != 2:
        sys.exit(f"usage : law_references.py {family} FICHIER")
    rng.seed(f"wsl::random/{family}")
    if check:
        check()
    out = [HEADER.format(family), "#pragma once\n\n", "namespace law_reference_data\n{\n"]
    total = 0
    for kind, rounder, ulp, fmt in (("double", g.round_double, g.ulp_double, g.hexd),
                                    ("float", g.round_float, g.ulp_float, g.hexf)):
        groups = {}
        listed = entries(kind)
        # tous les points de la famille dans un seul calcul parallèle
        tasks = [(i, x) for i, (_, _, _, xs, _) in enumerate(listed) for x in xs]
        values = iter(g.parallel_map(lambda i, v: listed[i][4](mpf(v)), tasks))
        for name, a, b, xs, f in listed:
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
    for struct, name, rows in tables():
        out.append(f"    // Quantiles exacts.\n    inline constexpr {struct} {name}[] = {{\n")
        out += rows
        out.append("    };\n\n")
        total += len(rows)
    out.append("}\n")
    output = Path(sys.argv[1])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(out))
    print(f"{total} points de référence ({family}) : {output.name}")


def normal_start(u):
    """Φ⁻¹(u) en double (départ des recherches de quantiles discrets)."""
    return float(-mp.sqrt(2) * g.erfcinv_exact(2 * mpf(u)))
