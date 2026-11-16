#!/usr/bin/env python3
"""Engendre les tables de wsl::special et les valeurs de référence de ses tests.

    uv run --with mpmath python tools/special/generate.py

Tout est calculé avec mpmath à 256 bits, puis arrondi au plus proche (pair) en
double ou en float : les constantes des en-têtes sont donc les arrondis corrects
des valeurs exactes, et les références des tests portent, en plus de l'arrondi
correct, la position exacte de la valeur entre deux flottants (en ulp), pour
mesurer l'erreur de nos fonctions au millième d'ulp près.

Le script est déterministe (graine fixe) : le relancer redonne les mêmes fichiers.

Fichiers engendrés (à ne pas modifier à la main) :
    include/wsl/special/detail/tables.hpp
    tests/data/special_elementary_reference.hpp
"""

import math
import random
import struct
from pathlib import Path

from mpmath import mp, mpf

mp.prec = 256
ROOT = Path(__file__).resolve().parents[2]
HEADER = "// ENGENDRÉ par tools/special/generate.py — ne pas modifier à la main.\n"


# ---------------------------------------------------------------------------
# Arrondis corrects
# ---------------------------------------------------------------------------

def round_double(v):
    """Arrondi au plus proche (pair) de v en double ; ±inf au-delà."""
    if v == 0:
        return 0.0
    limit = mpf(2) ** 1024 * (1 - mpf(2) ** -54)          # DBL_MAX + ulp/2
    if abs(v) >= limit:
        return math.copysign(math.inf, float(v))
    guess = float(v)
    best = None
    for c in (math.nextafter(guess, -math.inf), guess, math.nextafter(guess, math.inf)):
        if math.isinf(c):
            continue
        d = abs(mpf(c) - v)
        if best is None or d < best[0] or (d == best[0] and _even_double(c)):
            best = (d, c)
    return best[1]


def _even_double(c):
    return struct.unpack("<Q", struct.pack("<d", c))[0] % 2 == 0


def to_f32(x):
    return struct.unpack("<f", struct.pack("<f", x))[0]


def f32_bits(f):
    return struct.unpack("<I", struct.pack("<f", f))[0]


def f32_from_bits(b):
    return struct.unpack("<f", struct.pack("<I", b))[0]


def next_f32(f, up):
    """Float voisin de f (f fini), vers le haut ou vers le bas."""
    if f == 0:
        return f32_from_bits(1) if up else -f32_from_bits(1)
    b = f32_bits(f)
    if (f > 0) == up:
        b += 1
    else:
        b -= 1
    return f32_from_bits(b)


FLT_MAX = f32_from_bits(0x7F7FFFFF)


def round_float(v):
    """Arrondi au plus proche (pair) de v en float ; ±inf au-delà."""
    if v == 0:
        return 0.0
    limit = mpf(2) ** 128 * (1 - mpf(2) ** -25)           # FLT_MAX + ulp/2
    if abs(v) >= limit:
        return math.copysign(math.inf, float(v))
    guess = to_f32(float(v)) if abs(float(v)) <= FLT_MAX else math.copysign(FLT_MAX, float(v))
    best = None
    for c in (next_f32(guess, False), guess, next_f32(guess, True)):
        if math.isinf(c):
            continue
        d = abs(mpf(c) - v)
        if best is None or d < best[0] or (d == best[0] and f32_bits(c) % 2 == 0):
            best = (d, c)
    return best[1]


def ulp_double(c):
    return math.ulp(c)


def ulp_float(c):
    a = abs(c)
    return next_f32(a, True) - a


def double_double(v):
    hi = round_double(v)
    return hi, round_double(v - mpf(hi))


def truncated(v, bits_kept):
    """v arrondi à `bits_kept` bits significatifs (valeur exacte, en mpf)."""
    e = int(mp.floor(mp.log(abs(v), 2)))
    scale = mpf(2) ** (bits_kept - 1 - e)
    return mp.nint(v * scale) / scale


def hexd(x):
    if x == 0:
        return "-0.0" if math.copysign(1, x) < 0 else "0.0"
    return float.hex(x)


def hexf(x):
    return hexd(x) + "f" if x != 0 else ("-0.0f" if math.copysign(1, x) < 0 else "0.0f")


# ---------------------------------------------------------------------------
# Tables des en-têtes
# ---------------------------------------------------------------------------

EXP_BITS = 7                      # 2^7 = 128 sous-intervalles
EXP_N = 1 << EXP_BITS
LOG_BITS = 8                      # c = j / 256
LOG_FIRST, LOG_LAST = 181, 362    # z dans [√2/2, √2]  ⇒  j = round(256 z) dans [181, 362]


def tables():
    ln2 = mp.log(2)
    out = [HEADER, "#pragma once\n\n",
           "// Constantes et tables de exp et log (voir elementary.hpp pour leur usage).\n",
           "// Chaque valeur est l'arrondi correct de la valeur exacte ; les paires\n",
           "// (hi, lo) sont des double-double : hi + lo vaut la valeur à 2^-106 près.\n\n",
           "namespace wsl::special::detail\n{\n"]

    # exp : x = k ln2/N + r. k ln2/N doit être exact pour |k| < 2^18 : hi sur 35 bits.
    ln2_n_hi = truncated(ln2 / EXP_N, 35)
    out += [f"    inline constexpr int exp_table_bits = {EXP_BITS};\n",
            f"    inline constexpr double exp_inv_ln2_n = {hexd(round_double(EXP_N / ln2))};   // N / ln 2\n",
            f"    inline constexpr double exp_ln2_n_hi = {hexd(float(ln2_n_hi))};     // ln 2 / N, 35 bits\n",
            f"    inline constexpr double exp_ln2_n_lo = {hexd(round_double(ln2 / EXP_N - ln2_n_hi))};\n\n",
            "    // 2^(j/N) pour j = 0 … N-1.\n",
            "    inline constexpr double exp_table[][2] = {\n"]
    for j in range(EXP_N):
        hi, lo = double_double(mpf(2) ** (mpf(j) / EXP_N))
        out.append(f"        {{{hexd(hi)}, {hexd(lo)}}},\n")
    out.append("    };\n\n")

    # log : x = 2^e z, z = c (1 + r), c = j / 256. e ln2_hi doit être exact pour |e| < 2^11.
    ln2_hi = truncated(ln2, 42)
    out += [f"    inline constexpr int log_table_bits = {LOG_BITS};\n",
            f"    inline constexpr int log_table_first = {LOG_FIRST};\n",
            f"    inline constexpr double log_ln2_hi = {hexd(float(ln2_hi))};     // ln 2, 42 bits\n",
            f"    inline constexpr double log_ln2_lo = {hexd(round_double(ln2 - ln2_hi))};\n\n",
            f"    // Pour c = j / 256, j = {LOG_FIRST} … {LOG_LAST} : ln(c) en double-double, et 1/c arrondi.\n",
            "    inline constexpr double log_table[][3] = {\n"]
    for j in range(LOG_FIRST, LOG_LAST + 1):
        c = mpf(j) / 2 ** LOG_BITS
        hi, lo = double_double(mp.log(c))
        out.append(f"        {{{hexd(hi)}, {hexd(lo)}, {hexd(round_double(1 / c))}}},   // j = {j}\n")
    out.append("    };\n}\n")
    return "".join(out)


# ---------------------------------------------------------------------------
# Points de référence
# ---------------------------------------------------------------------------

rng = random.Random(20261116)


def random_positive_double():
    while True:
        b = rng.getrandbits(63)
        if 0 < b < 0x7FF0000000000000:
            return struct.unpack("<d", struct.pack("<Q", b))[0]


def random_positive_float():
    while True:
        b = rng.getrandbits(31)
        if 0 < b < 0x7F800000:
            return f32_from_bits(b)


def exp_points_double():
    xs = [rng.uniform(-745.13, 709.78) for _ in range(3000)]
    xs += [math.copysign(2.0 ** rng.uniform(-60, 9.4), rng.choice((-1, 1))) for _ in range(1500)]
    xs += [rng.uniform(-745.13, -708.4) for _ in range(500)]          # résultats sous-normaux
    xs += [rng.uniform(-1, 1) * 2 ** -10 for _ in range(500)]
    ln2 = math.log(2)
    for k in range(-20, 21):                                            # bords de sous-intervalles
        base = k * ln2 / EXP_N + ln2 / (2 * EXP_N)
        xs += [base, math.nextafter(base, math.inf), math.nextafter(base, -math.inf)]
    xs += [1.0, -1.0, 0.5, 2.0, 10.0, -10.0, 100.0, -100.0, 700.0, -700.0,
           float.fromhex("0x1.62e42fefa39efp+9"), float.fromhex("-0x1.6232bdd7abcd2p+9"),
           2.0 ** -1074, -(2.0 ** -1074), 2.0 ** -1022, 1e-300, -1e-300]
    return xs


def log_points_double():
    xs = [random_positive_double() for _ in range(3000)]
    xs += [1 + math.copysign(2.0 ** rng.uniform(-60, -1), rng.choice((-1, 1))) for _ in range(1500)]
    xs += [rng.uniform(0.5, 2.0) for _ in range(1000)]
    for j in range(LOG_FIRST, LOG_LAST + 1, 3):                        # bords de sous-intervalles
        base = (j + 0.5) / 2 ** LOG_BITS
        xs += [base, math.nextafter(base, 0), math.nextafter(base, 2)]
    sqrt2 = math.sqrt(2)
    xs += [sqrt2, math.nextafter(sqrt2, 0), math.nextafter(sqrt2, 2), sqrt2 / 2,
           math.nextafter(1.0, 2), math.nextafter(1.0, 0), 2.0, 0.5, math.e, 10.0,
           2.0 ** -1074, 2.0 ** -1060, 2.0 ** -1022, math.nextafter(2.0 ** -1022, 0),
           1.7976931348623157e308, 1e-300, 1e300]
    return [x for x in xs if x != 1.0]


def exp_points_float():
    xs = [to_f32(rng.uniform(-103.27, 88.72)) for _ in range(2000)]
    xs += [to_f32(math.copysign(2.0 ** rng.uniform(-30, 6.4), rng.choice((-1, 1)))) for _ in range(800)]
    xs += [to_f32(rng.uniform(-103.27, -87.33)) for _ in range(300)]
    xs += [to_f32(v) for v in (1.0, -1.0, 0.5, 10.0, -10.0, 88.72, -87.0)]
    return xs


def log_points_float():
    xs = [random_positive_float() for _ in range(2000)]
    xs += [to_f32(1 + math.copysign(2.0 ** rng.uniform(-24, -1), rng.choice((-1, 1)))) for _ in range(800)]
    xs += [to_f32(rng.uniform(0.5, 2.0)) for _ in range(500)]
    xs += [f32_from_bits(1), f32_from_bits(0x00800000), FLT_MAX, to_f32(math.sqrt(2)), to_f32(math.e)]
    return [x for x in xs if x != 1.0]


def references(name, xs, f, rounder, ulp, fmt, ctype):
    lines = [f"    // {name} : x, arrondi correct de f(x), position exacte de f(x) en ulp\n",
             f"    inline constexpr reference<{ctype}> {name}[] = {{\n"]
    kept = 0
    for x in sorted(set(xs)):
        exact = f(mpf(x))
        cr = rounder(exact)
        if math.isinf(cr) or cr == 0:
            continue
        offset = float((exact - mpf(cr)) / mpf(ulp(cr)))
        lines.append(f"        {{{fmt(x)}, {fmt(cr)}, {offset!r}}},\n")
        kept += 1
    lines.append("    };\n\n")
    return "".join(lines), kept


def reference_file():
    out = [HEADER, "#pragma once\n\n",
           "// Valeurs de référence de wsl::special, calculées par mpmath à 256 bits.\n\n",
           "namespace reference_data\n{\n",
           "    template <class T>\n    struct reference\n    {\n        T x, rounded;\n        double offset;\n    };\n\n"]
    counts = {}
    for name, xs, f, rounder, ulp, fmt, ctype in (
            ("exp_double", exp_points_double(), mp.exp, round_double, ulp_double, hexd, "double"),
            ("log_double", log_points_double(), mp.log, round_double, ulp_double, hexd, "double"),
            ("exp_float", exp_points_float(), mp.exp, round_float, ulp_float, hexf, "float"),
            ("log_float", log_points_float(), mp.log, round_float, ulp_float, hexf, "float")):
        text, counts[name] = references(name, xs, f, rounder, ulp, fmt, ctype)
        out.append(text)
    out.append("}\n")
    return "".join(out), counts


def main():
    (ROOT / "include/wsl/special/detail/tables.hpp").write_text(tables())
    text, counts = reference_file()
    (ROOT / "tests/data/special_elementary_reference.hpp").write_text(text)
    for name, n in counts.items():
        print(f"{name:12s} {n:5d} points")


if __name__ == "__main__":
    main()
