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
    include/wsl/special/detail/trig_tables.hpp
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
LOG_BITS = 9                      # c = j / 512
LOG_FIRST, LOG_LAST = 362, 724    # z dans [√2/2, √2)  ⇒  j = round(512 z) dans [362, 724]
LOG_QUANTUM = 42                  # log(c)_hi multiple de 2^-42, comme e·ln2_hi


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

    # log : x = 2^e z, z = c (1 + r), c = j / 512. e ln2_hi doit être exact pour |e| < 2^11,
    # et e ln2_hi + log(c)_hi aussi : les deux sont des multiples de 2^-42 (ln2_hi a 42
    # bits significatifs, d'exposant -1), et leur somme reste sous 2^10.
    ln2_hi = truncated(ln2, 42)
    assert ln2_hi * 2 ** LOG_QUANTUM == mp.nint(ln2_hi * 2 ** LOG_QUANTUM)
    out += [f"    inline constexpr int log_table_bits = {LOG_BITS};\n",
            f"    inline constexpr int log_table_first = {LOG_FIRST};\n",
            f"    inline constexpr double log_ln2_hi = {hexd(float(ln2_hi))};     // ln 2, 42 bits\n",
            f"    inline constexpr double log_ln2_lo = {hexd(round_double(ln2 - ln2_hi))};\n\n",
            f"    // Pour c = j / {2 ** LOG_BITS}, j = {LOG_FIRST} … {LOG_LAST} : ln(c) = hi + lo, hi multiple de\n",
            f"    // 2^-{LOG_QUANTUM} (d'où e ln2_hi + hi exact), lo l'arrondi du reste ; puis 1/c arrondi.\n",
            "    inline constexpr double log_table[][3] = {\n"]
    for j in range(LOG_FIRST, LOG_LAST + 1):
        c = mpf(j) / 2 ** LOG_BITS
        exact = mp.log(c)
        hi = float(mp.nint(exact * 2 ** LOG_QUANTUM) / 2 ** LOG_QUANTUM)
        lo = round_double(exact - mpf(hi))
        out.append(f"        {{{hexd(hi)}, {hexd(lo)}, {hexd(round_double(1 / c))}}},   // j = {j}\n")
    out.append("    };\n\n")

    # Constantes en double-double : ln 2 (exp2), 1/ln 2 (log2), 1/ln 10 (log10).
    out.append("    // Constantes en double-double (hi + lo).\n")
    for name, value in (("ln2", ln2), ("inv_ln2", 1 / ln2), ("inv_ln10", 1 / mp.log(10))):
        hi, lo = double_double(value)
        out.append(f"    inline constexpr double {name}_hi = {hexd(hi)};\n")
        out.append(f"    inline constexpr double {name}_lo = {hexd(lo)};\n")
    out.append("}\n")
    return "".join(out)


TRIG_BITS = 6                     # tables aux points c = i / 64
SINCOS_LAST = 52                  # |r| <= π/4 + ε  ⇒  i = round(64 |r|) <= 51
ATAN_LAST = 64                    # v dans [0, 1]
TWO_OVER_PI_WORDS = 21            # bits de 2/π nécessaires jusqu'aux plus grands doubles


def trig_tables():
    out = [HEADER, "#pragma once\n\n#include <cstdint>\n\n",
           "// Constantes et tables de sin, cos, tan, atan, atan2 (voir trigonometric.hpp).\n\n",
           "namespace wsl::special::detail\n{\n"]
    half_pi = mp.pi / 2
    # π/2 en cinq morceaux de 33 bits : k · morceau est exact pour |k| < 2^20.
    parts, rest = [], half_pi
    for _ in range(5):
        part = truncated(rest, 33)
        parts.append(part)
        rest -= part
    out.append("    // π/2 = pio2_1 + … + pio2_5 à 2^-165 près, chaque morceau sur 33 bits.\n")
    for n, v in enumerate(parts, start=1):
        out.append(f"    inline constexpr double pio2_{n} = {hexd(float(v))};\n")
    out.append(f"    inline constexpr double inv_pio2 = {hexd(round_double(2 / mp.pi))};   // 2/π\n\n")
    out.append("    // Constantes en double-double, et arrondis de π/4 et 3π/4.\n")
    for name, v in (("pi", mp.pi), ("pio2", half_pi)):
        hi, lo = double_double(v)
        out.append(f"    inline constexpr double {name}_hi = {hexd(hi)};\n    inline constexpr double {name}_lo = {hexd(lo)};\n")
    out.append(f"    inline constexpr double pio4 = {hexd(round_double(mp.pi / 4))};\n")
    out.append(f"    inline constexpr double three_pio4 = {hexd(round_double(3 * mp.pi / 4))};\n\n")

    # 2/π en mots de 64 bits, précédés d'un mot nul : le bit i (i >= 1) de 2/π est le
    # bit (i − 1 + 64) de la suite, compté depuis le poids fort du premier mot.
    # À calculer bien au-delà des 256 bits de travail : la table en compte 1344.
    with mp.workprec(64 * TWO_OVER_PI_WORDS + 256):
        scaled = int(mp.floor((2 / mp.pi) * mpf(2) ** (64 * TWO_OVER_PI_WORDS)))
    words = [0] + [(scaled >> (64 * (TWO_OVER_PI_WORDS - 1 - w))) & (2 ** 64 - 1) for w in range(TWO_OVER_PI_WORDS)]
    out.append(f"    // Bits de 2/π, {64 * TWO_OVER_PI_WORDS} après la virgule, précédés d'un mot nul (Payne–Hanek).\n")
    out.append("    inline constexpr std::uint64_t two_over_pi_bits[] = {\n")
    for w in words:
        out.append(f"        0x{w:016x},\n")
    out.append("    };\n\n")

    out.append(f"    // Pour c = i / {2 ** TRIG_BITS}, i = 0 … {SINCOS_LAST} : sin(c) = hi + lo, cos(c) = hi + lo.\n")
    out.append("    inline constexpr double sincos_table[][4] = {\n")
    for i in range(SINCOS_LAST + 1):
        c = mpf(i) / 2 ** TRIG_BITS
        sh, sl = double_double(mp.sin(c))
        ch, cl = double_double(mp.cos(c))
        out.append(f"        {{{hexd(sh)}, {hexd(sl)}, {hexd(ch)}, {hexd(cl)}}},   // i = {i}\n")
    out.append("    };\n\n")
    out.append(f"    // Pour c = i / {2 ** TRIG_BITS}, i = 0 … {ATAN_LAST} : atan(c) = hi + lo.\n")
    out.append("    inline constexpr double atan_table[][2] = {\n")
    for i in range(ATAN_LAST + 1):
        hi, lo = double_double(mp.atan(mpf(i) / 2 ** TRIG_BITS))
        out.append(f"        {{{hexd(hi)}, {hexd(lo)}}},   // i = {i}\n")
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


def exp2_points_double():
    xs = [rng.uniform(-1080, 1024) for _ in range(1500)]
    xs += [math.copysign(2.0 ** rng.uniform(-60, 3), rng.choice((-1, 1))) for _ in range(600)]
    xs += [rng.uniform(-1075, -1022) for _ in range(300)]
    xs += [float(n) for n in range(-1074, 1024, 37)] + [n + 0.5 for n in range(-50, 50, 7)]
    return xs


def expm1_points_double():
    xs = [rng.uniform(-38, 709.78) for _ in range(800)]
    xs += [math.copysign(2.0 ** rng.uniform(-60, 0), rng.choice((-1, 1))) for _ in range(1200)]
    xs += [rng.uniform(-1, 1) for _ in range(600)]
    ln2 = math.log(2)
    for k in range(-6, 7):                                              # bords de sous-intervalles
        base = k * ln2 / EXP_N + ln2 / (2 * EXP_N)
        xs += [base, math.nextafter(base, math.inf), math.nextafter(base, -math.inf)]
    xs += [2.0 ** -1074, -(2.0 ** -1074), 1e-300, -1e-300, -37.9, 709.0, 1.0, -1.0]
    return xs


def log1p_points_double():
    xs = [math.copysign(2.0 ** rng.uniform(-60, -1), rng.choice((-1, 1))) for _ in range(1200)]
    xs += [rng.uniform(-1, 10) for _ in range(600)]
    xs += [-1 + 2.0 ** -rng.uniform(1, 53) for _ in range(300)]                      # près de −1
    xs += [random_positive_double() for _ in range(500)]
    xs += [2.0 ** -1074, -(2.0 ** -1074), 1e-300, -0.5, 1.0, 1e300, 1.7976931348623157e308]
    return [x for x in xs if x > -1]


def exp2_points_float():
    xs = [to_f32(rng.uniform(-150, 128)) for _ in range(800)]
    xs += [to_f32(math.copysign(2.0 ** rng.uniform(-30, 3), rng.choice((-1, 1)))) for _ in range(400)]
    xs += [float(n) for n in range(-149, 128, 11)]
    return xs


def expm1_points_float():
    xs = [to_f32(rng.uniform(-18, 88.72)) for _ in range(500)]
    xs += [to_f32(math.copysign(2.0 ** rng.uniform(-30, 0), rng.choice((-1, 1)))) for _ in range(800)]
    return xs


def log1p_points_float():
    xs = [to_f32(math.copysign(2.0 ** rng.uniform(-30, -1), rng.choice((-1, 1)))) for _ in range(800)]
    xs += [to_f32(rng.uniform(-1, 10)) for _ in range(300)]
    xs += [random_positive_float() for _ in range(300)]
    return [x for x in xs if x > -1]


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


def near_multiples_of_half_pi(count, k_max):
    xs = []
    for _ in range(count):
        k = rng.randint(1, k_max)
        base = round_double(k * mp.pi / 2)
        xs += [base, math.nextafter(base, math.inf), math.nextafter(base, 0)]
    return xs


def trig_points_double():
    xs = [rng.uniform(-10, 10) for _ in range(1200)]
    xs += [math.copysign(2.0 ** rng.uniform(-40, 0), rng.choice((-1, 1))) for _ in range(400)]
    xs += [rng.uniform(-2.0 ** 20, 2.0 ** 20) for _ in range(400)]
    large = []
    while len(large) < 600:
        x = random_positive_double()
        if x >= 2.0 ** 20:
            large.append(math.copysign(x, rng.choice((-1, 1))))
    xs += large
    xs += near_multiples_of_half_pi(150, 2 ** 20)
    xs += near_multiples_of_half_pi(50, 2 ** 60)
    xs += [6381956970095103 * 2.0 ** 797, 1e22, 1e300, 1.7976931348623157e308, math.pi / 4,
           math.nextafter(math.pi / 4, 1), 2.0 ** 20, math.nextafter(2.0 ** 20, 0), 0.5, 1.0, 3.0]
    return xs


def trig_points_float():
    xs = [to_f32(rng.uniform(-10, 10)) for _ in range(600)]
    xs += [to_f32(math.copysign(2.0 ** rng.uniform(-20, 0), rng.choice((-1, 1)))) for _ in range(200)]
    xs += [math.copysign(random_positive_float(), rng.choice((-1, 1))) for _ in range(400)]
    return xs


def atan_points_double():
    xs = [math.copysign(2.0 ** rng.uniform(-60, 60), rng.choice((-1, 1))) for _ in range(1500)]
    xs += [rng.uniform(-2, 2) for _ in range(500)]
    xs += [2.0 ** -1074, 1e-300, 1.0, -1.0, 1e300, 1.7976931348623157e308, 2.0 ** 60, 2.0 ** 61]
    return xs


def atan_points_float():
    xs = [to_f32(math.copysign(2.0 ** rng.uniform(-30, 30), rng.choice((-1, 1)))) for _ in range(800)]
    xs += [to_f32(rng.uniform(-2, 2)) for _ in range(300)]
    return xs


def atan2_points_double():
    pts = []
    for _ in range(1500):
        y = math.copysign(2.0 ** rng.uniform(-60, 60) * rng.uniform(1, 2), rng.choice((-1, 1)))
        x = math.copysign(2.0 ** rng.uniform(-60, 60) * rng.uniform(1, 2), rng.choice((-1, 1)))
        pts.append((y, x))
    for _ in range(300):                                    # rapports extrêmes
        pts.append((math.copysign(random_positive_double(), rng.choice((-1, 1))),
                    math.copysign(random_positive_double(), rng.choice((-1, 1)))))
    for _ in range(200):                                    # |y| ≈ |x|
        x = math.copysign(rng.uniform(0.1, 10), rng.choice((-1, 1)))
        pts.append((math.copysign(abs(x) * (1 + rng.uniform(-1e-6, 1e-6)), rng.choice((-1, 1))), x))
    return pts


def atan2_points_float():
    pts = []
    for _ in range(800):
        y = math.copysign(to_f32(2.0 ** rng.uniform(-30, 30)), rng.choice((-1, 1)))
        x = math.copysign(to_f32(2.0 ** rng.uniform(-30, 30)), rng.choice((-1, 1)))
        pts.append((y, x))
    return pts


def pow_points_double():
    pts = []
    for _ in range(1500):                                   # x quelconque, |y log x| < 700
        x = random_positive_double()
        lx = abs(math.log(x))
        if lx > 0:
            pts.append((x, rng.uniform(-1, 1) * 700 / lx))
    for _ in range(800):                                    # x près de 1, y grand
        x = 1 + math.copysign(2.0 ** -rng.uniform(1, 52), rng.choice((-1, 1)))
        pts.append((x, rng.uniform(-1, 1) * 700 / abs(math.log(x))))
    for _ in range(600):                                    # y entier, x quelconque de signe quelconque
        x = math.copysign(math.exp(rng.uniform(-5, 5)), rng.choice((-1, 1)))
        n = rng.randint(-int(700 / abs(math.log(abs(x)) or 1)) - 1, int(700 / abs(math.log(abs(x)) or 1)) + 1)
        n = max(-1000, min(1000, n))
        if n != 0 and abs(n * math.log(abs(x))) < 700:
            pts.append((x, float(n)))
    for y in (0.5, -0.5, 1 / 3, 2.0, 3.0, -1.0, -2.0, 1.5, 0.1):
        for _ in range(60):
            pts.append((math.exp(rng.uniform(-200, 200)), y))
    for _ in range(300):                                    # résultats sous-normaux
        x = math.exp(rng.uniform(-50, 50))
        if abs(math.log(x)) > 1e-3:
            pts.append((x, rng.uniform(-745, -708.5) / math.log(x)))
    pts += [(2.0, float(n)) for n in range(-1074, 1024, 41)]
    pts += [(10.0, float(n)) for n in range(-20, 23)] + [(4.0, 0.5), (8.0, 1 / 3), (0.25, -0.5)]
    return pts


def pow_points_float():
    pts = []
    for _ in range(800):
        x = random_positive_float()
        lx = abs(math.log(x))
        if lx > 0:
            pts.append((x, to_f32(rng.uniform(-1, 1) * 85 / lx)))
    for _ in range(300):
        x = math.copysign(to_f32(math.exp(rng.uniform(-3, 3))), rng.choice((-1, 1)))
        pts.append((x, float(rng.randint(-20, 20) or 1)))
    for _ in range(300):
        x = to_f32(1 + math.copysign(2.0 ** -rng.uniform(1, 23), rng.choice((-1, 1))))
        if x != 1:
            pts.append((x, to_f32(rng.uniform(-1, 1) * 85 / abs(math.log(x)))))
    return pts


def pow_exact(x, y):
    mx, my = mpf(x), mpf(y)
    if mx < 0:                                              # y entier (garanti par les points)
        v = mp.power(-mx, my)
        return -v if int(y) % 2 else v
    return mp.power(mx, my)


def references2(name, pts, f, rounder, ulp, fmt, ctype):
    lines = [f"    // {name} : x, y, arrondi correct de f(x, y), position exacte en ulp\n",
             f"    inline constexpr reference2<{ctype}> {name}[] = {{\n"]
    kept = 0
    for x, y in sorted(set(pts)):
        exact = f(x, y)
        cr = rounder(exact)
        if math.isinf(cr) or cr == 0:
            continue
        offset = float((exact - mpf(cr)) / mpf(ulp(cr)))
        lines.append(f"        {{{fmt(x)}, {fmt(y)}, {fmt(cr)}, {offset!r}}},\n")
        kept += 1
    lines.append("    };\n\n")
    return "".join(lines), kept


def reference_file():
    out = [HEADER, "#pragma once\n\n",
           "// Valeurs de référence de wsl::special, calculées par mpmath à 256 bits.\n\n",
           "namespace reference_data\n{\n",
           "    template <class T>\n    struct reference\n    {\n        T x, rounded;\n        double offset;\n    };\n\n",
           "    template <class T>\n    struct reference2\n    {\n        T x, y, rounded;\n        double offset;\n    };\n\n"]
    counts = {}
    for name, xs, f, rounder, ulp, fmt, ctype in (
            ("exp_double", exp_points_double(), mp.exp, round_double, ulp_double, hexd, "double"),
            ("log_double", log_points_double(), mp.log, round_double, ulp_double, hexd, "double"),
            ("exp_float", exp_points_float(), mp.exp, round_float, ulp_float, hexf, "float"),
            ("log_float", log_points_float(), mp.log, round_float, ulp_float, hexf, "float"),
            ("exp2_double", exp2_points_double(), lambda x: mp.power(2, x), round_double, ulp_double, hexd, "double"),
            ("expm1_double", expm1_points_double(), mp.expm1, round_double, ulp_double, hexd, "double"),
            ("log2_double", log_points_double(), lambda x: mp.log(x, 2), round_double, ulp_double, hexd, "double"),
            ("log10_double", log_points_double(), mp.log10, round_double, ulp_double, hexd, "double"),
            ("log1p_double", log1p_points_double(), mp.log1p, round_double, ulp_double, hexd, "double"),
            ("exp2_float", exp2_points_float(), lambda x: mp.power(2, x), round_float, ulp_float, hexf, "float"),
            ("expm1_float", expm1_points_float(), mp.expm1, round_float, ulp_float, hexf, "float"),
            ("log2_float", log_points_float(), lambda x: mp.log(x, 2), round_float, ulp_float, hexf, "float"),
            ("log10_float", log_points_float(), mp.log10, round_float, ulp_float, hexf, "float"),
            ("log1p_float", log1p_points_float(), mp.log1p, round_float, ulp_float, hexf, "float"),
            ("sin_double", trig_points_double(), mp.sin, round_double, ulp_double, hexd, "double"),
            ("cos_double", trig_points_double(), mp.cos, round_double, ulp_double, hexd, "double"),
            ("tan_double", trig_points_double(), mp.tan, round_double, ulp_double, hexd, "double"),
            ("atan_double", atan_points_double(), mp.atan, round_double, ulp_double, hexd, "double"),
            ("sin_float", trig_points_float(), mp.sin, round_float, ulp_float, hexf, "float"),
            ("cos_float", trig_points_float(), mp.cos, round_float, ulp_float, hexf, "float"),
            ("tan_float", trig_points_float(), mp.tan, round_float, ulp_float, hexf, "float"),
            ("atan_float", atan_points_float(), mp.atan, round_float, ulp_float, hexf, "float")):
        text, counts[name] = references(name, xs, f, rounder, ulp, fmt, ctype)
        out.append(text)
    atan2_exact = lambda y, x: mp.atan2(mpf(y), mpf(x))
    for name, pts, f, rounder, ulp, fmt, ctype in (
            ("pow_double", pow_points_double(), pow_exact, round_double, ulp_double, hexd, "double"),
            ("pow_float", pow_points_float(), pow_exact, round_float, ulp_float, hexf, "float"),
            ("atan2_double", atan2_points_double(), atan2_exact, round_double, ulp_double, hexd, "double"),
            ("atan2_float", atan2_points_float(), atan2_exact, round_float, ulp_float, hexf, "float")):
        text, counts[name] = references2(name, pts, f, rounder, ulp, fmt, ctype)
        out.append(text)
    out.append("}\n")
    return "".join(out), counts


def main():
    (ROOT / "include/wsl/special/detail/tables.hpp").write_text(tables())
    (ROOT / "include/wsl/special/detail/trig_tables.hpp").write_text(trig_tables())
    text, counts = reference_file()
    (ROOT / "tests/data/special_elementary_reference.hpp").write_text(text)
    for name, n in counts.items():
        print(f"{name:12s} {n:5d} points")


if __name__ == "__main__":
    main()
