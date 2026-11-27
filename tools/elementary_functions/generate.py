#!/usr/bin/env python3
"""Engendre les tables de wsl::special et les valeurs de référence de ses tests.

    uv run --no-project --with mpmath python tools/elementary_functions/generate.py --tables
    uv run --no-project --with mpmath python tools/elementary_functions/generate.py --references FICHIER

--tables réécrit les tables des en-têtes (versionnées : la bibliothèque en a
besoin). --references écrit les valeurs de référence des tests dans FICHIER ;
elles ne sont pas versionnées : CMake les engendre dans le dossier de build quand
il construit les tests.

Tout est calculé avec mpmath à 256 bits, puis arrondi au plus proche (pair) en
double ou en float : les constantes des en-têtes sont donc les arrondis corrects
des valeurs exactes, et les références des tests portent, en plus de l'arrondi
correct, la position exacte de la valeur entre deux flottants (en ulp), pour
mesurer l'erreur de nos fonctions au millième d'ulp près.

Le script est déterministe (graine fixe) : le relancer redonne les mêmes fichiers.

Fichiers engendrés (à ne pas modifier à la main) :
    include/wsl/elementary_functions/detail/exponential_tables.hpp
    include/wsl/elementary_functions/detail/trigonometric_tables.hpp
    include/wsl/special_functions/detail/erf_tables.hpp
    include/wsl/special_functions/detail/gamma_tables.hpp
    include/wsl/special_functions/detail/digamma_tables.hpp
    include/wsl/special_functions/detail/incomplete_gamma_tables.hpp
"""

import math
import multiprocessing
import os
from fractions import Fraction
import random
import struct
from pathlib import Path

from mpmath import mp, mpf

mp.prec = 256
ROOT = Path(__file__).resolve().parents[2]
HEADER = "// ENGENDRÉ par tools/elementary_functions/generate.py — ne pas modifier à la main.\n"


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


def to_f32_safe(x):
    return to_f32(x) if abs(x) < 3e38 else math.inf


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
    out = [HEADER, "#pragma once\n\n#include <wsl/modes.hpp>\n\n",
           "// Constantes et tables de exp et log (voir elementary.hpp pour leur usage).\n",
           "// Chaque valeur est l'arrondi correct de la valeur exacte ; les paires\n",
           "// (hi, lo) sont des double-double : hi + lo vaut la valeur à 2^-106 près.\n\n",
           "namespace wsl::inline WSL_DEFAULTS::elementary_functions::detail\n{\n"]

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
    out = [HEADER, "#pragma once\n\n#include <wsl/modes.hpp>\n\n#include <cstdint>\n\n",
           "// Constantes et tables de sin, cos, tan, atan, atan2 (voir trigonometric.hpp).\n\n",
           "namespace wsl::inline WSL_DEFAULTS::elementary_functions::detail\n{\n"]
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


ERF_SMALL = 0.25                  # |x| < 1/4 : série ; au-delà, erfc(x) = e^{-x²} F(x)
ERF_SMALL_TERMS = 10              # termes de la série (reste < 2^-66)
ERFC_DEGREE = 15                  # degré de Taylor de F par morceau
ERFC_SPLIT = 4.0                  # morceaux de largeur 1/8 sous 4, 1/2 au-delà
ERFC_END = 28.0                   # erfc(x) < 2^-1100 au-delà : 0


def erfc_pieces():
    """Centres et demi-largeurs des morceaux de F, dans l'ordre de la table."""
    pieces = []
    x = ERF_SMALL
    while x < ERFC_SPLIT:
        pieces.append((x + 1 / 16, 1 / 16))
        x += 1 / 8
    while x < ERFC_END:
        pieces.append((x + 1 / 4, 1 / 4))
        x += 1 / 2
    return pieces


def erfc_factor(x):
    """F(x) = e^{x²} erfc(x), lisse et lentement variable."""
    return mp.exp(x * x) * mp.erfc(x)


def erfc_factor_taylor(c, degree):
    """Coefficients de Taylor de F en c, par F' = 2xF − 2/√π : (k+1) f_{k+1} = 2c f_k + 2 f_{k−1}.
    La récurrence amplifie l'erreur de f_0 d'un facteur e^{c²} (jusqu'à 2^1131) : on
    la conduit à 1 700 bits."""
    with mp.workprec(1700):
        c = mpf(c)
        f = [erfc_factor(c)]
        f.append(2 * c * f[0] - 2 / mp.sqrt(mp.pi))
        for k in range(1, degree):
            f.append((2 * c * f[k] + 2 * f[k - 1]) / (k + 1))
        return [+v for v in f]


def erf_tables():
    out = [HEADER, "#pragma once\n\n#include <wsl/modes.hpp>\n\n",
           "// Constantes et tables de erf et erfc (voir error_function.hpp).\n\n",
           "namespace wsl::inline WSL_DEFAULTS::special_functions::detail\n{\n"]
    hi, lo = double_double(2 / mp.sqrt(mp.pi))
    out.append(f"    inline constexpr double two_over_sqrt_pi_hi = {hexd(hi)};\n")
    out.append(f"    inline constexpr double two_over_sqrt_pi_lo = {hexd(lo)};\n\n")
    out.append(f"    // erf(x) = (2/√π) x (1 + x² Q(x²)) pour |x| < {ERF_SMALL} : Q(u) = Σ q_n u^(n−1),\n")
    out.append("    // q_n = (−1)^n / (n! (2n + 1)).\n")
    out.append("    inline constexpr double erf_small_q[] = {\n")
    for n in range(1, ERF_SMALL_TERMS + 1):
        q = mpf((-1) ** n) / (mp.factorial(n) * (2 * n + 1))
        out.append(f"        {hexd(round_double(q))},\n")
    out.append("    };\n")
    q1_hi, q1_lo = double_double(mpf(-1) / 3)
    out.append(f"    inline constexpr double erf_small_q1_lo = {hexd(q1_lo)};      // −1/3 − q_1\n\n")
    pieces = erfc_pieces()
    out.append(f"    inline constexpr int erfc_degree = {ERFC_DEGREE};\n")
    out.append(f"    inline constexpr int erfc_fine_pieces = {sum(1 for c, h in pieces if c < ERFC_SPLIT)};\n\n")
    out.append("    // F(x) = e^{x²} erfc(x) par morceaux : F(c + t) = Σ a_k t^k. Par ligne :\n")
    out.append("    // a_0 (hi, lo), a_1 (hi, lo), a_2 … a_15. Morceaux de largeur 1/8 de 1/4 à 4,\n")
    out.append("    // puis de largeur 1/2 jusqu'à 28 ; c au milieu.\n")
    out.append("    inline constexpr double erfc_table[][4 + erfc_degree - 1] = {\n")
    for c, half in pieces:
        a = erfc_factor_taylor(c, ERFC_DEGREE)
        # Reste de Taylor vérifié sur l'intervalle.
        worst = 0
        for step in range(-16, 17):
            t = mpf(half) * step / 16
            approx = sum(a[k] * t ** k for k in range(ERFC_DEGREE + 1))
            worst = max(worst, abs(approx / erfc_factor(mpf(c) + t) - 1))
        assert worst < mpf(2) ** -66, (c, worst)
        a0, a0l = double_double(a[0])
        a1, a1l = double_double(a[1])
        row = [a0, a0l, a1, a1l] + [round_double(v) for v in a[2:]]
        out.append("        {" + ", ".join(hexd(v) for v in row) + f"}},   // c = {c}\n")
    out.append("    };\n\n")
    out.append(erf_inverse_tables())
    out.append("}\n")
    return "".join(out)


ERFINV_CENTRAL_DEGREE = 14
ERFINV_TAIL_DEGREE = 12
ERFINV_TAIL_CUTS = [None, 3.2, 4.5, 7.0, 12.0, 20.0, None]   # bornes en s = √w ; extrêmes calculés


def erfc_inverse_exact(z):
    """erfc⁻¹(z), z dans (0, 1), par Newton sur log erfc (précis même pour z sous-normal)."""
    z = mpf(z)
    return mp.findroot(lambda x: mp.log(mp.erfc(x)) - mp.log(z), mp.sqrt(-mp.log(z)))


def erf_inverse_tables():
    """Approximation initiale de erf⁻¹(y) = y P, à 2^-35 près : P(w) au centre, P(√w) dans les queues,
    w = −log((1 − y)(1 + y)), chaque polynôme ajusté (Tchebychev) en u dans [−1, 1]."""
    out = []
    hi, lo = double_double(mp.sqrt(mp.pi) / 2)
    out.append(f"    inline constexpr double half_sqrt_pi_hi = {hexd(hi)};\n")
    out.append(f"    inline constexpr double half_sqrt_pi_lo = {hexd(lo)};\n\n")
    with mp.workprec(120):
        def central(u):
            w = (u + 1) * mpf(5) / 2
            y = mp.sqrt(-mp.expm1(-w))
            return mp.erfinv(y) / y if y != 0 else mp.sqrt(mp.pi) / 2
        coefficients, error = mp.chebyfit(central, [mpf(-1), mpf(1)], ERFINV_CENTRAL_DEGREE + 1, error=True)
        assert error < mpf(2) ** -34, error
        out.append("    // Centre : w dans [0, 5], u = (w − 5/2)/(5/2) ; coefficients du degré le plus haut.\n")
        out.append("    inline constexpr double erfinv_central[] = {\n")
        out += [f"        {hexd(round_double(c))},\n" for c in coefficients]
        out.append("    };\n\n")

        def tail(s):
            q = mp.exp(-s * s)
            z = q / (1 + mp.sqrt(1 - q))
            return erfc_inverse_exact(z) / (1 - z)
        cuts = list(ERFINV_TAIL_CUTS)
        cuts[0] = mp.sqrt(5)
        cuts[-1] = mp.sqrt(mp.log(2) * 1075)
        out.append("    // Queues : s = √w dans [√5, √(1075 ln 2)], six morceaux [a, b], u = (s − m)/h.\n")
        out.append("    inline constexpr double erfinv_tail_bounds[] = {")
        out.append(", ".join(hexd(round_double(c)) for c in cuts) + "};\n")
        out.append("    inline constexpr double erfinv_tail[][%d] = {\n" % (ERFINV_TAIL_DEGREE + 1))
        for a, b in zip(cuts, cuts[1:]):
            a, b = mpf(round_double(a)), mpf(round_double(b))
            mid, half = (a + b) / 2, (b - a) / 2
            coefficients, error = mp.chebyfit(lambda u: tail(mid + half * u), [mpf(-1), mpf(1)],
                                               ERFINV_TAIL_DEGREE + 1, error=True)
            assert error / abs(tail(a)) < mpf(2) ** -35, (a, b, error)
            out.append("        {" + ", ".join(hexd(round_double(c)) for c in coefficients) + "},\n")
        out.append("    };\n")
    return "".join(out)


LGAMMA_FIRST, LGAMMA_LAST = 4, 80   # morceaux centrés en c = k/8 : [0,4375 ; 10,0625)
LGAMMA_DEGREE = 26
LGAMMA_WIDE = {8: (7, 9), 16: (15, 17)}  # zéros 1 et 2 : leur morceau couvre aussi les voisins
STIRLING_TERMS = 10


def gamma_tables():
    out = [HEADER, "#pragma once\n\n#include <wsl/modes.hpp>\n\n",
           "// Constantes et tables de lgamma et tgamma (voir gamma.hpp).\n\n",
           "namespace wsl::inline WSL_DEFAULTS::special_functions::detail\n{\n"]
    for name, v in (("log_pi", mp.log(mp.pi)), ("half_log_two_pi", mp.log(2 * mp.pi) / 2)):
        hi, lo = double_double(v)
        out.append(f"    inline constexpr double {name}_hi = {hexd(hi)};\n    inline constexpr double {name}_lo = {hexd(lo)};\n")
    out.append("\n    // Stirling : lgamma(x) = (x − 1/2) log x − x + log(2π)/2 + Σ B_2k / (2k (2k − 1) x^(2k−1)).\n")
    out.append("    inline constexpr double stirling[] = {\n")
    for k in range(1, STIRLING_TERMS + 1):
        out.append(f"        {hexd(round_double(mp.bernoulli(2 * k) / (2 * k * (2 * k - 1))))},\n")
    out.append("    };\n\n")
    out.append(f"    inline constexpr int lgamma_first = {LGAMMA_FIRST};\n")
    out.append(f"    inline constexpr int lgamma_degree = {LGAMMA_DEGREE};\n\n")
    out.append("    // lgamma(c + t) = Σ a_k t^k pour c = k/8, |t| <= 1/16 : a_0 = lgamma(c), a_1 = ψ(c),\n")
    out.append("    // a_k = ψ^(k−1)(c)/k!. Par ligne : a_0, a_1, a_2 (hi, lo chacun), a_3 … a_26. Les centres\n")
    out.append("    // 1 et 2, zéros de lgamma, donnent a_0 = 0 exactement ; leurs morceaux s'étendent à\n")
    out.append("    // |t| <= 3/16 (k = 7 à 9 et 15 à 17) : près d'un zéro, partir d'un centre voisin\n")
    out.append("    // annulerait a_0 contre a_1 t et grossirait les erreurs d'arrondi.\n")
    out.append("    inline constexpr double lgamma_table[][6 + lgamma_degree - 2] = {\n")
    for k in range(LGAMMA_FIRST, LGAMMA_LAST + 1):
        c = mpf(k) / 8
        a = [mp.loggamma(c), mp.digamma(c)] + [mp.polygamma(n - 1, c) / mp.factorial(n) for n in range(2, LGAMMA_DEGREE + 1)]
        half = mpf(3) / 16 if k in LGAMMA_WIDE else mpf(1) / 16
        worst = 0
        for step in range(-16, 17):
            t = half * step / 16
            exact = mp.loggamma(c + t)
            if exact != 0:
                approx = sum(a[n] * t ** n for n in range(LGAMMA_DEGREE + 1))
                worst = max(worst, abs(approx - exact) / abs(exact))
        assert worst < mpf(2) ** -66, (k, worst)
        a0, a0l = double_double(a[0])
        a1, a1l = double_double(a[1])
        a2, a2l = double_double(a[2])
        row = [a0, a0l, a1, a1l, a2, a2l] + [round_double(v) for v in a[3:]]
        out.append("        {" + ", ".join(hexd(v) for v in row) + f"}},   // c = {k}/8\n")
    out.append("    };\n}\n")
    return "".join(out)


DIGAMMA_ROOT_DEGREE = 24
DIGAMMA_ROOT_HALF_WIDTH = mpf(1) / 8
DIGAMMA_ASYMPTOTIC_FROM = 10
DIGAMMA_ASYMPTOTIC_TERMS = 14


def digamma_tables():
    out = [HEADER, "#pragma once\n\n#include <wsl/modes.hpp>\n\n",
           "// Constantes et tables de digamma (voir gamma.hpp).\n\n",
           "namespace wsl::inline WSL_DEFAULTS::special_functions::detail\n{\n"]
    root = mp.findroot(mp.digamma, mpf("1.4616321449683623"))
    r1 = round_double(root)
    r2 = round_double(root - mpf(r1))
    r3 = round_double(root - mpf(r1) - mpf(r2))
    out.append("    // x₀, zéro positif de ψ, en trois parties.\n")
    out.append(f"    inline constexpr double digamma_root_hi = {hexd(r1)};\n")
    out.append(f"    inline constexpr double digamma_root_lo = {hexd(r2)};\n")
    out.append(f"    inline constexpr double digamma_root_lo2 = {hexd(r3)};\n\n")
    gamma_hi, gamma_lo = double_double(mp.euler)
    out.append(f"    inline constexpr double euler_gamma_hi = {hexd(gamma_hi)};\n")
    out.append(f"    inline constexpr double euler_gamma_lo = {hexd(gamma_lo)};\n\n")
    # Taylor autour de x₀ : ψ(x₀ + t) = Σ_{k>=1} c_k t^k, c_k = ψ^(k)(x₀)/k!.
    c = [None] + [mp.polygamma(k, root) / mp.factorial(k) for k in range(1, DIGAMMA_ROOT_DEGREE + 1)]
    worst = 0
    for step in range(-32, 33):
        t = DIGAMMA_ROOT_HALF_WIDTH * step / 32
        if t == 0:
            continue
        exact = mp.digamma(root + t)
        approx = sum(c[k] * t ** k for k in range(1, DIGAMMA_ROOT_DEGREE + 1))
        worst = max(worst, abs(approx - exact) / abs(exact))
    assert worst < mpf(2) ** -66, worst
    out.append(f"    inline constexpr int digamma_root_degree = {DIGAMMA_ROOT_DEGREE};\n")
    out.append(f"    inline constexpr double digamma_root_half_width = {hexd(float(DIGAMMA_ROOT_HALF_WIDTH))};\n")
    out.append("    // ψ(x₀ + t) = Σ c_k t^k, |t| <= 1/8 (reste < 2^-66 relativement) : c_1, c_2 en\n")
    out.append("    // double-double (hi, lo), puis c_3 … c_24.\n")
    row = list(double_double(c[1])) + list(double_double(c[2])) + [round_double(v) for v in c[3:]]
    out.append("    inline constexpr double digamma_root_taylor[] = {\n")
    out += [f"        {hexd(v)},\n" for v in row]
    out.append("    };\n\n")
    # Développement asymptotique : ψ(y) = log y − 1/(2y) − Σ_{k>=1} B_2k / (2k y^(2k)).
    y = mpf(DIGAMMA_ASYMPTOTIC_FROM)
    k = DIGAMMA_ASYMPTOTIC_TERMS + 1
    rest = abs(mp.bernoulli(2 * k) / (2 * k) / y ** (2 * k)) / mp.digamma(y)
    assert rest < mpf(2) ** -66, rest
    out.append(f"    inline constexpr double digamma_asymptotic_from = {hexd(float(y))};\n")
    out.append("    // B_2k / (2k), k = 1 … ; ψ(y) = log y − 1/(2y) − Σ c_k y^(−2k) pour y >= 10\n")
    out.append("    // (reste < 2^-66 relativement).\n")
    out.append("    inline constexpr double digamma_asymptotic[] = {\n")
    out += [f"        {hexd(round_double(mp.bernoulli(2 * j) / (2 * j)))},\n" for j in range(1, DIGAMMA_ASYMPTOTIC_TERMS + 1)]
    out.append("    };\n}\n")
    return "".join(out)


# Gamma incomplète : développement uniforme de Temme (DLMF 8.12) et lgamma(1 + a).
TEMME_FROM = 50                                 # a >= 50 et |η| <= 1
TEMME_ETA = 1
TEMME_TERMS = 12                                # c_0 … c_12 (reste < 2^-78 mesuré, a >= 50)
TEMME_SERIES_DEGREE = 60                        # degré des séries exactes en η
LGAMMA1P_TERMS = 40                             # lgamma(1 + a) = Σ_{k>=1} l_k a^k, 0 < a < 1/4


def _series_mul(a, b, n):
    r = [Fraction(0)] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]):
                r[i + j] += x * y
    return r


def _series_inv(a, n):
    r = [Fraction(0)] * (n + 1)
    r[0] = 1 / a[0]
    for k in range(1, n + 1):
        r[k] = -sum(a[j] * r[k - j] for j in range(1, min(k, len(a) - 1) + 1)) / a[0]
    return r


def temme_coefficients():
    """c_k(η) = Σ_n d_kn η^n, en rationnels exacts (DLMF 8.12.8–8.12.10)."""
    M = TEMME_SERIES_DEGREE
    # η = μ √(2 f(μ)/μ²), μ = λ − 1, f(μ) = μ − ln(1 + μ) ; série g(μ)/μ, puis réversion μ = η m(η).
    h = [Fraction(2 * (-1) ** j, j) for j in range(2, M + 3)]
    gm = [Fraction(0)] * (M + 1)
    gm[0] = Fraction(1)
    for k in range(1, M + 1):
        gm[k] = (h[k] - sum(gm[j] * gm[k - j] for j in range(1, k))) / 2
    m = [Fraction(1)] + [Fraction(0)] * M
    for _ in range(M + 1):
        mu = [Fraction(0)] + m[:M]
        comp = [Fraction(0)] * (M + 1)
        power = [Fraction(1)] + [Fraction(0)] * M
        for k in range(M + 1):
            if gm[k]:
                for i in range(M + 1):
                    comp[i] += gm[k] * power[i]
            power = _series_mul(power, mu, M)
        new = _series_inv(comp, M)
        if new == m:
            break
        m = new
    # g_k : Γ(z) ~ e^{−z} z^z √(2π/z) Σ g_k z^{−k}, Σ g_k t^k = exp(Σ_j B_2j t^(2j−1) / (2j (2j − 1))).
    K = TEMME_TERMS
    st = [Fraction(0)] * (K + 1)
    for j in range(1, K // 2 + 2):
        if 2 * j - 1 <= K:
            num, den = mp.bernfrac(2 * j)
            st[2 * j - 1] = Fraction(int(num), int(den)) / (2 * j * (2 * j - 1))
    g = [Fraction(1)] + [Fraction(0)] * K
    ds = [(k + 1) * st[k + 1] for k in range(K)]
    for k in range(1, K + 1):
        g[k] = sum(ds[j] * g[k - 1 - j] for j in range(k)) / k
    im = _series_inv(m, M)
    c = [[im[k + 1] for k in range(M)]]                     # c_0 = (1/m − 1)/η
    for k in range(1, K + 1):
        prev = c[-1]
        d = [(i + 1) * prev[i + 1] for i in range(len(prev) - 1)]
        num = [d[i] + (-1) ** k * g[k] * im[i] for i in range(len(d))]
        assert num[0] == 0, k
        c.append(num[1:])
    return c


def incomplete_gamma_tables():
    out = [HEADER, "#pragma once\n\n#include <wsl/modes.hpp>\n\n",
           "// Tables de la gamma incomplète (voir incomplete_gamma.hpp).\n\n",
           "namespace wsl::inline WSL_DEFAULTS::special_functions::detail\n{\n"]
    # lgamma(1 + a) = −γ a + Σ_{k>=2} (−1)^k ζ(k) a^k / k
    coefficients = [-mp.euler] + [(-1) ** k * mp.zeta(k) / k for k in range(2, LGAMMA1P_TERMS + 1)]
    for a in (mpf(1) / 4, mpf(1) / 1000, mpf(10) ** -30):
        approx = sum(cf * a ** (k + 1) for k, cf in enumerate(coefficients))
        assert abs(approx - mp.loggamma(1 + a)) / abs(mp.loggamma(1 + a)) < mpf(2) ** -70, a
    out.append("    // lgamma(1 + a) = Σ_{k>=1} l_k a^k pour 0 < a <= 1/4 (reste < 2^-70 relativement),\n")
    out.append("    // l_1 = −γ, l_k = (−1)^k ζ(k)/k ; en double-double (hi, lo).\n")
    out.append(f"    inline constexpr double lgamma1p_series[{LGAMMA1P_TERMS}][2] = {{\n")
    for cf in coefficients:
        hi, lo = double_double(cf)
        out.append(f"        {{{hexd(hi)}, {hexd(lo)}}},\n")
    out.append("    };\n\n")
    # Temme
    c = temme_coefficients()
    def value(fr):
        return mpf(fr.numerator) / fr.denominator
    degrees = []
    for k, ck in enumerate(c):
        # degré : reste Σ_{n>N} |d_kn| (|η| <= 1) · a^-k sous 2^-74 pour a >= TEMME_FROM
        scale = mpf(TEMME_FROM) ** -k
        n = len(ck) - 1
        tail = 0
        while n > 0 and (tail + abs(value(ck[n]))) * scale < mpf(2) ** -74:
            tail += abs(value(ck[n]))
            n -= 1
        assert n < len(ck) - 1 or k == 0 and False, (k, n)
        degrees.append(n)
    width = max(degrees) + 1
    out.append(f"    inline constexpr double temme_from = {hexd(float(TEMME_FROM))};\n")
    out.append(f"    inline constexpr double temme_eta = {hexd(float(TEMME_ETA))};\n")
    out.append(f"    inline constexpr int temme_terms = {TEMME_TERMS + 1};\n")
    out.append(f"    inline constexpr int temme_width = {width};\n")
    out.append("    // Développement uniforme de Temme (DLMF 8.12) :\n")
    out.append("    //     Q(a, x) = erfc(η √(a/2))/2 + e^{−aη²/2} / √(2πa) Σ_k c_k(η) a^{−k},\n")
    out.append("    //     η = signe(λ − 1) √(2 (λ − 1 − log λ)), λ = x/a,\n")
    out.append("    // c_k(η) = Σ_n d_kn η^n (rationnels exacts arrondis), tronqués pour |η| <= 1, a >= 50\n")
    out.append("    // (reste < 2^-74). Degrés par k, puis d_0n en double-double et d_kn (k >= 1).\n")
    out.append("    inline constexpr int temme_degree[] = {" + ", ".join(str(d) for d in degrees) + "};\n")
    out.append(f"    inline constexpr double temme_c0[{degrees[0] + 1}][2] = {{\n")
    for v in c[0][:degrees[0] + 1]:
        hi, lo = double_double(value(v))
        out.append(f"        {{{hexd(hi)}, {hexd(lo)}}},\n")
    out.append("    };\n")
    out.append(f"    inline constexpr double temme_c[{TEMME_TERMS}][{width}] = {{\n")
    for k in range(1, TEMME_TERMS + 1):
        row = [round_double(value(v)) for v in c[k][:degrees[k] + 1]] + [0.0] * (width - degrees[k] - 1)
        out.append("        {" + ", ".join(hexd(v) for v in row) + "},\n")
    out.append("    };\n}\n")
    return "".join(out)


# ---------------------------------------------------------------------------
# Points de référence
# ---------------------------------------------------------------------------

rng = random.Random(20261116)


# Bessel modifiée I_ν : développement uniforme de Debye (DLMF 10.41.3), ν >= BESSEL_DEBYE_FROM.
BESSEL_DEBYE_FROM = 50
BESSEL_DEBYE_TERMS = 28                         # u_0 … u_27 (reste < 2^-112 vérifié, ν >= 50)
BESSEL_HANKEL_FROM = 60                         # Hankel pour x >= max(60, ν²/2), ν < 50


def debye_polynomials(count):
    """u_k(t), DLMF 10.41.10 : u_{k+1} = ½ t²(1 − t²) u_k' + ⅛ ∫₀ᵗ (1 − 5s²) u_k(s) ds, en rationnels."""
    from fractions import Fraction

    def deriv(q):
        return {e - 1: c * e for e, c in q.items() if e > 0}

    def mul(a, b):
        r = {}
        for e1, c1 in a.items():
            for e2, c2 in b.items():
                r[e1 + e2] = r.get(e1 + e2, 0) + c1 * c2
        return r

    u = [{0: Fraction(1)}]
    while len(u) < count + 1:
        r = mul({2: Fraction(1, 2), 4: Fraction(-1, 2)}, deriv(u[-1]))
        for e, c in mul({0: Fraction(1), 2: Fraction(-5)}, u[-1]).items():
            r[e + 1] = r.get(e + 1, 0) + c / (e + 1) / 8
        u.append({e: c for e, c in r.items() if c != 0})
    return u


def bessel_tables():
    u = debye_polynomials(BESSEL_DEBYE_TERMS)
    worst = max(abs(sum(mpf(c.numerator) / c.denominator * mpf(t) ** e for e, c in u[BESSEL_DEBYE_TERMS].items()))
                for t in [i / mpf(400) for i in range(401)])
    rest = worst / mpf(BESSEL_DEBYE_FROM) ** BESSEL_DEBYE_TERMS
    assert rest < mpf(2) ** -112, rest
    out = [HEADER, "#pragma once\n\n#include <wsl/modes.hpp>\n\n",
           "// Tables de la fonction de Bessel modifiée I_ν (voir bessel.hpp).\n\n",
           "namespace wsl::inline WSL_DEFAULTS::special_functions::detail\n{\n",
           f"    inline constexpr double bessel_debye_from = {BESSEL_DEBYE_FROM};\n",
           f"    inline constexpr int bessel_debye_terms = {BESSEL_DEBYE_TERMS};\n",
           f"    inline constexpr double bessel_hankel_from = {BESSEL_HANKEL_FROM};\n\n",
           "    // Polynômes de Debye u_k(t) = Σ_j c_kj t^{k + 2j}, j = 0 … k, k = 0 … terms − 1 :\n",
           "    // c_kj (hi, lo) en k(k + 1)/2 + j. Reste après u_{terms − 1} < 2^-112 pour ν >= 50.\n",
           "    inline constexpr double bessel_debye[][2] = {\n"]
    for k in range(BESSEL_DEBYE_TERMS):
        for j in range(k + 1):
            c = u[k].get(k + 2 * j, 0)
            hi, lo = double_double(mpf(c.numerator) / c.denominator) if c != 0 else (0.0, 0.0)
            out.append(f"        {{{hexd(hi)}, {hexd(lo)}}},\n")
    out.append("    };\n}\n")
    return "".join(out)


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


_parallel_f = None


def _parallel_call(args):
    # Tuple brut : dépicler un mpf l'arrondirait à la précision courante (256 bits), or les
    # références à 320 ou 400 bits en ont besoin de tous.
    return _parallel_f(*args)._mpf_


def parallel_map(f, argsets):
    """[f(*a) for a in argsets], sur tous les cœurs (fork : f peut être une fonction locale).
    Résultats identiques au calcul séquentiel, dans le même ordre."""
    global _parallel_f
    argsets = list(argsets)
    if len(argsets) < 64:
        return [mpf(f(*a)) for a in argsets]
    _parallel_f = f
    with multiprocessing.get_context("fork").Pool() as pool:
        raw = pool.map(_parallel_call, argsets, chunksize=1)               # points triés : coûts groupés
    return [mp.make_mpf(v) for v in raw]


def references(name, xs, f, rounder, ulp, fmt, ctype):
    lines = [f"    // {name} : x, arrondi correct de f(x), position exacte de f(x) en ulp\n",
             f"    inline constexpr reference<{ctype}> {name}[] = {{\n"]
    kept = 0
    xs = sorted(set(xs))
    for x, exact in zip(xs, parallel_map(lambda v: f(mpf(v)), [(x,) for x in xs])):
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


def real_cbrt(x):
    """Racine cubique réelle (mp.cbrt rend la racine complexe principale pour x < 0)."""
    return mp.sign(x) * mp.cbrt(abs(x))


def cbrt_points_double():
    xs = [math.copysign(random_positive_double(), rng.choice((-1, 1))) for _ in range(1500)]
    xs += [rng.uniform(-10, 10) for _ in range(500)]
    xs += [float(n ** 3) for n in range(1, 200, 7)] + [2.0 ** -1074, 2.0 ** -1060, 1.7976931348623157e308]
    return [x for x in xs if x != 0]


def cbrt_points_float():
    xs = [math.copysign(random_positive_float(), rng.choice((-1, 1))) for _ in range(1000)]
    return xs


def hypot_points_double():
    pts = []
    for _ in range(1200):
        a = math.copysign(random_positive_double(), rng.choice((-1, 1)))
        pts.append((a, math.copysign(abs(a) * 2.0 ** rng.uniform(-60, 60), rng.choice((-1, 1)))))
    for _ in range(400):
        pts.append((rng.uniform(-10, 10), rng.uniform(-10, 10)))
    for _ in range(200):                                    # près du débordement et des sous-normaux
        pts.append((1.7976931348623157e308 * rng.uniform(0.3, 0.7), 1.7976931348623157e308 * rng.uniform(0.3, 0.7)))
        pts.append((2.0 ** -1074 * rng.randint(1, 2 ** 40), 2.0 ** -1074 * rng.randint(1, 2 ** 40)))
    pts += [(3.0, 4.0), (5.0, 12.0), (1.0, 1.0)]
    return [(x, y) for x, y in pts if not math.isinf(x) and not math.isinf(y) and abs(y) < 1.7976931348623157e308]


def hypot_points_float():
    pts = []
    for _ in range(800):
        a = math.copysign(random_positive_float(), rng.choice((-1, 1)))
        y = math.copysign(abs(a) * 2.0 ** rng.uniform(-30, 30), rng.choice((-1, 1)))
        if abs(y) < FLT_MAX:                                # sinon, pas un float
            pts.append((a, to_f32(y)))
    return pts


def erf_points_double():
    xs = [math.copysign(2.0 ** rng.uniform(-40, -1), rng.choice((-1, 1))) for _ in range(600)]
    xs += [rng.uniform(-6, 6) for _ in range(1500)]
    for c, half in erfc_pieces()[:46]:                      # bords de morceaux
        for b in (c - half, c + half):
            xs += [b, math.nextafter(b, 0), math.nextafter(b, 10), -b]
    xs += [2.0 ** -1074, -(2.0 ** -1074), 2.0 ** -1022, 1e-300, 0.25, math.nextafter(0.25, 0), 5.9, 6.0]
    return xs


def erfc_points_double():
    xs = [rng.uniform(-6, 28) for _ in range(1500)]
    xs += [rng.uniform(26, 27.3) for _ in range(300)]       # résultats sous-normaux
    xs += [rng.uniform(-0.5, 0.5) for _ in range(300)]
    xs += [math.copysign(2.0 ** rng.uniform(-40, -1), rng.choice((-1, 1))) for _ in range(300)]
    for c, half in erfc_pieces():
        for b in (c - half, c + half):
            xs += [b, math.nextafter(b, 0), math.nextafter(b, 100)]
    xs += [2.0 ** -1074, -6.0, -5.9, 0.25, -0.25, 27.0, 27.2]
    return xs


def erf_points_float():
    xs = [to_f32(math.copysign(2.0 ** rng.uniform(-20, -1), rng.choice((-1, 1)))) for _ in range(300)]
    xs += [to_f32(rng.uniform(-4.5, 4.5)) for _ in range(800)]
    return xs


def erfc_points_float():
    xs = [to_f32(rng.uniform(-4.5, 10.05)) for _ in range(1000)]
    xs += [to_f32(math.copysign(2.0 ** rng.uniform(-20, -1), rng.choice((-1, 1)))) for _ in range(200)]
    return xs


def log_abs_gamma(x):
    return mp.re(mp.loggamma(x))                            # log |Γ(x)|, réel aussi pour x < 0


def lgamma_positive_points_double():
    xs = [2.0 ** rng.uniform(-60, 1000) for _ in range(1000)]
    xs += [rng.uniform(0, 12) for _ in range(1200)]
    for base in (1.0, 2.0):
        xs += [base + math.copysign(2.0 ** -rng.uniform(1, 52), rng.choice((-1, 1))) for _ in range(150)]
    for k in range(LGAMMA_FIRST, LGAMMA_LAST + 1, 3):       # bords de morceaux
        b = (k + 0.5) / 8
        xs += [b, math.nextafter(b, 0), math.nextafter(b, 100)]
    xs += [2.0 ** -1074, 1e-300, 0.5, 1.5, 2.5, 10.0625, 1e300, 2.5599833278516383e305]
    return [x for x in xs if x not in (1.0, 2.0)]


def lgamma_negative_points_double():
    xs = []
    while len(xs) < 1000:
        x = rng.uniform(-30, 0)
        if x != round(x):
            xs.append(x)
    xs += [-(2.0 ** rng.uniform(-60, 50)) for _ in range(300)]
    # Au plus près des zéros de lgamma sur les négatifs : la réflexion y retranche
    # deux termes presque égaux (zone mesurée exprès, voir design/special.md).
    for start in (-2.457, -2.747, -3.143, -3.955, -4.039, -4.999, -5.0001, -6.0001, -9.0001):
        zero = mp.findroot(lambda t: mp.re(mp.loggamma(t)), mpf(start))
        base = round_double(zero)
        xs += [base, math.nextafter(base, 0), math.nextafter(base, -100)]
        xs += [base + 2.0 ** -e for e in (20, 30, 40, 45)] + [base - 2.0 ** -e for e in (20, 30, 40, 45)]
    return [x for x in xs if x != round(x)]


def tgamma_points_double():
    xs = [rng.uniform(0, 171.62) for _ in range(1200)]
    xs += [2.0 ** rng.uniform(-60, 0) for _ in range(300)]
    xs += [float(n) for n in range(1, 171)] + [n + 0.5 for n in range(0, 170, 3)]
    neg = []
    while len(neg) < 800:
        x = rng.uniform(-184, 0)
        if x != round(x):
            neg.append(x)
    xs += neg
    for n in range(1, 180, 7):                               # près des pôles
        xs += [-n + 2.0 ** -rng.uniform(1, 40), -n - 2.0 ** -rng.uniform(1, 40)]
    return xs


def lgamma_points_float():
    xs = [to_f32(2.0 ** rng.uniform(-30, 120)) for _ in range(600)]
    xs += [to_f32(rng.uniform(0, 12)) for _ in range(400)]
    return [x for x in xs if x not in (1.0, 2.0)]


digamma_rng = random.Random(20261118)                       # à part : ne décale pas les autres points


def digamma_points_double():
    r = digamma_rng
    root = float(mp.findroot(mp.digamma, mpf("1.4616321449683623")))
    xs = [r.uniform(0, 30) for _ in range(800)] + [2.0 ** r.uniform(-1000, 1000) for _ in range(400)]
    xs += [root + r.uniform(-0.2, 0.2) for _ in range(200)] + [root * (1 + r.choice([-1, 1]) * 2.0 ** -r.uniform(10, 53)) for _ in range(100)]
    xs += [10.0, math.nextafter(10.0, 0), 1.0, 2.0, 0.5, root, math.nextafter(root, 0), math.nextafter(root, 2)]
    return [x for x in xs if x > 0]


def digamma_negative_points_double():
    r = digamma_rng
    xs = []
    while len(xs) < 1000:
        x = r.uniform(-30, 0)
        if x != round(x):
            xs.append(x)
    xs += [-(2.0 ** r.uniform(-60, 50)) for _ in range(300)]
    for start in (-0.5, -1.6, -2.6, -3.6, -9.7):                      # près des zéros de ψ
        zero = mp.findroot(mp.digamma, mpf(start))
        base = round_double(zero)
        xs += [base, math.nextafter(base, 0), math.nextafter(base, -100)]
        xs += [base + 2.0 ** -e for e in (20, 30, 40, 45)] + [base - 2.0 ** -e for e in (20, 30, 40, 45)]
    return [x for x in xs if x != round(x) and x > -2.0 ** 52]


def digamma_points_float():
    r = digamma_rng
    xs = [to_f32(r.uniform(0, 30)) for _ in range(500)] + [to_f32(2.0 ** r.uniform(-120, 120)) for _ in range(300)]
    return [x for x in xs if x > 0]


igamma_rng = random.Random(20261119)


def igamma_points(count, to_type=lambda v: v):
    r = igamma_rng
    pts = []
    for i in range(count):
        kind = i % 6
        if kind == 0:                                                   # a quelconque, x autour de a
            a = 2.0 ** r.uniform(-30, 20)
            x = a * 2.0 ** r.uniform(-4, 4)
        elif kind == 1:                                                 # près de x = a (Temme, transition)
            a = 2.0 ** r.uniform(0, 20)
            x = a + r.uniform(-6, 6) * math.sqrt(a)
        elif kind == 2:                                                 # petits a
            a = 2.0 ** r.uniform(-60, -2)
            x = 2.0 ** r.uniform(-20, 3)
        elif kind == 3:                                                 # bord de Temme (a ≈ 50, |η| ≈ 1)
            a = r.uniform(40, 70)
            x = a * r.choice([r.uniform(0.25, 0.35), r.uniform(2.2, 2.5), r.uniform(0.9, 1.1)])
        elif kind == 4:                                                 # queues lointaines
            a = 2.0 ** r.uniform(-5, 10)
            x = a * 2.0 ** r.choice([r.uniform(-60, -5), r.uniform(3, 7)])
        else:                                                           # a modéré, x quelconque
            a = r.uniform(0.25, 60)
            x = r.uniform(0, 3 * a + 10)
        a, x = to_type(a), to_type(x)
        if a > 0 and x > 0:
            pts.append((a, x))
    return pts


_igamma_cache = {}


def igamma_reference(a, x):
    """(P, Q) à 320 bits : série pour P si x < a + 1, fraction continue de Legendre (Lentz)
    pour Q sinon ; l'autre vaut 1 − elle, sans perte notable à 320 bits pour a >= 2^-60.
    (mp.gammainc ne converge pas pour a de l'ordre de 10^6.)"""
    with mp.workprec(320):
        a, x = mpf(a), mpf(x)
        tiny = mpf(2) ** -330
        if x < a + 1:                                   # la fraction continue converge mal en deçà
            term = s = mpf(1)
            n = 1
            while True:
                term *= x / (a + n)
                s += term
                if term < s * tiny:
                    break
                n += 1
            p_ = mp.exp(a * mp.log(x) - x - mp.loggamma(a + 1)) * s
            return p_, 1 - p_
        f = c = x + 1 - a
        d = mpf(0)
        n = 1
        while True:
            an = -n * (n - a)
            bn = x + 2 * n + 1 - a
            d = 1 / (bn + an * d)
            c = bn + an / c
            delta = c * d
            f *= delta
            if abs(delta - 1) < tiny:
                break
            n += 1
        q_ = mp.exp(a * mp.log(x) - x - mp.loggamma(a)) / f
        return 1 - q_, q_


def igamma_exact(a, x, upper):
    key = (a, x)
    if key not in _igamma_cache:
        _igamma_cache[key] = igamma_reference(a, x)
    return _igamma_cache[key][1 if upper else 0]


def igamma_inverse_exact(a, t, upper):
    """x tel que P(a, x) = t (Q si upper), à 320 bits : bissection en u = log x, puis Newton
    (dL/du = x f(x)/F, f la densité)."""
    with mp.workprec(320):
        a, t = mpf(a), mpf(t)
        lt = mp.log(t)
        def residual(u):
            p_, q_ = igamma_reference(a, mp.exp(u))
            v = q_ if upper else p_
            return (mp.log(v) - lt) if v > 0 else -mpf(10) ** 9
        sign = -1 if upper else 1                                      # L croît avec u côté P
        # Départ : Wilson–Hilferty (a >= 1), sinon les queues ; encadrement élargi au besoin.
        pside = (1 - t) if upper else t
        z = -mp.sqrt(2) * erfcinv_exact(2 * pside)
        w = 1 - 1 / (9 * a) + z / (3 * mp.sqrt(a))
        if a >= 1 and w > 0:
            u0 = mp.log(a * w ** 3)
        elif not upper:
            u0 = (lt + mp.loggamma(a + 1)) / a
        else:
            u0 = mp.log(max(a, mpf(1)) - lt)
        lo, hi, width = u0 - 1, u0 + 1, mpf(2)
        while sign * residual(lo) > 0:
            lo -= width
            width *= 2
        width = mpf(2)
        while sign * residual(hi) < 0:
            hi += width
            width *= 2
        for _ in range(200):
            if hi - lo < mpf(2) ** -8:
                break
            mid = (lo + hi) / 2
            if sign * residual(mid) > 0:
                hi = mid
            else:
                lo = mid
        u = (lo + hi) / 2
        for _ in range(60):
            x = mp.exp(u)
            p_, q_ = igamma_reference(a, x)
            v = q_ if upper else p_
            g = mp.log(v) - lt
            density = mp.exp((a - 1) * mp.log(x) - x - mp.loggamma(a))
            step = g / (sign * x * density / v)
            u -= step
            if abs(step) < mpf(2) ** -300:
                break
        return mp.exp(u)


def igamma_inverse_points(count, to_type=lambda v: v):
    r = igamma_rng
    pts = []
    for i in range(count):
        a = to_type(2.0 ** r.uniform(-30, 20) if i % 3 else r.uniform(0.2, 80))
        u = r.random()
        t = r.random() if u < 0.4 else 2.0 ** -r.uniform(1, 1000) if u < 0.8 else 0.5 - 2.0 ** -r.uniform(2, 50)
        t = to_type(t)
        if a > 0 and 0 < t < 1:
            pts.append((a, t))
    return pts


ibeta_rng = random.Random(20261120)
_ibeta_cache = {}


def ibeta_reference(a, b, x):
    """(I, 1 − I) à 400 bits : fraction continue de Iₓ(a, b) (Lentz) si x < (a + 1)/(a + b + 2),
    sinon celle de I_{1−x}(b, a) ; l'autre vaut 1 − elle."""
    key = (a, b, x)
    if key in _ibeta_cache:
        return _ibeta_cache[key]
    with mp.workprec(400):
        a, b, x = mpf(a), mpf(b), mpf(x)
        y = 1 - x
        prefix = mp.exp(a * mp.log(x) + b * mp.log(y) - (mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)))
        def fraction(a, b, x):
            tiny, eps = mpf(2) ** -1000, mpf(2) ** -410
            c, d = mpf(1), 1 - (a + b) * x / (a + 1)
            d = 1 / (d if abs(d) > tiny else tiny)
            h = d
            m = 1
            while True:
                for num in (m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m)),
                            -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))):
                    d = 1 + num * d
                    d = 1 / (d if abs(d) > tiny else tiny)
                    c = 1 + num / c
                    c = c if abs(c) > tiny else tiny
                    delta = c * d
                    h *= delta
                if abs(delta - 1) < eps:
                    return h
                m += 1
        if x < (a + 1) / (a + b + 2):
            i = prefix * fraction(a, b, x) / a
            result = (i, 1 - i)
            derived = 1                                        # le côté obtenu par 1 − v
        else:
            j = prefix * fraction(b, a, y) / b
            result = (1 - j, j)
            derived = 0
        # Côté obtenu par 1 − v (paramètre minuscule) et trop petit pour 400 bits : directement.
        if result[derived] < mpf(2) ** -200:
            if derived == 0:
                result = (mp.betainc(a, b, 0, x, regularized=True), result[1])
            else:
                result = (result[0], mp.betainc(a, b, x, 1, regularized=True))
    _ibeta_cache[key] = result
    return result


def ibeta_points(count, to_type=lambda v: v):
    r = ibeta_rng
    pts = []
    for i in range(count):
        kind = i % 4
        if kind == 0:
            a, b = 2.0 ** r.uniform(-20, 20), 2.0 ** r.uniform(-20, 20)
        elif kind == 1:
            a, b = r.uniform(0.1, 60), r.uniform(0.1, 60)
        elif kind == 2:
            a, b = 2.0 ** r.uniform(-30, 0), 2.0 ** r.uniform(-5, 15)
        else:
            a = 2.0 ** r.uniform(3, 20)
            b = a * 2.0 ** r.uniform(-3, 3)
        if r.random() < 0.5:
            a, b = b, a
        u = r.random()
        mean = a / (a + b)
        if u < 0.3:
            x = r.random()
        elif u < 0.5:
            x = 2.0 ** -r.uniform(1, 60)
        elif u < 0.7:
            x = 1 - 2.0 ** -r.uniform(1, 53)
        else:
            sd = math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
            x = mean + r.uniform(-6, 6) * sd
        a, b, x = to_type(a), to_type(b), to_type(x)
        if a > 0 and b > 0 and 0 < x < 1:
            pts.append((a, b, x))
    return pts


def ibeta_inverse_exact(a, b, t, complement):
    """x tel que Iₓ(a, b) = t (1 − Iₓ(a, b) = t si complement), à 400 bits : variable u = x
    (ou 1 − x, paramètres échangés) dans (0, 1/2], côté minoritaire, bissection en log u puis
    Newton."""
    with mp.workprec(400):
        a, b, t = mpf(a), mpf(b), mpf(t)
        # Sans jamais former 1 − t (nul à 400 bits pour t < 2^-400) : on résout sur le côté
        # demandé, dans la variable où la solution est <= 1/2.
        half_i, half_j = ibeta_reference(a, b, mpf(1) / 2)
        in_x = (t >= half_j) if complement else (t <= half_i)
        aa, bb = (a, b) if in_x else (b, a)
        upper = complement if in_x else not complement         # côté 1 − Iᵤ(aa, bb)
        target = mp.log(t)
        sign = -1 if upper else 1
        def residual(v):
            i_, j_ = ibeta_reference(aa, bb, mp.exp(v))
            f = j_ if upper else i_
            return mp.log(f) - target if f > 0 else -mpf(10) ** 9
        lo, hi = mpf(-3000), mp.log(mpf(1) / 2)
        if sign * residual(hi) < 0:                            # (ne devrait pas arriver)
            return mpf(1) / 2 if in_x else mpf(1) / 2
        for _ in range(80):
            mid = (lo + hi) / 2
            if sign * residual(mid) > 0:
                hi = mid
            else:
                lo = mid
            if hi - lo < mpf(2) ** -10:
                break
        v = (lo + hi) / 2
        lb = mp.loggamma(aa) + mp.loggamma(bb) - mp.loggamma(aa + bb)
        for _ in range(60):
            u = mp.exp(v)
            i_, j_ = ibeta_reference(aa, bb, u)
            f = j_ if upper else i_
            g = mp.log(f) - target
            density_u = mp.exp(aa * mp.log(u) + (bb - 1) * mp.log(1 - u) - lb)   # u f(u)
            step = g / (sign * density_u / f)
            v -= step
            if abs(step) < mpf(2) ** -380:
                break
        u = mp.exp(v)
        return u if in_x else 1 - u


def ibeta_inverse_points(count, to_type=lambda v: v):
    r = ibeta_rng
    pts = []
    for a, b, _ in ibeta_points(count, to_type):
        if max(a, b) > 2.0 ** 12:                              # coût des références (fraction à 400 bits)
            continue
        u = r.random()
        t = r.random() if u < 0.4 else 2.0 ** -r.uniform(1, 600) if u < 0.8 else 0.5 + r.choice([-1, 1]) * 2.0 ** -r.uniform(2, 40)
        t = to_type(t)
        if 0 < t < 1:
            pts.append((a, b, t))
    return pts


def tgamma_points_float():
    xs = [to_f32(rng.uniform(0, 35.04)) for _ in range(800)]
    neg = []
    while len(neg) < 300:
        x = to_f32(rng.uniform(-40, 0))
        if x != round(x):
            neg.append(x)
    return xs + neg


def erfinv_points_double():
    ys = [rng.uniform(-1, 1) for _ in range(1500)]
    ys += [math.copysign(1 - 2.0 ** -k, rng.choice((-1, 1))) for k in range(1, 54)]
    ys += [math.copysign(2.0 ** rng.uniform(-60, -1), rng.choice((-1, 1))) for _ in range(300)]
    ys += [0.5, -0.5, math.nextafter(0.5, 1), math.nextafter(1.0, 0), 2.0 ** -1074, 1e-300]
    return [y for y in ys if abs(y) < 1 and y != 0]


def erfcinv_exact(z):
    z = mpf(z)
    if z < 1:
        return erfc_inverse_exact(z)
    if z == 1:
        return mpf(0)
    return -erfc_inverse_exact(2 - z)


def erfcinv_points_double():
    zs = [rng.uniform(0, 2) for _ in range(800)]
    zs += [2.0 ** -rng.uniform(1, 1074) for _ in range(800)]
    zs += [2 - 2.0 ** -k for k in range(1, 53)] + [2.0 ** -1074, 2.0 ** -1022, 0.5, 1.5, math.nextafter(2.0, 0)]
    return [z for z in zs if 0 < z < 2 and z != 1]


def erfinv_points_float():
    ys = [to_f32(rng.uniform(-1, 1)) for _ in range(800)]
    ys += [to_f32(math.copysign(1 - 2.0 ** -k, rng.choice((-1, 1)))) for k in range(1, 24)]
    return [y for y in ys if abs(y) < 1 and y != 0]


def erfcinv_points_float():
    zs = [to_f32(2.0 ** -rng.uniform(1, 149)) for _ in range(500)]
    zs += [to_f32(rng.uniform(0, 2)) for _ in range(500)]
    return [z for z in zs if 0 < z < 2 and z != 1]


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
    pts = sorted(set(pts))
    for (x, y), exact in zip(pts, parallel_map(f, pts)):
        cr = rounder(exact)
        if math.isinf(cr) or cr == 0:
            continue
        offset = float((exact - mpf(cr)) / mpf(ulp(cr)))
        lines.append(f"        {{{fmt(x)}, {fmt(y)}, {fmt(cr)}, {offset!r}}},\n")
        kept += 1
    lines.append("    };\n\n")
    return "".join(lines), kept


bessel_rng = random.Random(20261118)


def bessel_points(count, convert=float, low=-700.0, high=700.0):
    """(ν, x) : ν entiers, quelconques sous 50, autour de 50 (bascule vers Debye), grands ;
    x log-uniforme ou autour de ν ; résultat dans l'intervalle [e^low, e^high] (estimé)."""
    pts = []
    while len(pts) < count:
        u = bessel_rng.random()
        nu = (float(bessel_rng.randrange(0, 61)) if u < 0.3 else 2.0 ** bessel_rng.uniform(-20, math.log2(50)) if u < 0.6
              else bessel_rng.uniform(45, 55) if u < 0.75 else 2.0 ** bessel_rng.uniform(math.log2(50), math.log2(2e5)))
        x = 2.0 ** bessel_rng.uniform(-30, 11) if bessel_rng.random() < 0.6 or nu < 1 else nu * 2.0 ** bessel_rng.uniform(-3, 2)
        nu, x = convert(nu), convert(x)
        if not (x > 0 and nu >= 0) or math.isinf(x) or math.isinf(nu):
            continue
        if nu >= 1:
            z = x / nu
            q = math.sqrt(1 + z * z)
            estimate = nu * (q + math.log(z / (1 + q))) - 0.5 * math.log(2 * math.pi * nu) - 0.5 * math.log(q)
        else:
            estimate = max(x - 0.5 * math.log(2 * math.pi * x), nu * math.log(x / 2) - math.lgamma(nu + 1))
        if low < estimate < high:
            pts.append((nu, x))
    return pts


def bessel_exact(nu, x):
    return mp.besseli(mpf(nu), mpf(x), maxterms=10 ** 7)


def references3(name, pts, f, rounder, ulp, fmt, ctype):
    lines = [f"    // {name} : a, b, x, arrondi correct de f(a, b, x), position exacte en ulp\n",
             f"    inline constexpr reference3<{ctype}> {name}[] = {{\n"]
    kept = 0
    pts = sorted(set(pts))
    for (a, b, x), exact in zip(pts, parallel_map(f, pts)):
        cr = rounder(exact)
        if math.isinf(cr) or cr == 0:
            continue
        offset = float((exact - mpf(cr)) / mpf(ulp(cr)))
        lines.append(f"        {{{fmt(a)}, {fmt(b)}, {fmt(x)}, {fmt(cr)}, {offset!r}}},\n")
        kept += 1
    lines.append("    };\n\n")
    return "".join(lines), kept


def reference_file():
    out = [HEADER, "#pragma once\n\n",
           "// Valeurs de référence de wsl::special, calculées par mpmath à 256 bits.\n\n",
           "namespace reference_data\n{\n",
           "    template <class T>\n    struct reference\n    {\n        T x, rounded;\n        double offset;\n    };\n\n",
           "    template <class T>\n    struct reference2\n    {\n        T x, y, rounded;\n        double offset;\n    };\n\n",
           "    template <class T>\n    struct reference3\n    {\n        T a, b, x, rounded;\n        double offset;\n    };\n\n"]
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
            ("atan_float", atan_points_float(), mp.atan, round_float, ulp_float, hexf, "float"),
            ("cbrt_double", cbrt_points_double(), real_cbrt, round_double, ulp_double, hexd, "double"),
            ("cbrt_float", cbrt_points_float(), real_cbrt, round_float, ulp_float, hexf, "float"),
            ("erf_double", erf_points_double(), mp.erf, round_double, ulp_double, hexd, "double"),
            ("erfc_double", erfc_points_double(), mp.erfc, round_double, ulp_double, hexd, "double"),
            ("erf_float", erf_points_float(), mp.erf, round_float, ulp_float, hexf, "float"),
            ("erfc_float", erfc_points_float(), mp.erfc, round_float, ulp_float, hexf, "float"),
            ("lgamma_positive_double", lgamma_positive_points_double(), log_abs_gamma, round_double, ulp_double, hexd, "double"),
            ("lgamma_negative_double", lgamma_negative_points_double(), log_abs_gamma, round_double, ulp_double, hexd, "double"),
            ("tgamma_double", tgamma_points_double(), mp.gamma, round_double, ulp_double, hexd, "double"),
            ("lgamma_float", lgamma_points_float(), log_abs_gamma, round_float, ulp_float, hexf, "float"),
            ("tgamma_float", tgamma_points_float(), mp.gamma, round_float, ulp_float, hexf, "float"),
            ("digamma_double", digamma_points_double(), mp.digamma, round_double, ulp_double, hexd, "double"),
            ("digamma_negative_double", digamma_negative_points_double(), mp.digamma, round_double, ulp_double, hexd, "double"),
            ("digamma_float", digamma_points_float(), mp.digamma, round_float, ulp_float, hexf, "float"),
            ("erf_inv_double", erfinv_points_double(), mp.erfinv, round_double, ulp_double, hexd, "double"),
            ("erfc_inv_double", erfcinv_points_double(), erfcinv_exact, round_double, ulp_double, hexd, "double"),
            ("erf_inv_float", erfinv_points_float(), mp.erfinv, round_float, ulp_float, hexf, "float"),
            ("erfc_inv_float", erfcinv_points_float(), erfcinv_exact, round_float, ulp_float, hexf, "float")):
        text, counts[name] = references(name, xs, f, rounder, ulp, fmt, ctype)
        out.append(text)
    atan2_exact = lambda y, x: mp.atan2(mpf(y), mpf(x))
    igamma_double = igamma_points(3000)
    igamma_float = [(a, x) for a, x in igamma_points(800, to_f32_safe) if a < 3e38 and x < 3e38]
    igamma_inverse_double = igamma_inverse_points(500)
    beta_points = [(2.0 ** ibeta_rng.uniform(-40, 40), 2.0 ** ibeta_rng.uniform(-40, 40)) for _ in range(1500)]
    beta_points += [(ibeta_rng.uniform(0.01, 30), ibeta_rng.uniform(0.01, 30)) for _ in range(500)]
    igamma_inverse_float = igamma_inverse_points(200, to_f32_safe)
    for name, pts, f, rounder, ulp, fmt, ctype in (
            ("pow_double", pow_points_double(), pow_exact, round_double, ulp_double, hexd, "double"),
            ("pow_float", pow_points_float(), pow_exact, round_float, ulp_float, hexf, "float"),
            ("atan2_double", atan2_points_double(), atan2_exact, round_double, ulp_double, hexd, "double"),
            ("atan2_float", atan2_points_float(), atan2_exact, round_float, ulp_float, hexf, "float"),
            ("hypot_double", hypot_points_double(), lambda x, y: mp.hypot(mpf(x), mpf(y)), round_double, ulp_double, hexd, "double"),
            ("hypot_float", hypot_points_float(), lambda x, y: mp.hypot(mpf(x), mpf(y)), round_float, ulp_float, hexf, "float"),
            ("gamma_p_double", igamma_double, lambda a, x: igamma_exact(a, x, False), round_double, ulp_double, hexd, "double"),
            ("gamma_q_double", igamma_double, lambda a, x: igamma_exact(a, x, True), round_double, ulp_double, hexd, "double"),
            ("gamma_p_float", igamma_float, lambda a, x: igamma_exact(a, x, False), round_float, ulp_float, hexf, "float"),
            ("gamma_q_float", igamma_float, lambda a, x: igamma_exact(a, x, True), round_float, ulp_float, hexf, "float"),
            ("gamma_p_inv_double", igamma_inverse_double, lambda a, t: igamma_inverse_exact(a, t, False), round_double, ulp_double, hexd, "double"),
            ("gamma_q_inv_double", igamma_inverse_double, lambda a, t: igamma_inverse_exact(a, t, True), round_double, ulp_double, hexd, "double"),
            ("gamma_p_inv_float", igamma_inverse_float, lambda a, t: igamma_inverse_exact(a, t, False), round_float, ulp_float, hexf, "float"),
            ("lbeta_double", beta_points, lambda a, b: mp.loggamma(mpf(a)) + mp.loggamma(mpf(b)) - mp.loggamma(mpf(a) + mpf(b)), round_double, ulp_double, hexd, "double"),
            ("cyl_bessel_i_double", bessel_points(1500), bessel_exact, round_double, ulp_double, hexd, "double"),
            ("cyl_bessel_i_float", bessel_points(400, to_f32_safe, -100.0, 88.0), bessel_exact, round_float, ulp_float, hexf, "float"),
            ("beta_double", beta_points, lambda a, b: mp.exp(mp.loggamma(mpf(a)) + mp.loggamma(mpf(b)) - mp.loggamma(mpf(a) + mpf(b))), round_double, ulp_double, hexd, "double")):
        text, counts[name] = references2(name, pts, f, rounder, ulp, fmt, ctype)
        out.append(text)
    ibeta_double = ibeta_points(2400)
    ibeta_float = [(a, b, x) for a, b, x in ibeta_points(600, to_f32_safe) if max(a, b) < 3e38]
    ibeta_inverse_double = ibeta_inverse_points(400)
    ibeta_inverse_float = [(a, b, t) for a, b, t in ibeta_inverse_points(160, to_f32_safe) if max(a, b) < 3e38]
    for name, pts, f, rounder, ulp, fmt, ctype in (
            ("ibeta_double", ibeta_double, lambda a, b, x: ibeta_reference(a, b, x)[0], round_double, ulp_double, hexd, "double"),
            ("ibetac_double", ibeta_double, lambda a, b, x: ibeta_reference(a, b, x)[1], round_double, ulp_double, hexd, "double"),
            ("ibeta_float", ibeta_float, lambda a, b, x: ibeta_reference(a, b, x)[0], round_float, ulp_float, hexf, "float"),
            ("ibetac_float", ibeta_float, lambda a, b, x: ibeta_reference(a, b, x)[1], round_float, ulp_float, hexf, "float"),
            ("ibeta_inv_double", ibeta_inverse_double, lambda a, b, t: ibeta_inverse_exact(a, b, t, False), round_double, ulp_double, hexd, "double"),
            ("ibetac_inv_double", ibeta_inverse_double, lambda a, b, t: ibeta_inverse_exact(a, b, t, True), round_double, ulp_double, hexd, "double"),
            ("ibeta_inv_float", ibeta_inverse_float, lambda a, b, t: ibeta_inverse_exact(a, b, t, False), round_float, ulp_float, hexf, "float")):
        text, counts[name] = references3(name, pts, f, rounder, ulp, fmt, ctype)
        out.append(text)
    out.append("}\n")
    return "".join(out), counts


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tables", action="store_true", help="réécrire les tables des en-têtes")
    parser.add_argument("--references", metavar="FICHIER", help="écrire les valeurs de référence des tests")
    arguments = parser.parse_args()
    if not arguments.tables and not arguments.references:
        parser.error("préciser --tables et/ou --references FICHIER")
    if arguments.tables:
        (ROOT / "include/wsl/elementary_functions/detail/exponential_tables.hpp").write_text(tables())
        (ROOT / "include/wsl/elementary_functions/detail/trigonometric_tables.hpp").write_text(trig_tables())
        (ROOT / "include/wsl/special_functions/detail/erf_tables.hpp").write_text(erf_tables())
        (ROOT / "include/wsl/special_functions/detail/gamma_tables.hpp").write_text(gamma_tables())
        (ROOT / "include/wsl/special_functions/detail/digamma_tables.hpp").write_text(digamma_tables())
        (ROOT / "include/wsl/special_functions/detail/incomplete_gamma_tables.hpp").write_text(incomplete_gamma_tables())
        (ROOT / "include/wsl/special_functions/detail/bessel_tables.hpp").write_text(bessel_tables())
    if arguments.references:
        text, counts = reference_file()
        output = Path(arguments.references)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text)
        print(f"{sum(counts.values())} points de référence, {len(counts)} fonctions : {output}")


if __name__ == "__main__":
    main()
