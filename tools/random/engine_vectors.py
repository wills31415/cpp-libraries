#!/usr/bin/env python3
"""Vecteurs de test des moteurs de wsl::random, de sources indépendantes de son code.

    uv run --no-project --with numpy python tools/random/engine_vectors.py

- pcg64 (PCG XSL RR 128/64) et philox4x64 (Philox 4×64, 10 tours) : sorties de
  NumPy, dont on fixe l'état brut.
- xoshiro256++ et splitmix64 : absents de NumPy ; implémentation de référence
  ci-dessous, transcrite de l'algorithme publié (Blackman et Vigna), et deux
  contrôles indépendants : la première sortie de l'état {1, 2, 3, 4} se calcule à
  la main (rotl(1 + 4, 23) + 1 = 41943041), et les sauts (jump : 2^128 pas,
  long_jump : 2^192) sont vérifiés par la puissance de la matrice de transition
  dans GF(2), sans les constantes de saut.

Écrit tests/random_engine_vectors.hpp (petit, versionné).
"""

from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MASK = (1 << 64) - 1


def rotl(x, k):
    return ((x << k) | (x >> (64 - k))) & MASK


# --- splitmix64 ------------------------------------------------------------

def splitmix64(state, count):
    out = []
    for _ in range(count):
        state = (state + 0x9E3779B97F4A7C15) & MASK
        z = state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
        out.append(z ^ (z >> 31))
    return out


# --- xoshiro256++ ----------------------------------------------------------

def xoshiro_step(s):
    s0, s1, s2, s3 = s
    result = (rotl((s0 + s3) & MASK, 23) + s0) & MASK
    t = (s1 << 17) & MASK
    s2 ^= s0
    s3 ^= s1
    s1 ^= s2
    s0 ^= s3
    s2 ^= t
    s3 = rotl(s3, 45)
    return result, [s0, s1, s2, s3]


def xoshiro_outputs(s, count):
    out = []
    for _ in range(count):
        r, s = xoshiro_step(s)
        out.append(r)
    return out


JUMP = [0x180EC6D33CFD0ABA, 0xD5A61266F0C9392C, 0xA9582618E03FC9AA, 0x39ABDC4529B1661C]
LONG_JUMP = [0x76E15D3EFEFDCBBF, 0xC5004E441C522FB3, 0x77710069854EE241, 0x39109BB02ACBE635]


def xoshiro_jump(s, table):
    acc = [0, 0, 0, 0]
    for word in table:
        for b in range(64):
            if word >> b & 1:
                acc = [a ^ x for a, x in zip(acc, s)]
            _, s = xoshiro_step(s)
    return acc


def to_vector(s):
    return s[0] | s[1] << 64 | s[2] << 128 | s[3] << 192


def from_vector(v):
    return [(v >> (64 * i)) & MASK for i in range(4)]


def transition_columns():
    """Matrice de la transition (linéaire sur GF(2)) : la colonne i est l'image de e_i."""
    return [to_vector(xoshiro_step(from_vector(1 << i))[1]) for i in range(256)]


def apply(columns, v):
    r = 0
    i = 0
    while v:
        if v & 1:
            r ^= columns[i]
        v >>= 1
        i += 1
    return r


def square(columns):
    return [apply(columns, c) for c in columns]


def checked_jumps(state):
    """Sauts de 2^128 et 2^192 pas par puissances de la matrice, comparés aux tables."""
    m = transition_columns()
    for _ in range(128):
        m = square(m)
    by_matrix_128 = from_vector(apply(m, to_vector(state)))
    for _ in range(64):
        m = square(m)
    by_matrix_192 = from_vector(apply(m, to_vector(state)))
    jumped = xoshiro_jump(list(state), JUMP)
    long_jumped = xoshiro_jump(list(state), LONG_JUMP)
    assert jumped == by_matrix_128, "table JUMP fausse"
    assert long_jumped == by_matrix_192, "table LONG_JUMP fausse"
    return jumped, long_jumped


# --- pcg64 -----------------------------------------------------------------

PCG_MULTIPLIER = 0x2360ED051FC65DA44385DF649FCCF645


def numpy_pcg64(state, inc, count):
    g = np.random.PCG64()
    g.state = {"bit_generator": "PCG64", "state": {"state": state, "inc": inc}, "has_uint32": 0, "uinteger": 0}
    return [int(v) for v in g.random_raw(count)], g


def pcg64_seed(seed, stream):
    """Initialisation de pcg-cpp (setseq) : inc = 2 stream + 1, état 0, un pas, + seed, un pas."""
    m = (1 << 128) - 1
    inc = ((stream << 1) | 1) & m
    state = (0 * PCG_MULTIPLIER + inc) & m
    state = (state + seed) & m
    state = (state * PCG_MULTIPLIER + inc) & m
    return state, inc


# --- philox4x64 ------------------------------------------------------------

def numpy_philox(key, counter, count):
    """Sorties de NumPy pour notre convention : bloc(counter) d'abord, puis counter + 1.
    NumPy incrémente son compteur avant de calculer un bloc : on lui donne counter − 1."""
    g = np.random.Philox(key=np.array(key, dtype=np.uint64))
    previous = (counter[0] | counter[1] << 64 | counter[2] << 128 | counter[3] << 192) - 1
    words = [(previous >> (64 * i)) & MASK for i in range(4)]
    g.state = {"bit_generator": "Philox",
               "state": {"counter": np.array(words, dtype=np.uint64), "key": np.array(key, dtype=np.uint64)},
               "buffer": np.zeros(4, dtype=np.uint64), "buffer_pos": 4, "has_uint32": 0, "uinteger": 0}
    return [int(v) for v in g.random_raw(count)]


# --- en-tête ---------------------------------------------------------------

def cpp_array(name, values, width=64):
    body = ", ".join(f"0x{v:0{width // 4}x}" for v in values)
    return f"    inline constexpr std::uint64_t {name}[] = {{{body}}};\n"


def u128(name, v):
    return f"    inline constexpr std::uint64_t {name}[2] = {{0x{v >> 64:016x}, 0x{v & MASK:016x}}};   // poids fort, poids faible\n"


def main():
    lines = ["// ENGENDRÉ par tools/random/engine_vectors.py — ne pas modifier à la main.\n",
             "// Vecteurs de test des moteurs, de sources indépendantes (NumPy, matrices GF(2)).\n",
             "#pragma once\n\n#include <cstdint>\n\nnamespace engine_vectors\n{\n"]

    lines.append("    // splitmix64, graines 0 et 1234567.\n")
    lines.append(cpp_array("splitmix64_seed_0", splitmix64(0, 8)))
    lines.append(cpp_array("splitmix64_seed_1234567", splitmix64(1234567, 8)))

    s1234 = [1, 2, 3, 4]
    first = xoshiro_outputs(list(s1234), 8)
    assert first[0] == 41943041, "première sortie de xoshiro256++ fausse"
    lines.append("\n    // xoshiro256++ : état {1, 2, 3, 4} ; état initialisé par splitmix64(42) ; sauts.\n")
    lines.append(cpp_array("xoshiro_state_1234", first))
    seeded = splitmix64(42, 4)
    lines.append(cpp_array("xoshiro_seed_42_state", seeded))
    lines.append(cpp_array("xoshiro_seed_42", xoshiro_outputs(list(seeded), 8)))
    jumped, long_jumped = checked_jumps(s1234)
    lines.append(cpp_array("xoshiro_1234_jump_state", jumped))
    lines.append(cpp_array("xoshiro_1234_long_jump_state", long_jumped))

    lines.append("\n    // pcg64 : état et incrément bruts (ceux de NumPy pour PCG64(0)), puis avance de 10^20.\n")
    state, inc = np.random.PCG64(0).state["state"]["state"], np.random.PCG64(0).state["state"]["inc"]
    outputs, g = numpy_pcg64(state, inc, 8)
    lines.append(u128("pcg64_state", state) + u128("pcg64_inc", inc))
    lines.append(cpp_array("pcg64_outputs", outputs))
    g2 = np.random.PCG64()
    g2.state = {"bit_generator": "PCG64", "state": {"state": state, "inc": inc}, "has_uint32": 0, "uinteger": 0}
    g2.advance(10 ** 20)
    lines.append(cpp_array("pcg64_after_advance_1e20", [int(v) for v in g2.random_raw(4)]))
    seed_state, seed_inc = pcg64_seed(42, 54)
    lines.append("    // pcg64 initialisé par (seed = 42, stream = 54), à la manière de pcg-cpp.\n")
    lines.append(cpp_array("pcg64_seed_42_stream_54", numpy_pcg64(seed_state, seed_inc, 8)[0]))

    lines.append("\n    // philox4x64 : clé et compteur, sorties dans l'ordre (bloc du compteur, puis compteur + 1).\n")
    key = [0x0123456789ABCDEF, 0xFEDCBA9876543210]
    counter = [5, 0, 0, 0]
    lines.append(cpp_array("philox_key", key) + cpp_array("philox_counter", counter))
    lines.append(cpp_array("philox_outputs", numpy_philox(key, counter, 12)))
    carry = [MASK, MASK, 0, 0]
    lines.append("    // Compteur dont l'incrément propage une retenue sur deux mots.\n")
    lines.append(cpp_array("philox_carry_counter", carry))
    lines.append(cpp_array("philox_carry_outputs", numpy_philox(key, carry, 8)))
    lines.append("}\n")

    (ROOT / "tests/random_engine_vectors.hpp").write_text("".join(lines))
    print("tests/random_engine_vectors.hpp écrit ; sauts de xoshiro vérifiés par GF(2).")


if __name__ == "__main__":
    main()
