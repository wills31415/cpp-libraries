#!/usr/bin/env python3
"""Engendre les tables de la ziggourat normale de wsl::random (tirage rapide).

    uv run --no-project --with mpmath python tools/random/ziggurat.py

Écrit include/wsl/random/detail/ziggurat_tables.hpp (versionné, ne pas modifier
à la main). Ziggourat de Marsaglia et Tsang (2000) pour f(x) = e^{−x²/2}, N = 256
couches d'aire V : x₀ = V/f(R), x₁ = R, x_{i+1} = f⁻¹(f(x_i) + V/x_i), x_N = 0 ;
R est la racine qui ferme la récurrence exactement. Calcul à 256 bits, tables
arrondies au plus proche.
"""

from pathlib import Path

from mpmath import mp, mpf

mp.prec = 256
ROOT = Path(__file__).resolve().parents[2]
N = 256


def f(x):
    return mp.exp(-x * x / 2)


def f_inverse(y):
    return mp.sqrt(-2 * mp.log(y))


def tail_area(r):
    return mp.sqrt(mp.pi / 2) * mp.erfc(r / mp.sqrt(2))


def layers(r):
    """x₀ … x_N pour R = r, et l'écart à la fermeture (f(x_{N−1}) + V/x_{N−1} − 1)."""
    v = r * f(r) + tail_area(r)
    xs = [v / f(r), r]
    for _ in range(2, N):
        y = f(xs[-1]) + v / xs[-1]
        if y >= 1:
            return xs, v, y - 1
        xs.append(f_inverse(y))
    return xs, v, f(xs[-1]) + v / xs[-1] - 1


def main():
    r = mp.findroot(lambda r: layers(r)[2], mpf("3.6541528853610088"))
    xs, v, gap = layers(r)
    assert len(xs) == N and abs(gap) < mpf(2) ** -200, (len(xs), gap)
    xs.append(mpf(0))
    fs = [f(x) for x in xs]
    fs[0] = f(r)                                            # couche de base : f au bord du rectangle
    lines = ["// ENGENDRÉ par tools/random/ziggurat.py — ne pas modifier à la main.\n",
             "#pragma once\n\n",
             f"// Ziggourat normale, N = {N} couches, R = {mp.nstr(r, 20)}, V = {mp.nstr(v, 20)}.\n",
             "// x[i] : bord droit de la couche i (x[0] = V/f(R), x[1] = R, x[N] = 0) ;\n",
             "// f[i] = e^{−x[i]²/2} (f[0] = f(R)).\n\n",
             "namespace wsl::detail\n{\n",
             f"    inline constexpr double ziggurat_r = {float(r).hex()};\n",
             f"    inline constexpr double ziggurat_inv_r = {float(1 / r).hex()};\n\n",
             f"    inline constexpr double ziggurat_x[{N + 1}] = {{\n"]
    lines += [f"        {float(x).hex()},\n" for x in xs]
    lines.append("    };\n\n")
    lines.append(f"    inline constexpr double ziggurat_f[{N + 1}] = {{\n")
    lines += [f"        {float(y).hex()},\n" for y in fs]
    lines.append("    };\n}\n")
    path = ROOT / "include/wsl/random/detail/ziggurat_tables.hpp"
    path.write_text("".join(lines))
    print(f"R = {mp.nstr(r, 17)}, V = {mp.nstr(v, 17)} : {path}")


if __name__ == "__main__":
    main()
