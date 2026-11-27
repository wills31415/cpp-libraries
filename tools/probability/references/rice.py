"""Loi de Rice."""

from common import *                                         # noqa: F401,F403

RICE = [(1.0, 0.0), (1.0, 0.5), (2.0, 3.0), (1.0, 15.0), (0.5, 40.0), (1.0, 1e-3)]


def rice_parts(sigma, nu, x, prec=240):
    """(g, F, S) en x. Les I_k(z) de tous les ordres par la récurrence descendante de Miller
    (I_{k−1} = I_{k+1} + (2k/z) I_k, termes positifs, stable), normalisés par I₀ + 2 Σ I_k = e^z ;
    la série de Marcum de rapport < 1 sommée directement (termes positifs), l'autre par
    complément ; ν = 0 : Rayleigh. Départ N : I_N/e^z < 2^-(prec + 60)."""
    with mp.workprec(prec + 30):
        s_, n_ = mpf(sigma), mpf(nu)
        y, c = mpf(x) / s_, n_ / s_
        e = -(y * y + c * c) / 2
        if c == 0:
            return y * mp.exp(e) / s_, -mp.expm1(e), mp.exp(e)
        z = c * y
        bits = (prec + 60) * 0.7
        N = int(math.sqrt(2 * float(z) * bits) + bits + 2 * math.log1p(float(z)) + 20)
        lower = y < c
        r = y / c if lower else c / y
        two_over_z, inv_r = 2 / z, 1 / r
        above, here = mpf(0), mpf(1)                                    # I_{k+1}, I_k (non normés)
        norm, series, rk = mpf(0), mpf(0), r ** N
        for k in range(N, 0, -1):
            norm += 2 * here
            series += rk * here
            above, here = here, above + k * two_over_z * here           # I_{k−1}
            rk *= inv_r
        norm += here                                                    # I₀
        if not lower:
            series += here                                              # terme k = 0 de S
        scale = mp.exp(z + e) / norm
        v, dens = series * scale, y * here * scale / s_
        return (dens, v, 1 - v) if lower else (dens, 1 - v, v)


def rice_quantile(sigma, nu, t, upper):
    """x tel que F(x) = t (ou S(x) = t) : bissection en log x à 80 bits, puis Newton à 240 bits."""
    t = mpf(t)
    value = lambda x: rice_parts(sigma, nu, x, 80)[2 if upper else 1]
    lo, hi = mpf(2) ** -600, mpf(nu) + 60 * mpf(sigma)
    for _ in range(60):
        mid = mp.sqrt(lo * hi)
        if (value(mid) > t) != upper:
            hi = mid
        else:
            lo = mid
    x = mp.sqrt(lo * hi)
    with mp.workprec(240):
        for _ in range(12):
            f, F, S = rice_parts(sigma, nu, x, 240)
            step = ((S if upper else F) - t) / f
            x = x + step if upper else x - step
            if abs(step) < x * mpf(2) ** -180:
                break
        return +x


def rice_laws():
    out = []
    for sigma, nu in RICE:
        hi = nu + 25 * sigma
        xs = sorted({float(sigma * 2.0 ** rng.uniform(-20, math.log2(hi / sigma))) for _ in range(60)}
                    | {float(nu + sigma * rng.uniform(-3, 3)) for _ in range(30) if nu + sigma * 3 > 0})
        xs = [x for x in xs if x > 0]
        out.append(("rice_pdf", sigma, nu, xs, lambda x, sigma=sigma, nu=nu: rice_parts(sigma, nu, x)[0]))
        out.append(("rice_log_pdf", sigma, nu, xs, lambda x, sigma=sigma, nu=nu: mp.log(rice_parts(sigma, nu, x)[0])))
        out.append(("rice_cdf", sigma, nu, xs, lambda x, sigma=sigma, nu=nu: rice_parts(sigma, nu, x)[1]))
        out.append(("rice_ccdf", sigma, nu, xs, lambda x, sigma=sigma, nu=nu: rice_parts(sigma, nu, x)[2]))
        qs = probabilities(30)
        out.append(("rice_quantile", sigma, nu, qs, lambda p, sigma=sigma, nu=nu: rice_quantile(sigma, nu, p, False)))
        out.append(("rice_cquantile", sigma, nu, qs, lambda q, sigma=sigma, nu=nu: rice_quantile(sigma, nu, q, True)))

        def constants(which, sigma=sigma, nu=nu):
            with mp.workprec(300):
                s_, c = mpf(sigma), mpf(nu) / mpf(sigma)
                w = c * c / 4
                L = mp.exp(-w) * ((1 + 2 * w) * mp.besseli(0, w) + 2 * w * mp.besseli(1, w))
                m = s_ * mp.sqrt(mp.pi / 2) * L
                return m if which == 0 else 2 * s_ * s_ + mpf(nu) ** 2 - m * m
        out.append(("rice_mean", sigma, nu, [0.0], lambda _, f=constants: f(0)))
        out.append(("rice_variance", sigma, nu, [0.0], lambda _, f=constants: f(1)))
    return out


def entries(kind):
    return rice_laws() if kind == "double" else []


def main():
    run("rice", entries)
