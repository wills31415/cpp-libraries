"""Normale asymétrique."""

from common import *                                         # noqa: F401,F403


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


def entries(kind):
    return skew_normal_laws(kind)


def main():
    run('skew_normal', entries)
