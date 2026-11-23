"""Lois stables (intégrales de Nolan ; double seulement)."""

from common import *                                         # noqa: F401,F403


STABLE = [(0.5, 1.0), (0.8, -0.3), (1.5, 0.5), (1.0, 0.7), (1.9, 0.0), (1.2, -1.0), (0.3, 0.5), (1.0, 0.0), (2.0, 0.3)]


def stable_parts(alpha, beta, y, prec=120, which=None):
    """(f, F, S) de S1(α, β, 1, 0) en y par les intégrales de Nolan (variable φ = θ + θ₀)."""
    with mp.workprec(prec):
        a, b, y = mpf(alpha), mpf(beta), mpf(y)
        if a == 1:
            if b == 0:
                return 1 / (mp.pi * (1 + y * y)), mp.atan(y) / mp.pi + mpf(1) / 2, mp.atan(-y) / mp.pi + mpf(1) / 2
            if b < 0:
                f, F, S = stable_parts(alpha, -beta, -y, prec, None if which is None else (0, 2, 1)[which])
                return f, S, F
            top = +mp.pi
            def logg(phi):
                th = phi - mp.pi / 2
                lin = mp.pi / 2 + b * th
                return mp.log(2 / mp.pi) + mp.log(abs(lin)) - mp.log(abs(mp.cos(th))) + lin * mp.tan(th) / b - mp.pi * y / (2 * b)
            pref_f, c, swap = 1 / (2 * b), mpf(0), None
        else:
            if y < 0:
                f, F, S = stable_parts(alpha, -beta, -y, prec, None if which is None else (0, 2, 1)[which])
                return f, S, F
            t = b * mp.tan(mp.pi * a / 2)
            th0 = mp.atan(t) / a
            top = mp.pi / 2 + th0
            c = (mp.pi / 2 - th0) / mp.pi
            if y == 0:
                if abs(th0) >= mp.pi / 2 * (1 - mpf(2) ** -100):
                    return mpf(0), c, 1 - c
                f0 = mp.gamma(1 + 1 / a) * mp.cos(th0) * mp.cos(a * th0) ** (1 / a) / mp.pi
                return f0, c, 1 - c
            if top <= 0:
                return mpf(0), mpf(1), mpf(0)
            def logg(phi):
                th = phi - th0
                lv = (mp.log(abs(mp.cos(a * th0))) / (a - 1) + a / (a - 1) * (mp.log(abs(mp.cos(th))) - mp.log(abs(mp.sin(a * phi))))
                      + mp.log(abs(mp.cos(a * th0 + (a - 1) * th))) - mp.log(abs(mp.cos(th))))
                return a / (a - 1) * mp.log(y) + lv
            pref_f = a / (mp.pi * abs(a - 1) * y)
        eps = top * mpf(2) ** -100
        l0, l1 = logg(eps), logg(top - eps)
        cut = None
        if (l0 > 0) != (l1 > 0):
            lo, hi = eps, top - eps
            for _ in range(200):
                m = (lo + hi) / 2
                if (logg(m) > 0) == (l1 > 0):
                    hi = m
                else:
                    lo = m
            cut = lo
        def integral(h):
            """Morceaux géométriques (raison 4) de part et d'autre du pic, intégrande normalisé par
            morceau (mp.quad juge en absolu), arrêt quand le morceau devient négligeable."""
            def piece(u, v):
                lo_, hi_ = (u, v) if u < v else (v, u)
                scale = max(h(lo_ + (hi_ - lo_) * j / 4) for j in range(1, 4)) or mpf(1)
                return scale * mp.quad(lambda q: h(q) / scale, [lo_, (lo_ + hi_) / 2, hi_])
            if cut is None:
                return sum(piece(top * j / 8, top * (j + 1) / 8) for j in range(8))
            total = mpf(0)
            dh = top * mpf(2) ** -40
            slope = abs(logg(cut + dh) - logg(cut - dh)) / (2 * dh)
            width = 1 / slope if slope > 0 else top
            for end in (top, mpf(0)):
                d = end - cut
                edges, w = [cut], width / 8
                while w < abs(d):
                    edges.append(cut + (w if d > 0 else -w))
                    w *= 4
                edges.append(end)
                part = mpf(0)
                for u, v in zip(edges[:-1], edges[1:]):
                    q = piece(u, v)
                    part += q
                    if q < (part + total) * mpf(2) ** -110 and h(v) <= h(u):
                        break
                total += part
            return total
        def safe(fn, big=mpf(0)):
            def h(p):
                lg = logg(p)
                return fn(mp.exp(lg)) if lg < 3000 else big
            return h
        need_a = which is None or (which == 1 and (a == 1 or a < 1)) or (which == 2 and a > 1)
        need_b = which is None or (which == 2 and (a == 1 or a < 1)) or (which == 1 and a > 1)
        I_dens = integral(safe(lambda g: g * mp.exp(-g))) if which in (None, 0) else mpf(0)
        A = integral(safe(lambda g: mp.exp(-g))) / mp.pi if need_a else mpf(0)
        B = integral(safe(lambda g: -mp.expm1(-g), mpf(1))) / mp.pi if need_b else mpf(0)
        f = pref_f * I_dens
        if a == 1:
            return f, A, B
        return (f, c + A, B) if a < 1 else (f, c + B, A)


def stable_value(alpha, beta, y, which, prec=120):
    """Densité (0), répartition (1) ou survie (2) de S1(α, β, 1, 0) en y."""
    if alpha == 2:
        with mp.workprec(prec):
            s2 = mp.sqrt(2)
            return (mp.npdf(y, 0, s2), mp.ncdf(y, 0, s2), mp.ncdf(-y, 0, s2))[which]
    return stable_parts(alpha, beta, y, prec, which)[which]


def stable_quantile(alpha, beta, t, upper):
    """y tel que F(y) = t (ou S(y) = t), t <= ½ : bissection grossière (encadrement exponentiel), puis Newton."""
    with mp.workprec(120):
        t = mpf(t)
        V = lambda y, pr=120: stable_value(alpha, beta, y, 2 if upper else 1, pr)
        too_big = lambda y: (V(y, 60) > t) != upper                         # y au-delà de la solution (60 bits suffisent)
        lo, hi = mpf(-1), mpf(1)
        if alpha < 1 and abs(beta) == 1:
            lo, hi = (mpf(0), mpf(1)) if beta == 1 else (mpf(-1), mpf(0))
        while not too_big(hi):
            lo, hi = hi, hi * 4 if hi > 0 else (mpf(1) if hi == 0 else hi / 4)
        while too_big(lo):
            lo, hi = (lo * 4 if lo < 0 else (mpf(-1) if lo == 0 else lo / 4)), lo
        for _ in range(6):
            m = (lo + hi) / 2
            if too_big(m):
                hi = m
            else:
                lo = m
        y = (lo + hi) / 2
        for i in range(12):
            pr = 60 if i < 3 else 120
            v, f = V(y, pr), stable_value(alpha, beta, y, 0, pr)
            step = (v - t) / f
            y = y + step if upper else y - step
            if i >= 3 and abs(step) < mpf(2) ** -95 * (1 + abs(y)):
                break
        return y


def stable_laws():
    out = []
    for a, b in STABLE:
        ys = sorted({rng.uniform(-4, 4) for _ in range(7)} | {rng.choice([-1, 1]) * 2.0 ** rng.uniform(2, 12) for _ in range(5)} | {0.0})
        if a < 1 and abs(b) == 1:
            ys = [y for y in ys if (y > 0) == (b > 0)]
        out.append(("stable_pdf", a, b, ys, lambda y, a=a, b=b: stable_value(a, b, y, 0)))
        out.append(("stable_cdf", a, b, ys, lambda y, a=a, b=b: stable_value(a, b, y, 1)))
        out.append(("stable_ccdf", a, b, ys, lambda y, a=a, b=b: stable_value(a, b, y, 2)))
        # quantiles : pas de référence (Newton dans les queues lourdes, 5 à 20 s par point) ; vérifiés
        # par encadrement F(x⁻) <= p <= F(x⁺) dans test_random_laws, la répartition l'étant ici
    return out


def entries(kind):
    return stable_laws() if kind == "double" else []


def main():
    run('stable', entries)
