"""Von Mises et von Mises–Fisher."""

from common import *                                         # noqa: F401,F403


VON_MISES = [(0.0, 0.01), (0.5, 0.5), (-2.0, 4.0), (1.0, 30.0), (0.0, 1e3), (3.0, 1e5)]
VON_MISES_FISHER = [(0.5, 3), (10.0, 3), (1e3, 3), (3.5, 5), (50.0, 10), (1e4, 4), (200.0, 120)]


def circular_laws(kind):
    """von Mises (a = μ, b = κ) ; von Mises–Fisher (a = κ, b = p, μ = e_p, x = (√(1 − t²), 0, …, 0, t))."""
    if kind != "double":
        return []
    out = []
    count = 120
    bessel = lambda nu, k: mp.besseli(mpf(nu), mpf(k), maxterms=10 ** 7)
    for mu, kappa in VON_MISES:
        m, k = mpf(mu), mpf(kappa)
        xs = [float(mu + rng.uniform(-math.pi, math.pi) * (2.0 ** -rng.uniform(0, 12) if rng.random() < 0.4 else 1)) for _ in range(count)]
        with mp.workprec(400):
            i0, i1 = bessel(0, k), bessel(1, k)
            lnorm = mp.log(2 * mp.pi) + mp.log(i0)
        out.append(("von_mises_pdf", mu, kappa, xs, lambda x, m=m, k=k, lnorm=lnorm: mp.exp(k * mp.cos(mpf(x) - m) - lnorm)))
        out.append(("von_mises_log_pdf", mu, kappa, xs, lambda x, m=m, k=k, lnorm=lnorm: k * mp.cos(mpf(x) - m) - lnorm))
        out.append(("von_mises_resultant", mu, kappa, [0.0], lambda _, i0=i0, i1=i1: i1 / i0))
        out.append(("von_mises_circular_variance", mu, kappa, [0.0], lambda _, i0=i0, i1=i1: 1 - i1 / i0))
        out.append(("von_mises_entropy", mu, kappa, [0.0], lambda _, k=k, i0=i0, i1=i1, lnorm=lnorm: lnorm - k * i1 / i0))
    for kappa, dim in VON_MISES_FISHER:
        k, d = mpf(kappa), dim
        ts = sorted({1 - 2.0 ** -rng.uniform(1, 40) if rng.random() < 0.5 else rng.uniform(-1, 1) for _ in range(count)})
        with mp.workprec(400):
            nu = mpf(d) / 2 - 1
            lc = nu * mp.log(k) - (mpf(d) / 2) * mp.log(2 * mp.pi) - mp.log(bessel(nu, k))
            a = bessel(nu + 1, k) / bessel(nu, k)
        out.append(("von_mises_fisher_pdf", kappa, float(d), ts, lambda t, k=k, lc=lc: mp.exp(lc + k * mpf(t))))
        out.append(("von_mises_fisher_log_pdf", kappa, float(d), ts, lambda t, k=k, lc=lc: lc + k * mpf(t)))
        out.append(("von_mises_fisher_resultant", kappa, float(d), [0.0], lambda _, a=a: a))
        out.append(("von_mises_fisher_entropy", kappa, float(d), [0.0], lambda _, k=k, lc=lc, a=a: -lc - k * a))
    return out


def von_mises_norm(k):
    return 2 * mp.pi * mp.besseli(0, k, maxterms=10 ** 7) * mp.exp(-k) if k > 0 else 2 * mp.pi


def von_mises_tail(k, phi, prec):
    """(G(φ), g(φ)), G(φ) = ∫_φ^π e^{κ(cos t − 1)} dt / (2π e^{−κ} I₀(κ)). En w = 2κ (sin²(t/2) − sin²(φ/2)),
    dt = dw/(κ sin t), sin t = 2√(x(1 − x)), x = sin²(t/2) : sur [0, w(π)/2], coupures géométriques
    depuis l'échelle 2κ sin²(φ/2) de la quasi-singularité en w = −2κ sin²(φ/2) ; sur [w(π)/2, w(π)],
    en u = √(w(π) − w), qui lève la singularité en t = π (intégrande e^{−w} √2/√(κ x))."""
    with mp.workprec(prec):
        k, phi = mpf(k), mpf(phi)
        norm = von_mises_norm(k)
        s0, c0 = mp.sin(phi / 2) ** 2, mp.cos(phi / 2) ** 2
        density = mp.exp(-2 * k * s0) / norm
        if k == 0:
            return (mp.pi - phi) / norm, density
        top = 2 * k * c0                                        # w(π)
        cap = mpf(3 * prec)                                     # e^{−w} négligeable au-delà
        half = top / 2
        lower_end = min(half, cap)
        f = lambda w: mp.exp(-w) / (2 * k * mp.sqrt((s0 + w / (2 * k)) * (c0 - w / (2 * k))))
        pts = [mpf(0)]
        b = max(2 * k * s0, mpf(2) ** -prec)
        while b < lower_end:
            if b > pts[-1]:
                pts.append(b)
            b *= 2
        pts.append(lower_end)
        total = mp.quad(f, pts)
        if half < cap:
            g = lambda u: mp.exp(-(top - u * u)) * mp.sqrt(2) / mp.sqrt(k * (1 - u * u / (2 * k)))
            total += mp.quad(g, [0, mp.sqrt(half)])
        return mp.exp(-2 * k * s0) * total / norm, density


def von_mises_head(k, phi, prec):
    """D(φ) = ∫_0^φ e^{κ(cos t − 1)} dt / (2π e^{−κ} I₀(κ)) (φ petit : D = ½ − G sans compensation)."""
    with mp.workprec(prec):
        k, phi = mpf(k), mpf(phi)
        n = int(min(64, 2 + phi * mp.sqrt(k + 1)))
        pts = [phi * j / n for j in range(n + 1)]
        return mp.quad(lambda t: mp.exp(-2 * k * mp.sin(t / 2) ** 2), pts) / von_mises_norm(k)


def von_mises_solve(k, t, low_prec=53, prec=120):
    """φ tel que G(φ) = t (t <= 1/2) : par G pour t <= 1/4, par D(φ) = 1/2 − t au-delà."""
    t = mpf(t)
    head = t > mpf(1) / 4
    r = mpf(1) / 2 - t if head else t
    if r == 0:
        return mpf(0)
    value = (lambda phi, pr: von_mises_head(k, phi, pr)) if head else (lambda phi, pr: von_mises_tail(k, phi, pr)[0])
    lo, hi = mpf(0), +mp.pi
    for _ in range(25):
        mid = (lo + hi) / 2
        if (value(mid, low_prec) > r) != head:
            lo = mid
        else:
            hi = mid
    phi = (lo + hi) / 2
    with mp.workprec(prec):
        for _ in range(10):
            g = von_mises_tail(k, phi, prec)[1]
            step = (value(phi, prec) - r) / g
            phi = phi - step if head else phi + step
            if abs(step) < mpf(2) ** -100 * (phi + mpf(2) ** -1000):
                break
        return +phi


def von_mises_distribution(kind):
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 50 if kind == "double" else 20
    out = []
    for mu_, kappa_ in VON_MISES:
        m, k = mpf(to_type(mu_)), mpf(to_type(kappa_))
        thetas = [rng.uniform(-math.pi, math.pi) for _ in range(count)]
        thetas += [rng.choice([-1, 1]) * math.pi * (1 - 2.0 ** -rng.uniform(1, 45)) for _ in range(count // 4)]
        thetas += [rng.choice([-1, 1]) * 2.0 ** -rng.uniform(1, 40) for _ in range(count // 4)]
        xs = sorted({to_type(float(m) + th) for th in thetas})
        xs = [x for x in xs if abs(mpf(x) - m) < mp.pi]
        def cdf(x, m=m, k=k, upper=False):
            th = mpf(x) - m
            G = von_mises_tail(k, abs(th), 120)[0]
            return G if (th < 0) != upper else 1 - G
        out.append(("von_mises_cdf", mu_, kappa_, xs, cdf))
        out.append(("von_mises_ccdf", mu_, kappa_, xs, lambda x, cdf=cdf: cdf(x, upper=True)))
        qs = probabilities(count // 8)
        def quantile(p, m=m, k=k, upper=False):
            p = mpf(p)
            low = not upper
            if p > 0.5:
                p, low = 1 - p, not low
            phi = von_mises_solve(k, p)
            return m - phi if low else m + phi
        out.append(("von_mises_quantile", mu_, kappa_, qs, quantile))
        out.append(("von_mises_cquantile", mu_, kappa_, qs, lambda q, quantile=quantile: quantile(q, upper=True)))
    return out


def entries(kind):
    return circular_laws(kind) + von_mises_distribution(kind)


def main():
    run('circular', entries)
