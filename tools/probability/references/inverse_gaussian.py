"""Inverse-gaussienne."""

from common import *                                         # noqa: F401,F403


INVERSE_GAUSSIAN = [(1.0, 1.0), (1.0, 0.01), (2.0, 50.0), (1e-3, 1e3), (1e5, 3.0), (1.0, 1e-6)]


def inverse_gaussian_parts(mu, lam, x, prec=320):
    """(densité, répartition, survie) de l'inverse-gaussienne (survie compensée d'au plus 2^35 ici)."""
    with mp.workprec(prec):
        x, mu, lam = mpf(x), mpf(mu), mpf(lam)                          # 2λ/μ en flottant : e^{2λ/μ} faux de 2λ/μ·2^-53
        c = mp.sqrt(lam / x)
        a, b = c * (x / mu - 1), c * (x / mu + 1)
        tail = mp.exp(2 * lam / mu) * mp.erfc(b / mp.sqrt(2)) / 2
        lower = mp.erfc(-a / mp.sqrt(2)) / 2 + tail
        upper = mp.erfc(a / mp.sqrt(2)) / 2 - tail
        density = mp.sqrt(lam / (2 * mp.pi * x ** 3)) * mp.exp(-lam * (x - mu) ** 2 / (2 * mu ** 2 * x))
        return +density, +lower, +upper


def inverse_gaussian_quantile(mu, lam, t, upper):
    """x tel que F(x) = t (ou S(x) = t) : bissection en log x à 160 bits, puis Newton à 320 bits."""
    t = mpf(t)
    value = lambda x: inverse_gaussian_parts(mu, lam, x, 160)[2 if upper else 1]
    lo, hi = mpf(2) ** -1100, mpf(2) ** 1100
    for _ in range(50):
        mid = mp.sqrt(lo * hi)
        if (value(mid) > t) != upper:
            hi = mid
        else:
            lo = mid
    x = mp.sqrt(lo * hi)
    with mp.workprec(320):
        for _ in range(10):
            f, F, S = inverse_gaussian_parts(mu, lam, x, 320)
            step = ((S if upper else F) - t) / f
            x = x + step if upper else x - step
            if abs(step) < x * mpf(2) ** -200:
                break
        return +x


def inverse_gaussian_laws(kind):
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 80 if kind == "double" else 30
    out = []
    for mu_, lam_ in INVERSE_GAUSSIAN:
        mu, lam = mpf(to_type(mu_)), mpf(to_type(lam_))
        xs = [to_type(float(mu * log_uniform(2.0 ** -12, 2.0 ** 30))) for _ in range(count)]
        xs += [to_type(float(mu * (1 + rng.choice([-1, 1]) * log_uniform(2.0 ** -30, 0.5)))) for _ in range(count // 2)]
        xs = sorted({x for x in xs if 0 < x < math.inf})
        out.append(("inverse_gaussian_pdf", mu_, lam_, xs, lambda x, mu=mu, lam=lam: inverse_gaussian_parts(mu, lam, x)[0]))
        out.append(("inverse_gaussian_log_pdf", mu_, lam_, xs, lambda x, mu=mu, lam=lam: mp.log(inverse_gaussian_parts(mu, lam, x)[0])))
        out.append(("inverse_gaussian_cdf", mu_, lam_, xs, lambda x, mu=mu, lam=lam: inverse_gaussian_parts(mu, lam, x)[1]))
        out.append(("inverse_gaussian_ccdf", mu_, lam_, xs, lambda x, mu=mu, lam=lam: inverse_gaussian_parts(mu, lam, x)[2]))
        qs = [q for q in probabilities(count // 5) if q > 1e-300 or kind == "double"]
        def quantile(p, mu=mu, lam=lam):
            p = mpf(p)
            return inverse_gaussian_quantile(mu, lam, p, False) if p <= 0.5 else inverse_gaussian_quantile(mu, lam, 1 - p, True)
        def cquantile(q, mu=mu, lam=lam):
            q = mpf(q)
            return inverse_gaussian_quantile(mu, lam, q, True) if q <= 0.5 else inverse_gaussian_quantile(mu, lam, 1 - q, False)
        out.append(("inverse_gaussian_quantile", mu_, lam_, qs, quantile))
        out.append(("inverse_gaussian_cquantile", mu_, lam_, qs, cquantile))
        t = 3 * mu / (2 * lam)
        for c, v in [("variance", mu ** 3 / lam), ("mode", mu * (mp.sqrt(1 + t * t) - t)), ("skewness", 3 * mp.sqrt(mu / lam)), ("excess_kurtosis", 15 * mu / lam)]:
            out.append((f"inverse_gaussian_{c}", mu_, lam_, [0.0], lambda _, v=v: v))
    return out


def entries(kind):
    return inverse_gaussian_laws(kind)


def main():
    run('inverse_gaussian', entries)
