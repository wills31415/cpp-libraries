"""Gumbel, Fréchet, Rayleigh, triangulaire, arcsinus."""

from common import *                                         # noqa: F401,F403


GUMBEL = [(0.0, 1.0), (1.0, 2.0), (-1e5, 3.0), (0.0, 1e-10)]
FRECHET = [(3.0, 2.0, 0.5), (0.5, 1.0, 0.0), (1.5, 3.0, -2.0), (50.0, 1.0, 0.0), (2.5, 1e-5, 1.0)]
RAYLEIGH = [1.0, 2.0, 1e-10, 3e5]
TRIANGULAR = [(0.0, 1.0, 3.0), (-1.0, -1.0, 2.0), (0.0, 5.0, 5.0), (1e3, 1e3 + 0.25, 1e3 + 1.0)]
ARCSINE = [(0.0, 1.0), (-2.0, 3.0), (1e3, 1e3 + 1e-3)]


def complement_laws(kind):
    """Gumbel, Fréchet, Rayleigh, triangulaire, arcsinus (Erlang : la gamma, déjà mesurée)."""
    to_type = (lambda x: x) if kind == "double" else (lambda x: g.to_f32(x))
    count = 100 if kind == "double" else 40
    out = []
    qs = lambda: probabilities(count // 3)

    def standard(name, a, b, xs, cdf, ccdf, pdf, logpdf, quantile, cquantile, constants):
        out.append((f"{name}_pdf", a, b, xs, pdf))
        out.append((f"{name}_log_pdf", a, b, xs, logpdf))
        out.append((f"{name}_cdf", a, b, xs, cdf))
        out.append((f"{name}_ccdf", a, b, xs, ccdf))
        out.append((f"{name}_quantile", a, b, qs(), quantile))
        out.append((f"{name}_cquantile", a, b, qs(), cquantile))
        for c, f in constants:
            out.append((f"{name}_{c}", a, b, [0.0], lambda _, f=f: f))

    for mu, beta in GUMBEL:
        m, b = mpf(to_type(mu)), mpf(to_type(beta))
        xs = [to_type(float(m + b * rng.choice([-1, 1]) * log_uniform(2.0 ** -20, 40.0))) for _ in range(count)]
        z = lambda x, m=m, b=b: (mpf(x) - m) / b
        standard("gumbel", mu, beta, xs,
                 lambda x, z=z: mp.exp(-mp.exp(-z(x))), lambda x, z=z: -mp.expm1(-mp.exp(-z(x))),
                 lambda x, z=z, b=b: mp.exp(-z(x) - mp.exp(-z(x))) / b, lambda x, z=z, b=b: -z(x) - mp.exp(-z(x)) - mp.log(b),
                 lambda p, m=m, b=b: m - b * mp.log(-mp.log(mpf(p))), lambda q, m=m, b=b: m - b * mp.log(-mp.log1p(-mpf(q))),
                 [("mean", m + b * mp.euler), ("variance", b ** 2 * mp.pi ** 2 / 6), ("median", m - b * mp.log(mp.log(2))), ("entropy", mp.log(b) + mp.euler + 1)])
    for al, sc, loc in FRECHET:
        a, sv, lv = mpf(to_type(al)), mpf(to_type(sc)), mpf(to_type(loc))
        xs = [to_type(float(lv + sv * log_uniform(2.0 ** -1, 2.0 ** 20) ** (1 / a) * (log_uniform(0.2, 3.0) if rng.random() < 0.5 else 1))) for _ in range(count)]
        xs = [x for x in xs if x > float(lv) and x < math.inf]
        zz = lambda x, sv=sv, lv=lv: (mpf(x) - lv) / sv
        const = [("median", lv + sv * mp.log(2) ** (-1 / a)), ("mode", lv + sv * (a / (1 + a)) ** (1 / a)), ("entropy", 1 + mp.euler / a + mp.euler + mp.log(sv / a))]
        if al > 1:
            const.append(("mean", lv + sv * mp.gamma(1 - 1 / a)))
        if al > 2:
            const.append(("variance", sv ** 2 * (mp.gamma(1 - 2 / a) - mp.gamma(1 - 1 / a) ** 2)))
        standard("frechet", al, sc, xs,
                 lambda x, zz=zz, a=a: mp.exp(-zz(x) ** -a), lambda x, zz=zz, a=a: -mp.expm1(-zz(x) ** -a),
                 lambda x, zz=zz, a=a, sv=sv: a / sv * zz(x) ** (-1 - a) * mp.exp(-zz(x) ** -a),
                 lambda x, zz=zz, a=a, sv=sv: mp.log(a / sv) - (1 + a) * mp.log(zz(x)) - zz(x) ** -a,
                 lambda p, a=a, sv=sv, lv=lv: lv + sv * (-mp.log(mpf(p))) ** (-1 / a), lambda q, a=a, sv=sv, lv=lv: lv + sv * (-mp.log1p(-mpf(q))) ** (-1 / a), const)
    for sig in RAYLEIGH:
        sv = mpf(to_type(sig))
        if kind == "float" and not 2.0 ** -100 < sig < 2.0 ** 100:
            continue
        xs = [to_type(float(sv * log_uniform(2.0 ** -30, 30.0))) for _ in range(count)]
        h = lambda x, sv=sv: mpf(x) ** 2 / (2 * sv ** 2)
        standard("rayleigh", sig, 0.0, xs, lambda x, h=h: -mp.expm1(-h(x)), lambda x, h=h: mp.exp(-h(x)),
                 lambda x, h=h, sv=sv: mpf(x) / sv ** 2 * mp.exp(-h(x)), lambda x, h=h, sv=sv: mp.log(mpf(x) / sv ** 2) - h(x),
                 lambda p, sv=sv: sv * mp.sqrt(-2 * mp.log1p(-mpf(p))), lambda q, sv=sv: sv * mp.sqrt(-2 * mp.log(mpf(q))),
                 [("mean", sv * mp.sqrt(mp.pi / 2)), ("variance", (4 - mp.pi) / 2 * sv ** 2), ("median", sv * mp.sqrt(2 * mp.log(2))), ("entropy", 1 + mp.log(sv / mp.sqrt(2)) + mp.euler / 2)])
    for lo_, mo_, hi_ in TRIANGULAR:
        a, c, b = mpf(to_type(lo_)), mpf(to_type(mo_)), mpf(to_type(hi_))
        xs = [to_type(float(a + (b - a) * mpf(rng.random()))) for _ in range(count)] + [to_type(float(a + (c - a) * mpf(2.0 ** -rng.uniform(1, 40)))) for _ in range(count // 4)]
        xs = [x for x in xs if float(a) < x < float(b)]
        def cdf_t(x, a=a, b=b, c=c):
            x = mpf(x)
            return (x - a) ** 2 / ((b - a) * (c - a)) if x <= c else 1 - (b - x) ** 2 / ((b - a) * (b - c))
        def pdf_t(x, a=a, b=b, c=c):
            x = mpf(x)
            return 2 * (x - a) / ((b - a) * (c - a)) if x < c else 2 * (b - x) / ((b - a) * (b - c))
        def q_t(p, a=a, b=b, c=c):
            p = mpf(p)
            return a + mp.sqrt(p * (b - a) * (c - a)) if p * (b - a) <= c - a else b - mp.sqrt((1 - p) * (b - a) * (b - c))
        out.append(("triangular_pdf", lo_, hi_, [x for x in xs if float(x) != float(c)], pdf_t))
        out.append(("triangular_cdf", lo_, hi_, xs, cdf_t))
        out.append(("triangular_ccdf", lo_, hi_, xs, lambda x, a=a, b=b, c=c, cdf_t=cdf_t: (b - mpf(x)) ** 2 / ((b - a) * (b - c)) if mpf(x) > c else 1 - (mpf(x) - a) ** 2 / ((b - a) * (c - a))))
        out.append(("triangular_quantile", lo_, hi_, qs(), q_t))
        out.append(("triangular_cquantile", lo_, hi_, qs(), lambda q, q_t=q_t: q_t(1 - mpf(q))))
        out.append(("triangular_variance", lo_, hi_, [0.0], lambda _, a=a, b=b, c=c: ((a - b) ** 2 + (a - c) ** 2 + (b - c) ** 2) / 36))
        out.append(("triangular_mean", lo_, hi_, [0.0], lambda _, a=a, b=b, c=c: (a + b + c) / 3))
    for lo_, hi_ in ARCSINE:
        a, b = mpf(to_type(lo_)), mpf(to_type(hi_))
        xs = [to_type(float(a + (b - a) * mpf(rng.random()))) for _ in range(count)] + [to_type(float(a + (b - a) * mpf(2.0 ** -rng.uniform(1, 40)))) for _ in range(count // 4)] \
             + [to_type(float(b - (b - a) * mpf(2.0 ** -rng.uniform(1, 40)))) for _ in range(count // 4)]
        xs = [x for x in xs if float(a) < x < float(b)]
        u = lambda x, a=a, b=b: (mpf(x) - a) / (b - a)
        def cdf_a(x, u=u):
            v = u(x)
            return 2 / mp.pi * mp.atan(mp.sqrt(v / (1 - v)))
        def ccdf_a(x, u=u):
            v = u(x)
            return 2 / mp.pi * mp.atan(mp.sqrt((1 - v) / v))
        def q_a(p, a=a, b=b):
            p = mpf(p)
            return a + (b - a) * mp.sin(mp.pi * p / 2) ** 2 if p <= mpf(1) / 2 else b - (b - a) * mp.sin(mp.pi * (1 - p) / 2) ** 2
        def cq_a(q, a=a, b=b):
            q = mpf(q)
            return b - (b - a) * mp.sin(mp.pi * q / 2) ** 2 if q <= mpf(1) / 2 else a + (b - a) * mp.sin(mp.pi * (1 - q) / 2) ** 2
        standard("arcsine", lo_, hi_, xs, cdf_a, ccdf_a,
                 lambda x, a=a, b=b: 1 / (mp.pi * mp.sqrt((mpf(x) - a) * (b - mpf(x)))), lambda x, a=a, b=b: -mp.log(mp.pi) - mp.log((mpf(x) - a) * (b - mpf(x))) / 2,
                 q_a, cq_a, [("variance", (b - a) ** 2 / 8), ("entropy", mp.log(mp.pi * (b - a) / 4))])
    return out


def entries(kind):
    return complement_laws(kind)


def main():
    run('complement', entries)
