"""Loi de Zipf (finie)."""

from common import *                                         # noqa: F401,F403


ZIPF = [(1000, 1.2), (10 ** 6, 1.0), (2 ** 40, 2.0), (50, 0.0), (10 ** 4, 0.5), (100, 30.0), (2 ** 50, 1.0000001)]

def zipf_partial(n, s, a=1):
    """Σ_{j=a}^{n} j^{−s} (directe sur moins de 10⁵ termes, sinon zêta de Hurwitz, ψ pour s = 1), s double exact."""
    with mp.workprec(200):
        s = mpf(s)
        if n < a:
            return mpf(0)
        if n - a < 10 ** 5:                                                 # somme directe : le ζ de Hurwitz de
            return mp.fsum(mpf(j) ** -s for j in range(a, n + 1))           # mpmath n'est précis qu'en absolu
        if s == 1:
            return mp.digamma(n + 1) - mp.digamma(a)
        if s == 0:
            return mpf(n - a + 1)
        return mp.zeta(s, a) - mp.zeta(s, n + 1)


def zipf_cdf(n, s, k):
    if k < 1:
        return mpf(0)
    if k >= n:
        return mpf(1)
    return zipf_partial(k, s) / zipf_partial(n, s)


def zipf_laws():
    out = []
    for n, sz in ZIPF:
        ks = sorted({min(n, max(1, int(2.0 ** rng.uniform(0, math.log2(n))))) for _ in range(40)} | {1, 2, n, n - 1})
        H = lambda n=n, sz=sz: zipf_partial(n, sz)
        out.append(("zipf_pmf", float(n), sz, [float(k) for k in ks], lambda k, n=n, sz=sz: mpf(int(k)) ** -mpf(sz) / zipf_partial(n, sz)))
        out.append(("zipf_log_pmf", float(n), sz, [float(k) for k in ks if k > 1 or sz != 0 or n > 1], lambda k, n=n, sz=sz: -mpf(sz) * mp.log(int(k)) - mp.log(zipf_partial(n, sz))))
        inner = [float(k) for k in ks if k < n]
        out.append(("zipf_cdf", float(n), sz, inner, lambda k, n=n, sz=sz: zipf_cdf(n, sz, int(k))))
        out.append(("zipf_ccdf", float(n), sz, inner, lambda k, n=n, sz=sz: zipf_partial(n, sz, int(k) + 1) / zipf_partial(n, sz)))
        if n <= 10 ** 7:
            def zipf_moment(n=n, sz=sz, which=0):
                with mp.workprec(200):
                    h, s_ = zipf_partial(n, sz), mpf(sz)
                    a = zipf_partial(n, s_ - 1, 2) - zipf_partial(n, s_, 2)
                    if which == 0:
                        return 1 + a / h
                    b = zipf_partial(n, s_ - 2, 2) - 2 * zipf_partial(n, s_ - 1, 2) + zipf_partial(n, s_, 2)
                    return b / h - (a / h) ** 2
            out.append(("zipf_mean", float(n), sz, [0.0], lambda _, f=zipf_moment: f(which=0)))
            out.append(("zipf_variance", float(n), sz, [0.0], lambda _, f=zipf_moment: f(which=1)))
    return out


def entries(kind):
    return zipf_laws() if kind == "double" else []


def tables():
    return [("quantile_point2", "zipf_quantile", discrete_quantiles("zipf_quantile", [(float(n), sz) for n, sz in ZIPF],
                                                                 lambda n, sz, k: zipf_cdf(int(n), sz, int(k)), lambda n, sz, u: 1.0, 1))]


def main():
    run('zipf', entries, tables)
