"""Loi de Skellam."""

from common import *                                         # noqa: F401,F403


SKELLAM = [(3.0, 1.5), (0.3, 0.02), (15.0, 18.0), (400.0, 380.0), (0.01, 6.0), (2500.0, 0.5)]

_SKELLAM_TABLES = {}


def skellam_table(m1, m2):
    """(lo, termes P(lo … hi), sommes cumulées) à 200 bits ; la table s'étend depuis le mode
    jusqu'à des masses sous 2^-1130 (quantiles de probabilités sous-normales compris)."""
    key = (m1, m2)
    if key not in _SKELLAM_TABLES:
        with mp.workprec(200):
            a, b = mpf(m1), mpf(m2)
            z = 2 * mp.sqrt(a * b)
            mass = lambda j: mp.exp(-(a + b)) * (a / b) ** (mpf(j) / 2) * mp.besseli(abs(j), z)
            c = int(math.floor(m1 - m2))
            tiny = mpf(2) ** -1130
            lo = c
            while mass(lo - 1) > tiny or lo > c - 3:
                lo -= 1
            hi = c
            while mass(hi + 1) > tiny or hi < c + 3:
                hi += 1
            terms = [mass(j) for j in range(lo, hi + 1)]
            prefix, acc = [], mpf(0)
            for t in terms:
                acc += t
                prefix.append(acc)
            _SKELLAM_TABLES[key] = (lo, terms, prefix, acc)
    return _SKELLAM_TABLES[key]


def skellam_parts(m1, m2, k):
    """(P(k), F(k), 1 − F(k)) par sommes directes à 200 bits."""
    lo, terms, prefix, z = skellam_table(m1, m2)
    with mp.workprec(200):
        i = k - lo
        if i < 0:
            return mpf(0), mpf(0), mpf(1)
        if i >= len(terms):
            return mpf(0), mpf(1), mpf(0)
        lower = prefix[i]
        upper = mp.fsum(terms[i + 1:]) if lower > z / 2 else z - lower      # queue directe quand elle est petite
        return terms[i], lower, upper


def skellam_laws():
    out = []
    for m1, m2 in SKELLAM:
        sd, m = math.sqrt(m1 + m2), m1 - m2
        ks = sorted({int(round(m + sd * rng.uniform(-9, 9))) for _ in range(40)} | {int(m), int(m) - 1, int(m) + 1, 0}
                    | {int(round(m + sd * 2.0 ** rng.uniform(3, 4.5) * (1 if i % 2 else -1))) for i in range(6)})
        for name, idx in (("pmf", 0), ("cdf", 1), ("ccdf", 2)):
            out.append((f"skellam_{name}", m1, m2, [float(k) for k in ks], lambda k, m1=m1, m2=m2, idx=idx: skellam_parts(m1, m2, int(k))[idx]))
        out.append(("skellam_log_pmf", m1, m2, [float(k) for k in ks], lambda k, m1=m1, m2=m2: mp.log(skellam_parts(m1, m2, int(k))[0])))
        out.append(("skellam_skewness", m1, m2, [0.0], lambda _, m1=m1, m2=m2: (mpf(m1) - m2) / (mpf(m1) + m2) ** mpf(1.5)))
        out.append(("skellam_excess_kurtosis", m1, m2, [0.0], lambda _, m1=m1, m2=m2: 1 / (mpf(m1) + m2)))
    return out


def entries(kind):
    return skellam_laws() if kind == "double" else []


def tables():
    return [("quantile_point2", "skellam_quantile", discrete_quantiles("skellam_quantile", SKELLAM, lambda m1, m2, k: skellam_parts(m1, m2, int(k))[1],
                                                                       lambda m1, m2, u: m1 - m2 + math.sqrt(m1 + m2) * normal_start(u), -10 ** 6))]


def main():
    run('skellam', entries, tables)
