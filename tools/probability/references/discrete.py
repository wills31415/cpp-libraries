"""Binomiale négative, hypergéométrique, bêta-binomiale, logarithmique, gaussienne discrète."""

from common import *                                         # noqa: F401,F403


NEGATIVE_BINOMIAL = [(3.5, 0.4), (0.3, 0.05), (1e3, 0.999), (50.0, 0.5), (1e-3, 0.9), (2.0, 1e-4)]
HYPERGEOMETRIC = [(50, 20, 10), (1000, 300, 400), (10 ** 6, 10 ** 3, 5 * 10 ** 4), (10 ** 9, 5 * 10 ** 8, 10 ** 6), (30, 29, 15)]


def negative_binomial_log_pmf(r, pp, k):
    with mp.workprec(300):
        r, pp, k = mpf(r), mpf(pp), mpf(k)
        return mp.loggamma(k + r) - mp.loggamma(k + 1) - mp.loggamma(r) + r * mp.log(pp) + k * mp.log1p(-pp)


def hypergeometric_log_pmf(N, K, n, k):
    lc = lambda a, b: mp.loggamma(a + 1) - mp.loggamma(b + 1) - mp.loggamma(a - b + 1)
    with mp.workprec(200):
        return lc(mpf(K), mpf(k)) + lc(mpf(N - K), mpf(n - k)) - lc(mpf(N), mpf(n))


def hypergeometric_tail(N, K, n, k, upper):
    """Σ_{j<=k} P(j) (ou Σ_{j>k}) par le rapport de masses consécutives, à 320 bits."""
    lo, hi = max(0, n - (N - K)), min(n, K)
    with mp.workprec(200):
        j = k + 1 if upper else k
        if j < lo or j > hi:
            return mpf(0)
        base = mp.exp(hypergeometric_log_pmf(N, K, n, j))
        total, term = mpf(1), mpf(1)
        while (j < hi) if upper else (j > lo):
            if upper:
                term *= mpf(K - j) * (n - j) / (mpf(j + 1) * (N - K - n + j + 1))
                j += 1
            else:
                term *= mpf(j) * (N - K - n + j) / (mpf(K - j + 1) * (n - j + 1))
                j -= 1
            total += term
            if term < total * mpf(2) ** -130:
                break
        return base * total


HYPERGEOMETRIC_DRAWS = {(N, K): n for N, K, n in HYPERGEOMETRIC}

BETA_BINOMIAL = [(20, 2.0, 3.0), (1000, 0.5, 0.5), (10 ** 6, 3.0, 40.0), (50, 0.05, 2.0)]
BETA_BINOMIAL_BETA = {(n, a): b for n, a, b in BETA_BINOMIAL}
LOGARITHMIC = [0.9, 0.01, 0.5, 0.999, 1e-6]
DISCRETE_GAUSSIAN = [(0.3, 2.5), (0.3, 0.4), (-7.25, 30.0), (1e6 + 0.5, 1.0), (0.0, 0.05), (2.0, 1e3)]


def beta_binomial_log_pmf(n, a, b, k):
    with mp.workprec(200):
        a, b = mpf(a), mpf(b)
        return mp.log(mp.binomial(n, k)) + mp.log(mp.beta(k + a, n - k + b)) - mp.log(mp.beta(a, b))


_PREFIX_SUMS = {}


def prefix_sum(key, first, ratio, k):
    """Σ_{j<=k} t_j, t_0 = first(), t_{j+1} = t_j ratio(j), à 200 bits ; sommes cumulées gardées par
    processus et prolongées à la demande."""
    with mp.workprec(200):
        if key not in _PREFIX_SUMS:
            t = first()
            _PREFIX_SUMS[key] = ([t], [t])
        terms, sums = _PREFIX_SUMS[key]
        while len(sums) <= k:
            j = len(terms) - 1
            terms.append(terms[-1] * ratio(j))
            sums.append(sums[-1] + terms[-1])
        return sums[k]


def beta_binomial_cdf(n, a, b, k):
    """Σ_{j<=k} P(j) par le rapport de masses, à 200 bits."""
    if k >= n:
        return mpf(1)
    a_, b_ = mpf(a), mpf(b)
    return prefix_sum(("beta_binomial", n, a, b), lambda: mp.exp(beta_binomial_log_pmf(n, a, b, 0)),
                      lambda j: mpf(n - j) * (j + a_) / ((j + 1) * (n - j - 1 + b_)), k)


def logarithmic_cdf(pp, k):
    """Σ_{j=1}^{k} p^j/j / L, à 200 bits (terme d'indice i : p^{i+1}/((i + 1) L))."""
    if k < 1:
        return mpf(0)
    with mp.workprec(200):
        p_ = mpf(pp)
        L = -mp.log1p(-p_)
    return prefix_sum(("logarithmic", pp), lambda: p_ / L, lambda i: p_ * (i + 1) / (i + 2), k - 1)


_DISCRETE_GAUSSIAN_TABLES = {}


def discrete_gaussian_table(mu, sigma):
    """(c, span, termes f(c − span … c + span), sommes cumulées, Z) à 200 bits, une fois par processus ;
    au-delà de c ± (40σ + 40), les termes sont sous 2^-1100 (les quantiles de probabilités
    sous-normales tombent jusque vers 39σ)."""
    key = (mu, sigma)
    if key not in _DISCRETE_GAUSSIAN_TABLES:
        with mp.workprec(200):
            m, s_ = mpf(mu), mpf(sigma)
            c, span = int(mp.nint(m)), int(40 * sigma) + 40
            terms = [mp.exp(-(j - m) ** 2 / (2 * s_ ** 2)) for j in range(c - span, c + span + 1)]
            prefix, acc = [], mpf(0)
            for t in terms:
                acc += t
                prefix.append(acc)
            _DISCRETE_GAUSSIAN_TABLES[key] = (c, span, terms, prefix, acc)
    return _DISCRETE_GAUSSIAN_TABLES[key]


def discrete_gaussian_parts(mu, sigma, k):
    """(P(k), F(k), 1 − F(k)) par sommes directes à 200 bits."""
    c, span, terms, prefix, z = discrete_gaussian_table(mu, sigma)
    with mp.workprec(200):
        i = k - (c - span)
        if i < 0:
            return mpf(0), mpf(0), mpf(1)
        if i >= len(terms):
            return mpf(0), mpf(1), mpf(0)
        lower = prefix[i]
        upper = mp.fsum(terms[i + 1:]) if lower > z / 2 else z - lower      # queue directe quand elle est petite
        return terms[i] / z, lower / z, upper / z


def hypergeometric_cdf_reference(N, K, k):
    n = HYPERGEOMETRIC_DRAWS[(N, K)]
    if k >= min(n, K):
        return mpf(1)
    return hypergeometric_tail(N, K, n, k, False)


def hypergeometric_start(N, K, u):
    n = HYPERGEOMETRIC_DRAWS[(N, K)]
    mean, sd = n * K / N, math.sqrt(max(n * K * (N - K) * (N - n) / (N * N * (N - 1)), 0.25))
    return min(min(n, K), max(max(0, n - (N - K)), mean + sd * float(-mp.sqrt(2) * g.erfcinv_exact(2 * mpf(u)))))


def discrete_laws_batch():
    out = []
    for n, a, b in BETA_BINOMIAL:
        mean = n * a / (a + b); sd = math.sqrt(n * a * b * (a + b + n) / ((a + b) ** 2 * (a + b + 1)))
        ks = sorted({min(n, max(0, int(mean + sd * rng.uniform(-5, 5)))) for _ in range(40)} | {0, n})
        out.append(("beta_binomial_pmf", float(n), a, [float(k) for k in ks], lambda k, n=n, a=a, b=b: mp.exp(beta_binomial_log_pmf(n, a, b, int(k)))))
        out.append(("beta_binomial_log_pmf", float(n), a, [float(k) for k in ks], lambda k, n=n, a=a, b=b: beta_binomial_log_pmf(n, a, b, int(k))))
        if n <= 1000:
            inner = [float(k) for k in ks if k < n]
            out.append(("beta_binomial_cdf", float(n), a, inner, lambda k, n=n, a=a, b=b: beta_binomial_cdf(n, a, b, int(k))))
            out.append(("beta_binomial_ccdf", float(n), a, inner, lambda k, n=n, a=a, b=b: 1 - beta_binomial_cdf(n, a, b, int(k))))
        out.append(("beta_binomial_variance", float(n), a, [0.0], lambda _, n=n, a=a, b=b: mpf(n) * a * b * (a + b + n) / ((mpf(a) + b) ** 2 * (mpf(a) + b + 1))))
    for pp in LOGARITHMIC:
        top = int(min(3000, 40 / -math.log(pp)))
        ks = sorted({max(1, int(2.0 ** rng.uniform(0, math.log2(top + 1)))) for _ in range(40)} | {1, 2})
        L = -mp.log1p(-mpf(pp))
        out.append(("logarithmic_pmf", pp, 0.0, [float(k) for k in ks], lambda k, pp=pp, L=L: mpf(pp) ** int(k) / (int(k) * L)))
        out.append(("logarithmic_cdf", pp, 0.0, [float(k) for k in ks], lambda k, pp=pp: logarithmic_cdf(pp, int(k))))
        out.append(("logarithmic_ccdf", pp, 0.0, [float(k) for k in ks], lambda k, pp=pp: 1 - logarithmic_cdf(pp, int(k))))
        out.append(("logarithmic_mean", pp, 0.0, [0.0], lambda _, pp=pp, L=L: mpf(pp) / ((1 - mpf(pp)) * L)))
        out.append(("logarithmic_variance", pp, 0.0, [0.0], lambda _, pp=pp, L=L: mpf(pp) * (L - pp) / ((1 - mpf(pp)) ** 2 * L ** 2)))
    for mu, sigma in DISCRETE_GAUSSIAN:
        c = round(mu)
        ks = sorted({int(c + sigma * rng.uniform(-9, 9)) for _ in range(40)} | {c, c - 1, c + 1})
        for name, idx in (("pmf", 0), ("cdf", 1), ("ccdf", 2)):
            out.append((f"discrete_gaussian_{name}", mu, sigma, [float(k) for k in ks], lambda k, mu=mu, sigma=sigma, idx=idx: discrete_gaussian_parts(mu, sigma, int(k))[idx]))
        def moments(mu=mu, sigma=sigma, which=0):
            c, span, terms, prefix, z = discrete_gaussian_table(mu, sigma)
            with mp.workprec(200):
                js = range(c - span, c + span + 1)
                m1 = mp.fsum(j * t for j, t in zip(js, terms)) / z
                return m1 if which == 0 else mp.fsum((j - m1) ** 2 * t for j, t in zip(js, terms)) / z
        out.append(("discrete_gaussian_mean", mu, sigma, [0.0], lambda _, moments=moments: moments(which=0)))
        out.append(("discrete_gaussian_variance", mu, sigma, [0.0], lambda _, moments=moments: moments(which=1)))
    for r, pp in NEGATIVE_BINOMIAL:
        q = 1 - pp
        mean, sd = r * q / pp, math.sqrt(r * q) / pp
        ks = sorted({max(0, int(mean + sd * rng.uniform(-6, 12))) for _ in range(60)} | {0, 1, int(mean)} | {int(mean + sd * 2.0 ** rng.uniform(3, 6)) for _ in range(8)})
        out.append(("negative_binomial_pmf", r, pp, [float(k) for k in ks], lambda k, r=r, pp=pp: mp.exp(negative_binomial_log_pmf(r, pp, k))))
        out.append(("negative_binomial_log_pmf", r, pp, [float(k) for k in ks if k != 0 or r != 1], lambda k, r=r, pp=pp: negative_binomial_log_pmf(r, pp, k)))
        out.append(("negative_binomial_cdf", r, pp, [float(k) for k in ks], lambda k, r=r, pp=pp: g.ibeta_reference(r, float(k) + 1, pp)[0]))
        out.append(("negative_binomial_ccdf", r, pp, [float(k) for k in ks], lambda k, r=r, pp=pp: g.ibeta_reference(r, float(k) + 1, pp)[1]))
    for N, K, n in HYPERGEOMETRIC:
        lo, hi = max(0, n - (N - K)), min(n, K)
        mean = n * K / N
        sd = math.sqrt(max(n * K * (N - K) * (N - n) / (N * N * (N - 1)), 0.25))
        ks = sorted({min(hi, max(lo, int(mean + sd * rng.uniform(-9, 9)))) for _ in range(50)} | {lo, hi})
        out.append(("hypergeometric_pmf", float(N), float(K), [float(k) for k in ks], lambda k, N=N, K=K, n=n: mp.exp(hypergeometric_log_pmf(N, K, n, int(k)))))
        out.append(("hypergeometric_log_pmf", float(N), float(K), [float(k) for k in ks], lambda k, N=N, K=K, n=n: hypergeometric_log_pmf(N, K, n, int(k))))
        inner = [float(k) for k in ks if lo <= k < hi]
        out.append(("hypergeometric_cdf", float(N), float(K), inner, lambda k, N=N, K=K, n=n: hypergeometric_tail(N, K, n, int(k), False)))
        out.append(("hypergeometric_ccdf", float(N), float(K), inner, lambda k, N=N, K=K, n=n: hypergeometric_tail(N, K, n, int(k), True)))
        out.append(("hypergeometric_mean", float(N), float(K), [0.0], lambda _, N=N, K=K, n=n: mpf(n) * K / N))
        out.append(("hypergeometric_variance", float(N), float(K), [0.0], lambda _, N=N, K=K, n=n: mpf(n) * K * (N - K) * (N - n) / (mpf(N) ** 2 * (N - 1))))
    return out


def entries(kind):
    return discrete_laws_batch() if kind == "double" else []


def tables():
    return [("quantile_point2", name, discrete_quantiles(name, cases, cdf, start, *low)) for name, cases, cdf, start, *low in (
        ("negative_binomial_quantile", NEGATIVE_BINOMIAL, lambda r, pp, k: g.ibeta_reference(r, float(k) + 1, pp)[0],
         lambda r, pp, u: max(0.0, r * (1 - pp) / pp + math.sqrt(r * (1 - pp)) / pp * normal_start(u))),
        ("beta_binomial_quantile", [(float(n), a) for n, a, _ in BETA_BINOMIAL if n <= 1000],
         lambda n, a, k: beta_binomial_cdf(int(n), a, BETA_BINOMIAL_BETA[(int(n), a)], int(k)), lambda n, a, u: n * u),
        ("logarithmic_quantile", [(pp, 0.0) for pp in LOGARITHMIC], lambda pp, _, k: mpf(0) if k < 1 else logarithmic_cdf(pp, int(k)), lambda pp, _, u: 1.0),
        ("discrete_gaussian_quantile", [(mu, sigma) for mu, sigma in DISCRETE_GAUSSIAN if sigma <= 30],
         lambda mu, sigma, k: discrete_gaussian_parts(mu, sigma, int(k))[1], lambda mu, sigma, u: mu + sigma * normal_start(u), -10 ** 6),
        ("hypergeometric_quantile", [(float(N), float(K)) for N, K, _ in HYPERGEOMETRIC],
         lambda N, K, k: hypergeometric_cdf_reference(int(N), int(K), int(k)), lambda N, K, u: hypergeometric_start(int(N), int(K), u)))]


def main():
    run('discrete', entries, tables)
