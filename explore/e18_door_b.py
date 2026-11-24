"""E18 verification library for walking all the way through Door B (dropping L13 Gaussianity:
exact-covariance elliptical W_{ell,Sigma} and spherical W_{sph} / T75 (4), T75 (6)).

Exploration branch `explore/e18-door-b-spherical` (from `origin/main` `62e0390`, 24 Nov 2026).
Pure Python 3 stdlib, zero external dependencies.
"""
import math
import random


Z95 = 1.959963984540054


# =============================================================================
# E18-1: Exact-covariance symmetric 3-point break of T45 & Chebyshev-min check
# =============================================================================

def simulate_t45_vs_cheb_min_exact_cov(N, *, trials=10000, seed=0):
    """Draw N iid direct scalar errors from the symmetric 3-point law in W_{ell,Sigma}:
        P(e = +-1/sqrt(p)) = p/2,   P(e = 0) = 1 - p,
    with exact E[e] = 0 and exact Var(e) = 1, choosing p = 0.99 / (N * Z95^2).
    Returns coverage rates of Gaussian T45 and exact-variance Chebyshev-min."""
    p = 0.99 / (N * Z95 * Z95)
    val = 1.0 / math.sqrt(p)
    thresh_t45 = Z95 / math.sqrt(N)
    thresh_min = 20.0 ** (1.0 / (2.0 * N))
    rng = random.Random(seed)
    hit_t45 = hit_cheb = 0
    for _ in range(trials):
        es = []
        for _ in range(N):
            u = rng.random()
            if u < 0.5 * p:
                es.append(val)
            elif u < p:
                es.append(-val)
            else:
                es.append(0.0)
        if abs(math.fsum(es) / N) <= thresh_t45:
            hit_t45 += 1
        if min(abs(e) for e in es) <= thresh_min:
            hit_cheb += 1
    return dict(p=p, jump=val, t45_cov=hit_t45 / trials, cheb_min_cov=hit_cheb / trials)


def cheb_min_2scalar_intervals(x1, x2, sigma, *, unimodal=False):
    """Exact 95% kept set on 2 direct scalar readings with exact believed variance sigma^2:
      - arbitrary W_{cov,Sigma}:   t = 20^(1/4) = 2.11474
      - symmetric unimodal (Gauss 1823, P(|e|>=t*sigma) <= 4/(9*t^2)): t = (4/9)^(1/2) * 20^(1/4) = 1.40983.
    Returns merged list of intervals ((lo, hi), ...)."""
    t = (20.0 ** 0.25) * ((4.0 / 9.0) ** 0.5 if unimodal else 1.0)
    r = t * sigma
    lo1, hi1 = x1 - r, x1 + r
    lo2, hi2 = x2 - r, x2 + r
    segs = sorted([(lo1, hi1), (lo2, hi2)])
    if segs[1][0] <= segs[0][1]:
        return ((segs[0][0], max(segs[0][1], segs[1][1])),)
    return tuple(segs)


# =============================================================================
# E18-2 & E18-3: Exact 2D spherical-null tests (unweighted score, Cauchy score,
# spatial-sign score, positive-ray LR, and whole-line LR)
# =============================================================================

def _lr_line_2d(U, V):
    q_uu = math.fsum(u[0] * u[0] + u[1] * u[1] for u in U)
    q_uv = math.fsum(u[0] * v[0] + u[1] * v[1] for u, v in zip(U, V))
    q_vv = math.fsum(v[0] * v[0] + v[1] * v[1] for v in V)
    return max(0.0, 0.5 * (q_vv - q_uu + math.hypot(q_uu - q_vv, 2.0 * q_uv)))


def _lr_pos_2d(U, V, th0):
    q_uu = math.fsum(u[0] * u[0] + u[1] * u[1] for u in U)
    q_uv = math.fsum(u[0] * v[0] + u[1] * v[1] for u, v in zip(U, V))
    q_vv = math.fsum(v[0] * v[0] + v[1] * v[1] for v in V)
    mid = 0.5 * (q_uu + q_vv)
    A, B = 0.5 * (q_vv - q_uu), q_uv
    lo, hi = -2.0 * th0, math.pi - 2.0 * th0
    cands = [mid + A * math.cos(lo) - B * math.sin(lo), mid + A * math.cos(hi) - B * math.sin(hi)]
    R = math.hypot(A, B)
    if R > 0.0:
        phi = math.atan2(B, -A)
        for sh in (-2.0 * math.pi, 0.0, 2.0 * math.pi):
            if lo - 1e-12 <= phi + sh <= hi + 1e-12:
                cands.append(mid - R)
                break
    return max(0.0, q_vv - min(cands))


def sph_2d_family_pvalues(data, m0, *, draws=3000, seed=18):
    """Evaluate 5 exact spherical-null T75 (4) tests simultaneously on 2D isotropic unit-noise
    vector readings `data = [((fx, fy), (ax, ay)), ...]` at candidate mass `m0`:
      - unw_score:    |sum U_i . V_i|                      (Gaussian-LMP)
      - cauchy_score: |sum (1 / (1 + ||V_i||^2)) U_i . V_i| (Cauchy-LMP)
      - sign_score:   |sum U_i . (V_i / ||V_i||)|           (Laplace-LMP, spatial sign)
      - sph_lr_pos:   K(m0) - min_{m in [0, inf]} K(m)      (T75 (4), L1 positive ray)
      - sph_lr_line:  K(m0) - min_{m in R} K(m)             (T75 (4), whole projective line)."""
    th0 = 0.5 * math.pi if math.isinf(m0) else math.atan(m0)
    su, cu = math.sin(th0), math.cos(th0)
    Us = [(f[0] * su + a[0] * cu, f[1] * su + a[1] * cu) for f, a in data]
    Vs = [(f[0] * cu - a[0] * su, f[1] * cu - a[1] * su) for f, a in data]
    v_norms = [math.hypot(v[0], v[1]) for v in Vs]

    def s_unw(U, V):
        return abs(math.fsum(u[0] * v[0] + u[1] * v[1] for u, v in zip(U, V)))

    def s_cauchy(U, V):
        return abs(math.fsum((u[0] * v[0] + u[1] * v[1]) / (1.0 + vn * vn)
                             for u, v, vn in zip(U, V, v_norms)))

    def s_sign(U, V):
        return abs(math.fsum((u[0] * v[0] + u[1] * v[1]) / vn
                             for u, v, vn in zip(U, V, v_norms) if vn > 1e-12))

    fns = (
        ("unw_score", s_unw),
        ("cauchy_score", s_cauchy),
        ("sign_score", s_sign),
        ("sph_lr_pos", lambda U, V: _lr_pos_2d(U, V, th0)),
        ("sph_lr_line", _lr_line_2d),
    )
    obs = [fn(Us, Vs) for _, fn in fns]
    ge = [0] * len(fns)
    rng = random.Random(seed)
    for _ in range(draws):
        Vd = []
        for vn in v_norms:
            if vn > 0.0:
                ang = 2.0 * math.pi * rng.random()
                Vd.append((vn * math.cos(ang), vn * math.sin(ang)))
            else:
                Vd.append((0.0, 0.0))
        for i, (_, fn) in enumerate(fns):
            if fn(Us, Vd) >= obs[i] - 1e-11:
                ge[i] += 1
    return {name: (1 + ge[i]) / (draws + 1) for i, (name, _) in enumerate(fns)}


def argmin_k_2d_unit(data):
    """Exact T81 argmin K on [0, inf] for 2D isotropic unit-noise readings."""
    s_ff = math.fsum(f[0] * f[0] + f[1] * f[1] for f, a in data)
    s_fa = math.fsum(f[0] * a[0] + f[1] * a[1] for f, a in data)
    s_aa = math.fsum(a[0] * a[0] + a[1] * a[1] for f, a in data)
    if s_fa > 0.0:
        diff = s_ff - s_aa
        return (diff + math.hypot(diff, 2.0 * s_fa)) / (2.0 * s_fa)
    return 0.0 if s_ff <= s_aa else math.inf


# =============================================================================
# E18-3(c): Why t_joint under-covers under per-reading heavy tails despite
# Bakirov-Szekely (2005) / Ibragimov-Muller (2010) holding for equal-weight means
# =============================================================================

def compare_bakirov_equal_weight_vs_eiv_tjoint(N=5, *, nu=1.0, signal=0.5, trials=4000, seed=183):
    """Compare 95% Student-t_{N-1} coverage under independent per-reading heavy-tailed scales
    v_i ~ Inv-Chi2(nu) (Cauchy when nu=1):
      - Equal-weight direct scalar mean (Bakirov-Szekely 2005 / Ibragimov-Muller 2010): >= 0.95!
      - EIV t_joint (T82 (3), where reading i's leverage U_i grows with sqrt(v_i)): UNDER-COVERS!"""
    # Exact 95% two-sided t critical values for small df
    t_crit = {1: 12.7062047364, 2: 4.3026527299, 3: 3.1824463053, 4: 2.7764451052, 9: 2.2621571627}[N - 1]
    rng = random.Random(seed)
    hit_direct = hit_eiv = 0
    for _ in range(trials):
        xs = []
        us = []
        vs = []
        for _ in range(N):
            # Per-reading heavy-tailed scale v_i = nu / chi2_nu
            chi = 2.0 * rng.gammavariate(0.5 * nu, 1.0)
            sd = math.sqrt(nu / chi)
            # Direct scalar mean e_i ~ N(0, sd^2)
            xs.append(rng.gauss(0.0, sd))
            # EIV true candidate frame at m0: U_i = signal + N(0, sd^2), V_i = N(0, sd^2)
            us.append(signal + rng.gauss(0.0, sd))
            vs.append(rng.gauss(0.0, sd))
        # 1. Equal-weight Student-t_{N-1} on direct scalar errors xs
        xbar = math.fsum(xs) / N
        s2 = math.fsum((x - xbar) ** 2 for x in xs) / (N - 1)
        t_dir = abs(xbar) / math.sqrt(s2 / N) if s2 > 0.0 else 0.0
        if t_dir <= t_crit:
            hit_direct += 1
        # 2. EIV t_joint (T82 (3)): Z = (sum U_i V_i) / sqrt(sum U_i^2), W = sum V_i^2 - Z^2
        q_uu = math.fsum(u * u for u in us)
        q_uv = math.fsum(u * v for u, v in zip(us, vs))
        q_vv = math.fsum(v * v for v in vs)
        z = q_uv / math.sqrt(q_uu) if q_uu > 0.0 else 0.0
        w = max(0.0, q_vv - z * z)
        t_eiv = abs(z) / math.sqrt(w / (N - 1)) if w > 0.0 else math.inf
        if t_eiv <= t_crit:
            hit_eiv += 1
    return dict(direct_equal_weight_cov=hit_direct / trials, eiv_tjoint_cov=hit_eiv / trials)
