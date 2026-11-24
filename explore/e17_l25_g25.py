"""E17 verification library for L25 algebraic structure theorems and G25 non-Gaussian full-type
vector consistency of T81's argmin K (pure Python 3 stdlib, zero external dependencies).

Exploration branch `explore/e17-l25-and-g25` (from `origin/main` `62e0390`, 24 Nov 2026).
"""
import math
import random

from met.ours34 import Reading, State, _cholesky


# =============================================================================
# Non-Gaussian exact-covariance error generators (mean 0, exact Cov = I)
# =============================================================================

def sample_unit_exp(rng):
    """Skewed non-Gaussian: Exp(1) - 1 has exact mean 0, variance 1, skewness +2, kurtosis 9."""
    return -math.log(1.0 - rng.random()) - 1.0


def sample_unit_skew_discrete(rng):
    """Asymmetric 2-point discrete law: P(X = 2) = 1/5, P(X = -1/2) = 4/5.
    Exact mean = 0, exact variance = 1, skewness = +1.5, kurtosis = 3.25."""
    return 2.0 if rng.random() < 0.2 else -0.5


def sample_unit_t_vec(rng, dim, nu):
    """Spherical multivariate Student-t_nu standardized to EXACT Cov = I_dim (requires nu > 2).
    Let g ~ N(0, I_dim) and chi2_nu = 2 * Gamma(nu/2, 1). Since E[1/chi2_nu] = 1/(nu - 2),
    scaling by sqrt((nu - 2) / chi2_nu) gives exact covariance I_dim."""
    if nu <= 2.0:
        raise ValueError("finite covariance requires nu > 2")
    chi = 2.0 * rng.gammavariate(0.5 * nu, 1.0)
    sc = math.sqrt((nu - 2.0) / chi)
    return [rng.gauss(0.0, 1.0) * sc for _ in range(dim)]


def sample_correlated_error(rng, cov, kind="exp"):
    """Draw e in R^m with EXACT E[e] = 0 and EXACT Cov(e) = cov from non-Gaussian laws:
      - 'exp':      independent standardized Exp(1)-1 coordinates mapped by Cholesky(cov)
      - 'discrete': independent asymmetric 2-point coordinates mapped by Cholesky(cov)
      - 't5':       spherical standardized t_5 (finite 4th moment, heavy tails)
      - 't2.5':     spherical standardized t_2.5 (finite covariance, INFINITE 4th moment)."""
    m = len(cov)
    L = _cholesky(cov)
    if kind == "exp":
        z = [sample_unit_exp(rng) for _ in range(m)]
    elif kind == "discrete":
        z = [sample_unit_skew_discrete(rng) for _ in range(m)]
    elif kind == "t5":
        z = sample_unit_t_vec(rng, m, 5.0)
    elif kind == "t2.5":
        z = sample_unit_t_vec(rng, m, 2.5)
    elif kind == "gaussian":
        z = [rng.gauss(0.0, 1.0) for _ in range(m)]
    else:
        raise ValueError(f"unknown error kind: {kind}")
    return [math.fsum(L[i][j] * z[j] for j in range(i + 1)) for i in range(m)]


# =============================================================================
# Fast 1-instrument and 2-ratio exact global curve minimizers
# =============================================================================

def argmin_k_1instr_gram(s_ff, s_fa, s_aa, s=1.0):
    """Exact closed-form global minimizer of K(m) on [0, inf] for one isotropic instrument
    with noise ratio `s` given whitened sums (sum (F/s)^2, sum (F/s)*a, sum a^2)."""
    if s_fa > 0.0:
        diff = s_ff - s_aa
        r = (diff + math.hypot(diff, 2.0 * s_fa)) / (2.0 * s_fa)
        return s * r
    return 0.0 if s_ff <= s_aa else math.inf


def exact_2ratio_counterexample_trial(rng, N, r1_sq, m0=1.0):
    """Single draw of Theorem E17-G2's global-ray counterexample:
      - 1 reading at s1 = 1 with true signal energy r1_sq at m0 = 1
      - N pure-noise Gaussian readings at s_i = N (so W_noise = N / (1 + N^2)^2 * N^2 ~ 1/N -> 0).
    Returns dict with compact_m (on [0.2, 5.0]), k_compact, k_remote (min over {N/2, N, 2N, inf}),
    and escaped (bool: whether global minimum beats compact minimum at m >= N/2)."""
    r1 = math.sqrt(r1_sq)
    u0_f = m0 / math.hypot(m0, 1.0)
    u0_a = 1.0 / math.hypot(m0, 1.0)
    f1 = r1 * u0_f + rng.gauss(0.0, 1.0)
    a1 = r1 * u0_a + rng.gauss(0.0, 1.0)

    # Exact Wishart sufficient statistics for N iid standard Gaussian 2D pairs:
    # sum ef^2 ~ chi2_N, sum ea^2 ~ chi2_N; exact generation via Bartlett decomposition in O(1)!
    # Bartlett for W_2(N, I_2): L11 = sqrt(chi2_N), L21 = N(0, 1), L22 = sqrt(chi2_{N-1}),
    # n_ff = L11^2, n_fa = L11 * L21, n_aa = L21^2 + L22^2. Exact for EVERY N >= 2 in O(1) time!
    l11_sq = 2.0 * rng.gammavariate(0.5 * N, 1.0)
    l22_sq = 2.0 * rng.gammavariate(0.5 * (N - 1), 1.0)
    l21 = rng.gauss(0.0, 1.0)
    n_ff = l11_sq
    n_fa = math.sqrt(l11_sq) * l21
    n_aa = l21 * l21 + l22_sq

    def K_tot(m):
        if math.isinf(m):
            return a1 * a1 + n_aa
        if m == 0.0:
            return f1 * f1 + n_ff
        k1 = (f1 - m * a1) ** 2 / (1.0 + m * m)
        kn = (N * N * n_ff - 2.0 * N * m * n_fa + m * m * n_aa) / (N * N + m * m)
        return k1 + kn

    cmp_grid = [0.2 * (25.0 ** (j / 40.0)) for j in range(41)]
    k_cmp, m_cmp = min((K_tot(m), m) for m in cmp_grid)
    k_rem, m_rem = min((K_tot(m), m) for m in (0.5 * N, float(N), 2.0 * N, math.inf))
    return dict(m_compact=m_cmp, k_compact=k_cmp, m_remote=m_rem, k_remote=k_rem, escaped=(k_rem < k_cmp))


# =============================================================================
# Spherical null (T75 (4)) exact 1D and Monte Carlo 2D continuity checks (E17-L2)
# =============================================================================

def _lr_pos_1instr(U, V, theta0):
    q_uu = math.fsum(u * u for u in U)
    q_uv = math.fsum(u * v for u, v in zip(U, V))
    q_vv = math.fsum(v * v for v in V)
    mid = 0.5 * (q_uu + q_vv)
    A = 0.5 * (q_vv - q_uu)
    B = q_uv
    lo_2d, hi_2d = -2.0 * theta0, math.pi - 2.0 * theta0
    cands = [
        mid + A * math.cos(lo_2d) - B * math.sin(lo_2d),
        mid + A * math.cos(hi_2d) - B * math.sin(hi_2d),
    ]
    R = math.hypot(A, B)
    if R > 0.0:
        phi = math.atan2(B, -A)
        for sh in (-2.0 * math.pi, 0.0, 2.0 * math.pi):
            if lo_2d - 1e-12 <= phi + sh <= hi_2d + 1e-12:
                cands.append(mid - R)
                break
    return max(0.0, q_vv - min(cands))


def spherical_p_1d_exact(pairs, m0=1.0):
    """Exact 2^n sign-flip p-value of T75 (4) spherical null on 1D isotropic unit-noise readings."""
    th0 = math.atan(m0)
    su, cu = math.sin(th0), math.cos(th0)
    us = [f * su + a * cu for f, a in pairs]
    vs = [f * cu - a * su for f, a in pairs]
    obs = _lr_pos_1instr(us, vs, th0)
    n = len(pairs)
    ge = 0
    for mask in range(1 << n):
        v_flip = [vs[i] if ((mask >> i) & 1) else -vs[i] for i in range(n)]
        if _lr_pos_1instr(us, v_flip, th0) >= obs - 1e-14:
            ge += 1
    return ge / (1 << n), obs


def spherical_p_2d_mc(vec_pairs, m0=1.0, *, draws=4000, seed=0):
    """Exact-direction Monte Carlo p-value of T75 (4) spherical null on 2D isotropic unit-noise
    readings `vec_pairs = [((fx, fy), (ax, ay)), ...]`. Uses fixed RNG draw directions on the
    base readings so tiny appended readings perturb the coupled p-value continuously."""
    th0 = math.atan(m0)
    su, cu = math.sin(th0), math.cos(th0)
    # For 2D isotropic unit noise, U_i and V_i are 2D vectors; K_i(th) = ||V_i cos d - U_i sin d||^2
    # so the sum K(th) depends ONLY on the 3 scalar dot sums q_uu = sum ||U_i||^2, q_uv = sum U_i.V_i, q_vv = sum ||V_i||^2!
    Us = [(f[0] * su + a[0] * cu, f[1] * su + a[1] * cu) for f, a in vec_pairs]
    Vs = [(f[0] * cu - a[0] * su, f[1] * cu - a[1] * su) for f, a in vec_pairs]
    v_norms = [math.hypot(v[0], v[1]) for v in Vs]

    def lr_vec(U_list, V_list):
        q_uu = math.fsum(u[0] * u[0] + u[1] * u[1] for u in U_list)
        q_uv = math.fsum(u[0] * v[0] + u[1] * v[1] for u, v in zip(U_list, V_list))
        q_vv = math.fsum(v[0] * v[0] + v[1] * v[1] for v in V_list)
        mid = 0.5 * (q_uu + q_vv)
        A = 0.5 * (q_vv - q_uu)
        B = q_uv
        lo_2d, hi_2d = -2.0 * th0, math.pi - 2.0 * th0
        cands = [
            mid + A * math.cos(lo_2d) - B * math.sin(lo_2d),
            mid + A * math.cos(hi_2d) - B * math.sin(hi_2d),
        ]
        R = math.hypot(A, B)
        if R > 0.0:
            phi = math.atan2(B, -A)
            for sh in (-2.0 * math.pi, 0.0, 2.0 * math.pi):
                if lo_2d - 1e-12 <= phi + sh <= hi_2d + 1e-12:
                    cands.append(mid - R)
                    break
        return max(0.0, q_vv - min(cands))

    obs = lr_vec(Us, Vs)
    rng = random.Random(seed)
    ge = 0
    for _ in range(draws):
        V_draw = []
        for vn in v_norms:
            if vn == 0.0:
                V_draw.append((0.0, 0.0))
            else:
                ang = 2.0 * math.pi * rng.random()
                V_draw.append((vn * math.cos(ang), vn * math.sin(ang)))
        if lr_vec(Us, V_draw) >= obs - 1e-12:
            ge += 1
    return (1 + ge) / (draws + 1), obs


# =============================================================================
# Full-type correlated non-Gaussian series helper for E17-G1
# =============================================================================

from met.ours34 import _SpanPre, _solve_lower, _solve_spd_small, _orthonormal_basis, _dot


class _CorrelatedInstrument:
    """Precomputed whitened 2k-dim span for a fixed full-type correlated covariance and subspace.
    Evaluates exact sum_r K_r(m) = tr((I - P(m)) S) from the 2k x 2k scatter matrix S = sum zW zW^T
    in O(1) time regardless of the number of readings (identical to met.ours34.State.curve(m))."""

    def __init__(self, cov, subspace=None):
        self.cov = cov
        self.d = len(cov) // 2
        if subspace is None:
            subspace = tuple(tuple(1.0 if i == j else 0.0 for j in range(self.d)) for i in range(self.d))
        self.subspace = subspace
        self.k = len(subspace[0])
        self.L = _cholesky(cov)
        n_cols = [tuple(subspace[row][col] for row in range(self.d)) for col in range(self.k)]
        span_cols = [_solve_lower(self.L, list(col) + [0.0] * self.d) for col in n_cols]
        span_cols += [_solve_lower(self.L, [0.0] * self.d + list(col)) for col in n_cols]
        self.W = _orthonormal_basis(span_cols)  # 2k orthonormal vectors in R^{2d}
        self.bf = tuple(tuple(_dot(w, span_cols[j]) for w in self.W) for j in range(self.k))
        self.ba = tuple(tuple(_dot(w, span_cols[self.k + j]) for w in self.W) for j in range(self.k))
        k = self.k
        self.gff = tuple(tuple(_dot(self.bf[i], self.bf[j]) for j in range(k)) for i in range(k))
        self.gfa = tuple(tuple(_dot(self.bf[i], self.ba[j]) for j in range(k)) for i in range(k))
        self.gaa = tuple(tuple(_dot(self.ba[i], self.ba[j]) for j in range(k)) for i in range(k))
        self.gsym = tuple(tuple(self.gfa[i][j] + self.gfa[j][i] for j in range(k)) for i in range(k))

    def to_zW(self, force, accel):
        z = _solve_lower(self.L, list(force) + list(accel))
        return tuple(_dot(w, z) for w in self.W)

    def curve_from_scatter(self, S, m):
        """Exact sum_r K_r(m) = tr(S) - tr(P(m) S) where P(m) = B G^{-1} B^T in 2k-dim W-coordinates."""
        if math.isinf(m):
            p, q = 1.0, 0.0
        elif m == 0.0:
            p, q = 0.0, 1.0
        else:
            inv = 1.0 / math.hypot(m, 1.0)
            p, q = m * inv, inv
        k = self.k
        m2k = 2 * k
        pp, pq, qq = p * p, p * q, q * q
        G = [[pp * self.gff[i][j] + pq * self.gsym[i][j] + qq * self.gaa[i][j] for j in range(k)]
             for i in range(k)]
        # Columns of B (2k x k)
        B = [[p * self.bf[j][t] + q * self.ba[j][t] for j in range(k)] for t in range(m2k)]
        # Compute M = B^T S B (k x k)
        SB = [[math.fsum(S[t][u] * B[u][j] for u in range(m2k)) for j in range(k)] for t in range(m2k)]
        # tr(P S) = tr(G^{-1} (B^T S B))
        tr_S = math.fsum(S[t][t] for t in range(m2k))
        tr_PS = 0.0
        for col in range(k):
            rhs = [math.fsum(B[t][i] * SB[t][col] for t in range(m2k)) for i in range(k)]
            sol = _solve_spd_small(G, rhs)
            tr_PS += sol[col]
        return max(0.0, tr_S - tr_PS)


_COV3 = (
    (1.0, 0.3, -0.2, 0.25, 0.0, 0.1),
    (0.3, 1.4, 0.1, 0.0, -0.2, 0.05),
    (-0.2, 0.1, 0.8, 0.1, 0.05, 0.15),
    (0.25, 0.0, 0.1, 0.5, 0.1, -0.05),
    (0.0, -0.2, 0.05, 0.1, 0.7, 0.08),
    (0.1, 0.05, 0.15, -0.05, 0.08, 0.6),
)
_SUB3 = ((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))
_COV2 = (
    (1.2, 0.2, 0.15, -0.1),
    (0.2, 0.9, -0.05, 0.2),
    (0.15, -0.05, 0.8, 0.1),
    (-0.1, 0.2, 0.1, 1.1),
)
_COV1 = ((1.5, 0.3), (0.3, 0.7))

_INST3 = _CorrelatedInstrument(_COV3, _SUB3)
_INST2 = _CorrelatedInstrument(_COV2, None)
_INST1 = _CorrelatedInstrument(_COV1, None)


def _add_outer(S, v):
    n = len(v)
    for i in range(n):
        vi = v[i]
        row = S[i]
        for j in range(n):
            row[j] += vi * v[j]


def simulate_full_type_t81_error(n_blocks, amp_scale, *, m0=2.0, seed=0, return_readings=False):
    """Simulate `n_blocks` correlated full-type blocks (3D k=2 Exp(1)-1 + 2D skewed discrete +
    1D standardized t_5; total D = 5 * n_blocks). Uses exact 2k x 2k scatter matrices S_j so
    evaluating K(m) is O(1) per angle (identical to met.ours34.State.curve(m) to 1e-13)."""
    rng = random.Random(seed)
    S3 = [[0.0] * 4 for _ in range(4)]
    S2 = [[0.0] * 4 for _ in range(4)]
    S1 = [[0.0] * 2 for _ in range(2)]
    a3 = (amp_scale, -0.6 * amp_scale, 0.0)
    f3 = tuple(m0 * x for x in a3)
    a2 = (0.7 * amp_scale, 0.5 * amp_scale)
    f2 = tuple(m0 * x for x in a2)
    a1 = (0.8 * amp_scale,)
    f1 = (m0 * a1[0],)
    readings = [] if return_readings else None

    for _ in range(n_blocks):
        e3 = sample_correlated_error(rng, _COV3, kind="exp")
        y3f = tuple(f3[i] + e3[i] for i in range(3))
        y3a = tuple(a3[i] + e3[3 + i] for i in range(3))
        _add_outer(S3, _INST3.to_zW(y3f, y3a))
        e2 = sample_correlated_error(rng, _COV2, kind="discrete")
        y2f = tuple(f2[i] + e2[i] for i in range(2))
        y2a = tuple(a2[i] + e2[2 + i] for i in range(2))
        _add_outer(S2, _INST2.to_zW(y2f, y2a))
        e1 = sample_correlated_error(rng, _COV1, kind="t5")
        y1f = (f1[0] + e1[0],)
        y1a = (a1[0] + e1[1],)
        _add_outer(S1, _INST1.to_zW(y1f, y1a))
        if return_readings:
            readings.append(Reading(y3f, y3a, _COV3, subspace=_SUB3))
            readings.append(Reading(y2f, y2a, _COV2))
            readings.append(Reading(y1f, y1a, _COV1))

    def K_tot(phi):
        m = 0.0 if phi <= 0.0 else (math.inf if phi >= 0.5 * math.pi else math.tan(phi))
        return (_INST3.curve_from_scatter(S3, m)
                + _INST2.curve_from_scatter(S2, m)
                + _INST1.curve_from_scatter(S1, m))

    # Global minimization on compact angle arc [0, pi/2] + golden polish
    phis = [0.5 * math.pi * i / 120.0 for i in range(121)]
    vals = [K_tot(p) for p in phis]
    j = min(range(len(phis)), key=lambda i: vals[i])
    a_p, b_p = phis[max(0, j - 1)], phis[min(len(phis) - 1, j + 1)]
    gr = 0.5 * (math.sqrt(5.0) - 1.0)
    c_p, d_p = b_p - gr * (b_p - a_p), a_p + gr * (b_p - a_p)
    fc, fd = K_tot(c_p), K_tot(d_p)
    for _ in range(36):
        if fc < fd:
            b_p, d_p, fd = d_p, c_p, fc
            c_p = b_p - gr * (b_p - a_p)
            fc = K_tot(c_p)
        else:
            a_p, c_p, fc = c_p, d_p, fd
            d_p = a_p + gr * (b_p - a_p)
            fd = K_tot(d_p)
    phi_hat = c_p if fc < fd else d_p
    if vals[j] < min(fc, fd):
        phi_hat = phis[j]
    m_hat = 0.0 if phi_hat <= 1e-14 else (math.inf if phi_hat >= 0.5 * math.pi - 1e-14 else math.tan(phi_hat))

    z3_0 = _INST3.to_zW(f3, a3)
    z2_0 = _INST2.to_zW(f2, a2)
    z1_0 = _INST1.to_zW(f1, a1)
    delta_one = _dot(z3_0, z3_0) + _dot(z2_0, z2_0) + _dot(z1_0, z1_0)
    out = dict(
        m_hat=m_hat,
        ang_err=abs(phi_hat - math.atan(m0)),
        Delta=n_blocks * delta_one,
        D=5 * n_blocks,
        curve_fn=lambda m: (_INST3.curve_from_scatter(S3, m)
                            + _INST2.curve_from_scatter(S2, m)
                            + _INST1.curve_from_scatter(S1, m)),
    )
    if return_readings:
        out["readings"] = readings
    return out
