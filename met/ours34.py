"""ours34: a probe on the T80 / E13 / E14 output type (no starting law pi, no handed-over nu_eff).

NOT the definition. Built on exploration branch `explore/e14-ours34-walls` (22 Nov 2026,
Antigravity with Can) to turn E13 (T73-T80) and E14 (Theorems E14-1 to E14-7) into one concrete
executable equation so it can be thrown at adversarial scenes and broken (`GOAL.md` section 0).
Nothing in `definition/equation.md`, `definition/premises.md`, or `rulings.md` is changed.

## The equation

Each reading i in the full type (data y_i = (F_i; a_i) in R^{2d}, SPD covariance Sigma_i,
orthonormal candidate subspace N_i in R^{d x k_i}) is mapped to its whitened 2k_i-dimensional
span (T36). At every candidate mass m = tan(phi) on the closed noise-angle ray [0, inf]
(phi in [0, pi/2]):

  K_i(phi) = ||(I - P_i(phi)) z_i||^2                   (T36 / T73 unexplained residual)
  S_i(phi) = -(1/2) dK_i/dphi = V_i^T (P_i'(phi) U_i)   (T45 conditional score)
  C_i(phi) = ||P_i'(phi) U_i||^2 = Var(S_i | U_i)       (T45 conditional variance)

Additive carried state across independent blocks of one body (Wall 3, Theorem E14-3):

  State = (K(phi), C(phi), M_nz),   K = sum_i K_i,  C = sum_i C_i,  M_nz = sum_{y_i != 0} k_i.

Exact 0/0 readings have (K_i, C_i, k_{nz,i}) = (0, 0, 0) and are the exact additive identity (Z).

Point (Wall 1, Theorems E14-1 & E14-2):
  m_hat = argmin_{m in [0, inf]} K(m)   (uniquely forced on interior-minimum curves by L5 + L14
                                         + L3 + exact-branch concentration; closed ray [0, inf]
                                         forced by L5/L14 on one-sided boundary curves).
  On flat or tied curves, returns the declared covariance scale with `point_is_convention=True`
  (R-2026-11-22-2).

Compatible set at declared level 1 - alpha (Wall 2, Theorems E14-5, E14-6, E14-7):
  Raw Gaussian score test (T45, Door A):
    Z_T45(phi) = S(phi) / sqrt(C(phi)) ~ N(0, 1)  exact given U.
  Joint-spherical common-scale Student-t test (t_joint, exact under any common scale v > 0):
    v_hat(phi) = (K(phi) - Z_T45(phi)^2) / (M_nz - 1),   M_nz >= 2,
    t_joint(phi) = Z_T45(phi) / sqrt(v_hat(phi)) ~ t_{M_nz - 1}  exact given U.
  Floor set (`T45` union `t_joint`): keep m iff max(p_T45(m), p_tjoint(m)) > alpha, equivalently
    Z_T45(m)^2 <= max(z_{1-alpha/2}^2,  t_{M_nz-1, 1-alpha/2}^2 * v_hat(m)).
  Leftover K_min = min_{m in [0, inf]} K(m) is carried beside the set (T73 (5)).

Two-body composition (`lump` for L5, `common_force` for L14):
  Curve: T74 infimal convolution K_lump(M) = inf_{m1+m2=M} [K1(m1) + K2(m2)] (and reciprocal dual).
  Point: m1* + m2* (and harmonic dual) on [0, inf].
  Set:   `method="joint"` (default): sup_{0 < m1 < M} p_joint(m1, M - m1) combining the two
         bodies' independent floor p-values via chi2_2 (chi2_1 when one body has C_j = 0):
         zero Monte Carlo, ghost-free on nonlinear asymmetric curves, preserves discordant
         widening without dilution (at the cost of 2-dof 2.45 sigma vs 1.96 sigma at strong signal).
         `method="qpos"` also provided (1-dof linearized Gaussian score composition, Theorem E14-4,
         exposing the nonlinear fake-split break on asymmetric masses).

## Four explicit provisional picks (open forks, NOT forced)

  1. Point codomain: closed ray [0, inf] (vs extending R-2026-11-22-2 to all non-interior curves).
  2. Disagreement widening: common-scale floor `T45 U t_joint` (vs Door A `T45` alone, vs per-reading
     scale floor `T75 (6)`), with `T45` fallback at M_nz = 1 (`single_scalar_policy="t45"`).
  3. Score-test family (additive state, zero MC) vs conditional LR (kills hilltop/pinhole ghosts,
     requires raw readings + MC).
  4. Two-body set composition: `Q_joint` (2-dof conservative, discordant-safe, ghost-free) vs
     `Q_pos` (1-dof Gaussian, fake-split break) vs `T74 (6)` (1-dof MC LR on raw readings).

## Known breaks (surfaces for phase-two sorting)

  - Discontinuity at exact y_i = 0 through M_nz (proved unavoidable for any L11+T30 discount
    without opening a continuum, Theorem E14-6).
  - Stationary hilltop & right-angle ghosts where S(m) = -(1/2) K'(m) = 0 (anti-parallel pairs,
    bimodal multi-instrument conflicts, perpendicular readings).
  - Under-coverage under per-reading heavy-tailed noise scale v_i (0.86 at nu=2 vs T75 (6)'s 0.94).
  - 2-dof width tax (+25% at strong signal) in `Q_joint` two-body composition (vs fake-split break
    in linearized 1-dof `Q_pos`).
  - No residual dof at M_nz = 1 (single scalar reading): cannot detect disagreement on one scalar.

Zero external dependencies (pure Python 3 stdlib). Accepts both `met.full.Reading` (numpy arrays)
and pure-Python `met.ours34.Reading` or any duck-typed reading with `.force`, `.acceleration`,
`.cov`, and optional `.subspace` (d x k matrix with orthonormal columns).
"""
from dataclasses import dataclass
import math


# =============================================================================
# Special functions & deterministic quadrature (pure stdlib)
# =============================================================================

def norm_cdf(x):
    if math.isinf(x):
        return 1.0 if x > 0.0 else 0.0
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def norm_pdf(x):
    if math.isinf(x):
        return 0.0
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def inv_norm_cdf(p):
    """Inverse standard normal CDF Phi^{-1}(p) via Acklam's rational approximation + 1 Halley step
    (exact to ~1e-15 across (1e-300, 1 - 1e-16))."""
    if p <= 0.0:
        return -math.inf
    if p >= 1.0:
        return math.inf
    if p == 0.5:
        return 0.0
    # Acklam coefficients
    a = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
    b = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00)
    p_low = 0.02425
    if p < p_low:
        q = math.sqrt(-2.0 * math.log(p))
        x = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / (
            (((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
    elif p <= 1.0 - p_low:
        q = p - 0.5
        r = q * q
        x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / (
            ((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0)
    else:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        x = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / (
            (((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
    # One Halley polish step (skip extreme tail where pdf underflows)
    pdf = norm_pdf(x)
    if pdf > 1e-250:
        err = norm_cdf(x) - p
        u = err / pdf
        x -= u / (1.0 + 0.5 * x * u)
    return x


def _betacf(a, b, x):
    maxit, eps, fpmin = 200, 3e-15, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def reg_inc_beta(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(a * math.log(x) + b * math.log(1.0 - x) - ln_beta)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_sf_two_sided(t_abs, df):
    """Exact two-sided survival probability P(|T_df| >= t_abs)."""
    if t_abs <= 0.0:
        return 1.0
    if math.isinf(t_abs):
        return 0.0
    if df == 1:
        return 1.0 - (2.0 / math.pi) * math.atan(t_abs)
    if df == 2:
        return 1.0 - t_abs / math.sqrt(2.0 + t_abs * t_abs)
    x = df / (df + t_abs * t_abs)
    return reg_inc_beta(0.5 * df, 0.5, x)


def chi2_sf(x, df):
    """Exact upper tail P(chi2_df >= x) via regularized upper incomplete gamma Q(df/2, x/2)."""
    if x <= 0.0:
        return 1.0
    if math.isinf(x):
        return 0.0
    if df <= 0:
        return 0.0
    a = 0.5 * df
    z = 0.5 * x
    if z < a + 1.0:
        # Series expansion for P(a, z)
        ap = a
        term = 1.0 / a
        total = term
        for _ in range(300):
            ap += 1.0
            term *= z / ap
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        p_lower = total * math.exp(-z + a * math.log(z) - math.lgamma(a))
        return max(0.0, min(1.0, 1.0 - p_lower))
    # Lentz continued fraction for Q(a, z)
    fpmin = 1e-300
    b = z + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d
    for i in range(1, 300):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return max(0.0, min(1.0, math.exp(-z + a * math.log(z) - math.lgamma(a)) * h))


def _gauss_legendre_01(n=24):
    nodes, weights = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(25):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            x -= p1 / dp
        w = 2.0 / ((1.0 - x * x) * dp * dp)
        nodes.append(0.5 * (x + 1.0))
        weights.append(0.5 * w)
    return tuple(nodes), tuple(weights)


_GL_X, _GL_W = _gauss_legendre_01(24)
_GR = 0.5 * (math.sqrt(5.0) - 1.0)


# =============================================================================
# Fast linear algebra (exact for d in {1, 2, 3}, k in {1, 2, 3}; zero QR per eval)
# =============================================================================

def _dot(u, v):
    return math.fsum(a * b for a, b in zip(u, v))


def _cholesky(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = math.fsum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                val = A[i][i] - s
                if val <= 0.0:
                    raise ValueError("covariance matrix must be positive definite (L2)")
                L[i][j] = math.sqrt(val)
            else:
                L[i][j] = (A[i][j] - s) / L[j][j]
    return L


def _solve_lower(L, b):
    n = len(L)
    x = [0.0] * n
    for i in range(n):
        x[i] = (b[i] - math.fsum(L[i][j] * x[j] for j in range(i))) / L[i][i]
    return x


def _solve_spd_small(G, b):
    """Solve k x k SPD linear system for k in {1, 2, 3} with zero allocations overhead."""
    k = len(b)
    if k == 1:
        return (b[0] / G[0][0],)
    if k == 2:
        g00, g01, g11 = G[0][0], G[0][1], G[1][1]
        det = g00 * g11 - g01 * g01
        return ((g11 * b[0] - g01 * b[1]) / det, (-g01 * b[0] + g00 * b[1]) / det)
    # k == 3 (or general small k): Cholesky
    L = _cholesky(G)
    y = _solve_lower(L, b)
    x = [0.0] * k
    for i in range(k - 1, -1, -1):
        x[i] = (y[i] - math.fsum(L[j][i] * x[j] for j in range(i + 1, k))) / L[i][i]
    return tuple(x)


def _orthonormal_basis(cols):
    """Modified Gram-Schmidt (2 passes) for constructing the 2k-dim whitened span basis ONCE."""
    m = len(cols[0])
    Q = []
    for col in cols:
        v = list(col)
        for _ in range(2):
            for q in Q:
                r = _dot(q, v)
                v = [v[t] - r * q[t] for t in range(m)]
        nrm = math.sqrt(_dot(v, v))
        if nrm <= 1e-14:
            raise ValueError("degenerate whitened span basis")
        Q.append(tuple(x / nrm for x in v))
    return tuple(Q)


# =============================================================================
# Reading type & adapter (accepts met.full.Reading or met.ours34.Reading)
# =============================================================================

def _to_1d_tuple(vec):
    if hasattr(vec, "tolist"):
        vec = vec.tolist()
    if isinstance(vec, (int, float)):
        return (float(vec),)
    return tuple(float(x) for x in vec)


def _to_2d_tuple(mat):
    if hasattr(mat, "tolist"):
        mat = mat.tolist()
    return tuple(tuple(float(x) for x in row) for row in mat)


@dataclass(frozen=True)
class Reading:
    """Full-type reading matching `met.full.Reading` conventions (`subspace` is d x k with
    orthonormal columns, or None for I_d), without requiring numpy."""
    force: tuple
    acceleration: tuple
    cov: tuple
    subspace: tuple = None
    window: tuple = None

    def __post_init__(self):
        f = _to_1d_tuple(self.force)
        a = _to_1d_tuple(self.acceleration)
        d = len(f)
        if d not in (1, 2, 3) or len(a) != d:
            raise ValueError("force and acceleration must have equal length in {1, 2, 3}")
        if any(not math.isfinite(x) for x in f + a):
            raise ValueError("force and acceleration entries must be finite")
        cov = _to_2d_tuple(self.cov)
        if len(cov) != 2 * d or any(len(row) != 2 * d for row in cov):
            raise ValueError(f"cov must be {2 * d} x {2 * d}")
        for i in range(2 * d):
            for j in range(2 * d):
                if not math.isfinite(cov[i][j]):
                    raise ValueError("cov entries must be finite")
                if abs(cov[i][j] - cov[j][i]) > 1e-12 * max(1.0, abs(cov[i][j]), abs(cov[j][i])):
                    raise ValueError("cov must be symmetric")
        _cholesky(cov)  # validates positive definiteness (L2)
        if self.subspace is None:
            sub = tuple(tuple(1.0 if i == j else 0.0 for j in range(d)) for i in range(d))
        else:
            raw = self.subspace
            if hasattr(raw, "tolist"):
                raw = raw.tolist()
            if raw and isinstance(raw[0], (int, float)):
                sub = tuple((float(x),) for x in raw)
            else:
                sub = tuple(tuple(float(x) for x in row) for row in raw)
            k = len(sub[0]) if sub else 0
            if len(sub) != d or not (1 <= k <= d) or any(len(row) != k for row in sub):
                raise ValueError("subspace must be d x k with 1 <= k <= d")
            # Validate orthonormal columns
            for c1 in range(k):
                for c2 in range(c1, k):
                    ip = math.fsum(sub[r][c1] * sub[r][c2] for r in range(d))
                    target = 1.0 if c1 == c2 else 0.0
                    if abs(ip - target) > 1e-9:
                        raise ValueError("subspace columns must be orthonormal")
        object.__setattr__(self, "force", f)
        object.__setattr__(self, "acceleration", a)
        object.__setattr__(self, "cov", cov)
        object.__setattr__(self, "subspace", sub)

    @property
    def dimension(self):
        return len(self.force)

    @property
    def k(self):
        return len(self.subspace[0])

    @property
    def is_exact_zero(self):
        return all(x == 0.0 for x in self.force) and all(x == 0.0 for x in self.acceleration)

    def swapped(self):
        d = self.dimension
        perm = list(range(d, 2 * d)) + list(range(d))
        new_cov = tuple(tuple(self.cov[i][j] for j in perm) for i in perm)
        return Reading(self.acceleration, self.force, new_cov, self.subspace, self.window)

    def scaled(self, lam):
        d = self.dimension
        s = [lam] * d + [1.0] * d
        new_cov = tuple(tuple(self.cov[i][j] * s[i] * s[j] for j in range(2 * d)) for i in range(2 * d))
        return Reading(tuple(lam * x for x in self.force), self.acceleration, new_cov, self.subspace, self.window)


def scalar(force, acceleration, force_sd, acceleration_sd, corr=0.0, window=None):
    """Convenience constructor for a 1D reading."""
    cov = (
        (force_sd * force_sd, corr * force_sd * acceleration_sd),
        (corr * force_sd * acceleration_sd, acceleration_sd * acceleration_sd),
    )
    return Reading((force,), (acceleration,), cov, window=window)


def isotropic(force, acceleration, force_sd, acceleration_sd, subspace=None, window=None):
    """Convenience constructor for an isotropic d-dimensional reading."""
    f = _to_1d_tuple(force)
    a = _to_1d_tuple(acceleration)
    d = len(f)
    cov = tuple(
        tuple((force_sd * force_sd if i < d else acceleration_sd * acceleration_sd) if i == j else 0.0
              for j in range(2 * d))
        for i in range(2 * d)
    )
    return Reading(f, a, cov, subspace=subspace, window=window)


def adapt_reading(r):
    """Convert a `met.full.Reading`, `met.ours34.Reading`, or duck-typed reading into `Reading`."""
    if isinstance(r, Reading):
        return r
    return Reading(
        force=r.force,
        acceleration=r.acceleration,
        cov=r.cov,
        subspace=getattr(r, "subspace", None),
        window=getattr(r, "window", None),
    )


# =============================================================================
# Precomputed Whitened Span (zero-QR Gram-block evaluator on phi in [0, pi/2])
# =============================================================================

class _SpanPre:
    """Precomputed whitened span blocks for fast cancellation-free evaluation of
    (K_i, S_{phi,i}, C_{phi,i}) at any angle phi = atan(m) in [0, pi/2]."""

    def __init__(self, reading):
        self.reading = reading
        self.d = reading.dimension
        self.k = reading.k
        self.is_exact_zero = reading.is_exact_zero
        if self.is_exact_zero:
            return
        L = _cholesky(reading.cov)
        y = list(reading.force) + list(reading.acceleration)
        z = _solve_lower(L, y)
        # Columns of subspace N_i (d x k)
        n_cols = [tuple(reading.subspace[row][col] for row in range(self.d)) for col in range(self.k)]
        span_cols = []
        for col in n_cols:
            span_cols.append(_solve_lower(L, list(col) + [0.0] * self.d))
        for col in n_cols:
            span_cols.append(_solve_lower(L, [0.0] * self.d + list(col)))
        W = _orthonormal_basis(span_cols)  # 2k orthonormal vectors in R^{2d}
        self.zW = tuple(_dot(w, z) for w in W)
        self.bf = tuple(tuple(_dot(w, span_cols[j]) for w in W) for j in range(self.k))
        self.ba = tuple(tuple(_dot(w, span_cols[self.k + j]) for w in W) for j in range(self.k))
        # Precompute k x k Gram blocks and RHS vectors
        k = self.k
        self.gff = tuple(tuple(_dot(self.bf[i], self.bf[j]) for j in range(k)) for i in range(k))
        self.gfa = tuple(tuple(_dot(self.bf[i], self.ba[j]) for j in range(k)) for i in range(k))
        self.gaa = tuple(tuple(_dot(self.ba[i], self.ba[j]) for j in range(k)) for i in range(k))
        self.gsym = tuple(tuple(self.gfa[i][j] + self.gfa[j][i] for j in range(k)) for i in range(k))
        self.uf = tuple(_dot(self.bf[i], self.zW) for i in range(k))
        self.ua = tuple(_dot(self.ba[i], self.zW) for i in range(k))

    def eval_phi_pq(self, p, q):
        """Evaluate (K_i, S_{phi,i}, C_{phi,i}) at (p, q) = (sin phi, cos phi), m = p / q.
        Exact and cancellation-free for all m in [0, inf], including m = 0 and m = inf."""
        if self.is_exact_zero:
            return 0.0, 0.0, 0.0
        k = self.k
        pp, pq, qq = p * p, p * q, q * q
        if k == 1:
            g = pp * self.gff[0][0] + pq * self.gsym[0][0] + qq * self.gaa[0][0]
            b = p * self.uf[0] + q * self.ua[0]
            alpha = b / g
            # Derivative projection coefficient h = B^T (B' alpha)
            h = (pq * (self.gff[0][0] - self.gaa[0][0]) + (qq - pp) * self.gfa[0][0]) * alpha
            gamma = h / g
            bf0, ba0, zW = self.bf[0], self.ba[0], self.zW
            K = S = C = 0.0
            for t in range(2):
                b_t = p * bf0[t] + q * ba0[t]
                r_t = zW[t] - alpha * b_t
                w_t = alpha * (q * bf0[t] - p * ba0[t]) - gamma * b_t
                K += r_t * r_t
                S += r_t * w_t
                C += w_t * w_t
            return K, S, C

        G = [[pp * self.gff[i][j] + pq * self.gsym[i][j] + qq * self.gaa[i][j] for j in range(k)]
             for i in range(k)]
        b = [p * self.uf[i] + q * self.ua[i] for i in range(k)]
        alpha = _solve_spd_small(G, b)
        # h_i = sum_j (pq*(gff_ij - gaa_ij) + qq*gfa_ji - pp*gfa_ij) * alpha_j
        h = [
            math.fsum((pq * (self.gff[i][j] - self.gaa[i][j])
                       + qq * self.gfa[j][i] - pp * self.gfa[i][j]) * alpha[j]
                      for j in range(k))
            for i in range(k)
        ]
        gamma = _solve_spd_small(G, h)
        K = S = C = 0.0
        m2k = 2 * k
        for t in range(m2k):
            fit_t = 0.0
            tan_t = 0.0
            for j in range(k):
                bf_jt = self.bf[j][t]
                ba_jt = self.ba[j][t]
                b_jt = p * bf_jt + q * ba_jt
                fit_t += alpha[j] * b_jt
                tan_t += alpha[j] * (q * bf_jt - p * ba_jt) - gamma[j] * b_jt
            r_t = self.zW[t] - fit_t
            K += r_t * r_t
            S += r_t * tan_t
            C += tan_t * tan_t
        return K, S, C


def _mass_to_pq(m):
    if m < 0.0 or math.isnan(m):
        raise ValueError("candidate mass must lie in [0, inf]")
    if math.isinf(m):
        return 1.0, 0.0
    if m == 0.0:
        return 0.0, 1.0
    inv = 1.0 / math.hypot(m, 1.0)
    return m * inv, inv


# =============================================================================
# Additive State (K(phi), C(phi), M_nz)
# =============================================================================

class State:
    """Minimal additive state (K, C, M_nz) carried by ours34. Independent blocks of the
    same body combine by exact addition `state1 + state2`."""

    def __init__(self, readings=(), *, scale=None, single_scalar_policy="t45", _spans=None):
        if single_scalar_policy not in ("t45", "vacuous"):
            raise ValueError("single_scalar_policy must be 't45' or 'vacuous'")
        self.readings = tuple(adapt_reading(r) for r in readings)
        self.spans = tuple(_SpanPre(r) for r in self.readings) if _spans is None else tuple(_spans)
        self.single_scalar_policy = single_scalar_policy
        if scale is None:
            src = [r for r in self.readings if not r.is_exact_zero] or list(self.readings)
            logs = []
            for r in src:
                d = r.dimension
                tr_f = math.fsum(r.cov[i][i] for i in range(d))
                tr_a = math.fsum(r.cov[d + i][d + i] for i in range(d))
                logs.append(0.5 * math.log(tr_f / tr_a))
            self.scale = math.exp(math.fsum(logs) / len(logs)) if logs else 1.0
        else:
            self.scale = float(scale)
            if not (self.scale > 0.0 and math.isfinite(self.scale)):
                raise ValueError("scale must be positive and finite")
        self.M_total = sum(s.k for s in self.spans)
        self.M_nz = sum(s.k for s in self.spans if not s.is_exact_zero)

    def __add__(self, other):
        if not isinstance(other, State):
            return NotImplemented
        return State(
            self.readings + other.readings,
            single_scalar_policy=self.single_scalar_policy,
            _spans=self.spans + other.spans,
        )

    def ksc_phi(self, m):
        """Return (K(m), S_phi(m), C_phi(m)) at m in [0, inf] (phi = atan(m))."""
        p, q = _mass_to_pq(m)
        K = S = C = 0.0
        for sp in self.spans:
            ki, si, ci = sp.eval_phi_pq(p, q)
            K += ki
            S += si
            C += ci
        return K, S, C

    def ksc_m(self, m):
        """Return (K(m), S_m(m), C_m(m)) with respect to mass m (for linear L5 sums)."""
        p, q = _mass_to_pq(m)
        K, S_phi, C_phi = self.ksc_phi(m)
        q2 = q * q
        return K, q2 * S_phi, (q2 * q2) * C_phi

    def curve(self, m):
        return self.ksc_phi(m)[0]

    def mass_at_angle(self, u):
        if u <= 0.0:
            return 0.0
        if u >= 0.5 * math.pi:
            return math.inf
        return self.scale * math.tan(u)

    def angle_of(self, m):
        if math.isinf(m):
            return 0.5 * math.pi
        if m == 0.0:
            return 0.0
        return math.atan(m / self.scale)

    def p_value(self, m):
        """Evaluate exact T45, t_joint, and floor p-values at candidate m in [0, inf]."""
        K, S, C = self.ksc_phi(m)
        if C <= 1e-30 * max(1.0, K):
            p_t45 = 1.0
            z_t45 = 0.0
            t_val = 0.0
            v_hat = None
            df = max(0, self.M_nz - 1)
            p_tjoint = 1.0
        else:
            z_t45 = S / math.sqrt(C)
            p_t45 = math.erfc(abs(z_t45) / math.sqrt(2.0))
            df = self.M_nz - 1
            if df <= 0:
                t_val = 0.0
                v_hat = None
                p_tjoint = 1.0
            else:
                perp_sq = max(0.0, K - z_t45 * z_t45)
                v_hat = perp_sq / df
                if perp_sq <= 1e-28 * max(1.0, K):
                    if abs(z_t45) <= 1e-13:
                        t_val = 0.0
                        p_tjoint = 1.0
                    else:
                        t_val = math.inf
                        p_tjoint = 0.0
                else:
                    t_val = abs(z_t45) / math.sqrt(v_hat)
                    p_tjoint = t_sf_two_sided(t_val, df)
        if self.M_nz <= 1 and self.single_scalar_policy == "t45":
            p_floor = p_t45
        else:
            p_floor = max(p_t45, p_tjoint)
        return dict(
            p_value=p_floor,
            p_t45=p_t45,
            p_tjoint=p_tjoint,
            z_t45=z_t45,
            t_joint=t_val,
            v_hat=v_hat,
            df=max(0, self.M_nz - 1),
            K=K,
            S_phi=S,
            C_phi=C,
        )

    def _base_angles(self, grid=361, log_decades=10, pts_per_decade=12):
        """Scale-equivariant (L3) and swap-invariant (T5) angle mesh on [0, pi/2].
        Combines a uniform angle grid with a log-spaced mass grid spanning 10^(-D)..10^(+D)
        around `self.scale`, so narrow features at tiny sigma_a (e.g. scale=1e5, m=3) are
        resolved without ever using data to build the grid."""
        us = [0.5 * math.pi * i / (grid - 1) for i in range(grid)]
        n_half = log_decades * pts_per_decade
        for k in range(-n_half, n_half + 1):
            r = 10.0 ** (k / pts_per_decade)  # m / scale
            us.append(math.atan(r))
        return sorted(set(us))

    def argmin_K(self, grid=361, tol=1e-12):
        """Global minimum of K on the closed ray [0, inf], with S_phi = 0 root polishing."""
        us = self._base_angles(grid=grid)
        vals = [self.curve(self.mass_at_angle(u)) for u in us]
        vmin = min(vals)
        vmax = max(vals)
        if vmax - vmin <= tol * max(1.0, vmax):
            return dict(
                mass=self.scale, status="flat", convention=True,
                K_min=vmin, minima=(), K_0=vals[0], K_inf=vals[-1],
            )
        n_u = len(us)
        cand = []
        for j in range(n_u):
            left = vals[j - 1] if j > 0 else math.inf
            right = vals[j + 1] if j < n_u - 1 else math.inf
            if vals[j] <= left and vals[j] <= right and vals[j] <= vmin + 1e-2 * max(1.0, abs(vmin)) + 1e-6 * (vmax - vmin):
                cand.append(j)
        refined = []
        for j in cand:
            a = us[max(0, j - 1)]
            b = us[min(n_u - 1, j + 1)]
            fu, u_best = _golden_min(lambda u: self.curve(self.mass_at_angle(u)), a, b)
            if vals[j] < fu:
                fu, u_best = vals[j], us[j]
            if 1e-15 < u_best < 0.5 * math.pi - 1e-15:
                span_u = max(1e-12, min(1e-4, 0.5 * (b - a)))
                lo_u, hi_u = max(a, u_best - span_u), min(b, u_best + span_u)
                s_lo = self.ksc_phi(self.mass_at_angle(lo_u))[1]
                s_hi = self.ksc_phi(self.mass_at_angle(hi_u))[1]
                if s_lo > 0.0 and s_hi < 0.0:
                    for _ in range(56):
                        mid_u = 0.5 * (lo_u + hi_u)
                        if self.ksc_phi(self.mass_at_angle(mid_u))[1] > 0.0:
                            lo_u = mid_u
                        else:
                            hi_u = mid_u
                    u_pol = 0.5 * (lo_u + hi_u)
                    f_pol = self.curve(self.mass_at_angle(u_pol))
                    if f_pol <= fu + 1e-12 * max(1.0, abs(fu)):
                        fu, u_best = min(f_pol, fu), u_pol
            refined.append((fu, u_best))
        refined.append((vals[0], 0.0))
        refined.append((vals[-1], 0.5 * math.pi))
        best_val = min(x[0] for x in refined)
        tie_tol = 1e-9 * max(1.0, abs(best_val))
        global_us = []
        for fu, u in sorted(refined, key=lambda x: x[1]):
            if fu <= best_val + tie_tol:
                if not global_us or abs(u - global_us[-1]) > 1e-6 * max(1e-6, u, 0.5 * math.pi - u):
                    global_us.append(u)
        masses = tuple(
            math.inf if u >= 0.5 * math.pi - 1e-14 else (0.0 if u <= 1e-14 else self.scale * math.tan(u))
            for u in global_us
        )
        if len(masses) > 1:
            return dict(
                mass=self.scale, status="tied", convention=True,
                K_min=best_val, minima=masses, K_0=vals[0], K_inf=vals[-1],
            )
        m_out = masses[0]
        if m_out == 0.0:
            status = "boundary_zero"
        elif math.isinf(m_out):
            status = "boundary_inf"
        else:
            status = "interior"
        return dict(
            mass=m_out, status=status, convention=False,
            K_min=best_val, minima=masses, K_0=vals[0], K_inf=vals[-1],
        )

    def kept_intervals(self, *, level=0.95, kind="floor", grid=361):
        """Exact connected components of {m in [0, inf] : p(m) > 1 - level} on the angle arc,
        resolving narrow intervals around argmin_K and disconnected pieces (including hilltop ghosts)."""
        if not 0.0 < level < 1.0:
            raise ValueError("level must lie strictly in (0, 1)")
        alpha = 1.0 - level
        key = {"floor": "p_value", "t45": "p_t45", "tjoint": "p_tjoint"}[kind]

        def p_at_u(u):
            return self.p_value(self.mass_at_angle(u))[key]

        pt = self.argmin_K(grid=grid)
        u_pts = self._base_angles(grid=grid)
        # Local relative-mass stencil around every K minimizer (resolves widths down to 1e-7 * m*)
        for m_loc in pt["minima"]:
            if 0.0 < m_loc < math.inf:
                for exp in (-7, -6, -5, -4, -3, -2, -1):
                    step = 10.0 ** exp
                    for sgn in (-3, -2, -1, 0, 1, 2, 3):
                        mc = m_loc * (1.0 + sgn * step)
                        if mc > 0.0:
                            u_pts.append(self.angle_of(mc))
            else:
                u0 = self.angle_of(m_loc)
                for exp in (-7, -5, -3, -2):
                    for sgn in (-1, 1):
                        uu = u0 + sgn * (10.0 ** exp)
                        if 0.0 <= uu <= 0.5 * math.pi:
                            u_pts.append(uu)
        u_pts = sorted(set(u_pts))
        ps = [p_at_u(u) for u in u_pts]
        return _extract_intervals(u_pts, ps, p_at_u, alpha, self.mass_at_angle)


def _golden_min(fn, a, b, steps=48):
    c = b - _GR * (b - a)
    d = a + _GR * (b - a)
    fc, fd = fn(c), fn(d)
    for _ in range(steps):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - _GR * (b - a)
            fc = fn(c)
        else:
            a, c, fc = c, d, fd
            d = a + _GR * (b - a)
            fd = fn(d)
    return (fc, c) if fc < fd else (fd, d)


def _extract_intervals(u_pts, ps, p_at_u, alpha, u_to_m):
    """Extract exact connected kept intervals (lo, hi) on [0, inf] from sampled (u_pts, ps)
    by bisecting every crossing of p(u) = alpha to ~1e-13."""
    kept = [p > alpha for p in ps]
    if all(kept):
        return ((0.0, math.inf),)
    if not any(kept):
        return ()

    def bisect_root(u_lo, u_hi, k_lo):
        for _ in range(52):
            mid = 0.5 * (u_lo + u_hi)
            if (p_at_u(mid) > alpha) == k_lo:
                u_lo = mid
            else:
                u_hi = mid
        return 0.5 * (u_lo + u_hi)

    intervals = []
    in_seg = kept[0]
    seg_start = u_pts[0] if in_seg else None
    for i in range(len(u_pts) - 1):
        if kept[i] != kept[i + 1]:
            root = bisect_root(u_pts[i], u_pts[i + 1], kept[i])
            if not in_seg:
                seg_start = root
                in_seg = True
            else:
                intervals.append((u_to_m(seg_start), u_to_m(root)))
                in_seg = False
                seg_start = None
    if in_seg:
        intervals.append((u_to_m(seg_start), u_to_m(u_pts[-1])))
    # Merge intervals separated by negligible numerical sliver (only when both endpoints are finite!)
    merged = []
    for lo, hi in intervals:
        if (merged and math.isfinite(lo) and math.isfinite(merged[-1][1])
                and lo <= merged[-1][1] + 1e-11 * max(1.0, lo)):
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return tuple(merged)


def classify_shape(intervals):
    """Classify compatible set shape per T80 (3): 'empty', 'all', 'at_least', 'at_most',
    'interval', or 'disconnected'."""
    if not intervals:
        return "empty"
    if len(intervals) > 1:
        return "disconnected"
    lo, hi = intervals[0]
    if lo == 0.0 and math.isinf(hi):
        return "all"
    if lo == 0.0:
        return "at_most"
    if math.isinf(hi):
        return "at_least"
    return "interval"


# =============================================================================
# Public ours34 probe entry point
# =============================================================================

def ours34(readings, *, level=0.95, scale=None, single_scalar_policy="t45", grid=361):
    """Run the ours34 probe on a series of readings.

    Takes NO starting law pi and NO handed-over nu_eff.
    Returns a dict with:
      point, point_status, point_is_convention, minima,
      set (floor: T45 U t_joint), raw_set (T45 alone), t_joint_set, shape,
      leftover (K_min, M_nz, M_total, v_hat_at_point, chi2_bound),
      state (additive State for combining blocks or composing bodies), level.
    """
    st = State(readings, scale=scale, single_scalar_policy=single_scalar_policy)
    pt = st.argmin_K(grid=grid)
    floor_set = st.kept_intervals(level=level, kind="floor", grid=grid)
    raw_set = st.kept_intervals(level=level, kind="t45", grid=grid)
    tj_set = st.kept_intervals(level=level, kind="tjoint", grid=grid)
    at_pt = st.p_value(pt["mass"])
    leftover = dict(
        K_min=pt["K_min"],
        M_nz=st.M_nz,
        M_total=st.M_total,
        v_hat_at_point=at_pt["v_hat"],
        chi2_bound=chi2_sf(pt["K_min"], st.M_total) if st.M_total > 0 else 1.0,
    )
    return dict(
        point=pt["mass"],
        point_status=pt["status"],
        point_is_convention=pt["convention"],
        minima=pt["minima"],
        set=floor_set,
        raw_set=raw_set,
        t_joint_set=tj_set,
        shape=classify_shape(floor_set),
        raw_shape=classify_shape(raw_set),
        leftover=leftover,
        state=st,
        level=float(level),
    )


# =============================================================================
# Two-body composition: lump (L5) and common_force (L14)
# =============================================================================

def _extract_state(obj):
    if isinstance(obj, State):
        return obj
    if isinstance(obj, dict) and isinstance(obj.get("state"), State):
        return obj["state"]
    return State(obj)


def _p_to_chi1_z2(p):
    if p >= 1.0:
        return 0.0
    if p <= 0.0:
        return math.inf
    z = inv_norm_cdf(1.0 - 0.5 * p)
    return z * z


def _qpos_statistic(z1, z2, a1, a2, m1, m2):
    L2 = a1 * a1 + a2 * a2
    if L2 <= 0.0:
        return 0.0
    d = max(-m1, min(m2, (a1 * z1 - a2 * z2) / L2))
    c_lump = (z1 - a1 * d) ** 2 + (z2 + a2 * d) ** 2
    c1 = max(0.0, -a1 * m1 - z1) ** 2
    c2 = max(0.0, -a2 * m2 - z2) ** 2
    return max(0.0, c_lump - c1 - c2)


def _qpos_sf_exact(q, a1, a2, m1, m2):
    if q <= 0.0:
        return 1.0
    if q > 80.0:
        # Q_pos <= Z1^2 + Z2^2 ~ chi2_2 => P(Q_pos >= q) <= exp(-q/2) < 4.3e-18
        return 0.0
    r = math.sqrt(q)
    A1, A2, L2 = a1 * m1, a2 * m2, a1 * a1 + a2 * a2
    # Exact interior regime: when both positivity boundaries lie > 9 sigma beyond r, Q_pos == W_perp^2 ~ chi2_1
    L = math.sqrt(L2)
    if min(A1, A2, L * m1, L * m2) > r + 9.0:
        return math.erfc(r / math.sqrt(2.0))

    def p_z2_ge_q(z1):
        bps = sorted({-1e9, -A2, (a1 * z1 - L2 * m2) / a2, (a1 * z1 + L2 * m1) / a2, 1e9})
        c1 = max(0.0, -A1 - z1) ** 2
        prob = 0.0
        for lo, hi in zip(bps[:-1], bps[1:]):
            if hi - lo < 1e-14:
                continue
            mid = 0.5 * (max(lo, -1e6) + min(hi, 1e6))
            d_star = (a1 * z1 - a2 * mid) / L2
            if d_star <= -m1:
                al, be, ga = 1.0, -2.0 * a2 * m1, (z1 + A1) ** 2 + (a2 * m1) ** 2
            elif d_star >= m2:
                al, be, ga = 1.0, 2.0 * A2, (z1 - a1 * m2) ** 2 + A2 * A2
            else:
                al = a1 * a1 / L2
                be = 2.0 * a1 * a2 * z1 / L2
                ga = a2 * a2 * z1 * z1 / L2
            if mid < -A2:
                al -= 1.0
                be -= 2.0 * A2
                ga -= A2 * A2
            ga -= (c1 + q)
            if abs(al) < 1e-14:
                if abs(be) < 1e-14:
                    if ga >= 0.0:
                        prob += norm_cdf(hi) - norm_cdf(lo)
                else:
                    root = -ga / be
                    a_, b_ = (max(lo, root), hi) if be > 0.0 else (lo, min(hi, root))
                    if b_ > a_:
                        prob += norm_cdf(b_) - norm_cdf(a_)
            else:
                disc = be * be - 4.0 * al * ga
                if al > 0.0:
                    if disc <= 0.0:
                        prob += norm_cdf(hi) - norm_cdf(lo)
                    else:
                        s = math.sqrt(disc)
                        r1, r2 = (-be - s) / (2.0 * al), (-be + s) / (2.0 * al)
                        if r1 > lo:
                            prob += norm_cdf(min(hi, r1)) - norm_cdf(lo)
                        if r2 < hi:
                            prob += norm_cdf(hi) - norm_cdf(max(lo, r2))
                else:
                    if disc > 0.0:
                        s = math.sqrt(disc)
                        r1, r2 = (-be + s) / (2.0 * al), (-be - s) / (2.0 * al)
                        a_, b_ = max(lo, r1), min(hi, r2)
                        if b_ > a_:
                            prob += norm_cdf(b_) - norm_cdf(a_)
        return prob

    cand = [-9.0, 9.0, 0.0, -A1, r, -r, -A1 - r, -A1 + r, a1 * m2 - r, a1 * m2 + r]
    z1_bps = sorted({x for x in cand if -9.0 <= x <= 9.0})
    total = 0.0
    for a_, b_ in zip(z1_bps[:-1], z1_bps[1:]):
        if b_ - a_ < 1e-14:
            continue
        for x, w in zip(_GL_X, _GL_W):
            z1 = a_ + (b_ - a_) * x
            total += w * (b_ - a_) * norm_pdf(z1) * p_z2_ge_q(z1)
    return max(0.0, min(1.0, total))


def lump_curve(s1, s2, M, grid=201):
    """T74 infimal convolution K_lump(M) = inf_{m1 + m2 = M} [K1(m1) + K2(m2)]."""
    if math.isinf(M):
        return min(s1.curve(math.inf) + s2.argmin_K()["K_min"],
                   s1.argmin_K()["K_min"] + s2.curve(math.inf))
    if M == 0.0:
        return s1.curve(0.0) + s2.curve(0.0)
    def f(t):
        return s1.curve(t * M) + s2.curve((1.0 - t) * M)
    ts = [i / (grid - 1) for i in range(grid)]
    vs = [f(t) for t in ts]
    j = min(range(grid), key=lambda i: vs[i])
    best, _ = _golden_min(f, ts[max(0, j - 1)], ts[min(grid - 1, j + 1)], steps=36)
    return min(best, vs[j])


def _lump_p_split(s1, s2, M, m1, method):
    m2 = M - m1
    if method in ("joint", "t45_joint"):
        r1 = s1.p_value(m1)
        r2 = s2.p_value(m2)
        c1 = r1["C_phi"]
        c2 = r2["C_phi"]
        if c1 <= 1e-30 and c2 <= 1e-30:
            return 1.0
        p1 = r1["p_value"] if method == "joint" else r1["p_t45"]
        p2 = r2["p_value"] if method == "joint" else r2["p_t45"]
        if c2 <= 1e-30:
            return p1
        if c1 <= 1e-30:
            return p2
        return math.exp(-0.5 * (_p_to_chi1_z2(p1) + _p_to_chi1_z2(p2)))
    if method == "qpos":
        _, S1, C1 = s1.ksc_m(m1)
        _, S2, C2 = s2.ksc_m(m2)
        if C1 <= 0.0 and C2 <= 0.0:
            return 1.0
        if C2 <= 0.0:
            a1 = math.sqrt(C1)
            q = max(0.0, S1 / a1 - a1 * m2) ** 2
            return 1.0 if q == 0.0 else 1.0 - norm_cdf(a1 * m2 + math.sqrt(q))
        if C1 <= 0.0:
            a2 = math.sqrt(C2)
            q = max(0.0, S2 / a2 - a2 * m1) ** 2
            return 1.0 if q == 0.0 else 1.0 - norm_cdf(a2 * m1 + math.sqrt(q))
        a1, a2 = math.sqrt(C1), math.sqrt(C2)
        return _qpos_sf_exact(_qpos_statistic(S1 / a1, S2 / a2, a1, a2, m1, m2), a1, a2, m1, m2)
    raise ValueError("method must be 'joint', 't45_joint', or 'qpos'")


def _lump_p_total(s1, s2, M, pt1, pt2, method, n_splits=31):
    if M <= 0.0:
        return _lump_p_split(s1, s2, 0.0, 0.0, method) if method != "qpos" else 0.0
    if math.isinf(M):
        # At M = inf, one part is inf and the other can sit anywhere in (0, inf)
        if method == "qpos":
            return max(s1.p_value(math.inf)["p_t45"], s2.p_value(math.inf)["p_t45"])
        key = "p_value" if method == "joint" else "p_t45"
        p1_inf = s1.p_value(math.inf)[key]
        p2_inf = s2.p_value(math.inf)[key]
        p1_best = max([s1.p_value(m)[key] for m in pt1["minima"]] or [1.0])
        p2_best = max([s2.p_value(m)[key] for m in pt2["minima"]] or [1.0])
        c1_inf = s1.ksc_phi(math.inf)[2]
        c2_inf = s2.ksc_phi(math.inf)[2]
        cand = []
        cand.append(p1_inf if s2.M_nz == 0 else math.exp(-0.5 * (_p_to_chi1_z2(p1_inf) + _p_to_chi1_z2(p2_best))))
        cand.append(p2_inf if s1.M_nz == 0 else math.exp(-0.5 * (_p_to_chi1_z2(p1_best) + _p_to_chi1_z2(p2_inf))))
        if c1_inf <= 1e-30 and c2_inf <= 1e-30:
            cand.append(1.0)
        return max(cand)

    shares = [1e-6, 1e-4, 1e-2] + [i / n_splits for i in range(1, n_splits)] + [1.0 - 1e-2, 1.0 - 1e-4, 1.0 - 1e-6]
    pts = [t * M for t in shares]
    for m_loc in pt1["minima"]:
        if 0.0 < m_loc < M:
            pts.append(m_loc)
            for delta in (-1e-3, 1e-3, -1e-2, 1e-2):
                cand_m = m_loc * (1.0 + delta)
                if 0.0 < cand_m < M:
                    pts.append(cand_m)
    for m_loc in pt2["minima"]:
        if 0.0 < m_loc < M:
            m1_c = M - m_loc
            pts.append(m1_c)
            for delta in (-1e-3, 1e-3, -1e-2, 1e-2):
                cand_m = M - m_loc * (1.0 + delta)
                if 0.0 < cand_m < M:
                    pts.append(cand_m)
    if not pt1["convention"] and not pt2["convention"]:
        m1s, m2s = pt1["mass"], pt2["mass"]
        if 0.0 < m1s < math.inf and 0.0 < m2s < math.inf:
            pts.append(M * m1s / (m1s + m2s))
    vals = [_lump_p_split(s1, s2, M, x, method) for x in pts]
    best = max(vals)
    if best >= 0.999999:
        return best
    idx = max(range(len(pts)), key=lambda i: vals[i])
    x0 = pts[idx]
    lo = max(1e-9 * M, x0 - 0.05 * M)
    hi = min((1.0 - 1e-9) * M, x0 + 0.05 * M)
    if best > 1e-12:
        neg_best, _ = _golden_min(lambda x: -_lump_p_split(s1, s2, M, x, method), lo, hi, steps=18)
        best = max(best, -neg_best)
    elif method == "qpos":
        # Check for narrow nonlinear cancellation split where S1/C1 + S2/C2 crosses zero
        s_pts = sorted(set(pts))
        prev_g = None
        prev_x = None
        for x in s_pts:
            _, s1v, c1v = s1.ksc_m(x)
            _, s2v, c2v = s2.ksc_m(M - x)
            if c1v > 0.0 and c2v > 0.0:
                g = s1v / c1v + s2v / c2v
                if prev_g is not None and prev_g * g <= 0.0:
                    a_, b_, ga = prev_x, x, prev_g
                    for _ in range(32):
                        mid = 0.5 * (a_ + b_)
                        _, sm1, cm1 = s1.ksc_m(mid)
                        _, sm2, cm2 = s2.ksc_m(M - mid)
                        gm = sm1 / cm1 + sm2 / cm2
                        if gm * ga <= 0.0:
                            b_ = mid
                        else:
                            a_, ga = mid, gm
                    best = max(best, _lump_p_split(s1, s2, M, 0.5 * (a_ + b_), method))
                prev_g, prev_x = g, x
    return best


def lump(first, second, *, level=0.95, method="joint", grid=161, n_splits=31):
    """L5 two-body lumping M = m1 + m2 on ours34 states/results.

    `method="joint"` (default): 2-dof curve-pair composition combining the two bodies'
    independent floor p-values (zero MC, ghost-free, discordant-safe).
    `method="t45_joint"`: same using Gaussian T45 p-values.
    `method="qpos"`: 1-dof linearized Gaussian Q_pos (Theorem E14-4; sharp 1-dof on near-Gaussian
    curves, vulnerable to nonlinear fake-split ghosts on asymmetric curves).
    """
    s1 = _extract_state(first)
    s2 = _extract_state(second)
    pt1 = s1.argmin_K()
    pt2 = s2.argmin_K()
    scale = s1.scale + s2.scale

    if pt1["convention"] or pt2["convention"]:
        point = scale
        point_status = "convention"
        point_is_convention = True
    else:
        point = pt1["mass"] + pt2["mass"]
        point_is_convention = False
        if math.isinf(point):
            point_status = "boundary_inf"
        elif point == 0.0:
            point_status = "boundary_zero"
        else:
            point_status = "interior"

    def u_to_m(u):
        return 0.0 if u <= 0.0 else (math.inf if u >= 0.5 * math.pi else scale * math.tan(u))

    def p_value(M):
        return _lump_p_total(s1, s2, M, pt1, pt2, method, n_splits=n_splits)

    u_pts = [0.5 * math.pi * i / (grid - 1) for i in range(grid)]
    for k in range(-96, 97):
        u_pts.append(math.atan(10.0 ** (k / 12.0)))
    # Add local stencils around finite point AND around each body's finite minimizer (crucial for 5+/-0.5 + 0/0!)
    anchors = []
    if not point_is_convention and 0.0 < point < math.inf:
        anchors.append(point)
    for pt in (pt1, pt2):
        for m_loc in pt["minima"]:
            if 0.0 < m_loc < math.inf:
                anchors.append(m_loc)
    for anc in anchors:
        for exp in (-5, -4, -3, -2, -1):
            step = 10.0 ** exp
            for sgn in (-3, -2, -1, 0, 1, 2, 3):
                mc = anc * (1.0 + sgn * step)
                if mc > 0.0:
                    u_pts.append(math.atan(mc / scale))
    u_pts = sorted(set(u_pts))
    p_at_u = lambda u: p_value(u_to_m(u))
    ps = [p_at_u(u) for u in u_pts]
    kept = _extract_intervals(u_pts, ps, p_at_u, 1.0 - level, u_to_m)

    return dict(
        point=point,
        point_status=point_status,
        point_is_convention=point_is_convention,
        set=kept,
        shape=classify_shape(kept),
        curve=lambda M: lump_curve(s1, s2, M),
        p_value=p_value,
        method=method,
        level=float(level),
    )


def common_force(first, second, *, level=0.95, method="joint", grid=161, n_splits=31):
    """L14 common-force composition 1/M = 1/m1 + 1/m2 via exact T5 channel-swap duality."""
    s1 = _extract_state(first)
    s2 = _extract_state(second)
    s1_sw = State(tuple(r.swapped() for r in s1.readings), single_scalar_policy=s1.single_scalar_policy)
    s2_sw = State(tuple(r.swapped() for r in s2.readings), single_scalar_policy=s2.single_scalar_policy)
    sw = lump(s1_sw, s2_sw, level=level, method=method, grid=grid, n_splits=n_splits)

    def recip(x):
        return math.inf if x == 0.0 else (0.0 if math.isinf(x) else 1.0 / x)

    inv_set = tuple(sorted((recip(hi), recip(lo)) for lo, hi in sw["set"]))
    point = recip(sw["point"]) if not sw["point_is_convention"] else 1.0 / (1.0 / s1.scale + 1.0 / s2.scale)
    st = sw["point_status"]
    if st == "boundary_inf":
        st = "boundary_zero"
    elif st == "boundary_zero":
        st = "boundary_inf"
    return dict(
        point=point,
        point_status=st,
        point_is_convention=sw["point_is_convention"],
        set=inv_set,
        shape=classify_shape(inv_set),
        curve=lambda M: sw["curve"](recip(M)),
        p_value=lambda M: sw["p_value"](recip(M)),
        method=method,
        level=float(level),
    )
