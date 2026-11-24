"""ours35: Door A1 probe — pure Gaussian-believed unnormalized evidence on the full reading type.

NOT the definition. Taken across from exploration branch `explore/e16-ours35-door-a` (24 Nov 2026,
Antigravity with Can; R-2026-11-24-1; not yet checked by Claude) as a probe beside `ours33` and
`ours34`. Implements Door A1 from `explore/e15-door-a-uniqueness`:

  - Believed declared Gaussian error covariance Sigma_i (L13; no nu_eff, no set-widening discount)
  - No starting law pi (T80, R-2026-11-22-1)
  - Contradiction reported strictly in the leftover K_min beside the compatible set (T73 (5))
  - Pure additive carried state State = (K, C) with zero hidden reading count M or M_nz (T77 (6), L25)

## The equation

Each reading i in the full type (data y_i = (F_i; a_i) in R^{2d}, SPD covariance Sigma_i,
orthonormal candidate subspace N_i in R^{d x k_i}) is pre-whitened at its own natural scale
s_i = sqrt(tr(Sigma_{FF,i}) / tr(Sigma_{aa,i})) for full dynamic-range stability (even at
sigma_F / sigma_a = 1e+-100). At every candidate mass m on the closed ray [0, inf]:

  K_i(m) = ||(I - P_i(m)) z_i||^2                   (T36 / T73 unexplained residual)
  S_i(m) = -(1/2) dK_i/dtheta = V_i^T (P_i'(theta) U_i)   (T45 conditional score)
  C_i(m) = ||P_i'(theta) U_i||^2 = Var(S_i | U_i)   (T45 conditional variance)

Additive carried state across independent blocks of one body (T77 (6)):

  State = (K, C),    K = sum_i K_i,   C = sum_i C_i   (exact additive monoid; no M or M_nz).

Exact 0/0 readings have (K_i, C_i) = (0, 0) identically (exact Z, continuous at 0/0, and
exact L17 reading-splitting invariance).

Point (T81):
  m* = argmin_{m in [0, inf]} K(m). On strictly flat or tied curves, returns the declared
  covariance scale with `point_is_convention=True` (R-2026-11-22-2). Boundary minima at 0
  or inf are marked `boundary_zero` / `boundary_inf` (pending Can's ruling on T81 (5)).

Compatible set at declared level 1 - alpha (T45; uniquely UMPU/LMPU under Door A by E15-1 & E15-2):
  Z_T45(m) = S(m) / sqrt(C(m)) ~ N(0, 1) exact given U whenever C(m) > 0 (p = 1 at C = 0).
  Keep m iff |Z_T45(m)| < z_{1 - alpha/2} (equivalently p_T45(m) > alpha).
  Computed by exact S = 0 stationary-point discovery on a swap-symmetric two-half coordinate
  (r = m/scale in [0, 1] and t = scale/m in [0, 1]), so ultra-narrow strong-signal intervals,
  extreme scales, hilltop pieces, and exact single-reading pinholes are all resolved with
  exact L3 unit equivariance and T5 channel-swap symmetry.

Leftover beside the set (T73 (5)):
  leftover = K_min = min_{m in [0, inf]} K(m). Purely a function of K (no dof count M or M_nz),
  so every single output field is exact Z, continuous at 0/0, and L17 reading-split invariant.

Two-body composition (`lump` for L5, `common_force` for L14):
  Curve: T74 infimal convolution K_lump(M) = inf_{m1 + m2 = M} [K1(m1) + K2(m2)].
  Point: m1* + m2* (and harmonic dual) on [0, inf] (T81).
  Set:   Exact curve-pair `Q_joint` from (K, C) alone: sup_{m1 + m2 = M} p_joint(m1, m2), where
         p_joint = exp(-(Z1^2 + Z2^2)/2) (chi2_2 survival when both C_j > 0, collapsing cleanly
         to 1-dof p_T45 when either body has C_j = 0, e.g. 0/0). Zero Monte Carlo, ghost-free,
         uniformly exact coverage >= 1 - alpha everywhere. Known-mass shifts `lump_known` and
         `common_force_known` are exact 1-dof.

## Known provable costs of Door A1 (theorems, not bugs; surfaces for Can's ruling)

  1. No set widening on disagreement (Door A choice): wildly conflicting readings produce a
     narrow set around the compromise mass with a screaming leftover K_min beside it (E15 / T82).
  2. Hilltop ghosts, right-angle pinholes, and anti-aligned boundary-end rejection (Theorem E15-1
     & E15-3): unavoidable on (K, C) with exact size == alpha and continuity at 0/0. On anti-aligned
     (10, -3) & (10, -3.1), T81's point is inf (K_min = 18.61) while T45 rejects inf (p = 1.6e-5)
     and keeps only the hilltop around the global maximum m = 0.305 (K = 218.61).
  3. 2-dof width tax (+24.9% at strong signal, 2.45 sigma vs 1.96 sigma) in two-body `lump` /
     `common_force` from (K, C) alone (G27; vs 1-dof raw-reading MC T74 (6)).

Zero external dependencies (pure Python 3 stdlib). Duck-typed to accept `met.full.Reading`,
`met.ours34.Reading`, `met.ours35.Reading`, or any object with `.force`, `.acceleration`, `.cov`,
and optional `.subspace` / `.window`.
"""
from dataclasses import dataclass
import math


_GR = 0.5 * (math.sqrt(5.0) - 1.0)


# =============================================================================
# Pure-stdlib linear algebra (exact for d in {1, 2, 3}, k in {1, 2, 3})
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
    k = len(b)
    if k == 1:
        return (b[0] / G[0][0],)
    if k == 2:
        g00, g01, g11 = G[0][0], G[0][1], G[1][1]
        det = g00 * g11 - g01 * g01
        return ((g11 * b[0] - g01 * b[1]) / det, (-g01 * b[0] + g00 * b[1]) / det)
    L = _cholesky(G)
    y = _solve_lower(L, b)
    x = [0.0] * k
    for i in range(k - 1, -1, -1):
        x[i] = (y[i] - math.fsum(L[j][i] * x[j] for j in range(i + 1, k))) / L[i][i]
    return tuple(x)


def _orthonormal_basis(cols):
    """Two-pass modified Gram-Schmidt for the 2k-dim whitened span basis (computed once)."""
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
# Full-type Reading & exact branch (T1)
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
    """Full-type reading matching `met.full.Reading` without numpy."""
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
        # Validate positive definiteness on diagonally scaled cov so extreme scales (1e+-150) pass cleanly
        diag_sd = []
        for i in range(2 * d):
            if cov[i][i] <= 0.0:
                raise ValueError("covariance matrix must be positive definite (L2)")
            diag_sd.append(math.sqrt(cov[i][i]))
        corr_mat = [[cov[i][j] / (diag_sd[i] * diag_sd[j]) for j in range(2 * d)] for i in range(2 * d)]
        _cholesky(corr_mat)
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
        """Scale force units by `lam` (L3): m -> lam * m."""
        lam = float(lam)
        if not (lam > 0.0 and math.isfinite(lam)):
            raise ValueError("unit scale lam must be positive and finite")
        d = self.dimension
        s = [lam] * d + [1.0] * d
        new_cov = tuple(tuple(self.cov[i][j] * s[i] * s[j] for j in range(2 * d)) for i in range(2 * d))
        return Reading(tuple(lam * x for x in self.force), self.acceleration, new_cov, self.subspace, self.window)

    def rotated(self, R):
        """Rotate spatial vectors and subspace by orthogonal d x d matrix R (L4)."""
        d = self.dimension
        R = _to_2d_tuple(R)
        if len(R) != d or any(len(row) != d for row in R):
            raise ValueError(f"rotation must be {d} x {d}")
        rf = tuple(_dot(R[i], self.force) for i in range(d))
        ra = tuple(_dot(R[i], self.acceleration) for i in range(d))
        # Block-diagonal diag(R, R) @ cov @ diag(R, R)^T
        B = [[0.0] * (2 * d) for _ in range(2 * d)]
        for i in range(d):
            for j in range(d):
                B[i][j] = R[i][j]
                B[d + i][d + j] = R[i][j]
        tmp = [[_dot(B[i], [self.cov[t][j] for t in range(2 * d)]) for j in range(2 * d)] for i in range(2 * d)]
        new_cov = tuple(tuple(_dot(tmp[i], B[j]) for j in range(2 * d)) for i in range(2 * d))
        k = self.k
        new_sub = tuple(
            tuple(math.fsum(R[i][t] * self.subspace[t][j] for t in range(d)) for j in range(k))
            for i in range(d)
        )
        return Reading(rf, ra, new_cov, new_sub, self.window)


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
    """Accept `met.full.Reading`, `met.ours34.Reading`, `met.ours35.Reading`, or duck-typed reading."""
    if isinstance(r, Reading):
        return r
    return Reading(
        force=r.force,
        acceleration=r.acceleration,
        cov=r.cov,
        subspace=getattr(r, "subspace", None),
        window=getattr(r, "window", None),
    )


def exact_mass(force, acceleration, rtol=1e-12):
    """Exact branch (T1) for zero noise: returns ('one', m), ('none', reason), or ('every', None)."""
    f = _to_1d_tuple(force)
    a = _to_1d_tuple(acceleration)
    if len(f) != len(a) or len(f) not in (1, 2, 3):
        raise ValueError("force and acceleration must have equal length in {1, 2, 3}")
    ff, aa, fa = _dot(f, f), _dot(a, a), _dot(f, a)
    if ff == 0.0 and aa == 0.0:
        return ("every", None)
    if aa == 0.0:
        return ("none", "acceleration is zero and force is not: m2 is infinite")
    if ff == 0.0:
        return ("none", "force is zero and acceleration is not: m1 is zero")
    if fa * fa < ff * aa * (1.0 - rtol):
        return ("none", "not codirectional: the two pairings disagree (P < 1)")
    if fa <= 0.0:
        return ("none", "anti-aligned: the pairings agree on a negative value")
    return ("one", fa / aa)


def exact_series(pairs, rtol=1e-12):
    """Exact branch (T1) across a series of (force, acceleration) exact pairs."""
    found = None
    for f, a in pairs:
        tag, val = exact_mass(f, a, rtol=rtol)
        if tag == "none":
            return ("none", val)
        if tag == "one":
            if found is None:
                found = val
            elif abs(val - found) > rtol * max(val, found):
                return ("none", "two exact readings determine distinct positive masses")
    if found is None:
        return ("every", None)
    return ("one", found)


# =============================================================================
# Precomputed Scale-Normalized Whitened Span (zero-QR, overflow-free)
# =============================================================================

class _SpanPre:
    """Scale-normalized whitened span for one reading.

    Rescales (F, a, cov) internally by (f_scale, a_scale) = (sqrt(tr(Sigma_FF)/d), sqrt(tr(Sigma_aa)/d))
    so all Gram blocks are O(1) across arbitrary unit ratios s_i = f_scale / a_scale.
    Evaluates (K_i, S_{theta_i}, C_{theta_i}) at reading-local angle (p_i, q_i) where
    m = s_i * (p_i / q_i), with exact detection of physical right-angle pinholes (b = 0 => C = 0).
    """

    def __init__(self, reading):
        self.d = reading.dimension
        self.k = reading.k
        self.is_exact_zero = reading.is_exact_zero
        tr_f = math.fsum(reading.cov[i][i] for i in range(self.d))
        tr_a = math.fsum(reading.cov[self.d + i][self.d + i] for i in range(self.d))
        f_scale = math.sqrt(tr_f / self.d)
        a_scale = math.sqrt(tr_a / self.d)
        self.s_i = f_scale / a_scale
        if self.is_exact_zero:
            return
        scales = [f_scale] * self.d + [a_scale] * self.d
        cov_norm = [
            [reading.cov[i][j] / (scales[i] * scales[j]) for j in range(2 * self.d)]
            for i in range(2 * self.d)
        ]
        L = _cholesky(cov_norm)
        y = [reading.force[i] / f_scale for i in range(self.d)] + [
            reading.acceleration[i] / a_scale for i in range(self.d)
        ]
        z = _solve_lower(L, y)
        n_cols = [tuple(reading.subspace[row][col] for row in range(self.d)) for col in range(self.k)]
        span_cols = []
        for col in n_cols:
            span_cols.append(_solve_lower(L, list(col) + [0.0] * self.d))
        for col in n_cols:
            span_cols.append(_solve_lower(L, [0.0] * self.d + list(col)))
        W = _orthonormal_basis(span_cols)
        self.zW = tuple(_dot(w, z) for w in W)
        self.energy = _dot(self.zW, self.zW)
        self.bf = tuple(tuple(_dot(w, span_cols[j]) for w in W) for j in range(self.k))
        self.ba = tuple(tuple(_dot(w, span_cols[self.k + j]) for w in W) for j in range(self.k))
        k = self.k
        self.gff = tuple(tuple(_dot(self.bf[i], self.bf[j]) for j in range(k)) for i in range(k))
        self.gfa = tuple(tuple(_dot(self.bf[i], self.ba[j]) for j in range(k)) for i in range(k))
        self.gaa = tuple(tuple(_dot(self.ba[i], self.ba[j]) for j in range(k)) for i in range(k))
        self.gsym = tuple(tuple(self.gfa[i][j] + self.gfa[j][i] for j in range(k)) for i in range(k))
        self.uf = tuple(_dot(self.bf[i], self.zW) for i in range(k))
        self.ua = tuple(_dot(self.ba[i], self.zW) for i in range(k))

    def eval_local_pq(self, p, q):
        """Return (K_i, S_{theta_i}, C_{theta_i}) at unit (p, q) with m = s_i * (p / q)."""
        if self.is_exact_zero:
            return 0.0, 0.0, 0.0
        k = self.k
        pp, pq, qq = p * p, p * q, q * q
        if k == 1:
            b = p * self.uf[0] + q * self.ua[0]
            b_abs = abs(p * self.uf[0]) + abs(q * self.ua[0])
            if b == 0.0 or (b_abs > 0.0 and abs(b) <= 1e-14 * b_abs):
                # Exact physical right-angle pinhole (U = 0 => C = 0, S = 0, K = E)
                return self.energy, 0.0, 0.0
            g = pp * self.gff[0][0] + pq * self.gsym[0][0] + qq * self.gaa[0][0]
            alpha = b / g
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

        b = [p * self.uf[i] + q * self.ua[i] for i in range(k)]
        b_max = max(abs(x) for x in b)
        b_abs_max = max(abs(p * self.uf[i]) + abs(q * self.ua[i]) for i in range(k))
        if b_max == 0.0 or (b_abs_max > 0.0 and b_max <= 1e-14 * b_abs_max):
            return self.energy, 0.0, 0.0
        G = [[pp * self.gff[i][j] + pq * self.gsym[i][j] + qq * self.gaa[i][j] for j in range(k)]
             for i in range(k)]
        alpha = _solve_spd_small(G, b)
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

    def eval_in_scale(self, p_th, q_th, s_ref):
        """Return (K_i, S_{theta, i}, C_{theta, i}) in common series angle m = s_ref * (p_th / q_th)."""
        if self.is_exact_zero:
            return 0.0, 0.0, 0.0
        ratio = s_ref / self.s_i
        pr = p_th * ratio
        nrm = math.hypot(pr, q_th)
        p_i = pr / nrm
        q_i = q_th / nrm
        K, S_i, C_i = self.eval_local_pq(p_i, q_i)
        jac = ratio / (nrm * nrm)  # exact d(theta_i) / d(theta) > 0
        return K, jac * S_i, (jac * jac) * C_i


# =============================================================================
# Swap-symmetric two-half ray coordinates:
#   half == 0, x in [0, 1]  =>  m = scale * x     (x ascending: 0 -> 1)
#   half == 1, t in [0, 1]  =>  m = scale / t     (t descending: 1 -> 0, t=0 is inf)
# =============================================================================

def _half_to_pq(half, x):
    """Map two-half coordinate (half, x) to unit (p_th, q_th) with m = scale * p_th / q_th."""
    if half == 0:
        if x <= 0.0:
            return 0.0, 1.0
        inv = 1.0 / math.hypot(x, 1.0)
        return x * inv, inv
    if x <= 0.0:
        return 1.0, 0.0
    inv = 1.0 / math.hypot(1.0, x)
    return inv, x * inv


def _half_to_mass(half, x, scale):
    if half == 0:
        return 0.0 if x <= 0.0 else scale * x
    return math.inf if x <= 0.0 else scale / x


def _mass_to_half(m, scale):
    if m < 0.0 or math.isnan(m):
        raise ValueError("candidate mass must lie in [0, inf]")
    if math.isinf(m):
        return (1, 0.0)
    if m == 0.0:
        return (0, 0.0)
    if m <= scale:
        return (0, m / scale)
    return (1, scale / m)


def _mass_to_pq_in_scale(m, scale):
    half, x = _mass_to_half(m, scale)
    return _half_to_pq(half, x)


def _half_order_key(pt):
    """Monotonic sort key along [0, inf] from (0, 0.0) through (0, 1.0)==(1, 1.0) to (1, 0.0)."""
    half, x = pt
    return (0, x) if half == 0 else (1, -x)


def _bisect_scalar_root(fn, a, b, fa, fb, *, steps=60):
    """Bisect a sign change of `fn` on [a, b] (a <= b), using log-space narrowing if a == 0."""
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    lo, hi = a, b
    flo = fa
    # If left bracket is exact 0.0 and root is tiny (e.g. m ~ 1e-20 * scale), zoom geometrically first
    if lo == 0.0 and hi > 0.0:
        probe = hi * 1e-4
        while probe > 1e-300:
            fp = fn(probe)
            if fp == 0.0:
                return probe
            if (fp > 0.0) == (flo > 0.0):
                lo, flo = probe, fp
                break
            hi = probe
            probe *= 1e-4
    for _ in range(steps):
        if lo > 0.0 and hi / lo > 16.0:
            mid = math.sqrt(lo * hi)
        else:
            mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        fm = fn(mid)
        if fm == 0.0:
            return mid
        if (fm > 0.0) == (flo > 0.0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _bisect_bool_boundary(pred, a, b, pa, *, steps=60):
    """Bisect a True/False transition of `pred` on [a, b] (a <= b)."""
    lo, hi = a, b
    if lo == 0.0 and hi > 0.0:
        probe = hi * 1e-4
        while probe > 1e-300:
            if pred(probe) == pa:
                lo = probe
                break
            hi = probe
            probe *= 1e-4
    for _ in range(steps):
        if lo > 0.0 and hi / lo > 16.0:
            mid = math.sqrt(lo * hi)
        else:
            mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        if pred(mid) == pa:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


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


# =============================================================================
# Pure Additive State: (K, C) with zero hidden reading count
# =============================================================================

class State:
    """Pure additive state (K(m), C(m)) carried by ours35 (Door A1).

    Independent blocks of the same body combine by exact addition `state1 + state2`.
    Carries NO reading count M or M_nz, NO prior pi, and NO discount parameters.
    Exact 0/0 readings are stripped on construction (zero contribution to (K, C), exact Z).
    `scale` is purely a unit-equivariant coordinate anchor (and convention fallback on flat/tied curves).
    """

    def __init__(self, readings=(), *, scale=None, _spans=None, _all_scales=None):
        if _spans is None:
            adapted = tuple(adapt_reading(r) for r in readings)
            raw_spans = tuple(_SpanPre(r) for r in adapted)
            self._spans = tuple(sp for sp in raw_spans if not sp.is_exact_zero)
            all_scales = tuple(sp.s_i for sp in raw_spans)
        else:
            self._spans = tuple(sp for sp in _spans if not sp.is_exact_zero)
            all_scales = tuple(_all_scales) if _all_scales is not None else tuple(sp.s_i for sp in self._spans)
        self._all_scales = all_scales
        if scale is None:
            # Prefer non-zero spans so adding exact 0/0 does not even move the coordinate anchor (exact Z)
            src = tuple(sp.s_i for sp in self._spans) or all_scales
            if src:
                self.scale = math.exp(math.fsum(math.log(s) for s in src) / len(src))
            else:
                self.scale = 1.0
        else:
            self.scale = float(scale)
            if not (self.scale > 0.0 and math.isfinite(self.scale)):
                raise ValueError("scale must be positive and finite")

    @property
    def is_zero(self):
        """True iff (K, C) == (0, 0) identically (empty or all exact 0/0)."""
        return len(self._spans) == 0

    def __add__(self, other):
        if not isinstance(other, State):
            return NotImplemented
        return State(
            _spans=self._spans + other._spans,
            _all_scales=self._all_scales + other._all_scales,
        )

    def ksc_theta_pq(self, p_th, q_th):
        """Evaluate (K, S_theta, C_theta) at unit (p_th, q_th), where m = self.scale * (p_th / q_th)."""
        K = S = C = 0.0
        s_ref = self.scale
        for sp in self._spans:
            ki, si, ci = sp.eval_in_scale(p_th, q_th, s_ref)
            K += ki
            S += si
            C += ci
        return K, S, C

    def ksc_theta(self, m):
        """Return (K(m), S_theta(m), C_theta(m)) at m in [0, inf] (theta = atan(m / scale))."""
        p_th, q_th = _mass_to_pq_in_scale(m, self.scale)
        return self.ksc_theta_pq(p_th, q_th)

    def ksc_phi(self, m):
        """Return (K(m), S_phi(m), C_phi(m)) in raw angle phi = atan(m) (matches met.conditional_score)."""
        p_th, q_th = _mass_to_pq_in_scale(m, self.scale)
        K, S_th, C_th = self.ksc_theta_pq(p_th, q_th)
        s = self.scale
        dtheta_dphi = (q_th * q_th + (s * p_th) * (s * p_th)) / s
        return K, dtheta_dphi * S_th, (dtheta_dphi * dtheta_dphi) * C_th

    def ksc_m(self, m):
        """Return (K(m), S_m(m), C_m(m)) with respect to mass m (finite m >= 0)."""
        if math.isinf(m):
            K = self.curve(math.inf)
            return K, 0.0, 0.0
        p_th, q_th = _mass_to_pq_in_scale(m, self.scale)
        K, S_th, C_th = self.ksc_theta_pq(p_th, q_th)
        fac = (q_th * q_th) / self.scale
        return K, fac * S_th, (fac * fac) * C_th

    def curve(self, m):
        return self.ksc_theta(m)[0]

    def score_test(self, m):
        """Exact T45 conditional Gaussian score test at candidate mass m in [0, inf]."""
        p_th, q_th = _mass_to_pq_in_scale(m, self.scale)
        return self._score_test_pq(p_th, q_th)

    def _score_test_pq(self, p_th, q_th):
        K, S, C = self.ksc_theta_pq(p_th, q_th)
        if C <= 0.0:
            return dict(p_value=1.0, z=0.0, statistic=0.0, K=K, S_theta=S, C_theta=C, degenerate=True)
        z = S / math.sqrt(C)
        p_val = math.erfc(abs(z) / math.sqrt(2.0))
        return dict(p_value=p_val, z=z, statistic=z * z, K=K, S_theta=S, C_theta=C, degenerate=False)

    def p_value(self, m):
        return self.score_test(m)["p_value"]

    def _half_mesh(self, grid=241, log_decades=12, pts_per_decade=12):
        """Symmetric 1D mesh on [0, 1] used identically on both lower (m<=scale) and upper (m>=scale) halves."""
        xs = [i / (grid - 1) for i in range(grid)]
        for k in range(1, log_decades * pts_per_decade + 1):
            xs.append(10.0 ** (-k / pts_per_decade))
        return sorted(set(xs))

    def _eval_half(self, half, x):
        p_th, q_th = _half_to_pq(half, x)
        return self.ksc_theta_pq(p_th, q_th)

    def stationary_points(self, grid=241):
        """Find every interior stationary point (S_theta = 0) on (0, inf) via exact sign-change bisection
        on both swap-symmetric halves. Returns list of dicts:
          {'half', 'x', 'mass', 'K', 'C', 'kind': 'min'|'max'|'pinhole'|'saddle'}."""
        if self.is_zero:
            return []
        xs = self._half_mesh(grid=grid)
        found = []

        for half in (0, 1):
            # Evaluate K and S along x ascending from 0 to 1
            vals = [self._eval_half(half, x) for x in xs]
            k_min_h = min(v[0] for v in vals)
            k_max_h = max(v[0] for v in vals)
            if k_max_h <= 0.0 or k_max_h - k_min_h <= 1e-13 * k_max_h:
                continue

            # Collect candidate brackets for S = 0: every sign change of S, plus any K local extremum cell
            brackets = []
            for i in range(len(xs) - 1):
                s0, s1 = vals[i][1], vals[i + 1][1]
                if s0 == 0.0 and i > 0:
                    brackets.append((xs[i], xs[i], 0.0, 0.0))
                elif (s0 > 0.0 and s1 < 0.0) or (s0 < 0.0 and s1 > 0.0):
                    brackets.append((xs[i], xs[i + 1], s0, s1))
                elif i > 0 and i < len(xs) - 1:
                    # Guard against two close roots inside a single cell via K curvature check
                    k_prev, k_cur, k_next = vals[i - 1][0], vals[i][0], vals[i + 1][0]
                    if (k_cur < k_prev and k_cur < k_next) or (k_cur > k_prev and k_cur > k_next):
                        # Subdivide [xs[i-1], xs[i+1]] once on 9 points to catch double sign change
                        sub = [xs[i - 1] + (xs[i + 1] - xs[i - 1]) * j / 8.0 for j in range(9)]
                        ss = [self._eval_half(half, u)[1] for u in sub]
                        for j in range(8):
                            if (ss[j] > 0.0 and ss[j + 1] < 0.0) or (ss[j] < 0.0 and ss[j + 1] > 0.0):
                                brackets.append((sub[j], sub[j + 1], ss[j], ss[j + 1]))

            s_fn = lambda u, h=half: self._eval_half(h, u)[1]
            for a, b, sa, sb in brackets:
                if a == b:
                    xr = a
                else:
                    xr = _bisect_scalar_root(s_fn, a, b, sa, sb, steps=60)
                if xr <= 0.0:
                    continue
                if half == 1 and xr >= 1.0 - 1e-13:
                    # Boundary x=1 belongs to half=0 to avoid duplicating m = scale
                    continue
                K, S, C = self._eval_half(half, xr)
                # Classify via S sign to the left vs right in increasing mass m:
                # Recall S = -(1/2) dK/dtheta, so at a local minimum of K, S goes + -> - as m increases.
                eps = max(1e-12, 1e-8 * xr)
                x_lo = max(0.0, xr - eps)
                x_hi = min(1.0, xr + eps)
                if half == 0:
                    s_left = s_fn(x_lo)
                    s_right = s_fn(x_hi)
                else:
                    # On half == 1, increasing m corresponds to decreasing x!
                    s_left = s_fn(x_hi)
                    s_right = s_fn(x_lo)
                if C <= 0.0:
                    kind = "pinhole"
                elif s_left > 0.0 and s_right < 0.0:
                    kind = "min"
                elif s_left < 0.0 and s_right > 0.0:
                    kind = "max"
                else:
                    kind = "saddle"
                found.append(dict(
                    half=half, x=xr, mass=_half_to_mass(half, xr, self.scale),
                    K=K, C=C, kind=kind,
                ))

        # Also check exact midpoint m = scale ((0, 1.0)) if S changes sign across half 0 <-> half 1
        s_below = self._eval_half(0, 1.0 - 1e-8)[1]
        s_mid = self._eval_half(0, 1.0)[1]
        s_above = self._eval_half(1, 1.0 - 1e-8)[1]
        if abs(s_mid) <= 1e-13 * max(abs(s_below), abs(s_above)) or (s_below * s_above <= 0.0):
            K_m, _, C_m = self._eval_half(0, 1.0)
            kind = "pinhole" if C_m <= 0.0 else (
                "min" if (s_below > 0.0 and s_above < 0.0) else (
                    "max" if (s_below < 0.0 and s_above > 0.0) else "saddle"
                )
            )
            found.append(dict(half=0, x=1.0, mass=self.scale, K=K_m, C=C_m, kind=kind))

        # Sort monotonically along [0, inf] and deduplicate
        found.sort(key=lambda d: _half_order_key((d["half"], d["x"])))
        dedup = []
        for item in found:
            if not dedup:
                dedup.append(item)
                continue
            prev = dedup[-1]
            if prev["mass"] > 0.0 and item["mass"] > 0.0 and abs(item["mass"] - prev["mass"]) <= 1e-8 * max(prev["mass"], item["mass"]):
                if item["K"] < prev["K"] or (item["kind"] == "pinhole" and prev["kind"] != "pinhole"):
                    dedup[-1] = item
            else:
                dedup.append(item)
        return dedup

    def argmin_K(self, grid=241, tol=1e-12):
        """Global minimum of K on the closed ray [0, inf] (T81)."""
        if self.is_zero:
            return dict(
                mass=self.scale, status="flat", convention=True,
                K_min=0.0, minima=(), K_0=0.0, K_inf=0.0,
            )
        xs = self._half_mesh(grid=grid)
        k0 = self._eval_half(0, 0.0)[0]
        kinf = self._eval_half(1, 0.0)[0]
        all_k = [self._eval_half(0, x)[0] for x in xs] + [self._eval_half(1, x)[0] for x in xs]
        vmin, vmax = min(all_k), max(all_k)
        if vmax <= 0.0 or vmax - vmin <= tol * vmax:
            return dict(
                mass=self.scale, status="flat", convention=True,
                K_min=vmin, minima=(), K_0=k0, K_inf=kinf,
            )

        stats = self.stationary_points(grid=grid)
        cands = [(k0, (0, 0.0), 0.0), (kinf, (1, 0.0), math.inf)]
        for st in stats:
            cands.append((st["K"], (st["half"], st["x"]), st["mass"]))
        best_val = min(c[0] for c in cands)
        tie_tol = 1e-9 * max(vmax, abs(best_val))
        best_cands = [c for c in cands if c[0] <= best_val + tie_tol]
        best_cands.sort(key=lambda c: _half_order_key(c[1]))
        masses = []
        for _, _, m in best_cands:
            if not masses:
                masses.append(m)
            elif math.isinf(m) and math.isinf(masses[-1]):
                continue
            elif m == 0.0 and masses[-1] == 0.0:
                continue
            elif math.isfinite(m) and math.isfinite(masses[-1]) and abs(m - masses[-1]) <= 1e-7 * max(1.0, m, masses[-1]):
                continue
            else:
                masses.append(m)
        masses = tuple(masses)
        if len(masses) > 1:
            return dict(
                mass=self.scale, status="tied", convention=True,
                K_min=best_val, minima=masses, K_0=k0, K_inf=kinf,
            )
        m_out = masses[0]
        status = "boundary_zero" if m_out == 0.0 else ("boundary_inf" if math.isinf(m_out) else "interior")
        return dict(
            mass=m_out, status=status, convention=False,
            K_min=best_val, minima=masses, K_0=k0, K_inf=kinf,
        )

    def kept_intervals(self, *, level=0.95, grid=241):
        """Exact connected components of {m in [0, inf] : p_T45(m) > 1 - level}.

        Seeds the mesh with EVERY stationary point (local minima, maxima/hilltops, and pinholes)
        plus multi-scale local stencils around each, then bisects every p = alpha crossing to ~1e-14.
        Guarantees exact L3 unit equivariance and T5 swap invariance even on SNR = 1e10 narrow
        intervals and isolated single-reading pinholes.
        """
        if not 0.0 < level < 1.0:
            raise ValueError("level must lie strictly in (0, 1)")
        if self.is_zero:
            return ((0.0, math.inf),)
        alpha = 1.0 - level
        base = self._half_mesh(grid=grid)
        pts_half = {0: set(base), 1: set(base)}
        stats = self.stationary_points(grid=grid)
        for st in stats:
            h, x0 = st["half"], st["x"]
            pts_half[h].add(x0)
            for exp in (-12, -10, -8, -6, -4, -3, -2, -1):
                step = 10.0 ** exp
                for sgn in (-3, -1, 1, 3):
                    xn = x0 * (1.0 + sgn * step)
                    if 0.0 < xn <= 1.0:
                        pts_half[h].add(xn)
                    elif xn > 1.0:
                        # Crossed m = scale into the opposite half
                        pts_half[1 - h].add(1.0 / xn)

        def is_kept_half(h, x):
            p_th, q_th = _half_to_pq(h, x)
            res = self._score_test_pq(p_th, q_th)
            if x == 0.0 and res["degenerate"]:
                # Compatible set lives on the positive ray (0, inf) (L1): an isolated C=0 pinhole
                # strictly at the non-positive boundary 0 or inf is kept iff the one-sided limit
                # from (0, inf) is kept (fixes e.g. m=0 pinhole on (10, 0), whose set is at_least).
                p_lim, q_lim = _half_to_pq(h, 1e-50)
                res = self._score_test_pq(p_lim, q_lim)
            return res["p_value"] > alpha

        # Extract kept intervals separately on half 0 (x ascending) and half 1 (t descending), then stitch at m = scale
        segs = []
        for h in (0, 1):
            xs = sorted(pts_half[h])
            ks = [is_kept_half(h, x) for x in xs]
            if h == 1:
                xs = list(reversed(xs))
                ks = list(reversed(ks))
            in_seg = ks[0]
            start_m = _half_to_mass(h, xs[0], self.scale) if in_seg else None
            for i in range(len(xs) - 1):
                if ks[i] != ks[i + 1]:
                    if h == 0:
                        xr = _bisect_bool_boundary(lambda u: is_kept_half(0, u), xs[i], xs[i + 1], ks[i])
                    else:
                        # xs[i+1] < xs[i], with pred(xs[i+1]) == ks[i+1]
                        xr = _bisect_bool_boundary(lambda u: is_kept_half(1, u), xs[i + 1], xs[i], ks[i + 1])
                    mr = _half_to_mass(h, xr, self.scale)
                    if not in_seg:
                        start_m = mr
                        in_seg = True
                    else:
                        segs.append((start_m, mr))
                        in_seg = False
                        start_m = None
            if in_seg:
                segs.append((start_m, _half_to_mass(h, xs[-1], self.scale)))

        # Stitch intervals meeting at m = scale or separated by numerical sliver (keep genuine distinct intervals)
        merged = []
        for lo, hi in segs:
            if not merged:
                merged.append((lo, hi))
                continue
            plo, phi = merged[-1]
            if math.isfinite(phi) and math.isfinite(lo) and lo <= phi + 1e-11 * max(1e-300, phi, lo):
                merged[-1] = (plo, max(phi, hi))
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
# Public ours35 (Door A1) entry point
# =============================================================================

def ours35(readings, *, level=0.95, scale=None, grid=241):
    """Run the ours35 (Door A1) probe on a series of readings in the full type.

    Takes NO starting law pi, NO nu_eff, and NO discount parameter.
    Carries ONLY the additive state (K, C) (no M or M_nz).

    Returns a dict with:
      point, point_status, point_is_convention, minima,
      set (exact T45 compatible set at `level`), shape,
      leftover (pure K_min beside the set, Z-safe and L17-invariant),
      state (additive State for block combination and two-body composition),
      level.
    """
    st = State(readings, scale=scale)
    pt = st.argmin_K(grid=grid)
    kept = st.kept_intervals(level=level, grid=grid)
    return dict(
        point=pt["mass"],
        point_status=pt["status"],
        point_is_convention=pt["convention"],
        minima=pt["minima"],
        set=kept,
        shape=classify_shape(kept),
        leftover=pt["K_min"],
        state=st,
        level=float(level),
    )


# =============================================================================
# Two-body composition: lump (L5), common_force (L14), and known-mass shifts
# =============================================================================

def _extract_state(obj):
    if isinstance(obj, State):
        return obj
    if isinstance(obj, dict) and isinstance(obj.get("state"), State):
        return obj["state"]
    return State(obj)


def _recip(x):
    if x == 0.0:
        return math.inf
    if math.isinf(x):
        return 0.0
    return 1.0 / x


def lump_curve(s1, s2, M, grid=201):
    """T74 infimal convolution K_lump(M) = inf_{m1 + m2 = M} [K1(m1) + K2(m2)]."""
    s1 = _extract_state(s1)
    s2 = _extract_state(s2)
    if math.isinf(M):
        return min(
            s1.curve(math.inf) + s2.argmin_K()["K_min"],
            s1.argmin_K()["K_min"] + s2.curve(math.inf),
        )
    if M == 0.0:
        return s1.curve(0.0) + s2.curve(0.0)
    def f(t):
        return s1.curve(t * M) + s2.curve((1.0 - t) * M)
    ts = [i / (grid - 1) for i in range(grid)]
    for st in s1.stationary_points():
        if 0.0 < st["mass"] < M:
            ts.append(st["mass"] / M)
    for st in s2.stationary_points():
        if 0.0 < st["mass"] < M:
            ts.append(1.0 - st["mass"] / M)
    ts = sorted(set(ts))
    vs = [f(t) for t in ts]
    j = min(range(len(ts)), key=lambda i: vs[i])
    best, _ = _golden_min(f, ts[max(0, j - 1)], ts[min(len(ts) - 1, j + 1)], steps=36)
    return min(best, vs[j])


def _lump_p_split(s1, s2, m1, m2):
    """Exact Q_joint p-value at a candidate split (m1, m2):
      - C1 == 0 and C2 == 0 => 1.0
      - C2 == 0             => p_T45(s1, m1)   (1-dof exact, e.g. s2 is 0/0)
      - C1 == 0             => p_T45(s2, m2)   (1-dof exact, e.g. s1 is 0/0)
      - both C_j > 0        => exp(-(Z1^2 + Z2^2) / 2)   (exact chi2_2 survival)."""
    r1 = s1.score_test(m1)
    r2 = s2.score_test(m2)
    if r1["degenerate"] and r2["degenerate"]:
        return 1.0
    if r2["degenerate"]:
        return r1["p_value"]
    if r1["degenerate"]:
        return r2["p_value"]
    return math.exp(-0.5 * (r1["statistic"] + r2["statistic"]))


def _lump_p_total(s1, s2, M, stats1, stats2, pt1, pt2, n_splits=31):
    if M <= 0.0:
        return _lump_p_split(s1, s2, 0.0, 0.0)
    if math.isinf(M):
        # At M = inf, one body is inf and the other sits at any m in [0, inf]
        cands_1 = [0.0, math.inf] + [st["mass"] for st in stats1]
        cands_2 = [0.0, math.inf] + [st["mass"] for st in stats2]
        if s1.is_zero:
            cands_1.append(s1.scale)
        if s2.is_zero:
            cands_2.append(s2.scale)
        best = max(_lump_p_split(s1, s2, math.inf, m2) for m2 in cands_2)
        best = max(best, max(_lump_p_split(s1, s2, m1, math.inf) for m1 in cands_1))
        return best

    # Exact 1-dof open-interval supremum when one body is identically 0/0:
    # sup_{0 < m1 < M} p_T45(s1, m1) checked across endpoints, stationary points, and interior scan.
    if s1.is_zero and s2.is_zero:
        return 1.0
    shares = [1e-9, 1e-6, 1e-4, 1e-2] + [i / n_splits for i in range(1, n_splits)] + [1.0 - 1e-2, 1.0 - 1e-4, 1.0 - 1e-6, 1.0 - 1e-9]
    if s2.is_zero:
        cands = [0.0, M] + [t * M for t in shares] + [st["mass"] for st in stats1 if 0.0 < st["mass"] < M]
        return max(s1.p_value(x) for x in cands)
    if s1.is_zero:
        cands = [0.0, M] + [t * M for t in shares] + [st["mass"] for st in stats2 if 0.0 < st["mass"] < M]
        return max(s2.p_value(x) for x in cands)
    pts = [t * M for t in shares]
    for st in stats1:
        m_loc = st["mass"]
        if 0.0 < m_loc < M:
            pts.append(m_loc)
            for delta in (-1e-5, 1e-5, -1e-3, 1e-3, -1e-2, 1e-2):
                cm = m_loc * (1.0 + delta)
                if 0.0 < cm < M:
                    pts.append(cm)
    for st in stats2:
        m_loc = st["mass"]
        if 0.0 < m_loc < M:
            m1_c = M - m_loc
            pts.append(m1_c)
            for delta in (-1e-5, 1e-5, -1e-3, 1e-3, -1e-2, 1e-2):
                cm = M - m_loc * (1.0 + delta)
                if 0.0 < cm < M:
                    pts.append(cm)
    # Analytic local variance-weighted optimal split around every finite stationary-point pair
    for st1 in stats1:
        m1s = st1["mass"]
        if not (0.0 < m1s < math.inf):
            continue
        c1 = s1.ksc_m(m1s)[2]
        for st2 in stats2:
            m2s = st2["mass"]
            if not (0.0 < m2s < math.inf):
                continue
            c2 = s2.ksc_m(m2s)[2]
            if c1 > 0.0 and c2 > 0.0:
                w1, w2 = 1.0 / c1, 1.0 / c2
                dM = M - (m1s + m2s)
                m1_opt = m1s + dM * (w1 / (w1 + w2))
                if 0.0 < m1_opt < M:
                    pts.append(m1_opt)
                    for fac in (-0.5, 0.5, -2.0, 2.0):
                        cm = m1_opt + fac * max(abs(dM), 1e-9 * M)
                        if 0.0 < cm < M:
                            pts.append(cm)
    if not pt1["convention"] and not pt2["convention"]:
        m1s, m2s = pt1["mass"], pt2["mass"]
        if 0.0 < m1s < math.inf and 0.0 < m2s < math.inf:
            pts.append(M * m1s / (m1s + m2s))
    pts = sorted(set(pts))
    vals = [_lump_p_split(s1, s2, x, M - x) for x in pts]
    best = max(vals)
    if best >= 0.999999999:
        return best
    idx = max(range(len(pts)), key=lambda i: vals[i])
    lo = pts[idx - 1] if idx > 0 else 1e-12 * M
    hi = pts[idx + 1] if idx + 1 < len(pts) else (1.0 - 1e-12) * M
    if best > 1e-15 and hi > lo:
        neg_best, _ = _golden_min(lambda x: -_lump_p_split(s1, s2, x, M - x), lo, hi, steps=32)
        best = max(best, -neg_best)
    return best


def lump(first, second, *, level=0.95, grid=161, n_splits=31):
    """L5 two-body lumping M = m1 + m2 purely from Door A1 states (K, C).

    Zero Monte Carlo, exact coverage >= `level` uniformly everywhere, ghost-free on nonlinear
    asymmetric curves, collapses to 1-dof when one body is 0/0. Known cost: 2-dof width tax
    (+24.9% at strong signal when both bodies are informative; G27).
    """
    if not 0.0 < level < 1.0:
        raise ValueError("level must lie strictly in (0, 1)")
    s1 = _extract_state(first)
    s2 = _extract_state(second)
    pt1 = s1.argmin_K()
    pt2 = s2.argmin_K()
    stats1 = s1.stationary_points()
    stats2 = s2.stationary_points()
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

    def p_value(M):
        return _lump_p_total(s1, s2, M, stats1, stats2, pt1, pt2, n_splits=n_splits)

    base = s1._half_mesh(grid=grid, log_decades=10, pts_per_decade=10)
    pts_half = {0: set(base), 1: set(base)}
    anchors = []
    if not point_is_convention and 0.0 < point < math.inf:
        anchors.append(point)
    for st in stats1 + stats2:
        if 0.0 < st["mass"] < math.inf:
            anchors.append(st["mass"])
    for st1 in stats1:
        for st2 in stats2:
            s_m = st1["mass"] + st2["mass"]
            if 0.0 < s_m < math.inf:
                anchors.append(s_m)
    for anc in anchors:
        h, x0 = _mass_to_half(anc, scale)
        pts_half[h].add(x0)
        for exp in (-6, -4, -3, -2, -1):
            step = 10.0 ** exp
            for sgn in (-3, -1, 1, 3):
                mc = anc * (1.0 + sgn * step)
                if mc > 0.0:
                    hc, xc = _mass_to_half(mc, scale)
                    pts_half[hc].add(xc)

    alpha = 1.0 - level

    def is_kept_lump(h, x):
        # Compatible set lives on the positive ray (0, inf) (L1): evaluate boundary limits M -> 0+ and M -> inf
        # via x = 1e-30 so isolated C=0 boundary pinholes (e.g. m=0 on (10, 0)) do not create fake slivers.
        m_eval = _half_to_mass(h, 1e-30 if x == 0.0 else x, scale)
        return p_value(m_eval) > alpha

    segs = []
    for h in (0, 1):
        xs = sorted(pts_half[h])
        ks = [is_kept_lump(h, x) for x in xs]
        if h == 1:
            xs = list(reversed(xs))
            ks = list(reversed(ks))
        in_seg = ks[0]
        start_m = _half_to_mass(h, xs[0], scale) if in_seg else None
        for i in range(len(xs) - 1):
            if ks[i] != ks[i + 1]:
                if h == 0:
                    xr = _bisect_bool_boundary(lambda u: is_kept_lump(0, u), xs[i], xs[i + 1], ks[i], steps=48)
                else:
                    xr = _bisect_bool_boundary(lambda u: is_kept_lump(1, u), xs[i + 1], xs[i], ks[i + 1], steps=48)
                mr = _half_to_mass(h, xr, scale)
                if not in_seg:
                    start_m = mr
                    in_seg = True
                else:
                    segs.append((start_m, mr))
                    in_seg = False
                    start_m = None
        if in_seg:
            segs.append((start_m, _half_to_mass(h, xs[-1], scale)))

    merged = []
    for lo, hi in segs:
        if merged and math.isfinite(merged[-1][1]) and math.isfinite(lo) and lo <= merged[-1][1] + 1e-10 * max(1e-300, lo, merged[-1][1]):
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    kept = tuple(merged)

    return dict(
        point=point,
        point_status=point_status,
        point_is_convention=point_is_convention,
        set=kept,
        shape=classify_shape(kept),
        leftover=pt1["K_min"] + pt2["K_min"],
        curve=lambda M: lump_curve(s1, s2, M),
        p_value=p_value,
        level=float(level),
    )


def _swap_state(st):
    """Exact T5 channel-swap on a State (m <-> 1/m, scale <-> 1/scale)."""
    # Build swapped spans directly without storing raw readings!
    # Wait: swapping a Reading swaps (bf <-> ba, gff <-> gaa, gfa <-> gfa^T, uf <-> ua, s_i <-> 1/s_i).
    sw_spans = []
    for sp in st._spans:
        nsp = object.__new__(_SpanPre)
        nsp.d = sp.d
        nsp.k = sp.k
        nsp.is_exact_zero = False
        nsp.s_i = 1.0 / sp.s_i
        nsp.zW = sp.zW
        nsp.energy = sp.energy
        nsp.bf = sp.ba
        nsp.ba = sp.bf
        nsp.gff = sp.gaa
        nsp.gaa = sp.gff
        k = sp.k
        nsp.gfa = tuple(tuple(sp.gfa[j][i] for j in range(k)) for i in range(k))
        nsp.gsym = tuple(tuple(nsp.gfa[i][j] + nsp.gfa[j][i] for j in range(k)) for i in range(k))
        nsp.uf = sp.ua
        nsp.ua = sp.uf
        sw_spans.append(nsp)
    return State(
        scale=1.0 / st.scale,
        _spans=sw_spans,
        _all_scales=tuple(1.0 / s for s in st._all_scales),
    )


def common_force(first, second, *, level=0.95, grid=161, n_splits=31):
    """L14 common-force composition 1/M = 1/m1 + 1/m2 via exact T5 channel-swap duality."""
    s1 = _extract_state(first)
    s2 = _extract_state(second)
    sw = lump(_swap_state(s1), _swap_state(s2), level=level, grid=grid, n_splits=n_splits)
    inv_set = tuple(sorted((_recip(hi), _recip(lo)) for lo, hi in sw["set"]))
    point = _recip(sw["point"]) if not sw["point_is_convention"] else 1.0 / (1.0 / s1.scale + 1.0 / s2.scale)
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
        leftover=sw["leftover"],
        curve=lambda M: sw["curve"](_recip(M)),
        p_value=lambda M: sw["p_value"](_recip(M)),
        level=float(level),
    )


def lump_known(body, m_known, *, level=0.95, grid=241):
    """Exact 1-dof L5 composition of a measured body with an exact known mass m_known > 0: M = m + m_known."""
    m_known = float(m_known)
    if not (m_known > 0.0 and math.isfinite(m_known)):
        raise ValueError("m_known must be positive and finite")
    s = _extract_state(body)
    pt = s.argmin_K(grid=grid)
    kept_raw = s.kept_intervals(level=level, grid=grid)
    shifted = tuple((lo + m_known, hi + m_known) for lo, hi in kept_raw)
    point = (s.scale + m_known) if pt["convention"] else (pt["mass"] + m_known)
    return dict(
        point=point,
        point_status=("convention" if pt["convention"] else ("boundary_inf" if math.isinf(point) else "interior")),
        point_is_convention=pt["convention"],
        set=shifted,
        shape=classify_shape(shifted),
        leftover=pt["K_min"],
        curve=lambda M: math.inf if M < m_known else s.curve(M - m_known),
        p_value=lambda M: 0.0 if M < m_known else s.p_value(M - m_known),
        level=float(level),
    )


def common_force_known(body, m_known, *, level=0.95, grid=241):
    """Exact 1-dof L14 composition with an exact known mass m_known > 0: 1/M = 1/m + 1/m_known."""
    m_known = float(m_known)
    if not (m_known > 0.0 and math.isfinite(m_known)):
        raise ValueError("m_known must be positive and finite")
    s = _extract_state(body)
    sw = lump_known(_swap_state(s), 1.0 / m_known, level=level, grid=grid)
    inv_set = tuple(sorted((_recip(hi), _recip(lo)) for lo, hi in sw["set"]))
    point = _recip(sw["point"]) if not sw["point_is_convention"] else 1.0 / (1.0 / s.scale + 1.0 / m_known)
    st = sw["point_status"]
    if st == "boundary_inf":
        st = "boundary_zero"
    return dict(
        point=point,
        point_status=st,
        point_is_convention=sw["point_is_convention"],
        set=inv_set,
        shape=classify_shape(inv_set),
        leftover=sw["leftover"],
        curve=lambda M: sw["curve"](_recip(M)),
        p_value=lambda M: sw["p_value"](_recip(M)),
        level=float(level),
    )
