"""E14 prototype & verification library for the three walls of ours34 (pure Python stdlib).

Exploration branch `explore/e14-ours34-walls` (from `claude/intake-2026-11-22`, 22 Nov 2026).
Zero external dependencies so every check runs under plain `python3 -m unittest`.

Implements on the full reading type (arbitrary dimension d, arbitrary SPD 2d x 2d covariance
Sigma_i including cross-channel correlation, arbitrary k-dimensional orthonormal subspace N_i):
  - Exact T36/T45 per-reading geometry: K_i(m), S_i(m) = -0.5 * dK_i/dm, C_i(m) = Var(S_i | U_i)
  - Minimal additive state `EvidenceState(K, C, M_nz)` (Wall 3 & Wall 2B)
  - Exact tropical point `argmin_K` on the closed ray [0, inf] (Wall 1, Theorems E14-1 & E14-2)
  - Exact T45 Gaussian score test, exact closed-form joint-spherical Student-t test `t_joint`,
    and their union `floor_score` (Wall 2B, Theorem E14-7)
  - Exact curve-pair two-body composition `lump_qpos` (L5) and `common_force_qpos` (L14)
    using ONLY (K_1, C_1) and (K_2, C_2), with zero Monte Carlo (Wall 3, Theorem E14-4)
  - Exact tropical infimal convolutions (`lump_curve`, `common_force_curve`, `iterated_pinch`)
  - T10 noise-angle circular mean reference (Wall 1 / B28 item 2)
"""
from dataclasses import dataclass
import math


# =============================================================================
# Special functions & deterministic quadrature (zero external dependencies)
# =============================================================================

def norm_cdf(x):
    if math.isinf(x):
        return 1.0 if x > 0.0 else 0.0
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def norm_pdf(x):
    if math.isinf(x):
        return 0.0
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


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
    """Regularized incomplete beta I_x(a, b)."""
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
# Pure-Python linear algebra for full correlated covariance & subspace N_i
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
                    raise ValueError("covariance matrix must be positive definite")
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


def _solve_upper(R, b):
    n = len(R)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - math.fsum(R[i][j] * x[j] for j in range(i + 1, n))) / R[i][i]
    return x


def _qr_complete(cols):
    """Modified Gram-Schmidt (with reorthogonalization) + basis completion in R^m."""
    m = len(cols[0])
    k = len(cols)
    Q = []
    R = [[0.0] * k for _ in range(k)]
    for j in range(k):
        v = list(cols[j])
        for _ in range(2):
            for i in range(j):
                r = _dot(Q[i], v)
                R[i][j] += r
                v = [v[t] - r * Q[i][t] for t in range(m)]
        nrm = math.sqrt(_dot(v, v))
        if nrm <= 0.0:
            raise ValueError("linearly dependent columns in span basis")
        R[j][j] = nrm
        Q.append([x / nrm for x in v])
    Q_all = list(Q)
    for e_idx in range(m):
        if len(Q_all) == m:
            break
        v = [1.0 if t == e_idx else 0.0 for t in range(m)]
        for _ in range(2):
            for q in Q_all:
                r = _dot(q, v)
                v = [v[t] - r * q[t] for t in range(m)]
        nrm = math.sqrt(_dot(v, v))
        if nrm > 1e-8:
            Q_all.append([x / nrm for x in v])
    return Q_all[:k], Q_all[k:], R


# =============================================================================
# Full-type Reading & Whitened Span (T36 / T73 / T45)
# =============================================================================

@dataclass(frozen=True)
class Reading:
    """Full-type reading: force and acceleration in R^d, SPD 2d x 2d covariance,
    optional k orthonormal subspace columns in R^d (default full I_d)."""
    force: tuple
    acceleration: tuple
    cov: tuple
    subspace: tuple = None

    @staticmethod
    def scalar(F, a, sf, sa, rho=0.0):
        cov = ((sf * sf, rho * sf * sa), (rho * sf * sa, sa * sa))
        return Reading((float(F),), (float(a),), cov)

    @staticmethod
    def isotropic(F_vec, a_vec, sf, sa, subspace=None):
        d = len(F_vec)
        cov = tuple(
            tuple((sf * sf if i < d else sa * sa) if i == j else 0.0 for j in range(2 * d))
            for i in range(2 * d)
        )
        sub = None if subspace is None else tuple(tuple(float(x) for x in col) for col in subspace)
        return Reading(tuple(float(x) for x in F_vec), tuple(float(x) for x in a_vec), cov, sub)

    @property
    def d(self):
        return len(self.force)

    @property
    def basis_cols(self):
        if self.subspace is not None:
            return self.subspace
        d = self.d
        return tuple(tuple(1.0 if t == j else 0.0 for t in range(d)) for j in range(d))

    @property
    def k(self):
        return len(self.basis_cols)

    @property
    def is_exact_zero(self):
        return all(x == 0.0 for x in self.force) and all(x == 0.0 for x in self.acceleration)

    def swapped(self):
        """T5 channel swap (F <-> a) with matching covariance block swap."""
        d = self.d
        perm = list(range(d, 2 * d)) + list(range(d))
        new_cov = tuple(tuple(self.cov[i][j] for j in perm) for i in perm)
        return Reading(self.acceleration, self.force, new_cov, self.subspace)

    def scaled(self, lam):
        """L3 mass-unit scaling m -> lam * m (F -> lam * F, sf -> lam * sf)."""
        d = self.d
        s = [lam] * d + [1.0] * d
        new_cov = tuple(tuple(self.cov[i][j] * s[i] * s[j] for j in range(2 * d)) for i in range(2 * d))
        return Reading(tuple(lam * x for x in self.force), self.acceleration, new_cov, self.subspace)


class WhitenedSpan:
    """Precomputed whitened 2k-dimensional span W for one Reading (T36)."""

    def __init__(self, reading):
        self.reading = reading
        self.d = reading.d
        self.k = reading.k
        self.is_exact_zero = reading.is_exact_zero
        L = _cholesky(reading.cov)
        y = list(reading.force) + list(reading.acceleration)
        z = [0.0] * (2 * self.d) if self.is_exact_zero else _solve_lower(L, y)
        span_cols = []
        for col in reading.basis_cols:
            span_cols.append(_solve_lower(L, list(col) + [0.0] * self.d))
        for col in reading.basis_cols:
            span_cols.append(_solve_lower(L, [0.0] * self.d + list(col)))
        W_basis, _, _ = _qr_complete(span_cols)
        self.zW = [0.0] * (2 * self.k) if self.is_exact_zero else [_dot(w, z) for w in W_basis]
        self.bf = [[_dot(w, span_cols[j]) for w in W_basis] for j in range(self.k)]
        self.ba = [[_dot(w, span_cols[self.k + j]) for w in W_basis] for j in range(self.k)]
        self.energy = _dot(self.zW, self.zW)

    def eval_ksc(self, m):
        """Return (K_i(m), S_i(m), C_i(m), H_i(m)) at m in [0, inf].
        S_i = -0.5 * dK_i/dm, C_i = Var(S_i | U_i), H_i = ||U_i||^2."""
        if self.is_exact_zero:
            return 0.0, 0.0, 0.0, 0.0
        k = self.k
        if math.isinf(m):
            b_cols = self.bf
            Q_par, Q_perp, _ = _qr_complete(b_cols)
            U = [_dot(q, self.zW) for q in Q_par]
            V = [_dot(q, self.zW) for q in Q_perp]
            return _dot(V, V), 0.0, 0.0, _dot(U, U)
        b_cols = [[m * self.bf[j][t] + self.ba[j][t] for t in range(2 * k)] for j in range(k)]
        Q_par, Q_perp, R = _qr_complete(b_cols)
        U = [_dot(q, self.zW) for q in Q_par]
        V = [_dot(q, self.zW) for q in Q_perp]
        coeff = _solve_upper(R, U)
        db_coeff = [math.fsum(self.bf[j][t] * coeff[j] for j in range(k)) for t in range(2 * k)]
        tang = [_dot(q, db_coeff) for q in Q_perp]
        return _dot(V, V), _dot(V, tang), _dot(tang, tang), _dot(U, U)


# =============================================================================
# Minimal Additive Evidence State: (K(m), C(m), M_nz)
# =============================================================================

class EvidenceState:
    """Carries the minimal additive state (K(m), C(m), M_nz) of a series of readings.
    Independent blocks of the same body combine by exact addition: state1 + state2."""

    def __init__(self, readings=(), scale=None, _spans=None):
        self.readings = tuple(readings)
        self.spans = tuple(WhitenedSpan(r) for r in self.readings) if _spans is None else tuple(_spans)
        if scale is None:
            src = [r for r in self.readings if not r.is_exact_zero] or list(self.readings)
            logs = []
            for r in src:
                d = r.d
                tr_f = math.fsum(r.cov[i][i] for i in range(d))
                tr_a = math.fsum(r.cov[d + i][d + i] for i in range(d))
                logs.append(0.5 * math.log(tr_f / tr_a))
            self.scale = math.exp(math.fsum(logs) / len(logs)) if logs else 1.0
        else:
            self.scale = float(scale)
        self.M_total = sum(s.k for s in self.spans)
        self.M_nz = sum(s.k for s in self.spans if not s.is_exact_zero)

    def __add__(self, other):
        if not isinstance(other, EvidenceState):
            return NotImplemented
        return EvidenceState(self.readings + other.readings, _spans=self.spans + other.spans)

    def ksc(self, m):
        """Return (K(m), S_m(m), C_m(m)) at m in [0, inf]."""
        K = S = C = 0.0
        for sp in self.spans:
            ki, si, ci, _ = sp.eval_ksc(m)
            K += ki
            S += si
            C += ci
        return K, S, C

    def K(self, m):
        return self.ksc(m)[0]

    def at_angle(self, u):
        m = math.inf if u >= 0.5 * math.pi else (0.0 if u <= 0.0 else self.scale * math.tan(u))
        return self.K(m)

    def argmin_K(self, grid=721, tol=1e-12):
        """Find the global minimum of K on the closed ray [0, inf] (angle u in [0, pi/2]).
        Detects flat curves, boundary minima (0 or inf), and ties (multiple global minima)."""
        us = [0.5 * math.pi * i / (grid - 1) for i in range(grid)]
        vals = [self.at_angle(u) for u in us]
        vmin = min(vals)
        vmax = max(vals)
        if vmax - vmin <= tol * max(1.0, vmax):
            return dict(
                mass=None, angle=None, K_min=vmin, status="flat",
                minima=(), K_0=vals[0], K_inf=vals[-1],
            )
        # Refine every local minimum candidate on the grid
        cand_indices = []
        for j in range(grid):
            left = vals[j - 1] if j > 0 else math.inf
            right = vals[j + 1] if j < grid - 1 else math.inf
            if vals[j] <= left and vals[j] <= right and vals[j] <= vmin + 1e-3 * max(1.0, vmax - vmin):
                cand_indices.append(j)
        refined = []
        for j in cand_indices:
            a = us[max(0, j - 1)]
            b = us[min(grid - 1, j + 1)]
            fu, u_best = _golden_min(self.at_angle, a, b)
            if vals[j] < fu:
                fu, u_best = vals[j], us[j]
            # Polish interior minimum with exact S(m) = -0.5 * dK/dm = 0 bisection
            if 1e-12 < u_best < 0.5 * math.pi - 1e-12:
                lo_u, hi_u = max(a, u_best - 1e-4), min(b, u_best + 1e-4)
                s_lo = self.ksc(self.scale * math.tan(lo_u))[1]
                s_hi = self.ksc(self.scale * math.tan(hi_u))[1]
                if s_lo > 0.0 and s_hi < 0.0:
                    for _ in range(60):
                        mid_u = 0.5 * (lo_u + hi_u)
                        if self.ksc(self.scale * math.tan(mid_u))[1] > 0.0:
                            lo_u = mid_u
                        else:
                            hi_u = mid_u
                    u_pol = 0.5 * (lo_u + hi_u)
                    f_pol = self.at_angle(u_pol)
                    if f_pol <= fu + 1e-12 * max(1.0, abs(fu)):
                        fu, u_best = min(f_pol, fu), u_pol
            refined.append((fu, u_best))
        # Also explicitly check exact endpoints 0 and pi/2
        refined.append((vals[0], 0.0))
        refined.append((vals[-1], 0.5 * math.pi))
        best_val = min(x[0] for x in refined)
        # Deduplicate distinct global minima
        global_us = []
        for fu, u in sorted(refined, key=lambda x: x[1]):
            if fu <= best_val + 1e-9 * max(1.0, vmax):
                if not global_us or abs(u - global_us[-1]) > 1e-5:
                    global_us.append(u)
        masses = tuple(
            math.inf if u >= 0.5 * math.pi - 1e-12 else (0.0 if u <= 1e-12 else self.scale * math.tan(u))
            for u in global_us
        )
        if len(masses) > 1:
            status = "tied"
            mass_out = None
            ang_out = None
        else:
            mass_out = masses[0]
            ang_out = global_us[0]
            if mass_out == 0.0:
                status = "boundary_zero"
            elif math.isinf(mass_out):
                status = "boundary_inf"
            else:
                status = "interior"
        return dict(
            mass=mass_out, angle=ang_out, K_min=best_val, status=status,
            minima=masses, K_0=vals[0], K_inf=vals[-1],
        )

    def t45_test(self, m):
        """Exact T45 conditional Gaussian score test using ONLY (K, C)."""
        K, S, C = self.ksc(m)
        if C <= 0.0:
            return dict(z=0.0, statistic=0.0, p_value=1.0, K=K, S=S, C=C, degenerate=True)
        z = S / math.sqrt(C)
        return dict(z=z, statistic=z * z, p_value=math.erfc(abs(z) / math.sqrt(2.0)),
                    K=K, S=S, C=C, degenerate=False)

    def t_joint_test(self, m):
        """Exact closed-form joint-spherical Student-t test on (K, C, M_nz) (Theorem E14-7).
        Exact under ANY common unknown error scale v > 0 on the full reading type."""
        K, S, C = self.ksc(m)
        df = self.M_nz - 1
        if C <= 0.0 or df <= 0:
            return dict(t=0.0, df=max(0, df), p_value=1.0, v_hat=None, degenerate=True)
        z2 = (S * S) / C
        perp_sq = max(0.0, K - z2)
        v_hat = perp_sq / df
        if perp_sq <= 1e-30 * max(1.0, K):
            if abs(S) <= 1e-15 * math.sqrt(C):
                return dict(t=0.0, df=df, p_value=1.0, v_hat=0.0, degenerate=False)
            return dict(t=math.inf, df=df, p_value=0.0, v_hat=0.0, degenerate=False)
        t_val = abs(S) / math.sqrt(C * v_hat)
        return dict(t=t_val, df=df, p_value=t_sf_two_sided(t_val, df), v_hat=v_hat, degenerate=False)

    def floor_test(self, m, *, single_scalar_policy="t45"):
        """Union ("floor") of T45 and t_joint: keeps m when EITHER keeps m.
        `single_scalar_policy`:
          - 't45': at M_nz == 1 (df == 0, no residual dof to estimate disagreement), fall back to T45;
          - 'vacuous': keep strict v-validity at M_nz == 1 (p_t = 1 => keeps all m)."""
        g = self.t45_test(m)
        t = self.t_joint_test(m)
        if self.M_nz <= 1 and single_scalar_policy == "t45":
            p_floor = g["p_value"]
        else:
            p_floor = max(g["p_value"], t["p_value"])
        return dict(p_value=p_floor, p_t45=g["p_value"], p_tjoint=t["p_value"],
                    v_hat=t["v_hat"], df=t["df"])

    def summary(self, level=0.95):
        """Full ours34-style prototype summary: tropical point, leftover, and kept set endpoints."""
        pt = self.argmin_K()
        return dict(
            point=pt["mass"],
            point_status=pt["status"],
            minima=pt["minima"],
            K_min=pt["K_min"],
            M_nz=self.M_nz,
            M_total=self.M_total,
            scale=self.scale,
        )


def _golden_min(fn, a, b, steps=52):
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
# Wall 1: Tropical infimal convolutions & iterated pinch Q_n = (K_+^(n))_par^(n)
# =============================================================================

def lump_curve(k1_fn, k2_fn, M, grid=241):
    """T74 infimal convolution (K1 square_+ K2)(M) = inf_{m1+m2=M} [K1(m1) + K2(m2)]."""
    if math.isinf(M):
        min1 = _min_on_ray(k1_fn)
        min2 = _min_on_ray(k2_fn)
        return min(k1_fn(math.inf) + min2, min1 + k2_fn(math.inf))
    if M == 0.0:
        return k1_fn(0.0) + k2_fn(0.0)
    def f(t):
        return k1_fn(t * M) + k2_fn((1.0 - t) * M)
    ts = [i / (grid - 1) for i in range(grid)]
    vs = [f(t) for t in ts]
    j = min(range(grid), key=lambda i: vs[i])
    lo, hi = ts[max(0, j - 1)], ts[min(grid - 1, j + 1)]
    best, _ = _golden_min(f, lo, hi, steps=40)
    return min(best, vs[j])


def common_force_curve(k1_fn, k2_fn, M, grid=241):
    """T74 parallel infimal convolution (K1 square_par K2)(M) via T5 swap."""
    inv = math.inf if M == 0.0 else (0.0 if math.isinf(M) else 1.0 / M)
    def s1(w):
        return k1_fn(math.inf if w == 0.0 else (0.0 if math.isinf(w) else 1.0 / w))
    def s2(w):
        return k2_fn(math.inf if w == 0.0 else (0.0 if math.isinf(w) else 1.0 / w))
    return lump_curve(s1, s2, inv, grid=grid)


def _min_on_ray(k_fn, scale=1.0, grid=241):
    us = [0.5 * math.pi * i / (grid - 1) for i in range(grid)]
    def at_u(u):
        return k_fn(math.inf if u >= 0.5 * math.pi else (0.0 if u <= 0.0 else scale * math.tan(u)))
    vs = [at_u(u) for u in us]
    j = min(range(grid), key=lambda i: vs[i])
    best, _ = _golden_min(at_u, us[max(0, j - 1)], us[min(grid - 1, j + 1)], steps=36)
    return min(best, vs[j])


def t10_circular_mean_discrete(atoms, s):
    """T10 zero-padding endpoint on a discrete proper law [(m_j, p_j), ...]:
    s * E[sin theta] / E[cos theta] with theta = atan(m / s)."""
    num = math.fsum(p * math.sin(math.atan(m / s)) for m, p in atoms)
    den = math.fsum(p * math.cos(math.atan(m / s)) for m, p in atoms)
    return s * num / den


# =============================================================================
# Wall 3: Exact Curve-Pair Two-Body Composition Q_pos (Zero Monte Carlo!)
# =============================================================================

def q_pos_statistic(z1, z2, a1, a2, m1, m2):
    """Positive-quadrant profiled quadratic score statistic:
    min_{-m1<=d<=m2} [(z1 - a1 d)^2 + (z2 + a2 d)^2] - (-a1 m1 - z1)_+^2 - (-a2 m2 - z2)_+^2."""
    L2 = a1 * a1 + a2 * a2
    if L2 <= 0.0:
        return 0.0
    d = max(-m1, min(m2, (a1 * z1 - a2 * z2) / L2))
    c_lump = (z1 - a1 * d) ** 2 + (z2 + a2 * d) ** 2
    c1 = max(0.0, -a1 * m1 - z1) ** 2
    c2 = max(0.0, -a2 * m2 - z2) ** 2
    return max(0.0, c_lump - c1 - c2)


def q_pos_sf_exact(q, a1, a2, m1, m2):
    """Exact deterministic survival probability P(Q_pos(Z1, Z2) >= q) for Z1, Z2 iid N(0, 1).
    Uses exact analytical integration in Z2 and kink-split Gauss-Legendre in Z1."""
    if q <= 0.0:
        return 1.0
    r = math.sqrt(q)
    A1, A2, L2 = a1 * m1, a2 * m2, a1 * a1 + a2 * a2

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


def lump_qpos_split(state1, state2, M, m1):
    """Exact conditional p-value at split (m1, M - m1) using ONLY (K1, C1) and (K2, C2).
    Zero raw readings, zero counts, zero Monte Carlo."""
    m2 = M - m1
    _, S1, C1 = state1.ksc(m1)
    _, S2, C2 = state2.ksc(m2)
    if C1 <= 0.0 and C2 <= 0.0:
        return 1.0
    if C2 <= 0.0:
        a1 = math.sqrt(C1)
        Z1 = S1 / a1
        q = max(0.0, Z1 - a1 * m2) ** 2
        return 1.0 if q == 0.0 else 1.0 - norm_cdf(a1 * m2 + math.sqrt(q))
    if C1 <= 0.0:
        a2 = math.sqrt(C2)
        Z2 = S2 / a2
        q = max(0.0, Z2 - a2 * m1) ** 2
        return 1.0 if q == 0.0 else 1.0 - norm_cdf(a2 * m1 + math.sqrt(q))
    a1, a2 = math.sqrt(C1), math.sqrt(C2)
    Z1, Z2 = S1 / a1, S2 / a2
    q = q_pos_statistic(Z1, Z2, a1, a2, m1, m2)
    return q_pos_sf_exact(q, a1, a2, m1, m2)


def lump_qpos_total(state1, state2, M, n_splits=41):
    """Valid p-value for total mass M = m1 + m2 (L5) using ONLY (K1, C1) and (K2, C2)."""
    shares = [1e-6, 1e-4, 1e-2] + [i / n_splits for i in range(1, n_splits)] + [1.0 - 1e-2, 1.0 - 1e-4, 1.0 - 1e-6]
    return max(lump_qpos_split(state1, state2, M, t * M) for t in shares)


def common_force_qpos_total(state1, state2, M, n_splits=41):
    """Valid p-value for reduced mass 1/M = 1/m1 + 1/m2 (L14) via exact T5 swap duality."""
    s1_sw = EvidenceState(tuple(r.swapped() for r in state1.readings))
    s2_sw = EvidenceState(tuple(r.swapped() for r in state2.readings))
    return lump_qpos_total(s1_sw, s2_sw, 1.0 / M, n_splits=n_splits)
