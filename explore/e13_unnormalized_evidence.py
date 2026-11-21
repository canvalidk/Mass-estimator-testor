"""E13: Unnormalized evidence ('bell with sigma = inf') prototype and killer gallery.

Pure Python (standard library only) so the prototype and its checks run cleanly without
external packages. Nothing in met/ or definition/ is modified.

What this module computes for a series of full-covariance, arbitrary-subspace readings:
  1. T36 residual cost K(m) = sum_i K_i(m) and unnormalized evidence L(m) = exp(-K(m)/2),
     with no prior pi baked in. On all-zero 0/0, K(m) == 0 (L(m) == 1 everywhere: a bell
     with sigma = inf).
  2. Door A (strict Gaussian consequences):
       - T45 exact conditional score 95% compatibility set (passes exact-zero Z; can have
         hilltop ghosts on anti-parallel data).
       - Moreira (2003) exact conditional likelihood-ratio (CLR) 95% compatibility set on
         the common-ratio isotropic slice (removes hilltop ghosts; uses noise-dimension
         count M = sum k_i).
       - Separate leftover surprise readout K_min = min_{m>0} K(m) and upper-bound p-value
         P(chi^2_M >= K_min).
  3. Door B (declared heavy-tailed / fluke-allowing error model, unnormalized):
       - T48 Student-t residual cost sum_i (nu + k_i) * log(1 + K_i(m) / nu), with NO prior
         pi and NO series power/warp eta.
  4. Prior-free composition of unnormalized residual curves:
       - Same object: K_total(m) = K_1(m) + K_2(m) (so 0/0 adds 0 — exact identity).
       - Two-body lumping (L5, M = m1 + m2): infimal convolution
           K_lump(M) = inf_{m1 + m2 = M} [K_1(m1) + K_2(m2)].
       - Two-body common force (L14, 1/M = 1/m1 + 1/m2): infimal convolution
           K_reduced(M) = inf_{1/m1 + 1/m2 = 1/M} [K_1(m1) + K_2(m2)].

Run the killer gallery from the repo root:
    python3 explore/e13_unnormalized_evidence.py
"""
from dataclasses import dataclass
import math
import random

Z95 = 1.959963984540054
CHI2_95_1 = Z95 * Z95  # 3.841458820694124


# ---------------------------------------------------------------------------
# Minimal exact linear algebra (matrices as lists of rows)
# ---------------------------------------------------------------------------

def mat_t(A):
    return [list(col) for col in zip(*A)]


def mat_vec(A, x):
    return [sum(a * b for a, b in zip(row, x)) for row in A]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def cholesky(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                val = A[i][i] - s
                if val <= 0.0:
                    raise ValueError("Covariance matrix must be positive definite")
                L[i][j] = math.sqrt(val)
            else:
                L[i][j] = (A[i][j] - s) / L[j][j]
    return L


def solve_lower(L, b):
    n = len(L)
    x = [0.0] * n
    for i in range(n):
        s = sum(L[i][j] * x[j] for j in range(i))
        x[i] = (b[i] - s) / L[i][i]
    return x


def solve_upper(R, b):
    n = len(R)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(R[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (b[i] - s) / R[i][i]
    return x


def solve_spd(A, b):
    L = cholesky(A)
    y = solve_lower(L, b)
    return solve_upper(mat_t(L), y)


def solve_lower_mat(L, B):
    cols = [solve_lower(L, col) for col in mat_t(B)]
    return mat_t(cols)


def householder_qr_full(B):
    """Full QR of m x k matrix B (m >= k): returns Q (m x m), R_top (k x k)."""
    m = len(B)
    k = len(B[0])
    R = [row[:] for row in B]
    Q = [[1.0 if i == j else 0.0 for j in range(m)] for i in range(m)]
    for j in range(k):
        norm_x = math.sqrt(sum(R[i][j] ** 2 for i in range(j, m)))
        if norm_x == 0.0:
            continue
        sign = 1.0 if R[j][j] >= 0.0 else -1.0
        v = [0.0] * (m - j)
        v[0] = R[j][j] + sign * norm_x
        for i in range(1, m - j):
            v[i] = R[j + i][j]
        v_norm2 = sum(val * val for val in v)
        if v_norm2 == 0.0:
            continue
        beta = 2.0 / v_norm2
        for col in range(j, k):
            s = sum(v[i] * R[j + i][col] for i in range(m - j))
            for i in range(m - j):
                R[j + i][col] -= beta * v[i] * s
        for row in range(m):
            s = sum(Q[row][j + i] * v[i] for i in range(m - j))
            for i in range(m - j):
                Q[row][j + i] -= beta * s * v[i]
    R_top = [R[i][:k] for i in range(k)]
    return Q, R_top


# ---------------------------------------------------------------------------
# Special functions: chi-square survival & exact Moreira (2003) CLR p-value
# ---------------------------------------------------------------------------

def chi2_sf(x, df):
    """Exact upper-tail survival probability P(chi^2_df >= x)."""
    if x <= 0.0:
        return 1.0
    if df == 1:
        return math.erfc(math.sqrt(0.5 * x))
    if df == 2:
        return math.exp(-0.5 * x)
    a = 0.5 * df
    z = 0.5 * x
    if z < a + 1.0:
        ap = a
        term = 1.0 / a
        total = term
        for _ in range(500):
            ap += 1.0
            term *= z / ap
            total += term
            if abs(term) <= abs(total) * 1e-15:
                break
        p_low = total * math.exp(-z + a * math.log(z) - math.lgamma(a))
        return max(0.0, min(1.0, 1.0 - p_low))
    b = z + 1.0 - a
    c = 1.0 / 1e-300
    d = 1.0 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    q_val = math.exp(-z + a * math.log(z) - math.lgamma(a)) * h
    return max(0.0, min(1.0, q_val))


def moreira_clr_pvalue(lr, h, df_total):
    """Exact conditional p-value P(LR >= lr | H = h) for M = df_total noise dimensions.

    Uses Moreira's (2003) identity: LR <= c  <=>  T^2 <= c * (1 - W / (h + c)),
    where T ~ N(0, 1) and W ~ chi^2_{M-1} are independent.
    """
    if lr <= 0.0:
        return 1.0
    if df_total <= 1:
        return chi2_sf(lr, 1)
    if h <= 0.0:
        return chi2_sf(lr, df_total)
    k = df_total - 1
    lim = h + lr
    w_max = min(lim, k + 20.0 * math.sqrt(2.0 * k) + 80.0)
    half_k = 0.5 * k
    log_norm = half_k * math.log(2.0) + math.lgamma(half_k)
    sqrt_half_lr = math.sqrt(0.5 * lr)
    n_pts = 240
    dt = (0.5 * math.pi) / n_pts
    total = 0.0
    for j in range(n_pts):
        t = (j + 0.5) * dt
        s = math.sin(t)
        co = math.cos(t)
        w = w_max * s * s
        if w <= 0.0:
            continue
        log_fw = (half_k - 1.0) * math.log(w) - 0.5 * w - log_norm
        if log_fw < -700.0:
            continue
        ratio = max(0.0, 1.0 - w / lim)
        jac = 2.0 * w_max * s * co
        total += math.exp(log_fw) * math.erf(sqrt_half_lr * math.sqrt(ratio)) * jac * dt
    return max(0.0, min(1.0, 1.0 - total))


# ---------------------------------------------------------------------------
# Reading & core T36 / T45 / CLR / Student-t evaluations
# ---------------------------------------------------------------------------

@dataclass
class Reading:
    force: list
    acceleration: list
    cov: list
    subspace: list = None  # d x k matrix (rows = d, cols = k)

    def __post_init__(self):
        self.force = [float(x) for x in self.force]
        self.acceleration = [float(x) for x in self.acceleration]
        self.d = len(self.force)
        if len(self.acceleration) != self.d or self.d not in (1, 2, 3):
            raise ValueError("force and acceleration must have equal length 1..3")
        if self.subspace is None:
            self.subspace = [[1.0 if i == j else 0.0 for j in range(self.d)] for i in range(self.d)]
        self.k = len(self.subspace[0])
        self.chol = cholesky(self.cov)
        self.data = self.force + self.acceleration
        self.whitened = solve_lower(self.chol, self.data)
        z = [[0.0] * self.k for _ in range(self.d)]
        span = [self.subspace[i] + z[i] for i in range(self.d)] + [z[i] + self.subspace[i] for i in range(self.d)]
        span_w = solve_lower_mat(self.chol, span)
        Q_span, _ = householder_qr_full(span_w)
        proj_span = [dot([Q_span[r][c] for r in range(2 * self.d)], self.whitened) for c in range(2 * self.k)]
        self.energy = sum(v * v for v in proj_span)


def isotropic(F, a, sigma_f, sigma_a):
    if isinstance(F, (int, float)):
        F = [float(F)]
        a = [float(a)]
    d = len(F)
    cov = [[0.0] * (2 * d) for _ in range(2 * d)]
    for i in range(d):
        cov[i][i] = float(sigma_f) ** 2
        cov[d + i][d + i] = float(sigma_a) ** 2
    return Reading(F, a, cov)


def eval_reading_at_phi(reading, phi):
    """Return (H, K, score, cond_var) for one reading at candidate angle phi (m = tan(phi))."""
    s, co = math.sin(phi), math.cos(phi)
    n = reading.subspace
    d, k = reading.d, reading.k
    basis = [[s * n[i][j] for j in range(k)] for i in range(d)] + [
        [co * n[i][j] for j in range(k)] for i in range(d)
    ]
    deriv = [[co * n[i][j] for j in range(k)] for i in range(d)] + [
        [-s * n[i][j] for j in range(k)] for i in range(d)
    ]
    bw = solve_lower_mat(reading.chol, basis)
    dw = solve_lower_mat(reading.chol, deriv)
    Q, R = householder_qr_full(bw)
    Qt = mat_t(Q)
    q_par = Qt[:k]
    q_perp = Qt[k:]
    u = [dot(row, reading.whitened) for row in q_par]
    prec_y = solve_spd(reading.cov, reading.data)
    nt = mat_t(n)
    normal_rhs = [s * dot(nt[j], prec_y[:d]) + co * dot(nt[j], prec_y[d:]) for j in range(k)]
    if all(abs(val) <= 1e-15 * max(1.0, sum(abs(x) for x in reading.data)) for val in normal_rhs):
        if all(val == 0.0 for val in normal_rhs) or all(x == 0.0 for x in reading.data):
            u = [0.0] * k
    h = sum(v * v for v in u)
    k_res = max(0.0, reading.energy - h)
    r_perp = [dot(row, reading.whitened) for row in q_perp]
    alpha_hat = solve_upper(R, u)
    dw_alpha = mat_vec(dw, alpha_hat)
    tangent = [dot(row, dw_alpha) for row in q_perp]
    score = dot(r_perp, tangent)
    var = dot(tangent, tangent)
    return h, k_res, score, var


def evaluate_series(readings, phi):
    """Aggregate T36 & T45 quantities at candidate angle phi (m = tan(phi))."""
    H = K = S = C = 0.0
    per_K = []
    for r in readings:
        h, k_res, s, c = eval_reading_at_phi(r, phi)
        H += h
        K += k_res
        S += s
        C += c
        per_K.append((k_res, r.k))
    T = 0.0 if C <= 1e-300 else S / math.sqrt(C)
    return dict(H=H, K=K, S=S, C=C, T=T, T2=T * T, per_K=per_K)


def format_intervals(phis, mask):
    """Format connected intervals in m = tan(phi) from a boolean mask on open (0, pi/2)."""
    n = len(phis)
    if not any(mask):
        return "EMPTY"
    if all(mask):
        return "(0, inf)"
    runs = []
    i = 0
    while i < n:
        if mask[i]:
            j = i
            while j + 1 < n and mask[j + 1]:
                j += 1
            lo = "0" if i == 0 else f"{math.tan(phis[i]):.4g}"
            hi = "inf" if j == n - 1 else f"{math.tan(phis[j]):.4g}"
            runs.append(f"[{lo}, {hi}]" if (i > 0 and j < n - 1) else f"({lo}, {hi})")
            i = j + 1
        else:
            i += 1
    return " U ".join(runs)


def _minimize_1d(fn, a, b, steps=50):
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc, fd = fn(c), fn(d)
    for _ in range(steps):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = fn(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = fn(d)
    mid = 0.5 * (a + b)
    fm = fn(mid)
    return min((fc, c), (fd, d), (fm, mid), key=lambda x: x[0])


def analyze_case(readings, *, n_grid=1500, common_ratio=None, nu_door_b=10.0, scale_hint=1.0):
    """Evaluate the unnormalized evidence architecture on one series (no prior pi, no eta)."""
    eps = 1e-6
    us = [eps + (0.5 * math.pi - 2 * eps) * j / (n_grid - 1) for j in range(n_grid)]
    phis_coarse = [math.atan(scale_hint * math.tan(u)) for u in us]
    df_total = sum(r.k for r in readings)

    evals_c = [evaluate_series(readings, phi) for phi in phis_coarse]
    Ks_c = [e["K"] for e in evals_c]
    best_i = Ks_c.index(min(Ks_c))
    lo_b = phis_coarse[max(0, best_i - 2)]
    hi_b = phis_coarse[min(len(phis_coarse) - 1, best_i + 2)]
    K_min, phi_best = _minimize_1d(lambda ph: evaluate_series(readings, ph)["K"], lo_b, hi_b)

    extra_phis = [phi_best]
    for half_w in (0.05, 0.005, 0.0005):
        for j in range(-80, 81):
            p = phi_best + half_w * j / 80.0
            if eps < p < 0.5 * math.pi - eps:
                extra_phis.append(p)
    phis = sorted(set(phis_coarse + extra_phis))
    evals = [evaluate_series(readings, phi) for phi in phis]
    Ks = [e["K"] for e in evals]
    K_min = min(
        K_min,
        min(Ks),
        evaluate_series(readings, 0.0)["K"],
        evaluate_series(readings, 0.5 * math.pi)["K"],
    )
    K_max = max(Ks)
    m_hat = None if (K_max - K_min) < 1e-12 else math.tan(phis[Ks.index(min(Ks))])
    surprise_p = chi2_sf(K_min, df_total) if df_total > 0 else 1.0

    t45_mask = [e["T2"] <= CHI2_95_1 for e in evals]
    wilks_mask = [(k - K_min) <= CHI2_95_1 for k in Ks]

    clr_mask = None
    if common_ratio is not None:
        full_us = [-0.5 * math.pi + math.pi * j / 720 for j in range(720)]
        full_phis = [
            math.atan(scale_hint * math.tan(u)) if abs(math.cos(u)) > 1e-12 else 0.5 * math.pi
            for u in full_us
        ]
        full_ks = [evaluate_series(readings, ph)["K"] for ph in full_phis]
        wj = full_ks.index(min(full_ks))
        w_lo = full_phis[max(0, wj - 2)]
        w_hi = full_phis[min(len(full_phis) - 1, wj + 2)]
        if w_lo > w_hi:
            w_lo, w_hi = w_hi, w_lo
        k_whole_min, _ = _minimize_1d(lambda ph: evaluate_series(readings, ph)["K"], w_lo, w_hi)
        k_whole_min = min(k_whole_min, K_min)
        clr_mask = [
            moreira_clr_pvalue(max(0.0, e["K"] - k_whole_min), e["H"], df_total) >= 0.05
            for e in evals
        ]

    t_costs = [
        sum((nu_door_b + k_i) * math.log1p(k_res / nu_door_b) for k_res, k_i in e["per_K"])
        for e in evals
    ]
    t_min = min(t_costs)
    t_max = max(t_costs)
    t_best = None if (t_max - t_min) < 1e-12 else math.tan(phis[t_costs.index(t_min)])
    door_b_mask = [(tc - t_min) <= CHI2_95_1 for tc in t_costs]

    return dict(
        m_hat=m_hat,
        K_min=K_min,
        df_total=df_total,
        surprise_p=surprise_p,
        t45_set=format_intervals(phis, t45_mask),
        wilks_set=format_intervals(phis, wilks_mask),
        clr_set=None if clr_mask is None else format_intervals(phis, clr_mask),
        door_b_hat=t_best,
        door_b_set=format_intervals(phis, door_b_mask),
        phis=phis,
        Ks=Ks,
    )


# ---------------------------------------------------------------------------
# Prior-free composition: same object (K1 + K2), lumping (L5), common force (L14)
# ---------------------------------------------------------------------------

def compose_lump(K1_fn, K2_fn, M_values):
    """L5 lumping: K_lump(M) = inf_{0 < m1 < M} [K1(m1) + K2(M - m1)]."""
    out = []
    for M in M_values:
        n_scan = 200
        best = float("inf")
        best_u = 0.5
        for j in range(1, n_scan):
            u = j / n_scan
            val = K1_fn(u * M) + K2_fn((1.0 - u) * M)
            if val < best:
                best = val
                best_u = u
        a = max(1e-7, best_u - 1.5 / n_scan)
        b = min(1.0 - 1e-7, best_u + 1.5 / n_scan)
        val_ref, _ = _minimize_1d(lambda u: K1_fn(u * M) + K2_fn((1.0 - u) * M), a, b, steps=30)
        out.append(min(best, val_ref))
    return out


def compose_common_force(K1_fn, K2_fn, M_values):
    """L14 common force (reduced mass): K_red(M) = inf_{1/m1 + 1/m2 = 1/M} [K1(m1) + K2(m2)].

    Equivalent to lumping in reciprocal mass w = 1/m: w1 + w2 = 1/M.
    """
    return compose_lump(
        lambda w: K1_fn(1.0 / w),
        lambda w: K2_fn(1.0 / w),
        [1.0 / M for M in M_values],
    )


# ---------------------------------------------------------------------------
# Killer gallery CLI
# ---------------------------------------------------------------------------

def run_killer_gallery(reps=500, seed=2026):
    scenes = [
        ("1. All-zero 0/0 (3D, sigma_F=2, sigma_a=1)",
         [isotropic([0, 0, 0], [0, 0, 0], 2.0, 1.0)], 2.0, 2.0),
        ("2. Good aligned reading (F=10+-1, a=5+-1 => m=2)",
         [isotropic(10.0, 5.0, 1.0, 1.0)], 1.0, 2.0),
        ("3a. Can's Scene A: near pair (5.5, 5.6) +- 0.1, a=1+-0.002",
         [isotropic(5.5, 1.0, 0.1, 0.002), isotropic(5.6, 1.0, 0.1, 0.002)], 50.0, 5.5),
        ("3b. Can's Scene B: far pair (5.0, 6.0) +- 0.1 (10 sigma)",
         [isotropic(5.0, 1.0, 0.1, 0.002), isotropic(6.0, 1.0, 0.1, 0.002)], 50.0, 5.5),
        ("3c. Can's Scene B extreme: (5.0, 6.0) +- 0.02 (50 sigma)",
         [isotropic(5.0, 1.0, 0.02, 0.0004), isotropic(6.0, 1.0, 0.02, 0.0004)], 50.0, 5.5),
        ("4. Right-angle Scene C (G29: F=1,a=0 and F=0,a=1, sigma=0.05)",
         [isotropic(1.0, 0.0, 0.05, 0.05), isotropic(0.0, 1.0, 0.05, 0.05)], 1.0, 1.0),
        ("5. Anti-parallel pair ((2,-1) and (1,-2), sigma=0.1)",
         [isotropic(2.0, -1.0, 0.1, 0.1), isotropic(1.0, -2.0, 0.1, 0.1)], 1.0, 1.0),
        ("6. E12 3-reading triple (F=0.2, 1.0, 1.4, a=1, sigma=0.01)",
         [isotropic(f, 1.0, 0.01, 0.01) for f in (0.2, 1.0, 1.4)], 1.0, 1.0),
        ("7. 3D perpendicular single reading (F=(1,0,0), a=(0,1,0), sigma=0.1)",
         [isotropic([1, 0, 0], [0, 1, 0], 0.1, 0.1)], 1.0, 1.0),
    ]

    print("=== DETERMINISTIC KILLER SCENES (zero prior pi, zero eta) ===")
    for title, rds, cr, sh in scenes:
        res = analyze_case(rds, common_ratio=cr, scale_hint=sh)
        m_str = "NONE (flat)" if res["m_hat"] is None else f"{res['m_hat']:.4f}"
        tb_str = "NONE (flat)" if res["door_b_hat"] is None else f"{res['door_b_hat']:.4f}"
        print(f"\n{title}")
        print(f"  Peak m_hat={m_str} | K_min={res['K_min']:.3f} (df={res['df_total']}, p<={res['surprise_p']:.3g})")
        print(f"  T45 95%       : {res['t45_set']}")
        print(f"  Moreira CLR   : {res['clr_set']}")
        print(f"  Door B (t_10) : {res['door_b_set']}  (peak {tb_str})")

    rng = random.Random(seed)

    def rand_unit_vec(d):
        v = [rng.gauss(0.0, 1.0) for _ in range(d)]
        nrm = math.sqrt(dot(v, v))
        return [x / nrm for x in v]

    def make_world(d, N, r_len, m_true):
        th = math.atan(m_true)
        s, co = math.sin(th), math.cos(th)
        rds = []
        for _ in range(N):
            u = rand_unit_vec(d)
            F = [r_len * s * x + rng.gauss(0.0, 1.0) for x in u]
            a = [r_len * co * x + rng.gauss(0.0, 1.0) for x in u]
            rds.append(isotropic(F, a, 1.0, 1.0))
        return rds

    def sim_coverage(name, gen_series, m_true):
        cov_t45 = cov_clr = cov_wilks = 0
        unb_t45 = unb_clr = 0
        df_tot = None
        phi0 = math.atan(m_true)
        for _ in range(reps):
            rds = gen_series()
            if df_tot is None:
                df_tot = sum(r.k for r in rds)
            e0 = evaluate_series(rds, phi0)
            if e0["T2"] <= CHI2_95_1:
                cov_t45 += 1
            s_ff = sum(dot(r.force, r.force) for r in rds)
            s_aa = sum(dot(r.acceleration, r.acceleration) for r in rds)
            s_fa = sum(dot(r.force, r.acceleration) for r in rds)
            tr = s_ff + s_aa
            disc = math.hypot(s_ff - s_aa, 2.0 * s_fa)
            k_whole_min = max(0.0, 0.5 * (tr - disc))
            k_pos_min = k_whole_min if s_fa >= 0.0 else min(s_ff, s_aa)
            if moreira_clr_pvalue(max(0.0, e0["K"] - k_whole_min), e0["H"], df_tot) >= 0.05:
                cov_clr += 1
            if (e0["K"] - k_pos_min) <= CHI2_95_1:
                cov_wilks += 1
            eb = evaluate_series(rds, 0.5 * math.pi - 1e-6)
            if eb["T2"] <= CHI2_95_1:
                unb_t45 += 1
            if moreira_clr_pvalue(max(0.0, eb["K"] - k_whole_min), eb["H"], df_tot) >= 0.05:
                unb_clr += 1
        print(f"\n  {name} (reps={reps}, seed={seed}, true m={m_true})")
        print(f"    T45 95% cov={cov_t45/reps:5.1%} (unb {unb_t45/reps:5.1%}) | "
              f"CLR 95% cov={cov_clr/reps:5.1%} (unb {unb_clr/reps:5.1%}) | "
              f"Raw LR={cov_wilks/reps:5.1%}")

    print("\n=== MONTE CARLO COVERAGE ON GRAVEYARD WORLDS ===")
    sim_coverage("Pure noise (3D, N=20, signal=0, m=20)", lambda: make_world(3, 20, 0.0, 20.0), 20.0)
    sim_coverage("Weak signal (3D, N=20, signal=1, m=2)", lambda: make_world(3, 20, 1.0, 2.0), 2.0)

    def make_t43():
        th = math.atan(2.0)
        s, co = math.sin(th), math.cos(th)
        rds = [isotropic(20.0 * s + rng.gauss(0.0, 1.0), 20.0 * co + rng.gauss(0.0, 1.0), 1.0, 1.0)]
        for _ in range(99):
            rds.append(isotropic(rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0), 1.0, 1.0))
        return rds

    sim_coverage("T43 leverage (1D, 1x r=20 + 99 nulls, m=2)", make_t43, 2.0)
    sim_coverage("Strong signal (3D, N=10, signal=5, m=2)", lambda: make_world(3, 10, 5.0, 2.0), 2.0)


if __name__ == "__main__":
    run_killer_gallery()
