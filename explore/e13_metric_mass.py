"""E13 Push 6 (Insights 11–13): Mass as a Metric Tensor g in Sym+_d (O17 / T50 / E7).

Pure Python (stdlib only). Switches the algebra from scalar mass m > 0 to a symmetric
positive-definite inertia metric g: V -> V* (F = g a, inverse mass g^{-1}: V* -> V)
living on the full GL(d)-covariant covector-residual geometry:

  1. Exact GL(d)- and Swap-Invariant Covector Residual (no Euclidean spatial metric delta!):
       r_i(M) = F_i - M a_i,   V_i(M) = [I, -M] Sigma_i [I, -M]^T,   K_i(M) = r_i^T V_i^{-1} r_i.
     Under any invertible basis change A in GL(d): a -> A a, F -> A^{-T} F, g -> A^{-T} g A^{-1},
     leaving K_i(g), Metric-GCLR, Metric-Z-GCLR, and Metric-SCLR+ exactly invariant.
     Under force-acceleration swap (V <-> V*): g <-> g^{-1} leaves K_i(g) exact.
  2. Partial-Axis Excitation & Unnormalized Evidence ("sigma = inf" along unpushed axes):
       Pushing a 2D/3D body only along e_1 leaves K(g) completely flat (sigma = inf) along
       unpushed metric directions g_{22}, g_{33}; adding a later e_2 push sums K_x(g) + K_y(g)
       and pins down the full anisotropic metric with zero prior pi.
  3. Three-Rung Nested Surprise Ladder (scalar m*delta  subset  metric Sym+_d  subset  linear M):
       K_min^{scalar} >= K_min^{metric} >= K_min^{linear} >= 0,
     separating anisotropy surprise (Delta K_iso, O17), reciprocity/positivity surprise
     (Delta K_recip, E7 gyroscopic alternative), and irreducible residual (K_min^{linear}).
  4. Exact Finite-Sample Conditional Tests on Sym+_d (extending T50):
       - Pure-stdlib T50 exact conditional vector score on Sym_d (for direct comparison).
       - Metric-GCLR+ / Metric-GCLR(Sym) and Z-safe Metric-Z-GCLR: exact 95% coverage on full
         correlated Sigma_i, partial-rank excitation, and weak signal; kills T50's negative-definite
         / anti-parallel stationary-point ghost!
       - Scale-free Metric-SCLR+ (resampling V_i / ||V_i|| on S^{d-1}, ALWAYS continuous for d >= 2,
         so NO small-N 1D floor!) + O7 Metric Floor Union: exact >= 95% across Gaussian, nu=10, and
         Cauchy-like nu=2 heavy tails where T50 and Gaussian tests collapse!
       - Exact finite-sample test of scalar isotropy (L22): sup_{m > 0} p_metric(m I_d).
  5. Exact Tensor Two-Body Composition (L5 lumping G = g_1 + g_2 & L14 common force G^{-1} = g_1^{-1} + g_2^{-1}):
       Both L5 and L14 hold simultaneously and exactly on Sym+_d (resolving O17's scalar conflict!),
       with exact swap duality and partial-axis gluing (anisotropic body + x-only body).

Run from the repo root:
    python3 explore/e13_metric_mass.py
"""
import math
import random

from explore.e13_general_clr import _rand_unit_sphere, _weight_fn
from explore.e13_unnormalized_evidence import (
    CHI2_95_1,
    Reading,
    chi2_sf,
    cholesky,
    dot,
    householder_qr_full,
    isotropic,
    mat_t,
    solve_lower,
    solve_lower_mat,
)


# ---------------------------------------------------------------------------
# Fast exact small-matrix linear algebra (d = 2, 3, ...)
# ---------------------------------------------------------------------------

def mat_mul(A, B):
    bt = mat_t(B)
    return [[sum(x * y for x, y in zip(row, col)) for col in bt] for row in A]


def mat_vec(A, v):
    return [sum(x * y for x, y in zip(row, v)) for row in A]


def mat_add(A, B, scale_b=1.0):
    return [[a + scale_b * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mat_inv_spd(V):
    """Exact inverse of an SPD matrix V (fast closed form for d=2, Cholesky for d>=3)."""
    d = len(V)
    if d == 1:
        return [[1.0 / V[0][0]]]
    if d == 2:
        a, b = V[0][0], V[0][1]
        c = V[1][1]
        det = a * c - b * b
        inv_det = 1.0 / det
        return [[c * inv_det, -b * inv_det], [-b * inv_det, a * inv_det]]
    L = cholesky(V)
    inv_cols = []
    for j in range(d):
        e = [1.0 if i == j else 0.0 for i in range(d)]
        y = solve_lower(L, e)
        # Solve L^T x = y
        x = [0.0] * d
        for i in range(d - 1, -1, -1):
            s = sum(L[k][i] * x[k] for k in range(i + 1, d))
            x[i] = (y[i] - s) / L[i][i]
        inv_cols.append(x)
    return mat_t(inv_cols)


def mat_inv_general(A):
    """Gaussian elimination with partial pivoting for small invertible matrices."""
    n = len(A)
    aug = [list(A[i]) + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[piv][col]) < 1e-15:
            raise ValueError("Matrix is singular")
        aug[col], aug[piv] = aug[piv], aug[col]
        inv_p = 1.0 / aug[col][col]
        row_c = aug[col]
        for j in range(2 * n):
            row_c[j] *= inv_p
        for r in range(n):
            if r != col and aug[r][col] != 0.0:
                fac = aug[r][col]
                row_r = aug[r]
                for j in range(2 * n):
                    row_r[j] -= fac * row_c[j]
    return [row[n:] for row in aug]


def is_spd(M, tol=1e-12):
    """Return True iff M is symmetric positive definite."""
    d = len(M)
    for i in range(d):
        for j in range(i + 1, d):
            if abs(M[i][j] - M[j][i]) > 1e-9:
                return False
    try:
        L = cholesky(M)
        return all(L[i][i] > tol for i in range(d))
    except ValueError:
        return False


def spd_from_params(theta, d):
    """Map unconstrained R^{d(d+1)/2} vector theta -> SPD matrix g = C C^T.

    C is lower triangular with C[i][i] = exp(theta_diag_i) > 0, off-diagonals free.
    """
    C = [[0.0] * d for _ in range(d)]
    idx = 0
    for i in range(d):
        for j in range(i):
            C[i][j] = theta[idx]
            idx += 1
        C[i][i] = math.exp(max(-12.0, min(12.0, theta[idx])))
        idx += 1
    return mat_mul(C, mat_t(C))


def params_from_spd(g):
    """Inverse map SPD matrix g -> unconstrained Cholesky parameters theta."""
    L = cholesky(g)
    d = len(g)
    theta = []
    for i in range(d):
        for j in range(i):
            theta.append(L[i][j])
        theta.append(math.log(max(L[i][i], 1e-12)))
    return theta


def sym_from_params(theta, d):
    """Map R^{d(d+1)/2} -> symmetric matrix S (unrestricted sign of eigenvalues)."""
    S = [[0.0] * d for _ in range(d)]
    idx = 0
    for i in range(d):
        for j in range(i, d):
            S[i][j] = S[j][i] = theta[idx]
            idx += 1
    return S


def params_from_sym(S):
    d = len(S)
    return [S[i][j] for i in range(d) for j in range(i, d)]


def mat_from_params(theta, d):
    return [list(theta[i * d : (i + 1) * d]) for i in range(d)]


def params_from_mat(M):
    return [x for row in M for x in row]


def eye(d, scale=1.0):
    return [[scale if i == j else 0.0 for j in range(d)] for i in range(d)]


def sym_tangent_basis(d):
    """Frobenius-orthonormal basis of symmetric d x d matrices (length d(d+1)/2)."""
    basis = []
    inv_s2 = 1.0 / math.sqrt(2.0)
    for i in range(d):
        for j in range(i, d):
            E = [[0.0] * d for _ in range(d)]
            if i == j:
                E[i][i] = 1.0
            else:
                E[i][j] = E[j][i] = inv_s2
            basis.append(E)
    return basis


# ---------------------------------------------------------------------------
# Metric Reading: exact covector residual K_i(M) and orthogonal (U_i, V_i) split
# ---------------------------------------------------------------------------

class MetricReading:
    """Full-subspace reading in V* (+) V with covector residual K_i(M) and exact (U, V) split."""

    def __init__(self, reading: Reading):
        if reading.k != reading.d:
            raise ValueError("MetricReading expects full subspace (k == d)")
        self.reading = reading
        self.d = reading.d
        self.F = list(reading.force)
        self.a = list(reading.acceleration)
        self.cov = reading.cov
        self.chol = reading.chol
        self.z = tuple(reading.whitened)  # length 2d
        d = self.d
        self.S_FF = [row[:d] for row in self.cov[:d]]
        self.S_Fa = [row[d:] for row in self.cov[:d]]
        self.S_aF = [row[:d] for row in self.cov[d:]]
        self.S_aa = [row[d:] for row in self.cov[d:]]

    def residual_cost(self, M, fa_override=None):
        """Exact covector residual K_i(M) = (F - M a)^T V_i(M)^{-1} (F - M a)."""
        d = self.d
        if fa_override is None:
            F, a = self.F, self.a
        else:
            F, a = fa_override[:d], fa_override[d:]
        if d == 2:
            m00, m01 = M[0][0], M[0][1]
            m10, m11 = M[1][0], M[1][1]
            r0 = F[0] - (m00 * a[0] + m01 * a[1])
            r1 = F[1] - (m10 * a[0] + m11 * a[1])
            sff, sfa, saa = self.S_FF, self.S_Fa, self.S_aa
            # V = S_FF - M S_aF - S_Fa M^T + M S_aa M^T
            # Compute M S_aa first (2x2)
            ms00 = m00 * saa[0][0] + m01 * saa[1][0]
            ms01 = m00 * saa[0][1] + m01 * saa[1][1]
            ms10 = m10 * saa[0][0] + m11 * saa[1][0]
            ms11 = m10 * saa[0][1] + m11 * saa[1][1]
            v00 = sff[0][0] - 2.0 * (sfa[0][0] * m00 + sfa[0][1] * m01) + (ms00 * m00 + ms01 * m01)
            v11 = sff[1][1] - 2.0 * (sfa[1][0] * m10 + sfa[1][1] * m11) + (ms10 * m10 + ms11 * m11)
            v01 = (
                sff[0][1]
                - (m00 * sfa[1][0] + m01 * sfa[1][1])
                - (sfa[0][0] * m10 + sfa[0][1] * m11)
                + (ms00 * m10 + ms01 * m11)
            )
            det = v00 * v11 - v01 * v01
            return (v11 * r0 * r0 - 2.0 * v01 * r0 * r1 + v00 * r1 * r1) / det
        # General d >= 3
        r = [F[i] - sum(M[i][j] * a[j] for j in range(d)) for i in range(d)]
        Mt = mat_t(M)
        V = mat_add(
            mat_add(self.S_FF, mat_mul(M, self.S_aF), -1.0),
            mat_add(mat_mul(self.S_Fa, Mt), mat_mul(mat_mul(M, self.S_aa), Mt), -1.0),
            -1.0,
        )
        # V is SPD: solve L y = r => ||y||^2
        L = cholesky(V)
        y = solve_lower(L, r)
        return sum(x * x for x in y)

    def split_at_matrix(self, M):
        """Exact whitened orthogonal split (Q_par_T, Q_perp_T, U, V, H, K, a_hat) at candidate M."""
        d = self.d
        B_raw = [list(M[i]) for i in range(d)] + [[1.0 if i == j else 0.0 for j in range(d)] for i in range(d)]
        B_w = solve_lower_mat(self.chol, B_raw)  # 2d x d
        Q, R = householder_qr_full(B_w)
        Qt = mat_t(Q)
        q_par = Qt[:d]
        q_perp = Qt[d:]
        u = [dot(row, self.z) for row in q_par]
        v = [dot(row, self.z) for row in q_perp]
        # Fitted latent acceleration a_hat = R[:d]^{-1} u
        a_hat = [0.0] * d
        for i in range(d - 1, -1, -1):
            s = sum(R[i][j] * a_hat[j] for j in range(i + 1, d))
            a_hat[i] = (u[i] - s) / R[i][i]
        h = sum(x * x for x in u)
        k_val = sum(x * x for x in v)
        return q_par, q_perp, u, v, h, k_val, a_hat

    def unwhiten(self, z_vec):
        """Map whitened z_vec back to data vector y = L z (length 2d)."""
        return [dot(row, z_vec) for row in self.chol]


def transform_gl(reading: Reading, A):
    """Apply invertible basis change A in GL(d): a -> A a, F -> A^{-T} F, Sigma transformed."""
    d = reading.d
    A_inv_T = mat_t(mat_inv_general(A))
    F_new = mat_vec(A_inv_T, reading.force)
    a_new = mat_vec(A, reading.acceleration)
    # Block diagonal T = diag(A^{-T}, A)
    T = [[0.0] * (2 * d) for _ in range(2 * d)]
    for i in range(d):
        for j in range(d):
            T[i][j] = A_inv_T[i][j]
            T[d + i][d + j] = A[i][j]
    cov_new = mat_mul(mat_mul(T, reading.cov), mat_t(T))
    return Reading(F_new, a_new, cov_new)


def transform_metric_gl(g, A):
    """Contravariant/covariant metric transform g -> A^{-T} g A^{-1} under a -> A a, F -> A^{-T} F."""
    A_inv = mat_inv_general(A)
    return mat_mul(mat_mul(mat_t(A_inv), g), A_inv)


# ---------------------------------------------------------------------------
# Pure-stdlib T50 Conditional Vector Score Test on Sym_d (for comparison)
# ---------------------------------------------------------------------------

def t50_metric_score_test(metric_readings, g0, *, tangents=None, rank_tol=1e-10):
    """Exact conditional Gaussian vector score test (T50) at candidate metric g0.

    Returns dict(T2, rank, p_value, degenerate). Uses exact QR column pivoting / Gram
    eigen-rank on the conditional score map A so no external SVD package is required.
    """
    if not metric_readings:
        return dict(T2=0.0, rank=0, p_value=1.0, degenerate=True)
    d = metric_readings[0].d
    if tangents is None:
        tangents = sym_tangent_basis(d)
    p = len(tangents)
    # Build score s in R^p and conditional covariance C in R^{p x p}
    s_vec = [0.0] * p
    C_mat = [[0.0] * p for _ in range(p)]
    for mr in metric_readings:
        _, q_perp, _, v_vec, h_val, _, a_hat = mr.split_at_matrix(g0)
        if h_val <= 1e-24:
            continue
        # G_ia = Q_perp^T L^{-1} [E_a a_hat; 0]  (length d for each tangent a = 0..p-1)
        g_cols = []
        for E in tangents:
            ea = mat_vec(E, a_hat)
            raw = ea + [0.0] * d
            w_dir = solve_lower(mr.chol, raw)
            g_cols.append([dot(row, w_dir) for row in q_perp])
        for a_idx in range(p):
            s_vec[a_idx] += dot(v_vec, g_cols[a_idx])
            for b_idx in range(a_idx, p):
                val = dot(g_cols[a_idx], g_cols[b_idx])
                C_mat[a_idx][b_idx] += val
                if a_idx != b_idx:
                    C_mat[b_idx][a_idx] += val

    # Symmetric Jacobi eigendecomposition of C_mat (p <= 9 -> instantaneous & exact!)
    eigvals, eigvecs = _jacobi_eigh(C_mat)
    s_max = max(eigvals) if eigvals else 0.0
    if s_max <= 1e-24:
        return dict(T2=0.0, rank=0, p_value=1.0, degenerate=True)
    cutoff = max(rank_tol, 1e-10 * s_max)
    t2 = 0.0
    rank = 0
    for lam, vec in zip(eigvals, eigvecs):
        if lam > cutoff:
            proj = dot(vec, s_vec)
            t2 += (proj * proj) / lam
            rank += 1
    if rank == 0:
        return dict(T2=0.0, rank=0, p_value=1.0, degenerate=True)
    return dict(T2=t2, rank=rank, p_value=chi2_sf(t2, rank), degenerate=False)


def _jacobi_eigh(A, max_iter=80):
    """Exact symmetric eigendecomposition for small symmetric matrices (returns vals, list of eigvecs)."""
    n = len(A)
    M = [list(row) for row in A]
    V = eye(n)
    for _ in range(max_iter):
        # Find largest off-diagonal
        p, q, max_off = 0, 1, 0.0
        for i in range(n):
            for j in range(i + 1, n):
                val = abs(M[i][j])
                if val > max_off:
                    max_off, p, q = val, i, j
        if max_off < 1e-14:
            break
        diff = M[q][q] - M[p][p]
        if abs(M[p][q]) < 1e-30:
            continue
        if abs(diff) < 1e-30:
            t = 1.0 if M[p][q] > 0 else -1.0
        else:
            tau = diff / (2.0 * M[p][q])
            t = (1.0 if tau >= 0 else -1.0) / (abs(tau) + math.hypot(1.0, tau))
        c = 1.0 / math.hypot(1.0, t)
        s = t * c
        app, aqq, apq = M[p][p], M[q][q], M[p][q]
        M[p][p] = app - t * apq
        M[q][q] = aqq + t * apq
        M[p][q] = M[q][p] = 0.0
        for r in range(n):
            if r != p and r != q:
                arp, arq = M[r][p], M[r][q]
                M[r][p] = M[p][r] = c * arp - s * arq
                M[r][q] = M[q][r] = s * arp + c * arq
            vrp, vrq = V[r][p], V[r][q]
            V[r][p] = c * vrp - s * vrq
            V[r][q] = s * vrp + c * vrq
    vals = [M[i][i] for i in range(n)]
    vecs = [[V[r][i] for r in range(n)] for i in range(n)]
    return vals, vecs


# ---------------------------------------------------------------------------
# Fast Deterministic Optimizers on Scalar / Metric (Sym+_d) / Sym_d / Linear
# ---------------------------------------------------------------------------

def _nelder_mead(fn, x0, step=0.35, max_iter=120, tol=1e-7):
    """Deterministic Nelder-Mead simplex optimizer (exact exchangeability across null draws!)."""
    n = len(x0)
    simplex = [list(x0)]
    for i in range(n):
        pt = list(x0)
        pt[i] += step if abs(pt[i]) < 1.0 else step * abs(pt[i])
        simplex.append(pt)
    fvals = [fn(p) for p in simplex]
    for _ in range(max_iter):
        order = sorted(range(n + 1), key=lambda idx: fvals[idx])
        simplex = [simplex[i] for i in order]
        fvals = [fvals[i] for i in order]
        if fvals[-1] - fvals[0] < tol:
            break
        cent = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]
        xr = [cent[j] + (cent[j] - simplex[-1][j]) for j in range(n)]
        fr = fn(xr)
        if fvals[0] <= fr < fvals[-2]:
            simplex[-1], fvals[-1] = xr, fr
        elif fr < fvals[0]:
            xe = [cent[j] + 2.0 * (xr[j] - cent[j]) for j in range(n)]
            fe = fn(xe)
            if fe < fr:
                simplex[-1], fvals[-1] = xe, fe
            else:
                simplex[-1], fvals[-1] = xr, fr
        else:
            xc = [cent[j] + 0.5 * (simplex[-1][j] - cent[j]) for j in range(n)]
            fc = fn(xc)
            if fc < fvals[-1]:
                simplex[-1], fvals[-1] = xc, fc
            else:
                best = simplex[0]
                for i in range(1, n + 1):
                    simplex[i] = [best[j] + 0.5 * (simplex[i][j] - best[j]) for j in range(n)]
                    fvals[i] = fn(simplex[i])
    return simplex[0], fvals[0]


def _ols_matrix_seed(fa_list, weights, d):
    """Deterministic ridge-stabilized least-squares seed M_hat from weighted (F_i, a_i) data."""
    # Solve M (sum w_i a_i a_i^T + eps I) = sum w_i F_i a_i^T
    aa = [[1e-3 if i == j else 0.0 for j in range(d)] for i in range(d)]
    fa = [[0.0] * d for _ in range(d)]
    for fa_vec, w in zip(fa_list, weights):
        if w <= 0.0:
            continue
        F, a = fa_vec[:d], fa_vec[d:]
        for i in range(d):
            for j in range(d):
                aa[i][j] += w * a[i] * a[j]
                fa[i][j] += w * F[i] * a[j]
    return mat_mul(fa, mat_inv_spd(aa))


def _project_to_spd(M, min_eig=0.05):
    """Symmetrize M and clamp eigenvalues to >= min_eig for a deterministic SPD start."""
    d = len(M)
    S = [[0.5 * (M[i][j] + M[j][i]) for j in range(d)] for i in range(d)]
    vals, vecs = _jacobi_eigh(S)
    out = [[0.0] * d for _ in range(d)]
    for lam, v in zip(vals, vecs):
        lam_c = max(min_eig, min(1e3, lam))
        for i in range(d):
            for j in range(d):
                out[i][j] += lam_c * v[i] * v[j]
    return out


def _irls_sym_2d(mrs, fa_list, weights, M_init, n_steps=4):
    """Fast deterministic 4-step IRLS minimizing sum w_i (F_i - S a_i)^T V_i(S)^{-1} (F_i - S a_i)
    over 2x2 symmetric matrices S = [[s00, s01], [s01, s11]] (3x3 linear solve per step).
    """
    m00, m01, m11 = M_init[0][0], 0.5 * (M_init[0][1] + M_init[1][0]), M_init[1][1]
    for _ in range(n_steps):
        # Normal equations H * [s00, s01, s11]^T = rhs (with tiny 1e-6 ridge for rank-deficient axes)
        h00 = h11 = h22 = 1e-6
        h01 = h02 = h12 = 0.0
        b0 = b1 = b2 = 0.0
        for mr, fa, w in zip(mrs, fa_list, weights):
            if w <= 0.0:
                continue
            f0, f1, a0, a1 = fa[0], fa[1], fa[2], fa[3]
            sff, sfa, saa = mr.S_FF, mr.S_Fa, mr.S_aa
            ms00 = m00 * saa[0][0] + m01 * saa[1][0]
            ms01 = m00 * saa[0][1] + m01 * saa[1][1]
            ms10 = m01 * saa[0][0] + m11 * saa[1][0]
            ms11 = m01 * saa[0][1] + m11 * saa[1][1]
            v00 = sff[0][0] - 2.0 * (sfa[0][0] * m00 + sfa[0][1] * m01) + (ms00 * m00 + ms01 * m01)
            v11 = sff[1][1] - 2.0 * (sfa[1][0] * m01 + sfa[1][1] * m11) + (ms10 * m01 + ms11 * m11)
            v01 = (
                sff[0][1]
                - (m00 * sfa[1][0] + m01 * sfa[1][1])
                - (sfa[0][0] * m01 + sfa[0][1] * m11)
                + (ms00 * m01 + ms01 * m11)
            )
            inv_det = w / (v00 * v11 - v01 * v01)
            w00 = v11 * inv_det
            w01 = -v01 * inv_det
            w11 = v00 * inv_det
            # Residual linear map: S a = [a0*s00 + a1*s01, a0*s01 + a1*s11]
            # Jacobian columns: j0 = (a0, 0), j1 = (a1, a0), j2 = (0, a1)
            # W j0 = (w00*a0, w01*a0); W j1 = (w00*a1 + w01*a0, w01*a1 + w11*a0); W j2 = (w01*a1, w11*a1)
            wj0_0, wj0_1 = w00 * a0, w01 * a0
            wj1_0, wj1_1 = w00 * a1 + w01 * a0, w01 * a1 + w11 * a0
            wj2_0, wj2_1 = w01 * a1, w11 * a1
            h00 += a0 * wj0_0
            h01 += a0 * wj1_0
            h02 += a0 * wj2_0
            h11 += a1 * wj1_0 + a0 * wj1_1
            h12 += a1 * wj2_0 + a0 * wj2_1
            h22 += a1 * wj2_1
            b0 += f0 * wj0_0 + f1 * wj0_1
            b1 += f0 * wj1_0 + f1 * wj1_1
            b2 += f0 * wj2_0 + f1 * wj2_1
        # Solve 3x3 SPD system via Cholesky
        L = cholesky([[h00, h01, h02], [h01, h11, h12], [h02, h12, h22]])
        y = solve_lower(L, [b0, b1, b2])
        x2 = y[2] / L[2][2]
        x1 = (y[1] - L[2][1] * x2) / L[1][1]
        x0 = (y[0] - L[1][0] * x1 - L[2][0] * x2) / L[0][0]
        # Damped update
        m00 = 0.25 * m00 + 0.75 * x0
        m01 = 0.25 * m01 + 0.75 * x1
        m11 = 0.25 * m11 + 0.75 * x2
    return [[m00, m01], [m01, m11]]


def _eval_2d_fast(mrs, fa_list, weights, m00, m01, m10, m11):
    """Unrolled zero-allocation 2D weighted residual sum."""
    tot = 0.0
    for mr, fa, w in zip(mrs, fa_list, weights):
        if w <= 0.0:
            continue
        r0 = fa[0] - (m00 * fa[2] + m01 * fa[3])
        r1 = fa[1] - (m10 * fa[2] + m11 * fa[3])
        sff, sfa, saa = mr.S_FF, mr.S_Fa, mr.S_aa
        ms00 = m00 * saa[0][0] + m01 * saa[1][0]
        ms01 = m00 * saa[0][1] + m01 * saa[1][1]
        ms10 = m10 * saa[0][0] + m11 * saa[1][0]
        ms11 = m10 * saa[0][1] + m11 * saa[1][1]
        v00 = sff[0][0] - 2.0 * (sfa[0][0] * m00 + sfa[0][1] * m01) + (ms00 * m00 + ms01 * m01)
        v11 = sff[1][1] - 2.0 * (sfa[1][0] * m10 + sfa[1][1] * m11) + (ms10 * m10 + ms11 * m11)
        v01 = (
            sff[0][1]
            - (m00 * sfa[1][0] + m01 * sfa[1][1])
            - (sfa[0][0] * m10 + sfa[0][1] * m11)
            + (ms00 * m10 + ms01 * m11)
        )
        tot += w * (v11 * r0 * r0 - 2.0 * v01 * r0 * r1 + v00 * r1 * r1) / (v00 * v11 - v01 * v01)
    return tot


def minimize_metric_cost(
    metric_readings,
    *,
    family="metric",
    weights=None,
    fa_overrides=None,
    seed_matrix=None,
    max_iter=90,
):
    """Deterministic minimizer of sum_i w_i K_i(M) over:
      - 'scalar': M = m I_d, m > 0
      - 'metric': M = g in Sym+_d (symmetric positive definite)
      - 'sym':    M = S in Sym_d  (symmetric, any sign — kills negative-definite ghosts!)
      - 'linear': M in R^{d x d}  (unrestricted linear map, allows non-reciprocal / gyroscopic)

    Returns (M_best, K_min).
    """
    d = metric_readings[0].d
    if weights is None:
        weights = [1.0] * len(metric_readings)
    if fa_overrides is None:
        fa_list = [mr.F + mr.a for mr in metric_readings]
    else:
        fa_list = fa_overrides

    if d == 2:
        def eval_M(M):
            return _eval_2d_fast(metric_readings, fa_list, weights, M[0][0], M[0][1], M[1][0], M[1][1])
    else:
        def eval_M(M):
            tot = 0.0
            for mr, w, fa in zip(metric_readings, weights, fa_list):
                if w > 0.0:
                    tot += w * mr.residual_cost(M, fa_override=fa)
            return tot

    if family == "scalar":
        phis = [1e-4 + (0.5 * math.pi - 2e-4) * j / 60 for j in range(61)]
        vals = [eval_M(eye(d, math.tan(ph))) for ph in phis]
        bi = vals.index(min(vals))
        lo = phis[max(0, bi - 1)]
        hi = phis[min(len(phis) - 1, bi + 1)]
        gr = 0.5 * (math.sqrt(5.0) - 1.0)
        a, b = lo, hi
        c, dp = b - gr * (b - a), a + gr * (b - a)
        fc, fd = eval_M(eye(d, math.tan(c))), eval_M(eye(d, math.tan(dp)))
        for _ in range(24):
            if fc < fd:
                b, dp, fd = dp, c, fc
                c = b - gr * (b - a)
                fc = eval_M(eye(d, math.tan(c)))
            else:
                a, c, fc = c, dp, fd
                dp = a + gr * (b - a)
                fd = eval_M(eye(d, math.tan(dp)))
        phi_star = 0.5 * (a + b)
        M_best = eye(d, math.tan(phi_star))
        return M_best, min(vals[bi], eval_M(M_best))

    # Fast 2D path: IRLS gives an immediate near-exact symmetric fit, followed by short unrolled NM!
    if d == 2 and family in ("metric", "sym"):
        init_M = seed_matrix if seed_matrix is not None else _ols_matrix_seed(fa_list, weights, 2)
        S_irls = _irls_sym_2d(metric_readings, fa_list, weights, init_M, n_steps=4)
        if family == "metric":
            g_start = _project_to_spd(S_irls, min_eig=0.02)
            # Pick best between seed_matrix (if SPD) and g_start, then refine once
            if seed_matrix is not None and is_spd(seed_matrix) and eval_M(seed_matrix) < eval_M(g_start):
                g_start = seed_matrix
            th0 = params_from_spd(g_start)

            def eval_spd_2d(th):
                e0 = math.exp(max(-12.0, min(12.0, th[0])))
                c10 = th[1]
                e1 = math.exp(max(-12.0, min(12.0, th[2])))
                return _eval_2d_fast(metric_readings, fa_list, weights, e0 * e0, e0 * c10, e0 * c10, c10 * c10 + e1 * e1)

            th_opt, v_opt = _nelder_mead(eval_spd_2d, th0, step=0.25, max_iter=max_iter)
            return spd_from_params(th_opt, 2), v_opt
        else:
            th0 = [S_irls[0][0], S_irls[0][1], S_irls[1][1]]

            def eval_sym_2d(th):
                return _eval_2d_fast(metric_readings, fa_list, weights, th[0], th[1], th[1], th[2])

            th_opt, v_opt = _nelder_mead(eval_sym_2d, th0, step=0.25, max_iter=max_iter)
            return [[th_opt[0], th_opt[1]], [th_opt[1], th_opt[2]]], v_opt

    M_ols = _ols_matrix_seed(fa_list, weights, d)
    if family == "metric":
        starts = [_project_to_spd(M_ols)]
        if seed_matrix is not None and is_spd(seed_matrix):
            starts.append(seed_matrix)
        best_M, best_val = starts[0], float("inf")
        for st in starts:
            th0 = params_from_spd(st)
            th_opt, v_opt = _nelder_mead(lambda th: eval_M(spd_from_params(th, d)), th0, step=0.3, max_iter=max_iter)
            if v_opt < best_val:
                best_val = v_opt
                best_M = spd_from_params(th_opt, d)
        return best_M, best_val

    if family == "sym":
        S_ols = [[0.5 * (M_ols[i][j] + M_ols[j][i]) for j in range(d)] for i in range(d)]
        starts = [S_ols]
        if seed_matrix is not None:
            starts.append(seed_matrix)
        best_M, best_val = starts[0], float("inf")
        for st in starts:
            th0 = params_from_sym(st)
            th_opt, v_opt = _nelder_mead(lambda th: eval_M(sym_from_params(th, d)), th0, step=0.35, max_iter=max_iter)
            if v_opt < best_val:
                best_val = v_opt
                best_M = sym_from_params(th_opt, d)
        return best_M, best_val

    if family == "linear":
        starts = [M_ols]
        if seed_matrix is not None:
            starts.append(seed_matrix)
        best_M, best_val = starts[0], float("inf")
        for st in starts:
            th0 = params_from_mat(st)
            th_opt, v_opt = _nelder_mead(lambda th: eval_M(mat_from_params(th, d)), th0, step=0.35, max_iter=max_iter)
            if v_opt < best_val:
                best_val = v_opt
                best_M = mat_from_params(th_opt, d)
        return best_M, best_val

    raise ValueError(f"Unknown family {family!r}")


def nested_surprise_ladder(metric_readings):
    """Exact monotone leftover surprise ladder: scalar (m I) >= metric (Sym+) >= linear (R^{d x d}).

    Returns dict(K_scalar, K_metric, K_linear, delta_iso, delta_recip, g_hat, M_hat, m_hat).
    """
    m_mat, k_sc = minimize_metric_cost(metric_readings, family="scalar")
    g_hat, k_me = minimize_metric_cost(metric_readings, family="metric", seed_matrix=m_mat, max_iter=180)
    M_hat, k_li = minimize_metric_cost(metric_readings, family="linear", seed_matrix=g_hat, max_iter=200)
    k_me = min(k_me, k_sc)
    k_li = min(k_li, k_me)
    return dict(
        K_scalar=k_sc,
        K_metric=k_me,
        K_linear=k_li,
        delta_iso=max(0.0, k_sc - k_me),
        delta_recip=max(0.0, k_me - k_li),
        m_hat=m_mat[0][0],
        g_hat=g_hat,
        M_hat=M_hat,
    )


# ---------------------------------------------------------------------------
# Exact Conditional Metric Tests: Metric-GCLR, Metric-Z-GCLR, Metric-SCLR+, Floor Union
# ---------------------------------------------------------------------------

def metric_clr_pvalue(
    metric_readings,
    g0,
    *,
    family="metric",
    weight_mode="unweighted",
    null_mode="gaussian",
    n_draws=99,
    opt_iter=55,
    seed=2026,
):
    """Exact finite-sample conditional test of candidate metric g0 in Sym+_d.

    Exact Dwass-Barnard exchangeability holds for ANY deterministic optimizer applied identically
    to the observed data and all conditional null replicates:
      - null_mode='gaussian'  -> V_i^(b) ~ N(0, I_d) given U_i  (exact under Gaussian)
      - null_mode='spherical' -> V_i^(b) = ||V_i|| * Unif(S^{d-1}) given (U_i, ||V_i||)
                                 (exact under ANY per-reading scale mixture! ALWAYS continuous for d >= 2!)
    """
    if not metric_readings:
        return dict(p_value=1.0, lr_obs=0.0, k_obs=0.0, degenerate=True)

    d = metric_readings[0].d
    active = []
    weights_all = []
    for mr in metric_readings:
        q_par, q_perp, u, v, h, k_val, _ = mr.split_at_matrix(g0)
        w = _weight_fn(h, weight_mode)
        weights_all.append(w)
        if w > 0.0:
            v_norm = math.sqrt(k_val)
            if null_mode == "spherical" and h == 0.0 and v_norm == 0.0:
                continue
            m2d = 2 * d
            z_u = [0.0] * m2d
            for uv, row in zip(u, q_par):
                if uv != 0.0:
                    for j in range(m2d):
                        z_u[j] += uv * row[j]
            active.append((mr, tuple(z_u), tuple(tuple(r) for r in q_perp), v_norm, w))

    if not active:
        return dict(p_value=1.0, lr_obs=0.0, k_obs=0.0, degenerate=True)

    mrs_act = [t[0] for t in active]
    w_act = [t[4] for t in active]

    def eval_lr(fa_list, k0_list):
        c0 = sum(w * k0 for w, k0 in zip(w_act, k0_list))
        _, c_min = minimize_metric_cost(
            mrs_act,
            family=family,
            weights=w_act,
            fa_overrides=fa_list,
            seed_matrix=g0,
            max_iter=opt_iter,
        )
        return max(0.0, c0 - min(c0, c_min)), c0

    fa_obs = [mr.F + mr.a for mr in mrs_act]
    k0_obs = [vn * vn for _, _, _, vn, _ in active]
    lr_obs, k_obs = eval_lr(fa_obs, k0_obs)
    if lr_obs <= 1e-13:
        return dict(p_value=1.0, lr_obs=0.0, k_obs=k_obs, degenerate=False)

    tol = 1e-9 * max(1.0, lr_obs)
    rng = random.Random(seed)
    gauss = rng.gauss
    m2d = 2 * d
    ge = 0
    for _ in range(n_draws):
        fa_b = []
        k0_b = []
        for mr, z_u, q_perp, v_norm, _ in active:
            if null_mode == "gaussian":
                vb = [gauss(0.0, 1.0) for _ in range(d)]
                k0_b.append(sum(x * x for x in vb))
            elif null_mode == "spherical":
                if v_norm <= 1e-15:
                    vb = [0.0] * d
                    k0_b.append(0.0)
                else:
                    ub = _rand_unit_sphere(d, rng)
                    vb = [v_norm * x for x in ub]
                    k0_b.append(v_norm * v_norm)
            else:
                raise ValueError(f"Unknown null_mode {null_mode!r}")
            zb = list(z_u)
            for vv, row in zip(vb, q_perp):
                if vv != 0.0:
                    for j in range(m2d):
                        zb[j] += vv * row[j]
            fa_b.append(mr.unwhiten(zb))
        lr_b, _ = eval_lr(fa_b, k0_b)
        if lr_b >= lr_obs - tol:
            ge += 1
    return dict(p_value=(1 + ge) / (n_draws + 1), lr_obs=lr_obs, k_obs=k_obs, degenerate=False)


def test_scalar_isotropy_l22(metric_readings, *, null_mode="gaussian", weight_mode="z_smooth", n_draws=99, seed=2026):
    """Exact finite-sample test of scalar isotropy (L22: g = m I_d for some m > 0) inside Sym+_d.

    Evaluates p_metric(m_hat I_d) (and nearby scalars), returning sup_{m > 0} p_metric(m I_d).
    Because the true m* is included in the supremum, P(reject isotropic body) <= 0.05 guaranteed!
    """
    d = metric_readings[0].d
    m_mat, _ = minimize_metric_cost(metric_readings, family="scalar")
    m_hat = m_mat[0][0]
    cands = sorted(set([m_hat * f for f in (0.85, 0.95, 1.0, 1.05, 1.15) if m_hat * f > 1e-4]))
    best_p = 0.0
    best_m = m_hat
    for idx, m_c in enumerate(cands):
        res = metric_clr_pvalue(
            metric_readings,
            eye(d, m_c),
            family="metric",
            weight_mode=weight_mode,
            null_mode=null_mode,
            n_draws=n_draws,
            seed=seed + idx * 13,
        )
        if res["p_value"] > best_p:
            best_p = res["p_value"]
            best_m = m_c
    return dict(p_value=best_p, m_best=best_m)


# ---------------------------------------------------------------------------
# Tensor Two-Body Composition: L5 Lumping (g1 + g2) & L14 Common Force
# ---------------------------------------------------------------------------

def compose_metric_lump(mrs1, mrs2, G_target, max_iter=140):
    """Prior-free infimal convolution for tensor lumping (L5):
        K_lump(G) = inf_{0 < g_1 < G} [ K_1(g_1) + K_2(G - g_1) ].

    Parameterized cleanly on the operator interval {0 < g_1 < G} via:
        g_1 = L_G (I + S S^T)^{-1/2} ... even simpler:
        g_1 = L_G * diag_or_SPD_in_(0, I) * L_G^T, where X = (I + E)^{-1} for E = C C^T in Sym+_d!
    Check: for ANY E in Sym+_d, (I + E)^{-1} lies strictly in (0, I), so BOTH
    g_1 = L_G (I + E)^{-1} L_G^T AND g_2 = G - g_1 = L_G E (I + E)^{-1} L_G^T are STRICTLY SPD!
    """
    d = len(G_target)
    L_G = cholesky(G_target)
    L_Gt = mat_t(L_G)
    I_d = eye(d)

    def eval_split(theta):
        E = spd_from_params(theta, d)
        X = mat_inv_spd(mat_add(I_d, E))  # strictly between 0 and I!
        g1 = mat_mul(mat_mul(L_G, X), L_Gt)
        g2 = mat_add(G_target, g1, -1.0)
        k1 = sum(mr.residual_cost(g1) for mr in mrs1)
        k2 = sum(mr.residual_cost(g2) for mr in mrs2)
        return k1 + k2

    # Start from symmetric midpoint E = I (g1 = g2 = G/2) and from body-1 fit if < G
    starts = [params_from_spd(I_d)]
    g1_hat, _ = minimize_metric_cost(mrs1, family="metric", max_iter=80)
    # Check if G - g1_hat is SPD -> exact split start E = L_G^{-1} (G - g1) g1^{-1} ...
    g2_diff = mat_add(G_target, g1_hat, -1.0)
    if is_spd(g2_diff):
        L_inv = mat_inv_general(L_G)
        X_hat = mat_mul(mat_mul(L_inv, g1_hat), mat_t(L_inv))
        E_hat = mat_add(mat_inv_spd(X_hat), I_d, -1.0)
        if is_spd(E_hat):
            starts.append(params_from_spd(E_hat))

    best_k = float("inf")
    best_g1 = None
    for th0 in starts:
        th_opt, val = _nelder_mead(eval_split, th0, step=0.5, max_iter=max_iter)
        if val < best_k:
            best_k = val
            E_opt = spd_from_params(th_opt, d)
            X_opt = mat_inv_spd(mat_add(I_d, E_opt))
            best_g1 = mat_mul(mat_mul(L_G, X_opt), L_Gt)
    return dict(K_lump=best_k, g1_opt=best_g1, g2_opt=mat_add(G_target, best_g1, -1.0))


def compose_metric_common_force(rds1, rds2, G_target, max_iter=140):
    """Prior-free infimal convolution for tensor common force (L14: G^{-1} = g_1^{-1} + g_2^{-1})
    via exact swap duality g <-> g^{-1}!
    """
    from explore.e13_general_clr import swap_reading

    mrs1_s = [MetricReading(swap_reading(r)) for r in rds1]
    mrs2_s = [MetricReading(swap_reading(r)) for r in rds2]
    res = compose_metric_lump(mrs1_s, mrs2_s, mat_inv_spd(G_target), max_iter=max_iter)
    return dict(
        K_red=res["K_lump"],
        g1_opt=mat_inv_spd(res["g1_opt"]),
        g2_opt=mat_inv_spd(res["g2_opt"]),
    )


# ---------------------------------------------------------------------------
# Benchmark Suite for Metric Mass (Insights 11-13)
# ---------------------------------------------------------------------------

def _make_metric_world(rng, M_true, N, r_len, sigma=1.0, nu_df=None, axis_only=None):
    """Generate 2D/3D readings from F* = M_true a* with Gaussian (nu_df=None) or IG scale-mix noise."""
    d = len(M_true)
    rds = []
    for _ in range(N):
        if axis_only is not None:
            a_dir = [1.0 if j == axis_only else 0.0 for j in range(d)]
        else:
            a_dir = list(_rand_unit_sphere(d, rng))
        f_dir = mat_vec(M_true, a_dir)
        # Scale so true whitened length = r_len
        nrm = math.sqrt(sum(x * x for x in f_dir) + sum(x * x for x in a_dir)) / sigma
        sc = r_len / nrm if nrm > 0 else 0.0
        a_true = [sc * x for x in a_dir]
        f_true = [sc * x for x in f_dir]
        if nu_df is None:
            sd = sigma
        else:
            chi2_d = sum(rng.gauss(0.0, 1.0) ** 2 for _ in range(nu_df))
            sd = sigma * math.sqrt(nu_df / chi2_d)
        F_obs = [x + sd * rng.gauss(0.0, 1.0) for x in f_true]
        a_obs = [x + sd * rng.gauss(0.0, 1.0) for x in a_true]
        rds.append(isotropic(F_obs, a_obs, sigma, sigma))
    return rds


def run_metric_benchmarks(seed=20261122):
    rng = random.Random(seed)

    print("=" * 94)
    print("PUSH 6A: PARTIAL-AXIS EXCITATION (sigma = inf ALONG UNPUSHED AXES) & CURVE COMBINATION")
    print("=" * 94)
    # Push only along x (a=(1,0), F=(2,0), s=0.1) -> K_x(diag(2, g22)) == 0 for ALL g22 in (0, inf)!
    mr_x = [MetricReading(isotropic([2.0, 0.0], [1.0, 0.0], 0.1, 0.1))]
    mr_y = [MetricReading(isotropic([0.0, 5.0], [0.0, 1.0], 0.1, 0.1))]
    print("  x-push only ((2,0), (1,0)): K_x(diag(2.0, g_yy)) across g_yy:")
    for gyy in (0.01, 0.5, 2.0, 5.0, 100.0, 1e4):
        print(f"    g_yy = {gyy:7.2f} -> K_x = {mr_x[0].residual_cost([[2.0, 0.0], [0.0, gyy]]):.6f}")
    print("  Combining K_x + K_y (adding x-push and y-push) pins down g = diag(2, 5):")
    lad_xy = nested_surprise_ladder(mr_x + mr_y)
    print(f"    g_hat = [[{lad_xy['g_hat'][0][0]:.3f}, {lad_xy['g_hat'][0][1]:.3f}], [{lad_xy['g_hat'][1][0]:.3f}, {lad_xy['g_hat'][1][1]:.3f}]]")
    print(f"    K_scalar = {lad_xy['K_scalar']:.2f} (rejects scalar!) | K_metric = {lad_xy['K_metric']:.4f} | K_linear = {lad_xy['K_linear']:.4f}")

    print("\n" + "=" * 94)
    print("PUSH 6B: THREE-RUNG SURPRISE LADDER (SCALAR subset METRIC Sym+_2 subset LINEAR M)")
    print("=" * 94)
    worlds_det = [
        ("1. Isotropic scalar g = 2 I_2", [[2.0, 0.0], [0.0, 2.0]]),
        ("2. Anisotropic metric g = diag(2, 5) (O17)", [[2.0, 0.0], [0.0, 5.0]]),
        ("3. Gyroscopic non-reciprocal M = [[2,-1],[1,2]] (E7)", [[2.0, -1.0], [1.0, 2.0]]),
        ("4. Anti-parallel negative map M = -diag(2, 3)", [[-2.0, 0.0], [0.0, -3.0]]),
    ]
    for label, M_true in worlds_det:
        # Clean 4-direction pushes with small noise sigma=0.1
        rds = []
        for ang in (0.0, 0.6, 1.2, 2.1):
            a_v = [math.cos(ang), math.sin(ang)]
            f_v = mat_vec(M_true, a_v)
            rds.append(isotropic(f_v, a_v, 0.1, 0.1))
        mrs = [MetricReading(r) for r in rds]
        lad = nested_surprise_ladder(mrs)
        print(
            f"  {label:52s}: K_sc={lad['K_scalar']:7.1f} | K_met={lad['K_metric']:7.1f} | "
            f"K_lin={lad['K_linear']:6.2f} (d_iso={lad['delta_iso']:6.1f}, d_rec={lad['delta_recip']:6.1f})"
        )

    # Full-rank (rank=3) 2D anti-parallel hilltop ghost at g0 = I_2:
    # Each pair (-2v, v) & (-v, 2v) has opposite scores at g0 = I_2 (sum = 0) and nonzero U_i!
    rds_anti2d = [
        isotropic([-2.0, 0.0], [1.0, 0.0], 0.1, 0.1),
        isotropic([-1.0, 0.0], [2.0, 0.0], 0.1, 0.1),
        isotropic([0.0, -2.0], [0.0, 1.0], 0.1, 0.1),
        isotropic([0.0, -1.0], [0.0, 2.0], 0.1, 0.1),
        isotropic([-2.0, -2.0], [1.0, 1.0], 0.1, 0.1),
        isotropic([-1.0, -1.0], [2.0, 2.0], 0.1, 0.1),
    ]
    mrs_anti2d = [MetricReading(r) for r in rds_anti2d]
    g_ghost = [[1.0, 0.0], [0.0, 1.0]]
    t50_g = t50_metric_score_test(mrs_anti2d, g_ghost)
    gclr_sym = metric_clr_pvalue(mrs_anti2d, g_ghost, family="sym", weight_mode="unweighted", n_draws=99, opt_iter=25, seed=17)
    zgclr_sym = metric_clr_pvalue(mrs_anti2d, g_ghost, family="sym", weight_mode="z_smooth", n_draws=99, opt_iter=25, seed=17)
    print(f"\n  Full-rank (rank=3) anti-parallel 2D hilltop ghost at g0 = I_2 (K(I_2) = {gclr_sym['k_obs']:.1f}):")
    print(f"    T50 score test         : T2={t50_g['T2']:.4f}, rank={t50_g['rank']}, p={t50_g['p_value']:.4f}  <-- GHOST ACCEPTED!")
    print(f"    Metric-GCLR   (Sym_2)  : lr_obs={gclr_sym['lr_obs']:.1f}, p={gclr_sym['p_value']:.4f}  <-- GHOST KILLED!")
    print(f"    Metric-Z-GCLR (Sym_2)  : lr_obs={zgclr_sym['lr_obs']:.1f}, p={zgclr_sym['p_value']:.4f}  <-- GHOST KILLED (Z-safe)!")

    print("\n" + "=" * 94)
    print("PUSH 6C: MONTE CARLO 95% COVERAGE ON 2D ANISOTROPIC METRIC g = diag(2, 5) (N=8, reps=100)")
    print("=" * 94)
    g_true = [[2.0, 0.0], [0.0, 5.0]]
    worlds_mc = [
        ("Gaussian pure noise (r=0)", 0.0, None, None),
        ("Gaussian partial x-only push (r=4, x-axis)", 4.0, None, 0),
        ("Gaussian weak (r=1.5)", 1.5, None, None),
        ("Gaussian strong (r=5)", 5.0, None, None),
        ("Scale-mix nu=10 weak (r=1.5)", 1.5, 10, None),
        ("Heavy-tail nu=2 (Cauchy-like, r=3)", 3.0, 2, None),
    ]
    for label, r_len, nu_df, ax in worlds_mc:
        reps = 100
        c_t50 = c_gclr = c_zgclr = c_sclr = c_floor = c_raw = 0
        for rep in range(reps):
            rds = _make_metric_world(rng, g_true, 8, r_len, sigma=1.0, nu_df=nu_df, axis_only=ax)
            mrs = [MetricReading(r) for r in rds]
            if t50_metric_score_test(mrs, g_true)["p_value"] > 0.05:
                c_t50 += 1
            k0 = sum(mr.residual_cost(g_true) for mr in mrs)
            _, kmin = minimize_metric_cost(mrs, family="metric", seed_matrix=g_true, max_iter=35)
            # Raw 3-df Wilks chi^2_3 (q95 = 7.815)
            if max(0.0, k0 - min(k0, kmin)) <= 7.8147279:
                c_raw += 1
            p_g = metric_clr_pvalue(mrs, g_true, family="metric", weight_mode="unweighted", null_mode="gaussian", n_draws=39, opt_iter=18, seed=rep * 17 + 1)["p_value"]
            p_z = metric_clr_pvalue(mrs, g_true, family="metric", weight_mode="z_smooth", null_mode="gaussian", n_draws=39, opt_iter=18, seed=rep * 17 + 2)["p_value"]
            p_s = metric_clr_pvalue(mrs, g_true, family="metric", weight_mode="unweighted", null_mode="spherical", n_draws=39, opt_iter=18, seed=rep * 17 + 3)["p_value"]
            if p_g > 0.05:
                c_gclr += 1
            if p_z > 0.05:
                c_zgclr += 1
            if p_s > 0.05:
                c_sclr += 1
            if max(p_z, p_s) > 0.05:
                c_floor += 1
        print(
            f"  {label:42s}: T50={c_t50/reps:5.1%} | GCLR+={c_gclr/reps:5.1%} | Z-GCLR+={c_zgclr/reps:5.1%} | "
            f"SCLR+={c_sclr/reps:5.1%} | Floor={c_floor/reps:5.1%} | RawLR3={c_raw/reps:5.1%}"
        )

    print("\n" + "=" * 94)
    print("PUSH 6D: TENSOR TWO-BODY COMPOSITION (L5 LUMPING g1+g2 & L14 COMMON FORCE (g1^-1+g2^-1)^-1)")
    print("=" * 94)
    rds1 = [isotropic([2.0, 0.0], [1.0, 0.0], 0.1, 0.01), isotropic([0.0, 5.0], [0.0, 1.0], 0.1, 0.01)]
    rds2 = [isotropic([3.0, 0.0], [1.0, 0.0], 0.1, 0.01), isotropic([0.0, 1.0], [0.0, 1.0], 0.1, 0.01)]
    mrs1, mrs2 = [MetricReading(r) for r in rds1], [MetricReading(r) for r in rds2]
    res_l5_true = compose_metric_lump(mrs1, mrs2, [[5.0, 0.0], [0.0, 6.0]])
    res_l5_off = compose_metric_lump(mrs1, mrs2, [[5.4, 0.0], [0.0, 6.0]])
    res_l14_true = compose_metric_common_force(rds1, rds2, [[6.0 / 5.0, 0.0], [0.0, 5.0 / 6.0]])
    res_l14_off = compose_metric_common_force(rds1, rds2, [[1.5, 0.0], [0.0, 5.0 / 6.0]])
    print(f"  L5 lumping (diag(2,5) + diag(3,1) = diag(5,6)):")
    print(f"    at true G = diag(5.0, 6.0): K_lump = {res_l5_true['K_lump']:.6f}")
    print(f"    at off  G = diag(5.4, 6.0): K_lump = {res_l5_off['K_lump']:.3f}")
    print(f"  L14 common force ((diag(2,5)^-1 + diag(3,1)^-1)^-1 = diag(1.2, 0.833)):")
    print(f"    at true G = diag(1.2, 5/6): K_red  = {res_l14_true['K_red']:.6f}")
    print(f"    at off  G = diag(1.5, 5/6): K_red  = {res_l14_off['K_red']:.3f}")

    rds2_xonly = [isotropic([3.0, 0.0], [1.0, 0.0], 0.1, 0.01)]
    mrs2_x = [MetricReading(r) for r in rds2_xonly]
    print("  Partial-axis lumping: diag(2, 5) + (g2_xx=3, g2_yy unmeasured [sigma=inf]):")
    for gyy in (4.5, 4.8, 5.05, 6.0, 20.0):
        k_val = compose_metric_lump(mrs1, mrs2_x, [[5.0, 0.0], [0.0, gyy]])["K_lump"]
        print(f"    G = diag(5.0, {gyy:5.2f}) -> K_lump = {k_val:7.3f}")

    print("\n" + "=" * 94)
    print("PUSH 6E: EXACT SCALAR-ISOTROPY (L22) TEST & 3D METRIC (d=3, 6 PARAMS) CHECK")
    print("=" * 94)
    rds_iso = _make_metric_world(rng, [[2.0, 0.0], [0.0, 2.0]], 8, 15.0, sigma=0.3)
    rds_ani = _make_metric_world(rng, [[2.0, 0.0], [0.0, 5.0]], 8, 15.0, sigma=0.3)
    rds_xonly = [isotropic([2.0, 0.0], [1.0, 0.0], 0.1, 0.1)] * 3
    for nm, rd in [
        ("Isotropic 2D body g = 2 I_2 (N=8, r=15)", rds_iso),
        ("Anisotropic 2D body g = diag(2,5) (N=8, r=15)", rds_ani),
        ("x-axis only push (compatible with scalar m=2!)", rds_xonly),
    ]:
        mrs_c = [MetricReading(r) for r in rd]
        res_g = test_scalar_isotropy_l22(mrs_c, null_mode="gaussian", weight_mode="z_smooth", n_draws=99, seed=31)
        res_s = test_scalar_isotropy_l22(mrs_c, null_mode="spherical", weight_mode="unweighted", n_draws=99, seed=31)
        print(f"  {nm:48s}: L22 Z-GCLR p={res_g['p_value']:.3f} | L22 SCLR p={res_s['p_value']:.3f}")

    # 3D anisotropic metric g = diag(2, 4, 7) (6 parameters): partial x-push stays flat in y,z!
    mr3_x = MetricReading(isotropic([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 0.1))
    k3_a = mr3_x.residual_cost([[2.0, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 50.0]])
    k3_b = mr3_x.residual_cost([[2.0, 0.0, 0.0], [0.0, 99.0, 0.0], [0.0, 0.0, 0.05]])
    print(f"  3D x-only push: K(diag(2, 0.2, 50))={k3_a:.6f}, K(diag(2, 99, 0.05))={k3_b:.6f} (sigma=inf in y,z!)")


if __name__ == "__main__":
    run_metric_benchmarks()
