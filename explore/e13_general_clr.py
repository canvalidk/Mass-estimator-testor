"""E13 Push 1–3: General-Type Conditional LR, Exact Two-Body Composition, and Door B / Spherical CLR.

Pure Python (stdlib only). Extends explore/e13_unnormalized_evidence.py to:
  1. Push 1 — Full-type General Conditional Likelihood Ratio (GCLR+, GCLR, Z-GCLR+, Z-GCLR):
       Exact finite-sample conditional test given (U_1, ..., U_N) for ARBITRARY positive-definite
       correlated covariances Sigma_i, heterogeneous noise ratios, and declared subspaces N_i.
       Includes BOTH:
         - Unweighted GCLR+ / GCLR(RP^1) (w_i = 1; kills ALL ghosts including rank-1 anti-parallel,
           at the cost of failing both exact-zero Z and sigma -> inf via noise count M = sum k_i),
         - U-weighted Z-GCLR+ / Z-GCLR(RP^1) with w_i = H_i(m0) / (1 + H_i(m0)):
           EXACTLY bit-for-bit Z-invariant AND sigma -> inf invariant, exact finite-sample 95%
           coverage on the full type, kills the multi-direction anti-parallel ghost, and sharpens
           T43 leverage — with the honest theorem-forced trade-off of a rank-1 pinhole ghost where
           ALL readings satisfy U_i(m0) = 0 simultaneously.
  2. Push 2 — Exact Weak-Signal Two-Body Composition (L5 lumping M = m1 + m2 and L14 common force):
       Compares raw infimal convolution against exact conditional profiled-lump LR
       (unweighted and Z-weighted), lifting weak-signal coverage from ~89% to >= 95% while
       preserving 1-df chi^2_1 sharpness at strong signal, exact swap duality for L14, and exact
       [4.02, inf) when gluing 5 +- 0.5 kg onto stationary 0/0.
  3. Push 3 — Door B Stress Test, Scale-Free Spherical CLR (SCLR+), and Floor Hybrid (O7 floor row):
       - Exposes naive unnormalized Student-t (Wilks 3.841) Neyman-Scott undercoverage at weak signal
         (even on the exact nu=10 Student-t world it was built for: 62%).
       - Implements SCLR+ (conditioning on both U_i and ||V_i||, drawing V_i^(b)/||V_i|| uniform on
         S^{k_i-1}, or exact 2^N sign enumeration for small-N 1D): exact finite-sample >= 95%
         coverage under ANY per-reading scale mixture (Gaussian, nu=10, Cauchy-like nu=2), exact
         bit-for-bit Z invariance, and automatic discordance widening, while honestly exhibiting
         its 1D N <= 4 floor (1/2^N > 0.05) and scale-blindness on over-agreeing data.
       - Tests the zero-tuning O7 floor-row union C_floor = C_95(Z-GCLR+) U C_95(SCLR+): guaranteed
         >= 95% coverage under all v_i > 0, never narrower than declared sigma on over-agreeing data,
         and automatically widens on discordant data (finite from N >= 2 in 2D/3D, N >= 5 in 1D).

Run from the repo root:
    python3 explore/e13_general_clr.py
"""
import math
import random

from explore.e13_unnormalized_evidence import (
    CHI2_95_1,
    Reading,
    cholesky,
    compose_lump,
    dot,
    evaluate_series,
    format_intervals,
    householder_qr_full,
    isotropic,
    mat_t,
    solve_lower,
    solve_lower_mat,
)


# ---------------------------------------------------------------------------
# Fast within-span W_i representation for arbitrary correlated Sigma_i and N_i
# ---------------------------------------------------------------------------

class SpanReading:
    """Precomputed 2k-dimensional within-span representation of a general Reading."""

    def __init__(self, reading: Reading):
        self.reading = reading
        self.d = reading.d
        self.k = reading.k
        n = reading.subspace
        z = [[0.0] * self.k for _ in range(self.d)]
        span = [n[i] + z[i] for i in range(self.d)] + [z[i] + n[i] for i in range(self.d)]
        span_w = solve_lower_mat(reading.chol, span)  # 2d x 2k
        Q_full, _ = householder_qr_full(span_w)
        Qt = mat_t(Q_full)
        self.E_W_T = Qt[: 2 * self.k]  # 2k x 2d
        self.z = tuple(dot(row, reading.whitened) for row in self.E_W_T)
        bw_F = solve_lower_mat(
            reading.chol,
            [n[i] for i in range(self.d)] + [[0.0] * self.k for _ in range(self.d)],
        )
        bw_a = solve_lower_mat(
            reading.chol,
            [[0.0] * self.k for _ in range(self.d)] + [n[i] for i in range(self.d)],
        )
        self.BF = [
            [dot(self.E_W_T[r], [bw_F[p][c] for p in range(2 * self.d)]) for c in range(self.k)]
            for r in range(2 * self.k)
        ]
        self.Ba = [
            [dot(self.E_W_T[r], [bw_a[p][c] for p in range(2 * self.d)]) for c in range(self.k)]
            for r in range(2 * self.k)
        ]

    def split_at_phi(self, phi):
        """Return (Q_par_T, Q_perp_T) in W_i coordinates, each k x 2k with orthonormal rows."""
        s, co = math.sin(phi), math.cos(phi)
        m2k, k = 2 * self.k, self.k
        B = [[s * self.BF[r][c] + co * self.Ba[r][c] for c in range(k)] for r in range(m2k)]
        Q, _ = householder_qr_full(B)
        Qt = mat_t(Q)
        return Qt[:k], Qt[k:]

    def uv_at_phi(self, phi):
        """Return (Q_par_T, Q_perp_T, U, V, H, K) at angle phi."""
        q_par, q_perp = self.split_at_phi(phi)
        u = [dot(row, self.z) for row in q_par]
        v = [dot(row, self.z) for row in q_perp]
        h = sum(x * x for x in u)
        k_val = sum(x * x for x in v)
        return q_par, q_perp, u, v, h, k_val


def swap_reading(reading: Reading) -> Reading:
    """Swap force and acceleration channels (and covariance blocks) for L14 common-force duality."""
    d = reading.d
    idx = list(range(d, 2 * d)) + list(range(d))
    cov_s = [[reading.cov[idx[i]][idx[j]] for j in range(2 * d)] for i in range(2 * d)]
    return Reading(reading.acceleration, reading.force, cov_s, subspace=reading.subspace)


def _weight_fn(h, mode):
    if mode == "unweighted":
        return 1.0
    if mode == "z_smooth":
        return h / (1.0 + h)
    if mode == "z_quad":
        h2 = h * h
        return h2 / (1.0 + h2)
    if mode == "z_cap10":
        h2 = h * h
        return h2 / (100.0 + h2)
    if mode == "z_linear":
        return h
    if mode == "z_pow2":
        return h * h
    if mode == "z_exp":
        return 1.0 - math.exp(-0.5 * h)
    raise ValueError(f"Unknown weight mode {mode!r}")


def _k_of_z(q_perp_T, z_vec):
    """Squared norm ||Q_perp^T z||^2, unrolled for k = 1, 2, 3."""
    k = len(q_perp_T)
    if k == 1:
        r0 = q_perp_T[0]
        v0 = r0[0] * z_vec[0] + r0[1] * z_vec[1]
        return v0 * v0
    if k == 2:
        r0, r1 = q_perp_T
        z0, z1, z2, z3 = z_vec
        v0 = r0[0] * z0 + r0[1] * z1 + r0[2] * z2 + r0[3] * z3
        v1 = r1[0] * z0 + r1[1] * z1 + r1[2] * z2 + r1[3] * z3
        return v0 * v0 + v1 * v1
    if k == 3:
        r0, r1, r2 = q_perp_T
        z0, z1, z2, z3, z4, z5 = z_vec
        v0 = r0[0] * z0 + r0[1] * z1 + r0[2] * z2 + r0[3] * z3 + r0[4] * z4 + r0[5] * z5
        v1 = r1[0] * z0 + r1[1] * z1 + r1[2] * z2 + r1[3] * z3 + r1[4] * z4 + r1[5] * z5
        v2 = r2[0] * z0 + r2[1] * z1 + r2[2] * z2 + r2[3] * z3 + r2[4] * z4 + r2[5] * z5
        return v0 * v0 + v1 * v1 + v2 * v2
    total = 0.0
    for row in q_perp_T:
        val = sum(a * b for a, b in zip(row, z_vec))
        total += val * val
    return total


def _reconstruct_z(q_par_T, q_perp_T, u_vec, v_vec):
    """Compute z = Q_par U + Q_perp V (length 2k)."""
    k = len(u_vec)
    if k == 1:
        u0, v0 = u_vec[0], v_vec[0]
        p0, n0 = q_par_T[0], q_perp_T[0]
        return (u0 * p0[0] + v0 * n0[0], u0 * p0[1] + v0 * n0[1])
    m2k = 2 * k
    out = [0.0] * m2k
    for u_val, row in zip(u_vec, q_par_T):
        if u_val != 0.0:
            for j in range(m2k):
                out[j] += u_val * row[j]
    for v_val, row in zip(v_vec, q_perp_T):
        if v_val != 0.0:
            for j in range(m2k):
                out[j] += v_val * row[j]
    return tuple(out)


def _rand_unit_sphere(k, rng):
    """Exact uniform draw on S^{k-1} (for k=1: exact uniform +-1)."""
    if k == 1:
        return (1.0,) if rng.getrandbits(1) else (-1.0,)
    while True:
        v = [rng.gauss(0.0, 1.0) for _ in range(k)]
        n2 = sum(x * x for x in v)
        if n2 > 1e-30:
            inv = 1.0 / math.sqrt(n2)
            return tuple(x * inv for x in v)


def precompute_grid(span_readings, *, domain="positive", n_grid=41, scale_hint=1.0):
    if domain == "positive":
        us = [0.5 * math.pi * j / (n_grid - 1) for j in range(n_grid)]
        grid_phis = [
            math.atan(scale_hint * math.tan(u)) if j < n_grid - 1 else 0.5 * math.pi
            for j, u in enumerate(us)
        ]
    else:
        us = [-0.5 * math.pi + math.pi * j / n_grid for j in range(n_grid)]
        grid_phis = [
            math.atan(scale_hint * math.tan(u)) if abs(math.cos(u)) > 1e-12 else 0.5 * math.pi
            for u in us
        ]
    per_q_perp = [
        [tuple(tuple(r) for r in sr.split_at_phi(ph)[1]) for ph in grid_phis]
        for sr in span_readings
    ]
    return grid_phis, per_q_perp


def general_clr_pvalue(
    span_readings,
    phi0,
    *,
    domain="positive",
    weight_mode="unweighted",
    null_mode="gaussian",
    nu_student=None,
    n_grid=41,
    n_draws=199,
    seed=12345,
    scale_hint=1.0,
    precomputed_grid=None,
    exact_1d_spherical_max_n=8,
):
    """Exact finite-sample conditional p-value for candidate angle phi0 on the full input type.

    Args:
        span_readings: list of SpanReading.
        phi0: candidate angle in (0, pi/2) (m0 = tan(phi0)).
        domain: 'positive' (inf over m > 0: GCLR+) or 'whole' (inf over RP^1: GCLR).
        weight_mode: 'unweighted' (w_i = 1) or 'z_smooth' (w_i = H_i / (1 + H_i), exact Z!).
        null_mode:
            'gaussian'  -> draw V_i^(b) ~ N(0, I_{k_i}) conditional on U_i (exact under Gaussian).
            'spherical' -> draw V_i^(b) = ||V_i|| * Unif(S^{k_i-1}) conditional on (U_i, ||V_i||)
                           (exact under ANY per-reading scale mixture, exact Z!).
        nu_student: if float (Door B), use Student-t cost (nu + k_i) log(1 + K_i / nu).
        n_grid: number of grid angles for the infimum (exact exchangeability holds for ANY n_grid!).
        n_draws: number of conditional null replicates B (choose B+1 multiple of 20 for exact 0.05).
        exact_1d_spherical_max_n: for 1D spherical null with <= this many active nonzero-V readings,
            enumerate all 2^N sign configurations deterministically (exact, zero MC variance).

    Returns dict(p_value, lr_obs, k_obs, weights, degenerate).
    """
    if not span_readings:
        return dict(p_value=1.0, lr_obs=0.0, k_obs=0.0, weights=[], degenerate=True)

    if precomputed_grid is None:
        grid_phis, per_q_perp = precompute_grid(
            span_readings, domain=domain, n_grid=n_grid, scale_hint=scale_hint
        )
    else:
        grid_phis, per_q_perp = precomputed_grid

    splits = []
    weights = []
    for idx, sr in enumerate(span_readings):
        q_par, q_perp, u, v, h, k_val = sr.uv_at_phi(phi0)
        w = _weight_fn(h, weight_mode)
        weights.append(w)
        if w > 0.0:
            v_norm = math.sqrt(k_val)
            # In spherical mode, a reading with both U=0 and V=0 contributes identically 0
            # (exact bit-for-bit Z invariance without consuming RNG or sign bits).
            if null_mode == "spherical" and h == 0.0 and v_norm == 0.0:
                continue
            m2k = 2 * sr.k
            z_u = [0.0] * m2k
            for u_val, row in zip(u, q_par):
                if u_val != 0.0:
                    for j in range(m2k):
                        z_u[j] += u_val * row[j]
            splits.append((idx, sr.k, tuple(z_u), tuple(tuple(r) for r in q_perp), v_norm, w))

    if not splits:
        return dict(p_value=1.0, lr_obs=0.0, k_obs=0.0, weights=weights, degenerate=True)

    n_g = len(grid_phis)
    active_q = [[per_q_perp[idx][g] for g in range(n_g)] for idx, _, _, _, _, _ in splits]
    log1p = math.log1p

    if nu_student is None:
        def compute_lr(z_list, k0_list):
            c0 = 0.0
            for (_, _, _, _, _, w), k0 in zip(splits, k0_list):
                c0 += w * k0
            c_min = c0
            for g in range(n_g):
                cg = 0.0
                for s_i, (_, _, _, _, _, w) in enumerate(splits):
                    cg += w * _k_of_z(active_q[s_i][g], z_list[s_i])
                    if cg >= c_min:
                        break
                if cg < c_min:
                    c_min = cg
            return max(0.0, c0 - c_min), c0
    else:
        nu_val = float(nu_student)
        inv_nu = 1.0 / nu_val

        def compute_lr(z_list, k0_list):
            c0 = 0.0
            for (_, k_dim, _, _, _, w), k0 in zip(splits, k0_list):
                c0 += w * (nu_val + k_dim) * log1p(k0 * inv_nu)
            c_min = c0
            for g in range(n_g):
                cg = 0.0
                for s_i, (_, k_dim, _, _, _, w) in enumerate(splits):
                    kg = _k_of_z(active_q[s_i][g], z_list[s_i])
                    cg += w * (nu_val + k_dim) * log1p(kg * inv_nu)
                    if cg >= c_min:
                        break
                if cg < c_min:
                    c_min = cg
            return max(0.0, c0 - c_min), c0

    z_obs = [span_readings[idx].z for idx, _, _, _, _, _ in splits]
    k0_obs = [v_norm * v_norm for _, _, _, _, v_norm, _ in splits]
    lr_obs, k_obs = compute_lr(z_obs, k0_obs)

    if lr_obs <= 1e-15:
        return dict(p_value=1.0, lr_obs=0.0, k_obs=k_obs, weights=weights, degenerate=False)

    # Relative tolerance prevents floating-point roundoff from breaking exact ties at huge lr_obs.
    tol = 1e-10 * max(1.0, lr_obs)

    # Exact deterministic enumeration for small-N 1D spherical tests (eliminates MC jitter)
    if null_mode == "spherical" and all(k_dim == 1 for _, k_dim, _, _, _, _ in splits):
        active_flip_indices = [i for i, (_, _, _, _, vn, _) in enumerate(splits) if vn > 1e-15]
        n_flips = len(active_flip_indices)
        if n_flips == 0:
            return dict(p_value=1.0, lr_obs=lr_obs, k_obs=k_obs, weights=weights, degenerate=False)
        if n_flips <= exact_1d_spherical_max_n:
            n_total = 1 << n_flips
            ge_count = 0
            # Pre-extract plus/minus z vectors for each reading
            z_plus = []
            z_minus = []
            for _, _, z_u, q_perp, vn, _ in splits:
                n0 = q_perp[0]
                z_plus.append((z_u[0] + vn * n0[0], z_u[1] + vn * n0[1]))
                z_minus.append((z_u[0] - vn * n0[0], z_u[1] - vn * n0[1]))
            flip_pos = {s_i: bit for bit, s_i in enumerate(active_flip_indices)}
            n_sp = len(splits)
            for mask in range(n_total):
                z_draw = [
                    z_plus[s_i] if (s_i not in flip_pos or ((mask >> flip_pos[s_i]) & 1)) else z_minus[s_i]
                    for s_i in range(n_sp)
                ]
                lr_b, _ = compute_lr(z_draw, k0_obs)
                if lr_b >= lr_obs - tol:
                    ge_count += 1
            return dict(p_value=ge_count / n_total, lr_obs=lr_obs, k_obs=k_obs, weights=weights, degenerate=False)

    rng = random.Random(seed)
    gauss = rng.gauss
    ge_count = 0
    for _ in range(n_draws):
        z_draw = []
        k0_draw = []
        for _, k_dim, z_u, q_perp, v_norm, _ in splits:
            if null_mode == "gaussian":
                if k_dim == 1:
                    v0 = gauss(0.0, 1.0)
                    n0 = q_perp[0]
                    z_draw.append((z_u[0] + v0 * n0[0], z_u[1] + v0 * n0[1]))
                    k0_draw.append(v0 * v0)
                else:
                    v_b = [gauss(0.0, 1.0) for _ in range(k_dim)]
                    m2k = 2 * k_dim
                    z_vec = list(z_u)
                    k0_b = 0.0
                    for v_val, row in zip(v_b, q_perp):
                        k0_b += v_val * v_val
                        for j in range(m2k):
                            z_vec[j] += v_val * row[j]
                    z_draw.append(tuple(z_vec))
                    k0_draw.append(k0_b)
            elif null_mode == "spherical":
                if v_norm <= 1e-15:
                    z_draw.append(z_u)
                    k0_draw.append(0.0)
                elif k_dim == 1:
                    v0 = v_norm if rng.getrandbits(1) else -v_norm
                    n0 = q_perp[0]
                    z_draw.append((z_u[0] + v0 * n0[0], z_u[1] + v0 * n0[1]))
                    k0_draw.append(v_norm * v_norm)
                else:
                    dir_b = _rand_unit_sphere(k_dim, rng)
                    m2k = 2 * k_dim
                    z_vec = list(z_u)
                    for d_val, row in zip(dir_b, q_perp):
                        v_val = v_norm * d_val
                        for j in range(m2k):
                            z_vec[j] += v_val * row[j]
                    z_draw.append(tuple(z_vec))
                    k0_draw.append(v_norm * v_norm)
            else:
                raise ValueError(f"Unknown null_mode {null_mode!r}")
        lr_b, _ = compute_lr(z_draw, k0_draw)
        if lr_b >= lr_obs - tol:
            ge_count += 1

    p_val = (1 + ge_count) / (n_draws + 1)
    return dict(p_value=p_val, lr_obs=lr_obs, k_obs=k_obs, weights=weights, degenerate=False)


def _refine_phi_min(readings, phi_lo, phi_hi):
    """Golden-section refinement of K(phi) on [phi_lo, phi_hi] for scan centering."""
    gr = 0.5 * (math.sqrt(5.0) - 1.0)
    a, b = phi_lo, phi_hi
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = evaluate_series(readings, c)["K"]
    fd = evaluate_series(readings, d)["K"]
    for _ in range(28):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = evaluate_series(readings, c)["K"]
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = evaluate_series(readings, d)["K"]
    return 0.5 * (a + b)


def gclr_confidence_set(
    readings,
    *,
    domain="positive",
    weight_mode="unweighted",
    null_mode="gaussian",
    nu_student=None,
    floor_union=False,
    n_scan=64,
    n_grid=91,
    n_draws=199,
    seed=2026,
    scale_hint=1.0,
):
    """Scan (0, inf) with accurate peak zoom, returning 95% compatibility set.

    If floor_union=True, returns the O7 floor-row union C_95(Z-GCLR+) U C_95(SCLR+),
    which is guaranteed >= 95% under ANY per-reading scale mixture v_i > 0, never narrower
    than declared sigma on over-agreeing data, and automatically widens on discordance.
    """
    srs = [SpanReading(r) for r in readings]
    pre = precompute_grid(srs, domain=domain, n_grid=n_grid, scale_hint=scale_hint)
    eps = 1e-4
    us = [eps + (0.5 * math.pi - 2 * eps) * j / (n_scan - 1) for j in range(n_scan)]
    phis = [math.atan(scale_hint * math.tan(u)) for u in us]
    ks = [evaluate_series(readings, ph)["K"] for ph in phis]
    best_i = ks.index(min(ks))
    lo = phis[max(0, best_i - 1)]
    hi = phis[min(len(phis) - 1, best_i + 1)]
    phi_best = _refine_phi_min(readings, lo, hi)
    extra = [phi_best]
    for hw in (0.04, 0.008, 0.0012):
        for j in range(-12, 13):
            p = phi_best + hw * j / 12.0
            if eps < p < 0.5 * math.pi - eps:
                extra.append(p)
    phis = sorted(set(phis + extra))
    mask = []
    for idx, ph in enumerate(phis):
        if floor_union:
            p_g = general_clr_pvalue(
                srs, ph, domain=domain, weight_mode="z_smooth", null_mode="gaussian",
                n_grid=n_grid, n_draws=n_draws, seed=seed + idx * 17,
                scale_hint=scale_hint, precomputed_grid=pre,
            )["p_value"]
            p_s = general_clr_pvalue(
                srs, ph, domain=domain, weight_mode="unweighted", null_mode="spherical",
                nu_student=nu_student, n_grid=n_grid, n_draws=n_draws, seed=seed + idx * 17 + 1,
                scale_hint=scale_hint, precomputed_grid=pre,
            )["p_value"]
            mask.append(max(p_g, p_s) > 0.05)
        else:
            res = general_clr_pvalue(
                srs, ph, domain=domain, weight_mode=weight_mode, null_mode=null_mode,
                nu_student=nu_student, n_grid=n_grid, n_draws=n_draws, seed=seed + idx * 17,
                scale_hint=scale_hint, precomputed_grid=pre,
            )
            mask.append(res["p_value"] > 0.05)
    return format_intervals(phis, mask)


# ---------------------------------------------------------------------------
# Push 2: Exact Conditional Two-Body Composition (L5 lumping & L14 common force)
# ---------------------------------------------------------------------------

def two_body_lump_conditional_pvalue(
    sr1,
    sr2,
    M0,
    m1_candidate,
    *,
    weight_mode="unweighted",
    n_draws=199,
    n_grid=41,
    n_split=45,
    seed=777,
):
    """Exact conditional p-value for the profiled lump LR at split (m1, M0 - m1).

    Statistic: LR_lump(M0; z1, z2) = inf_{0 < m1 < M0} [K1_w(m1) + K2_w(M0 - m1)] - K1_w,min - K2_w,min.
    Conditional on (U1(m1_cand), U2(M0 - m1_cand)), at the true split (m1*, m2*),
    V1 ~ N(0, I) and V2 ~ N(0, I) are exact independent standard Gaussians and all U-based
    weights w_i = H_i / (1 + H_i) are fixed constants. Taking sup_{0 < m1 < M0} over m1 yields
    a guaranteed >= 95% exact confidence set for M0.
    With weight_mode='z_smooth', an unmeasured stationary body (0/0) drops out bit-for-bit!
    """
    m2_cand = M0 - m1_candidate
    if m1_candidate <= 0.0 or m2_cand <= 0.0:
        return 0.0

    # Active readings (with positive weight at candidate split)
    sp1 = []
    for sr in sr1:
        q_par, q_perp, u, v, h, k_val = sr.uv_at_phi(math.atan(m1_candidate))
        w = _weight_fn(h, weight_mode)
        if w > 0.0:
            sp1.append((sr, q_par, q_perp, u, w))
    sp2 = []
    for sr in sr2:
        q_par, q_perp, u, v, h, k_val = sr.uv_at_phi(math.atan(m2_cand))
        w = _weight_fn(h, weight_mode)
        if w > 0.0:
            sp2.append((sr, q_par, q_perp, u, w))

    if not sp1 and not sp2:
        return 1.0

    # Deterministic split grid anchored at m1_candidate (fixed across obs and all null draws)
    base_splits = [M0 * (j + 0.5) / n_split for j in range(n_split)]
    local_splits = [
        m1_candidate + M0 * frac
        for frac in (-0.06, -0.03, -0.015, -0.005, 0.0, 0.005, 0.015, 0.03, 0.06)
        if 1e-6 * M0 < m1_candidate + M0 * frac < (1.0 - 1e-6) * M0
    ]
    split_m1s = sorted(set(base_splits + local_splits))

    sr1_act = [t[0] for t in sp1]
    sr2_act = [t[0] for t in sp2]
    w1_act = [t[4] for t in sp1]
    w2_act = [t[4] for t in sp2]

    q_split_1 = [[tuple(tuple(r) for r in sr.split_at_phi(math.atan(m))[1]) for m in split_m1s] for sr in sr1_act]
    q_split_2 = [[tuple(tuple(r) for r in sr.split_at_phi(math.atan(M0 - m))[1]) for m in split_m1s] for sr in sr2_act]

    pre1 = precompute_grid(sr1_act, domain="positive", n_grid=n_grid, scale_hint=m1_candidate)[1] if sr1_act else []
    pre2 = precompute_grid(sr2_act, domain="positive", n_grid=n_grid, scale_hint=m2_cand)[1] if sr2_act else []
    n_sp = len(split_m1s)

    def lr_lump(z1_list, z2_list):
        if sr1_act:
            k1_min = min(
                sum(w1_act[i] * _k_of_z(pre1[i][g], z1_list[i]) for i in range(len(sr1_act)))
                for g in range(n_grid)
            )
        else:
            k1_min = 0.0
        if sr2_act:
            k2_min = min(
                sum(w2_act[i] * _k_of_z(pre2[i][g], z2_list[i]) for i in range(len(sr2_act)))
                for g in range(n_grid)
            )
        else:
            k2_min = 0.0
        best_sum = float("inf")
        for s_idx in range(n_sp):
            v1 = sum(w1_act[i] * _k_of_z(q_split_1[i][s_idx], z1_list[i]) for i in range(len(sr1_act)))
            v2 = sum(w2_act[i] * _k_of_z(q_split_2[i][s_idx], z2_list[i]) for i in range(len(sr2_act)))
            if v1 < k1_min:
                k1_min = v1
            if v2 < k2_min:
                k2_min = v2
            tot = v1 + v2
            if tot < best_sum:
                best_sum = tot
        return max(0.0, best_sum - k1_min - k2_min)

    obs = lr_lump([sr.z for sr in sr1_act], [sr.z for sr in sr2_act])
    if obs <= 1e-15:
        return 1.0
    tol = 1e-10 * max(1.0, obs)
    rng = random.Random(seed)
    ge = 0
    for _ in range(n_draws):
        z1_b = [
            _reconstruct_z(qp, qn, u, [rng.gauss(0.0, 1.0) for _ in range(sr.k)])
            for sr, qp, qn, u, _ in sp1
        ]
        z2_b = [
            _reconstruct_z(qp, qn, u, [rng.gauss(0.0, 1.0) for _ in range(sr.k)])
            for sr, qp, qn, u, _ in sp2
        ]
        if lr_lump(z1_b, z2_b) >= obs - tol:
            ge += 1
    return (1 + ge) / (n_draws + 1)


def two_body_common_force_conditional_pvalue(
    readings1,
    readings2,
    M_red_0,
    m1_candidate,
    *,
    weight_mode="unweighted",
    n_draws=199,
    n_grid=41,
    n_split=45,
    seed=777,
):
    """Exact conditional p-value for L14 common force (1/M = 1/m1 + 1/m2) via exact swap duality."""
    if M_red_0 <= 0.0 or m1_candidate <= M_red_0:
        return 0.0
    sr1_s = [SpanReading(swap_reading(r)) for r in readings1]
    sr2_s = [SpanReading(swap_reading(r)) for r in readings2]
    return two_body_lump_conditional_pvalue(
        sr1_s,
        sr2_s,
        1.0 / M_red_0,
        1.0 / m1_candidate,
        weight_mode=weight_mode,
        n_draws=n_draws,
        n_grid=n_grid,
        n_split=n_split,
        seed=seed,
    )


# ---------------------------------------------------------------------------
# Benchmark suite for Push 1, Push 2, and Push 3
# ---------------------------------------------------------------------------

def _make_correlated_world(rng, N, m_true, r_len, d=3, k=2, rho=0.6, ratio_spread=True):
    """Generate full-type readings: correlated cross-channel Sigma_i, heterogeneous ratios, subspace k <= d."""
    rds = []
    for i in range(N):
        sf = 1.0 + (0.8 * (i % 3) if ratio_spread else 0.0)
        sa = 1.0 + (0.5 * ((i + 1) % 2) if ratio_spread else 0.0)
        cov = [[0.0] * (2 * d) for _ in range(2 * d)]
        for j in range(d):
            cov[j][j] = sf * sf * (1.0 + 0.3 * j)
            cov[d + j][d + j] = sa * sa * (1.0 + 0.2 * j)
            cov[j][d + j] = rho * math.sqrt(cov[j][j] * cov[d + j][d + j])
            cov[d + j][j] = cov[j][d + j]
        if k == d:
            sub = [[1.0 if r == c else 0.0 for c in range(d)] for r in range(d)]
        else:
            raw = [[rng.gauss(0.0, 1.0) for _ in range(k)] for _ in range(d)]
            Q_sub, _ = householder_qr_full(raw)
            sub = [[Q_sub[r][c] for c in range(k)] for r in range(d)]
        L = cholesky(cov)
        u_dir = _rand_unit_sphere(k, rng)
        a_unit = [sum(sub[r][c] * u_dir[c] for c in range(k)) for r in range(d)]
        f_unit = [m_true * x for x in a_unit]
        w_unit = solve_lower(L, f_unit + a_unit)
        w_norm = math.sqrt(sum(x * x for x in w_unit))
        scale = r_len / w_norm if w_norm > 0 else 0.0
        mu = [scale * x for x in (f_unit + a_unit)]
        eps_w = [rng.gauss(0.0, 1.0) for _ in range(2 * d)]
        noise = [sum(L[r][c] * eps_w[c] for c in range(2 * d)) for r in range(2 * d)]
        y = [m_val + n_val for m_val, n_val in zip(mu, noise)]
        rds.append(Reading(y[:d], y[d:], cov, subspace=sub))
    return rds


def _make_scale_mixture_world(rng, d, N, r_len, m_true, nu_df=None):
    """Generate d-dim readings with Gaussian (nu_df=None) or IG(nu/2, nu/2) noise scales."""
    th = math.atan(m_true)
    s, co = math.sin(th), math.cos(th)
    rds = []
    for _ in range(N):
        if nu_df is None:
            sd = 1.0
        else:
            chi2_draw = sum(rng.gauss(0.0, 1.0) ** 2 for _ in range(nu_df))
            sd = math.sqrt(nu_df / chi2_draw)
        u = _rand_unit_sphere(d, rng)
        F = [r_len * s * x + sd * rng.gauss(0.0, 1.0) for x in u]
        a = [r_len * co * x + sd * rng.gauss(0.0, 1.0) for x in u]
        rds.append(isotropic(F, a, 1.0, 1.0))
    return rds


def run_push_benchmarks(seed=20261121):
    rng = random.Random(seed)

    print("=" * 92)
    print("PUSH 1A: DETERMINISTIC SCENES — UNWEIGHTED GCLR+ vs Z-SAFE U-WEIGHTED Z-GCLR+")
    print("=" * 92)
    scenes = [
        ("All-zero 0/0 (3D)", [isotropic([0, 0, 0], [0, 0, 0], 2.0, 1.0)], 2.0),
        ("Scene A: near (5.5, 5.6) +- 0.1", [isotropic(5.5, 1.0, 0.1, 0.002), isotropic(5.6, 1.0, 0.1, 0.002)], 5.5),
        ("Scene B: far (5.0, 6.0) +- 0.1 (10 sigma)", [isotropic(5.0, 1.0, 0.1, 0.002), isotropic(6.0, 1.0, 0.1, 0.002)], 5.5),
        ("Right-angle Scene C (G29: F=1,a=0 & F=0,a=1, s=0.05)", [isotropic(1.0, 0.0, 0.05, 0.05), isotropic(0.0, 1.0, 0.05, 0.05)], 1.0),
        ("Anti-parallel pair ((2,-1) & (1,-2), s=0.1)", [isotropic(2.0, -1.0, 0.1, 0.1), isotropic(1.0, -2.0, 0.1, 0.1)], 1.0),
    ]
    for name, rds, sh in scenes:
        s_gclr_p = gclr_confidence_set(rds, domain="positive", weight_mode="unweighted", scale_hint=sh, n_grid=121, n_draws=199)
        s_gclr_w = gclr_confidence_set(rds, domain="whole", weight_mode="unweighted", scale_hint=sh, n_grid=121, n_draws=199)
        s_zgclr_p = gclr_confidence_set(rds, domain="positive", weight_mode="z_smooth", scale_hint=sh, n_grid=121, n_draws=199)
        s_zgclr_w = gclr_confidence_set(rds, domain="whole", weight_mode="z_smooth", scale_hint=sh, n_grid=121, n_draws=199)
        print(f"\n{name}")
        print(f"  GCLR+   (unweighted, m>0)  : {s_gclr_p}")
        print(f"  GCLR    (unweighted, RP^1) : {s_gclr_w}")
        print(f"  Z-GCLR+ (Z-safe w,   m>0)  : {s_zgclr_p}")
        print(f"  Z-GCLR  (Z-safe w,   RP^1) : {s_zgclr_w}")

    print("\n" + "=" * 92)
    print("PUSH 1B: FULL-TYPE CORRELATED + HETEROGENEOUS + SUBSPACE COVERAGE (reps=200)")
    print("  d=3, subspace k=2, cross-channel rho=0.6, mixed noise ratios, N=10, true m=2.5")
    print("=" * 92)
    for label, r_len in [("Pure noise (r=0)", 0.0), ("Weak signal (r=1)", 1.0), ("Strong signal (r=5)", 5.0)]:
        reps = 200
        c_t45 = c_gp = c_gw = c_zgp = c_wilks = 0
        phi0 = math.atan(2.5)
        for rep in range(reps):
            rds = _make_correlated_world(rng, 10, 2.5, r_len, d=3, k=2, rho=0.6, ratio_spread=True)
            e0 = evaluate_series(rds, phi0)
            if e0["T2"] <= CHI2_95_1:
                c_t45 += 1
            srs = [SpanReading(r) for r in rds]
            pre_p = precompute_grid(srs, domain="positive", n_grid=31, scale_hint=2.5)
            k_min_g = min(
                sum(_k_of_z(pre_p[1][i][g], srs[i].z) for i in range(len(srs)))
                for g in range(len(pre_p[0]))
            )
            if (e0["K"] - min(k_min_g, e0["K"])) <= CHI2_95_1:
                c_wilks += 1
            if general_clr_pvalue(srs, phi0, domain="positive", weight_mode="unweighted", n_draws=99, seed=rep * 31 + 1, precomputed_grid=pre_p)["p_value"] > 0.05:
                c_gp += 1
            if general_clr_pvalue(srs, phi0, domain="whole", weight_mode="unweighted", n_draws=99, seed=rep * 31 + 2, n_grid=31, scale_hint=2.5)["p_value"] > 0.05:
                c_gw += 1
            if general_clr_pvalue(srs, phi0, domain="positive", weight_mode="z_smooth", n_draws=99, seed=rep * 31 + 3, precomputed_grid=pre_p)["p_value"] > 0.05:
                c_zgp += 1
        print(
            f"  {label:22s}: T45={c_t45/reps:5.1%} | GCLR+={c_gp/reps:5.1%} | "
            f"GCLR(RP1)={c_gw/reps:5.1%} | Z-GCLR+={c_zgp/reps:5.1%} | RawLR={c_wilks/reps:5.1%}"
        )

    print("\n" + "=" * 92)
    print("PUSH 2: TWO-BODY COMPOSITION COVERAGE (L5 lumping M = 2 + 3 = 5, 3D N=5 each, reps=150)")
    print("=" * 92)
    for label, r_len in [("Weak (signal=1)", 1.0), ("Moderate (signal=3)", 3.0), ("Strong (signal=6)", 6.0)]:
        reps = 150
        c_raw = c_cond_u = c_cond_z = 0
        for rep in range(reps):
            th1, th2 = math.atan(2.0), math.atan(3.0)
            r1 = [
                isotropic(
                    [r_len * math.sin(th1) * x + rng.gauss(0, 1) for x in u],
                    [r_len * math.cos(th1) * x + rng.gauss(0, 1) for x in u],
                    1.0, 1.0,
                )
                for u in (_rand_unit_sphere(3, rng) for _ in range(5))
            ]
            r2 = [
                isotropic(
                    [r_len * math.sin(th2) * x + rng.gauss(0, 1) for x in u],
                    [r_len * math.cos(th2) * x + rng.gauss(0, 1) for x in u],
                    1.0, 1.0,
                )
                for u in (_rand_unit_sphere(3, rng) for _ in range(5))
            ]
            s1_ff = sum(dot(r.force, r.force) for r in r1)
            s1_aa = sum(dot(r.acceleration, r.acceleration) for r in r1)
            s1_fa = sum(dot(r.force, r.acceleration) for r in r1)
            s2_ff = sum(dot(r.force, r.force) for r in r2)
            s2_aa = sum(dot(r.acceleration, r.acceleration) for r in r2)
            s2_fa = sum(dot(r.force, r.acceleration) for r in r2)
            K1 = lambda m: (s1_ff - 2.0 * m * s1_fa + m * m * s1_aa) / (1.0 + m * m)
            K2 = lambda m: (s2_ff - 2.0 * m * s2_fa + m * m * s2_aa) / (1.0 + m * m)
            k1_min = 0.5 * ((s1_ff + s1_aa) - math.hypot(s1_ff - s1_aa, 2.0 * s1_fa)) if s1_fa >= 0 else min(s1_ff, s1_aa)
            k2_min = 0.5 * ((s2_ff + s2_aa) - math.hypot(s2_ff - s2_aa, 2.0 * s2_fa)) if s2_fa >= 0 else min(s2_ff, s2_aa)
            if (compose_lump(K1, K2, [5.0])[0] - k1_min - k2_min) <= CHI2_95_1:
                c_raw += 1
            sr1 = [SpanReading(r) for r in r1]
            sr2 = [SpanReading(r) for r in r2]
            if two_body_lump_conditional_pvalue(sr1, sr2, 5.0, 2.0, weight_mode="unweighted", n_draws=99, n_grid=25, n_split=25, seed=rep * 53 + 7) > 0.05:
                c_cond_u += 1
            if two_body_lump_conditional_pvalue(sr1, sr2, 5.0, 2.0, weight_mode="z_smooth", n_draws=99, n_grid=25, n_split=25, seed=rep * 53 + 8) > 0.05:
                c_cond_z += 1
        print(
            f"  {label:20s}: Raw inf-conv={c_raw/reps:5.1%} | "
            f"Cond Lump LR={c_cond_u/reps:5.1%} | Z-Cond Lump LR={c_cond_z/reps:5.1%}"
        )

    print("\n" + "=" * 92)
    print("PUSH 3: DOOR B (Unnorm t_10 Wilks) vs GCLR+ vs SPHERICAL SCLR+ vs FLOOR UNION (3D, N=15, reps=150)")
    print("=" * 92)
    worlds = [
        ("Gaussian weak (r=1, m=2)", lambda: _make_scale_mixture_world(rng, 3, 15, 1.0, 2.0, None), 2.0),
        ("Gaussian strong (r=5, m=2)", lambda: _make_scale_mixture_world(rng, 3, 15, 5.0, 2.0, None), 2.0),
        ("Scale-mix nu=10 weak (r=1, m=2)", lambda: _make_scale_mixture_world(rng, 3, 15, 1.0, 2.0, 10), 2.0),
        ("Scale-mix nu=10 strong (r=5, m=2)", lambda: _make_scale_mixture_world(rng, 3, 15, 5.0, 2.0, 10), 2.0),
        ("Heavy-tail nu=2 (Cauchy-like, r=3, m=2)", lambda: _make_scale_mixture_world(rng, 3, 15, 3.0, 2.0, 2), 2.0),
    ]
    for label, gen, m_true in worlds:
        reps = 150
        c_door_b_wilks = c_gclr = c_sclr_g = c_sclr_t = c_floor = 0
        phi0 = math.atan(m_true)
        for rep in range(reps):
            rds = gen()
            srs = [SpanReading(r) for r in rds]
            pre = precompute_grid(srs, domain="positive", n_grid=31, scale_hint=m_true)
            e0 = evaluate_series(rds, phi0)
            t0 = sum((10.0 + k_i) * math.log1p(k_r / 10.0) for k_r, k_i in e0["per_K"])
            t_min = t0
            for g in range(len(pre[0])):
                tg = sum((10.0 + srs[i].k) * math.log1p(_k_of_z(pre[1][i][g], srs[i].z) / 10.0) for i in range(len(srs)))
                if tg < t_min:
                    t_min = tg
            if (t0 - t_min) <= CHI2_95_1:
                c_door_b_wilks += 1
            p_g = general_clr_pvalue(srs, phi0, domain="positive", weight_mode="unweighted", null_mode="gaussian", n_draws=99, seed=rep * 41 + 1, precomputed_grid=pre)["p_value"]
            p_sg = general_clr_pvalue(srs, phi0, domain="positive", weight_mode="unweighted", null_mode="spherical", n_draws=99, seed=rep * 41 + 2, precomputed_grid=pre)["p_value"]
            p_st = general_clr_pvalue(srs, phi0, domain="positive", weight_mode="unweighted", null_mode="spherical", nu_student=10.0, n_draws=99, seed=rep * 41 + 3, precomputed_grid=pre)["p_value"]
            if p_g > 0.05:
                c_gclr += 1
            if p_sg > 0.05:
                c_sclr_g += 1
            if p_st > 0.05:
                c_sclr_t += 1
            if max(p_g, p_sg) > 0.05:
                c_floor += 1
        print(
            f"  {label:40s}: DoorB-W={c_door_b_wilks/reps:5.1%} | GCLR+={c_gclr/reps:5.1%} | "
            f"SCLR+(G)={c_sclr_g/reps:5.1%} | SCLR+(t10)={c_sclr_t/reps:5.1%} | Floor={c_floor/reps:5.1%}"
        )

    print("\n" + "=" * 92)
    print("PUSH 3B: SPHERICAL SCLR+ & O7 FLOOR UNION ON DISCORDANT vs OVER-AGREEING SCENES")
    print("=" * 92)
    rds_far10 = [isotropic(5.0, 1.0, 0.1, 0.002)] * 5 + [isotropic(6.0, 1.0, 0.1, 0.002)] * 5
    rds_tight10 = [isotropic(5.5 + 0.002 * (i - 4.5), 1.0, 0.1, 0.002) for i in range(10)]
    rds_2d_far2 = [isotropic([5.0, 0.0], [1.0, 0.0], 0.1, 0.002), isotropic([6.0, 0.0], [1.0, 0.0], 0.1, 0.002)]
    rds_2d_aligned1 = [isotropic([5.0, 0.0], [1.0, 0.0], 0.1, 0.002)]
    rds_1d_n4 = [isotropic(2.0 + 0.005 * ((i % 3) - 1), 1.0, 0.01, 0.01) for i in range(4)]
    for label, rds, sh in [
        ("1D N=10 far (5x 5.0, 5x 6.0, sigma=0.1)", rds_far10, 5.5),
        ("1D N=10 over-agreeing (5.491..5.509, sigma=0.1)", rds_tight10, 5.5),
        ("2D N=2 far ((5,0) & (6,0), sigma=0.1)", rds_2d_far2, 5.5),
        ("2D N=1 perfectly aligned ((5,0), (1,0), sigma=0.1)", rds_2d_aligned1, 5.0),
        ("1D N=4 aligned (m=2.0, s=0.01) [exact 1D N<=4 floor]", rds_1d_n4, 2.0),
    ]:
        s_gclr = gclr_confidence_set(rds, domain="positive", null_mode="gaussian", scale_hint=sh, n_grid=121, n_draws=199)
        s_sclr = gclr_confidence_set(rds, domain="positive", null_mode="spherical", scale_hint=sh, n_grid=121, n_draws=199)
        s_sclr_t = gclr_confidence_set(rds, domain="positive", null_mode="spherical", nu_student=10.0, scale_hint=sh, n_grid=121, n_draws=199)
        s_floor = gclr_confidence_set(rds, domain="positive", floor_union=True, scale_hint=sh, n_grid=121, n_draws=199)
        print(f"\n  {label}")
        print(f"    Gaussian GCLR+  : {s_gclr}")
        print(f"    Spherical SCLR+ : {s_sclr}")
        print(f"    Spherical t_10  : {s_sclr_t}")
        print(f"    O7 Floor Union  : {s_floor}")


if __name__ == "__main__":
    run_push_benchmarks()
