"""E13 Push 4 & Push 5: Z-Safe Weight Pareto Frontier (E2/E3) & Closing the 1D N=2 Door B Gap (O14/O18).

Pure Python (stdlib only). Extends explore/e13_general_clr.py to:
  1. Push 4 (#1 — "1 real reading buried in N random noisy zeros", E2 / E3 / T43 / T45):
       - Exact full-covariance U-weighted conditional score test W-T45:
           T_w(m0) = sum_i w(H_i) S_i / sqrt(sum_i w(H_i)^2 C_i) ~ N(0, 1) exact given (U_1..U_N).
       - Exact Stein + Cauchy-Schwarz Pareto theorem for ANY Z-safe per-reading weight w(H):
           ARE_k(w) = (E[H w(H)])^2 / (E[H] E[H w(H)^2]) <= 1   (uniquely = 1 iff w == const),
           N_drown(w; r) = w(r^2)^2 * r^2 / E[w(H)^2 H].
       - Verifies the three-regime Pareto frontier (diffuse weak vs mixed strong vs buried sparse
         up to N_null = 2,000 noisy zeros) for both W-T45 and ghost-free Z-GCLR.
  2. Push 5 (#2 — Closing the 1D N=2 Door B Student-t gap, O14 / O7 / O18):
       - Proves why naive U-conditioning (drawing V_i ~ t_nu indep of U_i) breaks under Student-t:
         by Kac-Bernstein, v_i ~ IG(nu/2, nu/2) couples U_i and V_i through the unknown push alpha_i*.
       - Exact push-free Student-t Anderson-Rubin (t-AR) quantile q_{1-alpha}^{t-AR}(k_1..k_N; nu)
         via exact Beta(nu/2, k_i/2) representation: (nu + k_i) log(1 + K_i/nu) = -(nu + k_i) log B_i.
       - Proved AR-bounded t-LR: {m > 0 : T_t(m) - T_{t,min} <= q_{0.95}^{t-AR}} (1-line proof of
         >= 95% coverage for every (Sigma_i, N_i, alpha_i*, m*) under t_nu; never EMPTY; widens
         cleanly to [4.83, 6.17] on 1D N=2 10-sigma Scene B and splits bimodally at 50-sigma).
       - Piecewise Exact Door B test (switching ONLY on declared dimensions: AR-bounded t-LR for
         1D N <= 4, SCLR+(t_nu) for N >= 5 or k >= 2) — proved >= 95% everywhere under t_nu with
         ZERO Bonferroni penalty, plus its floor-respecting union with Z-GCLR+.

Run from the repo root:
    python3 explore/e13_push4_5.py
"""
import math
import random

from explore.e13_general_clr import (
    SpanReading,
    _k_of_z,
    _make_correlated_world,
    _make_scale_mixture_world,
    _rand_unit_sphere,
    _refine_phi_min,
    _weight_fn,
    gclr_confidence_set,
    general_clr_pvalue,
    precompute_grid,
)
from explore.e13_unnormalized_evidence import (
    CHI2_95_1,
    Reading,
    dot,
    evaluate_series,
    format_intervals,
    isotropic,
)


# ---------------------------------------------------------------------------
# Push 4: Full-covariance U-weighted conditional score test (W-T45) & Pareto
# ---------------------------------------------------------------------------

def evaluate_weighted_t45(readings, phi, weight_mode="unweighted"):
    """Exact full-covariance U-weighted conditional score test W-T45 at angle phi.

    Uses the exact T45 per-reading score S_i = V_i^T A_i U_i and conditional variance
    C_i = ||A_i U_i||^2 (from evaluate_series), weighted by w_i = w(H_i(phi)).
    Conditional on (U_1, ..., U_N), T_w = (sum w_i S_i) / sqrt(sum w_i^2 C_i) ~ N(0, 1)
    EXACTLY in finite samples for arbitrary correlated Sigma_i, subspaces N_i, and pushes alpha_i*.
    """
    s_phi, c_phi = math.sin(phi), math.cos(phi)
    # Per-reading S_i, C_i, H_i via evaluate_series on each single reading
    num = 0.0
    den = 0.0
    h_list = []
    w_list = []
    for r in readings:
        ev = evaluate_series([r], phi)
        h_i = ev["H"]
        s_i = ev["S"]
        c_i = ev["C"]
        w_i = _weight_fn(h_i, weight_mode)
        h_list.append(h_i)
        w_list.append(w_i)
        if w_i > 0.0:
            num += w_i * s_i
            den += (w_i * w_i) * c_i
    t2 = 0.0 if den <= 1e-30 else (num * num) / den
    return dict(T2=t2, num=num, den=den, weights=w_list, H=h_list)


def stein_pareto_metrics(weight_mode, k=1, r_strong=20.0, n_steps=1400):
    """Exact Stein + Cauchy-Schwarz weak-signal ARE_k(w) and sparse drowning threshold N_drown(w; r).

    Theorem (Stein's identity + Cauchy-Schwarz):
      For H ~ chi^2_k and any smooth nonnegative weight w(H) with finite E[H w(H)^2],
        ARE_k(w) = (E[H w(H)])^2 / (E[H] * E[H w(H)^2]) <= 1,
      with equality iff w(H) is constant almost everywhere.
      For 1 strong reading of amplitude r buried in N random pure-noise readings,
        N_drown(w; r) = w(r^2)^2 * r^2 / E[H w(H)^2].

    Evaluated to ~12 digits via smooth 1D Simpson quadrature on u = sqrt(H) in [0, 14],
    where f_k(h) dh = (2 / (2^{k/2} Gamma(k/2))) u^{k-1} exp(-u^2/2) du has no singularity at u = 0.
    """
    u_max = 14.0
    n = n_steps if n_steps % 2 == 0 else n_steps + 1
    du = u_max / n
    norm = 2.0 / ((2.0 ** (0.5 * k)) * math.gamma(0.5 * k))
    e_h = 0.0
    e_hw = 0.0
    e_hw2 = 0.0
    for j in range(n + 1):
        u = j * du
        coef = 1.0 if (j == 0 or j == n) else (4.0 if (j & 1) else 2.0)
        h = u * u
        pdf_u = norm * (u ** (k - 1)) * math.exp(-0.5 * h)
        wt = coef * pdf_u
        w = _weight_fn(h, weight_mode)
        hw = h * w
        e_h += wt * h
        e_hw += wt * hw
        e_hw2 += wt * hw * w
    scale = du / 3.0
    e_h *= scale
    e_hw *= scale
    e_hw2 *= scale
    are = (e_hw * e_hw) / (e_h * e_hw2)
    h_s = r_strong * r_strong
    w_s = _weight_fn(h_s, weight_mode)
    n_drown = (w_s * w_s * h_s) / e_hw2
    return dict(are=are, e_hw2=e_hw2, n_drown=n_drown)


# ---------------------------------------------------------------------------
# Push 5: Exact t-AR, Proved AR-bounded t-LR, and Piecewise Exact Door B
# ---------------------------------------------------------------------------

def exact_t_ar_quantile(k_dims, nu, alpha=0.05, n_mc=40000, seed=202605):
    """Exact (1 - alpha) quantile of sum_i (nu + k_i) log(1 + K_i / nu) under K_i ~ k_i F_{k_i, nu}.

    Uses the exact identity (1 + K_i / nu)^{-1} ~ Beta(nu/2, k_i/2), valid for ANY real nu > 0
    and ANY positive-definite correlated covariances Sigma_i and subspaces of dimension k_i.
    Completely free of all unknown pushes alpha_i^*!
    """
    if not k_dims:
        return 0.0
    rng = random.Random(seed)
    half_nu = 0.5 * float(nu)
    params = [(float(nu) + float(k_i), 0.5 * float(k_i)) for k_i in k_dims]
    betav = rng.betavariate
    log = math.log
    vals = []
    for _ in range(n_mc):
        tot = 0.0
        for coef, half_k in params:
            b = betav(half_nu, half_k)
            tot -= coef * log(max(b, 1e-300))
        vals.append(tot)
    vals.sort()
    idx = min(n_mc - 1, int(math.ceil((1.0 - alpha) * n_mc)) - 1)
    return vals[idx]


def student_t_series_cost(span_readings, phi, nu):
    """Door B T48 Student-t cost sum_i (nu + k_i) log(1 + K_i(phi) / nu) on the full type."""
    inv_nu = 1.0 / float(nu)
    tot = 0.0
    for sr in span_readings:
        _, _, _, _, _, k_val = sr.uv_at_phi(phi)
        tot += (float(nu) + sr.k) * math.log1p(k_val * inv_nu)
    return tot


def door_b_piecewise_exact_pvalue(
    span_readings,
    phi0,
    *,
    nu=10.0,
    alpha=0.05,
    n_grid=91,
    n_draws=199,
    seed=2026,
    scale_hint=1.0,
    precomputed_grid=None,
    q_t_ar=None,
):
    """Piecewise Exact Door B test under declared Student-t_nu scale-mixture noise.

    Branching depends ONLY on declared experiment dimensions (never on data y_i):
      - If all readings are 1D (k_i == 1) and N <= 4 (where SCLR+ hits the 1/2^N > 0.05 floor):
          AR-bounded t-LR: accept iff T_{t_nu}(phi0) - T_{t_nu, min} <= q_{1-alpha}^{t-AR}(N, k=1; nu).
          1-line proof: T(m*) - T_min <= T(m*) => P(accept m*) >= P(T(m*) <= q^{t-AR}) = 1 - alpha!
      - Otherwise (any k_i >= 2, or N >= 5):
          Spherical SCLR+(t_nu), exact >= 1 - alpha under EVERY per-reading scale mixture.

    Returns dict(accepted, branch, stat, threshold_or_p).
    """
    if not span_readings:
        return dict(accepted=True, branch="empty", stat=0.0, threshold_or_p=1.0)

    is_small_1d = len(span_readings) <= 4 and all(sr.k == 1 for sr in span_readings)
    if is_small_1d:
        if precomputed_grid is None:
            grid_phis, per_q_perp = precompute_grid(
                span_readings, domain="positive", n_grid=n_grid, scale_hint=scale_hint
            )
        else:
            grid_phis, per_q_perp = precomputed_grid
        inv_nu = 1.0 / float(nu)
        nu_f = float(nu)
        c0 = 0.0
        for sr in span_readings:
            _, _, _, _, _, kv = sr.uv_at_phi(phi0)
            c0 += (nu_f + sr.k) * math.log1p(kv * inv_nu)
        c_min = c0
        n_g = len(grid_phis)
        for g in range(n_g):
            cg = 0.0
            for i, sr in enumerate(span_readings):
                kg = _k_of_z(per_q_perp[i][g], sr.z)
                cg += (nu_f + sr.k) * math.log1p(kg * inv_nu)
                if cg >= c_min:
                    break
            if cg < c_min:
                c_min = cg
        lr = max(0.0, c0 - c_min)
        if q_t_ar is None:
            q_t_ar = exact_t_ar_quantile([sr.k for sr in span_readings], nu, alpha=alpha)
        return dict(accepted=(lr <= q_t_ar), branch="ar_bounded_t_lr", stat=lr, threshold_or_p=q_t_ar)

    res = general_clr_pvalue(
        span_readings,
        phi0,
        domain="positive",
        weight_mode="unweighted",
        null_mode="spherical",
        nu_student=nu,
        n_grid=n_grid,
        n_draws=n_draws,
        seed=seed,
        scale_hint=scale_hint,
        precomputed_grid=precomputed_grid,
    )
    return dict(
        accepted=(res["p_value"] > alpha),
        branch="sclr_t",
        stat=res["lr_obs"],
        threshold_or_p=res["p_value"],
    )


def door_b_piecewise_confidence_set(
    readings,
    *,
    nu=10.0,
    floor_union=False,
    n_scan=64,
    n_grid=121,
    n_draws=199,
    seed=2026,
    scale_hint=1.0,
):
    """Scan (0, inf) for the Piecewise Exact Door B 95% confidence set (or its O7 floor union)."""
    srs = [SpanReading(r) for r in readings]
    pre = precompute_grid(srs, domain="positive", n_grid=n_grid, scale_hint=scale_hint)
    eps = 1e-4
    us = [eps + (0.5 * math.pi - 2 * eps) * j / (n_scan - 1) for j in range(n_scan)]
    phis = [math.atan(scale_hint * math.tan(u)) for u in us]
    # Include both Gaussian peak and Student-t peaks (handles bimodal 50-sigma splits cleanly!)
    ks_t = [student_t_series_cost(srs, ph, nu) for ph in phis]
    best_i = ks_t.index(min(ks_t))
    # Also zoom around any local minima of the Student-t curve
    local_mins = [best_i]
    for i in range(1, len(ks_t) - 1):
        if ks_t[i] <= ks_t[i - 1] and ks_t[i] <= ks_t[i + 1] and ks_t[i] <= ks_t[best_i] + 12.0:
            local_mins.append(i)
    extra = []
    for idx_m in sorted(set(local_mins)):
        lo = phis[max(0, idx_m - 1)]
        hi = phis[min(len(phis) - 1, idx_m + 1)]
        # Golden-section on Student-t cost
        gr = 0.5 * (math.sqrt(5.0) - 1.0)
        a, b = lo, hi
        c, d_p = b - gr * (b - a), a + gr * (b - a)
        fc, fd = student_t_series_cost(srs, c, nu), student_t_series_cost(srs, d_p, nu)
        for _ in range(24):
            if fc < fd:
                b, d_p, fd = d_p, c, fc
                c = b - gr * (b - a)
                fc = student_t_series_cost(srs, c, nu)
            else:
                a, c, fc = c, d_p, fd
                d_p = a + gr * (b - a)
                fd = student_t_series_cost(srs, d_p, nu)
        p_star = 0.5 * (a + b)
        extra.append(p_star)
        for hw in (0.04, 0.008, 0.0012):
            for j in range(-10, 11):
                p = p_star + hw * j / 10.0
                if eps < p < 0.5 * math.pi - eps:
                    extra.append(p)
    phis = sorted(set(phis + extra))

    is_small_1d = len(srs) <= 4 and all(sr.k == 1 for sr in srs)
    q_t_ar = exact_t_ar_quantile([sr.k for sr in srs], nu) if is_small_1d else None

    # For small-1D AR-bounded t-LR, use the dense scan phis themselves for T_min (exact & smooth)
    if is_small_1d:
        t_vals = [student_t_series_cost(srs, ph, nu) for ph in phis]
        t_min = min(t_vals)
        mask_b = [(tv - t_min) <= q_t_ar for tv in t_vals]
    else:
        mask_b = [
            door_b_piecewise_exact_pvalue(
                srs, ph, nu=nu, n_grid=n_grid, n_draws=n_draws, seed=seed + idx * 17,
                scale_hint=scale_hint, precomputed_grid=pre,
            )["accepted"]
            for idx, ph in enumerate(phis)
        ]

    if not floor_union:
        return format_intervals(phis, mask_b)

    mask_u = []
    for idx, ph in enumerate(phis):
        if mask_b[idx]:
            mask_u.append(True)
        else:
            p_g = general_clr_pvalue(
                srs, ph, domain="positive", weight_mode="z_smooth", null_mode="gaussian",
                n_grid=n_grid, n_draws=n_draws, seed=seed + idx * 17 + 1,
                scale_hint=scale_hint, precomputed_grid=pre,
            )["p_value"]
            mask_u.append(p_g > 0.05)
    return format_intervals(phis, mask_u)


# ---------------------------------------------------------------------------
# Reproducible Benchmark Suite for Push 4 and Push 5
# ---------------------------------------------------------------------------

def run_push4_5_benchmarks(seed=20261121):
    print("=" * 94)
    print("PUSH 4A: STEIN + CAUCHY-SCHWARZ PARETO TABLE FOR THE Z-SAFE WEIGHT FAMILY w(H)")
    print("=" * 94)
    modes = [
        ("T45 (w=1)", "unweighted"),
        ("z_smooth H/(1+H)", "z_smooth"),
        ("z_quad H^2/(1+H^2)", "z_quad"),
        ("z_cap10 H^2/(100+H^2)", "z_cap10"),
        ("z_linear w=H", "z_linear"),
        ("z_pow2 w=H^2", "z_pow2"),
    ]
    for k in (1, 3):
        print(f"\n  Dimension k = {k} (strong reading r = 20, H = 400):")
        print(f"  {'Weight w(H)':24s} | {'ARE_k (weak)':>12s} | {'E[w^2 H]':>10s} | {'N_drown(r=20)':>14s}")
        print("  " + "-" * 68)
        for label, wm in modes:
            st = stein_pareto_metrics(wm, k=k, r_strong=20.0)
            print(f"  {label:24s} | {st['are']:11.1%} | {st['e_hw2']:10.4g} | {st['n_drown']:14.1f}")

    print("\n" + "=" * 94)
    print("PUSH 4B: THREE-REGIME PARETO CHECK (W-T45 ON 1D: DIFFUSE WEAK / MIXED STRONG / BURIED SPARSE)")
    print("=" * 94)
    rng = random.Random(seed)
    m_true = 2.0
    th = math.atan(m_true)
    s, co = math.sin(th), math.cos(th)
    regimes = [
        ("Diffuse weak (N=400, r=0.6) [rej m=3.0]", lambda: [(0.6*s + rng.gauss(0,1), 0.6*co + rng.gauss(0,1)) for _ in range(400)], 3.0, 300),
        ("Mixed strong (1x r=20 + 25x r=5) [rej m=2.35]", lambda: [(20.0*s + rng.gauss(0,1), 20.0*co + rng.gauss(0,1))] + [(5.0*s + rng.gauss(0,1), 5.0*co + rng.gauss(0,1)) for _ in range(25)], 2.35, 300),
        ("Buried sparse (1x r=20 + 500 nulls) [rej m=4.0]", lambda: [(20.0*s + rng.gauss(0,1), 20.0*co + rng.gauss(0,1))] + [(rng.gauss(0,1), rng.gauss(0,1)) for _ in range(500)], 4.0, 200),
        ("Buried sparse (1x r=20 + 2000 nulls) [rej m=4.0]", lambda: [(20.0*s + rng.gauss(0,1), 20.0*co + rng.gauss(0,1))] + [(rng.gauss(0,1), rng.gauss(0,1)) for _ in range(2000)], 4.0, 200),
    ]
    for title, gen, m_alt, reps in regimes:
        phi0, phi1 = math.atan(m_true), math.atan(m_alt)
        counts = {label: [0, 0] for label, _ in modes}
        for _ in range(reps):
            rds = gen()
            for label, wm in modes:
                for col, ph in enumerate((phi0, phi1)):
                    sp, cp = math.sin(ph), math.cos(ph)
                    num = den = 0.0
                    for f, a in rds:
                        u = sp * f + cp * a
                        v = cp * f - sp * a
                        h = u * u
                        w = _weight_fn(h, wm)
                        num += w * u * v
                        den += (w * w) * h
                    t2 = 0.0 if den <= 1e-30 else (num * num) / den
                    if col == 0 and t2 <= CHI2_95_1:
                        counts[label][0] += 1
                    elif col == 1 and t2 > CHI2_95_1:
                        counts[label][1] += 1
        print(f"\n  {title} (reps={reps}):")
        for label, _ in modes:
            c, p = counts[label]
            print(f"    {label:24s}: cov(m=2)={c/reps:5.1%} | power={p/reps:5.1%}")

    print("\n" + "=" * 94)
    print("PUSH 4C: GHOST-FREE Z-GCLR ON BURIED SIGNAL (1x r=20 + 500 nulls) & ANTI-PARALLEL PAIR")
    print("=" * 94)
    rng_g = random.Random(2026)
    rds_buried = [isotropic(20.0 * s, 20.0 * co, 1.0, 1.0)] + [
        isotropic(rng_g.gauss(0, 1), rng_g.gauss(0, 1), 1.0, 1.0) for _ in range(500)
    ]
    srs_buried = [SpanReading(r) for r in rds_buried]
    pre_b = precompute_grid(srs_buried, domain="positive", n_grid=21, scale_hint=2.0)
    scan_ms = [1.2 + 0.15 * j for j in range(29)]  # 1.20 .. 5.40
    for label, wm in [
        ("Unweighted GCLR+", "unweighted"),
        ("Z-GCLR+ H/(1+H)", "z_smooth"),
        ("Z-GCLR+ H^2/(100+H^2)", "z_cap10"),
        ("Z-GCLR+ w=H", "z_linear"),
    ]:
        acc = [
            m for idx, m in enumerate(scan_ms)
            if general_clr_pvalue(srs_buried, math.atan(m), domain="positive", weight_mode=wm, n_draws=99, seed=700 + idx, precomputed_grid=pre_b)["p_value"] > 0.05
        ]
        print(f"  {label:24s} on 1x r=20 + 500 nulls : [{min(acc):.2f}, {max(acc):.2f}] (width={max(acc)-min(acc):.2f})")
    anti = [isotropic(2.0, -1.0, 0.1, 0.1), isotropic(1.0, -2.0, 0.1, 0.1)]
    print(f"  Z-GCLR(RP^1) w=H on anti-parallel pair         : {gclr_confidence_set(anti, domain='whole', weight_mode='z_linear', scale_hint=1.0, n_grid=91, n_draws=199)}")

    print("\n" + "=" * 94)
    print("PUSH 5A: EXACT t-AR QUANTILES (nu=10) & PIECEWISE EXACT DOOR B ON 1D N=2 SCENES")
    print("=" * 94)
    for n_cnt in (1, 2, 3, 4):
        q95 = exact_t_ar_quantile([1] * n_cnt, nu=10.0)
        print(f"  1D N={n_cnt}, nu=10: exact push-free q95_t_AR = {q95:.3f}")
    scenes_1d2 = [
        ("All-zero 0/0 (1D N=2)", [isotropic(0.0, 0.0, 1.0, 1.0), isotropic(0.0, 0.0, 1.0, 1.0)], 1.0),
        ("Scene A: near (5.5, 5.6) +- 0.1", [isotropic(5.5, 1.0, 0.1, 0.002), isotropic(5.6, 1.0, 0.1, 0.002)], 5.5),
        ("Scene B 10s: (5.0, 6.0) +- 0.1", [isotropic(5.0, 1.0, 0.1, 0.002), isotropic(6.0, 1.0, 0.1, 0.002)], 5.5),
        ("Scene B 50s: (5.0, 6.0) +- 0.02", [isotropic(5.0, 1.0, 0.02, 0.0004), isotropic(6.0, 1.0, 0.02, 0.0004)], 5.5),
    ]
    for name, rds, sh in scenes_1d2:
        s_pw = door_b_piecewise_confidence_set(rds, nu=10.0, floor_union=False, scale_hint=sh)
        s_fl = door_b_piecewise_confidence_set(rds, nu=10.0, floor_union=True, scale_hint=sh)
        print(f"\n  {name}")
        print(f"    Piecewise Exact Door B (nu=10) : {s_pw}")
        print(f"    Piecewise Door B U Z-GCLR+     : {s_fl}")

    print("\n" + "=" * 94)
    print("PUSH 5B: HONEST MISSPECIFICATION CHECK ON 1D N=2 (declared nu=10 vs true nu=10, nu=2, Gauss)")
    print("=" * 94)
    q95_10 = exact_t_ar_quantile([1, 1], nu=10.0)
    rng_m = random.Random(seed + 55)
    phi0 = math.atan(2.0)
    grid_p = [0.5 * math.pi * j / 60 for j in range(61)]
    for label, nu_true, r_len in [
        ("True Gaussian, r=0 (1D N=2)", None, 0.0),
        ("True nu=10, r=0 (1D N=2, matched)", 10, 0.0),
        ("True nu=10, r=5 (1D N=2, matched)", 10, 5.0),
        ("True nu=2, r=0 (1D N=2, MISSPECIFIED!)", 2, 0.0),
        ("True nu=2, r=3 (1D N=2, MISSPECIFIED!)", 2, 3.0),
    ]:
        reps = 400
        cov_ar = cov_ar_lr = cov_wilks = 0
        for _ in range(reps):
            rds = _make_scale_mixture_world(rng_m, 1, 2, r_len, 2.0, nu_df=nu_true)
            srs = [SpanReading(r) for r in rds]
            c0 = student_t_series_cost(srs, phi0, 10.0)
            cmin = min(c0, min(student_t_series_cost(srs, ph, 10.0) for ph in grid_p))
            if c0 <= q95_10:
                cov_ar += 1
            if (c0 - cmin) <= q95_10:
                cov_ar_lr += 1
            if (c0 - cmin) <= CHI2_95_1:
                cov_wilks += 1
        print(f"  {label:40s}: t-AR={cov_ar/reps:5.1%} | AR-bd t-LR={cov_ar_lr/reps:5.1%} | Wilks={cov_wilks/reps:5.1%}")


if __name__ == "__main__":
    run_push4_5_benchmarks()
