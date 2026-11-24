"""E15 verification library for Door A uniqueness (`T45` vs `M_nz-GCLR` vs `Z-GCLR`) (pure stdlib).

Exploration branch `explore/e15-door-a-uniqueness` (from `origin/main` `62e0390`, 24 Nov 2026).
Zero external dependencies; runs under plain `python3 -m unittest`.

Provides:
  - Exact 1-instrument conditional likelihood-ratio statistics on the positive ray [0, inf] and
    the full projective line (`lr_positive_weighted`, `lr_line_weighted`)
  - Exact Monte Carlo conditional p-values for `M_nz-GCLR` and `Z-GCLR` (`gclr_pvalue_1instr`)
  - Exact `T45` score test and `psi_hill` hilltop-killing competitor (for Theorem E15-1 UMPU check)
  - Exact conservative chi-square band supremum (`chi2_band_sup_prob`, Theorem E15-1(c))
  - Door A1 full-type wrapper (`door_a1`) built directly on `met/ours34.py`'s `State` (zero discount)
"""
import math
import random

from met.ours34 import State, classify_shape, chi2_sf


Z95 = 1.959963984540054
C95_CENTER = 0.06270677794321385  # P(|N(0, 1)| < C95_CENTER) == 0.05


def lr_line_weighted(Us, Vs, ws):
    """Weighted LR drop W(m0) - min_{whole line} W(m) on one isotropic instrument in the
    candidate (u_0, v_0) frame: W(delta) = sum w_i (U_i sin delta - V_i cos delta)^2."""
    q_uu = math.fsum(w * u * u for w, u in zip(ws, Us))
    q_uv = math.fsum(w * u * v for w, u, v in zip(ws, Us, Vs))
    q_vv = math.fsum(w * v * v for w, v in zip(ws, Vs))
    return max(0.0, 0.5 * (q_vv - q_uu + math.hypot(q_uu - q_vv, 2.0 * q_uv)))


def lr_positive_weighted(Us, Vs, ws, theta0):
    """Exact weighted LR drop W(theta0) - min_{theta in [0, pi/2]} W(theta) on one isotropic
    instrument, where delta = theta - theta0 in [-theta0, pi/2 - theta0]."""
    q_uu = math.fsum(w * u * u for w, u in zip(ws, Us))
    q_uv = math.fsum(w * u * v for w, u, v in zip(ws, Us, Vs))
    q_vv = math.fsum(w * v * v for w, v in zip(ws, Vs))
    # W(delta) = 0.5*(q_uu + q_vv) + A*cos(2*delta) - B*sin(2*delta),
    # with A = 0.5*(q_vv - q_uu), B = q_uv, and 2*delta in [-2*theta0, pi - 2*theta0].
    mid = 0.5 * (q_uu + q_vv)
    A = 0.5 * (q_vv - q_uu)
    B = q_uv
    lo_2d = -2.0 * theta0
    hi_2d = math.pi - 2.0 * theta0
    cands = [
        mid + A * math.cos(lo_2d) - B * math.sin(lo_2d),
        mid + A * math.cos(hi_2d) - B * math.sin(hi_2d),
    ]
    R = math.hypot(A, B)
    if R > 0.0:
        # Global minimum of A*cos(phi) - B*sin(phi) is at (cos phi, sin phi) = (-A/R, B/R)
        phi_min = math.atan2(B, -A)
        for shift in (-2.0 * math.pi, 0.0, 2.0 * math.pi):
            p = phi_min + shift
            if lo_2d - 1e-12 <= p <= hi_2d + 1e-12:
                cands.append(mid - R)
                break
    return max(0.0, q_vv - min(cands))


def gclr_pvalue_1instr(pairs, m0, *, s=1.0, weight=None, drop_exact_zero=True,
                       domain="positive", draws=2000, seed=0):
    """Exact conditional Monte Carlo p-value on one isotropic scalar instrument (noise ratio s).
    `pairs`: list of (F/sf, a/sa) whitened scalar readings.
    `weight`: None for unweighted GCLR (w=1), or callable w(H_i) for Z-GCLR."""
    theta0 = 0.5 * math.pi if math.isinf(m0) else math.atan(m0 / s)
    su, cu = math.sin(theta0), math.cos(theta0)
    # Candidate line u = (sin theta0, cos theta0), leftover v = (cos theta0, -sin theta0)
    Us, Vs, ws = [], [], []
    for yf, ya in pairs:
        if drop_exact_zero and yf == 0.0 and ya == 0.0:
            continue
        u = yf * su + ya * cu
        v = yf * cu - ya * su
        w = 1.0 if weight is None else float(weight(u * u))
        if w > 0.0:
            Us.append(u)
            Vs.append(v)
            ws.append(w)
    if not Us:
        return dict(p_value=1.0, lr_obs=0.0, active=0)
    stat_fn = (lambda U, V, W: lr_positive_weighted(U, V, W, theta0)) if domain == "positive" else lr_line_weighted
    lr_obs = stat_fn(Us, Vs, ws)
    if lr_obs <= 1e-13:
        return dict(p_value=1.0, lr_obs=0.0, active=len(Us))
    rng = random.Random(seed)
    ge = 0
    tol = 1e-11 * max(1.0, lr_obs)
    n = len(Us)
    for _ in range(draws):
        v_draw = [rng.gauss(0.0, 1.0) for _ in range(n)]
        if stat_fn(Us, v_draw, ws) >= lr_obs - tol:
            ge += 1
    return dict(p_value=(1 + ge) / (draws + 1), lr_obs=lr_obs, active=n)


def chi2_band_sup_prob(a, b, k_max=None):
    """Compute sup_{k >= 1} P(a <= chi2_k <= b) exactly via chi2_sf(a, k) - chi2_sf(b, k).
    Demonstrates why exact size == alpha (not conservative <= alpha) is required in Theorem E15-1."""
    if k_max is None:
        k_max = max(30, int(b) + 40)
    best_p, best_k = 0.0, 1
    for k in range(1, k_max + 1):
        p = max(0.0, chi2_sf(a, k) - chi2_sf(b, k))
        if p > best_p:
            best_p, best_k = p, k
    return best_p, best_k


def simulate_t45_vs_psi_hill(M, signal_per_comp, delta, *, trials=8000, seed=0):
    """Compare exact rejection rates of T45 and the hilltop-killing exact-size competitor
    psi_hill(q, X, W) = 1{|X| > z95} if q >= W else 1{|X| < c95} on one isotropic instrument."""
    rng = random.Random(seed)
    c_d, s_d = math.cos(delta), math.sin(delta)
    rej_t45 = rej_hill = 0
    for _ in range(trials):
        q = uv = vv = 0.0
        for _ in range(M):
            u = signal_per_comp * c_d + rng.gauss(0.0, 1.0)
            v = signal_per_comp * s_d + rng.gauss(0.0, 1.0)
            q += u * u
            uv += u * v
            vv += v * v
        x = uv / math.sqrt(q)
        w = max(0.0, vv - x * x)
        if abs(x) > Z95:
            rej_t45 += 1
        if (abs(x) > Z95) if (q >= w) else (abs(x) < C95_CENTER):
            rej_hill += 1
    return dict(t45=rej_t45 / trials, psi_hill=rej_hill / trials)


def door_a1(readings, *, level=0.95, scale=None, grid=361):
    """Pure Door A1 equation on the full reading type: zero discount, additive state (K, C),
    point = argmin K on [0, inf] (T81), set = T45 exact score set, leftover = K_min beside it."""
    st = State(readings, scale=scale)
    pt = st.argmin_K(grid=grid)
    kept = st.kept_intervals(level=level, kind="t45", grid=grid)
    return dict(
        point=pt["mass"],
        point_status=pt["status"],
        point_is_convention=pt["convention"],
        minima=pt["minima"],
        set=kept,
        shape=classify_shape(kept),
        leftover=dict(
            K_min=pt["K_min"],
            M_total=st.M_total,
            chi2_bound=chi2_sf(pt["K_min"], st.M_total) if st.M_total > 0 else 1.0,
        ),
        state=st,
        level=float(level),
    )
