# E16: `met/ours35.py` — Executable Door A1 Probe & Premise Scorecard

**Branch:** `explore/e16-ours35-door-a` (off `origin/main` `62e0390`, 24 Nov 2026; Antigravity with Can).  
**Status:** Exploration probe only (`explore/README.md` rules). Nothing in `definition/`, `problems/`, `rulings.md`, or existing `met/` files is modified.  
**Builds directly on:** `explore/e15-door-a-uniqueness` (Theorems `E15-1`–`E15-3` and the `definition/proof.md` audit) and `main`'s `T45`, `T73`, `T74`, `T77 (6)`, `T80`, `T81`, `T82`.  
**Standard mathematics:** Lehmann & Romano (2005) §§4.2–4.3 (conditional score / UMPU); Moreira (2003); Anderson & Rubin (1949); Zehna (1966) (invariance of ML point under composition); Rockafellar (1970) (infimal convolution).

---

## 1. What `ours35` Is (Pure Door A1, Zero Free Knobs)

`met/ours35.py` strips away every `ours34` discount knob (`t_joint`, `M_nz`, `single_scalar_policy`, `chi2_bound(M_total)`) and implements the exact parameter-free Door A1 equation on the full reading type (arbitrary dimension $d \in \{1, 2, 3\}$, correlated SPD $\Sigma_i \in \mathbb R^{2d \times 2d}$, orthonormal candidate subspaces $N_i \in \mathbb R^{d \times k_i}$, mixed instruments, optional bookkeeping `window`, pure Python 3 stdlib, duck-typed to accept `met.full.Reading`):

- **Carried state (`T77 (6)`):** Pure additive monoid `State = (K(m), C(m))` with **zero hidden reading count** (`M` and `M_nz` do not exist on `State`). Independent blocks combine by exact addition `s1 + s2`. Exact `0/0` readings contribute `(0, 0)` and are the exact additive identity.
- **Curve (`T73`):** $K(m) = \sum_i \|(I - P_i(m)) z_i\|^2$ on the closed ray $m \in [0, \infty]$.
- **Point (`T81`):** $m^* = \arg\min_{m \in [0, \infty]} K(m)$. Returns `point_is_convention=True` (declared scale) on flat or tied curves (`R-2026-11-22-2`), and tags boundary minima cleanly (`boundary_zero`, `boundary_inf`) pending Can's ruling on `T81 (5)`.
- **Compatible set (`T45`, `E15-1`, `E15-2`):** $\{m \in [0, \infty] : |Z_{\rm T45}(m)| < z_{1-\alpha/2}\}$, where $Z_{\rm T45}(m) = S(m) / \sqrt{C(m)} \sim \mathcal N(0, 1)$ exact given $U$ (and $p = 1$ at $C = 0$). Closed-form, zero Monte Carlo, uniquely UMPU/LMPU on $(K, C)$ under Door A.
- **Leftover (`T73 (5)`):** `leftover = K_min` alone (pure function of $K$, so every single output field is exact `Z`, continuous at `0/0`, and `L17` reading-split invariant).
- **Exact zero-noise branch (`T1`):** `exact_mass` and `exact_series` for $\sigma = 0$.
- **Two-body composition (`L5` `lump`, `L14` `common_force`):** `T74` infimal convolution curve, `T81` exact point addition / harmonic dual, and exact `Q_joint` set composition directly from `(K, C)` via $\exp(-(Z_1^2 + Z_2^2)/2)$ ($\chi^2_2$ survival when both bodies have $C_j > 0$, collapsing cleanly to exact 1-dof $p_{\rm T45}$ when either body is `0/0`). Known-mass shifts `lump_known` / `common_force_known` are exact 1-dof.

---

## 2. Numerical Upgrades Over `met/ours34.py`

Building and torture-testing `ours35` caught and resolved four numerical/structural defects that affected `ours34`:

1. **Scale-normalized per-reading pre-whitening (`_SpanPre`):** Each reading is pre-whitened at its own natural scale $s_i = \sqrt{\operatorname{tr}\Sigma_{FF,i} / \operatorname{tr}\Sigma_{aa,i}}$ and mapped to the common series angle $\theta = \arctan(m / s_{\rm ref})$ via the exact analytic Jacobian $d\theta_i / d\theta = r_i / ((r_i \sin\theta)^2 + \cos^2\theta)$. Result: exact 14–15 digit `L3` unit equivariance across a **$10^{400}$ dynamic range** ($s = 10^{-200}$ to $s = 10^{+200}$, zero overflow or underflow).
2. **Swap-symmetric two-half ray parameterization (`(0, r)` and `(1, t)`):** Parameterizing $m = s_{\rm ref}\,r$ ($r \in [0, 1]$) and $m = s_{\rm ref}/t$ ($t \in [0, 1]$) preserves full 53-bit precision near $m = \infty$ and exact `T5` channel-swap symmetry.
3. **Exact `S = 0` stationary-point discovery + pinhole cancellation detection:** Bisecting every sign change of $S_\theta$ finds every local minimum, every hilltop maximum, and every physical right-angle pinhole ($b = 0 \implies U = 0 \implies C = 0$) regardless of grid alignment. Fixes `ours34`'s off-grid non-equivariance on single anti-parallel readings (`(10, -7.13)` now returns the exact singleton pinhole `{0.713}` and transforms equivariantly under `L3` scaling and `T5` swap) and resolves ultra-narrow `SNR = 1e10` intervals (`width = 1.23959e-9`, matching theory to 7+ digits) and ultra-narrow high-SNR hilltops (`(1000, -713) & (1000, -714)`).
4. **Energy-scale-free flat check + exact `0/0` and high-SNR two-body `lump`:** Tiny non-flat signals (`(3e-10, 1e-10)` at unit noise) are no longer falsely flagged flat; `lump(body, 0/0)` matches `body`'s 1-dof lower endpoint to 8+ digits; analytic local variance-weighted splits resolve `SNR = 1e6` two-body lumping to `< 1e-3` relative accuracy.
5. **Open-ray $(0, \infty)$ boundary pinhole limit (`L1`):** On aligned one-sided readings like `(10, 0)`, $m = 0$ is an isolated non-positive $C = 0$ pinhole surrounded by $p = 1.5\times 10^{-23}$. Evaluating degenerate endpoints `0` and `∞` via their one-sided positive-ray limit (`x = 1e-50` in `kept_intervals`, `x = 1e-30` in `lump`) keeps `(10, 0)` cleanly `at_least`, `(0, 10)` cleanly `at_most`, and `lump((10, 0), (10, 0))` cleanly `at_least` without fake zero-width boundary slivers.

---

## 3. Premise & Adversarial Scorecard (`tests/test_ours35.py`)

| Premise / Scene | `ours35` (Door A1) Verdict | Verified Evidence (`tests/test_ours35.py`) |
| --- | --- | --- |
| Full-type agreement with `met.ours34` (`met.evidence` / `met.conditional_score`) | **PASS** | 20 random mixed `d=1,2,3` correlated-$\Sigma_i$ subspace-$N_i$ series: `K`, `Z`, `argmin_K` agree to `1e-12`–`1e-15` |
| `T1` exact branch ($\sigma = 0$) | **PASS** | `exact_mass` & `exact_series` handle one mass, every mass, anti-aligned, non-codirectional, and inconsistent pairs |
| `L2` input validation | **PASS** | Non-SPD, asymmetric, non-orthonormal $N_i$, NaN/inf, invalid level/scale all raise `ValueError` |
| `L3` unit equivariance (including $s = 10^{\pm 200}$) | **PASS** | Exact 12–15 digit equivariance across $\lambda = 13.7$ and a $10^{400}$ ratio span |
| `T5` channel swap ($F \leftrightarrow a$, $m \leftrightarrow 1/m$) | **PASS** | Exact reciprocation of point & set, invariance of `leftover` (including off-grid pinholes) |
| `L4` spatial rotation invariance | **PASS** | Exact invariance on correlated 3D reading with $k = 2$ subspace |
| `L16` order invariance & `L18` window indifference | **PASS** | Identical to `1e-12` |
| `L20` exact `Z` (`0/0` identity) | **PASS** | Adding 50 exact `0/0` readings (any $\Sigma$) leaves `point`, `set`, `leftover`, `p_value` bit-for-bit identical |
| **Continuity at `0/0`** (kills `M_nz`-`GCLR` & `ours34` `t_joint`) | **PASS** | Adding 50×`(1e-9, 1e-9)` moves `p`, `point`, `set`, `leftover` by `< 1e-11` |
| **`L17` reading splitting** (kills `M_nz`-`GCLR`) | **PASS** | Splitting 1 reading into 10 or 50 equal parts keeps `point`, `set`, `leftover`, `p` identical to `1e-11` |
| `T77 (6)` pure additive `(K, C)` monoid | **PASS** | Associative `s1 + s2`, `0/0` identity, no `M` or `M_nz` anywhere on `State` |
| Uniform 95% exact frequentist coverage (full type) | **PASS** | Correlated 3D ($k=2$) + correlated 1D mixed-instrument series across null (`0`), weak (`0.7`), strong (`5.0`) + `0/0`: all within `[0.925, 0.975]` |
| `G29` right-angle `(10, 0) & (0, 10)` | **PASS** | Flat `K ≡ 100`, `point_is_convention=True`, `set = (0, ∞)`, `leftover = 100.0` |
| `L5` / `L14` two-body composition + `0/0` + asymmetric `100+10` | **PASS** | Exact `T81` point addition/harmonic dual, exact 1-dof collapse on `0/0`, ghost-free on `100+10`, exact 1-dof `lump_known`/`common_force_known`, coverage $\ge 95\%$ |
| **Door A cost #1:** No set widening on disagreement | **KNOWN COST (by design / `T82`)** | Can's `near (5.5, 5.6)` (`set = (5.4114, 5.6886)`, `leftover = 0.5`) vs `far (5, 6)` (`set = (5.3614, 5.6386)`, `leftover = 50.0`): identical width `0.2772`, `100×` leftover jump; `(100,10)&(10,100)` gives `(0.897, 1.114)` with `leftover = 8100` |
| **Door A cost #2:** Hilltop ghosts, pinholes, & anti-aligned boundary rejection | **KNOWN COST (`E15-1`, `E15-3`)** | On `(10,-3)&(10,-3.1)`: `T81` point is `∞` (`leftover = 18.61`), `T45` rejects `∞` ($p = 1.6\times 10^{-5}$) and keeps only the hilltop at `m = 0.305` ($K = 218.61$); aligned `(10, 0)` keeps `∞` (`p = 1.0`, `at_least`) |
| **Door A cost #3:** 2-dof width tax in curve-pair two-body `lump` / `common_force` | **KNOWN COST (`G27`)** | $+24.9\%$ width ($2.45\sigma$ vs $1.96\sigma$) at strong signal when both bodies are informative (verified analytically at `SNR = 1e6`); zero tax on `lump_known` or `body + 0/0` |

---

## 4. How to Run

```bash
python3 -m unittest -v tests/test_ours35.py
```

And on an environment with `numpy` + `scipy` + `pytest` (e.g., Claude's harvest check on `main`), `met.ours35.State` and `met.ours35.ours35` accept `met.full.Reading` directly via `adapt_reading`, so it drops straight into a `test_ours35_against_main.py` cross-check against `met.evidence` and `met.conditional_score`.
