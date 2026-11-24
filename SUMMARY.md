# E16: `met/ours35.py` — Executable Door A1 Probe & Premise Scorecard

Branch: `explore/e16-ours35-door-a` (from `origin/main` `62e0390`, 24 Nov 2026). Exploration only; nothing in `definition/`, `problems/`, `rulings.md`, or existing `met/` files modified.

**Core target:** Build and stress-test `met/ours35.py`, the zero-knob executable implementation of **Door A1** (`explore/e15-door-a-uniqueness`): believed declared Gaussian $\Sigma_i$ (`L13`), no starting law $\pi$ (`T80`), no $\nu_{\rm eff}$ or set-widening discount, pure additive state `State = (K, C)` with zero hidden reading count `M` or `M_nz` (`T77 (6)`), point `argmin K` on `[0, inf]` (`T81`), exact `T45` score set (`E15-1`/`E15-2`), pure `K_min` leftover (`T73 (5)`), exact zero-noise branch (`T1`), and curve-pair `L5`/`L14` two-body composition (`T74` + `Q_joint`, plus exact 1-dof `lump_known`/`common_force_known`).

## What passed (and what numerical defects over `ours34` were caught & fixed)

1. **Full-type cross-check & exact coverage:** Matches `met/ours34.py` (already cross-checked against `main`'s `met.evidence` & `met.conditional_score`) to `1e-12`–`1e-15` across 20 random mixed `d=1,2,3` correlated-$\Sigma_i$ subspace-$N_i$ series; holds exact 95% frequentist coverage across null, weak, and strong signal on correlated 3D (`k=2`) + 1D mixed-instrument series with `0/0` contamination.
2. **All Door-A invariances hold on every output field:** Exact `Z` (`L20`), **continuity at `0/0`** (immune to `50×1e-9` blowup that breaks `M_nz`-`GCLR` and `ours34`'s `t_joint`), **exact `L17` reading splitting** (`1` vs `50` equal parts), `L3` unit equivariance across a **$10^{400}$ scale span** ($s = 10^{\pm 200}$, via per-reading scale-normalized pre-whitening), exact `T5` channel swap (via swap-symmetric two-half coordinates), `L4` 3D rotation, `L16` order invariance, and `L18` window indifference.
3. **Adversarial & numerical fixes over `ours34`:** Exact `S = 0` stationary-point discovery resolves `SNR = 1e10` narrow intervals (`width = 1.24e-9` exact), ultra-narrow high-SNR hilltops, and off-grid single-reading pinholes (`(10, -7.13)` now `L3`/`T5`-equivariant); energy-scale-free flat check stops tiny signals (`1e-10`) being falsely flagged flat; open-ray $(0, \infty)$ boundary pinhole limit keeps `(10, 0)` cleanly `at_least` and `(0, 10)` cleanly `at_most`; `lump(body, 0/0)` matches the 1-dof bound to 8+ digits; `lump` resolves `SNR = 1e6` and ghost-free asymmetric `100 + 10`.

## Known provable Door A1 costs (theorems, not bugs; surfaced explicitly in tests)

- **No set widening on disagreement (Door A / `T82`):** Can's `near (5.5, 5.6)` and `far (5, 6)` keep identical width `0.2772` while `leftover` jumps `0.5 → 50.0`; `(100,10)&(10,100)` keeps `(0.897, 1.114)` with `leftover = 8100`.
- **Hilltop ghosts & anti-aligned boundary rejection (`E15-1`, `E15-3`):** On `(10,-3)&(10,-3.1)`, `T81` point is `inf` (`leftover = 18.61`), while `T45` rejects `inf` ($p = 1.6\times 10^{-5}$) and keeps only the hilltop at `m = 0.305` ($K = 218.61$).
- **2-dof width tax in curve-pair two-body `lump` / `common_force` (`G27`):** $+24.9\%$ at strong signal (`2.45σ` vs `1.96σ`; zero tax on `lump_known` or `body + 0/0`).

## Files & verification

- Pure-stdlib full-type Door A1 probe: `met/ours35.py` (duck-typed to accept `met.full.Reading`)
- Premise scorecard & design notes: `explore/e16-ours35-door-a.md`
- 12 executable checks: `python3 -m unittest -v tests/test_ours35.py`
