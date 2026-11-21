# E13: Unnormalized Evidence ("Bell with σ = ∞") and the Five Insights

**Branch:** `explore/e13-unnormalized-evidence` (from `main` `d2a3352`, 21 Nov 2026).  
**Status:** exploration input for Claude to harvest (`explore/README.md`, R-2026-11-07-2). Nothing in `definition/`, `problems/`, `rulings.md`, or `met/` is modified.  
**Code & checks:** `explore/e13_unnormalized_evidence.py`, `tests/test_e13_unnormalized_evidence.py` (pure stdlib; runs with `python3 -m unittest tests/test_e13_unnormalized_evidence.py` or `pytest`).

## Standard mathematics first (GOAL.md standing rule 6)

- **Kac–Bernstein theorem** (Kac 1939; Bernstein 1941): for independent $X_1, X_2$, the sample mean $\bar X$ is independent of the difference $X_1 - X_2$ *if and only if* $X_1, X_2$ are Gaussian.
- **Unnormalized / profile likelihood & infimal convolution** (Royall 1997; Barndorff-Nielsen; convex/epigraphical addition): combining independent evidence adds negative log-likelihoods $K_1 + K_2$; combining two parameters under $M = m_1 + m_2$ profiles out the split via $K_{\rm lump}(M) = \inf_{m_1+m_2=M}[K_1(m_1)+K_2(m_2)]$.
- **Weak-identification exact test inversion** (Moreira 2003, *Econometrica* 71:1027–1048, CLR and conditional score `T45`; Gleser & Hwang 1987; Dufour 1997): confidence sets with uniform finite-sample coverage near non-identification must allow unbounded sets $(0,\infty)$ with positive probability.

---

## The 5 Insights (countable summary for harvest)

### Insight 1 — `0/0` forces the identity of combination ($L \equiv 1$, "bell with $\sigma = \infty$"), *not* a normalized prior $\pi_c$ (`E5`, `O6`, `O11`, `O2`, `G28`, `B16`)

In `docs/vd-docs-context/why_mass_cannot_be_undefined_2026-09-29.md` §2.3, the density-route contradiction proves that a stationary `0/0` trace must act as the **identity element of combination** (combining `0/0` with `5 ± 0.5 kg` must return `5 ± 0.5 kg` untouched). The leap to §2.4 demanded a **proper area-1 probability law** $\pi_c$ and a finite point $c = \sigma_F/\sigma_a$ even when no data exist.

- An unnormalized Gaussian bell $e^{-(x-c)^2/(2\sigma^2)}$ as $\sigma \to \infty$ is literally $L(m) \equiv 1$ ($K(m) \equiv 0$) — its center $c$ vanishes and its inverse-variance combination weight is $1/\infty^2 = 0$ (matching R-2026-09-29-3).
- `B16` already proved that `T36`'s residual likelihood never drops to 0 as $m \to 0$ or $m \to \infty$ on any finite dataset, so its area is *always* infinite without a finite-width prior $\pi$.
- **Consequence:** Dropping the requirement to normalize area to 1 (returning the unnormalized residual curve $K(m)$ + a horizontal compatibility/confidence set that may be $(0,\infty)$) dissolves `O6` (shape of $\pi$), `O11` (multi-instrument center), `G6`/`O2` (forced point on `0/0`), and the Gleser–Hwang / null-floor impossibility (`O15`, `G28`), while keeping `0/0` total and never `undefined`. Any caller who wants a Bayesian posterior multiplies $e^{-K/2}$ by their own declared prior $\pi$ on top.

### Insight 2 — Exact Gaussian `L13` + believed $\sigma$ (`L18`) and discordance widening (`R-2026-11-17-3` / `T30`) contradict each other by Kac–Bernstein (`O14`, `G13`, `G18`, `G19`, `G22`)

Under independent exact Gaussian errors (`L13`), the mean $\bar X$ and gap $X_1 - X_2$ are statistically independent (Kac–Bernstein, verified in `tests/test_believed_less.py`). Therefore, if exact Gaussianity and exact $\sigma$ are truly believed (`R-2026-11-17-1`), Scene B (`5.0, 6.0 ± 0.1`) genuinely has the same conditional precision for $m$ as Scene A (`5.5, 5.6 ± 0.1`): `5.5 ± 0.14`. Human intuition expects widening because real errors have heavy tails (occasional flukes) or off-graph departures — which `explore/e12-scalar-mass-proof` (`e12-adding-gaussian-noise.md`) independently confirmed via Gaussian noise-degradation towers.

This splits `O14` cleanly into **two contradiction-free doors** (no ad-hoc $g$ or $\eta$):
- **Door A (Strict consequence of exact Gaussian inputs):** keep exact `T36` residual $K(m)$ and exact compatibility set (`T45` / Moreira CLR), drop automatic widening inside the mass law, and report the leftover $K_{\min} = \min_{m>0} K(m)$ (with exact upper-bound surprise $p$-value $P(\chi^2_M \ge K_{\min})$) beside the mass. Virtually all of `O14`'s holes vanish and the object is essentially unique.
- **Door B (Heavy-tailed / fluke-allowing error model):** drop exact Gaussian `L13` as universal and declare the error-tail shape (e.g., `T48` Student-$t$ with declared $\nu$) as an instrument input alongside $\sigma$ — evaluated *unnormalized* (no $\pi$) and **without $\eta$**.

### Insight 3 — Exact-zero `Z` (`y = 0` at finite $\sigma$) vs. truly uninformative data ($\sigma \to \infty$), and `T45` vs. Moreira CLR (`L20`, `O7`, `T45`, `T69`, `B24`)

An observed `(F, a) = (0, 0)` at finite $\sigma = 1$ carries real amplitude/noise-count information; a truly uninformative reading is $\sigma \to \infty$ ($\Sigma^{-1} \to 0$). Treating `y = 0` at finite $\sigma$ as `Z` blocks count-bearing exact tests like Moreira's CLR (`B24`) and maximal-invariant likelihood (`T69`).

On the killer gallery (table below), both `T45` and `Moreira CLR` achieve ~95% coverage with **zero $\pi$ and zero $\eta$** across pure noise ($m = 20$), weak-signal 3D $N=20$, `T43` leverage, and strong signal, and both handle right-angle Scene C (`G29`) cleanly by returning `(0, inf)` with $K_{\min} = 400$. Their exact trade-off shows up cleanly on the **anti-parallel pair `(2,-1), (1,-2)`**:
- `T45` passes exact-zero `Z`, **but** accepts a ghost interval `[0.954, 1.048]` around the *worst* positive mass (the hilltop where the slope vanishes).
- `Moreira CLR` eliminates the hilltop ghost (`EMPTY` + $K_{\min} = 500$, $p \le 2.7 \times 10^{-109}$), at the sole cost of using the noise-dimension count $M = \sum k_i$ (failing exact-zero `Z`, while passing $\sigma \to \infty$).

### Insight 4 — Prior-free two-body composition works by infimal convolution (`L5`, `L14`, `B26`)

Without area-1 normalization, combining two independent bodies' residual curves $K_1, K_2$ under lumping (`L5`, $M = m_1 + m_2$) or common force (`L14`, $1/M = 1/m_1 + 1/m_2$) uses **infimal convolution**:
$$
K_{\rm lump}(M) = \inf_{m_1 + m_2 = M} \bigl[K_1(m_1) + K_2(m_2)\bigr],
\qquad
K_{\rm red}(M) = \inf_{1/m_1 + 1/m_2 = 1/M} \bigl[K_1(m_1) + K_2(m_2)\bigr].
$$
This is exact on noiseless inputs, reproduces exact Gaussian variance addition ($\sigma_1^2 + \sigma_2^2$), satisfies exact swap duality between `L5` and `L14`, is immune to `B26`'s non-closure obstruction (no prior $\pi$ to close), and maps `5.0 ± 0.5 kg` glued to an unmeasured stationary body (`0/0`, $K_2 \equiv 0$) to `[4.05, inf)` (ruling out $< 4.02$ kg while leaving $> 5$ kg completely open).

### Insight 5 — Unnormalized Door B (`T48` Student-$t$, no $\pi$, no $\eta$) fixes both `ours33` breaks on discordant readings (`O18`, `O14`)

Evaluating `ours33`'s Student-$t$ residual factor (`T48`, $\nu = 10$) unnormalized (likelihood-ratio slice $\Delta(-2\ell) \le 3.841$, no $\pi$, no $\eta$):
- returns `(0, inf)` on `0/0` and on 3D perpendicular singletons;
- gives `[4.915, 6.090]` on Scene B at 10σ (eliminating `ours33`'s 79-decade double-discount blowup);
- does **not** crash at 50σ (where `ours33` raises `ValueError`), instead splitting cleanly into `[4.968, 5.041] ∪ [5.963, 6.035]` (*"one reading is a fluke: mass is near 5 or near 6"*).

---

## Reproducible Killer-Gallery Results (`python3 explore/e13_unnormalized_evidence.py`)

### Part A — Deterministic scenes (zero $\pi$, zero $\eta$)

| Scene | Peak $\hat m$ | $K_{\min}$ (surprise $p \le$) | `T45` 95% set | `Moreira CLR` 95% set | Door B ($t_{10}$) 95% set |
|---|---|---|---|---|---|
| 1. All-zero `0/0` (3D, $\sigma_F=2, \sigma_a=1$) | none (flat) | `0.0` (`1`) | `(0, inf)` | `(0, inf)` | `(0, inf)` |
| 2. Aligned (`F=10±1, a=5±1`) | `2.000` | `0.0` (`1`) | `[1.345, 3.379]` | `[1.345, 3.379]` | `[1.324, 3.478]` |
| 3a. Scene A (`5.5, 5.6 ± 0.1`) | `5.550` | `0.494` (`0.781`) | `[5.412, 5.688]` | `[5.412, 5.688]` | `[5.407, 5.694]` |
| 3b. Scene B (`5.0, 6.0 ± 0.1`, 10σ) | `5.501` | `49.40` (`1.87e-11`) | `[5.361, 5.639]` | `[5.361, 5.639]` | `[4.915, 6.09]` |
| 3c. Scene B extreme (`5.0, 6.0 ± 0.02`, 50σ) | `5.501` | `1235.1` (`6.48e-269`) | `[5.473, 5.528]` | `[5.473, 5.528]` | `[4.968, 5.041] U [5.963, 6.035]` |
| 4. Right-angle Scene C (`G29`, $\sigma=0.05$) | none (flat) | `400.0` (`1.38e-87`) | `(0, inf)` | `(0, inf)` | `(0, 0.1041) U (9.605, inf)` |
| 5. Anti-parallel pair (`(2,-1), (1,-2)`, $\sigma=0.1$) | boundary | `500.0` (`2.67e-109`) | `[0.9539, 1.048]` *(ghost)* | `EMPTY` | `(0, 0.08297) U (12.05, inf)` |
| 6. E12 triple (`F=0.2, 1.0, 1.4`, $\sigma=0.01$) | `1.000` | `4000.0` (`0`) | `[0.9823, 1.018]` | `[0.9839, 1.016]` | `[0.9777, 1.037] U [1.357, 1.425]` |
| 7. 3D perp singleton (`F=e_1, a=e_2`, $\sigma=0.1$) | none (flat) | `100.0` (`1.55e-21`) | `(0, inf)` | `(0, inf)` | `(0, inf)` |

### Part B — Monte Carlo 95% coverage on graveyard worlds (500 reps, seed 2026)

| World | `T45` cov (unb.) | `Moreira CLR` cov (unb.) | Raw Gaussian LR *(contrast)* |
|---|---|---|---|
| Pure noise (`3D, N=20, signal=0, m=20`) | `93.2%` (`93.6%`) | `96.2%` (`96.4%`) | `49.6%` |
| Weak signal (`3D, N=20, signal=1, m=2`) | `94.2%` (`83.2%`) | `96.2%` (`83.0%`) | `62.8%` |
| `T43` leverage (`1D, 1x r=20 + 99 nulls, m=2`) | `95.6%` (`0.0%`) | `95.6%` (`0.0%`) | `93.0%` |
| Strong signal (`3D, N=10, signal=5, m=2`) | `94.0%` (`0.0%`) | `94.0%` (`0.0%`) | `92.0%` |

---

## What this exploration does NOT solve (open for harvest / next steps)

1. **General-covariance CLR:** `moreira_clr_pvalue` here assumes the common-ratio isotropic slice (one noise ratio across readings); generalizing exact hilltop-free CLR to arbitrary heterogeneous $\Sigma_i$ and subspaces $N_i$ is open (`T45` already works on arbitrary $\Sigma_i, N_i$).
2. **Weak-signal calibration of composed curves:** Infimal convolution of raw profile $K$ covered 91.5% on `M = 2 + 3 = 5` (3D, $N=5$ each, signal=3, 400 reps) — inheriting raw profile LR's mild Neyman–Scott gap; calibrating composed sets at weak signal (and calibrating Door B Student-$t$ sets at weak signal) remains open.
3. **Random pure noise beside an informative block:** `explore/e2-null-combination`'s obstruction still applies to random noisy zero-signal blocks (distinct from exact observed `0/0`).
4. **The attitude fork (Can / Claude):** Door A (strict Gaussian consequences + separate $K_{\min}$ surprise) vs. Door B (declared heavy-tailed error model) is a choice about the goal, not a calculation.
