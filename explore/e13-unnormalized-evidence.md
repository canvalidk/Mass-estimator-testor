# E13: Unnormalized Evidence ("Bell with σ = ∞") — Eight Insights

**Branch:** `explore/e13-unnormalized-evidence` (from `main` `d2a3352`, 21 Nov 2026).  
**Status:** exploration input for Claude to harvest (`explore/README.md`, R-2026-11-07-2). Nothing in `definition/`, `problems/`, `rulings.md`, or `met/` is modified.  
**Code & checks (pure Python stdlib — zero external dependencies):**
- Initial prototype & killer gallery (Insights 1–5): `explore/e13_unnormalized_evidence.py`, `tests/test_e13_unnormalized_evidence.py`
- Push 1–3 general-type CLR, exact two-body composition, and Door B / spherical CLR / O7 floor union (Insights 6–8): `explore/e13_general_clr.py`, `tests/test_e13_general_clr.py`
- Run all 15 checks: `python3 -m unittest -v tests/test_e13_unnormalized_evidence.py tests/test_e13_general_clr.py` (also pytest-discoverable)

## Standard mathematics first (`GOAL.md` standing rule 6)

- **Kac–Bernstein theorem** (Kac 1939; Bernstein 1941): for independent $X_1, X_2$, the sample mean $\bar X$ is independent of the difference $X_1 - X_2$ *if and only if* $X_1, X_2$ are Gaussian.
- **Unnormalized / profile likelihood & infimal convolution** (Royall 1997; Barndorff-Nielsen; epigraphical addition): combining independent evidence adds negative log-likelihoods $K_1 + K_2$; combining two parameters under $M = m_1 + m_2$ profiles out the split via $K_{\rm lump}(M) = \inf_{m_1+m_2=M}[K_1(m_1)+K_2(m_2)]$, with Chernoff (1954) one-sided boundary mixture $\frac{1}{2}\delta_0 + \frac{1}{2}\chi^2_1$ when one body is unmeasured (`0/0`) under $m_2 > 0$.
- **Weak-identification exact conditional test inversion** (Moreira 2003, *Econometrica* 71:1027–1048; Andrews–Moreira–Stock 2006; Gleser & Hwang 1987; Dufour 1997; Dwass 1957 / Barnard 1963 exact Monte Carlo tests): although no single nuisance-free *conditional likelihood* exists because $U_i(m_0)$ depends on $m_0$ (Kalbfleisch & Sprott 1970, cited in `problems/O15.md`), **pointwise conditional test inversion** given $(U_1(m_0), \dots, U_N(m_0))$ is exact in finite samples on the full covariance/subspace type without a prior $\pi$ or warp $\eta$, provided unbounded sets $(0,\infty)$ are allowed.
- **Spherical / sign-flip randomization tests** (Fisher 1935; Lehmann & Romano 2005, §15.2): conditioning on both $U_i(m_0)$ and $\|V_i(m_0)\|$ makes $V_i / \|V_i\| \sim \mathrm{Unif}(S^{k_i-1})$ exact under *any* per-reading noise-scale mixture $v_i \Sigma_i$ (Gaussian, Student-$t$, Cauchy, adversarial scale spikes).

---

## The 8 Insights (countable summary for harvest)

### Insight 1 — `0/0` forces the identity of combination ($L \equiv 1$, "bell with $\sigma = \infty$"), *not* a normalized prior $\pi_c$ (`E5`, `O6`, `O11`, `O2`, `G28`, `B16`)

In `docs/vd-docs-context/why_mass_cannot_be_undefined_2026-09-29.md` §2.3, the density-route contradiction proves that a stationary `0/0` trace must act as the **identity element of combination** (combining `0/0` with `5 ± 0.5 kg` must return `5 ± 0.5 kg` untouched). The leap to §2.4 demanded a **proper area-1 probability law** $\pi_c$ and a finite point $c = \sigma_F/\sigma_a$ even when no data exist.

- An unnormalized Gaussian bell $e^{-(x-c)^2/(2\sigma^2)}$ as $\sigma \to \infty$ is literally $L(m) \equiv 1$ ($K(m) \equiv 0$) — its center $c$ vanishes and its inverse-variance combination weight is $1/\infty^2 = 0$ (matching R-2026-09-29-3).
- `B16` already proved that `T36`'s residual likelihood never drops to 0 as $m \to 0$ or $m \to \infty$ on any finite dataset, so its area is *always* infinite without a finite-width prior $\pi$.
- **Consequence:** Dropping the requirement to normalize area to 1 (returning the unnormalized residual curve $K(m)$ + a horizontal compatibility/confidence set that may be $(0,\infty)$) dissolves `O6` (shape of $\pi$), `O11` (multi-instrument center), `G6`/`O2` (forced point on `0/0`), and the Gleser–Hwang / null-floor impossibility (`O15`, `G28`), while keeping `0/0` total and never `undefined`. Any caller who wants a Bayesian posterior multiplies $e^{-K/2}$ by their own declared prior $\pi$ on top.

### Insight 2 — Exact Gaussian `L13` + believed $\sigma$ (`L18`) and discordance widening (`R-2026-11-17-3` / `T30`) contradict each other by Kac–Bernstein (`O14`, `G13`, `G18`, `G19`, `G22`)

Under independent exact Gaussian errors (`L13`), the mean $\bar X$ and gap $X_1 - X_2$ are statistically independent (Kac–Bernstein, verified in `tests/test_believed_less.py`). Therefore, if exact Gaussianity and exact $\sigma$ are truly believed (`R-2026-11-17-1`), Scene B (`5.0, 6.0 ± 0.1`) genuinely has the same conditional precision for $m$ as Scene A (`5.5, 5.6 ± 0.1`): `5.5 ± 0.14`. Human intuition expects widening because real errors have heavy tails (occasional flukes) or off-graph departures — which `explore/e12-scalar-mass-proof` (`e12-adding-gaussian-noise.md`) independently confirmed via Gaussian noise-degradation towers.

This splits `O14` cleanly into **two contradiction-free doors** (no ad-hoc $g$ or $\eta$):
- **Door A (Strict consequence of exact Gaussian inputs):** keep exact `T36` residual $K(m)$ and exact compatibility set (`T45` / `GCLR` / `Z-GCLR`), drop automatic widening inside the mass law, and report the leftover $K_{\min} = \min_{m>0} K(m)$ (with exact upper-bound surprise $p$-value $P(\chi^2_M \ge K_{\min})$) beside the mass.
- **Door B (Heavy-tailed / scale-uncertain error model):** drop exact fixed-scale Gaussian `L13` as universal — either declare a heavy-tailed law (`T48` Student-$t$, evaluated unnormalized without $\eta$, though see Insight 8 for its weak-signal break) or condition out per-reading noise scales via spherical `SCLR+` / the `O7` floor union (Insight 8).

### Insight 3 — Exact-zero `Z` (`y = 0` at finite $\sigma$) equals $\sigma \to \infty$ in whitened coordinates, and `T45` vs. Moreira CLR (`L20`, `O7`, `T45`, `T69`, `B24`)

In whitened coordinates $z_i = \Sigma_i^{-1/2} y_i$, an exact zero observation $(F_i, a_i) = (0, 0)$ at finite $\sigma$ and a fixed nonzero observation at $\sigma \to \infty$ both map to $z_i = 0$.
- The **unnormalized residual curve $K(m)$** itself and the **`T45` score test** pass both exact-zero `Z` and $\sigma \to \infty$.
- **Unweighted `Moreira CLR`** (and unweighted `GCLR` below) **fail both** exact-zero `Z` *and* $\sigma \to \infty$ (because a reading with $z_i = 0$ still contributes $k_i$ degrees of freedom to the simulated null draws $V_i^{(b)} \sim N(0, I_{k_i})$ and the noise count $M = \sum k_i$).
- On the common-ratio killer gallery (Table A1 below), both `T45` and `Moreira CLR` achieve ~95% coverage with zero $\pi$ and zero $\eta$ across pure noise ($m = 20$), weak-signal 3D $N=20$, `T43` leverage, and strong signal, and both return `(0, inf)` on right-angle Scene C (`G29`). On the anti-parallel pair `(2,-1), (1,-2)`, `T45` accepts a hilltop ghost `[0.954, 1.048]` around the *worst* positive mass, while `Moreira CLR` returns `EMPTY` ($K_{\min} = 500$, $p \le 2.7 \times 10^{-109}$).

### Insight 4 — Prior-free two-body composition works by infimal convolution (`L5`, `L14`, `B26`)

Without area-1 normalization, combining two independent bodies' residual curves $K_1, K_2$ under lumping (`L5`, $M = m_1 + m_2$) or common force (`L14`, $1/M = 1/m_1 + 1/m_2$) uses **infimal convolution**:
$$
K_{\rm lump}(M) = \inf_{m_1 + m_2 = M} \bigl[K_1(m_1) + K_2(m_2)\bigr],
\qquad
K_{\rm red}(M) = \inf_{1/m_1 + 1/m_2 = 1/M} \bigl[K_1(m_1) + K_2(m_2)\bigr].
$$
This is exact on noiseless inputs, reproduces exact Gaussian variance addition ($\sigma_1^2 + \sigma_2^2$), satisfies exact swap duality between `L5` and `L14`, is immune to `B26`'s non-closure obstruction (no prior $\pi$ to close), and maps `5.0 ± 0.5 kg` glued to an unmeasured stationary body (`0/0`, $K_2 \equiv 0$) to `[4.02, inf)` (ruling out $< 4.02$ kg while leaving $> 5$ kg completely open).

### Insight 5 — Unnormalized Door B (`T48` Student-$t$, no $\pi$, no $\eta$) fixes `ours33`'s discordant blowup and 50σ crash (`O18`, `O14`)

Evaluating `ours33`'s Student-$t$ residual factor (`T48`, $\nu = 10$) unnormalized (likelihood-ratio slice $\Delta(-2\ell) \le 3.841$, no $\pi$, no $\eta$):
- returns `(0, inf)` on `0/0` and on 3D perpendicular singletons;
- gives `[4.915, 6.090]` on Scene B at 10σ (eliminating `ours33`'s 79-decade double-discount blowup);
- does **not** crash at 50σ (where `ours33` raises `ValueError`), instead splitting cleanly into `[4.968, 5.041] ∪ [5.963, 6.035]` (*"one reading is a fluke: mass is near 5 or near 6"*). *(See Insight 8 for its weak-signal coverage break.)*

### Insight 6 (Push 1) — Full-type General Conditional LR (`GCLR` / `GCLR+`) and `Z`-safe `U`-weighted `Z-GCLR` extend `B24` past the common-ratio slice (`B24`, `T45`, `E2`, `E3`, `L20`, `T43`)

`problems/O15.md` (`B24`, 10 Nov) noted that heterogeneous covariances $\Sigma_i$ remove the common two-column rotation so *"the CLR formula below is not thereby extended."* Conditioning pointwise at candidate $\phi_0 = \arctan(m_0)$ on the whitened along-line projections $U_i = Q_{\parallel,i}(\phi_0)^\top z_i \in \mathbb{R}^{k_i}$ (inside each reading's $2k_i$-dimensional whitened span $W_i$) makes the orthogonal residuals $V_i = Q_{\perp,i}(\phi_0)^\top z_i \sim N(0, I_{k_i})$ **exact, independent across readings, independent of $U_i$, and free of all unknown pushes $\alpha_i^*$** on the **entire full type** (arbitrary correlated $\Sigma_i$, heterogeneous ratios, subspaces $N_i$). Drawing $V_i^{(b)} \sim N(0, I_{k_i})$ and evaluating the profile likelihood ratio on a fixed grid yields an exact finite-sample Dwass–Barnard conditional test (`GCLR+` on $m>0$, `GCLR` on $\mathbb{RP}^1$):
- **Full-type exact coverage (Table B2):** on correlated $d=3, k=2, \rho=0.6$ mixed-ratio worlds (`N=10`), raw profile LR collapses to `44.0%` (pure noise) and `67.0%` (weak signal), while `GCLR+`, `GCLR(RP^1)`, and `Z-GCLR+` hold exact ~95% (`95.0%–97.0%`) with zero $\pi$ and zero $\eta$.
- **Exact-`Z`-safe `U`-weighted `Z-GCLR` (`w_i = H_i(m_0)/(1 + H_i(m_0))`):** because $H_i(m_0) = \|U_i\|^2$ is fixed by conditioning, any invariant weight $w(H_i)$ with $w(0)=0, w(\infty)=1$ preserves **exact finite-sample 95% coverage** (`95.0%, 95.5%, 95.0%`), passes **units, swap, rotation, exact-zero `Z` (bit-for-bit!), and $\sigma \to \infty$**, eliminates the multi-direction anti-parallel hilltop ghost `(2,-1), (1,-2)` (`EMPTY` on $\mathbb{RP}^1$, boundary on $m>0$), returns `(0, inf)` on right-angle Scene C (`G29`), and **sharpens `T43` leverage width** (`[1.65, 2.65]`, width `1.00`, vs `1.10` for `T45` and unweighted `GCLR+`) by downweighting pure-noise readings.
- **Sharp theorem-forced trade-off (`E2`, `E3`):**
  1. $w(H)$ is a **family** of Z-safe weights (e.g., $H/(1+H)$ vs $1 - e^{-H/2}$), directly answering `E2`/`E3`'s question *"one member or a family?"* — nothing in the Gaussian type selects a unique member.
  2. **Rank-1 pinhole ghost:** under $U$-conditioning alone, exact `Z`-invariance forces ignoring any reading with $U_i(m_0) = 0$. Therefore on a *rank-1* anti-parallel dataset like a single reading `(2, -1)` at its orthogonal mass $m_0 = 0.5$ (where *every* reading has $U_i(0.5) = 0$ despite $K = 500$), `Z-GCLR` must give $p = 1.0$ (a pinhole ghost), whereas unweighted `GCLR` rejects even that pinhole ($p = 0.005$) at the price of failing both `Z` and $\sigma \to \infty$.

### Insight 7 (Push 2) — Exact weak-signal two-body composition (`L5` lumping and `L14` common force) lifts coverage to 95% and hits the exact Chernoff boundary (`L5`, `L14`, `B26`)

Raw infimal convolution (`Insight 4`) inherits raw profile LR's Neyman–Scott gap at weak signal (`89.3%` on 3D `M = 2 + 3 = 5`, $N=5$ each, signal=1). Conditioning on $(U_1(m_1), U_2(M_0 - m_1))$ calibrates the profiled lump LR $\mathrm{LR}_{\rm lump}(M_0) = \inf_{m_1+m_2=M_0}[K_1(m_1)+K_2(m_2)] - K_{1,\min} - K_{2,\min}$ exactly in finite samples (taking $\sup_{0 < m_1 < M_0}$ for a guaranteed $\ge 95\%$ confidence set for $M_0$):
- **Weak-signal calibration fixed (Table B3):** lifts weak-signal 3D lumping coverage from `89.3%` (`Raw infimal-conv`) to `97.3%` (unweighted) and **`96.0%` (`Z-Cond Lump LR`)**, with `95.3%` at moderate and strong signal.
- **1-df sharpness preserved at strong signal:** on strong `2 + 3 = 5` (combined $\sigma_M = 0.10$), $\sup_{m_1} p(m_1, M_0 - m_1)$ recovers the exact 1-df interval `[4.80, 5.20]` ($p(4.76)=0.018$, $p(4.80)=0.060$, $p(5.20)=0.068$, $p(5.24)=0.015$) with zero 2-df Bonferroni widening.
- **Exact `L14` common-force duality:** swapping $(F, a)$ turns `L14` ($1/M = 1/m_1 + 1/m_2$) literally into `L5` on $1/M$ (bit-identical p-values).
- **Chernoff one-sided boundary on `(5.0 ± 0.5 kg) + (0/0)`:** unweighted conditional lump LR lets a 3D `0/0` body inject 3 noise dims into the null distribution, over-widening the lower bound. With `Z`-safe weights (`z_smooth`), the `0/0` body drops out bit-for-bit and the test automatically recovers the exact **Chernoff $\frac{1}{2}\delta_0 + \frac{1}{2}\chi^2_1$ one-sided 95% bound** `[4.18, inf)` ($5.0 - 1.645 \times 0.5$: $p(4.00)=0.020$, $p(4.18)=0.051$, $p(4.35)=0.100$, $p(5.0)=p(50.0)=1.000$).
- **Honest trade-off:** exact weak-signal composition needs the readings (or their $2k$-dim within-span summaries), not merely the 1D curves $K_1(m), K_2(m)$.

### Insight 8 (Push 3) — Door B naive Wilks breaks at weak signal; Spherical Conditional LR (`SCLR+`) and the `O7` floor union deliver exact scale-free coverage (`O14`, `O7`, `T37`, `T48`, `O18`)

Stress-testing Door B (`Insight 5`, unnormalized `T48` Student-$t_{10}$ with Wilks `3.841` cutoff) on Monte Carlo worlds (`3D, N=15, 150 reps`, Table B4) reveals an **honest break**:
- **Door B naive Wilks fails at weak signal:** `71.3%` on Gaussian weak, **`62.0%` even on the exact $\nu = 10$ scale-mixture world it was built for**, and `74.0%` on heavy-tail $\nu = 2$ — Neyman–Scott weak-signal breakdown still hits Door B, just as it hit `ours31` (`O18`).
- **Spherical Conditional LR (`SCLR+`):** conditioning on **both** $U_i(m_0)$ and $\|V_i(m_0)\|$, and resampling $\hat V_i^{(b)} = \|V_i(m_0)\| \cdot \mathrm{Unif}(S^{k_i-1})$ (using exact deterministic $2^N$ sign flips for small-$N$ 1D), achieves **exact finite-sample $\ge 95\%$ coverage under *every* per-reading noise-scale mixture** (`94.7%` Gaussian weak, `94.0%` $\nu=10$ weak, **`96.0%` Cauchy-like $\nu=2$** where Gaussian `GCLR+` collapses to `40.7%`!), is **bit-for-bit `Z`-invariant**, needs **zero $\pi$, zero $\eta$, and zero $\nu$**, and automatically widens on discordant data (`[5.02, 5.95]` on both 1D $N=10$ and 2D $N=2$ `5.0` vs `6.0`).
- **Two honest limitations of `SCLR+` alone (directly mapping to `O7`):**
  1. **1D small-$N$ floor:** in 1D ($k=1$), $S^0 = \{\pm 1\}$ has exact minimum p-value $1/2^N$, so for $N \le 4$ ($1/16 = 0.0625 > 0.05$) 1D `SCLR+` *always* returns `(0, inf)` (finite starting at $N = 5$, $1/32 = 0.03125$). In 2D/3D ($k \ge 2$), $S^{k-1}$ is continuous so $N = 2$ is already finite.
  2. **Scale-blindness on over-agreeing / noiseless data:** `SCLR+` throws away the absolute scale of $\Sigma_i$ (`O7` row *"$\Sigma$ up to per-reading unknown scale"*). On 1D $N=10$ over-agreeing readings (`5.491..5.509` at $\sigma=0.1$), `SCLR+` shrinks to `[5.497, 5.503]` (14× narrower than $\sigma=0.1$ allows); on a single perfectly aligned 2D reading `((5,0), (1,0))`, `SCLR+` collapses to `[5, 5]`.
- **Zero-tuning `O7` Floor-Row Union ($C_{\rm floor} = C_{95}(\text{Z-GCLR+}) \cup C_{95}(\text{SCLR+})$):** because `SCLR+` alone already has exact 95% coverage for every $v_i > 0$, taking the union with `Z-GCLR+` needs **no Bonferroni penalty** and gives a **guaranteed $\ge 95\%$ finite-sample realization of `O7`'s floor row** (*"actual noise is at least as large as declared $\Sigma_i$"*): `96.0%–98.0%` across all five Push 3 worlds (Table B4), **exact `Z`-invariance**, **respects the $\sigma$ floor** on over-agreeing 1D (`[5.459, 5.542]`) and aligned 2D (`[4.832, 5.179]`), and **automatically widens** on discordant 1D $N=10$ (`[5.023, 5.928]`) and 2D $N=2$ (`[5.023, 5.950]`). Its honest trade-offs: slightly conservative (~97–98%) on exact Gaussian $v_i \equiv 1$, and `(0, inf)` on 1D $N \le 4$ (inherited from `SCLR+`).

---

## Reproducible Tables

### Part A — Initial Killer Gallery (`python3 explore/e13_unnormalized_evidence.py`)

#### Table A1 — Deterministic scenes (zero $\pi$, zero $\eta$)

| Scene | Peak $\hat m$ | $K_{\min}$ (surprise $p \le$) | `T45` 95% set | `Moreira CLR` 95% set | Door B ($t_{10}$ Wilks) 95% set |
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

#### Table A2 — Common-ratio Monte Carlo 95% coverage (500 reps, seed 2026)

| World | `T45` cov (unb.) | `Moreira CLR` cov (unb.) | Raw Gaussian LR *(contrast)* |
|---|---|---|---|
| Pure noise (`3D, N=20, signal=0, m=20`) | `93.2%` (`93.6%`) | `96.2%` (`96.4%`) | `49.6%` |
| Weak signal (`3D, N=20, signal=1, m=2`) | `94.2%` (`83.2%`) | `96.2%` (`83.0%`) | `62.8%` |
| `T43` leverage (`1D, 1x r=20 + 99 nulls, m=2`) | `95.6%` (`0.0%`) | `95.6%` (`0.0%`) | `93.0%` |
| Strong signal (`3D, N=10, signal=5, m=2`) | `94.0%` (`0.0%`) | `94.0%` (`0.0%`) | `92.0%` |

---

### Part B — Push 1–3 Benchmarks (`python3 explore/e13_general_clr.py`, seed 20261121)

#### Table B1 — Push 1A: Deterministic scenes (`GCLR+` vs `GCLR` vs `Z-GCLR+` vs `Z-GCLR`)

| Scene | `GCLR+` ($m>0$) | `GCLR` ($\mathbb{RP}^1$) | `Z-GCLR+` ($m>0$, `Z`-safe) | `Z-GCLR` ($\mathbb{RP}^1$, `Z`-safe) |
|---|---|---|---|---|
| All-zero `0/0` (3D) | `(0, inf)` | `(0, inf)` | `(0, inf)` | `(0, inf)` |
| Scene A (`5.5, 5.6 ± 0.1`) | `[5.426, 5.68]` | `[5.426, 5.658]` | `[5.426, 5.68]` | `[5.426, 5.658]` |
| Scene B (`5.0, 6.0 ± 0.1`, 10σ) | `[5.358, 5.628]` | `[5.358, 5.628]` | `[5.358, 5.628]` | `[5.358, 5.628]` |
| Right-angle Scene C (`G29`) | `(0, inf)` | `(0, inf)` | `(0, inf)` | `(0, inf)` |
| Anti-parallel pair `(2,-1), (1,-2)` | `(0, 0.0001) U (361.4, inf)` | `EMPTY` | `(0, 0.0001) U (361.4, inf)` | `EMPTY` |
| Single anti-parallel `(2,-1)` at $m=0.5$ | rejects (`p=0.005`) | rejects (`p=0.005`) | pinhole (`p=1.0`) | pinhole (`p=1.0`) |

#### Table B2 — Push 1B: Full-type correlated + heterogeneous + subspace coverage (200 reps)
*$d=3$, subspace $k=2$, cross-channel $\rho=0.6$, heterogeneous noise ratios, $N=10$, true $m=2.5$.*

| World | `T45` | `GCLR+` | `GCLR(RP^1)` | `Z-GCLR+` | Raw Profile LR *(contrast)* |
|---|---|---|---|---|---|
| Pure noise (`r = 0`) | `93.0%` | `97.0%` | `96.0%` | **`95.0%`** | `44.0%` |
| Weak signal (`r = 1`) | `96.5%` | `95.0%` | `95.0%` | **`95.5%`** | `67.0%` |
| Strong signal (`r = 5`) | `95.5%` | `96.0%` | `96.0%` | **`95.0%`** | `95.0%` |

#### Table B3 — Push 2: Two-body composition coverage (`L5` lumping `M = 2 + 3 = 5`, 3D `N = 5` each, 150 reps)

| Regime | Raw infimal-conv | Exact Cond Lump LR | `Z`-Cond Lump LR |
|---|---|---|---|
| Weak (`signal = 1`) | `89.3%` | `97.3%` | **`96.0%`** |
| Moderate (`signal = 3`) | `90.7%` | `93.3%` | **`95.3%`** |
| Strong (`signal = 6`) | `94.0%` | `94.7%` | **`95.3%`** |

#### Table B4 — Push 3: Door B stress test vs. `GCLR+`, Spherical `SCLR+`, and `O7` Floor Union (`3D, N=15`, 150 reps)

| World | Door B ($t_{10}$ Wilks) | Gaussian `GCLR+` | Spherical `SCLR+(G)` | Spherical `SCLR+(t10)` | `O7` Floor Union |
|---|---|---|---|---|---|
| Gaussian weak (`r=1, m=2`) | `71.3%` | `95.3%` | `94.7%` | `96.7%` | **`98.0%`** |
| Gaussian strong (`r=5, m=2`) | `93.3%` | `94.0%` | `93.3%` | `96.0%` | **`97.3%`** |
| Scale-mix $\nu=10$ weak (`r=1, m=2`) | **`62.0%`** | `80.0%` | `94.0%` | `94.0%` | **`96.0%`** |
| Scale-mix $\nu=10$ strong (`r=5, m=2`) | `94.7%` | `93.3%` | `98.0%` | `92.7%` | **`98.0%`** |
| Heavy-tail $\nu=2$ (`r=3, m=2`) | `74.0%` | **`40.7%`** | `96.0%` | `95.3%` | **`96.0%`** |

#### Table B5 — Push 3B: Discordant vs. over-agreeing scenes (`GCLR+` vs `SCLR+` vs `O7` Floor Union)

| Scene | Gaussian `GCLR+` | Spherical `SCLR+` | Spherical `t_10` | `O7` Floor Union |
|---|---|---|---|---|
| 1D $N=10$ discordant (`5x 5.0, 5x 6.0`, $\sigma=0.1$) | `[5.459, 5.543]` | `[5.023, 5.831] U [5.95, 5.95]` | `[5.023, 5.95]` | **`[5.023, 5.928]`** |
| 1D $N=10$ over-agreeing (`5.491..5.509`, $\sigma=0.1$) | `[5.459, 5.542]` | `[5.497, 5.503]` *(breaks floor)* | `[5.497, 5.503]` | **`[5.459, 5.542]`** |
| 2D $N=2$ discordant (`(5,0) & (6,0)`, $\sigma=0.1$) | `[5.365, 5.639]` | `[5.023, 5.95]` | `[5.023, 5.95]` | **`[5.023, 5.95]`** |
| 2D $N=1$ aligned (`(5,0), (1,0)`, $\sigma=0.1$) | `[4.832, 5.179]` | `[5, 5]` *(collapses)* | `[5, 5]` | **`[4.832, 5.179]`** |
| 1D $N=4$ aligned (`m=2.0`, $\sigma=0.01$) | `[1.979, 2.019]` | `(0, inf)` *($1/2^4 > 0.05$)* | `(0, inf)` | `(0, inf)` |

---

## Summary matrix: how the three `O7` rows map cleanly onto unnormalized evidence + exact sets

| `O7` declaration of $\Sigma_i$ | Compatibility curve | Exact 95% set (zero $\pi$, zero $\eta$) | Passes exact `Z`? | Automatic widening on discordance? | Honest limitation |
|---|---|---|---|---|---|
| **Row 1: $\Sigma_i$ exact Gaussian** | `T36` $K(m)$ + $K_{\min}$ surprise | `Z-GCLR` / `Z-GCLR+` (or unweighted `GCLR`) | `Z-GCLR`: **yes** (`GCLR`: no) | No (Kac–Bernstein; surprise lives in $K_{\min}$) | `Z-GCLR` has rank-1 pinhole ghost; $w(H)$ is a family |
| **Row 2: $\Sigma_i$ up to per-reading scale** | `T36` / `T48` shape | `SCLR+` (spherical conditional LR) | **yes** (bit-for-bit) | **yes** (exact $\ge 95\%$ for *all* $v_i > 0$) | Shrinks past $\sigma$ on over-agreeing data; 1D $N \le 4$ is `(0, inf)` |
| **Row 3: $\Sigma_i$ is a floor ($v_i \ge 1$)** | `T36` $K(m)$ + $K_{\min}$ | `C_floor = Z-GCLR+ ∪ SCLR+` | **yes** (bit-for-bit) | **yes** (and respects $\sigma$ floor!) | Slightly conservative (~97–98%) at $v_i \equiv 1$; 1D $N \le 4$ is `(0, inf)` |

## What this exploration does NOT solve (open for harvest / Can's decision)

1. **The output-type & attitude fork (Can / Claude):** changing the output from a normalized area-1 probability law + point (`L2`, `G28`) to an unnormalized residual curve $K(m)$ + exact compatibility set (allowed to return `(0, inf)`) plus a separate leftover surprise $K_{\min}$ is a ruling about the goal, not a calculation — and choosing among `O7`'s three rows (exact Gaussian vs. up-to-scale vs. floor) is an instrument declaration.
2. **Weight-function non-uniqueness (`E2`, `E3`):** `Z-GCLR`'s $w(H)$ is a whole family of valid exact `Z`-safe weights; whether anything narrows the family (e.g., minimax power or leverage optimality) is open.
3. **1D small-$N$ floor for exact scale-free tests ($N \le 4$):** any exact distribution-free sign-flip test in 1D cannot reject at $\alpha = 0.05$ with $N \le 4$ ($1/2^4 = 0.0625 > 0.05$); finite 1D widening at $N = 2$ requires either a parametric tail assumption (Door B with calibrated weak-signal critical values) or $k \ge 2$ spatial dimensions.
4. **Random pure noise beside an informative block:** `explore/e2-null-combination`'s obstruction still applies to random noisy zero-signal blocks (distinct from exact observed `0/0`), though `Z-GCLR`'s $H/(1+H)$ weights substantially dampen them compared to unweighted `GCLR`.
