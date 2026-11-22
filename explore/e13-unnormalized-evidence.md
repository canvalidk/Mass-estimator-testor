# E13: Unnormalized Evidence ("Bell with σ = ∞") — Thirteen Insights

**Branch:** `explore/e13-unnormalized-evidence` (from `main` `d2a3352`, 21–22 Nov 2026).  
**Status:** exploration input for Claude to harvest (`explore/README.md`, R-2026-11-07-2). Nothing in `definition/`, `problems/`, `rulings.md`, or `met/` is modified.  
**Code & checks (pure Python stdlib — zero external dependencies):**
- Initial prototype & killer gallery (Insights 1–5): `explore/e13_unnormalized_evidence.py`, `tests/test_e13_unnormalized_evidence.py`
- Push 1–3 general-type CLR, exact two-body composition, and Door B / spherical CLR / O7 floor union (Insights 6–8): `explore/e13_general_clr.py`, `tests/test_e13_general_clr.py`
- Push 4–5 Z-safe weight Pareto frontier & closing the 1D $N=2$ Door B gap (Insights 9–10): `explore/e13_push4_5.py`, `tests/test_e13_push4_5.py`
- Push 6 $GL(d)$-covariant metric mass $g \in \mathrm{Sym}^+_d$ (Insights 11–13): `explore/e13_metric_mass.py`, `tests/test_e13_metric_mass.py`
- Run all 24 checks: `python3 -m unittest -v tests/test_e13_unnormalized_evidence.py tests/test_e13_general_clr.py tests/test_e13_push4_5.py tests/test_e13_metric_mass.py` (also pytest-discoverable)

## Standard mathematics first (`GOAL.md` standing rule 6)

- **Kac–Bernstein theorem** (Kac 1939; Bernstein 1941): for independent $X_1, X_2$, the sample mean $\bar X$ is independent of the difference $X_1 - X_2$ *if and only if* $X_1, X_2$ are Gaussian. (Equivalently, orthogonal projections $U_i, V_i$ of an isotropic error vector are independent iff the noise is Gaussian — directly governing why $U$-conditioning works for Gaussian `GCLR` and breaks for Student-$t$.)
- **Unnormalized / profile likelihood & infimal convolution** (Royall 1997; Barndorff-Nielsen; epigraphical addition): combining independent evidence adds negative log-likelihoods $K_1 + K_2$; combining two parameters under $M = m_1 + m_2$ (or $G = g_1 + g_2$ on the SPD cone $\mathrm{Sym}^+_d$) profiles out the split via infimal convolution, with Chernoff (1954) one-sided boundary mixture $\frac{1}{2}\delta_0 + \frac{1}{2}\chi^2_1$ when one body is unmeasured (`0/0`) under positivity.
- **Weak-identification exact conditional test inversion** (Anderson & Rubin 1949; Moreira 2003, *Econometrica* 71:1027–1048; Andrews–Moreira–Stock 2006; Gleser & Hwang 1987; Dufour 1997; Dwass 1957 / Barnard 1963 exact Monte Carlo tests): although no single nuisance-free *conditional likelihood* exists because $U_i(m_0)$ depends on $m_0$ (Kalbfleisch & Sprott 1970, cited in `problems/O15.md`), **pointwise conditional test inversion** given $(U_1, \dots, U_N)$ is exact in finite samples on both the full scalar type and the metric cone $\mathrm{Sym}^+_d$ (extending `T50` / Markovsky & Van Huffel TLS) without a prior $\pi$ or warp $\eta$, provided unbounded sets are allowed.
- **Stein's identity & Cauchy–Schwarz** (Stein 1981): for $H = \|\xi\|^2 \sim \chi^2_k$, $E[k w(H) + 2 H w'(H)] = E[H w(H)]$, pinning down the exact diffuse weak-signal asymptotic relative efficiency $\text{ARE}_k(w) = (E[H w(H)])^2 / (E[H] E[H w(H)^2]) \le 1$ of every $U$-weighted exact test.
- **Spherical / sign-flip randomization tests** (Fisher 1935; Lehmann & Romano 2005, §15.2): conditioning on both $U_i$ and $\|V_i\|$ makes $V_i / \|V_i\| \sim \mathrm{Unif}(S^{k_i-1})$ exact under *any* per-reading noise-scale mixture $v_i \Sigma_i$ (Gaussian, Student-$t$, Cauchy, adversarial scale spikes).
- **Metric / covector duality on $V \oplus V^*$** (bilinear forms on $\mathrm{Sym}^+_d$; BFGS secant positivity; `O17` / `T50`): modeling mass as $g: V \to V^*$ ($\mathbf{F} = g\mathbf{a}$, $g = g^\top \succ 0$) makes the covector residual $K_i(g) = (\tilde{\mathbf{F}}_i - g\tilde{\mathbf{a}}_i)^\top V_i(g)^{-1}(\tilde{\mathbf{F}}_i - g\tilde{\mathbf{a}}_i)$ strictly invariant under all invertible coordinate changes $A \in GL(d)$ (no Euclidean spatial metric $\delta$) and under covector duality $g \leftrightarrow g^{-1}$.

---

## The 13 Insights (countable summary for harvest)

### Insight 1 — `0/0` forces the identity of combination ($L \equiv 1$, "bell with $\sigma = \infty$"), *not* a normalized prior $\pi_c$ (`E5`, `O6`, `O11`, `O2`, `G28`, `B16`)

In `docs/vd-docs-context/why_mass_cannot_be_undefined_2026-09-29.md` §2.3, the density-route contradiction proves that a stationary `0/0` trace must act as the **identity element of combination** (combining `0/0` with `5 ± 0.5 kg` must return `5 ± 0.5 kg` untouched). The leap to §2.4 demanded a **proper area-1 probability law** $\pi_c$ and a finite point $c = \sigma_F/\sigma_a$ even when no data exist.

- An unnormalized Gaussian bell $e^{-(x-c)^2/(2\sigma^2)}$ as $\sigma \to \infty$ is literally $L(m) \equiv 1$ ($K(m) \equiv 0$) — its center $c$ vanishes and its inverse-variance combination weight is $1/\infty^2 = 0$ (matching R-2026-09-29-3).
- Even with noiseless $\sigma_{\rm in} = 0$ excluded from the noisy domain (`type.md`), finite $\sigma_{\rm in} > 0$ **forces** infinite-tail likelihoods (`B16`) and unbounded sets `(0, ∞)` on totally ordinary nonzero measurements (e.g. `F = 10 ± 1, a = 1 ± 1` or perpendicular 3D singletons), because dividing by a noisy acceleration can never rule out arbitrarily small true $a^*$.
- **Consequence:** Dropping area-1 normalization (returning the unnormalized residual curve $K(m)$ + a horizontal compatibility/confidence set allowed to return $(0,\infty)$) dissolves `O6` (shape of $\pi$), `O11` (multi-instrument center), `G6`/`O2` (forced point on `0/0`), and the Gleser–Hwang / null-floor impossibility (`O15`, `G28`), while keeping `0/0` total and never `undefined`. Any caller who wants a Bayesian posterior multiplies $e^{-K/2}$ by their own declared prior $\pi$ on top.

### Insight 2 — Exact Gaussian `L13` + believed $\sigma$ (`L18`) and discordance widening (`R-2026-11-17-3` / `T30`) contradict each other by Kac–Bernstein (`O14`, `G13`, `G18`, `G19`, `G22`)

Under independent exact Gaussian errors (`L13`), the mean $\bar X$ and gap $X_1 - X_2$ are statistically independent (Kac–Bernstein, verified in `tests/test_believed_less.py`). Therefore, if exact Gaussianity and exact $\sigma$ are truly believed (`R-2026-11-17-1`), Scene B (`5.0, 6.0 ± 0.1`) genuinely has the same conditional precision for $m$ as Scene A (`5.5, 5.6 ± 0.1`): `5.5 ± 0.14`. Human intuition expects widening because real errors have heavy tails (occasional flukes) or off-graph departures — which `explore/e12-scalar-mass-proof` (`e12-adding-gaussian-noise.md`) independently confirmed via Gaussian noise-degradation towers.

This splits `O14` cleanly into **two contradiction-free doors** (no ad-hoc $g$ or $\eta$):
- **Door A (Strict consequence of exact Gaussian inputs):** keep exact `T36` residual $K(m)$ and exact compatibility set (`T45` / `GCLR` / `Z-GCLR`), drop automatic widening inside the mass law, and report the leftover $K_{\min} = \min_{m>0} K(m)$ (with exact upper-bound surprise $p$-value $P(\chi^2_M \ge K_{\min})$) beside the mass.
- **Door B (Heavy-tailed / scale-uncertain error model):** drop exact fixed-scale Gaussian `L13` as universal — either declare a heavy-tailed law (`T48` Student-$t$, calibrated rigorously via Insight 10) or condition out per-reading noise scales via spherical `SCLR+` / the `O7` floor union (Insight 8).

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
- does **not** crash at 50σ (where `ours33` raises `ValueError`), instead splitting cleanly into `[4.968, 5.041] ∪ [5.963, 6.035]` (*"one reading is a fluke: mass is near 5 or near 6"*). *(See Insights 8 & 10 for its weak-signal break and exact calibration.)*

### Insight 6 (Push 1) — Full-type General Conditional LR (`GCLR` / `GCLR+`) and `Z`-safe `U`-weighted `Z-GCLR` extend `B24` past the common-ratio slice (`B24`, `T45`, `E2`, `E3`, `L20`, `T43`)

`problems/O15.md` (`B24`, 10 Nov) noted that heterogeneous covariances $\Sigma_i$ remove the common two-column rotation so *"the CLR formula below is not thereby extended."* Conditioning pointwise at candidate $\phi_0 = \arctan(m_0)$ on the whitened along-line projections $U_i = Q_{\parallel,i}(\phi_0)^\top z_i \in \mathbb{R}^{k_i}$ (inside each reading's $2k_i$-dimensional whitened span $W_i$) makes the orthogonal residuals $V_i = Q_{\perp,i}(\phi_0)^\top z_i \sim N(0, I_{k_i})$ **exact, independent across readings, independent of $U_i$, and free of all unknown pushes $\alpha_i^*$** on the **entire full type** (arbitrary correlated $\Sigma_i$, heterogeneous ratios, subspaces $N_i$). Drawing $V_i^{(b)} \sim N(0, I_{k_i})$ and evaluating the profile likelihood ratio on a fixed grid yields an exact finite-sample Dwass–Barnard conditional test (`GCLR+` on $m>0$, `GCLR` on $\mathbb{RP}^1$):
- **Full-type exact coverage (Table B2):** on correlated $d=3, k=2, \rho=0.6$ mixed-ratio worlds (`N=10`), raw profile LR collapses to `44.0%` (pure noise) and `67.0%` (weak signal), while `GCLR+`, `GCLR(RP^1)`, and `Z-GCLR+` hold exact ~95% (`95.0%–97.0%`) with zero $\pi$ and zero $\eta$.
- **Exact-`Z`-safe `U`-weighted `Z-GCLR` (`w_i = w(H_i(m_0))`):** because $H_i(m_0) = \|U_i\|^2$ is fixed by conditioning, any invariant weight $w(H_i)$ with $w(0)=0$ preserves **exact finite-sample 95% coverage** (`95.0%, 95.5%, 95.0%`), passes **units, swap, rotation, exact-zero `Z` (bit-for-bit!), and $\sigma \to \infty$**, eliminates the multi-direction anti-parallel hilltop ghost `(2,-1), (1,-2)` (`EMPTY` on $\mathbb{RP}^1$, boundary on $m>0$), returns `(0, inf)` on right-angle Scene C (`G29`), and sharpens `T43` leverage width.
- **Rank-1 pinhole boundary:** under $U$-conditioning alone, exact `Z`-invariance forces ignoring any reading with $U_i(m_0) = 0$. Therefore on a *rank-1* anti-parallel dataset like a single reading `(2, -1)` at its orthogonal mass $m_0 = 0.5$ (where *every* reading has $U_i(0.5) = 0$ despite $K = 500$), `Z-GCLR` must give $p = 1.0$ (a pinhole ghost), whereas unweighted `GCLR` rejects even that pinhole ($p = 0.005$) at the price of failing both `Z` and $\sigma \to \infty$.

### Insight 7 (Push 2) — Exact weak-signal two-body composition (`L5` lumping and `L14` common force) lifts coverage to 95% and hits the exact Chernoff boundary (`L5`, `L14`, `B26`)

Raw infimal convolution (`Insight 4`) inherits raw profile LR's Neyman–Scott gap at weak signal (`89.3%` on 3D `M = 2 + 3 = 5`, $N=5$ each, signal=1). Conditioning on $(U_1(m_1), U_2(M_0 - m_1))$ calibrates the profiled lump LR $\mathrm{LR}_{\rm lump}(M_0) = \inf_{m_1+m_2=M_0}[K_1(m_1)+K_2(m_2)] - K_{1,\min} - K_{2,\min}$ exactly in finite samples (taking $\sup_{0 < m_1 < M_0}$ for a guaranteed $\ge 95\%$ confidence set for $M_0$):
- **Weak-signal calibration fixed (Table B3):** lifts weak-signal 3D lumping coverage from `89.3%` (`Raw infimal-conv`) to `97.3%` (unweighted) and **`96.0%` (`Z-Cond Lump LR`)**, with `95.3%` at moderate and strong signal.
- **1-df sharpness preserved at strong signal:** on strong `2 + 3 = 5` (combined $\sigma_M = 0.10$), $\sup_{m_1} p(m_1, M_0 - m_1)$ recovers the exact 1-df interval `[4.80, 5.20]` with zero 2-df Bonferroni widening.
- **Exact `L14` common-force duality:** swapping $(F, a)$ turns `L14` ($1/M = 1/m_1 + 1/m_2$) literally into `L5` on $1/M$ (bit-identical p-values).
- **Chernoff one-sided boundary on `(5.0 ± 0.5 kg) + (0/0)`:** with `Z`-safe weights (`z_smooth`), the `0/0` body drops out bit-for-bit and the test automatically recovers the exact **Chernoff $\frac{1}{2}\delta_0 + \frac{1}{2}\chi^2_1$ one-sided 95% bound** `[4.18, inf)` ($5.0 - 1.645 \times 0.5$: $p(4.00)=0.020$, $p(4.18)=0.051$, $p(4.35)=0.100$, $p(5.0)=p(50.0)=1.000$).
- **Honest trade-off:** exact weak-signal composition needs the readings (or their $2k$-dim within-span summaries), not merely the 1D curves $K_1(m), K_2(m)$.

### Insight 8 (Push 3) — Door B naive Wilks breaks at weak signal; Spherical Conditional LR (`SCLR+`) and the `O7` floor union deliver exact scale-free coverage (`O14`, `O7`, `T37`, `T48`, `O18`)

Stress-testing Door B (`Insight 5`, unnormalized `T48` Student-$t_{10}$ with Wilks `3.841` cutoff) on Monte Carlo worlds (`3D, N=15, 150 reps`, Table B4) reveals an **honest break**:
- **Door B naive Wilks fails at weak signal:** `71.3%` on Gaussian weak, **`62.0%` even on the exact $\nu = 10$ scale-mixture world it was built for**, and `74.0%` on heavy-tail $\nu = 2$ — Neyman–Scott weak-signal breakdown still hits Door B, just as it hit `ours31` (`O18`).
- **Spherical Conditional LR (`SCLR+`):** conditioning on **both** $U_i(m_0)$ and $\|V_i(m_0)\|$, and resampling $\hat V_i^{(b)} = \|V_i(m_0)\| \cdot \mathrm{Unif}(S^{k_i-1})$ (using exact deterministic $2^N$ sign flips for small-$N$ 1D), achieves **exact finite-sample $\ge 95\%$ coverage under *every* per-reading noise-scale mixture** (`94.7%` Gaussian weak, `94.0%` $\nu=10$ weak, **`96.0%` Cauchy-like $\nu=2$** where Gaussian `GCLR+` collapses to `40.7%`!), is **bit-for-bit `Z`-invariant**, needs **zero $\pi$, zero $\eta$, and zero $\nu$**, and automatically widens on discordant data (`[5.02, 5.95]` on both 1D $N=10$ and 2D $N=2$ `5.0` vs `6.0`).
- **Two honest limitations of `SCLR+` alone (directly mapping to `O7`):**
  1. **1D small-$N$ floor:** in 1D ($k=1$), $S^0 = \{\pm 1\}$ has exact minimum p-value $1/2^N$, so for $N \le 4$ ($1/16 = 0.0625 > 0.05$) 1D `SCLR+` *always* returns `(0, inf)` (finite starting at $N = 5$, $1/32 = 0.03125$). In 2D/3D ($k \ge 2$), $S^{k-1}$ is continuous so $N = 2$ is already finite.
  2. **Scale-blindness on over-agreeing / noiseless data:** `SCLR+` throws away the absolute scale of $\Sigma_i$ (`O7` row *"$\Sigma$ up to per-reading unknown scale"*). On 1D $N=10$ over-agreeing readings (`5.491..5.509` at $\sigma=0.1$), `SCLR+` shrinks to `[5.497, 5.503]` (14× narrower than $\sigma=0.1$ allows); on a single perfectly aligned 2D reading `((5,0), (1,0))`, `SCLR+` collapses to `[5, 5]`.
- **Zero-tuning `O7` Floor-Row Union ($C_{\rm floor} = C_{95}(\text{Z-GCLR+}) \cup C_{95}(\text{SCLR+})$):** because `SCLR+` alone already has exact 95% coverage for every $v_i > 0$, taking the union with `Z-GCLR+` needs **no Bonferroni penalty** and gives a **guaranteed $\ge 95\%$ finite-sample realization of `O7`'s floor row**: `96.0%–98.0%` across all five Push 3 worlds (Table B4), **exact `Z`-invariance**, **respects the $\sigma$ floor** on over-agreeing 1D (`[5.459, 5.542]`) and aligned 2D (`[4.832, 5.179]`), and **automatically widens** on discordant 1D $N=10$ (`[5.023, 5.928]`) and 2D $N=2$ (`[5.023, 5.950]`).

### Insight 9 (Push 4) — The `Z`-safe per-reading weight class $w(H_i)$ is a genuine three-regime Pareto family governed by Stein's lemma + Cauchy–Schwarz (`E2`, `E3`, `T43`, `T45`, `L10`, `L19`)

Answering `E2` and `E3`'s central question (*"a per-reading statistic that is `Z`-safe, order-free, borrows nothing, declares no size law, and makes a true statement on noise: one member or a family, and what selects a member?"*):

For **any** nonnegative per-reading weight $w_i = w(H_i(m_0))$ with $w(0) = 0$ (depending only on reading $i$'s own along-line energy $H_i = \|U_i\|^2$), both the weighted conditional score `W-T45` ($T_w = \sum w_i S_i / \sqrt{\sum w_i^2 C_i} \sim N(0, 1)$) and weighted `Z-GCLR` have **exact finite-sample 95% coverage**, pass units/swap/rotation/order/exact-`Z` bit-for-bit, borrow nothing (`L19`), declare no size law (`L10`), and preserve the `L11` rate $\Delta_N/\sqrt{N} \to \infty$.

By **Stein's identity** ($E[k w(H) + 2 H w'(H)] = E[H w(H)]$ for $H \sim \chi^2_k$) and **Cauchy–Schwarz**, every member $w(H)$ obeys two exact closed-form laws (Table C1):
$$
\text{ARE}_k(w) = \frac{\bigl(E[H w(H)]\bigr)^2}{E[H]\, E[H w(H)^2]} \le 1
\quad (\text{equality iff } w \equiv \text{const}),
\qquad
N_{\rm drown}(w; r) = \frac{w(r^2)^2 r^2}{E[H w(H)^2]}.
$$
**No single weight dominates all three signal geometries (Table C2):**
- **`w = 1` (`T45`):** wins **diffuse weak** (`ARE = 100%`, uniquely optimal) and **mixed strong** (`55.0%` power on `1x r=20 + 25x r=5`), **but drowns** when 1 real reading ($r=20$) is buried in 2,000 random noisy zeros (`43.0%` power, $N_{\rm drown} = 400$; matching `E2`'s null-erasure result).
- **`w = H` (parameter-free linear weight):** survives **4.27 million** random noisy zeros in 1D ($6.1 \times 10^5$ in 3D — **10,667× better than `T45`**!), keeps **`99.0%` power** under 2,000 noisy zeros (and restores ghost-free `Z-GCLR+` on `1x r=20 + 500 nulls` from width `1.65` right back to the oracle single-reading width **`0.90`** while still returning `EMPTY` on anti-parallel data!), and retains `60.0%` (1D) / `71.4%` (3D) diffuse weak efficiency — **but** pays a clear tax on **mixed strong** (`28.0%` vs `55.0%`) because unbounded $w=H$ overweights $r=20$ vs $r=5$ by $16\times$.
- **`w = H^2 / (100 + H^2)` (bounded power weight):** wins both **mixed strong** (`52.7%`) and **buried sparse** (`98.0%` power under 2,000 nulls; `Z-GCLR+` width `0.90` on 500 nulls), **but** cuts **diffuse weak** efficiency to `40.0%` (1D) and introduces a tunable scale $c=10$.
- **Verdict for `E2` / `E3`:** the class is an irreducible **Pareto family**, and selecting a member is a power preference over signal geometry (diffuse weak vs. mixed strong vs. sparse buried), not a property of Newtonian mass. *(Consistent with `E2`: for any fixed weight with $E[H w(H)^2] > 0$, sufficiently large $N \gg N_{\rm drown}$ eventually still drowns a fixed reading.)*

### Insight 10 (Push 5) — Why naive $U$-conditioning breaks under Student-$t_\nu$, and how Piecewise Exact Door B (`AR-bounded t-LR` + `SCLR+`) closes the 1D $N=2$ gap (`O14`, `O7`, `O18`, `T48`)

1. **Exact Kac–Bernstein barrier for Student-$t_\nu$ conditioning:** Under `T48`'s scale-mixture noise ($v_i \sim \mathrm{IG}(\nu/2, \nu/2)$), $U_i$ and $V_i$ are **coupled** by the shared scale $v_i$ via $v_i \mid (U_i, \alpha_i^*) \sim \mathrm{IG}((\nu+k)/2, (\nu + \|U_i - \alpha_i^*\|^2)/2)$, which depends on the **unknown** push $\alpha_i^*$. Consequently, naive $U$-conditioning (`indep-t-CLR`, drawing $V_i^{(b)} \sim t_\nu$ independent of $U_i$) **breaks even when $\nu$ is known exactly**: coverage drops to **`85.3%`** on 3D $N=15, \nu=3$ pure noise (`90.7%` weak signal) and **`89.0%`** on 1D $N=2, \nu=1$.
2. **Exact push-free `t-AR` and 1-line proved `AR-bounded t-LR`:** Unconditionally (not conditioning on $U_i$), $V_i(m^*) \sim t_\nu(0, I_{k_i})$ is **100% free of all unknown pushes $\alpha_i^*$** on the full correlated/subspace type, with exact representation $(1 + K_i(m^*)/\nu)^{-1} \sim \mathrm{Beta}(\nu/2, k_i/2)$.
   - **`t-AR`** ($T_{t_\nu}(m_0) \le q_{0.95}^{\text{t-AR}}$, where $q_{0.95}^{\text{t-AR}} = 6.939$ for 1D $N=2, \nu=10$) has exact 95.0% coverage everywhere, **but** returns **`EMPTY`** on 10σ and 50σ Scene B.
   - **`AR-bounded t-LR`** ($T_{t_\nu}(m_0) - T_{t_\nu,\min} \le q_{0.95}^{\text{t-AR}}$) has a **1-line mathematical proof of $\ge 95\%$ coverage** for every $(\Sigma_i, N_i, \alpha_i^*, m^*)$ under $t_\nu$ (since $T(m^*) - T_{\min} \le T(m^*)$), **never returns `EMPTY`**, and closes the 1D $N=2$ gap cleanly: **`[5.353, 5.734]`** on Scene A, **`[4.854, 6.154]`** on 10σ Scene B, **`[4.963, 5.067] ∪ [5.937, 6.041]`** on 50σ Scene B, and **`(0, inf)`** on `0/0` (Table C4).
3. **Piecewise Exact Door B (`door_b_piecewise_exact_pvalue`):** Switching **only on declared experiment dimensions** (never on data $y_i$) — using `AR-bounded t-LR` on 1D $N \le 4$ and spherical `SCLR+(t_\nu)` everywhere else ($N \ge 5$ or $k \ge 2$) — is **rigorously proved $\ge 95\%$ on every experiment under declared $t_\nu$ with zero Bonferroni penalty**, closing both 1D $N=2$ (`98.2%–99.2%`) and 3D $N=15$ weak signal (`94.0%`, vs naive Door B's `62.0%`).
4. **Honest limitations of the 1D $N \le 4$ branch (Table C5):** Unlike `SCLR+` (which is scale-free for $N \ge 5$ or $k \ge 2$), 1D $N \le 4$ `AR-bounded t-LR` genuinely relies on the declared $\nu$ (if true noise is $\nu=2$ while $\nu=10$ is declared, 1D $N=2$ coverage drops to `85.5%–89.0%`) and fails exact-zero `Z` (since $q_{0.95}^{\text{t-AR}}(N)$ grows with $N$: `4.442, 6.939, 8.971, 10.987`).

### Insight 11 (Push 6A–B) — Switching to metric mass $g \in \mathrm{Sym}^+_d$ ($\mathbf{F} = g\mathbf{a}$) makes "$\sigma = \infty$" unavoidable even without noise, is strictly $GL(d)$-covariant, and splits surprise into a 3-rung ladder (`O17`, `E7`, `L22`, `T36`, `T50`)

Reading acceleration as a vector $\mathbf{a} \in V$, force as a covector $\mathbf{F} \in V^*$, and mass as a symmetric positive-definite metric $g: V \to V^*$ ($F_i = g_{ij} a^j$, inverse mass $g^{-1}: V^* \to V$):
- **Full $GL(d)$ and swap invariance (no Euclidean spatial metric $\delta$):** under any invertible basis change $A \in GL(d)$ (unequal axis units, non-orthogonal shears) $\mathbf{a} \mapsto A\mathbf{a}$, $\mathbf{F} \mapsto A^{-\top}\mathbf{F}$, $g \mapsto A^{-\top} g A^{-1}$, the covector residual $K_i(g) = (\tilde{\mathbf{F}}_i - g\tilde{\mathbf{a}}_i)^\top V_i(g)^{-1}(\tilde{\mathbf{F}}_i - g\tilde{\mathbf{a}}_i)$ and all conditional tests built on it are **100% invariant** (verified in 2D and 3D to $10^{-10}$). Force–acceleration swap ($V \leftrightarrow V^*$) maps $g \leftrightarrow g^{-1}$ identically.
- **Unnormalized "$\sigma = \infty$" is forced on partial-axis pushes:** pushing a 2D or 3D body only along $e_1$ leaves $K_x(g)$ identically `0.000000` across every $g_{yy}, g_{zz} \in (0, \infty)$ — literally $\sigma = \infty$ along unpushed axes. Adding a later $e_2$ push simply sums $K_x(g) + K_y(g)$ and pins down `diag(2, 5)` with zero prior $\pi$.
- **Three-rung nested surprise ladder (`scalar ⊂ metric ⊂ linear`, Table D1):** $K_{\min}^{\rm scalar} \ge K_{\min}^{\rm metric} \ge K_{\min}^{\rm linear} \ge 0$ separates **anisotropy surprise** (`O17`: `d_iso = 59.9, d_rec = 0.0` on `diag(2, 5)`), **reciprocity / gyroscopic surprise** (`E7`: `d_rec = 66.2, K_lin = 0.0` on `[[2, -1], [1, 2]]`), and irreducible residual.

### Insight 12 (Push 6B–C & 6E) — Exact conditional tests on $\mathrm{Sym}^+_d$ (`Metric-GCLR`, `Metric-Z-GCLR`, `Metric-SCLR+`, and `L22` isotropy) kill `T50`'s full-rank ghost and eliminate the 1D small-$N$ floor (`T50`, `O17`, `L22`, `O7`)

- **Full-rank (`rank = 3`) 2D anti-parallel hilltop ghost in `T50`:** on 6 symmetric anti-parallel readings at $g_0 = I_2$ ($K(I_2) = 3600$), `T50`'s conditional vector score accepts with `T2 = 0.0000, rank = 3, p = 1.0000`; both `Metric-GCLR(Sym_2)` and bit-for-bit `Z`-safe `Metric-Z-GCLR(Sym_2)` kill it (`p = 0.01`).
- **Exact finite-sample coverage on 2D anisotropic `g = diag(2, 5)` (Table D2):** across pure noise ($r = 0$), partial $x$-axis-only excitation ($r = 4$), weak signal ($r = 1.5$), and strong signal ($r = 5$), `Metric-Z-GCLR+` and `Metric-SCLR+` hold ~95% (where raw 3-df Wilks LR drops to `65%`).
- **No small-$N$ discrete spherical floor in $d \ge 2$:** because $S^{d-1}$ ($S^1$ in 2D, $S^2$ in 3D) is continuous, `Metric-SCLR+` and the `O7` Metric Floor Union have **no $1/2^N$ floor** and hold **`93.0% / 94.0%`** even under Cauchy-like $\nu = 2$ heavy tails where `T50` (`41.0%`) and Gaussian `GCLR+` (`42.0%`) collapse!
- **Exact finite-sample test of scalar isotropy (`L22`):** $\sup_{m > 0} p_{\rm metric}(m I_d)$ accepts isotropic `2 I_2` (`p = 0.98 / 0.96`) and $x$-only pushes (`p = 1.00`, honestly admitting scalar compatibility!), while rejecting anisotropic `diag(2, 5)` (`p = 0.01 / 0.01`).

### Insight 13 (Push 6D) — Switching to the metric cone $\mathrm{Sym}^+_d$ resolves `O17`'s `L5`/`L14` conflict simultaneously (`O17`, `L5`, `L14`, `B26`)

- `O17` showed that no single scalarization of an anisotropic metric satisfies both `L5` (lumping) and `L14` (common force). Keeping the full unnormalized metric curve $K(g)$ on $\mathrm{Sym}^+_d$ makes **both** `L5` lumping ($G = g_1 + g_2$, infimal convolution over the operator interval $0 \prec g_1 \prec G$) and `L14` common force ($G^{-1} = g_1^{-1} + g_2^{-1}$) hold **simultaneously and exactly** with exact swap duality $g \leftrightarrow g^{-1}$ (`diag(2,5) + diag(3,1) -> diag(5,6)` and reduced `diag(6/5, 5/6)` both with $K = 0.000000$; Table D3).
- Gluing an anisotropic body `diag(2, 5)` onto a body pushed *only* along $x$ ($g_{2,xx} = 3$, $g_{2,yy}$ unmeasured = $\sigma_\infty$) rules out $G_{yy} < 5$ ($K = 20.79$ at $G_{yy} = 4.5$) while leaving $G_{yy} > 5$ identically flat ($K = 0.000$ at $5.05, 6, 20$).
- **Honest limitations of metric mass:** (1) requires a shared known frame across readings (not arbitrary per-reading spatial rotations); (2) any single reading with $\langle \mathbf{F}, \mathbf{a} \rangle > 0$ fits some $g \succ 0$ with $K_{\min}^{\rm metric} = 0$ (leftover surprise begins at $N \ge 2$ via reciprocity failure or $N > d$); (3) symmetry $g = g^\top$ imports reciprocity (parked by `R-2026-10-27-2` as a later-version commitment).

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

### Part C — Push 4 & Push 5 Benchmarks (`python3 explore/e13_push4_5.py`)

#### Table C1 — Push 4A: Exact Stein + Cauchy–Schwarz Pareto Metrics (12-digit quadrature, $r_{\rm strong} = 20$)

| Weight $w(H)$ | `ARE_1` (1D weak) | `N_drown` (1D, $r=20$) | `ARE_3` (3D weak) | `N_drown` (3D, $r=20$) |
|---|---:|---:|---:|---:|
| `T45` ($w \equiv 1$) | **`100.00%`** | `400.0` | **`100.00%`** | `133.3` |
| `z_smooth` ($H/(1+H)$) | `92.05%` | `852.2` | `97.59%` | `212.0` |
| `z_quad` ($H^2/(1+H^2)$) | `87.58%` | `654.2` | `97.26%` | `162.9` |
| `z_cap10` ($H^2/(100+H^2)$) | `39.95%` | `14,997.4` | `59.07%` | `1,868.7` |
| `z_linear` ($w = H$, **0 params**) | `60.00%` | **`4,266,666.7`** | `71.43%` | **`609,523.8`** |
| `z_pow2` ($w = H^2$) | `23.81%` | `1.08 × 10¹⁰` | `35.35%` | `9.85 × 10⁸` |

#### Table C2 — Push 4B: Three-Regime Pareto Check on 1D `W-T45` (all weights hold ~95% coverage)

| Weight $w(H)$ | Diffuse weak (`N=400, r=0.6`, rej `m=3`) | Mixed strong (`1x r=20 + 25x r=5`, rej `m=2.35`) | Buried 500 nulls (`1x r=20 + 500`, rej `m=4`) | Buried 2,000 nulls (`1x r=20 + 2000`, rej `m=4`) |
|---|---:|---:|---:|---:|
| `T45` ($w \equiv 1$) | **`13.3%`** | **`55.0%`** | `80.0%` | `43.0%` *(drowned)* |
| `z_smooth` ($H/(1+H)$) | **`13.7%`** | **`55.0%`** | `90.0%` | `64.5%` |
| `z_cap10` ($H^2/(100+H^2)$) | `8.0%` *(tax)* | **`52.7%`** | **`100.0%`** | **`98.0%`** |
| `z_linear` ($w = H$) | `9.7%` | `28.0%` *(tax)* | **`100.0%`** | **`99.0%`** |
| `z_pow2` ($w = H^2$) | `6.7%` | `23.7%` | `99.5%` | **`99.0%`** |

#### Table C3 — Push 4C: Ghost-Free `Z-GCLR+` on `1x r=20 (m=2) + 500 random noisy zeros`

| Method | 95% interval on `1x r=20 + 500 nulls` | Width | Anti-parallel `(2,-1), (1,-2)` (`RP^1`) |
|---|---|---:|---|
| Unweighted `GCLR+` | `[1.50, 3.15]` | `1.65` | `EMPTY` |
| `Z-GCLR+` ($H/(1+H)$) | `[1.65, 3.00]` | `1.35` | `EMPTY` |
| `Z-GCLR+` ($H^2/(100+H^2)$) | **`[1.65, 2.55]`** | **`0.90`** | `EMPTY` |
| `Z-GCLR+` ($w = H$, 0 params) | **`[1.65, 2.55]`** | **`0.90`** | `EMPTY` |

#### Table C4 — Push 5A: Closing the 1D $N=2$ Door B Gap ($\nu = 10$, exact push-free $q_{0.95}^{\text{t-AR}} = 6.939$)

| 1D $N=2$ Scene | `t-AR` (exact push-free) | **Piecewise Exact Door B (`AR-bounded t-LR`)** | Piecewise Door B $\cup$ `Z-GCLR+` | Naive Wilks `t-LR` *(undercovers)* |
|---|---|---|---|---|
| All-zero `0/0` | `(0, inf)` | **`(0, inf)`** | `(0, inf)` | `(0, inf)` |
| Scene A (`5.5, 5.6 ± 0.1`) | `[5.360, 5.740]` | **`[5.353, 5.734]`** | `[5.353, 5.734]` | `[5.407, 5.693]` |
| Scene B 10σ (`5.0, 6.0 ± 0.1`) | `EMPTY` | **`[4.854, 6.154]`** | `[4.854, 6.154]` | `[4.907, 6.094]` |
| Scene B 50σ (`5.0, 6.0 ± 0.02`) | `EMPTY` | **`[4.963, 5.067] U [5.937, 6.041]`** | `[4.963, 5.067] U [5.937, 6.041]` | `[4.967, 5.047] U [5.953, 6.034]` |

#### Table C5 — Push 5B: Honest Matched vs. Misspecified Coverage on 1D $N=2$ (declared $\nu = 10$, 400 reps)

| True World (1D $N=2$) | `t-AR` ($q_{95}=6.94$) | **AR-bounded `t-LR`** (proved $\ge 95\%$ if matched) | Naive Wilks `t-LR` (`3.841`) |
|---|---:|---:|---:|
| True Gaussian, $r=0$ | `99.0%` | **`99.8%`** | `96.2%` |
| True $\nu=10$, $r=0$ *(matched)* | `94.0%` | **`98.2%`** | `90.0%` |
| True $\nu=10$, $r=5$ *(matched)* | `94.5%` | **`99.2%`** | `94.5%` |
| True $\nu=2$, $r=0$ *(misspecified!)* | `76.0%` | `85.5%` *(relies on declared $\nu$)* | `74.5%` |
| True $\nu=2$, $r=3$ *(misspecified!)* | `75.5%` | `89.0%` *(relies on declared $\nu$)* | `76.8%` |

---

### Part D — Push 6: Metric Mass $g \in \mathrm{Sym}^+_d$ (`python3 explore/e13_metric_mass.py`, seed 20261122)

#### Table D1 — Push 6B: Three-Rung Nested Surprise Ladder (`scalar ⊂ metric Sym+_2 ⊂ linear M`)

| 2D World (4 clean pushes, $\sigma = 0.1$) | $K_{\min}^{\rm scalar}$ | $K_{\min}^{\rm metric}$ | $K_{\min}^{\rm linear}$ | $\Delta K_{\rm iso}$ (`O17`) | $\Delta K_{\rm recip}$ (`E7`) |
|---|---:|---:|---:|---:|---:|
| 1. Isotropic scalar $g = 2 I_2$ | `0.0` | `0.0` | `0.00` | `0.0` | `0.0` |
| 2. Anisotropic metric $g = \mathrm{diag}(2, 5)$ (`O17`) | `59.9` | **`0.0`** | `0.00` | **`59.9`** | `0.0` |
| 3. Gyroscopic non-reciprocal $M = \begin{pmatrix}2&-1\\1&2\end{pmatrix}$ (`E7`) | `68.6` | `66.2` | **`0.00`** | `2.5` | **`66.2`** |
| 4. Anti-parallel negative map $M = -\mathrm{diag}(2, 3)$ | `400.2` | `400.0` | `0.00` | `0.2` | `400.0` |

*Full-rank (`rank = 3`) 2D anti-parallel hilltop ghost at $g_0 = I_2$ ($K(I_2) = 3600$):* `T50` gives `T2 = 0.0000, p = 1.0000` *(ghost accepted)*; `Metric-GCLR(Sym_2)` gives `lr_obs = 3200.0, p = 0.0100`; `Metric-Z-GCLR(Sym_2)` gives `lr_obs = 3152.8, p = 0.0100` *(both kill the ghost)*.

#### Table D2 — Push 6C: Monte Carlo 95% Coverage on 2D Anisotropic Metric $g = \mathrm{diag}(2, 5)$ (`N = 8`, 100 reps)

| World | `T50` score | `Metric-GCLR+` | `Metric-Z-GCLR+` | `Metric-SCLR+` | `O7` Metric Floor | Raw 3-df Wilks *(contrast)* |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian pure noise (`r = 0`) | `96.0%` | `89.0%` | **`95.0%`** | **`96.0%`** | `100.0%` | `65.0%` |
| Gaussian partial $x$-axis only (`r = 4`) | `97.0%` | `94.0%` | **`94.0%`** | **`98.0%`** | `99.0%` | `88.0%` |
| Gaussian weak (`r = 1.5`) | `94.0%` | `94.0%` | **`97.0%`** | **`96.0%`** | `99.0%` | `84.0%` |
| Gaussian strong (`r = 5`) | `93.0%` | `91.0%` | **`93.0%`** | `90.0%` | `96.0%` | `92.0%` |
| Scale-mix $\nu = 10$ weak (`r = 1.5`) | `86.0%` | `84.0%` | `82.0%` | **`93.0%`** | **`96.0%`** | `70.0%` |
| Heavy-tail $\nu = 2$ (`r = 3`) | **`41.0%`** | **`42.0%`** | **`40.0%`** | **`93.0%`** | **`94.0%`** | `36.0%` |

#### Table D3 — Push 6D–E: Simultaneous `L5`/`L14` Tensor Composition & Exact `L22` Isotropy Test

- **Simultaneous `L5` & `L14`** (`g_1 = diag(2, 5)`, `g_2 = diag(3, 1)`):
  - `L5` lumping: $K_{\rm lump}(\mathrm{diag}(5.0, 6.0)) = 0.000000$; $K_{\rm lump}(\mathrm{diag}(5.4, 6.0)) = 7.439$.
  - `L14` common force: $K_{\rm red}(\mathrm{diag}(6/5, 5/6)) = 0.000000$; $K_{\rm red}(\mathrm{diag}(1.5, 5/6)) = 57.925$.
- **Partial-axis lumping** (`diag(2, 5)` + `g2_xx = 3`, $g_{2,yy}$ unmeasured = $\sigma_\infty$): $K_{\rm lump}(\mathrm{diag}(5.0, G_{yy}))$ equals `20.790` at `4.50`, `3.251` at `4.80`, and **`0.000`** at `5.05`, `6.00`, and `20.00`.
- **Exact `L22` scalar-isotropy test** (`N = 8, r = 15`): Isotropic `2 I_2` -> `Z-GCLR p = 0.980, SCLR p = 0.960`; Anisotropic `diag(2, 5)` -> `Z-GCLR p = 0.010, SCLR p = 0.010`; $x$-only push -> `p = 1.000` (both).

---

## Summary matrix: how the `O7` noise-scale declarations map onto unnormalized evidence + exact sets

| `O7` declaration of $\Sigma_i$ | Compatibility curve | Exact 95% set (zero $\pi$, zero $\eta$) | Passes exact `Z`? | Automatic widening on discordance? | Honest limitation |
|---|---|---|---|---|---|
| **Row 1: $\Sigma_i$ exact Gaussian** | `T36` $K(m)$ + $K_{\min}$ surprise | `Z-GCLR` / `Z-GCLR+` (or unweighted `GCLR`) | `Z-GCLR`: **yes** (`GCLR`: no) | No (Kac–Bernstein; surprise lives in $K_{\min}$) | `Z-GCLR` has rank-1 pinhole ghost; $w(H)$ is a 3-regime Pareto family |
| **Row 2: $\Sigma_i$ up to per-reading scale** | `T36` / `T48` shape | `SCLR+` (spherical conditional LR) | **yes** (bit-for-bit) | **yes** (exact $\ge 95\%$ for *all* $v_i > 0$) | Shrinks past $\sigma$ on over-agreeing data; 1D $N \le 4$ is `(0, inf)` |
| **Row 3: $\Sigma_i$ is a floor ($v_i \ge 1$)** | `T36` $K(m)$ + $K_{\min}$ | `C_floor = Z-GCLR+ ∪ SCLR+` | **yes** (bit-for-bit) | **yes** (and respects $\sigma$ floor!) | Slightly conservative (~97–98%) at $v_i \equiv 1$; 1D $N \le 4$ is `(0, inf)` |
| **Door B (declared $t_\nu$)** | `T48` Student-$t_\nu$ | Piecewise Exact Door B (`AR-bd t-LR` on 1D $N\le4$, `SCLR+(t_\nu)` else) | $N\ge5$ or $k\ge2$: **yes**; 1D $N\le4$: no | **yes** (even on 1D $N=2$ Scene B!) | 1D $N\le4$ branch depends on declared $\nu$ and count $N$ |

## What this exploration does NOT solve (where math ends and rulings begin)

1. **The output-type & attitude fork (Can / Claude):** changing the output from a normalized area-1 probability law + point (`L2`, `G28`) to an unnormalized residual curve $K(m)$ + exact compatibility set (allowed to return `(0, inf)`) plus a separate leftover surprise $K_{\min}$ is a ruling about the goal, not a calculation — and choosing among `O7`'s rows (exact Gaussian vs. up-to-scale vs. floor vs. declared $t_\nu$) is an instrument declaration.
2. **Selecting a weight inside the `Z`-safe Pareto family (`E2`, `E3`):** Insight 9 proves the family cannot collapse to one dominant weight without preferring one signal geometry over another.
3. **Minimal streaming state for correlated weak-signal composition:** Push 2 shows exact two-body composition needs more than the 1D curves $K_1(m), K_2(m)$ (on equal-ratio sensors the $2 \times 2$ Gram + count suffices; on arbitrary heterogeneous $\Sigma_i$, the $2k$-dim within-span reading summaries are kept).
