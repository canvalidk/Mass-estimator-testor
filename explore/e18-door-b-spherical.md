# E18: Walking All the Way Through Door B — Dropping Gaussianity, Spherical Non-Uniqueness, and the Door A vs. Door B Verdict

**Branch:** `explore/e18-door-b-spherical` (from `origin/main` `62e0390`, 24 Nov 2026; Antigravity with Can).  
**Status:** Exploration only (`explore/README.md` rules). Nothing in `definition/`, `problems/`, `rulings.md`, or `met/` is modified.  
**Goal:** Give **Door B** (dropping `L13` Gaussianity while keeping believed $\Sigma_i$ `L18`, so disagreement genuinely widens the compatible set without revising $\sigma$) the same full uniqueness and `definition/proof.md` audit that `explore/e15-door-a-uniqueness` gave Door A — and put Door A and Door B side by side with zero hand-waving.

---

## Standard Mathematics First (`GOAL.md` Rule 6)

- **Spherical / permutation randomization tests:** Fisher (1935), *The Design of Experiments*; Lehmann & Romano (2005), *Testing Statistical Hypotheses*, §15.2; Kariya & Eaton (1977), *Ann. Stat.* 5:206–215 (robustness of invariant tests under spherical symmetry).
- **Exact-variance / unimodal concentration bounds:** Chebyshev (1867); Gauss (1823), *Theoria Combinationis Observationum* ($P(|e| \ge t\sigma) \le \frac{4}{9 t^2}$ for every symmetric unimodal law with exact variance $\sigma^2$).
- **Self-normalized / heterogeneous-variance Student-$t$ robustness:** Bakirov & Székely (2005), *J. Math. Sci.* 131:5350–5356; Ibragimov & Müller (2010), *JBES* 28:453–468, Theorem 1 (exact conservativeness of equal-weight Student-$t_{N-1}$ under arbitrary heterogeneous Gaussian scales $\sigma_1, \dots, \sigma_N$).
- **Non-existence of UMPU tests in nonparametric location/spherical families:** Lehmann & Romano (2005), §§3.2, 4.3 (the locally most powerful directional score depends explicitly on the radial score $-g'(\rho^2)$).

---

## 1. Theorem E18-1 (The Non-Gaussian World Hierarchy & Universal Small-Scalar Impossibility)

Dropping `L13` (Gaussianity) while retaining believed $\Sigma_i$ (`L18`) yields a nested hierarchy of worlds around the declared covariances:
$$\mathcal W_{\rm Gauss}(\Sigma) \;\subset\; \mathcal W_{\rm ell,\Sigma}(\Sigma) \;\subset\; \mathcal W_{\rm sph}(\Sigma) \;\subset\; \mathcal W_{\rm cov}(\Sigma),$$
where $\mathcal W_{\rm ell,\Sigma}$ is independent $\Sigma_i$-spherically symmetric errors with **exact covariance $\operatorname{Cov}(e_i) = \Sigma_i$** (standardized $t_\nu$ for $\nu > 2$, Laplace, uniform ball, arbitrary $E[v_i]=1$ scale mixtures), $\mathcal W_{\rm sph}$ allows **arbitrary unknown radial laws per reading** (`T75 (4)`, including Cauchy $\nu = 1$ and arbitrary unknown per-reading scales $v_i > 0$), and $\mathcal W_{\rm cov}$ allows arbitrary zero-mean asymmetric laws with $\operatorname{Cov}(e_i) = \Sigma_i$.

1. **Exact-variance symmetric non-Gaussian noise ($\mathcal W_{\rm ell,\Sigma}$) immediately breaks Door A (`T45`):**
   On $N$ direct scalar readings with exact $\sigma = 1$, take the symmetric 3-point law $P(e = 0) = 1 - p$, $P(e = \pm 1/\sqrt p) = p/2$ ($E[e] = 0$, $\operatorname{Var}(e) = 1$ **exact** for every $p \in (0, 1]$). Setting $p = 0.99 / (N z_{1-\alpha/2}^2)$ makes any single nonzero draw push $|\bar e| > z_{1-\alpha/2}/\sqrt N$, so **`T45`'s 95% coverage collapses to `74.2%` ($N=1$), `76.8%` ($N=2$), `78.0%` ($N=10$), `78.5%` ($N=50$)** (verified in `tests/test_e18_door_b_spherical.py`).
2. **Universal Small-Scalar Impossibility in $\mathcal W_{\rm sph}$ (extends `T75 (7)` from LR to *every* test!):**
   Let $\psi \in \{0, 1\}$ be **any** non-randomized test with uniform size $\le \alpha$ over $\mathcal W_{\rm sph}$ on $N$ scalar readings. On the direct-reading limit ($\sigma_a \to 0$, $a_i = 1$, $F_i = m_0 + e_i$), for *any* fixed positive magnitudes $(|v_1|, \dots, |v_N|) \in (0, \infty)^N$, the discrete 2-point symmetric law $P(e_i = +|v_i|) = P(e_i = -|v_i|) = \frac12$ belongs to $\mathcal W_{\rm sph}$. Under this law $(e_1, \dots, e_N)$ is uniform over the $2^N$ sign patterns $(\pm|v_1|, \dots, \pm|v_N|)$, forcing
   $$\frac{1}{2^N} \sum_{s \in \{\pm 1\}^N} \psi(s_1|v_1|, \dots, s_N|v_N|) \le \alpha.$$
   - For **$N \le 4$** at $\alpha = 0.05$, $1/2^4 = 0.0625 > 0.05$, forcing $\psi \equiv 0$ identically!
   - If $\psi$ is two-sided symmetric ($\psi(-e) = \psi(e)$, e.g. `T5` swap-equivariant or unbiased for both sides), the sum is even, forcing $2/2^N \le \alpha \implies$ **$N \ge 6$** at $\alpha = 0.05$ (`2/2^5 = 0.0625 > 0.05`).
   - Thus Can's 2-reading scalar scenes (`near` and `far`) MUST return `(0, ∞)` under **every** exact non-randomized test over $\mathcal W_{\rm sph}$.
3. **What exact $\operatorname{Cov}(e_i) = \Sigma_i$ ($\mathcal W_{\rm ell,\Sigma}$) unlocks on 2 scalar readings:**
   By independence + Chebyshev (or Gauss's 1823 unimodal inequality $P(|e| \ge t\sigma) \le \frac{4}{9t^2}$), $P(\min(|e_1|, |e_2|) \ge t\sigma) \le t^{-4}$ (or $\frac{16}{81}t^{-4}$ unimodal), giving exact 95% thresholds $t = 20^{1/4} = 2.115$ (general) and $t = 1.410$ (unimodal):
   - On Can's `near` `(5.5, 5.6) ± 0.1`: single tight interval **`[5.289, 5.811]`** (`[5.359, 5.741]` unimodal).
   - On Can's `far` `(5.0, 6.0) ± 0.1`: **two disconnected islands `[4.789, 5.211] ∪ [5.789, 6.211]`** (`[4.859, 5.141] ∪ [5.859, 6.141]` unimodal), **rejecting the midpoint `5.5`** — the exact signature of per-reading heavy tails (`O14`, 23 Nov: one reading was an outlier, not both!).

---

## 2. Theorem E18-2 (Shattering of `F2` Uniqueness Inside Door B — Exact Counterpart to `E15-1`)

Condition on $\mathcal C = (U_1, \dots, U_n,\; r_1, \dots, r_n)$ where $r_i = \|V_i\|$ (`T75 (4)`). Under $H_0(m_0)$ across **all** of $\mathcal W_{\rm sph}$, the unit directions $\hat V_i = V_i / r_i$ are independent $\operatorname{Unif}(S^{k_i - 1})$ (exact for *every* spherically symmetric density $f_i(\|\varepsilon_i\|^2) = \exp(-g_i(\|\varepsilon_i\|^2))$, not only Gaussian scale mixtures).

**Exact Local Directional Score Derivation:**
Under a local mass shift $m_1 = m_0 + \delta$ ($\delta \to 0$), $\mu_i(m_1) = \mu_{U,i} + \delta\,P_i'(m_0)\mu_{U,i} + O(\delta^2)$, so at fixed $(U_i, r_i)$:
$$\|z_i - \mu_i(m_1)\|^2 = \underbrace{\|U_i - \mu_{U,i}\|^2 + r_i^2}_{=: \rho_i^2 \text{ (fixed by } \mathcal C\text{)}} - 2\delta\, V_i^\top (P_i'(m_0)\mu_{U,i}) + O(\delta^2).$$
Differentiating the conditional log-density $-\sum_i g_i(\|z_i - \mu_i(m_1)\|^2)$ at $\delta = 0$ yields:
$$\left.\frac{\partial}{\partial \delta}\right|_{\delta=0} \log p_\delta(\hat V_1, \dots, \hat V_n \mid U, r) \;=\; 2 \sum_{i=1}^n g_i'(\rho_i^2)\, V_i^\top \big(P_i'(m_0) \mu_{U,i}\big).$$

**Consequences (Direct Opposite of Door A's Theorem `E15-1`!):**
1. **No UMPU or LMP test exists in $\mathcal W_{\rm sph}$:** Because dropping `L13` leaves the radial log-density $g_i$ free, **every** positive radial weight $w(\|V_i\|)$ defines a distinct, exact, admissible spherical-null test $T_w(m_0) = \big|\sum_i w(\|V_i\|)\,V_i^\top (P_i'(m_0) U_i)\big|$ that is **uniquely locally most powerful** at its own radial family (at strong signal $\mu_{U,i} \approx U_i$):
   - Gaussian ($g' \equiv \frac12$) $\longleftrightarrow$ unweighted spherical score $w \equiv 1$ (`unw_score`, Fisher 1935);
   - Cauchy ($g'(t) \propto \frac{1}{1+t}$) $\longleftrightarrow$ redescending score $w(r) = \frac{1}{1 + r^2}$ (`cauchy_score`);
   - Laplace ($g'(t) \propto \frac{1}{\sqrt t}$) $\longleftrightarrow$ spatial-sign score $w(r) = \frac{1}{r}$ (`sign_score`; also discontinuous at exact residual $V_i = 0$!);
   - Plus the global spherical LR (`T75 (4)` `sph_lr_pos`).
2. **Severe real-data disagreement (breaks `F2`):** On **10 noisy agreeing 2D vector readings at $m = 5.0$ ($\pm 0.3$) + 1 wild outlier at $m = 50.0$** (`T81` point $\arg\min K = 13.78$, verified in `tests/test_e18_door_b_spherical.py`):

   | Candidate $m_0$ | `unw_score` ($w\equiv 1$) | `cauchy_score` ($1/(1+r^2)$) | `sign_score` ($1/r$) | `sph_lr_pos` (`T75 (4)`) |
   | --- | --- | --- | --- | --- |
   | **True $m = 5.0$** | `0.0645` (barely kept) | **`0.6406`** (peak) | **`0.6921`** (peak) | `0.0645` (barely kept) |
   | **$m = 6.5$** | **`0.3449` (KEPT)** | **`0.0007` (REJECTED)** | `0.2452` (KEPT) | **`0.3449` (KEPT)** |
   | **`T81` point $m^* = 13.78$** | **`1.0000` (KEPT)** | **`0.0342` (REJECTED)** | `0.2472` (KEPT) | **`1.0000` (KEPT)** |
   | **$m = 25.0$** | **`0.0170` (REJECTED)** | **`0.3877` (KEPT — island!)** | `0.2532` (KEPT) | **`0.0170` (REJECTED)** |

   `unw_score`/`sph_lr_pos` return one wide compromise around `13.78`; `cauchy_score` rejects `6.5` AND rejects `T81`'s point `13.78`, returning two disconnected islands around `5.0` and the outlier; `sign_score` keeps all four! Even forcing `T81` point-inclusion (`w ≡ 1`) leaves `unw_score` vs `sph_lr_pos` unpinned (and `T81`'s point itself breaks under $2 < \nu < 4$ heavy tails by `E17-G3`).

---

## 3. Theorem E18-3 (Three Structural Walls of Spherical `T75 (4)`, Floor `T75 (6)`, and `t_joint`)

1. **Spherical LR (`T75 (4)`) CANNOT kill anti-aligned hilltops — even on `domain="line"` (allowing rejection of Newton!):**
   Because conditioning on $\|V_i\|$ freezes $Q_{VV} = \sum \|V_i\|^2$ identically across all spherical null draws:
   - On anti-aligned 2D `(10, -10)` at the global maximum $m_0 = 1.0$: $U_i = 0$ for every reading, so `LR_line = Q_VV` is identical for every draw $\implies$ **`p_line(1.0) = 1.0000`**, while `p_line(0.5) = p_line(2.0) = 0.0003`! Even dropping `L1` ("allowing the empty set to reject Newton") leaves an isolated hilltop ghost at the worst-fitting mass $m = 1.0$!
   - On anti-aligned 2D `(10, -3) & (10, -3.1)` at the global maximum $m_0 = 0.305$ ($K = 218.6$), **`sph_lr_pos(0.305) = 0.5004`** (while `sph_lr_pos(1.0) = 0.0002`), and `T75 (6)` ("the floor") inherits this hilltop ghost through the union!
2. **Pure Spherical `T75 (4)` violates believed $\sigma$ (`L18`) on over-agreeing readings, and Floor `T75 (6)` re-imports `L13` without uniqueness:**
   Because $\mathcal W_{\rm sph}$ discards the absolute scale of $\Sigma_i$, on four 2D readings agreeing to `±0.002` when believed Gaussian $1\sigma \approx 0.28$, `sph_lr_pos` rejects `m0 = 5.505` (only `0.02σ` from the mean!) at **`p = 0.0013`**. Taking the union with Gaussian (`T75 (6)`) repairs this — **but** re-imports `L13` Gaussianity, over-covers (`0.97–0.995`), and opens a massive `F2` non-uniqueness (any valid test on $\mathcal W_{\rm Gauss}$ unioned with any valid test on $\mathcal W_{\rm sph}$ gives a valid floor).
3. **Why `t_joint` (`T82` / `ours34`) fails inside $\mathcal W_{\rm sph}$ despite Bakirov–Székely (2005) / Ibragimov–Müller (2010):**
   For equal-weight direct scalar means $X_i \sim \mathcal N(\mu, \sigma_i^2)$ with arbitrary unknown per-reading scales $\sigma_i$, Student-$t_{N-1}$ keeps exact $\ge 95\%$ coverage (Bakirov–Székely 2005; Ibragimov–Müller 2010 Theorem 1 — verified: `0.977` under Cauchy $\nu = 1$ scales at $N = 5$). In EIV (`T45`), reading $i$'s leverage $U_i / \sqrt{\sum U_j^2}$ scales *with* $\sqrt{v_i}$: when one reading $i^*$ draws a huge scale $v_{i^*} \gg 1$, both $U_{i^*}$ and $V_{i^*}$ blow up by $\sqrt{v_{i^*}}$, driving $|Z| \approx |V_{i^*}| \sim \sqrt{v_{i^*}}$ while $W = K - Z^2 = \sum_{j \neq i^*} V_j^2 = O(1)$ stays small — collapsing `t_joint` 95% coverage to **`~0.77`** at $\nu = 1, N = 5$!

---

## 4. Theorem E18-4 (Complete Side-by-Side `definition/proof.md` Verdict: Door A vs. Door B)

Walking all the way through **both** doors (`E15`, `E17`, `E18`) turns the project's central trade-off into a completely transparent ledger:

| Proof / Premise Criterion | **Door A (`E15` + `met/ours35.py`)** | **Door B — Spherical (`T75 (4)` / `(6)`)** | **Door B — Common Scale (`T82` / `met/ours34.py`)** |
| --- | --- | --- | --- |
| **`G1` (`L13` Gaussianity)** | **Required** (breaks to `74–78%` in $\mathcal W_{\rm ell,\Sigma}$) | **DISSOLVES** (exact on all of $\mathcal W_{\rm sph}$, unless Floor re-imports `L13`) | Required per-reading (breaks to `~0.77` on per-reading heavy tails) |
| **`F2` Uniqueness (`G26`)** | **PROVED** on `(K, C)` (`E15-1`/`E15-2`) | **SHATTERS** (`E18-2`: infinite LMP continuum disagreeing sharply on real data) | Proved only on common-scale Gaussian world (`T82 (4)`) |
| **`R-2026-11-17-3` disagreement widening & `T30`** | **FAILS** (puts contradiction in leftover $K_{\min}$ only) | **PASSES** on 2D/3D vectors (`T75 (6)`) | **PASSES** (at cost of widening 88% of honest pairs) |
| **Can's 2-scalar scenes (`near` vs `far`)** | `near`: `±0.14`, `far`: `±0.14` (`K_min: 0.5 → 50`) | **UNIVERSALLY VACUOUS `(0, ∞)`** on `< 6` scalars (`E18-1`) | `near`: `±0.64`, `far`: `(0, 11.85)` |
| **Believed $\sigma$ (`L18`) on over-agreement** | **PASSES** | Pure `T75 (4)` **FAILS** (`0.02σ` rejected); Floor `T75 (6)` passes via `L13` | **PASSES** (via `T45` floor) |
| **`L25 (b)` fixed-type additive state & `L17`** | **PASSES** (`State = (K, C)`, zero MC) | **FAILS** (needs growing bag of readings + MC; `E17-L2`) | **FAILS** (needs `M_nz`, breaks `0/0` continuity & `L17`) |
| **Continuity at `0/0`** | **PASSES** | **FAILS in 1D** (`0.0625 → 0.03125`); passes in 2D/3D (`E17-L2`) | **FAILS** (`M_nz` jump) |
| **Hilltop ghosts on anti-aligned `(10,-10)` / `(10,-3)`** | **FAILS** (`E15-1`, `E15-3`) | **STILL FAILS!** (`p_line(1.0) = 1.0`, `p_pos(0.305) = 0.50`; `E18-3`) | **FAILS** (`T45` floor keeps hilltops) |
| **`G27` two-body set composition** | 2-dof `+24.9%` tax from `(K, C)` | 1-dof exact via spherical `T74 (6)` (raw readings + MC) | 2-dof `Q_joint` |
| **`T81` point consistency (`G25`)** | **PROVED** globally for all finite-4th-moment laws (`E17-G1`) | `T81` **BREAKS** for $2 < \nu < 4$ (`E17-G3`); robust point breaks under skewness | Same as Door A |

**Bottom line:** Dropping Gaussianity (Door B) genuinely solves `G1` and delivers disagreement widening without doubting $\sigma$, **and simultaneously proves why no unique (`F2`), fixed-state (`L25 (b)`), small-scalar-capable confidence set exists once Gaussianity is dropped** — reinforcing the core finding from `E15`–`E17`: the clean physical object is the unnormalized evidence/misfit geometry, while any frequentist confidence-set wrapper forces a sharp choice between Door A's uniqueness (`F2` + `L25 (b)`) and Door B's non-Gaussian robustness (at the cost of an unpinned radial family `E18-2`, small-scalar vacuity `E18-1`, and raw-reading MC state).

## Verification

```bash
python3 -m unittest -v tests/test_e18_door_b_spherical.py
```
