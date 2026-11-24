# E18: Walking All the Way Through Door B — Dropping Gaussianity, Spherical Non-Uniqueness, and Door A vs. Door B Verdict

Branch: `explore/e18-door-b-spherical` (from `origin/main` `62e0390`, 24 Nov 2026). Exploration only; nothing in `definition/`, `problems/`, `rulings.md`, or `met/` modified.

**Core target:** Walk all the way through **Door B** (drop `L13` Gaussianity while keeping believed $\Sigma_i$ `L18`, testing $\mathcal W_{\rm ell,\Sigma}$ and $\mathcal W_{\rm sph}$ / `T75 (4)`–`(6)`) to determine whether Door B can become provably unique (`F2`) and what holes in `definition/proof.md` it closes vs. opens.

## Four theorems (for Can & Claude to harvest)

1. **`E18-1` (Non-Gaussian Hierarchy & Universal Small-Scalar Impossibility):**
   - Dropping `L13` to exact-variance symmetric laws ($\mathcal W_{\rm ell,\Sigma}$, exact $\sigma=1$ 3-point law) immediately collapses Door A `T45`'s 95% coverage to **`74.2%`–`78.5%`** for all $N \in \{1, 2, 10, 50\}$.
   - **Universal small-scalar impossibility on $\mathcal W_{\rm sph}$:** Using 2-point symmetric laws $P(e_i = \pm|v_i|)=1/2 \in \mathcal W_{\rm sph}$, **NO** non-randomized exact test on $\mathcal W_{\rm sph}$ (not just `T75 (4)` LR) can reject anything for $N \le 4$ direct scalar readings ($1/2^4 > 0.05$) or $N \le 5$ under two-sided symmetry ($2/2^5 > 0.05$).
   - Exact-variance $\mathcal W_{\rm ell,\Sigma}$ moment bounds (Chebyshev-min $2.115\sigma$ / Gauss-unimodal $1.410\sigma$) work on 2 scalars: tight single interval on Can's `near` (`5.55 ± 0.26`/`0.19`) and **two disconnected islands rejecting `5.5`** on `far` (`[4.79, 5.21] ∪ [5.79, 6.21]`).
2. **`E18-2` (Shattering of `F2` Uniqueness Inside Door B — Counterpart to `E15-1`):** Differentiating the conditional spherical density $\log p_\delta(\hat V \mid U, r)$ gives exact LMP score $2\sum g_i'(\rho_i^2)\, V_i^\top(P_i'\mu_{U,i})$. Because $g_i$ is free once `L13` is dropped, **no UMPU or LMP test exists in $\mathcal W_{\rm sph}$**: every radial weight $w(\|V_i\|)$ is exact, admissible, and uniquely LMP for its own radial law (Gaussian $w\equiv 1$, Cauchy $1/(1+r^2)$, Laplace spatial-sign $1/r$, plus `sph_LR`). On 10 noisy agreeing 2D readings at $m=5$ + 1 outlier at $m=50$ (`T81` $m^* = 13.78$): `unw_score`/`sph_LR` keep `6.5` (`p=0.345`) and `13.78` (`p=1.0`), reject `25` (`p=0.017`); `cauchy_score` centers at `5.0` (`p=0.641`), **rejects `6.5` (`p=0.0007`) and rejects `T81`'s `13.78` (`p=0.034`)**, and keeps `25` (`p=0.388`); `sign_score` keeps all four!
3. **`E18-3` (Three Structural Surprises of Spherical `T75 (4)`, Floor `T75 (6)`, & `t_joint`):**
   - **Spherical LR (`T75 (4)`) CANNOT kill anti-aligned hilltops — even on `domain="line"` (allowing rejection of Newton!):** Conditioning on $\|V_i\|$ freezes $Q_{VV}$ in every null draw, so on anti-aligned 2D `(10, -10)` at the hilltop $m=1.0$, **`p_line(1.0) = 1.000`** (while `p_line(2.0) = 0.0003`), and on `(10,-3)&(10,-3.1)` **`p_pos(0.305) = 0.500`** (inherited by Floor `T75 (6)`!).
   - **Pure spherical violates believed $\sigma$ (`L18`) on over-agreement:** rejects `0.02σ` offset at `p = 0.0013`; Floor `T75 (6)` repairs this only by re-importing `L13` Gaussianity without uniqueness.
   - **EIV leverage breaks Bakirov–Székely (2005) / Ibragimov–Müller (2010) for `t_joint`:** equal-weight Student-$t$ holds `>= 0.95` under Cauchy $\nu=1$ per-reading scales, while EIV `t_joint` drops below `0.85` because outlier leverage $U_i$ scales with $\sqrt{v_i}$.
4. **`E18-4` (Door A vs. Door B `definition/proof.md` Verdict):** Full side-by-side table showing the exact structural impossibility: Door A proves `F2` uniqueness + `L25 (b)` + small-sample sharpness at the cost of `L13` and no set widening; Door B dissolves `G1` (`L13`) and delivers vector disagreement widening at the provable cost of shattering `F2`, small-scalar vacuity, losing fixed-type `L25 (b)` state, and still keeping anti-aligned hilltops.

## Files & verification

- Full proofs & Door A vs. Door B audit: `explore/e18-door-b-spherical.md`
- Pure-stdlib verification library: `explore/e18_door_b.py`
- 6 executable checks: `python3 -m unittest -v tests/test_e18_door_b_spherical.py`
