# E17: Door-Independent Foundations — `L25` Algebraic Structure & `G25` Non-Gaussian Vector Consistency of `T81`

Branch: `explore/e17-l25-and-g25` (from `origin/main` `62e0390`, 24 Nov 2026). Exploration only; nothing in `definition/`, `problems/`, `rulings.md`, or `met/` modified.

**Core target:** Prove the two highest-leverage **door-independent** results — (I) the algebraic structure forced by order- and grouping-independent combination (`L25`), and (II) global full-type vector consistency (`G25` / `O19`) of `T81`'s $\arg\min_{[0,\infty]} K$ under **non-Gaussian** world errors.

## Five theorems (for Can & Claude to harvest)

1. **`E17-L1` (Aczél Curve-Addition Uniqueness):** Any continuous, strictly monotone, associative pointwise curve-combination rule with `0/0` identity (`L25(b) + L20`) is uniquely isomorphic (up to a positive scale) to **unnormalized addition** $K_{A\cup B} = K_A + K_B$ (`Aczél 1966`). Commutativity is automatic; normalizing (`ours30`–`ours33`) breaks `L25(b)`. Holds for both Door A ($K_i$) and non-Gaussian Door B curves ($\rho(K_i)$).
2. **`E17-L2` (`L25` Continuity-and-Quadraticity Wall & Exact Price of Door B):**
   - No continuous homomorphism of *any* continuous `L25(b) + L20` monoid can recover $N$, $D$, or $M_{\rm nz}$ (generalizing `T82 (2)` to all topological monoids).
   - Any $C^2$-at-0 additive state satisfying `L25(b) + L20 + L16` (frame sign-flip $\Phi(-y)=\Phi(y)$) is automatically quadratic near `0/0` ($D\Phi(0)=0$). Adding reading-splitting (`L17`: $n\,\Phi(y/\sqrt n)=\Phi(y)$) **forces the Gaussian Gram matrix / curve $K(m)$ WITHOUT assuming `L13`!**
   - Exact algebraic price of Door B: any disagreement-widening procedure must drop `L17` (non-quadratic additive curves, unpinned $\nu$), `0/0` continuity (`t_joint`), or fixed-type monoids (`T75 (4)` spherical — proved & verified **discontinuous at `0/0` in 1D** where `5×(5,1)+(1e-6,0)` drops `p: 0.0625 → 0.03125`, yet **continuous at `0/0` in 2D/3D**!).
3. **`E17-G1` (Global Non-Gaussian Full-Type Vector Consistency of `T81`):** Under independent arbitrary non-Gaussian errors with exact covariance $\Sigma_i$ (`E[K(m)-K(m0)] == G_n(m)` exact zero bias!) and uniform trace-free 4th-moment bound $M_4$ (all bounded-kurtosis laws in $d \le 3$; elliptical $\nu > 4$ in any $k_i$), with $\kappa(\Gamma_i) \le \kappa_0$, `T81`'s **global** $\arg\min_{[0,\infty]} K$ attains `B27`'s exact sharp rate $|\widehat\phi_n - \phi_0| = O_p(\sqrt{M_4 D_n + \Delta_n}/\Delta_n)$ whenever $\Delta_n^2/(D_n + \Delta_n) \to \infty$ (verified on correlated 3D $k=2$ + 2D + 1D series with skewed `Exp(1)-1`, asymmetric discrete, and $t_5$ errors). Closes `definition/proof.md` lines 89–90.
4. **`E17-G2` (Unbounded Ratios $s_i$: Layer-Cake/Doob Bound + Global Ray Counterexample):** 3-line Layer-Cake + Doob $L_q$ maximal inequality extends `B27`'s compact-interval $J^2/(W+J) \to \infty$ to non-Gaussian errors with zero log factors, **AND exposes a new global-ray counterexample:** on $m_0 = 1$, 1 reading ($s_1 = 1, r_1^2 = 8\log N$) + $N$ noise readings at $s_i = N$, `B27`'s $J^2/(W+J) \to \infty$ holds (compact accuracy $\to 1$), yet global $\arg\min_{[0,\infty]} K$ **escapes to $\ge N/2$ with probability $\to \ge 1/2$** (`0.217 → 0.387 → 0.470`). Gives exact global chord bound $\sin^2(\theta_i-\theta_{i,0}) \ge w_i^2\rho(m,m_0)^2$ and global sufficiency $J_n/\sqrt{M_4 D_n + \Delta_n} \to \infty$.
5. **`E17-G3` (Infinite 4th Moment $2 < \nu < 4$, Exact Finite $\Sigma_i$ — Link to Door B):** By von Bahr–Esseen (1965) + Doob $L_q$ ($q = \nu/2$), `T81`'s threshold degrades from $D^{1/2}$ to $D^{2/\nu}$. Under standardized $t_{2.5}$ ($\operatorname{Cov}=I_2$ exact, $2/\nu=0.8$) at $\Delta_N = 8 N^{0.65}$ ($\Delta/\sqrt N \to \infty$), single outliers hijack `T81` (`fail: 0.365 → 0.425`), while Door B's robust log-curve converges (`fail: 0.105 → 0.000`, at the cost of requiring symmetry where `T81` stays unbiased under skewness).

## Files & verification

- Full proofs & synthesis: `explore/e17-l25-and-g25.md`
- Pure-stdlib verification library: `explore/e17_l25_g25.py`
- 6 executable checks: `python3 -m unittest -v tests/test_e17_l25_and_g25.py`
