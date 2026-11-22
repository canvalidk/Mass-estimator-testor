# E14: The Three Walls of `ours34` — Tropical Point Uniqueness, Curve-Pair Composition, and the Disagreement Discount

**Branch:** `explore/e14-ours34-walls` (from `claude/intake-2026-11-22`, 22 Nov 2026; Antigravity, at Can's request).
**Scope:** Exploration only (per `explore/README.md`). Nothing in `definition/`, `problems/`, `rulings.md`, or `met/` is modified.
**Code & tests (pure Python 3.12 stdlib, zero dependencies):**
- Prototype library: `explore/e14_ours34.py`
- Verification suite: `python3 -m unittest -v tests/test_e14_ours34_walls.py`

---

## Standard mathematics cited first (`GOAL.md` standing rule 6)

- **Wall 1 (tropical uniqueness & boundaries):** Infimal convolution (epigraphical addition) and Legendre–Fenchel duality (Rockafellar 1970, *Convex Analysis*, §§5, 16); max-plus / tropical cumulants; compactified projective ray $[0, \infty] \cong [0, \pi/2]$ (`T73`, `T80`).
- **Wall 3 (minimal state & curve-pair composition):** Profile score identities for orthogonal projection families; conditional score sufficiency (Moreira 2003; `T45`; `T77 (4)`); positive-quadrant / boundary likelihood-ratio $\bar\chi^2$ (chi-bar-squared) distributions (Chernoff 1954; Shapiro 1988, *International Statistical Review* 56(1), 49–62).
- **Wall 2 (`T76` collapse, discount no-go, joint-spherical $t$-test):** Conditional sufficiency and scale-invariant $t$-tests (Student 1908; Fisher 1925; Lehmann & Romano 2005, §§4.3, 6.3, 15.2); Kac (1939) and Bernstein (1941) Gaussian independence; Gleser & Hwang (1987) / Dufour (1997) unbounded confidence sets (`T80`).

---

## 1. Wall 1 (`B28`): Tropical Composition Uniqueness & Boundary Obstructions

Let $\mathcal K_{\rm int}$ be the set of finite-series `T73` curves $K:[0,\infty]\to[0,\infty)$ whose global minimum $m^*(K)\in(0,\infty)$ is unique and interior, and let $\overline{\mathcal K}_{\rm int}$ be the closure of $\mathcal K_{\rm int}$ under `T74`'s infimal convolutions $\square_+$ (`L5`), $\square_\parallel$ (`L14`), and `L3` scaling $K_\lambda(m) = K(m/\lambda)$.

### Theorem E14-1 (Tropical uniqueness of $\arg\min_{m>0} K(m)$ on $\overline{\mathcal K}_{\rm int}$)

Suppose a point readout $R:\overline{\mathcal K}_{\rm int}\to(0,\infty)$ satisfies:
1. **`L5` (lumping under `T74`):** $R(K_1 \square_+ K_2) = R(K_1) + R(K_2)$;
2. **`L14` (common force under `T74`):** $1/R(K_1 \square_\parallel K_2) = 1/R(K_1) + 1/R(K_2)$;
3. **`L3` (units covariance):** $R(K_\lambda) = \lambda\,R(K)$ for all $\lambda > 0$;
4. **Exact-branch concentration (`L7`/`L11` limit):** if a sequence $Q_n \in \overline{\mathcal K}_{\rm int}$ satisfies $Q_n(m_0) = Q_{n,\min}$ and $\inf_{m\in[0,\infty]\setminus(m_0-\epsilon,m_0+\epsilon)} (Q_n(m) - Q_{n,\min}) \to +\infty$ for every $\epsilon > 0$, then $R(Q_n) \to m_0$.

Then **uniquely**, for every $K \in \overline{\mathcal K}_{\rm int}$,
$$R(K) = \arg\min_{m>0} K(m).$$
*(Note: neither `T5` swap nor `L20` constant-shift invariance is even needed — both follow as corollaries. And we do NOT use the circular reading "$K(m)=0\iff m=m_0 \implies R(K)=m_0$" on noisy zero-residual curves.)*

**Proof.** Fix $K \in \overline{\mathcal K}_{\rm int}$ with unique global minimum $m^* \in (0,\infty)$ and $K_{\min} = K(m^*)$. Since $K$ is continuous on the compact interval $[0,\infty]$ with unique minimum $m^*$, for every $\epsilon \in (0, m^*)$,
$$c(\epsilon) := \inf_{m\in[0,\infty]\setminus(m^*-\epsilon,m^*+\epsilon)} \big(K(m) - K_{\min}\big) > 0.$$
Define the $n$-fold scaled additive and parallel infimal convolutions:
$$K_+^{(n)}(m) := K^{\square_+ n}(n m) = \inf_{\frac1n\sum_{j=1}^n m_j = m} \sum_{j=1}^n K(m_j), \qquad
  L_n(m) := K_+^{(n)}(m) - n K_{\min}.$$
By (1) and (3), $R(K_+^{(n)}) = \frac1n R(K^{\square_+ n}) = R(K)$. We bound $L_n(m)$ uniformly on both sides of $m^*$:
- **Left side ($m \in [0, m^* - \epsilon]$):** Let $r = \#\{j : m_j \le m^* - \epsilon/2\}$. Since $m_j \ge 0$,
  $$n(m^*-\epsilon) \ge n m = \sum_{j=1}^n m_j \ge (n-r)(m^* - \epsilon/2) \implies r \ge \frac{\epsilon}{2m^*}\,n.$$
  Each of these $r$ indices pays at least $c(\epsilon/2)$, so
  $$L_n(m) \ge \frac{\epsilon\, c(\epsilon/2)}{2 m^*}\,n \qquad \forall m \in [0, m^* - \epsilon].$$
- **Right side ($m \in [m^* + \epsilon, \infty]$):** At least one $m_j$ satisfies $m_j \ge m^* + \epsilon$, so
  $$L_n(m) \ge c(\epsilon) > 0 \qquad \forall m \in [m^* + \epsilon, \infty].$$
  *(Crucial note: on bounded `T73` curves, $L_n(m) \le K(n m - (n-1) m^*) - K_{\min} \to K(\infty) - K_{\min} < \infty$ for $m > m^*$, so `L5` alone NEVER drives the right side to $+\infty$!)*

Now apply the $n$-fold scaled parallel infimal convolution to $K_+^{(n)}$:
$$Q_n(m) := \big(K_+^{(n)}\big)^{\square_\parallel n}(m / n) = \inf_{\frac1n\sum_{j=1}^n \frac{1}{m_j} = \frac1m} \sum_{j=1}^n K_+^{(n)}(m_j).$$
By (2) and (3), $R(Q_n) = n\, R\big((K_+^{(n)})^{\square_\parallel n}\big) = n \cdot \frac{1}{n} R(K_+^{(n)}) = R(K)$ for every $n \ge 1$.
Put $M_n(m) := Q_n(m) - n^2 K_{\min} = \inf_{\frac1n\sum 1/m_j = 1/m} \sum_{j=1}^n L_n(m_j)$. Then $M_n(m^*) = 0 = M_{n,\min}$, and for every $\epsilon \in (0, m^*)$:
- **If $m \in [0, m^* - \epsilon]$:** at least one $m_j \le m^* - \epsilon$, so $M_n(m) \ge \frac{\epsilon\,c(\epsilon/2)}{2m^*}\,n$.
- **If $m \in [m^* + \epsilon, \infty]$:** let $s = \#\{j : m_j \ge m^* + \epsilon/2\}$. Since $1/m_j \ge 0$,
  $$\frac{n}{m^* + \epsilon} \ge \frac{n}{m} = \sum_{j=1}^n \frac{1}{m_j} \ge \frac{n-s}{m^* + \epsilon/2} \implies s \ge \frac{\epsilon}{2(m^* + \epsilon)}\,n.$$
  Each of those $s$ indices pays $L_n(m_j) \ge c(\epsilon/2)$, so
  $$M_n(m) \ge \frac{\epsilon\, c(\epsilon/2)}{2(m^* + \epsilon)}\,n.$$

Combining both sides, for every $\epsilon \in (0, m^*)$,
$$\inf_{m\in[0,\infty]\setminus(m^*-\epsilon,m^*+\epsilon)} \big(Q_n(m) - Q_{n,\min}\big) \ge \frac{\epsilon\, c(\epsilon/2)}{2(m^* + \epsilon)}\,n \xrightarrow{n\to\infty} +\infty.$$
Hence $Q_n \in \overline{\mathcal K}_{\rm int}$ concentrates at $m^*$. By (4), $R(K) = \lim_{n\to\infty} R(Q_n) = m^*$. $\square$

**Why both `L5` and `L14` matter on `T73` curves (inversion of the 28 Oct Jensen conflict):**
- On **probability laws** (`O2`, 28 Oct), `L5` (arithmetic mean) and `L14` (harmonic mean) contradict each other on every spread law by Jensen's inequality.
- On **bounded `T73` curves**, `L5` and `L14` **cooperate and complete each other**: `L5` alone pinches only $m < m^*$ to $+\infty$ while leaving $m > m^*$ capped at the right plateau $K(\infty) - K_{\min}$; `L14` alone pinches only $m > m^*$ while leaving $m < m^*$ capped at $K(0) - K_{\min}$. Together they pinch both sides to $+\infty$.
- On smooth coercive curves (e.g., `T74 (2)` Gaussians), Legendre–Fenchel cumulants $c_r(K)$ (Taylor coefficients of $K^*(p) = \sup_m (pm - K(m))$ at $p=0$) add under $\square_+$ and scale as $\lambda^r$ under `L3`, so additivity + `L3` already kills every $r \ge 2$ and leaves $a_1 c_1 = a_1 m^*$.
- **Explicit failure of level-set midpoints** (verified in `tests/test_e14_ours34_walls.py`): on `T74 (2)` Gaussians $(\mu_1, s_1)=(4,1)$ and $(\mu_2, s_2)=(6,2)$, the geometric midpoint of the $\tau=1$ sublevel set gives $9.747$ for the lump against $3.873 + 5.657 = 9.530$ for the sum (fails `L5`); the arithmetic midpoint passes `L5` on Gaussians but fails `T5` swap ($4/15 = 0.2667 \neq 1/4 = 0.25$).

---

### Theorem E14-2 (Boundary obstruction on $(0,\infty)$ vs. exact closure on $[0,\infty]$)

Every valid measured series ($\Sigma_i \succ 0$, `L2`) whose `T73` curve $K$ does **not** have a unique interior minimum on $(0,\infty)$ falls into one of four classes (all constructed and verified in `tests/test_e14_ours34_walls.py`):
1. **Flat curves ($K \equiv c \ge 0$):** exact `0/0` ($K \equiv 0$), **and** non-zero isotropic-scatter series such as `{(1,0), (0,1)}` at $\Sigma = I_2$ ($K \equiv 1$). Every mass is equally compatible (`R-2026-11-22-2`).
2. **One-sided boundary minima ($0$ or $\infty$):** e.g., $F = 10\pm 1, a = 0\pm 1$ ($K(m) = 100/(1+m^2)$, unique minimum at $\infty$, `T80` set $[5.00, \infty)$), and its swap $F = 0\pm 1, a = 10\pm 1$ (unique minimum at $0$, `T80` set $(0, 0.20]$).
3. **Anti-aligned double-boundary ties ($\{0,\infty\}$):** e.g., $F = 1, a = -1$ at $\Sigma = I_2$ ($K(m) = 1 + 2m/(1+m^2)$, maximum at $1$, tied minima at $\{0, \infty\}$).
4. **Multi-instrument interior ties ($\{m_1^*, m_2^*\}$):** impossible on one instrument (by `T73 (4)`, $K(\theta)$ is a pure $\cos(2\theta)$ sinusoid on $[0,\pi/2]$ with at most one interior minimum), but **achievable on two instruments with $s_1 \neq s_2$**: e.g., $r_1 = (0.1, 1)$ at $(\sigma_F,\sigma_a)=(0.1,1)$ and its swap $r_2 = (1, 0.1)$ at $(1, 0.1)$ give $K(0.1) = K(10) = 0.980 < K(1) = 1.604 < K(0) = K(\infty) = 2.000$.

**Obstruction theorem (no finite extension to $(0,\infty)$):**
Let $R_{\rm ext} : \mathcal K_{\rm all} \to (0,\infty)$ be *any* extension agreeing with $\arg\min K$ on $\mathcal K_{\rm int}$. Then:
1. **$R_{\rm ext}$ CANNOT preserve `L14` on $[L,\infty)$ curves or `L5` on $(0,U]$ curves** (even without continuity, `L3`, or `T5`!).
   *4-line proof:* Let $K_{\rm inf}(m) = 100/(1+m^2)$ (unique minimum $0$ at $m=\infty$) and $K_1(m) = (3-m)^2/(1+m^2) \in \mathcal K_{\rm int}$ (unique minimum $0$ at $m_1^* = 3$). For any $m_1, m_2 \in [0,\infty]$ with $1/m_1 + 1/m_2 = 1/M$, $K_1(m_1) + K_{\rm inf}(m_2) \ge 0$ with equality **iff** $m_1 = 3$ and $m_2 = \infty$, i.e., $M = 3$. By compactness of $[0,\infty]^2$, $K_1 \square_\parallel K_{\rm inf} \in \mathcal K_{\rm int}$ has **unique interior global minimum $M^* = 3$** (verified numerically: `K_par(2.9)=0.00099, K_par(3.0)=0.00000, K_par(3.1)=0.00094`). Hence `L14` forces
   $$\frac{1}{3} = \frac{1}{R_{\rm ext}(K_1 \square_\parallel K_{\rm inf})} = \frac{1}{3} + \frac{1}{R_{\rm ext}(K_{\rm inf})} \implies R_{\rm ext}(K_{\rm inf}) = +\infty.$$
   Symmetrically, $K_1 \square_+ K_{\rm zero} \in \mathcal K_{\rm int}$ has unique interior minimum $3$, so `L5` forces $R_{\rm ext}(K_{\rm zero}) = 0$.
2. **$R_{\rm ext}$ CANNOT be continuous at one-sided boundaries:** $F=10, a=\epsilon \to 0^+$ has $K_\epsilon \to K_{\rm inf}$ uniformly on $[0,\infty]$ while $R(K_\epsilon) = 10/\epsilon \to +\infty$.
3. **$R_{\rm ext}$ CANNOT preserve `L5` or `L3` at flat curves ($K_0 \equiv 0$):** $K_0 \square_+ K_0 = K_0$ forces $R_{\rm ext}(K_0) = 2 R_{\rm ext}(K_0) \implies R_{\rm ext}(K_0) \in \{0, \infty\}$.

**Consequence for `ours34`:** The forced tropical point lives naturally on the **closed ray $[0,\infty]$** (equivalently the quarter-turn angle $\theta \in [0,\pi/2]$ already used by `T73` and `T80`), where `L5`, `L14`, and continuity hold cleanly across all unique-minimum curves (interior and one-sided). Only measure-zero ties and flat curves lack a unique minimizer (`R-2026-11-22-2`).

---

## 2. Wall 3: Minimal Sufficient State Past One Instrument & Exact Curve-Pair Composition

### Theorem E14-3 ($K(m)$ alone is insufficient past one instrument; $(K(m), C(m))$ is the minimal additive state)

`T77.md` (line 104) left open whether the unnormalised curve $K(m)$ alone suffices to combine independent blocks past one instrument. **It does not.**

1. **Exact full-type score identity:** For arbitrary dimension $d$, arbitrary correlated $\Sigma_i \succ 0$, and arbitrary subspace $N_i$ (`T36`), differentiating $K_i(m) = \|z_i\|^2 - z_i^\top P_i(m) z_i$ gives $z_i^\top P_i'(m) z_i = 2\,V_i^\top (P_i'(m) U_i) = 2 S_i(m)$, so
   $$S(m) = \sum_i S_i(m) = -\frac12 K'(m)$$
   holds identically on the full reading type (verified to relative error $1.16 \times 10^{-9}$ on a random correlated $6\times 6$ covariance with a 2D subspace in $\mathbb R^3$).
2. **Exact 4-line counterexample to $K$-only sufficiency past one instrument:**
   - **Series A** (instrument $s_1 = 1$, $\Sigma_1 = \mathrm{diag}(1,1)$): 20 copies of `{(1,0), (0,1)}`.
     $K_A(m) \equiv 20$, and $C_A^{(m)}(m) = \frac{20}{(1+m^2)^2} \implies C_A^{(m)}(1) = 5.0$.
   - **Series B** (instrument $s_2 = 2$, $\Sigma_2 = \mathrm{diag}(4,1)$): 20 copies of `{(2,0), (0,1)}`.
     $K_B(m) \equiv 20$ (**bit-for-bit identical** to $K_A(m)$ for every $m \in [0,\infty]$!), yet $C_B^{(m)}(m) = \frac{80}{(4+m^2)^2} \implies C_B^{(m)}(1) = 3.2 \neq 5.0$.
   - Append the **same** new reading $r_3 = (5, 1)$ at $\Sigma = I_2$ and test $m = 1$ under `T45`:
     - Series A + $r_3$: $T_{\rm T45} = 3.7895 \le 3.8415$ (**KEPT** at 95%).
     - Series B + $r_3$: $T_{\rm T45} = 4.6753 > 3.8415$ (**REJECTED** at 95%).
3. **Minimal additive state for exact Gaussian score combination:** The curve pair **$(K(m),\, C(m))$**, which adds componentwise across independent blocks ($(K_1+K_2,\, C_1+C_2)$), carries exact `0/0` identity $(0, 0)$ (`Z` / `L20`), determines $S(m) = -\frac12 K'(m)$, and reduces to $K(m)$ alone on one instrument (`T77 (4)`).

---

### Theorem E14-4 (Exact curve-pair two-body composition `Q_pos`, zero Monte Carlo)

`T74 (5)` showed that cutting $K_{\rm lump}(M) - K_{1,\min} - K_{2,\min}$ at $\chi^2_1 = 3.8415$ under-covers at weak signal (`0.840`), and `T74 (6)` repaired coverage using raw readings and Monte Carlo. **Neither raw readings nor Monte Carlo is needed:** $(K_1, C_1)$ and $(K_2, C_2)$ alone suffice for an exact finite-sample valid composite test!

Fix a total $M > 0$ and candidate split $(m_1, m_2)$ with $m_2 = M - m_1$. Put $a_j = \sqrt{C_j(m_j)}$ and $Z_j = S_j(m_j)/a_j$ ($j = 1, 2$). At the true split, given $(U_1, U_2)$, $(Z_1, Z_2) \stackrel{\text{i.i.d.}}{\sim} \mathcal N(0,1)$ (with $Z_j = 0$ if $C_j = 0$) independently of $(a_1, a_2)$, and under small mass shifts $(\delta_1, \delta_2)$ their means shift by $(a_1\delta_1, a_2\delta_2)$. Define the positive-quadrant profiled quadratic score:
$$Q_{\rm pos}(m_1; M) := \min_{-m_1 \le \delta \le m_2} \Big[(Z_1 - a_1\delta)^2 + (Z_2 + a_2\delta)^2\Big]
  - (-a_1 m_1 - Z_1)_+^2 - (-a_2 m_2 - Z_2)_+^2.$$
Because $(a_1, a_2, m_1, m_2)$ are fixed by conditioning, the exact null survival probability $p(m_1; M) = P(Q_{\rm pos}(Z_1, Z_2) \ge q_{\rm obs} \mid a_1, a_2)$ is a deterministic 1D kink-split Gauss–Legendre integral (analytic in $Z_2$, exact to $< 10^{-6}$). Taking $p_{\rm lump}(M) = \sup_{0 < m_1 < M} p(m_1; M)$ gives a valid level-$(1-\alpha)$ set for $M = m_1 + m_2$ (and via exact `T5` swap, for $1/M = 1/m_1 + 1/m_2$ under `L14`):

| Check | `Q_pos` (uses ONLY $(K_1,C_1), (K_2,C_2)$, zero MC) | Reference |
| --- | --- | --- |
| Interior exactness ($2.3^2$, strong) | `0.02144822` (err `6.9e-18`) | exact $\chi^2_1$: `0.02144822` |
| Corner $(A=0, B=0)$ at $3.8415$ | `0.0616241` (err `8.9e-7`) | exact Shapiro $\bar\chi^2 = \frac14\delta_0+\frac12\chi^2_1+\frac14\chi^2_2$: `0.0616250` |
| Strong signal $2.3\sigma$ below / above total | `0.02174` / `0.02174` (1-dof!) | `T74 (6)`: `0.022` (1-dof `0.0215`, 2-dof `0.071`) |
| $(5\pm 0.5) + (0/0)$ at $M = 4.17 / 4.18$ | `0.0485` (out) / `0.0505` (in) | exact Chernoff $1-\Phi((5-M)/0.5)$: `4.1776` |
| `T74 (5)` 3D weak signal (`sig = 0, 1, 3, 6`) | true-split coverage `0.970, 0.950, 0.960, 0.943` ($\sup \ge$ this) | raw $K_{\rm lump}$: `0.845, 0.875, 0.910` (fails!) |
| Naive compressed $(K_{\rm lump}, C_{\rm lump})$ | **`0.620, 0.785, 0.945`** at `sig = 0, 1, 3` | **fails at weak signal!** |

**Why naive compression fails (exact state boundary for Wall 3):** Compressing two bodies into a single curve pair $(K_{\rm lump}(M),\, (1/C_1(m_1^*) + 1/C_2(M-m_1^*))^{-1})$ evaluated at the data-dependent split $m_1^*(M)$ collapses coverage to `0.620` at pure noise and `0.785` at signal 1. Exact finite-sample composition requires retaining the parts' curve pairs $(K_1, C_1)$ and $(K_2, C_2)$ (and extends identically to combining a composite with a direct lump measurement $(K_3, C_3)$ for `R-2026-10-28-1` clause (c)).

---

## 3. Wall 2: `T76` Collapse, Curve-Only Discount No-Go, and Joint-Spherical Score Floor `t_joint`

### Theorem E14-5 (Collapse of `T76`'s weighted-score continuum to $w \equiv 1$)

`T76` left an entire continuum of exact conditional scores $T_w(m) = \sum w(H_i) S_i / \sqrt{\sum w(H_i)^2 C_i}$ unpinned by exactness alone. Inside `T76`, **the unweighted score $w \equiv 1$ (`T45`) is uniquely selected** by any one of three premises already on the project's list:
1. **Point-inside-set consistency with $\arg\min K(m)$ (Wall 1):** $T_w(m^*) = 0$ on every interior-minimum series $\iff w(H) \equiv \text{const}$.
   *Explicit counterexample (verified in `tests/test_e14_ours34_walls.py`):* 1 reading $(20, 2)$ plus 100 readings $(0.2, 2)$ at $\Sigma = I_2$ has unique minimum $m^* = \arg\min K(m) = 1.0000$.
   - $w \equiv 1$ (`T45`): $Z(m^*) = 0.000$, $p = 1.000$ (**KEPT**).
   - $w(H) = H/(1+H)$: $Z(m^*) = +3.003$, $p = 0.00267$ (**REJECTS $\arg\min K$ at 99%!**).
   - $w(H) = H$: $Z(m^*) = +12.600$, $p < 10^{-35}$ (**REJECTS $\arg\min K$!**).
2. **`L17` / curve sufficiency on one instrument (`T77 (4)`):** Splitting one reading $(\sqrt n F, \sqrt n a)$ into $n$ identical copies $(F, a)$ leaves $G$ and $K(m)$ bit-for-bit invariant, yet changes $T_w$ for every non-constant $w$.
3. **Curve-pair composability (Wall 3):** $S_w(m) = -\frac12 K'(m)$ holds $\iff w \equiv \text{const}$.

---

### Theorem E14-6 (The Curve-Only Discount No-Go Theorem)

**Statement.** No set rule $S = S(K, C)$ that depends only on $(K(m), C(m))$ (equivalently, on one instrument, only on the $2\times 2$ Gram matrix $G$) can simultaneously satisfy:
1. **`T30` / `R-2026-11-17-3` (bounded confidence on conflicting readings):** on fixed readings $X_1 = m_0 - \Delta/2$, $X_2 = m_0 + \Delta/2$, $\operatorname{width}(S) \ge w_0 > 0$ as $\sigma \to 0$; and
2. **Set-width consistency (`L11` for sets):** on $N$ honest Gaussian readings $X_i \stackrel{\text{i.i.d.}}{\sim} \mathcal N(m_0, \sigma_B^2)$ at fixed $\sigma_B = \Delta/2$, $\operatorname{width}(S) \to 0$ as $N \to \infty$.

**Proof.** Fix $m_0 = 5.5$, $\Delta = 1.0$, and any $N \ge 2$ (in the exact-acceleration limit $\sigma_a \to 0$).
- **Series A (2 conflicting readings at shrinking $\sigma_A$):** $X_1 = 5.0, X_2 = 6.0$ at $\sigma_A = \Delta / \sqrt{2(N-1)}$.
  $$K_A(m) = \frac{4(N-1)}{\Delta^2}(m - m_0)^2 + (N - 1), \qquad C_A^{(m)}(m) = \frac{4(N-1)}{\Delta^2}.$$
- **Series B ($N$ honest agreeing readings at fixed $\sigma_B = \frac{\Delta}{2}\sqrt{\frac{N}{N-1}} \approx 0.5$):** Any $N$-sample realization with sample mean $\bar X = m_0$ and exact sample variance $\frac{1}{N-1}\sum (X_i - m_0)^2 = \sigma_B^2$ (zero surprise: $K_{B,\min} = N - 1 = E[\chi^2_{N-1}]$!) has
  $$K_B(m) = \frac{N}{\sigma_B^2}(m - m_0)^2 + (N - 1) = K_A(m), \qquad C_B^{(m)}(m) = \frac{N}{\sigma_B^2} = C_A^{(m)}(m)$$
  **bit-for-bit identical for every $m$** (verified at $N = 50$ to $10^{-6}$ in `tests/test_e14_ours34_walls.py`). Any curve-only rule must return the exact same width on both: if it stays $\ge w_0 > 0$ on Series A as $N \to \infty$ ($\sigma_A \to 0$), it stays $\ge w_0 > 0$ on $N$ honest agreeing readings at $\sigma_B = 0.5$ and breaks `L11`. $\square$

**Corollary (the exact price of a disagreement discount):** To tell Series A ($M = 2, K_{\min} = N-1$) from Series B ($M = N, K_{\min} = N-1$), any widening rule must carry an effective reading count beyond $(K, C)$:
- Raw count $M = \sum k_i$: breaks `Z` (`L20`) on exact `0/0` readings.
- Non-zero count $M_{\rm nz} = \sum_{y_i \neq 0} k_i$: preserves **exact `Z`** and **exact additive state $(K, C, M_{\rm nz})$**, at the price of a discontinuity at $y_i = 0$ (the exact same discontinuity `T75 (4)`'s spherical null already has).
- Smooth zero-safe count $\sum h_b(\|z_i\|^2)$ (`T78`): continuous at $y_i = 0$, at the price of breaking `L17` and opening an uncountable continuum in $h_b$.

---

### Theorem E14-7 (Exact Closed-Form Joint-Spherical Score Floor `T45` $\cup$ `t_joint`, and the 2-Scalar-Reading Theorem)

On the full reading type (arbitrary $d$, correlated $\Sigma_i \succ 0$, subspaces $N_i$, mixed instruments), suppose the whitened errors share any unknown common scale $v > 0$: $z_i - \mu_i \sim \mathcal N(0, v I_{2k_i})$. At true mass $m_0$, given all $U_i$, the stacked residual vector $V = (V_1;\dots;V_N) \in \mathbb R^{M_{\rm nz}}$ satisfies $V \mid (U, v) \sim \mathcal N(0, v I_{M_{\rm nz}})$, with $\|V\|^2 = K(m_0)$, tangent $t = (P_i' U_i)$ fixed by $U$, $\|t\|^2 = C(m_0)$, and $V^\top t = S(m_0) = -\frac12 K'(m_0)$. Decomposing $V$ along $\hat t$ and its $(M_{\rm nz}-1)$-dimensional complement gives $V_\parallel = S/\sqrt{C} \sim \mathcal N(0, v)$ independently of $\|V_\perp\|^2 = K - S^2/C \sim v\,\chi^2_{M_{\rm nz}-1}$. Hence, whenever $C(m_0) > 0$ and $M_{\rm nz} \ge 2$,
$$t_{\rm joint}(m_0) := \sqrt{M_{\rm nz} - 1}\,\frac{S(m_0)}{\sqrt{K(m_0)\,C(m_0) - S(m_0)^2}} = \frac{Z_{\rm T45}(m_0)}{\sqrt{\hat v(m_0)}} \sim t_{M_{\rm nz}-1} \quad \text{exact},$$
where $\hat v(m_0) = \big(K(m_0) - Z_{\rm T45}(m_0)^2\big)/(M_{\rm nz}-1)$.

**Algebraic equivalence to `O14`'s "expected over observed leftover, capped at 1, first power":**
Joining `T45` and `t_joint` (`floor_score`: keep $m_0$ when either keeps it) is algebraically identical to
$$Z_{\rm T45}(m_0)^2 \le \max\!\Big(z_{1-\alpha/2}^2,\;\; t_{M_{\rm nz}-1,\,1-\alpha/2}^2 \cdot \frac{K(m_0) - Z_{\rm T45}(m_0)^2}{M_{\rm nz}-1}\Big),$$
i.e., discounting $Z_{\rm T45}^2$ by $\min\big(1,\; \frac{z^2}{t_{M_{\rm nz}-1}^2}\frac{M_{\rm nz}-1}{K - Z_{\rm T45}^2}\big)$. Requiring exact $t_{M_{\rm nz}-1}$ validity under unknown common scale $v > 0$ **uniquely pins the exponent to 1** (closing `T78 (5)`'s free exponent $\alpha$ under this premise) and uses **only the additive triple $(K(m), C(m), M_{\rm nz})$**, with **zero Monte Carlo**.

**Measured behavior (verified in `tests/test_e14_ours34_walls.py` and scratch checks):**

| Scene / World | `T45` (Gaussian) | `t_joint` (common scale $v$) | `floor` (`T45` $\cup$ `t_joint`) | Compare `T75 (6)` |
| --- | --- | --- | --- | --- |
| 3D Near `(5.5, 5.6)`, $\sigma = 0.1$ ($M_{\rm nz}=6$) | `[5.412, 5.688]` | `[5.493, 5.607]` | **`[5.412, 5.688]`** (no inflation!) | `5.55 ± 0.14` |
| 3D Far `(5.0, 6.0)`, $\sigma = 0.1$ | `[5.362, 5.638]` | `[4.926, 6.074]` | **`[4.926, 6.074]`** | `5.5 ± 0.45–0.49` |
| 3D Far `(5.0, 6.0)`, $\sigma = 0.01$ | `[5.487, 5.513]` | `[4.926, 6.074]` | **`[4.926, 6.074]`** ($\sigma$-invariant!) | `5.05–5.95` |
| 1D Near `(5.5, 5.6)`, $\sigma = 0.1$ ($M_{\rm nz}=2$) | `[5.416, 5.686]` | `[4.917, 6.181]` | **`[4.917, 6.181]`** (finite!) | **`(0, ∞)`** (vacuous $n \le 5$) |
| 1D Far `(5.0, 6.0)`, $\sigma = 0.1 \to 0.01$ | collapses | `(0, 11.85]` | **`(0, 11.85]`** ($\sigma$-invariant!) | **`(0, ∞)`** |
| `99@5.0 + 1@10.0`, $\sigma = 0.1$ | `[5.031, 5.069]` (misses 5) | `[4.951, 5.149]` | **`[4.951, 5.149]`** (covers 5) | Monte Carlo |
| `G29` right-angle `(2,0), (0,2)` | `(0, ∞)` | `(0, ∞)` | **`(0, ∞)`** (no explosion) | `(0, ∞)` |
| Full correlated 6×6 + mixed subspaces ($v=1$) | `0.952` | `0.944` | `0.980` | `0.97–0.995` |
| Same, common random scale $v \sim \mathrm{Inv}\chi^2(2)$ | `0.810` (fails) | **`0.938`** (exact) | **`0.966`** | `0.943` |
| Per-reading random scales $v_i \sim \mathrm{Inv}\chi^2(2)$ | `0.634` (fails) | **`0.856` (under-covers!)** | **`0.862` (under-covers!)** | **`0.943`** (exact) |
| Honest Gaussian limit ($N = 2, 5, 20, 100, 400$) | `0.95` | `0.95` | **`0.997 → 0.981 → 0.973 → 0.960 → 0.953`** | stays `> 0.97` |

**The 2-Scalar-Reading Theorem (exact classification for Can's scene):**
On $n = 2$ scalar Gaussian readings `(X_1, X_2)` with declared $\sigma$:
1. Under **per-reading unknown scales $(v_1, v_2)$** (`T75 (4)`): no non-trivial 95% equivariant set exists ($2/2^n > 0.05$ for all $n \le 5$); must return $(0,\infty)$ on both near and far pairs.
2. Under **common unknown scale $v > 0$** (`t_joint`): $(\bar X, |X_1 - X_2|)$ is complete sufficient, so the **unique** shortest equivariant 95% interval is Student's $t_1$ interval $\bar X \pm 6.353 |X_1 - X_2|$ (`± 0.635` near, `± 6.353` far — finite, widens on far, $\sigma$-invariant as $\sigma \to 0$).
3. Under **declared Gaussian $v = 1$ only**: any floor $\bar X \pm \max(1.96\sigma/\sqrt{2},\, c |X_1 - X_2|)$ keeps $\ge 95\%$ coverage (e.g., $c = 0.5$ leaves near untouched at `± 0.139` and gives `[5, 6]` on far), **but $c > 0$ is completely unforced** without the common-scale invariance principle.

**Asymptotic reconciliation with `T78 (6)`:**
- At **fixed $N$** (e.g., $N = 2$) with signal $\to \infty$: `floor` over-covers (`0.997` at $N=2$), obeying `T78 (6)`.
- As **$N \to \infty$** on any honest Gaussian identifying series: $\hat v(m_0) = \chi^2_{M-1}/(M-1) \xrightarrow{p} 1$ and $t_{M-1,\,1-\alpha/2} \to z_{1-\alpha/2}$, so `floor` converges to **exact 95.0% coverage** (`0.953` at $N=400$), while still widening persistently whenever readings systematically disagree ($\hat v > 1$).

---

## 4. `ours34` Probe (`met/ours34.py`): Equation, Adversarial Breaks, and `GOAL.md` §4 Report

**Executable probe (`met/ours34.py`, pure stdlib, zero external dependencies; verified by `tests/test_ours34.py`):**
Accepts both `met.full.Reading` (numpy arrays, $d \times k$ orthonormal subspace) and pure-Python `met.ours34.Reading`. Uses a zero-QR $(p, q) = (\sin\phi, \cos\phi)$ Gram-block engine so evaluation is exact and cancellation-free across the entire closed ray $[0, \infty]$.
- **State carried:** additive triple $(K(\phi),\, C(\phi),\, M_{\rm nz})$ (zero prior $\pi$, zero $\nu_{\rm eff}$, zero Monte Carlo).
- **Point:** $\widehat m = \arg\min_{m \in [0,\infty]} K(m)$ on the closed ray $[0,\infty]$ (provisional pick #1; flat/tied curves return declared covariance scale with `point_is_convention=True` per `R-2026-11-22-2`).
- **Compatible set:** exact connected components (all crossing roots on $[0,\pi/2]$) of `floor_score` (`T45` $\cup$ `t_joint`, provisional picks #2 & #3) beside `raw_set` (`T45` alone) and `leftover` ($K_{\min}, M_{\rm nz}, \hat v, \chi^2$ bound).
- **Two-body composition (`lump` for `L5`, `common_force` for `L14`):** `T74` infimal convolution for curve and exact point addition/harmonic addition on $[0,\infty]$; for the composite set, **`Q_joint`** (`method="joint"`, provisional pick #4):
  $$p_{\rm joint}(m_1; M) = \begin{cases} p_1(m_1) & C_2 = 0, \\ p_2(M-m_1) & C_1 = 0, \\ \exp\!\big(-\tfrac12(\tilde Z_1^2 + \tilde Z_2^2)\big) & C_1, C_2 > 0, \end{cases}\qquad \tilde Z_j := \Phi^{-1}\!\big(1 - \tfrac12 p_j^{\rm floor}(m_j)\big).$$
  (`method="qpos"` is also exposed to contrast 1-dof linearized `Q_pos` directly against `Q_joint`.)

### Concrete breaks discovered by throwing adversarial scenes at `ours34` (`tests/test_ours34.py`)

Building and stress-testing `met/ours34.py` immediately surfaced **five concrete breaks** (two of them brand-new discoveries beyond Theorems `E14-1`–`E14-7`):

1. **New Break A — Linearized `Q_pos` fake-split break & discordant collapse (vs. `Q_joint`'s 2-dof tax):**
   - On asymmetric nonlinear bodies (`Body 1 = 100`, `Body 2 = 10`, $\sigma_F = \sigma_a = 0.01$, true total `110`), our own `E14-4` linearized `Q_pos` **breaks catastrophically**: it accepts `M = 30` and `M = 40` with **`p = 1.000`** (kept set `[12.27, 112]`!) because the sub-quadratic curvature of $K_1(m_1)$ makes $S_1/C_1 + S_2/C_2 = 0$ at the fake split $(m_1, m_2) = (13.16, 16.84)$ even though $Z_1 = +500$ and $Z_2 = -500$.
   - On discordant `Body 1 = {5.0, 6.0}` (3D, $\sigma=0.01$, `[4.925, 6.075]`) + clean `Body 2 = 20@10.0` ($\sigma=0.01$), Gaussian `Q_pos` collapses to `[15.49, 15.51]` (30× narrower than Body 1 alone).
   - **`Q_joint` fixes BOTH breaks** with zero MC using only the two bodies' states (`[107.6, 112.5]` on `100+10`, `[14.68, 16.32]` $\sigma$-invariant without dilution on `Disc+Clean`, `[4.02, ∞)` on `5±0.5 + 0/0`), **at the exact price of a 2-dof width tax** at strong signal ($\sqrt{\chi^2_{2,0.95}}/z_{0.975} = 2.448/1.960$, +25% wider). This yields a sharp **3-way composition trade-off**: `Q_pos` (1-dof, curve-pair, breaks on nonlinear fake splits & discordant bodies) vs. `Q_joint` (2-dof, curve-pair, ghost-free & discordant-safe) vs. `T74 (6)` (1-dof, ghost-free, requires raw readings + MC).
2. **New Break B — `M_nz` discontinuity dilution break (`Theorem E14-6` in action):**
   - Appending 20 exact `0/0` readings to the 3D far pair `(5.0, 6.0)` leaves the set bit-for-bit at `[4.925, 6.075]`.
   - Appending 20 readings of tiny size $10^{-9}$ ($10^{-8}\sigma$, indistinguishable from `0/0`) jumps $M_{\rm nz}: 6 \to 66$, diluting $\hat v$ by $13\times$ and **collapsing** the set all the way back to `T45`'s `[5.361, 5.639]`.
   - By contrast, per-reading `T75 (6)` (`SCLR+`) is immune to tiny-$\epsilon$ dilution (since $\|V_i\| = 10^{-9}$ contributes nothing to its spherical draws), at the price of vacuous $(0,\infty)$ on $n \le 5$ scalar readings + MC.
3. **Break C — Stationary hilltop & right-angle pinhole ghosts (score-test limitation):**
   - Anti-parallel `(10, -10)` ($\sigma=0.1$) keeps a razor-thin ghost `[1, 1]` centered on the **maximum** of $K(m)$!
   - Bimodal 2-instrument conflict `(0.1 & 10)` keeps a third ghost island `[0.844, 1.185]` in `raw_set`.
   - One-sided `F=10, a=0` keeps an isolated right-angle pinhole `[0, 0]` alongside `[5.003, ∞)` (which propagates as `[0, 0]` into `common_force` beside `[1.703, 3.245]`).
   - 1D Far `(5.0, 6.0)` at tiny $\sigma_a = 10^{-8}$ keeps a high-mass hilltop ghost `[1.17e10, ∞)` beside `[0, 11.85]`.
4. **Break D — Per-reading heavy tails ($v_i$ varies per reading):**
   - Common-scale `t_joint` covers `>= 0.95` under shared $v \sim \mathrm{Inv}\chi^2(2)$ (`0.966`), but under-covers (`~0.86`, vs. `T75 (6)`'s `0.94`) when each reading draws an independent $v_i \sim \mathrm{Inv}\chi^2(2)$.
5. **Break E — Single-scalar $M_{\rm nz} = 1$:**
   - $df = M_{\rm nz} - 1 = 0$ leaves no residual dof for `t_joint`; requires either `T45` fallback (`single_scalar_policy="t45"`) or vacuous `(0, ∞)` (`"vacuous"`).

### What `ours34` passes where `ours30`–`ours33` broke (`tests/test_ours34.py`)
- **`T39` opposed pair (`ours32` bias & `ours33` >60-decade blowup):** `pt = 2.000`, finite set `[0.798, 34.02]` (log10 span `< 1.7`, zero Cauchy-tail blowup).
- **`G29` right-angle scene (`ours33` diverged to $\infty$):** clean `(0, ∞)`.
- **`R-2026-11-17-3` & `T30`:** 3D Near stays `[5.411, 5.689]` (zero inflation); 3D Far widens to `[4.925, 6.075]` and stays bit-for-bit $\sigma$-invariant across $100\times$ $\sigma$ reduction; `T78 (2)` `{0.2, 1.0, 1.4}` stays `[0.601, 1.664]`; `99@5.0 + 1@10.0` covers `5.0` (`[4.951, 5.149]`).
- **`T80` & `B26` & `O6`/`O11`:** no starting law $\pi$, no handed-over $\nu_{\rm eff}$, exact 4-shape unbounded sets, exact `L5`/`L14` tropical points (including boundary chaining `(F=10, a=0) || (3 ± 0.1) -> pt = 3.0`).

### Honest `GOAL.md` §4 matrix (Standing Rule 1: *passes / fails / unknown*, plus premise named)

| Premise | Status | Details / Premise named on failure |
| --- | --- | --- |
| L1 (one $m>0$, no gate) | **passes** | Never empty on $m \in [0,\infty]$; interior point always inside `T45`/`t_joint` set |
| L2 ($\sigma>0$, never undefined) | **fails on open $(0,\infty)$ / passes on $[0,\infty]$** | **Theorem E14-2:** one-sided sets force point at $0$ or $\infty$; names whether `L2`'s type for the point is open $(0,\infty)$ or closed $[0,\infty]$ (`T80`'s noise-angle arc) |
| L3 (units covariance) | **passes** | Exact on curve, point, `T45`, `t_joint`, `Q_joint`, and `Q_pos` |
| L4 ($\pi$ symmetry) | **vacuous** | No starting law $\pi$ (`O6`, `O11`, `B26` dissolved) |
| L5 & L14 (composition) | **passes (point & `Q_joint` valid coverage); fails 1-dof sharpness (`Q_joint` 2-dof tax)** | Point exact at finite noise (Theorem E14-1); `Q_joint` works for clean & discordant bodies from state alone (+25% 2-dof width tax); linearized 1-dof `Q_pos` fails on nonlinear fake splits & discordant bodies |
| L8, L9, L10, L13, L16 | **passes** | Built directly on `T36` whitened span |
| L11 (consistency at sharp rate) | **passes** | `argmin K` attains `B27` rate; `floor` converges to `T45` exact 95% as $N\to\infty$ |
| L17 (same curve, same mass) | **passes for point & `T45`; fails for `t_joint`** | **Theorem E14-6:** unavoidable for ANY discount keeping `T30` + `L11` (requires $M_{\rm nz}$ beside $(K, C)$) |
| L18 (believed $\sigma$, no handed-over $\nu$) | **passes (with caveat)** | No $\sigma$ revised, no $\nu_{\rm eff}$; **caveat:** `t_joint` imports common-scale $v$ validity on top of Kac–Bernstein (`T75 (2)`) |
| L19 (no borrowing) | **passes (same caveat as `T78`)** | Per-reading split at $m_0$ borrows nothing; pooled $\hat v(m_0)$ aggregates residuals across readings |
| L20 (`Z` on exact `0/0`) | **passes on exact `0/0`; fails continuity at `0/0`** | **Break B:** exact `0/0` changes nothing, but 20 readings of size $10^{-9}$ dilute $\hat v$ and collapse `[4.925, 6.075]` to `[5.361, 5.639]` |
| L21 & T30 (disagreement widening) | **passes for $M_{\rm nz} \ge 2$; fails at $M_{\rm nz} = 1$** | Exact $\sigma$-invariance on 1D ($N\ge 2$) and 2D/3D ($N\ge 1$); at $M_{\rm nz} = 1$, $df = 0$ leaves no dof for $v$ |
| Anti-parallel hilltop & right-angle ghost | **fails (Break C, score-test limitation)** | Anti-parallel `(10,-10)` keeps `[1,1]` at maximum; bimodal keeps ghost `[0.844, 1.185]`; `F=10, a=0` keeps pinhole `[0,0]` |
| Per-reading heavy tails ($v_i$ varies per reading) | **fails (Break D, `~0.86` at $\nu=2$)** | `t_joint` assumes common unknown scale $v$; `T75 (6)` assumes per-reading $v_i$ (`0.94`, at cost of vacuous $n \le 5$ scalar sets + MC) |
| Property N / `T80` (unbounded sets) | **passes** | Valid $\ge 1-\alpha$ for every $m>0$ and every signal; four shapes continuous in $\theta$ |
| L22, L23 | **passes** | Scalar mass, parametrisation-free |

### The four sharp forks left for Can (compressed from the entire E13/E14 space)

1. **Point codomain (Wall 1, Theorem E14-2):** On one-sided sets ($[L,\infty)$, $(0,U]$), `L5` + `L14` prove no finite point in $(0,\infty)$ works. Does the point type expand to the closed noise-angle arc $[0,\infty]$ (matching `T73`/`T80`), or does `R-2026-11-22-2` ("point value doesn't matter") extend to every curve with no interior minimum?
2. **Disagreement widening principle (Wall 2, Theorem E14-6 & E14-7, Breaks B & D):** Under Kac–Bernstein, widening on disagreement requires an explicit extra principle beyond believed Gaussian $\sigma$:
   - **Option A (Door A, `T45` only):** Exact $(K, C)$ state, continuous at `0/0`, `L17` intact, exact 95% everywhere, disagreement strictly in leftover $K_{\min}$ (contradicts 17 Nov reading of `R-2026-11-17-3`).
   - **Option B (`ours34` default: `T45` $\cup$ `t_joint`, common scale $v$):** Additive state $(K, C, M_{\rm nz})$, zero MC, uniquely pins exponent $1$, works on 2 scalar readings (`± 0.63` near, `± 6.35` far) and 2 3D readings (`± 0.138` near, `± 0.574` far), 95.0% limit as $N\to\infty$; costs: Break B ($10^{-9}$ readings dilute $\hat v$), Break D (per-reading $\nu=2$ under-covers `0.86`), Break E ($df=0$ at $M_{\rm nz}=1$).
   - **Option C (`T75 (6)`, per-reading scales $v_i$):** Immune to Break B and Break D (`0.94` at $\nu=2$), kills multi-direction hilltops; costs: vacuous $(0,\infty)$ on $n \le 5$ scalar readings, persistent over-coverage, requires raw readings + MC, `Z-GCLR`'s weight $w(H)$ breaks `L17` and can reject $\arg\min K$.
3. **Hilltop ghosts vs. additive closed-form state (Break C):** Score tests (`T45`, `t_joint`) have exact additive state $(K, C, M_{\rm nz})$ and zero MC, but keep stationary hilltops ($K'=0$) and right-angle pinholes ($U=0$); conditional LR ($M_{\rm nz}$-GCLR) kills both, at the cost of raw readings + MC.
4. **Two-body set composition (Break A):** `Q_joint` (default in `ours34`) uses only additive states, has zero MC, kills fake splits, and handles discordant bodies without dilution — at a +25% 2-dof width tax (`2.448σ` vs `1.960σ`). Recovering 1-dof sharpness without fake splits currently requires `T74 (6)` (raw readings + MC, Gaussian only).

