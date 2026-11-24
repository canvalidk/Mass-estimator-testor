# E15: Walking Through Door A — Laguerre–UMPU Uniqueness of `T45`, the Door A Set Dichotomy, and Proof-Hole Audit

**Branch:** `explore/e15-door-a-uniqueness` (from `origin/main` `62e0390`, 24 Nov 2026; Antigravity, at Can's request).
**Scope:** Exploration only (per `explore/README.md`). Nothing in `definition/`, `problems/`, `rulings.md`, or `met/` is modified.
**Code & tests (pure Python 3 stdlib, zero external dependencies):**
- Verification library: `explore/e15_door_a.py`
- Test suite (6 checks): `python3 -m unittest -v tests/test_e15_door_a_uniqueness.py`

---

## Standard mathematics cited first (`GOAL.md` standing rule 6)

- **Neyman structure, completeness, and UMPU tests:** Lehmann & Romano (2005), *Testing Statistical Hypotheses*, §§4.2–4.3 (Theorem 4.2.1 for 1-parameter exponential families / normal location; Theorem 4.3.1 for similarity via complete sufficient statistics).
- **Completeness of Laguerre polynomials / analytic Fourier uniqueness:** Szegő (1939), *Orthogonal Polynomials*, Theorem 5.7.1 (generalized Laguerre basis $\{L_j^{(\beta)}\}$ in $L^2(\mathbb R_+, x^\beta e^{-x}dx)$ for all $\beta > -1$); elementary analytic Fourier-transform proof given self-contained below.
- **Weak-instrument conditional tests:** Moreira (2003), *Econometrica* 71(4), 1027–1048 (conditional likelihood ratio `GCLR`); Anderson & Rubin (1949) (`AR`); Kleibergen (2002) (`T45` score test).
- **Functional invariance & boundary minima:** Zehna (1966) (`T81`); Rockafellar (1970) (`T74`); Gleser & Hwang (1987) / Dufour (1997) (`T80`).

---

## 0. What "Walking Through Door A" Means for the Equation

Walking through Door A takes `R-2026-11-17-1` + `L13` + `L18` literally: **no starting law $\pi$, no handed-over $\nu_{\rm eff}$, no discount $g$ or $\eta$** (retiring `R-2026-11-17-3`'s demand that the *compatible set itself* widen when readings disagree, and reporting the contradiction strictly in the leftover $K_{\min}$ beside the set).

The entire measured-branch equation collapses to **four lines**:
1. **Additive carried state (`T77 (6)`):** $(K(m),\, C(m))$, zero at exact `0/0`, adds across independent blocks.
2. **Evidence curve & leftover (`T73`):** $K(m) = \sum_i \|(I - P_i(m)) z_i\|^2$ on $[0,\infty]$, with $K_{\min} = \min_{m\in[0,\infty]} K(m)$ carried beside the result.
3. **Point (`T81`):** $\widehat m = \arg\min_{m\in[0,\infty]} K(m)$ on the closed arc $[0,\infty]$ (type-safety convention on flat/tied curves per `R-2026-11-22-2`).
4. **Compatible set at reporting level $1-\alpha$ (`T80`):** exact level-$(1-\alpha)$ set under the declared Gaussian inputs (`Door A1`: `T45` score set $S(m)^2/C(m) \le \chi^2_{1,\,1-\alpha}$, where $S(m) = -\frac12 K'(m)$; or `Door A2`: `M_nz`-`GCLR`).

**The central question:** Under Door A, is the compatible set itself **uniquely forced** (`F2`), or does a non-uniqueness hole survive inside Door A?

---

## 1. Theorem E15-1 (Exact Laguerre–UMPU Uniqueness of `T45` on One Instrument)

Work on one isotropic instrument (fixed noise ratio $s = \sigma_F / \sigma_a$, arbitrary dimension $d \ge 1$, arbitrary series with total noise components $M = \sum_{i=1}^N k_i \ge 1$). At candidate mass $m_0 = s\tan\theta_0$ ($\theta_0 \in (0, \pi/2)$), stack the whitened readings in the candidate's $(u_{\theta_0}, v_{\theta_0})$ frame as $U, V \in \mathbb R^M$. By `T73 (4)`, the entire unnormalised curve $K(\cdot)$ is bijective with the $2 \times 2$ Gram matrix $G = \begin{pmatrix} Q_{UU} & Q_{UV} \\ Q_{UV} & Q_{VV} \end{pmatrix}$ ($Q_{UU} = \|U\|^2$, $Q_{UV} = U^\top V$, $Q_{VV} = \|V\|^2$), and on $\{Q_{UU} > 0\}$ (a.s. for every $M \ge 1$) with the three coordinates
$$q := Q_{UU} = C_\theta(m_0) > 0,\qquad X := \frac{Q_{UV}}{\sqrt{Q_{UU}}} = Z_{\rm T45}(m_0) \in \mathbb R,\qquad W := Q_{VV} - X^2 = K(m_0) - Z_{\rm T45}(m_0)^2 \ge 0.$$

**Statement.** Let $\psi(\theta_0; K) \in [0, 1]$ be **any** measurable (possibly randomized) test of $m_0$ (allowed to depend arbitrarily on $\theta_0 = \arctan(m_0/s)$) that satisfies:
1. **Curve-only (`L17` / `T77 (4)` carried state):** depends on the data only through the unnormalised curve $K(\cdot)$, i.e., through $(q, X, W)$ (no reading count $M$);
2. **Exact size $\alpha \in (0, 1)$ under the declared Gaussian (`L13`, `L18`, `L10`, `T80`):** $E_{m_0, \mu}[\psi(\theta_0; K)] = \alpha$ for every $M \ge 1$ and every true push vector $\mu \in \mathbb R^M$;
3. **Locally unbiased:** $\frac{\partial}{\partial m_1}\big|_{m_1 = m_0} E_{m_1, \mu}[\psi(\theta_0; K)] = 0$ for every $M \ge 1$ and every $\mu \in \mathbb R^M$ (in particular, met by every unbiased test $E_{m_1, \mu}[\psi] \ge \alpha$; alternatively, by sign symmetry $\psi(\theta_0; q, -X, w) = \psi(\theta_0; q, X, w)$).

Then:
- **(a) Laguerre–Gamma cancellation of the leftover $W$:**
  $$E_{X \sim \mathcal N(0,1)}[\psi(\theta_0; q, X, w)] = \alpha \quad\text{and}\quad E_{X \sim \mathcal N(0,1)}[X\,\psi(\theta_0; q, X, w)] = 0$$
  for Lebesgue-almost every $(q, w) \in (0,\infty)^2$ (and at $w = 0$ for a.e. $q > 0$). Thus the leftover $W = K(m_0) - Z_{\rm T45}^2$ has **zero usable rejection information** for any curve-only exact test.
- **(b) Strict UMPU dominance and uniqueness of `T45`:**
  For **every** $M \ge 1$, **every** nonzero push $\mu \in \mathbb R^M \setminus \{0\}$, and **every** alternative mass $m_1 \neq m_0$,
  $$E_{m_1, \mu}[\psi(\theta_0; K)] \le E_{m_1, \mu}[\psi_{\rm T45}(\theta_0; K)],\qquad \psi_{\rm T45} := \mathbf 1\big\{|Z_{\rm T45}(m_0)| > z_{1-\alpha/2}\big\},$$
  with **strict inequality** unless $\psi = \psi_{\rm T45}$ almost everywhere. Hence `T45` is the **unique admissible (and unique UMPU)** curve-only exact locally-unbiased test under Door A.
- **(c) Exact size $= \alpha$ is strictly necessary:**
  Relaxing exact size $= \alpha$ to conservative size $\le \alpha$ breaks (a): rejecting on the narrow band $W \in [10.0, 10.5]$ satisfies $\sup_{k \ge 1} P(\chi^2_k \in [10.0, 10.5]) = 0.04434 < 0.05$ (attained at $k = 11$; verified in `tests/test_e15_door_a_uniqueness.py`). Such finite bands still cannot kill anti-parallel hilltops as $\sigma \to 0$ (where $W \to \infty$).
- **(d) Hilltop-ghost obstruction:**
  Even a curve-only exact-size test that explicitly switches rules on the hilltop side $q < w$,
  $$\psi_{\rm hill}(q, X, w) := \mathbf 1\{|X| > z_{1-\alpha/2}\}\,\mathbf 1\{q \ge w\} + \mathbf 1\{|X| < c_\alpha\}\,\mathbf 1\{q < w\}\qquad (P(|X| < c_\alpha) = \alpha,\; c_{0.05} = 0.06271),$$
  holds exact size $\alpha$, keeps `T81`'s interior point, and rejects the exact hilltop $m_0 = 1$ on `(10, -10)` — **but** is strictly dominated by `T45` at every alternative (`M=2`: `0.0761` vs `0.0700`; `M=5`: `0.0968` vs `0.0873`; `M=10`: `0.0788` vs `0.0713` at $\delta = 0.4$, `12,000` trials) and still accepts almost all of $(0, \infty)$ on `(10, -10)` (e.g. at $m = 0.95$ on `(10,-10) & (9.5,-10.5)`, $q < w$ and $|X| \gg c_\alpha$).

### Proof of Theorem E15-1

**Step 1 (Distribution under null and alternative).** Under true mass $m_1 = s\tan(\theta_0 + \delta)$ ($\delta \in (-\pi/2, \pi/2)$) and push vector $\mu \in \mathbb R^M$,
$$U \sim \mathcal N(\xi, I_M),\qquad V \sim \mathcal N(\eta, I_M),\qquad U \perp V,$$
with $\xi = (\cos\delta)\,\mu$ and $\eta = (\sin\delta)\,\mu$. Given $U \neq 0$, decompose $V$ along $\hat U = U/\|U\|$ and $\hat U^\perp$:
$$X = \hat U^\top V \sim \mathcal N(\hat U^\top \eta,\, 1),\qquad W = \|(I - \hat U\hat U^\top)V\|^2 \sim \chi^2_{M-1}(\|\eta_\perp\|^2)\quad (W \equiv 0 \text{ if } M = 1),$$
and **$X \perp W$ given $U$ for EVERY $\delta$ and EVERY $\mu$** (since $\operatorname{Cov}(V) = I_M$). At the null $\delta = 0$: $q \sim \chi^2_M(\Delta)$ ($\Delta = \|\mu\|^2$), $X \sim \mathcal N(0, 1)$, and $W \sim \chi^2_{M-1}$ are **mutually independent**.

**Step 2 (Self-contained analytic Fourier completeness lemma).** Let $f : (0, \infty) \to \mathbb R$ be bounded measurable and $\beta > -1$. If $\int_0^\infty f(x)\, x^{\beta + j} e^{-x/2}\,dx = 0$ for every integer $j \ge 0$, then $f(x) = 0$ for Lebesgue-almost every $x > 0$.
*Proof:* Put $\Phi(t) := \int_0^\infty f(x)\, x^\beta e^{-x/2} e^{i t x}\,dx$. Since $|f| \le C$ and $\beta > -1$, $\Phi(t)$ is holomorphic on the strip $|\operatorname{Im} t| < 1/2$, with derivatives $\Phi^{(j)}(0) = i^j \int_0^\infty f(x)\, x^{\beta+j} e^{-x/2}\,dx = 0$ for all $j \ge 0$. Hence $\Phi(t) \equiv 0$ for all $t \in \mathbb R$. Because $x \mapsto f(x) x^\beta e^{-x/2} \in L^1(0, \infty)$, uniqueness of the $L^1$ Fourier transform gives $f(x) = 0$ a.e. on $(0, \infty)$. *(Valid for all $M \ge 1$, including $M = 1$ where $\beta = -1/2 > -1$.)*

**Step 3 (Cancellation of $q$ and $W$ for size).** Fix $M \ge 1$ and put $h_M(q) := E_{X \sim \mathcal N(0,1),\, W \sim \chi^2_{M-1}}[\psi(\theta_0; q, X, W)] - \alpha \in [-\alpha, 1-\alpha]$. By condition 2 and mutual independence under $\delta = 0$,
$$0 = E_{q \sim \chi^2_M(\Delta)}[h_M(q)] = e^{-\Delta/2} \sum_{j=0}^\infty \frac{(\Delta/2)^j}{j!} \int_0^\infty h_M(q)\, g_{M+2j}(q)\,dq \qquad \forall \Delta \ge 0,$$
where $g_\nu(q) = \frac{1}{2^{\nu/2}\Gamma(\nu/2)} q^{\nu/2 - 1} e^{-q/2}$. Because $|h_M| \le 1$, the power series converges absolutely on all of $\mathbb C$ and vanishes on $[0, \infty)$, so every coefficient vanishes: $\int_0^\infty h_M(q)\, q^{M/2 - 1 + j} e^{-q/2}\,dq = 0$ for all $j \ge 0$. Applying Step 2 ($\beta = M/2 - 1 > -1$) yields $h_M(q) = 0$ for a.e. $q > 0$, for every $M \ge 1$.

Let $N_0 \subset (0, \infty)$ be the countable union of null sets where $h_M(q) \neq 0$ for some $M \ge 1$ ($\operatorname{Leb}(N_0) = 0$). For every $q \notin N_0$, set $f_q(w) := E_{X \sim \mathcal N(0,1)}[\psi(\theta_0; q, X, w)] - \alpha$. From $M = 1$ ($W \equiv 0$), $f_q(0) = 0$. From every odd $M = 2m + 1$ ($m = 1, 2, \dots$),
$$0 = h_{2m+1}(q) = E_{W \sim \chi^2_{2m}}[f_q(W)] = \frac{1}{2^m (m-1)!} \int_0^\infty f_q(w)\, w^{m-1} e^{-w/2}\,dw \qquad \forall m \ge 1.$$
Applying Step 2 ($\beta = 0$) gives $f_q(w) = 0$ for a.e. $w > 0$, so $E_X[\psi(\theta_0; q, X, w)] = \alpha$ for a.e. $(q, w)$ (and at $w = 0$ for a.e. $q > 0$).

**Step 4 (Cancellation of $q$ and $W$ for local unbiasedness).** Differentiate $E_{\delta, \mu}[\psi(\theta_0; q, X, W)]$ once at $\delta = 0$: since $W \sim \chi^2_{M-1}(\delta^2\|\mu_\perp\|^2)$ has noncentrality $O(\delta^2)$ (and $E_X[\psi(q, X, w)] \equiv \alpha$ is constant in $w$ anyway!), only the $X$-shift $\eta_\parallel = \delta\,(\hat U^\top \mu)$ enters at first order:
$$0 = \frac{\partial}{\partial\delta}\Big|_{\delta=0} E_{\delta, \mu}[\psi] = E_{U \sim \mathcal N(\mu, I_M)}\!\big[(\hat U^\top \mu)\, H_M(\|U\|^2)\big],\qquad H_M(q) := E_{X, W}[X\,\psi(\theta_0; q, X, W)].$$
Differentiating $E_{U \sim \mathcal N(\mu, I_M)}[r(\|U\|^2)]$ with respect to $\mu$ (where $r(q) := H_M(q)/\sqrt{q}$) gives the exact Gaussian identity
$$E_\mu\big[(\hat U^\top \mu)\, H_M(\|U\|^2)\big] = 2\Delta\, e^{-\Delta/2} \frac{\partial}{\partial \Delta}\!\left(e^{+\Delta/2} E_{q \sim \chi^2_M(\Delta)}\!\left[\frac{H_M(q)}{\sqrt{q}}\right]\right) = 0 \qquad \forall \Delta > 0.$$
Expanding $e^{+\Delta/2} E_\Delta[H_M(q)/\sqrt{q}] = \sum_{j=0}^\infty \frac{(\Delta/2)^j}{j!} \int_0^\infty \frac{H_M(q)}{\sqrt{q}} g_{M+2j}(q)\,dq$, vanishing of its derivative for all $\Delta > 0$ kills every term $j \ge 1$:
$$\int_0^\infty H_M(q)\, q^{(M+1)/2 - 1 + k} e^{-q/2}\,dq = 0 \qquad \forall k = j - 1 \ge 0.$$
Because $|H_M(q)| \le E|X| = \sqrt{2/\pi} < \infty$ is bounded, Step 2 ($\beta = (M-1)/2 \ge 0$) gives $H_M(q) = 0$ for a.e. $q > 0$, for every $M \ge 1$. Repeating Step 3 over $M = 1$ and $M = 2m + 1$ ($m \ge 1$) yields $E_{X \sim \mathcal N(0,1)}[X\,\psi(\theta_0; q, X, w)] = 0$ for a.e. $(q, w)$ (and at $w = 0$).

**Step 5 (Pointwise Neyman–Pearson / UMPU dominance).** Fix any $M \ge 1$, $\mu \neq 0$, and $m_1 \neq m_0$ ($\delta \neq 0$). Condition on $U$ and $W = w$. Because $X \perp W$ given $U$ even at $\delta \neq 0$, $X \mid (U, W = w) \sim \mathcal N(\theta, 1)$ with $\theta = (\sin\delta)\,\hat U^\top \mu$. By Steps 3–4, $\psi(\theta_0; q, \cdot, w)$ satisfies $E_0[\psi] = \alpha$ and $E_0[X\,\psi] = 0$ for a.e. $(q, w)$. By the UMPU theorem for the 1-parameter normal location family (Lehmann & Romano 2005, Theorem 4.2.1),
$$E_{X \sim \mathcal N(\theta, 1)}[\psi(\theta_0; q, X, w)] \le E_{X \sim \mathcal N(\theta, 1)}\big[\mathbf 1\{|X| > z_{1-\alpha/2}\}\big] \qquad \forall \theta \in \mathbb R,$$
with strict inequality for every $\theta \neq 0$ whenever $\psi(\theta_0; q, \cdot, w) \neq \psi_{\rm T45}$ on a set of positive measure. Since $P(\hat U^\top \mu \neq 0) = 1$ for every $\mu \neq 0$ and $(q, W)$ has strictly positive density on $(0,\infty)^2$ for every $M \ge 2$ (and $q$ on $(0,\infty)$ at $M = 1$), integrating over $(U, W)$ yields strict dominance $E_{m_1, \mu}[\psi] < E_{m_1, \mu}[\psi_{\rm T45}]$ unless $\psi = \psi_{\rm T45}$ almost everywhere. $\square$

---

## 2. Theorem E15-2 (Extension to the Full Reading Type & Multiple Instruments)

1. **Local minimal additive state $(K(m_0), S(m_0), C(m_0))$ on the full reading type (`T77 (6)`):**
   On the full type (arbitrary $d$, arbitrary correlated $\Sigma_i \succ 0$, arbitrary subspaces $N_i$, mixed instruments), `T82 (3)` / `E14-7` already establishes that at true mass $m_0$, given all $U_i$, $X = S(m_0)/\sqrt{C(m_0)} \sim \mathcal N(0, 1)$ and $W = K(m_0) - X^2 \sim \chi^2_{M-1}$ are independent ($M = \sum k_i \ge 1$). Because every triple $(C(m_0), X, W) \in (0,\infty) \times \mathbb R \times [0,\infty)$ is already realized on a single isotropic instrument (a subset of the full type), applying Theorem E15-1 on that slice immediately forces
   $$\psi\big(m_0;\; K(m_0), S(m_0), C(m_0)\big) = \mathbf 1\big\{|Z_{\rm T45}(m_0)| > z_{1-\alpha/2}\big\}$$
   **uniquely across the entire full reading type.**
2. **Global curve-pair state $(K(\cdot), C(\cdot))$ across $J$ distinct isotropic instruments ($s_1 < \dots < s_J$):**
   The global pair $(K(\cdot), C(\cdot))$ is bijective with the $J$ independent triples $\big((q_1, X_1, W_1), \dots, (q_J, X_J, W_J)\big)$. Because exact size $\alpha$ holds for every count tuple $(M_1, \dots, M_J) \in \mathbb Z_{\ge 1}^J$ and every push energy tuple $(\Delta_1, \dots, \Delta_J) \in [0,\infty)^J$, multivariate analytic Fourier completeness kills all $(W_1, \dots, W_J)$ simultaneously:
   $$E_{X_1,\dots,X_J \stackrel{\text{iid}}{\sim} \mathcal N(0,1)}[\psi(q, X, w)] = \alpha \qquad \text{a.e.}$$
   Under a local mass shift $m_1 = m_0 + \delta$, at strong signal ($\mu_j \approx U_j$) the mean vector in $\mathbb R^J$ is $\delta\,(\sqrt{C_1(m_0)}, \dots, \sqrt{C_J(m_0)})$, whose unique 1-dimensional sufficient statistic is $Z_{\rm T45}(m_0) = \frac{\sum \sqrt{C_j} X_j}{\sqrt{\sum C_j}} = \frac{\sum S_j}{\sqrt{\sum C_j}}$. Equivalently, by `T82 (1)` / `E14-5`, among linear combinations $\sum a_j(q) X_j$, only $a_j \propto \sqrt{C_j}$ keeps `T81`'s point $\arg\min K(m)$ inside the set. Hence `T45` is uniquely **strong-signal LMPU** (and unique among `T81`-consistent admissible linear tests) even on the global multi-instrument curve pair $(K(\cdot), C(\cdot))$. *(At weak unknown push allocations across distinct instruments, trading off instrument 1 vs instrument 2 reproduces `T76`'s Pareto frontier unless `T81` consistency or `L17` is imposed.)*

---

## 3. Theorem E15-3 (The Door A Set Dichotomy & `T81`'s Boundary-End Inclusion)

Even inside strict Door A (believed Gaussian $\Sigma_i$, zero discount, zero $\pi$, exact coverage $1-\alpha$), four candidate exact sets exist in the repo. Testing them against `L1`, `L20` (`Z`), continuity at `0/0`, `L17`, and `T81` kills two outright and leaves one razor-sharp dichotomy:

1. **`AR` ($K(m) \le \chi^2_{M,\,1-\alpha}$): DEAD.** Fails `L1` (empty whenever $K_{\min} > \chi^2_{M,\,1-\alpha}$, e.g. `(5, 6)` at $\sigma=0.1$), fails `Z` through $M$, loses power as $N$ grows.
2. **`Z-GCLR` (`T75 (3)`, $w(0) = 0$): DEAD under `T81`.** Extends `T82 (1)` from `T76` scores to `T75 (3)`'s weighted conditional LR: on `(20, 2)` beside `100` readings `(0.2, 2)` at unit noise (`T81` point $m^* = 1.000$ exact), `Z-GCLR` with $w(H) = H/(1+H)$ moves the weighted minimizer to $m = 2.142$ ($LR_{\rm obs} = 20.749$) and **rejects `T81`'s own point $m^* = 1.0$ with $p = 0.0007$** (verified in `tests/test_e15_door_a_uniqueness.py`). Also unpinned continuum of weights, breaks `L17`, and still has the right-angle pinhole.
3. **The Surviving Door A Dichotomy: `T45` (Door A1) vs. `M_nz`-`GCLR` (Door A2):**

| Property | **Door A1: `T45` (Conditional Score)** | **Door A2: `M_nz`-`GCLR` (Conditional LR, $w\equiv 1$, drop exact `0/0`)** |
| --- | --- | --- |
| Exact size $= 1-\alpha$ everywhere | **passes** | **passes** |
| Exact `Z` (`L20`) | **passes** | **passes** |
| **Continuous at `0/0`** | **passes** (`(5,1)` + 50×`1e-9`: `p = 0.00468` unchanged) | **FAILS** (`(5,1)` at $m_0=1$: `p = 0.0020 → 0.0455` at 10×`1e-9` `→ 0.9980` at 50×`1e-9` `→ 1.0000` at 150×`1e-9`!) |
| **`L17` reading splitting** (`1` reading split into `n` equal parts, identical curve $K$) | **passes** (`p = 0.00468` bit-for-bit at $n=1, 10, 50$) | **FAILS** (`p = 0.0053` at $n=1$ `→ 0.0306` at $n=10$ `→ 0.9847` at $n=50$!) |
| Minimal additive state $(K, C)$ (`T77 (6)`) | **passes** (uniquely UMPU/LMPU, Theorems E15-1 & E15-2; zero MC) | **FAILS** (needs $M_{\rm nz}$ + raw readings + MC) |
| **Hilltop ghosts & right-angle pinholes** | **FAILS** (provably unavoidable by Theorem E15-1) | **passes** (kills every hilltop & pinhole) |
| **Answers `T81` "What is not shown" #3: keeps `T81` boundary minimum (`0` or `∞`)?** | **FAILS on anti-aligned ends!** On `(10,-3) & (10,-3.1)`: `T81` point is $m^*=\infty$ ($K(\infty)=18.61$), yet `T45` **rejects $\infty$ ($p = 1.6\times 10^{-5}$)** and keeps ONLY the hilltop ghost at global maximum $m=0.305$ ($K=218.61, p=0.984$)! *(Passes aligned-zero ends like `(10,0)` where $K'(\pi/2)=0$.)* | **passes ALWAYS** ($T_{\rm LR}(m_{\rm end}) = K(m_{\rm end}) - K_{\min} = 0 \implies p = 1.0$ on `domain="positive"`) |
| Two-body set composition (`L5`/`L14`) | `Q_joint` from $(K, C)$ alone (zero MC, ghost-free, **+25% 2-dof tax**); `Q_pos` breaks on asymmetric `100+10` | `T74 (6)`: **1-dof exact**, ghost-free, needs raw readings + MC |

---

## 4. Complete Audit of `definition/proof.md` (`G1`–`G29`) Under Door A

Walking through Door A on `origin/main` (`62e0390`) transforms the proof ledger as follows:

- **Existence (`proof.md` line 92):** Flips from *"no function is known to meet the written requirements together"* to **PROVED** — `door_a1` (`explore/e15_door_a.py`, or `met/ours34.py` with `t_joint` off) already meets every Door-A premise simultaneously on the full reading type.
- **16 of 29 holes dissolve completely:**
  - Starting law: `G8`, `G9` (`O6`, `O11`, `B26`).
  - Discount / tempering / $\eta$ / $\nu_{\rm eff}$: `G3`, `G4`, `G5`, `G13`, `G18`, `G19`, `G20`, `G21`, `G22`, `G23`, `G24` (`O13`, `O14`).
  - Flat-push latent measure: `G2`, `G17`.
  - Right-angle divergence: `G29` (`door_a1` returns flat $K \equiv 400$, set `(0, ∞)`, $K_{\min} = 400$, $p \le 10^{-87}$).
- **5 already closed or provisionally sorted on `main`:** `G7`, `G12`, `G10` (`T71`), `G28` (`T80`), `G6` (`T81`).
- **What remains between Door A and a complete proof (only 4 concrete items):**
  1. **The Door A set fork (replaces `G26` / `O20` finding 5):** **Door A1 (`T45`)** — uniquely forced on $(K, C)$ by Theorems E15-1 & E15-2, continuous at `0/0`, closed-form, zero MC, *at the provable cost of hilltop ghosts and rejecting anti-aligned boundary minima* — **vs. Door A2 (`M_nz`-`GCLR`)** — kills hilltops and keeps boundary minima, *at the provable cost of breaking continuity at `0/0`, breaking `L17` reading splitting, losing $(K, C)$ state, and requiring MC*.
  2. **Two-body set composition (`G27` / `E14-4`):** 2-dof curve-pair `Q_joint` vs. 1-dof raw-reading MC `T74 (6)` (plus clause (c) third-measurement combination).
  3. **`T81` cleanups (`bounds.md`, `G25`):** Can's goal answer on whether `0` and `∞` count as returned masses (`T81 (5)`), stating `T81 (2)`'s concentration axiom `(iv)` explicitly, and writing `L11`'s vector rate (`G25`) directly for `T81`.
  4. **Door-independent foundations:** `G1` (`L13` Gaussianity still stipulated / conditionally Herschel–Maxwell), `G11` ($N_i$ input + VD-Newton rungs R3 unwritten, R5 provisional), `G15` (now reduced to the short Door A chain `T1, T2, T33, T36, T45, T73, T74, T77, T80, T81` + E15), and `G14` (later experimental check on the survivor).
