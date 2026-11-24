# E15: Walking Through Door A — Uniqueness of `T45`, the Set Dichotomy, and Proof-Hole Audit

Branch: `explore/e15-door-a-uniqueness` (from `origin/main` `62e0390`, 24 Nov 2026). Exploration only; nothing in `definition/`, `problems/`, `rulings.md`, or `met/` changed.

**Core target:** Walk all the way through Door A (believed Gaussian $\Sigma_i$, no $\pi$, no $\nu_{\rm eff}$, no set-widening discount; contradiction reported strictly in the leftover $K_{\min}$ beside the set) and test whether the full measured mass definition becomes **provably unique** (`F2`).

## Three theorems + complete `proof.md` audit (for Can & Claude to harvest)

1. **`E15-1` — Exact Laguerre–UMPU uniqueness of `T45` under Door A on one instrument** (`O20`, `G26`, `T45`, `T73`, `T77`, `T80`): any curve-only test $\psi(\theta_0; K) \in [0,1]$ (no reading count $M$, allowed to depend arbitrarily on $\theta_0 = \arctan(m_0/s)$) with **exact size $= \alpha$** for every $M \ge 1$ and every push $\mu$, and **local unbiasedness** ($\partial_{m_1}|_{m_0} E[\psi] = 0$), satisfies $E_X[\psi(q,X,w)]=\alpha$ and $E_X[X\psi(q,X,w)]=0$ a.e. (via noncentral-$\chi^2$ + Laguerre/analytic-Fourier completeness, valid for all $M \ge 1$). Therefore the leftover $W = K(m_0) - Z_{\rm T45}^2$ carries **zero usable rejection information**, and by Lehmann–Romano Theorem 4.2.1 **`T45` ($\mathbf 1\{|Z_{\rm T45}| > z_{1-\alpha/2}\}$) strictly dominates every competitor uniformly at every $(M, \mu, m_1 \neq m_0)$**. Shows exactness $= \alpha$ is necessary (band $W \in [10, 10.5]$ has $\sup_{k\ge 1} P(\chi^2_k \in [10,10.5]) = 0.0443 < 0.05$) and proves hilltop ghosts are unavoidable on $(K, C)$ (hilltop-killing exact competitor $\psi_{\rm hill}$ is strictly dominated everywhere and still keeps almost all of $(0,\infty)$ on `(10,-10)`).
2. **`E15-2` — Full-type and multi-instrument extension** (`T77 (6)`, `T82 (1)`): on the full reading type (arbitrary $d$, correlated $\Sigma_i$, subspaces $N_i$), restricting to one isotropic slice immediately forces `T45` uniquely over all local tests $\psi(m_0; K(m_0), S(m_0), C(m_0))$. On the global multi-instrument curve pair $(K(\cdot), C(\cdot))$, multivariate Laguerre completeness kills all $(W_1,\dots,W_J)$, and strong-signal LMPU (or `T81` point-inside-set + linear admissibility) uniquely selects `T45`.
3. **`E15-3` — The Door A set dichotomy (`T45` vs `M_nz`-`GCLR` vs `Z-GCLR`) & `T81`'s boundary-end inclusion**:
   - `AR` fails `L1` (can be empty) and `Z`.
   - **`Z-GCLR` ($w(H)=H/(1+H)$) REJECTS `T81`'s own point $\arg\min K$** on `(20,2) + 100×(0.2,2)` at $m^*=1$ ($LR_{\rm obs}=20.75$, **`p = 0.0007`**), extending `T82 (1)` from `T76` scores to `T75 (3)` weighted conditional LR.
   - **`M_nz`-`GCLR` kills hilltops and ALWAYS keeps `T81`'s boundary minima (`p = 1.0`),** *at the cost of breaking continuity at `0/0`* (`(5,1)` at $m_0=1$: `p = 0.0020 → 0.0455` with 10×`1e-9` `→ 0.9980` with 50×`1e-9`), *breaking `L17` reading splitting* (1 reading split into 50 equal parts with identical curve $K$: `p = 0.0053 → 0.9847`, vs `T45` bit-for-bit `0.00468`), and losing additive $(K, C)$ state.
   - **Answers `T81` "What is not shown" #3:** on anti-aligned `(10,-3) & (10,-3.1)` (`T81` point $m^* = \infty$, $K(\infty) = 18.61$), `M_nz`-`GCLR` keeps $m^* = \infty$ (`p = 1.0`), whereas **`T45` REJECTS $m^* = \infty$ (`p = 1.6e-5`)** and keeps ONLY the hilltop ghost at the global maximum $m = 0.305$ ($K = 218.61$, `p = 0.984`).
4. **Proof-hole audit (`definition/proof.md` `G1`–`G29`):** Under Door A, **existence is proved** (`door_a1` meets every Door-A premise simultaneously), **16 of 29 holes dissolve** (`G2–G5, G8–G9, G13, G17–G24, G29`), 5 are already closed/sorted (`G6, G7, G10, G12, G28`), leaving only the `T45` vs `M_nz`-`GCLR` set fork, `G27` (2-dof `Q_joint` vs MC `T74 (6)`), `T81`/`G25` cleanups, and the door-independent stipulations (`G1`, `G11`, `G14`).

## Files & verification

- Full proofs, tables, and `proof.md` audit: `explore/e15-door-a-uniqueness.md`
- Pure-stdlib verification library & Door A1 wrapper: `explore/e15_door_a.py`
- All 6 checks: `python3 -m unittest -v tests/test_e15_door_a_uniqueness.py`
