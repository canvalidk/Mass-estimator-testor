# Learning the push scale in noise units: the fork, a tempering identity, and coverage

**2026-09-27.** Record, filed at the user's request ("Great! Yes do those."),
written by Claude from the checks of that afternoon.

**Relations.** Tests the $\omega$ hypothesis and the units fork of this folder's
`README.md`. Supports the noise-unit reading in
`ours26_starting_law_2026-09-27.md` §6. Bears on the lemma record's L9 and L10
(`../mass_estimator_lemmas_2026-09-27.md`). Adds an option to ours26 and binds
nothing.

**Status.** Standalone simulation and algebra in Claude's sandbox, not yet the
tester. Intervals were checked; the readout under the $\omega$-law was not.

---

## 1. Two ways to give the pushes a finite scale

Per noise pair $k$ and per component, in that pair's noise units, a candidate
push is $v$ at angle $\theta_k=\arctan(m/s_k)$. Give the pushes a Gaussian of
width $w_k$. Integrating $v$ out, each component pair $(x,y)$ is
$N\big(0,\;I+w_k^2\,\mathbf u_k\mathbf u_k^\top\big)$ with
$\mathbf u_k=(\sin\theta_k,\cos\theta_k)$.

| reading of isotropy | width in noise units |
|---|---|
| in noise units | $w_k=\omega$, the same on every axis |
| in physical units (true acceleration $\sim N(0,\omega_a^2I)$) | $w_k=\omega_a/(\sigma_{a,k}\cos\theta_k)$ |

**One noise pair.** At fixed $\theta$ the two are the same family, since
$w=\omega_a/(\sigma_a\cos\theta)$ is just a rescaling of $\omega_a$. Under
scale-free priors ($dw/w=d\omega_a/\omega_a$) they are the same model. The fork
only exists with several noise pairs.

## 2. The tempering identity

For one component pair,
$\lVert\mathbf y\rVert^2=(\mathbf y\cdot\mathbf u)^2+(\mathbf y\cdot\mathbf u_\perp)^2$,
and the flat (ours26) likelihood is
$\exp\!\big[-\tfrac12(\mathbf y\cdot\mathbf u_\perp)^2\big]$. So, with the same
width $\omega$ on every component,

$$
\log L_\omega(\theta)=\text{const}(\omega)+\tfrac12\,c(\omega)\,H(\theta)^2,
\qquad c(\omega)=\frac{\omega^2}{1+\omega^2},
$$

while the flat law's log-likelihood is $\tfrac12H(\theta)^2+\text{const}$. With
several noise pairs, $H^2\to\sum_gH_g^2$. Hence

$$
\boxed{\;p_\omega(\theta)\;\propto\;\pi(\theta)\,L_{\rm flat}(\theta)^{\,c(\omega)}\;}
$$

- **The peak does not move.** The likelihood's maximum in $\theta$ is ours26's
  for every $\omega$. Checked: the argmax agrees to 8 digits at
  $\omega$ = 0.3, 1, 5 and 50.
- **The width grows** by $1/c=1+1/\omega^2$.
- The flat law is the limit $\omega\to\infty$, $c\to1$.

## 3. Learned $\omega$ gives the missing width factor

Maximising jointly over $\theta$ and $\omega$ gives

$$
\hat\omega^2=\frac{\lambda_{\max}}{M}-1,\qquad
\lambda_{\max}=\max_\theta H^2=\frac T2+2\lvert C\rvert,
$$

with $(C,T,M)$ the combination state. At the true angle
$E[H^2]=R^{\ast2}+M$, so $\hat\omega^2\approx R^{\ast2}/M$ and

$$
\frac1{c(\hat\omega)}\approx1+\frac{M}{R^{\ast2}},
$$

which is exactly the factor by which the flat law is too narrow (lemma T18).
It holds for any distribution of pushes, fixed-length ones included, because
$\hat\omega^2$ estimates the mean square push per component.

## 4. The fork, decided by consistency

L11 asks for convergence **whatever the pushes are**. Two axes, axis ratios 5
and 0.8 ($\sigma_F=1,1$; $\sigma_a=0.2,1.25$), 40,000 readings per cell,
maximum-likelihood mass:

| world | true $m$ | flat (ours26) | noise-unit $\omega$ | physical $\omega_a$ |
|---|---|---|---|---|
| pushes isotropic in physical units | 0.5 | 0.4966 | 0.4966 | 0.4969 |
| | 2 | 2.0006 | 2.0006 | 1.9993 |
| | 8 | 8.0024 | 8.0024 | 8.0021 |
| pushes isotropic in noise units | 0.5 | 0.4968 | 0.4968 | **1.0144** |
| | 2 | 2.0045 | 2.0045 | **2.8345** |
| | 8 | 7.9502 | 7.9502 | **8.2599** |

- The physical model converges only when its world is the true one. Its
  per-reading determinant $1+w_k^2$ depends on $m$: evidence when the model is
  right, drag when it is not.
- The noise-unit model always returns ours26's point (§2), so it converges
  wherever ours26 does.
- **So L11 selects isotropy in noise units.** That is also the premise the
  arc-length prior's Jeffreys route needs (equal pushes per axis in noise
  units).

## 5. Coverage

Plug-in $\hat\omega$ (floored at $\omega^2=10^{-3}$), uniform-angle prior, 95%
central interval of the law in $\theta$ (so equally in $m$), 2,000 series per
cell, one noise pair with $\sigma_F=\sigma_a=1$, new push per reading in a
random direction:

| cell | flat | learned $\omega$ | $1/c$ (median) | $1+d/r^{\ast2}$ |
|---|---|---|---|---|
| 3D, $N=21$, SNR 1, $m=s$ | 78.3% | **96.2%** | 2.45 | 2.50 |
| 3D, $N=21$, SNR 1, $m=s/4$ | 77.9% | **95.9%** | 3.55 | 3.82 |
| 1D, $N=21$, SNR 1, $m=s$ | 89.4% | **95.2%** | 1.49 | 1.50 |
| 3D, $N=6$, SNR 1, $m=s$ | 84.3% | 97.8% | 2.32 | 2.50 |
| 3D, $N=21$, SNR 3, $m=s$ | 93.4% | **95.1%** | 1.17 | 1.17 |
| 3D, $N=21$, SNR 1, $m=s$, push sizes spread | 83.2% | **95.5%** | 1.95 | — |

Pushes have fixed length (the tester's kind) except in the last row, where
sizes are half-normal with the same mean. The flat column reproduces the
tester's undercoverage (76.5% in `../the_combination_equation_is_exact_2026-09-26.md`
§4). Standard error about 0.5 points near 95%, 0.9 near 80%.

## 6. The prior on $\omega$

- **Scale-free does not work.** As $\omega\to0$ the likelihood tends to a
  positive constant that does not depend on the mass, so under $d\omega/\omega$
  the posterior is improper and the readings wash out. With five readings, the
  posterior ratio $p(m{=}2)/p(m{=}0.5)$ drifted from 0.709 to 0.833 as the lower
  cut-off went from $10^{-1}$ to $10^{-32}$, towards 1.
- The runs above used a plug-in instead. A proper prior on $\omega$ needs a
  scale; in noise units the noise itself is one, and it has no physical units,
  so L3 is not touched. Which prior: open.
- So L10's "no preferred push size before data" cannot be met literally by a
  learned $\omega$. What is available is "no preferred size apart from the
  noise".

## 7. What survives, what changes

- **Survives:** axis splitting (one $\omega$ shared by a reading and its axes);
  L3 ($\omega$ is a number in noise units); the swap ($H^2$ is symmetric under
  it); the $(C,T,M)$ state, since $\hat\omega$ comes from $T/2+2\lvert C\rvert$
  and $M$.
- **Changes:** L10 read as "$\omega$ unknown and learned", not
  "$\omega=\infty$"; the units fork settled on noise units, pending the tester;
  ours25's sandwich replaced by a width that comes out of the law.
- **New behaviour at weak signal:** when no push is detectable, $\hat\omega\to0$
  and the data are tempered towards the prior. Reasonable, but new in the
  Fieller corner.

## 8. Open

- The prior on $\omega$ (§6), and plug-in against integrating $\omega$ out.
- The readout under the $\omega$-law (its weight $A$ changes: the push given
  $\theta$ shrinks by $c$). Not checked.
- Several noise pairs: one common $\omega$ calibrates exactly only if the
  pushes are equally large in each pair's noise units.
- Slight over-coverage at small $N$ (97.8% at $N=6$).
- The tester: implementation, the presets above, and the axis-splitting check.

## 9. Reproduction

Scripts run in Claude's sandbox on 27 September, not kept: the fork table
(two-axis worlds, joint maximum likelihood over mass and scale by Nelder–Mead
from nine starts); the peak check; the hyperprior drift (quadrature over
$\ln\omega$); the coverage table (angle grid of 20,001 points). The tester
implementation is to follow.
