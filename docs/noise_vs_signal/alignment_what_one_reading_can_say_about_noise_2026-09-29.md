# Alignment: what one reading can say about noise

**2026-09-29 (night).** Record, written by Claude at the user's request
("And do number 2"), the second of the day's hints. It corrects a claim
Claude made in conversation the same afternoon (that in $d\ge2$ a
reading's misalignment supplies per-reading noise evidence). It does not.

**What it records.** Per reading, the law sees one complex number $C_i$
and one real number $T_i$. In 1D every reading is exactly aligned with
some mass ($|C_i|=T_i$). In $d\ge2$ a reading can be misaligned, but its
misalignment says almost nothing about whether it is a push or noise (2D:
median alignment 0.87 under noise, 0.88 under a push of one noise unit).
The only per-reading evidence of signal is the reading's size along its
best line, $(T_i+|C_i|)/2$, which is what ours28's per-reading scale reads.
Evidence that a *series* is noise is collective: $|\sum_iC_i|/\sum_iT_i$
sits at $1/\sqrt{dN}$ under noise and rises with signal that no single
reading can show. So the collision between a noise measure and L19 is
exact in every dimension, not only in 1D.

**Relations.** Scopes `null_series_tilt_and_property_N_2026-09-29.md` §3
(which deferred this) and §5, the view record §2, the lemma record (L19,
R1, T25) and `what_forces_the_ratio_of_means_readout.md` §7 (the old
gate's $\chi^2_{d-1}$). Binds nothing.

---

## 1. One reading is two numbers

Whitened reading $(\tilde{\mathbf F}_i,\tilde{\mathbf a}_i)\in\mathbb R^{2d}$,
identity noise. Write

$$
T_i=|\tilde{\mathbf F}_i|^2+|\tilde{\mathbf a}_i|^2,\qquad
C_i=\big(|\tilde{\mathbf a}_i|^2-|\tilde{\mathbf F}_i|^2\big)+2i\,\tilde{\mathbf F}_i\!\cdot\!\tilde{\mathbf a}_i .
$$

Then at $m=\tan\theta$ the explained energy is
$H_i(\theta)=\tfrac12T_i+\tfrac12|C_i|\cos(2\theta-\varphi_i)$, and

$$
\max_\theta H_i=\tfrac12\big(T_i+|C_i|\big),\qquad
\min_\theta Q_i=\tfrac12\big(T_i-|C_i|\big),
$$

the size along the reading's best line and the residual off it. (Checked
against a grid in 2D to $10^{-6}$.) The **alignment fraction**
$|C_i|/T_i\in[0,1]$ is the share of the reading's energy some mass can
explain. The series state is the sum, $C=\sum_iC_i$, $T=\sum_iT_i$ (the
combination record's $(C,T)$), and the law's tilt is $|C|$.

**In 1D, $|C_i|=T_i$ always.** $|C_i|^2=(\tilde a^2-\tilde F^2)^2+4\tilde F^2\tilde a^2=(\tilde a^2+\tilde F^2)^2$.
Every reading is exactly aligned with $m=\tilde F/\tilde a$ (or with a
boundary when the signs differ). One 1D reading has no residual to show.

## 2. In $d\ge2$, misalignment exists and says almost nothing

Median alignment fraction $|C_i|/T_i$ [10th, 90th percentile], 20 000
readings, push of $r$ noise units along a random direction, true $m=1$:

| $d$ | $r=0$ (noise) | $r=1$ | $r=3$ | $r=10$ |
|---|---|---|---|---|
| 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2 | 0.87 [0.44, 0.99] | 0.88 [0.46, 1.00] | 0.96 [0.75, 1.00] | 1.00 [0.97, 1.00] |
| 3 | 0.71 [0.32, 0.95] | 0.72 [0.32, 0.95] | 0.88 [0.62, 0.98] | 0.99 [0.96, 1.00] |

A noise reading in 2D is 87% aligned with *some* mass; a one-unit push is
88%. Misalignment separates strong pushes from noise (which the size
already does) and cannot separate weak pushes from noise, which is the
case that matters.

## 3. The residual is not noise evidence

$\min_\theta Q_i=\tfrac12(T_i-|C_i|)$, the old gate's statistic, in 2D:

| $r$ | mean | 90th pct |
|---|---|---|
| 0 | 0.43 | 1.17 |
| 1 | 0.60 | 1.62 |
| 3 | 0.95 | 2.56 |
| 10 | 1.00 | 2.69 |

Under noise the fit absorbs part of the residual (mean 0.43, below the
$\chi^2_1$ mean of 1); under a strong push the line is pinned and the
residual is the full $\chi^2_{d-1}$. So it is weakly informative, in the
direction opposite to a noise measure: a *larger* residual means a
*stronger* push. Under R1 (declared $\sigma$ correct) it checks nothing
the definition needs. The law discarding $Q$ is right, not a loss. If R1
were ever dropped, $Q$ would become a check on $\sigma$, still not on
signal.

## 4. The only per-reading evidence is size along the line

$\max_\theta H_i=\tfrac12(T_i+|C_i|)$, 2D, median: 2.96 at $r=0$, 4.56 at
$r=1$, 20.2 at $r=3$, 202 at $r=10$. Under noise it is a $\chi^2$-like
quantity of order $d$; a push adds $r^2$. This is what ours28's
per-reading scale $\beta_i$ reads through $e^{\beta H_i/2}(1-\beta)^{k/2}$,
and it is all a single reading can offer. It separates pushes of several
noise units from noise, and nothing weaker.

## 5. Evidence that a series is noise is collective

$|C|/T=\big|\sum_iC_i\big|/\sum_iT_i$, $N=100$, $d=2$, median over 500
series:

| $r$ per reading | $\lvert C\rvert/T$ |
|---|---|
| 0 | 0.083 ($1/\sqrt{dN}=0.071$) |
| 0.3 | 0.096 |
| 0.5 | 0.130 |
| 1 | 0.344 |

Under noise the per-reading $C_i$ point in random directions and their
sum grows like $\sqrt{dN}$ against $T\sim2dN$; under a common mass they
add coherently. A push of half a noise unit per reading is invisible to
every reading on its own (§§2, 4) and visible to the sum. That is the
noise measure of the null-series record, and it is a property of the
series.

## 6. The lemma

1. Per reading, the law sees $(T_i,C_i)$; its tilt from that reading is
   $|C_i|$, its residual $(T_i-|C_i|)/2$ is discarded.
2. In 1D, $|C_i|=T_i$: a reading is always aligned and shows nothing but
   its size.
3. In $d\ge2$, misalignment exists, is uninformative about weak signal,
   and the residual is (weakly) evidence *for* a push, not against.
4. The only per-reading signal evidence is $(T_i+|C_i|)/2$, size along the
   best line; ours28 already uses it; it detects pushes of a few noise
   units and no weaker.
5. Evidence that a series is noise is $|C|/T$ against $1/\sqrt{dN}$: a
   function of the series, in every dimension.

So a noise measure and L19 collide in every $d$. The null-series record's
§5 dichotomy (keep L19 and keep the $\sqrt N$, or share one scalar) is not
a 1D artefact. The afternoon's suggestion that $d\ge2$ might let a noise
measure be built inside L19 is withdrawn.

## 7. Checks

Sandbox, numpy, seed 3. Identity of §1: one aligned-sign 2D reading,
grid of 3201 log-spaced masses over $[e^{-8},e^{8}]$. §§2–4: 20 000
readings per cell, unit noise on both channels, $\tilde{\mathbf F}=r\mathbf u+\mathbf e$,
$\tilde{\mathbf a}=r\mathbf u+\mathbf e'$, $\mathbf u$ uniform on the
sphere. §5: 500 series of 100 such readings. Script not kept.

## Sources

- `null_series_tilt_and_property_N_2026-09-29.md` §§1, 3, 5.
- `the_noise_measure_view_2026-09-29.md` §2.
- `the_combination_equation_is_exact_2026-09-26.md` §1.
- `what_forces_the_ratio_of_means_readout.md` §7.
- `ours27/ours28_reference_implementation_and_readout_2026-09-28.md` (the
  per-reading scale).
- `mass_estimator_lemmas_2026-09-29.md`: L19, R1, T25.
