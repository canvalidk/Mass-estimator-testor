# ours28b: ours27 with one candidate-size scale per noise pair, shared over the series

**2026-09-28 (night).** Record, written by Claude at the user's request ("So
would you be able to construct a new equation? ours28b?", then "make a
ours28b folder and put the appropriate files that should be there there").
It names ours28b and records the checks made on it. It edits nothing and binds
nothing. **ours28b is a candidate**; §5 says why it is not yet the definition.

**Where it comes from.** `ours27_one_scale_per_noise_pair_2026-09-28.md` (moved
into this folder with it) derives the structure. If the law's width must equal
the spread of its own centre over repeats of the same experiment, whatever the
pushes, then:

- no rule with separate per-reading factors works;
- a shared scale works exactly when it pools only components of equal
  sensitivity;
- the precise-null limit picks classes by declared noise pair.

ours28b is that structure written as an equation, with the arcsine scale law.
The ours28 of `../ours27/ours28_candidate_push_scale_as_nuisance_2026-09-28.md`
is a different candidate (one scale per reading); the two share only the
Gaussian candidate family.

---

## 1. The equation

**Types.**

- $N\ge1$ readings in $d\in\{1,2,3\}$ dimensions: measured
  $\tilde{\mathbf F}_i,\tilde{\mathbf a}_i\in\mathbb R^d$, any values.
- Per reading, a declared error covariance
  $\Sigma_i=\operatorname{blockdiag}(\Sigma_{FF,i},\Sigma_{aa,i})$ whose two
  blocks **commute**. This covers per-axis instruments seen in any frame, and
  one isotropic pair. Cross-channel correlation and non-commuting blocks are
  outside the type (open, §5).
- A declared starting law $\pi$, currently the half-Cauchy $c/(c^2+m^2)$.
- A declared scale law $G$ on $\beta\in(0,1)$, currently arcsine,
  $\mathrm{Beta}(\tfrac12,\tfrac12)$.

**Components and classes.** Rotate each reading into the common eigenbasis
of $\Sigma_{FF,i}$ and $\Sigma_{aa,i}$. Each component $c$ then has a
declared noise pair $(\sigma_{F,c},\sigma_{a,c})$. A **class** $g$ is the set
of all components, over the whole series, with the same pair. For a class,
$s_g=\sigma_F/\sigma_a$, noise units $x_c=\tilde F_c/\sigma_F$,
$y_c=\tilde a_c/\sigma_a$, and

$$
\theta_g(m)=\arctan\frac{m}{s_g},\qquad
A_c=x_c\sin\theta_g+y_c\cos\theta_g,\qquad
H_g(m)=\sum_{c\in g}A_c^2,\qquad M_g=\lvert g\rvert .
$$

$H_g$ depends on the data only through
$P_g=\sum x_c^2$, $Q_g=\sum y_c^2$, $D_g=\sum x_cy_c$.

**The law of mass.**

$$
p(m)\;\propto\;\pi(m)\prod_g\int_0^1G(d\beta)\,(1-\beta)^{M_g/2}\exp\!\Big[\tfrac12\beta\,H_g(m)\Big],
$$

and with the arcsine law each factor is ${}_1F_1(\tfrac12;\,1+\tfrac{M_g}2;\,\tfrac{H_g}2)$
up to a constant.

**The candidates.** Given $m$ and the class scales, independently per
component, in noise units $v_c\sim N(\beta_gA_c,\beta_g)$, with $\beta_g$
drawn from its law given $m$:
$G(d\beta)(1-\beta)^{M_g/2}e^{\beta H_g/2}$, normalised. Then
$a^\ast_c=\sigma_{a,c}\cos\theta_g\,v_c$ and $F^\ast_c=m\,a^\ast_c$.

**The readout.** Unchanged in form:

$$
\widehat m=\frac{E\lVert\mathbf F^\ast\rVert}{E\lVert\mathbf a^\ast\rVert}
=\frac{\int m\,A(m)\,p(m)\,dm}{\int A(m)\,p(m)\,dm},
\qquad A(m)=E\big[\lVert\mathbf a^\ast\rVert\;\big|\;m\big],
$$

over the stacked candidates. $A$ is computed from the Laplace transform of
$\lVert\mathbf a^\ast\rVert^2$, which factorises over classes given $m$, and
within a class depends on the data only through $(H_g,M_g)$.

**Carried state.** $(P_g,Q_g,D_g,M_g)$ per noise pair, with the pair itself.
Combining series is addition (T11, T12 per pair).

**Reductions.**

- $G$ a point mass at $\beta=1$: ours27 on this input type (ours26 when the
  noise is diagonal in the recording axes, ours25b with one pair).
- One noise pair: a single class, law
  $\pi(m)\,{}_1F_1(\tfrac12;1+\tfrac M2;\tfrac{H(\theta)}2)$ with
  $H=P\sin^2\theta+Q\cos^2\theta+2D\sin\theta\cos\theta$.

## 2. Checks

Reference implementation `ours28b.py` and `ours28b_checks.py` in this folder.
The law is on an 801-point grid in $\log m$ (half-width 9), the scale on a
1601-point grid in $\operatorname{logit}\beta$, and $A(m)$ by a 721-point
trapezoid in $\log t$.

| check | data | result |
|---|---|---|
| flat limit against the ours25b closed form for $A(m)$ | one pair, two 3D readings | max relative error $1.4\times10^{-4}$ (far tails) |
| swap (channels exchanged, $c\to1/c$) | one pair, six 3D readings | readout product 1.000000, median product 1.000000 |
| units (force ×3) | same | readout ratio 3.000000 |
| each reading rotated in its own frame, with its covariance | same; and per-axis $\sigma_a=(0.2,0.5,1.25)$ | readout changes by $9\times10^{-14}$ and $4\times10^{-13}$ |
| axis splitting | six 3D readings as eighteen 1D readings | identical |
| zero experiment alone | three 3D zero readings | law = starting law (median 2.000000 at $c=2$); readout exactly $s$ |
| precise nulls tend to known-empty | 1D reading $F=3$, $a=1$ against the same with two 0/0 axes of noise $\varepsilon$ | $\varepsilon$ = 0.1, 0.01, 0.001: law equal to the 1D law; readout 2.1126, 2.1140, 2.1140 against 2.1140 |
| width calibration | from the derivation record §5 | 94–95% on the regular counterexample at 100–1600 repetitions; 94–98% one-pair |

Against ours27 on the same data:

| data | ours27 | ours28b |
|---|---|---|
| one pair, six weak 3D readings, $c=2$ | readout 4.194, median 4.276, 95% [1.61, 31.4] | readout 3.266, median 3.224, 95% [0.22, 51.9] |
| per-axis, three 3D readings, $c=2$ | readout 1.923, median 1.971, 95% [0.84, 3.93] | readout 1.909, median 1.945, 95% [0.58, 4.34] |

With weak data ours28b is wider and leans on the starting law; with stronger
data the two nearly agree.

## 3. What it keeps and what it gives up

**Kept:** T7 (the scale factor is mass-free), swap, units, rotation,
axis splitting, grouping within a pair (enough for the C3 derivation of
`../ours27/ours27_c3_from_betweenness_and_the_scale_fork_2026-09-28.md` with one
pair), consistency, the per-pair four-number state, and the weak
zero-reading property.

**Given up:** the strong zero-reading property. Appending one exact 3D zero
reading to the six-reading example moves the median 3.224 → 2.970 and the
readout 3.266 → 3.054. The derivation record §2 and §7 give the reason: no
rule in which a zero reading contributes a constant factor can have
calibrated widths.

## 4. Premises this rests on that are not yet ruled

- **Calibration:** the width equals the spread of the law's own centre over
  repeats of the same experiment, whatever the pushes, in the regular limit.
- **Dropping the strong zero-reading property** in favour of the weak one.

## 5. Why it is a candidate, not the definition

1. **A jump at equal noise pairs.** Classes by exact equality of pairs make
   the law discontinuous in the declared noise. With two 0/0 axes of noise
   0.999 beside the 1D reading above, the law is the 1D law (median 2.125). At
   exactly 1.0 they join the signal axis's class (median 1.822). Pooling by
   noise ratio would remove the jump but fails the precise-null limit. The
   likely repair is a joint law on the pairs' scales that couples them as the
   pairs converge. That is a new shape, not yet found.
2. **General covariance.** Cross-channel correlation and non-commuting noise
   blocks have no class structure yet. The candidate is the generalized
   eigendirections of the residual covariance $V$ and the candidate
   covariance $K$.
3. **Two shapes.** The starting law's and the scale law's. Both are symmetric
   and matter only at weak signal. Under the project's rule, each remaining
   freedom is a rival definition that agrees wherever the evidence identifies
   mass.
4. **Not done:** the readout point's own interval; the perturbation matrix
   rerun on ours28b; the tester build; the over-coverage at weak pushes with
   moderate $N$ (derivation record §5).

## 6. How we got here

The derivation record was filed the same night. The user asked whether a new
equation could be constructed and named it ours28b. Claude wrote the reference
implementation, validated its flat limit against the closed form, and ran the
checks of §2. The half-grid-step median bias noted in the perturbation matrix
record appeared again (median swap product 0.978 before correcting the
cumulative sum to midpoints) and was fixed before the table was read. The jump
of §5.1 was found while checking the precise-null limit. Asked whether this
is "it", Claude said: nearly, with one noise pair; not with several
instruments. The user then asked for this folder.
