# Calibration forces the size measure: one scale per noise pair, shared over the series

**2026-09-28 (night).** Record, written by Claude after the user's instruction
"continue to work on the solution", in the session that filed
`ours27_c3_from_betweenness_and_the_scale_fork_2026-09-28.md`. It edits nothing
and binds nothing. It uses the score and curvature identities of
`ours27_full_covariance_width_and_candidate_geometry_2026-09-28.md` and the
regular counterexample of `ours27_limits_of_scalar_and_tensor_size_repairs_2026-09-28.md`
as test cases.

**In one paragraph.** Ask the law of mass for one thing: its width should equal
the spread of its own centre over repeats of the same experiment, whatever the
pushes were, in the regular many-reading limit (the Bartlett identity at fixed
candidates). That requirement decides the size measure's structure. No rule
whose readings contribute separate factors can meet it: the per-reading
identity reduces to an Ornstein–Uhlenbeck equation whose only smooth radial
solutions are constants. Shared scales meet it exactly when each scale is
shared only by components of equal sensitivity to mass, for any push sizes
and any proper scale law. The precise-null limit (a measured-empty axis whose
noise shrinks to zero must approach the known-empty answer) then fixes the
classes as the **declared noise pairs**. So the candidate measure is one
Gaussian scale per noise pair, shared over the series. Simulations agree with
the limit theory to about a point everywhere tested. The per-pair sums
$(P,Q,D,M)$ survive as the carried state. The price is the strong zero-reading
property. The weak one (a zero experiment alone gives the starting law) holds.

---

## 1. The requirement

For a reading $i$ with candidate $a_i^*$, the flat score is
$U=\sum_ib_i^\top V_i^{-1}r_i$ with $b_i$ and $r_i$ independent, and exactly
$\operatorname{Var}U=J+G$ against expected curvature $J$ (the width record
§2). A rule that replaces the flat candidate measure changes both. Call a rule
**calibrated** if, for every sequence of fixed candidates along which the
law concentrates regularly, the variance of its score equals its expected
curvature to leading order. Then its law's variance equals the sampling
variance of its centre.

- The requirement uses repeats of one experiment, so it needs no mass behind
  the readings. It asks that the width the law reports be the spread the law's
  own centre shows.
- "Every sequence of fixed candidates" means no assumption about how the
  experimenter pushed: equal, mixed, weak or strong.
- It is a regular-limit requirement. At weak signal the Fieller obstruction
  stands (width record §6), and the shapes of $\pi$ and of the scale law
  matter there.

## 2. Rules with separate per-reading factors cannot be calibrated

Take one noise pair and a reading of $d$ components in noise units. At the
tested mass write $A\in\mathbb R^d$ for the components along the candidate
line and $B\in\mathbb R^d$ across it; $A\sim N(\mu,I)$, $B\sim N(0,I)$,
independent, with $\mu$ the candidate. Let the reading's evidence be
$e^{\phi(A)}$ for any candidate measure, and $w$ the sensitivity
$d\theta/d\log m$. The reading's score is $w\,\nabla\phi(A)\cdot B$, and

$$
V=w^2E|\nabla\phi|^2,\qquad
S=w^2E\big[A\cdot\nabla\phi-\Delta\phi\big].
$$

A per-reading rule is calibrated for every push sequence only if $V=S$ at
every $\mu$. Taking all readings with the same push shows this. The Gaussian
location family is complete, so

$$
|\nabla\phi|^2-A\cdot\nabla\phi+\Delta\phi=0 .
$$

With $u=e^{\phi}$ this is $\Delta u=A\cdot\nabla u$. Own frames (ruled
28 Sept) make $u$ radial, $u=u(\rho)$, and then

$$
\big(\rho^{\,d-1}e^{-\rho^2/2}u'\big)'=0
\quad\Rightarrow\quad
u'=C\rho^{\,1-d}e^{\rho^2/2}.
$$

Evidence from any candidate measure is smooth and even, so $u'(0)=0$ and
$C=0$. **Only the constant evidence, which carries no information, is
calibrated at every push size.**

Numerically, the per-reading ratio $V/S$ (one pair, 3D readings, push length
$r$):

| scale law per reading | $r=0.5$ | $r=1$ | $r=2$ | $r=3$ | $r=5$ | $r=8$ |
|---|---|---|---|---|---|---|
| arcsine $\mathrm{Beta}(\tfrac12,\tfrac12)$ | 5.40 | 1.90 | 1.08 | 0.97 | 0.97 | 0.99 |
| $\mathrm{Beta}(\tfrac14,\tfrac14)$ | 4.72 | 1.78 | 1.08 | 0.99 | 0.98 | 0.99 |
| $\mathrm{Beta}(2,2)$ | 6.04 | 2.02 | 1.02 | 0.89 | 0.89 | 0.95 |
| reference, about $\mathrm{Beta}(1,0)$ | 7.72 | 2.56 | 1.28 | 1.08 | 1.01 | 1.00 |
| $\mathrm{Beta}(\tfrac12,3)$ | 2.49 | 0.91 | 0.60 | 0.65 | 0.80 | 0.91 |
| $\mathrm{Beta}(3,\tfrac12)$ | 8.80 | 2.87 | 1.34 | 1.09 | 1.00 | 0.99 |

Every rule is close at strong pushes and fails at weak ones. A single weak
reading cannot tell how much of itself is push, so it cannot discount its own
noise correctly.

## 3. Shared scales are calibrated exactly when they do not mix sensitivities

Let a scale $\beta_g$ be shared by a class $g$ of $M_g$ components with total
candidate energy $R_g$ (noise units). The evidence of the class is

$$
F_g(H_g)=\int_0^1G(d\beta)\,(1-\beta)^{M_g/2}e^{\beta H_g/2}
$$

for a proper scale law $G$ with full support. As the class grows, the
posterior of $\beta_g$ concentrates at $c_g=R_g/(R_g+M_g)$, so the class's
score is $c_g$ times its flat score. If every component of the class has the
same sensitivity $w_g$:

$$
V_g=c_g^2\,w_g^2(R_g+M_g)=c_g\,w_g^2R_g=S_g .
$$

This holds **for any push sizes inside the class**: equal, mixed, or with
some pushes zero. If two sensitivities share one scale, the common $c$ is the
pooled energy fraction, and $V=S$ needs $c=\sum w_g^2R_g/\sum w_g^2(R_g+M_g)$,
which is a different number unless the $w_g$ or the $R_g/M_g$ agree. That is
the regular counterexample of the repair-limits record, reproduced here
exactly:

| rule, on that counterexample | variance ratio | SD ratio |
|---|---|---|
| flat (ours27) | 3.1287 | 1.769 |
| one scale for all components | 2.6072 | 1.615 |
| one scale per reading, across its axes | 2.3616 | 1.537 |
| one scale per component | 1.6031 | 1.266 |
| **one scale per noise pair, shared over the series** | **1.0000** | **1.000** |

The scale law's shape drops out of the limit, as the starting law's does.

**Why "over the series".** Any partition of a class into groups that each
grow is also calibrated in the limit. But readings are exchangeable: arrival
order carries nothing (T20), and nothing declared marks out blocks of
readings. The only partitions that don't depend on arbitrary grouping are
the whole class, or single readings, and single readings fail by §2.

## 4. Which classes: the declared noise pairs

Sensitivity depends only on the ratio $s=\sigma_F/\sigma_a$, so §3 alone
allows classes by ratio or by pair. The precise-null limit decides. Take an
active axis and two measured-empty axes whose force and acceleration noise
are both $\varepsilon$. That is the same ratio as the active axis, with a
different pair.

- **By ratio:** the empty axes join the active axis's class. They raise its
  $M$ to 3 at every $\varepsilon$, since noise units cannot see $\varepsilon$.
  The law never approaches the 1D answer.
- **By pair:** the empty axes form their own class. At exact zeros its
  evidence is the constant $F(0)$ for every $\varepsilon$, so the law is the
  1D law exactly. Their candidates' physical sizes are $\varepsilon$ times
  noise-unit sizes, so the readout tends to the 1D readout.

So measured-empty approaches known-empty only with classes by declared noise
pair. The same rule means precise empty axes never discount a signal axis
measured with a different instrument.

## 5. Finite samples

Coverage of the 95% law, half-Cauchy starting law at the true mass,
arcsine scale law.

**The regular counterexample** ($s=(1,10)$, $r^2=(0.1,9.9)$ per repetition):

| repetitions | flat | one scale for all | per reading | per component | **per noise pair, shared** |
|---|---|---|---|---|---|
| 20 (600 series) | 84.3 | 90.2 | 93.5 | 96.2 | 98.7 |
| 100 (600) | 73.3 | 78.3 | 79.0 | 88.0 | **94.3** |
| 400 (500) | 70.6 | 75.8 | 77.4 | 87.4 | **95.4** |
| 1600 (300) | 71.3 | 75.7 | 79.7 | 88.0 | **94.7** |
| limit (§3) | 73.2 | 77.6 | 79.8 | 87.8 | 95 |

**One noise pair**, 3D, 200 readings, fixed push lengths, random directions
(500 series):

| pushes | flat | one scale shared over series | per reading |
|---|---|---|---|
| all $r=1$ | 67.2 | 98.0 | 84.4 (limit 84.5) |
| alternating $r=3$, $r=0.5$ | 87.8 (limit 87.3) | 94.0 | 94.6 |
| all $r=3$ | 92.8 | 96.0 | 94.8 |

The shared rule over-covers at weak pushes with moderate $N$. The scale's own
uncertainty (about 17% relative at $r=1$, $N=200$) widens the law, and this
shrinks with $N$. Per reading is good when some pushes are strong and fails
when all are weak, as §2 predicts.

## 6. The equation this gives

With one noise pair, in noise units, $\theta=\arctan(m/s)$:

$$
p(m)\;\propto\;\pi(m)\int_0^1G(d\beta)\,(1-\beta)^{M/2}\exp\!\Big[\tfrac12\beta\,H(\theta)\Big],
\qquad
H=P\sin^2\theta+Q\cos^2\theta+2D\sin\theta\cos\theta,
$$

with the candidate components $v\mid\beta,m\sim N(\beta A,\beta)$ and the
readout unchanged in form, $E\lVert\mathbf F^\ast\rVert/E\lVert\mathbf a^\ast\rVert$
over the stacked candidates. For the arcsine law the integral is
${}_1F_1(\tfrac12;\,1+\tfrac M2;\,\tfrac H2)$. With diagonal per-axis noise,
$p(m)\propto\pi(m)\prod_gF_g(H_g)$ over the declared noise pairs, each with
its own $\theta_g=\arctan(m/s_g)$ and its own scale.

What carries over from the chain:

- **Sufficiency:** $(P_g,Q_g,D_g,M_g)$ per pair is the carried state, and
  combination is still addition (T11, T12).
- **T7:** the factor $(1-\beta)^{M/2}$ does not depend on mass.
- **Swap and units:** a pair maps to a pair; the scale in noise units is
  shared by both channels.
- **Axis splitting and grouping:** the class structure is grouping-free, so
  the law is invariant under any relisting or regrouping that keeps each
  component's pair. That is enough for step 2 of the C3 derivation in the
  betweenness record: its rotation closure used only regroupings of same-axis
  components across readings, which are within a pair. With one pair this is
  complete; with per-axis noise it also needs each reading's rotation, with
  its covariance, to carry the classes along, which waits on §8.
- **Kingman:** within a pair, components exchangeable and spherically
  symmetric under regrouping give exactly one Gaussian scale per pair. §3
  reaches the same structure from calibration alone.
- **Consistency:** each class's tempering tends to a positive constant.

## 7. The price: the strong zero-reading property

Appending an exact zero reading raises $M_g$ for its pairs and tempers them,
so the law changes. What survives is the weak form: a zero experiment alone
gives the starting law, since $H_g=0$ for every $m$ makes each $F_g$ constant.
§2 is why the strong form can't be kept alongside calibration. Rules in which
a zero reading contributes a constant factor are per-reading rules, and none
of them is calibrated. A noisy null does add pure noise to the score
($\operatorname{Var}AB=1$ at zero push, curvature 0), so a calibrated law has
to widen for it. The exact zero is the limit of such nulls.

## 8. Open

- **Full covariance.** Classes are declared pairs when the noise is diagonal
  in some frame. With a general $\Sigma_i$ the classes must be defined
  covariantly. The natural candidate is the generalized eigendirections of
  the residual covariance $V$ and the candidate covariance $K$, where
  sensitivity is uniform. Not worked out.
- **Shapes.** The starting law's shape and the scale law's shape, both
  weak-signal only. The joint law of several pairs' scales is also free at
  finite $N$ (independent is the natural choice).
- **Finite $N$.** Over-coverage at weak pushes with moderate $N$ (§5).
- **Not done:** the readout point's own interval under this rule; the
  perturbation matrix on it; the tester build.

## 9. Reproduction

Scripts in Claude's workspace, not kept; numpy and scipy.

- §2 and §3 ratios: $V$ and $S$ by Monte Carlo over $A$ (200,000 to 400,000
  draws) with $E[\beta\mid H]$ and its derivative from Kummer functions;
  plug-in $c_g$ for the shared rules.
- §5: the law on a 2401-point grid in $\log m$ over $[-6,6]$, shared scales
  by a 400-point uniform-$\varphi$ grid in log space, per-reading and
  per-component rules by ${}_1F_1$.

## 10. How we got here

After the betweenness record the user said four new records had dropped and
to keep working on the solution. The width record's exact $J+G$ against $J$
suggested asking which candidate measures restore the Bartlett identity for
fixed candidates. The plug-in argument of §3 and the counterexample's
numbers came first, then the §2 impossibility, then the coverage runs, then
the precise-null test that picks pairs over ratios.
