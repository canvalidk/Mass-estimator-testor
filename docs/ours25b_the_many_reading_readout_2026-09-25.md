# ours25b: the many-reading readout that axis splitting forces

**2026-09-25.** Record, filed at the user's request ("Will you now write that
record. You will have to present the new equation as well"), written by Claude
from the conversation of the afternoon and evening of 25 September. It follows
`axis_splitting_and_the_24_september_change_2026-09-25.md` in the same
discussion and pays the debt that record's §4 left open: "If axis splitting is
a law, the readout owes an account of this." It names a new equation, ours25b.
It edits nothing and binds nothing.

**Names used here.**

| name | what it is |
|---|---|
| ours21 | the 21 September equation (`mass_estimator_equation/`; tester `reference: flat`) |
| ours24 | the 24 September law (tester `reference: cartesian2`) with the plain-sum readout (contribution 11, eq. 28) |
| ours25 | ours24 with sandwich-calibrated widths (tester `calibrate: sandwich`) |
| **ours25b** | ours24's law with the stacked-length readout below (tester readout `ratio_of_means_rss`) |

ours25 and ours25b change different things (widths, and the point readout), so
they can be combined. Neither changes the law.

**In one paragraph.** The readout divides the average force size by the
average acceleration size. With many readings, ours24 took "the size of the
series" to be the plain sum of the per-reading sizes. Nobody argued for that.
Axis splitting, applied to the answer rather than only to the law, rules it
out: a 2D reading and its two axes read separately must have the same size,
and $|F_1|+|F_2|$ is not the length of $(F_1,F_2)$. Given one fixed rule for
joining sets of readings, which the plain sum also assumed, the only rule left
is the root sum of squares, i.e. the length of all the readings' true vectors
stacked into one. Nothing earlier breaks: both composition laws still hold, as
long as the pooling is done on each candidate before averaging. What it buys:
a 3D series and its axis readings now agree exactly (to $10^{-15}$), the
recording frame no longer matters, and N readings in d dimensions become
exactly one reading in $Nd$ dimensions, so the answer depends on the data only
through three sums. On accuracy it ties with ours24. It is adopted because it
is forced, as the 24 September change was.

---

## 1. The equation (ours25b)

**Types.**

- $N\ge 1$ readings of one object in $d\in\{1,2,3\}$ dimensions: measured net
  force $\tilde{\mathbf F}_i\in\mathbb R^d$ and measured acceleration
  $\tilde{\mathbf a}_i\in\mathbb R^d$, any values including zero.
- $\sigma_F>0$ and $\sigma_a>0$: one per channel for the whole series, the same
  on every axis and in every reading, with independent errors. The "$>0$" is
  part of the type.
- Output: a point $\widehat m\in(0,\infty)$ and a law over $(0,\infty)$. Never
  "undefined".

**Noise units and the three sums.**

$$
\mathbf x_i=\frac{\tilde{\mathbf F}_i}{\sigma_F},\qquad
\mathbf y_i=\frac{\tilde{\mathbf a}_i}{\sigma_a},\qquad
s=\frac{\sigma_F}{\sigma_a},\qquad
m=s\tan\theta,\quad \theta\in\left(0,\tfrac{\pi}{2}\right)
$$

$$
P=\sum_{i=1}^{N}|\mathbf x_i|^2,\qquad
Q=\sum_{i=1}^{N}|\mathbf y_i|^2,\qquad
D=\sum_{i=1}^{N}\mathbf x_i\cdot\mathbf y_i,\qquad
M=Nd .
$$

**The law of mass (unchanged from ours24).**

$$
p(m)\;\propto\;\frac{s}{s^2+m^2}\;
\exp\!\left[-\,\frac{\sum_{i}\big|\tilde{\mathbf F}_i-m\,\tilde{\mathbf a}_i\big|^2}{2\left(\sigma_F^2+m^2\sigma_a^2\right)}\right]
$$

or, in the angle,

$$
p(\theta)\;\propto\;\exp\!\left[-\tfrac12\left(P\cos^2\theta+Q\sin^2\theta-2D\sin\theta\cos\theta\right)\right].
$$

The first factor is the starting belief (equal belief for equal angles,
counted once); the exponent adds one residual per reading.

**The readout (new).** Stack the N readings' true force vectors into one
vector $\mathbf F^\ast\in\mathbb R^{M}$, and the true accelerations into
$\mathbf a^\ast\in\mathbb R^{M}$. Then

$$
\boxed{\;
\widehat m\;=\;\frac{E\big[\,|\mathbf F^\ast|\,\big]}{E\big[\,|\mathbf a^\ast|\,\big]}
\;=\;s\;\frac{\displaystyle\int_0^{\pi/2}\sin\theta\;\mu_M\big(H(\theta)\big)\,p(\theta)\,d\theta}
{\displaystyle\int_0^{\pi/2}\cos\theta\;\mu_M\big(H(\theta)\big)\,p(\theta)\,d\theta}
\;}
$$

with

$$
H(\theta)^2=P\sin^2\theta+Q\cos^2\theta+2D\sin\theta\cos\theta
$$

and $\mu_M(h)$ the mean length of a Gaussian vector in $M$ dimensions:

$$
\mu_M(h)=E\,|\mathbf v|,\quad \mathbf v\sim N(\mathbf w,I_M),\ |\mathbf w|=h,\qquad
\mu_M(h)=\sqrt2\,\frac{\Gamma\!\left(\frac{M+1}{2}\right)}{\Gamma\!\left(\frac{M}{2}\right)}\;
{}_1F_1\!\left(-\tfrac12;\tfrac{M}{2};-\tfrac{h^2}{2}\right).
$$

For $M=1,2,3$ this is the folded-normal, Rice and 3D closed forms already in
the tester. At fixed $\theta$ the stacked true push, in noise units, is
$N(\mathbf W(\theta),I_M)$ with
$\mathbf W=\mathbf X\sin\theta+\mathbf Y\cos\theta$ and
$|\mathbf W|=H(\theta)$, so $E[|\mathbf a^\ast|\mid\theta]=\sigma_a\cos\theta\,\mu_M(H)$
and $E[|\mathbf F^\ast|\mid\theta]=\sigma_F\sin\theta\,\mu_M(H)$.

**ours24, for comparison.** Everything above is the same except the weight:

$$
\text{ours24:}\quad \mu_M\big(H(\theta)\big)\;\longrightarrow\;\sum_{i=1}^{N}\mu_d\big(h_i(\theta)\big),\qquad
h_i(\theta)=|\mathbf x_i\sin\theta+\mathbf y_i\cos\theta| .
$$

For one reading the two are identical. The median and every interval come from
the law alone, so they are identical for any N.

**Where each piece comes from.**

| piece | what holds it |
|---|---|
| candidates, and averaging over them | Newton II with $m>0$; $\sigma>0$ makes the average always a mass |
| divide average force size by average acceleration size | any two of lumping, common force and force/acceleration symmetry (6 Sept theorem) |
| flat weighting over each push vector | axis splitting applied to the law (25 Sept record, §4) |
| equal belief for equal angles, counted once | centre and symmetry forced by the zero-reading result; shape open |
| readings multiply, nothing per reading depends on mass | the combination rule (25 Sept convergence argument) |
| **size of the series = length of the stacked vector** | **axis splitting applied to the answer, given one joining rule (this record, §3)** |

## 2. The gap this fills

With one reading, "force size" is the length of the true force vector. With N
readings the readout needs the size of the whole series. Contribution 11's
eq. 28 took $\sum_i|\mathbf F^\ast_i|$. That was the obvious default and was
never argued.

The 25 September record measured the consequence: one 3D reading and the same
reading as three 1D readings give laws that agree to $10^{-15}$ but readouts
that differ by 0.3–0.4% at the median and 1.1–1.4% at the 90th percentile. It
left that as owed.

## 3. What forces the root sum of squares

**Premises.**

1. **One reading (step 5 as it stands).** The size of one reading's true force
   is the length of the vector.
2. **One joining rule.** When two sets of readings of the same object are
   joined, the size of the whole is computed from the sizes of the two parts
   by one fixed rule $\oplus$. The plain sum assumed this too (adding is such a
   rule), so it is not new. It only rules out exotic rules.
3. **Axis splitting, applied to the answer.** A reading and the same reading
   taken as one-dimensional readings, one per axis, must give the same mass.

**Proof.** Take a 2D reading with true force $(F_1,F_2)$. By premise 1 its size
is $\sqrt{F_1^2+F_2^2}$. Taken as two 1D readings its size is
$|F_1|\oplus|F_2|$. Premise 3 says they are equal for every candidate, so

$$
p\oplus q=\sqrt{p^2+q^2}\qquad\text{for all } p,q\ge0 .
$$

By premise 2 the same rule joins any two sets, so for any series

$$
|\mathbf F^\ast|=\sqrt{\sum_i|\mathbf F^\ast_i|^2},
$$

the length of the stacked vector. The 3D case is then a check, not a further
assumption: $\big(|a|\oplus|b|\big)\oplus|c|=\sqrt{a^2+b^2+c^2}$. The plain sum
fails at the first step, since $|F_1|+|F_2|\neq\sqrt{F_1^2+F_2^2}$.

**Why premise 3 reaches the readout.** The 25 September record stated axis
splitting for the law, but its reason was "whether the channels are recorded
as one vector or as three numbers is bookkeeping. The mass should not depend
on it." The mass is the answer, i.e. the readout. The premise's own
justification already covers it.

**Strong and weak forms.** The proof uses the strong form: one candidate has
one size, however its axes are grouped. A weaker form asks only that the final
answers agree. That form is only sketched: at strong signal the gap between
two pooling rules comes down to how their expected size changes with the
angle, and matching that for every reading appears to force the vector length.
**Sketch, not checked.** It does not affect the choice between eq. 28 and
ours25b, because the plain sum fails both forms.

**Steps 7 and 11 are one move.** Axis splitting applied to the law gives the
flat weighting over push vectors (the ours21 → ours24 change). Applied to the
answer, it gives the stacked length (ours24 → ours25b). The second depends on
the first: under ours21's law the series does not collapse to one stacked
reading (§5.3), and the stacked readout has no clean form.

## 4. Nothing earlier breaks

**Composition.** Within one candidate every reading has the same mass, so
$\mathbf F^\ast_i=m\,\mathbf a^\ast_i$ for all $i$, and any size rule $S$ that
scales ($S(c\,\mathbf v)=c\,S(\mathbf v)$) gives $S(\mathbf F^\ast)=m\,S(\mathbf a^\ast)$.

- Lumping (common accelerations): $S(\mathbf F^\ast_A+\mathbf F^\ast_B)=(m_A+m_B)\,S(\mathbf a^\ast)$, so the readouts add.
- Common force: the inverse masses add, by the mirror argument.

So the composition laws cannot tell the plain sum from the stacked length.
Axis splitting can. The 6 September theorem then applies unchanged to the
pooled pair $\big(|\mathbf F^\ast|,|\mathbf a^\ast|\big)$.

**The pooling must happen on each candidate, before averaging.** Pooling the
averages instead, $S(E\mathbf f)/S(E\boldsymbol\alpha)$, breaks lumping,
because the per-reading ratios $E[f_i]/E[\alpha_i]$ differ across readings. On
constructed laws for two bodies over three readings (Monte Carlo, 400,000
candidates):

| pooling | $R(A{+}B)$ | $R(A)+R(B)$ | common force $1/R$ vs $1/R_A+1/R_B$ |
|---|---|---|---|
| stacked length, on each candidate | 2.04069 | 2.04069 | 3.26682 = 3.26682 |
| stacked length, on the averages | 1.98882 | 1.99362 | — |
| plain sum (linear, so either way) | 2.05801 | 2.05801 | 3.26703 = 3.26703 |

**Other properties kept.** Tested in the tester (`tests/test_stacked.py`):

- one reading: ours25b equals ours24 exactly;
- all readings zero: $\widehat m=s$, for any N;
- the force/acceleration swap ($\tilde{\mathbf F}\leftrightarrow\tilde{\mathbf a}$,
  $\sigma_F\leftrightarrow\sigma_a$) gives $1/\widehat m$;
- the two cartesian implementations (`cartesian1`, `cartesian2`) give the same
  readout.

From the formula, not separately tested:

- a change of force or acceleration units rescales $\widehat m$ as a mass
  ($s$ carries the units; $\mathbf x,\mathbf y$ do not);
- convergence to the true mass with many readings: the law is unchanged, and
  once it concentrates any positive smooth weight gives the same limit.

## 5. The symmetries it gains

**5.1 Axis splitting holds for the readout.** A $d$-dimensional series and the
same series as $Nd$ one-dimensional readings give the same $\widehat m$.

- Tester, `presets/ours25b_axis_splitting.yaml` (3D; $m=s/4,s,4s$; SNR 1 and 3;
  N = 1 and 5; 400 series per cell): ours25b's largest difference over all
  cells is $3.3\times10^{-15}$; the medians agree to $5\times10^{-15}$.
  ours24's readout differs by 0.006–0.42% at the median, up to 1.37% at the
  90th percentile, and up to 4.05% in the worst series.
- The 25 September record's worked reading ($\sigma_F=0.5$, $\sigma_a=0.1$,
  $\tilde{\mathbf F}=(0.3,0.2,-0.4)$, $\tilde{\mathbf a}=(0.15,-0.05,0.02)$):
  as a 3D vector, 3.8855; as three axes, 3.9051 under ours24 (+0.50%) and
  3.8855 under ours25b.

**5.2 The recording frame does not matter.** Under ours24, a 3D reading
recorded as three 1D readings gives an answer that depends on the orientation
of the lab axes: over 200 random rotations of one reading, it moved by 2.8%.
Under ours25b it does not move. The median does not move under either, since
it comes from the law.

**5.3 The series collapses to one reading.** Under the flat-in-the-vector
weighting, N readings in $d$ dimensions have exactly the law of one reading in
$Nd$ dimensions; the law never used the grouping. With the stacked length, the
readout does not either. So ours25b depends on the data **only through
$P$, $Q$, $D$ and $M$**. Not the order of the readings, not how they are grouped,
and not how the push was shared out among them.

- Checked: two different series of two 2D readings with the same three sums
  (per-reading push sizes $2.0,\,0.5$ against $1.67,\,1.21$) give
  $\widehat m=1.2146$ under ours25b for both; ours24 gives 1.2133 and 1.2145.
- Every single-reading intuition transfers: a series behaves as one reading
  in $M$ dimensions whose push, in noise units, has the three sums as its
  lengths and dot product.

**5.4 It is invariant under any rotation of the stacked space.** Because only
$P$, $Q$ and $D$ enter, ours25b is unchanged by any rotation of all $M$
components together. Rotating each reading's axes separately, rotating the
lab, permuting readings and regrouping axes are all special cases.

**5.5 The width formula in the same form.** The 25 September record's §6
gives, per reading at the true angle, $J=d+r^{\ast2}$ and $H=r^{\ast2}$ (its
symbols). Summed over readings with pushes $r_i^\ast$:

$$
\frac{\text{true variance}}{\text{posterior variance}}=\frac{J}{H}
=1+\frac{M}{R^{\ast2}},\qquad R^{\ast2}=\sum_i r_i^{\ast2},
$$

i.e. one plus the number of components over the total squared push. It is
$1+d/r^{\ast2}$ when the readings are equal. This is algebra on the record's
formula, a law-level result, and it matches the collapse of §5.3.

## 6. What it costs: nothing on accuracy

**Standalone study** (written for this question, independent of the tester;
72 cells: 1D and 3D; $m=s/4,s,4s$; equal pushes at SNR 1 and 3, mixed pushes
SNR $U[0,4]$, and half the readings coasting with zero push; N = 2, 10, 50;
400 series per cell):

- series by series, ours24 and ours25b differ by at most 0.56% at the median,
  1.4% at the 90th percentile, 4.8% in the worst series (N = 2, weak signal);
- the gap shrinks with N and with signal;
- within-×1.25 and within-×2 scores agree to within a point in every cell,
  inside the simulation noise (about ±2 points).

**Tester** (`presets/ours24_vs_ours25b.yaml`, 54 cells: 1D and 3D; three
masses; SNR 1, 3 and $U[0,4]$; N = 2, 10, 50; 400 series): largest median gap
0.43%, largest 90th-percentile gap 1.24%, worst series 3.31%; largest
within-factor difference 0.8 points.

**Why the gap is small.** It falls quickly with signal. One 3D reading against
its axes, $m=s$, one fixed noise draw: 0.29% at SNR 1, 0.08% at SNR 4, 0.007% at
SNR 16, 0.0008% at SNR 32.

**Is it better?** On the requirements, strictly: it keeps every property
ours24's readout has and gains §5.1–5.4. On accuracy it is a tie. And it is a
correction more than an improvement: given axis splitting and one joining
rule, ours24's readout breaks a premise its own law already relies on.

## 7. What it does not fix: padding with zero readings

Append readings that read exactly zero, such as a coasting interval. The law
does not change and neither does the median, but the readout drifts. Both
poolings drift to the same limit, $s\,E[\sin\theta]/E[\cos\theta]$ under the law.

| zero readings appended | 0 | 1 | 3 | 10 | 100 |
|---|---|---|---|---|---|
| strong 3D reading, ours24 | 1.1376 | 1.1354 | 1.1336 | 1.1320 | 1.1310 |
| strong 3D reading, ours25b | 1.1376 | 1.1359 | 1.1342 | 1.1325 | 1.1311 |
| weak 3D reading, ours24 | 1.0621 | 1.0565 | 1.0533 | 1.0512 | 1.0500 |
| weak 3D reading, ours25b | 1.0621 | 1.0567 | 1.0536 | 1.0513 | 1.0500 |

(median 1.1457 and 1.0717 throughout.)

- **Under ours25b this is forced.** A zero reading still has a noise-sized
  unknown push, and the stacked length counts it. Removing the drift would
  need a size rule that ignores some components, which premise 3 forbids.
- **It is the same effect** as the dimension dependence the 25 September
  record noted: a 1D reading is not the same as $(x,0,0)$ in 3D.
- **The honest statement:** $d$ is physics, the number of directions the push
  can take (a track fixes it at 1). Axis splitting is only about bookkeeping.

## 8. What changes, what does not, what is owed

**Unchanged.** The law; the median; every interval; the zero-reading answer;
the fake-reading reading of λ; the convergence argument; the width formula;
all single-reading results.

**Changes.** Only ratio-of-means numbers computed with eq. 28 on more than one
reading. The shift is small (§6), but earlier records' many-reading
ratio-of-means tables were computed with the plain sum and are worth rerunning.

**Proposed wording for step 11 of the equation note's chain.**

> 11. **Many-reading readout.** With many readings, "force size" means the size
> of the whole series: the square root of the sum of the squared per-reading
> sizes, worked out for each candidate before averaging. Forced by axis
> splitting applied to the answer, given one rule for joining sets of readings
> (the plain sum assumed that too). Consequence: N readings in d dimensions are
> exactly one reading in N·d dimensions, and the answer depends only on three
> sums. Costs nothing on accuracy. *(proved, given axis splitting)*

**Owed.**

- The weak form of premise 3 (§3), sketched only.
- **ours25's width correction depends on grouping.** It estimates the score
  spread from per-reading scores (`sandwich_weight` in the tester), so a 3D
  reading and its three axes give different widths, and one reading gives
  none. Found by reading the code, not measured.
- Padding (§7): no tester world generates coasting readings yet.
- Scope: one $\sigma_F$ and one $\sigma_a$ per series. If precision varied
  reading by reading with a common ratio, it is the stacked length in physical
  units that stays independent of grouping; that belongs with the open
  question of different noise per instrument.
- Rerun the earlier many-reading ratio-of-means results with ours25b; a second
  seed; real data.

## 9. How we got here

1. In the afternoon session the user asked for the moves up to ours24 to be
   justified as forced moves ("put your hand on the left wall and walk"). The
   chain drafted there marked step 11, the many-reading readout, as never
   argued, and proposed the root sum of squares as the account the 25
   September record owed.
2. He asked for it to be looked through rigorously and recorded, so the agent
   writing the equation note would not be slowed by it, then asked for it
   worked through first, before any record.
3. The argument (§3) and the checks (§4–7) were done in a standalone
   implementation of ours24, which reproduces the 25 September record's
   medians (3.72 and 3.46) and agrees with brute-force sampling.
4. He pushed the ours24/ours25 tester commits and asked for the readout to go
   into the tester under a proper name without touching the old equations:
   ours25b. He asked whether it is a strict improvement (§6), then for this
   record.

## 10. Reproduction

Tester `canvalidk/Mass-estimator-testor`: commit `e288110` on top of
`19d3b41`. It holds ours25b together with the later `stacked_state` and ours26
work; ours25b's own commit, made in Claude's sandbox, was not kept.

- Readout `ratio_of_means_rss` (`met/analyst.py`), the stacked acceleration
  weight `ReadingSet.stacked_cond_alpha` and $\mu_M$ as
  `stacked_mean_radius` (`met/law.py`). $\mu_M$ for $M>3$ uses the Poisson
  mixture of central chi means (exact to about $10^{-13}$, checked against
  arbitrary precision up to $M=1200$), tabulated by a spline in
  $\operatorname{asinh}h$ (about $5\times10^{-9}$).
- `reduce: split_axes` feeds a series in as its axis readings.
- With `scores.compare_to`, ours25b's `ratio_of_means_rss` is compared with a
  reference's `ratio_of_means`.
- `tests/test_stacked.py`: 40 checks, against closed forms, a two-dimensional
  quadrature, and an independent Gauss–Legendre quadrature in $\theta$. All
  722 tests pass.

```sh
python -m pytest -q tests/test_stacked.py
python -m met run presets/ours25b_axis_splitting.yaml
python -m met run presets/ours24_vs_ours25b.yaml
```

The standalone study of §6, the rotation, collapse, padding and composition
checks were run in Claude's workspace and are not kept; the tester's tests
cover the rotation, collapse and splitting properties.
