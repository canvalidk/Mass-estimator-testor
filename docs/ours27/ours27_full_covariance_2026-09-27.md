# ours27: the estimator with a full noise covariance

**2026-09-27.** Record, filed at the user's request ("Ok then lets just create
ours27. With the definite improvements. We'll deal with the rest later."),
written by Claude from the conversation of the evening of 27 September. It names
ours27 and records the checks made on it. It follows
`ours26_per_axis_noise_2026-09-26.md` and
`ours26_full_noise_covariance_2026-09-27.md` in this folder. It edits nothing
and binds nothing.

**Names used here.**

| name | what it is |
|---|---|
| ours25b | one noise pair for every component; stacked-length readout |
| ours26 | one noise SD per component, independent errors |
| **ours27** | ours26 with a declared full covariance of the errors, per reading |

**In one paragraph.** ours26 attaches the noise to the instrument's axes and
assumes the components' errors are independent. That model is not closed under
a rotation of the axes: rotate a reading and its per-axis noise together and the
noise becomes correlated, which ours26 cannot represent, so ours26 gives a
different answer in a different frame. Mass does not depend on how the axes were
drawn, so the noise model has to be one that rotates with the data. That is a
full covariance $\Sigma_i$ of the error pair $(\mathbf e_F,\mathbf e_a)$ for
each reading. Everything else in the equation is forced by the same premises as
before and is unchanged: the law is $\pi$ times one Mahalanobis misfit per
reading, the readout is the ratio of the mean stacked sizes. ours27 reduces to
ours26 when every $\Sigma_i$ is diagonal and to ours25b when they are all the
same isotropic pair. Checked: under a joint rotation of data and covariance the
readout is unchanged to eight digits where ours26's moves by 25%; a zero reading
leaves the law unchanged for any $\Sigma$; the swap and units hold under a
correlated $\Sigma$; on one-pair data ours27 agrees with ours25b's closed form to
$10^{-9}$. What ours27 does not change: the push scale ($\omega\to\infty$), the
starting law's shape and centre, C3, and the width problem. Those are dealt with
later, by the user's ruling.

---

## 1. The equation (ours27)

**Types.**

- $N\ge1$ readings of one object in $d\in\{1,2,3\}$ dimensions: measured force
  $\tilde{\mathbf F}_i\in\mathbb R^d$ and acceleration $\tilde{\mathbf a}_i\in\mathbb R^d$,
  any values including zero, in one set of axes.
- For each reading a declared, positive-definite covariance
  $\Sigma_i\in\mathbb R^{2d\times2d}$ of the error pair
  $(\mathbf e_{F,i},\mathbf e_{a,i})$, written in blocks
  $\Sigma_{FF},\Sigma_{Fa},\Sigma_{aF},\Sigma_{aa}$. Correlation between axes
  and between the two channels is allowed. Errors are Gaussian, zero-mean, and
  independent across readings.
- A declared starting law $\pi$ over $(0,\infty)$; currently
  $\pi(m)=c/(c^2+m^2)$ with a declared centre $c>0$.
- Output: a point $\widehat m\in(0,\infty)$ and a law over $(0,\infty)$.

**The residual and its covariance.** For reading $i$ at candidate mass $m$,

$$
\mathbf r_i(m)=\tilde{\mathbf F}_i-m\,\tilde{\mathbf a}_i,\qquad
V_i(m)=\Sigma_{FF,i}-m\big(\Sigma_{Fa,i}+\Sigma_{aF,i}\big)+m^2\,\Sigma_{aa,i},
$$

$$
Q_i(m)=\mathbf r_i(m)^\top V_i(m)^{-1}\,\mathbf r_i(m).
$$

**The law of mass.**

$$
p(m)\;\propto\;\pi(m)\,\exp\!\Big[-\tfrac12\sum_{i=1}^N Q_i(m)\Big].
$$

**The candidates at fixed mass.** Write $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)$
and $B_m=(mI_d;\,I_d)$, so that $\mathbf y_i=B_m\mathbf a^\ast_i+\mathbf e_i$.
With the push weighted flat, the true acceleration of reading $i$ given the data
and $m$ is Gaussian,

$$
\mathbf a^\ast_i\mid m\;\sim\;N(\mathbf b_i,K_i),\qquad
K_i=\big(B_m^\top\Sigma_i^{-1}B_m\big)^{-1},\qquad
\mathbf b_i=K_iB_m^\top\Sigma_i^{-1}\mathbf y_i ,
$$

independent across readings, and $\mathbf F^\ast_i=m\,\mathbf a^\ast_i$.

**The readout.**

$$
\boxed{\;
\widehat m=\frac{E\lVert\mathbf F^\ast\rVert}{E\lVert\mathbf a^\ast\rVert}
=\frac{\displaystyle\int_0^\infty m\,A(m)\,p(m)\,dm}{\displaystyle\int_0^\infty A(m)\,p(m)\,dm},
\qquad
A(m)=E\big[\lVert\mathbf a^\ast\rVert\,\big|\,m\big]\;}
$$

with $\mathbf a^\ast\in\mathbb R^{Nd}$ the stacked true accelerations, as in
ours25b and ours26.

**How $A$ is computed.** $X=\lVert\mathbf a^\ast\rVert^2$ is a quadratic form in
a Gaussian vector. Its Laplace transform is closed: with $\kappa_{ij}$ the
eigenvalues of $K_i$ and $\beta_{ij}$ the components of $\mathbf b_i$ in the
eigenbasis,

$$
E\,e^{-tX}=\prod_{i,j}\big(1+2t\kappa_{ij}\big)^{-1/2}
\exp\!\Big(-\frac{t\,\beta_{ij}^2}{1+2t\kappa_{ij}}\Big),
\qquad
A=\frac{1}{2\sqrt\pi}\int_0^\infty\big(1-E\,e^{-tX}\big)\,t^{-3/2}\,dt .
$$

One trapezoid sum in $\log t$, centred on $t=1/E[X]$ (§5).

**Reductions.**

- $\Sigma_i=\mathrm{diag}(\sigma_{F,c}^2,\sigma_{a,c}^2)$: $V_i$ is diagonal
  with entries $\sigma_{F,c}^2+m^2\sigma_{a,c}^2$, so $Q_i$ is ours26's sum of
  Beale terms; $K_i$ is diagonal with entries $\sigma_{a,c}^2\cos^2\theta_c$ and
  $\mathbf b_i$ is ours26's $\sigma_{a,c}\cos\theta_c\,w_c$. **ours27 = ours26.**
- All $\Sigma_i$ equal and isotropic: **ours27 = ours25b**, with the closed-form
  weight $\sigma_a\cos\theta\,\mu_M(H)$.

**Where each piece comes from.**

| piece | status |
|---|---|
| divide mean true force size by mean true acceleration size | forced, unchanged (composition, C3) |
| size of a series = length of the stacked true vector | forced, unchanged (axis splitting on the candidates, one joining rule) |
| flat push weight | forced, unchanged (isotropy, no push scale) |
| the weight cancels each reading's determinant, so each reading's evidence is $e^{-Q_i/2}$ | forced, unchanged (nothing per reading may depend on $m$; T7) |
| **the error model is a full covariance per reading** | **forced by coordinate independence (§2); new in ours27** |
| readings multiply, prior counted once | forced, unchanged |
| $\pi$: symmetric, declared before data, scale from the instruments | unchanged; shape and centre open |

## 2. Why the covariance is forced

A mass does not depend on how the axes were drawn. So the same readings
expressed in rotated axes, with the noise description rotated along with them,
must give the same law and the same point. ours26's noise description is one SD
per axis with independent errors. Rotating that description gives a covariance
with off-diagonal terms, which ours26 cannot hold. Two things follow:

- **ours26 is frame-dependent.** Applied in the rotated frame, with only the
  diagonal of the rotated covariance available to it, it gives a different
  answer (table below). This is the loss ours26's own record listed under
  "lost: rotation of the lab frame".
- **The covariance form is frame-covariant.** $Q_i$ is a Mahalanobis distance:
  rotate $\mathbf r_i$ and $V_i$ together and it is unchanged. $A$ depends on
  $\mathbf a^\ast$ only through its length, and $K_i$, $\mathbf b_i$ rotate
  with the data. So both the law and the readout are unchanged.

Checked, one 2D reading, $\tilde{\mathbf F}=(2,1)$, $\tilde{\mathbf a}=(0.9,0.3)$,
$\sigma_F=(1,1)$, $\sigma_a=(0.2,1.25)$ in the instrument's axes, $c=2$; the
reading and its covariance rotated together by $\varphi$:

| $\varphi$ (rad) | ours27 readout | ours27 median | ours26 on the rotated diagonal |
|---|---|---|---|
| 0 | 1.85529234 | 1.878668 | 1.855292 |
| 0.5 | 1.85529234 | 1.878668 | 2.247868 |
| 1.0 | 1.85529234 | 1.878668 | 2.339910 |
| 2.2 | 1.85529234 | 1.878668 | 2.298894 |

With an already-correlated $\Sigma$ (correlations up to 0.94 between channels)
rotated by 0.9 rad: 1.834154 before and after.

The same premise makes Galilean boosts and reflections trivial: a boost changes
neither force nor acceleration, and a reflection is a rotation for everything
here, since only lengths and Mahalanobis distances enter.

**Correlation between the channels is allowed.** The shared-quantity argument
requires that force and acceleration each be *obtainable* without the other. It
does not require their errors to be statistically independent; shared
calibration, timing or frame information may correlate them, and the record
that states the argument says so (`vd_shared_quantity_argument…` §5). ours27
represents that correlation; the channels stay independently determined.

## 3. What holds

All checked with the reference implementation of §5.

- **Zero reading leaves the law unchanged, for any $\Sigma$.** $\mathbf r_i=0$
  gives $Q_i=0$ at every $m$. Checked with a correlated $\Sigma$ at
  $m=0.3,1,5$: $Q=0$ exactly. Without the determinant term the reading would
  have contributed $\lvert B_m^\top\Sigma^{-1}B_m\rvert^{-1/2}$ = 0.0848, 0.0665,
  0.0345 at those masses, which is the term T7 forces the weight to cancel.
- **Swap.** Exchange $\tilde{\mathbf F}\leftrightarrow\tilde{\mathbf a}$,
  $\Sigma_{FF}\leftrightarrow\Sigma_{aa}$, $\Sigma_{Fa}\leftrightarrow\Sigma_{aF}$,
  $c\to1/c$: product of the two readouts 0.99999990 (correlated $\Sigma$;
  quadrature-limited).
- **Units.** Force readings and force covariance blocks scaled by 3, $c$ by 3:
  readout ratio 3.00000000.
- **Reduction to ours25b.** Three 2D readings, one pair $\sigma_F=0.7$,
  $\sigma_a=0.3$: ours27 4.0910703550, ours25b closed form on the same grid
  4.0910703515. $A(m)$ against the closed form $\sigma_a\cos\theta\,\mu_M(H)$:
  relative difference under $2.7\times10^{-6}$ everywhere the law has mass.
- **Reduction to ours26.** At $\varphi=0$ in the table above the two agree to
  the printed digits.
- **The weight against Monte Carlo.** Correlated $\Sigma$, $m=1.7$: quadrature
  1.45126, 400,000 draws 1.45032.
- **Totality.** $A(m)=O(1/m)$ still holds: as $m\to\infty$ every $K_i\to0$ and
  $\mathbf b_i\to0$ like $1/m$, so $mA$ stays bounded and the readout is finite
  for any proper $\pi$ (T13).
- **Combination.** Per reading nothing changes. Across readings that share one
  covariance, the record `ours26_full_noise_covariance_2026-09-27.md` gives the
  carried state $(S=\sum_i\mathbf y_i\mathbf y_i^\top,\,n)$ per covariance group,
  which reduces to $(P,Q,D,M)$ for one isotropic pair. Not built in the tester.

## 4. What ours27 does not change, and one thing it makes visible

**Not changed, by ruling (27 Sept: "we'll deal with the rest later").**

- The push scale. ours27 keeps $\omega\to\infty$. The many-reading width
  problem (T18; `README.md` failure (b)) stands. The learned-$\omega$ study
  stands as a hint (`ours26_learned_omega_is_a_hint_2026-09-27.md`). One point
  settled tonight and worth keeping: any push scale must live in **noise
  units**, because a scale in physical units would sit on one channel and break
  the swap, and no built-in scale is allowed. So the physical-units branch of
  the fork is closed by the premises, not by a test.
- $\pi$'s shape, and its centre when the instruments' ratios differ. Under a
  general $\Sigma$ the per-axis ratios $s_k$ are not even defined, so the
  arc-length candidate of `ours26_starting_law_2026-09-27.md` §5 needs
  restating before it can be used with ours27. Open.
- C3's justification (lemma record §3.2). Open.
- The axis-splitting premise as *stated* (a reading and its axes read as 1D
  readings give the same law) only makes sense when $\Sigma_i$ is diagonal in
  the chosen axes. In ours27 its two jobs are done by other premises: the law
  comes from the flat weight plus T7 (the determinant cancellation, which is
  what `ours26_full_noise_covariance` derives), and the stacked size comes from
  axis splitting applied to the candidates, which is noise-free and survives.
  The lemma record's L8 owes a restatement.

**Made visible.** With ours26 and ours27 an axis may be measured by a much
worse instrument than the others. A zero reading on such an axis leaves the law
unchanged and moves the point. One reading at SNR 3 on a good axis, $c=2$;
second axis added, reading zero, 1000 times noisier:

| second axis | readout | median |
|---|---|---|
| none | 2.907 | 2.899 |
| instrument ratio 2 (same as the good axis) | 2.881 | 2.899 |
| instrument ratio 10 | 3.669 | 2.899 |
| instrument ratio 0.4 | 2.479 | 2.899 |

At weak signal the moves are larger (2.69 → 3.90 or 1.67). The mechanism is
the null result's, applied to one axis: with no push scale, the only scale for
the unread push on that axis is its instrument's noise, so the stacked size is
dominated by that axis and the readout becomes the law's mean tamed at that
instrument's ratio. The shift saturates as the instrument gets worse; its limit
is set by the instrument's ratio and the law. Under the premises this is forced,
not a defect: composition, axis splitting, the swap and no built-in scale
together imply it, and a push scale in noise units would only rescale it.

What it leaves for a ruling is the premise the user stated on 27 Sept: "if you
get no additional information from the other dimensions, the answer must not
change." That holds for the law and cannot hold for the point under the other
premises. The user's reading on the night: the property is about the law; a
zero reading on a real instrument carries a little information (about the push),
and the point prices it in. Recorded as his position, not yet as a ruling.

## 5. Numerics of the reference implementation

Written in Claude's workspace on 27 Sept (`ours27_reference.py`, handed to
the user; not in this repo). The tester build is owed (§6).

- The law on a grid of 700 to 4000 points in $\log m$ over $\pm9$ around $c$.
- $A(m)$ by the trapezoid rule in $u=\log t$, step 0.01 over $[-26,26]$, with
  $t$ rescaled by $1/E[X]$ so the transition of $1-Ee^{-tX}$ sits at $u=0$.
  Without that rescaling the integral loses digits where $A$ is small (relative
  error $10^{-3}$ at $m\sim10^3c$); with it, $2.7\times10^{-6}$ everywhere
  tested.
- Accuracy of the whole readout on one-pair data against the closed form:
  $10^{-9}$ (§3).

## 6. What is owed

- **The tester build.** A per-reading covariance in `ReadingSet`; `Q` as a
  Mahalanobis form; the weight $A$ by the Laplace transform with a general $K$;
  worlds with correlated axis noise and with force–acceleration correlation;
  tests for the reductions, the rotation, the swap, units and the zero reading.
  The user's clone of `canvalidk/Mass-estimator-testor` is where it goes.
- **The width formula under correlation.** The per-component form found in the
  sweep, $1+\sum_c\theta_c'^2/\sum_c\theta_c'^2r_c^{\ast2}$, is for a diagonal
  $\Sigma$. Its covariance form is not written.
- **The starting law under a general $\Sigma$**: what replaces the per-axis
  ratios in the centre and in the arc-length candidate.
- **L8 restated** for the covariance form (§4).
- **Rerun** the ours26 preset and the coverage cells under ours27 with rotated
  and correlated worlds; expect no change on diagonal worlds by construction.
- **The ruling** on the no-information property (§4).

## 7. How we got here

1. The sweep of the flat records against the lemma list (afternoon, filed) left
   the equation as ours26 and listed what was open.
2. The user redirected the evening from "is the equation forced" to "what is
   mass": the definition must equal how mass is read from data, and progress
   comes from recognising properties of mass. Two he named: the answer cannot
   depend on the choice of axes, and dimensions that add no information must
   not change it.
3. Claude carved the equation from the properties and found that the first of
   them forces the full covariance and that the covariance form is covariant
   (§2). The second holds for the law and fails for the point (§4); Claude first
   called that a flaw and proposed a physical push scale; the user pointed out
   that instrument information is legitimate, which closed the physical-units
   branch (§4) and turned the effect into a forced consequence.
4. The user: "lets just create ours27. With the definite improvements. We'll
   deal with the rest later." This record.

## 8. Sources

- `ours26_per_axis_noise_2026-09-26.md` (the equation ours27 extends).
- `ours26_full_noise_covariance_2026-09-27.md` (the covariance form, the
  determinant cancellation, the per-covariance-group state).
- `ours25b_the_many_reading_readout_2026-09-25.md` (the stacked size).
- `../mass_estimator_lemmas_2026-09-27.md` (the premises; L8 and L13 owe
  restatements for ours27).
- `ours26_sweep_of_the_flat_records_2026-09-27.md` (the rotation failure of
  ours26, the per-component width formula).
- `../unsanctioned/vd_shared_quantity_argument_mass_estimation_2026-09-06.md`
  §5, on correlated channel errors (an `unsanctioned/` record, cited for that
  one statement only).
