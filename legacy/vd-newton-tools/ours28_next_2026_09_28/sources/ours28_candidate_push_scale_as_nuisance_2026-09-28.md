# ours28, a candidate: the push scale is a nuisance of the experiment, not a premise about mass

**2026-09-28 (evening).** Record, filed at the user's request ("Go ahead"), written
by Claude from the conversation of the evening of 28 September. It proposes a
change to ours27 (`ours27_full_covariance_2026-09-27.md`) and records the checks
made on it. It edits nothing and binds nothing. **ours28 is a candidate**, under
the standing rule of this folder: investigate without assuming we understand.

**The question it answers.** The user, 28 September: "We either made a mistake,
or we missed a property of mass... It could also be that the answer is in there
and we just need to stumble on to it." This record's answer: one mistake in the
premise table, and the answer already in two records that had not been put next
to each other.

---

## 0. In one paragraph

L10 ("no push scale", $\omega\to\infty$) is not a symmetry premise. The flat push
weight is a delta at infinite size (`ours27_candidate_size_and_null_readings_2026-09-28.md`
§4), so L10 is the assumption that every push is infinitely large in noise
units, and the width factor $1+M/R^{\ast2}$ (T18) is exactly that assumption's
error. The chain already forces the push weight into the Gaussian family
$e^{-\lVert\mathbf v\rVert^2/2\omega^2}$ (axis splitting and isotropy; note d2
step 7); L10 only picked the endpoint. Keep $\omega$ instead, as one more
nuisance per reading, averaged over like the push itself, under the rule that
already gave the starting law: uniform in the angle the noise defines,
$\tan\varphi=\omega$. That is ours28. Everything downstream is unchanged: the
readout, the stacked length, T7. Tested: the 95% intervals of the law of mass
cover 92–98% in the cells where ours27 covers 66–87%; a zero reading leaves the
law of mass unchanged to $0.0$; the median tracks ours27's and goes to the truth
with $N$. The price: T11's four-number state. What is still open: C3, and the
angle rule, which is now one choice instead of two.

## 1. The mistake in L10

Two records state one fact without meeting.

- **The flat size measure is a delta at infinite size.** Every finite range of
  push sizes gets zero share of the weight (candidate-size record §4, the
  user's "opposite of a Dirac delta").
- **The law's width is wrong by $1+M/R^{\ast2}$** (lemma record T18; the sweep
  §5 gives the per-component form). $R^{\ast2}$ is the total true push in noise
  units. The factor is 1 when $R^{\ast2}=\infty$ and nowhere else.

So L10 asserts $R^{\ast2}=\infty$ and T18 measures how far that is from true.
The many-reading undercoverage is not a defect that needs an added correction.
It is L10 failing at weak signal by its own margin, which is why the sandwich
multiplier was right to refuse (the user, 27 Sept: the fix should "drop out of
the chain rather than be added experimentally").

"No preferred push size" reads as scale invariance, but scale invariance in the
size would be $d^d\mathbf v/\lVert\mathbf v\rVert^d$, the log-flat measure
excluded on 2 September because it collapses zero readings onto the exact
boundary. The flat measure $d^d\mathbf v$ is not scale-free; it transforms by
$\lambda^d$. Neither endpoint is a symmetry. The push measure is not fixed by
symmetry at all.

## 2. What is already forced, and what was a choice

Note d2 step 7 (from `excitation_weighting_premise_2026-09-24.md` and
`axis_splitting_and_the_24_september_change_2026-09-25.md` §4):

- axis splitting (L8) makes the push weight a product over axes;
- isotropy (L9) makes it a function of the push's length;
- a product over axes that depends only on length is, for $d\ge2$, Gaussian,
  $e^{-\lVert\mathbf v\rVert^2/2\omega^2}$ (Herschel–Maxwell), with the flat
  weight as $\omega\to\infty$;
- L10 then "leaves only the flat one".

So L8 and L9 already force a one-parameter family, and the parameter $\omega$
is dimensionless in noise units, so L3 does not touch it (perturbation matrix
§5). $\omega$ is a fact about how hard the experimenter pushed, measured
against the instruments. It is not a property of the object. That is why no
property of mass constrains it: the matrix found the push-scale bend passing
every column. The residual freedom was never about mass. It was a nuisance of
the experiment left unhandled, and the way the chain handles every other
nuisance is step 4: average over it, do not pick it.

## 3. ours28

**Types.** As ours27, with nothing new declared. The push of reading $i$ has an
unknown scale $\omega_i>0$ in noise units.

**Per-reading evidence.** In noise units per axis, with $\mathbf u=(\sin\theta,\cos\theta)$
and $\mathbf z=(x,y)$ the reading's components on that axis, the Gaussian push
of scale $\omega$ gives the reading the marginal covariance $I+\omega^2\mathbf u\mathbf u^\top$,
so (Sherman–Morrison, $\beta=\omega^2/(1+\omega^2)$)

$$
\log L_{ik}(\theta,\beta)=-\tfrac12\log(1+\omega^2)-\tfrac12\lVert\mathbf z\rVert^2+\tfrac12\,\beta\,(\mathbf z\cdot\mathbf u)^2 .
$$

Summed over the reading's axes, with $h_i^2(\theta)=P_i\sin^2\theta+Q_i\cos^2\theta+2D_i\sin\theta\cos\theta$
the reading's own share of ours25b's three sums,

$$
L_i(\theta,\beta_i)\;\propto\;(1-\beta_i)^{d/2}\exp\!\big[\tfrac12\,\beta_i\,h_i^2(\theta)\big].
$$

$\beta_i=1$ is ours27's evidence. $\beta_i<1$ tempers it (the tempering identity,
`../ours26/ours26_learned_push_scale_2026-09-27.md` §2), and the factor
$(1-\beta_i)^{d/2}$ is the price of claiming a small push. **The determinant is
$\theta$-independent for every $\omega$** ($\det(I+\omega^2\mathbf u\mathbf u^\top)=1+\omega^2$),
so nothing attached to a reading depends on mass: T7 holds unchanged.

**The size-angle rule.** $\omega_i=\tan\varphi_i$, $\varphi_i$ uniform on
$(0,\pi/2)$, i.e. $\beta_i=\sin^2\varphi_i$, the arcsine law $\mathrm{Beta}(\tfrac12,\tfrac12)$
(candidate-size record §4). This is the rule that gives the half-Cauchy for
$\theta$: uniform in the angle the noise defines.

**The law of mass.**

$$
p(\theta)\;\propto\;\pi(\theta)\prod_{i=1}^N g_d\!\big(h_i^2(\theta)\big),\qquad
g_d(h^2)=\int_0^{\pi/2}(\cos\varphi)^d\exp\!\big[\tfrac12\sin^2\!\varphi\;h^2\big]\,d\varphi .
$$

$g_d$ is one fixed function per dimension. With a full covariance $\Sigma_i$
the same construction goes through with ours27's $Q_i$ in place of $h_i^2$
(not written out here).

**Candidates at fixed mass and scale.** ours27's conditional Gaussian with the
push prior included: $K_i=\big(B_m^\top\Sigma_i^{-1}B_m+\omega_i^{-2}I\big)^{-1}$,
$\mathbf b_i=K_iB_m^\top\Sigma_i^{-1}\mathbf y_i$. In noise units per axis this
is $v\mid\mathbf z\sim N(\beta\,\mathbf z\cdot\mathbf u,\;\beta)$: the push is
shrunk toward zero by its own $\beta$.

**The readout.** Unchanged in form: $E\lVert\mathbf F^\ast\rVert/E\lVert\mathbf a^\ast\rVert$
over the stacked true vectors, the average now running over $(\theta,\beta_1,\dots,\beta_N)$
and the pushes. Composition, C3 and the stacked length are untouched, since
they act on the law of the sizes, whatever produced it (matrix §2: law bends
pass composition by construction).

**Reductions.** $\beta_i\to1$ for every $i$: ours27. So ours25b, ours26 and
ours27 are all its limits at infinite push.

## 4. Shared or per reading: the zero-reading property decides

The scale could be one $\omega$ for the whole series, or one per reading. Both
calibrate (§5). They differ on three properties.

| | one shared $\omega$ | one $\omega_i$ per reading |
|---|---|---|
| zero reading leaves the law of mass unchanged | **no**: it raises $M$, shifts the $\beta$ weighting toward small pushes, and widens the law of mass (a null says "pushes are small" and discounts every other reading) | **yes**, exactly: the reading's factor is $g_d(0)$, a constant |
| sufficiency (T11) | kept: the joint law is $(1-\beta)^{M/2}e^{\beta H^2(\theta)/2}$, a function of $(P,Q,D,M)$ only | lost: each reading's $(P_i,Q_i,D_i)$ enters on its own |
| axis splitting (L8) | holds as before | holds if the scale sharing rides along with the split: a reading's axes read as 1D readings share its $\omega_i$ (that is what "bookkeeping" means; a per-component scale would break isotropy and rotation) |

Z ("a reading that adds no information from other dimensions must not change
the answer", the user's property; T14/T20's "zero readings leave the law
unchanged") picks per reading. ours28 is the per-reading version.

## 5. Tests

All on one noise pair, $d=3$, Gaussian pushes with per-component SD equal to
the stated SNR, $\pi$ the half-Cauchy at $c=s$, 600 series per cell (Monte
Carlo error about $\pm1$ point). "Coverage" is the 95% central interval of the
law of mass containing the true mass. "Median" is the law's median over the
true mass, cell median.

**Coverage and median.**

| cell | ours27 (flat) | shared $\omega$, integrated | **ours28** (per reading, integrated) | per reading, learned (profile) |
|---|---|---|---|---|
| $m=s$, SNR 1, $N=21$ | 81.5 / 0.99 | 96.2 / 0.99 | **93.3 / 0.99** | 90.8 / 0.98 |
| $m=s/4$, SNR 1, $N=21$ | 86.7 / 1.05 | 98.0 / 1.19 | **96.8 / 1.15** | 94.8 / 1.14 |
| $m=s$, SNR 1, $N=6$ | 85.2 / 1.08 | 98.8 / 1.05 | **97.5 / 1.03** | 95.8 / 1.05 |
| $m=s/4$, SNR 3, $N=21$ | 93.3 / 1.00 | 94.8 / 1.00 | **95.5 / 1.00** | 95.2 / 1.00 |
| $m=s$, SNR ½, $N=21$ | 65.7 / 1.04 | 99.3 / 1.02 | **92.3 / 1.05** | 90.3 / 1.04 |

- The width failure goes. Strong cells are untouched.
- The pull toward the starting law at weak signal ($m=s/4$, SNR 1: 1.05 →
  1.15) is smaller than the record's §7 rule gave (1.36 → 2.43), because that
  rule used one scale per series; per reading, the shrinkage is local.
- The learned (profile) version under-covers by a few points, as profiling a
  nuisance usually does. The chain's own rule is to average (step 4).

**Zero reading.** Two 3D readings; a 3D zero reading appended: maximum change
in the normalised law of mass $0.0$ (exactly, not to rounding: the appended
factor is a constant).

**Consistency.** Median over true mass, mean of six runs, ours27 and ours28 on
the same data (grid floor about 0.5%):

| world | cell | $N=50$ | $N=400$ | $N=1600$ |
|---|---|---|---|---|
| Gaussian pushes | $m=s$, SNR 1 | 1.018 / 1.095 | 1.000 / 0.992 | 1.020 / 1.011 |
| Gaussian pushes | $m=s/4$, SNR 1 | 0.767 / 0.825 | 1.007 / 1.021 | 0.950 / 0.940 |
| fixed-length pushes | $m=s$, SNR ½ | 2.958 / 1.880 | 1.230 / 1.232 | 0.944 / 0.959 |
| fixed-length pushes | $m=5s$, SNR 1.5 | 1.060 / 1.133 | 0.950 / 0.947 | 0.982 / 0.985 |

ours28 tracks ours27 to within the runs' own scatter in both push worlds and
goes to the truth with $N$. This is a look, not the matrix's quadrature test;
the latter is owed (§8). Note that the fixed-length world is the one where the
Gaussian push model is wrong, and ours28 is no worse there.

**Not tested.** The readout point under ours28 (only the law was run); the
swap, units and axis-splitting columns numerically (they hold by the algebra
of §3, but the matrix should say so); per-axis and full-covariance noise; the
interval of the point rather than of the law.

## 6. The price

T11 goes. Under ours27 the data enter through $(P,Q,D,M)$; under ours28 they
enter through each reading's $(P_i,Q_i,D_i)$, because each reading's $\beta_i$
is averaged separately. Combination is still exact, as a product of per-reading
factors, so persistence (G1) and T12's additivity survive in the form "the
log-law of a union is the sum of the log-laws". What is carried is the log-law
as a function of $\theta$ (a curve), not four numbers. The four-number state
returns only with a shared $\omega$, and that costs Z (§4). This record takes
Z.

The combination record's "adding a reading means carrying four numbers"
becomes "adding a reading means adding one term $\log g_d(h_i^2(\theta))$ to
a curve". The curve is determined by the readings' triples, so the state can
also be kept as the list of triples.

## 7. What is still open

- **C3.** Untouched by this record.
- **The angle rule.** Before this record there were two open measures, $\pi$'s
  shape and the size measure (candidate-size record §2: "one measure on
  candidate space, cut along a coordinate line"). ours28 applies one rule to
  both: uniform in the angle the noise defines. That leaves one choice, not
  two. It is still a choice.
- **A check that could turn the choice into a derivation, or break it.** The
  starting-law record §6 gives the Jeffreys route for $\theta$: with the push
  held fixed, Jeffreys is uniform in $\theta$, the half-Cauchy. Run the same
  route for $\omega$: for one reading $y\sim N(0,1+\omega^2)$ per component the
  Fisher information is $2\omega^2/(1+\omega^2)^2$ per component, so the
  Jeffreys density in $\varphi$ is $\propto\tan\varphi$, **not uniform**
  (Claude, 28 Sept, algebra only, unchecked). If that stands, "uniform in the
  noise-defined angle" and "Jeffreys" agree for the mass and disagree for the
  size, and the disagreement will say which principle is actually doing the
  work in step 8. If they must agree, the candidate is Jeffreys-in-$\varphi$
  for the size, and it should be run through §5.
- **Which $\omega$ is measured against.** Per reading in noise units, one scale
  per reading shared across its axes. With several noise pairs the units fork
  (`../ours26/README.md`) reappears here: whether $\omega_i$ is one number
  across axes of different noise ratios. Not examined.

## 8. What this changes in the lemma record (proposed, not made)

- **L10:** struck. Replaced by "the push scale is a per-reading nuisance,
  averaged under the size-angle rule". Its old content is the assumption
  $R^{\ast2}=\infty$.
- **L9 (isotropy)** and **L8 (axis splitting):** unchanged, now read as forcing
  the Gaussian *family* (they did already; d2 step 7).
- **T6 (flat pushes):** becomes the $\omega\to\infty$ limit, not the law.
- **T18 (width):** becomes the measure of L10's error, and is absent under
  ours28 to Monte Carlo precision (§5).
- **T11 (sufficiency):** marked [flat push], as T14 and the others are marked
  [one pair].
- **T7:** unchanged; strengthened, since the determinant is $\theta$-free for
  every $\omega$.
- **§3.3 (coverage, $\omega$ and the units of isotropy):** the coverage part
  closes if ours28 stands; the units fork stays.
- **§3.1 (the starting law):** the shape question and the size measure are one
  question with one current rule; the Jeffreys check of §7 is owed.

## 9. Owed

- the tester build of ours28 (`ratio_of_means_rss` with the per-reading
  $g_d$ law; a preset for the five cells of §5, with the tester's own push
  worlds);
- the perturbation matrix rerun on ours28: every column, plus a bend in the
  size-angle rule;
- the Jeffreys-in-$\varphi$ check (§7), and its coverage if it differs;
- the readout point under ours28, and the point's interval;
- the sweep's coverage cells (per-axis ratios 0.2/0.5/1.25 and 0.1/1/3), which
  need the several-pair form of $\omega_i$.

## 10. Reproduction

Scripts written in Claude's workspace on 28 September and not kept; numpy and
scipy only, a few minutes each. The law on a $\theta$ grid (400 or 800 points,
uniform in $\theta$, which is the half-Cauchy weight), $\beta$ on a $\varphi$
grid (80 or 160 points, uniform, which is the arcsine weight); the shared and
per-reading laws as written in §3 and §4; coverage from the cumulative law;
the profile version with $\beta_i=\max(0,1-d/h_i^2)$. The zero-reading check
compares normalised laws with and without an appended zero reading. The
consistency runs use six seeds per cell and a grid of 800.

## 11. How we got here

1. The user asked whether "the estimator tends to $F/a$ as uncertainties reduce
   and the pair is aligned" had ever been stated. It had not, beyond the
   1 August slogan. Claude checked it numerically on ours25b (aligned readings:
   yes, along every path; exact 0/0: path-dependent, every value reachable;
   one side zero or anti-aligned: diverges; acute misalignment: a finite limit
   the premises call impossible). It is a check every candidate passes, so it
   forces nothing; not filed separately.
2. The user: "We either made a mistake, or we missed a property of mass...
   It could also be that the answer is in there." Claude's first answer named
   L10 as the mistake (§1) and was told off for also calling a wall a live
   option. Dropped.
3. The shared-$\omega$ law was run first: coverage fixed, four-number state
   kept. Then the zero-reading tension was noticed (§4), and the per-reading
   version run, which passes it exactly. A broadcasting bug (`np.outer` on a
   2D array) produced a first per-reading table of zeros; fixed and rerun.
4. The user: "Go ahead."

