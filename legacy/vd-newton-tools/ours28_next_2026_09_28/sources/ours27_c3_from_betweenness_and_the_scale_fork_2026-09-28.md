# C3 from betweenness, the size measure, and a fork with ours28

**2026-09-28 (afternoon to evening).** Record, filed at the user's request
("Yes go ahead"), written by Claude from a session that ran in parallel with
the one that produced the three ours28 records in this folder. It edits
nothing and binds nothing. The equation throughout is ours27
(`ours27_full_covariance_2026-09-27.md`), one noise pair unless stated.

**Rulings it uses (the user, 28 Sept).**

- Betweenness, provisionally: "i think we can accept betweenness
  provisionally."
- Each reading keeps its own frame and sign convention: "Yes of course.
  Bookkeeping can't matter."
- Same candidate law, same mass, whichever instruments produced it: "Sure."
- The σ's are correct (morning ruling, recorded in
  `ours27_candidate_size_and_null_readings_2026-09-28.md` §6).

**In one paragraph.** C3 can be derived instead of assumed: betweenness makes
the readout a fixed per-candidate score, the bookkeeping symmetries make the
score depend only on $(X,Y)$, and composition leaves the ratio of means. The
derivation's second step needs the **grouping of components into readings to
be bookkeeping too**. The same premise, with properness and "readings not yet
taken can't change the measure over readings already taken", forces the size
measure to be one Gaussian scale per series (Kingman 1972). But one proper
scale per series fails the zero-reading property once other readings are
present, which is why the ours28 record chose one scale per reading. So there
is a fork. Grouping as bookkeeping, the zero-reading property, and a proper
size measure cannot all hold. ours27 drops properness, a shared proper scale
drops the zero-reading property, and ours28 drops grouping as bookkeeping,
which reopens a narrow part of C3. The ruling that settles it is the user's:
**is which components were measured together bookkeeping, or physics?**

---

## 1. C3 derived, given grouping as bookkeeping

### 1.1 The three steps

**Step 1: betweenness gives a fixed score.** Betweenness: if two candidate
laws give readouts $r_1\ne r_2$, their mixture gives a readout strictly
between; if they give the same $r$, so does the mixture. Take two single
candidates, one heavier than $r$ and one lighter. Exactly one mixing ratio
gives readout $r$: an exchange rate between them at level $r$. Betweenness
makes the set of laws with readout $r$ a flat slice (a hyperplane in the
space of laws), so the exchange rates are consistent and define a score
$\phi(c,r)$ with

$$
R(P)=r\iff E_P[\phi(c,r)]=0 .
$$

The exchange rates come from pairs of single hypothetical candidates, so the
score is fixed before any data. That is why "declared before data" cannot
exclude the 6 Sept direction weight: its rule is declared, it just reads the
law. Betweenness is what excludes it (§1.4). Examples: ratio of means,
$\phi=X-rY$; median, $\phi=\pm\tfrac12$ as $m\gtrless r$. The representation
is Dekel's theorem (1986), proved there for prizes on a compact interval with
a monotonicity axiom. Candidate space is neither, so this step is borrowed:
the flat-slice argument looks general, but it is not written out here.

**Step 2: bookkeeping makes the score depend only on $(X,Y)$.** The score
must be unchanged by every move that does not change the physics. For one
noise pair, in noise units, the candidate is a mass angle $\theta$ and a
stacked vector $V\in\mathbb R^{Nd}$. The moves generate this many independent
rotations of $V$ (for 3 readings in 3D; all rotations would be 36):

| moves allowed | rotations generated |
|---|---|
| one shared frame for all readings | 3 |
| + regrouping axes across readings (x of reading 1 with y of reading 2 ...) | 18 |
| + each 1D reading keeps its own sign convention | **36 (all)** |

Checked for up to 5 readings in 2D and 3D. With each reading in its own frame
(also ruled) it holds for every $N$ by an argument: any two axes of different
labels can be put in one reading and rotated, and rotations in neighbouring
planes generate the rest, including the plane of two same-label axes. With all rotations, the score sees
only $\theta$ and $\lVert V\rVert$, i.e. $(X,Y)$.

**Step 3: composition.** Among scores of $(X,Y)$, lumping and common force
leave $X-rY$, the ratio of means (the 6 Sept composition theorem, used as it
stands). The median dies here. A readout weight that is any function of
$(X,Y)$ breaks lumping (residual $3.5\times10^{-2}$ for $1+0.1\lVert V\rVert^2$,
$1.3\times10^{-2}$ for $1+0.1m$; $0$ for a constant).

**What survives without the last row of the table.** With one shared frame
and regrouping but no own sign conventions, the resultant $\lVert\sum_iv_i\rVert$
survives: how aligned the candidates of different readings are. Own sign
conventions kill it.

**The pure 1D case.** A 1D reading is a 3D reading whose other components are
known to be exactly zero, so its candidates are 3D candidates on a subspace
and inherit the 3D score. This is a statement about candidates, not about any
law, so it holds whatever the size measure.

### 1.2 What step 2 needs that the table hides

The middle row, regrouping axes across readings, is only a symmetry if
**which components were measured together is bookkeeping**. Under ours27's
flat size measure it is: the law cannot tell. Under ours28 it is not: a
reading's axes share one push scale, so a composite reading made from x of
reading 1 and y, z of reading 2 is not isotropic. Checked: with a scale per
reading, rotating such a composite reading changes its kurtosis 8.29 → 6.98
→ 5.63 at 0°, 22.5°, 45°; with one scale for the series it stays 8.29, 8.29,
8.25 (exact 8.155; 4 million draws, heavy tails).

Without cross-reading regrouping, the moves generate only each reading's own
rotations, and the score can also depend on the profile of per-reading sizes.
Example: $h=\sum_i\lVert v_i\rVert^4/(\sum_i\lVert v_i\rVert^2)^2$ built from
the per-reading acceleration sizes. It is degree-0 homogeneous, so a split of
mass or of acceleration that is common to a candidate's readings leaves it
unchanged. So it passes lumping, common force, the swap, units, calibration,
consistency, betweenness (it is linear) and same-law-same-mass. Under ours28,
C3 is needed to exclude it and nothing else does (algebra; not simulated).

### 1.3 Per-axis noise: "same candidate law, same mass" is needed

With per-axis noise the declared covariance has directions, so step 2 fails
unless the score cannot see the instrument. The survivor, data-free:

$$
h=1+\varepsilon\,\tfrac12\Big(\frac{\mathbf u^\top S_a\mathbf u}{\operatorname{tr}S_a}+\frac{\mathbf u^\top S_F\mathbf u}{\operatorname{tr}S_F}\Big),
$$

with $\mathbf u$ the stacked candidate direction and $S_a,S_F$ the stacked
declared noise blocks. Two 3D readings, $\sigma_F=1$,
$\sigma_a=(0.2,0.5,1.25)$, $c=2$:

- swap: $R\times R_{\rm swapped}=1.000000$; units (force ×3): 3.000000;
- each reading rotated in its own frame with its covariance:
  $-0.0002\pm0.0003$ (8 seeds);
- axis splitting (2 × 3D against 6 × 1D): identical to 8 digits;
- lumping, common force, calibration, consistency: exact by algebra (the
  weight depends only on direction);
- betweenness: 0 violations in 20,000 random mixtures;
- one isotropic pair: $g$ is constant ($1/Nd$), the weight vanishes exactly;
- size: $-0.52\%$ at $\varepsilon=1$, $-2.3\%$ at $\varepsilon=10$; with
  $\sigma_a=(0.05,0.5,3)$, $-1.0\%$ and $-4.5\%$. $\varepsilon$ is unbounded.

"Same candidate law, same mass, whichever instruments produced it" excludes
it, and step 2 then runs in physical units. ours27 obeys it. The starting law
may still depend on the instruments (the arc-length candidate does), because
it enters through the law.

### 1.4 The 6 Sept weight fails betweenness

$h=1+(\mathbf u\cdot\mathbf e)^2$, $\mathbf e=E_P[\mathbf u]$:

- random laws: 683 violations in 20,000 mixtures; worst: readouts 0.858 and
  0.851, mixture 1.033;
- laws ours27 produces (2D, two readings, the starting law a mixture of two
  declared starting laws): 7 in 7200; the worst re-run with 5 × 20,000 fresh
  candidates: readouts 2.8376 and 2.8338, mixture 2.8481, outside by
  $0.366\%\pm0.014\%$; the second $0.169\%\pm0.006\%$; the third was sampling
  noise;
- ratio of means and median: 0 violations anywhere.

## 2. Self-consistency: what coverage means without a planted mass

Coverage against a planted mass presupposes a mass behind the readings. The
non-realist version: feed the definition data generated by its own measure,
and its intervals should be calibrated on average. Any proper measure passes
this automatically, so it cannot choose among proper measures; it can only
exclude improper ones. One pair, 3D, 4000 series per row, nominal 95%, in the
world of the uniform-angle rule with one shared scale (mass angle and
$\varphi=\arctan\omega$ uniform):

| readings | uniform-angle rule | flat ours27, same data |
|---|---|---|
| 1 | 95.0 | 92.4 |
| 2 | 94.7 | 90.3 |
| 6 | 95.2 | 85.9 |
| 21 | 94.9 | 79.8 |

- **Flat is calibrated only at infinite scale.** At fixed $\omega$ = 0.5, 1,
  2, 5, 20, 100 (21 readings, 1500 series each): 70.7, 85.5, 92.6, 93.9,
  94.9, 94.9, matching T18 at large $\omega$.
- **At 21 readings flat's failure tracks the scale, not the mass:** 58.6% at
  $\omega<0.41$ up to 93.8% at $\omega>2.4$; 78–81% in every mass bin.
- **The planted-mass yardstick condemns every starting law.** One reading,
  masses below $0.27s$: flat 81%, the rule 86%; central masses 97–99%. Any
  proper starting law does this. So "cover every planted mass" cannot be
  what width means, and "calibrated in its own world" can.

## 3. The size measure

### 3.1 Forced, given grouping as bookkeeping

- **One scale per series.** The regrouping symmetry of §1.2 applied to the
  measure.
- **Gaussian given the scale.** Kingman (1972): an infinite sequence whose
  every finite section is spherically symmetric is, given one random
  variance, independent $N(0,\omega^2)$. The premise needed: readings not yet
  taken cannot change the measure over readings already taken. With
  properness, the measure is $V\mid\omega\sim N(0,\omega^2I)$, $\omega\sim G$.
  Flat is $G$ with all weight at $\omega=\infty$, outside the theorem only
  because it is improper.
- **The scale cannot depend on the mass.** A zero reading with no other data
  must leave the law of mass unchanged; its factor
  $\int G(d\omega\mid\theta)(1+\omega^2)^{-M/2}$ for every $M$ pins $G$ at
  each $\theta$, so $G$ cannot depend on $\theta$. A scale tilted with mass
  moves the median from 1 to 0.49 after one 3D zero reading. ours28's
  per-reading rule satisfies this too.

### 3.2 The fork

The ours28 record §4 found that one shared scale fails the zero-reading
property once other readings are present. Confirmed here, two 3D readings
plus appended 3D zero readings, uniform-angle scale law ("changes by" is the
largest pointwise change in the law relative to its peak):

| scale | +1 zero reading | +5 zero readings |
|---|---|---|
| one per series | law changes by 4.9%; median 0.9893 → 0.9930 | 10.4%; median → 0.9971; 95% log-width 6.741 → 6.547 |
| one per reading (ours28) | $2\times10^{-15}$ | $7\times10^{-15}$ |

In general, with a shared scale the zero-reading property holds for all data
only if $G$ is a point mass: $\omega$ fixed and declared, or flat. So three
things cannot all hold:

| | grouping is bookkeeping | zero reading leaves the law unchanged | proper size measure |
|---|---|---|---|
| ours27 (flat) | yes | yes | **no** |
| one shared proper scale | yes | **no** | yes |
| ours28 (one scale per reading) | **no** | yes | yes |
| a fixed declared $\omega_0$ | yes | yes | yes, but calibrated only at $\omega=\omega_0$ |

What rides on the first column: step 2 of §1.1 across readings (so C3's
per-reading size-profile hole, §1.2), and the one-scale-per-series half of
§3.1. The second column is the user's property Z. The third is self-consistency
(§2) and the widths.

### 3.3 Which scale law: Jeffreys and the others

This session computed the same Fisher matrix as the Jeffreys-check record
(`ours28_size_measure_jeffreys_check_2026-09-28.md` §1) and agrees with it:
$I_{\theta\theta}=\omega^4/(1+\omega^2)$, $I_{\omega\omega}=2\omega^2/(1+\omega^2)^2$,
$I_{\theta\omega}=0$, checked by finite differences to six digits at three
values of $\theta$ each for $\omega$ = 0.3, 1, 3. Added here:

- **A third route to uniform $\theta$.** First-order probability matching
  for $\theta$ with an orthogonal nuisance gives $\pi\propto I_{\theta\theta}^{1/2}g(\omega)$
  with $g$ arbitrary (Tibshirani 1989; Staicu and Reid 2008). $I_{\theta\theta}$
  does not depend on $\theta$, so the $\theta$-part is uniform: the
  half-Cauchy. Jeffreys, reference and matching agree on the mass.
- **Why the scale has no canonical law.** The only linear maps of the
  whitened plane that preserve the model family
  $N(0,I+\omega^2\mathbf u\mathbf u^\top)$ are rotations and reflections,
  which move $\theta$ and never $\omega$: the declared noise fixes the scale.
  So there is no group-invariant measure for $\omega$. Every invariance-based
  rule is improper at $\omega\to\infty$, in $\beta=\omega^2/(1+\omega^2)$:
  joint Jeffreys $\mathrm{Beta}(2,-\tfrac12)$, reference $\mathrm{Beta}(1,0)$,
  matching with $g=1$ $\mathrm{Beta}(\tfrac32,-1)$; uniform angle
  $\mathrm{Beta}(\tfrac12,\tfrac12)$ is proper. The Jeffreys-check record's
  L4′ ($\omega\leftrightarrow1/\omega$) supplies by premise what no group
  supplies.
- **The sharing changes the verdict on Jeffreys.** With one shared scale
  (this session, 1500 series per cell, 3D) joint Jeffreys covers
  92.6–95.5% and the uniform angle 94–96%. With one scale per reading (the
  Jeffreys-check record) joint Jeffreys covers 73–94%. A shared scale pools
  all $M$ components, so an improper tail costs little. Pull toward the
  starting law at $m=s/4$ with a shared scale, median over truth: $\omega=0.5$,
  21 readings: flat 1.64, Jeffreys 2.25, uniform angle 2.88; $\omega=1$,
  6 readings: 1.44, 1.73, 2.12.

## 4. Proposed for the lemma record (not made)

- **L6 (C3):** a theorem given betweenness, own frames and sign conventions,
  grouping as bookkeeping, same law same mass, and composition. If grouping
  is physics (ours28), C3 shrinks to "no weight on the per-reading size
  profile".
- **New premises, ruled 28 Sept:** betweenness (provisional); own frames and
  sign conventions; same candidate law, same mass.
- **Owed ruling:** grouping, bookkeeping or physics (§3.2).
- **L10:** as in the ours28 record, not a symmetry. Under grouping as
  bookkeeping, replaced by "one Gaussian scale per series with a proper law
  independent of mass", at the price of Z.

## 5. Not done

- The flat-slice step (§1.1, step 1) written out for our candidate space,
  with the continuity it needs (ratio of means is continuous where means
  converge; general readouts not checked) and the strict form of betweenness.
- Per-axis noise for the size measure: the symmetry of §1.1 acts only within
  a noise pair, which matches the known-empty record's "one scale per noise
  pair" candidate; the joint law of the pairs' scales is open.
- The size-profile weight of §1.2 simulated.
- The self-consistency test with per-reading scales.

## 6. Reproduction

Scripts in Claude's workspace on 28 September, not kept; numpy and scipy
only.

- Step 2: Lie-algebra closure of the generated rotations (numerical rank of
  iterated commutators) for $N\le5$, $d=2,3$.
- §1.3: ours27 candidates drawn from the conditional Gaussians on 401 nodes in
  $\log m$, 3000 per node, common random numbers.
- §1.4: exact readouts on random discrete laws; ours27 laws on 241 nodes with
  1500 candidates per node, worst cases re-run with 20,000.
- §2 and §3.3: the law on a 3001-point $\theta$ grid (refined around the
  mass), the scale on a 721-point $\log\omega$ grid.
- §3.2: 4001-point $\theta$ grid, 801-point $\varphi$ grid.

## 7. How we got here

1. The user asked to be caught up ("we are getting really close to resolving
   this thing"). Claude proposed that the three open items are one measure on
   candidate space and offered two leads: C3 with many readings, and
   self-consistency coverage.
2. Lead 1 found the table of §1.1 and the per-axis survivor of §1.3.
   Claude's first framing, that "declared before data" kills the 6 Sept
   weight, was wrong and was replaced by betweenness. The user accepted
   betweenness provisionally, then the sign-convention premise, then "same
   candidate law, same mass".
3. The user asked whether the candidates occupy a volume rather than a slice.
   They do: the slice is in the space of laws. With per-axis noise the
   candidate cloud leaves the plane of the two readings (4.9% of its length
   in one example; 0 with one pair).
4. The user: "Keep going." The size-measure moves of §3.1 and the Jeffreys
   computation followed. On filing, this session found the three ours28
   records written in parallel the same evening, confirmed their zero-reading
   finding, and wrote §3.2 as the fork between the two lines of work.

