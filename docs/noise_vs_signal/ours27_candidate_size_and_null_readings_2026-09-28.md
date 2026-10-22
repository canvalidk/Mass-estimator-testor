# Candidate size and null readings: seven short reminders

**2026-09-28.** Record, filed at the user's request ("Will you make a record of
this last ideas… Keep em brief. These are for both of us"). Written by Claude
from the conversation of the night of 27–28 September and the morning of the
28th. Each entry gives the idea, what it showed, and what it rules out or leaves
open. The numbers come from standalone scripts in Claude's workspace, which are
not kept. The equation throughout is ours27
(`ours27_full_covariance_2026-09-27.md`). It edits nothing and binds nothing.

**The standing rule these were worked under** (the user, 27 Sept): investigate
without assuming we already understand. Several candidate fixes were proposed
confidently that night and failed, so treat any "solution" below as a
candidate.

**Terms.**

- *Candidate*: one hypothetical aligned pair, $\mathbf F^\ast=m\,\mathbf a^\ast$,
  that the definition averages over. Neither slot is a true value that gets
  recovered.
- *Candidate size*: the candidate's length in noise units, $r$. In the records
  before 28 Sept it was called the "push" (earlier, the "excitation").
- *Size measure*: the weighting the definition puts on candidate size. ours27's
  is flat ($\omega=\infty$).

---

## 1. Has the point estimate changed, and what did rotation fix?

- **It has changed at every version:**
  - 24 Sept, law in 1D and 3D: a 5/0 reading went 6.2666 → 6.1370 in 1D;
  - 25 Sept, many-reading readout: median gap 0.56%, worst 4.8%;
  - ours26 and ours27: no change on old data.
- **ours27 made both the point and the law rotation-invariant.** Data and
  covariance rotated together: readout 1.85529234 at every angle, median
  1.878668. ours26 on the same data: 1.855, 2.248, 2.340.
- **The size-measure question is mostly about widths, not the point.** A finite
  scale leaves the law's peak where it is and moved the point by −0.01% (one
  pair) and −0.5% (per-axis) at the scale tested.

*Rules out:* the worry that ours27's rotation fix only touched the
uncertainties.

## 2. The starting law and the size measure as overlapping problems

- **Idea (the user's OP technique):** two questions that each block the other.
  Find a case where both have the same answer.
- **Shown:**
  - For one reading, a candidate is a point in the plane (candidate force,
    candidate acceleration), in noise units: $r(\sin\theta,\cos\theta)$. The
    angle carries the mass; $r$ is the size.
  - The starting law is the measure on the angle; the size measure is the
    measure on $r$. They are one measure on candidate space, cut along a
    coordinate line.
  - Each blocks the other:
    - the Jeffreys route to the starting law depends on what is held fixed
      about size (`../ours26/ours26_starting_law_2026-09-27.md` §6);
    - which size measure counts as "flat" depends on which coordinate the mass
      part is written in. ours21 is flat in the whole plane, which splits as
      uniform angle × $\lvert v\rvert$.
- **Hinge candidates, none tested:**
  - the starting law as an unseen reading (T16);
  - lumping at the level of the law;
  - the complex log of the candidate, which puts the angle and log-size on
    equal footing.

*Leaves open:* the hinges. The overlap was not checked with many readings, where
the angle is shared but the sizes are per reading.

## 3. Repeated 0/0 readings

- **The law stays exactly the starting law** for any number of nulls, any noise
  structure.
- **The point:**
  - with one noise pair: exactly $s$ for every number of nulls and every
    swap-symmetric starting-law shape;
  - per-axis (ratios 5, 2, 0.8, centre 2): it drifts with the number of nulls,
    and the size of the drift depends on the starting law's shape.

  | starting law | 1 null | 5 nulls | ∞ nulls |
  |---|---|---|---|
  | half-Cauchy | 1.759 | 1.730 | 1.719 |
  | log-normal, spread 1 | 1.830 | 1.807 | 1.799 |
  | log-normal, spread 0.3 | 1.975 | 1.972 | 1.970 |

- **The cause is in the readout, not the size measure.** At nulls any common
  size scale multiplies every component by the same $\sqrt\beta$ and cancels
  (exact algebra). The drift comes from the stacked length averaging
  $\sqrt{\sum_c\sigma_{a,c}^2\cos^2\theta_c Z_c^2}$ over components whose
  weights depend on the mass differently.
- **With no readings at all, the readout does not exist** (an empty stack is
  0/0). Only the law exists.

*Rules out:* fixing the null drift through the size measure. See §5 for why the
drift is not a flaw.

## 4. "How can 50 × 0/0 gain confidence?" The delta at infinity

- **The law does not gain confidence.** The 95% range stays
  [0.039, 25.45] × centre for any number of nulls. Only the point moves.
- **The flat size measure has no normalisable form.** It is the limit of wider
  and wider bells, and every finite range of sizes loses its share of the
  weight. It is a delta at infinite size: the user's "opposite of a Dirac
  delta". The opposite extreme, a delta at zero size ($\omega\to0$), makes every
  reading carry nothing.
- **The asymmetry between the two questions:**
  - the mass lives on a bounded angle (0–90°), so a uniform "know nothing" law
    exists: the half-Cauchy;
  - size is unbounded, so no uniform law exists.
- **Candidate from this (Claude's, tested in §7):** give the size scale the
  same rule as the mass, uniform in the angle $\tan\varphi=\omega$, which is
  the arcsine law Beta(½, ½) on $\beta=\omega^2/(1+\omega^2)$. One scale per
  series.

*Keeps in mind:* a flat size measure makes every null reading an infinitely
large surprise about size. That surprise cancels in the law and not in the
readout.

## 5. Rotate a non-zero 3D reading to 0/0 on one or two axes

- **Shown:**
  - Any two vectors lie in a plane, so every 3D reading is, by rotation, a 2D
    reading plus a 0/0 axis. A parallel reading is a 1D reading plus two.
  - ours27 gives a generic 3D reading exactly its embedded answer: 1.37849
    against a 2D answer of 1.38392. The law is identical across dimensions.
  - So "does a 0/0 axis count?" and "does the dimension count?" are one
    question.
- **Resolved: known-empty and measured-empty are different inputs.**
  - Recording in 1D or 2D declares the other components known to be zero, as
    with a constraint.
  - A 0/0 axis says they were measured as zero, with noise.
  - ours27 joins the two continuously as that noise shrinks. Parallel reading,
    1D answer 1.16161:

    | noise on the empty axes | readout |
    |---|---|
    | 0.001 | 1.16161 |
    | 0.1 | 1.15984 |
    | 1 | 1.12926 |
    | 10 | 1.10405 |
    | 1000 | 1.10319 (saturated) |

  - Non-parallel, 2D answer 1.38392: 1.38392 → 1.37849 → 1.34627.
- **The property as it holds:** an axis *known* to be empty does not change the
  answer. An axis *measured* as empty does not change the law, and changes the
  point by an amount that vanishes as the measurement becomes exact.
- **Why the point moves:** a measured null is information about candidate
  sizes, not about the mass. The point is a ratio of mean sizes (C3), and
  composition forces that. A point that heard only the law of mass would be a
  summary of that law, like the median, which fails lumping.

*Rules out:*

- `ours25b` §7's "d is physics". The dimension declares what is known.
- The null drift (§3) and the useless-axis shift
  (`ours27_full_covariance_2026-09-27.md` §4) as flaws.
- A readout that ignores 0/0 axes. Under rotation invariance it would jump when
  force and acceleration become exactly parallel.

*Caveat:* the clean version needs one noise pair. Per-axis, rotating makes the
hidden axis's noise correlated with the others.

## 6. Misaligned vectors with shrinking uncertainty

- **What ours27 does.** One pair, noise shrinking:
  - **10° apart:** the misfit grows (Q = 242 at noise 0.01, odds about
    10⁻⁵³), yet the law *narrows* around 2.018: [1.97, 2.06].
  - **Perpendicular:** the median runs to about 1/σ (257 at 0.01).
  - **Anti-parallel:** the median runs to 28,000.

  The law is normalised, so only relative fit between masses counts. The
  absolute badness of the best fit is discarded. That was the old gate's job,
  removed on 7 Sept.
- **The user's ruling (28 Sept): the σ's are correct.** Under the premises
  (correct σ, Newton II), such data essentially never occur. When they do, the
  equation's output is the exact conditional: the least-impossible explanation.
  For perpendicular readings that is "the smaller reading is all error", so the
  mass is at least about 1/σ, with the starting law's tail beyond.

*Rules out:* treating these cases as the equation responding incorrectly, and
inflating σ from the misfit (the Birge-ratio route), under the ruling.
*Points to:* the real failure on data the premises actually produce, the
many-reading width (§7).

## 7. Lower stakes for unrepresentative readings: is that ω?

- **No; they are two different levers.**
  - **Per-reading down-weighting** (a reading that fits badly counts less) means
    a heavy-tailed noise law. With the flat size measure it is inconsistent
    even in its own world: a t₅ model converges to 5.64 against a true 5
    (`ours27_perturbation_matrix_2026-09-27.md` §4). The Gaussian ruling excludes
    it anyway.
  - **ω** (a finite size scale, one per series) raises *every* reading's
    evidence to the same power β < 1. The discount is uniform, set by the whole
    series, and never picks out one reading.
- **Why it matters.** The width failure happens on typical data with no
  outliers, so per-reading down-weighting would not touch it.
- **The §4 rule tested on the width cells** (one pair, 3D, 600 series per cell;
  fixed-length and Gaussian push worlds give nearly the same numbers):

  | cell | flat coverage | new rule coverage |
  |---|---|---|
  | $m=s$, SNR 1, N = 21 | 79.3% | 97.5% |
  | $m=s/4$, SNR 1, N = 21 | 79.0% | 98.3% |
  | $m=s$, SNR 1, N = 6 | 85.0% | 99.7% |
  | $m=s/4$, SNR 3, N = 21 | 92.2% | 96.2% |

  - Across the 18 cells: 94–99.8% against flat's 79–96%. Strong cells are
    barely changed.
  - **Costs:**
    - over-coverage in the weakest cells (99.7–99.8%);
    - a stronger pull toward the starting law for light masses at weak signal.
      Median/truth at $m=s/4$, SNR 1: 1.36 → 2.43 (N = 21), 1.92 → 3.24
      (N = 6). At SNR 3 the two agree.
  - Not tested per-axis. The sweep showed one scale per series cannot calibrate
    several noise pairs, so it would need one angle per pair.

*Rules out:*

- per-reading robust weighting as the width fix;
- the fixed-length tester pushes as the cause of the undercoverage.

*Leaves open:* the trade between calibrated widths and a less prior-led point
at weak signal. The properties we have do not decide it. It is also measured
against a planted mass, a yardstick not yet settled under "mass is its
estimator".

---

## Status in one line

Off the list: the null drift and the dimension question (§5), and impossible
data (§6). Still open: the starting law's shape, C3, and the size measure (the
widths, §7).
