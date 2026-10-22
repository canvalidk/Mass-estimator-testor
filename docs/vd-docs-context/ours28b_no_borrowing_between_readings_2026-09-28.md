# No borrowing between readings: ours28b withdrawn, and what per-reading rules cannot do

**2026-09-28 (late night).** Record, filed at the user's request ("Yes
please"), written by Claude in the session that constructed ours28b. It edits
nothing and binds nothing. It closes the ours28b line under a ruling the user
made tonight, and records two results about the rules that survive it.

## 1. The ruling

The user, 28 September, in his words:

> "in general we actually don't want to learn by assuming anything about the
> upcoming new measurement. So we don't want to claim that the upcoming push
> will be estimatable."

> "info can come from any measurement domain in any combination. Like one F = 1,
> then F = 50000000, then F = 0.00001. So if you learn for example that there
> is some consistent direction or measurement uncertainty between guesses. We
> don't want that."

Stated as a premise: **the definition borrows nothing between readings.**
Readings may come from any regime in any combination. Nothing inferred from
one reading about pushes (their size, their direction, their spread) may be
used on another. Readings are joined only through the mass they share, which
is what Newton II says they have in common. What is known about a reading is
declared for that reading (its $\Sigma_i$, its candidate subspace), never
learned from the others.

## 2. What it excludes and what it keeps

**Excluded:**

- **ours28b** (`ours28b_one_scale_per_noise_pair_2026-09-28.md`): one scale
  per noise pair, shared over the series.
- **A learned population of per-reading scales** (tested tonight, called
  "middle"): each reading's scale drawn from a population whose parameters
  are learned from the series.
- **A shared direction tensor** (tested tonight): Codex's tensor, with a
  proper symmetric matrix-Beta law, shared over the series.
- Any one scale per series.

**Kept:**

- **ours27** (flat).
- **ours28** in the full-covariance form of
  `../ours27/ours28_full_covariance_form_and_declared_subspace_2026-09-28.md`:
  one scale per reading from a declared law, and the declared candidate
  subspace.

The ruling and the zero-reading property turn out to be one principle seen
from two sides. If nothing is shared, an exact zero reading cannot affect
anything, so the strong zero-reading property returns.

**What the tests of tonight showed** (one noise pair unless stated; 95%
coverage of the law; 400–500 series per cell; true mass at the starting law's
centre):

| world | flat | ours28 | middle | ours28b | shared tensor |
|---|---|---|---|---|---|
| one pair, every push $r=1$, 200 readings | 68.4 | 84.0 | 96.2 | 96.4 | – |
| one pair, pushes alternating $r=3$, $r=0.5$ | 90.4 | 94.8 | 96.0 | 95.8 | – |
| one pair, every push $r=3$ | 91.8 | 96.0 | 96.6 | 96.4 | – |
| Codex's heterogeneous example, direction undeclared, 400 repetitions | 74.5 | 79.5 | 80.0 | 96.0 | 93.8 |
| same, direction declared ($k=1$) | 95.3 | 96.3 | 95.8 | – | – |

Appending one exact zero reading to 20 weak readings moves the median 1.78 →
1.65 under the middle and 1.80 → 1.65 under ours28b. It does not move flat or
ours28 at all.

- **The rules that borrow do better in worlds whose pushes share a pattern**
  (all alike, or a fixed direction), because they learn it. The ruling says
  the definition must not assume such a pattern, so those worlds are not a
  standard it has to meet.
- **The heterogeneous example is handled the way the ruling allows.** The
  direction is declared per reading, and then every rule calibrates.

## 3. Result 1 (restated): no per-reading rule has the right width at every push size

From `ours27_one_scale_per_noise_pair_2026-09-28.md` §2, in this folder. For
a rule whose readings contribute separate factors, the width is right for
every sequence of pushes only if each reading's evidence
$e^{\phi(A)}$ satisfies $|\nabla\phi|^2-A\cdot\nabla\phi+\Delta\phi=0$. With
$u=e^\phi$ this is $\Delta u=A\cdot\nabla u$, whose smooth radial solutions
are constants. So under the ruling, some push sequences get widths that are
too narrow, and no choice of per-reading law avoids it. Tonight's example:
200 readings each with a push of one noise unit, ours28 covers 84%.

## 4. Result 2: no per-reading rule is never overconfident

The user asked whether some per-reading law is never overconfident, only ever
too cautious, which would pick the law's shape on honesty grounds. It is not
possible.

**Argument.** A per-reading score $\psi$ has mean zero at the tested mass
whatever that mass is. Differentiating that statement with respect to the
tested mass gives

$$
S=-E\big[\partial\psi\big]=E\big[\psi\;U_{\rm true}\big],
$$

where $S$ is the reading's contribution to the law's curvature and
$U_{\rm true}$ is the score of the data's actual distribution. For a reading
whose push is zero the data's distribution does not depend on the mass, so
$U_{\rm true}=0$ and $S=0$. Its score variance $V=E[\psi^2]$ is positive for
any rule that gives the reading non-zero weight, and a per-reading rule cannot
know that the reading is pure noise. So each such reading adds noise to the
score and nothing to the curvature, and enough of them make the law too
narrow, for every per-reading law.

**Numbers** (one pair, 3D readings, Monte Carlo with 300,000–400,000 draws).
$S$ at zero push was $\pm0.002$ or less for every law tried, zero within
sampling. $V/S$ by push length $r$ in noise units, and $V(0)$, the noise a
pure-noise reading injects (flat ours27: 3.00):

| scale law | $r=0.5$ | $r=1$ | $r=2$ | $r=3$ | $r=5$ | $r=8$ | $V(0)$ |
|---|---|---|---|---|---|---|---|
| $\mathrm{Beta}(2,2)$ | 6.18 | 2.03 | 1.02 | 0.89 | 0.89 | 0.95 | 0.640 |
| $\mathrm{Beta}(1,1)$ | 5.91 | 1.99 | 1.06 | 0.94 | 0.94 | 0.97 | 0.546 |
| **arcsine $\mathrm{Beta}(\tfrac12,\tfrac12)$ (ours28)** | 5.45 | 1.90 | 1.08 | 0.97 | 0.97 | 0.99 | 0.422 |
| $\mathrm{Beta}(\tfrac14,\tfrac14)$ | 4.82 | 1.78 | 1.08 | 0.99 | 0.98 | 0.99 | 0.291 |
| $\mathrm{Beta}(0.1,0.1)$ | 4.01 | 1.62 | 1.07 | 1.00 | 0.99 | 1.00 | 0.151 |
| $\mathrm{Beta}(0.01,0.01)$ | 2.75 | 1.32 | 1.03 | 1.00 | 1.00 | 1.00 | 0.019 |
| $\mathrm{Beta}(\tfrac15,20)$ | 0.14 | 0.04 | 0.02 | 0.02 | 0.08 | – | 0.000 |

Curvature credited to a reading relative to the best per-reading weighting
$r^4/(r^2+d)$, for the symmetric laws:

| $a$ in $\mathrm{Beta}(a,a)$ | $r=1$ | $r=2$ | $r=3$ | $r=5$ |
|---|---|---|---|---|
| 2 | 1.91 | 0.93 | 0.81 | 0.86 |
| 0.5 | 1.55 | 0.86 | 0.84 | 0.94 |
| 0.1 | 0.89 | 0.62 | 0.73 | 0.95 |
| 0.01 | 0.24 | 0.26 | 0.46 | 0.91 |
| 0.001 | 0.04 | 0.08 | 0.22 | 0.79 |

**What the shape trades.** A law more suspicious of weak readings (smaller
$a$) injects less noise from near-pure-noise readings and is less
overconfident at $r\approx1$. It pays by discounting moderately informative
readings. At the extreme it ignores nearly everything, which is "never
overconfident" only because it is uninformative (the $\mathrm{Beta}(\tfrac15,20)$
row). Every law is within a few percent for strong readings. So honesty does
not select a shape. The shape trades "too narrow when many readings are near
pure noise" against "wastes moderately weak readings", and **under the ruling
this trade is intrinsic: it can be placed, not removed.**

## 5. Where the definition stands after tonight

- **Form:** ours28 in the full-covariance form, each reading contributing its
  explained energy $H_i(m)$ on its declared candidate subspace, and the
  readout unchanged.
- **Premises added tonight:** no borrowing between readings (this record);
  betweenness (provisional), own frames and sign conventions, same candidate
  law same mass (`../ours27/ours27_c3_from_betweenness_and_the_scale_fork_2026-09-28.md`).
- **Declared shapes:** the starting law's and the per-reading scale law's.
  The current rule for both is "uniform in the noise-defined angle" (the
  half-Cauchy on $m/s$ and the arcsine on $\beta$). Relative to ours27 it cuts
  the noise injected by a pure-noise reading about seven-fold and keeps about
  85% of what a moderate reading says.
- **Width standard:** calibrated in the definition's own world (a theorem for
  a proper law). Worlds whose pushes carry an undeclared pattern are not the
  standard. Readings below about one noise unit make the law somewhat
  overconfident, a proven and bounded cost.
- **C3:** under "a reading is one push", the cross-reading regrouping that
  derived C3 in full is not a symmetry. The narrow form (no weight on the
  profile of per-reading sizes) stays a premise, as the full-covariance
  record proposes.

## 6. Status of ours28b

**Withdrawn under the ruling of §1.** It borrows push scale across readings.
The folder stays as the record of why a shared scale was tried, what it
achieved in patterned worlds, and why the definition does not use it. The
derivation record's two theorems (per-reading impossibility, and exact
calibration by per-class sharing) stand as mathematics. The first is Result 1
here.

## 7. Reproduction

Scripts in Claude's workspace, not kept; numpy and scipy.

- §2: the law on a 401-point grid in $\log m$; the per-reading scale on a
  161-point logit grid; the population over $(\mu,\tau)$ on a 15 × 14 grid
  with $\mu$ arcsine-spaced and $\tau$ log-spaced on $[0.25,2000]$.
- The tensor: eigenvalues on a 37-point logit grid and angle on 24 points,
  matrix-Beta(1,1) prior, with the evidence
  $\lvert I-W\rvert^{N/2}e^{\operatorname{tr}(WS_m)/2}$.
- §4: $E[\beta\mid H]$ and its derivative from Kummer functions, Monte Carlo
  over the reading.
