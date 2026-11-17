# State of the mass estimator

**Last updated:** 2026-11-09, by Claude (intake of the 9 Nov explorations E6, E9 and E10, left on their branches; "Where it stands" rewritten; hole G28). Update this page in
the same commit as any change it describes (see "Keeping this page true" at the end).

This is the entry point. A newcomer, human or AI, should be able to read this page,
open the files under "Read first", and be working within the hour.

**The living layer** (added 2026-10-24) is three folders and one file at the repo root:

| Place | Holds | Rule |
| --- | --- | --- |
| `definition/` | the deliverable: type, premises (F, L, R), one file per theorem, the equation | edit in place; git is the history |
| `problems/` | one file per open item: requirements, attempts, what killed each, questions for Can | a conversation's findings go here, not into a record |
| `rulings.md` | every ruling of Can's, dated, verbatim, append-only | a status may say "ruled" only if its date is in here |
| `tools/audit.py` | lists what is out of order across the three | run before committing |

Everything in `docs/` is now **archive**: the records that produced the living layer, kept
for citation and provenance. Read them when a theorem or problem file points at one; do
not reconstruct the state from them.

## What the project is trying to do

Find the function that is the definition of inertial mass in Newtonian mechanics.
Newton II has three corners. Two are total operations (F = m a, a = F/m). The third,
m from F and a, is not: one vector cannot be divided by another, and independently
measured force and acceleration never line up exactly. So mass has to be estimated,
and the project's claim (F1, ruled 27 Sept) is that **the estimator is the
definition**: there is no further mass behind the readout. F2 (ruled Sept) sets the
standard: any choice the premises leave open is an alternative Newtonian mechanics,
so the goal is uniqueness, not a good estimator.

**Where it stands (9 Nov).** No candidate meets the premise list, and the list is still being restated.
`definition/proof.md` opens with what is to be proved, part by part (9 Nov): existence and uniqueness without a named
function. Its chain from the type to ours30 (29 Oct, eighteen steps) is kept as the parts' arguments; what it proves is
conditional: on isotropic uncorrelated channels, any function meeting the premises as ruled, with R1 read as not
constraining the misfit weight, L11 read as fixed-N amplitude consistency, and L22, L23 and η's finite-N form taken as
given, is ours30 up to the declared inputs. ours30 fails two rows of its own list: R1 (a weight below 1 declares a
larger noise, G18) and L11 at fixed signal as N grows (the sign of F·a is in its power weight, G22). Each later probe
breaks: ours31 at weak signal, ours32 on T39's inputs, ours33 at the quoted width (O18). By part:

- **The input type:** settled, apart from errors shared between readings (O16).
- **What one reading contributes:** four candidate bodies (G19). One forcing argument is forming, L11 with T30
  against the power form (T46); the latent measure is unforced (B17, E5).
- **How much the series counts:** η's target is known (O14, 4 Nov); its finite-N form and its place (a power, a warp)
  are not (G3, G24, G26).
- **The starting law:** the family is forced and no member is; two routes to a member closed on 9 Nov (O6).
- **The point:** a family of location functionals of ln m; the median is the one in `met/` (G6).
- **The quoted width:** what it promises on noise is not written (G28).

Two rungs the proof leans on are not fixed in VD-Newton: R3 (persistence, lumping) is unwritten and R5 (Newton III)
is provisional (`definition/ladder.md`). 19 theorem files, ten on the path, still say "carried from the records"
(G15). The holes are typed in `definition/proof.md` §5, and that table is the queue.

**Since 9 Nov** (three exploration branches of Codex's; read by Claude and left on their branches, `../explore/README.md`):

- **The law cannot return to π on pure noise and keep L11's rate** (exploration E10; O15, O19). With reading order
  free, a law that washes out on noise is inconsistent on some identifying series (triangular arrays and one fixed
  history, the median included: both checked by Claude). So the bounded null tilt stays, and "we know nothing" is the quoted
  width's to carry; that statement is unwritten (G28), bracketed by this and B13's corner.
- **No principal member of π from the type or from the null limit** (E9, exploration E10; O6, O11). Noise ratios
  1, 2, 3 give a proper curve volume with no reciprocal centre (checked by Claude); no full-support π on positive
  mass is preserved on average by T42's pure-noise update (not checked). Sorted, provisional: this bars L3′ past one
  noise ratio, not L4; π is declared, centre and shape.
- **Shared signal structure** (E6; a pointer in O16). A joint law of mass and a shared calibrated direction fails T30
  and L19 as written; whether such evidence may enter is unsorted, and no hole depends on it.
- Briefs B21–B27 are open and not yet run; B24 (the width and η against the weak-instrument literature) and B25
  (which per-reading factors keep the truth) are first (`../briefs/README.md`).

**Since 3 Nov** (Codex's exploration and briefs B15–B19, on main since 7 Nov as delivered; **not checked by Claude**):

- **L11 as written cannot be met by any estimator** (T41): true signals $b/\sqrt i$ make the total energy diverge
  (like $\log N$) without ever identifying m. The errors-in-variables literature's condition for a consistent
  slope is that the total signal energy $\Delta_N$ grows faster than $\sqrt N$ (B19: Gallo 1982, Kukush–Van Huffel 2004,
  Shklyar 2018), the same threshold as T42's for the flat Bartlett law and T53's minimax one. So the strong premise
  was mis-stated; **restated 7 Nov (Claude, provisional)** with that rate, the old wording struck in its slot (G25, O19). A
  restatement, not a drop, so not evidence against the working assumption; open: the full type, and the phase-two bin.
- **ours33**, the sixteenth probe (`met/ours33.py`, O18): T48's exponent (ν+k)/2 with T47's normal-score warp in place
  of the power; point the raw median. Break: at small η the warped width is π's tail. **No width exists without π**
  (B16: the evidence alone is improper on every finite series, so π is in every quoted interval; ties the width to O6).
- **G19:** T48's latent measure is forced neither by T36 (B17) nor by the instrument's moments (E5). **G3:** η with T45's
  exact conditional score variance still misses its small-N target (B15). **Composition:** under ours33 the warp's order
  relative to the pushforward changes composite coverage, 87.75% against 91.75% (B18), as T51's non-commutation predicts (G27).
- Left on the exploration branches, not harvested: T54–T68, the joint, score-law and anchor line (`../explore/README.md`).

## The current equation: ours30 (since 29 Oct, R-2026-10-29-1)

Inputs per reading i: measured force and acceleration (any values, including zero),
an error covariance with every declared sigma > 0, a declared candidate subspace,
and (since 10 Oct) a window. Declared once: the starting law pi with centre c.

- **Law.** p(m) ∝ pi(m) · exp( eta · Σᵢ gᵢ Hᵢ(m) / 2 ), where Hᵢ is reading i's
  explained energy in noise units, gᵢ = min(1, κᵢ/Q₊ᵢ) weights it by its own misfit at its
  best **positive** mass Q₊ᵢ (κᵢ = kᵢ − ½), and eta = min(1, curvature / J) is the Bartlett
  power, with J the sum of the readings' squared scores, both taken at the series' best line
  on the whole projective line of slopes (eta = 1 when J = 0). It measures the series'
  own disagreement about m. ours29, the state before (30 Sept to 29 Oct), had gᵢ = min(1, ν/Qᵢ),
  ν = kᵢ − 1, on the whole-line misfit Qᵢ, which cannot see the sign of F·a.
- **Starting law.** Half-Cauchy, pi(m) = c / (c² + m²), with c = s = sigma_F / sigma_a
  by default (L3′, not yet ruled).
- **Point (since 28 Oct, provisional).** m̂ is the law's median (T35; L23 proposed;
  O2). Reported with the quoted interval (eta applied) and the raw law's interval beside it.
  The old readout, E[|F\*|] / E[|a\*|] over the candidates on the stacked length (T10),
  is kept as `readout` for comparison; it fails Z and reads the instrument's s.
- **Exact inputs** (sigma = 0) go to a separate three-way check: one mass, no mass,
  every mass. "Undefined" exists only there, at exact 0/0.

Stated in `../definition/equation.md`; the explainer for VD-docs is equation note d3
(`../outbox/VD-docs/`). ours29's defining record:
`noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md`. Premises and theorems:
`../definition/` (split from the lemma list d4).

**Run it.** In `met/` (isotropic or per-axis noise only), from the repo root:

```python
from met.analyst import estimate

LAW = {"reference": "cartesian2", "nuisance": "radius"}
OURS30 = {"rule": "likelihood_product",
          "prior": [{"uniform_angle": {"center": "first_reading"}}],
          "misfit": "kink", "positivity": True, "calibrate": "bartlett"}

# four 3D readings that disagree about m; sigma_F = 0.2 N, sigma_a = 0.1 m/s^2
F = [[4.0, 0.3, 0.0], [6.1, -0.2, 0.1], [1.0, 0.8, 0.0], [3.0, 0.1, -0.2]]
a = [[2.0, 0.0, 0.1], [2.0, 0.1, 0.0], [0.6, -0.4, 0.0], [1.0, 0.0, 0.0]]
estimate(F, a, 0.2, 0.1, ["median", "interval_95"], law=LAW, combine=OURS30)
# median 2.6451, interval_95 [2.0911, 3.4004]   (the third reading: g 0.079; eta 0.12)
# ours29 (no positivity): 2.6454, [2.0914, 3.4009] (g 0.064)
# without misfit and calibrate: 2.6262, [2.4203, 2.8549];
# the tempering widens the interval because the readings disagree
# all-zero readings return s = 0.2 / 0.1 = 2 (checked 2026-10-24 and 2026-10-29)
```

The returned dict also carries an `unresolved` flag.

**The full input type** (any covariance with cross-channel correlation, a declared subspace)
is `met/full.py`, for one series:

```python
from met.full import Reading, isotropic, ours30
r = ours30([Reading(F_i, a_i, Sigma_i) for ...], c=2.0)   # c is required: pi is an input
r["point"], r["interval"], r["raw_interval"], r["eta"], r["g"], r["Q_plus"]
# ours29(...) gives the state before 29 Oct; r["readout"] is the cloud readout, kept for comparison
```

The window is carried on each `Reading` but not used by the equation (T33). The batch
pipeline, presets and world still take isotropic and per-axis noise only.

**The forward model (2 Nov, Can's ask: "given F, a, … what distribution of noise would we get").**
`met/forward.py` holds the truth still (true pairs, Σ believed, optionally T37's ν_eff size) and
returns the draws of every noise measure the equations read, computed by met.full's own code:
E, H and K at the true mass, Q, Q₊, g29, g30, each reading's own line; per series K_total, η
(flat, ours29, ours30), curv, J, the series' line. The noise vector itself is returned split into
its three channels in the noise metric: latent (k), score (1, signed), misalignment (k − 1),
each exactly N(0, vI) at every nonzero push. For isotropic channels the laws depend only on
k, the push, the true line's angle atan(m₀/s), N and ν_eff (`invariants`). Exact reference
laws beside the draws (`reference`, `null_misfit`, `null_positive_misfit`,
`boundary_positive_misfit`); `tests/test_forward.py` checks them, B02's quadrature, and the two
engines against met.full.

```python
import math
from met.forward import forward, from_invariants, on_line, summarize
r = forward(from_invariants(3, push=4.0, angle=math.radians(45), readings=20), draws=4000, seed=1)
summarize(r)["Q_plus"]["mean"], summarize(r)["tempered30"], r["eta30"].mean()
r = forward(on_line(2.0, [[1, 0, 0], [0.5, 0.2, 0]], 0.2, 0.1), draws=4000, seed=1, nu_eff=10)
```

The table over that grid: `python -m met run presets/forward_noise_measures.yaml` (252 cells, 1 min).

**Pressure tools (Codex, 26 Oct).** `met/inspection/` inspects saved cases (any covariance,
declared subspaces, several readings) against independent geometry and justified checks,
with replayable reports: `python -m met inspect presets/inspection_cases.yaml`. Three studies
use it or its predecessors, all serial (`--workers 1`), no equation change:
`presets/misalignment_shrinking_noise.yaml`, `presets/misalignment_positive_boundaries.yaml`,
`presets/inspection_noise_paths.yaml` and `presets/inspection_prior_paths.yaml`. Their
findings live in `../problems/O13.md` (the sign of F·a, laws escaping to both boundaries
under a steady point) and `../problems/O2.md` (the perpendicular reading: law exactly π,
cloud readout R(c) ≠ c). T32 carries the positive-domain bounds.

**Real data.** `met/robot.py` turns VD-Newton's robot runs into `Reading`s: calibration
re-derived (matches the lab's `robot_cart.json` to 1e-12), both channels averaged over one
window each, IMU lag 8416 µs, budget = steady-state variance over the window's sample count
plus the calibration residual as a floor (Claude's declaration). `python tools/robot_run.py`
prints the table (ours30; `--equation ours29` for the table before 29 Oct). On 1-baseline
with 0.2 s windows ours30 gives median 0.9213, [0.9145, 0.9282] (ours29: 0.9213,
[0.9137, 0.9290]), against the static calibration mass 0.9323 (strict xfail under both);
static orientations held out give 0.9302 and 0.9320 and cover it. Why the dynamic runs
miss is O16. The calibration-side budget (autocorrelated white part, static rotation; B01, `--white`,
`--rotation`, `--matrix`) does not close it, so the miss points at error shared between readings.
The robot is a testbed for probes, not the experiment check (R-2026-10-31-1).

## Where the truth lives

| Question | Authority |
| --- | --- |
| What has been decided | Can's rulings, quoted verbatim with a date in the records. A dated ruling beats any summary of it, including this page. |
| What the premises are and what they force | The lemma list (`vd-docs-context/mass_estimator_lemmas_2026-10-10.md`, provisional). It lags one ruling: see Gaps. |
| What the code computes | `met/`, checked by `tests/`. |
| Whether a record's number is reproduced | `tests/test_reference_records_*.py`. A disagreement is a strict xfail naming both numbers (README rule 5). |
| How a study was run | Its preset in `presets/`, seed inside. |

The canonical copies of the decision records and the two notes live in VD-docs
(private). `vd-docs-context/` holds read-only copies. Only decisions go back, through
`../outbox/VD-docs/`.

## The breakthroughs still needed

### Closed 29 Oct: the sign of F·a (ours30)

ours29's misfit Q was taken on the whole projective line, so an anti-parallel pair had $Q=0$ and
full weight, and T30's bounded-confidence requirement failed there (`tests/test_sign_fork.py`,
strict xfail, kept as ours29's record). Sorted 27 Oct as a wrong premise (L8's count omits L1's
inequality); the misfit is the smallest error under the law, $m>0$ included (R-2026-09-07-1), so
it is $Q_+$; the numerator is O12's meaning read without a size, $k-\tfrac12$ (O13, 29 Oct).
**R-2026-10-29-1** switched the definition to ours30. The value $k-\tfrac12$ is proved (T38, 31 Oct); owed:
κ under a general covariance (unknown yet), only if G18 keeps O12's meaning. **31 Oct (B05):** the
fix costs consistency: a weight that sees the sign is correlated with the score, so ours30 fails L11 at
fixed per-reading signal as N grows; ours29's sign-blind weight does not (G22, a sort owed).

### O2: the readout (the median is provisional; 6 Nov: L23 does not force it, G6)

The three requirements were composition at the point (L5, L14), Z at the point (never ruled:
R-2026-10-29-2) and a value, c, at the all-zero experiment. Sorted 28 Oct (`../problems/O2.md`): under L19 the
point can only be a functional of the law of m; no law functional meets L5 and L14 at the
point on spread laws (Jensen; the law has no mean, T19); and on data no point meets either
at finite noise, the cloud readout and the median alike (lumping 0.27% strong, 21–23% weak;
common force 5.8% strong, 90% weak). So point-level composition at finite noise was the
mis-stated premise (TASK §3, outcome 3): the rungs commit exact values. **R-2026-10-28-1**
restates L5 and L14 at the level of the law in three clauses (exact inputs add; the
composite's law is computed from the parts' laws; a new measurement of the composite must
be consistent with it, as coverage), checked by `presets/composition_coverage.yaml`: both laws
cover at ~95% and overlap 99% wherever the parts do.

What is left for the point: Z-exact (every law functional is, T24), c at the all-zero
experiment for every c, the swap and units. Those leave a family of location functionals of
ln m. **L23** (proposed: the point commutes with every increasing reparametrisation, as T31
already asks of η) picks a quantile, and the swap makes it the **median** (T35, tested).
`met/full.py` returns it as `point`; the cloud readout (T10) stays as `readout` for
comparison, since it fails Z and reads the instrument's s.

The six-candidate table and the measurement-coarsening bridge (Codex, 28 Oct, in O2) are the
record of the candidates under the old point-level reading: on T4's full finite-law domain
every size rule fails composition or point Z; and fitting the summed measurements differs
from the pushforward of one joint law, which is why clause (c) asks for consistency, not
equality. **6 Nov (Claude, provisional; G6):** as a forcing argument L23 is mis-stated: it asks for
ordinal-scale invariance, and mass's commitments make it a ratio scale. Units and the swap leave the
reflection-equivariant location functionals of ln m, of which the robust ones remain; the point is a family, and
the median is the member in `met/`. Still open in O2: that sort, and whether the point must be continuous at the
all-zero experiment.

### O6: the shape of the starting law

Forced: symmetric under m/c ↔ c/m (L4), counted once (T7), declared before data
(T20), unit-covariant (L3). Free: any m⁻¹ g(ln(m/c)) with g even. The current
half-Cauchy is uniform in the angle the noise defines. Under ours28's model T28 found
the reference prior on mass to be uniform in that angle (the half-Cauchy, at c = s);
under ours29 it says only that pi's reference form does not depend on a size law, and
the Jeffreys-check record judged the agreement not to be evidence for the principle.
Only weak-signal answers depend on the shape. Under F2 it needs either a forcing
argument or a ruling that records it as a parameter of an alternative mechanics.

**9 Nov (E9, exploration E10; `../problems/O6.md`):** two routes to a forcing argument are closed. The geometry of the
declared instruments does not give a member of the family past one noise ratio, and the pure-noise limit selects
nothing. What is left is the question about the goal in `../definition/bounds.md`: whether a declared shape is a
branch cut or an alternative mechanics.

## Other open items

Queued for other agents (Codex), sorted, one deliverable each: `../briefs/README.md`.

| Item | Kind | One line |
| --- | --- | --- |
| O4 | proof | Gaussian noise forced? B05 (31 Oct): neither stated route forces it; L13 stipulated, conditionally forced by Herschel–Maxwell; the consistency route found G22 |
| O5 | proof | Write out the betweenness representation step behind T21 |
| O7 | ruling; part sorted as input (27 Oct) | Does Z cover zeros beside a signal? Separately: Z for mass and Z for noise split by whether Σ is declared exact or up to a scale |
| O9 | ruling | Is the window the criterion for one reading versus several? (proposed 10 Oct) |
| O11 | ruling | Does the default centre reference the instruments only (L3′), and is the swap-of-declarations premise accepted, which gives c = s for one noise pair? The several-instrument and anisotropic cases: sorted 9 Nov (Claude, provisional), declared, not derived |
| O13 | ruling; sign fork closed (ours30, 29 Oct) | Is anything more owed for readings that contradict each other? The sign of F·a: closed by R-2026-10-29-1 (g on Q₊, numerator k − ½; the value proved, T38) |
| O14 | proof | eta at small N: a correction that counts effective readings. The 26 Oct second noise channel is superseded by T36 (29 Oct). **30 Oct sort (provisional):** g and η are noise scales, per reading and per series, on an unwritten declaration that Σ is a floor; which levels a caller declares exact / floor / unknown is an input; G13 sorted |
| O15 | proof | Is the bounded null floor the best any law can do? **9 Nov (exploration E10):** an order-free law consistent at L11's rate cannot remove the null tilt; what the quoted width promises there is G28 |
| O17 | ruling and research | Mass as a metric: is isotropy of inertia (L22, proposed 27 Oct, read as a clause of L1) a premise; what it buys. On a tensor-mass world ours29 returns a confident pseudomass, not "no mass" (Codex, 27 Oct), so isotropy has to be declared |
| O19 | proof; premise sort | What growth of signal makes L11 attainable? T41: divergent total energy is not enough for any estimator; B19: the literature's rate is $\Delta_N/\sqrt N\to\infty$ (= T42, T53); L11 restated with it 7 Nov (Claude, provisional; G25). Open: the full type, the row's quantifiers, and (9 Nov) that the rate holds only beside a law that keeps its null tilt |
| O16 | fork (input, Claude's sort) | Errors shared across readings have no slot in the type. On the robot: the static rows show a real 1.15° force–IMU misalignment, but neither it nor a lever arm explains the miss, and the white budget is understated 2–19× by autocorrelated samples (28 Oct). Next: the budget from calibration-side data and the rerun (briefed), then the slot |
| L13, L18 | ruling | Confirm "the noise model is forced Gaussian, the data need not be"; is the subspace an input? |
| E10 | ruling | Does the VD entry still name the estimator without supplying it, now that composition forces the readout? (ours26 sweep §3.2) |

## Gaps (as of 2026-10-24; 27 Oct additions marked)

- (29 Oct) `tools/eta_sign_check.py` (T31) is a seeded script, not a preset kind of `met run`. (7 Nov) So are `tools/readout_gap.py` and `tools/cloud_readout_vs_median.py` (O2). (30 Oct) So is `tools/ours31_coverage.py` (O18), and (4 Nov) `tools/eta_vs_size_law.py` (O14).

- (27 Oct) L22 (isotropy of inertia) and T34 are Claude's, from the 26 Oct scratch session,
  whose two toy scripts (`comp.py`, `qnoise.py`) were not kept; numbers from them are marked
  "toy" in O14, O16, O17.

- The lemma list d4 still shows O1 open; the 29 Sept ruling closed it. The living layer
  records O1 as closed; d5 of the lemma list, when Can wants one, is generated from
  `definition/` (`python tools/audit.py --table`).
- Three of the 33 theorems have no full test (`python tools/audit.py`): T6 (the Herschel–Maxwell
  step is a derivation), T21 (a derivation; its representation step is O5) and T27 (needs a
  declared world of signals, which ours29 does not have). Their files say why.
- L3 (units covariance) is used by T14's argument but no theorem's From line cites it.
- The current equation note in VD-docs states ours25b, five equation states behind
  (ours26 to ours30), and its many-reading step rests on the strong form of axis
  splitting, ruled false on 29 Sept (L8). (29 Oct) Draft d3, for ours30, is in
  `../outbox/VD-docs/`, not yet carried to VD-docs.
- The ours29 record says "eta never moves the law's peak". With pi included the code
  moves it (ln m 1.1181 with eta 0.403, 1.2181 with eta 1). The likelihood part alone
  peaks at ln m 1.3181 for both, so the peak moves only because eta changes pi's
  relative weight: the record's sentence is wrong, not the code. Strict xfail.
- ours29 and ours30 have been run on the robot data only (Zenodo 10.5281/zenodo.11096791, CC-BY 4.0,
  in `../legacy/vd-newton-lab/raw/robot/`), not on the cart data. The robot's dynamic
  readings carry error beyond the declared budget (O16); no budget is ruled.
- The three strict rotation xfails in `tests/test_reference_records_ours29.py` sit near
  their 3e-8 bound (numpy 2.5.3 gives eta +3.3e-8, readout +6.0e-8, median +3.2e-8), so on
  another numpy or CPU one could pass, and so fail. The golden-output test no longer has
  this problem: it checks zero differences against 1e-7 (`ZERO_BOUND`), not their digits.
- ours25 (the sandwich, kept only for comparison) decides between the full law and the
  prior on rounding noise for some series. When a series' best line is at negative mass,
  its peak search on m > 0 ends at the grid's edge, where the law is flat to about 1e-8 and
  the finite-difference curvature H is the size of per-reading rounding (example: preset
  ours29_noise_vs_signal, cell c065, 4 replicates, series 3: H = 4e-9, so w lands on the
  1e-6 floor or on 1 depending on summation order). Found 24 Oct while making the batch law
  faster; the sandwich's search was left reading by reading so its recorded numbers stay
  reproducible (`sandwich_weight`'s docstring). Not a fork of the definition: ours25 is not
  the equation.
- `ratio_of_means_rss` in the batch pipeline depends on which series share a batch, at about
  1e-12 relative: `stacked_cond_alpha` picks the weighted route for the whole batch when any
  series has several noise pairs, so batch_elements can move a last digit. Older than the 24 Oct
  speed work (found by its review); the law itself is now chosen per series.
- The 10 Oct check script behind T32 is not yet kept (pending since 10 Oct).
- The reproducibility inventory still cites VD-Newton paths for scripts now in `../legacy/`.

## Read first

0. `../GOAL.md`: the goal in one page (target, how the work goes — probes that surface
   premises — the premise list, the probe history, the open holes), then `../TASK.md`: the
   task itself (find the function which is unique; what the readings cannot decide is
   input), with the logarithm example. Everything below assumes both. Before either, if you
   arrive with the textbook prior that estimators are choices:
   `existence_argument_audit_2026-10-29.md`, the argument for why the function exists,
   assembled and audited.
1. This page. Then `../definition/theorems/T36.md`: the frame everything is written in (the
   reading space, its noise metric, a mass as a subspace, one residual).
2. `../definition/bounds.md`: what Can's corrections have said the goal means, as Claude reads it, and the open
   questions about the goal (7 Nov). Then `../definition/ladder.md`: which commitment of the version removes which freedom, and
   what is left. Then `../definition/proof.md`: what is to be proved, part by part (its head, 9 Nov), then the chain from the type to ours30 with every
   step's status and the holes, typed; it is the queue. Then `../definition/equation.md`,
   `../definition/premises.md` and `../definition/theorems/README.md`: what is claimed and
   what it rests on.
3. `../problems/O2.md` and `../problems/O6.md`: the two items that need an idea.
4. `../rulings.md`: what Can has decided, in his words.
5. Then, only as the files above point at them: the lemma list d4
   (`vd-docs-context/mass_estimator_lemmas_2026-10-10.md`, the source the living layer was
   split from), the ours29 record, the undefined-mass record, and the 10 Oct units-and-axes
   record (the newest idea; aligned case only).

## Map of docs/

| Where | What | Tier |
| --- | --- | --- |
| `vd-docs-context/` | Copies of the 12 decision records and the current lemma and equation drafts | read first / reference |
| `noise_vs_signal/` | The work that produced ours29: misfit weight, Bartlett power, null series, ours28b | current work |
| `ours29/` | A pointer README only | pointer |
| `credit_and_originality_of_the_mass_estimator_2026-09-26.md` | What is borrowed (the law: Lindley–Beale, Liseo) and what is new (the selection principles) | reference |
| `mass_estimator_reproducibility_inventory_2026-10-01.md` | Which numbers can be regenerated, and from what | reference |
| `papers/` | Published papers, local only (not in git); `papers/README.md` has the citations | reference |
| `README_of_the_discussion_before_2026-10-22.md` | The history of this material in date order, with every relation | history |
| `ours26/`, `ours27/` | ours26 to ours28 working states | history |
| `mass_estimator_equation/` | The 21 Sept equation (ours21) with its spec, PDF and standalone code | history |
| flat `*_2026-09-*.md` | September work-state records (ours24, ours25b, axis splitting, priors, von Mises, combination) | history |
| `disagreement_rule_mass_estimator_trace.md`, `mass_estimation_zero_zero_two_cases.md`, `median_based_estimator_inertial_mass.md` | Records of 27 Aug to 5 Sept: the exact/measured split, the two zero-zero cases, and the median candidate | history |
| `literature_review/`, `unsanctioned/` | Output of one Codex session, 6–8 Sept, kept but never adopted | history |

## Keeping this page true

- Update this page in the same commit as the change: a new ruling, a new equation
  state, an open item closed, a change to the law in `met/`.
- When this page and a dated ruling disagree, the ruling wins: fix this page and say so
  in the commit message.
- Add one line to the log below at the end of any working session.
- A session writes into the living layer only: a theorem file, a problem file's attempts,
  a ruling appended to `rulings.md`, a test. It does not write a new record. If a result
  fits none of those slots, say so here under Gaps rather than inventing a fifth kind.
- Start a session by reading this page and running `python tools/audit.py` and the tests;
  say which mode the work is in (generation or pressure).

## Log

- 2026-11-10 (Codex, B21; not checked by Claude): T69 derives the maximal-invariant Wishart laws and profile curvature tending to Delta/(Delta+Nk); both laws fail Z and T30. O15 records the 2,400-series paired coverage study, group/MILE qualifications and unresolved L10/L19 permission; G3/G28 stay open. Full suite: 2,735 passed, 28 skipped, 25 xfailed; final focused checks 54 passed, 8 xfailed; audit unchanged. No production or premise change.

- 2026-11-09 (Claude, night): **`definition/proof.md`'s head restated**: what is to be proved no longer names ours30. A table by part (the type,
  the reading level, the series level, the starting law, consistency, the quoted width, composition, the point): the requirement, whether it is
  written with its quantifiers, what is known to meet it, and whether that is one body, a family or none. The consistency requirement written out
  for the scalar slice (sequences of experiments, the point; provisional). The chain and the holes table are unchanged.
- 2026-11-09 (Claude, night): **the fixed-history half of exploration E10 checked**, line by line, with its median corollary: no gap (O15).
  For each order-free procedure that returns to π on pure noise there is one series with Δ_N/√N → ∞ on which the median does not reach the
  truth. `../tests/test_washout_fixed_history.py` (20) holds the finite ingredients; the branch's 31 tests pass here. Unchecked still: the
  controls, and the averaged-null result behind O6's second closed route.
- 2026-11-09 (Claude, later): **queued, on Can's "ok add em":** briefs B24 (the quoted width and η against the weak-instrument literature; G28, G3),
  B25 (which per-reading factors keep the truth as the population optimum, and what bounded confidence adds; G19, G22), B26 (is π's family closed
  under lumping, and is that B18's weak-signal under-coverage; O6, G27) and B27 (L11's rate on vector readings; G25). B22 amended before any run:
  a sign-blind base beside ours30. Exploration E2 rewritten around G28's two walls and put first; E11 added (the sibling constant, Ohm's R = V/I:
  the phase-two bins and a venue for G14).
- 2026-11-09 (Claude): **intake of the 9 Nov explorations** (Codex: E6, E9 and exploration E10; read, left on their branches, nothing merged).
  E10: a law that returns to π on pure noise cannot keep L11's rate when reading order is free (the triangular-array bound checked by hand;
  O15, O19). E9 with E10: neither the type's geometry nor the pure-noise limit selects a member of π's family (ratios 1, 2, 3 checked by hand;
  O6); sorted, provisional: L3′ does not extend past one noise ratio (O11). E6: a pointer in O16. New hole G28 (what the quoted width promises
  on noise). "Where it stands" rewritten: it had called the equation in good shape since before G18 and G22. GOAL.md's ours32 row brought up to
  T39. No premise, equation or code changed. Suite 2,664 passed, 28 skipped, 17 xfailed; audit unchanged.
- 2026-11-07 (Claude, later): T70's misfit matched exactly to Hickey–Gour's geometric imaginarity; phase-two sort of the qubit picture into GOAL's bins
  added to T70; brief B23 (what makes mass real: Schur; parity in 2D) and exploration directions E6–E9 (coherent sum, gyroscopic alternative,
  Z-safe Gram estimates of the shrinkage, Fisher–Rao under the full type) queued.
- 2026-11-07 (Claude): **T70, the reading as a qubit** (Can's question: is misalignment's lost confidence the complex numbers?). Misalignment is the
  relative phase; explained energy is the Born rule; the whole-line misfit is imaginarity; the sign is not a phase (G22). The series' Gram cannot tell
  misalignment from disagreement, but T30 with Z keeps g per reading (T70 (6)). η_G = 1 − Q/H tracks η's target far better than the Bartlett η (O14); B22 queued.
- 2026-11-07 (Claude, afternoon): `claude/b21-brief` merged: brief B21 (Moreira's maximal invariant likelihood, O15/G3) queued;
  the 6 Nov G6 sort and the geometric-mean comparison readout (now on ours33 as well; `tests/test_geometric_mean.py`, 24). O2's
  B18 section moved into date order. No branch deleted.
- 2026-11-07 (Claude, afternoon): **the 5 Nov cloud-readout study committed** (`../tools/cloud_readout_vs_median.py`, O2): the
  cloud readout (T10) on ours30's and ours33's laws beside the median and the geometric mean; the 5 Nov numbers reproduce
  exactly with the same seeds. The cloud is a pull toward s on any law; whether a point that reads the sizes as well as the
  law is borrowing (L19, 28 Oct) is parked, not sorted. Tests: `tests/test_cloud_readout_vs_median.py` (6).
- 2026-11-07 (Claude, later): **the premises are deductions, Claude's to make** (R-2026-11-07-1, Can's words in `../rulings.md`):
  "ruled" is provenance, not a grade; grades say how much of the goal goes with a row; the hole type "rule (Can)" becomes
  *sort* or *goal* (a question to Can, asked as one); the Owed table is retyped as deductions Claude owes and questions about
  the goal. **L11 restated** (provisional, wrong premise): m̂ → m when Δ_N/√N → ∞, the literature's rate (B19; T42, T53); the
  "total push → ∞" wording struck in its slot (T41). New: `../definition/bounds.md`, Claude's running reading of the goal from
  the corrections, kept apart from TASK.md.
- 2026-11-07 (Claude): **intake.** Provisional protocol (Can, 7 Nov; `../briefs/README.md`, Intake): briefs merged as
  delivered, Claude owns the files, not a re-check; explorations are input. Harvested from `codex/ours32-consistency` by need:
  T39–T53, O19, G23–G27, the ours32 numerics fix and ours33 (squash and cherry-pick; each file marked). B15–B19 cherry-picked
  onto main as delivered (their branches sat on the exploration branch). E5's T48 section and tests harvested. T54–T68 left on
  the branches; `tools/audit.py` reads `../explore/README.md`'s table of them. Suite 2,619 passed, 28 skipped, 17 xfailed (main before: 2,189, 28, 16; the new xfail is T39's native-grid one); audit unchanged.
- 2026-11-05 (B19, Codex): O19 now compares functional EIV/TLS consistency conditions with T41–T53 and T61: the verified Gaussian one-slope rate is Delta_N/sqrt(N) -> infinity; dropping the multivariate eigenvalue-ratio condition does not weaken that scalar rate. Primary-source versions/pages and unverified originals are marked; no premise, input or estimator changed.
- 2026-11-05 (B18, Codex): extended the locked composition kind via explicit later-probe protocol; 400-replicate historical and q=c/4c worlds under ours30/32/33 show weak centred-sum undercoverage and an ours33 warp-order gap (87.75% pre versus 91.75% post-min), despite high overlap. O2 records all procedures, masses, MC errors and numerical scope; no warp order, eta combination, dependence model or equation adopted. 38 new checks pass. Full suite: 2,632 passed, 28 skipped, 17 xfailed; explicit 400-replicate ours29 reproduction passed; audit findings unchanged.
- 2026-11-05 (B15, Codex): T45 derives C33=sum[(nu+k)/(nu+k+2)]||RP'U||^2 under T48's specified residual model; the opt-in H/C study repairs the flat Gaussian leverage benchmark but still misses small-N signal fractions and keeps T39's pair at eta_C=0.0360. O14 records the 24-cell/three-leverage-size seeded study, cap rates and the fixed-signal conditioning limit. 27 new checks pass; defaults and the point unchanged. Final full suite: 2,621 passed, 28 skipped, 17 xfailed; audit findings unchanged. The standalone study replay reproduces every saved value exactly.
- 2026-11-04 (B17, Codex): T48/O18 refute the proposed T36 forcing of the latent scale measure: b=-1 also gives a normalized residual density; fixed-v normalization would force b=0 but is an additional requirement. G19 remains a sort/rule hole; stopped per brief, no exponent or premise changed. Counterexample tests: `tests/test_latent_measure_exponent.py` (20). Full suite: 2,614 passed, 28 skipped, 17 xfailed; audit findings unchanged.
- 2026-11-04 (Codex, B16): false-premise stop: a positive uniform bound proves ours33's evidence alone is improper in both mass and log mass for every finite valid series, including conflicting and strong readings; a finite nonnegative power does not cure it and T47 has no parent CDF. All proposed evidence-only interval cells are unresolved, so no module or coverage study is constructed. O18 holds the proof and scope; 21 independent projection/endpoint/integral/cutoff checks passed, as did mathematical review. Full default suite: 2611 passed, 28 skipped, 17 existing xfails; four Windows sandbox pipe failures passed on an isolated unrestricted rerun (six worker checks passed). Audit findings unchanged (T6/T21/T27, unused L24), diff checks clean. Production code, priors, points, output type and rulings unchanged.
- 2026-11-06 (Claude, evening): every equation's result dict now carries `geometric_mean` (exp E[ln m]) beside `median`; `point` unchanged.
  G6's fork measured (`tools/readout_gap.py`): the two centres differ by 5–9% on typical weak-signal series, 0.3% at signal 5.
- 2026-11-06 (Claude, later): **G6 sorted, provisionally** (O2): L23 is Stevens' ordinal-scale invariance, but mass's commitments (L1, L5, order)
  make it a ratio scale, so as a forcing argument L23 is mis-stated. Units and the swap leave the reflection-equivariant location functionals of ln m;
  L2 with O6 open removes the moment-based ones; the robust ones remain. G6 is now a sort, not a rule. Citations added to T5 (Horn, Umeyama) and T36 (TLS).
- 2026-11-06 (Claude): **brief B21 queued** (O15, G3). Moreira's (2009) maximal invariant likelihood is the one standard incidental-parameter
  cure not in O15's table and the only one acting across readings: on one instrument, $O(Nk)$ on the stacked rows collapses the $Nk$ latent
  pushes to one total $\Delta$ (noncentral Wishart, LIML/TLS point). Question for the brief: does its tempering reach η's target $\Delta/(\Delta+M)$?
- 2026-11-04 (Claude, evening): two boards. Briefs stay on `briefs/README.md` (branches `codex/<brief id>`, merged); exploration
  directions go on `explore/README.md` (branches `explore/<id>-<slug>`, never merged, each with a `SUMMARY.md`; E1–E4). B20 struck as moot.
- 2026-11-04 (Claude, on `codex/ours32-consistency`; harvested to main 7 Nov): **ours33** (`met/ours33.py`, O18): the probe Codex's T47 and T48 add up to and did not
  build — T48's exponent (ν+k)/2 on ours32's residual, η applied as the normal-score warp of the raw law instead of a power,
  point = raw median. Passes L2, Z, all-zero, units, swap, T30 and both T39 breaks (`tests/test_ours33.py`, 18). Coverage
  (`tools/ours33_coverage.py`, truth at c and at 4c): no better than ours32 at weak signal and 1.3–2.5× wider; the earlier
  cells' truth-at-c design rewarded the power form's pull toward c. **Break:** at small η the warped width is π's Cauchy tail
  (T39's pair: 78 decades around a raw interval of 1.5%). Names the output-type question: the quoted width is not a quantile of
  π·L. `met/ours32.py`: `bartlett_eta` takes a `terms` function (default unchanged). Equation unchanged.
- 2026-11-04 (Claude, later): **η's target, and what ours28 was doing** (O14; Can's standing hypothesis that the correct
  equation fudges back to ours28). On long series η → Δ/(Δ+M), and that is exactly the empirical-Bayes shrinkage
  β̂ = 1 − M/H² of ours28b's declared shared size law (Lindley–El-Sayyad): η is the size law's factor recovered from the
  Bartlett identity's failure instead of from a declaration. The cheat made precise: ours28 declared the law of the one
  quantity η must estimate, buying finite-N accuracy (N = 5, snr 0.5: η 0.90 against a target of 0.2–0.3) and, per
  reading, the weighting the leverage case (Codex's theorem 43, on its branch) needs. β̂ fails Z where η does not (a zero reading adds k to M); constraint
  on G3. `met/size_law.py`, `tools/eta_vs_size_law.py` (seed 2811), `tests/test_eta_target.py` (5). Equation unchanged.
- 2026-11-04 (Claude): **the brief queue is empty** (B01–B14 merged). ours32's reading level checked against G22: its population optimum is the
  truth (2.00 ± 0.01 at signal 1–2, 1D/3D) where ours30's power form gives 1.73/1.81, so the inconsistency is the power form's, not the sign's;
  a candidate forcing argument for G19 (O18; `tests/test_weight_consistency.py`). GOAL.md refreshed (ours32's row, G18/G19 standing, property N).
- 2026-11-04 (Claude): **B13 checked and merged.** Gleser–Hwang's corner holds for the type: any always-finite interval has infimum
  coverage exactly 0 at every finite N (direct total-variation proof; the 1987 PDF unavailable, a verified restatement used). Codex corrected
  the brief: at ce^{±3} the null coverage is positive (0.54/0.56). Claude's check: under ours30 that floor erodes with N (0.52 → 0.30 at
  N = 20 → 500 on pure noise; the law concentrates on c), G22's sign-in-the-weight mechanism in the corner.
- 2026-11-04 (Claude): **B14 checked and merged.** Barndorff-Nielsen's modified profile is not flat (Jacobian 1/|cos(θ − θ̂)|, exact for
  k = 1) but not a law: it has poles at positive masses on anti-parallel readings, and its Bartlett moments do not exist. O15's cure table
  is closed: no non-hierarchical latent measure does η's job, so η's job stays forced given L10 and L19.
- 2026-11-04 (Claude): **B10 checked and merged.** T5 (the swap, now for the median, not through T4), T12, T13, T22 written; T12's
  η = 32/225 example and T22's 13/10 vs 29/35 recomputed by hand. 19 theorem files still "carried", ten on the path (proof.md G15).
- 2026-11-04 (Claude): **B09 checked and merged.** T3, T4, T9–T11, T14, T15, T21, T29, L12, O5 and two ladder rows marked
  "comparison readout only, since 28 Oct" (still true under ours32, whose point is the median); nothing removed. `tools/audit.py` keeps
  every finding and tags T21 as not load-bearing; the theorem index gains a scope column (T37, T38 rows carried over).
- 2026-11-02 (later, Claude): **forward model on ours31/32** (O18, 2 Nov): `w31` and `eta32` added to `met/forward.py`
  (equal to met.ours32 to 1e-15) with `calibration` (weight against the true precision 1/v); preset
  `presets/forward_ours31_ours32.yaml`. In its own world ours31's weight is the true precision at strong signal; holes 2 and 3
  are ≤ 5%; so its weak-signal break is not the reading weight (confirms the 30 Oct re-attribution); ours32's η does not count the
  size spread twice, ours30's absorbs it (0.90 vs 0.95). Tests: `tests/test_forward.py` (36).
- 2026-11-02 (Claude): **the forward model** (`met/forward.py`, Can's "reverse function"): truth held still,
  the laws of E, K(m₀) and its channels, Q, Q₊, g, η, from met.full's code; preset kind `forward`
  (`presets/forward_noise_measures.yaml`, seed 20261102); `tests/test_forward.py` (28). It reproduces B02's null
  means, T38's k − ½ (only on the boundary; inside the half-line Q₊ → χ²ₖ₋₁) and B04's 22–29% tempered by the cap.
  T36 (4) refined: K(m₀) splits into a unit score and k − 1 misalignment, exactly, at every nonzero push. No equation change.
- 2026-11-01 (evening, Claude): **ours32, the fifteenth probe** (`met/ours32.py`, O18 1 Nov): ours31's reading level with η put back (Can: R-2026-11-01-1, which also leaves G18/G19 standing until forced; proof.md retyped them). Same-seed study `tools/ours32_coverage.py`: ours31's weak-signal break (0.57–0.71) goes to 0.82–0.87, ours30's level or slightly better; better than ours30 where the noise really varies (ν_eff = 10, 3D, N = 5). The remaining many-reading weak-signal shortfall is η's finite-N form and the latent measure (O14, O15). Tests: `tests/test_ours32.py` (8).
- 2026-11-01 (Claude): **B08 checked and merged.** π's family from L3, L4, properness: m⁻¹g(ln(m/c)), g even and normalized, or M = ce^X with X
  symmetric for general laws (δ_c included); T28 written. Claude's note: L11 adds a full-support clause (prior mass near every m), which
  excludes δ_c and gapped shapes. met/full.py takes only the half-Cauchy. 23 theorem files still "carried".
- 2026-10-31 (Claude): **B07 checked and merged.** T6's argument is in T06.md (Herschel–Maxwell, cited and proved; the Gaussian
  integral; T7's convention). Codex found my brief mis-sorted: L7 calibrates a point and cannot fix the coefficient ½. The remainder
  attaches to G18: by T36's homogeneity any τ ≠ 0 member is the flat factor under a declared (1 − 2τ)Σ (R1 with L24, candidate).
- 2026-10-31 (Claude): **B06 checked and merged.** T16–T19, T24, T25 written out; checked by hand (T25's correlated example gives H = 4/5 vs
  1/2 directly). Three statements lacked conditions and now carry scope lines: T16 (tilt gλ under ours30), T17 (fixed N only; false for
  ours30 through N, G22), T25 (zero axes inert only if uncorrelated; k = 1 through the conditional covariance). 25 files still "carried".
- 2026-10-31 (Claude): **B05 checked and merged; ours30 breaks L11.** Both routes to "Gaussian forced" fail as stated (Z admits
  non-Gaussian elliptical and independent-channel models; Codex's density checked to 9 digits): L13 is stipulated. The consistency
  counterexample confirmed with production code (1D, truth 2, signal 1, N = 3000: ours30 1.764 [1.597, 1.959]; ours29 2.002) and in
  2D/3D by population optimum (1.79, 1.81). Cause: g(Q₊) sees the sign of F·a, so it is correlated with the score; ours29's
  eigenvalue-only weight stays consistent. New hole G22 (sort); `tests/test_weight_consistency.py`. O15's ours30 N-study read as finite-N.
- 2026-10-31 (Claude): **B04 checked and merged.** The leverage identity E[J] = C(1 − 1/n_eff) holds for fixed nonnegative curvatures;
  ours30's observed curvatures are signed, so n_eff is not a count (G3). Strong-signal J/C (0.955, 0.943 in 2D, 3D) is E[g²]/E[g]: the cap tempers
  22–29% of readings the exact declared noise explains, T37 (3)'s efficiency cost seen from η (G18). η′ reported, not adopted.
- 2026-10-31 (Claude): **B03 checked and merged.** The untouched 28 Oct script and the new wrapper, both run at 400 replicates,
  give identical JSON and O2's printed table; the preset kind sits beside B01's robot_budget (merge conflicts in met/ resolved,
  one elif restored). The preset is locked to the historical ours29 run; composition under ours30 has not been checked.
- 2026-10-31 (Claude): **B02 checked and merged** (audit unchanged apart from T38; branch suite green). Proofs checked by hand;
  an independent Monte Carlo (four seeds, 4×10⁶ draws) agrees with the quadrature (|z| < 1.7). Written up as **T38**: sup E[Q₊] = k − ½,
  sup E[Q] = k − 1, neither attained; the brief's "only as both channels grow" was wrong. Load-bearing only if G18 keeps O12's meaning.
- 2026-10-31 (Codex, B14): O15 derives the exact regular scalar BN inverse Jacobian and the labelled Severini expected-score approximation: the isotropic correction is $|\cos(\theta-\widehat\theta)|^{-k}$, not constant. Its global continuation has nonintegrable positive-mass poles; restricting the fit to positive mass instead produces allowed inputs outside the regular interior-MLE domain. Eight tests check the Jacobian, expected covariance, properness obstruction, Z limits, coordinate invariance and private/shared constructions. The global-law/calibration branch stops under the false-premise protocol; no coverage of an undefined law is reported, no production code or premise changes, and no uniqueness claim for η follows. Full suite: 1,912 passed, 26 skipped, 16 xfailed; audit findings unchanged.
- 2026-10-31 (Claude): **B01 checked and merged** (audit unchanged; branch suite 1787 passed; default robot table byte-identical;
  one row reproduced from the CLI). The calibration-side budget does not close O16's gap: without the floor the misfit is 13–1100×
  what the declared noise produces and autocorrelation raises g only 2–5×, so the miss is error shared between readings (O16).
- 2026-10-31 (Codex, B13): O15 proves zero uniform coverage in the weak-identification corner for every finite reading count, even along nonzero-signal sequences, by a Gaussian total-variation argument. The original Gleser–Hwang PDF was unavailable; a published restatement and Dufour's author version were verified with that access limit explicit. The brief's fixed-$ce^{\pm3}$ zero-limit prediction is false: both limits are positive, with 0.5367/0.5600 observed in the single 300-series paired regression (about 5 s). The prior-averaged floor and regime (b) stay open/unchanged; no production code, prior, tempering, premise or PDF catalogue changed. Full suite: 1,845 passed, 26 skipped, 16 xfailed; audit findings unchanged.
- 2026-10-31 (Codex, B10): wrote T5's median-based swap proof (physical $1/m$ and calibrated $c^2/m$), T12's counted-once combination proof with the exact two-reading $\eta=32/225$ example, T13's projection/normalizer/median totality proof, and T22's covariance-congruence proof with $H(1)=13/10$ versus $29/35$ after induced correlations are discarded. Nine direct ours30 checks cover both swaps, the union, rotation and totality. G7's writing and B10's part of G15 are done, with maximizing-tie and compact-state conditions stated; no Statement, premise, equation default, `met/`, T6, T30 or T31 changed. Full suite: 1,853 passed, 26 skipped, 16 xfailed; audit findings unchanged.
- 2026-10-31 (Codex, B08): O6 derives the unit-covariant reciprocal family with both Jacobians, nonnegative normalization and converse; the general-law form is $M=ce^X$ for symmetric $X$, with the density formula its absolutely continuous case. T28 writes out the historical Fisher calculation and uniform-angle half-Cauchy member. The 36 new checks cover three null-prior families, medians, units, reciprocal symmetry, atomic/gapped controls and the marginal Fisher matrix; no general-prior API exists in `met/full.py`, so the three-member update tests are algebraic. G8's writing is done; O6 stays open, the default is unchanged, no member is selected and T36's measure question is parked. Full suite: 1,880 passed, 26 skipped, 16 xfailed; audit findings unchanged.
- 2026-10-30 (Codex, B07): wrote T06's conditional Herschel–Maxwell proof, exact Gaussian integral and four coefficient/width checks; G17, G13 and T36's measure question remain open. Found that literal L7 calibrates the point, not interval width: the archive's independent coefficient argument needs that extra condition, while L10 already states the flat factor; recorded the mismatch in T06 and proof step 3/G2 without changing the Statement, premises or `met/`. Full suite: 1,848 passed, 26 skipped, 16 xfailed; audit findings unchanged.
- 2026-10-30 (Codex, B06): wrote T16–T19, T24 and T25 proofs with all Statements unchanged. Added needed scope: T16 raw g=η=1 and negative-tilt attenuation; T17 fixed-N amplitude growth with positive g/η bounds, not generic limits 1; T18 unweighted Gaussian local asymptotic widths, since the exact mass-law mean is infinite; T25 matching covariance/candidate blocks or conditional Schur reduction. T19's tail and T24's zero-law proofs cover arbitrary covariance/subspace. All 25 new focused tests pass; full suite: 1,869 passed, 26 skipped, 16 xfailed; audit findings unchanged.
- 2026-10-30 (Codex, B05): T23 now derives the scope of both claimed Gaussian-forcing routes. Ours30's positive-misfit weights give a Gaussian-world counterexample to L11 with fixed weak nonzero signal and N growing; Z admits an explicit non-Gaussian independent-channel model. Fixed-N amplitude consistency is proved for Gaussian and bounded-mixture models. O4 stays open and proof step 2 is marked a gap; production code, premises and historical numerical assertions are unchanged. Six focused T23 checks pass; full suite: 1,847 passed, 26 skipped, 16 xfailed; audit findings unchanged.
- 2026-10-30 (night, Claude): **O15, the latent measure:** the flat law is a profile likelihood with N incidental parameters
  (Neyman–Scott); its information identity fails, and η = 𝓗/J is the standard quasi-likelihood correction. Every non-hierarchical
  cure coincides with the flat law here (Cox–Reid, reference prior: the push's information is m-free in whitened arc-length) or does
  not exist (conditional, marginal); the hierarchical cure, a size law, restores Bayesian calibration (0.96/0.93/0.94 vs the flat
  law's 0.70/0.85/0.91; `tests/test_latent_measure.py`) and is forbidden by L10. So η's job is forced; B12's noise test was the wrong
  test for it; L24 proposed narrowed to the reading channel. One cure unchecked: the modified profile likelihood, brief B14.
- 2026-10-30 (evening, Claude): **coverage sorted into two regimes** (O15). The corner (signal → 0) is Gleser–Hwang: no always-finite
  law covers uniformly; brief B13 checks the hypotheses. N growing at fixed signal is Neyman–Scott: the flat factor gets worse with N
  (0.84 → 0.72 at signal 1), ours30 reaches nominal (0.93–0.95 at N = 80–200; `tools/ours30_coverage_vs_n.py`). Property N as ruled is
  consistent with both; η has a real target, its form is the problem (L24).
- 2026-10-30 (evening, Claude): **ours31's break re-attributed** (O18): ours27's raw law and ours30 without η fail the same weak-signal
  cells as badly (0.575–0.71), so the cause is the flat latent measure (L10's endpoint, T18), which η patched, not ours31's size
  posterior. So η corrected the latent measure, not the noise (consistent with B12). The conflict is L10 + property N + L24 at weak signal.
- 2026-10-30 (evening, Claude): **ours31, the fourteenth probe** (`met/ours31.py`, O18), beside ours30, not replacing it: T37's exact
  cut factor per reading with ν_eff declared, no η. Premise report (`tests/test_ours31.py`): passes L2 on B11's inputs, Z, all-zero
  → c, T30 (anti-parallel included), units, swap, L24. **Breaks** at weak signal in 3D with N = 20 (coverage 0.57–0.71; ours30
  0.81–0.88), as predicted: the misfit's law given the noise size depends on the undeclared signal size (L10). Study:
  `tools/ours31_coverage.py` (a script, not a preset kind, as composition_coverage).
- 2026-10-30 (Codex, B04): exact fixed-curvature Gaussian score identity checked; signed fitted-curvature cancellation gives target expectation -infinity and J/C infinite second moment (O14). New `eta_effective` preset ran 72 cells × 1,000 paired series, retaining η/η′ coverage and explicitly empirical ratio errors; no correction adopted. Full suite: 1,767 passed, 26 skipped, 16 xfailed; 36 focused tests pass and audit findings unchanged.
- 2026-10-30 (Claude): **B11 checked and merged.** Codex found the brief's premise false: the integrated common-scale law
  π·K^−(n/2+a−1) with an improper scale prior is improper on valid measured inputs (an exact 1D fit puts a pole at its slope;
  boundary fits; the all-zero series), failing L2's totality, so it stopped before the coverage study, as the protocol asks.
  The cause is the improper prior on the scale: a proper one (a declared ν_eff, T37 G20) keeps every factor bounded.
- 2026-10-30 (Claude): **B12 checked and merged.** Audit unchanged; suite 1807 passed. Codex corrected the sketch (Bernstein gives a
  measure; a proper noise-size law also needs a normalisable inverse tilt) and found: every capped weight min(1, b/Q) fails Bernstein
  (a plateau forces a point mass), so ours29's and ours30's g fail L24 in k ≥ 2 under both normalisations; the smooth form has no
  proper law; the floor fails both normalisations; inverse-gamma controls pass in their own normalisation; **η fails too**: equal
  precisions give the same plateau, unequal ones make J depend on residual direction, and a shared noise size is a joint
  declaration the type does not make. Checked by hand: the kink jump, the improper smooth tilt, the η counterexample.
- 2026-10-30 (Codex, B12): corrected the Bernstein criterion with the proper-prior recovery gate and both mean conventions; `met/admissible.py` gives analytic verdicts for the requested forms and a finite-grid diagnostic that remains unknown. The 57-case isolated preset and 76 tests check the smooth form's improper inverse tilt, floor/inverse-gamma controls and the shared-noise-size obstruction for η at N>=2 (O14). No production equation, premise row or ruling changed. Full suite: 1,807 passed, 26 skipped, 16 xfailed; audit findings unchanged.
- 2026-10-30 (Codex, B11): the proposed common-scale density fails totality on valid 1D readings (aligned, boundary and all-zero examples); O14 derives the scale integral and its conditions, O7 records a proper null-padding example whose median moves, and `tests/test_scale_branch.py` checks them (19 pass). Stopped before the coverage study under the brief protocol; no production equation or ruling changed. Full suite: 1,750 passed, 26 skipped, 16 xfailed; audit findings unchanged.
- 2026-10-30 (Codex, B03): historical ours29 composition coverage is now an explicit preset kind; the old tool is a compatibility wrapper, and raw replicate outputs, job seeds and effective inputs are retained. Both the untouched original and the migrated 400-replicate runs reproduce all 21 recorded coverage/overlap entries and six L14 defect entries, and their full-precision legacy payloads are identical. Seed layout, numerical core and scientific records are preserved. Focused checks: 60 passed; slow record check: 1 passed; full regression: 1785 passed, 27 skipped, 16 xfailed; audit unchanged.
- 2026-10-30 (Codex, B01, pressure): declared the robot's autocorrelated white budget and optional static rotation from calibration-side data only; seeded preset `robot_budget.yaml` retains all 48 ours30 rows with ours29 medians, tau bootstrap spreads and rotation covariance. O16 contains the table and four predeclared answers: the gap remains, width compatibility is mixed, and no autocorrelated group covers the reference at all three widths. Defaults, existing xfails and equation unchanged. Full regression run: 1785 passed, 26 skipped, 16 xfailed; updated CLI guards: 8 passed; audit unchanged.
- 2026-10-30 (afternoon, Claude): **G20 sorted** (T37): α = ν_eff/2, the effective degrees of freedom of the stated uncertainty
  (GUM Annex G), on calibration certificates; T30 is the GUM's "no uncertainty is exact". Correction: practice normalises the
  precision (E[1/v] = 1), not the covariance; so ours30's κ = k − ½ is a σ declared with ν_eff = ½, and the cap holds exactly
  in 1D, failing in k ≥ 2 by ν/ν_eff. L24's row and B12 note both normalisations.
- 2026-10-30 (afternoon, Claude): `proof.md` brought up to T37 and L24. **The rulings at steps 9–10 conflict:** R1, R-2026-09-29-8
  and R-2026-10-30-1 (the σ's believed) against R-2026-09-30-10 (O12's cap) and R-2026-10-29-1 (κ's meaning); by homogeneity alone
  g < 1 reads misfit as evidence the σ's were too small. Sorted, provisionally: O12's meaning is the mis-stated one (G18, Can rules).
  What is proved is restated as conditional on R1 not constraining g. New holes G18–G21.
- 2026-10-30 (afternoon, Claude): **L24 candidate** (the law may not alter the declared noise law; a tempering is
  admissible iff it is the cut posterior of a mean-1 noise-size law), from Can's "any degree of freedom in the law
  that changes sigma as a banned one", R1 and R-2026-09-29-8. The test briefed to Sol as **B12** (Bernstein criterion;
  ours30's g, the smooth form, the floor, T37 controls, and η).
- 2026-10-30 (later, Claude): **R-2026-10-30-1 and T37.** Can: a misaligned reading is possible under the declared σ's,
  only unlikely; it conveys low information, direction kept. Formalised as Σ the covariance, each reading's noise size
  its shared state (L8), mean 1. Then T30 forces a power-law noise size, tail weight (ν + 2α)/Q with α > 1, so above
  k + 1; ours30's κ = k − ½ needs a noise size with no mean, and its cap g ≤ 1 holds only for the fixed-size Gaussian.
  A break of ours30, provisional; at strong signal it costs width (3D: 0.093 vs 0.084), not coverage. α is an input.
  `tests/test_noise_size.py`.
- 2026-10-30 (Codex, B02, pressure): O13 now proves the positive-domain supremum k - 1/2 and the sign-blind supremum k - 1, with finite nonattainment and the k=1 exception; independent compact quadrature checks a 256-cell grid per dimension with explicit exact-arithmetic discretization bounds (`tools/positivity_supremum.py`, `tests/test_positivity_supremum.py`). No equation change. Full suite: 1748 passed, 26 skipped, 16 xfailed; audit unchanged.
- 2026-10-30 (Claude, pressure): **what g and η are** (O14, provisional). Both are noise scales (tempering by w is
  declaring Σ/w); their caps read Σ as a floor, a declaration ours29/30 has run on unwritten (O7 gets the row).
  G13 sorted: the power is a scale read once (L21's cut, standard name), bounded confidence follows from the floor.
  On a common 4× inflation in 3D they over-correct together (g·η 0.165 vs 0.25); η alone gets it in 1D; g catches one
  bad reading. Declarations per level (exact / floor / unknown) are the input. `tests/test_scale_reading.py`.
- 2026-10-30 (Codex, B09): marked the nine theorems about the comparison cloud readout, O5 and the two series-size/weight questions as not load-bearing for ours30; the audit retains all three untested-theorem findings, labels T21's scope from its Status, and shows it in `--table` (checked by `tests/test_audit.py`). Full suite: 1707 passed, 26 skipped, 16 xfailed.
- 2026-10-29 (night, Claude): briefs checked against T36: B11 added (the raw law on O7's "up to a scale" branch,
  a Student-type law in the one residual; the scale prior's exponent sorted as an input); frame pointers in B04, B06,
  B07; B08 told to leave T36's measure item alone. T36 gets the declared-subspace count (2d − k per reading on that branch).
- 2026-10-29 (night, Claude): GOAL.md standing rule 6 and a line in CLAUDE.md, on Can's feedback: bring the standard
  mathematics in before working a fork (T36 was textbook and he had to supply it).
- 2026-10-29 (night, Claude): **T36 and R-2026-10-29-3.** On Can's "start over, with Einstein notation...
  have the metric emerge": a reading is a vector in the stacked force–acceleration space, whose only
  metric is the noise's (the axes share no unit; the 10 Oct record's parked §8.5); a mass is a
  subspace; Hᵢ is projection and Eᵢ − Hᵢ the residual; k residual dimensions per reading. η's and the
  Q channel are one residual; the 26 Oct pooled scale's bias is the incidental-parameter count (O14).
  Kept central on his ruling: pointers from equation.md, proof.md, GOAL.md and this page. Equation
  unchanged. `tests/test_reading_space.py`.
- 2026-10-29 (late, Claude): `GOAL.md` at the root and `docs/existence_argument_audit_2026-10-29.md`
  added (both also project docs on claude.ai); `definition/ladder.md` §2 gets a bin column: every
  measurement-constitution premise sorted as Newton's / measurement's / the dictionary's / input,
  on Can's "have they always been in Newton?" (Claude's sort, provisional). The sort needs four
  bins, not three: the dictionary's bin (L10, L17, L19, L23, no-gate) is the project's own
  contribution. Read-first updated.
- 2026-10-29 (night, Claude): **R-2026-10-29-2.** Z at the point was never ruled: the 27 Sept
  property is about the law, and the ours27 record had recorded his position as the opposite
  ("the point prices it in"), "not yet as a ruling". The lemma list's extension was carried into
  L20, O2's R-b and the pooled-scale verdict as if ruled. Corrected: L20's row, O2, O7, O14 (the
  second noise channel reopened; the bias objection stands), proof.md steps 13, 16, 17 and §3.
  The median's replacement of the cloud readout rests on L19 and the no-go, not on Z.
- 2026-10-29 (evening, Claude): `definition/proof.md` started: the chain from the type to ours30
  in eighteen steps, each marked ruled / forced / forced (argument in archive) / input /
  provisional / gap, and sixteen holes typed as write, derive, sort, rule or VD-Newton. Findings
  from writing it: the switch to the median took T3, T4, T9–T11, T14, T15, T21, T29 and L12 off
  the load-bearing path; T5 now rests on the law's swap symmetry, not T4; η's finite-N form is
  the one unsorted freedom inside the law; T30's bounded-confidence requirement has no named
  premise; the 1D candidate factor rests on an unstated premise. Briefs B05–B10 cut from the
  holes table (O4; six law theorems; T6's argument; π's family; the demoted theorems'
  bookkeeping; four short proofs). The three sorts (η at finite N, T30's requirement, the 1D
  premise) stay Claude's.
- 2026-10-29 (Claude): **the switch to ours30** (R-2026-10-29-1). `met/full.py`: `ours30()` is the
  definition, `ours29()` kept; the batch pipeline takes `positivity: true` (closed form for
  isotropic noise, a positive-slope search for per-axis noise, 1D included);
  `tools/robot_run.py` defaults to ours30. `definition/equation.md`, T8, T30, T32, L8, L21, the
  ladder, O13, TASK's table switched. Equation note d3 (the explainer) in `outbox/VD-docs/`.
  Robot baseline median unchanged to 1e-4, interval [0.9145, 0.9282]. `tests/test_ours30.py`
  (parity with the batch pipeline, 1D, the validator, the robot row). No ours29 number moved.
- 2026-10-29 (Claude, pressure): O13's misfit grounded in R-2026-09-07-1 (smallest error under the
  law, m > 0 included: the 7 Sept record's anti-aligned infimum 6 is Q₊; ours29's Q is 0). η's
  whole-line evaluation checked under the same ruling: not the same slip (forced by Z and the
  stationary point; negative-φ̂ series undercover as much without η, and η raises them).
  `tools/eta_sign_check.py`, T31.
- 2026-10-29 (Claude, pressure): O13's numerator sorted. g is a function of Q₊ (L21 + L1),
  which excludes the two-factor form; κ is O12's meaning read size-free (L10), the supremum over
  sizes: k − 1 for Q (the ruled ν) and k − ½ for Q₊ (seeded grid; proof queued as B02). Candidate
  `ours29(..., positivity=True)` in met/full.py, default unchanged; T30's tilt bounded on
  anti-parallel pairs in 1D–3D; acute pairs discounted by (k − ½)/(k − 1) less. Changes the ruled
  ν = k − 1: Can's check owed. `tests/test_sign_numerator.py`.
- 2026-10-29 (Claude): `briefs/` added, the queue for Codex runs: one sorted fork and a named
  deliverable per file, each run on its own branch from `main`, checked and merged by Claude
  (protocol in `briefs/README.md`; CLAUDE.md points there). Queued: B01 robot budget (O16),
  B02 exact S_k (O13), B03 composition coverage as a preset (O2, Gaps), B04 J's bias against
  effective readings (O14). Not queued: the T32 example numbers (inputs not recorded; where the
  script goes is a ruling owed) and the three rotation xfails (a call about ZERO_BOUND's rule).
- 2026-10-28 (night, Claude, pressure): O16 checked before briefing the slot. The floor is the
  static fit's zero-mean residual, not a shared offset; the static rows show a 1.15° force–IMU
  misalignment that does not close the 1.4% floor/no-floor gap; a constant offset (lever arm)
  fitted across runs does not either; the rest recordings are autocorrelated (tau up to 19
  samples), so the white budget is understated 2–19×. Sorted: fix the declared budget from
  calibration-side data first (briefed to Codex), then the slot. `tests/test_robot_budget_checks.py`.
- 2026-10-28 (merge, Claude): `sorts-2026-10-27` and the Codex chain of 27–28 Oct
  (`codex/composition-measurement-bridge`, six runs) merged into main. Conflicts were
  appended text only (O2, O17, T04, this log); O2 carries a merge note placing the Codex
  candidate and bridge results against R-2026-10-28-1. No number moved.
- 2026-10-28 (later, Claude): R-2026-10-28-1, composition restated at the level of the law.
  `tools/composition_coverage.py` (seed 2810, 400 replicates per size): computed and measured
  laws of the lump and of the reduced mass cover at ~95% and overlap 99% wherever the parts
  cover; L14 on data, cloud and median fail alike. `met/full.py` returns `point` = the median;
  equation.md states it. No recorded number moved. The study is a tool script, not yet a
  preset kind (Gaps).
- 2026-10-28 (Codex, pressure): O2's physical-composition bridge and T4's scope:
  shared-channel covariance, omitted Gaussian contrast, special-split equality and
  all-positive-noise fine/coarse counterexample, checked in `tests/test_composition_bridge.py`.
  The 28 Sept archived derivation is credited; no equation/type change or domain closure proof.
  Serial regression: 175 passed, one existing strict xfail (source orientation count);
  final 21 bridge checks pass after adding the common-force dual. Audit unchanged.
- 2026-10-28 (Claude, pressure): O2's no-go. Under L19 the point is a functional of the law;
  no law functional meets point-level L5 and L14; on data neither the cloud readout nor the
  median meets L5 at finite noise (0.273% vs 0.273% strong, 21% vs 23% weak). Sorted: L5/L14
  as finite-noise point constraints are mis-stated; the rungs commit exact values, delivered by
  T1 and L11. L23 (parametrisation-free point) proposed; T35: with the swap it forces the
  median. Candidate readout: the law's median. `tests/test_point_composition.py`. Equation in
  `met/` unchanged.
- 2026-10-28 (Codex, pressure): O2 perpendicular-first six-candidate table, exact Gaussian
  point-Z proof (global RMS 1 → √2) and propagated finite-law composition counterexamples
  (conditional mean 2 ≠ 10/3); deterministic serial preset `presets/readout_candidates.yaml`.
  Serial regression: 521 passed, three existing strict xfails; 22 study/replay checks
  passed after provenance review, including one additional real replay test. Audit unchanged.
- 2026-10-27 (evening, Claude, pressure): L12 checked against the ladder. The series join is
  second order (under 0.2% across l_1 to l_inf at weak signal) and fixed by no rung: unknown
  yet, parked. O2's Z failure located: null readings replace the readout's weight by
  (1 + m²/s²)^(-1/2), a function of Σ alone (1000 nulls: 2.3931 against 2.3930). Sorted
  provisionally as outcome 3, T21's domain against L17's reason; the three conditions a repair
  must meet are in O2. `tests/test_readout_weight.py` (one slow case under MET_SLOW).
- 2026-10-27 (Codex, pressure): O17's excited-eigenspace criterion, conservative scalar
  compatibility diagnostic and concentrating-law counterexamples; paired tensor worlds
  in `presets/anisotropic_world.yaml`, with independent geometry/angle/world checks.
  Serial regression: 460 passed, three existing strict xfails; audit unchanged.
  L22 remains proposed; production type and equation unchanged.
- 2026-10-27 (Codex, pressure): O16's shared rotation field, observability counterexamples,
  and noise-aware robot excitation; serial seeded recipe `presets/shared_sensor_tilt.yaml`.
  Baseline resolves at least two directions under both budgets; fits remain descriptive.
  Source count 23 versus 24 retained as a strict xfail. Equation and robot loader unchanged.
  Serial regression: 402 passed, three strict xfails (sign, robot interval, source count);
  audit findings unchanged.
- 2026-10-27 (Codex, pressure): O14's general Q bias curve, finite-N moments and
  scale-only null-point counterexample, linked in O7; serial paired recipe
  `presets/pooled_noise_scale.yaml` and independent reference/replay tests. Serial regression:
  360 passed, one existing strict xfail; audit unchanged. Equation unchanged.
- 2026-10-27 (Codex, pressure): exact positivity-null moments, a finite counterexample
  to the proposed universal numerators, and the acute-preservation obstruction in O13;
  candidate checks in `tests/test_positivity_candidate.py`. Serial regression: 320 passed,
  one existing strict xfail; audit findings unchanged. Equation unchanged.
- 2026-10-27 (Claude, pressure, branch `sorts-2026-10-27`): the sign fork sorted as a wrong
  premise (L8's count; O13, T30, T32; `tests/test_sign_fork.py` with a strict xfail on T30);
  the Z split sorted as an input (O7); L22 proposed; the 26 Oct scratch session folded into
  T34 (tested), O17 (new), and attempts on O2, O14, O16. CLAUDE.md: how to brief Codex runs.
  Later: R-2026-10-27-1 (moves are forced or "unknown yet"; Claude's sorts carry the status
  *provisional (from the goal)*, not a request for a ruling) and R-2026-10-27-2 (no energy in
  this version; the estimator's form is a theorem of the commitments after the headword;
  TASK.md §7). O17's two ruling questions answered; L22 read as a clause of L1. Then
  `definition/ladder.md` from VD-Newton `61afb9e`: rungs R0–R6, the measurement constitution as
  a second column, every open freedom's status; findings: L5's readout form may sit at R5 (Newton
  III), not R3; E10 d4's "estimation policy" slot is what this work fills; R3 (NML2.2) is
  unwritten and R5 provisional, so the body is provisional on both.
- 2026-10-27 (Claude, merge): the Codex branch merged into main. Its five study paragraphs
  above became one pointer block; the findings stay in O13, O2 and T32. The flat-CDF median
  defect is fixed: on a run of nodes where the CDF equals the probability to 1e-12,
  met.full returns the run's midpoint in log m, the one plateau point that commutes with the
  swap and units (`_plateau_centre`). Both strict xfails pass and lost their marks; the
  inspector's reciprocal-median check now passes, with a test that it still fails on a wrong
  median; a swap test at c = 0.3 and 3. No recorded number moved.
- 2026-10-26 (Codex, pressure priors): inspected twelve fixed-input cases over three
  declared prior centres, added shared independent CDF/cloud-point references, and
  isolated prior multiplication from data-only information and cloud weighting.
  Exact perpendicular readout and reciprocal opposite limits recorded in O2/O13;
  known median defect retained for the preset's exact covariance. Audit unchanged.
  Serial regression: 281 passed, 2 strict xfails for the existing median defect.
- 2026-10-26 (Codex, pressure paths): inspected 48 cases through the modular tools,
  added callable path transforms and an independent angle-probability reference,
  and recorded the information flow in O13. All results/replays retained; native
  tail warnings and the known median defect remain visible. No equation change.
  Serial regression: 247 passed, 1 existing strict xfail; 32 focused follow-up checks
  cover close path coordinates and widely shifted priors. Audit findings unchanged.
- 2026-10-26 (Codex, pressure tools): added the composable general case inspector,
  independent restricted-subspace fits, source conversion with isolated truth, a
  catalogue of justified checks and self-contained reports. Eight replayable cases;
  known scientific/numerical failures remain visible. Existing equation unchanged.
  Regression checks: 217 passed, 1 existing strict xfail, with worker-spawning tests
  omitted. Includes all eight frame transformations and a native-density corruption
  check. Follow-up checks also guard exact-symmetry applicability. Audit unchanged.
- 2026-10-26 (Codex, pressure follow-up): derived positive-domain residual bounds and
  normalization identity; added a 20-cell symmetric/boundary study and independent
  tests. Stable point and mathematical median can coexist with probability escaping
  every finite band. Kept the SD 0.1 numerical median discrepancy as a strict xfail;
  estimator unchanged. Targeted checks: 72 passed, 1 strict xfail. The exact findings
  and remaining physical fork are in O13; the existing audit findings are unchanged.
- 2026-10-26 (Codex, pressure): fixed-reading misalignment diagnostics, 32 cells from
  `presets/misalignment_shrinking_noise.yaml`, through the existing `met run` command.
  Independent residual-covariance projection, positive-domain boundary infima and Rice
  cloud mean; full-precision curves/corrections, optional self-contained report and PNG/SVG.
  Observations recorded in O13, without adopting a new criterion or changing the equation.
  Added 45 tests using independent least squares, analytic limits, angle quadrature and
  reproducible output checks. Suite: 1,236 passed, none failed, 13 xfailed, 26 skipped.
  Existing audit findings unchanged (T6, T21, T27 and unused L3).
- 2026-10-25 (speed, 2): `python -m met run` runs cells side by side (`--workers`, default half
  the cores, `max` for all); output bit for bit the same for any count, checked on all 22
  presets (tests/test_parallel.py). Making that true found a defect: ratio_of_means_rss's
  spline table for many readings was widened on demand, so a cell's value moved by ~1e-13 with
  what ran earlier in the process; it is now a fixed ladder of tables (values moved <= 8e-14).
  met/full.py's readout integral: window +-60 with 401 nodes, replacing +-40 with 801: A(m)
  went from 6e-13 to 8e-16 relative error, with half the nodes; robot table identical, 8 s.
  Suite: 1,191 pass, none fail, 13 xfail.
- 2026-10-24 (night, speed): no numbers changed beyond rounding. The batch law no longer
  computes each reading's E[alpha | m] unless the plain ratio_of_means readout asks for it
  (bit for bit the same), and with one s per series the likelihood product is summed from
  three sums instead of reading by reading (every preset's readouts within 8e-14 relative of
  before; tests/test_summed_law.py). met/full.py's readout integral is bit for bit the same
  and 2.5 times faster. The ours29 study at 50 replicates: 50 s (the same cells took over
  10 minutes before); the robot table: 28 s to 12 s; the suite: 225 s to about 135 s
  (reference scripts now run side by side; test_stacked's oracle density written out).
  Found on the way: ours25's fallback is decided by rounding (Gaps). An independent review
  confirmed the equivalences and found two faults, both fixed: the side-by-side reference
  scripts were relaunched and orphaned, and the fast path was chosen per batch, not per
  series. Suite: 1,173 pass, none fail, 13 xfail.
- 2026-10-24 (night, last): the golden-output test checks printed zero differences against
  1e-7 instead of their digits (the old failure was one of them); requirements-lock.txt pins
  the versions last run. Suite: 1,114 pass, none fail, 13 xfail.
- 2026-10-24 (night, later): tests/test_theorems.py (25 tests) for T2–T4, T6, T7, T15–T20, T23,
  T26, T28; each file names its test. Found on the way: T15's condition is on g-weighted sums;
  T23's Z route needs independent channels, which no premise states (O4); eta's spread on one
  series of 400 is about 0.09 (O14). Suite: 1,112 pass, 1 fails (the old golden output), 13 xfail.
- 2026-10-24 (night): met/robot.py, tools/robot_run.py, tests/test_robot.py (7 tests, one
  strict xfail: the baseline interval misses the calibration mass), tests/test_window.py (T33's
  example; T33 now names its tests). O16 opened from the robot run. Suite: 1,087 pass,
  1 fails (the old golden output), 13 xfail.
- 2026-10-24 (evening): met/full.py, the full input type, with tests/test_full.py (34 tests:
  agreement with the reference and the batch pipeline; T1, T5, T13, T14, T22, T24, T25, T29, T32).
  Suite: 1,078 pass, 1 fails (the old golden output), 12 xfail.
- 2026-10-24 (evening): TASK.md added, the principle as the task (R-2026-10-24-3); CLAUDE.md
  and this page now send readers there first.
- 2026-10-24 (later): the living layer added: definition/ split from lemma list d4, problems/
  (O2 and O6 written out), rulings.md (62 entries extracted from the records, 28 quotes
  re-checked), tools/audit.py. docs/ relabelled as archive.
- 2026-10-24: page created from a full map of the tester, VD-docs and VD-Newton, then
  corrected after an independent fact-check (c = s qualifier on the cloud readout; the
  value c derived, not ruled; eta's J; the example replaced by one where the tempering
  acts). Tests: 1,044 pass, 1 fails (golden output, see Gaps), 12 strict xfails.
