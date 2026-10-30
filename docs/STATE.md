# State of the mass estimator

**Last updated:** 2026-10-30 (later), by Claude (T37: the σ's believed; ours30's κ and cap break on it). Update this page in
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

**Where it stands.** The equation is in good shape; the proof is not yet written. The chain
from the type to ours30 is in `definition/proof.md` (29 Oct), eighteen steps: what is proved
today is that on isotropic uncorrelated channels any function meeting the premises as ruled,
with L22, L23 and η's finite-N form taken as given, is ours30 up to the declared inputs. Those
three, plus 31 theorem files whose proofs still say "carried from the records", are the
uniqueness half's work; the holes are typed there and that table is the queue. O2 (the readout) is provisionally
resolved as of 28 Oct, pending Can's check of L23; O6 still needs a real idea, not more
tests. The law's sign defect is closed by ours30 (29 Oct). Some premises are graded below bedrock or are
provisional (L12 medium; L14 and L15 provisional), and some choices are natural
rather than forced (eta's cap at 1 and its power form at finite N). The misfit weight's
numerator κ = k − ½ is fixed by meaning (O12, read size-free); its value awaits a proof (brief B02). The rest is rulings Can owes and proofs, some routine
(O4, O5) and some open research (O14, O15).

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
miss is O16: the declared budget is not yet right (28 Oct; brief B01).

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
**R-2026-10-29-1** switched the definition to ours30. Owed: the proof that $k-\tfrac12$ is the
supremum (brief B02); κ under a general covariance (unknown yet); ours30's coverage, not yet measured.

### O2: the readout (provisionally resolved 28 Oct; Can's check of L23 pending)

The three requirements were composition at the point (L5, L14), Z at the point (never ruled:
R-2026-10-29-2) and a value, c, at the all-zero experiment. Sorted 28 Oct (`../problems/O2.md`): under L19 the
point can only be a functional of the law of m; no law functional meets L5 and L14 at the
point on spread laws (Jensen; the law has no mean, T19); and on data no point meets either
at finite noise, the cloud readout and the median alike (lumping 0.27% strong, 21–23% weak;
common force 5.8% strong, 90% weak). So point-level composition at finite noise was the
mis-stated premise (TASK §3, outcome 3): the rungs commit exact values. **R-2026-10-28-1**
restates L5 and L14 at the level of the law in three clauses (exact inputs add; the
composite's law is computed from the parts' laws; a new measurement of the composite must
be consistent with it, as coverage), checked by `tools/composition_coverage.py`: both laws
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
equality. Still open in O2: Can's check of L23 against the goal, and whether the point must
be continuous at the all-zero experiment.

### O6: the shape of the starting law

Forced: symmetric under m/c ↔ c/m (L4), counted once (T7), declared before data
(T20), unit-covariant (L3). Free: any m⁻¹ g(ln(m/c)) with g even. The current
half-Cauchy is uniform in the angle the noise defines. Under ours28's model T28 found
the reference prior on mass to be uniform in that angle (the half-Cauchy, at c = s);
under ours29 it says only that pi's reference form does not depend on a size law, and
the Jeffreys-check record judged the agreement not to be evidence for the principle.
Only weak-signal answers depend on the shape. Under F2 it needs either a forcing
argument or a ruling that records it as a parameter of an alternative mechanics.

## Other open items

Queued for other agents (Codex), sorted, one deliverable each: `../briefs/README.md`.

| Item | Kind | One line |
| --- | --- | --- |
| O4 | proof | Gaussian noise forced under the tempered evidence |
| O5 | proof | Write out the betweenness representation step behind T21 |
| O7 | ruling; part sorted as input (27 Oct) | Does Z cover zeros beside a signal? Separately: Z for mass and Z for noise split by whether Σ is declared exact or up to a scale |
| O9 | ruling | Is the window the criterion for one reading versus several? (proposed 10 Oct) |
| O11 | ruling | Does the default centre reference the instruments only (L3′), and is the swap-of-declarations premise accepted, which gives c = s for one noise pair? Then the several-instrument and anisotropic cases |
| O13 | ruling; sign fork closed (ours30, 29 Oct) | Is anything more owed for readings that contradict each other? The sign of F·a: closed by R-2026-10-29-1 (g on Q₊, numerator k − ½; proof of the value: B02) |
| O14 | proof | eta at small N: a correction that counts effective readings. The 26 Oct second noise channel is superseded by T36 (29 Oct). **30 Oct sort (provisional):** g and η are noise scales, per reading and per series, on an unwritten declaration that Σ is a floor; which levels a caller declares exact / floor / unknown is an input; G13 sorted |
| O15 | proof | Is the bounded null floor the best any law can do? |
| O17 | ruling and research | Mass as a metric: is isotropy of inertia (L22, proposed 27 Oct, read as a clause of L1) a premise; what it buys. On a tensor-mass world ours29 returns a confident pseudomass, not "no mass" (Codex, 27 Oct), so isotropy has to be declared |
| O16 | fork (input, Claude's sort) | Errors shared across readings have no slot in the type. On the robot: the static rows show a real 1.15° force–IMU misalignment, but neither it nor a lever arm explains the miss, and the white budget is understated 2–19× by autocorrelated samples (28 Oct). Next: the budget from calibration-side data and the rerun (briefed), then the slot |
| L13, L18 | ruling | Confirm "the noise model is forced Gaussian, the data need not be"; is the subspace an input? |
| E10 | ruling | Does the VD entry still name the estimator without supplying it, now that composition forces the readout? (ours26 sweep §3.2) |

## Gaps (as of 2026-10-24; 27 Oct additions marked)

- (28 Oct) `tools/composition_coverage.py` is a seeded script, not a preset kind of `met run`
  (brief B03). (29 Oct) So is `tools/eta_sign_check.py` (T31). (30 Oct) So is `tools/ours31_coverage.py` (O18).

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
2. `../definition/ladder.md`: which commitment of the version removes which freedom, and
   what is left. Then `../definition/proof.md`: the chain from the type to ours30 with every
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

- 2026-10-30 (Codex, B06): wrote T16–T19, T24 and T25 proofs with all Statements unchanged. Added needed scope: T16 raw g=η=1 and negative-tilt attenuation; T17 fixed-N amplitude growth with positive g/η bounds, not generic limits 1; T18 unweighted Gaussian local asymptotic widths, since the exact mass-law mean is infinite; T25 matching covariance/candidate blocks or conditional Schur reduction. T19's tail and T24's zero-law proofs cover arbitrary covariance/subspace. All 25 new focused tests pass; full suite: 1,869 passed, 26 skipped, 16 xfailed; audit findings unchanged.
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
- 2026-10-30 (Claude, pressure): **what g and η are** (O14, provisional). Both are noise scales (tempering by w is
  declaring Σ/w); their caps read Σ as a floor, a declaration ours29/30 has run on unwritten (O7 gets the row).
  G13 sorted: the power is a scale read once (L21's cut, standard name), bounded confidence follows from the floor.
  On a common 4× inflation in 3D they over-correct together (g·η 0.165 vs 0.25); η alone gets it in 1D; g catches one
  bad reading. Declarations per level (exact / floor / unknown) are the input. `tests/test_scale_reading.py`.
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
