# State of the mass estimator

**Last updated:** 2026-10-24, by Claude in a session with Can. Update this page in
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

**Where it stands.** The problem is close to fully mapped: the premises, what they
force, and what is left are in the lemma list. Two open items need a real idea, not
more tests (O2 and O6 below). Some premises are graded below bedrock or are
provisional (L12 medium; L14 and L15 provisional), and some choices are natural
rather than forced (eta's cap at 1 and its power form at finite N; κ = ν is fixed by
meaning, not by a property). The rest is rulings Can owes and proofs, some routine
(O4, O5) and some open research (O14, O15).

## The current equation: ours29

Inputs per reading i: measured force and acceleration (any values, including zero),
an error covariance with every declared sigma > 0, a declared candidate subspace,
and (since 10 Oct) a window. Declared once: the starting law pi with centre c.

- **Law.** p(m) ∝ pi(m) · exp( eta · Σᵢ gᵢ Hᵢ(m) / 2 ), where Hᵢ is reading i's
  explained energy in noise units, gᵢ = min(1, ν/Qᵢ) weights it by its own off-line
  misfit Qᵢ (ν = kᵢ − 1), and eta = min(1, curvature / J) is the Bartlett power, with
  J the sum of the readings' squared scores, both taken at the series' best line on
  the whole projective line of slopes (eta = 1 when J = 0). It measures the series'
  own disagreement about m.
- **Starting law.** Half-Cauchy, pi(m) = c / (c² + m²), with c = s = sigma_F / sigma_a
  by default (L3′, not yet ruled).
- **Readout.** m̂ = E[|F\*|] / E[|a\*|], the ratio of the mean sizes of the true force
  and acceleration over the candidates, with the series' size taken as the stacked
  length. Reported with the median, the quoted interval (eta applied) and the raw
  law's interval beside it.
- **Exact inputs** (sigma = 0) go to a separate three-way check: one mass, no mass,
  every mass. "Undefined" exists only there, at exact 0/0.

Defining record: `noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md`.
Premises and theorems: `vd-docs-context/mass_estimator_lemmas_2026-10-10.md` (T8 is
the law, T10 the readout).

**Run it.** In `met/` (isotropic or per-axis noise only), from the repo root:

```python
from met.analyst import estimate

LAW = {"reference": "cartesian2", "nuisance": "radius"}
OURS29 = {"rule": "likelihood_product",
          "prior": [{"uniform_angle": {"center": "first_reading"}}],
          "misfit": "kink", "calibrate": "bartlett"}

# four 3D readings that disagree about m; sigma_F = 0.2 N, sigma_a = 0.1 m/s^2
F = [[4.0, 0.3, 0.0], [6.1, -0.2, 0.1], [1.0, 0.8, 0.0], [3.0, 0.1, -0.2]]
a = [[2.0, 0.0, 0.1], [2.0, 0.1, 0.0], [0.6, -0.4, 0.0], [1.0, 0.0, 0.0]]
estimate(F, a, 0.2, 0.1, ["ratio_of_means_rss", "median", "interval_95"],
         law=LAW, combine=OURS29)
# ratio_of_means_rss 2.6454, median 2.6454, interval_95 [2.0914, 3.4009]
# without misfit and calibrate: 2.6262, 2.6262, [2.4203, 2.8549];
# the tempering widens the interval because the readings disagree
# all-zero readings return s = 0.2 / 0.1 = 2 (checked 2026-10-24)
```

The returned dict also carries an `unresolved` flag. The full input type (any
covariance, a declared subspace) exists only in the ours29 reference implementation,
`noise_vs_signal/ours29_reference.py` (byte copy in `../reference/`). Nothing takes a
window yet.

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

### O2: the readout

Find a readout that satisfies all three of these at once:

1. **Composition.** Masses add when bodies are lumped (L5), and inverse masses add
   under a common force (L14, provisional). With C3 (T21), composition forces a
   ratio of means of the series' two sizes (T4). The open freedom is what "the
   series' size" is.
2. **Z at the point.** Appending a reading that carries no information leaves the
   answer unchanged (L20, ruled 27 Sept).
3. **A value at the all-zero experiment.** Ruled 29 Sept: on measured inputs mass
   is never undefined, so the readout must give a value there. That the value must
   be c is derived, not ruled: it follows from Z, units and the swap if the point
   keeps the swap property (`vd-docs-context/why_mass_cannot_be_undefined_2026-09-29.md`
   §4, offered for O2, not placed).

Each known candidate fails at least one:

| Candidate | Composition | Z | All-zero value c |
| --- | --- | --- | --- |
| Cloud readout (current, T10): sizes averaged over each reading's candidate cloud | yes | **no**: 50 appended null readings move it 2.2524 → 2.2406 (T29: ours28, 3D); under ours29 ten move it 3.5983 → 3.6017 (O1 row). It moves at c = s too | only when c = s: at c = 3s it gives 2.095 s (ours26 sweep, E1) |
| Best-guess repair: size of each reading's conditional-mean candidate | unknown | yes | **no**: 0/0 there; limits near zero depend on direction, 0.546 to 1.830 (1D, c = 1) |
| The law's median | **no**: it uses only the size law and meets C3–C5, so by the uniqueness theorem it fails lumping or common force | yes (the law is Z-exact, T24) | yes |

A breakthrough is either a readout that meets all three or a proof that none exists,
in which case a ruling picks the requirement to drop. Unanswered: must the readout
also be continuous at the all-zero experiment?

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

| Item | Kind | One line |
| --- | --- | --- |
| O4 | proof | Gaussian noise forced under the tempered evidence |
| O5 | proof | Write out the betweenness representation step behind T21 |
| O7 | ruling | Does Z cover zeros beside a signal, or only whole null readings? |
| O9 | ruling | Is the window the criterion for one reading versus several? (proposed 10 Oct) |
| O11 | ruling | Does the default centre reference the instruments only (L3′), and is the swap-of-declarations premise accepted, which gives c = s for one noise pair? Then the several-instrument and anisotropic cases |
| O13 | ruling | Is anything more owed for readings that contradict each other? |
| O14 | proof | eta at small N: a correction that counts effective readings |
| O15 | proof | Is the bounded null floor the best any law can do? |
| L13, L18 | ruling | Confirm "the noise model is forced Gaussian, the data need not be"; is the subspace an input? |
| E10 | ruling | Does the VD entry still name the estimator without supplying it, now that composition forces the readout? (ours26 sweep §3.2) |

## Gaps (as of 2026-10-24)

- The lemma list d4 still shows O1 open; the 29 Sept ruling closed it. The living layer
  records O1 as closed; d5 of the lemma list, when Can wants one, is generated from
  `definition/` (`python tools/audit.py --table`).
- 20 of the 33 theorems name no test, or only a partial one (`python tools/audit.py`).
  Several are derivations with nothing to test (T2, T3, T21); the rest are work.
- L3 (units covariance) is used by T14's argument but no theorem's From line cites it.
- The current equation note in VD-docs states ours25b, four equation states behind
  (ours26 to ours29), and its many-reading step rests on the strong form of axis
  splitting, ruled false on 29 Sept (L8).
- The ours29 record says "eta never moves the law's peak". With pi included the code
  moves it (ln m 1.1181 with eta 0.403, 1.2181 with eta 1). The likelihood part alone
  peaks at ln m 1.3181 for both, so the peak moves only because eta changes pi's
  relative weight: the record's sentence is wrong, not the code. Strict xfail.
- ours29 has never been run on real data. The robot data (Zenodo
  10.5281/zenodo.11096791, CC-BY 4.0) is in `../legacy/vd-newton-lab/raw/robot/`;
  `met/` has no file loader.
- One golden-output test fails under a newer numpy, on a number that is rounding noise
  around zero (`ours29_reference_checks.py`, token 103: +3.3e-08 against −3.8e-08).
  The three strict rotation xfails in `tests/test_reference_records_ours29.py` sit near
  their 3e-8 bound for the same reason and could pass, and so fail, on another numpy.
- The 10 Oct check script behind T32 is not yet kept (pending since 10 Oct).
- The reproducibility inventory still cites VD-Newton paths for scripts now in `../legacy/`.

## Read first

1. This page.
2. `../definition/equation.md`, then `../definition/premises.md` and
   `../definition/theorems/README.md`: what is claimed and what it rests on.
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

- 2026-10-24: page created from a full map of the tester, VD-docs and VD-Newton, then
  corrected after an independent fact-check (c = s qualifier on the cloud readout; the
  value c derived, not ruled; eta's J; the example replaced by one where the tempering
  acts). Tests: 1,044 pass, 1 fails (golden output, see Gaps), 12 strict xfails.
