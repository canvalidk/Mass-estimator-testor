# Mass estimator tester — design

**Status:** first build done, 2026-09-24; extended in two parallel lines and
merged on 2026-09-25 (see "Since the first build" below, and README.md for
what exists and what is not built yet). The open questions at the end were
settled provisionally as marked; each is Can's to overturn.

## What this program is for

One program that runs every mass-estimator test. You choose a test by
choosing settings. Tests that share a stage share that stage's code. A test
suite is a recorded preset: a settings file you can run again with one
command.

It replaces the VD-Newton lab (`03_trace_and_evaluator/mass_estimator/`),
which we are not porting. That code is a reference for the maths and for what
went wrong, not for its structure.

## The central split: the world versus the analyst

Every study has two sides, and the program keeps them apart.

- **The world** generates data. It knows the true mass, how each reading was
  excited, the true noise, and which readings share what.
- **The analyst** is the estimator. It sees only the readings and what it is
  told (supplied covariances, a declared model). It never sees the world's
  settings, except through an oracle declared by name in the preset
  (`known_excitation`, `known_beta`). An oracle is a benchmark, never an
  estimator.

An estimator is a **declared model plus a readout**. A study scores declared
models against worlds.

Wherever both sides have a version of the same thing, they are separate
settings:

| Thing | World setting | Analyst setting |
|---|---|---|
| Noise | true covariance | supplied covariance |
| Shared structure across readings | the design: same pair, or new excitation per reading | the declared joint model |
| Mass | true mass, or a distribution of masses | the prior factor |

Getting the model right means the two sides match. Getting it wrong on
purpose is a legitimate study, but it must be written in the preset, never
arrived at by accident. (The 23 Sept repetition study analysed same-pair data
with a new-excitation model, because the design was not a setting of the
generator.)

## The pipeline

```
world ──► readings ──► [reduce] ──► per-reading law ──► combine ──► readouts ──► scores
```

1. **World (generator).** Produces series of readings. Settings:
   - `dimension`: 1, 2 or 3.
   - `mass`: a value, a list (grid), or a distribution.
   - `design`: `same_pair` (one latent vector pair, measured N times) or
     `new_excitation` (each reading draws its own latent pair with the same
     mass), plus the excitation spec (magnitudes or SNRs, direction rule).
   - `readings_per_series`: N.
   - `noise`: the true covariance per reading. Readings may come from
     different instruments.
   - `calibration`: how the supplied covariance relates to the true one
     (exact, scaled, estimated from k samples).
2. **Reduce (optional, analyst).** Transforms the readings before any law is
   formed. The one planned: `pool_pair`, the inverse-covariance mean of the
   readings (contribution 11, eq. 5). This is only valid if the analyst
   declares that the readings share one pair.
3. **Per-reading law (analyst).** For each (reduced) reading, returns
   **evidence about mass as curves on one common log-mass grid**, not a point:
   - `log_like(u)`: the likelihood of mass, with this reading's nuisance
     (magnitude and direction) integrated under the declared per-reading
     nuisance measure;
   - `log_ref(u)`: the mass factor the single-reading reference measure
     implies (for the flat law, the half-Cauchy at σ_F/σ_a);
   - `cond_alpha(u)`: E[α | m], needed by the ratio-of-means readout.

   Keeping the likelihood and the reference factor separate lets the
   combination stage count a prior once or N times on purpose.
   Settings: `reference` (`flat`, `tube`, …) and `nuisance`
   (`radius`, `acceleration`, …).
4. **Combine (analyst).** Joins the per-reading curves into one joint law.
   Settings: `rule` (`single`, `likelihood_product` = count the prior once,
   `posterior_product` = multiply whole single-reading laws, in mass
   (eq. 18) or log mass (eq. 20)) and `prior` (a list of factors: flat in m,
   log m or 1/m, uniform angle, sech tilt, log-normal). The prior is never
   implicit. (As built; the draft named the rules symmetric / mass_matching /
   log_matching. Those are now particular settings, listed in the README.)
5. **Readouts (analyst).** Scalars and intervals from the joint law:
   `ratio_of_means` (eq. 28: ∫ m A p / ∫ A p), `median`, `geometric`,
   `reciprocal_root`, quantile intervals, log SD. **Direct rules** that skip
   the law (vector-length ratio, dot-product ratios, total least squares)
   take the readings directly and are listed alongside.
6. **Scores.** Per cell: within factor K, capped squared log error, interval
   coverage and width. Across rules: paired differences, regret against the
   best rule in the cell, worst-case regret.

An **estimator** in a preset is a named path through stages 2–5. All of its
readouts come from one evaluation of its joint law. (Sharing the per-reading
curves between estimators that use the same law is a planned speed-up; the
first build computes them once per estimator.)

## Settings with no default

A preset must state these, or the program refuses to run it:

- world `design` and `mass`;
- world `noise` and the analyst's supplied covariance (even when identical);
- for every law-based estimator: `reference`, `nuisance`, `rule`, `prior`;
- `replicates` and `seed`;
- the scores to compute.

Only numerical settings have defaults (grid spacing, quadrature order), and
those are recorded in the output.

## Presets

A preset is one YAML file in `presets/`. See
`presets/repetition_same_pair.yaml` for a complete one, and the README for
every setting. Any world setting given as a list becomes a grid axis; each
combination is a cell. `python -m met run presets/<file>.yaml` writes a
results folder holding a copy of the preset, the code version, every
readout for every replicate, and the per-cell scores.

## Tests

Two kinds, kept apart:

- **Tests** check that the code computes what it claims, exactly where
  possible. The built list is in the README. Planned and not yet written:
  the 1D law against Liseo's eq. (12), the known-direction truncated-normal
  means (known direction is not built), and the composition identities (L)
  and (R) on the readout.
- **Studies** are presets. They answer questions about estimators. A study is
  never a test.

## First build

1. The per-reading law for independent isotropic noise per channel, in 1D
   and 3D (the exact J(h) reduction). 2D and general covariance come later.
2. The world generator with both designs.
3. `pool_pair`, the three combination rules, the four readouts and an
   interval.
4. The runner and the scores.
5. Two presets: the repetition study run both ways (same pair, new
   excitation), and a single-reading operating-range grid.

## Since the first build

Two sessions extended the first build in parallel from the same commit
("Document the sequential rule"). Their lines were merged on 2026-09-25 with
every commit kept:

- **Line A:** the `vonmises_state` carry; a `cartesian` reference with its
  comparison presets; the `hierarchical` rule and the Gaussian-excitation
  world; the declared oracles and the bottom-regime presets.
- **Line B:** the `cartesian` reference as the 24 September law (the r^(d−2)
  reweighting carried into the nuisance measures, so the acceleration and
  force splits become d-dimensional volumes); `angle_power` priors;
  `calibrate: sandwich` ("ours25") and its presets.

Both lines built a `cartesian` reference. They agree for one reading and for
the `radius` nuisance, the only one Line A's presets use. They differ when
readings are combined under the `acceleration` or `force` split. The merge
keeps Line B's law, which follows the measure it states, and keeps Line A's
independent checks of it in `tests/test_cartesian.py`. Every preset of both
lines gives bit-identical values before and after the merge.

## Open questions, settled provisionally in the first build

- **No-default list:** as above. In addition, a prior can only be omitted
  (`prior: []`) where the law already carries its reference (`single`,
  `posterior_product`); `likelihood_product` refuses an empty prior.
- **Units:** physical. The world states true SDs and a physical mass; s is
  derived. Presets that set both SDs to 1 read in units of s.
- **Where the prior lives:** in the combine stage, as a list of factors.
  The per-reading law still reports its implied reference factor separately,
  so "count the reference once" and "count it N times" are both explicit
  choices (`likelihood_product` vs `posterior_product`).
