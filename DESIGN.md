# Mass estimator tester — design

**Status:** draft for review, 2026-09-24. Nothing below is built yet.

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
  settings.

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
   Settings: `rule` (`single`, `symmetric`, `mass_matching` = eq. 18,
   `log_matching` = eq. 20) and `prior` (`none`, a sech tilt with λ and m₀,
   log-normal …). A prior is always a named setting, including `none`.
5. **Readouts (analyst).** Scalars and intervals from the joint law:
   `ratio_of_means` (eq. 28: ∫ m A p / ∫ A p), `median`, `geometric`,
   `reciprocal_root`, quantile intervals, log SD. **Direct rules** that skip
   the law (vector-length ratio, dot-product ratios, total least squares)
   take the readings directly and are listed alongside.
6. **Scores.** Per cell: within factor K, capped squared log error, interval
   coverage and width. Across rules: paired differences, regret against the
   best rule in the cell, worst-case regret.

An **estimator** in a preset is a named path through stages 2–5. Estimators
that share a law share its computation: the curves are computed once per
reading and every readout reuses them.

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

A preset is one YAML file in `presets/`:

```yaml
name: repetition-same-pair
question: Does each combination rule recover the mass when one pair is measured 50 times?
world:
  dimension: 3
  design: same_pair
  excitation: {acceleration_snr: [1, 3], force_snr_from_mass: true}
  mass: [0.25, 1, 4]            # in units of s = σ_F/σ_a
  readings_per_series: 50
  noise: {force_sd: 1, acceleration_sd: 1}
  calibration: exact
estimators:
  - name: pooled
    reduce: pool_pair
    law: {reference: flat, nuisance: radius}
    combine: {rule: single, prior: none}
    readouts: [ratio_of_means, median, interval_95]
  - name: symmetric
    law: {reference: flat, nuisance: radius}
    combine: {rule: symmetric, prior: first_reading}
    readouts: [ratio_of_means, median, interval_95]
scores: [within_factor: [1.25, 1.5, 2], coverage: [0.8, 0.95]]
replicates: 2000
seed: 20260924
```

Any world setting given as a list becomes a grid axis; each combination is a
cell. `python -m met run presets/repetition-same-pair.yaml` writes a results
folder holding a copy of the preset, the code version, and the per-cell
scores.

## Tests

Two kinds, kept apart:

- **Tests** check that the code computes what it claims. They are exact where
  possible: the zero-reading null law (half-Cauchy, log SD π/2, the 95%
  interval [0.0393, 25.45]·s), reciprocal symmetry, unit covariance, the
  known-direction truncated-normal means, the 1D law against Liseo's eq. (12),
  the eq. (18)/(20) angle identities against the symmetric rule, and the
  composition identities (L) and (R) on the readout. One spot-check against
  the old lab's core law, on a handful of inputs, to catch a mistake on
  either side.
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

## Open, for Can to rule

- The no-default list above: anything to add or remove?
- Should masses and SNRs be stated in units of s = σ_F/σ_a (as the old lab
  did), or in physical units with s derived?
- Is a prior a property of the combination stage (as here), or its own stage?
