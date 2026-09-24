# Mass estimator tester

One program for every mass-estimator test. You choose a test by choosing
settings; a recorded study is a preset file; tests that share a stage share
its code. The design and its reasons are in [DESIGN.md](DESIGN.md).

```sh
python -m pip install -r requirements.txt
python -m met list                                   # presets and their questions
python -m met validate presets/repetition_same_pair.yaml
python -m met run presets/repetition_same_pair.yaml  # full run
python -m met run presets/repetition_same_pair.yaml --replicates 50   # quick look
python -m pytest -q                                  # the exact tests
```

A run writes `results/<preset>/<timestamp>/` (not committed):

- `summary.md`: per-cell tables and the worst-case regret table;
- `run.json`: every score, the world of every cell, timings, numerical
  settings, and the code version (git commit plus per-file hashes);
- `values.npz`: every readout for every replicate, keyed
  `cell/estimator/readout`, so new scores can be computed without rerunning;
- `preset.yaml`: the exact preset that was run.

A run is reproducible from its preset alone (see Numerics for the seeding).
All estimators in a cell see the same readings, so their scores are paired.

## The pipeline

```
world ──► readings ──► [reduce] ──► per-reading law ──► combine ──► readouts ──► scores
 world.py                 analyst.py     law.py          analyst.py   analyst.py    scores.py
```

**The world never meets the analyst.** The world knows the true mass, the
excitation, the true noise and which readings share a latent pair. An
estimator sees only the readings and the standard deviations it was
supplied.

| Module | Stage |
|---|---|
| `met/world.py` | World settings, validation, data generation (`same_pair`, `new_excitation`) |
| `met/law.py` | The per-reading law: closed-form radial kernels for d = 1, 2, 3, the three nuisance splits, E[α \| m], and `pool_pair` |
| `met/analyst.py` | Priors, combination rules, the adaptive log-mass grid, readouts, direct rules, the `Estimator` class |
| `met/scores.py` | Within-factor success, capped log error, interval coverage, regret |
| `met/preset.py` | Preset validation and grid expansion |
| `met/runner.py` | Runs a preset and writes the results folder |

## Writing a preset

Every setting below is required unless marked optional. There are no
defaults for anything that changes what is being modelled.

```yaml
name: my-study
question: The one question this study answers.
notes: optional free text
world:
  dimension: 3                 # 1, 2 or 3
  design: same_pair            # or new_excitation
  mass: [0.25, 1, 4]           # physical units; a list makes a grid axis
  acceleration_snr: [1, 3]     # |a*| / sigma_a; 0 = no excitation; or {uniform: [lo, hi]}
  direction: fixed             # or random (drawn per latent pair)
  readings: 50                 # per series
  noise: {force_sd: 1, acceleration_sd: 1}      # true per-coordinate SDs
  supplied_noise: exact        # or {force_scale: k, acceleration_scale: k}
estimators:
  - name: symmetric
    reduce: none               # or pool_pair
    law: {reference: flat, nuisance: radius}      # nuisance: radius | acceleration | force
    combine:
      rule: likelihood_product # single | likelihood_product | posterior_product
      prior: [{uniform_angle: {center: first_reading}}]
    readouts: [ratio_of_means, median, interval_95]
  - name: vector_length
    direct: norm_ratio         # norm_ratio | dot_acceleration | dot_force
    reduce: pool_pair          # direct rules act on one pair; this must be stated
scores:
  within_factor: [1.25, 1.5, 2]
replicates: 1000
seed: 1
```

Any world setting given as a list becomes a grid axis, including
`noise.force_sd`, `noise.acceleration_sd` and the `supplied_noise` scales.
Each combination of axis values is a cell.

`validate` (and `run`, before anything is computed) refuses a preset when,
in any cell:

- a law cannot be normalised. Whether it can depends only on the rule, the
  nuisance, the prior and the number of readings, never on the data: the
  radial kernel tends to a constant at both ends, so the log density's tail
  slopes in log m are fixed by the settings (`tail_slopes` in
  `met/analyst.py`, checked against the computed density in the tests);
- `rule: single` meets more than one reading;
- an estimator pools readings (`reduce: pool_pair`, which every direct rule
  does) in a `new_excitation` world without a `declared_mismatch: <why>`;
- a `new_excitation` world cannot vary (fixed direction, one SNR), which
  would silently be the same pair;
- a numeric or prior parameter is missing, mistyped or out of range.

### The per-reading law

For one reading the flat law is
`dP ∝ exp(−Q/2) df dα dΩ` (flat in the two positive magnitudes, uniform
in direction). The law module integrates direction and radius in closed
form (see the docstring of `met/law.py`) and returns, on a common log-mass
grid u = log m, three curves per reading:

- `log_like(u)`: the likelihood of mass with the reading's nuisance
  integrated;
- `log_ref(u)`: the mass factor the single-reading reference implies;
- `cond_alpha(u)`: E[α | m].

`nuisance` chooses how the single-reading law is split into likelihood and
reference factor. The product is the same for one reading. The split
matters only when readings are combined:

| nuisance | nuisance measure | implied reference factor on mass |
|---|---|---|
| `radius` | r dr dΩ | uniform angle: half-Cauchy at s = σ_F/σ_a |
| `acceleration` | α dα dΩ | flat in m |
| `force` | f df dΩ | flat in 1/m |

### Combination rules

| rule | joint log density | notes |
|---|---|---|
| `single` | log_like + log_ref (+ prior factors) | exactly one reading after reduce |
| `likelihood_product` | Σ log_like_i + prior | the prior is counted once and must be stated |
| `posterior_product` with `coordinate: mass` | Σ (log_like_i + log_ref_i) + (1 − N) u | contribution 11, eq. (18); the nuisance setting does not matter |
| `posterior_product` with `coordinate: log_mass` | Σ (log_like_i + log_ref_i) | contribution 11, eq. (20) |

Named combinations from the earlier work:

- contribution 12's **symmetric** rule: `likelihood_product`, nuisance
  `radius`, prior `[{uniform_angle: {center: first_reading}}]`;
- **eq. (18)** is also `likelihood_product` with nuisance `acceleration` and
  prior `[flat_mass]`, exactly (tested).

### Priors

A prior is a list of factors, multiplied together:
`flat_mass`, `flat_log_mass`, `flat_inverse_mass`,
`{uniform_angle: {center: C}}`, `{sech_tilt: {lambda: λ, center: C}}`,
`{lognormal: {center: C, width: w}}`. A centre `C` is a positive number or
`first_reading` (the instrument ratio s of the series' first reading).
`likelihood_product` refuses an empty prior. For `single` and
`posterior_product`, `[]` means no factor beyond the reference already in
the law.

A law that cannot be normalised (for example `likelihood_product` with
nuisance `radius` and prior `[flat_mass]`) is refused at validation. As a
safety net, a series whose law still runs off the search grid is reported as
unresolved: every readout is NaN and counts as a failure.

### Readouts

`ratio_of_means` is E[Σf]/E[Σα] under the joint law, computed as
∫ m A(m) p(m) dm / ∫ A(m) p(m) dm with A = Σ_i E[α_i | m] (contribution 11,
eq. 28). For one reading it is E[f]/E[α]. Also `median`, `geometric`
(exp E[log m]), `reciprocal_root` (E[√m]/E[1/√m]), `log_sd`, and equal-tail
intervals `interval_<pct>` (e.g. `interval_95`).

Direct rules (`norm_ratio`, `dot_acceleration`, `dot_force`) skip the law.
They are applied to the inverse-variance pooled pair of the series. An
invalid value (a non-positive dot product) is NaN and counts as a failure.

## Numerics

The joint law is found on a coarse log-mass grid (±60 around log s, step
0.05) and cut where it falls e⁻⁴⁰ below its peak. If the support spans fewer
than 100 coarse steps, the search is repeated on 401 points over it, until it
does (so very precise data are resolved: tested to SNR 10⁶). The law is then
evaluated on 2001 uniform points over the support. Moments use the trapezoid rule, which is exponentially
accurate for these smooth, fully decayed integrands. Quantiles use a
fourth-order CDF (trapezoid with the Euler–Maclaurin endpoint correction)
inverted by cubic Hermite interpolation. At the default grid, the zero-reading
95% interval is exact to 5×10⁻⁸ and the equation README example to
3×10⁻¹¹. Grid settings can be changed per estimator with
`numerics: {grid_points, cutoff, span, coarse_step, min_support_steps,
refine_points}`; unknown keys are refused, and the values used are recorded
in the output.

Each series has its own random stream,
`SeedSequence(seed, spawn_key=(cell_key, i))`, with `cell_key` a hash of
the cell's world settings. So `--replicates 50` sees exactly the first 50
series of the full run, and adding a value to a grid axis does not change
the data of the other cells.

## Tests

`tests/` holds exact checks of the code, never studies:

- the closed-form kernels and mean radius against direct numerical
  integration (d = 1, 2, 3);
- the three nuisance splits give the same single-reading law;
- the zero-reading law is the half-Cauchy at s (point, median, log SD π/2,
  95% interval) in every dimension;
- the equation README example (21 Sept) and contribution 11's single-trial,
  combined and scaled values;
- reciprocal symmetry, unit change, rotation invariance, grid convergence;
- all rules agree for one reading; eq. (18) equals the acceleration-nuisance
  product with a flat-mass prior; contribution 12's angle identities;
- zero-reading combinations (s/2 for eq. 18, s for eq. 20, the symmetric
  rule unchanged);
- contribution 12's tilted null log SDs;
- the static tail-slope table against the computed density, for every
  nuisance × rule × prior combination; improper laws refused at validation;
- very precise data (SNR up to 10⁶) resolved: log SD and interval against
  the normal limit;
- pooling against the N-reading same-pair law integrated by brute force;
- world designs and noise; per-series streams (a quick run is a prefix of the
  full run; grid position does not change a cell's data);
- presets refusing every missing required setting and each validation hole
  found in review; grid expansion; a reproducible end-to-end run.

## Not built yet

- Reference measures other than `flat` (the tube/Jeffreys measure).
- General or anisotropic covariance, correlated channels, and a known
  direction.
- Calibration estimated from samples, and calibration uncertainty.
- Readings from different instruments within one series (the law supports
  per-reading SDs; the world does not generate them yet).
- Paired-difference tests between estimators (the values are saved, so
  these can be added without rerunning).
