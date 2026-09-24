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
    law: {reference: flat, nuisance: radius}      # reference: flat | cartesian; nuisance: radius | acceleration | force
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

#### The `cartesian` reference (24 September)

`reference: cartesian` is the law of `excitation_weighting_premise_2026-09-24`:
flat in the true pair **as a vector**, d^d v dθ, with v the pair in noise
units and θ = atan(m/s) uniform. It is the flat measure times r^(d−2),
r² = (f/σ_F)² + (α/σ_a)², so it is identical to `flat` in 2D. Its radial
kernel is exp(h²/2) in every dimension, so a reading's law of mass is the
flat law's 2D form whatever d is (E[α | m], the noncentral chi mean, still
depends on d). The nuisance measures become d-dimensional volumes:

| nuisance | nuisance measure | implied reference factor on mass |
|---|---|---|
| `radius` | d^d v | uniform angle (as for `flat`) |
| `acceleration` | d^d a* = α^(d−1) dα dΩ | m · cos(θ)^(2−d) |
| `force` | d^d F* | (1/m) · sin(θ)^(2−d) |

So properness depends on d under this reference: `validate` checks it per
cell. The premise fixes the weighting only up to a factor in m; other
completions are this law times an `angle_power` prior (below), e.g.
`{sin: 0, cos: d−2}` for flat in m and `{sin: (d−2)/2, cos: (d−2)/2}` for the
σ-free (fα)^((d−2)/2) df dα dΩ.

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

### The sequential rule: old information in the prior slot

Under `likelihood_product` the exact law of N+1 readings factorises as
p_old(m) · L_new(m), and the readout's A(m) as A_old(m) + E[α_new | m]. So
the old readings' information enters exactly where a prior does. The
`sequential` rule treats the last reading of a series as new and the ones
before it as old, and `carry` says what of the old information is kept:

```yaml
combine:
  rule: sequential
  carry: exact            # exact | curve_only | lognormal_fit | tilt_fit
  old: {rule: likelihood_product, prior: [{uniform_angle: {center: first_reading}}]}
```

| carry | old information kept |
|---|---|
| `exact` | p_old and A_old (equals `old` applied to all readings; tested) |
| `curve_only` | p_old exactly; A_old dropped (median and intervals equal exact; tested) |
| `lognormal_fit` | two numbers: the mean and SD of log m under p_old, as a normal law in log m |
| `tilt_fit` | two numbers: m₀ and λ of sech(z)·exp(λ(sech z − 1)), z = log(m/m₀), matched to the same mean and SD |

The tilt family contains the zero-reading law (λ = 0, m₀ = s), and for large
λ it is close to normal with variance 1/(1+λ) in log m. So m₀ is the slot
for the old estimate and λ the slot for its precision. Add
`scores: {compare_to: exact}` to measure each carry against the exact law
series by series.

### Priors

A prior is a list of factors, multiplied together:
`flat_mass`, `flat_log_mass`, `flat_inverse_mass`,
`{uniform_angle: {center: C}}`, `{sech_tilt: {lambda: λ, center: C}}`,
`{lognormal: {center: C, width: w}}`,
`{angle_power: {sin: p, cos: q, center: C}}` (sin^p θ cos^q θ with
tan θ = m/C; `{sin: 1, cos: 1}` is `uniform_angle`). A centre `C` is a positive number or
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
- the `cartesian` reference (`tests/test_cartesian.py`): kernel and mean
  radius against radial quadrature; the whole law against brute force in the
  21 September coordinates with weight r^(d−2), and in 1D against brute force
  in the record's own variables (a*, m); the record's section 3 formulas term
  by term; 2D identical to `flat`; a 3D reading equal to three 1D readings;
  the σ-free completion by brute force; its d-dependent tail slopes;
- world designs and noise; per-series streams (a quick run is a prefix of the
  full run; grid position does not change a cell's data);
- presets refusing every missing required setting and each validation hole
  found in review; grid expansion; a reproducible end-to-end run.

## Not built yet

- Reference measures other than `flat` and `cartesian` (the tube/Jeffreys
  measure), and the hierarchical (finite-temperature) excitation prior.
- General or anisotropic covariance, correlated channels, and a known
  direction.
- Calibration estimated from samples, and calibration uncertainty.
- Readings from different instruments within one series (the law supports
  per-reading SDs; the world does not generate them yet).
- Paired-difference tests between estimators (the values are saved, so
  these can be added without rerunning).
