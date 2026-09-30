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
| `met/world.py` | World settings, validation, data generation (`same_pair`, `new_excitation`; fixed, uniform or Gaussian excitation) |
| `met/law.py` | The per-reading law: closed-form radial kernels for d = 1, 2, 3, the three nuisance splits, E[α \| m], and `pool_pair` |
| `met/analyst.py` | Priors, combination rules (including `sequential`, `hierarchical` and `calibrate: sandwich`), the adaptive log-mass grid, readouts, direct rules, declared oracles, the `Estimator` class |
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
  acceleration_snr: [1, 3]     # |a*| / sigma_a; 0 = no excitation; or {uniform: [lo, hi]},
                               # or {normal_rms: v} (a* ~ N(0, tau^2 I), RMS length v sigma_a; needs direction: random)
  direction: fixed             # or random (drawn per latent pair)
  readings: 50                 # per series
  noise: {force_sd: 1, acceleration_sd: 1}      # true per-coordinate SDs; either may be
                               # {per_axis: [s_1, ..., s_d]} (ours26: one SD per axis)
  supplied_noise: exact        # or {force_scale: k, acceleration_scale: k}
estimators:
  - name: symmetric
    reduce: none               # or pool_pair, or split_axes (each d-dim reading as d 1D readings)
    law: {reference: flat, nuisance: radius}      # reference: flat | cartesian1 | cartesian2; nuisance: radius | acceleration | force
    combine:
      rule: likelihood_product # single | likelihood_product | posterior_product | sequential | hierarchical
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

#### The two cartesian references: `cartesian1` and `cartesian2` (24 September)

Both implement the premise of `excitation_weighting_premise_2026-09-24`:
flat in the true pair **as a vector**, d^d v dθ, with v the pair in noise
units and θ = atan(m/s) uniform. That is the flat measure times r^(d−2),
r² = (f/σ_F)² + (α/σ_a)², so it is identical to `flat` in 2D. Its radial
kernel is exp(h²/2) in every dimension, so a reading's law of mass is the
flat law's 2D form whatever d is (E[α | m], the noncentral chi mean, still
depends on d).

Two sessions wrote this independently on the same day, and both versions are
kept so they can be checked against each other:

- `cartesian1`: Line A, commit 665a321 (18:13). Used by the `hierarchical`
  rule and Line A's presets.
- `cartesian2`: Line B, commit 7081a44 (23:27). Used by Line B's presets and
  the ours24/ours25 work.

The two give the same single-reading law, and the same combined law under the
`radius` nuisance. They differ only in how the `acceleration` and `force`
nuisances split a reading's law when readings are combined:

| nuisance | `cartesian1` nuisance measure | `cartesian2` nuisance measure | `cartesian2` reference factor on mass |
|---|---|---|---|
| `radius` | d^d v | d^d v | uniform angle (as for `flat`) |
| `acceleration` | α dα dΩ (as for `flat`) | d^d a* = α^(d−1) dα dΩ | m · cos(θ)^(2−d) |
| `force` | f df dΩ (as for `flat`) | d^d F* | (1/m) · sin(θ)^(2−d) |

Per reading, `cartesian1`'s likelihood is `cartesian2`'s times cos(θ)^(2−d)
(acceleration) or sin(θ)^(2−d) (force). `cartesian1` keeps `flat`'s reference
factors (flat in m, flat in 1/m). `cartesian2`'s splits are the integrals over
the nuisance under the stated measure, so its properness depends on d, and
`validate` checks it per cell.

The plain name `cartesian` is refused with a message naming both. Recorded
work from before the split (24–25 September) says `cartesian` and cites the
commit it ran at, which still reproduces it. Every study run so far used the
`radius` nuisance, and those studies give bit-identical values under either
reference.

The premise fixes the weighting only up to a factor in m. Other completions
are `cartesian2` times an `angle_power` prior (below): for example,
`{sin: 0, cos: d−2}` gives flat in m, and `{sin: (d−2)/2, cos: (d−2)/2}`
gives the σ-free (fα)^((d−2)/2) df dα dΩ.

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

### Calibrated width (`calibrate: sandwich`, "ours25")

With many readings sharing only the mass, the declared law is too narrow:
its curvature H understates the spread of its own peak, whose variance is
J/H², with J the variance of the summed per-reading score. For the
cartesian law (`radius` nuisance), per reading, J = d + r*² and H = r*² exactly (r* is the true
pair's noise-unit length), so the law is too narrow by 1 + d/r*² in variance.
`likelihood_product` accepts

```yaml
combine: {rule: likelihood_product, prior: [...], calibrate: sandwich}
```

which raises the summed likelihood to w = H/J, estimated per series from the
readings at the likelihood's peak (centred scores), clipped to [10⁻⁶, 1]
(never narrower than the declared law; with no curvature the law falls back
to the prior). The prior is not tempered. It needs at least two readings and
the `radius` nuisance. It is a large-N correction for each fixed, nonzero
excitation. No bounded interval can hold its coverage uniformly as the
excitation vanishes (Gleser & Hwang 1987); see
`presets/ours25_identification_limit.yaml`.

### ours29: the misfit weight and the Bartlett power (`misfit`, `calibrate: bartlett`, 30 September)

ours29 is the cartesian likelihood product (radius nuisance, prior counted once)
with two tempering steps, both computed from the readings alone:

```yaml
law: {reference: cartesian2, nuisance: radius}      # or cartesian1
combine:
  rule: likelihood_product
  prior: [{uniform_angle: {center: first_reading}}]
  misfit: kink          # or smooth; each reading's log-likelihood times g_i
  calibrate: bartlett   # the summed log-likelihood raised to eta
```

- **`misfit`** ("the readings are evidence, just less", 29 September): reading
  i's log-likelihood is multiplied by g_i = g(Q_i), where Q_i is its off-line
  misfit at its own best line (the smaller eigenvalue of its whitened Gram
  matrix, (T_i − |C_i|)/2) and ν = d − 1. `kink`: g = min(1, ν/Q);
  `smooth`: g = 1/(1 + Q/ν). g = 1 in 1D and for Q ≤ ν (kink). A reading no
  mass can explain is worth a bounded amount, ν a/(1 − a) nats with
  a = |C_i|/T_i, however small its noise.
- **`calibrate: bartlett`**: the summed (misfit-weighted) log-likelihood is
  raised to η = min(1, curv/J), curv = −Σ l_i'', J = Σ l_i'², at the series'
  best line on the **whole projective line** of slopes (negative and infinite m
  included, so every series has one). J is the raw sum: an exact zero reading
  has l_i' = l_i'' = 0 and leaves η alone (Z). ours25's `sandwich` differs in
  two places, both deliberate: its J carries the Bessel factor N/(N − 1), which
  counts readings and so moves when zero readings are appended, and it searches
  m > 0 only, so when a series' best line is at negative mass it falls back to
  the prior (w = 10⁻⁶; 3 of 40 pure-noise 3D series of 12 readings in one check). One
  reading, or no reading that disagrees (J = 0), gives η = 1. η never exceeds 1.

Neither step changes a reading's candidate cloud, so E[α | m] and the stacked
weight of `ratio_of_means_rss` are unchanged (shrinking the candidates by η
would infer signal size from the series). With one noise ratio per series
(isotropic noise) everything is closed form in C_i = (Q_i − P_i) + 2iD_i: the
best line is arg(Σ g_i C_i)/2, curv = |Σ g_i C_i|, l_i' = −½ g_i Im(C̄_i e^{2iθ̂});
otherwise (per-axis noise, mixed ratios) both are found numerically on the
projective line (golden section after a 721-point scan). The two paths agree to
about 10⁻⁸, and both agree with a 200 001-point brute-force grid and with the
reference implementation in VD-docs (`noise_vs_signal/ours29_reference.py`).

To report both the raw law and the quoted width (his ruling of 30 September),
run ours29 twice in one preset, with and without `calibrate: bartlett`; the
estimators see the same readings. See `presets/ours29_noise_vs_signal.yaml`.

### The sequential rule: old information in the prior slot

Under `likelihood_product` the exact law of N+1 readings factorises as
p_old(m) · L_new(m), and the readout's A(m) as A_old(m) + E[α_new | m]. So
the old readings' information enters exactly where a prior does. The
`sequential` rule treats the last reading of a series as new and the ones
before it as old, and `carry` says what of the old information is kept:

```yaml
combine:
  rule: sequential
  carry: exact            # exact | curve_only | lognormal_fit | tilt_fit | vonmises_state | stacked_state
  old: {rule: likelihood_product, prior: [{uniform_angle: {center: first_reading}}]}
```

| carry | old information kept |
|---|---|
| `exact` | p_old and A_old (equals `old` applied to all readings; tested) |
| `curve_only` | p_old exactly; A_old dropped (median and intervals equal exact; tested) |
| `lognormal_fit` | two numbers: the mean and SD of log m under p_old, as a normal law in log m |
| `tilt_fit` | two numbers: m₀ and λ of sech(z)·exp(λ(sech z − 1)), z = log(m/m₀), matched to the same mean and SD |
| `vonmises_state` | one doubled-angle vector: each old likelihood replaced by its von Mises part ((Q−P)/4) cos 2θ + (D/2) sin 2θ, m = s tan θ |
| `stacked_state` | ours25b's exact state, three running sums (C, T, M): the von Mises vector, the total energy T = Σ(P+Q), the component count M. Carries the law and `ratio_of_means_rss` exactly in every dimension (tested) |

The tilt family contains the zero-reading law (λ = 0, m₀ = s), and for large
λ it is close to normal with variance 1/(1+λ) in log m. So m₀ is the slot
for the old estimate and λ the slot for its precision. Add
`scores: {compare_to: exact}` to measure each carry against the exact law
series by series.

`vonmises_state` needs the old rule to be `likelihood_product` with nuisance
`radius`. With a common instrument ratio s, all the old readings reduce to one
vector Σ c_i, c_i = ¼ Σ_k (y_k + i x_k)². The state is exact under the
cartesian reference (either one) in every dimension and under `flat` in 2D. Under `flat`
in 1D and 3D it drops a factor of about h^(2−d) (tested both ways). Its peak
is the total-least-squares direction of the standardised scatter.

`stacked_state` (26 September) is the combination equation for ours25b. Under
either cartesian reference, with the radius nuisance, the prior counted once,
and one s and one σ_a per series, everything ours25b computes depends on the
data only through

    H(θ)² = T/2 + 2 Re[C e^(−2iθ)],   C = Σ ((Q−P) + 2iD)/4,   T = Σ (P+Q),

and the component count M = N·d: the log-likelihood is H²/2 and the readout's
weight is σ_a cos θ · μ_M(H). Combining sets of readings is adding their
(C, T, M), in any order and grouping (`met.law.StackedState`). It needs the
readout `ratio_of_means_rss`; it does not carry eq. 28's A_old, so the plain
`ratio_of_means` is refused with it. See `presets/old_info_stacked_state.yaml`.

### The hierarchical rule: a learned excitation scale

`rule: hierarchical` gives every reading's standardised latent vector a shared
prior N(0, ω² I) instead of a flat one. With n = N d and β = ω²/(1+ω²),
integrating every latent vector out leaves (1 − β)^(n/2) exp(β H(θ)/2), with
H = Σ_i h_i(θ)². The `cartesian1` reference is the β → 1 limit. β ~ Beta(1, b)
is integrated out exactly (a lower incomplete gamma function), so the
excitation scale is learned from the readings:

```yaml
law: {reference: cartesian1, nuisance: radius}    # required by this rule
combine:
  rule: hierarchical
  prior: [{uniform_angle: {center: first_reading}}]
  excitation: {beta_b: 1}                          # Beta(1, b) prior on beta, b > 0
```

E[α_i | m] averages over β's conditional law at 24 quantile midpoints. The
matching world is `acceleration_snr: {normal_rms: v}` with `direction: random`.
Its excitation law matches the rule's prior.

**With ours26 and the stacked readout (27 September).** The rule also takes
per-axis noise and the readout `ratio_of_means_rss`. Each component's push is
standardised in its own noise pair's units and all share one ω (isotropy in
noise units). Then, at fixed β, the law is ours26's likelihood raised to the
power β (the tempering identity): its peak is ours26's for every β, and it is
wider by 1/β. Learned, β ≈ R*²/(M + R*²), so 1/β ≈ 1 + M/R*², the factor by
which ours26's law is too narrow. Given θ and β the stacked push is
N(βW, βI), so the stacked readout's weight is E|a*| of that, averaged over
β | θ by tanh-sinh nodes in β's conditional quantile (about 1e-11 against
direct quadrature). The plain-sum `ratio_of_means` is not built for per-axis
noise, as for ours26. Tests: `tests/test_hierarchical_ours26.py` (the
tempering identity, the fixed peak, the law against brute-force integration
with per-axis noise, the shrunk weight against Monte Carlo, the β average
against quadrature, β → 1 gives ours26, axis splitting, the swap, units).
Study: `presets/ours26_learned_omega.yaml`.

### Declared oracles

An oracle is the one place an estimator is shown part of the truth. It is only
ever declared by name (`oracle: ...` on the estimator), so it cannot happen
by accident:

| oracle | told | role |
|---|---|---|
| `known_excitation` | every reading's true acceleration vector a*_i | a ceiling: the law is a normal in m truncated to m > 0. Takes exactly `name`, `oracle`, `readouts` (`ratio_of_means`, `median`, intervals). |
| `known_beta` | the true excitation scale β, from the world's settings | the `hierarchical` rule with β plugged in instead of learned. In the Gaussian-excitation world it is the best any method can do with the same readings and mass prior. |

See `presets/bottom_single_reading.yaml` and `presets/bottom_repeated.yaml`.

### Priors

A prior is a list of factors, multiplied together:
`flat_mass`, `flat_log_mass`, `flat_inverse_mass`,
`{uniform_angle: {center: C}}`, `{sech_tilt: {lambda: λ, center: C}}`,
`{lognormal: {center: C, width: w}}`,
`{angle_power: {sin: p, cos: q, center: C}}` (sin^p θ cos^q θ with
tan θ = m/C; `{sin: 1, cos: 1}` is `uniform_angle`). A centre `C` is a positive number or
`first_reading` (the instrument ratio s of the series' first reading), or
`first_reading_geometric` (the geometric mean of the first reading's axis ratios
σ_F,k/σ_a,k; the same as `first_reading` when the axes share one ratio). With
per-axis noise whose axes have different ratios, `first_reading` is refused: no
single ratio is singled out, so the centre is not forced and must be declared
(ours26, below). After `split_axes` the first reading is still the original
one, so both spellings mean the same centre whether or not the axes are split.
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
eq. 28). For one reading it is E[f]/E[α].

`ratio_of_means_rss` ("ours25b", 25 September) pools the readings' sizes as the
length of the stacked vector instead of their plain sum: E[|F*|]/E[|a*|], with
F*, a* the N readings' true vectors stacked into one of N·d components, so
A = E[|a*| | m] = σ_a cos θ · E|N(W, I_{Nd})| with |W|² = Σ_i h_i² (the
noncentral chi mean in N·d dimensions, `stacked_mean_radius` in `met/law.py`).
The pooling happens on each candidate, before the average; pooling the averages
instead would break lumping. Why this pooling: axis splitting applied to the
answer (a 2D reading and its two axes must have the same size), given one rule
for joining sets of readings, forces √(p² + q²); the plain sum fails. For one
reading it equals `ratio_of_means`; the medians and intervals are unchanged.
Consequences, all tested: a d-dimensional series and its N·d axis readings
give the same answer to rounding; the recording frame does not matter; and
under `cartesian1`/`cartesian2` N readings in d dimensions are exactly one
reading in N·d dimensions, so the answer depends on the data only through
Σ|x_i|², Σ|y_i|², Σ x_i·y_i (noise units) and N·d. It needs a cartesian
reference and rule `single`, `likelihood_product` or `posterior_product`. See
`presets/ours25b_axis_splitting.yaml` and `presets/ours24_vs_ours25b.yaml`.

**Per-axis noise ("ours26", 26 September).** A world may give each axis its own
SD (`noise: {acceleration_sd: {per_axis: [...]}}`), and a reading set may carry
SDs per axis, shape (B, N, d). ours26 is ours25b on such data, with nothing
changed where ours25b was defined (tested to 1e-12 on isotropic data):

- the law is the per-axis Beale sum
  Σ_i Σ_k −(F_ik − m a_ik)² / (2(σ_F,ik² + m² σ_a,ik²)), i.e. each axis as its own
  one-dimensional cartesian reading (axis splitting); only the radius nuisance
  with the likelihood product has this form, so other laws on per-axis readings
  are refused (split the axes first to use them);
- the readout keeps its definition, E|F*_stack|/E|a*_stack|; its weight is
  E√(Σ_c σ_a,c² cos²θ_c v_c²), v_c ~ N(w_c, 1), a weighted noncentral chi mean
  (`weighted_stacked_mean` in `met/law.py`: the Laplace transform of the
  weighted sum is closed, and the mean is one quadrature in log t, about 1e-10;
  components sharing a noise pair are grouped). On the posterior grid it is
  computed at nodes at most 0.02 apart in log m and carried by a cubic spline.
  The same weighted mean handles SDs that differ between readings of a series,
  which ours25b refused;
- the prior's centre is not forced when the axes' ratios differ, so it is
  declared (`first_reading_geometric` or a number);
- the answer depends on the data only through the three sums P, Q, D and the
  component count of each noise pair (tested), so the exact combination state
  becomes one (C, T, M) per noise pair; the `stacked_state` carry is not built
  for that yet, and neither is the plain `ratio_of_means` (ours24's sum of sizes);
- with zero readings the median is still the prior's centre, but the readout is
  not (the axes weigh differently): 1.76 against a centre of 2 for one zero
  reading with axis ratios 5, 2, 0.8, drifting with N.

When σ_F,k/σ_a,k is the same on every axis but the SDs are not, the angle is
common and `first_reading` is allowed again. See
`presets/ours26_per_axis_noise.yaml`.

With `scores.compare_to`, a ratio-of-means readout is compared with the
reference's same-named readout, or else with its other ratio-of-means readout
(so ours25b's `ratio_of_means_rss` is compared with ours24's `ratio_of_means`). Also `median`, `geometric`
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
- the cartesian references (`tests/test_cartesian.py`), each by the checks
  its own line wrote. `cartesian2`: kernel and mean radius against radial
  quadrature; the whole law against brute force in the 21 September
  coordinates with weight r^(d−2), and in 1D against brute force in the
  record's own variables (a*, m); the record's section 3 formulas term by
  term; 2D identical to `flat`; a 3D reading equal to three 1D readings; the
  σ-free completion by brute force; its d-dependent tail slopes.
  `cartesian1`: mean radius by 2D quadrature; the whole estimator identical to
  `flat` in 2D; the single-reading law by per-component brute force;
- `cartesian1` against `cartesian2`: mean radii, radius-nuisance terms,
  single-reading laws under every split and combined radius-nuisance laws
  agree to rounding. The acceleration and force splits differ by exactly
  cos(θ)^(2−d) and sin(θ)^(2−d), and not at all in 2D. The old name is refused;
- the sequential rule's `vonmises_state` carry: exact in 2D and under `cartesian1`
  in every dimension, approximate under `flat` in 1D and 3D, peak at the
  total-least-squares direction (`tests/test_sequential.py`,
  `tests/test_cartesian.py`);
- the `hierarchical` rule (`tests/test_hierarchical.py`): the β integral's
  closed form, the law against brute force, the E[α | m] quadrature, the
  strong-excitation limit approaching `cartesian1`, and the Gaussian-excitation
  world;
- `ratio_of_means_rss` (`tests/test_stacked.py`): the noncentral chi mean in
  k dimensions against its closed forms (k ≤ 3) and a two-dimensional quadrature
  (k up to 150); the readout against an independent Gauss–Legendre quadrature in
  θ; equal to `ratio_of_means` for one reading; s at zero readings; reciprocal
  symmetry; `cartesian1` = `cartesian2`; a 3D series equal to its axis readings
  (and `ratio_of_means` not); `split_axes` against splitting by hand; no
  dependence on the recording frame; only the three sums matter; its refusals;
- ours26, per-axis noise (`tests/test_ours26.py`): the weighted mean against the
  closed forms (one group), polar-coordinate integration and nested quadrature
  (two groups); the readout for a 2D reading with two ratios against an
  independent quadrature; the law against the Beale formula; isotropic data
  equal to ours25b; a vector series equal to its axes; only the per-noise-pair
  sums matter; units and reciprocal symmetry; per-axis pooling and the known-excitation oracle; the world's
  per-axis draws and its SNR reference; its refusals, world-level and
  data-level;
- the `stacked_state` carry (`tests/test_stacked_state.py`): equal to the rule
  on all readings (law, intervals and `ratio_of_means_rss`) in 1D, 2D and 3D
  under both cartesian references; order- and grouping-free; its H²/2 equals
  the summed per-reading likelihood and its weight the stacked one; states of
  different instruments refuse to add; its refusals;
- `calibrate: sandwich` (`tests/test_sandwich.py`): the weight against the
  exact large-N factor, the tempered density, and its refusals;
- ours29 (`tests/test_ours29.py`): η's closed form against the numerical path
  and a brute-force grid on the projective line (1D, 2D, 3D); g against the Gram
  eigenvalue rule; Z (exact zero readings leave η, the median and the interval
  unchanged) and the sandwich's Bessel factor moving where η does not; η = 1 for
  one reading; rotation of each reading, the swap and units; the impossible
  pair's tilt, 1.5 nats at every σ ≤ 1 against 1.5/σ² untempered; η against the
  sandwich up to N/(N − 1) where the peak is inside m > 0; the readout weight
  unchanged; its refusals;
- the oracles (`tests/test_oracles.py`): `known_excitation` against a
  truncated normal done by quadrature; β from the world; oracles only by name;
- world designs and noise; per-series streams (a quick run is a prefix of the
  full run; grid position does not change a cell's data);
- presets refusing every missing required setting and each validation hole
  found in review; grid expansion; a reproducible end-to-end run.

## Not built yet

- `calibrate: sandwich` (ours25) estimates J from per-reading scores, so it
  depends on how a series is grouped (a 3D reading against its three axes give
  different J, and one reading gives none). Found reading the code on 25
  September. Since the 29 September ruling that one reading is one completeness
  claim, this dependence is intended; `calibrate: bartlett` (ours29) keeps it
  and removes the two things that were not (the Bessel factor, and the m > 0
  search).
- The readout's drift when zero readings are appended (the law and median do
  not move; both poolings drift to the same limit). Forced under ours25b; no
  world generates coasting readings yet.
- Reference measures other than `flat`, `cartesian1` and `cartesian2` (the tube/Jeffreys
  measure).
- Correlated errors (between axes, or between the force and acceleration
  channels), and a known direction. Per-axis (diagonal) noise is ours26; a
  correlation that is the same for every reading can be whitened away first,
  but that is not built.
- ours26's prior centre: open, like the prior's shape. The preset measures how
  far the answer moves between centres; nothing yet chooses one.
- Calibration estimated from samples, and calibration uncertainty.
- Readings from different instruments within one series (the law and the
  ours26 readout support per-reading SDs; the world does not generate them yet).
- Paired-difference tests between estimators (the values are saved, so
  these can be added without rerunning).
