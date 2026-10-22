# First published-data trial: air-track mass estimation

Date: 2026-09-06. Status: local active analysis; awaiting transfer to VD-docs,
not submitted. This is an exploratory estimator comparison, not a calibrated
experimental validation or a claim of improved accuracy.

## Result

The full compatible-pair integral, specialized to the observed one-dimensional
motion, gives **1,516.952222 g**. The ordinary force/acceleration ratio on the
same processed inputs gives **1,516.952222 g**. No meaningful point-estimate
difference was resolved. Across the six uncertainty scenarios, their numerical
difference is below 0.000001 g; the actual observed maximum is recorded in JSON.

| Quantity | Mass |
|---|---:|
| Independently measured system mass reported in the paper | 1,502 g |
| Paper's reported mass from its acceleration-versus-load slope | 1,518 g |
| Our reconstruction from its published timings, ordinary ratio | 1,516.95 g |
| Our integral on exactly the same effective force and acceleration | 1,516.95 g |

The small difference from the published fitted mass is not a benefit of the
larger equation. Their displayed slope has only three decimal places; applying
their mass conversion to that rounded slope gives 1,517.885 g. Our reconstruction
retains more digits from the timing table. The integral does not move the result
closer to the independently measured mass.

## Experimental source and extraction

Hessel, R., Canola, S. R., and Vollet, D. R. (2013), *An experimental verification
of Newton's second law*, Revista Brasileira de Ensino de Fisica 35, 2504.
[DOI](https://doi.org/10.1590/S1806-11172013000200024),
[publisher PDF](https://sbfisica.org.br/rbef/pdf/352504.pdf).

Experiment 1 moves loads from an air-track glider to a hanging load while keeping
the combined mass fixed. Table 1 supplies five travel intervals per run. We use
its 12.00 cm post spacing, 1 kHz timer, and reported local g = 9.76 m/s^2.
The complete numerical transcription, source URLs, and extraction limitations
are in `03_trace_and_evaluator/mass_estimator/air_track_source.json`.

The PDF text supplied the timing table. Reconstructing and rounding the first
run's speeds reproduces Table 2. Reconstructing the force-line intercept also
closely matches the paper's displayed intercept. Web screenshots of the tables
were unavailable, so no visual source-table check is claimed.

The paper does not supply the full uncertainty budget needed by our equation.
The uncertainties below are our explicitly assumed scenarios.

## Mechanical boundary: a scalar reduction

The hanging load accelerates vertically while the glider accelerates
horizontally. Their combined Cartesian acceleration is not the glider's
acceleration. Simply calling the hanging weight a net force on a single rigid
object would be incorrect.

We adopt the standard scalar coordinate q = string travel. With an inextensible,
negligible-mass string, negligible pulley inertia, and approximately constant
resistance R, the reduced equation is

`h*g - R = M*q_ddot`.

Here h is hanging mass and M is the effective inertia in this coordinate. Under
the stated idealizations M equals glider-plus-hanging-load mass. Pulley inertia
would contribute effective inertia; varying drag and track tilt could violate
the constant-R approximation. This is a test of the scalar readout following
that mechanical reduction, not a test of the full Cartesian vector pipeline.

## Identical-input comparison

For each row, convert travel intervals to seconds, compute midpoint times and
average speeds, then fit speed against midpoint time with a free intercept:

`t_i = sum(intervals up to i) - interval_i/2`, `v_i = d/interval_i`.

The resulting accelerations in m/s^2 are:

`[0.30996668, 0.62531238, 0.94967099, 1.25927601, 1.60147324]`.

We use unrounded speeds. This is the paper's stated kinematic procedure; we do
not use the known system mass to infer acceleration or force.

To remove an unknown constant resistance without first estimating it from mass,
apply the same linear contrast to force and acceleration:

`w = [-2, -1, 0, 1, 2]/10`, `F_eff = sum(w_i*h_i*g)`,
`a_eff = sum(w_i*a_i)`.

Since sum(w)=0, the constant resistance cancels. These weights are proportional
to centered hanging loads and reproduce the unweighted regression slope per
50 g load increment. This gives

`F_eff = 0.488 N`, `a_eff = 0.3216976730 m/s^2`.

Then compare ordinary `F_eff/a_eff` against the integral's `E[f]/E[alpha]`
using that very same pair. The pair is a chosen summary of five runs. Applying
the flat latent measure to this summary is an analysis choice; it is not claimed
to equal a joint latent model of all five runs with one shared mass.

The inferred intercept corresponds to about 0.0242 N of constant resistance.
This is a diagnostic from the fit, not an independently measured force input.
Directly dividing each uncorrected hanging weight by its acceleration gives
apparent masses from about 1,524 to 1,574 g. Feeding those uncorrected inputs to
our integral leaves those values essentially unchanged. It cannot repair an
omitted force by itself. Per-run results are saved separately for that reason.

## Uncertainty scenarios

The reference code supports one spatial dimension with direction u = +1 or -1,
positive latent magnitudes, and the existing flat-magnitude measure. We use
independent Gaussian force/acceleration errors at the effective-pair level.

Two acceleration-error scenarios were examined:

1. **Timing and spacing model:** 100,000 seeded perturbations, with six
   independent uniform crossing-time errors of +/- 0.5 ms for each run, and
   one shared uniform spacing error of +/- 0.05 mm. The timing interval errors
   are correlated because adjacent intervals share a crossing. Each perturbation
   recomputes times, speeds, slopes and the contrast together. This gives an
   assumed acceleration standard uncertainty of 0.0023335 m/s^2. Clock calibration,
   individual post placement, and all other systematic errors are excluded.
   Display precision motivates this scenario; it does not establish that this
   was the experiment's actual error distribution.
2. **Between-run scatter model:** use the standard unweighted regression slope
   error, conditional on exact nominal loads and independent equal-variance
   acceleration residuals. This gives 0.0028060 m/s^2. Residuals could also contain
   model error, so this is an alternative sensitivity scenario, not a measured
   noise calibration and not an additional independent error to add in quadrature.

For each, assume force standard uncertainty of 0.1%, 1%, or 5%. These correspond
to a common multiplicative driving-force scale uncertainty, independent of the
timing errors. They are not independent errors on every load. Actual individual
load uncertainties and changes in resistance are not characterized.

Every scenario gives the same point estimate to numerical precision. With the
timing model and 1% force uncertainty, the induced latent mass-ratio distribution
has an equal-tail 95% interval of approximately **1.481 to 1.554 kg**. For 0.1%
force uncertainty it is about 1.495 to 1.539 kg; for 5% it is about 1.367 to
1.667 kg. These are conditional distribution intervals, not verified frequentist
confidence intervals or published experimental uncertainties. They must not be
used to claim a statistically established improvement over the source analysis.

## Why there is no visible point shift

In one dimension with independent errors, the positive-direction branch nearly
factorizes into two positive-truncated Gaussians. Here both measurements are
many standard uncertainties from zero. Removing their negative tails therefore
changes their means by a negligible amount. Thus `E[f]/E[alpha]` is effectively
the ordinary ratio even when the induced ratio distribution has visible width.

An independent analytic check includes both direction branches. For
`z_F=F/sigma_F`, `z_a=a/sigma_a`, and
`D=Phi(z_F)+Phi(z_a)-1`, the ratio is

`[F*D + sigma_F*phi(z_F)] / [a*D + sigma_a*phi(z_a)]`.

This is an evaluation of the current one-dimensional integral, not a replacement
estimator or an assumption that E[f/alpha] is finite.

As a separately labeled hypothetical exercise, enlarging acceleration uncertainty
to 50% while holding force uncertainty at 1% changes the point to approximately
1.476 kg. At 100% acceleration uncertainty it becomes approximately 1.178 kg.
These are not the measured error bars and do not demonstrate better accuracy.
They show where a point difference can become visible: near the positive-domain
boundary, not automatically whenever uncertainty is nonzero.

## Reproduction and checks

Run `air_track_experiment.py` with the existing Python/NumPy environment.
It writes `results/air_track.json` and `results/air_track_per_run.csv`.
It checks Table 2 reproduction, compares the full integral with the independent
analytic expression, and doubles the ratio quadrature starting at 8,192 intervals
until the endpoint quantile change is below 0.000005 kg. Point refinement and
quantile refinement are checked and recorded. The first 8,192/16,384 quantile
check failed; additional refinement resolved the narrowest distributions.
Source, analysis-script, and estimator SHA256 hashes are recorded in the result.
The existing estimator implementation was not changed.

Managed estimator provenance is inherited from the behavior laboratory:
`canvalidk/VD-docs` commit `743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6`,
`Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`.
Active local scope: estimator.py and README.md, plus this trial's new source,
script and outputs. No new managed-context synthesis or claim of exhaustive
document coverage was needed. No managed files were changed.

## What this trial establishes

The equation can be run on a scalar summary of published real measurements.
It reproduces the usual high-signal ratio on identical inputs. It does not show
an accuracy advantage here, nor does it test zeros, vector misalignment, the
adequacy of a full force inventory, or an independently calibrated uncertainty
law. A stronger follow-up needs measurements close enough to the uncertainty
boundary to distinguish estimators, with independently characterized force and
motion uncertainties and a withheld mass benchmark.
