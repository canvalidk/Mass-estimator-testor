# Mass estimator behavior laboratory

For the current understanding of the estimator, completed simulation evidence
and next question, start with the [investigation README](investigation/README.md).
Its [testing manual](investigation/comparison_guide.md) explains the reusable
simulation infrastructure. The material below documents the original numerical
reference and its dated behavior checks; its commands remain specific to that
reference suite.

**Current purpose, clarified by the user 2026-09-07:** investigate what follows
assuming Newton II is correct. For this conditional investigation, large error
scores describe unlikely measurement discrepancies; they do not terminate the
reconciliation. See the [steering record](https://github.com/canvalidk/VD-docs/blob/e22e150cc9cdb21d02932703e87eee5624fbf7e8/Newton-analysis/_dscn_mass_estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md).
Existing rejection/withholding tests reproduce the older working specification;
they are not acceptance requirements for the clarified purpose. The distinction
between likelihood-weighted averaging and selecting the minimum-error pair
remains explicit; the readout equation has not been replaced.

This is an executable test suite and numerical reference for the current VD
mass readout, `E[f] / E[alpha]`, dated 2026-09-06. It tests behavior; it does
not attempt another uniqueness proof or install a new production policy.

Read the [behavior walkthrough](https://github.com/canvalidk/VD-docs/blob/e22e150cc9cdb21d02932703e87eee5624fbf7e8/Newton-analysis/_dscn_mass_estimator/unsanctioned/mass_estimator_behavior_test_report_2026-09-06.md)
for the interpretation, counterexamples, and source record.

The [September 7 uncertainty walkthrough](https://github.com/canvalidk/VD-docs/blob/e22e150cc9cdb21d02932703e87eee5624fbf7e8/Newton-analysis/_dscn_mass_estimator/unsanctioned/mass_estimator_uncertainty_and_reuse_walkthrough_2026-09-07.md)
continues with the meaning of the null law, dependence when reusing the original
inputs, prediction under a newly specified acceleration, and continuous paths
through weak excitation and angular disagreement. Run
`python uncertainty_walkthrough.py` to regenerate
`results/uncertainty_walkthrough.json`. Its 17 additional tests are discovered
by `run_suite.py`; the original 40-test report remains a dated first-pass record.

## Run

For the compact formulas, convergence proofs, and a usable reporting recipe,
read [Practical limits and uncertainty](https://github.com/canvalidk/VD-docs/blob/e22e150cc9cdb21d02932703e87eee5624fbf7e8/Newton-analysis/_dscn_mass_estimator/unsanctioned/mass_estimator_practical_limits_and_uncertainty_2026-09-07.md).
`practical_formulas.py` adds a known-direction closed readout, a local aligned
uncertainty calculation with full covariance, and `summarize(estimate(...))`.
The summary supplies the existing mass point, an equal-tail probability
interval, and finite logarithmic spread; it does not apply a rejection gate.
`python practical_examples.py` regenerates its small example table in
`results/practical_examples.json`. Ten additional tests check the limits,
analytic formulas, covariance handling, and numerical refinement.

Python 3.10+ and NumPy 2.x are sufficient. No SciPy, web access, or private
repository access is needed to rerun the saved suite.

From this directory:

```powershell
python run_suite.py
python experiments.py
```

Or use the bundled runtime in this workspace:

```powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' '03_trace_and_evaluator/mass_estimator/run_suite.py'
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' '03_trace_and_evaluator/mass_estimator/experiments.py'
```

The first command runs the automated assertions and saves `results/tests.json`.
The second regenerates the named cases, vector-angle and excitation sweeps,
900,000 fit-calibration simulations, and 1,200 full estimator reruns in
`results/experiments.json`; its concise case table is `results/cases.csv`.
Results use seed 20260906. Test-specific simulations use their own fixed seeds.

## Files and supported scope

The first published-data comparison is `air_track_experiment.py`, with numerical
source transcription in `air_track_source.json`. Run it separately to regenerate
`results/air_track.json` and `results/air_track_per_run.csv`. The
[trial report](https://github.com/canvalidk/VD-docs/blob/e22e150cc9cdb21d02932703e87eee5624fbf7e8/Newton-analysis/_dscn_mass_estimator/unsanctioned/mass_estimator_air_track_trial_2026-09-06.md)
explains its scalar mechanical reduction, assumed uncertainty scenarios, and
the identical-input comparison. It is not part of the synthetic behavior suite.

- `estimator.py`: reference integration, exact branches, fixed-mass profile,
  isotropic free-mass minimum, full compatibility sets, and explicit gate.
- `test_estimator.py`: independent oracles, behavioral cases, numerical
  refinement, seeded Monte Carlo, and input failure checks.
- `experiments.py`: reproducible study rather than pass/fail tests.
- `run_suite.py`: writes a machine-readable verification record and hashes.
- `results/`: saved outputs, including failed tests if a future run regresses.

The integrator supports 1D, 2D and 3D vectors, uniform unknown common direction
or an explicitly supplied known direction, and a full positive definite Gaussian
covariance, including anisotropy and force/acceleration cross-correlation. The
latent measure is **flat in the two positive magnitudes**, as in the working
specification. It is not uniform Cartesian volume, flat mass, or flat log mass.

`profile` supports full covariance and an unknown latent vector direction.
`isotropic_qmin` and `compatibility_intervals` are specialized to independent
isotropic channels and unknown direction. Do not use them as the gates for a
known-direction or general-covariance experiment. `gate` only applies a supplied
threshold; it does not choose or calibrate that threshold.

Exact recovery compares the supplied floating representations without a tolerance
repair. It is a branch fixture for exact representable examples, not a symbolic
algebra engine. An empirically rounded zero is not an exact input.

## Numerical method and error handling

For numerical scaling let `sf` and `sa` be channel RMS standard uncertainties.
Use `f=sf*r*sin(theta)` and `alpha=sa*r*cos(theta)`, with positive radius and
`theta` between zero and pi/2. The Jacobian is proportional to `r`. For each
common direction, the Gaussian exponent is quadratic in `r`; radial zeroth and
first-magnitude integrals reduce to `J1` and `J2`. This reduces the original
integral without replacing the measure or truncating the mass domain.

The radial integrals use stable scaled formulas and a negative-tail quadrature.
The outer integrals use a periodic direction grid in 2D, Gauss-Legendre/polar
quadrature in 3D, and composite Simpson integration over the ratio angle.
Quantiles use an independently normalized trapezoidal CDF on that angle grid.
No mass cutoff is used to calculate the point or quantiles. The displayed
log-mass curves have a finite viewing window; the saved experiment records its
approximate displayed probability separately.

Saved cases compare 128 directions / 1024 ratio intervals with 256 / 2048.
This is an observed refinement error, not a certified bound for arbitrary inputs.
Highly concentrated signals, extreme covariance condition numbers, and extreme
tails need additional refinement. The caller chooses quadrature resolution;
the library does not silently promise convergence. The test suite independently
checks a full-vector integral in the original magnitude/direction coordinates.

The full empirical pipeline's object/time/frame and provenance checks, hard
bounds, non-Gaussian policies, common-mass posterior integration, and production
admission classifications are outside this original reference implementation.
Common-mass profile compatibility has a test here. The later investigation
develops [joint inference and confidence constructions](investigation/joint_inference_and_uncertainty.md)
under declared models; those results do not install a new production policy.

## Robot and cart real-data trials (2026-09-06)

Run `python robot_cart_experiments.py`, then `python robot_vector_check.py`.
The former reproduces the robot calibration, evaluates scalar projections and
actual weak components, and diagnoses the cart exports. The latter evaluates
the full unknown-direction 3D integral at the first weak-z robot sample and
compares ordinary vector least squares and a same-covariance profile fit.

Raw inputs and metadata are preserved under `raw/`; checksums and results are
in `results/robot_cart.json`, `results/robot_cart_examples.csv`, and
`results/robot_full_vector.json`. The report is
[mass_estimator_robot_cart_trials_2026-09-06.md](https://github.com/canvalidk/VD-docs/blob/e22e150cc9cdb21d02932703e87eee5624fbf7e8/Newton-analysis/_dscn_mass_estimator/unsanctioned/mass_estimator_robot_cart_trials_2026-09-06.md).
Uncertainty budgets are provisional. A finite weak-signal readout does not
establish mass identification, and these trials do not show improved accuracy.

## Source identity

Managed inventory: `canvalidk/VD-docs/catalog.yaml` at
`743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6`, resolved from `main` on 2026-09-06.
Construction: `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`.
The report records the additional managed and local sources. No managed sources
were copied into this production directory or changed.
