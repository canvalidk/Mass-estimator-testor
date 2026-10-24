# Refined mass-estimator comparison

This test compares **11 isotropic or 16 diagonal-noise point estimators on the same noisy observations**, using
several explicitly different meanings of accuracy. It also scores the two
available probability laws. The original full-covariance estimator is unchanged.
The benchmark supports independent Gaussian errors with equal or unequal
coordinate SDs, with the common latent 3D direction unknown to the estimators.
Newton II holds in every simulated trial; observations are never rejected for
failing an alignment or compatibility test.

## What is scored

Let the independently specified true mass be \(m\), the estimate be
\(\widehat m\), \(r=\widehat m/m\), and \(e=\log r\). All losses below are
smaller-is-better. An estimated probability law is always evaluated against
simulation truth, not against a value drawn from that same estimated law.

| Criterion | Formula or output | What it measures |
|---|---|---|
| Absolute log error | \(|e|\) | Typical multiplicative error; treats a factor two and one-half equally. |
| Squared log error | \(e^2\) | Multiplicative error with more weight on large log errors. |
| Reciprocal square-root loss | \(\sqrt r+1/\sqrt r-2\) | Reciprocal symmetry, with stronger penalties far out in the tails; approximately \(e^2/4\) near truth. |
| Absolute relative error | \(|r-1|\) | Conventional percentage error after multiplication by 100; over- and underestimation have different multiplicative penalties. |
| Squared relative error | \((r-1)^2\) | Squared mass error divided by true mass squared. |
| Absolute and squared mass error | \(|\widehat m-m|\), \((\widehat m-m)^2\) | Error in mass units, useful if absolute physical error is the objective. |
| Capped squared log error | \(\min(e^2,\log^2 K)\), \(K=2,5,10,20\) | Bounded alternatives whose failure penalty is stated in advance. |
| Tolerance success | \(P(1/K\le r\le K)\) | Reports eight thresholds and a 401-point curve from factor 1 to 20, retaining the old factor-two score. |
| Error diagnostics | Mean signed log error, log variance; 50th, 90th, 95th, 99th absolute log-error quantiles | Separates centring, variability and exceptional failures. |
| Joint magnitude loss | \(f/\sqrt{\widehat m}+\alpha\sqrt{\widehat m}-2\sqrt{f\alpha}\) | The original readout's loss; equals \(\sqrt{f\alpha}\) times reciprocal square-root loss when \(m=f/\alpha\). |

Log-ratio error and squared log error are supported by
[Tofallis (2015), §3](https://arxiv.org/pdf/2105.05249); absolute log error
applies ordinary absolute loss on that scale. Ordinary absolute,
squared and percentage losses, and the dependence of an optimal point on its
loss, are discussed by [Gneiting, §§2–3](https://arxiv.org/html/0912.0902v2).
The robotic mass-error comparison in
[Nadeau, Giamou and Kelly, §V-B Eq.10](https://arxiv.org/html/2203.00830v3#S5.SS2)
uses absolute percentage mass error. The reciprocal square-root and joint
losses come from our existing managed derivation and its normalization;
they are not attributed to these papers.

These are complementary objectives. They are not interchangeable votes for a
single overall winner. Changing a score after seeing which one favors an
estimator is not evidence that the estimator is generally better.

## Added equations and retained comparators

The vectors below are standardized by their component noise SDs. Define
\(p=\|\mathbf F\|^2\), \(q=\|\mathbf a\|^2\), and
\(d=\mathbf F\cdot\mathbf a\). Multiply a returned point by
\(\sigma_F/\sigma_a\) to obtain physical mass units.
Expectations in the first four rows use our existing flat-magnitude joint law
\(P\), with the same Gaussian likelihood and unknown common direction.

| Code name | Equation | Short description |
|---|---|---|
| `flat_joint` | \(E_P[f]/E_P[\alpha]\) | Our existing equation; ratio of two magnitude means. |
| `flat_geometric` | \(\exp(E_P[\log M])\) | Added: point minimizing expected squared log loss under the same law. |
| `flat_median` | \(\operatorname{median}_P(M)\) | Added: point minimizing expected absolute log loss under the same law. |
| `flat_reciprocal_root` | \(E_P[\sqrt M]/E_P[1/\sqrt M]\) | Added: point minimizing expected normalized reciprocal square-root loss. There is **no outer square**. |
| `tube_joint` | \(E_T[f]/E_T[\alpha]\) | Existing sensitivity comparator: changes reference volume to \((f^2+\alpha^2)df\,d\alpha\,d\Omega\) in standardized units. |
| `norm_ratio` | \(\sqrt{p/q}\) | Existing vector-length ratio; ignores magnitude noise bias. |
| `norm_floor` | \(\max(\sqrt p,c)/\max(\sqrt q,c)\), \(c=\sqrt{2/\pi}\) | Existing heuristic with an explicit noise floor. The floor is not the expected norm of 3D noise. |
| `positive_profile` | \(\arg\min_{m>0}(p-2dm+qm^2)/(1+m^2)\) | Existing Gaussian profile / total least-squares fit; uses errors in both channels, choosing the best latent fit rather than averaging it. |
| `forward_ols` | \(d/q\) | Added: forward least squares, treating measured acceleration as exact. |
| `reverse_ols` | \(p/d\) | Added: fit acceleration against force, treating measured force as exact, then invert the slope. |
| `noise_corrected_ols` | \(d/(q-3)\), when \(q>3\) | Added: subtract expected 3D acceleration-noise energy before taking a moment ratio. Nonpositive corrected energy is unavailable. |

The OLS and corrected-moment approaches are grounded in
[Kelly (2007), §3 and §7.1](https://arxiv.org/pdf/0705.2774).
The correction follows \(E[d]=m\|\mathbf a^*\|^2\) and
\(E[q-3]=\|\mathbf a^*\|^2\). An unbiased estimating equation does not make
its one-pair ratio unbiased. This is particularly unstable near \(q=3\).
The role of TLS when both channels contain measurement error also appears in
[Kubus et al. (2008), §IV-D Eq.23](https://www.researchgate.net/publication/224339569_On-Line_Estimation_of_Inertial_Parameters_Using_a_Recursive_Total_Least-Squares_Approach).
Our formulas are simple through-origin, single-pair specializations; they do
not reproduce Kelly's latent-population regression or a recursive robot solver.
The two OLS exact-channel assumptions are deliberately violated in this test,
so they serve as informative baselines, not fully matched competitors.

The three added summaries preserve our probability law while changing its
decision rule. Posterior optimality for a chosen loss does not guarantee
smaller repeated-experiment error at every fixed true mass. A posterior
arithmetic mean of \(M\) is deliberately absent: it diverges under the flat
law even though \(E[f]/E[\alpha]\) and the displayed log/half moments exist.

## Probability-law and prediction scores

For the flat and tube laws, the test reports:

- Central 50%, 80% and 95% interval **coverage**, mean width on the log-mass
  scale, and interval score. Coverage is checked separately at every true
  signal pair. Narrowness alone is not a success criterion.
- **CRPS of log mass**: \(E|Z-y|-\tfrac12E|Z-Z'|\), where
  \(Z=\log M\) and \(y=\log m\). This evaluates the whole distribution.
- **Negative log density of log mass**, \(-\log p_Z(y)\), with the density's
  reference coordinate fixed in advance.

The proper scoring rules and interval score come from
[Gneiting and Raftery (2007), §§4.2 and 6.2](https://sites.stat.washington.edu/people/raftery/Research/PDF/Gneiting2007jasa.pdf).
Applying them to log mass is our explicit choice for this positive, heavy-tailed
problem. The interval score at level \(1-\gamma\) is
\((u-l)+2(l-y)1_{y<l}/\gamma+2(y-u)1_{y>u}/\gamma\).
All probability-law scores use the smaller-is-better convention.
The other point equations are not assigned invented uncertainty intervals.
The three added flat-law points all share the same flat-law uncertainty score.

Two secondary tests predict a new measured force at twice the training
experiment's true excitation, using either a supplied exact acceleration or
an independently noisy acceleration. Both use new Gaussian force noise and
normalize squared vector residual by squared true force magnitude. This is a
specified transfer experiment, inspired by the held-out validation in
[Wensing et al. (2017), §VI](https://arxiv.org/pdf/1701.04395v3), not a robot
validation result. With noisy acceleration, prediction loss also penalizes
amplified input noise and can reward shrinkage; it must not be relabelled as
pure mass accuracy. The exact-input version closely tracks squared mass error.

## Fairness, failures and numerical checks

The original design estimates one mass from one vector pair. All methods receive the
same observations and the same declared noise information. The original 16
signal pairs are retained, using the original training seed. A further 30
cells hold true mass fixed at 0.25, 1, 4, 16 or 32 while total true signal
\(\sqrt{f^2+\alpha^2}\) varies over 0.5, 1, 2, 4, 8 and 16.
This separates changing excitation from the earlier changing-mass path.
No global average over these arbitrarily selected cells is presented.

Nonpositive, NaN or infinite point outputs remain failures: they fail every
finite tolerance, receive the maximum explicitly capped loss, and receive
infinite loss for unbounded criteria. This is an explicit extension of the
positive-mass losses to invalid outputs. A null summary means an unbounded or
undefined result, never zero error. No valid-only average is substituted.
The Gaussian profile and regression baselines have positive-probability invalid
outputs even where a finite simulation happens to observe none. A finite
empirical mean is not a proof of finite population risk.

The norm ratio has an additional subtlety in 3D: its squared mass error has a
finite mean but infinite variance. Its empirical squared-loss Monte Carlo
standard error is therefore not justified by the usual finite-variance
argument. The same caveat applies to its squared relative error, squared
held-out residuals and the corresponding paired differences. Reverse and
corrected regression also have denominator poles.
Bounded tolerance and capped-loss comparisons avoid these moment issues.

Every score difference is computed on paired trials. Negative
`flat_joint minus comparator` differences favor our equation. The test gives
empirical Monte Carlo standard errors, not a correction for searching across
many scenarios and scores. Errors are shared across scenarios too; cells
must not be pooled as independent experiments.

The probability calculation uses the exact angular reduction of the existing
3D model. The numerical implementation checks analytic zero-data identities,
channel reciprocity, agreement with the earlier angular calculation, and
refinement of points, quantiles and distribution scores. Each benchmark cell
also repeats its first 32 observations at twice the integration order.
Numerical approximation error and Monte Carlo error are reported separately.

## Running and extending the test

### Unequal coordinate noise

Select **Independent Gaussian noise · different coordinate SDs** in the
experiment form to supply three force SDs, three acceleration SDs, and the
true signal direction. The direction is used only to generate trials. Every
estimator still has to infer mass without knowing the true direction. Errors
remain independent across coordinates and channels; this is diagonal
covariance, not a correlated-noise experiment.

Noise profiles in the form preserve each channel's RMS SD. The scenario's
force and acceleration SNRs mean latent magnitude divided by the respective
channel RMS SD. Estimates and truth are reported in physical force/acceleration
units; fixed-mass designs hold physical mass fixed even when the channel RMS
SDs differ. Equal RMS noise does **not** mean equal directional information:
rotating a signal relative to precise/noisy axes is an explicit experimental
variable. The initial validation covers within-channel SD ratios up to 16.
Diagonal runs currently allow 8,192 trials per cell and 262,144 total trials.

The original flat law remains `df d(alpha) dOmega` in physical magnitudes.
Only scalar channel standardization is used. Independently whitening the two
vectors coordinate by coordinate and calling the old isotropic equation would
change their shared-direction relation and is not done. The tube sensitivity
law uses the explicit extra weight
`(f/RMS_force_SD)^2 + (alpha/RMS_acceleration_SD)^2`.

The diagonal calculation integrates the common direction through Gaussian
radius moments. In RMS-scaled coordinates, write
`(f, alpha) = r (sin(theta), cos(theta))`, `x = r u`. The flat measure becomes
`dtheta d^3x / ||x||`. At each angle the likelihood in `x` is Gaussian with
diagonal precision. Thus its ordinary Gaussian integral times `E[1/||x||]`
gives the flat probability kernel; multiplication by `E[||x||]` gives the tube
kernel. Flat first-magnitude integrals use the Gaussian integral itself;
tube first moments additionally use `E[||x||^2]`. A positive one-dimensional
Laplace integral evaluates the radius moments. This removes the expensive
direction grid while preserving the existing model. The log-mass integration
and scoring machinery are shared with the isotropic engine.

Alongside the original methods the diagonal run adds:

| Code name | Equation or role |
|---|---|
| `covariance_profile` | Minimize `sum_i (F_i - m a_i)^2 / (sigma_Fi^2 + m^2 sigma_ai^2)` over positive mass, comparing all numerical stationary points and both boundaries. Uses all supplied diagonal noise information. |
| `weighted_forward_ols` | `sum(F_i a_i / sigma_Fi^2) / sum(a_i^2 / sigma_Fi^2)`; force-noise weights but still treats acceleration as exact. |
| `weighted_reverse_ols` | `sum(F_i^2 / sigma_ai^2) / sum(F_i a_i / sigma_ai^2)`; acceleration-noise weights but still treats force as exact. |
| `weighted_corrected_ols` | Weighted forward ratio with `sum(sigma_ai^2 / sigma_Fi^2)` subtracted from its denominator; retains unavailable/nonpositive outputs. |
| `flat_joint_rms` | Our original isotropic equation using each channel's RMS SD. A control showing the cost or benefit of discarding coordinate-specific noise information. |

The covariance profile is the through-origin, independent-error specialization
of the weighted errors-in-variables criterion, described in
[Mikkonen et al. (2019), Appendix A4](https://acp.copernicus.org/articles/19/12531/2019/).
Our solver enumerates the degree-at-most-ten stationary polynomial in mass and
reciprocal mass, polishes roots against the rational derivative, and compares
their objectives with both boundaries. It does not assume one local minimum.
This is numerical root finding, not certified symbolic root isolation.

The legacy profile and norm floor use RMS channel SDs; their labels identify
that approximation. Ordinary noise-corrected OLS subtracts the full known
acceleration-noise trace. All positive-boundary and invalid-output policies
remain unchanged. In particular the covariance profile also has infinite
unbounded-loss risk under this failure policy, so bounded tolerance or capped
loss comparisons are necessary to compare its finite performance fairly.

Validation includes the independent `estimator.py` full-covariance reference,
analytic Gaussian-radius and zero-data laws, equal-noise reproduction,
coordinate permutation and sign changes, channel swapping and physical-unit
changes, and truth-independent point estimates. Every generated cell repeats
its first 32 observations at twice the quadrature order. This checks numerical
precision separately from Monte Carlo sampling error.

### Noise calibration and repeated readings

The experiment form separates the **true measurement noise**, used to generate
readings, from the **supplied noise information**, used to estimate mass.
An optional `measurement` configuration adds repeated readings and calibration:

```json
{
  "measurement": {
    "repeats": 10,
    "calibration": {
      "mode": "estimated",
      "samples": 30,
      "seed": 2026091923
    }
  }
}
```

Calibration has three modes. `known` supplies the true coordinate SDs and is
an oracle control. `scaled` multiplies the true SDs by declared three-coordinate
`force_scale` and `acceleration_scale` factors; for example, a factor of 0.75
understates an SD by 25%, rather than understating its variance by 25%.
`estimated` supplies an independently estimated SD for each coordinate in
each Monte Carlo trial. With `k` calibration readings, it draws
`s_hat = s_true * sqrt(chi_squared(k-1)/(k-1))`, exactly the sampling law of
the ordinary sample SD from independent Gaussian calibration measurements
after subtracting their sample mean. The calibration seed differs from both
measurement and held-out seeds. This simulates fresh calibration in each
trial, so Monte Carlo uncertainty includes calibration variability. It does
not describe many measurements sharing one fixed calibration error.

Every noise-aware estimator receives the same supplied SDs. Noise-blind
methods continue to ignore them. The posterior treats supplied variances as
fixed: this is **plug-in calibration**, not inference that integrates over
uncertainty in the calibrated variances. Interval coverage measures the
consequence of that assumption. Finite calibration can also create unequal
estimated coordinate SDs when the true sensor is isotropic.

Deterministically specified coordinate SD ratios remain limited to 16.
Estimated calibration supports realized ratios up to 32, covered by additional
numerical stress checks. A draw beyond the supported range stops the whole
run visibly; no trial is clipped, resampled, or silently omitted. In these
runs numerical refinement checks the first 32 trials plus extreme supplied
noise contrasts, standardized reading strengths, and channel RMS ratios;
the checked indices are saved. This is sampled convergence evidence, not
a proof of numerical accuracy for every possible observation.

Repeated readings keep the latent force vector and acceleration vector fixed.
Their errors are independent between repeats. For `n` readings of a vector
`x`, the Gaussian likelihood satisfies
`sum_j ||x_j - mu||^2_Cinv = sum_j ||x_j - x_bar||^2_Cinv + n ||x_bar - mu||^2_Cinv`.
The first term is independent of the latent vector. Therefore the sample mean
and covariance `C/n` give exactly the same posterior for the latent pair as
the full repeated-reading likelihood, conditional on the supplied covariance.
The engine samples this sufficient mean directly, with true SDs divided by
`sqrt(n)`, and supplies calibrated SDs divided by the same factor. It combines
observations before estimating mass; it does not average estimated masses.

The force and acceleration SNRs in the scenario remain **per-reading** SNRs.
Effective sample-mean SNRs are `sqrt(n)` times larger and are recorded
separately. Physical truth and held-out measurement noise are unchanged.
Shared offsets, drift, correlations between repeats, and changing excitation
are outside this repeated-reading model. Combining different latent force
and acceleration pairs with a common mass requires a different joint model.

The unchanged primary comparison outputs include all point losses, failures,
paired Monte Carlo errors, and posterior coverage/width/proper scores.
Saved trial arrays include supplied noise information so that fairness and
the separation of calibration from simulation truth can be checked.

### Local results browser

Start the browser server from the workspace root, then open
<http://127.0.0.1:8766>:

```powershell
python -B 03_trace_and_evaluator/mass_estimator/investigation/results_app.py
```

This is an offline local application with no frontend build or web-service
dependency. It requires the same Python and NumPy installation as the runner.
Stop the server with Ctrl+C; leave it running while an experiment is executing.
The library reads the saved full comparison, its validation run, and the
September 18 boundary studies from their existing `.tools/` locations. Original
files are not moved or rewritten. A fresh checkout without those local result
files starts with an empty library and can generate new runs through the form.

Choose a run, a path through the signal conditions, a score, and a comparator.
The plots distinguish individual mean losses from **paired differences on the
same trials**. Negative differences favor `flat_joint` for loss scores. Error
bars, where justified, are pointwise approximate 95% Monte Carlo intervals
(mean ± 1.96 MCSE); they are not mass uncertainty intervals, simultaneous bands,
or exact crossover boundaries. Lines interpolate sampled cells and do not
establish behavior between them. Cells reuse noise and must not be pooled as
independent trials. Coverage is compared with its nominal target; interval
width alone is not a performance ranking.

Historical boundary runs saved fewer quantities than the full comparison.
The browser shows their available paired scores without inventing missing
individual scores or standard errors. Antithetic studies retain their paired
Monte Carlo error calculation. Failed estimates, infinite population risk, and
the norm ratio's descriptive-only squared-loss MCSE remain explicit.

The experiment form saves one JSON configuration per run, including the
scenario, trial count, training and validation seeds, integration order,
and trial-array option. Presets, fixed-force paths, fixed-mass paths, and
explicit signal pairs all call the existing comparison engine. All supported
estimators and implemented scores are calculated together; choosing a
displayed score does not refit the estimates. New outputs go to
`.tools/mass_estimator_runs/<run-id>/`, with `config.json`, `status.json`, and
`results.json`. The results contain source hashes and numerical-refinement
checks. Failed runs retain their configuration and error instead of appearing
as successful comparisons.

Export a configuration to repeat an experiment using either the form or CLI:

```powershell
python -B 03_trace_and_evaluator/mass_estimator/investigation/compare_estimators.py --config path/to/config.json --output .tools/mass_estimator_runs/my-repeat
```

Reproducibility requires both the saved configuration and the recorded code
version; configuration alone does not freeze future code changes. The local
server accepts one calculation at a time and binds to the loopback interface.
There is one Python server process and, only while computing, one fresh Python
worker process. The worker loads the current engine for each experiment;
configuration and source checks prevent publishing a result if either changes
during execution. Closing the browser tab does not stop a run. Stopping the
server normally records interruption and stops its worker.

`results_app.py` serves the interface and manages jobs; `results_catalog.py`
adapts saved outputs for display; `experiment_config.py` validates configurations
and feeds the same simulation runner. The frontend is three static files in
`results_ui/`. New parameter sweeps should be configurations, not new Python
scripts. Diagonal configurations add `noise_model: "diagonal_gaussian_3d"`,
`noise: {"force_sd": [sx, sy, sz], "acceleration_sd": [sx, sy, sz]}` and
`direction: [ux, uy, uz]`. The direction is normalized when saved. Legacy
isotropic configurations retain their original schema and numerical behavior.

### Command-line comparison

Run from the workspace root with Python and NumPy available:

```powershell
python -B -m unittest discover -s 03_trace_and_evaluator/mass_estimator/investigation -p 'test_comparison_*.py'
python -B 03_trace_and_evaluator/mass_estimator/investigation/compare_estimators.py --samples 32768 --scenario-set all --output .tools/mass_comparison_20260917/full
python -B 03_trace_and_evaluator/mass_estimator/investigation/comparison_report.py .tools/mass_comparison_20260917/full/results.json
```

`comparison_scores.py` owns the point losses and interval score;
`comparison_estimators.py` owns the six direct comparators;
`comparison_posterior.py` owns the probability laws and their summaries;
`compare_estimators.py` owns simulation, paired evaluation and output;
`comparison_report.py` makes readable comparison tables without refitting.
Results, optional per-trial arrays (`--save-trials`), and generated summaries
belong in the specified scratch directory. The trial arrays let us add another score
without changing or refitting the estimator outputs. Training and validation
seeds, source hashes and numerical-refinement checks accompany `results.json`.
The completed 32,768-trial, 46-cell run passed every refinement check. All
39 focused tests passed; its retained 16 original cells reproduced the earlier
factor-two success fractions exactly, and squared-log means to within
\(1.62\times10^{-11}\). This run saved aggregate results; per-trial arrays
can be regenerated with the same seeds and the optional flag.

The comparison does not yet establish behavior for correlated noise,
calibration uncertainty integrated into inference, repeated pairs with varying
excitation, shared sensor offsets, non-Gaussian noise, or real robot data.
Those are changes to the experiment, not alternative scores for the present one.

## Context provenance

Managed derivation and steering remain pinned to
[VD-docs 8782a293297330f7a354d87ebe62792018fa8866](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866),
the coherent snapshot used in this conversation. Relevant managed sources are
`Newton-analysis/_dscn_mass_estimator/mass_estimator_additional_derivations_2026-09-05.md`
and `mass_estimator_law_assumed_correct_steering_2026-09-07.md` in the same folder.
Active local sources are `crossover_vector_checks.py`, `crossover_scalar_checks.py`,
`crossover_regime.md`, `operational_warrant.md`, and the original `estimator.py`.
This is a scoped implementation and source recheck, not a new exhaustive
managed-collection audit. No managed source has been edited.
