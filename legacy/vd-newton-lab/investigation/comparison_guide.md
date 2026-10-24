# Mass-estimator testing manual

Read the [investigation README](README.md) first for the current argument,
completed-study findings and next question. This guide owns the estimator
equations, experimental assumptions, scores and reproduction instructions.

The tester compares point estimates and probability laws on matched simulated
observations, using explicitly different meanings of accuracy. Available
methods depend on the experiment: the base isotropic mode has 11 point rules
and two laws; diagonal mode has 16 point rules; pooled-calibration mode has
22 point rules and four laws, including a true-noise oracle. The original
full-covariance reference estimator remains separate.
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

For each implemented probability law, the test reports the following. The
base laws are flat and tube; pooled-calibration experiments additionally
include calibration-integrated flat and known-noise oracle flat laws.

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

The reported intervals describe a **conditional mass law**: the measurement
likelihood normalized against the declared physical reference measure
`df d(alpha) dOmega` for the flat law. A central 95% interval contains 95% of
that conditional law. Its coverage at a fixed true mass is a separate
repeated-experiment property, measured by the simulations. Neither the
reference measure nor propagation of noise-calibration uncertainty promises
95% coverage at every mass and excitation. Coverage must be read alongside
both miss directions, interval width and proper scores.

Two different zero cases explain this distinction. At **exactly observed**
zero force and acceleration with known isotropic SDs, the flat-law mass divided
by `sigma_F/sigma_a` has a half-Cauchy distribution. This is a conditional
consequence of the reference measure and measurement scales. At **true zero
excitation**, the observations are random, but their sampling distribution
does not depend on the declared positive mass. For any procedure returning
`0 < lower < upper < infinity` almost surely, coverage then tends to zero as
the true mass tends to zero or infinity: the coverage indicator tends to zero
and is bounded by one. This is a direct argument for our null model, related
to the need for unbounded confidence sets in ratio problems discussed by
[von Luxburg and Franz (2009)](https://arxiv.org/abs/0711.0198).
Wider intervals from uncertain calibration cannot recover information absent
from that experiment.

Computational validation has a different target again.
[Talts et al. (2020), §§2–4](https://arxiv.org/html/1804.06788v2) formulate
simulation-based calibration (SBC) by drawing parameters from a proper prior,
simulating observations and comparing posterior ranks with their known
distribution. The current flat magnitude reference measure is improper,
although its conditional law is proper; it does not define a population from
which to draw a prior-predictive ensemble. Consequently our fixed-truth
coverage and posterior-CDF diagnostics are not SBC and need not be uniform.
A future SBC check would need an explicitly proper test prior used by both
the simulator and inference. It would assess computation under that joint
model, alongside the fixed-truth performance assessment.

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

The original design estimates one mass from one vector pair. Competing methods
receive the same observations and declared noise information, except for the
explicitly named true-noise oracle in pooled-calibration experiments. That
oracle has additional information; plug-in versus integrated inference is the
comparison using the same observed calibration statistics. The original 16
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

Within an experiment, method differences are computed on paired trials.
Negative `flat_joint minus comparator` differences favor our equation for
smaller-is-better losses. Coverage is judged against its nominal target;
neither smaller width nor a smaller miss fraction alone ranks uncertainty
procedures. Comparisons across independently seeded experiments are not paired.
The test gives
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
Diagonal runs currently allow 10,000 trials per cell and 262,144 total trials.

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

The existing coordinate-calibration modes are `known`, `scaled` and
`estimated`. The pooled isotropic extension is described below.
`known` supplies the true coordinate SDs and is
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

Within these coordinate-calibration modes, every noise-aware estimator
receives the same supplied SDs. Noise-blind methods continue to ignore them.
The posterior treats supplied variances as fixed: this is **plug-in
calibration**, not inference that integrates over uncertainty in the
calibrated variances. Interval coverage measures the consequence of that
assumption. Finite coordinate calibration can also create unequal estimated
coordinate SDs when the true sensor is isotropic.

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

#### Pooled isotropic calibration uncertainty

The `pooled_isotropic` calibration model assumes one unknown scalar variance
for each channel, constant across its three coordinates, with independent
Gaussian errors. The force and acceleration variances are separate unknowns.
For `k >= 4` independent calibration vectors per channel, subtract that
channel's estimated calibration mean vector and pool the residual squared
lengths into `S`. Its sampling law is `S / sigma^2 ~ chi_squared(nu)`, where
`nu = 3 * (k - 1)`. This differs from estimating three independent coordinate
variances, each with `k - 1` degrees of freedom. The lower bound on `k` is
the implemented scope, not a claim that the probability model is undefined
for every smaller calibration sample.

The inference input is the stable measurement **mean plus external
calibration scatter**. The within-experiment scatter of repeated measurements
is not included. With unknown variance that scatter would contain additional
information, so this interface does not claim to reproduce inference from
all raw repeated readings. For each channel the conditional mean likelihood
is `Normal_3(latent_vector, sigma^2 / n * I)`.

The new variance integration uses the explicit scale-invariant reference
`d(sigma^2) / sigma^2` independently for the two channels, while holding the
physical magnitude/direction reference `df d(alpha) dOmega` fixed. This is a
declared extension of the model, not a uniquely required prior for the full
mass problem. Calibration alone gives
`sigma^2 | S ~ InverseGamma(nu/2, S/2)`, using density proportional to
`v^(-shape-1) exp(-scale/v)` for variance `v`. Integrating that variance out
of the three-dimensional Gaussian mean likelihood gives a multivariate
Student-t likelihood with `nu` degrees of freedom and scale matrix
`S / (nu * n) * I`. Its covariance is `S / ((nu - 2) * n) * I`, so scale and
variance must not be interchanged. The normal–gamma integration is described
by [Murphy (2007), §4](https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf);
the pooling and constrained-mass application here are our model-specific
derivation.

Each channel has **one** multivariate t likelihood. Its three coordinates
share the integrated variance; a product of independent scalar t densities
would be a different model. The mass law integrates the likelihood against
the variance distribution before normalization over compatible latent
pairs. Averaging already normalized fixed-variance mass laws using only
calibration weights would omit their evidence weights. Scalar channel
standardization must also retain the factors implied by the fixed physical
reference measure when those scales vary during integration.

The comparison separates a known-SD **oracle**, which receives the true
channel variances, a **plug-in** method using `S/nu` as fixed variance, and
the **integrated** method using the same observed means and calibration
statistics with variance uncertainty propagated. The oracle is a benchmark
with additional information. The plug-in/integrated pair isolates how the
same calibration information is used. The integrated point retains the
ratio-of-magnitude-means decision rule; its interval is a conditional mass
interval, whose fixed-truth coverage remains an empirical assessment.

Calibration is redrawn independently for each simulated dataset. Sharing
one calibration across many future measurements is a separate extension,
not implemented by this model. With `B` independent calibration batches and
`J` measurements per batch, a future study would assess conditional behavior
within batches and use the standard deviation of batch mean scores divided
by `sqrt(B)` for the overall MCSE. Paired method differences require the same
batch-level calculation; treating all `B * J` scores as independent would
understate uncertainty when the shared calibration induces dependence.

#### Calibration-integration pilot

The fixed pilot in
[`studies/calibration_integration_pilot.json`](studies/calibration_integration_pilot.json)
compares pooled calibration with `k=4` and `k=30`, four stable readings,
three true masses and four acceleration SNRs including zero. Each cell has
512 independent deployment replicates. All methods within a cell share the
same observations; the two calibration-size experiments use independent
seeds. The protocol records the assumptions and Monte Carlo precision targets.

Run it from the workspace root with:

```powershell
python -B 03_trace_and_evaluator/mass_estimator/investigation/study_runner.py --study 03_trace_and_evaluator/mass_estimator/investigation/studies/calibration_integration_pilot.json --output .tools/mass_estimator_studies/my-calibration-study
```

The calibrated method and the oracle are available throughout the archive,
results browser and completed-study analysis. Paired uncertainty summaries
include coverage, both miss directions and log-width, as well as proper
scores. These use differences on the same trials, not differences between
independently estimated standard errors. The null diagnostic discovers all
archived probability laws and requires complete interval endpoints for each.

Numerical verification includes an independent two-precision mixture,
direct radial integration, analytic zero-data moments, channel exchange,
large changes of units, the known-noise limit and separate refinement of the
mass and calibration quadratures. The integrator works in calibration-relative
mass units and restores physical units afterward. This avoids a numerical
integration window that changes its meaning when units change.

The completed evidence is
`.tools/mass_estimator_studies/assessment-calibration-integrated-20260920-verified`.
The earlier directory without `-verified` omitted paired coverage summaries;
the corrected rerun reproduced all archived arrays exactly. Use the verified
run for analysis and do not count the two directories as independent evidence.
See the [current synthesis](README.md#what-the-completed-simulations-add) for
results and their limits.

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

That completed 46-cell comparison did not assess correlated noise,
calibration uncertainty integrated into inference, repeated pairs with varying
excitation, shared sensor offsets, non-Gaussian noise, or real robot data.
Those require changes to the experiment, rather than alternative scores for
its saved outputs. The pooled-calibration extension above defines one such
change separately.

## Designed studies and reusable trial outputs

The first assessment protocol is
[`studies/mass_excitation_pilot.json`](studies/mass_excitation_pilot.json).
It specifies aims, data generation, estimands, methods, performance measures,
sources and interpretation before execution. This is an exploratory pilot:
three masses crossed with four excitation levels and two stable-reading
designs, with 512 independent experiments per cell. The fixed replication
count is separate from the number of readings within each experiment.
Primary bounded losses and interval coverage have stated Monte Carlo precision
targets; the pilot estimates a later study's required size and never stops or
extends itself based on observed results.

### Explicit mass and excitation, including no excitation

The new scenario type keeps physical mass independent of excitation:

```json
{"scenario": {"type": "mass_excitation", "masses": [0.25, 1, 4], "acceleration_snrs": [0, 0.5, 2, 8]}}
```

Every mass is crossed with every acceleration SNR. Acceleration magnitude is
the configured SNR times the true per-reading acceleration RMS noise SD; force
is mass times acceleration. Supplied calibration and the number of repeated
readings do not redefine these physical truths. At zero excitation both true
vectors are zero, while the stipulated positive mass remains available for
assessment. The data distribution then cannot distinguish masses. This is a
simulation of noisy observations from a true null, distinct from the numerical
fixture that supplies exactly zero measured vectors.

Mass scores and interval coverage still apply. Relative force-prediction error
divides by zero true-force energy and is explicitly `not_applicable`. The
magnitude-weighted joint loss is also marked inapplicable for ranking mass at
this boundary: its magnitude weight is zero for every candidate. Neither case
is scored as perfect performance or counted as a numerical failure. Existing
positive-signal scenario types retain their previous meaning and results.
The form offers this design as **Mass and acceleration excitation**.

Coverage summaries now also include the fractions with truth below and above
the reported interval. Those two fractions plus coverage equal one; each has
its own Monte Carlo error. The results browser exposes them separately.

### Run, checkpoint and resume a study

From the workspace root:

```powershell
python -B 03_trace_and_evaluator/mass_estimator/investigation/study_runner.py --study 03_trace_and_evaluator/mass_estimator/investigation/studies/mass_excitation_pilot.json --validate
python -B 03_trace_and_evaluator/mass_estimator/investigation/study_runner.py --study 03_trace_and_evaluator/mass_estimator/investigation/studies/mass_excitation_pilot.json --output .tools/mass_estimator_studies/assessment-pilot-20260920
```

Use a fresh output directory for a new study. To continue the same study, use
the same command with `--resume`. The checkpoint locks the normalized protocol,
Python/NumPy versions and production-source hashes. It verifies every saved
file of completed experiments before reusing them. A changed protocol or code
requires a new study directory. Failures remain explicit and are not retried
silently.

Checkpointing currently operates at **completed-experiment** granularity. An
interrupted experiment starts again in a new numbered attempt directory;
earlier attempts remain available. It does not yet resume individual inference
batches or append Monte Carlo replicates. Completed studies can be resumed as
a verification operation without fitting again.

The study writes `study_checkpoint.json` and `precision_summary.json` beside
its experiment directories. The latter reports actual MCSEs, pilot-variance
sample-size suggestions and conservative bounded-variance requirements. It
flags requirements beyond the engine's current sample limits. These are
planning quantities, not acceptable-error thresholds, stopping rules or
simultaneous confidence guarantees. Zero observed pilot variance is explicitly
flagged rather than taken as proof that no more simulation is needed.

Completed, unchanged experiment results under `.tools/mass_estimator_studies/`
appear in the existing results browser. No result is pooled across scenario
cells, which share noise.

### Rescore saved points without running inference

With `save_trials: true`, new runs save a versioned `trial_manifest.json` and
one checksummed NPZ archive per scenario. Archives preserve ordered replicate
IDs, training and held-out observations, supplied and true noise scales,
points and selected posterior summaries. The executed inference source files
are saved under `source/`; the manifest includes runtime and source provenance.
The loader verifies hashes, shapes and numeric types with pickle disabled.

For example, after the pilot:

```powershell
python -B 03_trace_and_evaluator/mass_estimator/investigation/trial_archive.py --run .tools/mass_estimator_studies/assessment-pilot-20260920/experiments/stable_n1/attempt-001 --output .tools/mass_estimator_studies/assessment-pilot-20260920/rescored-factor-3.json --factor-tolerance 3 --cap-factor 5
```

This computes original point losses and additional configurable bounded losses
and paired comparisons from saved estimates. It never refits a method, drops
invalid estimates or overwrites an existing output. Existing population-risk
and MCSE caveats travel with the saved evidence. The new scoring definition and
its source hashes are recorded separately from the original experiment.

The current rescoring command covers points. Although held-out readings are
saved for future prediction checks, arbitrary distributional rescoring or new
interval widths still require more than the saved posterior summaries: full
probability curves are not archived. Older aggregate-only and unversioned
archives are not silently interpreted as the new format.

This workflow applies the already researched ADEMP framework of
[Morris, White and Crowther](https://doi.org/10.1002/sim.8086), the experimental
design guidance of [Chipman and Bingham](https://arxiv.org/abs/2111.13737), and
the failure-handling distinctions of [Pawel et al.](https://arxiv.org/abs/2409.18527).
The protocol's choices are explicit study decisions, not consequences forced
by those sources or a claim that one estimator is generally superior.

### Comparing whole error curves

The completed 21 September archive analysis asks whether performance
differences persist as the allowed error and tail penalty change. Its current
[results](../../../.tools/loss_sensitivity_20260921/main-reviewed/analysis.json)
and [figure](../../../.tools/loss_sensitivity_20260921/main-reviewed/tolerance-comparison.png)
are scratch evidence; the [README](README.md#what-the-full-error-curves-show)
holds the interpretation. This exploratory analysis reuses the main study's
trials and adds no independent simulations.

It compares `flat_joint` against `flat_geometric`, `flat_median`, `norm_ratio`,
`norm_floor` and `positive_profile` at all 24 conditions. Invalid outputs
remain failures. The eight declared factors are 1.1, 1.25, 1.5, 2, 3, 5, 10
and 20; each is used both as a tolerance and as a squared-log-loss cap.
The latter is divided by \(\log^2 K\) to put bounded losses on a common
0–1 scale. This normalization preserves rankings at a fixed cap; it does not
turn different caps into the same criterion.

For absolute log error \(E=|\log(\widehat m/m)|\), the survival curve
\(S(t)=\Pr(E>t)\) gives failure probability at factor \(K=e^t\).
The analysis examines every empirical error breakpoint for survival-curve
extrema, including the finite-threshold failure mass of invalid estimates
represented by \(E=+\infty\). Plotting uses a finite display grid; the
crossing check does not depend on that grid. Absence of a statistically
resolved crossing does not establish population dominance, particularly in
sparsely observed tails.

Paired MCSE remains the local precision measure. For a simultaneous check,
the Dvoretzky–Kiefer–Wolfowitz–Massart inequality bounds each marginal CDF
uniformly in its threshold; see [Reeve (2024)](https://arxiv.org/html/2403.16651v1).
Applying a union bound over **all archived main-study methods**, including
those omitted from this comparison, gives

$$
\epsilon=\sqrt{\frac{\log(2L/.05)}{2R}},\qquad
L=12\times11+12\times16=324,\quad R=10{,}000.
$$

Each survival-difference curve is then within \(2\epsilon=0.04351924\)
of its population curve, simultaneously over this family and every finite
threshold, with probability at least 95%. Shared observations between methods
or conditions do not invalidate the union bound; independent replicates
within each marginal experiment are required. The family includes all
archived methods to cover selecting a subset after seeing earlier results.
Formally, a monotone transform such as \(2\arctan(E)/\pi\) maps the invalid
atom to 1, permitting the ordinary real-variable inequality to be applied.

The same bound also controls normalized capped squared-log differences
uniformly over every cap, because for \(c=\log K>0\),

$$
E[\min(E^2/c^2,1)]=\int_0^c\frac{2t}{c^2}S(t)\,dt.
$$

The positive weights integrate to one. This integral consequence and the
family correction are our application of the inequality. A raw capped-loss
bound multiplies \(2\epsilon\) by \(c^2\). These conservative bounds quantify
Monte Carlo uncertainty for the fixed implemented procedures; they do not
cover quadrature bias, model error, other truth conditions, or application
adequacy. A genuine ordering of the population absolute-log-error curves
would support increasing losses of that error, not arbitrary asymmetric or
ordinary mass-error objectives.

Reproduce the numeric analysis from the workspace root, choosing a new output
directory:

```powershell
python -B .tools/loss_sensitivity_20260921/analyze.py --study .tools/mass_estimator_studies/assessment-main-20260920 --output .tools/loss_sensitivity_20260921/my-repeat
```

The script uses the existing strict study/archive loaders and production
scoring functions. It verifies all 1,440 overlapping saved paired summaries,
checks the survival-integral identity independently on finite and invalid
examples, and records source hashes and its source snapshot. The reviewed
result uses the actual 324-CDF inventory; the preliminary `main/` scratch
output used only the six selected methods in its family correction and is
not the result linked by the synthesis. Both use identical saved estimates.
No simulation or inference module was changed for this analysis.

## Main assessment and calibration sensitivity

The completed main fixed-size protocol is
[`studies/mass_excitation_main.json`](studies/mass_excitation_main.json):
the same 24 mass/excitation/reading-count conditions as the pilot, with
10,000 independent simulated experiments per condition and fresh seeds.
The study phase is `main_assessment`; this describes the study's role and
does not imply a confirmatory test or a universal estimator ranking.

The primary point criteria remain factor-two failure and capped squared log
loss. They express a multiplicative error preference, not an application
requirement that the user has supplied. Tighter tolerance curves remain
available. Comparisons use the same data within each condition and preserve
the directly calculated paired Monte Carlo standard error (MCSE).

The replication count has a conservative rationale independent of favorable
pilot variance estimates. A probability estimated from 10,000 independent
Bernoulli observations has population MCSE at most 0.005; a paired difference
of two probabilities has MCSE at most 0.01. The corresponding bounds for
individual and paired capped squared log loss are approximately 0.002403 and
0.004805. These are pointwise Monte Carlo precision bounds, not estimator
accuracy thresholds or simultaneous confidence guarantees. Empirical MCSE
uses a sample variance and may differ slightly from its population bound.

Coverage is examined together with log-interval width and both miss
directions, at 50% as well as 95%. The zero-excitation diagnostic reuses the
saved null observations and interval endpoints at masses
`[1/64, 1/16, 1/4, 1, 4, 16, 64]`. This is legitimate because the stipulated
null observation law does not depend on mass. The diagnostic verifies equal
observations and estimates across the original null cells. Its truth points
reuse the same trials and do not add independent simulations. Broad coverage near
one mass scale does not imply mass identification or calibrated intervals
across the positive mass domain.

The separate
[`studies/calibration_sensitivity_pilot.json`](studies/calibration_sensitivity_pilot.json)
holds the physical law and true Gaussian noise fixed while varying supplied
noise information. It has five calibration regimes: known SDs, both channels
at 0.75 times the true SD, both at 1.25, unequal channel factors, and estimated
coordinate SDs from 30 independent calibration readings. Each regime has
four stable readings, the same 12 truth conditions, and 1,024 Monte Carlo
replicates per condition. Across regimes the seeds are independent, so those
comparisons are not paired. Estimated calibration is redrawn per replicate;
inference plugs it in without propagating its uncertainty. Correlated noise,
non-Gaussian data, changing excitation and one calibration shared by many
measurement experiments are still outside these protocols.

Run either study with `study_runner.py --study <protocol> --output <new-dir>`.
The completed local output directories are
`.tools/mass_estimator_studies/assessment-main-20260920` and
`.tools/mass_estimator_studies/assessment-calibration-20260920`.
The expanded capacity retains the 262,144 total-trial guard for diagonal and
measurement experiments. The RMS control now respects `batch_size`, including
per-replicate calibration scales, rather than allocating a whole-cell
quadrature operation.

The analysis tools consume completed results independently of the currently
installed inference code. They verify saved evidence and record their own
analysis provenance, so old studies can be analyzed after code development.
Their plots show MCSE explicitly; missing metrics and invalid methods remain
visible. The study resume operation still requires the original execution-source
fingerprint because it can execute unfinished experiments. This fingerprint
covers the declared inference/configuration/scoring dependencies, the study
runner and the production estimator dependency. Analysis and browser source
changes are tracked separately and do not block study resumption.

The reusable analysis commands are:

```powershell
python -B 03_trace_and_evaluator/mass_estimator/investigation/assessment_summary.py --study .tools/mass_estimator_studies/assessment-main-20260920 --output .tools/assessment_main_20260920/main-analysis
python -B 03_trace_and_evaluator/mass_estimator/investigation/null_truth_sweep.py --study .tools/mass_estimator_studies/assessment-main-20260920 --output .tools/assessment_main_20260920/main-null
```

Choose a new output directory when an analysis already exists. Both tools
verify saved evidence without refitting. The first exports JSON and PNG/SVG
figures; `--no-figures` produces JSON without Matplotlib. The null tool exports
JSON and a coverage plot, and accepts an explicit increasing positive mass
grid through `--masses`. Only true-zero-excitation archives qualify for that
reuse: arbitrary nonzero truth changes would require a new simulated dataset.

Matplotlib is an optional plotting dependency, not an inference dependency.
The local bundled runtime uses the workspace installation at
`.tools/assessment_plot_dependencies`. Append that directory to `sys.path`
after its existing entries when exporting figures, preserving the inference
runtime's existing NumPy version. Do not prepend this directory as a replacement
simulation environment. Analysis outputs record their own runtime versions.

### Operating-range assessment

[`studies/operating_range_followup.json`](studies/operating_range_followup.json)
specifies a fresh follow-up to the archive sensitivity analysis. It crosses
effective force and acceleration SNR values `[1,2,3,4,6,8,12,16]`, retaining
all 64 cells and 10,000 replicates per cell. Both channels have independent
known isotropic Gaussian noise with coordinate SD 1 and a common latent
direction unknown to inference. The physical mass in these units is the
force SNR divided by acceleration SNR. Off-diagonal comparisons change mass
as well as signal balance; constant-ratio rays hold mass fixed.

The experiment uses one effective reading. For independent stable Gaussian
repeats with known SDs, the sufficient mean has SNR equal to the per-reading
SNR times the square root of the reading count. This lets the point-estimation
map inform that repeated-reading experiment. It does not reproduce a different
held-out prediction task, shared calibration, or changing excitation.

Absolute-accuracy benchmarks are 90% and 95% success within factors 1.25,
1.5 and 2. These are explicit reference targets, not application requirements.
All archived methods enter the simultaneous DKW family described above;
with 64 cells and 11 methods, the marginal halfwidth is about 2.263 pp and
the paired-difference halfwidth about 4.527 pp. A target is:

- **supported** when the simultaneous lower success bound reaches it;
- **below target** when the upper success bound is below it;
- **unresolved** when the bounds straddle it.

The study retains same-law median and geometric readouts as well as the norm
ratio, separating a benefit shared by the conditional law from an advantage
of the particular ratio-of-means summary. All failures remain in denominators.
No interpolation or monotonic extension from passing grid cells is warranted.
True-SNR maps do not supply an acceptance rule for observed readings.

Run a new copy of the study, or analyze the completed local run, with:

```powershell
python -B 03_trace_and_evaluator/mass_estimator/investigation/study_runner.py --study 03_trace_and_evaluator/mass_estimator/investigation/studies/operating_range_followup.json --output .tools/mass_estimator_studies/my-operating-study
python -B 03_trace_and_evaluator/mass_estimator/investigation/assessment_summary.py --study .tools/mass_estimator_studies/operating-range-20260921 --output .tools/operating_range_20260921/my-accuracy-analysis --accuracy-map
```

The reusable accuracy summary reads verified saved results, adds no fitting,
and preserves missing scores as unavailable. Its JSON retains each truth
condition, method, paired MCSE, simultaneous bounds and target classification.
Changing the requested factors or success targets changes the analysis rather
than the fitted estimates; supported factors must exist in the saved scores.

### Prior factors, switch rules, dominance and sequential updating

These 23 September modules are documented, with results, in
[working contribution 12](prior_factors_and_sequential_updating.md). They read
verified trial archives and never refit an archived method:

- **`paired_dominance.py`** finds cells where a focal rule beats every other
  saved rule. It uses paired differences and a two-sided Bonferroni normal
  bound over cells × factors × rules.
- **`switch_rules.py`** rescores two-branch pretest rules and their true-cell
  oracle.
- **`prior_weight_study.py`** recomputes ratio-of-means points from saved
  observations under declared mass factors (`sech`, `lognormal_factor`,
  `lognormal_prior`), then scores, fits (mean over factors of the largest cell
  regret) and summarizes them. Options: `--half first|second` for disjoint
  replicate halves, `--cells` for subsets, and `--min-snr` in `fit` and
  `summarize`.
- **`sequential_update.py`** splits a saved cell into consecutive series and
  compares the symmetric, eq. (18) and eq. (20) rules for combining readings.

`weighted_posterior.py` supplies the weighted law and its continuous
summaries; `archive_rescoring.py` supplies the shared loading and paired
statistics. The fresh studies are
[`studies/prior_update_validation_grid.json`](studies/prior_update_validation_grid.json),
[`studies/prior_update_offgrid_shifted.json`](studies/prior_update_offgrid_shifted.json)
and [`studies/prior_update_mass_ladder.json`](studies/prior_update_mass_ladder.json).

```powershell
$S = "03_trace_and_evaluator/mass_estimator/investigation"
$RUN = ".tools/mass_estimator_studies/prior-update-validation-20260923/experiments/validation_snr_grid/attempt-001"
python -B $S/paired_dominance.py --run $RUN --exclude tube_joint --output .tools/my-analysis/dominance.json
python -B $S/switch_rules.py --run $RUN --output .tools/my-analysis/switch.json
python -B $S/prior_weight_study.py sweep --run $RUN --family sech --lams 0 0.35 --output .tools/my-analysis/sweep.json
python -B $S/prior_weight_study.py summarize --sweep .tools/my-analysis/sweep.json --output .tools/my-analysis/summary.json
python -B $S/sequential_update.py --run <attempt directory> --cells 07_mass_excitation --output .tools/my-analysis/series.json
```

Outputs refuse to overwrite existing files. Sweep, dominance, switch and
sequential outputs record source hashes, runtime and archive manifests; fit
and summarize outputs record their input sweeps and source hashes. A sweep over 64 cells of 10,000 readings takes about a
minute; a sequential analysis about 14 s per cell. Every result is conditional
on the isotropic known-noise Gaussian model, and the modules reject
anisotropic supplied SDs.

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
