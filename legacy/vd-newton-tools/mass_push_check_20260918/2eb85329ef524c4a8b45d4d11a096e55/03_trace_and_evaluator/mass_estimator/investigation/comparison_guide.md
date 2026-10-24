# Refined mass-estimator comparison

This test compares **11 point estimators on the same noisy observations**, using
several explicitly different meanings of accuracy. It also scores the two
available probability laws. The original full-covariance estimator is unchanged.
This benchmark is the independent, isotropic, unknown-direction 3D specialization.
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

Each trial estimates one mass from one vector pair. All methods receive the
same observations and the same declared unit noise scales. The original 16
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

The comparison does not yet establish behavior for correlated/anisotropic
noise, learned error scales, repeated-pair estimators, non-Gaussian noise or
real robot data. Those are changes to the experiment, not alternative scores
for the present one.

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
