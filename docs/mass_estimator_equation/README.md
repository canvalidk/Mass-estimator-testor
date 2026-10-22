# The mass estimator equation

**Specification date: 21 September 2026.** This folder specifies the current
single-pair, Gaussian, flat-magnitude estimator and provides a standalone
implementation. The point is the **ratio of the mean compatible magnitudes**.

Open **[mass_estimator.pdf](mass_estimator.pdf)** for the presentation.
The compact equation is also available as **[SVG](equation.svg)** and
**[PNG](equation.png)** for inserting into other documents.

## The equation

For measured net force and acceleration, form the column vectors

$$
\mathbf y=\begin{bmatrix}\widetilde{\mathbf F}\\\widetilde{\mathbf a}\end{bmatrix},
\qquad
\mathbf x(f,\alpha,\mathbf u)=\begin{bmatrix}f\mathbf u\\\alpha\mathbf u\end{bmatrix},
\qquad
Q=(\mathbf y-\mathbf x)^T\Sigma^{-1}(\mathbf y-\mathbf x).
$$

The full compact formula is

$$
\boxed{
\widehat m(\widetilde{\mathbf F},\widetilde{\mathbf a},\Sigma)
=\frac{E_P[f]}{E_P[\alpha]}
=\frac{\displaystyle\int_{\mathcal C} f\,e^{-Q/2}\,d\nu}
        {\displaystyle\int_{\mathcal C} \alpha\,e^{-Q/2}\,d\nu}
}
$$

with the domain, measure and probability law fixed by

$$
\mathcal C=(0,\infty)^2\times S^{d-1},\qquad
d\nu=df\,d\alpha\,d\Omega(\mathbf u),\qquad
dP=\frac{e^{-Q/2}}{Z}\,d\nu,\qquad
Z=\int_{\mathcal C}e^{-Q/2}\,d\nu.
$$

Here $S^{d-1}$ is the unit sphere in $d$ dimensions and $d\Omega$ is its
uniform probability measure. In one dimension the two directions, $+1$ and
$-1$, each have weight $1/2$. Flat means flat in the two **physical positive
magnitudes** $f$ and $\alpha$. It does not mean flat in mass, log mass or
Cartesian vector volume. The unnormalised reference measure has infinite
total volume, but $Z$ and both magnitude means are finite and positive for
finite observations and positive-definite Gaussian covariance.

The operation averages possible force and acceleration magnitudes under one
joint law, then divides. Direction influences the weight of each pair. It
does not divide vectors, take the ratio of measured vector lengths, or compute
$E_P[f/\alpha]$.

## Types and inputs

$$
\operatorname{estimate\_mass}:
\mathbb R_F^d\times\mathbb R_a^d\times\mathrm{SPD}(2d)
\longrightarrow\mathrm{MassEstimate},\qquad d\in\{1,2,3\}.
$$

| Symbol / Python input | Type and meaning |
|---|---|
| $\widetilde{\mathbf F}$ / `force` | Finite measured **net force vector**, shape `(d,)`; SI unit N. Components may have either sign or be zero. |
| $\widetilde{\mathbf a}$ / `acceleration` | Finite measured acceleration vector, shape `(d,)`, in the same axes; SI unit m/s². |
| $\Sigma$ / `covariance` | Symmetric positive-definite matrix, shape `(2*d, 2*d)`, for the measurement errors in `[F1,...,Fd,a1,...,ad]` order. |
| $f,\alpha$ | Integration variables: positive **latent magnitudes**, with force and acceleration units. These are inferred, not supplied. |
| $\mathbf u$ | Integration variable: a common dimensionless unit direction. |
| $\widehat m$ / `result.mass` | A strictly positive scalar in force-unit / acceleration-unit; N divided by m/s² gives kg. |
| `probability` | Probability content of the central equal-tail interval; default `0.95`. |

The force block of $\Sigma$ has units N², the acceleration block (m/s²)²,
and the cross-block N·m/s². Diagonal entries are **variances**, not standard
deviations. For independent isotropic channels in 3D,

$$
\Sigma=\operatorname{diag}(\sigma_F^2,\sigma_F^2,\sigma_F^2,
                          \sigma_a^2,\sigma_a^2,\sigma_a^2).
$$

Each $\sigma$ is a **per-coordinate** standard uncertainty, not the uncertainty
of the vector norm. Supply the full matrix when coordinate errors or channels
are correlated. Independent provenance does not imply independent errors.

## Basic assumption and declared model

**Newton II holds for the true quantities of this same object:**
$\mathbf F^*=m\mathbf a^*$ with $m>0$. The force and acceleration must concern
the same object, time or interval, and inertial frame. Their determination
must not already assume the mass being sought.

The implementation additionally declares zero-mean jointly Gaussian
measurement errors with supplied covariance $\Sigma$, a uniform unknown
common direction, and flat positive magnitudes. These statistical choices
are part of the equation's specification; Newton II alone does not select
them. Observed nonalignment is handled through measurement uncertainty.
The conditional calculation has no empirical rejection threshold.

## Uncertainty: exactly what is returned

Use the **same** law $P$ to define the random scalar $M=f/\alpha$. Let
$q_p=\inf\{t>0:P(M\le t)\ge p\}$. For requested content $c$, return

$$
I_c=[q_{(1-c)/2},q_{(1+c)/2}],\qquad
s_{\log M}=\sqrt{\operatorname{Var}_P[\log(M/m_0)]},\quad m_0>0.
$$

The log spread is dimensionless and independent of the chosen reference mass
$m_0$. The code uses $m_0=s_F/s_a$, where $s_F,s_a$ are the square roots of
the mean diagonal variance in their respective covariance blocks.

`result.interval` and `result.log_standard_deviation` are these two outputs.
A 95% interval contains 95% of the **declared conditional mass law**. It is not
automatically a 95% repeated-experiment confidence interval. The point is a
different summary and need not equal the interval midpoint or median, or lie
inside an arbitrarily narrow requested interval.

Under this finite-noise Gaussian/flat model, the ordinary mean $E_P[M]$ and
second moment diverge. Consequently an ordinary mass SD is undefined; a
symmetric `mass +/- SD` is not an exact uncertainty description. The central
quantiles and log spread are finite. Log spread does not assert a lognormal law.

Inverse mass is derived from the same result:
$\widehat\mu=1/\widehat m$ and $I_\mu=[1/I_c^{\rm upper},1/I_c^{\rm lower}]$.
These are `result.inverse_mass` and `result.inverse_mass_interval`.

## How to use the Python function

Use Python 3.10 or newer. Install the dependency and run the example from this
folder:

```sh
python -m pip install -r requirements.txt
python example.py
```

```python
import numpy as np
from mass_estimator import estimate_mass

force = [10.0, 0.0, 0.0]       # N
acceleration = [2.0, 0.1, 0.0] # m/s^2
covariance = np.diag([0.5**2] * 3 + [0.1**2] * 3)

result = estimate_mass(force, acceleration, covariance)
print(result.mass)
print(result.interval)
print(result.log_standard_deviation)
```

For this example, the rounded outputs are **4.993754 kg**, a **95% conditional
interval [4.344146, 5.740363] kg**, and **log SD 0.071035**. At the default
settings, the final scaled numerical refinement change is about `6.61e-5`.
That diagnostic measures numerical stability, not measurement uncertainty.

1. Establish the object, frame, interval, measurement units and covariance.
2. Pass the measured vector components directly, retaining their directions.
3. Report the point together with the interval, its probability, and the model.
4. Retain the original readings, covariance and their provenance with the result.

**Numerics.** The function integrates radius analytically, then integrates
direction and ratio angle numerically. Independent isotropic 3D errors have
an exact direction reduction as well. It covers the full positive mass
domain. Quantiles invert a continuous piecewise-linear angle density; endpoint
transformations resolve the finite log moments. No positive mass cutoff is
introduced. The negative radial-tail quadrature uses the existing negligible
scaled-tail approximation documented in the code.

Every successful call compares at least two grid resolutions. Starting
`direction_order=24` and `ratio_order=512` are doubled as needed, up to
`max_refinements=3`; `rtol=2e-4` applies to all returned numerical summaries.
`numerical_change`, `direction_order` and `ratio_order` report the last check
and final resolution. The direction order is unused in the analytic 3D,
known-direction and 1D direction calculations. Refinement is an observed
stability check, not a certified error bound. Highly concentrated or badly
conditioned cases need particular care and independent resolution checks.
Invalid inputs raise `ValueError`; unresolved quadrature or numerical range
failures raise `ArithmeticError` rather than returning an unchecked number.

## Scope and special cases

- **Known direction:** if externally established, pass an oriented unit vector
  as `known_direction`. This replaces $d\Omega$ by a point mass at that
  direction. The function integrates only the two magnitudes. Do not infer
  this direction from the same noisy readings and treat it as exact.
- **Measured zeros:** are accepted with positive-definite covariance. At two
  measured zero vectors and independent isotropic errors, $\widehat m=
  \sigma_F/\sigma_a$; $M/(\sigma_F/\sigma_a)$ is half-Cauchy. Its central
  95% interval is approximately `[0.039290, 25.451700]` times that scale and
  its log SD is $\pi/2$. This apparatus-dependent output does not establish
  the object's mass. At true zero excitation the experiment cannot identify mass.
- **Exact quantities:** are a separate checked operation, outside this
  function. Exact nonzero codirectional vectors identify their magnitude
  ratio; exact zero-zero identifies no mass; other exact pairs contradict
  positive-mass Newton II. Zero covariance is not an empirical input.
- **Stable repeats:** independent readings of one unchanged latent vector
  pair with the same known covariance can be replaced by their vector mean
  and covariance $\Sigma/n$. Different excitations sharing a mass require
  a separate joint model. Averaging separately estimated masses is not that model.
- **Calibration uncertainty:** this interface conditions on the supplied
  covariance. It does not propagate uncertainty in that covariance. The recent
  pooled-calibration extension keeps the ratio-of-means readout but supplies a
  different likelihood and extra calibration inputs; it is not silently
  substituted here. Hard bounds and non-Gaussian models are also outside scope.

## Version check and provenance

**Selected form:** the compatible-pair ratio-of-means estimator, under the
7 September law-assumed-correct steering, still identified as `flat_joint`
in the local investigation synthesis dated 21 September. The earlier
best-point/stationary-root estimator and the later median proposal are the
two prior stages checked. Neither is the function delivered here. New
comparison summaries, alternative measures and calibration models do not
constitute a replacement of the selected point formula.

Managed `main` was resolved once on 21 September 2026 to
[`canvalidk/VD-docs@8782a293297330f7a354d87ebe62792018fa8866`](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866).
All managed reads used that snapshot. Selection began with
[`catalog.yaml`](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/catalog.yaml)
and the mass-estimator discussion guidance, then followed formula, readout,
measure, uncertainty and steering references. Principal sources:

- [Full working specification, especially sections 7-12 and 22-23](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md): current compatible-pair equation and its relation to the two earlier estimators.
- [Law-assumed-correct steering, 7 September](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md): governs the conditional purpose over the older specification's rejection gate.
- [Readout audit](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md): ratio of means versus induced median; point and interval are different functionals.
- [Additional derivations, 5 September](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_estimator_additional_derivations_2026-09-05.md): null law, scale behavior and uncertainty.
- [Earlier median proposal](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/median_based_estimator_inertial_mass.md): inspected for historical identification, not adopted.

Active local sources are in `03_trace_and_evaluator/mass_estimator/` within
VD-Newton. Local HEAD was `560d0ba2ccb5b44ec2a76f4f52cb8095b35aefec`, but several
investigation files are modified or untracked; that commit alone does **not**
identify the active contents. SHA-256 prefixes below identify the bytes read:

| Local source, relative to the laboratory | SHA-256 prefix | Use |
|---|---|---|
| `estimator.py` | `950786e8b6a71216` | General-covariance integration reference. |
| `practical_formulas.py` | `766bb8cca3458f2d` | Quantiles and endpoint-safe log spread. |
| `investigation/README.md` | `fc100798bab4a8ab` | Current synthesis; moment divergence and unchanged point. |
| `investigation/comparison_guide.md` | `5fe3de058bd63faf` | Current method identity, covariance and calibration scope. |
| `investigation/comparison_posterior.py` | `1a5d1c896c1f4add` | Recent analytic 3D reduction and independent numerical comparison. |
| `investigation/gaussian_conditioning.md` | `3cbff40230095e02` | Why Gaussianity does not alone select the reference measure. |
| `investigation/why_flat_conditioning.md` | `c2294d5b1a5e1b82` | Flat law has a conditioning interpretation under added assumptions. |

The local developments extend the managed mathematical analysis; they do not
replace the core readout. This was a targeted version check, not an exhaustive
inbox or literature review. Unsanctioned reports and literature-review files
were not used as current discussion authority.

**Verification:** nine checks passed: analytic null cases in 1D/2D/3D,
correlated-null quantiles and log spread, known-direction truncated-normal
means, agreement with the latest isotropic backend and the original full
covariance implementation, an independent integral in the original
magnitude/direction coordinates, rotations, unit changes, channel reciprocity,
anisotropic 3D inputs, invalid inputs and explicit numerical nonconvergence.
These verify calculation of the stated model, not experimental accuracy or
confidence coverage. Verification scripts and rendering intermediates stay
in workspace scratch space, outside this presentation folder.

**Submission status:** created at the user's request in
`pass to VD-docs/Mass estimator equation/`; awaiting transfer. Saving or
publishing this folder in VD-Newton does not itself submit it to VD-docs or
edit the managed collection. The catalogued
[8 September receipt](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/outbox/VD-Newton/PASS_RECEIPT_2026-09-08.md)
was checked; its closure is already documented in the workspace guidance.
