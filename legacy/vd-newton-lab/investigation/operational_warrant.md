# Choosing the mismatch and understanding its mass uncertainty

**Working contribution 7 — 15 September 2026.** Follows the composition
extension of [why_flat_conditioning.md](why_flat_conditioning.md). The question
is whether the remaining mismatch requirements follow from the estimator's
purpose, and what choosing them does to mass uncertainty.

## Findings

1. The stated mass-composition and Gaussian-error premises do not select a
   unique conditioning rule. The counterexamples are explicit; an additional
   operational requirement is needed to choose one.
2. Preserving additive balance error motivates the original area rule.
   Scaling discrepancy by local Gaussian error motivates the tube rule.
   Using exact candidate sampling covariance gives another rule. Each has a
   mathematical warrant; none thereby gains a general calibration guarantee.
3. We can describe a substantial uncertainty difference exactly. At equal
   standardized observed vector lengths, all three rules return the same
   mass point. Positive alignment makes the tube interval narrower than the
   original; negative alignment makes it wider. For perpendicular vectors
   their mass laws coincide exactly with the zero-reading half-Cauchy law.
4. A one-dimensional integration formula reproduces these results independently
   of the earlier radial-and-angle integration. No baseline policy is changed.

All calculations assume Newton II and the declared measurement-error model.
No across-trial population of excitations is introduced.

## 1. What an operational requirement would add

The previous contribution established a sufficient route:

- Preserve force additivity while eliminating mass from the balance error.
- Require rotation covariance, continuity, zero discrepancy on compatible
  pairs, and equal discrepancy magnitude under standardized channel exchange.
- Shrink a constant tolerance in that discrepancy norm.

These requirements give cross-product area and hence flat compatible
magnitudes. The first step extends the additive structure of Newton's balance
equation to errors away from compatibility; the last chooses what it means
to hold those errors equally small.

The managed specification explicitly declares uniform common direction and
flat positive magnitudes as its current policy. Its Gaussian likelihood and
the readout composition laws do not assert these stronger mismatch premises.
Conversely, the readout audit explicitly discusses noise-cone volume as a
possible alternative. The source documents therefore do not close this gap
in favor of fixed area.

This is an identified freedom, not a missing algebraic step in the current
proof. The conditional derivation remains useful: it states what adopting
the original rule commits us to.

## 2. Three precise discrepancy scales

For this comparison use independent isotropic three-dimensional channels,
standardized by their known component uncertainties. Write possible vectors
as \(X,Y\), their compatible magnitudes as \(x,y>0\),
\(r^2=x^2+y^2\), and area as \(\mathcal A=\|X\times Y\|\).
The positive-dot condition remains part of every limiting construction.

| Discrepancy held below a shrinking constant | Magnitude reference weight | Warrant |
|---|---|---|
| \(\mathcal A\) | \(1\) | Additive algebraic mismatch with fixed norm tolerance |
| \(\mathcal A/\sqrt{\|X\|^2+\|Y\|^2}\) | \(r^2\) | Local Gaussian error scale; noise-metric tube |
| \(\mathcal A/\sqrt{\|X\|^2+\|Y\|^2+2}\) | \(r^2+2\) | Candidate sampling-covariance field described below |

All three posterior laws are proper. Their nonconstant radial weights cannot
be dropped as if they were overall normalizing constants.

### Local Gaussian error scale

At compatible candidate means \(X=xu,Y=yu\), perturb both channels by
independent unit Gaussian errors \(\xi,\eta\). To first order,

$$
\delta B=\xi\times(yu)+(xu)\times\eta,
\qquad \operatorname{Cov}(\delta B)=r^2(I-uu^T).
$$

Its two transverse components therefore have standard deviation \(r\).
Measuring discrepancy in those local noise units gives \(\mathcal A/r\).
The constant-width region in the noise metric retains the surface factor
\(r^2\), producing the tube measure. This agrees with the earlier joint
Jeffreys-volume calculation.

This construction is invariant under regular re-expressions of the local
constraints when their covariance is transformed too. For independent
constraint coordinates \(c\), covariance \(V\), and nonsingular matrix
\(H\), the quadratic form is unchanged under \(c'=Hc,V'=HVH^T\).
That supplies a structural reason to scale residuals by their covariance.
It does not make a finite-noise residual Gaussian or select a universally
optimal mass prior.

### Exact candidate sampling covariance is a different scale

At arbitrary fixed means \(X,Y\), direct expansion of the noisy cross product
gives

$$
C(X,Y)=\operatorname{Cov}[(X+\xi)\times(Y+\eta)]
=(\|X\|^2+\|Y\|^2+2)I-XX^T-YY^T.
$$

The \(2I\) term comes from \(\xi\times\eta\). This matrix is positive
definite. The candidate residual \(B=X\times Y\) is perpendicular to both
means, so

$$
B^TC(X,Y)^{-1}B
=\frac{\mathcal A^2}{\|X\|^2+\|Y\|^2+2}.
$$

Using this candidate-dependent norm defines the third row. At compatibility
its transverse scale is \(\sqrt{r^2+2}\), so the limiting weight is \(r^2+2\).
A polar-cap bound times the Gaussian envelope proves the limit globally,
including the axes. This is a concrete alternative, not an adopted tuning.

**The candidate qualifier matters.** The covariance of the cross product
under the supplied Gaussian uncertainty cloud, centered at observed vectors
\(U,V\), is instead the single fixed matrix \(C(U,V)\). Using that fixed
matrix defines another construction, generally with a direction-dependent
weight. The field \(C(X,Y)\) uses the sampling model at each possible pair;
it is not the covariance of one fixed input cloud at every point.

Nor does this imply that observed cross products at a fixed true mean can be
fully whitened by dividing their norm by \(\sqrt{r^2+2}\). The covariance
at a compatible true pair has transverse variance \(r^2+2\) and longitudinal
variance 2. The simplification above applies to the candidate residual with
its own candidate covariance. No exact chi-square pivot or coverage claim
follows from it.

## 3. An exact reduction that reveals how the mass laws differ

Let \(U,V\) now denote the observed standardized vectors. Put

$$
x=r\sin\theta,\quad y=r\cos\theta,\quad
0<\theta<\pi/2,\quad M/s=\tan\theta,
\qquad s=\sigma_F/\sigma_a.
$$

The observations enter the radial integral through

$$
h(\theta)^2=\|\sin\theta\,U+\cos\theta\,V\|^2.
$$

After integrating the shared direction, the likelihood is proportional to
\(e^{-r^2/2}\sinh(rh)/(rh)\); the omitted observed-data constant is common
to all parameters. Define

$$
J_1(h)=\int_0^\infty r e^{-r^2/2}\frac{\sinh(rh)}{rh}\,dr
=\sqrt{\frac\pi2}\frac{e^{h^2/2}\operatorname{erf}(h/\sqrt2)}h,
\qquad J_1(0)=1.
$$

The corresponding radial-cubic integral is exactly

$$
J_3(h)=(1+h^2)J_1(h)+1.
$$

For example, writing \(A(h)=hJ_1(h)\) gives \(A'=hA+1\), and the
radial-cubic integral equals \(A''/h\). Thus the normalized mass-angle
densities are proportional to:

| Reference | Unnormalized \(\theta\)-density |
|---|---|
| Flat | \(J_1(h)\) |
| Tube | \(J_3(h)\) |
| Candidate covariance | \(J_3(h)+2J_1(h)\) |

In particular,

$$
\frac{p_{\rm tube}(\theta)}{p_{\rm flat}(\theta)}
\ \propto\ K(h)=\frac{J_3(h)}{J_1(h)}
=E_{\rm flat}[r^2\mid\theta]
=1+h^2+\frac1{J_1(h)}.
$$

The tube law reweights the same data by the latent total squared magnitude
associated with each possible mass angle. It does not access additional
observations. The intermediate law is exactly the mixture

$$
P_{r^2+2}=\frac{Z_{\rm tube}P_{\rm tube}+2Z_{\rm flat}P_{\rm flat}}
{Z_{\rm tube}+2Z_{\rm flat}},
$$

where both normalizers use the same underlying likelihood constants.

For point readouts, two further exact integrals remove radial quadrature:
\(J_2=\sqrt{\pi/2}e^{h^2/2}\) and \(J_4=(h^2+3)J_2\).
The ratio of means is the ratio of the \(\sin\theta\)-weighted and
\(\cos\theta\)-weighted integrals of \(J_2\), \(J_4\), or \(J_4+2J_2\),
respectively. These are point readouts; the mass-ratio mean still diverges.

## 4. Same mass point, exactly ordered uncertainty

Suppose \(\|U\|^2=\|V\|^2=a\), and set \(d=U\cdot V\).
For \(\ell=\log(M/s)\),

$$
h^2=a+d\operatorname{sech}\ell.
$$

All three distributions are symmetric under \(\ell\mapsto-\ell\), and
their ratio-of-means points equal \(s\). The function \(K(h)\) is strictly
increasing for \(h>0\). One proof differentiates its radial expectation:
\(K'(h)=\operatorname{Cov}(r^2,r\coth(hr)-1/h)>0\), since both functions
increase with \(r\) under a nondegenerate radial law.

It follows that:

- **Positive observed alignment:** the tube-to-flat density ratio decreases
  with \(|\ell|\). Tube gives stochastically smaller absolute log mass and
  narrower central log intervals.
- **Negative observed alignment:** the ordering reverses. Tube gives
  stochastically larger absolute log mass and wider central log intervals.
- **Perpendicular observations:** \(h\) is constant, so all three laws are
  exactly the same. Their mass ratio is half-Cauchy, just as at zero readings.

The \(r^2+2\) law lies between the other two in these equal-norm orderings.
This does not assert a general variance ordering for unequal norms.

### Numerical example: both observed standardized lengths equal 3

Every mass point below equals \(s\). Entries are central 95% conditional
intervals for \(M/s\), not repeated-sampling confidence guarantees.

| Observed relation | Original flat | Tube | Candidate covariance |
|---|---|---|---|
| Same direction | [0.284, 3.523] | [0.319, 3.134] | [0.315, 3.171] |
| Perpendicular | [0.0393, 25.452] | [0.0393, 25.452] | [0.0393, 25.452] |
| Opposite directions | [0.00721, 138.67] | [0.00540, 185.11] | [0.00570, 175.32] |

A procedure given only the two observed magnitudes and their component error
scales receives identical inputs in all three rows. The full-vector construction
uses the dot product to distinguish their mass uncertainty. This is a concrete
example of information used beyond the magnitudes; it does not quantify an
information advantage or prove better calibration.

The perpendicular result holds at **any common observed length**, not just 3.
Indeed, any proper radial reference \(\rho(r)\) leaves a uniform
\(\theta\)-law in this configuration. The invariance concerns uncertainty
about the mass ratio; the radial information can change substantially.

Opposed observed centers remain permitted inputs under the Gaussian
uncertainty model. Their broad mass law is obtained while retaining Newton II
as a premise, not by rejecting that premise.

### Weak observed signals: the first departure from the null law

Write \(P=\|U\|^2,Q=\|V\|^2,D=U\cdot V\), and let both observed vectors
scale as \(O(\varepsilon)\). All remainders below are \(O(\varepsilon^4)\):

| Reference | \(E[\ell]\) | \(\operatorname{Var}(\ell)\) | \(\widehat m/s\) |
|---|---|---|---|
| Flat | \((P-Q)/6\) | \(\pi^2/4-\pi D/9\) | \(1+(P-Q)/6\) |
| Tube | \((P-Q)/3\) | \(\pi^2/4-2\pi D/9\) | \(1+5(P-Q)/18\) |
| Candidate covariance | \((P-Q)/4\) | \(\pi^2/4-\pi D/6\) | \(1+7(P-Q)/30\) |

The norm imbalance moves the center; the dot product changes the spread at
this order. Tube reacts twice as strongly as flat in the log-mass mean and
variance corrections. That is a sensitivity difference, not twice the data
or twice the information. All three retain \(M^{-2}\) upper density tails
and an infinite ordinary mass-ratio mean for finite observed inputs.

## 5. A useful operational test: does mere bookkeeping change the answer?

Unobserved subdivision must be distinguished from additional measurements.
If one total force is split into positive parts using a normalized conditional
split distribution, integrating those unobserved parts restores the original
law. This works with flat, tube, or angular reference choices.

Assigning fresh independent flat reference measures to the parts is a model
change. With one shared acceleration and direction,

$$
\prod_{i=1}^N df_i\quad\xrightarrow{\ f=\sum_i f_i\ }\quad
\frac{f^{N-1}}{(N-1)!}\,df.
$$

For unchanged zero-centered unit-error total-vector data, this changes the
mass readout from 1 to \(\pi/2\) to 2 as the number of unobserved force parts
changes from one to two to three. The difference is decomposition volume,
not additional measured information. The statement concerns locally finite
reference measures and their resulting proper posteriors, not normalized
independent flat population distributions.

Consistent bookkeeping is therefore a meaningful requirement, but it does
not uniquely select flat magnitudes. Nor does final readout additivity require
a fresh fit to the total to equal propagation from every possible part model.

## 6. Where this leaves the choice

The current stated premises permit several constructions. To select one,
we need an additional declared meaning for discrepancy, or a performance
criterion tied to how the mass estimate will be used.

There are now concrete comparisons available:

- **Structural interpretation:** additive balance discrepancy versus local
  Gaussian distinguishability, with the tolerance rule specified.
- **Representation consistency:** preserve the same probability law through
  unit changes, transformations, and unobserved bookkeeping refinements.
- **Mass performance:** compare point error and interval coverage at fixed
  true masses and excitation strengths, using the same simulated data.
  Report interval width alongside coverage, rather than rewarding width alone.

The earlier [calibration study](weighting_and_calibration.md) already shows
that credible interval width and fixed-truth coverage can disagree. Its weak
excitation impossibility result still applies: tuning an always-finite interval
cannot create uniform 95% coverage over unrestricted masses and arbitrarily
weak excitation. The present exact results explain part of the finite-data
behavior without selecting a replacement measure.

## Verification and provenance

[operational_checks.py](operational_checks.py) is a separate one-dimensional
implementation. It checks mass points and CDFs against the existing
two-dimensional integrators on 12 fixed and seeded random observation pairs,
including zero, oblique, opposed, and strong aligned cases. Largest checked
differences were \(3.30\times10^{-13}\) relative for points and
\(2.46\times10^{-14}\) absolute for CDFs. These are checked-case differences,
not uniform numerical error bounds.

It also checks channel reciprocity, the intermediate measure, equal-norm
intervals, log-moment quadrature refinement, the exact perpendicular law, and
the weak-data expansion at three scales. The second-order expansion remainder
decreases by approximately a factor of 16 when the observed vectors are halved.
Quantiles are bracketed within log-ratio [-36, 36]; log moments use the same
finite range, with an extension to [-44, 44] checked on the equal-norm cases.
Fixed quadrature orders and these brackets are not guarantees for arbitrary
signal sizes. The candidate-law CDF is checked independently using the
normalizer-weighted mixture of the earlier flat and tube integrations.
The baseline estimator and its existing checks are unchanged.

```powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/operational_checks.py' --output '.tools/mass_equation_20260915/operational_warrant/results.json'
```

Independent derivation and scope reviews, and machine results, are in
`.tools/mass_equation_20260915/operational_warrant/`. The exact interval ordering
is a theorem; the table illustrates it. No new sampling-coverage claim is made
for the candidate-covariance rule.

Managed context remains the coherent VD-docs snapshot
[`8782a293297330f7a354d87ebe62792018fa8866`](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866).
Selection follows the [catalog](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/catalog.yaml),
the [working specification, especially Section 9](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md),
the [readout audit's composition and measure discussion](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md),
and the [law-assumed-correct steering](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md).
Active local scope is the preceding six investigation contributions and their
relevant integrators. This is a targeted mathematical investigation, not an
exhaustive context or literature review. No managed document was edited.
