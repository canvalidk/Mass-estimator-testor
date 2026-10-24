# Operational warrant for a compatibility neighborhood

Independent bounded scratch audit, 2026-09-15. No production or managed edits.

## Conclusions

1. The inspected project commitments do not select fixed cross-product-area
   tolerance over fixed distance in noise units. They require a Gaussian
   measurement account, preservation of the original errors, symmetric use of
   both channels, and a declared latent-pair measure. Flat magnitudes are the
   current declared baseline. The extra mismatch axioms in the new note are
   mathematically effective but need their own operational warrant.
2. An unobserved partition that is genuinely only a change of representation
   should not change total-mass inference. That requirement can be satisfied
   for any of the candidate total laws by adding a normalized conditional
   distribution of the split. It does not select area, tube, or flat weighting.
3. Assigning fresh independent flat measures to each new constituent is a
   substantive change of reference, not an innocuous refinement. An explicit
   null-data calculation below changes the readout from 1 to pi/2 to 2 as the
   unobserved number of force parts changes from one to two to three.
4. Consistent independent constituent reference families exist, including a
   gamma convolution family. They require an additive allocation parameter
   across partitions. Neither that parameter nor independence of unobserved
   parts is currently supplied. This construction shows possibility, not a
   reason to adopt a population model.

## 1. What the operational sources actually support

The full working specification, Section 9, explicitly distinguishes the
Gaussian factor exp(-Q/2) from the declared latent measure rho. Its current
implementation uses uniform direction and flat positive magnitudes, while
stating that this is not uniquely dictated by Newton II. Section 33 retains
the measure among policy choices. Sections 29-30 require preservation of
uncertainty, discrepancy, identification standing, and provenance.

The later law-assumed-correct steering governs the present investigation over
the old rejection-gate passages. It requires explaining even opposing observed
vectors through measurement error, measuring corrections relative to their
declared uncertainties, and treating the channels symmetrically as information.
It expressly warns that the current integrated ratio-of-means readout is not
generally the least-error reconciliation. It does not authorize replacing it
with an optimizer.

These requirements support the same Gaussian likelihood factor under every
candidate reference. They do not say that equal coordinate volumes, equal
angular volumes, equal cross-product residual volumes, or equal noise-distance
volumes must have equal prior support. In particular, the words 'smallest
error' cannot supply a reference-measure theorem after the source itself has
distinguished likelihood optimization from the agreed equation.

The readout-composition source explicitly treats a lump's law as a pushforward
of the parts' joint law and distinguishes it from a fresh fit with a new flat
measure. It also discusses noise-metric cone volume as another plausible
policy. Its favorable language about that alternative is not an established
optimality result; it does show that fixed area is not an inherited exclusive
operational commitment. The rigorous uniqueness theorem accepts a supplied
probability law and does not select its construction.

## 2. What must be added to favor one neighborhood

### Fixed area

The current note establishes a sufficient bridge: require a continuous,
rotation-covariant oriented discrepancy, force-side additivity off the
compatible set, symmetry of its norm on channel exchange, and a fixed
Euclidean tolerance in that norm. Then its magnitude is proportional to
||X cross Y||, and the area-conditioning limit has flat magnitudes.

The missing operational premises are specific:

- Why must an oriented discrepancy on incompatible inputs be recoverable by
  adding discrepancies of arbitrary force components? Physical mass addition
  on exactly compatible states sees only zero discrepancies.
- Why is this discrepancy an operationally conserved quantity rather than a
  useful algebraic representation of incompatibility?
- Why must the residual tolerance be independent of the candidate signal
  magnitudes and direction? State-dependent widths reinstate measure freedom.
- Why is norm exchange, a property of the discrepancy, required beyond
  reciprocal relabeling of the final inference?

An explicit requirement to transmit and sum residual vectors without retaining
the original vectors could motivate the first premise. The current result type
instead preserves the original measurements. Recomputing a nonlinear residual
from their sum therefore violates no stated information-storage requirement.

### Noise distance

Equal normal displacements in standardized measurement coordinates have equal
local Gaussian significance. Requiring equal transverse neighborhood widths
in those coordinates leads to metric-tube weighting. This is a clear possible
operational criterion, but it is also an extra criterion for the reference
volume. Using Mahalanobis likelihoods does not, by itself, require it.

The fixed-area rule and this local equal-significance criterion make different
predictions. At a regular compatible state with standardized combined magnitude
r, a small orthonormal normal displacement n has area r||n|| to first order.
A fixed area threshold therefore permits normal distance approximately epsilon/r.
It becomes tighter in noise units as r increases. A fixed noise-distance tube
permits approximately epsilon at every r.

This difference is not a reason to reject area conditioning without first
adopting an equal-significance-volume requirement. Conversely, physical
composition alone cannot override it and establish that area is preferable.
The exact sampling identity E[A_observed^2]=(d-1)(r_true^2+d), under independent
unit Gaussian channel errors at compatible means, confirms that a fixed area
threshold is not a fixed sampling-significance threshold. Normalizing a
discrepancy by a candidate-dependent sampling covariance would be yet another
explicit rule; it is not already the fixed input Gaussian covariance.

## 3. Bookkeeping refinement: the precise invariance that is available

Let z=(s,alpha,u) denote the total compatible pair, and let Pi(dz) be any
chosen reference measure. Introduce positive fractions t=(t_1,...,t_N) on the
unit simplex and a conditional probability kernel K_z(dt), so that

    integral K_z(dt) = 1,
    f_i = s t_i.

Define the augmented reference as Pi(dz) K_z(dt). If the data likelihood
L(data|z) depends only on the total, integrating out t returns exactly
L(data|z) Pi(dz). Thus the posterior of the total, the induced mass law, and
the ratio-of-means point remain unchanged. This works for flat, tube, angular,
or another normalizable total reference.

The same statement holds before taking a conditioning limit. Augment a
specified ambient probability law by a normalized kernel, and define the
compatibility event on the same total vectors. Marginalizing the split commutes
with conditioning for every positive-probability tolerance. Any existing limit
is consequently unchanged.

For example, a uniform conditional distribution on the simplex has coordinate
density (N-1)!. If the total-force reference is flat ds, its induced density
in constituent coordinates is

    (N-1)! / (f_1+...+f_N)^(N-1) times product_i df_i.

It is not a product of flat constituent measures. The normalization is exactly
what prevents an unobserved number of decompositions from supplying extra
weight to large totals.

This representation invariance is appropriate when a split is explicitly only
bookkeeping. Adding genuine separate measurements or imposing new physical
assumptions on constituents can change the inference. Requiring each part's
vector to satisfy a separately imposed finite tolerance is also a change of
conditioning procedure unless its relationship to the total event is specified.
Forces of real parts and independent measurements of them are not automatically
equivalent to adding an unobserved fraction variable.

## 4. Fresh flat constituents fail that bookkeeping test

Keep one shared acceleration and one shared direction. If the reference is

    (product_i df_i) d(alpha) dOmega,

then s=sum_i f_i and simplex fractions have Jacobian s^(N-1). Integrating the
simplex gives the total reference

    s^(N-1)/(N-1)! ds d(alpha) dOmega.                         (1)

This is an identity of locally finite improper measures, not an assertion
that normalized posteriors are uniform or that independent constituent priors
are physically appropriate. For a likelihood depending only on the total,
(1) is exactly the extra weighting. Part-specific likelihoods would also need
integration over the split.

A transparent counterexample uses the same zero-centered independent isotropic
unit Gaussian cloud for the total force and acceleration for every N. On the
compatible set its likelihood is proportional to exp[-(s^2+alpha^2)/2]. Then

    E[s] = sqrt(2) Gamma((N+1)/2) / Gamma(N/2),
    E[alpha] = sqrt(2/pi),
    R_N = sqrt(pi) Gamma((N+1)/2) / Gamma(N/2).

| Unobserved force parts N | Readout in standardized mass units |
| ---: | ---: |
| 1 | 1 |
| 2 | pi/2 = 1.5707963268 |
| 3 | 2 |
| 4 | 3pi/4 = 2.3561944902 |

The data and physical total are identical. Only a newly assigned reference on
unobserved decompositions changed. In physical units multiply the table by
sigma_F/sigma_a. These are measured zero centers with finite noise, not an
absence of measurements and not evidence for any identified true mass.

There is no failure of the readout's C1 axiom: for each resulting joint law,
linearity of expectation still gives R(sum f_i,alpha)=sum R(f_i,alpha).
The failure is equating freshly constructed laws across different reference
assignments. The argument is therefore directly relevant to a proposed
stronger composition principle, but is not a refutation of the existing theorem.

## 5. A consistent convolution family exists, without choosing a population model

For a>0, define the positive-magnitude reference

    nu_a(df) = f^(a-1) df / Gamma(a).

The beta integral proves

    nu_a * nu_b = nu_(a+b),

because integrating x^(a-1)(s-x)^(b-1) over 0<x<s gives
s^(a+b-1) Gamma(a)Gamma(b)/Gamma(a+b). These measures are locally finite but
improper on the whole positive half-line.

Fix a total shape A>0. Any finite partition with positive shapes a_i satisfying
sum_i a_i=A can use product_i nu_(a_i); its total is nu_A independently of the
partition. Conditional fractions have the normalized Dirichlet density with
those shapes. Refining one cell requires splitting its shape additively.

For A=1 the total reference is exactly flat ds. An equal N-way refinement uses
a_i=1/N, not a_i=1. The singularity of each f_i^(1/N-1) at zero is integrable.
The naive flat-parts construction instead gives every new part shape 1, so it
changes the total shape to N and reproduces (1).

There is no selection of A=1 in this convolution identity. Every A>0 is
consistent. To make this a physical partition rule would additionally require
an additive allocation of shape to pieces, an interpretation of that allocation,
and a reason for conditional independence of their force references. None is
provided by mass addition itself. If both channels use the same positive total
shape, reciprocal/channel symmetry still leaves a continuum of choices.

A proper common-rate gamma family supplies the analogous convolution identity
for probability laws, but adds a rate scale and a population interpretation.
No such distribution is adopted here. Only finite-partition reference identities
are needed; this audit does not assert an infinite improper stochastic process.

## 6. A falsifiable comparison criterion that preserves the original task

Two questions should be tested separately.

**Representation test.** Hold the supplied total-vector Gaussian law and the
compatibility event fixed. Add two and three unobserved parts by normalized
conditional split kernels. Require all total-mass probabilities and readouts
to remain unchanged, up to numerical tolerance. This is exact and falsifiable;
it catches resetting the reference after refinement. Every candidate can pass
it, so it cannot establish area's superiority.

**Statistical performance test.** Compare area/flat and noise-distance/tube
under identical fixed-truth, single-pair Gaussian measurement experiments.
Use a predeclared grid of positive true masses and standardized signal
strengths, including weak, moderate, and strong excitation. Do not generate
latent excitation vectors from a new population distribution. Each grid point
is a separate fixed physical state whose measurement error is repeatedly
simulated.

For each point, report the following for both policies:

- coverage of the induced-ratio 95% interval for the fixed true mass;
- mean interval length in log-mass units;
- reciprocal-symmetric point loss E[(log(m_hat/m_true))^2];
- the existing discrepancy/identification diagnostics.

Keep interval construction, covariance, and the ratio-of-means readout fixed
so the reference choice is the only intervention. A concrete uncertainty
criterion could require coverage within a predeclared tolerance throughout
the operating grid before rewarding shorter intervals. Monte Carlo error must
be smaller than that tolerance. Neither shorter intervals alone nor unchanged
central readouts establish that more information has been used.

This would empirically distinguish uncertainty behavior under the declared
measurement experiment. It would not prove universal optimality or select a
unique reference without an operating region and loss criterion. At true zero
excitation mass is not sampling-identified, so that boundary must be reported
as a separate diagnostic, not used to manufacture a uniquely correct finite
mass interval. No rejection of Newton II is part of this comparison.

## Sources and scope

Managed files were inspected at the single pinned snapshot
`canvalidk/VD-docs@8782a293297330f7a354d87ebe62792018fa8866`:

- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`,
  especially Sections 9, 29-30, and 32-35, from the connector cache.
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md`,
  reread through the connector for this audit. Its historical custody wording
  is not a statement that the present managed source is unsubmitted.
- `Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md`,
  especially Sections 4-5, from the connector cache.
- `Newton-analysis/_dscn_mass_estimator/mass_readout_composition_uniqueness_2026-09-06.md`,
  the previously inspected rigorous conditional-on-law theorem.

Local sources were the current `investigation/why_flat_conditioning.md` and the
earlier scratch `composition_bridge/physical_warrant_audit.md`. This is a
bounded operational and mathematical audit, not an exhaustive corpus or
literature search. The gamma identity and partition counterexample are derived
above; their four elementary numerical readouts were checked with Python.
