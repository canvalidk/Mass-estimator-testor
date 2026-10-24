# What the Gaussian vector assumptions determine

**Second contribution to the mass-equation investigation — 15 September 2026.**
This follows the user's question about Gaussian, unskewed, spherical vector
uncertainty. It develops the construction of the probability law supplied to
the readout. Newton II remains a premise.

## Result

**Even a completely specified joint Gaussian distribution over the two
vectors does not, by itself, select the current distribution over exactly
compatible pairs.** In dimensions greater than one, exact positive
collinearity has probability zero. Two explicit ways of approaching that
constraint give different limiting distributions.

This freedom is a choice in the inference construction. Gaussianity does
not force a particular adjustable parameter, and it does not license choosing
whichever construction produces the narrowest interval.

At measured zero-zero, in three dimensions, both constructions return the
same mass point, yet one has a finite mass variance and the other has a
divergent mass mean. Away from zero-zero, their mass points also differ.

## 1. Three things to specify

1. **The joint measurement model.** Gaussian, spherical uncertainty of each
   vector does not specify cross-channel dependence. Under a jointly Gaussian,
   simultaneously rotation-invariant error model, the covariance can be

   $$
   \Sigma=\begin{pmatrix}
   \sigma_F^2 I&c\sigma_F\sigma_a I\\
   c\sigma_F\sigma_a I&\sigma_a^2 I
   \end{pmatrix},\qquad -1<c<1.
   $$

   Independence fixes \(c=0\), as in the current isotropic laboratory case.
   Separate provenance alone does not imply statistical independence.
   Gaussian marginals alone need not even imply a jointly Gaussian law.

2. **What the Gaussian describes.** A Gaussian likelihood
   \(p(\widetilde{\mathbf F},\widetilde{\mathbf a}\mid\mathbf F,\mathbf a)\)
   is not automatically a probability distribution over the unknown vectors.
   A posterior needs an inference construction, such as a prior. Alternatively,
   the vector uncertainty distributions may themselves be explicitly supplied
   as the inputs. The analysis below grants that stronger input: a fully
   specified ambient joint Gaussian probability law.

3. **How compatibility is imposed.** A density gives probabilities only
   together with its reference volume. Evaluating the Gaussian on
   \((f\mathbf u,\alpha\mathbf u)\) does not specify that volume.

The working specification already separates these commitments: its §9 lists
Gaussian uncertainty, uniform common direction, and flat positive magnitudes
as the current policy, and explicitly identifies the latent measure as an
additional input. The present result explains why this separation matters.

## 2. Two explicit constructions from the same Gaussian

Let \(g(\mathbf F,\mathbf a)\) be a nonsingular Gaussian density with respect
to Cartesian volume. The compatible set is

$$
\mathcal C=\{(f\mathbf u,\alpha\mathbf u):f,\alpha>0,
\ \mathbf u\in S^{d-1}\}.
$$

It has dimension \(d+1\) inside a \(2d\)-dimensional ambient space. For
\(d>1\), ordinary division by \(P(\mathcal C)\) cannot define conditioning.

### A. Shrink the allowed angular mismatch

Condition on \(\angle(\mathbf F,\mathbf a)<\epsilon\), then let
\(\epsilon\downarrow0\). Cartesian polar volume is

$$
d\mathbf F\,d\mathbf a=
f^{d-1}\alpha^{d-1}\,df\,d\alpha\,d\Omega(\mathbf u)\,d\Omega(\mathbf v).
$$

The small cap around \(\mathbf v=\mathbf u\) has the same area for every
\(\mathbf u\). Dividing numerator and denominator of the conditional
expectation by that cap area gives

$$
dP_{\rm angle}\ \propto
g(f\mathbf u,\alpha\mathbf u)(f\alpha)^{d-1}
\,df\,d\alpha\,d\Omega. \tag{1}
$$

A Gaussian envelope times the displayed polynomial is integrable, so
dominated convergence justifies this limit, including small magnitudes.

### B. Shrink the distance to the compatible set

Choose the dimensionless Cartesian coordinates
\(\mathbf x=\mathbf F/\sigma_F\),
\(\mathbf y=\mathbf a/\sigma_a\). Condition on their Euclidean distance to
the compatible set being less than \(\epsilon\), then let
\(\epsilon\downarrow0\). This is a tube around the set, using a specified
metric; force and acceleration cannot be combined in an unscaled physical
Euclidean distance.

Put \(x=f/\sigma_F\), \(y=\alpha/\sigma_a\), \(r=\sqrt{x^2+y^2}\).
For the embedding \(\Phi(x,y,\mathbf u)=(x\mathbf u,y\mathbf u)\),

$$
\|d\Phi\|^2=dx^2+dy^2+r^2\|d\mathbf u\|^2.
$$

Consequently its surface volume is \(r^{d-1}dx\,dy\,d\Omega\), and the
thin-tube limit is

$$
dP_{\rm tube}\ \propto
g(f\mathbf u,\alpha\mathbf u)
\left(\frac{f^2}{\sigma_F^2}+\frac{\alpha^2}{\sigma_a^2}\right)^{(d-1)/2}
\,df\,d\alpha\,d\Omega. \tag{2}
$$

Appendix A checks the axes and singular origin, which a purely local surface
argument would omit. Constants independent of the latent state cancel in
both expressions.

These are different conditioning events approaching the same set. They are
not a disagreement produced by expressing one probability law in different
coordinates. Both respect simultaneous rotations and interchange of the
standardized channels. Those symmetries do not select between them.

The construction issue belongs to the established mathematics of singular
conditioning; [Bungert and Wacker, §§3–4](https://arxiv.org/html/2009.04778v2#S3)
distinguish a geometric surface measure from conditioning through a specified
random variable. A specified metric supplies a geometric convention; it does
not follow solely from a probability-zero set. Equations (1)–(2) above are
derived here for this estimator.

In \(d=1\), compatibility is the positive-probability event of equal signs;
both weighting powers vanish and these constructions agree. The comparison
here concerns an unknown shared direction. A direction fixed by the experiment
requires its own specified construction.

## 3. Exact zero-reading comparison in three dimensions

Take independent isotropic Gaussian inputs centered at zero. Define
\(s=\sigma_F/\sigma_a\), \(R=M/s=x/y\), and
\(x=r\sin\theta,\ y=r\cos\theta\). The Gaussian factor is
\(e^{-r^2/2}\).

| Compatible-pair measure | Radial density, up to normalization | Ratio-angle density |
|---|---|---|
| Current flat magnitudes | \(r e^{-r^2/2}\) | Uniform |
| Standardized distance tube | \(r^3 e^{-r^2/2}\) | Uniform |
| Angular tolerance | \(r^5 e^{-r^2/2}\) | \((4/\pi)\sin^2(2\theta)\) |

Under angular conditioning the two standardized magnitudes are independent
\(\chi_3\) variables. Under tube conditioning they are coupled, even though
the original ambient vectors were independent. All three have exchange
symmetry, so **their ratio-of-means points and mass medians equal \(s\)**.

| Quantity | Flat magnitudes | Distance tube | Angular tolerance |
|---|---:|---:|---:|
| \(\widehat m/s\) | 1 | 1 | 1 |
| Central 95% interval for \(M/s\) | [0.03929, 25.45170] | [0.03929, 25.45170] | [0.25450, 3.92927] |
| SD of \(\log(M/s)\) | 1.57080 | 1.57080 | 0.68367 |
| \(E[f]/\sigma_F=E[\alpha]/\sigma_a\) | 0.79788 | 1.19683 | 1.59577 |
| \(E[M]/s\) | Diverges | Diverges | \(4/\pi\) |

These are conditional probability intervals under the stated constructions,
not yet calibrated repeated-experiment confidence intervals.

The exact ratio densities are

$$
p_{R,\rm flat}(z)=p_{R,\rm tube}(z)=\frac{2}{\pi(1+z^2)},
\qquad
p_{R,\rm angle}(z)=\frac{16}{\pi}\frac{z^2}{(1+z^2)^3},\quad z>0.
$$

For the angular law,
\(\operatorname{Var}(\log R)=\pi^2/4-2\),
\(E[R^2]=3\), and \(\operatorname{Var}(R)=3-16/\pi^2\).
The other two laws have divergent first moments; ordinary variance about
their mean is therefore undefined. Their finite log spread remains useful.

Flat magnitudes and the tube agree on the **entire** null mass distribution,
yet their different mean magnitudes show that the full latent laws differ.
Agreement on a mass distribution alone therefore cannot establish equivalence
for reuse with the original force or acceleration.

**The tail distinction persists at every finite Gaussian input.** Angular
conditioning gives density orders \(m^{d-1}\) near zero and \(m^{-d-1}\)
at infinity, so \(E[M^t]\) exists exactly when \(-d<t<d\). Flat magnitudes
and distance tubes retain a positive density at zero and an \(m^{-2}\)
high-mass tail. For a tube, its radial weight is nonzero at either nonzero
axis; it does not suppress the small-denominator configurations responsible
for that tail. Thus the difference in existence of mass variance is structural.

## 4. The choice also changes nonzero mass points

The following use \(d=3\), independent channels, unknown shared direction,
and standardized measured vectors \(\widetilde{\mathbf F}/\sigma_F\),
\(\widetilde{\mathbf a}/\sigma_a\). Each column uses the **same** ambient
Gaussian density and the **same** readout \(E[f]/E[\alpha]\).

| Standardized force | Standardized acceleration | Flat \(\widehat m/s\) | Tube \(\widehat m/s\) | Angular \(\widehat m/s\) |
|---|---|---:|---:|---:|
| (3,0,0) | (0,0,0) | 3.17913 | 3.64008 | 2.08877 |
| (4,0,0) | (1,1,0) | 2.85627 | 2.95175 | 1.97118 |
| (4,0,0) | (−1,0,0) | 6.76502 | 7.50604 | 3.42463 |
| (8,0,0) | (4,0,0) | 1.99991 | 1.99993 | 1.84134 |

The comparison does not rank the constructions. It establishes that choosing
how to form the compatible-pair distribution affects both the central mass
and its uncertainty, including for displayed aligned readings.

## 5. Consequence for the investigation

The original direction-weighting freedom concerns the final readout applied
to an already supplied law. The freedom established here is upstream:
**which law is supplied?** These are distinct questions.

The current flat-magnitude measure remains a mathematically valid declared
policy. It has not been derived from spherical Gaussian vector uncertainty
alone. Neither of the two conditioning constructions above automatically
replaces it.

The next substantive question is operational: what measurement or inference
procedure is the compatible-pair law meant to represent? A choice justified
there can constrain its measure before empirical tuning begins. Assessment
of information use and preferable estimators then needs a specified sampling
experiment, uncertainty calibration, and performance criterion. A narrower
conditional interval alone cannot establish that more information was used.

## Appendix A. Boundary control for the tube limit

At an interior point \((x\mathbf u,y\mathbf u)\), let \(r=\sqrt{x^2+y^2}\)
and \(\mathbf w\perp\mathbf u\). Normal coordinates are

$$
\mathbf X=x\mathbf u+(y/r)\mathbf w,\qquad
\mathbf Y=y\mathbf u-(x/r)\mathbf w.
$$

For \(\|\mathbf w\|<r\), the closest compatible pair is
\((x\mathbf u,y\mathbf u)\), and the distance is \(\|\mathbf w\|\).
The two singular values of the matrix with columns \(\mathbf X,\mathbf Y\)
are \(r\) and \(\|\mathbf w\|\). Its leading component has positive
channel coefficients. Conversely, every pair with positive dot product has
this representation, apart from irrelevant degeneracy sets.

The volume Jacobian relative to \(dx\,dy\,d\Omega\,d\mathbf w\) is

$$
J=r^{d-1}\left(1-\frac{\|\mathbf w\|^2}{r^2}\right).
$$

One determinant check uses \(x=r\cos\eta,y=r\sin\eta\): after an orthogonal
channel rotation, the nontrivial block is
\(\left(\begin{smallmatrix}rI&\mathbf w\\-\mathbf w^T&-r\end{smallmatrix}\right)\).
Its absolute determinant is \(r^d(1-\|\mathbf w\|^2/r^2)\); divide by
\(dx\,dy=r\,dr\,d\eta\).

Substitute \(\mathbf w=\epsilon\mathbf v\) in the tube integral and divide
by the normal-ball volume \(\kappa_{d-1}\epsilon^{d-1}\). The integrand tends
to \(g\) on the compatible set times \(r^{d-1}\). The bound
\(0\le J\le r^{d-1}\), together with Gaussian decay, permits dominated
convergence even near \(r=0\).

If \(\mathbf X\cdot\mathbf Y\le0\), the closest point in the closed
positive cone lies on an axis, and the distance is
\(\min(\|\mathbf X\|,\|\mathbf Y\|)\). This part of the tube has Gaussian
probability \(O(\epsilon^d)\), negligible relative to
\(\epsilon^{d-1}\). Therefore neither the axes nor the origin acquire a
limiting atom. Distance to the open cone equals distance to its closure.

## Appendix B. Correlation is a separate, quantifiable input

Retain the current flat magnitude measure, measured zero-zero, and the
jointly Gaussian covariance in §1. This changes the joint measurement model;
it is separate from comparing conditioning geometries at fixed independence.
For \(R=M/s\), let \(\beta=\pi/2+\arcsin c\). Integrating the standardized
radial coordinate gives

$$
p_R(z)=\frac{\sqrt{1-c^2}}{\beta}\frac1{z^2-2cz+1},\qquad z>0.
$$

Channel symmetry preserves the point \(\widehat m/s=1\), and reciprocal
symmetry preserves the median 1. Nevertheless,

$$
E[R^t]=\frac{\pi\sin(\beta t)}{\beta\sin(\pi t)},\quad |t|<1,
\qquad
\operatorname{Var}(\log R)=\frac{\pi^2-\beta^2}{3}.
$$

The Mellin transform at \(t=0\) takes its continuous value 1. It follows by
integrating the displayed ratio density; differentiating its logarithm twice
at zero yields the log variance. At \(c=0\) this recovers the half-Cauchy law.
As \(c\to1\) it concentrates at 1; as \(c\to-1\) its density tends to
\(1/(1+z)^2\). These endpoint limits do not extend the nonsingular assumption
to \(c=\pm1\). The correlation is to be justified from the measurement
process, not selected to obtain a desired uncertainty width.

## Reproduction and provenance

[gaussian_conditioning.py](gaussian_conditioning.py) reproduces both numerical tables.
It integrates common direction analytically in three dimensions, then uses
numerical quadrature over the two magnitudes. It checks polar against Cartesian
integration, doubles order from 128 to 256, extends the radial cutoff, checks
reciprocity and exact null moments, and compares flat-policy points with the
existing estimator. The largest relative difference among quadrature variants
was \(2.56\times10^{-14}\) on the tabulated cases. All assertions passed.
Separate log-density quadrature checks Appendix B at \(c=-0.8,0,0.8\).
These numerical comparisons do not replace the tail and conditioning proofs.

~~~powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/gaussian_conditioning.py'
~~~

The geometric derivations and null uncertainty received independent parallel
reviews. An additional finite-tube check remains exploratory scratch work.
This contribution uses Python and NumPy; it does not change the baseline
estimator or select an estimator policy.

Managed main was re-resolved to the same snapshot,
[8782a293297330f7a354d87ebe62792018fa8866](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866).
The catalog was consulted; this pass reread the relevant sections of the
[working specification](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md),
[original vector-division analysis](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/vector_division_inertial_mass_measurement_d1.md),
and [chain audit](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_chain_to_the_estimator.md).
The [main investigation](README.md) documents the broader prior-pass sources
and governing user steering. Local scope was this investigation and the active
laboratory integrator. Selection targeted Gaussian, covariance, latent measure,
and null behavior. The one external probability paper supports the general
conditioning distinction; this is not a literature or priority review.
