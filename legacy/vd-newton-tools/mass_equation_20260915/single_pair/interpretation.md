# Single-pair Gaussian measurement law and the compatible-pair measure

Scratch interpretation audit, 2026-09-15. Original single-pair scope only.
This is not a managed record, new estimator policy, or production edit.

## 1. What is already declared

The working specification's construction distinguishes a Gaussian measurement
likelihood from the latent-pair measure. Its current reference is

    dPi_flat = df d(alpha) dOmega(u),
    f>0, alpha>0, ||u||=1,

and the normalized compatible-pair law is

    dP_y proportional to L(y | f u, alpha u) dPi_flat.

The source already declares this measure; the investigation is asking what
interpretation or further requirement could justify it. Discovering another
way to derive this same law does not establish that the earlier specification
contained a contradiction or that its declared construction was undefined.

## 2. Two distinct probabilistic readings can give the same posterior formula

### A. A reference prior directly on compatible latent pairs

Specify a shared uniformly distributed direction and a flat reference measure
over the two positive magnitudes, then apply the Gaussian observation model.
This is a coherent generalized-Bayesian posterior construction.

The two flat half-line magnitude measures are improper. They are not
normalized independent uniform random variables, and therefore do not by
themselves define a complete prior-predictive generative distribution.

A proper approximation makes that distinction concrete:

    f ~ Uniform(0,F_max),
    alpha ~ Uniform(0,A_max),
    u ~ uniform sphere,
    y | f,alpha,u ~ the specified Gaussian.

For finite cutoffs this is a proper generative model. As both cutoffs tend
to infinity, the posterior converges to the flat-magnitude posterior: the
cutoff indicators increase to one and the Gaussian-compatible likelihood
is integrable. The same limit holds along non-monotone cutoff sequences
that both diverge, by dominated convergence. This is a construction of the
posterior limit, not an assertion that real experiments generate uniformly
distributed latent magnitudes.

The improper limit does not supply a normalized prior-predictive law or an
absolute evidence normalization. The ordinary posterior calculations can
nevertheless be well defined.

### B. Conditioning a supplied ambient Gaussian law

Alternatively, begin with an explicitly supplied probability density
g(F,a) with respect to full Cartesian volume, and impose compatibility by
a specified limiting conditioning event. A Gaussian likelihood alone is
not automatically this probability law over unknown vectors.

For dimension d>1, exact positive collinearity has ambient probability zero.
There is no unique ordinary event-conditional probability until the
construction is specified. Angular, Cartesian-tube, and wedge tolerances
can approach the same set and give different limiting measures.

For the proposed wedge construction, use standardized vectors X,Y and

    X dot Y > 0,   ||X wedge Y|| < epsilon.

At fixed positive magnitudes x,y, the wedge norm is x y sin(theta).
The small angular cap has volume proportional to

    (epsilon/(x y))^(d-1).

This cancels the Cartesian polar factor (x y)^(d-1), leaving the flat
positive-magnitude reference after normalization. The other agent's proof
handles the small-product region and the limit rigorously. The key
interpretation is that this is an **absolute wedge-defect tolerance**.
Its allowed angular mismatch grows when either magnitude becomes small.

Thus this precisely defined conditioning rule can lead to exactly the
existing flat-magnitude compatible law. It provides one mathematical
interpretation of that law, rather than proving it follows from Gaussianity
alone. The wedge construction fixes the additional choice.

The compatible posterior in A and the limiting conditional distribution in
B can coincide, although their preconditioning models and their statements
about generative sampling are different.

## 3. Why Cartesian-flat and magnitude-flat are different

For two unconstrained vectors, Cartesian volume is

    dF da = f^(d-1) alpha^(d-1)
            df d(alpha) dOmega(u) dOmega(v).

It has two directions and radial weights. It is not the same measure as
flat positive magnitudes with one already-shared direction. Merely setting
u=v in the written density does not establish a conditional measure.

In particular:

* Angular tolerance retains the radial factor (f alpha)^(d-1).
* An absolute wedge-defect tolerance cancels it in the limiting construction.
* A specified metric tube produces the corresponding surface-volume factor.

Calling a prior "Cartesian flat on the compatible set" is ambiguous: the
set has zero ambient Cartesian volume for d>1. A surface measure or a
limiting construction still has to be identified.

The ambiguity is not repaired by changing coordinates while forgetting the
Jacobian. A genuine coordinate change leaves an already specified measure
unchanged.

## 4. Which axioms could choose flat magnitudes?

### Symmetries alone do not choose it

In d>=2, rotation invariance with no preferred direction selects uniform
conditional direction measure, but it leaves a magnitude density
rho(f,alpha) unspecified. In d=1 the corresponding claim needs sign-exchange
symmetry, since proper rotations alone do not exchange the two directions.
Channel exchange, unit covariance, and consistent reciprocal
reparameterization do not remove that freedom.

For example,

    rho_k(f,alpha)=(f alpha)^k,  k>-1,

gives an entire family. Reference-unit constants cancel in the normalized
posterior. Every member is symmetric in the magnitude channels, compatible
with their independent unit rescalings, and rotation invariant. The Gaussian
posterior is proper: k>-1 gives integrability at either axis, and Gaussian
decay controls large magnitudes. The flat choice is k=0, but those symmetry
requirements do not single it out.

Even a stronger relative-invariance rule under independent magnitude
rescalings generally admits power laws. Channel exchange equates their
exponents and still leaves the exponent free. Ordinary unit covariance is
weaker still: every correctly transformed measure is covariant.

Literal invariance with no allowed overall rescaling is a different axiom.
The scale-Haar candidate df/f d(alpha)/alpha is singular at the axes and
does not give a proper Gaussian posterior without additional boundary
treatment. It should not be conflated with the harmless unit covariance of
the current construction.

### An additional local translation axiom would choose it

Within the usual class of locally integrable reference densities, impose
independent local translation invariance:

    rho(f+h,alpha+k)=rho(f,alpha)

whenever the translated locations remain inside the positive quadrant.
Together with rotational invariance of the direction, this requires rho to
be constant almost everywhere. It thus selects df d(alpha) dOmega up to
an irrelevant global factor. Equivalently, equal absolute increments of
each magnitude receive equal reference weight everywhere in the interior.

This is a sufficient **additional modeling axiom**. It is not implied by
Newton II, a Gaussian observation kernel, or spherical error contours.
Positive magnitudes possess a distinguished boundary at zero, and the
experiment has not automatically supplied a translation symmetry on their
latent reference weights. A local translation rule can be proposed and
examined; it must not be smuggled into the word "flat" or "uninformative."

No proper uniform distribution exists on the unbounded positive half-line.
An appeal to maximum entropy without bounds or other constraints therefore
does not select a normalized version of this prior either.

## 5. Mass-coordinate descriptions require care

With m=f/alpha,

    df d(alpha) = alpha dm d(alpha).

With b=sqrt(f alpha) and z=log(m),

    df d(alpha) = b db dz.

The same two-dimensional measure can therefore look flat in m when the
nuisance coordinate is alpha, or flat in log mass when the nuisance
coordinate is b. These are not contradictory marginal mass priors: the
nuisance coordinate changed too, and its unbounded reference volume cannot
be integrated into a normalized mass marginal. This is another reason to
state the joint reference rather than merely describe it as flat in mass.

## 6. Operational questions for the investigation

These are questions to resolve through analysis, not requests for more user
data before work can continue:

1. Is the Gaussian supplied as an observation likelihood, or as a full
   uncertainty distribution over unconstrained latent vectors?
2. If compatibility is imposed by conditioning, which measurable defect is
   being made small: angle, absolute wedge defect, or distance in a declared
   measurement metric? Equivalent zero sets need not give equivalent limits.
3. If the measure is specified directly, what operation warrants equal
   absolute magnitude increments receiving equal reference weight?
4. Does the wedge construction represent an intended experimental or
   inferential procedure, or serve as one exact mathematical representation
   of an already selected posterior?
5. Which single-pair behavior should discriminate the choices: changes of
   measurement precision, weak-excitation uncertainty, nonalignment, or
   subsequent reuse? Narrower uncertainty alone is not a selection rule.

The next production contribution should keep these distinctions attached to
the original single-pair construction. No further repeated-trial or random
excitation model is needed to explain this ambiguity.

## Source scope

The construction in section 1 is the already-read working specification at
VD-docs commit 8782a293297330f7a354d87ebe62792018fa8866, section 9, as recorded
in the existing investigation's source map. This audit adds no fresh managed
source claim. The wedge interpretation is conditional on the parallel
agent's proof; the reference-prior and symmetry analysis above is an
independent mathematical interpretation audit.
