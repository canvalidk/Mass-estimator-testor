# What could justify the original flat weighting?

**Working contribution 6 — 15 September 2026.** Returns to the original
single-pair equation following the clarification about Gaussian measurement
uncertainty. No new distribution of excitations across trials is assumed here.
Newton II remains a premise. Extended on the same date to examine whether mass
composition supplies the proposed mismatch requirements.

## Main result

**The original flat-magnitude law has an exact Gaussian-conditioning
interpretation.** Start with the Gaussian uncertainty cloud over the two
vectors. Retain possible pairs that point in the same general direction and
have a sufficiently small cross-product magnitude. As that tolerance tends
to zero, the resulting law is precisely

$$
dP_{\rm flat}\propto
g(xu,yu)\,dx\,dy\,d\Omega(u),\qquad x,y>0.
$$

Thus the flat weighting is a coherent choice with a concrete interpretation.
Gaussianity alone still does not choose it: shrinking an angular mismatch or
a distance to the compatible set gives different weights. The choice can now
be described as **which mismatch is held small**, rather than only as an
unspecified reference measure.

A second result identifies a possible structural warrant. In three dimensions,
a nonzero continuous mismatch that is additive in each vector and transforms
appropriately under rotations must be proportional to the cross product.
If those properties are required of the mismatch, constant tolerance in that
mismatch selects the original flat law. These are additional assumptions;
they are not supplied by the existing readout-composition theorem.

The composition follow-up strengthens this result: force-side additivity plus
exchange symmetry of the mismatch magnitude is enough, with the other
premises in Section 4.2 below. Separate acceleration-side additivity is
unnecessary for selecting the area-conditioning event.

## 1. Precisely what is being conditioned

Use dimensionless vector coordinates, for example force divided by
\(\sigma_F\) and acceleration divided by \(\sigma_a\) in the independent
isotropic case. Grant a fully specified nonsingular joint Gaussian density
\(g(X,Y)\) over the uncertain possible vectors, with the observations as its
centers. This grants more than a Gaussian sampling likelihood alone.

The compatible set is

$$
\mathcal C=\{(xu,yu):x,y>0,\ u\in S^{d-1}\}.
$$

For \(d>1\), this set has probability zero under the ambient Gaussian.
Define the parallelogram area

$$
\mathcal A(X,Y)=\|X\wedge Y\|
=\sqrt{\|X\|^2\|Y\|^2-(X\cdot Y)^2}.
$$

In three dimensions it is \(\|X\times Y\|\); in two dimensions it is the
absolute determinant. Consider the positive-probability events

$$
E_\varepsilon={X\cdot Y>0,\ \mathcal A(X,Y)<\varepsilon\}. \tag{1}
$$

The sign condition selects positive collinearity rather than antiparallel
pairs. It is applied **inside the uncertainty cloud of possible latent pairs**.
It does not reject an experiment whose observed vector centers point in
opposite directions. Such centers still give positive Gaussian density on
the compatible set and are included in the numerical checks below.

## 2. Why area conditioning produces a flat magnitude measure

Except on a null set, write

$$
X=xu,\qquad Y=a u+w,\qquad w\perp u,
$$

where \(x>0\). Before taking the limit, \(a\) is the component of \(Y\)
parallel to \(X\), not its magnitude. The identities are exact:

$$
X\cdot Y=xa,\qquad \mathcal A=x\|w\|,
\qquad dX\,dY=x^{d-1}dx\,d\Omega\,da\,dw.
$$

On (1), \(a>0\) and \(\|w\|<\varepsilon/x\). Substitute
\(w=(\varepsilon/x)z\), with \(z\) in the unit ball of \(u^\perp\). The
Jacobian factors cancel:

$$
\boxed{x^{d-1}\,dw
=x^{d-1}(\varepsilon/x)^{d-1}\,dz
=\varepsilon^{d-1}\,dz.} \tag{2}
$$

After dividing by \(\varepsilon^{d-1}\), the unnormalized conditional
integral is therefore

$$
\int_{x,a>0}\int_{S^{d-1}}\int_{\|z\|<1}
g\!\left(xu,a u+\frac\varepsilon x z\right)
\,dz\,d\Omega\,dx\,da.
$$

As \(\varepsilon\downarrow0\), this tends to

$$
\kappa_{d-1}\int g(xu,a u)\,dx\,da\,d\Omega,
$$

where \(\kappa_{d-1}\) is the unit-ball volume. Normalization cancels that
constant, leaving the original flat magnitude measure.

This limit is rigorous even near zero magnitudes. For finite means and a
positive-definite covariance, constants \(C,c>0\) exist such that

$$
g\!\left(xu,a u+\frac\varepsilon x z\right)
\le C\exp[-c(x^2+a^2+\varepsilon^2\|z\|^2/x^2)]
\le C\exp[-c(x^2+a^2)].
$$

The final expression is integrable over the displayed domain. Dominated
convergence applies without excluding the axes or adding an atom at the
origin. It also applies with bounded continuous observables, proving weak
convergence of the ambient conditional laws to the supported compatible law.
Magnitude first moments converge as well: the extra perpendicular term is
controlled by the Gaussian decay in \(\varepsilon\|z\|/x\).

The argument works for arbitrary nonsingular joint Gaussian covariance;
independence and spherical covariance are not needed for this limit theorem.
The Euclidean vector coordinates defining area must be specified. Returning
from dimensionless to physical magnitudes multiplies the reference by constants
that cancel in normalization.

In \(d=1\), the area is identically zero. Conditioning on equal signs already
gives the flat law directly; no shrinking-tolerance argument is needed.

### An intuitive version of the cancellation

In polar coordinates, Cartesian volume contains the radial factor
\((xy)^{d-1}\). An area tolerance permits an angular mismatch of approximately
\(\varepsilon/(xy)\), whose small cap volume is proportional to
\((\varepsilon/(xy))^{d-1}\). These factors cancel.

For example, in three dimensions, doubling both magnitudes makes the radial
volume factor 16 times larger. The allowed angle becomes four times smaller,
making its small cap area 16 times smaller. Equal absolute magnitude increments
are left with equal reference weight.

## 3. The remaining choice is visible in the mismatch

At a compatible pair \((xu,yu)\), write \(r=\sqrt{x^2+y^2}\). A small normal
displacement \(n\) obeys \(\mathcal A\simeq r\|n\|\). In suitable normal
coordinates the identity is exact. The compatible surface volume is
\(r^{d-1}dx\,dy\,d\Omega\).

| Quantity whose tolerance shrinks | Local allowed normal radius | Resulting magnitude weight \(\rho\) |
|---|---|---|
| Area \(\mathcal A\) | \(\varepsilon/r\) | \(1\): original flat measure |
| Angular mismatch, equivalently \(\mathcal A/(xy)\) near compatibility | \(\varepsilon xy/r\) | \((xy)^{d-1}\) |
| Distance to the compatible set | \(\varepsilon\) | \(r^{d-1}\): metric-tube measure |

All three approach the same positive-collinearity set. They differ in how much
nearby ambient volume they retain at each compatible pair. Rotational symmetry,
channel exchange, and unit covariance do not distinguish them.

These are different conditioning constructions, not different coordinate
expressions of one fixed measure. Once a construction is fixed, ordinary
coordinate changes must preserve its measure through their Jacobians.
The general distinction is discussed in
[Bungert and Wacker, Sections 3–4](https://arxiv.org/html/2009.04778v2#S3).
The area limit above is derived specifically for the mass construction.

The earlier uncertainty differences consequently retain their meaning. For
example, at zero observations in three dimensions, flat and tube give the
same half-Cauchy mass-ratio law, while angular conditioning gives a different,
narrower law with finite second mass moment. Recovering the flat law from a
precise limit does not establish that its probability intervals have uniform
repeated-experiment coverage.

## 4. A possible structural reason to select the cross product

Suppose the mismatch is a map \(B:\mathbb R^3\times\mathbb R^3\to\mathbb R^3\)
with these declared properties:

1. It is continuous and additive in each argument separately.
2. For every proper rotation \(R\), \(B(RX,RY)=R B(X,Y)\).
3. It is not identically zero.

Then

$$
\boxed{B(X,Y)=c\,X\times Y\quad\text{for some }c\ne0.} \tag{3}
$$

An elementary proof is short. Continuity and additivity give bilinearity.
Rotations fixing \(e_1\) force \(B(e_1,e_1)\) to lie along \(e_1\). A half-turn
about \(e_2\), combined with bilinearity, forces that value to vanish.
Rotating and scaling gives \(B(X,X)=0\), hence antisymmetry. A half-turn about
\(e_3\) forces \(B(e_1,e_2)=c e_3\); cyclic rotations determine the remaining
basis values. Bilinearity now gives (3).

Under reflections the cross product transforms as an **axial vector**:
\(B(RX,RY)=\det(R)R B(X,Y)\). Requiring ordinary polar-vector covariance under
all reflections instead would force \(B=0\). Its norm is reflection invariant.

If a constant Euclidean-ball tolerance is then imposed on \(B\), the factor
\(|c|\) merely rescales \(\varepsilon\). The limiting law is the original flat
law. In this sense, a bilinear oriented mismatch supplies a specific structural
route to the measure.

**This is an additional warrant to assess.** The existing composition theorem
governs the final readout on an already supplied compatible-pair law. It does
not require an off-constraint mismatch to be additive in each uncertain vector.
Nor does Gaussianity require that property. The result makes a possible
justification explicit without declaring it compulsory.

### 4.1 What mass additivity already establishes

For one propagated joint law of positive compatible magnitudes,

$$
\frac{E[f_1+f_2]}{E[\alpha]}
=\frac{E[f_1]}{E[\alpha]}+\frac{E[f_2]}{E[\alpha]}.
$$

The common acceleration variable is retained. The dual rule adds inverse
masses when force is common. These identities hold for any such probability
law; they do not select how the law was constructed. In particular, they hold
with a law constructed using the tube measure as well as the flat measure.
They concern scalar magnitudes along a shared direction, not arbitrary sums
of force magnitudes with different directions.

The proposed additional requirement concerns an oriented mismatch before
compatibility has been imposed:

$$
B(X_1+X_2,Y)=B(X_1,Y)+B(X_2,Y). \tag{3a}
$$

It says that, at a common acceleration vector, force contributions have
additive vector discrepancies. This is a plausible extension of composition,
but it is stronger: it specifies behavior on incompatible pairs, where the
original readout theorem makes no such assertion. On exactly compatible
pairs all these discrepancies vanish, so their addition is uninformative.

There is a constructive connection to Newton's balance equation. In fixed
standardized coordinates let \(e(X,Y,t)=X-tY\), where \(t\) is the mass
coefficient in those units. At common \(Y\),

$$
e(X_1+X_2,Y,t_1+t_2)=e(X_1,Y,t_1)+e(X_2,Y,t_2).
$$

Choose a linear elimination map \(L_Y\) with \(L_YY=0\). Then
\(B(X,Y)=L_Ye=L_YX\) removes the unknown mass while preserving force-side
additivity. This explains how a discrepancy can inherit the balance equation's
additive structure. Requiring elimination to be linear on arbitrary errors is
the extra commitment; physical mass additivity alone does not require this
choice. The map need not be an idempotent projection. The positive-alignment
condition is still supplied separately, as in (1).

### 4.2 A weaker requirement than additivity in both inputs suffices

Use fixed common reference scales for each channel. Summands are expressed
in the same units; they are not individually divided by different uncertainty
scales before addition. Exchange below means exchanging the standardized
channels and their roles, not physically equating force and acceleration.

Let \(B:\mathbb R^3\times\mathbb R^3\to\mathbb R^3\) satisfy:

1. Continuity and force-side additivity (3a).
2. Proper-rotation covariance: \(B(RX,RY)=RB(X,Y)\).
3. Zero mismatch on positive collinear pairs.
4. Exchange symmetry of magnitude: \(\|B(X,Y)\|=\|B(Y,X)\|\).
5. \(B\) is not identically zero.

The map has no additional preferred direction or tensor argument. Then,
for some constant \(k>0\),

$$
\boxed{\|B(X,Y)\|=k\,\|X\times Y\|.} \tag{3b}
$$

**Proof.** Continuity and (3a) make \(B\) linear in \(X\). Fix
\(Y=yu\), \(y>0\), \(\|u\|=1\). Vanishing on positive multiples of \(u\)
and linearity make the entire parallel component irrelevant. Rotations about
\(u\) force its action on the perpendicular plane to be a scalar multiple
of the identity plus a scalar multiple of the quarter-turn. Thus

$$
B(X,yu)=p(y)[X-(X\cdot u)u]+q(y)X\times u.
$$

The two displayed vectors are perpendicular and have equal norms. Writing
\(s(y)=\sqrt{p(y)^2+q(y)^2}\), the mismatch norm is
\(s(y)x|\sin\theta|\), where \(x=\|X\|\) and \(\theta\) is the angle
between the inputs. Exchanging any nonparallel inputs gives
\(s(y)x=s(x)y\). Hence \(s(y)=ky\) for all \(y>0\), which proves (3b).
Continuity covers zero inputs; nonzero \(B\) forces \(k>0\).

Consequently, a fixed norm tolerance gives exactly the area events of (1),
with their tolerance rescaled by \(k\). The original flat law follows.
No separate acceleration-additivity axiom is needed for this result.

The vector mismatch itself need not be bilinear. For example,

$$
B_\perp(X,Y)=\|Y\|X-\frac{X\cdot Y}{\|Y\|}Y\quad(Y\ne0),
\qquad B_\perp(X,0)=0
$$

satisfies all five requirements and has cross-product magnitude. It is even
in \(Y\), so it cannot be additive in \(Y\). Continuity at zero follows from
\(\|B_\perp(X,Y)\|\le\|X\|\|Y\|\). This shows why selecting the
conditioning event requires less than selecting a unique vector formula.

Exchange symmetry here is stronger than the final identity
\(\widehat\mu=1/\widehat m\). It concerns the mismatch for every input pair.
Reciprocity of a final ratio does not imply it. The norm theorem allows both
polar and axial vector realizations under reflections; their magnitudes agree.
If signed exchange antisymmetry \(B(X,Y)=-B(Y,X)\) is required instead,
force-side additivity also implies acceleration-side additivity, and the
stronger vector theorem (3) applies.

### 4.3 What prevents this from being an unconditional derivation

**Mass composition alone leaves alternatives.** The continuous mismatch

$$
B_{\rm tube}(X,Y)=\frac{X\times Y}{\sqrt{\|X\|^2+\|Y\|^2}},
\qquad B_{\rm tube}(0,0)=0,
$$

has the same collinearity zeros, proper-rotation covariance, and exchange
symmetry of magnitude. It is not force-additive: its values at
\((2e_1,e_2)\) and twice \((e_1,e_2)\) have magnitudes \(2/\sqrt5\)
and \(\sqrt2\). Small tolerances give the tube weight in Section 3.
Using ratio-of-means readouts on propagated laws still preserves both mass
composition identities. This is an explicit counterexample to deriving flat
weighting from those identities alone.

**Force-side mismatch additivity alone also leaves freedom.** The smooth map

$$
B_1(X,Y)=\frac{X\times Y}{\sqrt{1+\|Y\|^2}}
$$

is force-additive, rotation-covariant, and has the same collinearity zeros.
It fails exchange symmetry of magnitude. In three dimensions its small-norm
conditioning law has \(\rho(x,y)=1+y^2\). One can prove the limit globally by
writing \(Y=yu\), \(X=bu+w\): the allowed perpendicular radius is
\(\varepsilon\sqrt{1+y^2}/y\), whose two-dimensional volume cancels the
Cartesian radial factor \(y^2\), leaving \(1+y^2\). Gaussian decay times
this polynomial supplies an integrable bound, including near \(y=0\).

For unit independent Gaussian errors and zero observed centers, \(x\) remains
half-normal. The \(y\) density is proportional to
\((1+y^2)e^{-y^2/2}\). With \(c_0=\sqrt{2/\pi}\), its mean is
\((c_0+2c_0)/2=3c_0/2\), so the standardized mass readout is exactly
\(2/3\), versus \(1\) for the original flat law. Thus the missing symmetry
can affect the actual mass value even when force-side mismatch additivity holds.

The corresponding standardized ratio \(R=x/y\) has density
\(p_R(r)=(r^2+3)/[\pi(1+r^2)^2]\), \(r>0\). Its uncertainty law changes
too, but its \(r^{-2}\) tail still makes \(E[R]\) infinite. This counterexample
isolates the remaining choice; it is not an uncertainty improvement claim.

**A fixed tolerance is also a requirement.** Even after choosing the cross
product, an event
\(\mathcal A(X,Y)<\varepsilon h(X,Y)\) with state-dependent width has local
magnitude weight \(h(xu,yu)^{d-1}\), when the resulting limit is integrable
and no boundary contribution dominates. Width \(h=r\) gives the tube measure;
width \(h=xy\) gives angular conditioning. Choosing the algebraic mismatch
does not by itself prohibit these widths.

Conversely, additivity is sufficient in the stated theorem, not necessary for
the flat law. Squared area has the same flat limit because
\(\mathcal A^2<\varepsilon\) is exactly \(\mathcal A<\sqrt\varepsilon\).
A global strictly increasing transformation of area that is continuous and
zero at zero simply relabels sufficiently small tolerances.

There is therefore a concrete possible bridge: **require additive oriented
force discrepancy, symmetric treatment of the two standardized channels, and
constant tolerance in its norm.** The result follows if these are warranted.
The original readout axioms do not already establish that bridge.

## 5. Area is not the same as Gaussian significance

To examine the error scale, consider a separate sampling experiment with
compatible true standardized means \(x_0u,y_0u\), independent unit Gaussian
measurement errors, and \(r_0^2=x_0^2+y_0^2\). Before any filtering of the
observations, the measured squared area has the exact distribution

$$
\mathcal A_{\rm obs}^2=UV,\qquad
U\sim\chi_d^2(r_0^2),\quad V\sim\chi_{d-1}^2,
\quad U\perp V. \tag{4}
$$

To see this, rotate the two channel columns so their means become \(r_0u\)
and zero. Independent unit Gaussian errors and wedge magnitude are preserved.
Conditioning on the first column, the perpendicular part of the second has
squared length \(\chi^2_{d-1}\), independently of the first column's norm.

In particular,

$$
E[\mathcal A_{\rm obs}^2]=(d-1)(r_0^2+d). \tag{5}
$$

Equal area tolerances therefore do not correspond to equal error probabilities
at different signal strengths. At strong signal, area divided by \(r_0\)
approaches a \(\chi_{d-1}\) residual. At weak signal, the additional noise term
\(d(d-1)\) matters. Dividing by an observed radius does not make this an exact
finite-signal pivot.

Thus the bilinear-area warrant and a requirement of equal statistical
significance are different commitments. This explains why Gaussian error
geometry can motivate the tube measure while algebraic mismatch can motivate
the flat measure. Neither description alone proves superior mass uncertainty.

## 6. Finite-width checks recover the actual equation

[area_conditioning_checks.py](area_conditioning_checks.py) integrates the
finite event (1) under independent unit Gaussian channels in three dimensions.
Observed centers are collinear for this numerical reduction; the proof does
not have that restriction.

Given \(X=xu\), the parallel part of \(Y\) is a one-dimensional Gaussian and
its perpendicular part must lie in a disk of radius \(\varepsilon/x\).
Integrating that disk gives the finite-event probability and conditional
\(E\|X\|\). Swapping the channels gives \(E\|Y\|\), with the same event
probability. This is an ambient-Gaussian calculation independent of the
existing compatible-surface integrator.

The table shows the finite-event norm ratio, which tends to the original
ratio-of-means readout. Finite tolerances are convergence diagnostics, not
proposed new estimators or an allowance for Newton II to fail.

| Observed centers, along one axis | \(\varepsilon=1\) | \(\varepsilon=0.1\) | \(\varepsilon=0.003\) | Original flat limit |
|---|---:|---:|---:|---:|
| \(0,0\) | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| \(3,0\) | 2.969531 | 3.173191 | 3.179124 | 3.179135 |
| \(4,1\) | 2.932227 | 2.991648 | 2.992752 | 2.992754 |
| \(4,-1\) | 6.070370 | 6.747053 | 6.764986 | 6.765017 |

All quantities are standardized; physical mass is multiplied by
\(\sigma_F/\sigma_a\). Both conditional magnitudes and the event normalization
also converge to the independent flat-law calculations.

At zero Gaussian centers, area has a \(\operatorname{Gamma}(d-1,1)\) law, and
the positive-dot sign has conditional probability one half. In three dimensions,

$$
P(E_\varepsilon)=\tfrac12[1-(1+\varepsilon)e^{-\varepsilon}]
\sim\varepsilon^2/4.
$$

This supplies an exact finite-width normalization check, including the region
near zero vectors. At \(\varepsilon=.003\), the numerical scaled probability
is 0.2495005621, tending to 0.25. Conditional mean magnitudes approach
\(\sqrt{2/\pi}\), not the unconditioned three-dimensional Gaussian norm mean.

## 7. What is now justified, and what remains a choice

The original single-pair equation can be read coherently in either of two ways:

- **Supplied Gaussian uncertainty:** impose positive compatibility through
  vanishing area tolerance. The flat compatible law follows from the limit.
- **Gaussian measurement likelihood:** use a shared uniform direction and a
  flat reference over the positive magnitudes. This is an improper prior
  reference; proper finite uniform magnitude cutoffs have the same posterior
  limit as the cutoffs grow. It is not a normalized population model for true
  force and acceleration magnitudes.

These readings should not be conflated or multiplied together as if each
required another prior. The construction depends on what probability input
has actually been supplied.

The answer to the user's clarification is therefore sharper: Gaussian
uncertainty supplies the weight of possible vector values; the original
equation also embodies a particular meaning of approaching compatibility.
That meaning can be stated exactly as an area limit. A sufficient structural
justification needs force-side mismatch additivity and exchange symmetry of
its magnitude, together with the other premises of Section 4.2 and a fixed
tolerance. Full bilinearity is a stronger sufficient route. Neither follows
from Gaussianity or from mass-readout composition alone.

The next unresolved issue is the operational warrant for those extra requirements:
should force discrepancies add before conditioning, should their magnitudes
be invariant under channel exchange, and should tolerance be constant in that
norm? These requirements now have an exact consequence. A noise-metric
distance implements a different choice. The readout theorem does not answer
these questions. This contribution
does not adopt a replacement measure, change the baseline estimator, or add
an across-trial excitation model.

## Verification and provenance

From the workspace root:

```powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/area_conditioning_checks.py' --output '.tools/mass_equation_20260915/single_pair/results.json'
```

The run checks four observation cases at six positive tolerances, quadrature
orders 64 and 96, an extended radial cutoff, channel exchange, and the exact
null event probability. The largest relative quadrature/refinement difference
was \(4.54\times10^{-11}\); channel-swapped event probabilities agreed within
\(4.02\times10^{-14}\). These are checked-case errors, not uniform numerical
error bounds. The perpendicular-disk tail replacement has probability error
at most \(e^{-72}\); the code bounds its supported Bessel arguments explicitly.

A separate 1.2-million-draw check across dimensions 2 and 3 verifies the sampling
moment (5) and the null Gamma event probabilities. The six area-squared sample
means were all within 1.07 Monte Carlo standard errors of their exact values.
Proofs, interpretation audits, and machine outputs are retained separately in
`.tools/mass_equation_20260915/single_pair/`. The limiting law and its structural
warrant were independently reviewed. The one-sided classification, weaker
norm-exchange theorem, and separation from physical composition received
further independent proof and scope reviews; their scratch derivations are in
`.tools/mass_equation_20260915/composition_bridge/`. These are analytical results;
the numerical checks in Section 6 continue to verify the same area limit.

Managed VD-docs `main` was resolved to
[`8782a293297330f7a354d87ebe62792018fa8866`](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866).
At that snapshot, the [catalog](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/catalog.yaml),
[working specification](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md)
(especially Section 9), and
[readout audit](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md)
were checked. The readout audit's common-law composition requirements and its
explicit distinction from a fresh fit were rechecked for the composition
follow-up. Local scope was the earlier Gaussian-conditioning contribution,
the active flat-law integrator, and the new single-pair derivations and checks.
The external probability paper supports the singular-conditioning distinction;
this is not a priority or exhaustive context review. No managed file was edited.
