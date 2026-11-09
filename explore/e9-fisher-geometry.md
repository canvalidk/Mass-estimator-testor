# E9: the full-type mass geometry and the limit of the angle prior

**Codex, 9 November 2026.** `explore/e9-fisher-geometry`, from current
`origin/main` at `80cb9c1`. E6 was pushed through `635a2d3` before this
separate branch; no E6 implementation or premise is imported here.
Exploration input, not adoption. **Checks:** `tests/test_mass_geometry.py`;
research helpers: `met/mass_geometry.py`.

The full type supplies a well-defined curve of noise-whitened hypothesis
subspaces. It does not generally supply one fixed two-dimensional inner
product on the Newtonian rays `(m,1)`. In an explicit structural Gaussian
model, Fisher orthogonality survives full covariance, but normalized
mass volume is independent of arbitrary relative signal sizes **iff the
instruments' metric profiles are proportional**. Orthogonality alone
does not establish that conclusion.

Even the normalized canonical product-curve volume can fail the written
prior symmetry L4 for an ordinary allowed experiment: three independent
scalar noise ratios `1,2,3` have no reciprocal centre at all. This refutes
that proposed general principal-prior rule, not the existence goal. The
permitted input family and production starting law remain unchanged.

## 1. Standard mathematics, and what was already known

The relevant theory is Gaussian Fisher information with nuisance parameters,
canonical Grassmann geometry, and pullback metrics on a one-dimensional
curve. [Edelman, Arias and Smith (1998), §§2.3.2 and 2.5](https://math.mit.edu/~edelman/publications/geometry_of_algorithms.pdf)
supplies the canonical subspace geometry used in T49. This is not a new
metric principle.

[Wootters (1981), *Statistical distance and Hilbert space*](https://doi.org/10.1103/PhysRevD.23.357)
starts with normalized outcome probabilities and distinguishing
fluctuations. The angular result is a property of that experiment. In
the local convention `ds_W^2=I/4`, binary probabilities `sin^2(theta)`
and `cos^2(theta)` have `I=4`; a Gaussian mean `r u(theta)` instead has
`I_theta,theta=r^2`. T70's encoding of observed data is not by itself an
identification of either outcome model with the mass sampling experiment.

[Berger, Bernardo and Sun (2009), introduction, §1.2 and §3.1 Definition 6](https://arxiv.org/pdf/0904.0156)
distinguishes scalar regular reference analysis from sequential
multiparameter analysis and specifies the full observation model. Calling
a Schur-complement square root "the reference prior" would skip the
nuisance construction. None of the calculations below claims to finish
that reference analysis.

Already in main:

- T28/B08 derives orthogonality and uniform angle in its specified
  noise-radius population model and distinguishes the fixed-signal model.
- T49 gives the general graph metric, correlated scalar curve volume,
  proper scalar half-Cauchy special case, latent-volume distinction,
  ambient-restriction obstruction and a heterogeneous zero-padding warning.
- The September 15 legacy Fisher notes already derive the joint
  mean-manifold volume and its incidental-signal coupling. They are read
  for provenance, not run or altered.

The new results here are the full-covariance structural Fisher link,
the exact allocation-independence criterion, the spectral reduction of
anisotropy, and classifications of one-linear-ray form and reciprocal
symmetry. The conditional-prior versus reference-prior distinction is
made concrete with two sections of the same nuisance family.

## 2. The geometric curve under the full declarations

For one reading in ambient dimension `2d`, use its fixed SPD covariance
`Sigma`, declared subspace basis `N`, and

$$C_m=\begin{pmatrix}mN\\N\end{pmatrix},\qquad
D=\begin{pmatrix}N\\0\end{pmatrix},\quad\Omega=\Sigma^{-1},$$
$$G_m=C_m^T\Omega C_m,\qquad
E_m=D^T\Omega D-D^T\Omega C_mG_m^{-1}C_m^T\Omega D.$$

`P_m` denotes the orthogonal projector onto `Sigma^-1/2 C_m`. T49's
canonical squared speed is

$$g(m)=\tfrac12\operatorname{tr}(P_m'^2)
=\operatorname{tr}(G_m^{-1}E_m).$$

Both Gram matrices are positive definite for a nonempty full-rank `N`.
The helper computes them by whitening and QR rather than subtracting
large precision terms. For a proper subspace, **restrict precision**:
use `diag(N,N)^T Sigma^-1 diag(N,N)`, not the inverse of a marginal
projected covariance. Outside-mean observations can carry correlated
noise information in the declared full experiment.

For independent readings, the canonical metric of their direct-sum graph
is `sum_i g_i(m)`. Equal weights follow once that direct-sum construction
is specified. Using its volume as the starting probability law is a
separate inference principle; arbitrary reliability weights do not follow
from the word "canonical."

### Independent channel anisotropy is several effective scalar ratios

If `Sigma=diag(Sigma_F,Sigma_a)`, put

$$A=N^T\Sigma_F^{-1}N,\qquad B=N^T\Sigma_a^{-1}N.$$

Let `lambda_j` be the positive generalized eigenvalues of `(A,B)`.
Congruence by `B^-1/2` and diagonalization give

$$g(m)=\sum_j\frac{\lambda_j}{(1+\lambda_jm^2)^2}
=\sum_j\frac{s_j^2}{(m^2+s_j^2)^2},\qquad s_j=\lambda_j^{-1/2}.$$

The channel covariance matrices need not commute. Consequently even one
anisotropic vector instrument can have the same geometric obstruction
as several scalar instruments with different ratios.

### The geometric volume stays proper, but is not a precision weight

The curve extends smoothly with a strictly positive finite metric to
`m=0`. Channel swap gives `g_D(m)=m^-4 g_(Dswap)(1/m)`, where the
declaration is also swapped. Thus `g_D(m)` has a strictly positive
`m^-4` leading term at infinity. Normalized intrinsic curve volume exists
for every fixed finite nonempty full-type Gaussian declaration. It has
positive finite density at zero and a tail proportional to `m^-2`.

Multiplying a covariance by a positive scalar leaves the graph and `g`
unchanged. At the same physical latent signal, the functional efficient
Fisher information `alpha^T E_m alpha` instead scales inversely with that
scalar. Graph geometry cannot silently supply statistical precision weights.

## 3. Fisher orthogonality survives; size-independent shape generally does not

To connect the geometry to an actual normalized sampling model, expressly
introduce the structural comparison model in whitened coordinates:

$$Z_i\sim N(0,V_i),\qquad V_i=I+\omega_i^2 P_i(m),\qquad\omega_i>0.$$

This makes the physical latent acceleration a Gaussian population with

$$\alpha_i\mid(m,\omega_i)\sim
N(0,\omega_i^2G_i(m)^{-1}).$$

It is isotropy in noise-metric arclength, not generally isotropy with a
mass-independent covariance in physical acceleration coordinates. It is
an explicit conditional comparison, **not** an input licensed by L10/L19.

For a normalized centered Gaussian, the covariance Fisher blocks are
`I_ab=tr(V^-1 V_a V^-1 V_b)/2`. Differentiating `P^2=P` makes `P'`
off-diagonal relative to `range(P)` and its orthogonal complement. With
`k_i=rank(P_i)`, substitution gives

$$I_{mm}=\sum_i w(\omega_i)g_i(m),\quad
w(\omega)=\frac{\omega^4}{1+\omega^2},$$
$$I_{m\omega_i}=0,\qquad
I_{\omega_i\omega_j}=\delta_{ij}
\frac{2k_i\omega_i^2}{(1+\omega_i^2)^2}.$$

Thus the efficient mass block is already `I_mm`. At fixed positive sizes,
its normalized scalar Fisher volume is

$$p_{\boldsymbol\omega}(m)=
\frac{\sqrt{\sum_iw(\omega_i)g_i(m)}}
{\int_0^\infty\sqrt{\sum_iw(\omega_i)g_i(u)}\,du}.$$

This formula is **not** being identified with a completed nuisance-reference
prior. It is the conditional volume whose size independence is in question.

**Allocation theorem.** For positive metric profiles with finite positive
curve lengths, this normalized volume is independent of **every** positive
relative weight vector iff all the profiles `g_i` are proportional.

Proof: if the common density is `p`, then

$$\sum_iw_i g_i(m)=Z(\mathbf w)^2p(m)^2.$$

Subtract two strictly positive weight vectors differing only at coordinate
`i`; this forces `g_i` to be proportional to `p^2`. The converse is
immediate. No zero-size degeneracy is used. The map `omega -> w(omega)`
is increasing from zero to infinity, so all positive weights are reachable.

For scalar correlated channels the profiles are

$$g_{s,\rho}(m)=\frac{s^2(1-\rho^2)}
{(m^2-2\rho s m+s^2)^2}.$$

Two such profiles are proportional iff `(s,rho)` agrees. For example,
uncorrelated ratios `1,4`, with weights `(9,1)` versus `(1,9)`, give
normalized-volume medians **1.154417611233** versus **3.464950604598**.
The cross Fisher blocks are zero in both cases. Orthogonality has not
removed the nuisance dependence of the shape.

**Positive control:** one expressly shared `omega` factors out of the sum,
giving the normalized direct-sum Grassmann volume even for heterogeneous
instruments. Assuming a shared population size is additional structure;
it is not inferred from the instruments or adopted into the definition.
More generally, proportional vector-instrument profiles can suffice even
when each profile itself is not a single scalar-ray metric.

## 4. Holding a nuisance fixed is not coordinate-free marginalization

The functional model already supplies a normalized experiment
`Y~N(C_m alpha,Sigma)` without a signal population. Its efficient mass
information `alpha^T E_m alpha` is T49's result. At one unit scalar
instrument, it is `alpha^2/(1+m^2)`. Evaluating this efficient-information
field along a fixed physical `alpha!=0` section gives a square-root tail
proportional to `1/m`, of infinite length. Along a fixed noise-radius
section `r=sqrt(1+m^2) alpha` it is `r^2/(1+m^2)^2`, of finite length
and half-Cauchy angular form. These trace different signal sections, not
different information values at the same physical point. They do not
declare the nuisance known: in the genuinely known-alpha unit-noise
submodel, the full mass Fisher information would instead be `alpha^2`.

The same issue occurs **within the structural Gaussian family** of §3.
For one unit instrument write the physical acceleration size as `tau`,
so `omega=tau sqrt(1+m^2)`. The full covariance family can be written

$$V=I+\tau^2\begin{pmatrix}m\\1\end{pmatrix}
\begin{pmatrix}m&1\end{pmatrix}.$$

The efficient mass information at the same point is unchanged by this
nuisance reparameterization:

$$I_{m\cdot\tau}=
\frac{\tau^4}{1+\tau^2(1+m^2)}
=w(\omega)\frac1{(1+m^2)^2}.$$

But evaluating that efficient block along a fixed `tau` section gives an
improper `1/m` square-root tail, while a fixed `omega` section gives the
proper angular form. Again, this is not the full Fisher metric of a
subexperiment that declares `tau` known. Fisher
orthogonality is expressed in the `(m,omega)` chart; the `(m,tau)` cross
block is not zero. Transforming the **full** Fisher matrix preserves its
volume and efficient information at each point, as the tests verify.

This is why taking a Schur complement and then declaring one nuisance
section fixed is not automatically an invariant reference-prior analysis.
Nor does the dependence prove no reference construction can work without
a physical population. Reference theory can treat fixed unknown nuisances;
its inference principle, grouping, normalization and regularity would
still need to be specified and checked. That route remains open here.

## 5. Exactly when one fixed two-dimensional mass-ray geometry survives

Consider any finite positive weighted sum of scalar channel metrics in
§3. It has the angular form of the original Newtonian rays `span(m,1)`
under **one fixed SPD two-dimensional inner product**, up to an overall
constant metric factor, iff all `(s_i,rho_i)` agree.

Proof: one fixed inner product gives a constant divided by the square
of one positive-definite quadratic in `m`. A distinct channel quadratic
has its own pair of nonreal roots. At its upper root
`z=s rho+i s sqrt(1-rho^2)`, the leading double-pole coefficient of
`g_(s,rho)` is `-1/4`; a positive channel weight makes it `-w_i/4`.
Thus no distinct pole pair cancels. A one-quadratic form is possible only
with one distinct quadratic; sufficiency follows by summing its weights.

The constant-factor qualification allows identical repetitions, which
multiply the canonical metric but not normalized curve volume. This is
not a claim that heterogeneous instruments lack a one-dimensional curve.
Every positive one-dimensional metric can be reparameterized by arclength
and represented as an arc. Such a nonlinear re-encoding is not one fixed
inner product on the original physical mass rays, and does not by itself
force the corresponding prior rule.

## 6. Reciprocal symmetry is weaker, and sometimes absent for every centre

For the reflection `T_c(m)=c^2/m`, direct substitution gives

$$g_{s,\rho}(T_c(m))|T_c'(m)|^2=g_{c^2/s,\rho}(m).$$

Therefore the normalized square-root volume of a finite positive weighted
sum is self-reciprocal about `c` iff its **aggregated weighted multiset**
of `(s,rho)` is invariant under

$$(s,\rho)\longmapsto(c^2/s,\rho).$$

Sufficiency follows term by term. Necessity follows from the same distinct
poles and their coefficients. Equivalently compare the squared density
identity after its reciprocal Jacobian; normalization introduces no free
constant because pullback preserves total curve length.

Two equal-weight unequal ratios `s_1,s_2` thus have centre
`sqrt(s_1 s_2)` when their correlations agree, despite failing §5's
one-ray test. All ratios equal to `c` remain reciprocal even with different
correlations. Neither fact implies allocation independence.

**A full-type obstruction to the proposed principal prior.** Take three
independent scalar instruments with uncorrelated ratios `1,2,3`, equally
weighted in the canonical direct-sum metric. Reflection must reverse
their ordered pole scales. The extremes require `c^2=1*3=3`; the middle
requires `c^2=2^2=4`. No centre exists.

The curve volume is nevertheless smooth, positive and proper. Independent
quadrature gives length `2.846513208821` and median `1.817989432560`.
Even matching the endpoint density constants forces the different necessary
centre

$$c=\left(\frac{\sum s_i^2}{\sum s_i^{-2}}\right)^{1/4}
=1.790847537606.$$

At `m=1`, the density is `.248411558041`, whereas its reciprocal image
about this necessary centre is `.256055611015`. The exact pole argument,
not this numerical discrepancy, proves failure for every centre.

### Physical swap covariance is not the written within-prior symmetry

Let `J(F,a)=(a,F)` and `Dswap` carry `J Sigma J^T`. Since
`J C_m=m C_(1/m)`,

$$g_D(m)=m^{-4}g_{Dswap}(1/m),\qquad
p_D(m)=m^{-2}p_{Dswap}(1/m).$$

The complete transformed declaration remains covariant, including in the
three-instrument counterexample. Unit changes likewise transform the
metric by the correct coordinate Jacobian.

However L4 explicitly requires **one** prior to be symmetric under
`m -> c^2/m`; O6's fixed-declaration family has an even log-mass shape.
Covariance across transformed declarations does not establish this extra
within-prior requirement. Under the current literal L4, the three-ratio
curve-volume candidate is inadmissible. We do not reinterpret L4 to rescue
it. Auditing whether L4 is the right general physical requirement would
be a separate premise argument, not a conclusion forced by this example.

## 7. Consequence for the uniqueness goal

| Proposed implication | Result |
| --- | --- |
| Full declarations define whitened hypothesis subspaces and canonical curve speeds | yes, conditional on the stated known Gaussian/noise geometry |
| The full type always has one fixed inner product on the linear mass rays | false; §5 gives the exact scalar-channel boundary |
| Orthogonal mass/size Fisher blocks imply size-independent normalized mass volume | false; §3 gives the exact proportional-profile boundary |
| Canonical product-curve volume always lies in the written L4/O6 family | false for the ordinary ratios `1,2,3` |
| These failures prove no instrument-based reference-prior construction exists | not shown; a complete nuisance-reference construction is a different task |
| Geometry alone now selects a new production starting law | not shown; no prior is adopted |

**Provisional sort, Codex:** the single-instrument angular agreement is a
conditional positive result, not an unrestricted uniqueness step. Extending
that step across the current type requires additional justification of
the statistical model, nuisance handling and prior-generating rule. Merely
declaring a population to make the factor cancel would evade L10/L19 rather
than supply that justification. Keeping `pi` as a declared input remains
consistent with TASK; O6's principal-member question stays open.

This does not repeat T49's zero-padding counterexample: that already shows
why internally recomputing a geometric prior is not Z at fixed supplied
`pi`. The new obstruction holds even before data are supplied: the proposed
canonical prior can fail L4 for the declarations themselves.

The conditional geometric objects are usable diagnostics. They do not
change the latent measure in T48, the shrinkage in G3, consistency claims,
or any ruling. No new input is proposed merely to save this construction.

## Validation

All **28 focused E9 checks pass**. With graph geometry, starting-law family,
reading-space, theorem, full-type and audit regressions: **196 passed**.
Independent mathematical/source reviews found no blocker after clarifying
efficient-information sections versus known-nuisance submodels and the
reference-prior citation. The full default suite was not rerun for this
isolated research module; production paths are unchanged.

`tools/audit.py` retains its existing T6/T21/T27 and unused-L24 findings;
`git diff --check` passes. Deterministic quadratures check explicit
mathematical witnesses, not a coverage study or a choice by test score.
Helpers are research-only and do not certify extreme floating point
condition numbers.
