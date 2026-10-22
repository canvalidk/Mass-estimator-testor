# Limits of scalar-size and tensor-reference repairs

**Record — 28 September 2026.** Written by Codex at the user's request.

**Findings.** The shared scalar variance law with reference measure $dt/t$ does not generally repair ours27's width under heterogeneous instrument geometry. A natural tensor reference law is improper with too few readings, and two rotation-invariant ways of taking its boundary limit give different mass laws.

These are explicit obstructions to particular repair proposals. They do not prove that every learned-size or tensor construction is impossible, and they do not justify announcing an ours28.

## 1. Why a shared scalar size does not suffice

Whiten the observations and let $P_m$ be the projector onto their aligned-candidate subspace, of rank $M$. Put

$$
H(m)=z^\top P_mz,\qquad
Q(m)=\|z\|^2-H(m).
$$

The scalar variance construction in the [zero-reading record](ours27_zero_readings_and_shared_candidate_information_2026-09-28.md) has evidence

$$
e^{-\|z\|^2/2}F_M(H),\qquad
F_M(H)=\int_0^1(1-\beta)^{M/2-1}e^{\beta H/2}\,d\beta.
$$

It is strictly increasing in $H$. It therefore preserves the locations and ordering of the flat evidence's profile maxima. At a stationary point, its curvature is the flat evidence curvature multiplied by

$$
\kappa(H)=2\frac{d}{dH}\log F_M(H)
=E[\beta\mid H].
$$

This is an exact statement about evidence. A fixed smooth starting density is negligible at the leading order in the regular limit below.

### A regular heterogeneous counterexample

Repeat the following two-axis experiment independently $N$ times:

$$
m_0=1,\qquad
\sigma_F=(1,10),\qquad
\sigma_a=(1,1),
$$

with all measurement errors independent. Let the squared standardized candidate lengths on the two axes be

$$
r_1^2=\frac1{10},\qquad r_2^2=\frac{99}{10}.
$$

These specify fixed admissible candidates, not a Gaussian distribution of generating pushes. For example, their physical acceleration components can be taken as $a_j^*=r_j/\sqrt{1+1/s_j^2}$, where $s=(1,10)$; their force components are equal because $m_0=1$.

Use $u=\log m$. The candidate-line angle on axis $j$ is $\theta_j(u)=\arctan(e^u/s_j)$, whose derivative at $u=0$ is

$$
w_1=\frac12,\qquad w_2=\frac{10}{101}.
$$

Resolve each noisy pair into parallel and perpendicular coordinates at $u=0$:

$$
A_j=r_j+\xi_j,\qquad B_j=\eta_j,
$$

with independent standard Gaussian errors. The flat score per repetition is $\sum_jw_jA_jB_j$. The exact score-variance and expected-curvature quantities from the [width record](ours27_full_covariance_width_and_candidate_geometry_2026-09-28.md) become

$$
J=\sum_jr_j^2w_j^2
=\frac{49801}{408040}
=0.12204930889128517,
$$

$$
G=\sum_jw_j^2
=\frac{10601}{40804}
=0.2598029604940692.
$$

The candidate subspaces identify $m_0$ uniquely. With independent repetition the profile estimate is consistent, the usual local expansion applies, and $NJ\to\infty$.

Here $M=2N$ and $H(\widehat u)/N\to12$. The scalar size law concentrates at

$$
\kappa\longrightarrow1-\frac MH=\frac56.
$$

To match sampling variance, however, the required multiplier would be

$$
\kappa_{\mathrm{match}}=\frac{J}{J+G}
=\frac{49801}{155811}
=0.31962441676133263.
$$

The resulting discrepancy is

$$
V_{\mathrm{sampling}}\sim\frac{J+G}{NJ^2},
\qquad
V_{\mathrm{law}}\sim\frac1{N\kappa J},
$$

$$
\frac{V_{\mathrm{sampling}}}{V_{\mathrm{law}}}
\longrightarrow
\kappa\frac{J+G}{J}
=\frac{259685}{99602}
=2.607226762514809.
$$

The standard-deviation ratio is about 1.61469. Both variances tend to zero, but their relative mismatch persists. This is not a vanishing-information or boundary example.

The shared scalar state learns total standardized energy. The directions carrying that energy have different sensitivities to mass. In this example most signal lies on the less sensitive direction, and one common shrinkage factor does not recover the relevant uncertainty.

This disproves general width calibration for this particular scalar reference law. It is not a theorem that every possible scalar correction fails. It motivates retaining candidate shape as well as size, subject to a justified joint model.

## 2. A tensor law retains directional structure

Consider $N$ independent $d$-component readings sharing the same positive-definite measurement covariance. Let $z_i$ be their whitened data and $U_m$ an orthonormal frame for the whitened aligned-candidate subspace. Set

$$
A_i=U_m^\top z_i,\qquad
S_m=\sum_{i=1}^N A_iA_i^\top.
$$

Introduce a shared positive-definite candidate covariance $R$ in these coordinates:

$$
v_i\mid R\sim N(0,R).
$$

The whitened observation covariance is $I+U_mRU_m^\top$. Define the matrix shrinkage

$$
W=R(I+R)^{-1},\qquad 0<W<I.
$$

Apart from a data-only factor, the likelihood and current candidate laws are

$$
L(m,R)\propto
|I+R|^{-N/2}
\exp\left(\frac12\operatorname{tr}(WS_m)\right),
$$

$$
v_i\mid W,m,D\sim N(WA_i,W).
$$

The tensor can be carried through an axis recoding with its full sharing structure. The [axis-splitting record](ours27_axis_splitting_does_not_force_independent_candidates_2026-09-28.md) explains why resetting it to independent scalar models is not a valid representation test.

That resolves one objection to investigating this family. It does not decide which readings should share $R$, how frames should be related across differing instrument designs, or which measure on $R$ should be used.

## 3. The natural conditional tensor measure is improper for small samples

The Fisher volume for the conditional Gaussian covariance experiment suggests

$$
q(R)\,dR=|I+R|^{-(d+1)/2}\,dR,
$$

where $dR$ is Lebesgue measure on real symmetric matrices. This is a conditional nuisance measure within the chosen family, not a uniquely derived joint mass prior.

Under $W=R(I+R)^{-1}$,

$$
dR=|I-W|^{-(d+1)}\,dW.
$$

The resulting conditional hyperposterior is

$$
|I-W|^{(N-d-1)/2}
\exp\left(\frac12\operatorname{tr}(WS_m)\right)\,dW,
\qquad 0<W<I.
$$

For finite data, the exponential factor is bounded above and below by positive constants. It cannot repair boundary integrability. The matrix-beta integral is finite precisely when

$$
N>d-1.
$$

For necessity, let one eigenvalue $\lambda$ of $R$ tend to infinity while the others remain bounded and separated. The eigenvalue Jacobian contributes $\lambda^{d-1}$, giving a tail proportional to

$$
\lambda^{(d-3-N)/2}\,d\lambda,
$$

which is nonintegrable when $N\leq d-1$. For sufficiency, the matrix-beta parameters are $(d+1)/2$ and $N/2$; both exceed $(d-1)/2$ exactly under the displayed condition.

In particular, a single three-dimensional reading has no normalized posterior under this measure. A boundary extension would be an additional construction.

## 4. Two invariant cutoffs give different boundary laws

Take $d=3$, $N=1$, initially at $S_m=0$. Put $X=I-W$ and order its eigenvalues as $0<x_1<x<y<1$. Up to orientation and normalization constants, the density is

$$
x_1^{-3/2}(xy)^{-3/2}
(x-x_1)(y-x_1)(y-x)\,dx_1\,dx\,dy.
$$

Orientation is uniform under rotations.

### Spectral ceiling

Truncate $R<LI$ and let $L\to\infty$. This requires $x_1>\delta=1/(L+1)$. Integrating over $x_1$ gives leading term

$$
2\delta^{-1/2}(xy)^{-1/2}(y-x).
$$

The normalized limiting density of the two nonzero deficits is

$$
p_{\mathrm{spectral}}(x,y)
=\frac32(xy)^{-1/2}(y-x),
\qquad 0<x<y<1.
$$

Equivalently, $y\sim\operatorname{Beta}(2,1)$ and $x/y\sim\operatorname{Beta}(1/2,2)$ independently. Hence

$$
E[\operatorname{tr}X]=\frac45,\qquad
E[\operatorname{tr}W]=\frac{11}{5},\qquad
E[W]=\frac{11}{15}I.
$$

### Determinant ceiling

Instead truncate $|I+R|<L$. This requires $x_1xy>1/L$. The lower bound $x_1>1/(Lxy)$ supplies a leading factor $2\sqrt{Lxy}$ in the $x_1$ integration. The normalized remaining density is

$$
p_{\mathrm{det}}(x,y)=6(y-x),\qquad 0<x<y<1.
$$

Now $y\sim\operatorname{Beta}(3,1)$ and $x/y\sim\operatorname{Beta}(1,2)$ independently, so

$$
E[\operatorname{tr}X]=1,\qquad
E[\operatorname{tr}W]=2,\qquad
E[W]=\frac23I.
$$

Both cutoffs are dimensionless, invariant under rotations, nested and exhaustive of the same positive-definite cone. Their limits differ. In both, $X$ has rank two: only one eigenvalue of $W$ reaches one.

These are genuine normalized limits, not just tail heuristics. After division by $2\delta^{-1/2}$ or $2\sqrt L$, respectively, the integrated densities are bounded by the integrable kernels $(xy)^{-1/2}(y-x)$ or $y-x$. Dominated convergence gives the limits. A fixed finite-data exponential tilt is bounded and continuous on the compact matrix interval, so it can subsequently be restored and normalized.

## 5. The ambiguity reaches the mass law

For a single reading let $A=\sqrt H\,n$, with $\|n\|=1$. Write the normalized limiting evidence factor as

$$
F(H)=E\left[\exp\left(\frac H2 n^\top Wn\right)\right],
\qquad F(0)=1.
$$

Rotation invariance and the expectations just calculated imply

$$
F_{\mathrm{spectral}}'(0)=\frac{11}{30},
\qquad
F_{\mathrm{det}}'(0)=\frac13.
$$

Thus their difference cannot be removed by a mass-independent normalization constant. When $H$ varies with mass, the two limits generally give different mass laws.

There is also a proper singular matrix-beta construction. Let independent Gaussian Wishart matrices be

$$
G_1\sim\operatorname{Wishart}_d(d+1,I),
\qquad
G_2\sim\operatorname{Wishart}_d(N,I),
$$

and set

$$
W=(G_1+G_2)^{-1/2}G_1(G_1+G_2)^{-1/2}.
$$

For integer $N\geq d$, this agrees with the ordinary interior matrix-beta law above, before the data tilt. For $d=3,N=1$, it instead has $\operatorname{rank}(I-W)=1$ and $E[W]=(4/5)I$. It is a third proper continuation, distinct from both cutoff limits. The singular-beta construction is supported by [Díaz-García and Gutiérrez-Jáimez, Definition 2.1 and equations (1), (3)](https://arxiv.org/pdf/0904.2138); the comparison with the two cutoffs is the present calculation.

For an explicit mass example, take unit channel noise, $x=(2,0,0)$, $y=0$, $m=\tan\theta$, and a uniform starting law in $0<\theta<\pi/2$. Deterministic integration gives:

| Continuation | $P(m<1\mid D)$ |
|---|---:|
| Spectral ceiling | 0.2705582155 |
| Determinant ceiling | 0.2859423924 |
| Singular Wishart construction | 0.2557286651 |

Rotation invariance does not select among them.

## 6. A proper current candidate law does not restore ordinary sequential updating

All these limits live on $0\leq W\leq I$. The conditional current-candidate law $N(WA_i,W)$ remains a proper Gaussian, with a degenerate interpretation if needed. Divergence of $R$ does not alone prevent a proper law for the current candidates.

However, the next reading's likelihood in these coordinates contains

$$
|I-W|^{1/2}
\exp\left(\frac12 A_{\mathrm{new}}^\top WA_{\mathrm{new}}\right).
$$

This is zero on the boundary supports just obtained. Multiplying a normalized boundary law by it produces zero total mass, not an updated probability law. Finite-cutoff models can be updated before taking the cutoff limit; that operation need not agree with first taking and normalizing the boundary limit.

The singular Wishart continuation likewise changes support rank with $N$, which ordinary likelihood multiplication cannot create from the previous boundary support. Thus these continuations do not automatically supply a sequential Bayesian rule over all sample sizes.

This is a limitation if ordinary sequential updating is required. It does not prove that every possible batch-defined operational rule is meaningless, or that no proper tensor measure can work.

## 7. Consequence for the next equation

The results narrow the remaining task:

- The scalar $dt/t$ law cannot be adopted as a general calibration repair.
- Retaining a shared tensor is compatible with representation invariance when the actual joint model is carried through the recoding.
- The conditional tensor reference measure needs an additional small-sample construction, and rotation invariance does not make that construction unique.
- A proper continuation must also be checked against whatever sequential-use and persistence requirements the project actually adopts.

The unresolved choice is substantive: different admissible-looking measures already give different finite-data mass laws. Instruments-only specification, uniqueness and persistence still require arguments.

## 8. Provenance and status

Managed context was read at **canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8**, selected through its catalog: the [full-covariance ours27 formulation](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours27/ours27_full_covariance_2026-09-27.md), [learned-push-scale proposal](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours26/ours26_learned_push_scale_2026-09-27.md), and [its status as a hint](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours26/ours26_learned_omega_is_a_hint_2026-09-27.md).

Exact scalar arithmetic and tensor limits were checked in .tools/ours27_direction_2026_09_28/geometry/shared_scale_heterogeneous_width.py and geometry/tensor_boundary_audit.py. The finite-data mass example was checked in measure/tensor_boundary_audit.py within the same scratch directory. Their saved results support the calculations reproduced here.

This local record has not been submitted to VD-docs. No production estimator or managed definition has been changed.
