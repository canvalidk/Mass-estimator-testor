# E6: what may be shared between readings?

**Codex, 9 November 2026. Exploration, on `explore/e6-declared-structure`,
from main `80cb9c1`.** Proofs and deterministic checks of specified Gaussian
experiments; the proposed interpretation of L19 is provisional. Tests:
`tests/test_declared_structure.py`. Research helpers:
`met/declared_structure.py`. No signal population is supplied or fitted in
the new comparisons. The production type, equations and premises are unchanged.

## Result and scope

The current independent-reading type cannot exploit a physical relation it
does not declare. Merely putting readings in one coordinate frame adds no
such relation. Declaring a known direction, an unknown common line, equality
of vectors, or identity of an observation makes four different experiments.
Their Gaussian profiles show precisely where their information differs.

Two implications sometimes suggested by the history fail as general claims:

1. Z does not imply absence of shared-signal inference. A profile with a
   common unknown line is Z-safe despite jointly fitting that line.
2. An observed zero with flat standalone mass evidence need not be neutral
   in an experiment with declared shared nuisance quantities. An independent
   zero observation of a shared vector changes its mass evidence.

These are counterexamples to those implications, not counterexamples to the
goal or permissions to use these constructions under current L19. The main
positive result is an exact factorization criterion for fixed linear signal
constraints, including mass-relative likelihoods (§4).

## 1. Standard mathematics and the historical correction

This is constrained functional errors-in-variables estimation: fixed unknown
signals, Gaussian errors and specified relationships among their means.
Profiling a linear constraint is generalized least squares; an unknown common
line leads to best rank-one approximation (Eckart–Young). A covariance model
and a constraint on the means are separate parts of an experiment.

- [Eckart and Young (1936), original publisher record](https://doi.org/10.1007/BF02288367),
  *The approximation of one matrix by another of lower rank*, 1:211–218.
  The publisher abstract and bibliographic record were checked; the full
  original paper was not read. The elementary rank-one proof is below.
- [JCGM 100:2008](https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf),
  §5.2 and F.1.2: uncertainty propagation must retain input correlations;
  common measurement sources can create them. This supports retaining a
  declared acquisition model, not any particular mass likelihood.
- [VIM, §2.48](https://www.iso.org/sites/JCGM/VIM/JCGM_200e_FILES/MAIN_JCGM_200e/02_e.html):
  a measurement model relates the quantities involved in a measurement.
  These sources do not establish that a new declaration passes TASK §5.

**Ours28b was the immediate origin of explicit L19.** Its shared quantity
was a Gaussian signal-size scale per equal-noise class, not one actual push
or a coherent sum of vector readings. Its factor was

$$L_g(m)=\int G(d\beta)(1-\beta)^{M_g/2}\exp[\beta H_g(m)/2].$$

Appending zero leaves H fixed but increases M: it updates the shared scale
and thereby changes relative mass evidence. The formula's zero-padding
effect is independently checked with its exact hypergeometric expression.
The historical readout/grid numbers are not reconstructed or edited.

The 28 September withdrawal also excluded a learned population of scales
and a shared direction tensor, while explicitly keeping a **known direction
declared per reading**. Ours28 was initially kept; on 30 September its
independent draws from a common size law were also rejected as a population
of setups. Thus L10's no-population restriction is stronger than merely
forbidding one fitted shared parameter. A learned-omega proposal on
27 September was an earlier precursor, rejected as unforced and score-selected.

Sources: [withdrawal §§1–2](../docs/vd-docs-context/ours28b_no_borrowing_between_readings_2026-09-28.md),
[ours28b §§1,3](../docs/noise_vs_signal/ours28b/ours28b_one_scale_per_noise_pair_2026-09-28.md),
[30 September reason and scope §2](../docs/noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md),
[27 September precursor](../docs/noise_vs_signal/ours26_learned_omega_is_a_hint_2026-09-27.md).
The 30 September statement explicitly allows reconsidering L19 if its reason
is upheld: the definition must not rely on a repeated experimental setup.

## 2. What a common frame does and does not supply

First fix independent unit isotropic Gaussian errors in both channels.
Work in noise units, so physical mass is the instrument ratio times m below.
Write s=m/sqrt(1+m²), c=1/sqrt(1+m²), and

$$b_i(m)=sF_i+ca_i,\quad r_i(m)=cF_i-sa_i,\quad
E=\sum_i(\|F_i\|^2+\|a_i\|^2).$$

At a candidate mass the change from (F,a) to (b,r) is orthogonal.
Every model has profile log likelihood -(E-H)/2, where:

| Declared true accelerations | H(m) |
| --- | --- |
| arbitrary vector at each reading | $\sum_i\lVert b_i\rVert^2$ |
| arbitrary signed multiple of known unit n at each reading | sum_i (n·b_i)² |
| arbitrary signed multiples of one unknown unit n | largest eigenvalue of sum_i b_i b_iᵀ |
| one vector equal across N independently observed pairs | $\lVert\sum_i b_i\rVert^2/N$ |

**Proof.** The perpendicular residual r is unavoidable. For free vectors,
fit every parallel b exactly. For known n, project each b onto n. For an
unknown n, minimizing sum_i ||b_i-t_i n||² first gives t_i=n·b_i, then
maximizes nᵀ(sum_i b_i b_iᵀ)n over unit n. For equality, fit all b_i by
their common least-squares mean. No prior on n, amplitudes or signal sizes
enters these profile calculations. This does not select profiling as the
definition's nuisance-elimination rule or prove calibration.

The first row holds even when all readings already use one laboratory
frame. Cross-reading dot products are then unnecessary for this profile.
A common frame makes them meaningful; it does not make a common signal true.
The known-line row is the already proposed per-reading input N_i=n. The
unknown-line row is different: the line's value is learned jointly, requiring
a new relation beyond the current independent N_i declarations.

For free isotropic readings each pair's Gram is a complete invariant under
its own O(d) frame changes: isometries between the two-vector spans extend
to orthogonal maps of the whole space. Cross-reading dot products do not
survive independent frame changes without a relative-frame declaration.
If physical frame maps T_i are declared, a common-line model is
a_i*=t_i T_i n. Changing recording frame by R_i also changes T_i to R_i T_i.
Transporting both preserves the experiment; holding T_i fixed changes it.

## 3. A common unknown line supplies different information, while passing Z

Take two 2D datasets with unit errors:

| | first (F,a) | second (F,a) |
| --- | --- | --- |
| P | (2e1,e1) | (0.5e1,e1) |
| O | (2e1,e1) | (0.5e2,e2) |

The per-reading Grams are identical. Free-vector evidence is identical:
H_free=(2s+c)²+(0.5s+c)², with best positive m=(9+sqrt(481))/20,
approximately 1.54658561. With a common unknown line declared, P retains
that profile; O has H_line=max((2s+c)²,(0.5s+c)²), whose best m is 2.
The difference depends on m, so normalization cannot remove it.

Both observed datasets are reachable in the constrained Gaussian experiment:
its positive error density does not require the observed vectors to obey
the true-signal constraint. A finite neighborhood preserves the profile
difference. No noise-free rank test is being applied to observations.

Appending an exact-zero pair adds zero to sum_i b_i b_iᵀ. Consequently this
shared-line profile is exactly Z-safe. Multiplying exp[-K_line(m)/2] by
any fixed proper mass prior gives a proper mass law because 0<=K_line<=E.
This proves the law-level Z counterexample, not only point neutrality.
It **fails strict current L19** by jointly learning n, even though no
signal-size population is declared. Z alone cannot certify L19.

The shared-line statistic still sums own-reading outer products. Its top
eigenvalue depends on relative spatial orientation (for example through
tr(M²)=sum_ij(b_i·b_j)²). It is not the signed coherent sum in the equality
row. No quantum phase or interference postulate has been derived.

## 4. When does a Gaussian profile factor over readings?

Let the stacked true acceleration lie in a **fixed declared linear
subspace** S of V=direct_sum_i V_i. Allow arbitrary SPD covariance within
each reading, independent between readings. Whiten blockwise and define

$$A_m=\operatorname{blockdiag}_i[\Sigma_i^{-1/2}C_i(m)],\qquad
C_i(m)v=(mv,v),\qquad P_m=\operatorname{proj}(A_m S).$$

The full profile cost is K_m(z)=||z||²-zᵀP_m z.

**Factorization criterion.** K_m is a sum of own-reading functions for
every observation z iff S is a Cartesian product over the partition.
Moreover, for two distinct finite masses m and t, K_m-K_t factors that
way iff S is a Cartesian product. Thus allowing a mass-independent mixed
data constant does not evade the criterion for this graph model.

**Proof.** Additivity of a quadratic form is equivalent to all its mixed
blocks vanishing, as its mixed finite difference at u_i,v_j is twice the
corresponding bilinear form. P_m is block diagonal iff its range is a
direct sum of subspaces supported on individual recording blocks. Each
block of A_m is injective onto its graph plane, so pulling back that
decomposition gives S=direct_sum_i S_i. The converse is immediate.

For relative costs, additivity says the mixed blocks of P_m and P_t
coincide. The range of their common (i,j) block lies in both whitened
graph planes at reading i. These planes have zero intersection for m!=t:
(mv,v)=(tw,w) implies v=w=0, and the same invertible whitening preserves
this fact. Their common mixed block is therefore zero. Apply the first
part. This proves the criterion for relative evidence too.

This classifies Gaussian **profiling**, not arbitrary integrated likelihoods,
Student factors, g, eta or quoted widths. Cross-reading noise correlations
invalidate the block-independence assumption. The unknown-line set is a
nonlinear union and is outside this linear theorem. Shared observations
require a unique-acquisition model (§6), not independent SPD blocks.

## 5. Equal-vector acquisitions: coherent averaging and the scope of Z

If every independent acquisition observes the same fixed acceleration,
the equality row has genuine cross terms b_i·b_j. Orthogonal transformation
of the reading index separates the normalized sum from N-1 contrasts.
The contrasts have zero mean in both channels, independent of m. Exactly,

$$K_{equal}(m)=\sum_i\|F_i-\bar F\|^2+
\sum_i\|a_i-\bar a\|^2+
K_{one}(m;\bar F,\bar a,\Sigma=I_{2d}/N).$$

Hence averaging each channel separately, retaining covariance I/N, loses no
mass information **in this declared model**. Force is not adjusted toward
acceleration before E10. In the free-vector model the discarded contrasts
can themselves identify mass: (F1,a1)=(2,1), (F2,a2)=(-2,-1) have zero
means but a nonconstant profile minimized at 2. Averaging is then lossy.
Splitting a record while retaining all constraints differs from this
coarsening. The earlier [axis-splitting audit](../docs/ours27/ours27_axis_splitting_does_not_force_independent_candidates_2026-09-28.md)
already made that general distinction; the equations here specialize it
without a shared population law.

Append a new independent pair whose observed value is zero, with the
same equality declaration. E is unchanged but H_equal becomes N/(N+1)
times its old value. For one pair (2,1), the profile-cost difference
K(1)-K(2)=1/2 becomes 1/4 after zero-padding. A standalone zero's flat
mass profile did not imply neutrality in the shared-vector experiment.
The added observation informs the shared vector. Requiring universal Z
after this type extension would reject the ordinary Gaussian likelihood
for that extension. This does not contradict Z on the current free type.

## 6. A shared observation is not a new observation

If acceleration A is acquired once and reused alongside N independent
force acquisitions of the same body at that acceleration, unit-noise
unique data are (F1,...,FN,A), with means (m u,...,m u,u). Then

$$H_{once}(m)=\frac{\|m\sum_iF_i+A\|^2}{Nm^2+1}.$$

Treating the N copies of A as independent instead gives
H_copy=||m sum_i F_i+N A||²/[N(m²+1)]. This changes precision.
For F=(1,3), A=1, the correct cost is 2+(2-m)²/(1/2+m²);
the copied-independent cost is 2+(2-m)²/(1/2+m²/2).
Copying a report with its full covariance has singular duplicate blocks
and must retain their deterministic support. It cannot be passed to the
current SPD Reading constructor as fresh independent information.
This agrees with the shared-observation issue already found in O2/T54;
the E6 scalar witness does not replace that joint-composition work.

## 7. What is forced, what is an additional declaration, what remains open

**Forced given a specified Gaussian experiment:** preserve its covariance,
mean constraints, frame transports and observation identity under recoding.
For any invertible T, transforming y, Sigma and X together preserves
min_beta ||y-X beta||²_(Sigma^-1). The determinant density change is
independent of m. Reordering or relisting the same experiment is such a
recoding; replacing its constraint set or dropping contrasts is not.

**Already available:** a known per-reading direction via N_i (proposed L18).
**Additional conditional models here:** a common unknown line, equal true
vectors, relative frame maps, or a single reused acquisition. Neither a
common instrument nor observed alignment establishes these facts. Merely
declaring them before data is insufficient: the apparatus, synchronization,
physical model or calibration must warrant them. Real deviations and
uncertainty in these relationships would have to be retained. No arbitrary
grouping knob is licensed by this exploration.

**Proposed reading of the no-cheating reason (Codex, provisional):** generic
mass determination must work without a repeated setup or a population of
signals. That leaves a distinct question whether independently warranted
relations can enter a richer experimental determination, or must be reduced
to the existing reading type before its mass function is called. The equality
example supplies one exact separate-channel reduction; an unknown common
line generally needs joint internal nuisance elimination. These routes
should not be silently identified.

Thus this exploration does not restate L19. Its strict wording excludes the
unknown-line comparison; its historical reason does not by itself prove
that every richer, warranted experiment must discard its extra relations.
Whether a common-but-unknown line passes TASK §5 is still unsorted. It would
be a restriction on allowable signals, not a population law, but that
mathematical distinction alone is not a permission argument.

T70's independent-frame picture is faithful to the current generic type.
Frame bookkeeping alone is not a physical decoherence mechanism. The
existence and authority of a shared nuisance relation, its uncertainty,
and its preservation through the prediction trace are the relevant questions.
No coverage, universal consistency, unique finite law or new probe is claimed.

## 8. Continuation: exact reduction for a fixed linear mean space

This is Gaussian sufficiency and generalized least squares. The minimal-
sufficiency test used in §10 is the likelihood-ratio criterion, stated in
[Balakrishnan, 36-705 Lecture 12, Theorem 12.1](https://stat.cmu.edu/~siva/teaching/705/lec12.pdf).
The specialization and square-completion proof here are explicit so the
result does not rely on an unverified appeal to a statistical theorem.

Let B be a fixed p-by-q full-column-rank matrix and L=diag(B,B). Suppose
both channel means lie in col(B), independently parameterized as
theta=(u_F,u_a), and y~N(L theta,Sigma) with fixed known Sigma SPD.
Sigma may correlate different acquisitions and both channels. Define

$$V=(L^T\Sigma^{-1}L)^{-1},\quad T=VL^T\Sigma^{-1},\quad
\widehat\theta=Ty,\quad e=y-L\widehat\theta,\quad R=e^T\Sigma^{-1}e.$$

Then, for every coefficient pair theta,

$$(y-L\theta)^T\Sigma^{-1}(y-L\theta)
=R+(\widehat\theta-\theta)^TV^{-1}(\widehat\theta-\theta).$$

**Proof.** L^T Sigma^-1 e=0 kills the cross term. Moreover theta_hat has
mean theta and covariance V; e has mean zero and covariance Sigma-LVL^T.
Their cross-covariance is zero. Joint Gaussianity makes them independent,
and the residual distribution is parameter-free. Thus the full density
factors into h(y) times the normalized N(theta,V) density of theta_hat.

Only now impose theta=(m t,t). The identity proves, at every candidate,

$$K_{original}(m)=R+K_{compressed}(m;\widehat u_F,\widehat u_a,V).$$

It preserves the entire likelihood in (m,t), not just its maximizing mass.
Integrating it preserves an integrated likelihood if the **same nuisance
measure** is retained or pushed forward consistently. It does not select
that measure. In particular the original and compressed declared-metric
graph arclengths agree: their Gram is the same
C_m^T L^T Sigma^-1 L C_m=C_m^T V^-1 C_m. Profile equality does not in
general authorize substituting a newly chosen flat reference measure.

The coefficient dimension q can exceed three and need not be a physical
spatial coordinate. A sufficient coefficient pair is not automatically
an instance of the current `Reading` type. Nor does statistical sufficiency
alone establish that g, eta, counts or a separately calibrated report must
ignore the discarded residual; those rules require their own argument.

## 9. Channel locality and the provenance of the reduction

Separate-channel compression holds exactly when T's two off-diagonal
channel blocks vanish. It holds if Sigma is block diagonal by channel.
It also holds for Sigma=C tensor R_index, even with correlated channels:

$$T=I_2\otimes[(B^TR_{index}^{-1}B)^{-1}B^TR_{index}^{-1}],\quad
V=C\otimes(B^TR_{index}^{-1}B)^{-1}.$$

It need not hold for the currently allowed within-reading correlations.
Consider two scalar acquisitions of one common force and one common
acceleration, with independent errors between acquisitions and
Sigma_1=[[1,rho],[rho,1]], Sigma_2=[[1,-rho],[-rho,1]]. Direct GLS gives

$$\widehat u_F=\bar F-\rho(A_1-A_2)/2,\qquad
\widehat u_a=\bar A-\rho(F_1-F_2)/2,\qquad
V=(1-\rho^2)I_2/2.$$

At rho=1/2, compare F=(2,2), A=(1,1) with F=(2,2), A=(2,0).
The separate channel means and their covariance agree. The exact reduced
pairs are (2,1) and (1.5,1), so the best masses are 2 and 1.5. Opposite-
channel contrasts contain information through declared noise correlations.
For A=(1+delta,1-delta), the full profile is exactly

$$K(m)=2\delta^2+
\frac{(2-\rho\delta-m)^2}{[(1-\rho^2)/2](1+m^2)}.$$

The compression never fits a mass or enforces u_F=m u_a. **But this alone
does not make it permissible upstream of E10.** It already assumes that
both force and acceleration obey B. If only acceleration's constraint is
independently supplied, the force constraint may itself come from Newton II
and persistence. Removing m from the algebra does not remove that premise.
A wholly acceleration-only declaration has the larger ambient mean basis
diag(I_p,B); it must not silently be replaced by diag(B,B).

This qualifies the separate-channel reduction in §5. It is an exact
conditional calculation; its location in a VD trace depends on the source
of *both* mean constraints. A joint statistic using both observations can
be evaluated inside E10, where both are available. Substituting it for an
independently obtainable outward force or acceleration input needs the
provenance argument in [why reconciliation cannot be upstream](../docs/vd-docs-context/why_reconciliation_cannot_be_upstream_of_E10.md).
The present theorem supplies no loophole around that argument.

## 10. Why a common unknown line has no equivalent fixed linear compression

Let p=Nd and consider the full experiment
Y~N((m vec(A),vec(A)),Sigma), m>0, rank(A)<=1, with fixed SPD Sigma.
This is a nonlinear Gaussian mean family, with singular strata at zero or
rank-deficient signals. Although its dimension is small, its linear span
is all R^(2p).

**Proof.** Elementary one-entry matrices are rank one and span R^(N d).
For each such a, two distinct positive masses give (m1 a,a) and (m2 a,a).
Their linear combinations span (a,0) and (0,a). These span both channels.
For y and y', their Gaussian density ratio is parameter-independent iff
(y-y')^T Sigma^-1 mu is constant over every allowed mean mu. Zero is an
allowed mean, so the constant must be zero. Full span and invertible Sigma
then imply y=y'. The likelihood-ratio criterion makes the raw observation
minimal sufficient for the **full (m,A) experiment**. Every noninjective
fixed linear compression loses information for that experiment.

This does not forbid smaller statistics for an explicitly selected mass
profile, nuisance integration or invariance reduction. Indeed §2's profile
can be stored in three d-by-d scatter matrices; profiling intentionally
eliminates the signals. Nor does this claim rule out pathological injective
encodings into fewer real coordinates. The substantive distinction is
loss of information, not the number of scalar labels in an encoding.

## 11. A supplied direction with uncertainty is still a shared nuisance

External supply of a direction is not equivalent to its being exact. A
finite-calibration example shows why this matters even at infinite signal.
Let the true direction be n0=e1, true mass m0=2, and means
F*=m0 A_signal e1, a*=A_signal e1. Independent channel covariances are
Sigma_F=diag(1,1/4), Sigma_a=I2. Independently calibrate the direction with
Z~N(n0,tau² I2), fixed tau>0.

Suppose n=Z/||Z|| is installed as an exact N_i. Write W_F=Sigma_F^-1,
W_a=Sigma_a^-1. At increasing signal the fixed-line fit has limit

$$m_{plug}=m_0\frac{n^TW_F n_0/(n^TW_F n)}
{n^TW_a n_0/(n^TW_a n)}.$$

**Proof.** At fixed n, fit the force and acceleration amplitudes separately
in their respective metrics. Where their signs agree, their ratio is a
positive mass achieving the joint optimum. Their large-signal limits are
the two metric projections displayed above. At n=(e1+e2)/sqrt(2) the ratio
is 2/5, so the plug-in result is 0.8 while the truth is 2. A neighborhood
of such calibration outcomes has positive probability for every fixed
tau>0. Measurement noise of order one does not remove this limit as
A_signal grows. Proportional channel metrics make the ratio cancel, so
the isotropic slice can hide this defect.

The correctly specified profile retains the calibration likelihood:
K_signal(m,n)+||Z-n||²/tau², with unit n and the same signal constraints.
The finite calibration cost is O_p(1); any direction bounded away from the
true line costs order A_signal² in the signal data. On the true line both
nonzero channels fix m0. The profile therefore recovers the line and mass
in this amplification limit. This is a profile consistency statement, not
a proper-law or interval-calibration theorem.

This is not a counterexample to correctness of an actually exact N_i. The
plug-in procedure falsely declares a noisy calibration to be exact. The
lesson for TASK §5 is narrower: an externally warranted relation needs its
uncertainty and provenance, not just an externally supplied numerical value.
Retaining finite calibration then introduces a legitimate statistical
shared nuisance in the specified experiment; its admission into the mass
definition's input type remains unsorted.

## 12. Limit of the Gaussian reduction: a shared unknown noise scale

The ancillary claim in §8 requires fixed known Gaussian covariance. With
one additional shared scale v, conditional GLS still separates the residual
R from the compressed mass cost K, but both conditional densities depend
on v. For v~IG(nu/2,nu/2), integrating the normalized Gaussian density
against the same declared-metric graph arclength, independent of v, gives

$$L_{full}(m)\propto(\nu+R+K(m))^{-(\nu+r_0+q)/2},$$

where r0=2p-2q is the ancillary dimension and q the compressed graph's
residual dimension. Dropping R and using only the compressed pair would
instead give (nu+K(m))^(-(nu+q)/2). Relative mass evidence differs.
The residual is ancillary for the mean **conditional on v**, but informs
v after mixing. This is the outside-span warning of T48 §2 applied to the
proposed measurement reduction. It is not an asserted common scale in the
current per-reading finite-nu model. With separate latent scales even the
conditional GLS transformation can depend on their values.

Consequently E6 has not produced a universal preprocessing adapter. It has
identified a lossless fixed-Gaussian reduction, its precise provenance
requirements, a curved-family obstruction to the same route, and a finite
calibration failure if uncertainty is erased. These are additional tests
of any proposed input permission; they do not settle that permission by
declaring a relationship to be an input.

## 13. Retain calibration as a joint state: a conditional construction

The standard mathematics here is nuisance marginalization, sequential joint
updating and propagation of a joint distribution through a deterministic
quantity map. The mixture-versus-reset distinction is already explicit in
the [axis-splitting audit](../docs/ours27/ours27_axis_splitting_does_not_force_independent_candidates_2026-09-28.md).
The following specialization uses a calibrated physical direction rather
than a population of push sizes. It does not select this construction as
the mass definition.

Let q0(dn) be a supplied proper direction reference and c(n) an independent
calibration likelihood with a finite positive integral. Write
q_C(dn)=c(n)q0(dn)/integral(c dq0). For the conditional experiment, each
reading has its own unrestricted signed amplitude t_i, a common direction
n, and fixed known Gaussian covariance Sigma_i independent of m,n.
The Gaussian errors are independent between readings and independent of
the external calibration observations.
Use the full ambient profile cost

$$K_i(m,n)=\min_{t_i}\|y_i-(m n t_i,n t_i)\|^2_{\Sigma_i^{-1}},\quad
f_i(m,n)=\exp[-K_i(m,n)/2].$$

Consider the **specified working joint law**

$$Q_D(dm,dn)\propto\pi(dm)q_C(dn)\prod_{i\in D}f_i(m,n).$$

The same factors also result from integrating the normalized Gaussian
density against each line's declared-metric arclength. That is an explicit
measure convention, not a new forcing argument for T48/G19. A proper q0
is uncertainty about one fixture parameter; no distribution of the t_i
is declared. Whether such fixture information belongs to the definition's
input remains a separate question. No default q0 or new production prior
is proposed.

**Propriety and Z.** Zero amplitude is permitted, so
0<=K_i<=E_i=y_i^T Sigma_i^-1 y_i. Hence
exp[-sum_i E_i/2]<=product_i f_i<=1. Every proper pi and q_C give a
strictly positive finite normalizer. The bounds also preserve divergent
nonnegative moments of pi at finite data; proper does not imply finite
mass variance. An exact-zero observed reading has K_i=0 at every (m,n).
It therefore leaves the full joint law, both marginals and all their
pushforwards unchanged. The empty experiment gives pi times q_C.

This Z statement needs the stated covariance and amplitude domain. A
direction-dependent noise law or constraints excluding zero amplitude
need a separate calculation. It does not apply automatically to every
calibrated experiment.

## 14. A mass-neutral history can still change a future mass update

Keeping Q_D permits exact associative updates:
Q_(D+B) proportional to Q_D times product_(i in B) f_i. Marginalizing first
instead requires a conditional direction law:

$$p_{D+B}(m)\propto p_D(m)\,
E_{Q_D(n\mid m)}[\prod_{i\in B}f_i(m,n)].$$

The missing object cannot in general be recovered from p_D.

**An exact Gaussian witness.** Let n take e1 or e2 with calibrated weights
w,1-w, 0<w<1, and use unit isotropic errors. The two-reading history
(F,a)=(e1,e1), (e1,-e1) has product factor e^-1 at e1 and e^-2 at e2,
independent of mass. To check, its explained energies are (s+c)^2 and
(s-c)^2 along e1, summing to 2; total observed energy is 4. Along e2
the explained energy is zero. The mass marginal is exactly pi, just as
for empty history, but the direction weight becomes

$$w'=\frac{w e}{w e+1-w}.$$

A later reading (F,a)=(2e1,e1) has factors
A(m)=exp[-(2-m)^2/(2(1+m²))] and b=exp(-5/2). Its mass evidence is
w A(m)+(1-w)b after empty history, and w' A(m)+(1-w')b after the block.
Since A is nonconstant and w'!=w, these two functions are not proportional.
For any prior supporting masses with different A, the normalized laws
therefore differ. Identical mass marginals require different future updates.

This is a nonzero mass-neutral block, not an appended exact-zero reading;
§13's exact-zero neutrality still holds. It refutes universal updating from
the mass marginal alone once the shared fixture is admitted. It does not
contradict L17's same-law/same-point readout rule or T12's use of additional
carried state. No stronger posterior-only closure axiom is adopted.

## 15. Resetting a direction and repeating calibration change the experiment

For two blocks A and B the common-fixture mass evidence is
E_(q_C)[f_A f_B]. Separately marginalizing and multiplying gives
E_(q_C)[f_A] E_(q_C)[f_B]. The difference is Cov_(q_C)(f_A,f_B).
The latter calculation assigns independent directions to the blocks.

For directions e1,e2 with weights 3/4,1/4 and the unit-noise pairs
(2e1,e1), (0.5e2,e2), common/reset evidences are respectively
0.331924444751 / 0.384720474049 at m=1 and
0.417832622069 / 0.463141157284 at m=2. Their ratio depends on mass;
normalizing the mass law cannot repair the reset. The same marginal
direction uncertainty can describe either one persistent fixture or two
independent fixtures. Its numerical probabilities do not encode identity.

A different error occurs if already updated joint block laws are multiplied
while retaining a common n. With densities in the same reference, on the
common prior support (where the denominator is nonzero),

$$Q_{A+B}\propto Q_AQ_B/(\pi q_C).$$

Dividing out only pi repeats the same calibration. For uniform q0 and
calibration odds 3, the repeated likelihood gives odds 9 (direction weight
0.9 instead of 0.75) before signal data. A genuinely new independent
calibration contributes another likelihood; copying the same calibration
record does not. These two mistakes are distinct: resetting the shared
parameter and counting its evidence twice.

The finite-state tests are exact supplied-reference experiments, not
quadrature approximations or a proposed discrete default. Gaussian
calibration can produce these weights: under uniform direction reference,
unit calibration variance and observation Z, c(e1)/c(e2)=exp(Z1-Z2).
The mismatches persist under small perturbations of observations and under
concentrated smooth proper references.

## 16. Composition can inherit dependence from geometry alone

Two bodies can have independent mass priors and independent measurement
errors while sharing an uncertain fixture direction. Their joint law is

$$Q(dm_1,dm_2,dn)\propto
\pi_1(dm_1)\pi_2(dm_2)q_C(dn)f_1(m_1,n)f_2(m_2,n).$$

Given n, the masses are independent. After marginalizing n they generally
are not. For bounded functions g1,g2 the law of total covariance gives
Cov(g1,g2)=Cov_(n|D)(E[g1|n,D],E[g2|n,D]). Use bounded functions here:
the finite-law mass moments need not exist under a half-Cauchy reference.

**Finite witness.** Each supplied mass prior is uniform on {1,2}; the
direction reference gives e1,e2 equal weight. Body 1 observes (4e1,2e1),
body 2 (4e2,2e2), with unit errors. Put a(1)=e^-1, a(2)=1. The likelihood
of the direction perpendicular to each observed pair is e^-10, so the
joint mass weights are proportional to a(m1)+a(m2). Normalization gives
P(1,1)=0.134470710685, P(1,2)=P(2,1)=1/4 and P(2,2)=0.365529289315.
Thus P(m1+m2=3)=0.5. Multiplying even the marginals of this same full
joint instead gives 0.473305966621. The missing calibration dependence
changes the actual sum distribution.

This is a pushforward of persistent body masses under an explicit joint
reference, not a new common-acceleration acquisition or a coverage result
for L5(c). Independent errors alone do not supply the copula of inferred
masses. It instantiates T51's existing obligation with shared geometry
rather than copied acceleration data. Both sum and parallel sum use the
retained joint law; neither operation manufactures independence.

## 17. What this continuation decides, and what it does not

Full ambient evidence must travel while n is uncertain. At fixed n,
K_full(m,n)=K_inside(m,n)+B_out(n), where B_out is the squared noise-metric
distance outside the two-channel span of n. B_out is mass-independent
conditional on n but updates n. Dropping it before direction marginalization
changes mass evidence, just as dropping §12's residual changed a shared
noise scale. The tests contain an explicit Gaussian witness.

These are compatibility results for bookkeeping, not a survivor equation:

| Obligation | Result for the specified joint construction |
| --- | --- |
| proper joint mass/direction law at finite data | passes under the supplied proper references |
| exact-zero law neutrality | passes under §13's Gaussian/noise/amplitude conditions |
| combining independent data blocks | passes carrying joint (m,n), with calibration once |
| universal update from mass marginal alone | fails, §14 |
| independent multiplication of body marginals for composition | fails, §16 |
| bounded confidence on fixed incompatible data (T30) | fails in general for these Gaussian factors |
| literal current L19 | fails: the signals still update a shared direction |
| unique nuisance measure, finite-law selection, general coverage | unproved |

For the T30 failure, take one fixed pair F=(2,0), a=(1,1), directions
e1/e2 equally weighted and covariance sigma² I4. The log evidence contrast
between m=2 and m=1 grows as 1/(4 sigma²): it is 6.25, 25 and 100 for
sigma=0.2,0.1,0.05. Retaining proper calibration does not cure the
fixed-Gaussian confidence problem or replace the current probe line.

**Proposed input sort, Codex, provisional.** If a richer experiment admits
fixture calibration, it needs the identity and persistence of the calibrated
parameter, its uncertainty, and the identities of the calibration observations.
Whether two readings refer to one fixture or two cannot be deduced from
their matching marginal direction laws. This information has an experimental
source and is a candidate measurement input, unlike inventing a population
of pushes from a pattern in the readings. But exact joint bookkeeping still
learns shared signal direction and does not satisfy literal L19. The next
goal question is therefore specific: whether the definition may consume such
addressed calibration evidence, or requires a determination in the existing
type which preserves all information it is entitled to use. A direction
value or its marginal error bar alone does not supply that determination.

## Born-reconstruction continuation

[E6 and Born reconstruction](e6-born-bridge.md) audits the operational
connection to T70. The joint state has an exact commutative Born/Lüders
representation, but that representation admits any bounded positive
kernel. Positive mass rays do not resolve a POVM; real mass probes do not
identify a unique complex state. The missing bridge is an operationally
justified experiment and nuisance reduction, not a new chosen kernel.
Independent checks: `tests/test_born_operational_bridge.py`.

## Validation

The 32 E6 checks pass. The default full suite gave 2,691 passed, 28 skipped,
17 expected failures and five Windows sandbox pipe-creation denials. All five
affected worker tests passed when rerun outside the sandbox: 2,696 distinct
passing tests across those runs, with no numerical failure. The focused E6
suite also passed after the final covariance-tolerance edit. Slow opt-in tests
were not enabled. `tools/audit.py` retains its existing T6/T21/T27 and unused
L24 findings; `git diff --check` passes. Independent mathematical review
found no blocker in the stated results or their scope.

Continuation: 46 focused E6 checks now pass (14 new cases), including full
Gaussian density factorization, correlated-channel witnesses, full-span
checks, noisy direction calibration and independent noise-mixture quadrature.
The focused suite plus the full-type, residual-mixture and latent-measure
regressions pass: 125 tests. Independent review found no mathematical blocker.
The original full-suite result above belongs to the first exploration
commit; that full suite was not repeated for this research-only extension.

Calibrated-state continuation: 53 focused E6 checks pass (seven new cases).
They check joint Z, sequential combination, mass-neutral history, direction
reset, duplicated calibration, the composition witness, ambient evidence
and the explicit T30 failure. These plus the composition-bridge, short-law
and full-type regressions pass: 126 tests. Independent review found no
mathematical blocker; audit findings are unchanged. The full default suite
was not repeated for this isolated research extension.
