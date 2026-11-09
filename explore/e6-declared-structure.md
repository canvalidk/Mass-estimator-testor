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

## Validation

The 32 E6 checks pass. The default full suite gave 2,691 passed, 28 skipped,
17 expected failures and five Windows sandbox pipe-creation denials. All five
affected worker tests passed when rerun outside the sandbox: 2,696 distinct
passing tests across those runs, with no numerical failure. The focused E6
suite also passed after the final covariance-tolerance edit. Slow opt-in tests
were not enabled. `tools/audit.py` retains its existing T6/T21/T27 and unused
L24 findings; `git diff --check` passes. Independent mathematical review
found no blocker in the stated results or their scope.
