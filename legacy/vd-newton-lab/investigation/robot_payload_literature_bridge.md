# What the mass equation could contribute to robot payload identification

**Contribution 10 - 16 September 2026.** Primary-source literature review and
mechanical reduction, with independently checked examples. This is a working
research contribution, not a claim of robotic validation or priority.

## Finding

There is a legitimate connection, but the experiment must be specified more
carefully than "robots move slowly, so their acceleration denominator is near
zero."

**The clearest current test case is fixed-orientation translation, using a
stationary force reference measured with the same payload attached.** After
subtraction, the force change and acceleration satisfy our mass relationship.
Both can be partially resolved. The reference uncertainty then correlates all
subsequent measurements and must be retained.

The reviewed robotics work can supply important parts our equation assumes:
kinematic filtering, force-sensor bias handling, useful excitation, and the
physical domain for a full inertia model. Our existing equation supplies a
particular conditional mass law and continuous point readout. A practical
accuracy, speed, or uncertainty-calibration advantage remains to be shown.

## 1. What was actually read

Selection started from the STARS group's 2022 paper: its direct comparators,
physical-consistency foundations, follow-up work, and the particularly relevant
statistical-consistency literature. This is a targeted surrounding review,
not an exhaustive search of robot identification or a novelty determination.
Page references below use author-manuscript numbering unless indicated.

| Source | Reading coverage and role |
|---|---|
| [Nadeau, Giamou, Kelly, 2022: Fast Object Inertial Parameter Identification for Collaborative Robots](https://arxiv.org/pdf/2203.00830v3) | All eight pages, including figures, tables and observability appendix; main application paper. |
| [Kubus, Kröger, Wahl, 2008: On-Line Estimation of Inertial Parameters Using a Recursive Total Least-Squares Approach](https://www.researchgate.net/publication/224339569_On-Line_Estimation_of_Inertial_Parameters_Using_a_Recursive_Total_Least-Squares_Approach) | Author-uploaded full text: model, offsets, trajectory design, algorithms and experiments. Extracted equations inspected; PDF matrices not visually verified. |
| [Farsoni et al., 2018: Real-Time Identification of Robot Payload Using a Multirate Quaternion-Based Kalman Filter and Recursive Total Least-Squares](https://sfera.unife.it/handle/11392/2395584) | **Abstract only.** Institutional archive restricts both listed PDFs. DOI: 10.1109/ICRA.2018.8461167. |
| [Farsoni et al., 2017: Compensation of Load Dynamics for Admittance Controlled Interactive Industrial Robots Using a Quaternion-Based Kalman Filter](https://sfera.unife.it/bitstream/11392/2363687/1/RAL_kalquat_final.pdf) | Open eight-page manuscript, relevant filtering, compensation and experiment sections. A companion, not a substitute for the 2018 text. |
| [Wensing, Kim, Slotine: Linear Matrix Inequalities for Physically-Consistent Inertial Parameter Identification](https://arxiv.org/pdf/1701.04395v3) | Complete eight-page author manuscript; physical moment geometry. DOI: 10.1109/LRA.2017.2729659. |
| [Lee, Wensing, Park, 2020: Geometric Robot Dynamic Identification](https://www.researchgate.net/publication/337753338_Geometric_Robot_Dynamic_Identification_A_Convex_Programming_Approach) | Full author preprint available; principal formulation, covariance, regularization and appendix inspected; no experiment replication. DOI: 10.1109/TRO.2019.2926491. |
| [Janot, Wensing, 2021: Sequential Semidefinite Optimization for Physically and Statistically Consistent Robot Identification](https://www.researchgate.net/publication/346610534_Sequential_Semidefinite_Optimization_For_Physically_and_Statistically_Consistent_Robot_Identification) | Targeted full-text reading of statistical assumptions, consistency propositions and covariance analysis; not an audit of every proof. DOI: 10.1016/j.conengprac.2020.104699. |
| [Nadeau, Giamou, Kelly, 2023: The Sum of Its Parts: Visual Part Segmentation for Inertial Parameter Identification of Manipulated Objects](https://arxiv.org/abs/2302.06685) | Full seven-page paper, project and relevant public implementation; visual shape/part extension. |
| [Nadeau et al., 2024: Automated Continuous Force-Torque Sensor Bias Estimation](https://arxiv.org/abs/2403.01068) | Full seven-page technical report, including uncertainty propagation and Kalman equations. |
| [Nadeau, author-hosted thesis](https://starslab.ca/wp-content/papercite-data/pdf/2025_nadeau_autonomous.pdf) | Chapter 3, printed pp.27-54 / PDF pp.35-62, and relevant conclusions pp.120-122. The URL says 2025 but the downloaded cover says copyright 2026; this is the inspected version, not an asserted unchanged 2025 text. |

### What these sources establish

**2022 STARS paper.** It estimates mass, centre of mass and rotational inertia,
using point masses constrained by object shape. It already blends reduced
gravity-dominated dynamics and full dynamics smoothly, with a motion-based
weight (Section IV). This is a different transition from our noise-dependent
mass readout. Its appendix shows the reduced model cannot identify arbitrary
rotational inertia. Table IV's real experiments report lower mass error for
OLS than PMD, although PMD substantially improves other parameters. Therefore
its headline improvement should not be read as universal mass superiority.
Its noise, shape, motion and offset problems need to be separated before
matching our equation to them. [Full paper](https://arxiv.org/html/2203.00830v3)

**RTLS and filtering.** Kubus already models errors in both the regression
matrix and measured wrench, rather than only force error (Eq.23). Its
recursive SVD method uses diagonal weights (Eqs.24-36); offsets can be
corrected or estimated jointly (Section II). Farsoni 2018's abstract describes
multirate sensor fusion followed by RTLS, but its exact updates remain unread.
The accessible 2017 companion supplies filtering and compensation details
(Eqs.14-33) while explicitly assuming the inertial parameters known. Those
details cannot be attributed to the unavailable 2018 paper.
[Kubus](https://doi.org/10.1109/IROS.2008.4650672),
[Farsoni 2018](https://doi.org/10.1109/ICRA.2018.8461167),
[Farsoni 2017](https://doi.org/10.1109/LRA.2017.2651393)

**Physical geometry.** Wensing et al. characterize physically realizable
inertias through a positive moment matrix and supply convex shape constraints
(Definition 7, Theorems 3-4, pp.5-6). Its spatial covariance describes where
material lies; it is not sensor-error covariance. Lee et al. additionally use
measurement-residual covariance (Eq.3, p.2), explain a Gaussian likelihood
conditional on the observed regressor (p.4), and introduce geometric
regularization (pp.5-7). Thus uncertainty weighting is present, but it is not
our latent force/acceleration integration.
[Wensing](https://arxiv.org/pdf/1701.04395v3),
[Lee](https://doi.org/10.1109/TRO.2019.2926491)

**Statistical consistency.** Janot and Wensing address bias from noisy
kinematics and feedback, and derive parameter-covariance estimates under
asymptotic/linearized and residual-noise assumptions (Eqs.55-56). Their
convergence claim concerns identifiable base parameters
and depends on correct model form, sufficient excitation, error assumptions,
and an iteration-matrix convergence assumption whose formal justification is
deferred (III-E; Propositions 2-3). Physical feasibility, asymptotic consistency
and finite-data calibration are distinct achievements. This paper rules out
positioning ours as introducing uncertainty-aware robot identification.
[Author manuscript](https://www.researchgate.net/publication/346610534_Sequential_Semidefinite_Optimization_For_Physically_and_Statistically_Consistent_Robot_Identification)

**Later STARS work.** The 2023 follow-up estimates up to four homogeneous-part
masses from stopped measurements with sufficiently informative geometry and
poses; known part geometry then supplies inertia (Eqs.6-7). It reports benefits
for centre of mass and inertia, but not uniform mass-error superiority
(Table IV). The thesis chapter discusses coplanar-part degeneracy and future
uncertainty-guided motion selection (pp.53-54). These are valuable model and
experiment-design ideas; added shape/homogeneity information must be credited
in any comparison. [Follow-up](https://arxiv.org/abs/2302.06685),
[thesis](https://starslab.ca/wp-content/papercite-data/pdf/2025_nadeau_autonomous.pdf)

**Bias estimation.** The 2024 report propagates regressor uncertainty into
wrench-residual covariance and tracks sensor bias and drift with a Kalman
filter (Eqs.38-55). It assumes the relevant load/gripper inertial parameters
known. It does not establish simultaneous unknown-payload and bias
identification or provide an experimental calibration study. Applying its
known-load residual model to an unknown added payload would require a justified
calibration/update protocol or a joint model.
[Technical report](https://arxiv.org/abs/2403.01068)

## 2. Our mechanical reduction: which force and which acceleration?

The rest of this note derives consequences for our estimator. These are our
calculations, not claims that the reviewed authors proposed our equation.

Consider a rigid payload attached to a sensor at point S. Let c be the vector
from S to the centre of mass C, and express all vectors in the same axes at
the instant considered. Let g denote known uniform gravitational acceleration. Assume the
payload experiences only gravity and the measured attachment force; human
contact or other unmeasured forces would require additional terms. Calibrate
the sensor sign so F_s is the force exerted by the robot on the payload.
Any intervening gripper/tool load must also be removed or included in the
body whose mass is being estimated.

Rigid-body kinematics and Newton II give

$$
a_C=a_S+\dot\omega\times c+\omega\times(\omega\times c),
\qquad F_s=m(a_C-g). \tag{1}
$$

Equivalently, with h=mc and [v]_x the cross-product matrix,

$$
F_s=m(a_S-g)+\bigl([\dot\omega]_x+[\omega]_x^2\bigr)h. \tag{2}
$$

The original compatible-pair construction can be applied to (F_s,b), where
b=a_C-g, after specifying its joint measurement covariance and adopting its
flat positive-magnitude weighting for this experiment. The mechanical
reduction alone does not justify that weighting. Here the magnitude in the
denominator is beta=||b||. Equation (1) does not license pairing wrist force
with raw acceleration a_S.

### Stationary does not mean an unresolved mass denominator

At rest, b=-g and F_s=-mg. A 2 kg object gives 19.62 N and 9.81 m/s² in the
corresponding direction. As a deliberately incorrect-input demonstration,
using the known-direction formula h_sigma(F)/h_sigma(a) with raw a=0,
sigma_F=0.05 N and sigma_a=0.05 m/s² gives 491.8005 kg. Using b=9.81 instead
gives 2 kg. Continuity cannot repair a wrongly specified force balance.

An accelerometer's ideal specific-force channel represents a-g, whereas
acceleration obtained from positions represents a; calibration, location and
frame must be checked before either is used. If the sensor and centre of mass
are separated and the body rotates, Equation (2)'s extra terms matter.

For a checked example with m=2 kg, c=(0.2,0,0) m, a_S=0,
omega=(0,0,3) rad/s and dot(omega)=(0,2,0) rad/s², the contact force is
(-3.6,0,18.82) N. It is not parallel to a_S-g. Correcting gravity alone is
insufficient; the centre-of-mass acceleration gives the correct pair.

## 3. A concrete experiment that reaches our crossover regime

Hold the same payload at a fixed orientation, take a stationary reference,
then translate without rotation. Suppose sensor bias is constant over the
experiment. Writing its noisy reading as

$$
F_i^{obs}=m(a_i-g)+b_F+\epsilon_{F,i},
\qquad F_0^{obs}=-mg+b_F+\epsilon_{F,0},
$$

subtraction gives

$$
\Delta F_i^{obs}=m a_i+\epsilon_{F,i}-\epsilon_{F,0}. \tag{3}
$$

This is a **loaded stationary reference**, not merely an empty-gripper tare.
An empty-gripper correction leaves the payload's weight in the measurement.
If the reference acceleration is itself uncertain, retain it and use
Delta F=m Delta a instead of declaring it exactly zero.

Equation (3) gives an exact mechanical candidate for applying the mass
construction, including its partially resolved regime, under the stated
assumptions. Its probability law still has to be declared. It does not
establish that subtraction is preferable when calibrated weight already
provides better mass information. Its attraction is cancellation of a
constant unknown offset and of gravity at the fixed orientation.

### Shared reference noise is not fresh independent noise

Let z_i=(F_i,a_i) have independent raw measurement errors with covariance
Sigma_i. For differences d_i=z_i-z_0,

$$
\operatorname{Cov}(d_i,d_j)=\delta_{ij}\Sigma_i+\Sigma_0. \tag{4}
$$

The off-diagonal blocks are nonzero. With equal covariance Sigma for every
raw reading, averaging N differences has covariance (1+1/N)Sigma. Treating
the differences as independent would instead give 2Sigma/N. At N=20 that
understates variance by a factor of 10.5. An averaged reference reduces its
Sigma_0; it does not make its uncertainty independent between subsequent
samples. Equation (4) differences both channels; if only force is subtracted
and reference acceleration is treated as known, reference covariance enters
only the subtracted channels. Drift, filtering and correlated raw readings
require further terms.

The existing single-pair estimator can be tested on one properly modelled
difference pair. A sequence requires one shared mass, the shared reference,
and any latent accelerations in a joint construction. Averaging separately
fitted mass points is not that construction. Likewise, feeding our result
back into a robot fit as an independent prior would double-count evidence
if both used the same force measurements.

### Why bias creates a real weak-excitation problem

For this subsection only, take effective accelerations b_i as exactly known
and independent force noise with covariance sigma_F² I. Estimate m and an
unknown constant three-component sensor bias from F_i=m b_i+b_F+epsilon_i.
Use sensor coordinates if it is sensor-frame bias that is assumed constant;
it generally rotates when expressed in world coordinates.
The information for mass after accounting for the bias is

$$
\mathcal I_{m\mid b_F}
=\frac1{\sigma_F^2}\sum_i\|b_i-\bar b\|^2. \tag{5}
$$

Proof: stack the design rows [b_i,I_3] and take the Schur complement of its
bias block in X^T X/sigma_F². The full design has rank four exactly when at
least two b_i differ. With constant b, any mass change can be cancelled by
a bias change; gravity alone no longer separates them. Different static
orientations can change sensor-frame gravity and restore information even
without dynamic acceleration. Noise in b_i requires an errors-in-variables
extension; Equation (5) is not its general information formula.

This identifies an actual need our smooth point might help summarize: weak
changes in the force balance after unknown bias is accounted for. It cannot
create missing information at exactly constant b with unrestricted bias.

## 4. Which problems can each side help solve?

| Problem | What our present equation contributes | What robotics supplies / what remains |
|---|---|---|
| Partly resolved translational force and acceleration | A specified compatible-pair law and finite, positive, continuous readout | A candidate test, not proven better accuracy or calibrated intervals. |
| Very light payload with resolved gravity | A possible numerator-noise transition after sensor bias is controlled | Must distinguish small weight from offset/drift; direct weighing is a comparator. |
| Noisy regressors | Joint Gaussian pair covariance can represent declared channel errors | RTLS and simulated-dynamics identification already address noisy inputs; fit and averaging must be compared fairly. |
| Force offset, lag, multirate sensing | Assumed corrected inputs or an explicitly enlarged model | Filtering and bias models provide upstream machinery; shared uncertainty must survive preprocessing. |
| Centre of mass and rotational inertia | Not supplied by the present scalar equation | Shape/part models and physical moment constraints offer a route to an extension. |
| Unexcited directions or parameters | A finite point is a summary under a chosen measure | Useful motion, shape information or prior knowledge is still required for identification. |
| Choosing the reduced/full dynamics blend | No direct rule follows from our mass point | A noise-calibrated approximation-error calculation is a possible research project; scalar smoothness alone does not determine a blend weight. |

For the inertia limitation, two centred uniform spheres of mass 2 kg and
radii 0.05 m and 0.10 m have identical gravitational force and torque at every
static orientation, but principal inertias 0.002 and 0.008 kg m². A mass
readout cannot distinguish them from those data. Geometric knowledge could;
that would be additional information rather than recovery from the same
uninformative force channel.

## 5. A constructive route from the scalar readout to physical inertia

The moment geometry suggests one algebraically valid extension architecture.
It is not yet an estimator or a claim of novelty.

Let S be the raw second spatial mass moment, h=mc, and

$$
J=\begin{pmatrix}S&h\\h^T&m\end{pmatrix},
\qquad I_O=\operatorname{tr}(S)I_3-S. \tag{6}
$$

For a nonnegative mass distribution J is positive semidefinite; with m>0
its Schur complement S-hh^T/m is the central second moment. The strict
positive-definite interior excludes ideal plates, rods and points. Positive
rotational inertia alone is weaker: I_C=diag(1,1,3) implies the impossible
central second moment diag(1.5,1.5,-0.5).

**Our deduction.** Suppose a joint law over admissible J and a compatible
pair already exists, with f=m beta, beta>=0, 0<E[beta]<infinity, and
E[beta ||J||]<infinity. Finite E[f] alone does not guarantee those inertia
moments. Define

$$
\widehat J=\frac{E[\beta J]}{E[\beta]}. \tag{7}
$$

Its mass entry is E[f]/E[beta], exactly the original readout for that pair.
For every v,

$$
v^T\widehat Jv=\frac{E[\beta v^TJv]}{E[\beta]}\ge0.
$$

It therefore preserves physical moment feasibility. Fixed linear shape
constraints survive the same positive averaging, and a rigid frame change
commutes with it by matrix congruence. Addition commutes **only with the same
weights and common frame**; separately fitted parts need not share those
weights. Rotating body parts need not even share acceleration magnitudes.

This gives a useful compatibility result between the two research threads.
It does not construct the required probability law, select its reference
measure, establish finite moments for that law, or justify this weighting for
every rotational task. With attachment force, beta must be ||a_C-g||; with
the loaded-reference experiment it is the relevant change magnitude. Torque
and angular acceleration are generally not parallel, so replacing scalar
mass by a tensor in the original parallel-vector parameterization would not
suffice.

## 6. The next experiment has a precise contract

The most informative first comparison is a small translational benchmark,
before attempting all ten inertial parameters:

1. **Measurements:** the same force and kinematic streams for every method;
   known payload mass for evaluation; a separately specified reference period;
   fixed orientation and a checked absence of other contact forces.
2. **Model:** document which channels are measured, their covariance and
   filtering, baseline uncertainty, bias/drift and one common mass. Vary
   signal strength across the previously identified crossover region.
3. **Comparators:** calibrated weight when available, a covariance-aware
   errors-in-variables fit, and positive regularized alternatives with any
   prior information made explicit. Compare a single-pair experiment first;
   any sequence version must derive its shared-mass law.
4. **Outcomes:** mass error under declared loss, failure/extreme-output rate,
   interval coverage and width, and computation time. Smoothness and finite
   outputs are properties to record, not substitutes for these outcomes.
5. **Stress cases:** incorrect reference covariance, modest drift, correlated
   channels, and small unintended rotations. Report model sensitivity
   separately from nominal-model performance.

The earlier crossover simulations are relevant preparation; they are not
robot experiments and do not meet this contract yet. An eventual full robot
comparison should also assess centre-of-mass/inertia estimates and downstream
prediction or control, using equivalent shape information.

## 7. Reproduction, access and provenance

[robot_payload_checks.py](robot_payload_checks.py) reproduces seven groups of
mechanics and algebra examples: stationary gravity, rotational correction,
bias identifiability/information, shared-reference covariance, inertia
non-identifiability, an impossible positive inertia, and the weighted moment
readout. All assertions passed on 16 September using Python and NumPy. These
checks validate the displayed examples, not a published solver or robot
performance. The general claims rest on the derivations above.

~~~powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/robot_payload_checks.py'
~~~

The 2022 [project](https://papers.starslab.ca/fast-inertial-identification/)
links publicly to [test-object model generation](https://github.com/utiasSTARS/pyb-sim-models/tree/main/pbsm/models/CompositeTestObject).
That inspected link is not evidence that its full solver and raw benchmark
streams are released. The 2023 [project](https://papers.starslab.ca/part-segmentation-for-inertial-identification/)
links a public [implementation and dataset repository](https://github.com/utiasSTARS/inertial-identification-with-part-segmentation).
Its relevant READMEs and identification code were inspected via public main
URLs; this pass did not run the pipeline or verify every dataset download.

**Open access sufficed for the central review.** Farsoni 2018 is a specific
remaining full-text gap for which authorized Imperial access could help.
Kubus's extracted equations need PDF-level visual verification before an
implementation-level reproduction. No author contact or institutional
sign-in was performed. No source paper's reported performance was independently
replicated.

Managed context uses
[VD-docs commit 8782a293297330f7a354d87ebe62792018fa8866](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866),
selected through catalog.yaml and the mass-discussion README in this ongoing
investigation. The pertinent canonical inputs are the
[working specification](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md),
[composition theorem](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_readout_composition_uniqueness_2026-09-06.md),
and [law-assumed-correct steering](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md).
The local investigation, especially contributions 4, 5 and 9, supplies active
analysis of common-mass inference and the crossover behavior. Managed files
were not edited; historical literature-review files were not treated as
discussion authority.

Research notes, downloaded source material and machine output remain under
`.tools/robot_pmd_20260916/`. The thesis download's SHA-256 is
`685ad6a20b724e0aeeb0e6807e4b41f8accd9f844a8687160f939c0df7ab5ce5`.
The public-paper search used exact titles/DOIs, the 2022 bibliography and
STARS project/author links; this was not an exhaustive managed-inbox or
unsubmitted-local-material audit.
