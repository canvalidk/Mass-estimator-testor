# STARS follow-up: shape assumptions, low motion, and bias uncertainty

Independent research contribution, 16 September 2026. External sources read
only; no production edits, code execution from the repositories, or outbound
communications. This is a bounded follow-up review, not a field-wide novelty
or completeness claim.

## 1. Reading coverage and source identity

| Source | Coverage actually read | Precise locators |
|---|---|---|
| Nadeau, Giamou, Kelly, *The Sum of Its Parts* (2023) | Full seven-page arXiv v2 text, including tables, captions and references | [PDF](https://arxiv.org/pdf/2302.06685v2): p1 Assumption 1; p3 Eqs.1-7 and footnote 2; pp5-6 Tables III-IV; p6 Sections VI-VII |
| [HPS project](https://papers.starslab.ca/part-segmentation-for-inertial-identification/) | Entire project page and its repository links | Dataset and Source Code sections |
| [Public HPS repository](https://github.com/utiasSTARS/inertial-identification-with-part-segmentation) | Root README, identification README, dataset README, HPS.py, InertialParametersEstimation.py; directory/commit listings | [HPS implementation](https://github.com/utiasSTARS/inertial-identification-with-part-segmentation/blob/main/inertial-parameters-identification/HPS.py); [data processing](https://github.com/utiasSTARS/inertial-identification-with-part-segmentation/blob/main/inertial-parameters-identification/InertialParametersEstimation.py) |
| Nadeau, *Autonomous Object Handling in Collaborative Robotics* | Chapter 3 in full: printed pp27-54 = PDF pp35-62; relevant concluding material printed pp120-122 = PDF pp128-130; title/contents and limited introductory provenance | [Supplied author PDF](https://starslab.ca/wp-content/papercite-data/pdf/2025_nadeau_autonomous.pdf): pp40-41; pp44-45 Eqs.3.18-3.26; pp51-54; pp121-122 |
| Nadeau, Rogel Garcia, Wise, Kelly, *Automated Continuous Force-Torque Sensor Bias Estimation* | Full seven-page technical report, including references; revision 1.2, March 5, 2024 | [PDF](https://arxiv.org/pdf/2403.01068): pp1-3 Sections1.1-1.3; pp3-5 Section2; pp5-7 Section3, Eqs.31-55 |

**No full-thesis reading is claimed.** Thesis chapter text is saved as
`thesis_chapter3_text.txt`; equation page45 and future-work page54 were also
rendered and visually checked. The thesis URL/website bibliography says 2025,
but the downloaded cover says copyright 2026. Preserve that discrepancy when
citing this version. Local PDF SHA-256:
`685ad6a20b724e0aeeb0e6807e4b41f8accd9f844a8687160f939c0df7ab5ce5`.

The GitHub history listing showed latest commit
`48928088e4877c872ea5aaed92f2ef315190153a` (February 6, 2023); reads used public
`main` URLs. The commit detail fetch failed, so this is a listing observation,
not a fresh Git/API resolution of all source blobs at that commit.

## 2. Factual synopses

### HPS paper (under 180 words)

HPS adds segmented object geometry and homogeneous density within each part.
Known part geometry makes every inertial parameter a function of part masses.
For stopped measurements, it estimates up to four nonnegative masses by
least squares; nondegenerate geometry/poses remain necessary. The paper reports
20 tools across four noise settings and a hammer-balancing demonstration.
Table IV shows gains principally for centre of mass and inertia under noise,
not mass alone: high-noise mass errors are HPS 2.79%, OLS 0.69%, GEO 0.64%.
GEO receives true mass in its prior. HPS selects low-acceleration samples;
OLS/GEO use all samples. The physical demonstration uses the estimated centre
of mass. Discussion identifies incorrect heterogeneous-part merging and
coplanar part centroids as limitations; future work proposes material-density
choices and improved segmentation. The reported evaluations are point-error
and physical-consistency measures, not calibrated parameter-interval coverage.

### Author-hosted thesis (under 180 words)

Chapter 3 retains the separate PMD and HPS approaches. PMD's motion weighting
uses a robot-tuned parameter; pp40-41 propose selecting it using noise-power
estimates. HPS discussion on p53 says geometry errors can affect inertia more
strongly than mass or centre of mass, and retains the coplanar-centroid issue.
Chapter future work on p54 proposes optimizing trajectories for shape and
inertial uncertainty, continuously refining models during handling, and using
visual material recognition for density priors. The general future-work
section on pp121-122 again proposes integrating identification into ordinary
task motion. The reviewed sections do not supply a calibrated posterior or
frequentist coverage analysis for the identified inertial parameters.

### HPS project page (under 180 words)

The page links a public repository for segmentation, identification, and data.
It describes 20 CAD-based workshop tools with watertight meshes, synthetic
coloured point clouds, ground-truth inertial parameters, part labels, and
grasping frames. It also links a demonstration video. The page's statement
about moving the repository after review is stale relative to its live GitHub
links; availability was checked directly rather than inferred from that text.

### Repository inspection (under 180 words)

The root identifies synthetic stop-and-go dynamics data and segmentation
outputs. Dataset README defines CAD-derived ground truth and material notes;
this is not a collection of independently calibrated physical-tool measurements.
Identification README provides HPS, GEO and OLS reproduction instructions,
requiring Python, NumPy, Open3D, SpatialMath, CVXPY, MOSEK, ROS and
With-Respect-To; it describes 240 algorithm executions. HPS.py enforces at most
four parts, computes volume-weighted centroids, and solves nonnegative squared
residual minimization using MOSEK. It returns point inertial parameters.
InertialParametersEstimation.py filters HPS updates by both measured linear
and angular acceleration norms being below 1 and constructs GEO's prior using
the true mass. MIT licensing is shown. Source availability is verified;
installation, numerical reproduction, downloadable completeness and release
quality are not certified by this read-only inspection.

### 2024 bias report (under 140 words)

The report estimates joint kinematics with a jerk-noise Kalman model, maps
them to sensor kinematics, then filters wrench bias and its drift. It assumes
known load inertial parameters; the abstract specifically says known gripper
parameters with certainty. Section3 includes regressor uncertainty:
`B=theta^T kron I`, with measurement covariance
`R=B Sigma_D B^T+Sigma_W` (p6, Eqs.38-49). Bias/drift state covariance and
time-extrapolated wrench correction follow on p7, Eqs.50-55. Thus covariance
handling and ongoing bias estimation are already part of this group's work.
The seven-page report contains a formulation rather than an experimental
validation or interval-coverage study. It does not present joint estimation
of unknown payload inertia and sensor bias.

## 3. Our interpretation and comparisons

### The low-motion gap is conditionally addressed, not simply left open

HPS must be included in any claim about what remains after PMD2022. It uses
additional information to reduce the number of independent unknowns. A method
that knows a shape and homogeneous partition is solving a more constrained
problem than a method given only noisy force and acceleration.

Here is the structural reason, independently derived from the static moment
relations. Let part centroids be `r_j`, and write

`[total_mass; first_moment] = B_part m_parts`,
`B_part = [1 ... 1; r_1 ... r_N]`.

Even perfect static wrench observations constrain at most these four moment
components. Full column rank of `B_part`, in addition to sufficiently varied
gravity poses, can identify up to four part masses. Four coplanar centroids
make this matrix rank-deficient. A strictly positive solution then normally
has nearby equally fitting mass allocations along its nullspace. A scalar
mass estimator cannot create missing second-moment information or repair this
structural degeneracy. Positivity and homogeneity are different assumptions.

### Sensor force is not automatically net force

Our mass equation needs the physical relation behind its inputs established.
For a supported object, wrist force contains gravity and tool contributions.
Ignoring rotation for illustration, `F_support=m(a_COM-g)` up to the sensor's
sign convention. At rest the support-force magnitude is approximately `m|g|`,
although kinematic acceleration is zero. Pairing that support force directly
with zero kinematic acceleration would pose the wrong Newton-II inference.

One can instead use the appropriate effective acceleration, or model all
external forces jointly. A controlled difference of fixed-orientation loaded
readings can remove constant gravity contributions, but tool inertia, bias,
orientation changes and the stochastic baseline still require treatment.
Rotating experiments additionally couple acceleration to the unknown centre
of mass and inertia tensor. HPS is therefore a relevant practical baseline,
not a drop-in scalar comparator on unprocessed wrist data.

### What remains a credible contribution target

Our current theory could support an errors-in-both-channels treatment of a
carefully isolated translational mass subproblem and a useful uncertainty
readout. It has not demonstrated improved full inertial identification, nor
does it supersede HPS's geometry information. A convincing comparison should
hold the available information fixed, distinguish mass from inertia-tensor
error, include realistic correlated uncertainty, and test prediction on data
not used for identification. It should include bias/calibration failure modes
and misspecified shape/homogeneity when those priors are used.

The absence of a parameter-coverage analysis in the inspected HPS material is
a bounded opening for evaluation, not evidence that the entire literature
lacks uncertainty methods. The 2024 report makes that caution concrete.

### Shared taring and the 2024 bias model

The bias report's known-inertia premise matters. If an unmodelled payload
contributes to the wrench, subtracting only a known-gripper prediction leaves
both payload signal and bias in the residual. Letting a bias filter absorb
that entire residual can remove part of the signal we intend to use. An
unknown-payload experiment therefore needs a justified calibration/update
protocol or a joint model; the known-gripper report is not automatically that
joint solution.

For a simple independent-noise example, define `Delta W_i=W_i-W_0` using one
shared reference measurement. Then

`Cov(Delta W_i,Delta W_j)=Sigma_0` for distinct i,j,
`Var(Delta W_i)=Sigma_i+Sigma_0`.

These are our direct covariance calculations. Reusing the reference does not
create independent force observations. Drift adds temporal structure; using
the same measurements to estimate the correction introduces cross-covariance
terms too. If corrections use a shared bias/drift state with covariance `P_b`
and propagation matrices `H_i`, an independent calibration-state model adds
`H_i P_b H_j^T` between corrected readings. The 2024 machinery is a sensible
starting point, while those dependence and unknown-payload questions still
need to be handled for our intended experiment.

## 4. Completion limits

Read-only source and data availability checks are complete. No repository was
executed or installed, no full dataset was downloaded, and no authors were
contacted. The thesis download succeeded after network escalation; local
working materials remain in this scratch directory. All proposed extensions
above are our interpretations, not claims attributed to the source authors.
