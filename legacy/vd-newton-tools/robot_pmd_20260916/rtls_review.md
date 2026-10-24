# RTLS and robot payload identification: targeted primary-source review

**16 September 2026.** Scratch contribution to the mass-estimator investigation.
Scope: Kubus–Kröger–Wahl 2008, Farsoni et al. 2018, and the accessible 2017
Farsoni companion needed to understand its filtering context. Paper synopses
are deliberately short; the subsequent mathematical comparison is our own
analysis and does not assert new results established by these papers.

## 1. Kubus, Kröger and Wahl (2008)

**Title:** *On-Line Estimation of Inertial Parameters Using a Recursive Total
Least-Squares Approach*. IROS 2008, pp. 3845–3852.

**DOI:** https://doi.org/10.1109/IROS.2008.4650672

**Primary-source routes:**

- Institutional metadata: https://publikationen.bibliothek.kit.edu/1000083705
- Author-uploaded full text, attributed to Daniel Kubus:
  https://www.researchgate.net/publication/224339569_On-Line_Estimation_of_Inertial_Parameters_Using_a_Recursive_Total_Least-Squares_Approach
- Author's former laboratory bibliography:
  https://www.rob.cs.tu-bs.de/www.rob.cs.tu-bs.de/bibliob407.html?f%5Bauthor%5D=232&page=8

**Reading coverage:** Author-uploaded article's extracted full-text material:
model and offsets; trajectory design; RLS/RIV/RTLS sections; update formulas;
experimental description, results discussion and conclusion. Equation references
below were checked against the extracted text. The PDF download endpoint failed;
this is not a visual verification of every matrix element or figure ordinate.

### Factual synopsis — 160 words

The paper estimates ten rigid-load inertial parameters from wrist wrench and
kinematics. Newton–Euler equations (1–5), pp. 3846–3847, include gravity,
rotation and centre-of-mass terms. Sensor offsets are handled by known-orientation
zeroing with gravity correction (7–9), or by expanding to sixteen estimated
parameters. Drift is neglected over the short experiment. Trajectories are
optimised for conditioning and limited mechanical excitation (Section III).
The RTLS error model corrects both the data matrix and wrench (23), unlike
the stated RLS model (16). Incremental SVD with diagonal weighting matrices
$D,T$ supplies recursive estimates (24–36), pp. 3849–3850. Equation (36)
requires a nonzero singular-vector denominator. The inspected RTLS presentation
does not specify a Gaussian probability law or a calibrated full residual
covariance. Trials used a Stäubli RX60, JR3 sensing, 1000 Hz estimation and
1.5-second duration. RTLS outperformed the tested RLS/RIV configurations;
combining measured linear acceleration with filtered angular acceleration
gave the best signal-source results in that setup (Section V). These are
experimental comparisons, not universal statistical dominance claims.

## 2. Farsoni et al. (2018)

**Title:** *Real-Time Identification of Robot Payload Using a Multirate
Quaternion-Based Kalman Filter and Recursive Total Least-Squares*.
ICRA 2018, pp. 2103–2109.

**DOI:** https://doi.org/10.1109/ICRA.2018.8461167

**Primary-source routes and access findings:**

- University of Modena/Reggio Emilia metadata and abstract:
  https://iris.unimore.it/handle/11380/1174438
- University of Ferrara metadata and abstract:
  https://sfera.unife.it/handle/11392/2395584
  This archive explicitly lists both `icra2018.pdf` (preprint, 1.33 MB) and
  `ICRA.2018.8461167 (1).pdf` (publisher text, 2.47 MB) as private/restricted.
- Publisher: https://ieeexplore.ieee.org/document/8461167
  No usable article text returned by this tool.
- ResearchGate:
  https://www.researchgate.net/publication/327807970_Real-Time_Identification_of_Robot_Payload_Using_a_Multirate_Quaternion-Based_Kalman_Filter_and_Recursive_Total_Least-Squares
  Explicitly says no full text available. No copy request was sent.

**Reading coverage:** Metadata and author/institution abstract only. Searches
by exact title, DOI, authors, `icra2018.pdf`, and title plus `filetype:pdf`
found no legitimate open full text. Institutional access remains a useful
next route. The archived files' access restriction was respected.

### Factual synopsis — 80 words

The institutional abstract describes real-time reconstruction of a rigid
end-effector load's inertial parameters. A multirate quaternion Kalman filter
combines robot-kinematic and inertial-sensor measurements, potentially sampled
at different frequencies, to estimate linear acceleration and angular
velocity/acceleration. RTLS then identifies load parameters without offline
post-processing. The accessible abstract does not provide covariance choices,
offset handling, update equations, experimental error figures or convergence
times. Those details remain unverified; neither later papers' summaries nor
the 2017 companion can substitute for the 2018 full text.

## 3. Essential accessible companion: Farsoni et al. (2017)

**Title:** *Compensation of Load Dynamics for Admittance Controlled Interactive
Industrial Robots Using a Quaternion-Based Kalman Filter*.
IEEE Robotics and Automation Letters 2(2), pp. 672–679.

**DOI:** https://doi.org/10.1109/LRA.2017.2651393

**Institutional metadata:** https://sfera.unife.it/handle/11392/2363687

**Open institutional author manuscript:**
https://sfera.unife.it/bitstream/11392/2363687/1/RAL_kalquat_final.pdf

**Reading coverage:** Eight-page author manuscript available. Inspected
Sections II–V covering kinematics, Kalman/EKF construction, compensation,
offsets, thresholds, covariance tuning and experiments. Page numbers below
are manuscript pages 1–8, not journal pagination. No details here are silently
attributed to the unavailable 2018 paper.

### Factual synopsis — 153 words

The companion estimates load kinematics to remove non-contact effects from
wrist force/torque signals. It explicitly assumes load mass, centre of mass
and inertia are known (Section IV, p. 5). A linear Kalman filter estimates
translation derivatives; an EKF handles quaternion orientation and angular
derivatives (14–27), pp. 4–5. Covariance matrices $Q,R$ are introduced in (14)
and tuned for response/smoothness; numerical diagonal choices and the 500 Hz
experimental rate appear on p. 6. Newton–Euler compensation includes gravity
and rotation (28–29). Offset zeroing and initial-orientation correction enter
(30), with adaptive thresholds (31–33), pp. 5–6. Experiments report mean
filter computation time of 408 microseconds and acceptable compensation
quality at sample periods up to 10 milliseconds. Differentiation plus a
20 Hz FIR filter produced undesirable lag in their comparison; the tuned
Kalman output had negligible lag. Residual oscillations are associated with
possible structural flexibility (p. 7). This supports the importance of
filtering and model error, not mass-estimation superiority.

## 4. Our mathematical comparison: what problem is actually shared?

The present estimator uses, for one candidate force–acceleration pair,

$$dP\propto\exp[-Q(f,\alpha,u)/2]\,df\,d\alpha\,d\Omega,
\qquad f,\alpha>0,\qquad \widehat m=E_P[f]/E_P[\alpha].$$

Its uncertainty model and its point readout are distinct ingredients. RTLS
is a method for fitting a relationship when both the regression matrix and
observed response contain errors. Accounting for uncertainty on both sides
is therefore shared territory, not an advantage established merely by naming
our two noisy channels.

### 4.1 Translational reduction and the gravity trap

Our reduction of the Newton–Euler model, with rotational contributions absent,
is

$$F_s=m(a_s-g).$$

Here $F_s$ is the offset-corrected sensor force; it is not automatically the
net force used in $F_{\rm net}=ma_s$. Defining $b=a_s-g$ gives a mass ratio
from $F_s$ and $b$. A stationary load with $a_s\simeq0$ still has
$b\simeq-g$: it need not be a weak-regressor mass experiment.

If the sensor was tared at a known stationary initial orientation, the
corresponding translational relation becomes

$$F_{\rm tare}=m(a_s-g+g_{\rm init}).$$

For fixed orientation, $g=g_{\rm init}$, so the relation reduces to
$F_{\rm tare}=ma_s$. **A tared, fixed-orientation translational experiment is
therefore a clean candidate for testing our weak-acceleration scalar method.**
That is our experimental deduction, not a result asserted by either paper.
The common initial tare uncertainty must still be propagated; subtraction
can correlate every subsequent sample through the same baseline error.

With rotation and unknown centre of mass, isolating the acceleration at the
centre of mass requires additional unknown geometry. A single-pair mass
formula cannot replace simultaneous inertial-parameter identification without
an explicit treatment of those nuisance parameters. More output parameters
alone do not make methods incomparable: mass is a component of the broader
problem, but the data reduction must be legitimate.

### 4.2 A like-for-like errors-in-variables comparator

For independent Gaussian scalar measurements $x_i,y_i$ with known variances,
an alternative *fit* under the same law minimises over common mass and latent
accelerations:

$$\sum_i\left[(x_i-ma_i)^2/s_{F,i}^2
 +(y_i-a_i)^2/s_{a,i}^2\right].$$

For a fixed mass, when the unconstrained minimiser has $a_i>0$, profiling out
$a_i$ leaves

$$D(m)=\sum_i\frac{(x_i-my_i)^2}{s_{F,i}^2+m^2s_{a,i}^2}.$$

This algebra follows by completing the square. If the minimiser lies outside
the positive domain, the appropriate constrained boundary solution is needed.
The vector analogue uses the residual covariance
$\Sigma_{FF}+m^2\Sigma_{aa}-m(\Sigma_{Fa}+\Sigma_{aF})$.

That fit and our posterior/reference-measure averaging generally return
different points. Our construction integrates candidate pairs using its
declared measure; the fit selects corrections of least weighted squared size.
Neither operation is validated by confusing it with the other. A claim of
superior use requires comparison on common inputs, shared nuisance parameters
and declared loss/coverage criteria.

### 4.3 A filter covariance is not automatically the covariance we need

Suppose a kinematic filter produces an acceleration estimate. Feeding that
estimate to the mass equation requires its uncertainty jointly with the
force estimate, plus any relevant cross-time correlations. Filtering raw
independent sensor readings can create temporal dependence. Frame rotations,
products of angular rates, gravity subtraction and common tare corrections
can create cross-channel dependence or non-Gaussian errors.

Therefore our ability to accept a general Gaussian covariance matrix is
useful only after that matrix and Gaussian approximation have been justified.
It does not itself estimate unknown offsets, align multirate streams, choose
filter bandwidth, remove lag, or model structural flexibility. Conversely,
the presence of diagonal weights in a fitting method does not establish that
all possible uncertainty models have already been solved there.

## 5. What has been established about our added value?

**Established from our construction:** for its specified Gaussian pair model
and flat positive magnitude/common-direction measure, the readout supplies a
positive, finite conditional point, with the continuity and composition
properties previously derived in the local investigation. Those are
mathematical properties of that defined operation.

**Not established by this literature comparison:** better robotic mass
accuracy, faster convergence, better calibrated confidence, lower control
error, or a unique useful response to poor excitation. A robot's need for
noise-aware mass estimation is established surrounding territory. Whether
our particular averaging rule adds a useful capability remains open.

The practical contribution to test is consequently specific: insert the
declared probability construction and readout into a legitimate translational
data model; propagate filtering/tare/gravity uncertainty; compare against
covariance-aware fitting on the same time series. Maintain the fixed common
mass across samples. Do not replace that joint problem with independent
per-sample estimates and assume averaging them preserves its information.

## 6. Remaining evidence needed

1. Read the actual Farsoni 2018 full text through authorized institutional
   access, checking the multirate update, RTLS weighting, offsets, covariance
   propagation and experimental trajectories.
2. Verify the exact matrix weighting/normalisation of Kubus 2008 visually
   from a usable PDF before reproducing its algorithm.
3. Choose the simplest genuine robot experiment whose measured quantities
   satisfy our model; fixed-orientation tared translation is one candidate.
4. Compare shared-model joint inference before claiming a practical advantage.

No contact with authors, library, publisher or other people was made. No
restricted PDF was obtained and no production file was modified.
