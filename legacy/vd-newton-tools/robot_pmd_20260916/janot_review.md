# Janot–Wensing 2021: consistency and uncertainty are already robotics topics

**16 September 2026.** Bounded primary-source review supporting the robotics
comparison. No benchmark replication or production change.

## Source and reading coverage

Alexandre Janot and Patrick M. Wensing, *Sequential Semidefinite Optimization
For Physically and Statistically Consistent Robot Identification*, Control
Engineering Practice 107 (2021), article 104699.

- DOI: https://doi.org/10.1016/j.conengprac.2020.104699
- Author-uploaded article, attributed to Alexandre Janot, uploaded 3 December
  2020:
  https://www.researchgate.net/publication/346610534_Sequential_Semidefinite_Optimization_For_Physically_and_Statistically_Consistent_Robot_Identification

**Coverage:** The author-uploaded extracted full text was accessible. Inspected
Sections II-C/D, III-E/F, IV-B/C/D/E and the covariance Monte Carlo discussion
in V-E; checked the listed definitions, propositions and equation references.
This is a targeted full-text review, not an independent proof audit of all
appendices or visual inspection of every table/matrix. Page references below
are the author manuscript's printed pages.

## Factual synopsis — under 170 words

PC-DIDIM combines simulated closed-loop dynamics with physical constraints.
Statistical consistency means convergence in probability of identifiable
base-parameter estimates (Definition 1, p. 3), not exact finite-sample accuracy.
Noisy kinematics enter regressors nonlinearly; feedback correlates regressors
and residuals. Physical constraints alone do not remove that bias:
PC-IDIM-LS requires $E[X^T\varepsilon]=0$ (Proposition 2, (44), p. 7).
PC-DIDIM uses simulated kinematics. Its consistency proposition requires no
modeling error, giving $X_s=X_{nf}$, together with zero-mean finite-variance
errors (Proposition 3, (52), p. 8). The analysis assumes sufficient excitation
and convergence of an iteration matrix; a formal justification of the latter
is deferred (III-E, p. 6). The residual model is Gaussian with serially
uncorrelated, joint-dependent variances (7–8); whiteness is checked (IV-C).
Parameter covariance is explicitly addressed: (55) incorporates constraint
effects, (56) transforms to base parameters, and (53) is a rough alternative
ignoring constraints (pp. 9–10). Monte Carlo covariance comparisons are
reported in V-E. The work is therefore neither uncertainty-free nor based on
an assumption that measured kinematics are exact.

## Our comparison and interpretation

The following points are our deductions for the present investigation, not
additional claims reported by Janot and Wensing.

### 1. Distinguish three kinds of correctness

For our positive-mass construction, these questions must remain separate:

1. **Physical admissibility:** does an output satisfy the imposed physical
   restrictions, such as positive mass?
2. **Statistical consistency:** does a sequence of estimates approach the
   identified true quantity as suitable data accumulate?
3. **Finite-data performance:** how accurate, variable and well calibrated is
   the method on the actual experiment?

A finite positive point at measured zero-zero answers part of the first
question. It cannot by itself answer the second or third. Conversely, an
asymptotic consistency result does not promise a useful finite-data mass
determination in a nearly unexcited experiment.

### 2. Measurement covariance and parameter uncertainty are different objects

Our Gaussian covariance $\Sigma$ describes the supplied force/acceleration
measurement errors. The induced law of $M=f/\alpha$ follows only after
choosing the compatibility model and reference measure. The sampling
distribution of the reported point $E[f]/E[\alpha]$ is another object again.

It would therefore be misleading to compare our **input** measurement
covariance with another method's **output** parameter covariance and conclude
that one supplies uncertainty while the other does not. A valid comparison
must specify whether it targets conditional mass probabilities, repeated-data
coverage, parameter-estimate covariance, or downstream prediction error.

### 3. Known Gaussian covariance does not automatically solve feedback

If a controller changes force in response to a noisy kinematic reading, the
resulting force and acceleration evidence may share disturbances. A correct
known joint covariance can represent certain contemporaneous linear Gaussian
dependencies. It does not automatically capture unknown controller errors,
time dependence, nonlinear transformations, filtering bias or model mismatch.

Accordingly, our framework should state whether it conditions on an imposed
input, models a feedback-generated input, or receives corrected measurements
from a separate model. Calling both observed channels uncertain is not a
substitute for this step.

### 4. Identifiable quantities must match

A mechanism can constrain combinations of inertial parameters while leaving
individual link parameters ambiguous. Our scalar mass target assumes the
body, net force and relevant acceleration have already been identified.
Comparing the methods requires matching that target and the nuisance
parameters; a scalar mass ratio is not automatically a substitute for a
complete robot identification problem.

### 5. What this rules out in our positioning

We should not describe the surrounding field as lacking:

- procedures that account for noisy kinematics;
- methods that address feedback-induced estimation bias;
- physically constrained identification;
- parameter covariance outputs or statistical model checks.

Our distinctive candidate is the **particular compatibility measure and
ratio-of-magnitude-means readout**, with its established algebraic properties.
Whether this yields a practical advantage under a legitimate robot
measurement model remains a separate question. This review supplies no
evidence of superiority for our construction.

### 6. A precise next comparison

For any proposed robot use, document the observed quantities, controller,
calibration and nuisance parameters; then derive the common measurement
model. Compare estimates for the same identifiable target under the same
excitation and error process. Evaluate repeated-data error and uncertainty
calibration separately from physical feasibility. If our method is meant to
remain useful near a failure of excitation, define what useful information
means there—such as a one-sided bound—before comparing finite point values.

The conditional assumptions in a competitor's guarantee are not a defect to
be used as evidence for ours. They specify a standard of explicitness that
our own claims must meet.
