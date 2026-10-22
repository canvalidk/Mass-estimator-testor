# Arrow drag candidate: data audit

Date: 2026-09-06. Status: local record awaiting transfer to VD-docs; not
submitted. Outcome: the methods are suitable in principle, but no admissible
measured force/acceleration pair with characterized uncertainty was obtained.
No mass-estimator run or new mass result is claimed for this candidate.

## Sources inspected

1. Okawa et al., *Free flight and wind tunnel measurements of the drag exerted
   on an archery arrow*, Procedia Engineering 60 (2013), 67-72.
   [DOI](https://doi.org/10.1016/j.proeng.2013.07.017).
   [Full paper inspected](https://archerycoachesguildblog2.wordpress.com/wp-content/uploads/2015/06/tl_arfd4.pdf).
   Read sections 2-4 and references, particularly the magnetic force measurement
   and two-camera flight analysis. Numerical output is mainly plotted drag
   coefficients; no per-shot camera table and matching force covariance found.
2. Miyazaki et al., *Aerodynamic properties of an archery arrow*, Sports
   Engineering 16 (2013), 43-54, DOI 10.1007/s12283-012-0102-y.
   [Publisher page and accessible appendix](https://link.springer.com/article/10.1007/s12283-012-0102-y).
   [Author-manuscript copy inspected](https://archerycoachesguildblog2.wordpress.com/wp-content/uploads/2015/06/tl_arfd2.pdf).
   Read methods, results, interpretation, and appendix. The manuscript is not
   asserted to be identical to the final version. The main camera-decay equations
   were checked against the final publisher's accessible appendix.

Searched both DOIs with data/supplementary/dataset and the companion title with
data/repository. No downloadable raw dataset or supplement was located in this
bounded search. This is not proof that none exists. Web PDF screenshots failed;
no graph digitization was performed. A related 2020 paper's landing page returned
HTTP 429 and was not used as evidence.

## Finding

These experiments offer independently measured aerodynamic force and motion in
principle. Published coefficients and qualitative agreement are insufficient to
reconstruct a specific matching pair and its joint uncertainty reliably. The
two physical setups also need condition matching: attitude, rotation, vibration,
configuration, air properties, and flow regime.

The companion manuscript attributes some large scatter to laminar/turbulent
switching rather than increased measurement error. Therefore a broad Gaussian
error bar around one fixed latent force would need justification. Pairing a
force from one flow regime with acceleration from another could produce an
apparent mass shift with either estimator.

## Recoverable quantities and circularity

The appendix first infers a velocity-decay parameter from motion:

`Dhat = -(2/s) ln(u2/u1)`.

Under its gravity-plus-drag/no-lift model and constant decay-parameter
approximation, the horizontal acceleration is

`a_x = -(Dhat/2) u sqrt(u^2+w^2)`.

The velocity-decay parameter does not require the unknown inertial mass.
But the published flight drag coefficient is subsequently formed as

`C_D,flight = 4 m_used Dhat / (rho*pi*D^2)`.

Here D is shaft diameter, distinct from Dhat. Thus the plotted flight
coefficient already includes a mass used by the authors. It is not an
independently measured force.

This does NOT mean that any reuse of such coefficients is necessarily circular:
if the exact normalization mass, air properties, diameter, and associated
uncertainties were supplied per observation, the normalization could be undone
to recover the mass-independent decay measurement. Combining that recovered
motion with a separate wind-tunnel force would still be a meaningful test.
However, using the same flight-derived drag force and acceleration to recover
mass would simply return the mass embedded in the force calculation.

A ratio of tunnel and flight drag coefficients could provide a normalized
mass-consistency check under matched aerodynamic conditions. Approximate prose
summaries of coefficient levels are not matched observations. We did not turn
those summaries into an apparently precise mass estimate.

## Minimum data for a next run

- Corrected per-shot horizontal/vertical velocity measurements at both camera
  stations, geometry, and their uncertainties/covariances; preferably camera
  position/time records and calibration information.
- Independently calibrated aerodynamic force components for the corresponding
  configuration, flow speed, air properties, attitude and flow regime, with
  uncertainty and calibration provenance.
- An explicit mapping from tunnel conditions to the free-flight observations,
  including discrepancy uncertainty rather than assuming perfect transfer.
- The mass of the actual flying configuration as a withheld benchmark, with
  any normalization masses and their use recorded separately.

Use matched horizontal drag and horizontal acceleration if other horizontal
forces are negligible. For a full Cartesian comparison, use total force and
center-of-mass acceleration; drag alone is not total force in a gravitational
field. Gravity depends on the target mass, so adding an independently known
weight as though it were mass-free evidence would also need a clear measurement
model. Direct force components are preferable to a scalar coefficient when
testing vector alignment.

## Implication for the estimator project

This candidate establishes data requirements and an important model-mismatch
case, not a new empirical result. The air-track comparison remains the completed
published-data run. Large physical scatter is not automatically uncertainty
about one compatible pair, and object size is not a substitute for weak signal
relative to measurement uncertainty.

Local scope: existing `air_track_experiment.py` read for comparison; new audit
only. Existing estimator, experimental outputs, and managed context unchanged.
No new managed-context synthesis was performed; estimator provenance remains
the behavior laboratory's pinned VD-docs commit
`743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6`. No authors were contacted.
