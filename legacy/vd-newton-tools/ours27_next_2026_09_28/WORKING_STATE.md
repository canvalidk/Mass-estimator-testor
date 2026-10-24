# Internal research state — 28 September continuation

Scratch only. No new equation has been adopted or named ours28. No curated
record, handoff, production edit or repository write was requested or made.

Managed source snapshot: canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.
The main managed sources are the ours27 full covariance record, the September
28 candidate-size/null-reading record, and the September 27 lemma list under
Newton-analysis/_dscn_mass_estimator/. The theorem's source-status audit is in
regrouping/null_neutrality_width_obstruction.md.

## Strongest findings

1. Exact conditional score under full declared covariance:
   r=F-ma, V=Var(r), b=a-(Sigma_aF-m Sigma_aa)V^-1 r.
   At a specified aligned candidate experiment, b and r are independent and
   Z=(sum b_i' V_i^-1 r_i)/sqrt(sum b_i' V_i^-1 b_i) is exactly standard
   normal conditional on the fitted components, when the denominator is
   positive. See full_covariance_conditional_score.py and its results.
   This is an uncertainty test, not a normalized mass law or candidate law.

2. In the common-noise-pair experiment, exact conditional LR inversion also
   eliminates unknown candidate sizes. The numerical six-cell check gives
   94.85–95.40% coverage at nominal 95%. See exact/conditional_lr.py.
   Its confidence sets may be disconnected/unbounded. One cannot simply
   normalize the confidence curve and claim a calibrated probability law.

3. Exact appended-null mass-law neutrality obstructs universal regular width
   calibration for ordinary mass-independent, isotropic Gaussian candidate
   integration. This includes improper size measures, and does NOT assume
   candidate-prior projectivity. Null neutrality derives projectivity for the
   proper zero-data tilted candidate laws. Their Gaussian-mixture
   representation forces a common asymptotic evidence slope b; calibration
   at R^2/M -> lambda requires b=lambda/(1+lambda). See the full proof in
   regrouping/null_neutrality_width_obstruction.md.
   Universal repeated-experiment calibration is explicitly NOT settled by
   the managed sources. Full stacked O(M) isotropy is also stronger than
   ordinary physical rotation. Do not promote this into a universal no-go.

## Constructed and rejected routes

- Tensor candidate covariance repairs regular anisotropic width behavior,
  but a 2D reading and two scalar readings have different profile evidence
  unless the original dependence metadata are retained. See
  regrouping/check_regrouping.py.
- The joint Grassmann graph volume is proper and invariant under rotations,
  units and channel swap, but fails known-empty continuity. On F=2,a=0 with
  unit noise the point is 1.989891273. Adding a zero axis whose two noise SDs
  approach zero gives 1.838380705, not the original point. It also fails prior
  marginal compatibility and arbitrary flattening. See measure/
  graph_volume_checks.py and graph_audit/graph_empty_axis.py.
- Proper projective inverse-Wishart families and a sigma-finite boundary
  measure repair small-N normalization, but leave shape/scale choices. See
  measure/measure_limits.py and projective_counterfamily.py.
- Uniform radius gives an exact hypergeometric candidate equation, but is
  an additional measure choice and fails M=2 to M=1 prior projection. See
  exact/radial_reference.py.
- Add fresh Gaussian noise to observations then project: the local angular
  variance is (R^2+3M)/R^4 instead of (R^2+M)/R^4. Positive finite projection
  also fails on a positive-probability set because the optimum is at 0 or
  infinity. Conditional residual REPLACEMENT repairs the leading width, but
  fitted-center simulation is only a bootstrap. See exact/noise_projection.py.
- Inverting the exact score pivot does not automatically assign a law: its
  observed range is bounded by sqrt(lambda_max)-sqrt(lambda_min), so a fresh
  standard normal draw has positive probability of no inverse candidate.
  Multiple roots and absent roots need an allocation rule.

## What remains unresolved

The ratio-of-mean-sizes readout remains conditional on C3 and on a joint
candidate law. The exact conditional noise experiment supplies a benchmark
without assigning that law. No uniqueness argument for the full size/mass
measure has been established. Instrument-dependent global candidate metrics
can escape some continuity obstructions, but their functional form is not
selected by current symmetries. Do not confuse an escape family with a
derived next equation.

The earlier continuation is in ../ours27_repair_2026_09_28/: full covariance
J/G width diagnosis, the overly broad R^2->infinity consistency counterexample,
and the first tensor construction. Its positive tensor result has now been
restricted by the regrouping failure above.
