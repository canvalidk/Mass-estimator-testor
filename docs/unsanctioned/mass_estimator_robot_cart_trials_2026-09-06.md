# Mass estimator: robot and cart data trials

Status: local analysis, awaiting submission to VD-docs. No managed documents changed.

## Outcome

Both public datasets were downloaded and the estimator was evaluated on their actual readings. Strong scalar readings give essentially the ordinary force/acceleration ratio. Weak readings can give finite positive readouts where that ratio is undefined or negative, accompanied by very broad compatible-pair mass laws. This is a behaviour demonstration, not evidence of improved mass accuracy.

The robot provides the better experimental setting. The cart data expose an inadequate measurement/model combination during short collisions. Neither source provides an independently weighed reference mass for the runs examined.

## Sources and reproducibility

Estimator specification: `canvalidk/VD-docs`, commit `743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6`, `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`. This is the relevant pinned managed source used by the existing local implementation; this trial is not an exhaustive review of VD context.

External sources:

- [Skrede's Zenodo dataset, version 2](https://doi.org/10.5281/zenodo.11096791). All 12 CSVs were downloaded, and their MD5 checksums match the record metadata.
- [Skrede, A Linear Discrete Kalman Filter to Estimate the Contact Wrench of an Unknown Robot End Effector](https://www.researchgate.net/publication/384034117_A_Linear_Discrete_Kalman_Filter_to_Estimate_the_Contact_Wrench_of_an_Unknown_Robot_End_Effector), RCAR 2024, DOI 10.1109/RCAR61438.2024.10671273. Author-uploaded full text, especially equations 23–25 and the limitations discussion.
- [University of Tennessee EF105 physics cart page](https://efcms.engr.utk.edu/ef105-2021-01/grav/datasets/physics-cart), with its [cart_data.zip](https://efcms.engr.utk.edu/ef105-2021-01/grav/datasets/physics-cart/cart_data.zip). The two extracted CSVs contain 128 flat-track and 158 incline-track rows.

Production surface: `03_trace_and_evaluator/mass_estimator/`. Raw files are retained in `raw/robot/` and `raw/cart/`. SHA-256 hashes for every downloaded/extracted input are recorded in `results/robot_cart.json`.

Run `robot_cart_experiments.py`, then `robot_vector_check.py` from that directory. Outputs: `results/robot_cart.json`, `results/robot_cart_examples.csv`, the three calibrated robot z-channel CSVs, and `results/robot_full_vector.json`. These scripts import the existing estimator rather than replace its equation.

## Robot: calibration and strong signals

The source study estimates end-effector mass and uses it to compensate force measurements for motion and gravity. It reports **0.932 kg**. Its authors explicitly leave independent accuracy validation to future work; this published number is an estimate, not ground truth.

Our ordinary calibration fit across the 24 recorded orientations is

$$
\mathbf F_i=\mathbf b_F+m\mathbf g_i+\mathbf e_i,
$$

which gives **0.932308057 kg**, reproducing the reported rounded value. The fitted force offset is `[9.076516, -1.018148, 9.984755] N`.

Accelerometer readings have magnitudes near one at rest. We use the supplied frame rotation and the calibration gravity magnitude 9.82085 m/s² to express them in force-sensor coordinates and units. The force/acceleration signs are matched consistently to the calibration records. The fitted IMU offset in its original coordinates is `[-0.00366873, 0.00886992, 0.07710856]` in the raw approximately-g units. This conversion also reproduces the paper's stationary acceleration variances to rounding.

For a separate strong-signal calculation, offsets are fitted on the first 12 orientations and applied to the remaining 12. Each held-out measurement is projected along its recorded gravity direction, and the signed projections are averaged before taking the ratio. The ordinary readout is **0.930776365 kg**. The extended scalar integral is **0.930776365 kg**, agreeing to numerical precision. Individual held-out ratios range from approximately 0.9252 to 0.9364 kg.

These are support-force/specific-acceleration measurements. A stationary supported object has a nonzero accelerometer signal. We have not called that signal its inertial acceleration, or treated its static gravity/support pair as a net-force zero/zero measurement.

## Robot: weak components and model failures

For dynamic runs, subtract the calibrated offsets, rotate/sign-convert the IMU data, subtract the documented 8416 microseconds from IMU timestamps, and linearly interpolate force samples to overlapping IMU timestamps. No extrapolation is used. The raw/bias wording in the dataset description is inconsistent; the numerical offsets and calibration records support treating these CSVs as biased raw readings.

The exploratory uncertainty budget adds stationary variance to between-orientation calibration residual variance, component by component. This deliberately retains a model/calibration floor instead of dividing all uncertainty by sample count. For z, it gives sigma_F = **0.0948313 N** and sigma_a = **0.131585 m/s²**. Cross-channel covariance is assumed zero. These are our uncertainty assumptions, not a complete published calibration budget; interpolation correlations and systematic rotation effects are not fully modeled.

Define a weak z reading by both absolute readings being no more than twice these standard uncertainties. There are **123 of 1591** overlapping baseline samples meeting that criterion. We select the first such sample, without selecting on the resulting mass:

| Readout from that same scalar measurement | Result |
|---|---:|
| Ordinary signed Fz/az | −1.78798 kg |
| Extended positive-mass readout | 2.07996 kg |
| Central 95% interval of the induced compatible-pair mass law | 0.07177–68.50059 kg |

The finite positive point does not identify a two-kilogram tool. The extremely broad law expresses weak information. Altering only the assumed force uncertainty also changes the readout; the sensitivity runs are retained in the output.

The other coordinates of this same measurement are strong: force norm about **9.266 N**, acceleration norm about **9.760 m/s²**. The separate `robot_vector_check.py` calculation evaluates the full unknown-direction three-dimensional integral, including those coordinates, and gives **0.949607834 kg**, with an induced central interval **0.915108–0.985066 kg**. Thus a near-zero component is not a zero vector. Exact numerical values and quadrature refinements are in `robot_full_vector.json`.

For this SAME full-vector reading, ordinary vector least squares gives **0.948955609 kg**. A conventional errors-in-variables profile fit, using the SAME covariance as the extended equation and minimizing Q(m), gives **0.949670869 kg**. The extended readout is **0.00664% lower** than that matched-covariance fit. This is a resolved numerical difference between estimation rules, not a demonstrated accuracy gain. Comparing with the paper's 0.932 kg additionally changes the observations and procedure, so that difference must not be attributed solely to our equation.

Dynamic measurements also expose model limitations. Using the calibration mass and the stated diagonal uncertainty budget, 338/1591 baseline z readings have Q(m)>9; the contact run has 748/1589. This is a discrepancy diagnostic, not a calibrated rejection rate, because observations are correlated and the budget is provisional. During contact the wrist sensor does not independently measure all the net force needed by a bare F=ma inversion. Rotational acceleration differences between the IMU and centre of mass also remain. We retain these failures rather than interpret every returned scalar as a valid mass.

## Cart: actual zero readings, poorly sampled impacts

The two exports contain both an acceleration column and separate x/y/z accelerometer channels. We use the x accelerometer channel with the bumper-force channel. During incline coasting the longitudinal accelerometer stays near its offset while the other acceleration column is about 0.65 m/s²; this is consistent with specific acceleration and gravitational acceleration being different quantities. The instructor page's broad account of derived channels does not establish the precise filtering or timing of each exported channel.

The first 30 stationary flat-track rows set force offset **0.700333 N** and x-acceleration offset **−0.100000 m/s²**. The latter is quantized at 0.1 m/s²; identical exported values do not demonstrate zero sensor uncertainty. Our provisional Gaussian budget combines observed baseline scatter and a uniform-rounding approximation using export steps 0.01 N and 0.1 m/s². It gives sigma_F = **0.00343188 N**, sigma_a = **0.0288675 m/s²**. This is only a minimal numerical sensitivity scenario, not a sensor accuracy specification. Reusing the offsets for the incline assumes unchanged sensor calibration.

At held-out stationary time **1.50 s**, the offset-corrected acceleration is zero and force is approximately −0.000333 N. Ordinary division is undefined. The extended scalar readout is **0.119444 kg**, with induced central 95% interval **0.004693–3.040054 kg**. The cart's actual mass has not thereby been measured.

A separately labelled, constructed exact-zero control gives **0.118884 kg = sigma_F/sigma_a**. That row is synthetic; its uncertainty scales come from the real-data scenario. It makes the scale dependence explicit.

Only negative-direction bumper impacts are used for impact diagnostics, excluding the initial hand launch and opposite-end reversal. The data spacing is 0.05 s. The eight qualifying samples yield:

| Recording | Times (s) | Ordinary instantaneous ratios (kg) |
|---|---|---|
| Flat | 3.85, 3.90 | 0.26988, 0.09407 |
| Incline | 1.95, 2.00 | 0.43887, 0.05407 |
| Incline | 4.85, 4.90 | 0.37754, 0.21740 |
| Incline | 7.40, 7.45 | 0.32034, 0.16331 |

The extended equation agrees with each of these strong scalar ratios to numerical precision. Their disagreement within a single impact is much larger than the provisional rounding/noise law allows. Narrow law intervals calculated under that budget are therefore not credible total mass uncertainties. The recordings have only one or two substantial bumper-force samples per collision, and their acceleration channels differ sharply during impacts. Sampling, channel response, and omitted mechanical effects require investigation before any mass-accuracy comparison. An instantaneous uncertainty-aware equation cannot by itself reconstruct a missed impulse or supply unmeasured forces.

## Interpretation and checks

All displayed intervals are quantiles of the estimator's specified induced compatible-pair law. They are not automatically frequentist confidence intervals, and their coverage has not been established by these trials. In weak cases they can have very long tails even when the ratio-of-means readout is finite.

The scalar examples are restrictions of the full equation to a line. The robot vector check separately integrates over unknown common direction on the sphere. It rotates coordinates and the complete covariance together to resolve the angular peak efficiently; it does not fix the direction. Angular and ratio quadratures are independently refined. Scalar quadrature is doubled until both the point and the two interval endpoints change by less than 0.02%.

All 40 existing estimator tests pass. Raw robot checksums pass, calibration reproduces the published rounded estimate, and numerical convergence is checked in the experiment scripts. These checks establish reproducibility and implementation consistency, not superiority over conventional mass estimators.

The useful next experimental target is a known, independently weighed mass with synchronized force and acceleration channels, a verified complete force model, repeated weak excitation, and a calibration/covariance budget. These data have already served a useful purpose: they distinguish a finite readout, actual information about mass, and adequacy of the measurement model.
