# Independent audit: transition study

2026-09-16. Scope: `crossover_vector_checks.py`, comparison semantics, and
independently derived scalar/3D-axis results. No production edits.

## Confirmed formulas and implementation

For independent isotropic unit channel errors in 3D, integration over a
uniform common direction contributes `sinh(r h)/(r h)`, where
`h²=p sin²(theta)+q cos²(theta)+2 dot sin(theta) cos(theta)`.
Flat magnitude measure gives radial power one; tube weighting adds two powers.
Writing `J_k=integral_0^infinity r^k exp(-r²/2)sinh(rh)/(rh) dr`, the identities are

- `J1=exp(h²/2)*sqrt(pi/2)*erf(h/sqrt(2))/h`, with continuous value 1 at h=0;
- `J2=sqrt(pi/2)*exp(h²/2)`;
- `J3=(h²+1)*J1+1`;
- `J4=(h²+3)*J2`.

These yield the script's densities and first-moment integrands. Direct radial
Gaussian quadrature for seven h values from 0 to 20 and all four powers agreed
to 9.66e-15 relative error. Stable common scaling uses the correct maximum
over the positive theta quadrant, including negative observed dot products.

The finite positive profile stationary point solves
`dot*m²+(q-p)*m-dot=0`. The script uses its stable positive root when dot>0.
For dot<=0 the optimum is generally an unattained zero or infinite boundary.
Three hundred random-vector checks verified interior stationarity to 1.64e-14
and objective values no greater than either boundary infimum.

The factor-two success comparison includes all trials and uses paired
differences on the same observed data. The posterior central-95% coverage event
is correctly implemented as `.025 <= CDF(true_mass) <= .975`.

## Issues raised and fixed by lead

1. The boundary-valued positive-profile rule has infinite population squared-log
   risk under every finite Gaussian signal: opposite/nonpositive observed dot
   products have positive probability. A finite sample missing those events
   cannot establish finite risk. Lead now reports its population log risk as
   infinite and its log-MSE field as null in every setting.
2. If p=q and dot=0, every positive mass is a profile optimum. In particular
   zero/zero observations cannot truthfully be labeled an unattained infinite
   boundary. Lead now returns NaN to mark this nonidentification tie.

## Additional analytical contributions reviewed

### Scalar log sensitivity

For independent signed scalar channels, the point is proportional to H(x)/H(y),
with H the positive-truncated Gaussian mean in noise units. H'=Var is positive
and less than 1. The truncated normal has increasing hazard, so its normalized
survival S satisfies S(x+y)<=S(x)S(y). Double integration gives
`E[X²]/2 <= E[X]²`, hence `Var <= H²`. Thus

`0 < d log H/dt <= min(H,1/H) <= 1`,

and the claimed global coordinatewise log-Lipschitz bound follows. This is a
data-sensitivity guarantee under that declared scalar model; it is not a
uniform bound on estimation error relative to unknown truth.

### Exact unknown-direction 3D axis

When observed force is zero and observed acceleration has standardized norm t,
the force magnitude remains half-normal with mean `sqrt(2/pi)*sigma_F`.
Define `D(t)=exp(-t²/2)*integral_0^t exp(s²/2) ds`. The acceleration mean is

`sigma_a * erf(t/sqrt(2))/D(t)`.

The proof follows by integrating `sinh(alpha*t)/(alpha*t)`: its normalization
is `sqrt(pi/2)/t * integral_0^t exp(s²/2) ds`, and its first moment is
`sqrt(pi/2)*exp(t²/2)*erf(t/sqrt(2))/t`. At zero, take the continuous limit
sqrt(2/pi). The resulting point agrees with the production angular integrator
at eight t values from .001 to 20 to 1.24e-12 relative error.

## Comparison and scope qualifications

- `norm_floor` is one explicitly fixed heuristic, not an optimized family of
  norm regularizers. Superiority to this rule cannot establish superiority to
  all practical regularization alternatives.
- The known-direction point in observed curves uses additional information.
  It can explain directional effects but is not a same-information competitor
  for the unknown-direction sampling experiment.
- Tube weighting should remain labeled as a change of reference measure, not
  another implementation of the user's flat law.
- Central 95% interval coverage must not be conflated with the earlier 95%
  one-sided lower-bound guarantee.
- The paired Monte Carlo standard errors quantify sampling uncertainty for
  specified contrasts. A finite grid does not identify exact crossover curves,
  rule out untested performance reversals, or prove global dominance.
- Refinement of 128 initial runs per scenario is a useful checked-case
  comparison, not a uniform numerical bound over all possible Gaussian draws.
- The code's finite quadrature can fail at arbitrary extreme signals or invalid
  ratio inputs. It is currently a scoped research script, not a validated
  general-purpose numerical API.

## Completed manuscript review

Reviewed `crossover_regime.md` against the scalar/vector JSON outputs and the
completed independent mathematical note. Scalar and vector performance tables,
paired errors, posterior coverage, axis values, sufficient thresholds, and
reported numerical-check magnitudes agree with the output.

Three corrections were requested:

1. Observed-angle points at 0/90/180 degrees are `2.9927539169`,
   `4.4560376562`, `6.7650170275`; rounded to three decimals these are
   `2.993`, `4.456`, `6.765`. The manuscript's earlier values were inaccurate.
2. Derivative-by-derivative domination alone proves smoothness, not real
   analyticity. The Gaussian tail also dominates the entire exponential
   kernel on every compact complex neighborhood of the readings. Stating
   this supplies the holomorphic-integral argument and therefore real
   analyticity. Positivity at a real reading makes the denominator nonzero
   in some neighborhood, so the ratio is locally analytic as claimed.
3. The opening resolved-ratio limit needs the qualification that the two
   nonzero observed vectors resolve compatibly/codirectionally. Fixed
   incompatible readings can have very different small-noise limits; the
   independent mathematical note correctly explains this distinction.

Minor rounding note: scalar raw-invalid percentage at true signals (8,2)
is exactly 2.255% in the sample, conventionally displayed as 2.26% rather
than 2.25% at two decimal places.

The candidate application in partial signal is responsibly bounded. The
manuscript distinguishes criteria, individual settings, reference-measure
changes, additional information, conditional intervals and sampling coverage.
It does not establish a general dominance or deployment claim. The positive
and negative results can be presented together without further qualification
beyond the corrections above.
