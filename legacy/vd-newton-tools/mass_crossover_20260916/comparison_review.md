# Scalar crossover: smoothness, accuracy, and competing readouts

**Independent investigation contribution — 16 September 2026.**
This comparison examines a proposed practical use of the mass estimator: a
continuous transition from instrument-scale values to the resolved force / 
acceleration ratio. It supports the parent investigation and does not alter
the production estimator.

## What this establishes

The ratio-of-means equation gives a well-behaved positive readout and a useful
transition to the familiar ratio. In the specified scalar model it avoids
negative estimates and division by an observed near-zero acceleration. These
are real benefits over the raw signed ratio.

**Those properties do not establish superior estimation accuracy.** The
comparisons below identify regions where it helps and regions where simple
alternatives or another readout of the same probability model do better.
In particular, the weak-acceleration / strong-force regime can make its
regularization toward the instrument ratio suppress genuinely large masses.

The relevant positive finding is narrower: this equation is one defensible,
smooth way to regularize the crossover, with good performance in portions of
that crossover. Its composition justification and the choice of an estimation
loss remain separate questions.

## 1. The exact scope

The physical direction is known and fixed. Signed scalar observations satisfy

$$\widetilde F\sim N(f,s_F^2),\quad
\widetilde a\sim N(a,s_a^2),\quad f,a>0,$$

independently, with known positive standard uncertainties. Define

$$X=\widetilde F/s_F,\quad Y=\widetilde a/s_a,\quad
\theta_F=f/s_F,\quad\theta_a=a/s_a,$$

so true mass is

$$m=(s_F/s_a)(\theta_F/\theta_a).$$

No correlation, unknown direction, censoring, estimated standard uncertainty,
or repeated-trial pooling is modeled here. Signed readings can be negative.
Observed vector norms cannot be substituted into these scalar calculations.

With the existing flat-positive magnitude law, put

$$h(t)=t+\phi(t)/\Phi(t).$$

The ratio-of-means estimate is exactly

$$\widehat m_h=(s_F/s_a)h(X)/h(Y).$$

## 2. Comparators chosen before examining results

All methods use identical observations and uncertainty scales. None receives
the true mass or additional measurements.

| Method | Standardized positive-channel transform / mass readout | Purpose |
|---|---|---|
| Ratio of means | $h(t)$; ratio of channel values | Existing readout |
| Continuous hard floor | $g(t)=\max(t,c)$, $c=\sqrt{2/\pi}$; ratio | Simple positive regularization, same value at zero and same large-positive asymptote |
| Smooth algebraic | $g(t)=(t+\sqrt{t^2+4c^2})/2$; ratio | Smooth positive comparator with the same zero and large-positive limits |
| Posterior geometric | $\exp\{E(\log f\mid X)-E(\log a\mid Y)\}$ | Different readout of exactly the same flat-positive posterior law |
| Raw signed ratio | $(s_F/s_a)X/Y$ | Familiar unregularized operation, retaining all its outputs |

The hard floor is continuous, but has a derivative discontinuity at $t=c$.
The smooth algebraic comparator has no fitted tuning parameter. It was chosen
for endpoint matching, not for optimising results on this grid. These simple
comparators do not claim to satisfy the composition theorem's readout demands.

The posterior geometric decision minimises the posterior expected squared
logarithmic error: differentiating $E[(\log d-\log M)^2\mid X,Y]$ gives
$\log d=E[\log M\mid X,Y]$. This is not a claim that it uniformly minimises
sampling error at every fixed true mass. It is deliberately labelled as a
different readout rather than a version of the same mass equation.

## 3. Reproducible design and metrics

The complete grid is

$$\theta_F,\theta_a\in\{.25,.5,1,2,3,5,8\}.$$

There are 49 cells, each with **500,000 IID Gaussian measurement pairs**, seed
20260916. The same error pairs are reused across methods and signal settings;
within each cell the trials are IID. This enables paired comparisons with much
smaller Monte Carlo uncertainty than independent simulation of each method.
No preferred distribution over grid cells is assumed and no pooled score is
used to select a winner.

For every method we report:

- Probability that its estimate is within a factor of 2 of the true mass.
- Median and 5th/95th percentiles of estimate / true mass.
- Positive-finite output rate.
- Mean squared log error for the four structurally positive methods.
- A capped squared-log error, $\min\{|\log(\widehat m/m)|,\log10\}^2$,
  assigning the maximum loss to nonpositive or nonfinite outputs.

The cap is a declared comparison convention, not a statistical theorem. The
uncapped squared-log risks of the four positive methods allow comparison
without that convention. The raw ratio's positive-subset log summaries are
separately labelled as conditional and are never substituted for all-trial
performance.

**Raw-ratio invalid outputs are not discarded.** They are factor-of-2 failures
and receive maximum capped log loss. Both-negative observations can give a
positive raw mass; their frequency is separately retained in the JSON. Its
exact nonpositive probability is

$$\Phi(-\theta_F)\Phi(\theta_a)+\Phi(\theta_F)\Phi(-\theta_a)>0.$$

Thus even a high-signal simulation with no observed invalid outputs does not
establish universal positivity. No ordinary mean squared error is computed
for the raw ratio: a Gaussian denominator has nonzero density at zero, and its
squared reciprocal has a divergent integral there.

There are two different moment questions. The **sampling distribution of the
scalar point** $h(X)/h(Y)$ has finite moments of every positive order under
fixed Gaussian measurement means: $h(t)\sim t$ for large positive $t$ and
$h(t)\sim1/|t|$ for large negative $t$, so the readout and its reciprocal grow
at most polynomially in Gaussian readings. However, the **conditional latent
mass** $M=f/a$ under the flat-positive joint law has a divergent mean because
the acceleration density remains positive at zero. The finite sampling risk
of a point readout and the infinite conditional mean of the underlying ratio
are compatible statements about different distributions. The logarithmic
moments used for the posterior geometric readout remain finite.

Each factor-of-2 probability's maximum Monte Carlo standard error is .000707
(about .071 percentage point). The JSON also reports the standard error of
each paired difference. These are simulation error estimates, not uncertainty
about whether the experimental assumptions are correct.

## 4. Factor-of-2 recovery results

Entries are percentages of **all trials**. The last column is the raw ratio's
observed nonpositive/nonfinite percentage.

| True $(\theta_F,\theta_a)$ | Ratio of means | Hard floor | Smooth algebraic | Posterior geometric | Raw ratio | Raw invalid |
|---|---:|---:|---:|---:|---:|---:|
| (.25,.25) | 72.31 | 86.18 | 61.16 | 62.29 | 21.71 | 48.14 |
| (1,1) | 71.23 | 70.26 | 66.08 | 59.90 | 36.39 | 26.77 |
| (2,2) | 76.18 | 67.45 | 77.70 | 65.55 | 62.79 | 4.46 |
| (3,3) | 85.23 | 82.10 | 88.32 | 78.36 | 81.96 | .27 |
| (8,.25) | 8.99 | 0.00 | 15.85 | 35.10 | 14.94 | 40.12 |
| (8,.5) | 48.51 | 65.91 | 44.54 | 68.08 | 28.95 | 30.88 |
| (8,1) | 79.34 | 83.39 | 70.86 | 74.74 | 52.51 | 15.89 |
| (8,2) | 89.45 | 80.82 | 89.72 | 79.98 | 80.26 | 2.25 |
| (8,3) | 94.18 | 92.17 | 95.94 | 89.25 | 92.17 | .13 |
| (8,8) | 99.965 | 99.965 | 99.975 | 99.957 | 99.965 | 0 observed |

The hard floor's zero entry at (8,.25) is a finite-simulation result, not a
proof of exactly zero recovery probability. Its denominator floor makes a
large enough output extremely unlikely in that cell.

### Findings with operational meaning

1. **Replacing the raw ratio helps substantially in much of the crossover.**
   At (1,1), ratio-of-means factor-of-2 recovery is 71.23%, compared with 36.39%
   for the raw ratio. At (8,2), recovery rises from 80.26% to 89.45%.
2. **The ratio of means is not the best positive regularizer everywhere.**
   At (3,3), the smooth algebraic comparator improves recovery from 85.23% to
   88.32%; at (8,1), the floor improves it from 79.34% to 83.39%.
3. **Strong force with very weak acceleration is an important limitation.**
   At (8,.25), true mass is 32 instrument-ratio units. The median ratio-of-means
   estimate is only .27645 of the true mass. The posterior geometric readout's
   median is .40170, and its factor-of-2 recovery is 35.10% versus 8.99%.
   Neither resolves the mass well, but the choice of readout matters greatly.
4. **Different criteria can favour different methods on the same trials.**
   At (8,2), algebraic smoothing has slightly higher factor-of-2 recovery,
   while the ratio of means has lower squared-log error. At (8,.25), the raw
   ratio beats the ratio of means on factor-of-2 success but is much worse on
   capped log loss: 3.23990 versus 1.90256, with 40.12% invalid outputs.
5. **The weak equal-signal diagonal is especially favourable to shrinkage.**
   When $\theta_F=\theta_a$, true mass equals the instrument ratio. The floor
   often returns exactly that ratio and consequently wins at (.25,.25).
   That result cannot be advertised as successful recovery of arbitrary
   masses from two weak signals. Unequal low-signal cells are essential.

## 5. Squared-log risk confirms the trade-off

These are unconditional mean squared log errors for structurally positive
methods. Lower is better. They do not use the cap or omit negative readings.

| True $(\theta_F,\theta_a)$ | Ratio of means | Hard floor | Smooth algebraic | Posterior geometric |
|---|---:|---:|---:|---:|
| (.25,.25) | .39534 | .17727 | .60739 | .59857 |
| (3,3) | .23254 | .28703 | .20050 | .33289 |
| (8,.25) | 1.91492 | 1.85132 | 1.99823 | 1.21414 |
| (8,.5) | .74103 | .61776 | .85518 | .48200 |
| (8,1) | .28485 | .23513 | .36999 | .33683 |
| (8,2) | .18704 | .24932 | .19033 | .31071 |

The posterior-log action's conditional decision optimality therefore has no
uniform fixed-parameter sampling dominance: it helps at (8,.25) and (8,.5)
but loses clearly at (8,2).

## 6. What smoothness does and does not prove

The exact scalar transform is positive and smooth,

$$h(0)=\sqrt{2/\pi},\qquad h(t)-t=\phi(t)/\Phi(t)\to0
\quad(t\to+\infty).$$

At exact measured zero-zero the mass point is $s_F/s_a$. At well-resolved
positive readings the formula approaches their ratio. Some channel values:

| Reading / sigma | $h(t)$ | Hard floor | Smooth algebraic |
|---|---:|---:|---:|
| 0 | .797885 | .797885 | .797885 |
| .5 | 1.009160 | .797885 | 1.086134 |
| 1 | 1.287600 | 1 | 1.441605 |
| 2 | 2.055248 | 2 | 2.279304 |
| 3 | 3.004438 | 3 | 3.199006 |
| 5 | 5.000001 | 5 | 5.124237 |

The ratio of means reaches the resolved positive channel reading much faster
than this algebraic comparator. The comparator can nevertheless have smaller
sampling risk: its residual regularization can reduce random spread enough
to compensate for bias. Therefore closeness to raw division on noiseless or
central readings is not by itself an accuracy criterion.

For positive central readings, the weak channel receives the larger
proportional correction. Consequently the mass ratio is pulled toward the
instrument ratio. This is useful when that pull reduces noise, and costly
when true mass is far from that ratio. The effect is not an error in the
equation; it is a trade-off that a practical claim must expose.

If both true signals approach zero while their ratio remains fixed, the
sampling distributions of the observations approach the same zero-mean
Gaussian distribution for many different masses. No smooth crossover can
create the missing information. The zero-zero *point* also must not be
confused with the estimator's sampling distribution under zero signals:
noisy readings still vary.

## 7. Numerical verification and reproducibility

Run from the workspace root:

```powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '.tools/mass_crossover_20260916/scalar_comparison.py' --trials 500000
```

The JSON contains all 49 cells, paired standard errors, conditional versus
unconditional metric labels, seed, method definitions, quadrature checks,
and the script SHA-256.

- Direct $h(t)$ evaluation agrees with the existing scalar implementation to
  maximum absolute error $1.28\times10^{-13}$ on the validation range.
- The posterior geometric comparator uses quadrature in $u=\log a$ to remove
  the logarithmic endpoint singularity, then a precomputed interpolation
  table independent of the trial results. Against order-768 quadrature,
  maximum checked log-expectation error is $8.61\times10^{-8}$.
- Its zero-reading log expectation also matches the exact half-normal log
  moment $-(\gamma+\log2)/2$ to the asserted tolerance $2\times10^{-12}$.

This is a simulation study within a specified model, with direct analytical
checks. It does not establish a best estimator for every operating condition,
select a reference measure, or transfer the scalar conclusions to an unknown
vector direction.

## Suggested contribution to the larger argument

The practical crossover claim can be stated precisely: the equation supplies
a smooth positive mass readout that automatically regularizes uncertain
channels and approaches the resolved force / acceleration ratio. A domain
claim should additionally specify signal-to-noise, acceptable multiplicative
error, and whether preserving the composition readout is a requirement.

**Smooth crossover is a useful behaviour; the present comparison does not
make it a unique advantage or a general accuracy optimum.** These findings
give a concrete starting point for identifying use cases where the existing
readout's other guarantees justify its particular bias–spread trade-off.
