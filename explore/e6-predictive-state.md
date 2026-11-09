# E6: which state preserves the allowed future experiments?

**Codex, 9 November 2026; continuation after `05b9309` on
`explore/e6-declared-structure`.** Pressure/derivation, not a new premise
or adopted estimator. **Check:** `tests/test_predictive_state.py`.

The operational fork now has a conditional answer. Preserving the mass
law and all future experiments that depend only on mass requires only
the mass marginal. If an addressed fixture can be calibrated again,
preserving its event probabilities and subsequent mass laws can require
the entire joint mass–fixture law. Even both separate marginals can fail.
This requirement is fixed by the admitted operations, not by Born notation.

A second result limits the bridge: at fixed finite data, a single
parameter-independent reduction whose entire normalized output law is
independent of arbitrary push amplitudes cannot retain mass information
when every hidden mass–fixture state shares one zero-signal law. This
rules out that proposed route to a nuisance-free experiment, not mass
inference or the project.

## 1. Standard mathematics and the exact operational criterion

This is **predictive sufficiency relative to declared future protocols**.
[Shalizi and Crutchfield, *Computational Mechanics*](https://arxiv.org/html/cond-mat/9907176v2),
§IV.1 Definition 5 and §V Theorem 2/Remark 2, groups histories by equality
of the distributions of futures of every length; one-step equality alone
does not suffice. [Jaeger, *Observable operator models*](https://www.ai.rug.nl/minds/uploads/oom_neco00.pdf),
§5 and Proposition 4, uses predictive spaces and linear observation
operators to obtain minimal-dimensional realizations. The finite static
specialization below is elementary. No stationary signal process or
population of pushes is imported from those time-series models.

Fix a finite persistent hidden state `h`, initially distributed as `Q`.
For each admitted setting `s`, specify normalized outcome probabilities

$$e_{a|s}(h)\ge0,\qquad \sum_a e_{a|s}(h)=1.$$

Fresh observation errors are independent conditional on `h`. Settings,
fixture identities, and which observations are genuinely fresh must be
part of this declaration. For an outcome string `w`, its effect and
probability are

$$e_w(h)=\prod_r e_{a_r|s_r}(h),\qquad P_Q(w)=Qe_w.$$

Also retain requested readouts `g(h)`, including `1`. For a mass law on
a finite mass grid these include every mass indicator `g_i(h)=1_{m=m_i}`.
Define

$$V=\operatorname{span}\{e_w g:\ w\text{ any admitted finite string},
\ g\text{ a requested readout}\}.$$

**Criterion.** Two beliefs are operationally equivalent precisely when
`Qv=Q'v` for every `v in V`, if the requirement is to preserve outcome
probabilities and all requested conditional readouts. An adequate retained
state must not merge inequivalent beliefs; it may retain more information.
For `P_Q(w)>0`,

$$E_Q[g\mid w]=\frac{Q(e_wg)}{Qe_w}.$$

Equality on `V` preserves numerator and denominator. Conversely, equality
of all those event probabilities and conditional readouts preserves their
products, hence every generating expectation. Zero-probability histories
need no conditional prediction. Adaptive settings introduce no new issue
when each path is an admitted string.

The space is closed under multiplication by every admitted effect whenever
the corresponding repetitions/continuations are allowed. This is what
makes updates possible from its coordinates:

$$Q^a(v)=\frac{Q(e_a v)}{Qe_a}.$$

For unrestricted repetitions, starting with `1` and multiplying effects
generates a finite-dimensional commutative algebra. It consists of all
functions constant on hidden states having identical effect-value vectors.
Proof: a generic linear combination separates the finitely many distinct
vectors, and Lagrange interpolation gives each class indicator. If mass
readouts are also required, refine the classes by mass identity. The exact
state is the probability of each such class, up to an equivalent encoding.
This is not a lower bound on bytes or on pathological nonlinear encodings.

**Why repetition matters.** A persistent Bernoulli bias `p in {.2,.5,.8}`
has beliefs `Q=delta_.5` and `Q'=.5 delta_.2+.5 delta_.8`. Both predict
next-yes probability `.5`. Their probabilities of two yes outcomes are
`.25` and `.34`; after one yes their next-yes probabilities are `.5` and
`.68`. The one-step space omitted `p^2`; closure recovers all three hidden
weights. Redrawing the bias at each step changes the experiment.

Without allowed combinations, even separating pure hidden states is not
enough to identify mixtures. Separately observing either bit of a hidden
two-bit table leaves its association unknown: equal diagonal versus equal
off-diagonal weights give the same two marginal observation laws. Observing
both bits on the same persistent state distinguishes them. Thus the
one-step expectation space must not silently be replaced by its algebra.

## 2. A small normalized fixture experiment forces the joint association

Here is a concrete family whose closure can be calculated completely.
Let `h=(m_i,[n_j])`, with `M` distinct masses and two unoriented axes
`[n_1]=[e_1]`, `[n_2]=[e_2]`. Admit a future fixture-calibration readout

$$Z\mid[n]\sim N(n_x^2,1),$$

whose noise is fresh and independent of prior records conditional on the
same persistent fixture. The sign-invariant squared component avoids an
orientation convention. This is an explicitly specified research
experiment, not a claim that the current input type already permits it.
There is no signal amplitude or population in this channel.

The exhaustive outcomes `B={Z>1/2}` and its complement have effects

$$e=(a,b),\quad a=\Phi(1/2),\quad b=\Phi(-1/2)=1-a,\quad a\ne b.$$

For each mass, preserving the current mass law and the mass numerators
after this outcome retains

$$s_i=Q(g_i)=q_{i1}+q_{i2},\qquad
t_i=Q(eg_i)=a q_{i1}+b q_{i2}.$$

These recover every joint cell:

$$q_{i1}=\frac{t_i-bs_i}{a-b},\qquad
q_{i2}=\frac{as_i-t_i}{a-b}.$$

Thus `{g_i,eg_i}` has rank `2M`; a normalized state has `2M-1` independent
probabilities. Closure stops immediately because `e^2=e-ab`. An exact
state preserving this experiment's probabilities **and** conditional mass
laws must distinguish any two different joint laws on this grid. This is a
claim on the class of all joint beliefs; a restricted reachable subclass
could have a shorter parameterization. One useful equivalent encoding is
the mass law together with the conditional fixture law given mass, on
the mass law's support. The fixture marginal alone does not supply that
conditional law.

The failure of separate marginals is explicit:

$$Q=\begin{pmatrix}.4&.1\\.1&.4\end{pmatrix},\qquad
Q'=\begin{pmatrix}.1&.4\\.4&.1\end{pmatrix}.$$

Both marginals agree and both assign probability `.5` to `B`. After `B`
their mass laws are `(.8a+.2b,.2a+.8b)` and its reverse, approximately
`(.61488,.38512)` and `(.38512,.61488)`. A calibration likelihood with no
direct mass dependence can change mass through an existing association.
It is not a mass-neutral update in this example.

**Two controls prevent overclaiming.** If all admitted effects and requested
readouts depend only on mass, their entire closure depends only on mass:
the mass marginal is then an exact state. Conversely, preserving only
conditional mass laws under calibration alone is weaker than preserving
the operational experiment. Every product state `p(m)q(n)` with a fixed
`p` keeps that mass law under all calibration outcomes, even when the
outcome probabilities differ. The full-joint conclusion above explicitly
requires both. Truly duplicate axis labels may also be merged within
each mass if all their allowed experiments agree.

## 3. A normalized raw-reading witness, with no amplitude population

Use the existing E6 mass-neutral history: masses `.5,1,2` have weights
`(1/4,1/2,1/4)`, axes are `e_1,e_2`, and initial direction weight is
`w=1/4`. The two observed pairs `(F,a)=(e_1,e_1),(e_1,-e_1)` multiply
the axes' profile weights by `exp(-1),exp(-2)`, independently of mass.
They leave the mass law unchanged but give `w'=e/(e+3)`. The retained
joint here is parent §13's **specified working joint law**; the normalized
future experiment below does not retroactively derive a Bayesian
posterior from an undeclared population of past amplitudes.

Now specify the normalized conditional future experiment

$$(F,A)\mid(m,n,t)\sim N_4((mt n,t n),I_4).$$

Here `t` is an arbitrary fixed next acceleration amplitude. It is not
drawn from a population, inferred from earlier push sizes, or integrated
out. The fixed raw-data event `B={|A_x|>1}` has probabilities

$$p_1(t)=\Phi(t-1)+\Phi(-t-1),\qquad p_2=2\Phi(-1)$$

at the two axes. These define a normalized binary experiment. `p_1` is
even, `p_1(0)=p_2`, and for `t>0`,

$$p_1'(t)=\varphi(t-1)-\varphi(t+1)>0.$$

Thus `P_t(B)=w p_1(t)+(1-w)p_2` changes after the history for **every
nonzero** `t`, with no choice of event depending on `t` or mass. At the
slice `t=1` it changes from `.368670413884` to `.414969702311`.
Mass neutrality therefore does not mean predictive irrelevance for this
family. These are genuine raw-event probabilities, unlike the future
profile-factor expectations used in the previous continuation.

**What is preserved:** a family `t -> P_t`, not a single unconditional
forecast. Even the full `Q(m,n)` does not specify the next amplitude.
A known setting or a supplied conditional amplitude law could specify
one forecast; the present type supplies neither in general. Requiring
one would be an extra obligation, not something already implied by
returning a mass law.

An exact-zero *observed* reading remains neutral in the fixed-known-noise
profile construction: it leaves the entire `Q`, and therefore this whole
conditional forecast family, unchanged. That statement is different from
the next signal's unknown amplitude being zero.

### A richer raw-event family distinguishes every joint weight

For any fixed `t!=0`, write `b_h=(m_i n_j,n_j)` and use the common known
Gaussian covariance. If two finite mixtures have equal probabilities for
all raw events, their Fourier transforms give

$$\sum_h(Q_h-Q'_h)\exp(i s^Tt b_h)=0\quad\text{for every }s.$$

The Gaussian characteristic factor has no zeros and cancels. After
aggregating coincident means, choose a projection separating the finitely
many means. Differentiating at zero gives a Vandermonde system, so every
aggregate weight agrees. This is the finite known-location Gaussian
mixture argument, not a claim of stable estimation from one observation.

For unoriented axes, explicitly restrict to inversion-symmetric events
`B=-B`. Their laws can be represented mathematically by
`S_h=.5 N(t b_h,Sigma)+.5 N(-t b_h,Sigma)`, whose characteristic factor
is `cos(s^Tt b_h)`. Expanding the cosines into exponentials gives the same
independence for distinct pairs `{+b_h,-b_h}`. Positive mass and distinct
unoriented axes give distinct pairs. The symmetric-event family therefore
distinguishes every finite joint law at each nonzero amplitude slice.
This symmetrization is a representation of a restricted event domain;
it does not assert a physical sign population. Unknown signed amplitude
alone is not a declaration to discard all asymmetric events.

Restricted sensors can preserve less. The tests exhibit normalized
symmetric projection events with full rank on the six-cell E6 grid,
and a rank-one collapse at zero amplitude. Full-rank exact separation
does not give a useful numerical conditioning bound near that limit.

## 4. A common null rules out the strongest nuisance-free reduction

State the proposed bridge requirement precisely. Let `P_(h,t)` be a
fixed finite raw-data experiment, with a common law `P_0` at zero signal
for every hidden mass/fixture state. Suppose a **single** Markov kernel
`K(y,dz)`, independent of `h,t`, produces a normalized reduced experiment
whose **entire law** is independent of amplitude:

$$P_{h,t}K=R_h\quad\text{for every allowed }t.$$

If `t=0` is permitted, immediately `R_h=P_0K` for every `h`: the reduced
experiment is uninformative about `h`. If zero is excluded but can be
approached, the same conclusion follows from total-variation continuity:

$$\|P_{h,t}K-P_0K\|_{TV}\le\|P_{h,t}-P_0\|_{TV}\longrightarrow0.$$

For equal-covariance Gaussian shifts the right side is exactly
`2 Phi(||mu_(h,t)||_(Sigma^-1)/2)-1`, which tends to zero. Independent
randomization cannot evade the argument because it is included in `K`.
The reasoning also covers a finite series with independently free
amplitudes, by sending every amplitude to zero.

An approximate version makes the quantifier visible. If a common kernel
is within `epsilon` in TV of amplitude-independent laws `R_h`, uniformly
down to zero signal, then `TV(R_h,R_h')<=2 epsilon`. Any substantial
separation must therefore give up such uniform amplitude independence.
The common-null idea is already used in T45/B13/O15's coverage arguments;
this is its application to a proposed normalized experiment reduction.

**Scope.** This does not rule out computing an estimator, a profile or
integrated likelihood, mass-dependent residual pivots, exact confidence
sets, conditional inference, or approximations in an identifying regime.
For example `(F-mA)/sqrt(1+m^2)` is nuisance-free at its matching true
candidate in the scalar unit-noise model. Changing `m` changes the
data-to-output map, so these pivots are not one kernel `K` on the data.
T45 already constructs exact candidate-wise tests without a population.
A positive signal floor, changing noise laws, or an infinite-data limit
requires its own analysis; none is resolved by this short proof.

The mass definition is not required here to forecast an unspecified next
push. Imposing an informative, normalized, exactly amplitude-free raw-data
experiment would add an incompatible bridge requirement on this slice.
The conclusion is to reject that proposed forcing step, not L10 or the
existence goal.

## 5. Positive control: GLS is an actual experiment reduction

The parent exploration's fixed-linear result has a stronger operational
description. For `Y~N(L theta,Sigma)`, with fixed full-rank `L` and known
positive-definite `Sigma`, put

$$V=(L^T\Sigma^{-1}L)^{-1},\qquad
S=VL^T\Sigma^{-1}Y.$$

Then `S~N(theta,V)`. Independently generate

$$\xi\sim N(0,\Sigma-LVL^T),\qquad Y^*=LS+\xi.$$

The residual covariance is positive semidefinite, possibly singular.
`Y*` has exactly the original mean and covariance, hence the original
normalized Gaussian law for every `theta`. This reconstruction kernel
depends on the fixed instrument and mean space, not the unknown parameter.
Thus every raw-data event probability is recoverable from the reduced
experiment. Its retained law still depends on the unknown signal
coordinates: it does not contradict §4 by eliminating all amplitudes.

No new Newton-II loophole follows: the mean-space declaration's physical
provenance remains the issue in parent §9. Unknown noise scales can make
the discarded residual informative (parent §12); its realized scale
cannot then be silently omitted from this reconstruction kernel.

## 6. What this closes, and the remaining physical fork

**Forced, conditionally:** once a normalized experiment family, its shared
identities and its requested readouts are specified, the closed expectation
space above determines the exact distinctions a retained state must keep.
For the two-axis calibration family it keeps the full joint association;
for mass-only operations it keeps only mass. Arbitrary extra fixture memory
is not required by predictive sufficiency when no admitted operation can use it.

**Refuted proposed bridge:** an informative normalized raw-data reduction
which removes unrestricted amplitudes exactly through one common kernel
cannot coexist with this common null. A candidate-residual density does
not establish such a reduction. Born notation supplies neither a missing
amplitude law nor permission to substitute a different sampling experiment.

**Still open:** whether addressed calibration uncertainty and later
recalibration belong to the physical mass-definition type, and how they
interact with literal L19. The current type does not supply this channel;
the explicitly declared fixture example is conditional input to that sort,
not permission inferred from its mathematical convenience. The Gaussian
profile joint also still fails T30. No production law or premise changes.
This does not settle the nuisance measure in T48/G19 or the consistency
rate and coverage questions in the probe line.

The next concrete physical task is to distinguish a fresh mass-only route,
a fresh fixture observation, and reused evidence in the typed combination
operation. The result here supplies a criterion for that task: an allowed
continuation must act identically on retained states declared equivalent.
It does not require a quantum effect algebra or a signal population.

## Validation

All **24 new checks pass**. Together with the Born bridge, T70, E6,
graph geometry, normalized residual mixture, conditional-score and audit
regressions: **188 passed**. Independent mathematical/source reviews
found no blocker after clarifying operational equivalence, the common
null law and the working-state status of the past-data joint.

The tests check the normalized witnesses, closure ranks and reconstruction
identities; the general claims are proved above, not established by finite
test counts. `tools/audit.py` retains its existing T6/T21/T27 and unused-L24
findings. `git diff --check` passes. The full default suite was not repeated
for this documentation/test-only continuation; the earlier full-suite
result in the parent exploration belongs to its first commit.
