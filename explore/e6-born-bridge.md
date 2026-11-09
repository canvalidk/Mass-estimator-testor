# E6 and Born reconstruction: what would the operational bridge need?

**Codex, 9 November 2026.** Continuation on `explore/e6-declared-structure`,
after `c5591d3`; input for later checking/harvesting, not an adopted premise.
**Standard mathematics:** Gleason-type probability assignments on effects;
commutative probability algebras; statistical sufficiency and equivalence of
experiments. Primary reconstruction sources are identified in §1.
**Checks:** `tests/test_born_operational_bridge.py`, together with the existing
T70 and E6 tests. No estimator, production default or ruling changes.

The connection to E6 is direct, but does not yet force a mass equation.
On a finite grid, E6's joint state has an exact diagonal Born representation;
its marginalization
and conditioning have partial-trace and diagonal Lüders representations.
Every bounded positive update kernel has that representation. T70 supplies
a different, two-dimensional squared-projection identity. Neither result
establishes the operational hypotheses of a Born reconstruction. Below are
specific obstructions to identifying them, and a narrower next question.

## 1. What the reconstruction theorem actually assumes

[Busch, *Quantum states and generalized observables*](https://arxiv.org/html/quant-ph/9909073v3),
conditions P1–P3 and the unnumbered theorem: a probability assignment on
**all** effects `0 <= E <= I`, normalized at `I` and countably additive when effects
sum to at most `I`, is `v(E) = tr(rho E)`. This includes dimension two.
The ordinary projector theorem requires complex dimension at least three.

[Caves, Fuchs, Manne and Renes, *Gleason-Type Derivations*](https://arxiv.org/html/quant-ph/0306179v1),
§I and §III, Eq. (3)/Theorem 1: Hilbert-space effects and measurement
noncontextuality are assumed; normalization on every POVM yields the trace
rule. Their qubit result for continuous frame functions on rotated trines
(§IV.1–IV.2.1) also shows that the full effect set is sufficient, not the
only possible sufficient domain. One fixed trine is not that hypothesis.

For the mass problem, the missing work is to identify experimental outcomes
with such effects and justify the available outcome partitions, refinements
and operational equivalences. The following are different requirements:

| Existing requirement | What it does not by itself establish |
| --- | --- |
| L3, units; L16, physical frame covariance | availability of every qubit unitary or effect |
| L8, complete readings | every operator-valued outcome partition is an experiment |
| L17, equal laws give equal readouts | equality of effects across distinct measurement procedures |
| L19, no borrowing | an additive probability assignment on a noncommuting effect domain |
| L20, exact-zero neutrality | normalization of exhaustive raw-data outcome events |
| E6, shared calibration/evidence identity | noncommutativity, coherence or quantum interference |

This is a gap in the proposed bridge, not a contradiction in those premises.
In particular T70's Hilbert dimension is **two**, even for force and
acceleration in three-dimensional physical space. The two dimensions
cannot be exchanged to invoke the ordinary Gleason theorem.

## 2. Positive mass hypotheses are not the outcomes of one qubit POVM

On T70's one-instrument isotropic slice put

$$v_\phi=(\sin\phi,\cos\phi)^T,\quad
P_\phi=v_\phi v_\phi^T,\quad m=c\tan\phi,\quad 0<\phi<\pi/2.$$

Every `P_phi` has strictly positive off-diagonal entry. Thus no nonnegative
scalar weighting supported on these positive, finite mass rays can sum or
integrate to `I_2`: the trace condition gives total weight two, while the
off-diagonal entry of the sum is strictly positive. For example,

$$\int_0^{\pi/2}\frac4\pi P_\phi\,d\phi
=\begin{pmatrix}1&2/\pi\\2/\pi&1\end{pmatrix}\ne I_2.$$

The complement `I-P_phi` is the orthogonal ray `P_(phi+pi/2)`, corresponding
to signed slope `m'=-c^2/m`. Adding this effect is not merely selecting a
different positive mass from the current hypothesis set. Allowing the two
boundary axes would give the elementary two-outcome identity resolution;
it would not turn the interior positive mass alternatives into one POVM.

Consequently `H(m)/E = tr(rho P_m)` is a valid overlap for each test, but
is not already an outcome distribution over positive mass. A prior and
normalization over mass define an inference law, not this missing POVM.
This obstruction concerns this direct encoding; it does not exclude other
encodings or additional experimental effects.

## 3. The mass tests do not identify a pure complex state

Let `G` be T70's real, positive-definite Gram matrix and `E=tr(G)>0`.
The whole family

$$\rho_t=\frac1E\begin{pmatrix}
G_{11}&G_{12}+it\\G_{12}-it&G_{22}
\end{pmatrix},\qquad |t|\le\sqrt{\det G},$$

is positive semidefinite with trace one. For every real mass ray,

$$\operatorname{tr}(\rho_t P_\phi)=v_\phi^T Gv_\phi/E,$$

independent of `t`. The endpoints are pure complex-conjugate completions;
`t=0` is a real mixed state. For singular nonzero `G` the interval collapses;
the ambiguity is generic, not asserted for every reading.

For perpendicular equal-length vectors, `G=I`, both `rho_0=I/2` and the
pure state `psi=(1,i)/sqrt(2)` give probability `1/2` to every real ray.
A projector onto this complex `psi` distinguishes them, but is an extra
test absent from the present mass family.

T70's pure-completion identity `Q/E = geometric imaginarity(psi)` remains
exact. It does not make that resource operationally identifiable from the
mass tests: their equivalent real mixed state has zero imaginarity. A
physical pure-state interpretation needs an additional operational reason
to prefer that completion. This is not a refutation of T70's identities.

## 4. Binary normalization and covariance leave non-Born assignments

For `alpha>0` define

$$g_\alpha(p)=\frac{p^\alpha}{p^\alpha+(1-p)^\alpha},\qquad
q_\alpha(\rho,P)=g_\alpha(\operatorname{tr}(\rho P)).$$

On qubit rank-one projectors these assignments are continuous, monotone,
jointly unitarily covariant and context independent, with eigenstate
certainty and `q(P)+q(I-P)=1`. For `alpha != 1` this family contains
assignments with no representing Born density.
For a pure input, certainty fixes any possible representing density to
that pure state. At projector overlap `1/4`, `q_2=1/10`, rather than `1/4`.
These are countermodels to that short list, not rival mass definitions
satisfying all project premises.

Two ways to exclude this family require operations not established above:

- **Preparation mixing:** for `rho=.25 P_0+.75 P_1` with orthogonal pure
  states, `q_2(rho,P_0)=.1`, whereas mixing their assigned probabilities
  gives `.25`. Adding readings multiplies evidence; it is not by itself
  random preparation mixing of a system represented by the data qubit.
- **Effect refinement:** the displayed extension
  `F_alpha(E)=g_alpha(tr(rho E))` normalizes every binary POVM, but fails
  general refinements. Three real trine effects `E_j=(2/3)P_j` sum to `I`.
  For `rho=I/2`, `sum_j F_2(E_j)=3/5`; for a pure state along one of the
  trine rays it is `57/65`, rather than one. The maximally mixed
  *projector* assignment itself has an ordinary Born extension: this
  example refutes the displayed nonlinear effect extension, not every
  extension of those particular projector values.

Thus a square cannot be forced here merely by binary normalization plus
symmetry. A richer operational condition might force it; establishing
that condition for this measurement problem is the substantive task.

There is also missing information after normalization of the data state.
Readings `y=R(2,1)` have the same normalized T70 qubit for every `R>0`.
In the unit-Gaussian flat law, however, the log evidence contrast of mass
two against mass one is `R^2/4`: `.25`, `1`, `25` at `R=1,2,10`.
The law uses `exp(H/2)`, not `H/E`. Carrying energy separately repairs this
particular information loss; it does not select that response kernel,
the tempering, nuisance measure, prior or readout. At an exact-zero reading
`E=0` the normalized data qubit is undefined; the neutral evidence factor
is still well-defined. Normalized-state geometry alone does not supply
the estimator's null extension.

## 5. E6 has an exact commutative Born representation

For a finite mass/direction grid the joint E6 state `Q_ij` defines

$$\rho_Q=\sum_{i,j}Q_{ij}|m_i\rangle\langle m_i|
\otimes|n_j\rangle\langle n_j|.$$

For any declared decision probabilities `0<=r_ij<=1`, set
`E_r=diag(r_ij)`. Then `tr(rho_Q E_r)=sum Q_ij r_ij`: ordinary expectation.
Disjoint event refinement is additive. Partial trace over direction is
the mass marginal. The state is a tensor product iff mass and direction
are independent. All these states are separable and all these effects
commute, including E6's statistically dependent composition witness.

For **any** chosen kernel `0<L_ij<=1`,

$$\frac{\sqrt{E_L}\rho_Q\sqrt{E_L}}
{\operatorname{tr}(\rho_QE_L)}
=\operatorname{diag}\!\left(\frac{Q_{ij}L_{ij}}{\sum Q_{ij}L_{ij}}\right).$$

This diagonal Lüders formula represents ordinary Bayesian reweighting.
It holds for `L`, `sqrt(L)` and `L^2`, among many other choices. Any bounded
positive kernel can be scaled into this interval without changing the
update. Its Born
form cannot select `L`. It also supplies no reason to identify this state
with T70's qubit: T70 encodes observed data and tests a mass ray; `rho_Q`
encodes an already chosen joint law over unknowns.

Partial trace preserves every mass-local effect `E_M tensor I`. It need
not preserve a future fixture-dependent effect. The mass-neutral E6
history in [§14](e6-declared-structure.md#14-a-mass-neutral-history-can-still-change-a-future-mass-update)
has the same mass marginal as no history but a different direction law;
a later reading therefore produces different mass evidence. The new
checks verify both the preserved local expectations and a lost joint
expectation. That expectation of a chosen profile factor is not yet a
predictive probability of an allowed raw-data event. E6 established a
classical reason to carry addressed joint
state; labeling this a partial trace does not add interference.

## 6. A measurement experiment needs effects for raw-data events

The formula in §5 makes any bounded kernel one valid binary effect with
its complement. It does **not** make the collection of all possible
observed-data kernels a normalized measurement of that data.

For a finite grid of full hypotheses `h=(m,n,t_1,...,t_N)`, the specified
normalized Gaussian sampling model does supply such a measurement:

$$E_A=\operatorname{diag}_h\{\Pr_h(Y\in A)\}.$$

These commuting effects are assigned before data. For every exhaustive
disjoint raw-data partition `{A_j}`, `sum E_(A_j)=I`, with additivity under
coarse-graining. This follows pointwise in the full hypothesis and needs
no probability population of signals. A probability law on that full
hypothesis would be a separate commitment. Continuous data require
events/bins or an operator-valued measure; point densities need not lie
between zero and one.

Profiling a signal is not automatically a reduction preserving those
event probabilities. In the whitened scalar example, write

$$u_m=(m,1)/\sqrt{1+m^2},\quad v_m=(1,-m)/\sqrt{1+m^2},
\quad y=t u_m+r v_m.$$

Then `K_m(y)=r^2`. The profile factor `exp(-K_m/2)` is constant along the
signal line. Its raw-data Lebesgue integral over `|t|<=A` is
`2A sqrt(2pi)`, diverging with `A`. It is not a normalized raw-data density
of the original experiment. T48 already distinguishes a normalized
candidate-residual density from a normalized full-data experiment, and
E5 already exhibits nuisance-measure ambiguity. Moreover, `r_m(y)` changes
its data-to-residual map with the candidate: a residual bin is not one
common raw-data event as `m` varies. This is the operational
connection to those findings, not a new general impossibility theorem.

**Important limit of that argument.** Nonintegrability does not exclude
every normalized family having these likelihood ratios. For a fixed-rank
candidate projector `P_h` in dimension `D`, multiply the profile by the
data-only weight `a(y)=exp(-||y||^2/(2 tau^2))`, `tau>0`. Then

$$q_h(y)=C\,a(y)\exp[-y^T(I-P_h)y/2]$$

is a normalized zero-mean Gaussian, with the **same** `C` for every `h`
and covariance

$$\tau^2P_h+\frac{\tau^2}{1+\tau^2}(I-P_h).$$

It has the same relative likelihoods as the profile. But its residual
variance is less than the declared unit noise variance, so it is a
different sampling experiment. This construction, checked independently,
shows why the declared probability model must do the work: neither
nonintegrability nor a formal Born representation settles it.

The E6 fixed-linear GLS compression (§8 of the parent exploration) is a
positive conditional example of a sufficient statistic for a declared
Gaussian experiment. Its constraint provenance remains to be justified.
The unknown-line sufficiency result (§10) concerns the full mean experiment;
it does not forbid every mass-specific or nonlinear reduction. Likewise
Born normalization on the full `h` does not decide how to eliminate `t`,
what measure to integrate against, or which calibration may be learned
under L19.

## 7. Sort and next direction

**Provisional assessment, Codex.** The live T70 work is relevant to the
estimator's geometry. E6 is relevant to the identity and scope of the
state carried between operations. A Born uniqueness argument remains
unestablished for this type. The direct positive-mass POVM identification,
unique pure completion, and sufficiency of binary normalization each
have explicit counterexamples above. None is a proof that the sought
mass definition does not exist.

The useful next fork is: **which reductions of a declared raw-data
experiment preserve all probabilities relevant to later allowed
measurements, without installing a signal-size population?** First state
the event domain and allowed operations; then seek an exact reduction or
a counterexample. Two histories count as equivalent only if every allowed
future experiment has the same relevant probabilities; seek a state
that preserves this equivalence. The candidate-dependent residual map
must not silently stand for a common physical outcome. Fixed-linear GLS
and E6's mass-neutral history provide
positive and negative controls. The answer must also account for T30,
Z and L19; the Gaussian joint construction is still not their survivor.
An effect algebra or a nuisance population is not automatically a
permissible new input under TASK §5 merely because a reconstruction
would be easier with it.

This does not redo the uncompleted B22 Gram-shrinkage coverage work,
B23 real/scalar-inertia argument, or E7–E9 alternatives and geometry.
Those remain separate tests of which ingredients the physical type
supplies. This continuation makes no claim that any of them has been
completed or rendered unnecessary.

## Validation

The 18 new checks pass. Together with T70's reading-qubit checks, E6's
declared-structure checks, graph geometry, normalized residual-mixture
and audit regressions: **138 passed**. Independent mathematical/source
review found no blocker after the finite-grid, countable-additivity and
special-state scope clarifications. The new tests check the explicit
calculations and counterexamples; they do not numerically prove the cited
reconstruction theorems or test all possible operational models.

`tools/audit.py` retains its existing T6/T21/T27 and unused-L24 findings;
`git diff --check` passes. The full default suite was not repeated for
this documentation/test-only continuation; the earlier full-suite result
in the parent exploration belongs to its first commit.
