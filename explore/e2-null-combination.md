# E2: null evidence beside an independent mass determination

**Status:** exploration, 17 Nov 2026, Codex; branch `explore/e2-null-combination`,
from main `08d36f0`. Conditional obstructions and a score-combination benchmark,
not a new premise, estimator, or interpretation adopted for property N.
**Checks:** `tests/test_null_combination.py`, `tests/test_null_combination_score.py`,
and the existing `tests/test_washout_fixed_history.py` for the E10 history bound.

## 1. The fork and its standard mathematics

This concerns **combining independent evidence about the same object's mass**.
It is not L5/L14's addition or parallel addition of different masses. The
external experiment may represent another measurement route; the proof needs
only its distribution and the baseline result, not a density/volume model.

R-2026-09-29-3 says the 0/0 result should carry zero weight beside an independent
density result. R-2026-09-29-4 asks for a true uncertainty statement on noise;
R-2026-09-30-3 permits that statement beside the law. The candidate sharpening
tested here is that **adding a growing pure-noise block changes the external
answer by an amount tending to zero in probability**. That sharpening is not
already a ruling. In particular, distinguish:

1. Exact observed zeros: Z's algebraic invariance.
2. True zero signal declared known: its likelihood is independent of mass.
3. True zero signal in the investigator's world, not disclosed to the procedure:
   random Gaussian observations. This is the experiment below.

The standard method is the **constrained-risk inequality**, applied to departure
from a random external baseline. Duchi and Ruan, *A constrained risk inequality
for general losses* (2020), [§2, Theorem 1 and Corollary 1; §3.1, Proposition 1](https://arxiv.org/pdf/1804.08116),
give the likelihood-ratio second-moment method. E10 already applies it to a
returned law. Here the external result is kept random throughout the proof.

For the constructive comparison, T45 is **conditional score-test inversion**;
see Moreira (2003), [§2.2, Theorem 1 and Example 2, pp. 1031–1032](https://drphilipshaw.com/Moreira-2003.pdf).
Sampling calibration is not a likelihood interpretation. Even for genuine
confidence distributions, normalized multiplication is not generally calibrated:
Singh, Xie and Strawderman (2005), [§4, Theorem 4.1 and Example 4.1, electronic-reprint pp. 13–14](https://arxiv.org/pdf/math/0504507).
T45's two-sided p-value curve is not itself a CDF in mass.

## 2. A precise conditional contract

Fix one positive mass m, a declared prior and the scalar common-instrument
Gaussian slice. Independently observe

\[
 X\sim R_m,\qquad Y_i=r_i u_m+Z_i,\quad Z_i\sim N(0,I_2),\qquad
 u_m=(\sin\theta_m,\cos\theta_m),\quad\theta_m=\arctan(m/s).
\]

The signed signals r_i are fixed, not draws from a declared population. The
external experiment X has fixed precision as N grows. Let B=b(X) be its baseline
positive point and U_N=U_N(X,Y_1,\ldots,Y_N) the combined point. Neither must be
a median. Assume for this m that for some epsilon>0

\[
 \delta=R_m\{|B-m|>3\epsilon\}>0. \tag{1}
\]

This is the ordinary nonexact external case. For example B=exp(X), with
X Gaussian about log m, is positive and has nonzero error. More generally, if
the external distributions at two distinct masses are equivalent, a baseline
cannot be exactly right almost surely at both; (1) holds at at least one mass.

The two candidate requirements are:

- **W_ext:** under R_m times the pure-null Y law, U_N-B tends to zero in probability.
- **C_ext:** for every deterministic sequence of signal experiments with
  Delta_N=sum r_i^2 and Delta_N^2/(N+Delta_N) tending to infinity,
  U_N tends to m in probability under the joint experiment.

C_ext is the proposed extension of L11 to a combination with an independent
route. It is **additional to** the currently written scalar L11 about the
estimator's own point. It would follow from requiring that adding a fixed
ordinary measurement cannot spoil an identifying determination. The proof
does not silently put external observations into the current input type.

**Internal pooling corollary.** If X is instead a fixed finite block of the
same scalar reading type and U_N is the estimator on the pooled readings,
no extension of L11 is needed. Adding the fixed count and fixed energy of X
preserves the identifying condition. Its finite Gaussian laws are equivalent
at different masses, so for any baseline, (1) holds at at least one of those
masses. Thus an L11-consistent
point cannot also leave every fixed informative block unchanged in probability
under growing random true-null padding. Z still concerns exact observed zeros.

## 3. Finite bound and the slowly identifying counterexample

Set A_N={|U_N-B|>epsilon} and alpha_N=(R_m times P_0)(A_N). For a proof mixture
give every reading squared radius t=Delta/N and average independent fair signs.
The Y likelihood ratio relative to pure noise is

\[
 L(Y)=\prod_i e^{-t/2}\cosh(\sqrt t\,u_m^TY_i),\qquad
 E_0L^2=\cosh(t)^N\le e^{\Delta^2/(2N)}. \tag{2}
\]

The external distribution is **the same** in reference and alternative, so
the joint likelihood ratio is still L(Y). Cauchy–Schwarz gives

\[
 P_{\rm mix}(A_N)\le\sqrt{\alpha_N}\,e^{\Delta^2/(4N)}. \tag{3}
\]

This applies even when A_N depends on both experiments. On |B-m|>3epsilon
and A_N's complement, |U_N-m|>2epsilon. Consequently

\[
 P_{\rm mix}(|U_N-m|>2\epsilon)
 \ge \delta-\sqrt{\alpha_N}\,e^{\Delta^2/(4N)}. \tag{4}
\]

W_ext gives alpha_N to zero. When 0<alpha_N<1 choose

\[
 d_N^2=\min\{\log(1/\alpha_N),\sqrt N\},\qquad
 \Delta_N=\sqrt N\,d_N.
\]

Then d_N diverges, d_N/sqrt N tends to zero and
Delta_N^2/(N+Delta_N)=d_N^2/(1+d_N/sqrt N) diverges. Nevertheless the
departure bound (3) is at most alpha_N^(1/4), tending to zero. At indices
alpha_N=0, absolute continuity makes departure impossible under every finite
alternative; use d_N=N^(1/4) there. Initial indices do not matter.

At each N choose a deterministic sign vector whose departure probability is
no larger than the mixture's. Those vectors form one identifying triangular
array on which the combined point's error probability has liminf at least
delta. Thus **W_ext and C_ext cannot both hold**. No law concentration,
particular estimator, Z, no-borrowing rule, or external error moments were used.
The bad array depends on the proposed procedure; it is not one array defeating
all procedures. The argument gives failure along one array, not merely a failure
of a uniform error bound over experiments.

### One fixed infinite history, with reading-order invariance

The fixed-history strengthening of E10 carries over, provided U_N is invariant
to permuting the Y readings while X remains distinguished. For completeness:
write a_i=r_i^2 and average independent signs and permutations of a prefix.
For a_i<=j/sqrt i the likelihood-ratio second moment obeys

\[
 \log E_0L_{\rm sym}^2\le
 2j^2+\frac{e^{j^2/2}j^4(1+\log N)^2}{8N}. \tag{5}
\]

The proof and its negative-association source are recorded in O15's 9 Nov
fixed-history intake and on [the E10 branch](https://github.com/canvalidk/Mass-estimator-testor/blob/explore/e10-null-limit-prior/explore/e10-null-washout-consistency.md#4-a-fixed-history-obstruction-when-readings-are-permutation-invariant).
Its finite ingredients already have independent tests on main.

Choose increasing endpoints N_j so that the remainder in (5) is at most 1,
alpha_(N_j) exp(2j^2+1)<=j^(-4), and N_j>=(j+1)^4. Define a_i=j/sqrt i
for N_(j-1)<i<=N_j. The previous-stage coefficients are at most j, so (5)
applies at N_j. Permutation invariance makes the fixed-order sign-mixture
departure probability equal to its symmetrized version, hence at most j^(-2).

This single amplitude history has a_i to zero and
Delta_N/sqrt N >= (3/4+o(1)) times the stage index at ceil(N/4), which diverges.
Its mixture error at the endpoints is at least delta-j^(-2). If every fixed
infinite signed version were consistent, dominated convergence over **one
product measure of infinite signs**, and the independent X law, would make that
average error tend to zero. Therefore some fixed signed identifying history
fails too. The permutation mixtures need not be nested; order invariance
identifies them with the one fixed sign mixture separately at each prefix.

### Stronger and weaker notions

- At finite N, exact almost-sure neutrality under random true-null data transfers
  by Gaussian absolute continuity to every signal vector. It discards the block
  entirely. This says nothing against exact *observed*-zero invariance.
- If the combined law converges in total variation to the random external law,
  and the latter has an almost-sure unique strict-crossing median B, its median
  is neutral too. For any epsilon, the baseline CDF's gaps at B±epsilon are
  positive almost surely; first bound them below by gamma, then let gamma fall
  to zero. Thus using laws does not evade the result for that readout.
- Neutrality of just an auxiliary field, with no consequence for the combined
  point, is not excluded. Neither is an averaged/unbiased combination promise.
- An external determination whose own precision grows, an exact oracle, or
  consistency only above an additional signal-growth margin is outside this
  fixed-external/full-rate obstruction.

## 4. A separate inverse-weight field faces the same detection boundary

Suppose the Y block reports a point t_N and a nonnegative finite precision w_N,
both computed without X. Let the fixed external report have a nonatomic point
B and fixed weight w_0>0. Consider the declared combination rule

\[
 U_N=\frac{w_0B+w_Nt_N}{w_0+w_N}=a_NB+(1-a_N)t_N,
 \qquad a_N=\frac{w_0}{w_0+w_N}. \tag{6}
\]

This is a mathematical consequence for inverse weighting, not a declaration
that any proposed w is the true inverse sampling variance. A positive continuous
external B suffices; a Gaussian raw mass observation is not required.

Let Q_B(h)=sup_x P(B in [x,x+h]) be the concentration function. Nonatomicity
implies Q_B(h) to zero as h falls to zero. Conditional on Y, on a_N>=epsilon,

\[
 P(|U_N-m|\le d\mid Y)\le Q_B(2d/\epsilon).
\]

Therefore

\[
 P(|U_N-m|>d)\ge
 P(a_N\ge\epsilon)[1-Q_B(2d/\epsilon)]. \tag{7}
\]

Point consistency forces a_N to zero and hence w_N to infinity along every
identifying experiment. But w_N to zero under random null data and w_N to
infinity along **every** identifying experiment are incompatible: apply (3)
to A_N={w_N>w_0} and choose the same slow alternatives. The fixed-history
version holds if w_N is invariant to reading order.

Also, w_N to zero alone does not ensure U_N-B to zero: the product w_N t_N
may survive. For example w_N=1/N and t_N=N gives U_N-B tending to 1/w_0.
An operational zero-influence requirement concerns both weighted numerator
and denominator. The detection obstruction above does not need this stronger
null-neutrality assumption; the precision claims already conflict.

This does **not** prove that moving uncertainty beside the law was wrong.
It excludes one added promise for that field combined with (6) and C_ext.

### A positive control for a weaker null promise

There is a simple compatibility witness for a **tight, nonvanishing** null
weight. In the scalar common-instrument slice let G=sum Y_iY_i^T, E=tr G,
H=lambda_max(G)-lambda_min(G), and w=H^2/E, with w=0 at E=0. For finite data
0<=w<=E. Zero padding and permutations leave it unchanged.

For arbitrary deterministic signals, Gaussian moments give

\[
 EG=NI_2+\Delta uu^T,\quad E\|G-EG\|_F^2=6(N+\Delta),\quad
 E[E]=2N+\Delta,\quad\operatorname{Var}(E)=4(N+\Delta).
\]

Under Delta^2/(N+Delta) tending to infinity, eigenvalue perturbation gives
H/Delta to one, while E/(2N+Delta) tends to one. Therefore
w is asymptotic to Delta^2/(2N+Delta), and diverges. At pure null, the
two-dimensional scatter CLT in (9) gives w converging in distribution to
2 chi-square_2. The null weight is tight, but does not tend to zero.

Combine this weight with a positive finite consistent point, for example
T42's flat Bartlett median in its stated Gaussian scope, and a fixed external
positive finite point. Formula (6) is then consistent: the external coefficient
tends to zero, and both point errors are tight. This is an attainability
control for the weakened weight promise only. A common change of mass units
can assign w/c^2 and w_0/c^2 the same inverse-mass-squared dimension using the
declared prior centre c; no data-derived physical variance is asserted.
There is no calibration, L19, T30 or uniqueness claim and no proposed default.

## 5. What the exact score construction does when combined

T45 supplies, for each candidate mass, a signed score S and its conditional
variance C. For independent blocks, express these in one common parameter
coordinate, retain the curves and add:

\[
 (S,C)_{A\cup B}=(S_A+S_B,C_A+C_B),\qquad
 T_{A\cup B}^2=\frac{(S_A+S_B)^2}{C_A+C_B}. \tag{8}
\]

Under the fixed candidate's true declared Gaussian experiment, conditioning
on all parallel projections leaves independent normal scores with these
variances. Hence (8) preserves exact pointwise coverage, is associative,
order-independent and Z-safe. It uses each reading's own projection; it
does not learn a signal population. Combining another measurement route
would additionally require that route to supply the analogous justified
score and variance. A mass interval alone does not do so.

**Even the entire p-value curve is insufficient for this merger.** On the
scalar isotropic slice take two possible blocks with Grams I_2 and 3 I_2.
Both have S(theta)=0 at every theta, so their p-value curves are identically
one and their confidence sets are the whole domain at every level. Their
flat Gaussian profile factors are constant in mass too, so with the same
prior their normalized raw laws agree. Add the
same reading y=(3,-1). At m=1 its score is 4 and variance 2. The combined
squared statistics are 16/3 and 16/5: one rejects and the other accepts at
95%. The state lost by the set is consequential, even though the original
sets agree at *all* levels. This proves insufficiency for (8), not a universal
minimal-state theorem or a prohibition on other set-combination rules.

**Exactness and these invariances leave a family.** At each candidate let
h_i=h(||U_i||^2)>0 be a bounded deterministic function of that reading's own
whitened parallel projection. Conditional on all projections,

\[
 \frac{\sum_i h_i S_i}{\sqrt{\sum_i h_i^2 C_i}}
\]

is again exactly standard normal at truth. Choosing h=1 or h(t)=1/(1+t)
preserves Z, reading order, frame/units covariance and candidate-coordinate
invariance, and does not borrow another reading's signal. They give different
sets on reachable data: the tests exhibit squared pivots 0.1 and 8.5 at the
same candidate. On the scalar unit-noise angle slice C_i=||U_i||^2, but the
weight's full-type definition uses the projection norm, not the
coordinate-dependent C_i. This establishes nonuniqueness under those stated
conditions only; it makes no T30, L11, proper-law or full property-N claim.

### A growing true-null block overwhelms a fixed informative block

For common unit isotropic noise put G_N=sum Y_iY_i^T, u=(sin theta,cos theta),
v=(cos theta,-sin theta). The signed statistic is S/sqrt C with
S=u^TG_Nv and C=u^TG_Nu. Under pure noise,

\[
 G_N=NI_2+\sqrt N A_N,\qquad
 \left(\frac{G_{11}-G_{22}}{2\sqrt N},\frac{G_{12}}{\sqrt N}\right)
 \Rightarrow (Z_1,Z_2),\quad Z_1,Z_2\text{ independent }N(0,1).
\]

The independent fixed block contributes finite S_X(theta), C_X(theta).
Uniformly on the closed angle interval [0,pi/2], its contributions vanish
after dividing score by sqrt N and variance by N. The combined signed-score
curve converges in distribution to

\[
 Z(theta)=Z_1\sin(2theta)+Z_2\cos(2theta), \tag{9}
\]

independently of X. This is convergence of continuous statistic curves;
no stronger random-set topology claim is needed. At any fixed candidate
the limiting acceptance probability at nominal p is p, and it has forgotten
X. If X alone accepts that candidate with probability q_X, the asymptotic
probability that the combined and X-only decisions differ is

\[
 p(1-q_X)+(1-p)q_X\ge\min(p,1-p)>0. \tag{10}
\]

At a true candidate with an ordinary nonempty T45 external block, q_X=p,
so disagreement tends to 2p(1-p), or 9.5% at p=0.95. Coverage remains correct;
neutrality fails. This is an independence limit, not an assertion of exact
finite-sample independence of the two decisions.

The limit accepts the **whole positive mass domain** iff Z_1^2+Z_2^2<=z_p^2,
where z_p=Phi^(-1)((1+p)/2). Its probability is

\[
 1-e^{-z_p^2/2}, \tag{11}
\]

about 0.8535 for p=0.95, not one. Thus exact coverage does not make the
reported set become wholly uninformative with growing pure-null data.
The limit set is nonempty almost surely because the sinusoid has a zero
inside the positive angle interval. Fixed-candidate coverage, whole-domain
acceptance, and zero influence on another result are different statements.

## 6. What this exploration decides

| Candidate requirement or construction | Finding | Scope |
| --- | --- | --- |
| Random-null zero influence plus full-rate combined point consistency | Incompatible | Fixed nonexact external result; scalar Gaussian arrays; fixed histories with order invariance |
| Separate inverse weight vanishes at null and dominates a fixed external result on every identifying series | Incompatible | Independent nonatomic external point, fixed positive external weight; no moment assumption on the block point |
| Tight, nonvanishing null weight and full-rate identifying dominance | Compatible | Scalar Gram control H^2/tr G; not a selected uncertainty rule |
| T45 score/variance curves carried and added | Exactly calibrated, associative and Z-safe | Declared independent Gaussian worlds; candidate-specific conditional test, not a proper mass law |
| T45 set or p-value curve alone as state for that merger | Insufficient | Explicit two-Gram counterexample |
| Exact conditional coverage, Z and the stated invariances select the score set | No: a family remains | Own-projection weighted score controls; no claim about the whole premise list |
| T45 merger gives a null block zero influence | Fails | Null block asymptotically erases the fixed informative block in the score statistic |

The provisional sort is a **candidate wrong-premise boundary**, not an input
or a selected repair: interpreting the 0/0 combination ruling as W_ext over
random true-null histories cannot coexist with C_ext at the full identifying
rate. The existing text does not force that interpretation or C_ext's
cross-route extension, so no recorded premise is struck.

G28 now has a sharper question: which operational promise is intended for
the output beside the law—finite-error control of its influence, a
signal-margin-qualified vanishing influence, or another statement? Exact
coverage alone does not answer it. A proposed promise must specify the
combination operation and its retained state as well as its null limit.
G3's finite-N form and the rest of the production equation remain open.

## Validation

On 17 Nov the final two new test files passed **40 tests**. Before the last
three weighted-score checks were added, the focused run with T45 and E10's
existing history checks passed 83 tests. Independent agent reviews of the
obstruction, inverse-weight argument and score limit found no mathematical
blocker; finite tests do not replace the written asymptotic arguments.

The full suite, which collected the initial 37 new tests, finished with
2,847 passed, 28 skipped, 25 expected failures and five failures. All five
were Windows `WinError 5` when creating worker-process communication pipes,
before numerical comparisons. The two composition worker tests then passed
with process permission, and all five tests in `tests/test_parallel.py`
passed with that permission (including the three affected cases). The final
40-test run includes the three subsequently added checks. No unresolved
test failure remains. The structural audit retains its four existing
findings: T6, comparison-only T21, T27, and unused L24. `git diff --check` passes.
