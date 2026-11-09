# E10 continuation: null washout versus slow identification

**Status:** exploration, 9 Nov 2026, Codex. A conditional compatibility
boundary, not a change of L11, Z, the prior, or a production equation.
**Checks:** `../tests/test_null_washout.py`; research helpers in
`../met/null_washout.py`. This follows `e10-null-limit-prior.md`.

## 1. Name the proposed requirements

Use the scalar Gaussian slice of T41--T53: independent
Y_i=r_i u_theta+Z_i, Z_i~N(0,I2), u_theta=(sin theta,cos theta),
theta in (0,pi/2). All covariances, reading windows and other declarations
are identical, with one fixed full-support prior pi. Signs and radii are
fixed nuisance quantities, not a declared population. Write

\[
\Delta_N=\sum_{i\le N}r_i^2,\qquad
\mathcal I_N=\frac{\Delta_N^2}{N+\Delta_N}.
\]

T42's consistency threshold I_N->infinity is equivalent, as N grows,
to Delta_N/sqrt(N)->infinity. T53 establishes its energy-class meaning.

**W (candidate null washout):** under r_i=0, the returned law Q_N tends to
the supplied pi in total variation in probability. This implies asymptotic
averaged neutrality E0 Q_N->pi and is *stronger* than that averaged limit
requirement (E10 tested an identity for its limiting random law, not exact
finite-N neutrality). Neither is already implied by
exact measured-zero invariance Z. The lower bound needs only one event
probability Q_N(B)->pi(B), not full total variation.

**C (law concentration):** Q_N gives probability tending to one to every
fixed neighbourhood of the true positive mass whenever I_N diverges.
We separately specify triangular arrays and fixed infinite histories below.
Current L11 is a **point** consistency statement, not automatically C.
Section 5 gives the implication for the current median readout.

**P (explicit scope assumption):** reordering whole readings with identical
declarations leaves the output unchanged. This holds for the summed/product
probes used here. We use it explicitly, not as a new ruling attributed to
L19. The type lists reading membership and windows; it does not separately
declare a chronological weighting protocol. Whether P is forced by that
type is a premise interpretation, not a statistical theorem.

## 2. Standard mathematics

This is a **constrained-risk/contiguity argument**: making a statistic settle
at one experiment limits its response to nearby experiments. Duchi and Ruan,
*A constrained risk inequality for general losses* (2020), section 2,
Theorem 1 and Corollary 1, use the likelihood-ratio second moment and
Cauchy--Schwarz; section 3.1, Proposition 1 shows the corresponding slow-local-
alternative effect in Gaussian mean estimation:
[primary paper](https://arxiv.org/pdf/1804.08116).
We apply that method to a **reported-law functional**, not to a supposedly
identified true mass at the null. The elementary model-specific proof is
below. Duchi's [2021 notes](https://stanford.edu/class/stats300b/Notes/contiguity-and-asymptotics.pdf),
section 1.2, define contiguity and its transfer of events of vanishing probability.

T41/T53 already supply the Gaussian sign-mixture identities. The fixed-history
extension additionally uses **negative association of sampling without
replacement** to control a permutation mixture. The source and finite bound
are recorded in section 4. These are proof devices; no physical distribution
of signal sizes, signs, order or setup is passed to the estimator.

## 3. The triangular-array obstruction and a finite bound

Choose an interior true angle and a neighbourhood B with pi(B)<1. Let

\[
A_N=\{Q_N(B)>(1+\pi(B))/2\},\qquad \alpha_N=P_0(A_N).
\]

W implies alpha_N->0. C would imply P_theta,r(A_N)->1 for identifying
arrays. For the proof, average independent fair signs at homogeneous
squared radius t=Delta_N/N. Every mixed signal vector has exactly that
energy and the same true mass. T41's likelihood-ratio identity gives

\[
E_0 L_N^2=\cosh(t)^N\le e^{\Delta_N^2/(2N)},\qquad
P_{\rm mix}(A_N)\le\sqrt{\alpha_N}\,e^{\Delta_N^2/(4N)}.
\]

The bound is independent of the estimator, Z, and any proposed tempering.
For 0<alpha_N<1 choose

\[
d_N^2=\min\{\log(1/\alpha_N),\sqrt N\},\qquad
\Delta_N=\sqrt N\,d_N.
\]

Both entries in the minimum diverge, so d_N->infinity. Also d_N<=N^(1/4),
and I_N=d_N^2/(1+d_N/sqrt N)->infinity. Nevertheless
P_mix(A_N)<=alpha_N^(1/4)->0. If alpha_N=0, finite Gaussian alternatives
are absolutely continuous with respect to the null, so their event
probabilities are zero too; take d_N=N^(1/4) at those indices.

At each N some deterministic sign vector has event probability at most
the mixture average. Thus W and C cannot both hold for **every identifying
triangular array**. This supplies one bad array for each proposed procedure,
not one fixed array defeating every procedure. It is not by itself a
fixed-history proof, nor an assertion that L11 already has these quantifiers.

## 4. A fixed-history obstruction when readings are permutation invariant

The additional assumption P lets us prove the stronger statement: for
each procedure with W and P, there is one fixed signed infinite history
with Delta_N/sqrt N->infinity on which C fails.

Put a_i=r_i^2. Average the fixed amplitudes over independent fair signs and
a uniform permutation of their N positions. Relative to pure noise,
expanding the product of two mixture likelihood ratios gives exactly

\[
E_0 L_{\rm sym}^2
=E_\sigma\prod_{i=1}^N\cosh(\sqrt{a_i a_{\sigma(i)}})
\le E_\sigma\exp\left(\tfrac12\sum_i a_i a_{\sigma(i)}\right).
\]

For the random permutation of a fixed nonnegative population, negative
association yields, with increasing nonnegative f_i(x)=exp(a_i x/2),

\[
E_0 L_{\rm sym}^2
\le\prod_i\left[\frac1N\sum_k e^{a_i a_k/2}\right].
\]

The primary verification is Karayel, *Negatively Associated Random Variables*,
[AFP proof document](https://isa-afp.org/browser_info/current/AFP/Negative_Association/document.pdf):
section 6, printed p.62, covers uniform permutations of populations with
possible repeated values; section 2, Property P2, printed pp.17--18, gives
the product inequality for nonnegative integrable monotone functions.
Finite populations make our exponentials integrable. The formalized proof
follows Dubhashi's FKG argument; see also Karayel et al. (2025), section 6,
printed pp.34:12--34:13, in their
[primary paper](https://drops.dagstuhl.de/storage/00lipics/lipics-vol352-itp2025/LIPIcs.ITP.2025.34/LIPIcs.ITP.2025.34.pdf).

For a_i<=b_i=j/sqrt i, monotonicity and
exp(x)<=1+x+exp(j^2/2)x^2/2 for 0<=x<=j^2/2 give

\[
\log E_0 L_{\rm sym}^2
\le\frac{(\sum_i b_i)^2}{2N}
 +\frac{e^{j^2/2}(\sum_i b_i^2)^2}{8N}
\le 2j^2+\frac{e^{j^2/2}j^4(1+\log N)^2}{8N}.
\tag{*}
\]

For each fixed j the last remainder vanishes. Now choose increasing endpoints
N_j so large that the remainder in (*) is at most 1 and
alpha_(N_j) exp(2j^2+1)<=j^(-4). Such endpoints exist by W. Also require
N_j>=(j+1)^4. Define one infinite deterministic amplitude history by

\[
a_i=\frac{j}{\sqrt i}\quad\text{for }N_{j-1}<i\le N_j,
\qquad N_0=0.
\]

At the endpoint N_j its entire prefix satisfies a_i<=j/sqrt i. P makes
the probability of A_(N_j) under the fixed-order sign mixture equal to
its probability under the symmetrized mixture. Cauchy--Schwarz and (*)
therefore bound it by j^(-2).

The stage coefficient kappa_i=j is nondecreasing and diverges. Consequently

\[
\frac{\Delta_N}{\sqrt N}
\ge (3/4+o(1))\,\kappa_{\lceil N/4\rceil}\longrightarrow\infty.
\]

Moreover a_i->0, by the endpoint growth requirement: all the individual
signals eventually become weak, while total identification grows.

Finally use **one fixed product measure on infinite sign streams**. If C
held on every deterministic signed version of this amplitude history,
bounded convergence over this sign measure would give P_sign(A_N)->1.
The subsequence bound j^(-2) contradicts that. Thus some one fixed sign
history fails. The permutation mixtures themselves need not be consistent
across N: P supplies their equality to the fixed sign mixture separately
at each N. This is the step that makes the fixed-history quantifiers valid.

No exact-zero premise or no-borrowing restriction was needed for this
obstruction. Adding those restrictions cannot rescue W+C+P. This still
does not make the newly proposed W an existing project requirement.

## 5. Connection to the point statement and why P matters

For a prior with a unique median c and strict CDF crossings on either side
of c (including the half-Cauchy), W forces every returned median to tend
to c at the null. Choose a true mass m0 different from c and use
A_N={|median(Q_N)-c|>|m0-c|/3}. Null probability tends to zero and point
consistency at m0 would make its probability tend to one. The preceding
proofs apply unchanged. Thus the permutation-invariant fixed-history
obstruction reaches the **current median's point consistency**, not merely
the stronger law-concentration promise. An arbitrary discontinuous readout
is outside this implication.

Without P, a fixed-history positive control exists. Define

\[
C_N=\sum_{i=1}^N\frac{Y_iY_i^T-I_2}{\sqrt i},\quad
\ell_N=\log(N+1),\quad
Q_N^{\rm ord}(d\theta)\propto\pi_\theta(d\theta)
\exp[u_\theta^TC_Nu_\theta/\ell_N].
\]

Every finite output is a proper law. Write
A_N^*=sum r_i^2/sqrt i, h_N=sum 1/i, B_N^*=sum r_i^2/i<=A_N^*.
Gaussian moments and independence give

\[
EC_N=A_N^*u_0u_0^T,\qquad
E\|C_N-EC_N\|_F^2=6(h_N+B_N^*).
\]

At null C_N/ell_N->0 in L2, hence TV washout. For a fixed identifying
history, summation by parts gives

\[
A_N^*=\Delta_N/\sqrt N+
\sum_{i<N}\Delta_i(1/\sqrt i-1/\sqrt{i+1}),
\qquad A_N^*/\ell_N\to\infty.
\]

The error C_N-EC_N is o_p(A_N^*). Its leading line tends to the true
line and its scaled eigenvalue gap diverges, so full support of pi gives
concentration at every fixed interior mass. Thus W and fixed-history C
are statistically compatible if chronological weights are permitted.
This control fails P and exact-zero padding: appending zero changes ell_N
(its additional scalar identity term cancels, but the denominator does not).
It is a countercontrol for the quantifiers, not a licensed estimator.

## 6. A zero-padding-safe control with an explicit additional margin

Let E=tr(sum Y_iY_i^T), and let K be T42's flat Bartlett precision
min(H,H^2/J), with its stated fallback. As a comparison only, define

\[
b(E)=\log(e+E),\qquad K^*=\frac{K^2}{K+b(E)},\qquad
Q_N^*(d\theta)\propto\pi_\theta(d\theta)
\exp[-K^*\sin^2(\theta-\hat\theta)/2].
\]

This continuous attenuation uses no true amplitudes and no declared signal
population. It is **unforced** and its use of total observed energy needs
a separate L19 justification. It does not reject any reading or claim to
certify the Newtonian model. It is not proposed as a production replacement.

At null K=O_p(1), E/N->2 and K*->0, giving W. If instead

\[
\frac{\mathcal I_N}{\log(e+N+\Delta_N)}\to\infty,
\]

T42 gives H/Delta->1 and J=O_p(N+Delta), while
E/(2N+Delta)->1. It follows that K/b(E)->infinity, K*/K->1,
and K*->infinity. The fitted line is consistent and the law concentrates.
This includes repeated fixed nonzero signals and fixed-N increasing
amplitudes. It retains Bartlett's known calibration limitations; relative
precision agreement is not a new coverage theorem.

Padding exact zeros leaves E,H,J and K unchanged. Reordering readings
also leaves them unchanged. Since 0<=K*<=K<=H, the law is continuous even
at eigenvalue ties. But the logarithmic schedule has no forcing argument;
other diverging schedules create other tradeoffs.

The sacrificed regime is explicit: take homogeneous triangular signals
with Delta_N=sqrt(N)(log N)^(1/4). Then I_N~sqrt(log N)->infinity,
H/Delta_N->1 and J_N/N->1 (the homogeneous weak-signal Gaussian
polynomial-moment law of large numbers, uniformly in angle). Hence
K~sqrt(log N), b(E)~log N and **K*->1**. The fitted angle still converges,
but the returned law approaches the broad law
pi_theta(dtheta) exp[-sin^2(theta-theta0)/2], normalized. The additional
margin cannot be omitted. Section 4 shows that an order-invariant rule
with W must also miss some fixed histories, not just this triangular example.

## 7. Project consequence

The new candidate W cannot be freely added to unrestricted slow-signal
consistency and order-independent combination. This is a constraint on
what "no information" should mean asymptotically; it does not select pi.
Exact observed-zero invariance remains a different requirement.

| requirement combination | result |
| --- | --- |
| W plus concentration on every identifying triangular array | impossible, without any permutation or Z assumption |
| W plus concentration on every fixed identifying history plus P | impossible by the fixed-history construction |
| W plus current-median point consistency on those histories plus P | impossible for a fixed prior with the stated unique-median property |
| W plus fixed-history concentration without P/Z | compatible: chronological control above |
| W plus P/Z and consistency with the explicit additional growth margin | compatible on this scalar slice: continuous attenuation control; L19 and full premise set not established |
| Merely E0 Q_N->pi plus consistency | not decided by this washout proof; E10's earlier kernel obstruction has its own narrower scope |

No current premise is restated. O6, O7, O15 and O19 remain open. The live
question is whether the goal needs law-by-law null washout, an averaged
property, or a weaker uncertainty statement. The theorem excludes treating
those as interchangeable while promising arbitrary slow identification.

## Validation

All 17 new tests passed, including an independent reviewer's run. They
check Gaussian likelihood-ratio moments by quadrature; the permutation
identity by direct Gaussian moment-generating-function enumeration; the
product and Taylor bounds; attenuation identities and invariances; and the
ordered control's moment identity and deliberate invariance failures.
Finite checks witness the proof ingredients, not the infinite diagonal
quantifiers. No stochastic simulation study was needed.

On 9 Nov, the focused regression covering `test_null_washout.py`,
`test_null_limit.py`, `test_gaussian_information.py`, `test_energy_minimax.py`,
`test_signal_identifiability.py`, `test_starting_law_family.py` and
`test_audit.py` passed **117 tests**. The full production suite was not rerun
for these isolated research helpers. Existing audit findings remain unchanged.
