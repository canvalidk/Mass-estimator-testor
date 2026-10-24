# Data-based profile power: exact Bayesian boundary

Scratch analysis, 2026-09-28. This identifies the inference rule changed by
the candidate, rather than treating successful width calibration as a proof
that the candidate is or is not the desired equation.

## Proposed law

For one common noise pair, whiten and stack the observations z=(F,a), define
U_theta v=(sin(theta) v, cos(theta) v), and set

    A_theta = U_theta' z,
    H_theta = ||A_theta||^2,
    T = ||z||^2,
    Q_theta = T-H_theta.

Let h>=q be the eigenvalues of the 2 by 2 Gram matrix of F and a. The proposed
coefficient is beta(z)=(h-q)/h when h>0. The proposed law is

    p_beta(theta | z) proportional to
       pi(theta) exp[-beta(z) Q_theta/2].

The zero-data convention for beta is irrelevant to this angle law because
all Q_theta are then zero. A possible full joint candidate law retains

    v | theta,z ~ N(A_theta,I_M).

This construction is not ruled out merely because separate batch mass laws
do not multiply. Hierarchical Bayesian methods can legitimately require
additional carried state. The exclusions below use the full declared
Gaussian observation likelihood and are stronger than a batching complaint.

## 1. Direct obstruction for the specified full joint law

The declared whitened Gaussian likelihood factors as

    L(z | theta,v) proportional to
      exp[-Q_theta/2] exp[-||v-A_theta||^2/2].

Divide the proposed full joint density by L. A candidate measure that could
produce it by Bayes' rule would have to be, up to a data-only constant,

    mu_z(dtheta,dv) proportional to
      pi(theta) exp[(1-beta(z)) Q_theta(z)/2] dtheta dv.

In particular, keeping the exact Gaussian conditional v|theta,z forces
flatness in v. The remaining angle factor is data-dependent.

For a direct contradiction avoiding the exact-zero data point, choose
anisotropic nonzero data with 0<beta(z)<1 and compare z with c z, c!=1.
The coefficient beta is unchanged because both Gram eigenvalues scale by c^2,
whereas Q_theta scales by c^2. The ratio of the two required measures is

    exp[(1-beta(z))(c^2-1) Q_theta(z)/2],

which varies with theta. It cannot be absorbed into a data-only normalizing
constant. No fixed positive candidate measure gives this full joint law for
both data sets under the original Gaussian likelihood.

Adding hidden fixed-prior variables does not repair this: integrating them
out supplies a fixed effective measure on (theta,v), to which the same
argument applies. A latent variable that changes the observation likelihood
instead changes the declared measurement model and is outside this argument.

## 2. Obstruction for the mass marginal alone, under isotropy

Even if the conditional law of v is changed, any fixed isotropic candidate
measure at each angle gives an integrated likelihood of the form

    L_mu(theta;z) = exp(-T/2) F_theta(H_theta).

The function F_theta can depend on angle; no common measure or flatness in
mass is assumed. Posterior odds between two fixed angles therefore depend
on their two projected energies, not separately on the mutual angle of
the projected vectors.

The proposed law does depend on that mutual angle through beta(z). Here is
an explicit interior-angle example for M>=2. Choose theta_1=atan(1/2) and
theta_2=atan(2), so

    t1=A_theta1=(F+2a)/sqrt(5),
    t2=A_theta2=(2F+a)/sqrt(5).

Fix ||t1||^2=1 and ||t2||^2=4, and vary c=t1't2 between 0 and 1. Both choices
are feasible already in two spatial dimensions. The inverse transform is

    F=sqrt(5)(2t2-t1)/3,
    a=sqrt(5)(2t1-t2)/3.

The Gram trace and determinant of (F,a) are

    T=5(25-8c)/9,
    determinant=25(4-c^2)/9.

Hence

    beta at c=0 = 2 sqrt(481)/(25+sqrt(481)),
    beta at c=1 = 2 sqrt(181)/(17+sqrt(181)),

which differ. Any fixed isotropic candidate-measure construction gives
unchanged odds F_theta1(1)/F_theta2(4). The proposed law gives odds

    [pi(theta1)/pi(theta2)] exp[-3 beta(z)/2],

which change. Thus no fixed isotropic angle-dependent candidate measure
under the declared Gaussian likelihood produces even this mass marginal.
This statement does not exclude arbitrary anisotropic measures, changed
measurement likelihoods, or an explicitly non-Bayesian law construction.

## 3. Consequence for composition and coarsening

Bayesian posterior coarsening under the declared generative model obeys a
tower identity using that model's predictive law for the discarded readings.
The proposed full joint kernel cannot have that interpretation for any fixed
candidate measure with the original Gaussian likelihood. Therefore one
cannot claim that it inherits that entire Bayesian update/coarsening package.

This does not mean every possible operational coherence rule excludes the
candidate. Its data Gram matrix is additive, so storing the underlying sums
and recomputing the law can still be exactly order-independent and invariant
to batching. That computational fact is different from being a conditional
probability law under the declared fixed model. A new operational update or
uncertainty principle could be proposed, but it must be stated and justified.

At each fixed data set, the full candidate law is normalized and the
ratio-of-means readout still satisfies C1/C2 on that same propagated law.
The readout theorem alone therefore does not exclude the correction.

## 4. Precise source commitments and open physical premise

Managed reads: installed GitHub connector, pinned throughout to
canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.

`Newton-analysis/_dscn_mass_estimator/ours27/ours27_full_covariance_2026-09-27.md`
declares Gaussian errors independent across readings (lines 47-52), gives
the exp[-sum Q_i/2] law (lines 68-71), and gives the Gaussian conditional
candidates (lines 74-85). Its status table says flat push weighting and
determinant cancellation are forced (lines 129-130), readings multiply and
the prior is counted once (line 132), and pi is symmetric, declared before
data, with its scale from instruments (line 133).

`Newton-analysis/_dscn_mass_estimator/ours26/ours26_starting_law_2026-09-27.md`,
lines 94-95, says the starting law is part of the specification and fixed
before data; changing it changes the specification rather than updating it.

Importantly, the proposal can keep the explicit pi fixed. Thus the literal
before-data-pi rule ALONE does not refute it. The conflict is with obtaining
the stated joint law from the declared Gaussian experiment and fixed
candidate weighting. Reinterpreting it as Bayes with the same Gaussian
likelihood reveals a data-dependent effective candidate measure; describing
it as a direct posterior correction instead openly changes the inference
operation.

`Newton-analysis/_dscn_mass_estimator/ours27/ours27_candidate_size_and_null_readings_2026-09-28.md`,
lines 223-226, leaves universal planted-mass width calibration unsettled.
Consequently good repeated-experiment calibration does not itself authorize
replacing the current posterior semantics, and current semantics do not
automatically prove that every calibrated operational law is inadmissible.

The missing physical premise is a rule specifying what the reported law
means and how its probabilities must update: a conditional candidate law of
the stated Gaussian experiment, or another operational uncertainty object.
The choice beta=(h-q)/h obtains the desired regular local width, but local
asymptotic calibration and symmetry alone do not uniquely select its entire
finite-data form, nor uniquely justify keeping v|theta,z unchanged.
