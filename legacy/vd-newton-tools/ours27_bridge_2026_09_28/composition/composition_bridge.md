# Composition of candidate laws: what is selected, and what is not

Scratch analysis, 2026-09-28. No managed or production writes.

## Source boundary

Managed reads used the installed GitHub connector at
canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.
The catalog lists the following as the canonical home sources; the uniqueness
proof also informs the equation bundle through a copy placed September 26.

1. `Newton-analysis/_dscn_mass_estimator/mass_readout_composition_uniqueness_2026-09-06.md`:
   C1/C2 are explicitly on one propagated joint law. Line 37 says they do not
   equate that law with a fresh fit using new evidence or a newly assigned
   measure. Lines 17 and 140 say the theorem does not construct the input law.
   Line 142 separates repeated observations of one candidate pair, different
   trials sharing one mass, and regrouping an existing law.
2. `Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md`:
   section 4, lines 98-109, identifies the composite support as the pushforward
   of the parts' joint support. Line 136 grounds composition in estimator-is-
   mass; line 138 expressly excludes equality to a fresh total-force fit under
   a fresh flat measure. Its older many-reading proposal in section 6 is
   historical and is not adopted here.
3. `Newton-analysis/_dscn_mass_estimator/what_forces_the_chain_to_the_estimator.md`:
   S10 treats composition as a premise whose readout consequence is derived.
4. `Newton-analysis/_dscn_mass_estimator/mass_estimator_additional_derivations_2026-09-05.md`:
   section 3 proves acceleration-weighted averaging when one law is
   partitioned into cases, and excludes arbitrary pooling of separately
   fitted experiments.
5. `Newton-analysis/_dscn_mass_estimator/ours27/ours27_candidate_size_and_null_readings_2026-09-28.md`:
   section 2 lists law-level lumping as an untested hinge; section 3 gives
   the one-noise-pair measured-null point sigma_F/sigma_a. No-readings 0/0 is
   separately undefined and is not used in the example below.

## 1. A constructive joint experiment

For two bodies certified to have the same acceleration candidate a, use

    F1* = m1 a,   F2* = m2 a,   m1,m2 > 0,
    M = m1+m2,    Ftotal* = M a.

Let the observed vector be y=(F1,F2,A) with known Gaussian covariance Sigma.
Choose a positive joint candidate measure Pi(dm1,dm2,da) with finite required
Gaussian integrals and first moments. This remains the unresolved measure
choice; writing the joint model does not make it unique.

The joint candidate law is

    P_y(d m1,d m2,d a) proportional to
      phi_Sigma(y - (m1 a,m2 a,a)) Pi(d m1,d m2,d a).

On this same law, let alpha=||a||. The readouts obey exactly

    r1 = E[m1 alpha]/E[alpha],
    r2 = E[m2 alpha]/E[alpha],
    rtotal = E[(m1+m2) alpha]/E[alpha] = r1+r2.

No independent-fitting or independent-mass assumption is used. The observed
acceleration likelihood is included once. Computing two separate posteriors
and multiplying them generally repeats this evidence and can impose two
different acceleration laws instead of one shared acceleration.

The common-force construction is the exact dual: use one force candidate b,
inverse masses mu_i=1/m_i, accelerations a_i*=mu_i b, and relative acceleration
a_rel*=(mu_1+mu_2)b. On that propagated law, inverse readouts add. Summing the
acceleration observations can again discard a contrast.

For persistent objects in several trials, use one addressed m1 and m2 across
trials, trial-specific common accelerations a_j when justified, and M=m1+m2
throughout. This formalizes a persistence/composition premise; it does not
derive persistence from a noise model or select Pi.

## 2. Coarsening is not an invertible change of representation

Define L y=(S,A), S=F1+F2. Then the correct coarse measurement covariance is

    Sigma_coarse = L Sigma L',

and the coarse candidate map is q(m1,m2,a)=(M,a). A coherent starting measure
must be Pi_coarse=q#Pi, not a separately reset default measure.

Even with this prior compatibility, fitting only (S,A) usually differs from
pushing forward the posterior fitted to all (F1,F2,A). The latter has more
evidence. It retains a contrast that summation erases.

For transparent algebra, take independent force errors v1 I_d and v2 I_d,
independent acceleration error, and put

    vT=v1+v2,  w=v1/vT,
    C=F1-w S,  vC=v1 v2/vT,
    delta=m1-w M.

S and C have independent noise. Their candidate means are M a and delta a,
and the positive-mass condition is

    -w M < delta < (1-w) M.

Disintegrate the measure as Pi_coarse(dM,da) K(d delta | M,a). The total
candidate posterior from the fine data is

    P_fine(dM,da | S,C,A) proportional to
      L_coarse(S,A | M,a) Pi_coarse(dM,da) Kbar(C | M,a),

    Kbar(C | M,a) = integral phi_(vC I)(C-delta a)
                             K(d delta | M,a).

The coarse posterior omits Kbar. Keeping (S,C,A) is an invertible
representation change and preserves the entire inference. Keeping only
(S,A) discards evidence. With general correlated Sigma, the corresponding
discarded factor is the conditional Gaussian likelihood obtained by the
usual block covariance / Schur complement decomposition.

## 3. Exact fine/coarse posterior equality forces a special split

Suppose equality of the full coarse candidate laws is demanded for every
data vector and every (M,a), with the same positive coarse measure. This
requires Kbar(C | M,a) to be parameter-independent. Gaussian convolution is
injective, so the distribution of delta a itself must be independent of M,a.

Fix any nonzero a and let M decrease to zero. The allowed interval for delta
shrinks to zero, so the common distribution can only be a point mass at zero.
Therefore

    delta=0,  m1=w M,  m2=(1-w) M

almost surely. Conversely this split makes C pure noise and ancillary.

Thus a universal exact fine/coarse equality requirement, in this model, fixes
the mass split by the force instruments' error variances. That is extra
object-specific structure, not a consequence of Newton II. The statement
assumes proper conditional split kernels and a domain extending to all M>0;
restricted candidate domains can evade the shrinking-support step.

This proves an obstruction for full (M,a)-law equality. Equality of only one
point or one mass marginal is weaker and is not excluded by this argument
alone; accidental cancellations are possible. The example next shows that
even mass-point equality fails in valid joint models.

## 4. Exact counterexample with the same summed observation

Take a=1 known exactly and force error variances v1=v2=1. Give equal prior
weight to two candidate states:

    state 1: (m1,m2)=(1,1), M=2, expected (S,C)=(2,0);
    state 2: (m1,m2)=(3,1), M=4, expected (S,C)=(4,1).

At observed S=3 the coarse likelihoods are equal, so its mass readout is 3.
But the contrast has variance 1/2, and the fine posterior odds of state 2
against state 1 are exp(2C-1). Hence

    r_fine = 2 + 2/(1+exp(1-2C)).

At C=0 this is about 2.53788; at C=1 it is about 3.46212. Both data sets have
the same (S,A,Sigma_coarse). Each fine fit satisfies rtotal=r1+r2 exactly.
No function of only the summed pair and its covariance can reproduce both.
This example uses a declared finite candidate measure to expose the lost
information; it is not a recommendation to use that measure as the default.

## 5. The correct coarsening equation

When Pi_coarse is the pushforward of Pi, Gaussian marginalization gives

    P_coarse(dM,da | S,A)
      = integral P_fine(dM,da | S,C,A) p(dC | S,A).

The density p(dC | S,A) is the predictive law supplied by this same model.
This is the posterior tower identity. It makes discarding information
coherent without pretending that the discarded contrast was never useful.

Let A_C=E[||a|| | S,C,A] and r_C be the fine total readout. Then

    r_coarse = integral A_C r_C p(dC | S,A)
               / integral A_C p(dC | S,A).

This is exactly the already-proved acceleration-weighted case-aggregation
law. The natural bridge is therefore: transport one candidate law through
physical composition and average consistently when evidence is actually
forgotten. Ordinary conditional-probability inference satisfies the bridge
for every admissible Pi. The bridge does not select a primitive Pi.

The carried state cannot generally be only the marginal mass law: the readout
also needs the conditional acceleration-size weights. A persistent shared
mass address fixes candidate identity across trials; it does not make those
weights recoverable from p(m) alone.

## 6. What compatibility does constrain about the measure

In common direction and common acceleration coordinates, suppose the fine
measure happens to be df1 df2 times a measure on acceleration and direction.
Under f=f1+f2, integrating the split gives

    (sum)#(df1 df2) = f df,
    (sum)#(product_i df_i) = f^(k-1)/(k-1)! df.

Thus resetting a flat df measure on the composite does not preserve the
starting measure. This is an elementary concrete law-level constraint. The
flat-force example is illustrative; it is not asserted to equal ours27's
current mass-angle/radius measure.

There are many coherent alternatives. Independent gamma force measures with
one common rate and positive shapes a_i push forward to a gamma law with
shape sum_i a_i; the conditional split is Dirichlet. Rates and shapes remain
free. A projectively coherent object-partition system would require an
additive base assignment of shapes to physical parts, which is additional
structure rather than a quantity supplied by instrument covariance alone.

## 7. A diagnostic for an over-strong bridge applied to ours27

Take two measured 0/0 part pairs, with independent force errors of equal SD
sigma_F and one shared acceleration observation of SD sigma_a. A separate
ours27 fit gives each part the point s=sigma_F/sigma_a. If a proposed bridge
preserves both separately fitted magnitude laws as marginals on a common
acceleration law, C1 forces the composite point to be 2s.

A fresh fit to the summed pair has force SD sqrt(2) sigma_F and therefore
gives sqrt(2) s at measured null. These differ. Independence of the candidate
masses is not needed: C1 is valid for every common-acceleration coupling.

This refutes the simultaneous demands of separate-fit marginal fidelity,
same-law composition, and equality to a fresh summed-pair default fit. It
does not refute C1, and it does not by itself refute ours27. A genuine joint
fit can change the part marginals; a coarse fit has discarded the contrast;
and a propagated composite measure need not equal the fresh default.

The unresolved physical choice is which of those distinct operations the
user's definition intends to identify. Existing sources deliberately keep
them distinct. Simply asserting that mass composes cannot silently add all
three requirements.
