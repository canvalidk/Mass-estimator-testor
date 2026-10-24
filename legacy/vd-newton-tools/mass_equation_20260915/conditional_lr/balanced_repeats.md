# Repeated readings of each latent vector

Bounded independent design audit, 2026-09-15. Scratch contribution only.

## Experiment and sufficient reduction

There are N groups, each with one fixed unknown acceleration vector v_i. In
group i, J independent readings observe the same pair (m v_i,v_i). The common
mass m is fixed and positive. Errors are independent Gaussian, with known
homogeneous per-reading variances sigma_F^2 and sigma_a^2. There is no drift,
within-group latent variation, or shared calibration error in this model.

Combine the J readings in each group **before fitting**:

\[
\overline F_i\sim N(mv_i,\sigma_F^2 I_d/J),\qquad
\overline a_i\sim N(v_i,\sigma_a^2 I_d/J).
\]

The group vector means are sufficient for the mean parameters. The within-
group residual scatter is ancillary when the covariance is known. Retaining
the raw readings instead of these sufficient means does not supply additional
mass likelihood information under this model.

The standardized group data are

\[
X_i=\sqrt J\,\overline F_i/\sigma_F,\qquad
Y_i=\sqrt J\,\overline a_i/\sigma_a.
\]

The instrument ratio s=sigma_F/sigma_a stays unchanged, while each group's
signal-to-noise ratios increase by sqrt(J). There are K=Nd nuisance vector
coordinates after reduction, not NJd.

Balanced replication also has a clean prior consequence. For fixed N and the
same parameter vector (m,v_1,...,v_N), its Fisher matrix is J times the
single-reading Fisher matrix. The joint Jeffreys volume therefore changes by
the constant `J^((Nd+1)/2)`. Recomputing that reference after balanced repeats
does not change its shape, and ordinary updating with the fixed reference
agrees. Adding new latent groups, or changing relative channel precisions, is
a different operation and need not preserve reference shape.

## Positive-signal coverage gap

Put

\[
e_N={1\over N}\sum_i\|v_i\|^2,\qquad
E_N=Ne_N,\qquad V_0=\sigma_F^2+m^2\sigma_a^2.
\]

Here e_N is the **average squared excitation magnitude**, with acceleration-
squared units. It is not the square of an average vector, which can cancel
when directions differ.

The squared standardized excitation in the reduced experiment is

\[
L_J={J E_N\over\sigma_a^2}.
\]

At the true r=m/s, the profile equation H=R dot S has

\[
\operatorname{Var}(H)=K+(1+r^2)L_J,\qquad E[H'(r)]=-L_J.
\]

Consequently the leading normal approximation for a locally consistent root
has mass variance

\[
\boxed{V_{\mathrm{sampling,local}}
={V_0\over JNe_N}
+{d\sigma_F^2\sigma_a^2\over J^2N e_N^2}.}
\]

The joint Jeffreys posterior's leading local variance is

\[
V_{\mathrm{posterior,local}}={V_0\over JNe_N}.
\]

Their ratio is therefore

\[
\boxed{\kappa_J
=1+{d\sigma_F^2\sigma_a^2\over J V_0 e_N}.}
\]

These are local asymptotic normal variances, not a claim that every finite-
sample root estimator has a finite exact second moment. The weak-identification
boundary and possible tail/branch pathologies remain separate questions.

For N tending to infinity with fixed J and positive limiting e_N=e, the
limiting coverage of its nominal 95% central posterior interval is

\[
2\Phi\left({z_{.975}\over\sqrt{\kappa_J}}\right)-1.
\]

In the existing example m=2, sigma_F=sigma_a=1, d=3, e=1,
`kappa_J=1+0.6/J`:

| Repeats J per fixed vector | Variance ratio | Limiting nominal-95% coverage as N grows |
|---:|---:|---:|
| 1 | 1.600000 | 87.8735% |
| 2 | 1.300000 | 91.4386% |
| 4 | 1.150000 | 93.2402% |
| 8 | 1.075000 | 94.1290% |
| 16 | 1.037500 | 94.5673% |
| 64 | 1.009375 | 94.8924% |

The table is an asymptotic calculation, not a finite-N simulation or a
universal recommendation for the number of repeats.

## Different limiting regimes

1. **N fixed, J grows, and E_N>0:** the nuisance dimension stays fixed; the
   group means become precise. The model is regular and identified, and the
   incidental-parameter coverage gap vanishes. Some groups may have v_i=0;
   at least one nonzero vector is sufficient for local mass identification.
2. **J fixed, N grows, with e_N tending to e>0:** the point is consistent,
   but the variance ratio stays above one. More noisy new nuisance vectors
   do not cure this posterior calibration gap.
3. **Both grow with positive average excitation:** within-group information
   grows, K/L_J tends to zero, and this particular local variance gap vanishes.
   The exact one-dimensional marginal still needs its stated regular regime;
   this is not a uniform weak-signal coverage theorem.
4. **Actual zero excitation in every group:** none of these replication
   arguments creates mass information. The next section is exact.

A useful distinction for weakening excitations is between direction learning
per nuisance coordinate and total identification. The ratio K/L_J depends
on `J e_N`; consistency of the profile root needs its signal L_J to dominate
the order-sqrt(K) random profile fluctuations. Thus a sequence can have a
consistent point while its posterior intervals remain badly miscalibrated.
The positive constant-e calculation should not be extrapolated uniformly to
all sequences approaching zero.

## Actual zero: repeating does not repair the pathology

If every v_i=0, then for **each J** the standardized sufficient means are
independent standard Gaussian vectors:

\[
X_i,Y_i\stackrel{ind}{\sim}N(0,I_d).
\]

Their distribution is exactly independent of J. The joint Jeffreys posterior
and its readout, as functions of this reduced standardized data, therefore
have exactly the same sampling law for every J at fixed N. Their coverage at
any fixed true mass is unchanged by within-group replication.

If N grows, the growing-N true-null behavior consequently persists for any
sequence J=J_N. Repetition shrinks the dimensional readings and their errors
together but supplies no signal relating force to acceleration. In the exact
zero-noise boundary, Newton II still permits every positive mass.

This is compatible with the separate observed-zero identity: when all
*observed* group means happen to equal zero, the posterior mass law is
half-Cauchy for every N,J. Random observations at a true null are noisy;
they need not produce that observed-zero posterior.

## Design implication at a fixed reading budget

If B=NJ readings are available and the average squared excitation e is held
comparable across designs, then

\[
I_{m\cdot V}={B e\over V_0}
\]

is independent of the allocation between new groups and repeated readings.
But the extra local sampling variance is

\[
{d\sigma_F^2\sigma_a^2\over J B e^2},
\]

so it decreases with J. Repeatedly measuring a genuinely stable vector can
improve nuisance estimation without sacrificing this local mass-information
budget. Treating every replicate as a new latent vector discards the known
within-group equality and leaves unnecessary nuisance noise.

Equivalently, let rho_F,J and rho_a,J be the force and acceleration SNRs of
each group's mean. Then

\[
\kappa_J
=1+{d\over \operatorname{average}(\rho_{F,J}^2+\rho_{a,J}^2)}.
\]

This concerns learning the nuisance direction. The absolute precision of mass
also has the distinct weaker-channel bottleneck
`I_log(m)^(-1)=rho_F^(-2)+rho_a^(-2)`. A well-learned direction does not by
itself make the mass precise when the other channel remains weak.

The design comparison is conditional on fixed vectors within groups, stable
mass, independent known Gaussian errors, homogeneous scales, and comparable
excitation. Correlated errors or a calibration floor prevent the simple 1/J
variance reduction. If the acceleration actually changes within a group,
averaging it as one unchanged latent vector is a model change, not a free
information gain. This calculation also does not replace any need for diverse
conditions to investigate systematic effects outside the stated model.

## Fair comparisons for the investigation

- Give compared algorithms the same raw Gaussian experiment or its same
  sufficient group means and correctly reduced covariance.
- Specify the same true-state sharing structure: one vector per reading,
  one vector per group, and one common vector are different models.
- Separate a fixed-observation-budget design comparison from an algorithm
  comparison. Changing excitation strength or the number of readings changes
  the available information.
- Compare conditional posterior probability with conditional posterior
  probability, or compare repeated-sampling uncertainty at a common coverage
  requirement. Nominal 95% labels alone do not equate those guarantees.
- For anisotropic covariance, replace squared excitation norms by the
  corresponding covariance-weighted information; changing direction can then
  change the mass-information budget.

No production code, policy, or managed source was edited in this audit.
