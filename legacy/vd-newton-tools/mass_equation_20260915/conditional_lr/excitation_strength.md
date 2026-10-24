# Excitation strength, repeated measurements, and the null LR limit

Scratch derivation for the mass-equation investigation, 2026-09-15. This is
not a managed record or a replacement estimator policy. No production file
is changed. The calculations below are derived from the explicit models.

## 1. Three experiments must remain distinct

1. **Repeated observations of one latent vector v.** The latent parameter is
   `(m,v)` throughout. Independent equal-covariance observations can be
   aggregated into their full vector mean with covariance `Sigma/N`; the
   remaining within-sample scatter has no likelihood information about the
   mean parameters. Nuisance dimension stays fixed.
2. **One new fixed unknown v_i per trial.** The latent parameter is
   `(m,v_1,...,v_N)`. Nuisance dimension increases with N. A joint reference
   prior can give a proper posterior without guaranteeing its fixed-design
   sampling calibration.
3. **Independent random excitations.** Declaring
   `v_i iid N(0,tau^2 I_d)` supplies an excitation model. It is not a
   consequence of Gaussian measurement error. If the same v is reused,
   integrating that shared latent variable produces cross-reading
   dependence, not the independent marginal likelihood below.

The random model describes acceleration excitation. A force-controlled
protocol may call for a different model. Reciprocal reparameterization must
transform the same joint distribution rather than reassign a fixed
excitation variance to the swapped channel.

## 2. Exact random-excitation marginal

Let

    Y_i | m,v_i ~ N(A_m v_i, Sigma),  A_m=[m I_d; I_d],
    v_i iid N(0,tau^2 I_d).

Then independently over trials,

    Y_i | m,tau ~ N(0, Sigma + tau^2 A_m A_m^T).

For independent isotropic channel errors, put

    A=sigma_F^2, B=sigma_a^2, k=tau^2,
    V=A+m^2 B, D=AB+kV, n=Nd.

The per-coordinate covariance is

    C(m,k) = [[A+m^2 k, mk], [mk, B+k]],   det C=D.

With S_F=sum_i ||F_i||^2, S_a=sum_i ||a_i||^2, and
S_Fa=sum_i F_i dot a_i, its log likelihood is

    ell(m,k) = -n log(D)/2
      - [(k+B) S_F + (m^2 k+A) S_a - 2mk S_Fa]/(2D) + constant.

These are second moments about the explicitly specified zero mean; they
are not sample covariances formed by silently subtracting fitted means.
Known noise and k>0 identify m and k. At k=0 the marginal law does not
depend on m.

## 3. Whitened eigenvalue form with unknown strength

Write s=sigma_F/sigma_a, r=m/s, and standardize observed coordinate pairs
as X=F/sigma_F and Y=a/sigma_a. Let S be their empirical 2-by-2 second-moment
matrix divided by n. Define

    u_r=(r,1)/sqrt(1+r^2),
    gamma=(k/B)(1+r^2),
    s_r=u_r^T S u_r.

The covariance is I_2+gamma u_r u_r^T. Its log likelihood is

    ell(r,gamma) = -n/2 [log(1+gamma) + tr(S)
                           - gamma s_r/(1+gamma)] + constant.

For each candidate r, profiling gamma>=0 gives

    gamma_hat(r)=max(s_r-1,0),

and the log-likelihood gain above the null covariance I_2 is

    G(r) = n/2 [s_r-1-log(s_r)] if s_r>1, and 0 otherwise.          (1)

This is monotone in s_r above 1. In an informative interior case it selects
the same direction as fixed-vector profile/TLS, which maximizes s_r.
If every allowed s_r<=1, the random-excitation fit selects gamma=0 and no
mass is identified. Fixed-vector profile can still select a noise direction.

### Positive mass restricts the eigenvalue optimization

For r>0 the two components of u_r have the same sign. Therefore

    sup_{r>0} s_r = lambda_max(S),              if S_12>=0,
                   max(S_11,S_22),             if S_12<0.           (2)

In the second case the supremum is at r->0 or r->infinity. Using the
unrestricted largest eigenvalue when S_12<0 instead fits an opposite-sign
force/acceleration direction, outside the positive-mass model. If signed
mass is deliberately included as an enlarged alternative, say so.

## 4. Learned strength supplies the extra sampling-variance term

At an interior maximum with s_r>1 and s_r'=0, differentiating (1) gives

    -ell_profile_strength''
       = [(s_r-1)/s_r] [-ell_profile_fixed_vectors''].               (3)

The distinction is curvature, even where the mass point agrees.

Direct covariance Fisher calculation, per coordinate, gives

    I_mm = k^2 (D+2m^2 B^2)/D^2,
    I_mk = kmBV/D^2,
    I_kk = V^2/(2D^2).

Eliminating unknown k gives

    I_{m.k}=k^2/D,
    Var(m_hat) approximately D/(Nd k^2).                           (4)

For fixed nuisance vectors let S_v=sum_i ||v_i||^2. The fixed-vector
profile score admits the decomposition

    U_i=(v_i+zeta_i) dot r_i / V,
    r_i=epsilon_F,i-m epsilon_a,i,
    zeta_i=epsilon_a,i+(mB/V)r_i.

The Gaussian vectors r_i and zeta_i are independent with component
variances V and AB/V. Hence expected curvature and score variance are

    H=S_v/V,
    J=S_v/V+Nd AB/V^2,

and the fixed-design sandwich variance is

    J/H^2=V/S_v+Nd AB/S_v^2.                                       (5)

When S_v approximately equals Nd k, (4) and (5) coincide. This needs positive
limiting average excitation and regularity/no-dominant-trial conditions.
It does not make the Gaussian excitation distribution a deduced physical
fact. The covariance likelihood also deserves examination as a working
second-moment procedure beyond a literal Gaussian excitation model.

For d=3, m=2, A=B=1, and ||v_i||=1, k=S_v/(Nd)=1/3 and V=5. The
curvature multiplier in (3) tends to 5/8. Inverse fixed-vector curvature
is 5/N; learned-strength curvature gives 8/N, matching (5).

By contrast, if k is declared known incorrectly when actual excitation
energy is k_0, the expected mass score at the true m is

    (k_0-k) kmBV/D^2

per coordinate. This is generally nonzero. An incorrectly fixed excitation
strength can therefore bias mass even with unlimited observations.

## 5. Null behavior and the reason ordinary Wilks calibration fails

At true zero excitation, the whitened observations are independent N(0,I_2)
coordinate pairs, and

    sqrt(n)(S-I_2) => G,

where G is a 2-by-2 GOE matrix: its independent diagonal entries have
variance 2 and its off-diagonal entry has variance 1. Uniformly over r,

    s_r=1+u_r^T G u_r/sqrt(n)+o_p(n^-1/2).

The fixed-vector profile has varying log-likelihood terms of order sqrt(n).
It can acquire an increasingly sharp noise-selected direction. Equation
(1) instead gives

    G_profile(r) => (u_r^T G u_r)_+^2/4.                            (6)

Its relative log-likelihood evidence stays O_p(1); mass remains unidentified.
A regular proper prior with suitable support near k=0 may avoid the previous
explosive concentration, but (6) is not a proof of Bayesian coverage or of
returning exactly to the original mass prior. The null is a boundary with an
unidentified mass; regular BvM/Wilks claims cannot be imported there.

For a fixed candidate direction u and a chosen set of alternative directions,
the null limit of twice the profiled log-likelihood ratio is

    T_limit = 1/2 [(g_max)_+^2 - (u^T G u)_+^2],                    (7)

where g_max is the supremum of the quadratic form over those alternatives.

* For unrestricted orientation, g_max=lambda_max(G).
* For positive mass, g_max is lambda_max(G) if G_12>=0 and is
  max(G_11,G_22) otherwise, by (2).

### An atom at zero disproves a chi-square-1 limiting distribution

For unrestricted orientation, write

    t=(G_11+G_22)/2, d=(G_11-G_22)/2, b=G_12.

These are iid standard normal, and lambda_max=t+sqrt(d^2+b^2).
Thus P(lambda_max<=0) is the solid-angle probability of a cone with half
angle pi/4:

    P(T_limit=0)=(1-1/sqrt(2))/2 = 0.1464466094.                    (8)

For an interior positive-mass candidate, with a positive-mass alternative,
the corresponding zero atom is

    P(G_12<0,G_11<=0,G_22<=0)
      + P(G_12>=0,lambda_max<=0)
      = 1/8 + (1-1/sqrt(2))/4 = 0.1982233047.                      (9)

Other exact equality with the maximum has probability zero for a fixed
interior candidate. Chi-square-1 has no atom at zero, so neither null LR
limit is chi-square-1. This alone does not prove a particular cutoff cannot
coincidentally have 95% coverage; the numerical check below examines it.

## 6. Bounded GOE simulation of the nominal chi-square cutoff

Independent check on 2026-09-15: 10,000,000 GOE draws, NumPy default_rng seed
763915, in 50 batches of 200,000; threshold 3.841458820694124. The diagonal
draws were sqrt(2) times standard normal and the off-diagonal draw standard
normal. The same draws were used for all candidate ratios.

| Alternative | Candidate r | Limit-statistic rejection fraction |
|---|---:|---:|
| Positive mass | 0.01 | 0.0368286 |
| Positive mass | 0.1 | 0.0342638 |
| Positive mass | 1 | 0.0266750 |
| Positive mass | 2 | 0.0272425 |
| Positive mass | 16 | 0.0354337 |
| Unrestricted orientation | 1 | 0.0534699 |

The positive-mass zero-atom frequency was 0.1981430, consistent with (9).
At a 5% rejection probability the Monte Carlo standard error is about
0.000069 (0.0069 percentage points). These are simulation estimates of the
derived limiting law, not exact quantiles or a finite-N calibration proof.
The unrestricted distribution does not depend on candidate direction by GOE
rotation invariance. The positive-mass distribution generally does depend
on the candidate's position inside the allowed angular sector.

Thus near-95% finite-N coverage with an unrestricted eigenvalue statistic is
compatible with its approximately 94.65% limiting coverage at that cutoff.
It must not be described as exact or asymptotically chi-square calibration.
The main investigation's separately derived conditional LR construction can
retain its exact guarantee without relying on this approximation.

## Scope and next use

The derivation explains why estimating a shared excitation strength changes
the nuisance treatment and uncertainty curvature. It does not recommend an
excitation model without experimental warrant. A next comparison should say
which quantities are physically repeated, whether an excitation distribution
is supplied or merely a working model, and whether the uncertainty claim is
posterior probability, exact fixed-design coverage, or an asymptotic result.
