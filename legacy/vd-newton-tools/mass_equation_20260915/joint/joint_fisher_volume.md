# One mass across several latent trials: Fisher volume and its marginal

Scratch derivation, 2026-09-15. No production or managed edits. The inputs are
the explicitly stated sampling model below; no new managed snapshot was read.

## 1. Model and Fisher determinant

For i=1,...,N, independently observe

    F_i ~ N(m v_i, sf_i^2 I_d),
    A_i ~ N(v_i, sa_i^2 I_d),

with independent channels, known positive standard deviations, one m>0, and
nuisance vectors v_i in R^d. Here F_i,A_i denote the observed vectors. Put

    a_i=sf_i^(-2), b_i=sa_i^(-2),
    c_i(m)=a_i m^2+b_i,
    D_i(m)=sf_i^2+m^2 sa_i^2,
    S(m,V)=sum_i ||v_i||^2/D_i(m).

The parameter dimension is 1+Nd. The Fisher matrix has blocks

    I_mm = sum_i a_i ||v_i||^2,
    I_m,vi = a_i m v_i^T,
    I_vi,vi = c_i I_d,
    I_vi,vj = 0 for i != j.

Taking the Schur complement of the nuisance block gives

    I_mm - sum_i a_i^2 m^2 ||v_i||^2/c_i
       = sum_i (a_i b_i/c_i) ||v_i||^2 = S.

Consequently the full joint Jeffreys volume is

    dPi_J proportional to
       sqrt(S(m,V)) product_i c_i(m)^(d/2) dm product_i dv_i.             (1)

The proposed determinant formula is exact. N=1 reduces to
||v||/(sf sa) * (m^2/sf^2+1/sa^2)^((d-1)/2), as in the single-trial result.

The mass information after accounting for the nuisance vectors is S: the
individual contributions add. The remaining product is the nuisance-volume
factor. Independence of the likelihood does not require the joint reference
prior to factor; the square root of a sum in (1) couples the nuisance norms.

If some v_i vanish and at least one does not, S>0 and the full Fisher matrix
still has full rank. At V=0, S=0 and mass is unidentifiable. The image of all
(m,0,...,0) is one point, the singular vertex; (1) places no atom there. A
model with a literal prior atom on that physical null needs separate treatment.

## 2. Integrating all nuisance vectors exactly

Complete each likelihood square using

    mu_i(m)=(a_i m F_i+b_i A_i)/c_i,
    eta_i(m)=sqrt(c_i) mu_i,
    Q_i(m)=||F_i-m A_i||^2/D_i,
    w_i(m)=sf_i sa_i/D_i.

Then

    a_i||F_i-mv_i||^2+b_i||A_i-v_i||^2
       = c_i||v_i-mu_i||^2+Q_i.

Set z_i=sqrt(c_i)(v_i-mu_i). The nuisance Jacobians cancel the product of
c_i^(d/2) in (1), while

    S = sum_i w_i^2 ||z_i+eta_i||^2.

Thus, up to a positive constant independent of m, the exact marginal is

    q_J(m) = exp[-sum_i Q_i(m)/2]
             E sqrt(sum_i w_i(m)^2 ||Z_i+eta_i(m)||^2),                 (2)

where the Z_i are independent standard d-dimensional Gaussian vectors.
The expectation in (2) is a residual finite-dimensional Gaussian norm
expectation, not an approximation or a fitted mass prior.

### Properness for every finite dataset

The eta_i are bounded, since

    ||eta_i(m)|| <= sqrt(||F_i||^2/sf_i^2+||A_i||^2/sa_i^2).

At m down to zero,

    eta_i -> A_i/sa_i,
    w_i -> sa_i/sf_i,
    Q_i -> ||F_i||^2/sf_i^2.

At m to infinity,

    eta_i -> F_i/sf_i,
    m^2 w_i -> sf_i/sa_i,
    Q_i -> ||A_i||^2/sa_i^2.

Gaussian integrability bounds the norm expectations uniformly for these
limits. Each limiting expectation is finite and strictly positive. Therefore

    q_J(m) -> C_0 in (0,infinity) as m down to zero,
    q_J(m) ~ C_infinity/m^2 with C_infinity>0 as m to infinity.          (3)

On every compact interval of positive m the function is continuous and finite.
Equation (3) proves a proper posterior for every finite dataset and every
finite N>=1. Exactly the power moments -1<p<1 are finite. In particular,
ordinary mass and inverse-mass means, and ordinary mass variance, diverge
even with many trials. All fixed-order absolute log-mass moments are finite.
These tail statements do not determine the widths of central intervals.

The reference prior itself is improper. Posterior properness is a consequence
of the Gaussian likelihood and the full joint volume, not a normalization of
the prior over the unbounded parameter space.

## 3. Common instrument ratio: a compact marginal

Suppose sf_i/sa_i=s for every i; their absolute scales may differ. Put r=m/s
and stack standardized observations X_i=F_i/sf_i, Y_i=A_i/sa_i into two
K=Nd component vectors X,Y. Define

    R_r=(X-rY)/sqrt(1+r^2),
    S_r=(rX+Y)/sqrt(1+r^2),
    h_K(t)=E||Z+t e_1||,  Z~N(0,I_K).

In (2), sum Q_i=||R_r||^2, the stacked eta vector equals S_r, and every
w_i=1/[s(1+r^2)]. Multiplying by dm=s dr gives

    p_J(r | data) proportional to
       exp[-||R_r||^2/2] h_K(||S_r||)/(1+r^2).                        (4)

The function h_K is the mean of a noncentral chi distribution. Equation (4)
is exact. In theta=atan(r) coordinates,

    p_J(theta | data) proportional to
       exp[-||R_theta||^2/2] h_K(||S_theta||), 0<theta<pi/2.

This stacking does not require the physical directions of different trials
to agree. It works because every standardized latent force component equals
r times its corresponding standardized acceleration component.

A bounded one-dimensional integral for the noncentral chi mean is

    h_K(t)=sqrt(2/pi) integral_0^(pi/2)
       [1-cos(beta)^K exp(-t^2 sin(beta)^2/2)]/sin(beta)^2 d beta.

It follows from sqrt(q)=(2 sqrt(pi))^(-1) integral_0^infinity
(1-exp(-u q))u^(-3/2)du and u=tan(beta)^2/2. The apparent zero-endpoint
singularity is removable, with integrand limit (K+t^2)/2. Numerically the
numerator can be evaluated as -expm1(K log(cos beta)-t^2 sin(beta)^2/2).

### The corresponding stacked-norm ratio-of-means point

For homogeneous channel standard deviations across trials, let beta be the
standardized latent acceleration stack. After square completion, use
beta=cos(theta)(S_theta+Z). The full conditional Jeffreys posterior of Z is
the standard Gaussian tilted by ||S_theta+Z||; it is not the unweighted
Gaussian square-completion law. Therefore the extra norm in the readout
turns the nuisance expectation into E||S_theta+Z||^2=K+||S_theta||^2:

    m_hat_stack = s *
       integral e^(-||R_theta||^2/2) sin(theta)(K+||S_theta||^2) dtheta
       / integral e^(-||R_theta||^2/2) cos(theta)(K+||S_theta||^2) dtheta.

This is E[m||V||]/E||V|| for the concatenated physical acceleration vector V.
Choosing that aggregate operation is an additional application decision;
physical addition of vectors across different trials is not assumed.

If only the instrument ratios are common while sa_i vary, this simple formula
still applies to a standardized-stack norm. A physical-stack norm instead
requires the expectation

    E[||S+Z|| sqrt(sum_i sa_i^2 ||S_i+Z_i||^2)]

in place of K+||S||^2 in both angular integrals. Homogeneous sa_i reduce that
expectation to a constant multiple of K+||S||^2, and the constant cancels.

### All observed vectors zero

Then R_r=S_r=0 for every r, and h_K(0) is constant. Thus for every N and d,

    p_J(m | all observations zero)=2s/[pi(s^2+m^2)],
    p_J(log(m/s)=ell)=1/[pi cosh(ell)].                               (5)

Adding further observed-zero trials with the same instrument ratio does not
narrow this conditional mass law. The ratio-of-means readout formed from any
one trial is s here, as is a pooled first-moment readout.

For unequal ratios s_i=sf_i/sa_i, the null marginal is instead

    q_J(m) proportional to
       E sqrt(sum_i [s_i/(s_i^2+m^2)]^2 Zeta_i^2),

with independent Zeta_i~chi_d. Its shape depends on the declared instruments,
but (3) still applies.

Observed zeros and true zero excitation remain distinct. At actual V=0 the
sampling distribution is identical for all m; the experiment cannot identify
mass there. A finite conditional distribution does not change that fact.

## 4. Why multiplying single-trial priors is a different construction

In common (m,V) coordinates, naively multiplying the N single-trial Jeffreys
density functions and retaining only one dm gives

    pi_product(m,V) proportional to
       product_i [sqrt(||v_i||^2/D_i) c_i^(d/2)].

Its ratio to the actual joint density is

    pi_J/pi_product = sqrt(sum_i s_i*) / product_i sqrt(s_i*),
    s_i*=||v_i||^2/D_i,

where this comparison is away from individual zero vectors. At a vanishing
single trial, the product density vanishes even if other trials identify m;
the joint density need not. Their different behavior reflects different
reference assignments, not a determinant approximation.

There is also a coordinate warning. Each single-trial Jeffreys measure is
invariant, but multiplying their densities and identifying their mass
coordinates is not automatically an invariant operation. Repeating that
prescription in inverse-mass/force-vector coordinates and transforming back
multiplies it by m^(2N-2). The true joint determinant transforms correctly.

Indeed under mu=1/m and W_i=m v_i, the full coordinate Jacobian is
|dmu product dW_i/(dm product dv_i)|=m^(Nd-2). The true joint density in the
new coordinates picks up m^(2-Nd), which cancels. Multiplying N separate
single-trial density transformations instead counts a mass Jacobian N times
while the common-coordinate volume contains it only once.

## 5. Products of null mass posteriors give two explicit different answers

Each single-trial null mass posterior with the same instrument ratio is
the half-Cauchy (5). Multiplying posterior densities is not generally the
same as multiplying likelihoods: a common prior and nuisance model must be
specified before such a combination is licensed.

### Product in mass coordinates

Multiplying the densities relative to dm gives

    p_product(m) proportional to (s^2+m^2)^(-N).

For r=m/s its normalized density is

    [2 Gamma(N)/(sqrt(pi) Gamma(N-1/2))] (1+r^2)^(-N).

For N>1 it concentrates toward the boundary m=0 as N grows. Its scale is
O(s/sqrt(N)); it is not reciprocal-symmetric about s.

### Product in log-mass coordinates

Multiplying the densities of ell=log(m/s) instead gives

    p_product(ell) =
       Gamma((N+1)/2)/(sqrt(pi) Gamma(N/2)) sech(ell)^N.

Its exact moment-generating function and log variance are

    E exp(t ell)=Gamma((N+t)/2)Gamma((N-t)/2)/Gamma(N/2)^2, |t|<N,
    Var(ell)=(1/2) trigamma(N/2) ~ 1/N.

It stays centered at m=s while becoming arbitrarily narrow. The corresponding
mass density is proportional to m^(N-1)/(s^2+m^2)^N. It is distinct from both
the dm-product and the actual joint Jeffreys posterior (5).

Thus an unchanged set of null likelihoods can acquire shrinking apparent
uncertainty merely through a pooling/reference convention. That narrowing is
not new information about m at the actual zero-excitation state. A hierarchical
model for excitation could deliberately supply extra assumptions; it must be
stated rather than inferred from a product of already-normalized posteriors.

## 6. Independent check of the proposed score confidence set

For the homogeneous standardized model above, at the true r the residual
R_r is N(0,I_K) and independent of S_r. Therefore

    Z(r)=R_r dot S_r / ||S_r||

is exactly standard normal conditional on S_r. The event S_r=0 has probability
zero for fixed r. Inverting |Z(r)|<=z_(.975) over positive r gives an exact
95% confidence set under that model. It need not be a bounded connected
interval, and it is not an induced posterior probability statement.

Put A=||X||^2, B=||Y||^2, c=X dot Y. The score numerator is

    c(1-r^2)+(A-B)r.

For c<0 its positive zero is the profile maximum rather than its minimum.
This does not invalidate the pivot, but it prevents interpreting an isolated
stationary root as the best-fit mass. For K>=2 the event that S_r vanishes
at any positive r requires exact antiparallel X,Y, a probability-zero data
set under the nonsingular model. A robust implementation still needs an
explicit rule for such degenerate supplied inputs.

## Verification status

The determinant, marginalization and tail results are analytical derivations.
Numerical checks in check_joint.py verify block determinants and reciprocal
coordinate transformation for several N,d; independent quadrature checks the
null pooling normalizations and log variances. These checks supplement the
proofs rather than establishing general sampling calibration.

## 7. True zero excitation with growing nuisance dimension

This is different from the exactly zero *observed* dataset in (5). Let every
true v_i be zero, let K=Nd tend to infinity, and fix any true r_0 in (0,infinity).
The standardized observations are independent X,Y~N(0,I_K); their distribution
does not depend on r_0. The posterior for every finite K remains the proper
joint Jeffreys posterior (4), but its central 95% credible intervals have
coverage tending to zero.

### The random angular likelihood

Put

    A=||X||^2, B=||Y||^2, C=X dot Y,
    U=(A+B)/2, G=(A-B)/2, rho=sqrt(G^2+C^2).

At theta=atan(r),

    T(theta)=||R_theta||^2=U+G cos(2theta)-C sin(2theta),
    ||S_theta||^2=U-G cos(2theta)+C sin(2theta).

The ordinary Gaussian central limit theorem gives

    (G,C)/sqrt(K) => N(0,I_2),
    U/K -> 1 in probability.

In particular rho is of order sqrt(K). The random angular variation in the
profile log likelihood is also of that order, even though the true experiment
contains no mass information at V=0.

### Why the norm factor does not prevent concentration

Uniformly over theta, ||S_theta||^2 lies between U-rho and U+rho. With
probability tending to one these bounds are positive constants times K.
The noncentral norm mean h_K(t)=E||Z+t e_1|| satisfies h_K(t)>=t by Jensen,
and its derivative has absolute value at most 1 by the norm's Lipschitz bound.
Therefore, for q>0,

    |d log h_K(sqrt(q))/dq| <= 1/(2q).

Across the angular domain, q varies by at most 2rho=O_p(sqrt(K)), while its
minimum is of order K. The entire angular range of log h_K(||S_theta||) is
therefore O_p(K^(-1/2)), tending uniformly to zero. Consequently

    log p_J(theta | X,Y)
       = a data-dependent constant
         + (1/2)[-G cos(2theta)+C sin(2theta)] + o_p(1),                (6)

uniformly in theta. The density is with respect to dtheta, including the
Jacobian in (4). These statements concern ratios of angular densities;
the data-dependent normalizing constant is absorbed above.

### Mode locations and local rates

The posterior density in theta is a strictly increasing function of
||S_theta||^2: both exp(||S||^2/2) and h_K(||S||) increase with ||S||.
Its mode therefore maximizes -G cos(2theta)+C sin(2theta).

If C>0, the mode is interior, with

    2 theta_* = atan2(C,-G) in (0,pi).

Writing delta=theta-theta_*, the angular term is rho cos(2delta). From (6),
the local log density equals its maximum minus rho delta^2 to leading order.
Away from boundary-transition events, the local posterior variance is
1/(2rho), and its width is of order K^(-1/4).

If C<0, the preferred boundary is theta=0 when G<0, and theta=pi/2 when G>0.
At theta=0, for example,

    (1/2)[-G cos(2theta)+C sin(2theta)]
       = -G/2 + C theta + G theta^2 + higher terms.

Since C<0 has magnitude of order sqrt(K), the inward distance has leading
exponential scale 1/|C|, of order K^(-1/2). The same holds at pi/2. The inferior
boundary is downweighted by exp(-|G|), so it does not keep substantial posterior
mass unless the boundaries are nearly tied. Events with coefficients so small
that these rates cross over have probabilities tending to zero. They can be
noticeable at finite K and must not be omitted from numerical checks.

The limiting angular mode law is

    (1/2) Uniform(0,pi/2) + (1/4) delta_0 + (1/4) delta_(pi/2).        (7)

This follows from the rotationally symmetric limit of (G,C). In fact the
isotropic Gaussian observation model makes the orientation of (G,C) uniform
already at finite K: rotating the two observation channels rotates this pair
by twice the channel angle. Thus the sign probabilities and conditional
interior-mode uniformity are not generated by genuine mass information.

### Credible coverage tends to zero

The posterior concentrates about its random preferred angle or boundary, and
its central equal-tailed 95% angular interval has width tending to zero in
probability. For any fixed interior theta_0=atan(r_0), the limiting mode law
(7) has no atom at theta_0. Fix a small neighborhood of theta_0. Outside an
event whose probability tends to zero, a credible interval containing theta_0
must have its preferred angle in that neighborhood. Letting the neighborhood
shrink proves

    P_(r_0,V=0){r_0 lies in the central 95% Jeffreys posterior interval}
       -> 0.

Monotonicity of r=tan(theta) makes angular and mass-interval membership
equivalent. This is a fixed positive true mass; zero or infinite mass are not
included as possible truths in this statement.

There is no contradiction with posterior properness, with the surviving
m^(-2) tails at every finite K, or with the unchanged posterior at the exact
observed-zero dataset. Under true zero excitation, an exactly observed-zero
dataset has probability zero. There is also no contradiction with regular
Bayesian asymptotics: the true null is nonidentified/singular for mass and
the number of unknown nuisance coordinates grows with K.

The exact confidence constructions retain their 95% sampling coverage at this
null. The collapse established here concerns the equal-tailed intervals of
this joint Jeffreys posterior, not all possible uncertainty procedures.
