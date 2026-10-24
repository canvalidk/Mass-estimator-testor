# Fisher volume and the tube weight

Working derivation, 2026-09-15. Scratch only; no production or managed edits.
No new managed source snapshot was resolved for this task.

## Result and scope

For the sampling model

    y | theta ~ N(Phi(theta), Sigma),
    Phi(f,alpha,u)=(f u,alpha u),
    f,alpha>0, u in S^(d-1), Sigma known, positive definite and theta-independent,

the expected Fisher-information metric is the pullback of the Cartesian
Mahalanobis metric Sigma^(-1). Its volume is therefore exactly the reference
weight obtained from a shrinking distance tube in that same metric.

This is the **joint Jeffreys-rule prior** for all d+1 latent parameters. It is
not automatically a mass-only reference prior, a unique uncertainty law, or
an optimal estimator. A tube defined with a different metric has a different
volume and need not be the sampling model's Fisher volume.

## 1. Direct Fisher derivation

Use local sphere coordinates q and let J=D_theta Phi. The log likelihood is

    ell(theta;y)=constant - (y-Phi(theta))^T Sigma^(-1)(y-Phi(theta))/2.

Its score and expected score product are

    score = J^T Sigma^(-1)(y-Phi),
    I(theta) = E_theta[score score^T] = J^T Sigma^(-1)J.

Thus

    dPi_J(theta) proportional to sqrt(det I(theta)) dtheta.

Under theta=h(eta), I_eta=(Dh)^T I_theta Dh. Its determinant square root picks
up exactly |det Dh|, so this construction defines one coordinate-independent
measure. The observed likelihood Hessian is a different object: for a nonlinear
Phi it also contains residual-times-second-derivative terms. Substituting that
data-dependent Hessian would not be the derivation here.

If Sigma depends on theta, the Fisher matrix additionally contains

    (1/2) tr(Sigma^(-1) Sigma_,i Sigma^(-1) Sigma_,j).

Unknown covariance parameters or fitted uncertainty scales therefore require
a separate model; the fixed-Sigma result cannot simply be asserted for them.

### Independent isotropic channels

Write Sigma=diag(sF^2 I_d,sa^2 I_d), H(q)=(D_q u)^T(D_q u), and

    K=f^2/sF^2+alpha^2/sa^2.

Because u is perpendicular to its derivatives,

    I = diag(sF^(-2), sa^(-2), K H(q)),
    det I = (sF sa)^(-2) K^(d-1) det H(q).

Since dOmega=sqrt(det H) dq,

    dPi_J = constant * (sF sa)^(-1)
            K^((d-1)/2) df d alpha dOmega.

This is exactly the earlier standardized-tube weight. More generally take an
orthonormal tangent frame E at u, so the derivative columns are
[(u,0),(0,u),(f E,alpha E)]. The weight relative to df d alpha dOmega is the
square root of their Gram determinant in Sigma^(-1). This includes anisotropy
and force/acceleration covariance without assuming a uniform direction weight.

## 2. What local distinguishability contributes

For two model parameters with the same known covariance, Gaussian KL divergence
is exactly

    KL(p_theta || p_theta')
      = (Phi(theta')-Phi(theta))^T Sigma^(-1)
        (Phi(theta')-Phi(theta))/2.

Locally this is (delta theta)^T I(theta)(delta theta)/2, up to higher orders.
The coordinate volume of a small Fisher ellipsoid

    (delta theta)^T I(theta)(delta theta) <= epsilon^2

is kappa_(d+1) epsilon^(d+1)/sqrt(det I(theta)). Jeffreys density cancels that
factor. To leading order it allocates equal mass to equally sized local cells
of distinguishable sampling distributions.

This is a defensible interpretation of the weighting: it uses the declared
experiment's resolving geometry to count parameter possibilities. It does not
mean the observed data contain extra information, or that posterior uncertainty
is automatically calibrated for the intended use.

In the isotropic model, K is the information about a small common direction
rotation. Either large force or large acceleration can make that rotation
distinguishable. The K^((d-1)/2) factor reflects the volume of these direction
possibilities. At f=alpha=0 the direction becomes unidentifiable and this
volume vanishes for d>1. At either nonzero axis the other channel can still
identify direction, so the weight is positive there. This helps explain why
the angular-tolerance weight (f alpha)^(d-1) has different boundary behavior.

## 3. Flat magnitudes and the change to (mass, acceleration vector)

The existing flat-magnitude measure is

    dPi_flat = df d alpha dOmega.

It gives equal reference weight to equal intervals of each positive magnitude
and uses uniform angular measure. It is distinct from Cartesian volume and
from Fisher surface volume. It is nonetheless a legitimate measure that can
be transformed into any coordinates. The objection is not that it somehow
fails change-of-variables calculus: declaring a new flat density after every
coordinate change would replace the measure.

Set m=f/alpha and v=alpha u, with alpha=||v||. Since

    df d alpha = alpha dm d alpha,
    dv = alpha^(d-1) d alpha dOmega,

the two joint measures become

    dPi_flat = ||v||^(2-d) dm dv,

    dPi_J proportional to
       ||v|| (m^2/sF^2+1/sa^2)^((d-1)/2) dm dv.

The omitted Jeffreys factor (sF sa)^(-1) is constant. In d=3, the flat density
is 1/||v||; the joint Jeffreys density is
||v||(m^2/sF^2+1/sa^2). Neither is a flat density on (m,v).

A direct determinant check in these coordinates uses Phi(m,v)=(mv,v):

    I = [ ||v||^2/sF^2        m v^T/sF^2 ]
        [ m v/sF^2    (m^2/sF^2+1/sa^2) I_d ],

whose determinant is

    ||v||^2/(sF^2 sa^2) * (m^2/sF^2+1/sa^2)^(d-1).

Both priors are improper on the complete unbounded latent space; Gaussian
likelihood restriction yields proper posteriors here. Integrating their joint
prior densities over all nuisance v gives divergent integrals, so the formal
m-dependent factors must not be presented as already normalized marginal
priors for mass. Choosing mass as the inferential target while treating v as
nuisance is a further problem, not solved by naming the joint Jeffreys rule.

## 4. Measurement design and sequential updating

First consider independent repeated observations of the SAME latent pair:

    y_i | theta ~ N(Phi(theta),Sigma_i).

For known fixed designs, information and precision add:

    G_total=sum_i Sigma_i^(-1),
    I_total=J^T G_total J.

The combined likelihood is equivalent, up to a theta-independent factor, to
one precision-weighted observation with covariance G_total^(-1). A Jeffreys
rule applied to that combined design gives

    pi_J,total(theta) proportional to sqrt(det[J^T G_total J]).

### A fixed prior remains sequentially coherent

If pi_0 is selected once, updating first by L1 and then by L2 gives exactly
L1 L2 pi_0, just as a batch calculation does. This includes a Jeffreys prior
selected for the entire planned design before observation. No inconsistency
in Bayes' rule has been found.

If instead the first analysis uses pi_J,1 and the later batch analysis replaces
it with pi_J,total, the resulting posteriors differ by pi_J,total/pi_J,1. That
is an added change of reference weighting, rather than multiplication by the
new measurement likelihood alone.

### When recomputation changes only a constant

For n identical designs, I_n=n I_1, so
pi_J,n=n^((d+1)/2) pi_J,1: the normalized prior shape and posterior are unchanged
by this recomputation. The same holds when every Sigma_i is a scalar multiple
of one fixed covariance matrix. It also holds for a common rescaling of both
isotropic channel uncertainties.

Changing their relative precision generally changes the shape. In d=3, let
the first standardized experiment have channel precisions (1,1). Its prior
weight is f^2+alpha^2. Add an experiment with precisions (9,1). The combined
weight is, apart from a constant,

    10 f^2+2 alpha^2.

The ratio to the first weight varies from 2 to 10 as f/alpha varies, so this is
not an irrelevant normalization. At f=alpha the ratio is 6; at f=3 alpha it
is 9.2. This example depends only on changed apparatus precision, not observed
values.

There is an exact readout contrast if both experiments observe zero vectors.
The combined likelihood is exp[-(10 f^2+2 alpha^2)/2]. With a general old
quadratic weight a f^2+b alpha^2 and combined precisions A,B, half-normal
moment identities give

    R_old = sqrt(B/A) * (2a/A+b/B)/(a/A+2b/B).

Keeping the original prior (a=b=1) gives R_old=(7/11)/sqrt(5), approximately
0.284590. Recomputing combined-design Jeffreys (a=A=10,b=B=2) gives
R_new=1/sqrt(5), approximately 0.447214. Both use the same combined likelihood;
the extra change comes entirely from the prior redesign.

### What an operational policy must decide

An investigation can choose one reference design and hold its prior fixed,
choose the full planned experiment's Jeffreys prior once, or explicitly report
a design-dependent reference analysis. Those choices have different meanings.
One should not claim ordinary prior-to-posterior updating while silently
replacing the prior after each design change.

Adaptive designs can make the selected design itself depend on earlier data;
stopping rules can also alter expected Fisher information. Such rules require
a declared sequential experiment. Merely plugging the latest fitted covariance
into a new prior is not an established sequential model. Conversely, for a
fixed prior, an ignorable design-selection rule does not license adding a new
prior factor. Different trials that share mass but have different latent
accelerations constitute another parameter model; the same-pair derivation
above does not cover them.

## 5. Relationship to the ambient-Gaussian conditioning construction

If a flat Cartesian ambient prior is used with the fixed-covariance Gaussian
measurement likelihood, the ambient posterior is the normalized Gaussian
expression in latent Cartesian coordinates. Conditioning that posterior by a
shrinking Sigma^(-1)-distance tube produces likelihood times Fisher volume.
Thus it agrees with the constrained-model joint Jeffreys posterior derived
directly above.

The direct Fisher derivation does not require adopting an unconstrained ambient
prior or a probability-zero conditioning event. It starts with the constrained
sampling model. Equality of the final constructions is a mathematical result;
the justification for adopting one of them is still an inferential choice.

If the ambient Gaussian posterior incorporates an informative prior, its
covariance need not equal the sampling covariance Sigma. A tube in that posterior
covariance metric then does not automatically yield the sampling model's joint
Jeffreys prior. Keep the likelihood covariance, posterior covariance, and chosen
tube metric explicit.

## Sources and status

The identities above are direct score, determinant, Jacobian, KL and
likelihood-product derivations. A targeted primary-source check used Sun and
Berger, *Objective Bayesian analysis under sequential experimentation* (2008),
[published-paper record](https://arxiv.org/abs/0805.3064) and
[author-hosted text](https://www2.stat.duke.edu/~berger/papers/OB-sequential.pdf),
especially section 2.1. It gives the joint Jeffreys determinant rule and shows
that stopping-rule information can modify default priors; it also distinguishes
joint Jeffreys and reference-prior approaches in multiparameter problems.
These are supporting definitions and context, not an attribution of this mass
derivation or evidence that the Jeffreys choice is uniquely appropriate here.

The 1946 Jeffreys publisher text was sought but inaccessible through the browser;
it was not relied on. No literature-priority or exhaustive survey claim is made.
