# Conditioning an ambient Gaussian onto positive collinearity

Working mathematical contribution, 2026-09-15. No production or managed edits.
This derivation starts from an explicitly supplied *probability law* in
Cartesian force/acceleration coordinates. A Gaussian measurement likelihood
alone does not supply that law; see the final section.

## 1. Set-up and results

Let z=(F,a) belong to R^d x R^d and have a nonsingular Gaussian density g(z)
with respect to Cartesian Lebesgue volume. Its covariance may be correlated;
the simple tube metric below is separately specified. For d>=2 the positive
collinearity manifold

    C = {(f u, alpha u): f,alpha>0, u in S^(d-1)}

has dimension d+1 and ambient probability zero. Two regularizations give two
different conditional limits:

1. Condition on angle(F,a)<epsilon, with epsilon tending to zero:

       dP_angle proportional to
       g(fu,alpha u) (f alpha)^(d-1) df d alpha dOmega(u).

2. Choose positive scalar scales sF,sa and condition on Euclidean distance
   less than epsilon to C in coordinates x=F/sF, y=a/sa:

       dP_tube proportional to
       g(fu,alpha u)
       (f^2/sF^2 + alpha^2/sa^2)^((d-1)/2)
       df d alpha dOmega(u).

Constants independent of (f,alpha,u) cancel on normalization. Both limits are
proper because their polynomial weights are Gaussian-integrable. They are not
two coordinate expressions for one conditional probability: the shrinking
events differ.

## 2. Angular tolerance: direct polar derivation

Write F=f u and a=alpha v, where f,alpha>0 and u,v lie on the unit sphere.
The ambient volume is

    dF da = f^(d-1) alpha^(d-1) df d alpha dOmega(u) dOmega(v).

For fixed u the cap angle(u,v)<epsilon has area

    kappa_(d-1) epsilon^(d-1) [1+O(epsilon^2)],

where kappa_j is the volume of the j-dimensional unit ball. Its normalized
average of any continuous test integrand tends to its value at v=u. Divide
the numerator and denominator of the conditional expectation by this cap
area. The displayed angular limit follows.

Dominated convergence is legitimate: a nonsingular, possibly noncentral
Gaussian admits a bound C exp[-c(f^2+alpha^2)] uniformly in u,v, after reducing
c if necessary. This bound times (f alpha)^(d-1) is integrable. The vector-zero
sets have probability zero, so defining their angle arbitrarily does not
change the finite-epsilon events or limit. There is no limiting atom on either
axis: the stated smooth limiting magnitude density already normalizes to one.

## 3. Fixed standardized-distance tube: geometry and boundary control

Put r=f/sF, s=alpha/sa, R=sqrt(r^2+s^2). In standardized Cartesian coordinates
the embedding is Phi(r,s,u)=(r u,s u). Since u dot du=0,

    |dPhi|^2 = dr^2 + ds^2 + (r^2+s^2)|du|^2.

Its induced (d+1)-dimensional surface volume is

    dH_(d+1) = R^(d-1) dr ds dOmega(u).

This explains the proposed tube weight. The following calculation also
controls the missing axes and the singular vertex rather than assuming they
make no contribution.

### Interior normal coordinates

For w perpendicular to u define

    x = r u + (s/R) w,
    y = s u - (r/R) w.

When |w|<R and r,s>0, the closest point on the closed positive-collinearity cone
is (r u,s u), and the distance is |w|. To see this, treat (x,y) as a d-by-2
matrix. It has orthogonal left factors u,w and orthogonal right factors
(r,s),(s,-r); its two singular values are R and |w|. The larger rank-one
component has positive right coefficients. Also

    x dot y = r s (1-|w|^2/R^2)>0.

Conversely, if x dot y>0, the largest eigenvector of the 2-by-2 Gram matrix has
two coefficients of the same sign. The best rank-one approximation therefore
belongs to the interior of the positive cone, uniquely except on a null
degeneracy set, and has these normal coordinates.

Their exact volume Jacobian, relative to dr ds dOmega(u) dw, is

    J = R^(d-1) (1-|w|^2/R^2).

One elementary determinant check: set r=R cos(theta), s=R sin(theta), use a
local orthonormal tangent frame at u, and rotate the two output channels by
the angle theta. At the chosen point the variations are

    dR:       (u,0),
    dtheta:   (w,-R u),
    dt_i:     (R e_i,-w_i u),
    dw_i:     (0,e_i).

The only nontrivial determinant block is

    [ R I   w ]
    [ -w^T -R ],

whose absolute determinant is R^d(1-|w|^2/R^2). Dividing by the planar polar
Jacobian dr ds=R dR dtheta gives J above.

For the interior tube |w|<epsilon, rescale w=epsilon v. Divide its integral by
kappa_(d-1) epsilon^(d-1). For fixed r,s>0, the integrand converges to the
ambient density at (r u,s u) times R^(d-1). The indicator |w|<R converges to one.
Because |(x,y)|^2=R^2+|w|^2 and 0<J<=R^(d-1), a Gaussian envelope supplies the
integrable dominating function C R^(d-1) exp(-cR^2). Thus the interior integral
has exactly the surface-density limit claimed above. This also controls
arbitrarily small R; no assumption of uniform curvature near the vertex is
being made.

### Axis-projected part of the tube

If x dot y<=0, a rank-one fit with nonnegative right coefficients maximizes

    |x|^2 c1^2 + |y|^2 c2^2 + 2(x dot y)c1 c2,
    c1,c2>=0, c1^2+c2^2=1.

The maximum is attained on an axis and equals max(|x|^2,|y|^2). Consequently
the distance to the closed positive cone is min(|x|,|y|). Its epsilon-tube
portion is contained in {|x|<epsilon} union {|y|<epsilon}; the Gaussian mass
there is O(epsilon^d). This is subleading relative to epsilon^(d-1), so the
limiting conditional law has no axis atom. The origin lies in this same bound;
its radius-O(epsilon) ball separately has mass O(epsilon^(2d)). Distance to
the open cone and its closure is identical, so using either in the tube event
does not alter this result.

### Restoring physical units

The standardized density is g_s(x,y)=sF^d sa^d g(sF x,sa y), while
dr ds=df d alpha/(sF sa). All remaining scale factors are constant, leaving
the weight (f^2/sF^2+alpha^2/sa^2)^((d-1)/2).

## 4. General metric: do not silently substitute full whitening

For a fixed positive-definite Cartesian metric G, replace Euclidean distance
by distance after a linear map L with L^T L=G. Let E be a d-by-(d-1)
orthonormal tangent frame at u and let

    D Phi = [ (u,0), (0,u), (f E,alpha E) ].

The induced tube weight, apart from a constant ambient determinant, is

    sqrt(det[(D Phi)^T G (D Phi)]).

For G=diag(sF^(-2) I,sa^(-2) I) this is

    (sF sa)^(-1)
    (f^2/sF^2+alpha^2/sa^2)^((d-1)/2),

as above. If G has isotropic blocks [[a I,c I],[c I,b I]], its nonconstant
factor is (a f^2+2c f alpha+b alpha^2)^((d-1)/2). A general anisotropic G also
introduces direction dependence through the Gram determinant.

Therefore, when 'standardized' means whitening the *full* joint covariance,
the simple two-scale formula applies only to the block-isotropic special
case. Separate scalar rescalings preserve angles within each vector; general
whitening can mix or distort them. Physical positive collinearity must be
transformed along with the coordinates.

## 5. Null uncertainty sharply distinguishes the two limits

For an independent isotropic ambient Gaussian centered at (0,0), let s0=sF/sa.

Angular conditioning gives independent standardized magnitudes U,V~chi_d.
For t=M/s0=U/V,

    p_angle(t) = [2 Gamma(d)/Gamma(d/2)^2]
                 t^(d-1)/(1+t^2)^d,       t>0,

    Var[log(M/s0)] = (1/2) trigamma(d/2).

The point readout is still s0 by channel symmetry. At d=2 the central 95%
ratio interval is [sqrt(.025/.975),sqrt(.975/.025)] rather than the
half-Cauchy interval.

For tube conditioning, use r=R cos(theta), s=R sin(theta). The null density
relative to dR dtheta is proportional to R^d exp(-R^2/2), with theta uniform
on (0,pi/2). The radial variable is chi_(d+1), but the ratio is half-Cauchy:

    p_tube(t)=2/[pi(1+t^2)],
    Var[log(M/s0)]=pi^2/4.

Thus the tube changes the latent radial distribution relative to the flat
positive-magnitude model while preserving its entire null mass-ratio law.

The tail distinction persists at arbitrary finite Gaussian observations. For
angular conditioning, integration over direction gives positive axis limits
for the Gaussian factor, and change of variables f=m alpha yields

    p_M(m) ~ C0 m^(d-1)        as m downarrow 0,
    p_M(m) ~ Cinf m^(-d-1)     as m -> infinity,

with finite positive coefficients. Exactly the power moments with -d<q<d
are finite. In d=3 the ordinary mass variance is therefore finite. The tube
weight remains positive at either nonzero axis, so its induced mass density
instead retains a positive limit at zero and a positive m^(-2) tail. Its
ordinary mass and inverse-mass means diverge, despite the same null point.

For d=1 the collinearity condition is simply equal signs, with positive
probability; the two weight exponents are zero and the two constructions agree.

## 6. Interpretation: what has and has not been selected

- A Gaussian likelihood p(observation | F,a) is a density in observation
  space. Using the same Gaussian expression as a probability law over latent
  Cartesian (F,a) is an additional inference construction. With fixed known
  covariance it is obtained, for example, from a flat Cartesian prior and
  normalization. It does not follow from Gaussian error alone.
- If the ambient latent law is explicitly supplied, both conditioning limits
  above follow mathematically. Neither limit is selected merely by saying
  'condition on exact Newtonian collinearity': that is a probability-zero
  event for d>=2, and its limiting interpretation must be specified.
- Different tolerance variables encode different experimental/operational
  questions. A relative/angular tolerance and an absolute standardized
  distance tube need not agree. Their difference is not a parameterization
  error.
- Both results distinguish positive collinearity from allowing the two
  opposite directions. Conditioning on absolute collinearity would add the
  antiparallel branch and would define a different problem.
- There are no boundary atoms under these nonsingular Gaussian limits. A
  model with an actual atom at zero excitation or a singular covariance
  requires a separate analysis.

The results are direct polar, tangent-metric, Jacobian, and dominated-limit
calculations, rather than a literature-priority claim. No claim of a uniquely
correct conditioning convention or empirically preferred estimator follows.
