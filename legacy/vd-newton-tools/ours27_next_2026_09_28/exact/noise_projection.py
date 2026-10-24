"""Audit of Gaussian-noise inversion followed by nearest law projection.

Scratch only. Managed context VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.
One noise pair suffices to exhibit the obstructions. Standardized data form
an M x 2 matrix Y. The compatible means are v u_theta^T, with
u_theta=(sin theta,cos theta), theta in the OPEN interval (0,pi/2).

Nearest projection and totality
------------------------------
For fixed theta, nearest v is Y u_theta and the squared residual is
tr(S)-u_theta^T S u_theta, S=Y^T Y. Write S=[[P,D],[D,Q]]. If D<0, the
maximum of P sin^2(theta)+Q cos^2(theta)+2D sin(theta)cos(theta) over the
closed positive quadrant is at theta=0 or pi/2. The corresponding infimum
over the OPEN quadrant is not attained. Thus no nearest candidate with
positive finite mass exists. D<0 occurs with positive Gaussian probability
for every finite observed Y in a fresh-noise inversion Y*=Y+E. At Y=0 it
has probability exactly 1/2. Closing the cone gives forbidden mass atoms at
zero and infinity. If D=0 and P=Q, every angle ties; probability zero does
not remove an every-input totality requirement.

Why ordinary perturb-and-TLS overstates conditional widths
---------------------------------------------------------
Rotate the TWO standardized channels to observed principal axes a,b, where
a'b=0, h=||a||^2>q=||b||^2. Fresh independent N(0,I_M) noise xi,eta gives

 C*=(a+xi)'(b+eta)=a'eta+b'xi+xi'eta,
 E[C*|Y]=0, Var(C*|Y)=h+q+M, EXACTLY.

Every cross-covariance of the three displayed terms vanishes. The mean
scatter gap remains h-q because both diagonals gain M. First-order TLS
angular variance is (h+q+M)/(h-q)^2. In regular sequences R^2=lambda*M,
h/M -> lambda+1 and q/M ->1, so this becomes
(R^2+3M)/R^4. The actual repeated-experiment TLS variance is
(R^2+M)/R^4. Their ratio is (lambda+3)/(lambda+1), not one.
Removing the quadratic fresh-noise term xi'eta is insufficient: h+q is
still asymptotically R^2+2M. Noise scale tuning could match one formula,
but would be a new data-dependent choice, not inversion of declared noise.

A constructive conditional replacement
---------------------------------------
Holding a fixed and REPLACING the perpendicular residual by fresh eta,
rather than adding noise to both observed channels, gives H*=h,
Q*=||eta||^2, C*=a'eta. The exact conditional variance is h and expected
gap h-M. Thus its first-order angle variance h/(h-M)^2 has the correct
(R^2+M)/R^4 limit. This is not arbitrary score tuning: it follows the
known-noise conditional experiment used by conditional_lr.py.

But centering this simulation at the fitted theta is a plug-in conditional
bootstrap, not an exact fiducial law. For finite-sample inversion the
conditioning statistic h(theta)=||Y u_theta||^2 changes with the TESTED
theta. conditional_lr.py retains that dependence and gives exact sets.
A single pushforward at fitted theta drops it. Boundary atoms still occur.

The exact conditional score C(theta)/sqrt(H(theta)) is standard normal at
the tested angle, but as a function of theta it is generally nonmonotone,
bounded, and can have missing/multiple roots for sampled normal values.
More precisely, for observed scatter eigenvalues lambda1>=lambda2,
 C(theta)^2/H(theta)=(lambda1-H)*(H-lambda2)/H
                   <=(sqrt(lambda1)-sqrt(lambda2))^2.
The maximum is at H=sqrt(lambda1*lambda2). A fresh standard normal draw
therefore has no inverse on even the WHOLE projective circle with
probability at least 2*Phi(-(sqrt(lambda1)-sqrt(lambda2)))>0. Restricting
to positive angles only increases failure. Isotropic scatter makes the
score identically zero, so almost every fresh normal draw has no root.
Turning those inversions into a total proper positive mass distribution
requires root weights, branch selection, or conditioning. None is forced
merely by knowing the Gaussian error covariance.

Candidate geometry supplies another operation, conditioning a Gaussian
ambient law on the rank-one manifold using Hausdorff surface volume. Its
Jacobian is ||v||: d(surface)=||v|| dv dtheta, because derivatives with
respect to v have identity Gram and the theta derivative has norm||v||.
This is NOT the pushforward of nearest projection. Relative to the flat
candidate measure it adds one radial power, an O(1) log-likelihood term;
under R^2 proportional to M it retains leading angular width 1/R rather
than the needed sqrt(R^2+M)/R^2. Thus that canonical surface choice also
does not settle the many-reading problem.
"""
from pathlib import Path
import json
import numpy as np


def full_perturbation(rng, dimension, h, q, trials):
    """Exact scatter of fixed orthogonal columns plus Gaussian noise."""
    xi = rng.normal(size=(trials, 2))
    eta = rng.normal(size=(trials, 2))
    n = dimension - 2
    w11 = rng.chisquare(n, trials)
    z = rng.normal(size=trials)
    w12 = np.sqrt(w11) * z
    w22 = z*z + rng.chisquare(n-1, trials)
    hs = h + 2*np.sqrt(h)*xi[:, 0] + np.sum(xi*xi, axis=1) + w11
    qs = q + 2*np.sqrt(q)*eta[:, 1] + np.sum(eta*eta, axis=1) + w22
    cs = np.sqrt(h)*eta[:, 0] + np.sqrt(q)*xi[:, 1] + np.sum(xi*eta, axis=1) + w12
    delta = .5*np.arctan2(2*cs, hs-qs)
    return cs, delta


def conditional_replacement(rng, dimension, h, trials):
    z = rng.normal(size=trials)
    q = z*z + rng.chisquare(dimension-1, trials)
    c = np.sqrt(h)*z
    delta = .5*np.arctan2(2*c, h-q)
    return c, delta


def main():
    rng = np.random.default_rng(280929)
    rows=[]
    trials=50000
    for dimension in (500, 5000):
        for strength in (.5, 1., 3.):
            energy=strength*dimension
            h=energy+dimension
            q=float(dimension)
            c,delta=full_perturbation(rng,dimension,h,q,trials)
            cc,dd=conditional_replacement(rng,dimension,h,trials)
            full_cross_ratio=float(np.var(c)/(h+q+dimension))
            conditional_cross_ratio=float(np.var(cc)/h)
            assert abs(full_cross_ratio-1)<.03
            assert abs(conditional_cross_ratio-1)<.03
            rows.append({
                "M":dimension,"R_squared_over_M":strength,"trials":trials,
                "full_exact_cross_variance":h+q+dimension,
                "full_empirical_cross_variance_over_exact":full_cross_ratio,
                "full_angular_variance":float(np.var(delta)),
                "full_first_order_angular_variance":(h+q+dimension)/(h-q)**2,
                "conditional_empirical_cross_variance_over_exact":conditional_cross_ratio,
                "conditional_angular_variance":float(np.var(dd)),
                "conditional_first_order_angular_variance":h/(h-dimension)**2,
                "required_repeated_experiment_angular_variance":(energy+dimension)/energy**2,
                "full_asymptotic_overdispersion_factor":(strength+3)/(strength+1),
            })
    boundary=[]
    for dimension in (1, 3, 21):
        x=rng.normal(size=(trials,dimension))
        y=rng.normal(size=(trials,dimension))
        p=np.sum(x*x,axis=1);q=np.sum(y*y,axis=1);d=np.sum(x*y,axis=1)
        zero=(d<0)&(q>p)
        infinity=(d<0)&(p>q)
        boundary.append({"M":dimension,"trials":trials,
                         "no_positive_finite_projection_fraction":float(np.mean(d<0)),
                         "closed_projection_mass_zero_fraction":float(np.mean(zero)),
                         "closed_projection_mass_infinity_fraction":float(np.mean(infinity)),
                         "exact_no_finite_projection_probability":.5})
    out={"source_commit":"b4bfc2040031ef010416334d0c5f4aa2ad05cbd8","seed":280929,
         "variance_checks":rows,"zero_data_projection_boundary_checks":boundary,
         "conclusion":"Raw perturb-and-project overdisperses and fails positive finite candidate totality; conditional replacement fixes leading width only, not exact proper-law construction."}
    Path(__file__).with_name("noise_projection_results.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
