"""Constructive Gaussian/coarea bridge and conditional-score-law audit.

Scratch, not a final proposed mass definition. Managed snapshot:
VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.

ONE COMPONENT: AN EXACT JOINT LAW
Standardize x=F/sigma_F,y=a/sigma_a and s=sigma_F/sigma_a. Gaussian
location inversion supplies independent X*~N(x,1),Y*~N(y,1). Conditional
on X*Y*>0, write X*=v sin(theta),Y*=v cos(theta), 0<theta<pi/2,v real.
The ordinary change-of-variables Jacobian is |v|. Thus
 q(theta,v) proportional to |v| exp[-B(theta)^2/2-(v-A(theta))^2/2],
 A=x sin(theta)+y cos(theta), B=x cos(theta)-y sin(theta).
The mass-angle law is phi(B)*mu_1(A), normalized over positive angles,
where mu_1(A)=E|N(A,1)|. The ratio-of-means mass is exactly
 s*integral sin(theta)phi(B)(1+A^2)dtheta
     /integral cos(theta)phi(B)(1+A^2)dtheta.
All finite data, including x=y=0, give a proper joint law and finite
positive point. The ordinary conditioning event has strictly positive
probability, so this operation has no missing-root ambiguity.

GENERAL GAUSSIAN COAREA BRIDGE
For M>1, compatibility is a measure-zero manifold in 2M-dimensional
observation space. The declared noise metric supplies a precise extension:
use Euclidean Hausdorff volume after whitening. In one noise pair its
surface element is ||v|| dv dtheta, since the derivative Gram is
diag(I_M,||v||^2). The mass law becomes exp(-Q/2)*mu_M(A). This is also
the Euclidean generalized-fiducial density for the additive Gaussian data
equation. It supplies a joint law but retains the leading ours27 width
error: one radial power changes only an O(1) log-likelihood term when
R^2 is proportional to M. It is not a sampling-exact probability bridge.

For arbitrary SPD Sigma, native coordinates (m,a) give
 G=B'W B, K=G^-1, H=E'(W-W B K B'W)E, W=Sigma^-1, E=dB/dm.
The whitened manifold Jacobian is
 J(m,a)=sqrt(det G)*sqrt(a'H a).
Integrating a yields p(m) proportional to
 exp(-Q(m)/2) E_[N(b,K)] sqrt(a'H a).
The ratio-of-means numerator/denominator use
 exp(-Q/2) E_[N(b,K)] [||a||sqrt(a'H a)] times m / 1.
The joint law is total and proper: as m->infinity, K=O(m^-2),
b=O(m^-1),H=O(m^-2), so its marginal mass weight is O(m^-2) and the
readout weights are O(m^-3) and O(m^-2). Near zero they are finite.
Swap, units and joint coordinate rotations act geometrically on the same
whitened manifold; composition holds for the ratio-of-means functional on
a shared joint law. Fixed-coordinate regrouping does not change it.
However, the companion coarea_known_empty.py rejects this candidate on a
stronger continuity test: adding an observed-zero axis and letting both
of its SDs tend to zero does not recover the original lower-dimensional
readout. For F=2,a=0 with unit noise, 2.5279113099 becomes 2.4780610395
in that limit. Whitening keeps the disappearing axis's standardized
candidate variance in the surface Jacobian. Properness and fixed-M
coordinate invariance therefore do not suffice.

WHAT COAREA DOES AND DOES NOT FORCE
The above is unique once the additive Gaussian data equation, whitening,
Euclidean distance and Hausdorff conditioning are specified. Positivity
and composition alone do not choose these ingredients. Tubes of radius
epsilon*w(p) around each manifold point p have limiting measure
w(p)^(M-1) times Hausdorff measure; every positive smooth w gives the same
zero set and a different conditional law. The noise metric motivates
w=constant; it does not make this law's many-reading widths correct.
Generalized fiducial inference also depends on its data equation and norm:
Murph, Hannig and Williams (2023), theorem1.1 and surrounding discussion,
https://arxiv.org/html/2302.14598 . Its formula gives exactly our Jacobian
for the additive Gaussian data equation and l2 norm. It explicitly does
not promise sampling-exact coverage in general.

CONDITIONAL NORMAL SCORE AS A DENSITY
An attractive alternative is p(theta) proportional to phi(Z)|Z'|, with
H=||A||^2,C=A'B,Z=C/sqrt(H), and
 Z'=(det(S)-H^2)/H^(3/2).
Near a well-identified principal line its curvature/width is correct:
1/(Z')^2=lambda1/(lambda1-lambda2)^2. But the all-branch density also
hears the MINOR line. On the full projective circle, every interior
normal pivot value has two inverse branches; normalizing this expression
is an equal-branch allocation plus truncation to the observed pivot range.
The integral is exactly 2*(2*Phi(sqrt(lambda1)-sqrt(lambda2))-1).

Mandatory positive mass does not fix the allocation: with negative cross
product the positive quadrant contains the minor line. The score density
then concentrates at the WORST-fitting positive angle rather than the
least-impossible explanation. Concrete S=[[100,-80],[-80,100]]: at theta
pi/4 it favors H=20,Q=180, whereas endpoint limits give H=100,Q=100.
The likelihood ratio against the favored angle is exp(40).

At S=cI, both Z and Z' vanish at every angle and the density is undefined.
There is no continuous unique extension: approaching with diagonal
anisotropy gives density proportional to |cos(2theta)|; approaching with
off-diagonal anisotropy gives density proportional to sin(2theta).

Even M1 exposes the mismatch: score density phi(B)|A| differs from the
proper Gaussian-inversion law phi(B)E|N(A,1)|. Replacing the former factor
by the latter repairs M1 totality by restoring candidate-size uncertainty;
the corresponding general-M coarea extension is precisely the law above,
whose many-reading width remains wrong. This is the concrete bridge and
the precise place where it falls short.
"""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss


def phi(x):
    return np.exp(-np.asarray(x)**2/2)/np.sqrt(2*np.pi)


def Phi(x):
    return (1+math.erf(x/np.sqrt(2)))/2


def folded_mean(x):
    x=np.asarray(x)
    erfs=np.array([math.erf(v/np.sqrt(2)) for v in x.ravel()]).reshape(x.shape)
    return 2*phi(x)+x*erfs


def grid(a,b,n=800):
    z,w=leggauss(n)
    return a+(b-a)*(z+1)/2,w*(b-a)/2


def one_component_checks():
    theta,w=grid(0,np.pi/2)
    sn,cs=np.sin(theta),np.cos(theta)
    rows=[]
    for x,y in ((0.,0.),(1.,0.),(2.,1.),(3.,-1.),(-2.,-1.)):
        a=x*sn+y*cs;b=x*cs-y*sn
        angular=phi(b)*folded_mean(a)
        probability=float(w@angular)
        point=float(w@(sn*phi(b)*(1+a*a))/(w@(cs*phi(b)*(1+a*a))))
        probability_exact=Phi(x)*Phi(y)+Phi(-x)*Phi(-y)
        xplus=x*Phi(x)+float(phi(x));xminus=float(phi(x))-x*Phi(-x)
        yplus=y*Phi(y)+float(phi(y));yminus=float(phi(y))-y*Phi(-y)
        x_size=xplus*Phi(y)+xminus*Phi(-y)
        y_size=yplus*Phi(x)+yminus*Phi(-x)
        exact_point=x_size/y_size
        assert abs(probability-probability_exact)<1e-12
        assert abs(point-exact_point)<1e-11
        rows.append({"x":x,"y":y,"compatibility_probability":probability,
                     "ratio_mean_point":point,"point_closed_form":exact_point,
                     "probability_error":abs(probability-probability_exact),
                     "point_error":abs(point-exact_point)})
    return rows


def jacobian_checks():
    rng=np.random.default_rng(290928)
    errors=[]
    for dimension in (1,2,3,6):
        for _ in range(6):
            raw=rng.normal(size=(2*dimension,2*dimension))
            sigma=raw@raw.T+np.eye(2*dimension)
            precision=np.linalg.inv(sigma)
            mass=float(np.exp(rng.normal()))
            acc=rng.normal(size=dimension)
            b=np.vstack((mass*np.eye(dimension),np.eye(dimension)))
            e=np.vstack((np.eye(dimension),np.zeros((dimension,dimension))))
            g=b.T@precision@b;k=np.linalg.inv(g)
            hh=e.T@(precision-precision@b@k@b.T@precision)@e
            j=np.sqrt(np.linalg.det(g)*(acc@hh@acc))
            derivative=np.column_stack((b,e@acc))
            direct=np.sqrt(np.linalg.det(derivative.T@precision@derivative))
            errors.append(float(abs(j/direct-1)))
    assert max(errors)<1e-11
    return max(errors)


def score_density(theta,p,q,d):
    sn,cs=np.sin(theta),np.cos(theta)
    h=p*sn*sn+q*cs*cs+2*d*sn*cs
    c=(p-q)*sn*cs+d*(cs*cs-sn*sn)
    z=c/np.sqrt(h)
    zp=(p*q-d*d-h*h)/h**1.5
    return phi(z)*abs(zp),h,z


def score_checks():
    theta,w=grid(0,np.pi,2000)
    normalizers=[]
    for l1,l2 in ((4.,1.),(25.,4.),(180.,20.)):
        density,_,_=score_density(theta,l1,l2,0.)
        numerical=float(w@density)
        exact=2*(2*Phi(np.sqrt(l1)-np.sqrt(l2))-1)
        normalizers.append({"lambda1":l1,"lambda2":l2,"numerical_normalizer":numerical,
                            "exact_normalizer":exact,"absolute_error":abs(numerical-exact)})
    assert max(r["absolute_error"] for r in normalizers)<1e-5
    theta,w=grid(0,np.pi/2,2000)
    density,h,z=score_density(theta,100.,100.,-80.)
    probability=density*w;probability/=sum(probability)
    central=float(sum(probability[abs(theta-np.pi/4)<.1]))
    assert central>.99
    eps=1e-5
    diagonal,_,_=score_density(theta,1+eps,1-eps,0.)
    offdiag,_,_=score_density(theta,1.,1.,eps)
    diagonal/=w@diagonal;offdiag/=w@offdiag
    region=abs(theta-np.pi/4)<.1
    path1=float(sum(w[region]*diagonal[region]))
    path2=float(sum(w[region]*offdiag[region]))
    assert abs(path1-path2)>.15
    return {"full_circle_normalizers":normalizers,
            "negative_cross_example":{"P":100.,"Q":100.,"D":-80.,
              "score_law_probability_within_point_one_rad_of_worst_angle":central,
              "score_favored_angle":np.pi/4,"H_at_favored_angle":20.,"residual_Q_at_favored_angle":180.,
              "residual_Q_at_positive_domain_endpoint":100.,"endpoint_likelihood_ratio_vs_favored_angle":math.exp(40)},
            "isotropic_limit_path_dependence":{"epsilon":eps,
              "diagonal_path_probability_near_pi_over_4":path1,
              "offdiagonal_path_probability_near_pi_over_4":path2}}


def main():
    out={"source_commit":"b4bfc2040031ef010416334d0c5f4aa2ad05cbd8",
         "one_component_Gaussian_conditioning":one_component_checks(),
         "full_covariance_coarea_Jacobian_relative_error":jacobian_checks(),
         "conditional_score_as_density":score_checks(),
         "status":"Gaussian/coarea is a precise total joint law but fails known-empty continuity and many-reading widths; normalized pivot-Jacobian law also fails as a general bridge."}
    Path(__file__).with_name("bridge_candidates_results.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
