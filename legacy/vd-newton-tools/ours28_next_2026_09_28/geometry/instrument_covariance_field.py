"""Constructive, continuous instrument-covariant size/shape field.

Candidate only, not an adopted equation. Sources read: all four files in
../sources, VD-docs@592d144531fa4a30d1c2229a4c3b528e84043d77. This construction
keeps an explicit proper starting law and a genuine joint Gaussian likelihood.
Its remaining premise is the sharing/correlation law specified below.

INSTRUMENT DESCRIPTORS (NO DATA)
--------------------------------
For reading i with full SPD 2d-by-2d covariance Sigma_i, W=Sigma_i^-1, set
 c_i^2=tr(Sigma_FF)/tr(Sigma_aa), B_c=[c_i I;I], E_c=[c_i I;0],
 K_c=(B_c' W B_c)^-1,
 V_c=K_c^(1/2) E_c' [W-W B_c K_c B_c' W] E_c K_c^(1/2),
 D_i=V_c(I+V_c)^-1,
 A_i=log(c_i*K_c/k_ref), ell_i=log(c_i).
Here k_ref is any common unit of force-times-acceleration variance. Its
numerical value cancels from the law, because only differences of the log
eigenvalues appear. D_i is the bounded Grassmann velocity metric of the
candidate subspace at the instrument-only reference c_i. No eigenvector is
selected as physical: all matrix functions use their spectral projectors.

A PROPER SHARED GAUSSIAN MATRIX FIELD
-------------------------------------
For each real lambda put H_i(lambda)=A_i+lambda*D_i, with spectral projectors
P_ia(lambda) and eigenvalues h_ia(lambda). Define jointly centered symmetric
Gaussian matrices Z_i by the covariance kernel

 Cov(tr(U Z_i),tr(V Z_j))
 = integral phi(lambda) sum_ab
     exp(-[(h_ia-h_jb)^2+(ell_i-ell_j)^2]/2)
     tr(U P_ia) tr(V P_jb) d lambda.

This is positive semidefinite: at each lambda it is the spectral functional
calculus pushforward of one shared stationary scalar Gaussian process on R^2,
and the integral is a positive mixture. It defines consistent finite-dimensional
Gaussian laws for any number of readings. Repeated identical instruments use
the SAME latent matrix, not a new draw. Own coordinate frames act separately
on each matrix by conjugation; no alignment of arbitrarily chosen eigenvectors
between readings is used. Degenerate eigenvalues pose no ambiguity because
their equal kernel entries combine their projectors.

Let Phi be the standard-normal CDF, interpreted as a matrix function, and set
 Beta_i=sin^2[(pi/2) Phi(Z_i)], R_i=Beta_i(I-Beta_i)^-1.
For every finite Z, 0<Beta<I and R>0. The law is proper. Z and -Z have the
same law, so the JOINT law is invariant under all R_i -> R_i^-1, equivalently
all Beta_i -> I-Beta_i. For an isotropic noise pair Z=z I with z~N(0,1), so
Phi(z) is uniform and the one scalar Beta is EXACTLY the ours28 arcsine law.

EVIDENCE, CONDITIONAL CANDIDATES, AND TRANSPORT IN MASS
------------------------------------------------------
For each tested m set B_m=[m I;I], K_m=(B_m' W B_m)^-1,
 b_m=K_m B_m' W y, a_bar=K_m^-1/2 b_m, T=y' W y.
The candidate coordinate is v=K_m^-1/2 a, using the symmetric positive square
root in the original physical component frame. Declare v|R_i~N(0,R_i).
Then the exact joint law is

 p(m,dR|y) proportional pi(m) P_Sigma(dR)
   product_i |I-Beta_i|^(1/2) exp[-T_i/2+a_bar_i' Beta_i a_bar_i/2].

The conditional physical candidate acceleration is
 a_i|m,Beta_i,y ~ N(K_m^(1/2) Beta_i a_bar_i,
                    K_m^(1/2) Beta_i K_m^(1/2)), F_i=m a_i.
This is the explicit frame transport; a tensor has NOT silently been kept
fixed in a moving arbitrary whitened frame. The ratio-of-stacked-means readout
can be applied directly to this proper joint law.

The likelihood factor is between 0 and 1 after the displayed Gaussian data
constant is retained, since 0<Beta<I and ||a_bar||^2<=T. Hence any proper pi
gives a positive finite normalizer for every finite datum. Conditional means
and covariances are bounded by their flat-candidate counterparts in the needed
quadratic norm. For a finite fixed instrument set and a proper positive pi,
the usual endpoint bounds K_m=O(1) at m=0, K_m=O(m^-2) at infinity give finite
positive numerator and denominator for the ratio-of-means point.

EXACT GEOMETRIC PROPERTIES
---------------------------
* Independent changes of each reading's frame: all its descriptors, Z_i,
  Beta_i and physical candidates conjugate by that reading's orthogonal O_i.
  The joint field kernel is covariant under every collection {O_i}.
* Force units multiplied by f, acceleration units by a: c_i multiplies f/a,
  K_ci by a^2, so A_i shifts by log(f*a) I and ell_i by log(f/a), both COMMON
  shifts. The stationary field law is unchanged, and mass rescales f/a.
* Channel swap: c_i ->1/c_i, K_ci ->c_i^2 K_ci, so c_i K_ci and D_i are
  unchanged while ell_i changes sign. The kernel is invariant. At arbitrary
  mass, K'_(1/m)=m^2 K_m and b'_(1/m)=m b_m, so v and the candidate tensor
  retain precisely the same physical-coordinate interpretation.
* Coincident eigenvalues: matrix functional calculus and its Gaussian covariance
  are continuous; there is no discrete choice of an eigenbasis or noise group.
* Independent precise-null axes: append a covariance block epsilon^2 Sigma_0
  and zero observations. Its log-product eigenvalues go to -infinity, while D
  stays bounded. Its cross-field covariance with all finite-precision blocks
  tends to zero for every fixed lambda, hence also after lambda integration.
  Gaussian independence follows. The active descriptor block tends to its
  original value. The empty block's evidence becomes an m-independent factor,
  so the normalized mass law reduces to the lower-dimensional law. Its physical
  candidate covariance is O(epsilon^2), so its readout contribution vanishes.
  This statement concerns independent appended blocks; persistent standardized
  error correlations can carry information and require the corresponding
  conditional-covariance experiment instead of simple deletion.

LEARNING RATHER THAN RESETTING
------------------------------
For finitely many distinct diagonal instrument types, the induced Gaussian
covariance of their scalar Z values is strictly positive definite (unless two
types have identical descriptors, when they share a scale). Thus their Beta
prior has a positive smooth density on (0,1)^G. Repeated readings learn these
group scales; the finite-dimensional prior correlation washes out asymptotically.
If group g has average squared standardized candidate size C_g, its learned
Beta is C_g/(1+C_g). With log-mass sensitivity w_g,
  variance(score)=sum M_g Beta_g^2(C_g+1)w_g^2
                =sum M_g Beta_g C_g w_g^2=expected curvature.
This recovers regular widths in the distinct diagonal-group case, including
the earlier scalar-scale counterexample, and allows deterministic bounded
push sequences with the indicated second-moment limits. It does not establish
arbitrary-full-covariance calibration or growing-number-of-types asymptotics.

REMAINING PREMISES / LIMITS
---------------------------
The Gaussian correlation kernel, its unit log-width, Gaussian lambda weighting,
reference c prescription, and sharing across readings are concrete added choices,
not consequences of the current axioms. The correlation kernel can be changed
while retaining the proved geometric properties. Sharing fixes the demonstrated
regular-width mechanism but means a finite-noise all-zero new reading can update
other readings' scales; unconditional literal-null neutrality is not claimed.
Splitting/repacking a reading must carry its original instrument block and the
same joint field variables; resetting their sharing changes the experiment.
The field therefore supplies a workable continuous extension and an explicit
remaining definition question, not uniqueness. The finite lambda quadrature
below defines proper approximating fields and preserves the listed symmetries.
These size/shape results presume the mass starting law pi itself transforms
under units/swap and has the required precise-null limit. They do not choose
or repair that separately specified mass marginal.
"""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss


def symfun(a,fun):
    val,vec=np.linalg.eigh((a+a.T)/2)
    return (vec*fun(val))@vec.T


def descriptor(sigma):
    d=len(sigma)//2;eye=np.eye(d);W=np.linalg.inv(sigma)
    c=float(np.sqrt(np.trace(sigma[:d,:d])/np.trace(sigma[d:,d:])))
    B=np.vstack((c*eye,eye));E=np.vstack((c*eye,np.zeros((d,d))))
    K=np.linalg.inv(B.T@W@B);root=symfun(K,np.sqrt)
    residual=W-W@B@K@B.T@W
    velocity=root@E.T@residual@E@root
    velocity=(velocity+velocity.T)/2
    D=symfun(velocity,lambda x:x/(1+x))
    return {'d':d,'c':c,'K':K,'V':velocity,'A':symfun(c*K,np.log),'D':D,'ell':math.log(c)}


def cross_covariance(desc1,desc2,order=17):
    nodes,w=hermgauss(order);nodes=nodes*np.sqrt(2);w=w/np.sqrt(np.pi)
    out=np.zeros((desc1['d']**2,desc2['d']**2))
    for t,weight in zip(nodes,w):
        a,u=np.linalg.eigh(desc1['A']+t*desc1['D'])
        b,v=np.linalg.eigh(desc2['A']+t*desc2['D'])
        pu=np.stack([np.outer(u[:,i],u[:,i]).reshape(-1) for i in range(len(a))],axis=1)
        pv=np.stack([np.outer(v[:,i],v[:,i]).reshape(-1) for i in range(len(b))],axis=1)
        kernel=np.exp(-.5*((a[:,None]-b[None,:])**2+(desc1['ell']-desc2['ell'])**2))
        out+=weight*(pu@kernel@pv.T)
    return out


def diagonal_covariance(sf,sa):
    return np.diag(np.r_[np.array(sf)**2,np.array(sa)**2])


def invariant_checks():
    rng=np.random.default_rng(280940);rows=[]
    for d in (2,3):
        b=rng.normal(size=(2*d,2*d));sigma=b@b.T+np.eye(2*d)*.5
        base=descriptor(sigma);gamma=cross_covariance(base,base)
        o,_=np.linalg.qr(rng.normal(size=(d,d)))
        big=np.zeros((2*d,2*d));big[:d,:d]=o;big[d:,d:]=o
        rotated=descriptor(big@sigma@big.T)
        tensor=np.kron(o,o)
        rot_error=float(np.max(abs(cross_covariance(rotated,rotated)-tensor@gamma@tensor.T)))
        units=np.diag(np.r_[np.full(d,3.7),np.full(d,.4)])
        scaled=descriptor(units@sigma@units)
        unit_error=float(np.max(abs(cross_covariance(scaled,scaled)-gamma)))
        idx=np.r_[np.arange(d,2*d),np.arange(d)]
        swapped=descriptor(sigma[np.ix_(idx,idx)])
        swap_error=float(np.max(abs(cross_covariance(swapped,swapped)-gamma)))
        assert max(rot_error,unit_error,swap_error)<1e-10
        ev=np.linalg.eigvalsh(gamma)
        rows.append({'d':d,'rotation_error':rot_error,'units_error':unit_error,
                     'swap_error':swap_error,'positive_covariance_eigenvalues':ev[ev>1e-11].tolist(),
                     'symmetric_matrix_dimension':d*(d+1)//2})
    # Cross-reading covariance also has independent-frame and unit/swap covariance.
    sigma2=diagonal_covariance([2.,.3,1.1],[.7,1.3,.5])
    d1=descriptor(sigma);d2=descriptor(sigma2);cross=cross_covariance(d1,d2)
    o2,_=np.linalg.qr(rng.normal(size=(3,3)));big2=np.zeros((6,6));big2[:3,:3]=o2;big2[3:,3:]=o2
    rot2=descriptor(big2@sigma2@big2.T)
    error=float(np.max(abs(cross_covariance(rotated,rot2)-np.kron(o,o)@cross@np.kron(o2,o2).T)))
    assert error<1e-10
    return rows,error


def degeneracy_checks():
    base=descriptor(diagonal_covariance([1.,1.],[1.,1.]))
    target=np.outer(np.eye(2).reshape(-1),np.eye(2).reshape(-1))
    assert np.max(abs(cross_covariance(base,base)-target))<1e-12
    rows=[]
    for delta in (.1,.01,.001,.0001):
        desc=descriptor(diagonal_covariance([1+delta,1-delta],[1.,1.]))
        rows.append({'delta':delta,'covariance_distance_to_isotropic_limit':float(np.linalg.norm(cross_covariance(desc,desc)-target))})
    return rows


def phi_normal(x):
    return .5*(1+np.vectorize(math.erf)(x/np.sqrt(2)))


def known_empty_law(order=64):
    nodes,w=hermgauss(order);nodes=nodes*np.sqrt(2);w=w/np.sqrt(np.pi)
    z1=nodes[:,None];z2=nodes[None,:];ww=w[:,None]*w[None,:]
    theta,wt=leggauss(240);theta=(theta+1)*np.pi/4;wt=wt*np.pi/4
    H=4*np.sin(theta)**2
    beta=np.sin(np.pi*phi_normal(nodes)/2)**2
    base=np.array([np.dot(w,np.sqrt(1-beta)*np.exp(.5*beta*h)) for h in H])
    base=base/np.dot(wt,base)
    rows=[]
    for eps in (1.,.5,.3,.1,.03,.01):
        # D is the same scalar in both axes, so the lambda integral cancels here.
        rho=math.exp(-.5*math.log(eps*eps)**2)
        b1=np.sin(np.pi*phi_normal(z1)/2)**2
        b2=np.sin(np.pi*phi_normal(rho*z1+math.sqrt(max(0,1-rho*rho))*z2)/2)**2
        p=ww*np.sqrt((1-b1)*(1-b2))
        curve=np.array([np.sum(p*np.exp(.5*b1*h)) for h in H])
        curve=curve/np.dot(wt,curve)
        rows.append({'empty_axis_sd':eps,'latent_scale_correlation':rho,
                     'mass_law_TV_from_1D':float(.5*np.dot(wt,abs(curve-base)))})
    assert rows[-1]['mass_law_TV_from_1D']<1e-12
    return rows


def diagonal_learning_check():
    # Declared as two independent instrument types; repeated readings share fields.
    first=descriptor(diagonal_covariance([1.],[1.]))
    second=descriptor(diagonal_covariance([10.],[1.]))
    rho=float(cross_covariance(first,second)[0,0])
    C=np.array([.1,9.9]);w=np.array([.5,10/101]);beta=C/(1+C)
    variance=float(np.sum(beta*beta*(1+C)*w*w))
    curvature=float(np.sum(beta*C*w*w))
    assert abs(variance-curvature)<1e-15
    assert 0<rho<.01
    return {'prior_scale_Gaussian_correlation':rho,
            'learned_group_beta':beta.tolist(),
            'sampling_score_variance_per_pair':variance,
            'expected_curvature_per_pair':curvature,
            'asymptotic_variance_ratio':variance/curvature}


def main():
    invariant,cross=invariant_checks()
    rng=np.random.default_rng(280941);q=rng.normal(size=(6,6))
    desc=descriptor(q@q.T+np.eye(6)*.5)
    lambda_error=float(np.max(abs(cross_covariance(desc,desc,17)-cross_covariance(desc,desc,33))))
    null_rows=known_empty_law(64);refined_null=known_empty_law(96)
    null_error=max(abs(x['mass_law_TV_from_1D']-y['mass_law_TV_from_1D']) for x,y in zip(null_rows,refined_null))
    out={'source_commit':'592d144531fa4a30d1c2229a4c3b528e84043d77',
         'invariance_checks':invariant,'independent_reading_frame_cross_covariance_error':cross,
         'degeneracy_continuity':degeneracy_checks(),
         'independent_precise_null_reduction':null_rows,
         'lambda_quadrature_17_to_33_max_covariance_change':lambda_error,
         'normal_quadrature_64_to_96_max_law_TV_change':null_error,
         'diagonal_repeated_group_width_identity':diagonal_learning_check(),
         'status':'Constructive proper continuous tensor-field candidate; sharing/kernel prescriptions remain additional premises; general full-covariance calibration unproved.'}
    Path(__file__).with_name('instrument_covariance_field_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
