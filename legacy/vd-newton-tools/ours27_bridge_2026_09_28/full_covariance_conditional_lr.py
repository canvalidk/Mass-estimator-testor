"""Global conditional noise test with arbitrary declared SPD covariance.

Not a posterior, candidate probability law, or next adopted mass definition.
Stack y=(F,a), whiten z=Sigma^(-1/2)y. P(theta) is the orthogonal projector
onto whitened B_theta=(sin(theta)I;cos(theta)I), theta in [0,pi/2]. Endpoints
represent limits of positive masses; no endpoint point estimate is output.

At tested theta0, t=U0'z and r=W0'z are independent, r~N(0,I) for EVERY fixed
aligned candidate. Holding t fixed and replacing r by independent Gaussian
noise therefore gives the exact conditional experiment without a size prior.

Lambda_0(z)=sup_theta z'P(theta)z-z'P(theta0)z is the global profile LR
statistic. Its conditional tail probability is a valid test, including at
nonquadratic peaks and multiple modes. It conditions on the model and uses
relative profile fit; it is not an absolute model-rejection gate.

This implementation uses a finite, predeclared grid including the tested
theta. The resulting GRID statistic remains an exact conditional test in
principle; Monte Carlo approximates its critical values. Grid density affects
power, not the conditional validity argument. For the quartic demonstration,
additional local grid points depend only on the conditioned fitted component,
so validity is retained. Monte Carlo comparisons use independent calibration
and assessment samples. Atoms at LR=0 are accepted, making the test conservative.

LIMIT: this global conditional test does not preserve known-empty-axis
continuity. In exact/coarea_known_empty.py and exact/conditional_lr_known_empty.py
the distinct geometry/calibration failures are recorded. For the LR test,
F=2,a=0,noise=1,tested m=1 gives p=.1572992071; adding a measured-zero axis
with both noise SDs epsilon gives p=.2391605842 for every epsilon>0. The
calibration experiment changes even as physical extra-axis noise tends to zero.
Do not infer a complete estimator satisfying the project's representation
requirements from the conditional-validity proof.

Source context: VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.
"""
from pathlib import Path
import json
import numpy as np


def whitening(sigma):
    e,u=np.linalg.eigh(sigma)
    return (u/np.sqrt(e))@u.T


def frame(white,theta):
    M=len(white)//2
    B=np.vstack([np.sin(theta)*np.eye(M),np.cos(theta)*np.eye(M)])
    qr,_=np.linalg.qr(white@B,mode='complete')
    U=qr[:,:M];W=qr[:,M:]
    return U,W,U@U.T


def grid_projectors(white,theta0,t=None):
    angles=np.r_[np.linspace(0,np.pi/2,301),theta0]
    if t is not None:
        local=theta0+np.linspace(-8.,8.,161)/np.sqrt(1+np.sum(t*t))
        angles=np.r_[angles,local[(local>0)&(local<np.pi/2)]]
    angles=np.unique(angles)
    return np.array([frame(white,a)[2] for a in angles])


def statistic(z,projectors,P0):
    z=np.atleast_2d(z)
    out=np.empty(len(z))
    flat=projectors.reshape(len(projectors),-1)
    for start in range(0,len(z),1024):
        zz=z[start:start+1024]
        products=(zz[:,:,None]*zz[:,None,:]).reshape(len(zz),-1)
        captured=products@flat.T
        fitted=products@P0.ravel()
        out[start:start+len(zz)]=np.maximum(captured.max(axis=1)-fitted,0.)
    return out


def conditional_draws(rng,U,W,t,n):
    return U@t+rng.normal(size=(n,len(t)))@W.T


def main():
    rng=np.random.default_rng(28092832)
    raw=np.array([[1.,.4,-.2,.1],[.2,1.3,.3,-.1],[.1,-.2,.8,.35],[.2,.1,.1,1.1]])
    sigma=raw@raw.T+.2*np.eye(4)
    white=whitening(sigma);theta0=.637
    U,W,P0=frame(white,theta0)
    projectors=grid_projectors(white,theta0)
    cells=[]
    for t in [np.zeros(2),np.array([1.,-.5]),np.array([8.,-3.])]:
        ref=statistic(conditional_draws(rng,U,W,t,20000),projectors,P0)
        critical=float(np.quantile(ref,.95,method='higher'))
        obs=statistic(conditional_draws(rng,U,W,t,20000),projectors,P0)
        coverage=float(np.mean(obs<=critical))
        assert .94<coverage<.96,(t,coverage)
        cells.append({'conditioned_component':t.tolist(),'reference_draws':len(ref),
            'independent_assessment_draws':len(obs),'conditional_95_critical_value':critical,
            'assessment_coverage':coverage,'grid_alternatives':len(projectors)})

    # The whole-profile operation retains the quartic peak that j/Vc erases.
    s=np.array([np.sqrt(2)+1,np.sqrt(2)-1])
    sigmaq=np.diag(np.r_[s*s,np.ones(2)]);whiteq=whitening(sigmaq)
    quartic=[]
    for k in [10.,100.]:
        physical=np.r_[s*k,np.full(2,k)];z=whiteq@physical
        for u0 in [0., .75/np.sqrt(k),1.5/np.sqrt(k),3./np.sqrt(k)]:
            tested=np.arctan(np.exp(u0))
            U,W,P0=frame(whiteq,tested);t=U.T@z
            ps=grid_projectors(whiteq,tested,t)
            obs=float(statistic(z,ps,P0)[0])
            simulated=statistic(conditional_draws(rng,U,W,t,12000),ps,P0)
            tail=float((1+np.sum(simulated>=obs))/(1+len(simulated)))
            quartic.append({'k':k,'tested_log_mass':float(u0),'profile_LR':obs,
                'conditional_tail_MC':tail,'conditional_reference_draws':len(simulated)})
    out={'status':'Global finite-grid conditional noise test; not a mass probability law.',
        'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8',
        'conditional_calibration_cells':cells,'quartic_profile_test':quartic,
        'warning':'Do not normalize these p-values as a posterior or infer candidate-size weights from them.'}
    Path(__file__).with_name('full_covariance_conditional_lr_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
