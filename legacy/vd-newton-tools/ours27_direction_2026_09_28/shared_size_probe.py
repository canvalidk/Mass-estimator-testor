"""Audit the unadopted reference-size law; no fitted hyperparameters.

An all-null mass marginal can be unchanged while its shared size state learns.
The prior dt/t on t=1+omega^2>=1 is an added reference-model premise.
These checks are not a claim of uniqueness or universal interval coverage.
"""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss


def gl(n):
    x, w = leggauss(n)
    return (x+1)/2, w/2


def folded_mean(h):
    h = np.asarray(h)
    erf = np.fromiter((math.erf(float(x)/math.sqrt(2)) for x in h.flat), float).reshape(h.shape)
    return np.sqrt(2/np.pi)*np.exp(-h*h/2)+h*erf


def radial_mean(h, M):
    if M == 1:
        return folded_mean(h)
    assert M == 3
    h = np.asarray(h)
    erf = np.fromiter((math.erf(float(x)/math.sqrt(2)) for x in h.flat), float).reshape(h.shape)
    return np.sqrt(2/np.pi)*np.exp(-h*h/2)+(h+1/h)*erf


def law(P, Q, D, M, ntheta=8001, nscale=200, point_dimension=None):
    theta = np.linspace(0, np.pi/2, ntheta)
    H = P*np.sin(theta)**2+Q*np.cos(theta)**2+D*np.sin(2*theta)
    u, w = gl(nscale)
    beta = 1-u*u
    # u=1/sqrt(t): dt/t times t^(-M/2) becomes 2 u^(M-1) du.
    logweights = np.log(2*w)+(M-1)*np.log(u)
    jointlog = logweights[:,None]+.5*beta[:,None]*H[None,:]
    shift = np.max(jointlog)
    joint = np.exp(jointlog-shift)
    marginal = np.sum(joint, axis=0)
    density = marginal/np.trapezoid(marginal, theta)
    cdf = np.r_[0., np.cumsum((density[:-1]+density[1:])*np.diff(theta)/2)]
    qs = np.tan(np.interp([.025,.5,.975], cdf, theta))
    meanbeta = np.sum(beta[:,None]*joint, axis=0)/marginal
    out = {'M': M, 'interval95_and_median': qs.tolist(),
           'beta_at_peak': float(meanbeta[np.argmax(H)]),
           'mass_cdf_at_one': float(np.interp(np.pi/4,theta,cdf))}
    if point_dimension is not None:
        # Aligned observations, all physical size in the first component.
        # point_dimension=1 with M=3 is the physical epsilon-to-zero limit:
        # transverse physical candidates vanish while scale inference retains M.
        h = np.sqrt(beta[:,None]*H[None,:])
        size = np.sqrt(beta[:,None])*radial_mean(h, point_dimension)
        size_marginal = np.sum(joint*size, axis=0)
        out['ratio_of_mean_sizes'] = float(np.trapezoid(np.sin(theta)*size_marginal,theta)/np.trapezoid(np.cos(theta)*size_marginal,theta))
    return out, theta, density


def mixture_factors(H, M, n=400):
    u,w = gl(n); beta=1-u*u
    z=2*w*u**(M-1)*np.exp(.5*beta*H)
    return float(z.sum()), float(z@beta/z.sum())


def main():
    rows=[]
    for M,point_dimension in [(1,1),(3,3),(3,1),(7,None),(63,None)]:
        row,theta,density=law(9.,1.,3.,M,point_dimension=point_dimension)
        fine,_,density_fine=law(9.,1.,3.,M,nscale=400,point_dimension=point_dimension)
        err=float(np.max(abs(density-density_fine)))
        assert err<2e-9,err
        row['physical_point_dimension']=point_dimension
        row['refinement_max_density_error']=err
        rows.append(row)

    # Exact sequential update: n null axes alone leave angular marginal flat,
    # but change t law by t^(-n/2), thus changing the next nonzero evidence.
    H=7.; initial_M=1; null_axes=2
    u,w=gl(400); beta=1-u*u
    joint=2*w*u**(initial_M-1)*np.exp(.5*beta*H)
    joint*=u**null_axes
    one_shot, one_beta=mixture_factors(H,initial_M+null_axes)
    assert abs(joint.sum()/one_shot-1)<1e-13
    assert abs(joint@beta/joint.sum()-one_beta)<1e-13
    nulls=[]
    for M in [1,3,63]:
        row,_,density=law(0.,0.,0.,M)
        assert np.max(abs(density-2/np.pi))<1e-12
        assert abs(row['beta_at_peak']-2/(M+2))<1e-12
        nulls.append({'M':M,'mean_beta':row['beta_at_peak'],'mass_law':'unchanged uniform angle'})

    # Mixed signal sizes: posterior curvature factor approaches the regular
    # lambda/(1+lambda) benchmark without selecting a scale from scores.
    regular=[]
    for M in [63,400,1600]:
        for lam in [1/3,1.,3.]:
            H=M*(1+lam)
            u,w=gl(500); beta=1-u*u
            lw=np.log(2*w)+(M-1)*np.log(u)+.5*beta*H
            z=np.exp(lw-lw.max())
            regular.append({'M':M,'lambda':lam,'posterior_mean_beta':float(z@beta/z.sum()),'regular_limit':lam/(1+lam)})

    out={'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8',
         'status':'Reopened reference-size candidate, not an adopted equation.',
         'aligned_force':3.,'aligned_acceleration':1.,'noise_SD_on_active_axis':1.,
         'aligned_null_padding':rows,'all_null_experiments':nulls,
         'sequential_shared_state_check':'passed; resetting to mass marginal fails',
         'regular_curvature_checks':regular}
    Path(__file__).with_name('shared_size_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
