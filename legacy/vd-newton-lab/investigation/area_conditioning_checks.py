"""Single-pair flat measure as vanishing positive wedge-area conditioning.

Finite-band quadrature is for 3D independent unit Gaussian channels with
collinear observed centers. It integrates the original six-dimensional ambient
Gaussian event using a radial/directional reduction, not the limiting surface
law. The reported finite-band norm ratio is only a convergence diagnostic.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from gaussian_conditioning import SOURCE_COMMIT, summary


@lru_cache(None)
def unit_nodes(order):
    z,w=np.polynomial.legendre.leggauss(order)
    return (z+1)/2,w/2


def normal_cdf(x):
    x=np.asarray(x)
    return np.array([.5*math.erfc(-float(v)/math.sqrt(2)) for v in x.flat]).reshape(x.shape)


def disk_probability_over_radius_squared(radius,offset,order):
    """P{|N(mu,I2)|<b}/b², where |mu|=offset and b=radius.

    If b>=|mu|+12, replace its probability by one. The missing probability is
    bounded by P{|N(0,I2)|>12}=exp(-72). Otherwise polar disk quadrature is stable
    for the observed-center magnitudes in this study.
    """
    radius,offset=np.broadcast_arrays(radius,offset)
    result=np.empty_like(radius)
    large=radius>=offset+12
    result[large]=1/radius[large]**2
    b,h=radius[~large],offset[~large]
    q,w=unit_nodes(order)
    z=b[:,None]*q
    argument=h[:,None]*z
    if argument.size and np.max(argument)>500:
        raise ValueError('Observed center outside this disk-integrator validation range')
    integrand=q*np.exp(-.5*(z*z+h[:,None]**2))*np.i0(argument)
    result[~large]=integrand@w
    return result


def finite_band_one_side(force,acceleration,epsilon,order=64,padding=12.):
    """Return P(event)/epsilon² and E[|X| | event].

    X~N(force*e3,I3), Y~N(acceleration*e3,I3).
    Event: X dot Y>0 and |X cross Y|<epsilon.
    Given X=x*u, Y_parallel~N(acceleration*u3,1) and Y_perpendicular is an
    independent 2D Gaussian. Its allowed disk radius is epsilon/x.
    """
    if epsilon<=0: raise ValueError('Positive tolerance required')
    upper=abs(force)+padding
    cuts=sorted(set([0.,upper]+[min(upper,v) for v in
        (epsilon/(abs(acceleration)+12),epsilon/(abs(acceleration)+1),
         epsilon,math.sqrt(epsilon))]))
    q,w=unit_nodes(order)
    x=np.concatenate([lo+(hi-lo)*q for lo,hi in zip(cuts[:-1],cuts[1:])])
    wx=np.concatenate([(hi-lo)*w for lo,hi in zip(cuts[:-1],cuts[1:])])
    cosine=2*q-1
    wt=2*w
    offset=abs(acceleration)*np.sqrt(1-cosine*cosine)
    radius=epsilon/x[:,None]
    disk=disk_probability_over_radius_squared(radius,offset[None,:],order)
    gaussian=np.exp(-.5*(x[:,None]**2+force*force-2*x[:,None]*force*cosine))
    positive=normal_cdf(acceleration*cosine)
    weighted=wx[:,None]*wt*gaussian*positive*disk/math.sqrt(2*math.pi)
    scaled_probability=float(weighted.sum())
    mean_x=float((weighted*x[:,None]).sum()/scaled_probability)
    return {'probability_over_epsilon_squared':scaled_probability,'mean_x':mean_x}


def finite_band(force,acceleration,epsilon,order=64,padding=12.):
    forward=finite_band_one_side(force,acceleration,epsilon,order,padding)
    reverse=finite_band_one_side(acceleration,force,epsilon,order,padding)
    discrepancy=abs(forward['probability_over_epsilon_squared']/reverse['probability_over_epsilon_squared']-1)
    return {**forward,'mean_y':reverse['mean_x'],
            'finite_band_norm_ratio':forward['mean_x']/reverse['mean_x'],
            'channel_swap_normalization_error':discrepancy}


def noise_checks(draws=200000,seed=20260918):
    """Verify area² moments and the exact null Gamma(d-1,1) law."""
    rng=np.random.default_rng(seed)
    rows=[]
    for d in (2,3):
        for strength in (0.,1.,4.):
            x,y=rng.normal(size=(2,draws,d))
            x[:,0]+=.6*strength
            y[:,0]+=.8*strength
            xx,yy,xy=(x*x).sum(-1),(y*y).sum(-1),(x*y).sum(-1)
            area_squared=np.maximum(0.,xx*yy-xy*xy)
            expected=(d-1)*(strength*strength+d)
            eu2=(d+strength*strength)**2+2*(d+2*strength*strength)
            variance=eu2*(d-1)*(d+1)-expected**2
            mcse=math.sqrt(variance/draws)
            observed=float(area_squared.mean())
            assert abs(observed-expected)<6*mcse
            row={'dimension':d,'combined_true_strength':strength,
                 'area_squared_mean':observed,'exact_mean':expected,
                 'monte_carlo_standard_error':mcse}
            if strength==0:
                probabilities=[]
                for epsilon in (.5,1.,2.):
                    gamma_cdf=-math.expm1(-epsilon)
                    if d==3: gamma_cdf-=epsilon*math.exp(-epsilon)
                    target=gamma_cdf/2
                    actual=float(np.mean((xy>0)&(area_squared<epsilon**2)))
                    assert abs(actual-target)<6*math.sqrt(target*(1-target)/draws)
                    probabilities.append({'epsilon':epsilon,'probability':actual,'exact':target})
                row['positive_area_event']=probabilities
            rows.append(row)
    return {'draws_per_case':draws,'seed':seed,'cases':rows}


def run():
    results={}
    worst_swap=0.
    worst_refinement=0.
    for name,force,acceleration in (
        ('zero_zero',0.,0.),('force_only',3.,0.),
        ('aligned',4.,1.),('opposed',4.,-1.)):
        target=summary([0.,0.,force],[0.,0.,acceleration],'flat',256)
        rows=[]
        for epsilon in (1.,.3,.1,.03,.01,.003):
            low=finite_band(force,acceleration,epsilon)
            high=finite_band(force,acceleration,epsilon,96)
            for key in ('probability_over_epsilon_squared','mean_x','mean_y','finite_band_norm_ratio'):
                error=abs(low[key]/high[key]-1)
                worst_refinement=max(worst_refinement,error)
                assert error<3e-9,(name,epsilon,key,error)
            worst_swap=max(worst_swap,high['channel_swap_normalization_error'])
            assert high['channel_swap_normalization_error']<3e-9
            if name=='zero_zero':
                exact=(-math.expm1(-epsilon)-epsilon*math.exp(-epsilon))/(2*epsilon**2)
                assert abs(high['probability_over_epsilon_squared']/exact-1)<3e-9
            rows.append({'epsilon':epsilon,**high})
        extended=finite_band(force,acceleration,.003,96,16.)
        for key in ('probability_over_epsilon_squared','mean_x','mean_y'):
            error=abs(extended[key]/rows[-1][key]-1)
            worst_refinement=max(worst_refinement,error)
            assert error<3e-9
        limit={'probability_over_epsilon_squared':target['normalization']/(2*math.pi),
               'mean_x':target['mean_x'],'mean_y':target['mean_y'],
               'mass_over_s':target['mass_over_s']}
        assert abs(rows[-1]['finite_band_norm_ratio']/limit['mass_over_s']-1)<.01
        results[name]={'limit':limit,'finite_tolerances':rows}
        print(json.dumps({'case':name,'flat_limit':limit,
                          'smallest_tolerance':rows[-1]}),flush=True)
    return {'finite_band_cases':results,'noise_sampling':noise_checks(),
            'worst_relative_quadrature_refinement':worst_refinement,
            'worst_relative_channel_swap_normalization_error':worst_swap,
            'perpendicular_tail_probability_bound':math.exp(-72)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=run()
    result.update({'source_commit':SOURCE_COMMIT,
                   'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'limiting_integrator_sha256':hashlib.sha256(Path(__file__).with_name('gaussian_conditioning.py').read_bytes()).hexdigest()})
    print('All assertions passed.',flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
