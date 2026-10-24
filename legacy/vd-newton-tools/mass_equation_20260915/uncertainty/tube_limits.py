"""Independent d=3 finite-tube check in original magnitude coordinates.

Two independent N(0,I_3) vectors have iid chi_3 norms and uniform relative
direction cosine. Its exact integration leaves the two magnitude integrals.
Only scratch verification output is written.
"""
from functools import lru_cache
from pathlib import Path
import json
import math
import numpy as np

C = math.sqrt(2/math.pi)
MAXIMUM = 14.


@lru_cache(None)
def nodes(order):
    return np.polynomial.legendre.leggauss(order)


def grid(edges, order):
    z, w = nodes(order)
    edges = sorted(set(float(x) for x in edges if 0 <= x <= MAXIMUM))
    points, weights = [], []
    for left,right in zip(edges, edges[1:]):
        t = (z+1)/2
        if left > 0:
            # Removes the square-root cusp at the tube's norm threshold.
            points.append(left + (right-left)*t*t)
            weights.append((right-left)*w*t)
        else:
            points.append(left + (right-left)*t)
            weights.append((right-left)*w/2)
    return np.concatenate(points), np.concatenate(weights)


def angular_event_scaled(x, y, eps):
    """P(distance to positive collinearity cone <= eps | x,y)/eps^2.

    If either norm <= eps, every relative orientation lies in the tube.
    Otherwise the uniform cosine must exceed sqrt((1-eps²/x²)(1-eps²/y²)).
    Rationalization avoids cancellation for thin tubes.
    """
    ex, ey = (eps/x)**2, (eps/y)**2
    interior = (x>eps) & (y>eps)
    product = np.maximum(0., (1-ex)*(1-ey))
    resolved = .5*(1/x**2+1/y**2-eps**2/(x*y)**2)/(1+np.sqrt(product))
    return np.where(interior, resolved, 1/eps**2)


def direct_moments(eps, order=256):
    x,w = grid([0,eps,MAXIMUM],order)
    chi_weights = w*C*x*x*np.exp(-x*x/2)
    xx, yy = x[:,None], x[None,:]
    weight = chi_weights[:,None]*chi_weights[None,:]*angular_event_scaled(xx,yy,eps)
    z = float(weight.sum())
    log_ratio = np.log(xx)-np.log(yy)
    # Minimum-norm<=eps is precisely the part having some boundary projections;
    # exactly half of those orientations have negative dot product.
    boundary_projection = float((weight*((xx<=eps)|(yy<=eps))).sum()/(2*z))
    return dict(normalizer_over_epsilon_squared=z,
                point_ratio=float((weight*xx).sum()/(weight*yy).sum()),
                mean_norm=float((weight*xx).sum()/z),
                log_sd=math.sqrt(float((weight*log_ratio**2).sum()/z)),
                ratio_mean=float((weight*xx/yy).sum()/z),
                ratio_second_moment=float((weight*(xx/yy)**2).sum()/z),
                projected_boundary_probability=boundary_projection)


def cdf(eps, ratio, order=128):
    """Direct radial integral with x <= ratio*y, avoiding ratio-bin error."""
    y, wy = grid([0,eps,eps/ratio,MAXIMUM],order)
    pyw = wy*C*y*y*np.exp(-y*y/2)
    x_upper = np.minimum(MAXIMUM,ratio*y)
    z,w = nodes(order)
    numerator = np.zeros_like(y)
    for panel,(lower, upper) in enumerate(((np.zeros_like(y),np.minimum(eps,x_upper)),
                                         (np.minimum(eps,x_upper),x_upper))):
        t = (z[None,:]+1)/2
        x = lower[:,None]+(upper-lower)[:,None]*(t if panel==0 else t*t)
        xw = (upper-lower)[:,None]*w[None,:]*(.5 if panel==0 else t)
        # The second interval is zero-length where x_upper<eps. Its positive
        # x endpoint still lets the event formula remain numerically well-defined.
        values = C*x*x*np.exp(-x*x/2)*angular_event_scaled(x,y[:,None],eps)
        numerator += (xw*values).sum(axis=1)
    # Denominator uses the matching quadrature order to avoid mixing resolutions.
    return float(pyw @ numerator)/direct_moments(eps,order)['normalizer_over_epsilon_squared']


def lower_quantile(eps, order=128):
    lo, hi = 0., 1.
    for _ in range(32):
        mid=(lo+hi)/2
        if cdf(eps,mid,order)<.025:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def chi3_ratio_cdf(r):
    theta=math.atan(r)
    return 2*theta/math.pi-math.sin(4*theta)/(2*math.pi)


def main():
    rows=[]
    for eps in (.5,.2,.1,.05,.02,.01,.005,.001,.0005):
        row=direct_moments(eps)
        high=direct_moments(eps,512)
        refinement=max(abs(row[k]/high[k]-1) for k in row)
        assert refinement < 2e-5
        q=lower_quantile(eps)
        row.update(epsilon=eps,central_95=[q,1/q],maximum_moment_refinement=refinement)
        rows.append(row)
    lo,hi=0.,1.
    for _ in range(64):
        mid=(lo+hi)/2
        if chi3_ratio_cdf(mid)<.025:
            lo=mid
        else:
            hi=mid
    q=(lo+hi)/2
    result=dict(angular_limit=dict(point=1.,mean_norm=2*C,
                                   central_95=[q,1/q],log_sd=math.sqrt(math.pi**2/4-2),
                                   ratio_mean=4/math.pi,ratio_second_moment=3.,
                                   ratio_variance=3-16/math.pi**2),
                tube_limit=dict(point=1.,mean_norm=1.5*C,
                                central_95=[math.tan(math.pi*.025/2),math.tan(math.pi*.975/2)],
                                log_sd=math.pi/2,ratio_mean='divergent',
                                normalizer_over_epsilon_squared=.5),finite_tubes=rows)
    assert abs(rows[-1]['normalizer_over_epsilon_squared']/.5-1) < .01
    assert abs(rows[-1]['log_sd']/(math.pi/2)-1) < .02
    out=Path(__file__).with_name('tube_verification.json')
    out.write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
