"""Full S^2 check at the first weak-z baseline sample (strong overall vector)."""
import json
import time

import numpy as np

from estimator import estimate, profile
from robot_cart_experiments import OUT, R, G, load_robot, fit_calibration


def profile_mass(f, a, cov):
    # Same-covariance conventional errors-in-variables profile fit. A coarse
    # logarithmic scan brackets the minimum before golden-section refinement.
    grid=np.linspace(-8,8,401)
    q=np.array([profile(f,a,cov,np.exp(x)) for x in grid])
    i=int(np.argmin(q)); assert 0<i<len(grid)-1
    lo,hi=grid[i-1],grid[i+1]
    for _ in range(90):
        x1=hi-(hi-lo)/((1+np.sqrt(5))/2)
        x2=lo+(hi-lo)/((1+np.sqrt(5))/2)
        if profile(f,a,cov,np.exp(x1)) < profile(f,a,cov,np.exp(x2)): hi=x2
        else: lo=x1
    return float(np.exp((lo+hi)/2))


def run():
    result = json.loads((OUT/'robot_cart.json').read_text())
    m, fb, ab, _, _ = fit_calibration(load_robot('0-calibration_fts-accel'))
    aa, ff = load_robot('1-baseline_accel'), load_robot('1-baseline_wrench')
    row = next(e for e in result['robot']['examples'] if e['label']=='1-baseline_first_weak')
    idx = row['imu_csv_data_row']-1
    a = -(aa[idx,1:4]-ab) @ R.T*G
    f = np.array([np.interp(aa[idx,0]-8416, ff[:,0], ff[:,j+1]) for j in range(3)])-fb
    sf = np.array(result['robot']['total_force_sigma_N'])
    sa = np.array(result['robot']['total_acceleration_sigma_m_s2'])
    cov = np.diag(np.r_[sf**2, sa**2])
    # Rotate coordinates to put the narrow angular peak near a quadrature pole.
    # Rotate BOTH vectors and the full covariance; the law and flat solid-angle
    # measure are unchanged. No direction is declared known or fixed.
    z = a/np.linalg.norm(a)
    x = np.cross(z, np.eye(3)[np.argmin(abs(z))]); x /= np.linalg.norm(x)
    y = np.cross(z,x)
    rot = np.array([x,y,z]); transform = np.zeros((6,6))
    transform[:3,:3]=rot; transform[3:,3:]=rot
    assert np.allclose(rot@rot.T,np.eye(3))
    previous = None; checks=[]
    for n in [32,64,128,256]:
        start = time.monotonic()
        e = estimate(rot@f, rot@a, transform@cov@transform.T,
                     direction_order=n, ratio_order=1024)
        values=np.array([e.mass,e.quantile(.025),e.quantile(.975)])
        checks.append(dict(direction_order=n, ratio_order=1024, mass_kg=e.mass,
                           quantiles_kg=values[1:].tolist(), seconds=time.monotonic()-start))
        print(checks[-1], flush=True)
        if previous is not None and np.max(abs(values/previous-1)) < 2e-4:
            break
        previous=values
    else:
        raise AssertionError('Angular quadrature not converged')
    # Independently refine the ratio-angle axis at the converged angular order.
    fine=estimate(rot@f,rot@a,transform@cov@transform.T,
                  direction_order=n,ratio_order=2048)
    fine_values=np.array([fine.mass,fine.quantile(.025),fine.quantile(.975)])
    assert np.max(abs(fine_values/values-1)) < 2e-4
    out=dict(force_vector_N=f.tolist(), acceleration_force_sign_m_s2=a.tolist(),
             covariance=cov.tolist(), component_example=row,
             ordinary_vector_OLS_mass_kg=float(f@a/(a@a)),
             conventional_same_covariance_profile_mass_kg=profile_mass(f,a,cov),
             extended_mass_kg=fine.mass, law_q025_kg=fine_values[1], law_q975_kg=fine_values[2],
             q_at_calibration_mass=profile(f,a,cov,m), angular_convergence=checks,
             ratio_refinement_relative_change=(fine_values/values-1).tolist(),
             caveat='Specific-acceleration/support-force convention; dynamic model and covariance provisional.')
    (OUT/'robot_full_vector.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__': run()
