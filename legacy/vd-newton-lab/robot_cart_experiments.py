"""Reproducible scalar restrictions of the compatible-pair mass estimator.

Real measurements, exploratory uncertainty budgets; no accuracy-validation claim.
Robot projections and cart longitudinal data test the S^0 integral, not the full
unknown-direction S^2 integral. Raw data are preserved alongside this script.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from estimator import estimate, isotropic_covariance, profile

HERE = Path(__file__).resolve().parent
RAW = HERE / 'raw'
OUT = HERE / 'results'
G = 9.82085
R = np.array([[0., -1., 0.], [0., 0., 1.], [-1., 0., 0.]])


def load_robot(name):
    return np.genfromtxt(RAW / 'robot' / (name + '.csv'), delimiter=',', skip_header=1)


def fit_calibration(c):
    f, a, g = c[:, :3], c[:, 6:9], c[:, 9:12]
    design = np.column_stack([g.ravel(), np.tile(np.eye(3), (len(c), 1))])
    assert np.linalg.matrix_rank(design) == 4
    fit = np.linalg.lstsq(design, f.ravel(), rcond=None)[0]
    ab = np.mean(a + g @ R / G, axis=0)
    return fit[0], fit[1:], ab, f-fit[1:]-fit[0]*g, -(a-ab) @ R.T*G-g


def readout(label, f, a, sf, sa, reference=None):
    # Adaptive quadrature checks mass and BOTH tail quantiles. No grid shifts
    # are interpreted as differences in the physical readout.
    cov = isotropic_covariance(1, sf, sa)
    previous = None
    for order in [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072]:
        e = estimate([f], [a], cov, ratio_order=order)
        values = np.array([e.mass, e.quantile(.025), e.quantile(.975)])
        if previous is not None and np.max(np.abs(values/previous-1)) < 2e-4:
            break
        previous = values
    else:
        raise AssertionError('Quadrature did not converge: ' + label)
    result = dict(label=label, force_N=float(f), acceleration_m_s2=float(a),
                  sigma_force_N=float(sf), sigma_acceleration_m_s2=float(sa),
                  ordinary_ratio_kg=float(f/a) if a != 0 else None,
                  extended_mass_kg=e.mass, law_q025_kg=values[1], law_q975_kg=values[2],
                  quadrature_order=order,
                  q_at_reference=profile([f], [a], cov, reference) if reference else None)
    return result


def robot():
    c = load_robot('0-calibration_fts-accel')
    mass, fb, ab, fr, ar = fit_calibration(c)
    fn = load_robot('0-steady-state_wrench')[:, :3]
    an = -load_robot('0-steady-state_accel') @ R.T * G
    # Deliberately retain between-orientation residuals as a systematic/model
    # floor; do not shrink this floor as if every dynamic sample were independent.
    sf = np.sqrt(fn.var(0, ddof=1) + fr.var(0, ddof=1))
    sa = np.sqrt(an.var(0, ddof=1) + ar.var(0, ddof=1))
    _, train_fb, train_ab, _, _ = fit_calibration(c[:12])
    held = c[12:]
    u = held[:, 9:12] / np.linalg.norm(held[:, 9:12], axis=1)[:, None]
    pf = np.sum((held[:, :3] - train_fb)*u, axis=1)
    pa = np.sum((-(held[:, 6:9]-train_ab) @ R.T*G)*u, axis=1)
    # Scalar average of signed projections, not average of individual ratios.
    # Conservative per-reading uncertainty, retained for the mean as a sensitivity
    # budget because calibration/systematic covariance is not supplied.
    examples = [readout('robot_heldout_static_projection_mean', pf.mean(), pa.mean(),
                        np.sqrt(np.mean((u*sf)**2)*3),
                        np.sqrt(np.mean((u*sa)**2)*3), mass)]
    summary = {}
    for run in ['1-baseline', '2-vibrations', '3-vibrations-contact']:
        aa, ff = load_robot(run+'_accel'), load_robot(run+'_wrench')
        oo = load_robot(run+'_orientations')
        origin = min(aa[0, 0], ff[0, 0], oo[0, 0])
        ta = (aa[:, 0]-origin-8416)/1e6
        tf = (ff[:, 0]-origin)/1e6
        valid = (ta >= tf[0]) & (ta <= tf[-1])
        ta, aa = ta[valid], aa[valid]
        forces = np.column_stack([np.interp(ta, tf, ff[:, j+1]) for j in range(3)])-fb
        accels = -(aa[:, 1:4]-ab) @ R.T*G
        f, a = forces[:, 2], accels[:, 2]
        weak = (np.abs(f) <= 2*sf[2]) & (np.abs(a) <= 2*sa[2])
        strong = (np.abs(f) >= 10*sf[2]) & (np.abs(a) >= 10*sa[2]) & (f*a > 0)
        q = (f-mass*a)**2/(sf[2]**2+mass**2*sa[2]**2)
        summary[run] = dict(n=int(len(f)), weak_count=int(weak.sum()),
                            strong_count=int(strong.sum()),
                            q_above_9_count=int((q>9).sum()),
                            strong_ratio_median=float(np.median(f[strong]/a[strong])))
        # Deterministic illustrative selection: first weak, largest |a|, first
        # opposite-sign weak. Includes model failures without suppressing them.
        selected = [('largest_acceleration', int(np.argmax(abs(a))))]
        if weak.any(): selected.append(('first_weak', int(np.flatnonzero(weak)[0])))
        anti = weak & (f*a < 0)
        if anti.any(): selected.append(('first_weak_opposite_sign', int(np.flatnonzero(anti)[0])))
        if run == '1-baseline':
            selected.append(('middle_time', len(a)//2))
        for label, idx in selected:
            e = readout(run+'_'+label, f[idx], a[idx], sf[2], sa[2], mass)
            e.update(time_s=float(ta[idx]), imu_csv_data_row=int(np.flatnonzero(valid)[idx]+1))
            examples.append(e)
        if run == '1-baseline':
            idx = int(np.flatnonzero(weak)[0])
            # Change the uncertainty RATIO while retaining the SAME real reading.
            for factor in [.5, 2.]:
                examples.append(readout('robot_weak_force_sigma_x'+str(factor),
                                        f[idx], a[idx], sf[2]*factor, sa[2], mass))
        np.savetxt(OUT/(run+'_calibrated_z.csv'), np.column_stack([ta, f, a, q, weak, strong]),
                   delimiter=',', header='time_s,Fz_N,az_force_sign_m_s2,Q_at_calibration_mass,weak,strong', comments='')
    return dict(reference_fit_mass_kg=float(mass), published_estimate_kg=.932,
                force_bias_N=fb.tolist(), imu_bias_g=ab.tolist(),
                stationary_force_sigma_N=fn.std(0, ddof=1).tolist(),
                stationary_acceleration_sigma_m_s2=an.std(0, ddof=1).tolist(),
                total_force_sigma_N=sf.tolist(), total_acceleration_sigma_m_s2=sa.tolist(),
                heldout_projection_individual_ratios_kg=(pf/pa).tolist(),
                summary=summary, examples=examples)


def cart():
    data = {n: np.genfromtxt(RAW/'cart'/('cart_'+n+'.csv'), delimiter=',',
                            skip_header=2, usecols=range(1,13)) for n in ['flat','incline']}
    base = data['flat'][:30]
    fb, ab = base[:, 4].mean(), base[:, 5].mean()
    # Export resolution is NOT a manufacturer accuracy specification. Gaussian
    # approximation to uniform rounding, plus observed baseline scatter; include
    # uncertainty of offset mean. x-accelerometer is quantized at 0.1 m/s^2.
    sf = np.sqrt(base[:,4].var(ddof=1)*(1+1/len(base)) + .01**2/12)
    sa = np.sqrt(base[:,5].var(ddof=1)*(1+1/len(base)) + .1**2/12)
    examples, summaries = [], {}
    for name, d in data.items():
        f, a = d[:,4]-fb, d[:,5]-ab
        # Recorded negative bumper pulses: avoids the hand launch and opposite-end
        # collision, which do not pass through the instrumented spring bumper.
        hit = (f < -.2) & (a < -1.)
        for idx in np.flatnonzero(hit):
            e = readout('cart_'+name+'_bumper_t'+str(d[idx,0]), f[idx], a[idx], sf, sa)
            e.update(time_s=float(d[idx,0]), csv_data_row=int(idx+1),
                     exported_kinematic_acceleration_m_s2=float(d[idx,3]),
                     status='sampling/model diagnostic; not validated mass')
            examples.append(e)
        summaries[name] = dict(rows=int(len(d)), sample_period_s=float(np.median(np.diff(d[:,0]))),
                               bumper_rows=int(hit.sum()), bumper_ratios_kg=(f[hit]/a[hit]).tolist())
    # Held-out stationary reading and coasting reading: these deliberately have
    # essentially no mass information. They are not pooled into a mass estimate.
    for t in [1.5, 3.0]:
        d = data['flat']; idx = int(np.argmin(abs(d[:,0]-t)))
        examples.append(readout('cart_flat_weak_t'+str(t), d[idx,4]-fb,
                                float(np.round(d[idx,5]-ab,12)), sf, sa))
    # Explicit constructed exact-zero control with the measured/assumed noise
    # scales. This row is not represented as an observed physical measurement.
    examples.append(readout('cart_constructed_zero_control', 0., 0., sf, sa))
    return dict(force_offset_N=float(fb), acceleration_offset_m_s2=float(ab),
                sigma_force_N=float(sf), sigma_acceleration_m_s2=float(sa),
                zero_limit_mass_scale_kg=float(sf/sa), summaries=summaries, examples=examples)


def run():
    OUT.mkdir(exist_ok=True)
    record = json.loads((RAW/'robot/record.json').read_text())
    for item in record['files']:
        assert 'md5:'+hashlib.md5((RAW/'robot'/item['key']).read_bytes()).hexdigest() == item['checksum']
    manifest = {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in RAW.rglob('*') if p.is_file()}
    result = dict(scope='Scalar restrictions only; exploratory uncertainty budgets, no independent ground truth.',
                  estimator_spec_commit='743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6',
                  raw_sha256=manifest, robot=robot(), cart=cart())
    (OUT/'robot_cart.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    rows = result['robot']['examples'] + result['cart']['examples']
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/'robot_cart_examples.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(json.dumps({k: {n:v for n,v in result[k].items() if n!='examples'} for k in ['robot','cart']}, indent=2))
    for row in rows:
        print(row['label'], 'ordinary=', row['ordinary_ratio_kg'], 'extended=', row['extended_mass_kg'],
              'law95=', [row['law_q025_kg'], row['law_q975_kg']], 'Q=', row['q_at_reference'])


if __name__ == '__main__':
    run()
