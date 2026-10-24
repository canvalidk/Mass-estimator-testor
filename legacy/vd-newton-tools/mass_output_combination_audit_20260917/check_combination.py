"""Scratch verification of declared output-pooling rules, not production policy.

The dm-product rule uses common-mass reference dm * product(alpha_i d alpha_i
dOmega_i). The log-product rule uses m**(N-1) times that reference. Both are
explicit choices, not consequences of the original single-trial estimator.
"""
from functools import lru_cache
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / '03_trace_and_evaluator/mass_estimator'))
from estimator import estimate, isotropic_covariance


@lru_cache(None)
def angle_rule(order):
    x, w = np.polynomial.legendre.leggauss(order)
    return (x + 1) * math.pi / 4, w * math.pi / 4


def log_j1(h2):
    h2 = np.maximum(h2, 0)
    h = np.sqrt(h2)
    out = np.empty_like(h)
    small = h < 1e-6
    out[small] = np.log1p(h2[small] / 3 + h2[small] ** 2 / 15)
    erf = np.array([math.erf(float(v) / math.sqrt(2)) for v in h[~small]])
    out[~small] = h2[~small] / 2 + .5 * math.log(math.pi / 2) + np.log(erf) - np.log(h[~small])
    return out


class Trial:
    def __init__(self, force, acceleration, sf, sa):
        self.force, self.acceleration = np.array(force), np.array(acceleration)
        self.sf, self.sa, self.s = sf, sa, sf / sa
        u, v = self.force / sf, self.acceleration / sa
        self.A, self.B, self.C = u @ u, v @ v, u @ v

    def log_angle_kernel(self, theta):
        sn, cs = np.sin(theta), np.cos(theta)
        h2 = self.A * sn**2 + 2*self.C*sn*cs + self.B*cs**2
        return log_j1(h2), h2

    def curves(self, mass, order):
        theta, weights = angle_rule(order)
        kernel, _ = self.log_angle_kernel(theta)
        peak = kernel.max()
        log_norm = peak + math.log(float(weights @ np.exp(kernel - peak)))
        local_theta = np.arctan(mass / self.s)
        kernel, h2 = self.log_angle_kernel(local_theta)
        log_q = kernel - log_norm + math.log(self.s) - np.log(mass**2 + self.s**2)
        # E[alpha | m] requires the radial first moment, not just the mass density.
        mean_alpha = self.sa * np.cos(local_theta) * np.exp(.5*math.log(math.pi/2) + h2/2 - kernel)
        return log_q, mean_alpha


def normalize(log_q, weights):
    peak = log_q.max()
    log_z = peak + math.log(float(weights @ np.exp(log_q - peak)))
    return log_q - log_z


def readout(log_q, a, mass, weights):
    p = np.exp(log_q)
    return float(weights @ (mass*a*p) / (weights @ (a*p)))


def run(order):
    theta, wt = angle_rule(order)
    mass, weights = np.tan(theta), wt / np.cos(theta)**2
    trials = [Trial([2,.2,0], [1,0,0], .5,.4),
              Trial([5,0,.1], [2.1,.2,0], 1,.6),
              Trial([1.5,-.1,.2], [.5,.1,0], .8,.35)]
    curves = [t.curves(mass, order) for t in trials]
    out = {'order': order, 'trials': []}
    for t, (log_q, a) in zip(trials, curves):
        ref = estimate(t.force, t.acceleration, isotropic_covariance(3,t.sf,t.sa),
                       direction_order=64, ratio_order=2048)
        point = readout(log_q, a, mass, weights)
        out['trials'].append({'mass': point, 'reference_mass': ref.mass,
                              'mass_absolute_difference': abs(point-ref.mass),
                              'normalization': float(weights @ np.exp(log_q))})

    for name, power in [('mass_coordinate',0), ('log_coordinate',len(trials)-1)]:
        direct = normalize(sum(q for q,a in curves) + power*np.log(mass), weights)
        pair = normalize(curves[0][0] + curves[1][0] + (power!=0)*np.log(mass), weights)
        grouped = normalize(pair + curves[2][0] + (power!=0)*np.log(mass), weights)
        a_sum = sum(a for q,a in curves)
        out[name] = {'max_grouping_log_density_difference': float(np.max(np.abs(direct-grouped))),
                     'mass_readout': readout(direct,a_sum,mass,weights),
                     'mean_mass': float(weights @ (mass*np.exp(direct)))}

    # Analytic null fixtures (s=1), including the original readout weighting.
    null = Trial([0,0,0], [0,0,0], 1,1)
    q, a = null.curves(mass, order)
    pooled = normalize(2*q,weights)
    pooled_log = normalize(2*q+np.log(mass),weights)
    exact = 4/math.pi/(1+mass**2)**2
    exact_log = 2*mass/(1+mass**2)**2
    out['null'] = {
        'mass_product_max_relative_difference': float(np.max(np.abs(np.exp(pooled)/exact-1))),
        'log_product_max_relative_difference': float(np.max(np.abs(np.exp(pooled_log)/exact_log-1))),
        'mass_product_readout': readout(pooled,2*a,mass,weights),
        'log_product_readout': readout(pooled_log,2*a,mass,weights),
        'mass_product_mean': float(weights @ (mass*np.exp(pooled))),
    }

    # Same mass curve and point under a common physical scaling of both channels,
    # but a changed conditional magnitude curve changes a sum-magnitude readout.
    original = trials[0]
    scaled = Trial(10*original.force,10*original.acceleration,10*original.sf,10*original.sa)
    qs, ass = scaled.curves(mass,order)
    pair = normalize(curves[0][0]+curves[1][0],weights)
    out['scale_witness'] = {
        'max_single_log_density_change': float(np.max(np.abs(qs-curves[0][0]))),
        'single_point_change': abs(readout(qs,ass,mass,weights)-readout(*curves[0],mass,weights)),
        'combined_point_original': readout(pair,curves[0][1]+curves[1][1],mass,weights),
        'combined_point_scaled': readout(pair,ass+curves[1][1],mass,weights),
    }
    return out


if __name__ == '__main__':
    results = [run(n) for n in (256,512)]
    for r in results:
        assert max(t['mass_absolute_difference'] for t in r['trials']) < 1e-8
        assert abs(r['null']['mass_product_readout']-.5) < 1e-10
        assert abs(r['null']['log_product_readout']-1) < 1e-10
        assert abs(r['null']['mass_product_mean']-2/math.pi) < 1e-10
        assert r['mass_coordinate']['max_grouping_log_density_difference'] < 1e-11
        assert r['log_coordinate']['max_grouping_log_density_difference'] < 1e-11
    output = Path(__file__).with_name('checks.json')
    output.write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps(results,indent=2))
