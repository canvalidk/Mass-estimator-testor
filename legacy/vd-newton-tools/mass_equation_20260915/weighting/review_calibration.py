"""Independent Cartesian-wedge check of selected CDF calculations; no full rerun."""
from pathlib import Path
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from calibration_probe import infer_batch
from gaussian_conditioning import sphere_average


def wedge_cdfs(force, acceleration, mass, order=256):
    force, acceleration = np.asarray(force), np.asarray(acceleration)
    upper = np.sqrt(force@force+acceleration@acceleration)+16.
    t, w = np.polynomial.legendre.leggauss(order)
    t, w = (t+1)/2, w/2

    def integrate(a0, a1, limit):
        alpha = a0+(a1-a0)*t[:, None]
        fmax = limit(alpha)
        f = fmax*t[None, :]
        weight = ((a1-a0)*w[:, None]*fmax*w[None, :]
                  *sphere_average(f, alpha, force, acceleration))
        return np.array([np.sum(weight),
                         np.sum(weight*(f*f+alpha*alpha)),
                         np.sum(weight*(f*alpha)**2)])

    denominator = integrate(0., upper, lambda a: upper+0*a)
    split = min(upper, upper/mass)
    numerator = integrate(0., split, lambda a: mass*a)
    if split < upper:
        numerator += integrate(split, upper, lambda a: upper+0*a)
    return numerator/denominator


def run():
    noise = np.random.default_rng(20260915).normal(size=(2048, 2, 3))
    indices = (0, int(np.argmin(noise[:,1,0])), int(np.argmax(noise[:,1,0])))
    worst = 0.
    for mass, excitation in ((4., 0.), (16., .25), (16., 4.)):
        for idx in indices:
            force, acceleration = noise[idx].copy()
            force[0] += mass*excitation
            acceleration[0] += excitation
            _, original = infer_batch([force], [acceleration], mass)
            cartesian = wedge_cdfs(force, acceleration, mass)
            error = float(np.max(abs(original[0]-cartesian)))
            worst = max(worst, error)
            assert error < 2e-8, (mass, excitation, idx, original, cartesian)
            assert np.array_equal((original[0]>=.025)&(original[0]<=.975),
                                  (cartesian>=.025)&(cartesian<=.975))
            print(mass, excitation, idx, error)
    print('Worst absolute CDF discrepancy:', worst)
    print('All selected coverage classifications agree.')


if __name__ == '__main__':
    run()
