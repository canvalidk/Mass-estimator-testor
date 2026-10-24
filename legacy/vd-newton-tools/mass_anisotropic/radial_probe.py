import json
import math
import numpy as np


def radial(mu, variance, order, scale='mean', mapping='u'):
    mu = np.atleast_2d(mu)
    variance = np.broadcast_to(np.atleast_2d(variance), mu.shape)
    mean2 = np.sum(mu * mu, axis=-1)
    if scale == 'mean':
        scale2 = (np.sum(variance, axis=-1) + mean2) / 3
    elif scale == 'geom':
        scale2 = np.prod(variance, axis=-1) ** (1 / 3) + mean2 / 3
    elif scale == 'min':
        scale2 = np.min(variance, axis=-1) + mean2 / 3
    elif scale == 'geom_only':
        scale2 = np.prod(variance, axis=-1) ** (1 / 3)
    else:
        raise ValueError(scale)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (nodes + 1) / 2, weights / 2
    if mapping == 'angle':
        angles = nodes * math.pi / 2
        nodes, weights = np.sin(angles), weights * (math.pi / 2) * np.cos(angles)
    inv, first = np.zeros(len(mu)), np.zeros(len(mu))
    for u, weight in zip(nodes, weights):
        one_minus = 1 - u * u
        den = one_minus + variance / scale2[:, None] * u * u
        log_l = -0.5 * (np.sum(np.log(den), axis=-1) +
                       u * u / scale2 * np.sum(mu * mu / den, axis=-1))
        integrand = np.exp(log_l)
        moment = np.sum(variance * one_minus / den + mu * mu * one_minus**2 / den**2, axis=-1)
        inv += weight * integrand
        first += weight * integrand * moment
    factor = math.sqrt(2 / math.pi) / np.sqrt(scale2)
    return inv * factor, first * factor


if __name__ == '__main__':
    rng = np.random.default_rng(20260919)
    profiles = np.array([[.25, 1, 4], [.25, .25, 4], [.25, 4, 4], [.5, 1, 2], [1,1,1]])
    means, variances = [], []
    for sd in profiles:
        sd /= np.sqrt(np.mean(sd * sd))
        for strength in (0, .1, 1, 3, 8, 25, 64):
            directions = np.vstack([np.eye(3), np.ones(3)/np.sqrt(3), rng.normal(size=(96,3))])
            directions /= np.linalg.norm(directions,axis=1)[:,None]
            means.extend(directions * strength)
            variances.extend([sd*sd] * len(directions))
    mu, c = np.array(means), np.array(variances)
    reference = radial(mu,c,512,'geom')
    for scale in ('mean','geom','geom_only'):
        for order in (16,24,32,48,64,96):
            result = radial(mu,c,order,scale,'angle')
            errors = np.array([np.abs(got/ref-1) for got,ref in zip(result,reference)])
            i,j=np.unravel_index(np.argmax(errors), errors.shape)
            print(json.dumps({'scale':scale,'order':order,'max_relative':np.max(errors,axis=1).tolist(),'worst_mu':mu[j].tolist(),'worst_sd':np.sqrt(c[j]).tolist()}))
