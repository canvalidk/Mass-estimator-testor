"""Scratch checks of a proposed ours26 continuation, not production code.

Source snapshot: VD-docs c877b1df131bc3a809ef25bc270160eb5b11b6d5.
No simulation-performance or uniqueness claim is made by these checks.
Requires numpy only. Run this file to reproduce results.json beside it.
"""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss


def gl(a, b, n):
    z, w = leggauss(n)
    return (a + b) / 2 + (b - a) * z / 2, w * (b - a) / 2


LOGT, LOGW = gl(-36, 36, 220)


def normal_norm_spectrum(variances, mean_squares, multiplicities=None):
    """Mean Euclidean norm from Gaussian eigenvalues and squared means.

Multiplicity is the number of independent copies of each eigenvalue;
mean_squares already sums the squared means over those copies.
"""
    v, z = np.asarray(variances), np.asarray(mean_squares)
    n = np.ones_like(v) if multiplicities is None else np.asarray(multiplicities)
    ex = np.sum(n * v + z)
    t = np.exp(LOGT) / ex
    den = 1 + 2 * t[:, None] * v
    lp = -.5 * (np.log1p(2 * t[:, None] * v) @ n)
    lp -= (t[:, None] * z / den).sum(axis=1)
    value = np.dot(LOGW, -np.expm1(lp) * np.exp(-LOGT / 2))
    # Leading analytic tails; neglected contributions at +/-36 are negligible.
    value += 2 * np.exp(-18) + 2 * np.exp(-18)
    return float(np.sqrt(ex) * value / (2 * np.sqrt(np.pi)))


def norm_mean(b, k):
    eig, u = np.linalg.eigh(k)
    return normal_norm_spectrum(eig, (u.T @ b) ** 2)


def block_diag(matrices):
    sizes = [len(a) for a in matrices]
    out = np.zeros((sum(sizes), sum(sizes)))
    start = 0
    for a, size in zip(matrices, sizes):
        out[start:start + size, start:start + size] = a
        start += size
    return out


def conditional(y, sigma, mass):
    """Precision formula for latent a; independent residual formula for Q."""
    d = len(y) // 2
    bmat = np.vstack((mass * np.eye(d), np.eye(d)))
    wb = np.linalg.solve(sigma, bmat)
    k = np.linalg.inv(bmat.T @ wb)
    j = k @ wb.T
    b = j @ y
    residual = y[:d] - mass * y[d:]
    v = sigma[:d, :d] - mass * (sigma[:d, d:] + sigma[d:, :d])
    v += mass * mass * sigma[d:, d:]
    q = residual @ np.linalg.solve(v, residual)
    return float(q), b, k, j


def direct_terms(y, sigma, masses):
    q, size = [], []
    for mass in masses:
        qi, b, k, _ = conditional(y, sigma, mass)
        q.append(qi)
        size.append(norm_mean(b, k))
    return np.asarray(q), np.asarray(size)


def state_terms(groups, masses):
    """Each group is (Sigma, S=sum y y^T, n). Exact additive state."""
    qs, sizes = [], []
    for mass in masses:
        qsum, variances, mean_squares, multiplicities = 0., [], [], []
        for sigma, scatter, count in groups:
            d = len(sigma) // 2
            _, _, k, j = conditional(np.zeros(2*d), sigma, mass)
            bmat = np.vstack((mass * np.eye(d), np.eye(d)))
            w = np.linalg.inv(sigma)
            r = w - w @ bmat @ j
            qsum += np.trace(r @ scatter)
            eig, u = np.linalg.eigh(k)
            means_scatter = j @ scatter @ j.T
            variances.extend(eig)
            mean_squares.extend(np.diag(u.T @ means_scatter @ u))
            multiplicities.extend([count] * d)
        qs.append(qsum)
        sizes.append(normal_norm_spectrum(variances, mean_squares, multiplicities))
    return np.asarray(qs), np.asarray(sizes)


def halfcauchy_grid(center, n=160):
    theta, w = gl(0, np.pi / 2, n)
    return center * np.tan(theta), w * 2 / np.pi


def point(q, size, masses, prior_weights):
    w = prior_weights * np.exp(-.5 * (q - np.min(q)))
    return float(np.dot(w, masses * size) / np.dot(w, size))


def infer(y, sigma, center=1., order=160):
    m, w = halfcauchy_grid(center, order)
    q, size = direct_terms(y, sigma, m)
    return point(q, size, m, w)


def check_near(name, actual, expected, tol, results):
    error = float(np.max(np.abs(np.asarray(actual) - expected)))
    if not error <= tol:
        raise AssertionError((name, error, tol, actual, expected))
    results[name] = error


def main():
    rng = np.random.default_rng(20260927)
    results = {"source_commit": "c877b1df131bc3a809ef25bc270160eb5b11b6d5"}

    # Norm integral checked against independent exact scalar and central formulas.
    errors = []
    for sd in (.2, 1., 4.):
        for mu in (0., .4, 4., -2.):
            h = abs(mu) / sd
            exact = sd * np.sqrt(2 / np.pi) * np.exp(-h*h/2)
            exact += abs(mu) * math.erf(h / np.sqrt(2))
            errors.append(abs(norm_mean(np.array([mu]), np.array([[sd*sd]])) - exact))
    check_near("folded_normal_max_error", errors, 0., 2e-8, results)
    for d in (2, 3, 8):
        exact = np.sqrt(2) * math.gamma((d+1)/2) / math.gamma(d/2)
        check_near("central_norm_d%d" % d, norm_mean(np.zeros(d), np.eye(d)), exact, 2e-8, results)

    # Residual/covariance and precision/least-squares forms checked independently.
    identity_errors = []
    for d in (1, 2, 3):
        for case in range(5):
            matrix = rng.normal(size=(2*d, 2*d))
            sigma = matrix @ matrix.T + .4 * np.eye(2*d)
            y = rng.normal(size=2*d)
            for mass in (.07, .7, 7.):
                q, b, k, _ = conditional(y, sigma, mass)
                bmat = np.vstack((mass*np.eye(d), np.eye(d)))
                e = y - bmat @ b
                identity_errors.append(abs(q - e @ np.linalg.solve(sigma, e)))
                cross = sigma[d:, :d] - mass*sigma[d:, d:]
                v = sigma[:d, :d] - mass*(sigma[:d, d:]+sigma[d:, :d]) + mass*mass*sigma[d:, d:]
                b2 = y[d:] - cross @ np.linalg.solve(v, y[:d]-mass*y[d:])
                k2 = sigma[d:, d:] - cross @ np.linalg.solve(v, cross.T)
                identity_errors.extend(np.abs(b-b2))
                identity_errors.extend(np.abs(k-k2).flat)
    check_near("conditional_identity_max_error", identity_errors, 0., 2e-11, results)

    # Fully correlated examples: channel swap, spatial rotation and force units.
    d = 3
    matrix = rng.normal(size=(2*d, 2*d))
    sigma = matrix @ matrix.T + np.eye(2*d)
    y = np.array([2., -1., 3., .4, 1., -.2])
    center = 1.3
    value = infer(y, sigma, center)
    results["correlated_example_mass"] = value
    check_near("mass_quadrature_refinement", infer(y, sigma, center, 320), value, 2e-8, results)
    swap = np.r_[np.arange(d, 2*d), np.arange(d)]
    inverse = infer(y[swap], sigma[np.ix_(swap, swap)], 1/center)
    check_near("reciprocity_product_error", value * inverse, 1., 2e-9, results)
    rot, _ = np.linalg.qr(rng.normal(size=(d, d)))
    transform = block_diag([rot, rot])
    check_near("rotation_error", infer(transform@y, transform@sigma@transform.T, center), value, 2e-9, results)
    scale = np.diag([7.]*d + [1.]*d)
    check_near("force_units_error", infer(scale@y, scale@sigma@scale.T, 7*center), 7*value, 2e-8, results)

    # Independent diagonal component formula (ours26) against full matrix formula.
    sf, sa = np.array([1., 2., .7]), np.array([.2, .8, 1.3])
    diag = np.diag(np.r_[sf*sf, sa*sa])
    diag_errors = []
    for mass in (.2, 1., 4.):
        q, b, k, _ = conditional(y, diag, mass)
        v = sf*sf + mass*mass*sa*sa
        expected_b = (mass*sa*sa*y[:d] + sf*sf*y[d:]) / v
        expected_k = sf*sf*sa*sa / v
        diag_errors.extend(np.abs(b-expected_b))
        diag_errors.extend(np.abs(k-np.diag(expected_k)).flat)
        diag_errors.append(abs(q - np.sum((y[:d]-mass*y[d:])**2 / v)))
    check_near("ours26_reduction_error", diag_errors, 0., 1e-11, results)

    # Exact state for full covariance blocks compared with all original readings.
    groups, readings, covariances, order = [], [], [], []
    for group, (d, n) in enumerate(((2, 3), (2, 2))):
        matrix = rng.normal(size=(2*d, 2*d))
        sig = matrix @ matrix.T + .6*np.eye(2*d)
        ys = rng.normal(size=(n, 2*d))
        scatter = ys.T @ ys
        groups.append((sig, scatter, n))
        for yi in ys:
            offset = sum(len(yy) for yy in readings)
            order.append((list(range(offset, offset+d)), list(range(offset+d, offset+2*d))))
            readings.append(yi)
            covariances.append(sig)
    permutation = np.array(sum((i[0] for i in order), []) + sum((i[1] for i in order), []))
    yy = np.concatenate(readings)[permutation]
    ss = block_diag(covariances)[np.ix_(permutation, permutation)]
    masses, weights = halfcauchy_grid(1.4, 180)
    qa, aa = direct_terms(yy, ss, masses)
    qs, ast = state_terms(groups, masses)
    check_near("state_Q_error", qs, qa, 2e-11, results)
    check_near("state_norm_error", ast, aa, 2e-10, results)
    check_near("state_readout_error", point(qs, ast, masses, weights), point(qa, aa, masses, weights), 2e-11, results)

    # Priors with the same symmetry/centre give genuinely different answers.
    y1, sigma1 = np.array([1., 0.]), np.eye(2)
    prior_examples = {"half_cauchy": infer(y1, sigma1, order=320)}
    for width in (.3, 1., 2.):
        z, w = hermgauss(160)
        masses, weights = np.exp(np.sqrt(2)*width*z), w/np.sqrt(np.pi)
        q, size = direct_terms(y1, sigma1, masses)
        prior_examples["lognormal_logsd_%g" % width] = point(q, size, masses, weights)
    results["same_center_different_prior"] = prior_examples

    # Optional geometric default: uniform Euclidean arc length of the curve of
    # whitened compatible subspaces. This is an EXTRA premise, not forced here.
    # Check trace(K H) expression against finite differences of the projector.
    matrix = rng.normal(size=(4, 4))
    sg = matrix @ matrix.T + np.eye(4)
    ev, eu = np.linalg.eigh(sg)
    root_inv = (eu * (1 / np.sqrt(ev))) @ eu.T
    w = np.linalg.inv(sg)
    ee = np.vstack((np.eye(2), np.zeros((2, 2))))
    def projector(mass):
        bm = np.vstack((mass * np.eye(2), np.eye(2)))
        bw = root_inv @ bm
        return bw @ np.linalg.inv(bw.T @ bw) @ bw.T
    speed_errors = []
    for mass in (.2, 1., 4.):
        _, _, k, j = conditional(np.zeros(4), sg, mass)
        bm = np.vstack((mass * np.eye(2), np.eye(2)))
        rr = w - w @ bm @ j
        hh = ee.T @ rr @ ee
        derivative = (projector(mass + 1e-5) - projector(mass - 1e-5)) / 2e-5
        speed_errors.append(abs(np.trace(k @ hh) - .5*np.sum(derivative*derivative)))
    check_near("geometry_speed_squared_error", speed_errors, 0., 2e-9, results)
    # Distinct invariant norms give different valid diagonal geometry priors.
    geometric_examples = {}
    yy = np.array([2., 0., .2, 0.])
    sg = np.diag([1., 16., 1., 1.])
    masses, base_w = halfcauchy_grid(2., 360)
    q, size = direct_terms(yy, sg, masses)
    base_density = (2 / np.pi) * 2 / (4 + masses*masses)
    j = np.array([1., 4.])[None, :] / (np.array([1., 16.])[None, :] + masses[:, None]**2)
    for power in (1, 2, 4):
        density = np.sum(j**power, axis=1)**(1/power)
        geometric_examples["angular_speed_L%d" % power] = point(q, size, masses, base_w*density/base_density)
    results["geometry_nonuniqueness_example"] = geometric_examples

    # Prior-once updating with mixed instruments: all-zero factors equal one,
    # so a first-instrument prior remains its original (different) law.
    zero_data = np.zeros(4)
    sigma_zero = np.diag([1., 16., 1., 1.])
    results["zero_data_first_instrument_order_example"] = {
        "prior_medians": [1., 4.],
        "mass_first_ratio_1": infer(zero_data, sigma_zero, 1.),
        "mass_first_ratio_4": infer(zero_data, sigma_zero, 4.),
    }

    # Fixed-mass Q has chi-square M distribution even with cross-correlations.
    dim, trials = 6, 100000
    matrix = rng.normal(size=(2*dim, 2*dim))
    sig = matrix @ matrix.T + np.eye(2*dim)
    true_mass = 2.3
    true_acc = rng.normal(size=dim)
    errors = rng.multivariate_normal(np.zeros(2*dim), sig, size=trials)
    data = errors + np.r_[true_mass*true_acc, true_acc]
    rr = data[:, :dim] - true_mass*data[:, dim:]
    vv = sig[:dim, :dim] - true_mass*(sig[:dim, dim:]+sig[dim:, :dim]) + true_mass**2*sig[dim:, dim:]
    qq = np.sum(rr * np.linalg.solve(vv, rr.T).T, axis=1)
    # Exact integer df6 quantile solves 1-exp(-x/2)*(1+x/2+(x/2)^2/2)=.95.
    lower, upper = 0., 50.
    for _ in range(70):
        middle = (lower+upper)/2
        cdf = 1-np.exp(-middle/2)*(1+middle/2+middle*middle/8)
        if cdf < .95: lower = middle
        else: upper = middle
    coverage = float(np.mean(qq <= (lower+upper)/2))
    results["confidence_set_check"] = {"trials": trials, "dimension": dim, "coverage": coverage, "MCSE": np.sqrt(.95*.05/trials), "Q_mean": float(np.mean(qq))}
    assert abs(coverage-.95) < 5*np.sqrt(.95*.05/trials)

    path = Path(__file__).with_name("results.json")
    path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
