"""Independent oracles and invariants for joint calibration inference."""
import math
import unittest

import numpy as np

from calibration_posterior import _log_radial, infer_batch
from comparison_posterior import infer_batch as known_infer


def _gamma_expectation_nodes(shape, order):
    # Independent direct precision-mixture oracle, using Laguerre rather
    # than the production collapse/Beta quadrature.
    j = np.arange(order,dtype=float)
    off = np.sqrt(j[1:]*(j[1:]+shape-1))
    nodes,vectors = np.linalg.eigh(np.diag(2*j+shape)+np.diag(off,1)+np.diag(off,-1))
    return nodes,vectors[0]**2


def _mixture_oracle(force, acceleration, sf, sa, degrees, truth):
    kf,ka = np.asarray(degrees)/2
    lf,wf = _gamma_expectation_nodes(kf,96)
    la,wa = _gamma_expectation_nodes(ka,96)
    lf,la = np.meshgrid(lf/(kf*sf**2),la/(ka*sa**2),indexing='ij')
    lf,la = lf.ravel(),la.ravel()
    weights = np.outer(wf,wa).ravel()
    p,q = lf*np.dot(force,force),la*np.dot(acceleration,acceleration)
    dot = np.sqrt(lf*la)*np.dot(force,acceleration)
    tx,tw = np.polynomial.legendre.leggauss(160)
    theta = (tx+1)*math.pi/4
    sn,cs = np.sin(theta),np.cos(theta)
    h2 = p[:,None]*sn**2+q[:,None]*cs**2+2*dot[:,None]*sn*cs
    expo = np.exp((h2-p[:,None]-q[:,None])/2)
    h = np.sqrt(np.maximum(h2,0))
    erf = np.fromiter((math.erf(float(v)/math.sqrt(2)) for v in h.flat),float,h.size).reshape(h.shape)
    radial = expo*np.divide(math.sqrt(math.pi/2)*erf,h,out=np.ones_like(h),where=h>1e-12)
    norm = np.sum(weights*lf*la*np.sum(tw*radial,axis=1))*math.pi/4
    nf = np.sum(weights*np.sqrt(lf)*la*np.sum(tw*sn*expo,axis=1))
    na = np.sum(weights*lf*np.sqrt(la)*np.sum(tw*cs*expo,axis=1))
    yz = math.log(truth)-.5*np.log(la/lf)
    ys,yc = np.sin(np.arctan(np.exp(yz))),np.cos(np.arctan(np.exp(yz)))
    yh2 = p*ys**2+q*yc**2+2*dot*ys*yc
    yh = np.sqrt(np.maximum(yh2,0))
    ye = np.fromiter((math.erf(float(v)/math.sqrt(2)) for v in yh),float,len(yh))
    yr = np.exp((yh2-p-q)/2)*np.divide(math.sqrt(math.pi/2)*ye,yh,out=np.ones_like(yh),where=yh>1e-12)
    density = np.sum(weights*lf*la*yr/(2*np.cosh(yz)))/norm
    return nf/na,-math.log(density)


class CalibrationPosteriorTests(unittest.TestCase):
    def setUp(self):
        self.force = np.array([[2.,.5,0.],[3.,-.4,.2],[0.,0.,0.]])
        self.acceleration = np.array([[1.,-.2,.1],[-2.,.2,.4],[0.,0.,0.]])
        self.options = dict(force_sd=1.2,acceleration_sd=.7,degrees=(12.,18.))

    def test_radial_against_direct_integral(self):
        nodes,weights = np.polynomial.legendre.leggauss(512)
        nodes,weights = (nodes+1)/2,weights/2
        h2 = np.array([0.,1e-14,1e-6,.1,2.,100.,10000.])
        residual = np.array([2.,0.,1.,10.,0.,2.,50.])
        for power in (3.,11.,89.,1001.):
            scale = -power*np.log1p(residual/2)
            reference = scale+np.log(np.sum(weights*np.exp(-power*np.log1p(h2[:,None]*nodes**2/(2+residual[:,None]))),axis=1))
            np.testing.assert_allclose(_log_radial(h2,residual,power),reference,rtol=0,atol=2e-11)

    def test_direct_two_precision_mixture_oracle(self):
        force,acceleration = self.force[0],self.acceleration[0]
        reference_point,reference_score = _mixture_oracle(force,acceleration,1.2,.7,(12.,18.),1.)
        result = infer_batch(force,acceleration,calibration_order=64,order=192,**self.options)
        self.assertAlmostEqual(result['points']['calibrated_joint'][0],reference_point,places=8)
        self.assertAlmostEqual(result['distributions']['calibrated_joint']['log_density_score'][0],reference_score,places=8)

    def test_zero_data_analytic_physical_first_moments(self):
        def mean(sd,nu):
            k = nu/2
            return sd*math.sqrt(nu/math.pi)*math.exp(math.lgamma(k+.5)-math.lgamma(k+1))
        result = infer_batch([0.,0.,0.],[0.,0.,0.],force_sd=1.4,acceleration_sd=.8,
                             degrees=(9.,15.),calibration_order=64)
        self.assertAlmostEqual(result['points']['calibrated_joint'][0],mean(1.4,9.)/mean(.8,15.),places=9)

    def test_channel_exchange(self):
        forward = infer_batch(self.force,self.acceleration,**self.options)
        reverse = infer_batch(self.acceleration,self.force,force_sd=.7,acceleration_sd=1.2,degrees=(18.,12.))
        for method in forward['points']:
            np.testing.assert_allclose(forward['points'][method]*reverse['points'][method],1.,rtol=3e-9)
        left,right = forward['distributions']['calibrated_joint'],reverse['distributions']['calibrated_joint']
        np.testing.assert_allclose(left['log_lower_95'],-right['log_upper_95'],atol=2e-8)

    def test_change_of_units(self):
        forward = infer_batch(self.force,self.acceleration,**self.options)
        scaled = infer_batch(self.force*2,self.acceleration/3,truth=6.,force_sd=2.4,
                            acceleration_sd=.7/3,degrees=(12.,18.))
        for name in forward['points']:
            np.testing.assert_allclose(scaled['points'][name],6*forward['points'][name],rtol=3e-8)
        left,right = forward['distributions']['calibrated_joint'],scaled['distributions']['calibrated_joint']
        for metric in ('cdf_truth','log_density_score','log_crps','log_sd'):
            np.testing.assert_allclose(left[metric],right[metric],atol=3e-8)
        np.testing.assert_allclose(right['log_lower_95'],left['log_lower_95']+math.log(6),atol=3e-8)

    def test_truth_only_scores(self):
        one = infer_batch(self.force,self.acceleration,truth=1.,**self.options)
        two = infer_batch(self.force,self.acceleration,truth=100.,**self.options)
        for method in one['points']:
            np.testing.assert_array_equal(one['points'][method],two['points'][method])
        np.testing.assert_array_equal(one['distributions']['calibrated_joint']['log_lower_95'],
                                      two['distributions']['calibrated_joint']['log_lower_95'])
        self.assertTrue(np.all(one['distributions']['calibrated_joint']['cdf_truth'] <
                               two['distributions']['calibrated_joint']['cdf_truth']))

    def test_extreme_change_of_units_preserves_relative_grid(self):
        baseline = infer_batch(self.force,self.acceleration,**self.options)
        for factor in (1e-16,1e16):
            changed = infer_batch(self.force*factor,self.acceleration,truth=factor,
                                 force_sd=1.2*factor,acceleration_sd=.7,degrees=(12.,18.))
            for method in baseline['points']:
                np.testing.assert_allclose(changed['points'][method]/factor,baseline['points'][method],rtol=2e-13)
            before,after = baseline['distributions']['calibrated_joint'],changed['distributions']['calibrated_joint']
            for metric in ('cdf_truth','log_density_score','log_crps','log_sd'):
                np.testing.assert_allclose(after[metric],before[metric],atol=2e-13)
            np.testing.assert_allclose(after['log_lower_95']-math.log(factor),before['log_lower_95'],atol=2e-13)

    def test_known_noise_limit(self):
        force = np.array([self.force[0],[32.,.2,0.]])
        acceleration = np.array([self.acceleration[0],[8.,-.2,.1]])
        truth = np.array([1.,4.])
        result = infer_batch(force,acceleration,truth=truth,force_sd=1.2,acceleration_sd=.7,degrees=100000.,calibration_order=24)
        known = known_infer(force/1.2,acceleration/.7,truth=truth*.7/1.2)
        np.testing.assert_allclose(result['points']['calibrated_joint'],known['points']['flat_joint']*1.2/.7,rtol=2e-6)
        calibrated,reference = result['distributions']['calibrated_joint'],known['distributions']['flat_joint']
        for metric in ('cdf_truth','log_sd','log_crps','log_density_score'):
            np.testing.assert_allclose(calibrated[metric],reference[metric],atol=2e-5)
        np.testing.assert_allclose(calibrated['log_lower_95'],reference['log_lower_95']+math.log(1.2/.7),atol=3e-5)

    def test_independent_resolution_refinement(self):
        options = dict(force_sd=1.2,acceleration_sd=.7,degrees=9.)
        ordinary = infer_batch(self.force,self.acceleration,**options)
        refined = infer_batch(self.force,self.acceleration,order=192,calibration_order=48,**options)
        for method in ordinary['points']:
            np.testing.assert_allclose(ordinary['points'][method],refined['points'][method],rtol=3e-6)
        for metric in ('log_lower_95','log_upper_95','cdf_truth','log_crps'):
            np.testing.assert_allclose(ordinary['distributions']['calibrated_joint'][metric],
                                       refined['distributions']['calibrated_joint'][metric],atol=3e-6)

    def test_batch_matches_individual_means(self):
        sf,sa = np.array([1.2,.8,2.]),np.array([.7,1.1,.9])
        batch = infer_batch(self.force,self.acceleration,force_sd=sf,acceleration_sd=sa,degrees=9.)
        for i in range(3):
            single = infer_batch(self.force[i],self.acceleration[i],force_sd=sf[i],acceleration_sd=sa[i],degrees=9.)
            np.testing.assert_allclose(batch['points']['calibrated_joint'][i],single['points']['calibrated_joint'][0],rtol=2e-14)

    def test_strong_signal_with_uncertain_calibration_refines(self):
        force = np.array([[32.,.2,-.4],[31.5,-.3,.1],[8.2,.4,.1]])
        acceleration = np.array([[8.,-.2,.1],[7.6,.3,.1],[7.8,.1,-.3]])
        for degrees in (9.,87.):
            options = dict(force_sd=np.array([.3,.8,.45]),acceleration_sd=np.array([.7,.3,.5]),degrees=degrees)
            ordinary = infer_batch(force,acceleration,**options)
            refined = infer_batch(force,acceleration,order=192,calibration_order=48,**options)
            np.testing.assert_allclose(ordinary['points']['calibrated_joint'],refined['points']['calibrated_joint'],rtol=2e-6)
            np.testing.assert_allclose(ordinary['distributions']['calibrated_joint']['log_lower_95'],
                                       refined['distributions']['calibrated_joint']['log_lower_95'],atol=2e-6)

    def test_rejects_invalid_calibration_inputs(self):
        for update in ({'force_sd':0.},{'acceleration_sd':float('nan')},{'degrees':0.},
                       {'degrees':(3.,-1.)},{'degrees':[3.,4.,5.]},{'calibration_order':7},
                       {'force_sd':np.ones((3,1))}):
            options = dict(self.options)
            options.update(update)
            with self.assertRaises(ValueError):
                infer_batch(self.force,self.acceleration,**options)


if __name__ == '__main__':
    unittest.main()
