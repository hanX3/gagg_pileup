"""Numerical and future-independence checks for isolated prefix inference."""
import unittest
import numpy as np
from scipy.linalg import cholesky
from scipy.optimize import lsq_linear

from .alpha_gamma_short_window import (
    PrefixCalibration, PrefixConfig, fit_prefix, find_trigger, linear_profile,
)
from .alpha_gamma_waveform_fit import sampled_component
from dataclasses import asdict


def synthetic_calibration():
    readout = dict(rise_ns=2.,decay_ns=20.,lowpass_ns=2.,aperture_ns=.5)
    grid = np.arange(-200.,3000.1,2.)
    covariance_grid = grid[grid<=1250.]
    means = {}
    classes = {}
    for particle,w in [('alpha',.7),('gamma',.21)]:
        curve = 13.618501389*((1-w)*sampled_component(grid,0.,95.,readout)+
                             w*sampled_component(grid,0.,351.,readout))
        means[particle+'_means'] = np.array([curve,curve])
        classes[particle] = {'charge_centers_pe':[500.,10000.]}
    arrays = dict(grid_ns=grid,covariance_grid_ns=covariance_grid,
                  shot_covariance_per_pe=np.eye(len(covariance_grid))*1e-5,**means)
    cfg = asdict(PrefixConfig(onset_tolerance_ns=.0001))
    meta = dict(readout=readout,area_per_pe_mV_ns=13.618501389,baseline_mV=10.,
                electronics_sigma_mV=.1,prefix_config=cfg,classes=classes)
    return PrefixCalibration(meta,arrays)


class PrefixTests(unittest.TestCase):
    def test_active_set_against_independent_bounded_solver(self):
        rng=np.random.default_rng(261091001)
        X=rng.normal(size=(70,2))
        cov=.2**2*.7**abs(np.arange(70)[:,None]-np.arange(70)[None,:])
        L=cholesky(cov,lower=True)
        for coefficients in ([3.,2.],[3.,-2.],[-2.,-1.]):
            y=10.+X@coefficients+.15+rng.normal(0.,.03,70)
            actual=linear_profile(X,y,10.,.2,L)
            A=np.vstack([np.linalg.solve(L,np.c_[X,np.ones(70)]),[0.,0.,5.]])
            b=np.r_[np.linalg.solve(L,y-10.),0.]
            expected=lsq_linear(A,b,bounds=([0.,0.,-np.inf],[np.inf]*3),tol=1e-12)
            np.testing.assert_allclose(actual['amplitudes_pe'],expected.x[:2],rtol=1e-8,atol=1e-8)
            self.assertAlmostEqual(actual['baseline_mV'],10.+expected.x[2],places=8)
            self.assertAlmostEqual(actual['objective'],float(np.sum((A@expected.x-b)**2)),places=7)

    def test_short_window_covariance_is_not_inverse_block(self):
        cal=synthetic_calibration()
        x=np.arange(-20.,100.1,2.)
        indices=np.rint((x+200.)/2).astype(int)
        original=cal.arrays['shot_covariance_per_pe']
        original[:]=.0001*.95**abs(np.arange(len(original))[:,None]-np.arange(len(original))[None,:])
        cov=cal.covariance(x,0.,1000.)
        expected=original[np.ix_(indices,indices)]*1000.+np.eye(len(x))*.01
        np.testing.assert_allclose(cov,expected)
        full=original*1000.+np.eye(len(original))*.01
        self.assertGreater(np.max(abs(np.linalg.inv(cov)-np.linalg.inv(full)[np.ix_(indices,indices)])),.01)

    def test_noiseless_parameter_recovery(self):
        cal=synthetic_calibration();t=np.arange(.25,6500.,2.)
        onset=5000.37; q=3000.
        y=10.+q*cal.meta['area_per_pe_mV_ns']*(.3*sampled_component(t,onset,95.,cal.readout)+
                                             .7*sampled_component(t,onset,351.,cal.readout))
        result=fit_prefix(t,y,cal,pre_ns=50.,post_ns=150.)
        self.assertTrue(result['detected'])
        for name,fit in result['methods'].items():
            self.assertLess(abs(fit['charge_pe']/q-1),.003,name)
            if name.startswith('template'):
                self.assertEqual(fit['particle'],'alpha')
        self.assertLess(abs(result['methods']['physics_ols']['slow_fraction']-.7),.001)

    def test_future_replacement_leaves_trigger_and_all_fits_identical(self):
        cal=synthetic_calibration();t=np.arange(.25,6500.,2.)
        y=10.+3000.*cal.template(t-5000.37,'gamma',3000.)
        a=fit_prefix(t,y,cal,pre_ns=20.,post_ns=100.)
        rng=np.random.default_rng(261091002)
        changed=y.copy();changed[t>a['last_observed_ns']]=rng.normal(1e8,1e6,np.sum(t>a['last_observed_ns']))
        b=fit_prefix(t,changed,cal,pre_ns=20.,post_ns=100.)
        self.assertEqual(a,b)
        trigger=a['trigger']; j=trigger['index']
        # Trigger threshold depends only on the declared twenty-ns lookback.
        local=find_trigger(t[j-10:j+1],y[j-10:j+1],20.,2.)
        self.assertEqual(local['time_ns'],trigger['time_ns'])


if __name__=='__main__':
    unittest.main()
