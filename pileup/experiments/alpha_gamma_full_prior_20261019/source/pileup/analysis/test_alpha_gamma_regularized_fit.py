"""Independent algebraic and physical controls for the regularized fitter."""
import json
from pathlib import Path
import unittest
import numpy as np
from .alpha_gamma_regularized_fit import (CalibratedBank,FrozenWhitening,RegularizedConfig,
                                          fit_regularized,positive_lasso)


class RegularizedChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root=Path(__file__).resolve().parents[2]
        readout=json.loads((root/'g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json').read_text())['readout']
        cls.cal=dict(readout=readout,area_per_pe_mV_ns=13.618501389,baseline_mV=10.,
            classes={p:dict(slow_intercept=f,slow_log_slope=0.,delay_ns=0.,charge_range_pe=[50.,20000.])
                     for p,f in (('alpha',.70),('gamma',.21))},
            noise=dict(baseline_variance_mV2=.01,variance_per_mV=.02,shot_rho_at_2ns=.9))
        cls.bank=CalibratedBank(cls.cal);cls.t=np.arange(.25,12000.,2.)

    def test_positive_lasso_known_orthogonal_solution(self):
        target=np.array([10.,-2.,3.,1.]);x,details=positive_lasso(np.eye(4),target,2.)
        np.testing.assert_allclose(x,[8.,0.,1.,0.],atol=1e-7)
        self.assertLess(details['projected_gradient_max'],1e-6)

    def test_whitening_inverts_heteroscedastic_covariance(self):
        w=FrozenWhitening(np.arange(1.,9.),2.,self.cal['noise'])
        correlation=np.eye(8)
        for i in range(8):
            for j in range(i+1,8):correlation[i,j]=correlation[j,i]=np.prod(w.rho[i:j])
        cov=correlation*w.sigma[:,None]*w.sigma[None,:]
        transform=w.apply(np.eye(8))
        np.testing.assert_allclose(transform@cov@transform.T,np.eye(8),atol=1e-12)

    def test_weak_component_kept_but_not_reported(self):
        times=[3000.,3400.,6000.];types=[0,1,1];charges=[3000.,180.,1200.]
        y=10.+self.bank.matrix(self.t,times,types,charges)@charges
        fit=fit_regularized(self.t,y,self.cal,RegularizedConfig(selection_penalty=3.))
        self.assertEqual(fit['particle_name'],['alpha','gamma'])
        np.testing.assert_allclose(fit['charge_pe'],[3000.,1200.],rtol=.01)
        np.testing.assert_allclose(fit['times_ns'],[3000.,6000.],atol=.2)
        self.assertGreaterEqual(fit['nuisance_component_count'],1)
        self.assertLess(fit['rms_mV'],.05)
        history=np.asarray(fit['weighted_cost_history'])
        self.assertTrue(np.all(np.diff(history)<=1e-7))

    def test_empty_electronics_noise_has_no_reported_pulse(self):
        y=10.+np.random.default_rng(261030).normal(0.,.1,len(self.t))
        fit=fit_regularized(self.t,y,self.cal)
        self.assertEqual(fit['times_ns'],[])

    def test_close_pair_unknown_order(self):
        times=[5000.,5040.];q=[3000.,1200.]
        for kinds in ([0,1],[1,0]):
            y=10.+self.bank.matrix(self.t,times,kinds,q)@q
            fit=fit_regularized(self.t,y,self.cal,RegularizedConfig(selection_penalty=3.))
            self.assertEqual(fit['particle_name'],self.bank.names[kinds].tolist())
            np.testing.assert_allclose(fit['times_ns'],times,atol=2.)
            np.testing.assert_allclose(fit['charge_pe'],q,rtol=.1)


if __name__=='__main__':unittest.main()
