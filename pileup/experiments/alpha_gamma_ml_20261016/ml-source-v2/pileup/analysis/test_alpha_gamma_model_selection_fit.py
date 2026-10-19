"""Regression checks for missed second pulses and approximate sparse solutions."""
import json
from pathlib import Path
import unittest
import numpy as np
from .alpha_gamma_model_selection_fit import strict_positive_lasso,fit_model_selection,ModelSelectionConfig
from .alpha_gamma_regularized_fit import CalibratedBank


class JointSelectionChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root=Path(__file__).resolve().parents[2]
        readout=json.loads((root/'g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json').read_text())['readout']
        cls.cal=dict(readout=readout,area_per_pe_mV_ns=13.618501389,baseline_mV=10.,
            classes={p:dict(slow_intercept=f,slow_log_slope=s,delay_ns=d,charge_range_pe=[250.,20000.])
                     for p,f,s,d in (('alpha',.681,.0324,-1.91),('gamma',.243,.0047,-.397))},
            noise=dict(baseline_variance_mV2=.01,variance_per_mV=.614,shot_rho_at_2ns=.9))
        cls.bank=CalibratedBank(cls.cal);cls.t=np.arange(.25,12000.,2.)

    def test_sparse_kkt_with_known_nonnegative_optimum(self):
        rng=np.random.default_rng(261041001)
        a=rng.uniform(0,1,(160,30));a/=np.linalg.norm(a,axis=0)
        expected=np.zeros(30);expected[::3]=rng.uniform(10,100,10)
        multiplier=np.where(expected>0,0.,rng.uniform(.1,2,30))
        rhs=a.T@(a@expected)+3.-multiplier
        target=np.linalg.lstsq(a.T,rhs,rcond=None)[0]
        actual,status=strict_positive_lasso(a,target,3.)
        np.testing.assert_allclose(actual,expected,atol=1e-6,rtol=1e-8)
        self.assertTrue(status['success']);self.assertLessEqual(status['projected_gradient_max'],1e-5)

    def check_pair(self,kinds,gap):
        times=[5000.,5000.+gap];charges=[3000.,1000.]
        y=10.+self.bank.matrix(self.t,times,kinds,charges)@charges
        y+=np.random.default_rng(261041002).normal(0,.1,len(self.t))
        fit=fit_model_selection(self.t,y,self.cal,ModelSelectionConfig(addition_cost=1.))
        self.assertEqual(fit['particle_name'],self.bank.names[kinds].tolist())
        np.testing.assert_allclose(fit['times_ns'],times,atol=2.)
        np.testing.assert_allclose(fit['charge_pe'],charges,rtol=.1)
        self.assertTrue(np.all(np.diff(fit['weighted_cost_history'])<=1e-6))
        self.assertEqual(fit['fit_truth_inputs'],[])

    def test_ten_ns_second_pulse_unknown_order(self):
        for kinds in ([0,1],[1,0]):self.check_pair(kinds,10.)

    def test_two_hundred_ns_gamma_alpha_class_trap(self):
        self.check_pair([1,0],200.)

    def test_single_and_empty_do_not_gain_reported_ghosts(self):
        for kinds in ([],[0],[1]):
            times=[5000.] if kinds else [];charges=[3000.] if kinds else []
            y=10.+self.bank.matrix(self.t,times,kinds,charges)@np.asarray(charges)
            y+=np.random.default_rng(261041003).normal(0,.1,len(self.t))
            fit=fit_model_selection(self.t,y,self.cal,ModelSelectionConfig(addition_cost=1.))
            self.assertEqual(fit['particle_name'],self.bank.names[kinds].tolist())


if __name__=='__main__':unittest.main()
