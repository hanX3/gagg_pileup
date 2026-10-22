"""End-to-end numerical controls; these do not measure native SiPM accuracy."""
import json
from pathlib import Path
import unittest
import numpy as np
from .alpha_gamma_regularized_fit import CalibratedBank
from .alpha_gamma_model_selection_fit import fit_model_selection, ModelSelectionConfig
from .alpha_gamma_profile_fit import fit_profile, ProfileConfig
from .alpha_gamma_stable_profile import fit_stable, StableConfig
from .alpha_gamma_shape_profile import fit_shape_profile, ShapeConfig, shape_matrix, fraction_center


class ProfileVariantChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root=Path(__file__).resolve().parents[2]
        readout=json.loads((root/'g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json').read_text())['readout']
        cls.cal=dict(readout=readout,area_per_pe_mV_ns=13.618501389,baseline_mV=10.,
            classes={p:dict(slow_intercept=f,slow_log_slope=s,delay_ns=d,charge_range_pe=[250.,20000.],rms_fraction_residual=r)
                for p,f,s,d,r in (('alpha',.681,.0324,-1.91,.0202),('gamma',.243,.0047,-.397,.0487))},
            noise=dict(baseline_variance_mV2=.01,variance_per_mV=.614,shot_rho_at_2ns=.9))
        cls.bank=CalibratedBank(cls.cal);cls.t=np.arange(.25,12000.,2.)

    def test_profile_and_stable_preserve_blind_controlled_pair(self):
        times=[5000.,5100.];kinds=[0,1];charges=[3000.,1000.]
        y=10.+self.bank.matrix(self.t,times,kinds,charges)@charges
        y+=np.random.default_rng(261051001).normal(0,.1,len(self.t))
        initial=fit_model_selection(self.t,y,self.cal,ModelSelectionConfig(addition_cost=1.))
        methods=[(fit_profile,ProfileConfig(passes=2)),(fit_stable,StableConfig(passes=2))]
        for fitter,config in methods:
            with self.subTest(method=fitter.__name__):
                fit=fitter(self.t,y,self.cal,config,initial=initial)
                self.assertEqual(fit['particle_name'],['alpha','gamma'])
                np.testing.assert_allclose(fit['charge_pe'],charges,rtol=.1)
                np.testing.assert_allclose(fit['times_ns'],times,atol=2.)
                history=fit.get('objective_history',fit.get('weighted_cost_history'))
                self.assertTrue(np.all(np.diff(history)<=1e-5))
                self.assertEqual(fit['fit_truth_inputs'],[])

    def test_shape_variation_preserves_charge_and_finite_class_bounds(self):
        for kind,offset in ((0,.012),(1,-.025)):
            with self.subTest(particle=kind):
                q=3000.;f=fraction_center(self.cal,kind,q)+offset
                y=10.+shape_matrix(self.bank,self.t,[5000.],[kind],[f])[:,0]*q
                y+=np.random.default_rng(261051010+kind).normal(0,.1,len(self.t))
                fit=fit_shape_profile(self.t,y,self.cal,ShapeConfig())
                self.assertEqual(fit['particle_name'],[self.bank.names[kind]])
                np.testing.assert_allclose(fit['charge_pe'],[q],rtol=.1)
                self.assertTrue(np.all(np.diff(fit['objective_history'])<=1e-5))
                c=fit['all_components'];ks=[0 if p=='alpha' else 1 for p in c['particle_name']]
                rebuilt=10.+shape_matrix(self.bank,self.t,c['times_ns'],ks,c['slow_fraction'])@c['charge_pe']
                np.testing.assert_allclose(rebuilt,fit['model_mV'],atol=1e-10)
                for k,s,v in zip(ks,c['shape_charge_pe'],c['slow_fraction']):
                    center=fraction_center(self.cal,k,s);r=self.cal['classes'][self.bank.names[k]]['rms_fraction_residual']
                    self.assertLessEqual(abs(v-center),3*r+1e-12)

    def test_empty_voltage_does_not_gain_reported_events(self):
        y=10.+np.random.default_rng(261051020).normal(0,.1,len(self.t))
        initial=fit_model_selection(self.t,y,self.cal,ModelSelectionConfig(addition_cost=1.))
        for fitter in (fit_profile,fit_stable,fit_shape_profile):
            with self.subTest(method=fitter.__name__):
                fit=fitter(self.t,y,self.cal,initial=initial)
                self.assertEqual(fit['times_ns'],[])
                self.assertEqual(fit['fit_truth_inputs'],[])


if __name__=='__main__':unittest.main()
