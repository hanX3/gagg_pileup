"""Numerical controls for causal history conditioning, not physical accuracy."""
import copy
import unittest
import numpy as np
from . import test_alpha_gamma_profile_variants as fixtures
from .alpha_gamma_context_calibration import preceding_voltage
from .alpha_gamma_context_fit import shape_parameters, fit_context
from .alpha_gamma_model_selection_fit import fit_model_selection, ModelSelectionConfig
from .alpha_gamma_shape_profile import fit_shape_profile, ShapeConfig


class ContextChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixtures.ProfileVariantChecks.setUpClass()
        cls.cal=copy.deepcopy(fixtures.ProfileVariantChecks.cal)
        cls.t=fixtures.ProfileVariantChecks.t
        cls.bank=fixtures.ProfileVariantChecks.bank
        cls.cal['context_classes']={p:dict(v,slow_per_100mV=.03,history_range_mV=[0.,200.],
            fraction_variance_coeff=[.0001,.0004]) for p,v in cls.cal['classes'].items()}

    def test_history_is_causal_and_clamped(self):
        y=np.full(len(self.t),30.)
        before=preceding_voltage(self.t,y,5000.)
        y[self.t>=4996.]=1e9
        self.assertEqual(before,20.)
        self.assertEqual(preceding_voltage(self.t,y,5000.),before)
        centers,sigmas,history=shape_parameters(self.t,y,np.array([5000.]),np.array([1000.]),self.cal)
        np.testing.assert_allclose(centers[0],[.687,.249])
        np.testing.assert_allclose(sigmas,np.sqrt(.0005))
        self.assertEqual(history.tolist(),[20.])

    def test_default_and_explicit_parameters_are_identical(self):
        q=[3000.,1000.];times=[5000.,5100.];kinds=[0,1]
        y=10.+self.bank.matrix(self.t,times,kinds,q)@q
        initial=fit_model_selection(self.t,y,self.cal,ModelSelectionConfig(addition_cost=1.))
        c=initial['all_components'];n=len(c['times_ns'])
        from .alpha_gamma_shape_profile import fraction_center
        centers=np.array([[fraction_center(self.cal,k,q) for k in (0,1)] for q in c['shape_charge_pe']])
        sigmas=np.broadcast_to([self.cal['classes'][p]['rms_fraction_residual'] for p in ('alpha','gamma')],(n,2))
        a=fit_shape_profile(self.t,y,self.cal,initial=initial)
        b=fit_shape_profile(self.t,y,self.cal,initial=initial,center_grid=centers,sigma_grid=sigmas)
        np.testing.assert_array_equal(a['model_mV'],b['model_mV'])
        for key in ('charge_pe','particle_name','times_ns','objective_history'):self.assertEqual(a[key],b[key])
        fitted=fit_context(self.t,y,self.cal,ShapeConfig(sigma_scale=.5),initial=initial)
        self.assertEqual(fitted['fit_truth_inputs'],[])
        self.assertTrue(np.all(np.diff(fitted['objective_history'])<=1e-7))
        self.assertEqual(fitted['particle_name'],['alpha','gamma'])
        np.testing.assert_allclose(fitted['charge_pe'],q,rtol=.1)


if __name__=='__main__':unittest.main()
