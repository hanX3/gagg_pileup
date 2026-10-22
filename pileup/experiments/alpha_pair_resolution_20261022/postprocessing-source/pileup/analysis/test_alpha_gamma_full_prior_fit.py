"""Numerical controls for separated penalties in the full hybrid optimizer."""
from dataclasses import replace
from pathlib import Path
import json
import unittest
import numpy as np
from .alpha_gamma_full_prior_fit import FullPriorConfig, fit_full_prior
from .alpha_gamma_ml_hybrid import MixedHybridConfig, fit_hybrid
from .alpha_gamma_regularized_fit import CalibratedBank


class FullPriorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        cls.cal = json.loads((root/'pileup/experiments/alpha_gamma_ml_20261016/calibration.json').read_text())
        cls.t = np.arange(2048)*2.
        cls.bank = CalibratedBank(cls.cal)

    def inputs(self, missing=False):
        truth_times, charges = [1500.6, 1750.4], [1200., 800.]
        voltage = self.bank.matrix(self.t, truth_times, [0, 1], charges) @ charges+self.cal['baseline_mV']
        count = 1 if missing else 2
        candidates = dict(times_ns=[1497., 1754.][:count],
                          particle_name=['gamma', 'alpha'][:count], charge_pe=charges[:count])
        neural = dict(all_components=candidates)
        baseline = dict(all_components=candidates, model_mV=voltage)
        dense = dict(event_score=np.tile([[.5], [.5]], (1, len(self.t))))
        return voltage, neural, dense, baseline

    def test_equal_strengths_reproduce_both_original_endpoints(self):
        voltage, neural, dense, baseline = self.inputs()
        for strength in (0., 10.):
            with self.subTest(strength=strength):
                old = fit_hybrid(self.t, voltage, self.cal, neural, dense, baseline,
                                 MixedHybridConfig(strength=strength, passes=1, max_additions=1))
                new = fit_full_prior(self.t, voltage, self.cal, neural, dense, baseline,
                                    FullPriorConfig(charge_strength=strength, class_strength=strength,
                                                    passes=1, max_additions=1))
                for key in ('times_ns', 'particle_name', 'charge_pe', 'all_components',
                            'weighted_cost_history', 'pair_trials', 'residual_additions'):
                    self.assertEqual(new[key], old[key], key)
                np.testing.assert_array_equal(new['model_mV'], old['model_mV'])

    def test_arrival_search_and_residual_addition_remain_active(self):
        voltage, neural, dense, baseline = self.inputs(missing=True)
        result = fit_full_prior(self.t, voltage, self.cal, neural, dense, baseline,
                                FullPriorConfig(passes=2, max_additions=2))
        self.assertAlmostEqual(result['times_ns'][0], 1500.6, delta=1.)
        self.assertGreater(len(result['residual_additions']), 0)
        self.assertTrue(any(abs(t-1750.4) < 10. for t in result['times_ns']))
        self.assertEqual(result['fit_truth_inputs'], [])
        self.assertEqual(result['all_components']['shape_charge_pe'][0], 1200.)

    def test_disabled_charge_penalty_and_invalid_alignment(self):
        voltage, neural, dense, baseline = self.inputs()
        cfg = FullPriorConfig(passes=1, max_additions=0)
        result = fit_full_prior(self.t, voltage, self.cal, neural, dense, baseline, cfg)
        self.assertEqual(result['config']['charge_strength'], 0.)
        self.assertEqual(result['config']['class_strength'], 10.)
        self.assertTrue(np.all(np.diff(result['weighted_cost_history']) <= 1e-7))
        for bad in (replace(cfg, charge_strength=-1.), replace(cfg, class_strength=float('nan'))):
            with self.assertRaises(ValueError):
                fit_full_prior(self.t, voltage, self.cal, neural, dense, baseline, bad)
        with self.assertRaises(ValueError):
            fit_full_prior(self.t*2., voltage, self.cal, neural, dense, baseline, cfg)


if __name__ == '__main__':
    unittest.main()
