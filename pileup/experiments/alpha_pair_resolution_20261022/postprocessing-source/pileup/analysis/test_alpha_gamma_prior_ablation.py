"""Numerical controls for a fixed-candidate diagnostic, not blind accuracy."""
from itertools import product
from pathlib import Path
import json
import unittest
import numpy as np
from scipy.optimize import nnls
from .alpha_gamma_prior_ablation import prepare, fit_fixed, decode_cached, PriorConfig, ARMS
from .alpha_gamma_regularized_fit import CalibratedBank


class PriorAblationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        cls.cal = json.loads((root/'pileup/experiments/alpha_gamma_ml_20261016/calibration.json').read_text())
        cls.t = np.arange(2048)*2.
        cls.bank = CalibratedBank(cls.cal)

    def common(self, count=2):
        times = [1500., 1540.][:count]
        q = np.array([1200., 800.][:count])
        y = self.bank.matrix(self.t, times, [0, 1][:count], np.ones(count)*1000.) @ q+self.cal['baseline_mV']
        fitted = dict(times_ns=times, particle_name=['gamma', 'alpha'][:count], charge_pe=[1000.]*count)
        neural = dict(all_components=fitted)
        dense = dict(event_score=np.tile([[.7], [.3]], (1, len(self.t))))
        baseline = dict(all_components=fitted, model_mV=y)
        return prepare(self.t, y, self.cal, neural, dense, baseline)

    def test_two_candidate_search_matches_exhaustive_discrete_minimum(self):
        common = self.common()
        original = common['R'].copy()
        for name, config in ARMS.items():
            with self.subTest(arm=name):
                fitted = fit_fixed(common, config)
                weights = config.charge_strength*common['unit_weight']
                cc = config.class_strength*common['unit_class']
                objectives = []
                for kinds in product((0, 1), repeat=2):
                    H = common['R'][:, 2*np.arange(2)+kinds]
                    a = np.vstack([H, np.diag(np.sqrt(weights))])
                    b = np.r_[common['target'], np.sqrt(weights)*common['prior']]
                    q = nnls(a, b)[0]
                    objectives.append(float(np.linalg.norm(a @ q-b)**2+cc[np.arange(2), kinds].sum()+common['constant']))
                self.assertAlmostEqual(fitted['objective_history'][-1], min(objectives), places=6)
                self.assertEqual(fitted['all_components']['times_ns'], common['times'].tolist())
                self.assertEqual(fitted['fit_truth_inputs'], [])
                self.assertTrue(np.all(np.diff(fitted['objective_history']) <= 1e-6))
        np.testing.assert_array_equal(common['R'], original)

    def test_prior_switches_change_only_the_declared_penalties(self):
        common = self.common(count=1)
        outputs = {name: fit_fixed(common, config) for name, config in ARMS.items()}
        self.assertEqual(len({r['candidate_sha256'] for r in outputs.values()}), 1)
        for name in ('none', 'class'):
            self.assertEqual(outputs[name]['charge_penalty'], 0.)
        for name in ('none', 'charge'):
            self.assertEqual(outputs[name]['class_penalty'], 0.)
        self.assertGreater(outputs['class']['class_penalty'], 0.)
        for r in outputs.values():
            self.assertEqual(r['all_components']['shape_charge_pe'], [1000.])
        self.assertNotEqual(outputs['charge']['all_components']['charge_pe'], outputs['none']['all_components']['charge_pe'])

    def test_empty_candidates_and_negative_strength(self):
        y = np.full_like(self.t, self.cal['baseline_mV'])
        empty = dict(times_ns=[], particle_name=[], charge_pe=[])
        common = prepare(self.t, y, self.cal, dict(all_components=empty),
                         dict(event_score=np.ones((2, len(y)))*.01), dict(all_components=empty, model_mV=y))
        fitted = fit_fixed(common)
        self.assertEqual(fitted['times_ns'], [])
        np.testing.assert_array_equal(fitted['model_mV'], y)
        with self.assertRaises(ValueError):
            fit_fixed(common, PriorConfig(class_strength=-1.))

    def test_cached_decoder_preserves_offsets_classes_and_thresholds(self):
        dense = dict(event_score=np.zeros((2, len(self.t))), charge_pe=np.full((2, len(self.t)), 1000.),
                     time_offset_ns=np.zeros((2, len(self.t))))
        dense['event_score'][0, 200] = .8
        dense['event_score'][1, 250] = .7
        dense['event_score'][1, 300] = .3
        dense['time_offset_ns'][0, 200] = .6
        decoded = decode_cached(self.t, dense, [.5, .5])
        self.assertEqual(decoded['times_ns'], [400.6, 500.])
        self.assertEqual(decoded['particle_name'], ['alpha', 'gamma'])
        with self.assertRaises(ValueError):
            prepare(self.t*2, np.zeros_like(self.t), self.cal, {}, {}, {})


if __name__ == '__main__':
    unittest.main()
