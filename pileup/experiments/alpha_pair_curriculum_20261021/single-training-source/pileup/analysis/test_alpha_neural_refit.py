"""Physical checks for neural proposal guards and covariance-preserving refit."""
from dataclasses import replace
import json
from pathlib import Path
import unittest

import numpy as np

from .alpha_correlated_fit import CorrelatedConfig, InnovationsMetric, ResidualCalibration
from .alpha_full_signal_study import make_template
from .alpha_neural_refit import NeuralProposalConfig, fit_neural_wls, propose_neural_times
from .alpha_pileup_fit import FitConfig, score


class NeuralRefitChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        settings = json.loads((root / 'jupyter/artifacts/alpha_correlated_global_20261020/baseline-settings.json').read_text())
        cls.template = make_template(settings)
        cls.base = FitConfig(**settings['fit_config'])
        cls.t = np.arange(0., 4000., 2.)
        cls.calibration = ResidualCalibration(.1, 0., .7976401192, .9819124013)
        cls.physical = CorrelatedConfig(diagonal=True, passes=1, residual_budget=4,
                                      penalty_per_pe=.0001082891147, penalty_multiplier=.5)

    def dense(self):
        return dict(probability=np.zeros(len(self.t)), charge_pe=np.full(len(self.t), 3000.),
                    time_offset_ns=np.zeros(len(self.t)))

    def test_record_edge_outlier_is_not_a_physical_proposal(self):
        d = self.dense()
        d['probability'][[-10, 600]] = .99
        d['charge_pe'][-10] = 1.e9
        times, _, info = propose_neural_times(self.t, d)
        np.testing.assert_array_equal(times, [1200.])
        self.assertEqual(info['excluded_edge_peaks'], 1)

    def test_neural_charge_is_not_an_amplitude_prior_and_covariance_is_unchanged(self):
        times = np.array([1200., 1240.])
        charge = np.array([4000., 2000.])
        scale = self.base.area_per_pe_mV_ns * self.template.unit_area_peak
        mu = self.template(self.t[:, None] - times) @ (charge * scale)
        # A one-candidate observed pilot deliberately lacks the second event.
        pilot_mu = self.template(self.t - times[0]) * (6000. * scale)
        pilot = dict(status='fitted', times_ns=[1200.], peak_mV=[6000. * scale],
                     charge_pe=[6000.], model_mV=10. + pilot_mu)
        d = self.dense()
        d['probability'][[600, 620]] = .99
        d['charge_pe'][[600, 620]] = [1.e8, 1.e8]
        fit = fit_neural_wls(self.t, 10. + mu, d, self.template, self.base,
                             self.calibration, self.physical, pilot)
        truth = dict(first_active_ns=times, emission_ns=times, charge_pe=charge)
        self.assertEqual(score(fit, truth, [500., 3000.])['recovered_count'], 2)
        expected = InnovationsMetric(pilot_mu, self.calibration, diagonal=True)
        self.assertAlmostEqual(fit['logdet_covariance'], expected.logdet, places=8)
        self.assertFalse(fit['neural_charge_prior'])
        self.assertEqual(fit['neural_proposals']['added_neural_components'], 1)

    def test_empty_record_has_no_hallucinated_events(self):
        pilot = dict(status='fitted', times_ns=[], peak_mV=[], charge_pe=[],
                     model_mV=np.full(len(self.t), 10.))
        fit = fit_neural_wls(self.t, np.full(len(self.t), 10.), self.dense(), self.template,
                             self.base, self.calibration, self.physical, pilot)
        self.assertEqual(fit['times_ns'], [])
        np.testing.assert_array_equal(fit['model_mV'], pilot['model_mV'])


if __name__ == '__main__':
    unittest.main()
