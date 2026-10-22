"""Independent numerical references for the new conditional-prefix models."""
from dataclasses import asdict
import math
import unittest

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import multivariate_normal

from .alpha_gamma_short_window import PrefixConfig
from .alpha_gamma_waveform_fit import sampled_component
from .alpha_gamma_advanced_fit import (
    AdvancedCalibration, AdvancedConfig, bootstrap_means, curvature_interval,
    fit_advanced_prefix, gaussian_profile, lowrank_covariance, project_covariance,
)


def synthetic_calibration(shot_variance=0., shrinkage=.5):
    readout = dict(rise_ns=2., decay_ns=20., lowpass_ns=2., aperture_ns=.5)
    grid = np.arange(-200., 3000.1, 2.)
    covariance_grid = grid[grid <= 1250.]
    arrays = {'grid_ns': grid, 'covariance_grid_ns': covariance_grid,
              'shot_covariance_per_pe': np.eye(len(covariance_grid))*shot_variance}
    classes = {}
    centers = [500., 2000., 10000.]
    for name, ws in (('alpha', .7), ('gamma', .21)):
        means = []
        for q in centers:
            fraction = ws + .02*math.log(q/2000.)
            means.append(13.618501389*((1-fraction)*sampled_component(grid, 0., 95., readout)
                                      + fraction*sampled_component(grid, 0., 351., readout)))
        arrays[name+'_means'] = np.array(means)
        arrays[name+'_bootstrap_deviations'] = np.zeros((3, 32, len(grid)))
        classes[name] = {'charge_centers_pe': centers, 'charge_range_pe': [250., 15000.]}
    for rank in (8, 16):
        c, modes, diagonal = lowrank_covariance(arrays['shot_covariance_per_pe'], rank)
        arrays[f'pca{rank}_shot_covariance'] = c
        arrays[f'pca{rank}_modes'] = modes
        arrays[f'pca{rank}_diagonal'] = diagonal
    meta = {'prefix_config': asdict(PrefixConfig(covariance_shrinkage=shrinkage)),
            'readout': readout, 'classes': classes, 'area_per_pe_mV_ns': 13.618501389,
            'baseline_mV': 10., 'electronics_sigma_mV': .1}
    return AdvancedCalibration(meta, arrays)


class AdvancedFitTests(unittest.TestCase):
    def test_covariance_projection_against_explicit_matrix(self):
        rng = np.random.default_rng(262020601)
        grid = np.arange(-4., 6.1, 2.)
        points = np.array([-5., -3.3, -1.1, 0., 4.9, 6., 8.])
        matrix = rng.normal(size=(len(grid), 9))
        covariance = matrix@matrix.T
        # Independent construction using interpolated basis vectors.
        projection = np.array([np.interp(points, grid, basis, left=0., right=0.)
                               for basis in np.eye(len(grid))]).T
        actual = project_covariance(covariance, grid, points)
        np.testing.assert_allclose(actual, projection@covariance@projection.T, atol=1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(actual)[0], -1e-12)

    def test_white_variance_does_not_depend_on_subsample_onset(self):
        cal = synthetic_calibration(shot_variance=1e-4, shrinkage=1.)
        points = np.arange(2., 50., 2.)
        a = cal.conditional_covariance(points, 1000., 'alpha', 'conditional_mle')
        b = cal.conditional_covariance(points+.93, 1000., 'alpha', 'conditional_mle')
        np.testing.assert_allclose(a, b, atol=1e-14)
        np.testing.assert_allclose(a, np.eye(len(points))*.11, atol=1e-14)

    def test_profiled_gaussian_against_independent_density_and_scalar_optimizer(self):
        rng = np.random.default_rng(262020602)
        matrix = rng.normal(size=(11, 11)); covariance = matrix@matrix.T + np.eye(11)*.2
        mean = rng.normal(size=11); values = 10.3+mean+rng.normal(size=11)
        def direct(b):
            return (-2*multivariate_normal.logpdf(values, mean=mean+b, cov=covariance)
                    +(b-10.)**2/.2**2-len(mean)*math.log(2*math.pi))
        expected = minimize_scalar(direct, bounds=(9., 11.), method='bounded',
                                   options={'xatol': 1e-12})
        actual = gaussian_profile(values, mean, covariance, 10., .2)
        self.assertAlmostEqual(actual['baseline_mV'], expected.x, places=7)
        self.assertAlmostEqual(actual['objective'], expected.fun, places=9)

    def test_logdet_gives_correct_variance_mle(self):
        values = np.array([1., -2., 4., -3.])
        def objective(log_variance):
            return gaussian_profile(values, np.zeros(4), np.eye(4)*math.exp(log_variance),
                                    0., 1e-10)['objective']
        result = minimize_scalar(objective, bounds=(-5., 5.), method='bounded',
                                 options={'xatol': 1e-12})
        self.assertAlmostEqual(math.exp(result.x), float(np.mean(values**2)), places=5)

    def test_lowrank_preserves_variances_and_psd(self):
        rng = np.random.default_rng(262020603)
        M = rng.normal(size=(27, 11)); covariance = M@M.T+.4*np.eye(27)
        approx, modes, diagonal = lowrank_covariance(covariance, 4)
        np.testing.assert_allclose(np.diag(approx), np.diag(covariance), atol=1e-12)
        np.testing.assert_allclose(approx, modes@modes.T+np.diag(diagonal), atol=1e-12)
        self.assertGreater(np.linalg.eigvalsh(approx)[0], 0.)

    def test_bootstrap_mean_variance_against_analytic_sampling_variance(self):
        curves = np.array([1., 2., 4., 9.])[:, None]*np.array([[1., 2.]])
        deviations = bootstrap_means(curves, np.ones(4), 12000,
                                     np.random.default_rng(262020604))
        actual = deviations.T@deviations/(len(deviations)-1)
        expected = np.cov(curves, rowvar=False, ddof=1)/len(curves)
        np.testing.assert_allclose(actual, expected, rtol=.025)

    def test_observed_curvature_scaling_against_exact_quadratic(self):
        information = np.array([[12., 1.], [1., 3.]])
        optimum = np.array([math.log(3000.), -3.])
        objective = lambda x: float((x-optimum)@information@(x-optimum))
        sigma, interval = curvature_interval(objective, optimum, [(0., 12.), (-20., 2.)], AdvancedConfig())
        expected_log_sigma = math.sqrt(np.linalg.inv(information)[0, 0])
        self.assertAlmostEqual(sigma, 3000.*expected_log_sigma, places=5)
        self.assertAlmostEqual(interval[0], 3000.*math.exp(-1.95996398454*expected_log_sigma), places=5)

    def test_blind_noiseless_Q_dependent_dual_shape_recovery(self):
        cal = synthetic_calibration()
        t = np.arange(.25, 6200., 2.)
        q, onset = 3000., 5000.37
        y = 10.+q*cal.template(t-onset, 'gamma', q)
        methods = ('conditional_mle', 'pca8_mle', 'pca16_mle', 'bootstrap_mle', 'pca8_bootstrap_mle')
        fitted = fit_advanced_prefix(t, y, cal, pre_ns=20., post_ns=150., methods=methods)
        self.assertTrue(fitted['detected'])
        for method, result in fitted['methods'].items():
            self.assertEqual(result['particle'], 'gamma', method)
            self.assertLess(abs(result['charge_pe']/q-1), .001, method)
            self.assertLess(abs(fitted['trigger']['time_ns']+result['onset_relative_ns']-onset), .05, method)

    def test_future_nan_and_global_time_shift(self):
        cal = synthetic_calibration()
        t = np.arange(.25, 6200., 2.)
        y = 10.+2000.*cal.template(t-5000.37, 'alpha', 2000.)
        methods = ('conditional_mle', 'bootstrap_mle')
        a = fit_advanced_prefix(t, y, cal, methods=methods, pre_ns=20., post_ns=100.)
        changed = y.copy(); changed[t > a['last_observed_ns']] = np.nan
        b = fit_advanced_prefix(t, changed, cal, methods=methods, pre_ns=20., post_ns=100.)
        self.assertEqual(a, b)
        shifted = fit_advanced_prefix(t+1234., y, cal, methods=methods, pre_ns=20., post_ns=100.)
        for method in methods:
            self.assertAlmostEqual(a['methods'][method]['charge_pe'], shifted['methods'][method]['charge_pe'], places=4)
            self.assertAlmostEqual(a['methods'][method]['onset_relative_ns'], shifted['methods'][method]['onset_relative_ns'], places=5)


if __name__ == '__main__':
    unittest.main()
