"""Independent numerical references for shared-lifetime variable projection."""
import unittest
import numpy as np
from scipy.optimize import lsq_linear

from .alpha_gamma_shared_tau import (
    ProfileProblem, fit_problem, fit_joint_shared, fit_shared_prefix, prior_precision,
)
from .alpha_gamma_short_window import PrefixCalibration, PrefixConfig
from dataclasses import asdict

READOUT = dict(rise_ns=2., decay_ns=20., lowpass_ns=2., aperture_ns=.5)
SHARED = dict(taus_ns=[95., 351.], log_tau_covariance=[[.002, -.0003], [-.0003, .003]])


def problem(post=1000., taus=(95., 351.), onset=-6.3, amplitudes=(1000., 2000.), noise=0., seed=8):
    x = np.arange(-20., post+.1, 2.)
    p = ProfileProblem(x, np.zeros(len(x)), READOUT, 13.618501389, 10., .2, .1)
    y = 10.+p.columns(taus, onset)@np.array(amplitudes)
    y += np.random.default_rng(seed).normal(0., noise, len(y))
    return ProfileProblem(x, y, READOUT, 13.618501389, 10., .2, .1)


class SharedTauTests(unittest.TestCase):
    def test_profile_against_augmented_bounded_linear_least_squares(self):
        p = problem(noise=.12)
        X = p.columns((95., 351.), -5.8)
        design = np.vstack([p.whiten(np.c_[X, np.ones(len(X))]), [0., 0., 1/.2]])
        target = np.r_[p.whiten(p.y-10.), 0.]
        independent = lsq_linear(design, target, bounds=([0., 0., -np.inf], [np.inf]*3), tol=1e-12)
        fit = p.profile((95., 351.), -5.8, True)
        np.testing.assert_allclose(fit['amplitudes_pe'], independent.x[:2], rtol=1e-9)
        self.assertAlmostEqual(fit['objective'], np.sum((design@independent.x-target)**2), places=7)

    def test_negative_amplitude_active_set(self):
        p = problem(amplitudes=(2000., -200.))
        result = p.profile((95., 351.), -6.3, True)
        self.assertEqual(result['amplitudes_pe'][1], 0.)

    def test_fixed_shared_identity(self):
        p = problem(post=150., noise=.1)
        a = fit_problem(p, 'fixed_tau', SHARED)
        b = fit_problem(p, 'shared_tau', SHARED)
        self.assertAlmostEqual(a['charge_pe'], b['charge_pe'], places=7)
        self.assertAlmostEqual(a['objective'], b['objective'], places=8)

    def test_free_noiseless_lifetime_and_charge_recovery(self):
        p = problem(taus=(110., 420.))
        result = fit_problem(p, 'free_tau', SHARED)
        np.testing.assert_allclose([result['tau_fast_ns'], result['tau_slow_ns']], [110., 420.], rtol=2e-4)
        self.assertLess(abs(result['charge_pe']/3000.-1), 1e-5)

    def test_prior_inflation_scales_precision(self):
        np.testing.assert_allclose(prior_precision(SHARED, 3.), prior_precision(SHARED, 1.)/9.)
        self.assertTrue(np.linalg.eigvalsh(prior_precision(SHARED, 1.)).min() > 0)

    def test_strong_prior_converges_to_shared_tau(self):
        p = problem(post=100., taus=(110., 420.), noise=.1)
        result = fit_problem(p, 'constrained_tau', SHARED, inflation=.001)
        np.testing.assert_allclose([result['tau_fast_ns'], result['tau_slow_ns']], SHARED['taus_ns'], rtol=1e-4)

    def test_joint_calibration_with_event_local_amplitudes_and_onsets(self):
        ps = [problem(taus=(100., 380.), amplitudes=a, onset=t, noise=.01, seed=i)
              for i, (a, t) in enumerate([((1000., 2000.), -6.), ((2000., 500.), -8.),
                                          ((700., 2800.), -7.), ((3000., 800.), -5.)])]
        fit = fit_joint_shared(ps)
        np.testing.assert_allclose(fit['taus_ns'], [100., 380.], rtol=.015)
        self.assertTrue(np.linalg.eigvalsh(fit['log_tau_covariance']).min() > 0)
        self.assertTrue(all(b['objective'] <= a['objective']+1e-7
                            for a, b in zip(fit['iterations'], fit['iterations'][1:])))

    def test_charge_curvature_matches_independent_fisher(self):
        p = problem(noise=.1)
        fit = fit_problem(p, 'fixed_tau', SHARED)
        a = np.array(fit['amplitudes_pe']); t = fit['onset_relative_ns']
        X = p.columns((95., 351.), t)
        derivative = (p.columns((95., 351.), t+1e-3)-p.columns((95., 351.), t-1e-3))@a/.002
        J = p.whiten(np.c_[X, np.ones(len(X)), derivative])
        info = J.T@J; info[2, 2] += 25.
        c = np.array([1., 1., 0., 0.]); sigma = np.sqrt(c@np.linalg.inv(info)@c)
        self.assertAlmostEqual((fit['charge_interval95_pe'][1]-fit['charge_interval95_pe'][0])/3.92/sigma, 1., places=4)

    def test_causal_prefix_future_nan_and_time_translation(self):
        grid = np.arange(-200., 1250.1, 2.)
        shape = ProfileProblem(grid, np.zeros(len(grid)), READOUT, 13.618501389, 10., .2, .1)
        cols = shape.columns((95., 351.), 0.)
        arrays = dict(grid_ns=grid, covariance_grid_ns=grid,
                      shot_covariance_per_pe=np.eye(len(grid))*1e-5,
                      alpha_means=np.tile(cols@np.array([.3, .7]), (3, 1)),
                      gamma_means=np.tile(cols@np.array([.79, .21]), (3, 1)))
        meta = dict(prefix_config=asdict(PrefixConfig()), readout=READOUT,
                    baseline_mV=10., area_per_pe_mV_ns=13.618501389, electronics_sigma_mV=.1,
                    classes={c: dict(charge_centers_pe=[500., 2000., 10000.]) for c in ('alpha', 'gamma')})
        cal = PrefixCalibration(meta, arrays)
        time = np.arange(0., 2000., 2.)
        curve = ProfileProblem(time, np.zeros(len(time)), READOUT, 13.618501389, 10., .2, .1)
        y = 10.+curve.columns((95., 351.), 500.)@np.array([1000., 2000.])
        fit = fit_shared_prefix(time, y, cal, SHARED, post_ns=150.)
        modified = y.copy(); modified[time > fit['last_observed_ns']] = np.nan
        repeat = fit_shared_prefix(time+2346., modified, cal, SHARED, post_ns=150.)
        for method in fit['methods']:
            self.assertAlmostEqual(fit['methods'][method]['charge_pe'], repeat['methods'][method]['charge_pe'], places=6)


if __name__ == '__main__':
    unittest.main()
