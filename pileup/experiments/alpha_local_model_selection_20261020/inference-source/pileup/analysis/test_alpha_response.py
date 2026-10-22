"""Checks for timing boundaries, normalization, calibration units and test isolation."""
from pathlib import Path
import unittest

import numpy as np

from pileup.analysis.alpha_response import (ResponseConfig, fit_charge, histogram_ps,
                                          invert_charge, invert_yield, make_split,
                                          rebin_counts)


class AlphaResponseChecks(unittest.TestCase):
    def test_integer_time_boundaries_and_conservation(self):
        actual = histogram_ps(np.array([0, 999, 1000, 1999, 4999], dtype=np.uint64), 5)
        np.testing.assert_array_equal(actual, [2, 2, 0, 0, 1])
        self.assertEqual(int(rebin_counts(actual, 5)[0]), 5)
        with self.assertRaises(ValueError):
            histogram_ps([5000], 5)  # The record endpoint is excluded.
        with self.assertRaises(ValueError):
            histogram_ps([-1], 5)
        np.testing.assert_array_equal(histogram_ps([], 5), np.zeros(5))

    def test_monotone_energy_inverse_and_range(self):
        e = np.array([.01, .1, .5, 1, 2, 4])
        q = 6000*e + 500*e*e
        coefficient = fit_charge(e, q)
        np.testing.assert_allclose(coefficient, [6000, 500], rtol=1e-10)
        np.testing.assert_allclose(invert_charge(q, coefficient), e, rtol=1e-10)
        self.assertTrue(np.isnan(invert_charge([1e8], coefficient, [.01, 4])[0]))
        np.testing.assert_allclose(invert_charge([0, 4], [0, 1]), [0, 2])
        np.testing.assert_allclose(invert_charge([0, 4], [2, 0]), [0, 2])

    def test_frozen_yield_inversion_does_not_clip_outliers(self):
        e = np.array([0, 1, 2, 3])
        light = np.array([0, 100, 400, 900])
        got = invert_yield([0, 50, 200, 450, 500], .5, e, light, [1, 3])
        np.testing.assert_allclose(got, [np.nan, 1, 2, 3, np.nan], equal_nan=True)

    def test_whole_window_split_is_repeatable_and_disjoint(self):
        project = Path(__file__).resolve().parents[2]
        if not (project / "jupyter/run_records/root_catalog.json").exists():
            self.skipTest("Local production catalog is not distributed in Git.")
        a, b = make_split(project, ResponseConfig()), make_split(project, ResponseConfig())
        self.assertEqual(a, b)
        self.assertEqual(len({r['seed'] for r in a}), 400)
        self.assertEqual(len({r['root_file'] for r in a}), 400)
        for rate in (100000, 1000000, 10000000, 100000000):
            for role, count in [('calibration',30),('validation',20),('test',50)]:
                self.assertEqual(sum(r['rate_hz']==rate and r['role']==role for r in a),count)


if __name__ == '__main__':
    unittest.main()
