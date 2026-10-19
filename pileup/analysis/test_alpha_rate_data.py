"""Regressions for photon ancestry collisions in dense pileup records."""
import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np

from pileup.analysis.alpha_rate_data import photon_key, photon_full_key, export_window


class PhotonKeyTests(unittest.TestCase):
    def test_same_time_and_boundary_x_remain_distinct(self):
        # The first pair is taken from the 20 MHz regression; the third also
        # shares y, so energy must distinguish it. Input photon order varies.
        t = np.array([8184989, 8184989, 8184989])
        x = np.array([5., 5., 5.], dtype=np.float32)
        y = np.array([.9336937069892883, 1.0416803359985352, 1.0416803359985352])
        energy = np.array([2.2829649448394775, 2.291916847229004, 2.31])
        self.assertEqual(len(np.unique(photon_key(t, x))), 1)
        keys = photon_full_key(t, x, y, energy)
        self.assertEqual(len(np.unique(keys)), 3)
        permutation = np.array([2, 0, 1])
        order = np.argsort(keys[permutation])
        located = np.searchsorted(keys[permutation][order], keys)
        np.testing.assert_array_equal(permutation[order[located]], np.arange(3))

    def test_long_window_times_do_not_wrap(self):
        keys = photon_full_key([0, 2**32], [1, 1], [1, 1], [2, 2])
        self.assertEqual(len(np.unique(keys)), 2)

    def test_identical_records_remain_ambiguous(self):
        keys = photon_full_key([1, 1], [2, 2], [3, 3], [4, 4])
        self.assertEqual(len(np.unique(keys)), 1)

    def test_mixed_cache_cannot_bypass_alpha_only_policy(self):
        with TemporaryDirectory() as directory:
            output=Path(directory)
            for name in ('metadata','truth'):(output/name).mkdir()
            (output/'metadata/window.json').write_text(json.dumps({'seed':1,'rate_hz':1e6}))
            np.savez(output/'truth/window.npz',particle_name=['alpha','gamma'])
            with self.assertRaisesRegex(ValueError,'Unexpected cached'):
                export_window('window.root',output)
            self.assertEqual(export_window('window.root',output,allowed_particles=('alpha','gamma'))['seed'],1)

    def test_cached_seed_still_checked(self):
        with TemporaryDirectory() as directory:
            output=Path(directory)
            for name in ('metadata','truth'):(output/name).mkdir()
            (output/'metadata/window.json').write_text(json.dumps({'seed':1,'rate_hz':1e6}))
            np.savez(output/'truth/window.npz',particle_name=['alpha'])
            with self.assertRaises(AssertionError):
                export_window('window.root',output,expected={'seed':2,'rate_hz':1e6})


if __name__ == '__main__':
    unittest.main()
