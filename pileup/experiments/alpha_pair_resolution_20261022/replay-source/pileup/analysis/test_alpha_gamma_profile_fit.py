"""Numerical profile checks that differ from the implementation path."""
import unittest
import numpy as np
from scipy.optimize import nnls
from .alpha_gamma_profile_fit import ProfileSystem


class ProfileChecks(unittest.TestCase):
    def test_subsets_match_original_space_nnls(self):
        rng=np.random.default_rng(109)
        matrix=rng.normal(size=(110,12));matrix[:,1]=.98*matrix[:,0]+.02*matrix[:,1]
        target=matrix[:,[0,4,7]]@np.array([2.,1.,3.])+rng.normal(0,.04,110)
        profile=ProfileSystem(matrix,target)
        for subset in ([0,2,4,6,8,10],[1,3,4,7,9,11],[0,1,4,7]):
            q,cost=profile.solve(np.array(subset))
            direct=nnls(matrix[:,subset],target,maxiter=1000)[0]
            np.testing.assert_allclose(q,direct,atol=1e-8,rtol=1e-8)
            self.assertAlmostEqual(cost,float(np.sum((matrix[:,subset]@direct-target)**2)),places=9)

    def test_rank_deficient_subset_preserves_fitted_signal(self):
        rng=np.random.default_rng(110);base=rng.normal(size=(70,3))
        matrix=np.c_[base[:,0],base[:,0],base[:,1:]];target=base@np.array([2.,3.,1.])
        profile=ProfileSystem(matrix,target);q,cost=profile.solve(np.arange(4))
        np.testing.assert_allclose(matrix@q,target,atol=1e-8)
        self.assertLess(abs(cost),1e-7)


if __name__=='__main__':unittest.main()
