"""Numerical equivalence and rank-deficiency regression checks."""
import unittest
import numpy as np
from scipy.optimize import nnls
from .alpha_gamma_nnls import reduced_nnls


class ReducedSolverChecks(unittest.TestCase):
    def test_tall_correlated_pulses(self):
        rng=np.random.default_rng(261027001)
        t=np.arange(600.)
        a=np.column_stack([np.exp(-np.maximum(t-x,0)/80)*(t>=x) for x in np.arange(20.)*7])
        charge=rng.uniform(0,3000,20);charge[::3]=0
        b=a@charge+rng.normal(0,.1,len(t))
        expected,_=nnls(a,b);actual,residual=reduced_nnls(a,b)
        np.testing.assert_allclose(actual,expected,atol=1e-7,rtol=1e-10)
        self.assertAlmostEqual(residual,np.linalg.norm(a@actual-b),places=10)

    def test_rank_deficient_duplicate_columns(self):
        rng=np.random.default_rng(261027002)
        base=rng.uniform(0,1,(400,16));a=np.c_[base,base[:,0],base[:,1]]
        b=base@rng.uniform(0,10,16)
        expected,_=nnls(a,b);actual,_=reduced_nnls(a,b)
        np.testing.assert_allclose(a@actual,a@expected,atol=1e-8,rtol=1e-10)
        self.assertTrue(np.all(actual>=0))

    def test_negative_signal_has_zero_nonnegative_solution(self):
        a=np.random.default_rng(261027003).uniform(0,1,(400,20))
        actual,_=reduced_nnls(a,-np.ones(400))
        np.testing.assert_array_equal(actual,np.zeros(20))


if __name__=='__main__':unittest.main()
