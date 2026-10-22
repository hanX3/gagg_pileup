"""Validate the actual NNLS-work budget and returned feasible best candidate."""
import unittest
from unittest.mock import patch

import numpy as np

from .alpha_correlated_fit import InnovationsMetric, ResidualCalibration, projected
from .alpha_optimizer_budget import search


class Template:
    grid=np.arange(-4.,6000.001,.25)
    shape=np.where(grid>=0,(1-np.exp(-grid/4))*np.exp(-grid/150),0.)
    def __call__(self,t):return np.interp(t,self.grid,self.shape,left=0.,right=0.)


class BudgetTests(unittest.TestCase):
    def test_equal_actual_work_and_feasible_return(self):
        t=np.arange(0.,3000.,12.);template=Template()
        values=template(t[:,None]-np.array([1001.3,1021.7]))@np.array([3.,2.])
        metric=InnovationsMetric(np.zeros(len(t)),ResidualCalibration(1.,0.,1.,.3))
        for population in (False,True):
            with patch('pileup.analysis.alpha_optimizer_budget.projected',wraps=projected) as calls:
                result=search(t,values,template,metric,2,seed=1701,budget=500,population=population)
            self.assertEqual(calls.call_count,500)
            self.assertEqual(result['profile_evaluations'],500)
            self.assertTrue(np.all(np.diff(result['times_ns'])>=3.-1e-5))
            cost,_,amplitude=projected(t,values,result['times_ns'],template,metric)
            self.assertAlmostEqual(cost,result['cost'],places=12)
            np.testing.assert_array_equal(amplitude,result['amplitude'])
            repeat=search(t,values,template,metric,2,seed=1701,budget=500,population=population)
            self.assertEqual(result['times_ns'],repeat['times_ns'])
            self.assertEqual(result['cost'],repeat['cost'])


if __name__=='__main__':unittest.main()
