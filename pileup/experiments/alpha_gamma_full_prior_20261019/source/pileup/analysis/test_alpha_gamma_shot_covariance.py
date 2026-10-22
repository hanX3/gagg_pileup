import unittest
import numpy as np
from scipy.linalg import cholesky_banded, solve_banded
from .alpha_gamma_shot_covariance import covariance_bands


class ShotCovarianceChecks(unittest.TestCase):
    def test_banded_whitening_matches_dense_covariance(self):
        t=np.arange(160)*4.
        cal=dict(readout=dict(rise_ns=2.,decay_ns=20.,lowpass_ns=2.,aperture_ns=.5),
            area_per_pe_mV_ns=13.618501389,noise=dict(baseline_variance_mV2=.01))
        intensity=.002+3.*np.exp(-np.maximum(t-40,0)/95.)*(t>=40)
        band=covariance_bands(t,intensity,cal,support_ns=120.)
        dense=np.zeros((len(t),len(t)))
        for j in range(len(t)):
            for lag in range(min(len(band),len(t)-j)):
                dense[j+lag,j]=dense[j,j+lag]=band[lag,j]
        lower=cholesky_banded(band,lower=True)
        transform=solve_banded((len(band)-1,0),lower,np.eye(len(t)))
        np.testing.assert_allclose(transform@dense@transform.T,np.eye(len(t)),atol=1e-9)
        dense_lower=np.linalg.cholesky(dense)
        np.testing.assert_allclose(transform,np.linalg.inv(dense_lower),atol=1e-9)

    def test_quadrature_and_tail_converge(self):
        t=np.arange(180)*4.
        cal=dict(readout=dict(rise_ns=2.,decay_ns=20.,lowpass_ns=2.,aperture_ns=.5),
            area_per_pe_mV_ns=13.618501389,noise=dict(baseline_variance_mV2=.01))
        intensity=np.full(len(t),2.)
        a=covariance_bands(t,intensity,cal,quadrature=12,support_ns=240.)
        b=covariance_bands(t,intensity,cal,quadrature=24,support_ns=320.)
        self.assertLess(np.max(abs(a-b[:len(a)]))/np.max(b),.001)


if __name__=='__main__':unittest.main()
