"""Independent covariance and numerical checks for correlated pile-up fitting."""
import unittest
import numpy as np
from scipy.linalg import cholesky, solve_triangular

from .alpha_correlated_fit import InnovationsMetric, ResidualCalibration, CorrelatedConfig, projected, optimize,fit_correlated
from .alpha_pileup_fit import FitConfig
from .alpha_local_fit_study import settings
from .alpha_full_signal_study import make_template


class CovarianceChecks(unittest.TestCase):
    def test_empty_blind_record(self):
        d=settings();template=make_template(d);t=np.arange(0.,4000.,2.)
        y=10+np.random.default_rng(21).normal(0,.1,len(t))
        fitted=fit_correlated(t,y,template,FitConfig(**d['fit_config']),ResidualCalibration(.1,0.,.8,.98))
        self.assertEqual(fitted['times_ns'],[])
        self.assertEqual(fitted['fit_truth_inputs'],[])
        np.testing.assert_allclose(fitted['model_mV'],10.,atol=1e-12)

    def test_nonstationary_precision_and_determinant(self):
        rng=np.random.default_rng(19);mu=rng.uniform(0,60,23)
        c=ResidualCalibration(.2,.07,.31,.93);w=InnovationsMetric(mu,c)
        covariance=np.outer(w.s,w.s)*c.rho**abs(np.arange(len(mu))[:,None]-np.arange(len(mu)))+c.electronic_sigma_mV**2*np.eye(len(mu))
        lower=cholesky(covariance,lower=True)
        x=rng.normal(size=(len(mu),4))
        reference=solve_triangular(lower,x,lower=True)
        np.testing.assert_allclose(w.apply(x),reference,atol=2e-12,rtol=2e-12)
        self.assertAlmostEqual(w.logdet,np.linalg.slogdet(covariance)[1],places=11)

    def test_diagonal_limit_and_empty_signal(self):
        mu=np.array([0.,1.,4.,50.]);x=np.arange(4.)+.3;c=ResidualCalibration(.1,.05,.2,.98)
        w=InnovationsMetric(mu,c,diagonal=True)
        np.testing.assert_allclose(w.apply(x),x/np.sqrt(.01+.05+.2*mu),atol=1e-12)
        w=InnovationsMetric(np.zeros(4),ResidualCalibration(.1,0.,.8,.99))
        np.testing.assert_allclose(w.apply(x),x/.1,atol=1e-12)

    def test_stationary_whitening_on_independent_draws(self):
        rng=np.random.default_rng(20);n=15;draws=30000
        w=InnovationsMetric(np.linspace(0,20,n),ResidualCalibration(.3,.04,.2,.85))
        latent=np.empty((n,draws));latent[0]=rng.normal(size=draws)
        for i in range(1,n):latent[i]=.85*latent[i-1]+np.sqrt(1-.85**2)*rng.normal(size=draws)
        samples=w.s[:,None]*latent+.3*rng.normal(size=(n,draws))
        actual=np.cov(w.apply(samples))
        np.testing.assert_allclose(actual,np.eye(n),atol=.025)

    def test_variable_projection_gradient(self):
        template=make_template(settings());t=np.arange(0.,1200.,2.)+.137
        times=np.array([150.173,188.617]);true=np.array([151.101,191.221])
        y=template(t[:,None]-true)@np.array([80.,40.])
        w=InnovationsMetric(y,ResidualCalibration(.1,.02,.2,.9))
        cost,g,_=projected(t,y,times,template,w)
        step=1e-5;numeric=[]
        for j in range(2):
            shift=np.zeros(2);shift[j]=step
            numeric.append((projected(t,y,times+shift,template,w)[0]-projected(t,y,times-shift,template,w)[0])/(2*step))
        np.testing.assert_allclose(g,numeric,rtol=5e-6,atol=2e-6)

    def test_known_pair_optimizer_and_time_translation(self):
        template=make_template(settings());t=np.arange(0.,1000.,2.)+.137
        true=np.array([150.173,178.617]);amp=np.array([80.,40.]);y=template(t[:,None]-true)@amp
        w=InnovationsMetric(y,ResidualCalibration(.1,.02,.2,.9));cfg=CorrelatedConfig(max_iterations=80,global_iterations=4)
        result,trials=optimize(t,y,[[148.1,181.1]],template,130.,200.,w,cfg,global_rng=np.random.default_rng(3))
        np.testing.assert_allclose(result['times'],true,atol=.005)
        np.testing.assert_allclose(result['amplitude'],amp,atol=.01)
        moved,_=optimize(t+713.,y,[[861.1,894.1]],template,843.,913.,w,cfg)
        np.testing.assert_allclose(moved['times']-713.,result['times'],atol=.005)
        large,_=optimize(t+1e6,y,[[1e6+148.1,1e6+181.1]],template,1e6+130.,1e6+200.,w,cfg,global_rng=np.random.default_rng(3))
        np.testing.assert_allclose(large['times']-1e6,true,atol=.005)
        self.assertTrue(any(r['kind']=='differential_evolution' for r in trials))


if __name__=='__main__':unittest.main()
