"""Numerical checks for local multiplicity searches (software controls only)."""
import unittest
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from .alpha_local_fit import LocalConfig, fit_local, projected_value_gradient
from .alpha_pileup_fit import AlphaTemplate, FitConfig, fit_waveform, score


class LocalFitChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        project = Path(__file__).resolve().parents[2]
        d = json.loads((project/'jupyter/artifacts/alpha_rate_fit_20261013_211541/analysis_config.json').read_text())
        p = d['template']; cls.template = AlphaTemplate(p['readout'], p['slow_fraction'], p['tau_fast_ns'], p['tau_slow_ns'])
        cls.baseline = FitConfig(**d['fit_config'])
        cls.t = np.arange(0., 4000., 2.)+.137

    def waveform(self, times, charges, seed=1):
        p = self.baseline.area_per_pe_mV_ns*self.template.unit_area_peak
        return 10+self.template(self.t[:, None]-np.asarray(times)[None, :])@(np.asarray(charges)*p)+np.random.default_rng(seed).normal(0,.1,len(self.t))

    def test_projected_gradient(self):
        times = np.array([1000.303, 1022.714])
        target = self.waveform([1000.4, 1020.9], [4000.,2000.])-10
        _, analytic, _ = projected_value_gradient(self.t,target,times,self.template)
        finite = []
        for j in range(2):
            plus=times.copy();minus=times.copy();plus[j]+=.00001;minus[j]-=.00001
            finite.append((projected_value_gradient(self.t,target,plus,self.template)[0]-projected_value_gradient(self.t,target,minus,self.template)[0])/.00002)
        np.testing.assert_allclose(analytic, finite, rtol=2e-5, atol=1e-4)

    def test_empty_and_single_do_not_split(self):
        for times, charges in [([],[]),([1000.3],[4000.])]:
            y=self.waveform(times,charges)
            fit=fit_local(self.t,y,self.template,self.baseline)
            self.assertEqual(len(fit['times_ns']),len(times))

    def test_hidden_double(self):
        times=np.array([1000.3,1012.3]);charges=np.array([4000.,2000.])
        y=self.waveform(times,charges)
        fitted=fit_local(self.t,y,self.template,self.baseline)
        result=score(fitted,dict(first_active_ns=times,emission_ns=times,charge_pe=charges),[500.,3000.])
        self.assertEqual(result['recovered_count'],2)
        self.assertEqual(result['unmatched_predictions'],0)

    def test_warm_start_and_model_identity(self):
        y=self.waveform([1000.3,1080.7],[4000.,2000.])
        initial=fit_waveform(self.t,y,self.template,self.baseline)
        warm=fit_local(self.t,y,self.template,self.baseline,initial=initial)
        cold=fit_local(self.t,y,self.template,self.baseline)
        np.testing.assert_allclose(warm['times_ns'],cold['times_ns'],atol=1e-8)
        np.testing.assert_allclose(warm['charge_pe'],cold['charge_pe'],atol=1e-8)
        c=warm['all_components']
        reconstructed=10+self.template(self.t[:,None]-np.array(c['times_ns'])[None,:])@np.array(c['peak_mV'])
        np.testing.assert_allclose(reconstructed,warm['model_mV'],atol=1e-10)

    def test_local_updates_above_old_joint_limit(self):
        # Conditional numerical fixture, not a blind event-recovery benchmark.
        times=500.3+np.arange(70)*40.
        charges=np.full(70,1000.)
        y=self.waveform(times,charges)
        peak=charges*self.baseline.area_per_pe_mV_ns*self.template.unit_area_peak
        initial=dict(status='fitted',times_ns=(times+.6).tolist(),peak_mV=peak.tolist())
        before=np.mean((y-10-self.template(self.t[:,None]-times[None,:]-.6)@peak)**2)
        cfg=replace(LocalConfig(),passes=1,max_additions=0,residual_budget=0)
        fitted=fit_local(self.t,y,self.template,self.baseline,cfg,initial=initial)
        self.assertGreater(fitted['trials']['same_k'],64)
        self.assertLess(fitted['rms_mV']**2,before)


if __name__=='__main__':
    unittest.main()
