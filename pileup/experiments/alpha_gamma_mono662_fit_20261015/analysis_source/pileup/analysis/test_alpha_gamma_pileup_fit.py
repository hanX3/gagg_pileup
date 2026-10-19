"""Checks of mixed reconstruction limits and class-independent evaluation."""
import json
from pathlib import Path
import unittest
import numpy as np

from pileup.analysis.alpha_gamma_pileup_fit import TemplateBank, MixedFitConfig, fit_waveform, score_mixed


class MixedFittingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root=Path(__file__).resolve().parents[2]
        readout=json.loads((root/'g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json').read_text())['readout']
        cls.bank=TemplateBank(readout);cls.config=MixedFitConfig();cls.t=np.arange(.25,12000.,2.)

    def test_dark_electronics_do_not_create_signal(self):
        y=10.+np.random.default_rng(701).normal(0,.1,len(self.t))
        fit=fit_waveform(self.t,y,self.bank,self.config)
        self.assertEqual(fit['times_ns'],[])

    def test_unknown_number_and_type_on_separated_controls(self):
        times=np.array([3000.,4500.,6500.]);types=[1,0,1];charges=np.array([900.,4000.,2200.])
        y=10.+self.bank.matrix(self.t,times,types,self.config)@charges
        fit=fit_waveform(self.t,y,self.bank,self.config)
        np.testing.assert_allclose(fit['times_ns'],times,atol=.03)
        np.testing.assert_allclose(fit['charge_pe'],charges,rtol=1e-4)
        self.assertEqual(fit['particle_name'],['gamma','alpha','gamma'])

    def test_close_pair_order_is_not_fixed_by_initialization(self):
        for types,charges in (([0,1],[4000.,2000.]),([1,0],[2000.,4000.])):
            times=[5000.,5020.];y=10.+self.bank.matrix(self.t,times,types,self.config)@charges
            fit=fit_waveform(self.t,y,self.bank,self.config)
            self.assertEqual(fit['particle_name'],self.bank.names[types].tolist())
            np.testing.assert_allclose(fit['times_ns'],times,atol=.03)
            np.testing.assert_allclose(fit['charge_pe'],charges,rtol=1e-4)

    def test_matching_does_not_use_particle_labels_or_charge(self):
        truth=dict(first_active_ns=np.array([5000.,5008.]),emission_ns=np.array([5000.,5008.]),
                   charge_pe=np.array([1000.,4000.]),particle_name=np.array(['alpha','gamma']))
        fit=dict(times_ns=[5000.,5008.],charge_pe=[4000.,1000.],particle_name=['gamma','alpha'])
        scored=score_mixed(fit,truth,[4000.,6000.])
        self.assertEqual(scored['matched_pairs'],[(0,0),(1,1)])
        self.assertEqual(scored['correct_type_count'],0)
        self.assertEqual(scored['charge_and_type_recovered_count'],0)
        self.assertEqual(scored['matched_count'],2)

    def test_guard_match_and_zero_signal_gamma_are_not_false_misses(self):
        truth=dict(first_active_ns=np.array([3999.,5000.,np.nan]),emission_ns=np.array([3990.,4990.,5200.]),
                   charge_pe=np.array([1000.,2000.,0.]),particle_name=np.array(['gamma','alpha','gamma']))
        fit=dict(times_ns=[4001.,5001.],charge_pe=[1000.,2000.],particle_name=['gamma','alpha'])
        scored=score_mixed(fit,truth,[4000.,6000.])
        self.assertEqual(scored['target_count'],1)
        self.assertEqual(scored['unmatched_predictions'],0)
        self.assertEqual(scored['by_particle']['gamma']['central_injected'],1)
        self.assertEqual(scored['by_particle']['gamma']['missed_count'],0)


if __name__=='__main__':unittest.main()
