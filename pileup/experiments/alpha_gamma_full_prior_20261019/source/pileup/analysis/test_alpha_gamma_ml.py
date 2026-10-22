"""Checks for class labels, inference seams, missed-class scores and hybrid swaps."""
from pathlib import Path
import json
import unittest
import numpy as np
import torch
from torch import nn
from .alpha_gamma_ml import (NetworkConfig,MixedEventNet,supervised_targets,
    loss_function,predict_features,decode,voltage_features)
from .alpha_gamma_ml_study import summarize_scores
from .alpha_gamma_ml_hybrid import fit_hybrid,MixedHybridConfig
from .alpha_gamma_pileup_fit import score_mixed
from .alpha_gamma_regularized_fit import CalibratedBank


def truth(times,classes,charges):
    return dict(first_active_ns=np.array(times,float),emission_ns=np.array(times,float),
                particle_name=np.array(classes),charge_pe=np.array(charges,float))


class MixedMLTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)
        cls.cal=json.loads(Path('pileup/experiments/alpha_gamma_ml_20261016/calibration.json').read_text())

    def test_subsample_and_two_classes_are_distinct_targets(self):
        t=np.arange(512)*2.;tr=truth([100.6,100.6,201.1],['alpha','gamma','gamma'],[1300,700,900])
        y,c=supervised_targets(t,tr)
        self.assertEqual(y[0,0,50],1.);self.assertEqual(y[1,0,50],1.)
        self.assertAlmostEqual(float(4*y[1,2,50]+t[50]),100.6,places=5)
        self.assertEqual(c,dict(alpha=0,gamma=0))

    def test_same_bin_collision_is_reported(self):
        _,c=supervised_targets(np.arange(512)*2.,truth([100.1,100.4],['gamma','gamma'],[600,900]))
        self.assertEqual(c['gamma'],1)

    def test_exact_encoded_events_recover_time_charge_and_class(self):
        t=np.arange(512)*2.;tr=truth([100.6,120.4,300.9],['alpha','gamma','alpha'],[1300,700,900])
        labels,_=supervised_targets(t,tr)
        dense=dict(event_score=labels[:,0],charge_pe=1000*np.exp(labels[:,1]),time_offset_ns=4*labels[:,2])
        fitted=decode(t,dense,[.5,.5]);s=score_mixed(fitted,tr,[0,1024])
        self.assertEqual(s['charge_and_type_recovered_count'],3)
        self.assertEqual(s['unmatched_predictions'],0)

    def test_unknown_empty_count_not_forced_positive(self):
        t=np.arange(512)*2.;dense=dict(event_score=np.full((2,512),.01),
                                     charge_pe=np.ones((2,512))*1000,time_offset_ns=np.zeros((2,512)))
        self.assertEqual(decode(t,dense)['times_ns'],[])

    def test_padding_and_core_seams_preserve_each_sample(self):
        class PointModel(nn.Module):
            def forward(self,x):
                return x[:,:6].reshape(x.shape[0],2,3,x.shape[-1])
        x=np.random.default_rng(7).normal(0,.1,(8,4387)).astype(np.float32)
        dense=predict_features(PointModel(),x,batch_size=3)
        np.testing.assert_allclose(dense['event_score'][0],1/(1+np.exp(-x[0])),rtol=1e-6)
        np.testing.assert_allclose(dense['time_offset_ns'][1],4*np.tanh(x[5]),rtol=1e-6)

    def test_both_architectures_have_finite_class_balanced_gradients(self):
        cfg=NetworkConfig();t=np.arange(cfg.crop_samples)*2.
        y,_=supervised_targets(t,truth([1701.1,1750.8],['alpha','gamma'],[800,600]))
        for kind in ('cnn','unet'):
            net=MixedEventNet(kind);z=net(torch.zeros(1,8,cfg.crop_samples))
            loss,_=loss_function(z,torch.from_numpy(y[None]));loss.backward()
            self.assertTrue(torch.isfinite(loss))
            self.assertTrue(all(torch.isfinite(p.grad).all() for p in net.parameters()))

    def test_empty_gamma_predictions_are_penalized_in_macro(self):
        tr=truth([100,200],['alpha','gamma'],[1000,1000]);fit=dict(times_ns=[100.],particle_name=['alpha'],charge_pe=[1000.])
        s=score_mixed(fit,tr,[0,500]);m=summarize_scores([dict(row=dict(gamma_model='brems',rate_hz=1e6),
            score=s,fit=fit,truth=tr,gate_ns=[0,500])])
        self.assertEqual(m['objective'],.5)
        # Real NumPy truth charges must produce JSON-serializable validation
        # counts; a positive match previously left a NumPy int64 here.
        self.assertEqual(json.loads(json.dumps(m))['objective'],.5)

    def test_features_validate_sampling_and_are_finite(self):
        t=np.arange(2048)*2.;v=np.full_like(t,self.cal['baseline_mV'])
        x=voltage_features(t,v,self.cal);self.assertEqual(x.shape,(8,2048));self.assertTrue(np.isfinite(x).all())
        with self.assertRaises(ValueError):voltage_features(t*2,v,self.cal)

    def test_hybrid_can_change_initial_classes_and_use_nn_only_candidate(self):
        # Known constructed candidates test optimizer mechanics, not blind accuracy.
        t=np.arange(2048)*2.;bank=CalibratedBank(self.cal);times=[1500.6,1550.4];q=[1200.,800.]
        voltage=bank.matrix(t,times,[0,1],q)@q+self.cal['baseline_mV']
        baseline=dict(times_ns=times[:1],particle_name=['gamma'],charge_pe=q[:1],
                      all_components=dict(times_ns=times[:1],particle_name=['gamma'],charge_pe=q[:1]),
                      model_mV=voltage)
        tr=truth(times,['alpha','gamma'],q);labels,_=supervised_targets(t,tr)
        dense=dict(event_score=labels[:,0],charge_pe=1000*np.exp(labels[:,1]),time_offset_ns=4*labels[:,2])
        neural=decode(t,dense)
        fitted=fit_hybrid(t,voltage,self.cal,neural,dense,baseline,MixedHybridConfig(strength=10.,passes=1,max_additions=0))
        s=score_mixed(fitted,tr,[1000,2500])
        self.assertEqual(s['charge_and_type_recovered_count'],2)
        self.assertEqual(fitted['fit_truth_inputs'],[])


if __name__=='__main__':unittest.main()
