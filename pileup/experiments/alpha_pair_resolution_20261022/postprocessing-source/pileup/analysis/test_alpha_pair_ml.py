"""Meaningful boundary/assignment/conservation checks for the new event slots."""
import json
from pathlib import Path
import unittest

import numpy as np
import torch

from .alpha_pair_ml import ComponentNet,blind_predict,set_loss


class PairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)
        path=Path(__file__).resolve().parents[2]/'pileup/experiments/alpha_component_ml_20261021/calibration-settings.json'
        cls.settings=json.loads(path.read_text())
        cls.settings['component_fixed_norm']=1e7

    def test_two_slot_assignment_is_permutation_invariant_and_has_finite_gradients(self):
        torch.manual_seed(9)
        model=ComponentNet(self.settings,head='direct')
        y=torch.ones((2,1024))*10.
        y[:,100:400]+=30.
        pred=model(y)
        target=dict(active=torch.tensor([[1.,1.],[1.,0.]]),
                    charge_pe=torch.tensor([[4000.,1000.],[2500.,1.]]),
                    times_ns=torch.tensor([[180.,184.],[200.,300.]]),
                    components_mV=torch.stack((torch.stack((y[0]-10,(y[0]-10)*.25)),torch.stack((y[0]-10,torch.zeros(1024))))))
        reverse={k:v[:,[1,0]] for k,v in target.items()}
        for norm in ('fixed','relative'):
            a,_=set_loss(pred,target,self.settings,normalization=norm)
            b,_=set_loss(pred,reverse,self.settings,normalization=norm)
            torch.testing.assert_close(a,b)
        a.backward()
        assert all(torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None)

    def test_direct_masks_conserve_positive_observed_voltage(self):
        model=ComponentNet(self.settings,head='direct')
        y=torch.rand((3,1024))*100.+9.8
        p=model(y)
        self.assertTrue(torch.all(p['components_mV']>=0))
        self.assertTrue(torch.all(p['components_mV'].sum(1)<=(y-10.).clamp_min(0.)+1e-4))
        self.assertTrue(torch.isfinite(p['charge_pe']).all())

    def test_inference_retains_two_close_event_slots_without_truth(self):
        class Fixed(torch.nn.Module):
            slots=2
            def forward(self,x):
                return dict(presence_logits=torch.tensor([[8.,8.]]),
                            times_ns=torch.tensor([[100.,104.]]),charge_pe=torch.tensor([[4000.,1000.]]))
        fit,_=blind_predict(Fixed(),np.full(1024,10.),threshold=.5)
        self.assertEqual(fit['times_ns'],[100.,104.])
        self.assertEqual(fit['truth_inputs'],[])


if __name__=='__main__':unittest.main()
