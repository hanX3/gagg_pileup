"""Numerical checks for the physical training loss and the blind boundary."""
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
import torch

from .alpha_component_ml import (PhysicalDecoder, OFFSETS, component_loss,
                                 blind_predict, model_for, batch_components)


class PhysicalTrainingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)
        cls.settings = json.loads(Path('pileup/experiments/alpha_ml_comparison_20261014/calibration-settings.json').read_text())

    def test_decoder_interpolation_gradient(self):
        decoder = PhysicalDecoder(self.settings)
        q = torch.tensor([2000.], requires_grad=True)
        dt = torch.tensor([1.37], requires_grad=True)
        reference = decoder.components(q.detach(), torch.zeros(1)).detach()
        value = (decoder.components(q, dt)-reference).square().sum()
        value.backward()
        def energy(t):
            return float((decoder.components(q.detach(), torch.tensor([t]))-reference).square().sum())
        numeric = (energy(1.38)-energy(1.36))/.02
        self.assertAlmostEqual(float(dt.grad)/numeric, 1., delta=.025)
        self.assertGreater(float(q.grad), 0.)

    def test_exact_components_empty_loss_and_error(self):
        decoder = PhysicalDecoder(self.settings)
        q = torch.tensor([1500.,4000.]); offsets = torch.tensor([.7,-1.2])
        curves = decoder.components(q, offsets).detach()
        target = dict(batch=torch.tensor([0,0]), anchor=torch.tensor([300,330]),
                      curve=curves, norm=curves.square().sum(1))
        z = torch.zeros(1,3,1024)
        z[0,1,target['anchor']] = torch.log(q/3000.)
        z[0,2,target['anchor']] = torch.atanh(offsets/4.)
        a,b = component_loss(z,target,decoder)
        self.assertLess(float(a+b),1e-10)
        z[0,1,300] += .3
        z.requires_grad_()
        a,b = component_loss(z,target,decoder)
        (a+.1*b).backward()
        self.assertGreater(float(a.detach()), .01)
        self.assertGreater(float(z.grad[0,1,300]),0)
        empty = dict(batch=torch.empty(0,dtype=torch.long))
        a,b = component_loss(z,empty,decoder)
        self.assertEqual(float((a+b).detach()),0.)

    def test_inference_does_not_read_arrays_or_truth(self):
        model = model_for('gru_base').eval()
        t = np.arange(10000)*2.+.25
        y = np.full_like(t,10.)
        with patch('numpy.load',side_effect=AssertionError('Truth/array file opened during inference')):
            pred,dense = blind_predict(model,t,y,self.settings,.95)
        self.assertEqual(len(dense['probability']),10000)
        self.assertEqual(pred['input_policy'],'complete_sampled_voltage')

    def test_paired_initializations_match(self):
        for seed in (262041001,262041002):
            for name in ('gru','unetgru'):
                torch.manual_seed(seed);a=model_for(name+'_base').state_dict()
                torch.manual_seed(seed);b=model_for(name+'_component').state_dict()
                self.assertTrue(all(torch.equal(a[k],b[k]) for k in a))


if __name__ == '__main__':
    unittest.main()
