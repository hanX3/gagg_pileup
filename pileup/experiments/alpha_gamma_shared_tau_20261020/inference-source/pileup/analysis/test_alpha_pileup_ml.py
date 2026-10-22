"""Physical label and event-coordinate checks; CPU torch is optional."""
import unittest
import numpy as np

try:
    from .alpha_pileup_ml import supervised_targets, decode_prediction, EventNet, loss_function
    import torch
except ImportError:
    torch=None


@unittest.skipIf(torch is None,'Optional PyTorch environment not installed')
class EventModelChecks(unittest.TestCase):
    def test_charge_cut_and_first_photon_time(self):
        t=np.arange(100)*2.+.25
        truth=dict(first_active_ns=np.array([40.25,80.25,np.nan]),
                   emission_ns=np.array([12.,20.,30.]),charge_pe=np.array([4000.,499.,8000.]))
        labels=supervised_targets(t,truth)
        self.assertEqual(np.argmax(labels[0]),20)
        self.assertAlmostEqual(labels[1,20],np.log(4000/3000),places=6)
        self.assertAlmostEqual(labels[2,20],0.)
        self.assertLess(labels[0,40],1e-20)

    def test_unknown_multiplicity_and_subthreshold_filter(self):
        t=np.arange(100)*2.+.25
        dense=dict(probability=np.zeros(100),charge_pe=np.ones(100)*2000.,time_offset_ns=np.zeros(100))
        dense['probability'][[10,30,60,90]]=.9
        dense['charge_pe'][30]=499.
        dense['time_offset_ns'][60]=.75
        result=decode_prediction(t,dense,.5)
        self.assertEqual(len(result['times_ns']),3)
        self.assertEqual(result['times_ns'],[20.25,121.,180.25])
        self.assertEqual(decode_prediction(t,dense,.95)['times_ns'],[])

    def test_training_guards_cannot_change_loss(self):
        torch.manual_seed(17)
        z=torch.randn(2,3,1024,requires_grad=True)
        labels=torch.zeros(2,4,1024)
        loss,_=loss_function(z,labels)
        other=labels.clone();other[:,:,:262]=100.;other[:,:,762:]=100.
        changed,_=loss_function(z,other)
        self.assertEqual(float(loss.detach()),float(changed.detach()))
        loss.backward()
        self.assertEqual(float(z.grad[:,:,:262].abs().sum()),0.)
        self.assertEqual(float(z.grad[:,:,762:].abs().sum()),0.)

    def test_arbitrary_record_length_and_no_truth_argument(self):
        torch.set_num_threads(1)
        model=EventNet().eval()
        with torch.inference_mode():
            for n in (1024,10000):
                prediction=model(torch.zeros(1,5,n))
                self.assertEqual(tuple(prediction.shape),(1,3,n))
                self.assertTrue(bool(torch.isfinite(prediction).all()))


if __name__=='__main__':
    unittest.main()
