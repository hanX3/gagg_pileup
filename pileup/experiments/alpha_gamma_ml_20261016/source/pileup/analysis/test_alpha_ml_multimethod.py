"""Window-coordinate and mathematical physical-fit controls (optional torch)."""
import unittest
import numpy as np
from scipy.optimize import nnls

try:
    import torch
    from .alpha_ml_models import make_model,windowed_prediction
    from .alpha_pileup_hybrid import HybridConfig,fit_candidates,amplitude_fit
    from .alpha_full_signal_study import make_template
    from .alpha_pileup_fit import FitConfig
except ImportError:
    torch=None


@unittest.skipIf(torch is None,'Optional CPU PyTorch dependency unavailable')
class AlternativeModelChecks(unittest.TestCase):
    def test_guarded_windows_preserve_coordinates_at_seams(self):
        class Identity(torch.nn.Module):
            def forward(self,x):return x[:,:3]
        for n in (1031,10000):
            rng=np.random.default_rng(18);x=rng.normal(0,.5,(5,n)).astype(np.float32)
            predicted=windowed_prediction(Identity(),x)
            np.testing.assert_allclose(predicted['probability'],1/(1+np.exp(-x[0])),rtol=1e-6)
            np.testing.assert_allclose(predicted['charge_pe'],3000*np.exp(x[1]),rtol=1e-6)
            np.testing.assert_allclose(predicted['time_offset_ns'],4*np.tanh(x[2]),rtol=1e-6)

    def test_both_architectures_have_all_three_gradients(self):
        torch.set_num_threads(1);torch.manual_seed(13)
        for kind in ('unet','transformer'):
            model=make_model(kind)
            prediction=model(torch.randn(2,5,1024))
            self.assertEqual(tuple(prediction.shape),(2,3,1024))
            prediction.square().mean().backward()
            gradients=[p.grad for p in model.parameters() if p.requires_grad]
            self.assertTrue(all(g is not None and bool(torch.isfinite(g).all()) for g in gradients))
            self.assertGreater(sum(float(g.abs().sum()) for g in gradients),0.)

    def test_large_bounded_solver_agrees_with_augmented_nnls(self):
        rng=np.random.default_rng(21);H=rng.uniform(0,1,(400,80))
        prior=rng.uniform(.1,2,80);values=H@rng.uniform(0,2,80)
        strength=.1
        predicted,info=amplitude_fit(H,values,prior,prior,strength,400)
        penalty=np.sqrt(strength*np.sum(H*H,axis=0))
        expected=nnls(np.vstack((H,np.diag(penalty))),np.r_[values,penalty*prior],maxiter=2000)[0]
        np.testing.assert_allclose(predicted,expected,rtol=2e-4,atol=2e-5)
        self.assertLess(info['kkt_relative'],1e-6)

    def test_physical_amplitude_refit_can_correct_wrong_neural_charge(self):
        import json
        from pathlib import Path
        settings=json.loads(Path('pileup/experiments/alpha_ml_multimethod_20261014/calibration-settings.json').read_text())
        cfg=FitConfig(**settings['fit_config']);template=make_template(settings)
        t=np.arange(4000)*2.+.25
        times=np.array([2200.3,2600.7]);charges=np.array([4000.,2000.])
        voltage=cfg.baseline_mV+template(t[:,None]-times)@(charges*cfg.area_per_pe_mV_ns*template.unit_area_peak)
        # Known candidate coordinates are deliberately used ONLY to test the solver,
        # not as blind event-reconstruction performance evidence.
        fitted=fit_candidates(t,voltage,times,[1200.,6000.],settings,HybridConfig(time_refine_ns=0,time_passes=0))
        np.testing.assert_allclose(fitted['charge_pe'],charges,rtol=1e-8)
        self.assertLess(fitted['rms_mV'],1e-8)

    def test_empty_physical_candidate_set(self):
        import json
        from pathlib import Path
        settings=json.loads(Path('pileup/experiments/alpha_ml_multimethod_20261014/calibration-settings.json').read_text())
        t=np.arange(100)*2.+.25
        result=fit_candidates(t,np.full_like(t,10.),[],[],settings,HybridConfig())
        self.assertEqual(result['times_ns'],[]);self.assertEqual(result['charge_pe'],[])
        self.assertEqual(result['rms_mV'],0.)


if __name__=='__main__':unittest.main()
