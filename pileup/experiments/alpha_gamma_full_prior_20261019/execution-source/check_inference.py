from pathlib import Path
import importlib.util
import json
import sys
import numpy as np
import torch

record=Path(__file__).resolve().parent
work=Path('/tmp/gagg-full-prior-20261019')
sys.path.insert(0,str(work/'source'));sys.path.insert(0,str(work/'dependency-source'))
from pileup.analysis.alpha_gamma_ml_study import load_model
from pileup.analysis.alpha_gamma_ml import predict_features,voltage_features,decode
torch.set_num_threads(1)
spec=importlib.util.spec_from_file_location('frozen_full_test',work/'test-source/alpha_gamma_full_prior_test.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
stem=sorted((work/'cases').glob('*.json'))[0].stem
case=json.loads((work/'cases'/f'{stem}.json').read_text())
wave=dict(np.load(work/case['inputs']['waveform']['path']))
cal=json.loads((work/'calibration.json').read_text())
model,cfg,frozen=load_model(work/'models/unet')
fits,dense,seconds=module.infer(wave['time_ns'],wave['analog_mV'],cal,model,cfg,frozen['selected']['thresholds'])
cached=dict(np.load(work/case['inputs']['dense']['path']))
maxima={}
for name in cached:
    np.testing.assert_allclose(dense[name],cached[name],rtol=1e-6,atol=1e-5)
    maxima[name]=float(np.max(abs(dense[name]-cached[name])))
references={
    'baseline':json.loads((work/case['inputs']['baseline']['path']).read_text())['fit'],
    'unet':json.loads((work/case['inputs']['neural']['path']).read_text())['fit'],
    **{m:json.loads((work/'outputs/fits'/old/f'{stem}.json').read_text())['fit'] for m,old in
       [('hybrid_none','replay_none'),('hybrid_both','replay_both'),('class_only','class_only')]}}
fit_differences={}
for method,reference in references.items():
    fit=fits[method];assert fit['particle_name']==reference['particle_name']
    diff=0.
    for key in ('times_ns','charge_pe'):
        np.testing.assert_allclose(fit[key],reference[key],rtol=1e-8,atol=1e-6)
        if len(fit[key]):diff=max(diff,float(np.max(abs(np.asarray(fit[key])-reference[key]))))
    fit_differences[method]=diff
result=dict(stem=stem,torch=torch.__version__,frozen_model_sha256=frozen['model_sha256'],
    dense_max_abs_differences=maxima,full_end_to_end_fit_max_abs_differences=fit_differences,
    all_five_methods_reproduced=True,truth_opened=False,seconds=seconds)
(record/'inference-preflight.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
