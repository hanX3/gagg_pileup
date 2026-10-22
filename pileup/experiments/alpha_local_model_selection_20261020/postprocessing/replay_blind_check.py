"""Cold replay of one real test waveform; never open scoring truth arrays."""
import json
from pathlib import Path
import sys
import numpy as np

project=Path.cwd();exp=project/'pileup/experiments/alpha_local_model_selection_20261020'
out=project/'jupyter/artifacts/alpha_local_model_selection_20261020';rec=project/'jupyter/run_records/alpha_local_model_selection_20261020'
sys.path.insert(0,str(exp/'inference-source'))
from pileup.analysis.alpha_local_fit import LocalConfig,fit_local
from pileup.analysis.alpha_pileup_fit import FitConfig,fit_waveform
from pileup.analysis.alpha_full_signal_study import make_template
from pileup.analysis.alpha_rate_campaign import digest,save

frozen=json.loads((out/'freeze.json').read_text())
for name,expected in frozen['frozen_files'].items():assert digest(project/name)==expected,name
settings=json.loads((out/'baseline-settings.json').read_text());template=make_template(settings)
base=FitConfig(**settings['fit_config']);cfg=LocalConfig(**frozen['selected']['config'])
stem='20261020_17h34m43s'
with np.load(out/'inputs/waveforms'/f'{stem}.npz') as z:t,y=z['time_ns'],z['analog_mV']
initial=fit_waveform(t,y,template,base);replay=fit_local(t,y,template,base,cfg,initial=initial)
old=json.loads((out/'test/selected'/f'{stem}.json').read_text())['methods']['local']['fit']
assert len(replay['times_ns'])==len(old['times_ns'])
np.testing.assert_allclose(replay['times_ns'],old['times_ns'],atol=1e-3,rtol=0)
np.testing.assert_allclose(replay['charge_pe'],old['charge_pe'],atol=.01,rtol=0)
with np.load(out/'test/selected'/f'{stem}_local.npz') as z:
    np.testing.assert_allclose(replay['model_mV'],z['model_mV'],atol=1e-4,rtol=0)
    delta=float(np.max(abs(replay['model_mV']-z['model_mV'])))
save(rec/'blind-replay-validation.json',dict(passed=True,root=f'g4/data/{stem}.root',
    scoring_truth_arrays_opened=False,inference_inputs=['time_ns','analog_mV','independent fixed template/config','observed-voltage baseline warm start'],
    parameters_received_from_truth=[],max_model_difference_mV=delta,
    time_tolerance_ns=.001,charge_tolerance_pe=.01,model_tolerance_mV=.0001,
    reported_components=len(replay['times_ns'])))
print('BLIND_REPLAY_PASSED',len(replay['times_ns']),'max model difference',delta,flush=True)
