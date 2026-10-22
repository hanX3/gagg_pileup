"""Check one complete optimizer objective from internal numerical state."""
from pathlib import Path
import json
import sys
import numpy as np

record=Path(__file__).resolve().parent;work=Path('/tmp/gagg-full-prior-20261019')
sys.path.insert(0,str(work/'source'));sys.path.insert(0,str(work/'dependency-source'))
from pileup.analysis.alpha_gamma_full_prior_fit import fit_full_prior,FullPriorConfig
stem=sorted((work/'cases').glob('*.json'))[0].stem
case=json.loads((work/'cases'/f'{stem}.json').read_text())
inputs={key:work/value['path'] for key,value in case['inputs'].items()}
wave=dict(np.load(inputs['waveform']));dense=dict(np.load(inputs['dense']))
baseline=json.loads(inputs['baseline'].read_text())['fit']
baseline['model_mV']=np.load(inputs['baseline_model'])['model_mV']
neural=json.loads(inputs['neural'].read_text())['fit']
cal=json.loads((work/'calibration.json').read_text())
checked={}

def profile(frame,event,result):
    if frame.f_code is fit_full_prior.__code__ and event=='return':
        values=frame.f_locals;inverse=np.argsort(values['order'])
        q=values['q'][inverse];kinds=values['types'][inverse]
        error=values['W'].apply((values['full_y']-values['model'])[::values['config'].fit_stride])
        sse=float(error@error)
        charge=float(np.sum(values['weight']*(q-values['prior'])**2))
        classes=float(values['class_cost'][np.arange(len(kinds)),kinds].sum())
        expected=values['costs'][-1]
        assert charge==0.
        np.testing.assert_allclose(sse+charge+classes,expected,rtol=1e-12,atol=1e-8)
        checked.update(stem=stem,weighted_voltage_sse=sse,explicit_charge_penalty=charge,
            explicit_class_penalty=classes,recorded_total=expected,absolute_difference=abs(sse+charge+classes-expected),
            candidates_including_zero_charge=len(q),zero_charge_candidates=int(np.count_nonzero(q<=1e-6)),
            interpretation='State trace of unchanged frozen function, including penalties on zero-Q candidates')

sys.setprofile(profile)
try:
    fitted=fit_full_prior(wave['time_ns'],wave['analog_mV'],cal,neural,dense,baseline,FullPriorConfig())
finally:
    sys.setprofile(None)
reference=json.loads((work/'outputs/fits/class_only'/f'{stem}.json').read_text())['fit']
assert fitted['times_ns']==reference['times_ns'] and fitted['charge_pe']==reference['charge_pe']
assert checked
(record/'objective-state-check.json').write_text(json.dumps(checked,indent=2)+'\n')
print(json.dumps(checked),flush=True)
