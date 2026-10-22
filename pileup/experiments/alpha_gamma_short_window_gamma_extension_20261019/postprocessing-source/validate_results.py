"""Validate new ROOT provenance and independently rescore prefix outputs."""
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys

import awkward as ak
import numpy as np

project=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(project))
from analysis.alpha_gamma_short_window import PrefixCalibration, fit_prefix
from analysis.alpha_gamma_short_window_study import digest, save
from pileup.analysis.alpha_rate_data import open_root

record=Path(__file__).resolve().parent
output=project/'jupyter/artifacts/alpha_gamma_short_window_implementation_20261019'
settings=json.loads((output/'frozen-settings.json').read_text())
freeze=json.loads((output/'freeze-complete.json').read_text())
frozen_at=datetime.fromisoformat(freeze['created_utc'])
for path,expected in settings['source_sha256'].items():assert digest(project/path)==expected
assert digest(output/'calibration/calibration-arrays.npz')==settings['calibration_arrays_sha256']
assert digest(output/'calibration/calibration.json')==settings['calibration_sha256']
experiments=[project/'pileup/experiments'/name for name in (
    'alpha_gamma_short_window_20261019','alpha_gamma_short_window_gamma_extension_20261019')]
runs=[]
for folder in experiments:
    data=json.loads((folder/'manifest.json').read_text())
    assert len(data['runs'])==len(data['planned_runs'])
    assert all(r['status']=='completed' and r.get('validation') for r in data['runs'])
    for path,expected in data['source_sha256'].items():assert digest(folder/path)==expected
    runs.extend(data['runs'])
assert len(runs)==174
assert len({r['root'] for r in runs})==len({r['seed'] for r in runs})==174
assert sum(r['particle']=='alpha' for r in runs)==30
assert sum(r['particle']=='gamma' for r in runs)==144
times=sorted(datetime.fromisoformat(r['started']) for r in runs)
assert times[0]>=frozen_at and min((b-a).total_seconds() for a,b in zip(times,times[1:]))>=1.
source_config=None;data_summary=[]
for i,row in enumerate(runs):
    root=project/row['root'];stem=root.stem
    assert re.fullmatch(r'\d{8}_\d{2}h\d{2}m\d{2}s',stem)
    assert digest(root)==row['root_sha256']
    macro=project/row['macro'];assert digest(macro)==row['macro_sha256']
    metadata=json.loads((output/'native-inputs/metadata'/f'{stem}.json').read_text())
    assert metadata['macro_contents']==macro.read_text() and metadata['primaries']==1
    assert metadata['root_sha256']==row['root_sha256']
    config=metadata['config']
    if source_config is None:source_config=config
    assert config==source_config
    with open_root(root) as f:
        assert all(f[k].num_entries==1 for k in ('RunInfo','WaveformEvent','PrimaryPhoton','SiPMEvent'))
        info=ak.to_list(f['RunInfo'].arrays(['random_seed','macro_contents','sipm_avalanche_base_seed',
                 'sipm_readout_base_seed','t_length_ps','sipm_config_json']))[0]
        assert info['random_seed']==row['seed'] and info['sipm_avalanche_base_seed']==row['sipm_seed']
        assert info['sipm_readout_base_seed']==row['readout_seed']
        assert info['macro_contents']==macro.read_text() and info['t_length_ps']==20_000_000
        assert json.loads(info['sipm_config_json'])==config
        primary=ak.to_list(f['PrimaryPhoton'].arrays(['primary_particle_name','primary_energy_MeV']))[0]
        assert primary['primary_particle_name']==row['particle']
        assert abs(primary['primary_energy_MeV']-row['energy_MeV'])<1e-7
    with np.load(output/'native-inputs/waveforms'/f'{stem}.npz') as wave, np.load(output/'native-inputs/truth'/f'{stem}.npz') as truth:
        t,y=wave['time_ns'],wave['analog_mV'];assert len(t)==len(y)==10000
        np.testing.assert_allclose(np.diff(t),2.,atol=0,rtol=0);assert np.isfinite(y).all()
        assert len(truth['charge_pe'])==1
        q=float(truth['charge_pe'][0]);q0=config['sensor']['gain_electrons']*1.602176634e-19
        recalculated=float(np.sum(truth['avalanche_charge_C'][truth['avalanche_owner']>=0])/q0)
        np.testing.assert_allclose(q,recalculated,rtol=1e-11,atol=1e-7)
        assert abs(q+metadata['dark_charge_pe']-metadata['charge_pe'])<1e-6
    result=json.loads((output/'test'/f'{stem}.json').read_text())
    assert result['fitter_sha256']==settings['source_sha256']['analysis/alpha_gamma_short_window.py']
    assert result['study_sha256']==settings['source_sha256']['analysis/alpha_gamma_short_window_study.py']
    assert result['calibration_sha256']==settings['calibration_arrays_sha256']
    assert result['inference_truth_inputs']==[] and len(result['prefixes'])==44
    assert abs(result['truth_for_evaluation_only']['charge_pe']-q)<1e-9
    data_summary.append({'root':row['root'],'particle':row['particle'],'charge_pe':q,
                         'waveform_sha256':digest(output/'native-inputs/waveforms'/f'{stem}.npz'),
                         'truth_sha256':digest(output/'native-inputs/truth'/f'{stem}.npz')})
    if i%25==0:print('VALIDATED_ROOT',i+1,'/174',flush=True)

summary=json.loads((output/'test-summary.json').read_text())
results={r['root']:json.loads((output/'test'/f"{Path(r['root']).stem}.json").read_text()) for r in runs}
for group in summary['groups']:
    low,high={'all_above_500':(500.,np.inf),'500_1000':(500.,1000.),
              '1000_3000':(1000.,3000.),'above_3000':(3000.,np.inf)}[group['band']]
    selected=[r for r in data_summary if r['particle']==group['particle'] and low<=r['charge_pe']<high]
    successes=0;matched=0;types=0;errors=[]
    for data in selected:
        event=results[data['root']];truth=event['truth_for_evaluation_only']
        fit=next(f for f in event['prefixes'] if f['pre_ns']==group['pre_ns'] and f['post_ns']==group['post_ns'])
        onset=truth['onset_first_active_ns']
        good=bool(fit['detected'] and onset is not None and abs(fit['trigger']['time_ns']-onset)<=60.)
        assert fit['evaluation']['matched_signal']==good
        if not good:continue
        matched+=1;r=fit['methods'][group['method']]
        error=(r['charge_pe']-data['charge_pe'])/data['charge_pe'];errors.append(error)
        successes+=abs(error)<=.05;types+=r['particle']==group['particle']
    assert group['targets']==len(selected) and group['charge_5pct_success']==successes and group['detected']==matched
    if errors:np.testing.assert_allclose(group['charge_bias'],np.mean(errors),atol=1e-13)
    if group['correct_type'] is not None:assert group['correct_type']==types

cal=PrefixCalibration.load(output/'calibration')
future_checks=[];replays=[]
for particle in ('alpha','gamma'):
    data=next(r for r in sorted(data_summary,key=lambda x:x['root']) if r['particle']==particle and r['charge_pe']>=500)
    stem=Path(data['root']).stem
    with np.load(output/'native-inputs/waveforms'/f'{stem}.npz') as wave:
        t,y=wave['time_ns'],wave['analog_mV']
        for pre,post in ((20.,100.),(50.,300.)):
            actual=fit_prefix(t,y,cal,pre_ns=pre,post_ns=post)
            future=y.copy();future[t>actual['last_observed_ns']]=np.nan
            assert actual==fit_prefix(t,future,cal,pre_ns=pre,post_ns=post)
            shifted=fit_prefix(t+2345.,y,cal,pre_ns=pre,post_ns=post)
            assert actual['methods']==shifted['methods']
            stored=next(f for f in results[data['root']]['prefixes'] if f['pre_ns']==pre and f['post_ns']==post)
            for name in actual['methods']:
                previous={k:v for k,v in stored['methods'][name].items() if k!='evaluation'}
                assert actual['methods'][name]==previous
            future_checks.append({'root':data['root'],'pre_ns':pre,'post_ns':post,'future_replaced_with_nan':True,'global_time_shift_invariance':True})
        destination=record/f'frozen-replay-{particle}.json'
        if not destination.exists():
            subprocess.run([sys.executable,str(record/'frozen_replay.py'),'--experiment',str(experiments[0]),
                '--root',data['root'],'--pre','20','--post','100','--output',str(destination)],cwd=project,check=True)
        replay=json.loads(destination.read_text());original=fit_prefix(t,y,cal,pre_ns=20.,post_ns=100.)
        assert replay==original
        replays.append(data['root'])
controls=json.loads((output/'controlled-doubles.json').read_text());assert len(controls['rows'])==72
save(record/'native-data-validation.json',{'created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed',
    'new_native_roots':174,'independent_groups_rescored':len(summary['groups']),
    'full_runinfo_macro_seed_and_config_checks':174,'charge_ancestry_checks':174,
    'preserved_raw_root_hashes':174,'future_checks':future_checks,'frozen_source_replays':replays,
    'controlled_prefixes':72,'fitter_sha256':settings['source_sha256']['analysis/alpha_gamma_short_window.py'],
    'source_frozen_before_native_generation':True,'input_data':data_summary})
print('NATIVE_PREFIX_VALIDATION_PASSED',flush=True)
