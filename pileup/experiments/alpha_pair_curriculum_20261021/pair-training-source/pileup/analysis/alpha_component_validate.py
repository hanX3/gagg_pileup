"""Recompute stored results, native provenance and frozen inference checks."""
from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path
import re
from unittest.mock import patch

import awkward as ak
import numpy as np

from .alpha_component_ml import PROJECT,EXPERIMENT,OUTPUT,VARIANTS,SEEDS,decode_prediction
from .alpha_component_study import load,RECORD
from .alpha_component_evaluate import check_freeze,infer,networks,configuration
from .alpha_full_signal_study import make_template
from .alpha_pileup_fit import score
from .alpha_rate_campaign import save,digest
from .alpha_rate_data_hashed import open_root


def compatible(value):return json.loads(json.dumps(value))


def numerical_fit(value):
    """Exclude only measured wall time from an exact reconstruction replay."""
    return {key:item for key,item in value.items() if key!='refinement_seconds'}


def reconstructed_model(time,times,amplitudes,baseline,template):
    """Identical per-sample sum without full-record interpolation temporaries."""
    output=np.empty_like(time,dtype=float)
    for start in range(0,len(time),128):
        output[start:start+128]=baseline+template(time[start:start+128,None]-times[None,:])@amplitudes
    return output


def validate():
    frozen=check_freeze();manifest=load(EXPERIMENT/'manifest.json');summary=load(OUTPUT/'test-summary.json')
    plans=load(EXPERIMENT/'assets/test-plan.json')['planned_runs']
    assert manifest['planned_runs']==plans
    assert len(manifest['runs'])==len(plans)==len(summary['windows'])
    assert {r['seed'] for r in plans}=={r['seed'] for r in manifest['runs']}
    assert all(r['status']=='completed' for r in manifest['runs'])
    split=load(EXPERIMENT/'data-split.json')['windows']
    train={r['root_sha256'] for r in split if r['split']=='train'}
    valid={r['root_sha256'] for r in split if r['split']=='validation'}
    test={r['root_sha256'] for r in manifest['runs']}
    assert len(train)==49 and len(valid)==26 and len(test)==len(plans)
    assert not (train&valid or train&test or valid&test)
    assert not {r['seed'] for r in split}&{r['seed'] for r in plans}
    root_checks=[]
    fields=['random_seed','sipm_avalanche_base_seed','sipm_readout_base_seed','alpha_p11b_rate_hz',
            'brems_enabled','c12_capture_enabled','t_length_ps']
    for r in manifest['runs']:
        assert datetime.fromisoformat(r['started'])>datetime.fromisoformat(frozen['frozen_at'])
        root=PROJECT/r['root'];assert root.parent==PROJECT/'g4/data'
        assert re.fullmatch(r'\d{8}_\d{2}h\d{2}m\d{2}s\.root',root.name)
        assert digest(PROJECT/r['macro'])==r['macro_sha256']
        meta=load(PROJECT/r['export_metadata'])
        assert meta['root_sha256']==r['root_sha256'] and meta['seed']==r['seed']
        assert meta['macro_contents'].strip()==(PROJECT/r['macro']).read_text().strip()
        assert meta['duration_ns']==20000. and meta['gate_ns']==[5000.,15000.]
        assert meta['source_rates_hz']['brems']==meta['source_rates_hz']['c12_capture']==0
        with open_root(root) as f:
            assert all(f[k].num_entries==1 for k in ('RunInfo','WaveformEvent','SiPMEvent'))
            info=ak.to_list(f['RunInfo'].arrays(fields))[0]
        assert info['random_seed']==r['seed']
        assert info['sipm_avalanche_base_seed']==r['sipm_seed']==r['seed']+100000
        assert info['sipm_readout_base_seed']==r['readout_seed']==r['seed']+200000
        assert info['alpha_p11b_rate_hz']==r['rate_hz'] and info['t_length_ps']==20000000
        assert not info['brems_enabled'] and not info['c12_capture_enabled']
        root_checks.append(dict(root=r['root'],seed=r['seed'],actual_runinfo_verified=True))
    settings,_,_=configuration();template=make_template(settings)
    score_count=dense_count=model_count=0;max_model_error=0.
    for r in summary['windows']:
        with np.load(PROJECT/r['truth']) as f:
            truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe')}
        with np.load(PROJECT/r['waveform']) as f:t=f['time_ns']
        assert digest(PROJECT/r['waveform'])==r['waveform_sha256']
        assert digest(PROJECT/r['truth'])==r['truth_sha256']
        stem=Path(r['root']).stem
        for name,m in r['methods'].items():
            s=compatible(score(m['fit'],truth,r['gate_ns']))
            assert s==m['score'],(stem,name)
            assert s['target_count']==s['recovered_count']+(s['matched_count']-s['recovered_count'])+s['missed_count']
            score_count+=1
            if name in ('sse','wls'):
                fit=m['fit'];full=fit.get('all_components',fit)
                tt=np.asarray(full['times_ns']);aa=np.asarray(full['peak_mV'])
                actual=reconstructed_model(t,tt,aa,settings['fit_config']['baseline_mV'],template)
                with np.load(OUTPUT/'test/models'/f'{stem}_{name}.npz') as f:expected=f['model_mV']
                error=float(np.max(abs(actual-expected)));max_model_error=max(error,max_model_error)
                assert error<1e-8;model_count+=1
            else:
                base_name=name.removesuffix('_matched_extra')
                with np.load(OUTPUT/'test/dense'/f'{stem}_{base_name}.npz') as f:dense={k:f[k] for k in f.files}
                assert compatible(decode_prediction(t,dense,m['fit']['threshold']))==m['fit']
                dense_count+=1
    caches=load(OUTPUT/'cache-validation.json')['windows']
    assert len(caches)==75
    for r in caches:assert digest(PROJECT/r['cache'])==r['cache_sha256']
    # Same architecture/seed objectives have paired initial weights and crops.
    for architecture in ('gru','unetgru'):
        for seed in SEEDS:
            a=load(OUTPUT/'models'/f'{architecture}_base_{seed}'/'config.json')
            b=load(OUTPUT/'models'/f'{architecture}_component_{seed}'/'config.json')
            for key in ('seed','steps','batch_size','split_sha256','calibration_sha256','parameters','source_sha256'):
                assert a[key]==b[key]
    first=min((r for r in summary['windows'] if r['rate_hz']==1000000),key=lambda r:r['seed'])
    with np.load(PROJECT/first['waveform']) as f:t,y=f['time_ns'].copy(),f['analog_mV'].copy()
    nets=networks()
    with patch('numpy.load',side_effect=AssertionError('Array/truth file opened inside blind boundary')):
        replay,_,_=infer(t,y,nets)
    for name,fit in replay.items():
        fit.pop('model_mV',None)
        assert numerical_fit(compatible(fit))==numerical_fit(first['methods'][name]['fit']),name
    galleries=load(OUTPUT/'gallery/index.json')['figures']
    assert len(galleries)==8
    for g in galleries:
        with np.load(OUTPUT/'gallery'/f'rate_{g["rate_hz"]/1e6:g}MHz'/'components.npz') as f:
            assert len(f['time_ns'])==10000
            assert np.max(abs(f['predicted_components_mV'].sum(0)+settings['fit_config']['baseline_mV']-f['neural_total_mV']))<1e-8
        assert g['native_additivity_max_error_mV']<1e-8
        assert g['ancestry_charge_max_error_pe']<1e-7
    result=dict(status='passed',frozen_files_checked=len(frozen['frozen_files']),
        fresh_roots_checked=len(root_checks),recomputed_scores=score_count,
        independently_redecoded_dense_predictions=dense_count,physical_models_reconstructed=model_count,
        physical_model_max_error_mV=max_model_error,cache_files_checked=75,
        paired_training_configurations=4,inference_with_array_loading_forbidden=True,
        cold_same_input_exact_replay_methods=len(replay),full_record_galleries_checked=8,
        replay_exempt_profiling_fields=['refinement_seconds'],
        test_events_excluded=0,scoring_inputs_to_inference=[],source_rate_inputs_to_inference=[],
        root_provenance=root_checks)
    save(OUTPUT/'validation.json',result)
    print('VALIDATED',score_count,'scores',model_count,'physical models,',len(root_checks),'new roots',flush=True)


if __name__=='__main__':validate()
