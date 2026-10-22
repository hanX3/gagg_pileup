"""Recompute every saved score and validate provenance/inference boundaries."""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
import subprocess
import sys

import numpy as np

from .alpha_pair_data import PROJECT,EXPERIMENT,OUTPUT,RECORD,load,read_batch
from .alpha_pair_native import check_freeze
from .alpha_pair_study import batch_scores,SEEDS
from .alpha_pair_photons import one_case
from .alpha_rate_campaign import digest,save
from .alpha_pileup_fit import score


def validate():
    frozen=check_freeze()
    split=load(EXPERIMENT/'data-split.json')
    fresh=load(OUTPUT/'native-records.json')['test_single']
    groups=[{r['root'] for r in group} for group in (split['train'],split['validation'],fresh)]
    assert all(not groups[i].intersection(groups[j]) for i in range(3) for j in range(i+1,3))
    assert len(fresh)==24
    for r in fresh:
        meta=load(PROJECT/r['metadata'])
        assert digest(PROJECT/r['root'])==meta['root_sha256']==r['root_original_sha256']
    manifest=load(EXPERIMENT/'manifest.json')
    for r in manifest['runs']:
        assert r['status']=='completed' and r['returncode']==0
        assert r['started']>frozen['frozen_at']
        meta=load(OUTPUT/'native-inputs/metadata'/f"{Path(r['root']).stem}.json")
        assert meta['macro_contents']==(PROJECT/r['macro']).read_text()
        assert r['seed']==meta['seed']
        assert r['sipm_seed']==r['seed']+100000 and r['readout_seed']==r['seed']+200000
        q0=meta['config']['sensor']['gain_electrons']*1.602176634e-19
        area=q0*meta['config']['readout']['transimpedance_ohm']*1e12
        assert abs(area-load(EXPERIMENT/'settings.json')['fit_config']['area_per_pe_mV_ns'])<1e-8
    score_count=0;decoded=0;direct_count=0
    for kind in ('native','mean','linear','shared'):
        batch=read_batch(OUTPUT/f'assets/test-{kind}.npz')
        result=load(OUTPUT/'test'/kind/'results.json')
        assert digest(OUTPUT/f'assets/test-{kind}.npz')==result['input_sha256']
        names=list(result['rows'][0]['methods'])
        for name in names:
            fits=[r['methods'][name]['fit'] for r in result['rows']]
            computed=batch_scores(fits,batch)
            for i,s in enumerate(computed):
                saved=result['rows'][i]['methods'][name]
                assert json.dumps(s,sort_keys=True)==json.dumps(saved['score'],sort_keys=True)
                active=batch['active'][i].astype(bool);times=batch['times_ns'][i,active]
                true=dict(first_active_ns=times,emission_ns=times,charge_pe=batch['charge_pe'][i,active])
                tolerance=min(10.,float(np.diff(np.sort(times))[0])/4.) if len(times)==2 else 10.
                strict=score(fits[i],true,[0.,2048.],tolerance_ns=tolerance)
                assert json.dumps(strict,sort_keys=True)==json.dumps(saved['strict_score'],sort_keys=True)
                score_count+=2
            if name.startswith('pair_') and not name.endswith('_refit'):
                with np.load(OUTPUT/'test'/kind/f'{name}-raw.npz') as f:
                    arrays={k:f[k] for k in f.files}
                threshold=load(OUTPUT/'models'/name/'selection.json')['threshold']
                if name.startswith('pair_ordered'):
                    from .alpha_pair_ordered import decode
                    expected,_=decode(arrays,threshold)
                    assert expected==fits
                    decoded+=len(fits)
                    if 'direct' in name:
                        assert np.all(arrays['components_mV']>=0.)
                        assert np.all(arrays['components_mV'].sum(1)<=np.maximum(batch['voltage_mV']-10.,0.)+.003)
                        direct_count+=len(fits)
                    continue
                probability=1/(1+np.exp(-np.clip(arrays['presence_logits'],-40,40)))
                for i,fit in enumerate(fits):
                    ids=np.flatnonzero((probability[i]>=threshold)&(arrays['charge_pe'][i]>=500.))
                    ids=ids[np.argsort(arrays['times_ns'][i,ids])]
                    assert fit['slot_ids']==ids.tolist()
                    assert fit['times_ns']==arrays['times_ns'][i,ids].astype(float).tolist()
                    assert fit['charge_pe']==arrays['charge_pe'][i,ids].astype(float).tolist()
                    assert fit['truth_inputs']==[]
                    decoded+=1
                if name.startswith('pair_direct'):
                    assert np.all(arrays['components_mV']>=0.)
                    positive=np.maximum(batch['voltage_mV']-10.,0.)
                    assert np.all(arrays['components_mV'].sum(1)<=positive+.003)
                    direct_count+=len(fits)
        print('VALIDATED_SCORES',kind,score_count,flush=True)
    # Exact cold replays of three cases in every digitization split. Transport
    # donors are immutable files; new RNG seeds and shared microcells are replayed.
    replays=0
    for role,records in (('train',split['train']),('validation',split['validation']),('test',fresh)):
        paths=sorted((OUTPUT/'photons'/role).glob('*.npz'))
        donors=[]
        for p in paths:
            with np.load(p) as f:donors.append({k:f[k].copy() for k in f.files})
        metadata=load(OUTPUT/f'{role}-shared-pair.json')['cases']
        batch=read_batch(OUTPUT/f'{role}-shared-pair.npz')
        config=load(PROJECT/records[0]['metadata'])['config']
        for i in (0,7,13):
            reconstructed,info=one_case(donors,config,metadata[i]['seed'])
            for key,value in reconstructed.items():
                np.testing.assert_allclose(np.asarray(value,np.float32),batch[key][i],rtol=0.,atol=0.)
            assert info['parents']==metadata[i]['parents']
            replays+=1
    # Independent process, frozen source, no truth file or rate/count metadata.
    b=read_batch(OUTPUT/'assets/test-shared.npz')
    cold_record=RECORD/('validation-'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    cold_record.mkdir()
    voltage_path=cold_record/'inference-voltage-only.npz'
    np.savez_compressed(voltage_path,voltage_mV=b['voltage_mV'][[0,17]])
    cold_runs=[]
    for name in ('pair_direct_262050001','pair_physical_262050001','pair_ordered_direct_262050001','pair_ordered_physical_262050001'):
        destination=cold_record/f'cold-{name}'
        code='''from pathlib import Path
import os,sys
import numpy as np
project=Path(sys.argv[1]);waveform=Path(sys.argv[2]);model=sys.argv[3];output=sys.argv[4]
os.environ['GAGG_PAIR_PROJECT']=str(project)
sys.path.insert(0,str(project/'pileup/experiments/alpha_pair_curriculum_20261021/inference-source'))
original=np.load
def observed_only(path,*args,**kwargs):
    if Path(path).resolve()!=waveform.resolve():
        raise AssertionError('Inference attempted an unapproved array read: '+str(path))
    return original(path,*args,**kwargs)
np.load=observed_only
from pileup.analysis.alpha_pair_predict import main
sys.argv=['frozen-inference','--model',model,'--waveforms',str(waveform),'--output',output]
main()
'''
        result=subprocess.run([sys.executable,'-c',code,str(PROJECT),str(voltage_path),name,str(destination)],
                              cwd=PROJECT,capture_output=True,text=True,timeout=60)
        assert result.returncode==0,result.stderr
        saved=load(OUTPUT/'test/shared/results.json')['rows']
        actual=load(destination/'events.json')['fits']
        for pred,i in zip(actual,(0,17)):
            old=saved[i]['methods'][name]['fit']
            assert len(pred['times_ns'])==len(old['times_ns'])
            np.testing.assert_allclose(pred['times_ns'],old['times_ns'],rtol=1e-5,atol=.001)
            np.testing.assert_allclose(pred['charge_pe'],old['charge_pe'],rtol=1e-5,atol=.05)
        cold_runs.append(dict(model=name,source='Frozen inference-source; only waveform NPZ reads permitted',cases=[0,17],passed=True))
    # Known-single scores are separate and do not turn into unknown-count claims.
    batch=read_batch(OUTPUT/'assets/test-native-single.npz')
    for r in load(OUTPUT/'single-test.json')['rows']:
        assert r['known_K']==1 and r['targets']==len(batch['active'])
        assert json.dumps(batch_scores(r['fits'],batch),sort_keys=True)==json.dumps(r['scores'],sort_keys=True)
    report=dict(status='passed',scoring_rows_recomputed=score_count,decoded_raw_predictions=decoded,
        direct_mask_conservation_cases=direct_count,shared_digitization_exact_replays=replays,
        cold_blind_inference=cold_runs,grouped_roots_disjoint=True,new_root_count=24,
        all_native_macro_seeds_and_hashes=True,frozen_files_checked=len(frozen['files']),
        known_K_diagnostic_separate=True,strict_gap_metric_evaluation_only=True,
        new_test_tuning=False,unit_tests='3 passed before production; engineering16-case overfit passed',
        limitations='Local0/1/2 task; native isolated vs offline shared controls distinct; not20us high-rate certification')
    save(OUTPUT/'validation.json',report)
    print('VALIDATION_PASSED',score_count,decoded,replays,flush=True)


if __name__=='__main__':validate()
