"""Blind comparison driver. Open event truth only after all fits are complete."""
import argparse
import importlib
import json
from pathlib import Path
import time
import numpy as np
from .alpha_rate_campaign import digest, save
from .alpha_gamma_model_selection_fit import fit_model_selection, ModelSelectionConfig
from .alpha_gamma_context_fit import fit_context
from .alpha_gamma_shape_profile import ShapeConfig
from .alpha_gamma_pileup_fit import score_mixed


def evaluate(folder,inputs,calibration,old_calibration,output,role,methods,shard=(0,1),baseline_cache=None):
    folder,inputs,calibration,old_calibration,output=map(Path,(folder,inputs,calibration,old_calibration,output))
    rows=[r for r in json.loads((folder/'manifest.json').read_text())['runs'] if r['role']==role]
    assert rows and all(r['status']=='completed' for r in rows)
    cal=json.loads(calibration.read_text());old=json.loads(old_calibration.read_text())
    hashes={name:digest(Path(importlib.import_module('.'+name[:-3],__package__).__file__)) for name in ('alpha_gamma_context_study.py','alpha_gamma_context_fit.py',
        'alpha_gamma_context_calibration.py','alpha_gamma_shape_profile.py','alpha_gamma_model_selection_fit.py',
        'alpha_gamma_regularized_fit.py','alpha_gamma_nnls.py')}
    for row in rows[shard[0]::shard[1]]:
        stem=Path(row['root']).stem;wavepath=inputs/'waveforms'/f'{stem}.npz'
        wave=dict(np.load(wavepath));fits={};times={};cache_info={}
        wanted=list(methods)
        if any(m.startswith('context_') for m in wanted) and 'expanded' not in wanted:wanted.insert(0,'expanded')
        for method in wanted:
            destination=output/method/f'{stem}.json'
            if destination.exists():
                record=json.loads(destination.read_text())
                assert record['calibration_sha256']==digest(old_calibration if method=='baseline' else calibration)
                assert record['waveform_sha256']==digest(wavepath) and record['source_sha256']==hashes
                fits[method]=dict(record['fit'],model_mV=np.load(destination.with_suffix('.npz'))['model_mV'])
                times[method]=record['fit_seconds'];continue
            start=time.perf_counter()
            if method=='baseline' and baseline_cache:
                path=Path(baseline_cache)/f'{stem}.json';r=json.loads(path.read_text())
                assert r['calibration_sha256']==digest(old_calibration) and r['waveform_sha256']==digest(wavepath)
                assert r['fitter_sha256']==hashes['alpha_gamma_model_selection_fit.py']
                fits[method]=dict(r['fit'],model_mV=np.load(path.with_suffix('.npz'))['model_mV'])
                times[method]=r['fit_seconds'];cache_info[method]=dict(path=str(path),sha256=digest(path))
            elif method in ('baseline','expanded'):
                fits[method]=fit_model_selection(wave['time_ns'],wave['analog_mV'],old if method=='baseline' else cal,
                    ModelSelectionConfig(addition_cost=1.))
                times[method]=time.perf_counter()-start
            else:
                scale=float(method.split('_')[1])
                fits[method]=fit_context(wave['time_ns'],wave['analog_mV'],cal,ShapeConfig(sigma_scale=scale),initial=fits['expanded'])
                times[method]=time.perf_counter()-start
        # Label access occurs after voltage-only inference (or immutable cache).
        truthpath=inputs/'truth'/f'{stem}.npz';truth=dict(np.load(truthpath))
        meta=json.loads((inputs/'metadata'/f'{stem}.json').read_text())
        for method in wanted:
            destination=output/method/f'{stem}.json'
            if destination.exists():continue
            fit=dict(fits[method]);model=fit.pop('model_mV');score=score_mixed(fit,truth,meta['gate_ns'])
            destination.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(destination.with_suffix('.npz'),time_ns=wave['time_ns'],model_mV=model)
            save(destination,dict(row=row,method=method,fit=fit,score=score,fit_seconds=times[method],
                pipeline_seconds=times[method]+(times['expanded'] if method.startswith('context_') else 0.),
                initialization_cache=cache_info.get(method),source_sha256=hashes,
                calibration_sha256=digest(old_calibration if method=='baseline' else calibration),
                waveform_sha256=digest(wavepath),truth_sha256=digest(truthpath),gate_ns=meta['gate_ns'],fit_truth_inputs=[]))
            print('FIT',row['label'],method,{p:(v['recovered_count'],v['target_count']) for p,v in score['by_particle'].items()},
                  'extra',score['unmatched_predictions'],'seconds',round(times[method],2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('folder','inputs','calibration','old-calibration','output'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--role',required=True);p.add_argument('--methods',nargs='+',default=['baseline','expanded','context_0.5','context_1'])
    p.add_argument('--shard',nargs=2,type=int,default=[0,1]);p.add_argument('--baseline-cache',type=Path)
    a=p.parse_args();evaluate(a.folder,a.inputs,a.calibration,a.old_calibration,a.output,a.role,a.methods,a.shard,a.baseline_cache)
