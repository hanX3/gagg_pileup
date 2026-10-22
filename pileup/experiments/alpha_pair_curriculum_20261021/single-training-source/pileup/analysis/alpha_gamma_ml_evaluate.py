"""Voltage-only inference, separate scoring, validation choices and test freeze."""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
from datetime import datetime
import numpy as np
import torch
from .alpha_rate_campaign import save,digest
from .alpha_gamma_ml import voltage_features,predict_features,decode
from .alpha_gamma_ml_study import project_root,load_model,summarize_scores
from .alpha_gamma_ml_hybrid import fit_hybrid,MixedHybridConfig
from .alpha_gamma_model_selection_fit import fit_model_selection,ModelSelectionConfig
from .alpha_gamma_pileup_fit import score_mixed
from .alpha_gamma_regularized_fit import CalibratedBank


def rows_for(folder,output,role):
    root=project_root();folder,output=Path(folder).resolve(),Path(output).resolve()
    if role=='test':
        assert (folder/'frozen-fit-settings.json').exists()
        m=json.loads((folder/'manifest.json').read_text())
        rows=[r for r in m['runs'] if r['role']=='test_mixed']
        assert len(rows)==16 and all(r['status']=='completed' for r in rows)
        prepared=[]
        for r in rows:
            stem=Path(r['root']).stem
            prepared.append(dict(r,split='test',**{k:str((output/'inputs'/sub/f'{stem}.{ext}').relative_to(root))
                for k,sub,ext in [('waveform','waveforms','npz'),('truth','truth','npz'),('metadata','metadata','json')]}))
        return prepared
    path=output/'data-split.json'
    if not path.exists():path=folder/'development-split.json'
    return [r for r in json.loads(path.read_text())['windows'] if r['split']==role]


def algorithm_hashes(method):
    names=['alpha_gamma_pileup_fit.py','alpha_pileup_fit.py']
    if method=='baseline' or method.startswith('hybrid'):
        names+=['alpha_gamma_model_selection_fit.py','alpha_gamma_regularized_fit.py','alpha_gamma_nnls.py']
    if method!='baseline':names+=['alpha_gamma_ml.py','alpha_gamma_ml_study.py']
    if method.startswith('hybrid'):names+=['alpha_gamma_ml_hybrid.py']
    return {n:digest(Path(__file__).with_name(n)) for n in names}


def evaluate(folder,output,role,methods,shard=(0,1)):
    root=project_root();folder,output=Path(folder).resolve(),Path(output).resolve();torch.set_num_threads(1)
    cal=json.loads((folder/'calibration.json').read_text());cal_hash=digest(folder/'calibration.json')
    networks={};requested=list(methods)
    for method in methods:
        kind=method.split('_')[1] if method.startswith('hybrid') else method
        if kind in ('cnn','unet') and kind not in networks:networks[kind]=load_model(output/kind)
    if any(m.startswith('hybrid') for m in methods) and 'baseline' not in requested:requested.insert(0,'baseline')
    for row in rows_for(folder,output,role)[shard[0]::shard[1]]:
        stem=Path(row['root']).stem;wavepath=root/row['waveform'];wave=dict(np.load(wavepath))
        t,y=wave['time_ns'],wave['analog_mV'];wave_hash=digest(wavepath)
        fits={};seconds={};dense={};neural={};provenance={};features=None
        for kind,(model,cfg,frozen) in networks.items():
            begin=time.perf_counter()
            if features is None:features=voltage_features(t,y,cal)
            dense[kind]=predict_features(model,features,cfg)
            neural[kind]=decode(t,dense[kind],frozen['selected']['thresholds'],cfg)
            seconds[kind]=time.perf_counter()-begin;fits[kind]=neural[kind]
            densepath=output/'dense'/role/kind/f'{stem}.npz';densepath.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(densepath,**dense[kind])
        for method in requested:
            destination=output/'results'/role/method/f'{stem}.json'
            if destination.exists():
                cache=json.loads(destination.read_text())
                assert cache['waveform_sha256']==wave_hash and cache['calibration_sha256']==cal_hash
                assert cache['algorithm_sha256']==algorithm_hashes(method)
                fits[method]=cache['fit'];seconds[method]=cache['fit_seconds']
                if destination.with_suffix('.npz').exists():
                    fits[method]=dict(fits[method],model_mV=np.load(destination.with_suffix('.npz'))['model_mV'])
                continue
            start=time.perf_counter()
            if method=='baseline':
                # Existing expanded fits were produced from voltage only, and
                # are reused only after waveform/calibration/fitter hash checks.
                candidates=[root/f'jupyter/artifacts/alpha_gamma_context_calibration_20261016/{part}/expanded/{stem}.json'
                            for part in ('development','test_mixed')]
                cache_path=next((p for p in candidates if p.exists()),None)
                cached=json.loads(cache_path.read_text()) if cache_path else None
                if (cached and cached['waveform_sha256']==wave_hash and cached['calibration_sha256']==cal_hash and
                    cached['source_sha256']['alpha_gamma_model_selection_fit.py']==algorithm_hashes(method)['alpha_gamma_model_selection_fit.py']):
                    fits[method]=dict(cached['fit'],model_mV=np.load(cache_path.with_suffix('.npz'))['model_mV'])
                    seconds[method]=cached['fit_seconds'];provenance[method]=dict(cache=str(cache_path.relative_to(root)),sha256=digest(cache_path))
                else:
                    fits[method]=fit_model_selection(t,y,cal,ModelSelectionConfig(addition_cost=1.))
                    seconds[method]=time.perf_counter()-start
            elif method.startswith('hybrid'):
                _,kind,strength=method.split('_')
                fits[method]=fit_hybrid(t,y,cal,neural[kind],dense[kind],fits['baseline'],MixedHybridConfig(strength=float(strength)))
                seconds[method]=time.perf_counter()-start
            elif method not in networks:raise ValueError(method)
        # Access labels only after all requested blind methods are finished.
        truth=dict(np.load(root/row['truth']));meta=json.loads((root/row['metadata']).read_text())
        for method in requested:
            destination=output/'results'/role/method/f'{stem}.json'
            if destination.exists():continue
            destination.parent.mkdir(parents=True,exist_ok=True)
            fit=dict(fits[method]);model_voltage=fit.pop('model_mV',None)
            if model_voltage is not None:np.savez_compressed(destination.with_suffix('.npz'),model_mV=model_voltage,time_ns=t)
            s=score_mixed(fit,truth,meta['gate_ns'])
            pipeline=seconds[method]
            if method.startswith('hybrid'):pipeline+=seconds['baseline']+seconds[method.split('_')[1]]
            save(destination,dict(row=row,method=method,fit=fit,score=s,gate_ns=meta['gate_ns'],
                fit_seconds=seconds[method],pipeline_seconds=pipeline,calibration_sha256=cal_hash,
                waveform_sha256=wave_hash,truth_sha256=digest(root/row['truth']),algorithm_sha256=algorithm_hashes(method),
                model_hashes={k:f['model_sha256'] for k,(_,_,f) in networks.items()},
                initialization_cache=provenance.get(method),fit_truth_inputs=[],
                inference_label_access='All reconstructions finished before opening truth file'))
            print('FIT',role,row['label'],method,{p:(v['recovered_count'],v['target_count']) for p,v in s['by_particle'].items()},
                  'extra',s['unmatched_predictions'],'seconds',round(seconds[method],2),flush=True)


def select_and_freeze(folder,output):
    root=project_root();folder,output=Path(folder).resolve(),Path(output).resolve()
    assert not (folder/'frozen-fit-settings.json').exists()
    neural={k:json.loads((output/k/'frozen-model.json').read_text()) for k in ('cnn','unet')}
    selected_kind=max(neural,key=lambda k:(neural[k]['selected']['objective'],k=='cnn'))
    choices=[]
    for strength in (0,10,100):
        method=f'hybrid_{selected_kind}_{strength}'
        paths=sorted((output/'results/validation'/method).glob('*.json'));assert len(paths)==20
        rows=[]
        for p in paths:
            r=json.loads(p.read_text());r['truth']=dict(np.load(root/r['row']['truth']));rows.append(r)
        choices.append(dict(method=method,strength=strength,**summarize_scores(rows),
            result_sha256={str(p.relative_to(root)):digest(p) for p in paths}))
    chosen=max(choices,key=lambda x:(x['objective'],-x['strength']))
    saved=dict(created=datetime.now().astimezone().isoformat(),state='frozen_before_new_test_generation',
        neural_models=neural,hybrid_kind=selected_kind,hybrid_method=chosen['method'],
        hybrid_config=asdict_config(chosen['strength']),hybrid_choices=choices,
        common_methods=['baseline','cnn','unet',chosen['method']],
        calibration_sha256=digest(folder/'calibration.json'),data_split_sha256=digest(output/'data-split.json'),
        validation_rule='20 disclosed ROOTs; alpha/gamma and four source conditions equal-weight full recovery F1',
        test_labels_seen=False,regression_used_for_selection=False,
        implementation_manifest_sha256=digest(folder/'ml-source-manifest.json'))
    save(folder/'frozen-fit-settings.json',saved)
    print('ALL METHODS FROZEN',saved['common_methods'],[(c['strength'],c['objective']) for c in choices],flush=True)


def asdict_config(strength):
    from dataclasses import asdict
    return asdict(MixedHybridConfig(strength=strength))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=('evaluate','freeze'))
    p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--role',choices=('validation','test','regression','control'),default='validation')
    p.add_argument('--methods',nargs='+',default=['baseline']);p.add_argument('--shard',nargs=2,type=int,default=[0,1])
    a=p.parse_args()
    if a.stage=='freeze':select_and_freeze(a.folder,a.output)
    else:evaluate(a.folder,a.output,a.role,a.methods,a.shard)


if __name__=='__main__':main()
