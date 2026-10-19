"""Fit completed fresh windows with frozen code; score only after both fits."""
from pathlib import Path
from datetime import datetime
import json,sys,time
FOLDER=Path(__file__).resolve().parent
PROJECT=FOLDER.parents[2]
sys.path.insert(0,str(FOLDER/'source'))
import numpy as np
from pileup.analysis.alpha_gamma_model_selection_fit import fit_model_selection,ModelSelectionConfig
from pileup.analysis.alpha_gamma_shape_profile import fit_shape_profile,ShapeConfig
from pileup.analysis.alpha_gamma_pileup_fit import score_mixed
from pileup.analysis.alpha_rate_campaign import save,digest


def main():
    artifacts=PROJECT/'jupyter/artifacts/alpha_gamma_profile_fit_20261016';inputs=artifacts/'inputs'
    frozen=json.loads((FOLDER/'frozen-fit-settings.json').read_text());calpath=FOLDER/'calibration.json';cal=json.loads(calpath.read_text())
    assert digest(calpath)==frozen['calibration_sha256']
    for name,key in (('alpha_gamma_model_selection_fit.py','initializer_fitter_sha256'),('alpha_gamma_shape_profile.py','fitter_sha256')):
        assert digest(FOLDER/'source/pileup/analysis'/name)==frozen[key]
    save(FOLDER/'test-driver.json',dict(created=datetime.now().astimezone().isoformat(),source_sha256=digest(Path(__file__)),
        selected_config=frozen['config'],calibration_sha256=digest(calpath),truth_opened_after_both_fits=True))
    while True:
        manifest=json.loads((FOLDER/'manifest.json').read_text());done=0
        assert all(r['status'] in ('running','completed') for r in manifest['runs'])
        for row in manifest['runs']:
            if row['status']!='completed' or 'validation' not in row:continue
            stem=Path(row['root']).stem
            paths=[artifacts/'test_mixed'/m/f'{stem}.json' for m in ('model_selection_1','shape_0.5')]
            if all(p.exists() for p in paths):done+=1;continue
            assert not any(p.exists() for p in paths),'Inspect incomplete saved pair before resuming'
            wp=inputs/'waveforms'/f'{stem}.npz';wave=dict(np.load(wp));start=time.perf_counter()
            initial=fit_model_selection(wave['time_ns'],wave['analog_mV'],cal,ModelSelectionConfig(addition_cost=1.))
            initial_seconds=time.perf_counter()-start;start=time.perf_counter()
            new=fit_shape_profile(wave['time_ns'],wave['analog_mV'],cal,ShapeConfig(**frozen['config']),initial=initial)
            new_seconds=time.perf_counter()-start
            # Scores and simulation truth are deliberately accessed after inference.
            tp=inputs/'truth'/f'{stem}.npz';truth=dict(np.load(tp));meta=json.loads((inputs/'metadata'/f'{stem}.json').read_text())
            for path,fit,seconds,key in zip(paths,(initial,new),(initial_seconds,new_seconds),('initializer_fitter_sha256','fitter_sha256')):
                fit=dict(fit);model=fit.pop('model_mV');score=score_mixed(fit,truth,meta['gate_ns'])
                path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path.with_suffix('.npz'),time_ns=wave['time_ns'],model_mV=model)
                record=dict(row=row,method=path.parent.name,fit=fit,score=score,fit_seconds=seconds,
                    fitter_sha256=frozen[key],calibration_sha256=digest(calpath),waveform_sha256=digest(wp),truth_sha256=digest(tp),
                    fit_truth_inputs=[],gate_ns=meta['gate_ns'],pipeline_driver_sha256=digest(Path(__file__)))
                if path.parent.name=='shape_0.5':record['initializer_sha256']=digest(paths[0])
                save(path,record)
                print('TEST',row['label'],path.parent.name,{p:(s['recovered_count'],s['target_count']) for p,s in score['by_particle'].items()},
                    'extra',score['unmatched_predictions'],'s',round(seconds,2),flush=True)
            done+=1
        if done==len(manifest['planned_runs']):break
        time.sleep(2.)


if __name__=='__main__':main()
