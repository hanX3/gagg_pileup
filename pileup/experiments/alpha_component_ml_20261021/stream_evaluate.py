"""Evaluate completed predeclared windows while remaining G4 jobs run.

Only orchestration/serialization is added after the freeze. Predictions use the
unchanged frozen infer() and model assets; every planned window is retained.
The original evaluate stage assembles/rechecks the full dataset afterwards.
"""
from pathlib import Path
import json
import os
import sys
import time

project=Path(__file__).resolve().parents[3]
os.environ['GAGG_COMPONENT_PROJECT']=str(project)
sys.path.insert(0,str(Path(__file__).resolve().parent/'inference-source'))
import numpy as np
from pileup.analysis import alpha_component_evaluate as evaluation
from pileup.analysis.alpha_rate_campaign import save,digest
from pileup.analysis.alpha_pileup_fit import score


def main():
    evaluation.check_freeze();nets=evaluation.networks()
    out=evaluation.OUTPUT;experiment=evaluation.EXPERIMENT
    for name in ('windows','models','dense'):(out/'test'/name).mkdir(parents=True,exist_ok=True)
    plan=json.loads((experiment/'assets/test-plan.json').read_text())['planned_runs']
    expected={r['label'] for r in plan};done=set()
    source_sha=digest(evaluation.__file__);driver_sha=digest(__file__)
    while done!=expected:
        try:manifest=json.loads((experiment/'manifest.json').read_text())
        except json.JSONDecodeError:
            time.sleep(2.);continue
        assert manifest['planned_runs']==plan
        if any(r['status'] in ('failed','interrupted') for r in manifest['runs']):
            raise RuntimeError('Native production has a failed/interrupted window; preserve it for diagnosis')
        pending=[r for r in manifest['runs'] if r['status']=='completed' and r['label'] not in done]
        for row in pending:
            stem=Path(row['root']).stem;destination=out/'test/windows'/(stem+'.json')
            if destination.exists():
                saved=json.loads(destination.read_text());assert saved['root_sha256']==row['root_sha256']
                done.add(row['label']);continue
            wp=out/'inputs/waveforms'/(stem+'.npz');tp=out/'inputs/truth'/(stem+'.npz')
            with np.load(wp) as f:t,y=f['time_ns'].copy(),f['analog_mV'].copy()
            fits,seconds,denses=evaluation.infer(t,y,nets)
            # Separate scoring export is opened only after all blind predictions.
            with np.load(tp) as f:truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe')}
            methods={}
            for name,fit in fits.items():
                if 'model_mV' in fit:
                    np.savez_compressed(out/'test/models'/f'{stem}_{name}.npz',model_mV=fit.pop('model_mV'))
                methods[name]=dict(fit=fit,score=score(fit,truth,[5000.,15000.]),seconds=seconds[name])
            for name,dense in denses.items():np.savez_compressed(out/'test/dense'/f'{stem}_{name}.npz',**dense)
            save(destination,dict(root=row['root'],root_sha256=row['root_sha256'],seed=row['seed'],rate_hz=row['rate_hz'],
                waveform=str(wp.relative_to(project)),truth=str(tp.relative_to(project)),gate_ns=[5000.,15000.],
                waveform_sha256=digest(wp),truth_sha256=digest(tp),methods=methods,
                inference_truth_inputs=[],rate_input_at_inference=False,
                frozen_inference_source_sha256=source_sha,streaming_driver_sha256=driver_sha))
            done.add(row['label'])
            print('STREAMED',len(done),'/',len(plan),row['rate_hz'],row['seed'],flush=True)
        if done!=expected:time.sleep(5.)
    save(project/'jupyter/run_records/alpha_component_ml_20261021/streaming-evaluation.json',dict(
        status='completed',windows=len(done),frozen_inference_source_sha256=source_sha,
        driver_sha256=driver_sha,parameter_or_seed_changes=0,scoring_after_inference=True))
    print('STREAMING_EVALUATION_COMPLETE',flush=True)


if __name__=='__main__':main()
