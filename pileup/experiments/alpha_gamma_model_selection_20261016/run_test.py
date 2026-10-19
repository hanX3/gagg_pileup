"""Compare three frozen fitters using at most four independent subprocesses."""
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
import os
import subprocess
import sys
import time

for variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[variable]='1'
FOLDER=Path(__file__).resolve().parent
sys.path.insert(0,str(FOLDER/'source'))
sys.path.insert(1,str(FOLDER.parents[2]))


def run_one(label,method,inputs,output):
    import numpy as np
    from pileup.analysis.alpha_rate_campaign import digest,save
    from pileup.analysis.alpha_gamma_model_selection_fit import fit_model_selection,ModelSelectionConfig
    from pileup.analysis.alpha_gamma_regularized_fit import fit_regularized,RegularizedConfig
    from pileup.analysis.alpha_gamma_study import calibration as old_calibration
    from pileup.analysis import alpha_gamma_pileup_fit as original
    from pileup.analysis.alpha_gamma_nnls import reduced_nnls
    frozen=json.loads((FOLDER/'frozen-fit-settings.json').read_text())
    row=next(r for r in json.loads((FOLDER/'manifest.json').read_text())['runs'] if r['label']==label)
    assert row['status']=='completed' and row['role']=='test_mixed'
    inputs,output=Path(inputs),Path(output);stem=Path(row['root']).stem
    source_name={'baseline':'alpha_gamma_pileup_fit.py','regularized_3':'alpha_gamma_regularized_fit.py',
        'model_selection':'alpha_gamma_model_selection_fit.py'}[method]
    source=FOLDER/'source/pileup/analysis'/source_name
    assert digest(source)==frozen['source_sha256'][str(source.relative_to(FOLDER))]
    calpath=FOLDER/'calibration.json';assert digest(calpath)==frozen['calibration_sha256']
    name=method if method!='model_selection' else f"model_selection_{frozen['config']['addition_cost']:g}"
    wp=inputs/'waveforms'/f'{stem}.npz';tp=inputs/'truth'/f'{stem}.npz'
    dest=output/name/f'{stem}.json'
    if dest.exists():
        previous=json.loads(dest.read_text());assert previous['fitter_sha256']==digest(source)
        assert previous['waveform_sha256']==digest(wp)
        print('REUSE',label,name,flush=True);return
    cal=json.loads(calpath.read_text());wave=dict(np.load(wp));start=time.perf_counter()
    if method=='model_selection':
        fit=fit_model_selection(wave['time_ns'],wave['analog_mV'],cal,ModelSelectionConfig(**frozen['config']))
    elif method=='regularized_3':
        fit=fit_regularized(wave['time_ns'],wave['analog_mV'],cal,RegularizedConfig(**frozen['reference_config']))
    else:
        original.nnls=reduced_nnls
        _,bank,config=old_calibration(FOLDER/'baseline-settings.json')
        fit=original.fit_waveform(wave['time_ns'],wave['analog_mV'],bank,config)
    seconds=time.perf_counter()-start
    # Test truth is opened after reconstruction and never passed to a fitter.
    truth=dict(np.load(tp));meta=json.loads((inputs/'metadata'/f'{stem}.json').read_text())
    score=original.score_mixed(fit,truth,meta['gate_ns']);model=fit.pop('model_mV')
    dest.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(dest.with_suffix('.npz'),time_ns=wave['time_ns'],model_mV=model)
    save(dest,dict(row=row,method=name,fit=fit,score=score,fit_seconds=seconds,
        fitter_sha256=digest(source),calibration_sha256=digest(calpath),
        baseline_settings_sha256=digest(FOLDER/'baseline-settings.json'),
        waveform_sha256=digest(wp),truth_sha256=digest(tp),gate_ns=meta['gate_ns'],fit_truth_inputs=[],
        launcher_sha256=digest(Path(__file__))))
    print('FIT',label,name,{p:(s['recovered_count'],s['target_count']) for p,s in score['by_particle'].items()},
        'extra',score['unmatched_predictions'],'seconds',round(seconds,2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,choices=range(1,5),default=4)
    parser.add_argument('--method',choices=('all','baseline','regularized_3','model_selection'),default='all')
    parser.add_argument('--label');args=parser.parse_args()
    if args.label:
        assert args.method!='all';run_one(args.label,args.method,args.inputs,args.output);return
    rows=json.loads((FOLDER/'manifest.json').read_text())['runs']
    assert len(rows)==8 and all(r['status']=='completed' for r in rows)
    methods=('baseline','regularized_3','model_selection') if args.method=='all' else (args.method,)
    def launch(label,method):
        command=[sys.executable,str(Path(__file__)), '--label',label,'--method',method,
                 '--inputs',str(args.inputs),'--output',str(args.output)]
        done=subprocess.run(command,capture_output=True,text=True)
        if done.returncode:raise RuntimeError(done.stdout+done.stderr)
        return done.stdout.strip()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs=[pool.submit(launch,r['label'],m) for r in rows for m in methods]
        for future in as_completed(jobs):print(future.result(),flush=True)


if __name__=='__main__':main()
