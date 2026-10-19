"""Controls and development validation for joint alpha/gamma reconstruction."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from .alpha_gamma_model_selection_fit import fit_model_selection,ModelSelectionConfig
from .alpha_gamma_pileup_fit import score_mixed
from .alpha_rate_campaign import save,digest


def evaluate(folder,inputs,calibration,output,role,penalties):
    folder,inputs,calibration,output=map(Path,(folder,inputs,calibration,output))
    cal=json.loads(calibration.read_text());manifest=json.loads((folder/'manifest.json').read_text())
    rows=[r for r in manifest['runs'] if r['role']==role]
    source=Path(__file__).with_name('alpha_gamma_model_selection_fit.py');code_hash=digest(source)
    for row in rows:
        stem=Path(row['root']).stem;wave=dict(np.load(inputs/'waveforms'/f'{stem}.npz'))
        for penalty in penalties:
            name=f'model_selection_{penalty:g}';dest=output/name/f'{stem}.json'
            if dest.exists():
                old=json.loads(dest.read_text());assert old['fitter_sha256']==code_hash
                continue
            start=time.perf_counter()
            fit=fit_model_selection(wave['time_ns'],wave['analog_mV'],cal,ModelSelectionConfig(addition_cost=penalty))
            seconds=time.perf_counter()-start
            truth=dict(np.load(inputs/'truth'/f'{stem}.npz'))
            meta=json.loads((inputs/'metadata'/f'{stem}.json').read_text())
            score=score_mixed(fit,truth,meta['gate_ns']);model=fit.pop('model_mV')
            dest.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(dest.with_suffix('.npz'),time_ns=wave['time_ns'],model_mV=model)
            save(dest,dict(row=row,method=name,fit=fit,score=score,fit_seconds=seconds,
                fitter_sha256=code_hash,calibration_sha256=digest(calibration),
                waveform_sha256=digest(inputs/'waveforms'/f'{stem}.npz'),truth_sha256=digest(inputs/'truth'/f'{stem}.npz'),
                gate_ns=meta['gate_ns'],fit_truth_inputs=[]))
            print('FIT',row['label'],name,{p:(s['recovered_count'],s['target_count']) for p,s in score['by_particle'].items()},
                'extra',score['unmatched_predictions'],'seconds',round(seconds,2),'KKT',(fit.get('selection') or {}).get('projected_gradient_max'),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('folder','inputs','calibration','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--role',default='validation_mixed');p.add_argument('--penalties',type=float,nargs='+',default=[1.,9.])
    a=p.parse_args();evaluate(a.folder,a.inputs,a.calibration,a.output,a.role,a.penalties)
