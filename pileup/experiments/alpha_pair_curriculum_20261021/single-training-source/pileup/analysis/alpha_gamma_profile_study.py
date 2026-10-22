"""Evaluate covariance/profile variants on disclosed development windows."""
from pathlib import Path
import argparse
from dataclasses import asdict
import json
import time
import numpy as np
from .alpha_gamma_profile_fit import fit_profile, ProfileConfig
from .alpha_gamma_pileup_fit import score_mixed
from .alpha_rate_campaign import digest, save


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('initial','inputs','calibration','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--label')
    parser.add_argument('--weighting',choices=('ar1','diagonal','uniform','shot'),default='ar1')
    parser.add_argument('--time-maxiter',type=int,default=60)
    parser.add_argument('--passes',type=int,default=3)
    args=parser.parse_args()
    cal=json.loads(args.calibration.read_text());config=ProfileConfig(weighting=args.weighting,time_maxiter=args.time_maxiter,passes=args.passes)
    method='profile_'+args.weighting+('_fixed' if args.time_maxiter==0 else '')
    source=Path(__file__).with_name('alpha_gamma_profile_fit.py')
    initial_source=Path(__file__).with_name('alpha_gamma_model_selection_fit.py')
    for path in sorted(args.initial.glob('*.json')):
        reference=json.loads(path.read_text());row=reference['row']
        if args.label and row['label']!=args.label:continue
        assert reference['fitter_sha256']==digest(initial_source)
        assert reference['calibration_sha256']==cal.get('parent_calibration_sha256',digest(args.calibration))
        stem=Path(row['root']).stem;wp=args.inputs/'waveforms'/f'{stem}.npz';tp=args.inputs/'truth'/f'{stem}.npz'
        assert digest(wp)==reference['waveform_sha256']
        dest=args.output/method/f'{stem}.json'
        if dest.exists():
            previous=json.loads(dest.read_text());assert previous['fitter_sha256']==digest(source)
            print('REUSE',dest,flush=True);continue
        wave=dict(np.load(wp));started=time.perf_counter()
        # Only the voltage-derived fit is passed, never its evaluation score.
        initial=dict(reference['fit']);initial['model_mV']=np.load(path.with_suffix('.npz'))['model_mV']
        fit=fit_profile(wave['time_ns'],wave['analog_mV'],cal,config,initial=initial)
        elapsed=time.perf_counter()-started
        truth=dict(np.load(tp));score=score_mixed(fit,truth,reference['gate_ns'])
        model=fit.pop('model_mV');dest.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(dest.with_suffix('.npz'),time_ns=wave['time_ns'],model_mV=model)
        save(dest,dict(row=row,method=method,fit=fit,score=score,fit_seconds=elapsed,
            fitter_sha256=digest(source),calibration_sha256=digest(args.calibration),
            initializer_sha256=digest(path),initializer_fitter_sha256=digest(initial_source),
            initializer_config=reference['fit'].get('config'),waveform_sha256=digest(wp),truth_sha256=digest(tp),
            fit_truth_inputs=[],gate_ns=reference['gate_ns'],config=asdict(config)))
        print('PROFILE',row['label'],args.weighting,{p:(s['recovered_count'],s['target_count']) for p,s in score['by_particle'].items()},
            'extra',score['unmatched_predictions'],'s',round(elapsed,1),flush=True)


if __name__=='__main__':main()
