"""Evaluate leading-edge stabilization with truth opened only after fitting."""
from pathlib import Path
import argparse,json,time
import numpy as np
from .alpha_gamma_stable_profile import fit_stable,StableConfig
from .alpha_gamma_pileup_fit import score_mixed
from .alpha_rate_campaign import digest,save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('initial','inputs','calibration','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--scale',type=float,required=True);p.add_argument('--label');args=p.parse_args()
    cal=json.loads(args.calibration.read_text());source=Path(__file__).with_name('alpha_gamma_stable_profile.py')
    for path in sorted(args.initial.glob('*.json')):
        r=json.loads(path.read_text());row=r['row']
        if args.label and row['label']!=args.label:continue
        assert r['calibration_sha256']==digest(args.calibration)
        assert r['fitter_sha256']==digest(Path(__file__).with_name('alpha_gamma_model_selection_fit.py'))
        wp=args.inputs/'waveforms'/f'{path.stem}.npz';tp=args.inputs/'truth'/f'{path.stem}.npz'
        assert digest(wp)==r['waveform_sha256']
        method=f'stable_{args.scale:g}';dest=args.output/method/path.name
        if dest.exists():raise FileExistsError(dest)
        wave=dict(np.load(wp));start=time.perf_counter()
        initial=dict(r['fit']);initial['model_mV']=np.load(path.with_suffix('.npz'))['model_mV']
        fit=fit_stable(wave['time_ns'],wave['analog_mV'],cal,StableConfig(relative_jump_scale=args.scale),initial=initial)
        seconds=time.perf_counter()-start;truth=dict(np.load(tp));score=score_mixed(fit,truth,r['gate_ns']);model=fit.pop('model_mV')
        dest.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(dest.with_suffix('.npz'),time_ns=wave['time_ns'],model_mV=model)
        save(dest,dict(row=row,method=method,fit=fit,score=score,fit_seconds=seconds,fitter_sha256=digest(source),
            profile_solver_sha256=digest(Path(__file__).with_name('alpha_gamma_profile_fit.py')),
            calibration_sha256=digest(args.calibration),initializer_sha256=digest(path),
            waveform_sha256=digest(wp),truth_sha256=digest(tp),fit_truth_inputs=[],gate_ns=r['gate_ns']))
        print('STABLE',row['label'],args.scale,{p:(s['recovered_count'],s['target_count']) for p,s in score['by_particle'].items()},
            'extra',score['unmatched_predictions'],'s',round(seconds,1),flush=True)


if __name__=='__main__':main()
