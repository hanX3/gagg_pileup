"""Development evaluation of regularized within-class shape variation."""
from pathlib import Path
import argparse,json,time
import numpy as np
from .alpha_gamma_shape_profile import fit_shape_profile,ShapeConfig
from .alpha_gamma_pileup_fit import score_mixed
from .alpha_rate_campaign import digest,save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('initial','inputs','calibration','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--scale',type=float,required=True);p.add_argument('--label');a=p.parse_args()
    cal=json.loads(a.calibration.read_text());source=Path(__file__).with_name('alpha_gamma_shape_profile.py')
    for path in sorted(a.initial.glob('*.json')):
        ref=json.loads(path.read_text());row=ref['row']
        if a.label and row['label']!=a.label:continue
        assert ref['calibration_sha256']==digest(a.calibration)
        assert ref['fitter_sha256']==digest(Path(__file__).with_name('alpha_gamma_model_selection_fit.py'))
        wp=a.inputs/'waveforms'/f'{path.stem}.npz';tp=a.inputs/'truth'/f'{path.stem}.npz'
        assert digest(wp)==ref['waveform_sha256']
        method=f'shape_{a.scale:g}';dest=a.output/method/path.name
        if dest.exists():
            old=json.loads(dest.read_text())
            assert old['fitter_sha256']==digest(source) and old['initializer_sha256']==digest(path)
            assert old['calibration_sha256']==digest(a.calibration) and old['waveform_sha256']==digest(wp)
            continue
        wave=dict(np.load(wp));start=time.perf_counter()
        initial=dict(ref['fit'])
        # JSON omits the saved waveform model. Empty fits return their pilot,
        # so restore this voltage-derived array before calling the fitter.
        initial['model_mV']=np.load(path.with_suffix('.npz'))['model_mV']
        fit=fit_shape_profile(wave['time_ns'],wave['analog_mV'],cal,ShapeConfig(sigma_scale=a.scale),initial=initial)
        seconds=time.perf_counter()-start
        truth=dict(np.load(tp));score=score_mixed(fit,truth,ref['gate_ns']);model=fit.pop('model_mV')
        dest.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(dest.with_suffix('.npz'),time_ns=wave['time_ns'],model_mV=model)
        save(dest,dict(row=row,method=method,fit=fit,score=score,fit_seconds=seconds,fitter_sha256=digest(source),
            calibration_sha256=digest(a.calibration),initializer_sha256=digest(path),waveform_sha256=digest(wp),
            truth_sha256=digest(tp),fit_truth_inputs=[],gate_ns=ref['gate_ns']))
        print('SHAPE',row['label'],a.scale,{p:(s['recovered_count'],s['target_count']) for p,s in score['by_particle'].items()},
            'extra',score['unmatched_predictions'],'s',round(seconds,1),flush=True)


if __name__=='__main__':main()
