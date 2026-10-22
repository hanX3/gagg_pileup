"""Explicitly truth-assisted covariance diagnosis on disclosed development."""
from pathlib import Path
import argparse
import json
import numpy as np
from .alpha_gamma_regularized_fit import CalibratedBank
from .alpha_gamma_shot_covariance import ShotWhitening
from .alpha_gamma_nnls import reduced_nnls
from .alpha_gamma_pileup_fit import score_mixed
from .alpha_rate_campaign import digest, save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('folder','inputs','calibration','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    cal=json.loads(args.calibration.read_text());bank=CalibratedBank(cal);rows=[]
    for row in json.loads((args.folder/'manifest.json').read_text())['runs']:
        if row['role']!='validation_mixed':continue
        stem=Path(row['root']).stem;wp=args.inputs/'waveforms'/f'{stem}.npz';tp=args.inputs/'truth'/f'{stem}.npz'
        wave=dict(np.load(wp));truth=dict(np.load(tp));keep=(truth['charge_pe']>0)&np.isfinite(truth['first_active_ns'])
        tt=truth['first_active_ns'][keep];q=truth['charge_pe'][keep];k=(truth['particle_name'][keep]=='gamma').astype(int)
        t=wave['time_ns'][::2];v=wave['analog_mV'][::2]-cal['baseline_mV'];W=ShotWhitening(t,tt,k,q,q,cal);y=W.apply(v)
        cases=[]
        for known in (True,False):
            times,kinds,shapes=(tt,k,q) if known else (np.repeat(tt,2),np.tile([0,1],len(tt)),np.repeat(q,2))
            H=W.apply(bank.matrix(t,times,kinds,shapes));amps=reduced_nnls(H,y,maxiter=10000)[0]
            if known:charges,classes=amps,k
            else:
                dual=amps.reshape(-1,2);charges=dual.sum(axis=1);classes=np.argmax(dual,axis=1)
            report=charges>=500
            fit=dict(times_ns=tt[report].tolist(),charge_pe=charges[report].tolist(),particle_name=bank.names[classes[report]].tolist())
            score=score_mixed(fit,truth,[5000.,15000.]);gradient=H.T@(H@amps-y)
            violation=max(float(np.max(np.maximum(-gradient,0.))),float(np.max(abs(gradient[amps>1e-6]))))/max(1.,float(np.max(abs(H.T@y))))
            assert violation<1e-7
            cases.append(dict(known_types=known,fit=fit,score=score,nnls_relative_kkt=violation))
            print('PHYSICAL ORACLE',row['label'],known,{p:(s['recovered_count'],s['target_count']) for p,s in score['by_particle'].items()},flush=True)
        rows.append(dict(row=row,cases=cases,waveform_sha256=digest(wp),truth_sha256=digest(tp)))
    save(args.output,dict(rows=rows,calibration_sha256=digest(args.calibration),source_sha256=digest(Path(__file__)),
        covariance_sha256=digest(Path(__file__).with_name('alpha_gamma_shot_covariance.py')),
        truth_inputs=['all nonzero event times/counts','true charge for template shape and covariance','true classes for covariance, also fit when known_types=True'],
        interpretation='Development diagnosis only; not blind performance and not a rigorous physical information bound'))


if __name__=='__main__':main()
