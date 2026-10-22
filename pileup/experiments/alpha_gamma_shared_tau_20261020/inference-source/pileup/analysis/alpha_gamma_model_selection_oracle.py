"""Truth-assisted development diagnosis, explicitly not blind performance."""
from pathlib import Path
import argparse
import json
import numpy as np
from .alpha_gamma_regularized_fit import CalibratedBank,FrozenWhitening
from .alpha_gamma_nnls import reduced_nnls
from .alpha_gamma_pileup_fit import score_mixed
from .alpha_rate_campaign import digest,save


def diagnose(folder,inputs,calibration,output):
    folder,inputs,calibration,output=map(Path,(folder,inputs,calibration,output))
    cal=json.loads(calibration.read_text());bank=CalibratedBank(cal);rows=[]
    for row in json.loads((folder/'manifest.json').read_text())['runs']:
        if row['role']!='validation_mixed':continue
        stem=Path(row['root']).stem;wp=inputs/'waveforms'/f'{stem}.npz';tp=inputs/'truth'/f'{stem}.npz'
        wave=dict(np.load(wp));truth=dict(np.load(tp))
        keep=(truth['charge_pe']>0)&np.isfinite(truth['first_active_ns'])
        tt=truth['first_active_ns'][keep];q=truth['charge_pe'][keep];typ=(truth['particle_name'][keep]=='gamma').astype(int)
        t=wave['time_ns'][::2];v=wave['analog_mV'][::2]-10.
        pilot=bank.matrix(t,tt,typ,q)@q;W=FrozenWhitening(pilot,4.,cal['noise'],smooth_ns=0.)
        target=W.apply(v);cases=[]
        for known in (True,False):
            if known:times,kinds,shapes=tt,typ,q
            else:times,kinds,shapes=np.repeat(tt,2),np.tile([0,1],len(tt)),np.repeat(q,2)
            H=W.apply(bank.matrix(t,times,kinds,shapes));amps=reduced_nnls(H,target,maxiter=10000)[0]
            if known:charge=amps;classes=typ
            else:dual=amps.reshape(-1,2);charge=dual.sum(axis=1);classes=np.argmax(dual,axis=1)
            report=charge>=500
            fit=dict(times_ns=tt[report].tolist(),charge_pe=charge[report].tolist(),particle_name=bank.names[classes[report]].tolist())
            score=score_mixed(fit,truth,[5000.,15000.])
            gradient=H.T@(H@amps-target);scale=max(1.,float(np.max(abs(H.T@target))))
            violation=max(float(np.max(np.maximum(-gradient,0.))),float(np.max(abs(gradient[amps>1e-6]))))/scale
            assert violation<1e-7
            cases.append(dict(known_types=known,score=score,fit=fit,weighted_cost=float(np.sum((target-H@amps)**2)),
                nnls_relative_kkt_violation=violation))
            print('ORACLE',row['label'],known,{p:(s['recovered_count'],s['target_count']) for p,s in score['by_particle'].items()},flush=True)
        rows.append(dict(row=row,cases=cases,waveform_sha256=digest(wp),truth_sha256=digest(tp)))
    result=dict(rows=rows,truth_inputs=['all nonzero counts and first_active times',
        'true charge for template shape and pilot covariance','true type for covariance; also reconstruction types when known_types=True'],
        purpose='Deliberately truth-assisted development diagnosis; never blind performance or a rigorous physical upper bound',
        calibration_sha256=digest(calibration),source_sha256=digest(Path(__file__)))
    save(output,result);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('folder','inputs','calibration','output'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();diagnose(args.folder,args.inputs,args.calibration,args.output)
