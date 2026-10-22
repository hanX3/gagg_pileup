"""Truth-only descriptive pulse-shape check; never a blind fit or resolution bound."""
from pathlib import Path
import argparse
import json
import numpy as np
from .alpha_rate_campaign import digest, save


def diagnose(folder, inputs, calibration, output):
    folder,inputs,calibration,output=map(Path,(folder,inputs,calibration,output))
    if output.exists():raise FileExistsError(output)
    cal=json.loads(calibration.read_text());records=[];hashes={}
    for row in json.loads((folder/'manifest.json').read_text())['runs']:
        if row['role'] not in ('calibration','validation_mixed'):continue
        path=inputs/'truth'/f"{Path(row['root']).stem}.npz";truth=dict(np.load(path));hashes[str(path)]=digest(path)
        at=truth['avalanche_time_ns'];aq=truth['avalanche_charge_C'];owner=truth['avalanche_owner']
        for i,q in enumerate(truth['charge_pe']):
            if q<500 or not 5000<=truth['emission_ns'][i]<15000:continue
            mask=owner==i;t=at[mask]-truth['first_active_ns'][i];charge=aq[mask]
            # Moment identity of a 95/351 ns two-exponential mixture. This is
            # descriptive and retains values outside [0,1] if fluctuations occur.
            mean=float(np.average(t,weights=charge));f_moment=(mean-95.)/(351.-95.)
            p=cal['classes'][str(truth['particle_name'][i])]
            sq=np.clip(q,*p['charge_range_pe']);pred=p['slow_intercept']+p['slow_log_slope']*np.log(sq/1000.)
            records.append(dict(root=row['root'],role=row['role'],particle=str(truth['particle_name'][i]),
                truth_index=i,charge_pe=float(q),charge_weighted_mean_ns=mean,slow_fraction_moment=f_moment,
                calibrated_slow=float(pred),moment_residual=float(f_moment-pred),
                gamma_model=row.get('gamma_model'),alpha_rate_hz=row.get('rate_hz'),
                interpretation='Uses avalanche ancestry and true onset/charge; population diagnosis only'))
    groups=[]
    for role in ('calibration','validation_mixed'):
        for particle in ('alpha','gamma'):
            group=[r for r in records if r['role']==role and r['particle']==particle]
            f=np.array([r['slow_fraction_moment'] for r in group]);d=np.array([r['moment_residual'] for r in group])
            groups.append(dict(role=role,particle=particle,n=len(group),
                slow_p16_p50_p84=np.quantile(f,[.16,.5,.84]).tolist(),
                residual_p16_p50_p84=np.quantile(d,[.16,.5,.84]).tolist()))
    result=dict(records=records,groups=groups,input_hashes=hashes,calibration_sha256=digest(calibration),
        source_sha256=digest(Path(__file__)),truth_inputs=['avalanche ancestry, charge and time','event true first-photon onset, charge and class'],
        limitations=['Moment estimate is not the voltage least-squares estimator',
            'Finite tail, onset delay, transport and occupancy can affect the moment',
            'Only threshold targets in disclosed calibration/development data; not blind performance'])
    save(output,result)
    for row in groups:print('SHAPE',row,flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('folder','inputs','calibration','output'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();diagnose(args.folder,args.inputs,args.calibration,args.output)
