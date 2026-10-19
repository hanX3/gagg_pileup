"""Evaluation-only purity, threshold migration and failure accounting."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from .alpha_rate_campaign import save
from .alpha_gamma_pileup_fit import score_mixed


def diagnose(output,inputs):
    output,inputs=Path(output),Path(inputs)
    rows=json.loads((output/'window-results.json').read_text());extras=[];selection=[];failures=[];signal_balance=[]
    for row in rows:
        stem=Path(row['row']['root']).stem;truth=dict(np.load(inputs/'truth'/f'{stem}.npz'))
        fit=row['fit'];score=row['score'];gate=[5000.,15000.]
        by_pred={j:i for i,j in score['all_matched_pairs']}
        by_truth={i:j for i,j in score['matched_pairs']}
        for j in score['unmatched_prediction_indices']:
            distance=abs(truth['first_active_ns']-fit['times_ns'][j])
            candidates=np.flatnonzero(np.isfinite(distance)&(distance<=10.)&(truth['charge_pe']>0))
            above=candidates[truth['charge_pe'][candidates]>=500.]
            tag='near_above_threshold_truth' if len(above) else ('near_subthreshold_truth' if len(candidates) else 'no_truth_within_10_ns')
            nearest=int(candidates[np.argmin(distance[candidates])]) if len(candidates) else None
            extras.append(dict(root=row['row']['root'],gamma_model=row['row']['gamma_model'],
                alpha_injection_MHz=row['row']['rate_hz']/1e6,prediction_index=j,predicted_type=fit['particle_name'][j],
                predicted_time_ns=fit['times_ns'][j],predicted_charge_pe=fit['charge_pe'][j],proximity_flag=tag,
                nearest_truth_type=str(truth['particle_name'][nearest]) if nearest is not None else None,
                nearest_truth_charge_pe=float(truth['charge_pe'][nearest]) if nearest is not None else None,
                nearest_truth_dt_ns=float(distance[nearest]) if nearest is not None else None))
        for particle in ('alpha','gamma'):
            central=(truth['emission_ns']>=gate[0])&(truth['emission_ns']<gate[1])&(truth['particle_name']==particle)
            weak=central&(truth['charge_pe']<500.)
            strong=central&(truth['charge_pe']>=500.)
            signal_balance.append(dict(gamma_model=row['row']['gamma_model'],alpha_injection_MHz=row['row']['rate_hz']/1e6,
                particle=particle,total_charge_pe=float(truth['charge_pe'][central].sum()),
                subthreshold_charge_pe=float(truth['charge_pe'][weak].sum()),
                subthreshold_nonzero_count=int(np.sum(weak&(truth['charge_pe']>0))),
                target_charge_pe=truth['charge_pe'][strong].tolist(),
                target_edep_MeV=truth['edep_MeV'][strong].tolist()))
            predicted=[j for j,t in enumerate(fit['times_ns']) if gate[0]<=t<gate[1] and fit['particle_name'][j]==particle]
            correct=[j for j in predicted if j in by_pred and truth['particle_name'][by_pred[j]]==particle]
            wrong=[j for j in predicted if j in by_pred and truth['particle_name'][by_pred[j]]!=particle]
            selection.append(dict(root=row['row']['root'],gamma_model=row['row']['gamma_model'],
                alpha_injection_MHz=row['row']['rate_hz']/1e6,particle=particle,predicted=len(predicted),
                correct_type_matches=len(correct),wrong_type_matches=len(wrong),unmatched=len(predicted)-len(correct)-len(wrong)))
        for i in score['target_indices']:
            j=by_truth.get(i);reason='missed'
            if j is not None:
                type_ok=fit['particle_name'][j]==truth['particle_name'][i]
                charge_ok=abs(fit['charge_pe'][j]/truth['charge_pe'][i]-1)<=.2
                reason='success' if type_ok and charge_ok else ('type_and_charge' if not type_ok and not charge_ok else 'type_only' if not type_ok else 'charge_only')
            failures.append(dict(root=row['row']['root'],gamma_model=row['row']['gamma_model'],alpha_injection_MHz=row['row']['rate_hz']/1e6,
                particle=str(truth['particle_name'][i]),truth_index=i,reason=reason,truth_charge_pe=float(truth['charge_pe'][i])))
    groups=[]
    for model in ('brems','mono662'):
        for rate in (1.,10.):
            for particle in ('alpha','gamma'):
                s=[v for v in selection if (v['gamma_model'],v['alpha_injection_MHz'],v['particle'])==(model,rate,particle)]
                ex=[v for v in extras if (v['gamma_model'],v['alpha_injection_MHz'],v['predicted_type'])==(model,rate,particle)]
                fail=[v for v in failures if (v['gamma_model'],v['alpha_injection_MHz'],v['particle'])==(model,rate,particle)]
                balance=[v for v in signal_balance if (v['gamma_model'],v['alpha_injection_MHz'],v['particle'])==(model,rate,particle)]
                group=dict(gamma_model=model,alpha_injection_MHz=rate,particle=particle,
                    **{k:sum(v[k] for v in s) for k in ('predicted','correct_type_matches','wrong_type_matches','unmatched')})
                group['selection_purity']=group['correct_type_matches']/group['predicted'] if group['predicted'] else None
                group['extra_proximity']={k:sum(v['proximity_flag']==k for v in ex) for k in ('near_above_threshold_truth','near_subthreshold_truth','no_truth_within_10_ns')}
                group['failure_partition']={k:sum(v['reason']==k for v in fail) for k in ('success','missed','type_only','charge_only','type_and_charge')}
                group['signal_balance']={k:sum(v[k] for v in balance) for k in ('total_charge_pe','subthreshold_charge_pe','subthreshold_nonzero_count')}
                total=group['signal_balance']['total_charge_pe']
                group['signal_balance']['subthreshold_charge_fraction']=group['signal_balance']['subthreshold_charge_pe']/total if total else None
                for source,destination in (('target_charge_pe','target_charge_quantiles_pe'),
                                           ('target_edep_MeV','target_edep_quantiles_MeV')):
                    values=[x for v in balance for x in v[source]]
                    group['signal_balance'][destination]=np.quantile(values,[.16,.5,.84]).tolist() if values else None
                groups.append(group)
    result=dict(groups=groups,extra_candidates=extras,selection_windows=selection,target_outcomes=failures,
        purity_definition='Correct particle class among predictions with fitted time in [5,15) us; truths Q>=500 pe matched across full window',
        proximity_note='10 ns proximity flags are diagnostics, not proof of causal threshold migration or splitting')
    save(output/'selection-diagnostics.json',result);return result


def shifted_time_control(output,inputs):
    """Keep reconstructed charges/classes/counts; destroy central time alignment.

    This is a descriptive accidental-match check, not a calibrated null model
    or a p value. Real guard predictions are retained at their original times.
    """
    output,inputs=Path(output),Path(inputs)
    windows=json.loads((output/'window-results.json').read_text());rows=[]
    for window in windows:
        truth=dict(np.load(inputs/'truth'/f'{Path(window['row']['root']).stem}.npz'))
        original=window['fit'];times=np.array(original['times_ns']);central=(times>=5000)&(times<15000)
        for shift in (1250.,2750.,4250.,5750.,7250.):
            changed=times.copy();changed[central]=5000.+(changed[central]-5000.+shift)%10000.
            fit=dict(times_ns=changed.tolist(),charge_pe=original['charge_pe'],particle_name=original['particle_name'])
            scored=score_mixed(fit,truth,[5000.,15000.])
            rows.append(dict(root=window['row']['root'],gamma_model=window['row']['gamma_model'],
                alpha_injection_MHz=window['row']['rate_hz']/1e6,shift_ns=shift,
                by_particle={p:{k:v[k] for k in ('target_count','matched_count','correct_type_count','recovered_count')}
                             for p,v in scored['by_particle'].items()}))
    groups=[]
    for model in ('brems','mono662'):
        for rate in (1.,10.):
            for particle in ('alpha','gamma'):
                selected=[r for r in rows if (r['gamma_model'],r['alpha_injection_MHz'])==(model,rate)]
                denominator=sum(r['by_particle'][particle]['target_count'] for r in selected)
                groups.append(dict(gamma_model=model,alpha_injection_MHz=rate,particle=particle,
                    shifted_detection_efficiency=sum(r['by_particle'][particle]['matched_count'] for r in selected)/denominator if denominator else None,
                    shifted_full_recovery_efficiency=sum(r['by_particle'][particle]['recovered_count'] for r in selected)/denominator if denominator else None))
    result=dict(groups=groups,rows=rows,definition=shifted_time_control.__doc__,shifts_ns=[1250.,2750.,4250.,5750.,7250.])
    save(output/'shifted-time-control.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True)
    args=p.parse_args();diagnose(args.output,args.inputs);shifted_time_control(args.output,args.inputs)
