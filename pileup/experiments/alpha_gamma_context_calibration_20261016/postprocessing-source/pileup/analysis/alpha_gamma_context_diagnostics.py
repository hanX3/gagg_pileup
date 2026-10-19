"""Post-fit failure stratification; all truth-derived features are evaluation only."""
from pathlib import Path
import argparse,json
import numpy as np
from .alpha_rate_campaign import save,digest
from .alpha_gamma_context_calibration import preceding_voltage


def diagnose(folder,artifacts):
    folder,artifacts=Path(folder),Path(artifacts)
    selected=json.loads((folder/'frozen-fit-settings.json').read_text())['selected_method']
    methods=['baseline',selected];events=[]
    for method in methods:
        paths=sorted((artifacts/'test_mixed'/method).glob('*.json'));assert len(paths)==16
        for path in paths:
            row=json.loads(path.read_text());stem=Path(row['row']['root']).stem
            truth=dict(np.load(artifacts/'inputs/truth'/f'{stem}.npz'))
            wave=dict(np.load(artifacts/'inputs/waveforms'/f'{stem}.npz'))
            fit=row['fit'];pairs={i:j for i,j in row['score']['all_matched_pairs']}
            eligible=(truth['charge_pe']>=50)&np.isfinite(truth['first_active_ns'])
            targets=np.flatnonzero((truth['charge_pe']>=500)&(truth['emission_ns']>=5000)&(truth['emission_ns']<15000))
            for i in targets:
                q=float(truth['charge_pe'][i]);t0=float(truth['first_active_ns'][i]);particle=str(truth['particle_name'][i])
                others=eligible.copy();others[i]=False
                gap=float(np.min(abs(truth['first_active_ns'][others]-t0))) if others.any() else None
                history=preceding_voltage(wave['time_ns'],wave['analog_mV'],t0)
                if i not in pairs:outcome='missed'
                elif fit['particle_name'][pairs[i]]!=particle:outcome='wrong_type'
                elif abs(fit['charge_pe'][pairs[i]]/q-1)>.2:outcome='bad_charge'
                else:outcome='recovered'
                events.append(dict(root=row['row']['root'],method=method,truth_index=int(i),particle=particle,
                    charge_pe=q,nearest_other_q50_ns=gap,preceding_voltage_mV=history,outcome=outcome,
                    charge_bin='500–1000' if q<1000 else '1000–2000' if q<2000 else '≥2000',
                    spacing_bin='<50' if gap is not None and gap<50 else '50–150' if gap is not None and gap<150 else '150–500' if gap is not None and gap<500 else '≥500',
                    history_bin='<50' if history<50 else '50–200' if history<200 else '≥200'))
            for p in ('alpha','gamma'):
                used=[e for e in events if e['root']==row['row']['root'] and e['method']==method and e['particle']==p]
                assert sum(e['outcome']=='recovered' for e in used)==row['score']['by_particle'][p]['recovered_count']
    groups=[];tables={}
    for feature,bins,unit in [('charge_bin',['500–1000','1000–2000','≥2000'],'p.e.'),
                              ('spacing_bin',['<50','50–150','150–500','≥500'],'ns'),
                              ('history_bin',['<50','50–200','≥200'],'mV')]:
        lines=[f'| 类别 | 分组 / {unit} | 原V3恢复 | 新标定恢复 | 漏检：旧/新 | 错分型：旧/新 | 仅Q失败：旧/新 |',
               '| --- | --- | --- | --- | --- | --- | --- |']
        for p in ('alpha','gamma'):
            for b in bins:
                pair=[]
                for method in methods:
                    subset=[e for e in events if e['method']==method and e['particle']==p and e[feature]==b]
                    r=dict(feature=feature,bin=b,unit=unit,method=method,particle=p,n=len(subset),
                        **{k:sum(e['outcome']==k for e in subset) for k in ('recovered','missed','wrong_type','bad_charge')})
                    groups.append(r);pair.append(r)
                a,c=pair;assert a['n']==c['n']
                lines.append(f"| {p} | {b} | {a['recovered']}/{a['n']} | {c['recovered']}/{c['n']} | {a['missed']}/{c['missed']} | {a['wrong_type']}/{c['wrong_type']} | {a['bad_charge']}/{c['bad_charge']} |")
        tables[feature]='\n'.join(lines)
    result=dict(groups=groups,events=events,tables=tables,source_sha256=digest(Path(__file__)),
        interpretation='Descriptive post-fit associations, not independent causal effects or resolution bounds; no retuning',
        truth_inputs='True charge/class/onset for bins and scores; nearest other real pulse with Q>=50 pe across full window',
        background_feature='Observed voltage in [-24,-4) ns relative to true onset, evaluation only',
        frozen_fit_unchanged=True)
    save(artifacts/'test-failure-bins.json',result)
    print(tables['charge_bin']);print(tables['spacing_bin'])
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder',type=Path,required=True);p.add_argument('--artifacts',type=Path,required=True)
    a=p.parse_args();diagnose(a.folder,a.artifacts)
