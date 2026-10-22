"""Common alpha/gamma metrics and paired-window uncertainty, without tuning."""
from pathlib import Path
import argparse,json
from collections import defaultdict
import numpy as np
from .alpha_rate_campaign import save,digest


def root_path():
    return next(p for p in [Path.cwd(),*Path.cwd().parents] if (p/'g4/data').is_dir())


def metrics(rows,particle,root):
    keys=('target_count','matched_count','correct_type_count','recovered_count','missed_count','wrong_type_count','central_injected','central_nonzero_charge')
    counts={k:sum(r['score']['by_particle'][particle][k] for r in rows) for k in keys}
    predictions=pure=precise=0;dq=[];dt=[];correct_dq=[]
    for r in rows:
        s=r['score']['by_particle'][particle];dq+=s['relative_charge_errors'];dt+=s['time_errors_ns']
        truth=dict(np.load(root/r['row']['truth']));fit=r['fit'];gate=r['gate_ns']
        mask=np.array([name==particle and gate[0]<=t<gate[1] for name,t in zip(fit['particle_name'],fit['times_ns'])],bool)
        predictions+=int(mask.sum())
        for i,j in r['score']['all_matched_pairs']:
            if mask[j] and truth['particle_name'][i]==particle:
                pure+=1;precise+=abs(fit['charge_pe'][j]/truth['charge_pe'][i]-1)<=.2
        for i,j in r['score']['matched_pairs']:
            if truth['particle_name'][i]==particle and fit['particle_name'][j]==particle:
                correct_dq.append(fit['charge_pe'][j]/truth['charge_pe'][i]-1)
    recall=counts['recovered_count']/counts['target_count'] if counts['target_count'] else None
    precision=precise/predictions if predictions else 0.
    f1=2*precision*recall/(precision+recall) if recall is not None and precision+recall else 0.
    def moments(x):
        return dict(n=len(x),bias=float(np.mean(x)) if len(x) else None,
            std=float(np.std(x,ddof=1)) if len(x)>1 else None,
            p16_p50_p84=np.quantile(x,[.16,.5,.84]).tolist() if len(x) else None)
    return dict(**counts,prediction_count=predictions,pure_prediction_count=pure,precise_prediction_count=int(precise),
        full_recovery=recall,type_purity=pure/predictions if predictions else None,full_precision=precision,full_f1=f1,
        charge_error_all_time_matches=moments(dq),charge_error_correct_type_time_matches=moments(correct_dq),
        time_rmse_ns=float(np.sqrt(np.mean(np.square(dt)))) if dt else None,
        unmatched_predictions=sum(r['score']['unmatched_predictions_by_class'][particle] for r in rows))


def summarize(folder,output,role='test'):
    root=root_path();folder,output=Path(folder),Path(output)
    frozen=json.loads((folder/'frozen-fit-settings.json').read_text());methods=frozen['common_methods']
    all_rows={m:[json.loads(p.read_text()) for p in sorted((output/'results'/role/m).glob('*.json'))] for m in methods}
    expected=16 if role in ('test','regression') else 43
    assert all(len(rows)==expected for rows in all_rows.values())
    groups=[];pooled=[];diagnostics={};timing={}
    for method,rows in all_rows.items():
        assert len({r['row']['root'] for r in rows})==expected
        pooled.extend(dict(method=method,particle=p,**metrics(rows,p,root)) for p in ('alpha','gamma'))
        if role!='control':
            for model in ('brems','mono662'):
                for rate in (1000000,10000000):
                    group=[r for r in rows if r['row']['gamma_model']==model and r['row']['rate_hz']==rate]
                    assert len(group)==4
                    for p in ('alpha','gamma'):
                        groups.append(dict(method=method,gamma_model=model,alpha_MHz=rate/1e6,particle=p,
                            windows=len(group),**metrics(group,p,root)))
        diagnostics[method]={k:sum(r['score'][k] for r in rows) for k in
            ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags')}
        timing[method]=dict(sum_pipeline_seconds=sum(r['pipeline_seconds'] for r in rows),
            median_pipeline_seconds=float(np.median([r['pipeline_seconds'] for r in rows])),
            maximum_pipeline_seconds=max(r['pipeline_seconds'] for r in rows))
    paired={}
    if role!='control':
        base={r['row']['root']:r for r in all_rows['baseline']}
        strata=[[key for key,r in base.items() if r['row']['gamma_model']==m and r['row']['rate_hz']==rate]
                for m in ('brems','mono662') for rate in (1000000,10000000)]
        for method in methods[1:]:
            other={r['row']['root']:r for r in all_rows[method]};rng=np.random.default_rng(261084000)
            draws={p:[] for p in ('alpha','gamma')}
            for _ in range(2000):
                roots=[str(k) for s in strata for k in rng.choice(s,len(s),replace=True)]
                for p in draws:
                    n=sum(base[k]['score']['by_particle'][p]['target_count'] for k in roots)
                    if n:
                        draws[p].append(sum(other[k]['score']['by_particle'][p]['recovered_count']-
                            base[k]['score']['by_particle'][p]['recovered_count'] for k in roots)/n)
            paired[method]={p:dict(valid_draws=len(x),delta_recovery_p025_p50_p975=np.quantile(x,[.025,.5,.975]).tolist()) for p,x in draws.items()}
    result=dict(role=role,methods=methods,windows=expected,groups=groups,pooled=pooled,diagnostics=diagnostics,timing=timing,
        paired_window_bootstrap=paired,bootstrap_seed=261084000,bootstrap_draws=2000,
        target_definition='Q>=500 pe; truth emission [5,15) us; time-only matching <=10 ns; correct class and |dQ/Q|<=20%',
        timing_scope='Offline complete 20 us records; shared feature computation charged to first neural method; concurrent CPU timings are diagnostic only',
        error_scope='Charge errors conditional on time matching, separately with/without correct class; not deposited-energy resolution',
        limits=['Four independent windows per source condition; limited gamma statistics',
                'Bootstrap describes window variability, not model or calibration systematics',
                'No real detector data or calibrated classification probabilities',
                'ML output min spacing 4 ns is a representation choice, not demonstrated detector resolution'],
        algorithm_frozen_before_new_test=True,source_sha256=digest(Path(__file__)))
    save(output/f'{role}-summary.json',result)
    print('SUMMARY',role,[(r['method'],r['particle'],r['recovered_count'],r['target_count'],r['full_precision']) for r in pooled],flush=True)
    return result


def diagnostics(folder,output,role='test'):
    root=root_path();folder,output=Path(folder),Path(output)
    methods=json.loads((folder/'frozen-fit-settings.json').read_text())['common_methods']
    bins=defaultdict(lambda:dict(n=0,recovered=0,missed=0,wrong_type=0,bad_charge=0));events=[]
    for method in methods:
        for path in sorted((output/'results'/role/method).glob('*.json')):
            r=json.loads(path.read_text());truth=dict(np.load(root/r['row']['truth']));wave=dict(np.load(root/r['row']['waveform']))
            pair=dict(r['score']['matched_pairs']);tt=truth['first_active_ns'];q=truth['charge_pe']
            physical=np.flatnonzero(np.isfinite(tt)&(q>=50.))
            for i in r['score']['target_indices']:
                near=physical[physical!=i];gap=float(np.min(abs(tt[near]-tt[i]))) if len(near) else float('inf')
                mask=(wave['time_ns']>=tt[i]-24)&(wave['time_ns']<tt[i]-4)
                history=max(0.,float(wave['analog_mV'][mask].mean()-10.)) if mask.any() else 0.
                p=str(truth['particle_name'][i]);outcome='missed'
                if i in pair:
                    j=pair[i]
                    outcome=('wrong_type' if r['fit']['particle_name'][j]!=p else
                             'bad_charge' if abs(r['fit']['charge_pe'][j]/q[i]-1)>.2 else 'recovered')
                descriptors=dict(charge='500–1000' if q[i]<1000 else '1000–2000' if q[i]<2000 else '≥2000',
                    spacing='<50' if gap<50 else '50–150' if gap<150 else '150–500' if gap<500 else '≥500',
                    history='<50' if history<50 else '50–200' if history<200 else '≥200')
                events.append(dict(root=r['row']['root'],method=method,particle=p,truth_index=i,charge_pe=float(q[i]),
                    nearest_other_q50_ns=gap if np.isfinite(gap) else None,preceding_voltage_mV=history,outcome=outcome))
                for feature,interval in descriptors.items():
                    count=bins[(method,p,feature,interval)];count['n']+=1;count[outcome]+=1
    result=dict(groups=[dict(method=m,particle=p,feature=f,bin=b,**v) for (m,p,f,b),v in sorted(bins.items())],events=events,
        definition='Post-fit diagnostic: true onsets/classes/Q for strata; observed voltage [-24,-4) ns before truth onset minus 10 mV',
        interpretation='Overlapping descriptive strata, not causal effects; these truth-indexed features are not inference inputs')
    save(output/f'{role}-failure-bins.json',result)
    return result


def figures(output,role='test'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output=Path(output);s=json.loads((output/f'{role}-summary.json').read_text())
    colors=['#666666','#2878b5','#ee8534','#3b9b5f'];fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
    for col,p in enumerate(('alpha','gamma')):
        for row,metric in enumerate(('full_recovery','full_precision')):
            ax=axes[row,col]
            for k,method in enumerate(s['methods']):
                vals=[next(g for g in s['groups'] if (g['method'],g['gamma_model'],g['alpha_MHz'],g['particle'])==(method,m,r,p))
                      for m in ('brems','mono662') for r in (1.,10.)]
                x=np.arange(4)+(k-1.5)*.19;y=[100*(v[metric] or 0.) for v in vals]
                ax.bar(x,y,.18,label=method,color=colors[k])
            ax.set_ylim(0,105);ax.set_ylabel('Percent');ax.set_title(p+' '+metric.replace('_',' '));ax.grid(axis='y',alpha=.2)
            ax.set_xticks(range(4),['Brems\n1 MHz','Brems\n10 MHz','662 keV\n1 MHz','662 keV\n10 MHz'])
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='upper center',ncol=4)
    fig.suptitle(f'{role}: gamma source 10 MHz; alpha source shown below',y=.94,fontsize=11)
    fig.tight_layout(rect=(0,0,1,.91))
    for ext in ('png','pdf'):fig.savefig(output/f'{role}-comparison.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,3.5))
    for kind in ('cnn','unet'):
        h=json.loads((output/kind/'training-history.json').read_text())
        axes[0].plot([r['step'] for r in h],[r['loss'] for r in h],marker='o',label=kind)
        axes[1].plot([r['step'] for r in h],[r['objective'] for r in h],marker='o',label=kind)
    axes[0].set(xlabel='Training step',ylabel='Training batch loss');axes[1].set(xlabel='Training step',ylabel='Validation macro full F1',ylim=(0,1))
    for ax in axes:ax.grid(alpha=.2);ax.legend()
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(output/f'training-history.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--role',choices=('test','regression','control'),default='test');a=p.parse_args()
    summarize(a.folder,a.output,a.role)
    if a.role!='control':diagnostics(a.folder,a.output,a.role);figures(a.output,a.role)
