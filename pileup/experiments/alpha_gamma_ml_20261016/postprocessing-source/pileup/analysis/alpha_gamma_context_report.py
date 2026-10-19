"""Reports retain all conditions/classes and paired window uncertainty."""
from pathlib import Path
import argparse,json
import numpy as np
from .alpha_rate_campaign import digest,save
from .alpha_gamma_regularized_report import metrics


def read_rows(path):
    return [json.loads(p.read_text()) for p in sorted(Path(path).glob('*.json'))]


def summarize(folder,artifacts,development_inputs):
    folder,artifacts,development_inputs=map(Path,(folder,artifacts,development_inputs))
    frozen=json.loads((folder/'frozen-fit-settings.json').read_text());selected=frozen['selected_method']
    methods=['baseline',selected];groups=[];pooled=[];diagnostics={};seconds={};test_rows={}
    for method in methods:
        rows=read_rows(artifacts/'test_mixed'/method);assert len(rows)==16;test_rows[method]=rows
        for model in ('brems','mono662'):
            for rate in (1e6,1e7):
                subset=[r for r in rows if r['row']['gamma_model']==model and r['row']['rate_hz']==rate];assert len(subset)==4
                for particle in ('alpha','gamma'):
                    groups.append(dict(method=method,gamma_model=model,alpha_MHz=rate/1e6,particle=particle,windows=4,
                                       **metrics(subset,particle,artifacts/'inputs')))
        pooled.extend(dict(method=method,particle=p,**metrics(rows,p,artifacts/'inputs')) for p in ('alpha','gamma'))
        diagnostics[method]={k:sum(r['score'][k] for r in rows) for k in
            ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags')}
        seconds[method]=sum(r['pipeline_seconds'] for r in rows)
    development=[]
    for c in frozen['validation_choices']:
        rows=read_rows(artifacts/'development'/c['method'])
        development.append(dict(method=c['method'],macro_full_f1=c['macro_full_f1'],
            by_particle={p:metrics(rows,p,development_inputs) for p in ('alpha','gamma')},
            extras=sum(r['score']['unmatched_predictions'] for r in rows)))
    # Paired bootstrap of complete independent windows, stratified by condition.
    # This retains within-window correlations and the fixed condition counts.
    rng=np.random.default_rng(261069999);draws={p:[] for p in ('alpha','gamma')}
    base={r['row']['root']:r for r in test_rows['baseline']};new={r['row']['root']:r for r in test_rows[selected]}
    strata=[[root for root,r in base.items() if r['row']['gamma_model']==m and r['row']['rate_hz']==rate]
            for m in ('brems','mono662') for rate in (1e6,1e7)]
    for _ in range(4000):
        roots=[root for s in strata for root in rng.choice(s,len(s),replace=True)]
        for p in draws:
            n=sum(base[root]['score']['by_particle'][p]['target_count'] for root in roots)
            if n:
                delta=sum(new[root]['score']['by_particle'][p]['recovered_count']-base[root]['score']['by_particle'][p]['recovered_count'] for root in roots)/n
                draws[p].append(delta)
    bootstrap={p:dict(valid_draws=len(v),delta_recovery_p025_p50_p975=np.quantile(v,[.025,.5,.975]).tolist()) for p,v in draws.items()}
    result=dict(selected_method=selected,methods=methods,groups=groups,pooled=pooled,diagnostics=diagnostics,
        pipeline_seconds=seconds,development=development,paired_window_bootstrap=bootstrap,
        bootstrap_seed=261069999,bootstrap_draws=4000,
        interpretation='Small fixed simulated sample; bootstrap covers window variability, not model/systematic uncertainty',
        target_definition='Emission [5,15) us, realized Q>=500 pe; time-only one-to-one match <=10 ns; correct class and |Qfit/Qtruth-1|<=20%',
        fit_truth_inputs=[],test_tuning=False,source_sha256=digest(Path(__file__)))
    save(artifacts/'summary.json',result)
    tables={}
    lines=['| 开发方案 | 等权完整恢复F1 | α完整恢复 | γ完整恢复 | 额外预测 |','| --- | ---: | ---: | ---: | ---: |']
    for d in development:
        a,g=(d['by_particle'][p] for p in ('alpha','gamma'))
        lines.append(f"| {d['method']} | {d['macro_full_f1']:.5f} | {a['recovered_count']}/{a['target_count']} | {g['recovered_count']}/{g['target_count']} | {d['extras']} |")
    tables['development']='\n'.join(lines)
    lines=['| γ来源 | α源/MHz | 类别 | 原V3 | 新标定候选 |','| --- | ---: | --- | --- | --- |']
    for model in ('brems','mono662'):
        for rate in (1.,10.):
            for p in ('alpha','gamma'):
                cells=[next(g for g in groups if (g['method'],g['gamma_model'],g['alpha_MHz'],g['particle'])==(m,model,rate,p)) for m in methods]
                lines.append(f'| {model} | {rate:g} | {p} | '+' | '.join(f"{g['recovered_count']}/{g['target_count']} ({g['full_recovery']:.1%})" for g in cells)+' |')
    tables['test']='\n'.join(lines)
    lines=['| 方法 | 类别 | 完整恢复 | 类型纯度 | 完整预测比例 | Q偏差/标准差 | 漏检 | 错分型 |','| --- | --- | --- | --- | --- | --- | ---: | ---: |']
    for g in pooled:
        lines.append(f"| {g['method']} | {g['particle']} | {g['recovered_count']}/{g['target_count']} | {g['pure_prediction_count']}/{g['prediction_count']} | {g['precise_prediction_count']}/{g['prediction_count']} | {g['charge_bias']:.1%} / {g['charge_std']:.1%} | {g['missed_count']} | {g['wrong_type_count']} |")
    tables['pooled']='\n'.join(lines)
    lines=['| 方法 | 目标 | 时间匹配 | 漏检 | 额外预测 | 误拆分标志 | 误合并标志 |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for m in methods:lines.append('| '+m+' | '+' | '.join(str(v) for v in diagnostics[m].values())+' |')
    tables['failures']='\n'.join(lines)
    save(artifacts/'tables.json',tables)
    return result


def figures(folder,artifacts):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder,artifacts=Path(folder),Path(artifacts)
    data=json.loads((artifacts/'calibration-observations.json').read_text())['rows']
    cal=json.loads((folder/'calibration.json').read_text());old=json.loads((folder/'old-calibration.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for k,p in enumerate(('alpha','gamma')):
        isolated=[r for r in data if r['used'] and r['particle']==p and r['role']=='calibration_single']
        ax=axes[0,k]
        for pos,color,label in [([0,0,0],'#2474a4','center'),([4,0,0],'#dc8b23','edge'),([4,4,0],'#9e5588','corner')]:
            group=[r for r in isolated if r['position_mm']==pos]
            ax.scatter([r['charge_pe'] for r in group],[r['slow_fraction'] for r in group],s=24,c=color,label=label,alpha=.8)
        q=np.geomspace(50,20000,300)
        for c,style,label in [(old,'--','old mean'),(cal,'-','expanded mean')]:
            cc=c['classes'][p];qq=np.clip(q,*cc['charge_range_pe'])
            ax.plot(q,np.clip(cc['slow_intercept']+cc['slow_log_slope']*np.log(qq/1000.),.02,.98),style,color='black',label=label)
        ax.axvline(500,color='grey',ls=':',lw=1);ax.set_xscale('log');ax.set_ylim(0,1)
        ax.set_title(p+': isolated native calibration');ax.set_xlabel('Realized charge / p.e.');ax.set_ylabel('Fitted slow fraction')
        ax.legend(fontsize=7,ncol=2);ax.grid(alpha=.2)
        ax=axes[1,k];cc=cal['context_classes'][p]
        history=[r for r in data if r['used'] and r['is_probe'] and r['particle']==p and r['role']=='calibration_history']
        for pos,color,label in [([0,0,0],'#2474a4','probe center'),([4,0,0],'#dc8b23','probe edge')]:
            group=[r for r in history if r['position_mm']==pos]
            residual=[r['slow_fraction']-cc['slow_intercept']-cc['slow_log_slope']*np.log(np.clip(r['charge_pe'],*cc['charge_range_pe'])/1000.) for r in group]
            ax.scatter([r['preceding_voltage_mV'] for r in group],residual,s=28,c=color,label=label)
        h=np.linspace(0,cc['history_range_mV'][1],100);ax.plot(h,cc['slow_per_100mV']*h/100.,color='black')
        ax.set_title(p+': controlled-history calibration');ax.set_xlabel('Observed pre-onset level / mV');ax.set_ylabel('Fraction after charge trend removal')
        ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Simulation calibration; ancestry used only to isolate labelled component shapes',fontsize=11)
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(artifacts/f'calibration-coverage.{ext}',dpi=150,bbox_inches='tight')
    plt.close(fig)
    s=json.loads((artifacts/'summary.json').read_text());fig,axes=plt.subplots(2,2,figsize=(11,7.5))
    for i,metric in enumerate(('full_recovery','full_precision')):
        for j,p in enumerate(('alpha','gamma')):
            ax=axes[i,j]
            for k,m in enumerate(s['methods']):
                gs=[next(g for g in s['groups'] if (g['method'],g['particle'],g['gamma_model'],g['alpha_MHz'])==(m,p,model,rate))
                    for model in ('brems','mono662') for rate in (1.,10.)]
                x=np.arange(4)+(k-.5)*.32;ys=[100*(g[metric] or 0.) for g in gs]
                ax.bar(x,ys,.30,label=m,color=['#2474a4','#dc8b23'][k])
                for xx,y,g in zip(x,ys,gs):
                    num=g['recovered_count'] if i==0 else g['precise_prediction_count'];den=g['target_count'] if i==0 else g['prediction_count']
                    ax.text(xx,y+2+4*k,f'{num}/{den}',ha='center',fontsize=7)
            ax.set_xticks(range(4),['Brems\n1 MHz','Brems\n10 MHz','662 keV\n1 MHz','662 keV\n10 MHz'])
            ax.set_ylim(0,115);ax.set_yticks([0,25,50,75,100]);ax.set_ylabel('Percent');ax.grid(axis='y',alpha=.2)
            ax.set_title(p+(' full recovery' if i==0 else ' full precision'))
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='upper center',ncol=2);fig.tight_layout(rect=[0,0,1,.95])
    for ext in ('png','pdf'):fig.savefig(artifacts/f'independent-test.{ext}',dpi=150,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ('folder','artifacts','development-inputs'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();s=summarize(a.folder,a.artifacts,a.development_inputs);figures(a.folder,a.artifacts)
    print('SUMMARY',s['selected_method'],[(r['method'],r['particle'],r['recovered_count'],r['target_count']) for r in s['pooled']],flush=True)
