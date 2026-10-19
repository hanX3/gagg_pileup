"""Compare all disclosed profile trials and the frozen shape-fit test."""
from pathlib import Path
import argparse,json
import numpy as np
from .alpha_gamma_regularized_report import metrics
from .alpha_rate_campaign import save


def summarize(folder,artifacts,previous,development_inputs):
    folder,artifacts,previous,development_inputs=map(Path,(folder,artifacts,previous,development_inputs))
    frozen=json.loads((folder/'frozen-fit-settings.json').read_text());selected=frozen['selected_method']
    development=[]
    for choice in frozen['validation_choices']:
        method=choice['method'];directory=previous/'validation_mixed'/method if method=='model_selection_1' else artifacts/'development'/method
        rows=[json.loads(p.read_text()) for p in sorted(directory.glob('*.json'))];assert len(rows)==4
        development.append(dict(method=method,macro_full_f1=choice['macro_full_f1'],
            by_particle={p:metrics(rows,p,development_inputs) for p in ('alpha','gamma')},
            diagnostics={k:sum(r['score'][k] for r in rows) for k in ('missed_count','unmatched_predictions','split_flags','merge_flags')},
            postprocess_seconds=sum(r['fit_seconds'] for r in rows)))
    methods=['model_selection_1',selected];all_rows=[];groups=[];pooled=[]
    for method in methods:
        rows=[json.loads(p.read_text()) for p in sorted((artifacts/'test_mixed'/method).glob('*.json'))];assert len(rows)==8,(method,len(rows))
        all_rows+=rows
        for model in ('brems','mono662'):
            for rate in (1.,10.):
                subset=[r for r in rows if r['row']['gamma_model']==model and r['row']['rate_hz']==rate*1e6];assert len(subset)==2
                for particle in ('alpha','gamma'):
                    groups.append(dict(method=method,gamma_model=model,alpha_MHz=rate,particle=particle,windows=2,**metrics(subset,particle,artifacts/'inputs')))
        for particle in ('alpha','gamma'):pooled.append(dict(method=method,particle=particle,**metrics(rows,particle,artifacts/'inputs')))
    totals={m:{k:sum(r['score'][k] for r in all_rows if r['method']==m) for k in ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags')} for m in methods}
    seconds={m:sum(r['fit_seconds'] for r in all_rows if r['method']==m) for m in methods}
    result=dict(methods=methods,selected_method=selected,development=development,groups=groups,pooled=pooled,diagnostics=totals,
        seconds_recorded=seconds,selected_pipeline_seconds=seconds[methods[0]]+seconds[selected],
        runtime_note='Selected method adds postprocessing to the full V3 initializer; costs must be added, not compared as independent total runtimes',
        target_definition='Truth emission [5,15) us, realized Q>=500 pe; class-blind matching <=10 ns; correct class and |Qfit/Qtruth-1|<=20%',
        uncertainty='Eight independent windows, two per condition. Small gamma counts and tiny development improvement; no robust gain or maximum-rate claim from small changes.',
        class_priority='Alpha and gamma both matter; class-conditional recovery and prediction purity are reported separately',
        truth_use='Truth opened for scoring after inference; separate oracle and shape diagnostic explicitly use truth',
        test_tuning=False)
    save(artifacts/'summary.json',result);save(artifacts/'test-window-results.json',all_rows)
    pct=lambda x:'—' if x is None else f'{100*x:.1f}%'
    names={'model_selection_1':'上一版联合拟合',selected:'受约束可变形状'}
    dev=['| 开发候选 | 等权完整恢复F1 | α完整恢复 | γ完整恢复 | 额外预测 |','| --- | ---: | ---: | ---: | ---: |']
    for d in development:
        a=d['by_particle']['alpha'];g=d['by_particle']['gamma']
        dev.append(f"| {d['method']} | {d['macro_full_f1']:.5f} | {a['recovered_count']}/{a['target_count']} | {g['recovered_count']}/{g['target_count']} | {d['diagnostics']['unmatched_predictions']} |")
    main=['| γ来源 | α源/MHz | 类别 | 上一版联合拟合 | 受约束可变形状 |','| --- | ---: | --- | --- | --- |']
    detail=['| γ / α MHz | 方法 | 类别 | 正确分型/预测 | 电荷偏差/标准差 | 电荷P16/P50/P84 | 时间RMSE/ns |','| --- | --- | --- | --- | --- | --- | ---: |']
    fail=['| γ / α MHz | 方法 | 类别 | 完整恢复 | 漏检 | 错分型 | 仅电荷不合格 |','| --- | --- | --- | ---: | ---: | ---: | ---: |']
    for model in ('brems','mono662'):
        for rate in (1.,10.):
            for particle in ('alpha','gamma'):
                gs=[next(g for g in groups if (g['method'],g['gamma_model'],g['alpha_MHz'],g['particle'])==(m,model,rate,particle)) for m in methods]
                main.append(f'| {model} | {rate:g} | {particle} | '+' | '.join(f"{g['recovered_count']}/{g['target_count']}（{pct(g['full_recovery'])}）" for g in gs)+' |')
                for g in gs:
                    quantiles=g['charge_p16_p50_p84'];timing='—' if g['time_rmse_ns'] is None else f"{g['time_rmse_ns']:.2f}"
                    detail.append(f"| {model} / {rate:g} | {names[g['method']]} | {particle} | {g['pure_prediction_count']}/{g['prediction_count']}（{pct(g['type_purity'])}） | {pct(g['charge_bias'])}/{pct(g['charge_std'])} | "+(' / '.join(map(pct,quantiles)) if quantiles else '—')+f' | {timing} |')
                    qbad=g['correct_type_count']-g['recovered_count']
                    assert g['target_count']==g['recovered_count']+g['missed_count']+g['wrong_type_count']+qbad
                    fail.append(f"| {model} / {rate:g} | {names[g['method']]} | {particle} | {g['recovered_count']} | {g['missed_count']} | {g['wrong_type_count']} | {qbad} |")
    counts=['| 方法 | 目标 | 时间匹配 | 漏检 | 额外预测 | 误拆分标志 | 误合并标志 |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for m in methods:counts.append('| '+names[m]+' | '+' | '.join(str(v) for v in totals[m].values())+' |')
    save(artifacts/'tables.json',dict(development='\n'.join(dev),main='\n'.join(main),details='\n'.join(detail),failures='\n'.join(fail),counts='\n'.join(counts)))
    return result


def figures(artifacts):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    artifacts=Path(artifacts);s=json.loads((artifacts/'summary.json').read_text())
    ds=s['development'];colors=['#c76512' if d['method']==s['selected_method'] else '#236daa' if d['method']=='model_selection_1' else '#8d969e' for d in ds]
    fig,ax=plt.subplots(figsize=(10,6))
    ax.barh([d['method'] for d in ds],[d['macro_full_f1'] for d in ds],color=colors)
    ax.axvline(ds[0]['macro_full_f1'],color='#236daa',ls='--',lw=1)
    for i,d in enumerate(ds):ax.text(d['macro_full_f1']+.004,i,f"{d['macro_full_f1']:.5f}",va='center',fontsize=9)
    ax.invert_yaxis();ax.set_xlim(0,.72);ax.set_xlabel('Equal class / condition full-recovery F1')
    ax.set_title('Development only: 12 variants, four previously disclosed windows');fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(artifacts/f'development-comparison.{ext}',dpi=160,bbox_inches='tight')
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,7.5),sharey=True)
    conditions=[(g,r) for g in ('brems','mono662') for r in (1.,10.)]
    for i,metric in enumerate(('full_recovery','type_purity')):
        for j,particle in enumerate(('alpha','gamma')):
            ax=axes[i,j]
            for k,(method,color,label) in enumerate(zip(s['methods'],['#236daa','#c76512'],['Previous joint fit','Bounded shape fit'])):
                gs=[next(g for g in s['groups'] if (g['method'],g['particle'],g['gamma_model'],g['alpha_MHz'])==(method,particle,*condition)) for condition in conditions]
                x=np.arange(4)+(k-.5)*.32;ys=[100*g[metric] if g[metric] is not None else 0. for g in gs]
                ax.bar(x,ys,width=.30,color=color,label=label)
                for xp,y,g in zip(x,ys,gs):
                    n=g['target_count'] if i==0 else g['prediction_count'];v=g['recovered_count'] if i==0 else g['pure_prediction_count']
                    ax.text(xp,y+2+4*k,f'{v}/{n}',ha='center',fontsize=8)
            ax.set_ylim(0,116);ax.set_yticks([0,25,50,75,100]);ax.grid(axis='y',alpha=.2)
            ax.set_xticks(np.arange(4),['Brems\n1 MHz','Brems\n10 MHz','662 keV\n1 MHz','662 keV\n10 MHz'])
            ax.set_title(particle+(': full recovery' if i==0 else ': type purity'))
            if j==0:ax.set_ylabel('Percent')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,ncol=2,loc='upper center',bbox_to_anchor=(.5,.96))
    fig.suptitle('Frozen parameters; fresh native tests, two windows per condition',y=1.01);fig.tight_layout(rect=[0,0,1,.93])
    fig.supxlabel('Alpha source injection rate (gamma source fixed at 10 MHz)',y=-.015,fontsize=11)
    for ext in ('png','pdf'):fig.savefig(artifacts/f'heldout-shape-comparison.{ext}',dpi=160,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('folder','artifacts','previous','development-inputs'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();s=summarize(a.folder,a.artifacts,a.previous,a.development_inputs);figures(a.artifacts)
    print(json.dumps(dict(pooled=s['pooled'],diagnostics=s['diagnostics'],seconds=s['seconds_recorded']),ensure_ascii=False),flush=True)
