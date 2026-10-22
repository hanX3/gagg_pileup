"""Report frozen three-way alpha/gamma tests and their numerical limits."""
from pathlib import Path
import argparse
import json
import numpy as np
from .alpha_gamma_regularized_report import metrics
from .alpha_gamma_regularized_fit import CalibratedBank
from .alpha_rate_campaign import save,digest


def summarize(folder,artifacts):
    folder,artifacts=Path(folder),Path(artifacts)
    frozen=json.loads((folder/'frozen-fit-settings.json').read_text())
    methods=['baseline','regularized_3',f"model_selection_{frozen['config']['addition_cost']:g}"]
    windows=[];groups=[]
    for method in methods:
        rows=[json.loads(p.read_text()) for p in sorted((artifacts/'test_mixed'/method).glob('*.json'))]
        assert len(rows)==8,(method,len(rows))
        windows.extend(rows)
        for gamma in ('brems','mono662'):
            for rate in (1.,10.):
                selected=[r for r in rows if r['row']['gamma_model']==gamma and r['row']['rate_hz']==rate*1e6]
                assert len(selected)==2
                for particle in ('alpha','gamma'):
                    groups.append(dict(method=method,gamma_model=gamma,alpha_MHz=rate,particle=particle,
                        windows=2,**metrics(selected,particle,artifacts/'inputs')))
    diagnostics={method:{k:sum(r['score'][k] for r in windows if r['method']==method)
        for k in ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags')}
        for method in methods}
    result=dict(groups=groups,methods=methods,diagnostics=diagnostics,
        fit_seconds={m:sum(r['fit_seconds'] for r in windows if r['method']==m) for m in methods},
        target_gate='Truth emission in [5,15) us; Q>=500 pe; type-blind maximum time matching <=10 ns; correct class and Q error<=20%',
        purity_gate='Fitted time in [5,15) us; correct class matched to any threshold truth including guards',
        uncertainty='Two windows/condition only; no precise maximum-rate or gamma-efficiency claim',
        calibration_sha256=frozen['calibration_sha256'],new_fitter_sha256=frozen['fitter_sha256'],
        equal_priority='Alpha and gamma are evaluated separately; neither class can be hidden in a total accuracy',
        truth_use='No truth enters any blind fitter; separate oracle is development diagnosis only')
    save(artifacts/'summary.json',result);save(artifacts/'test-window-results.json',windows)
    names={'baseline':'原始拟合','regularized_3':'上一版稀疏',methods[-1]:'本轮联合'}
    main=['| γ来源 | α源/MHz | 类别 | 原始拟合 | 上一版稀疏 | 本轮联合 |','| --- | ---: | --- | --- | --- | --- |']
    details=['| γ / α MHz | 方法 | 类别 | 时间匹配/目标 | 纯度分子/预测 | 电荷偏差/标准差 | 电荷误差P16/P50/P84 | 时间RMSE/ns |',
        '| --- | --- | --- | --- | --- | --- | --- | --- |']
    pct=lambda x:'—' if x is None else f'{100*x:.1f}%'
    for gamma in ('brems','mono662'):
        for rate in (1.,10.):
            for particle in ('alpha','gamma'):
                gs=[next(g for g in groups if (g['gamma_model'],g['alpha_MHz'],g['particle'],g['method'])==
                      (gamma,rate,particle,m)) for m in methods]
                main.append(f'| {gamma} | {rate:g} | {particle} | '+' | '.join(f"{g['recovered_count']}/{g['target_count']}（{pct(g['full_recovery'])}）" for g in gs)+' |')
                for g in gs:
                    quantiles=g['charge_p16_p50_p84'];timing='—' if g['time_rmse_ns'] is None else f"{g['time_rmse_ns']:.2f}"
                    details.append(f"| {gamma} / {rate:g} | {names[g['method']]} | {particle} | {g['matched_count']}/{g['target_count']} | "
                        f"{g['pure_prediction_count']}/{g['prediction_count']}（{pct(g['type_purity'])}） | {pct(g['charge_bias'])}/{pct(g['charge_std'])} | "+
                        (' / '.join(map(pct,quantiles)) if quantiles else '—')+f' | {timing} |')
    diagnostic=['| 方法 | 目标 | 时间匹配 | 漏检 | 额外预测 | 误拆分标志 | 误合并标志 |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for m in methods:diagnostic.append('| '+names[m]+' | '+' | '.join(str(v) for v in diagnostics[m].values())+' |')
    failure=['| γ / α MHz | 方法 | 类别 | 完整恢复 | 漏检 | 错类型 | 类型正确但电荷不合格 |',
        '| --- | --- | --- | ---: | ---: | ---: | ---: |']
    for g in groups:
        charge_fail=g['correct_type_count']-g['recovered_count']
        assert g['target_count']==g['recovered_count']+g['missed_count']+g['wrong_type_count']+charge_fail
        failure.append(f"| {g['gamma_model']} / {g['alpha_MHz']:g} | {names[g['method']]} | {g['particle']} | "
            f"{g['recovered_count']} | {g['missed_count']} | {g['wrong_type_count']} | {charge_fail} |")
    table=dict(main='\n'.join(main),details='\n'.join(details),diagnostics='\n'.join(diagnostic),failures='\n'.join(failure))
    save(artifacts/'tables.json',table)
    return result


def figures(folder,artifacts,previous):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder,artifacts,previous=Path(folder),Path(artifacts),Path(previous)
    summary=json.loads((artifacts/'summary.json').read_text());methods=summary['methods'];conditions=[(g,r) for g in ('brems','mono662') for r in (1.,10.)]
    colors=['#7b858f','#236daa','#c76512'];labels=['Original fit','Previous sparse fit','Joint fit + residual proposals']
    fig,axes=plt.subplots(2,2,figsize=(13,7.5),sharex=True,sharey=True)
    for row,metric in enumerate(('full_recovery','type_purity')):
        for col,particle in enumerate(('alpha','gamma')):
            ax=axes[row,col]
            for j,(method,color,label) in enumerate(zip(methods,colors,labels)):
                gs=[next(g for g in summary['groups'] if (g['method'],g['particle'],g['gamma_model'],g['alpha_MHz'])==
                         (method,particle,*condition)) for condition in conditions]
                x=np.arange(4)+(j-1)*.24
                values=[100*g[metric] if g[metric] is not None else 0. for g in gs]
                ax.bar(x,values,width=.23,color=color,label=label)
                for xp,value,g in zip(x,values,gs):
                    n=g['target_count'] if metric=='full_recovery' else g['prediction_count']
                    k=g['recovered_count'] if metric=='full_recovery' else g['pure_prediction_count']
                    ax.text(xp,value+2+j*5,f'{k}/{n}',ha='center',fontsize=7)
            ax.set_ylim(0,118);ax.set_yticks([0,25,50,75,100]);ax.grid(axis='y',alpha=.2)
            ax.set_title(particle+(': full recovery' if row==0 else ': type purity'))
            ax.set_xticks(np.arange(4),['Brems\n1 MHz','Brems\n10 MHz','662 keV\n1 MHz','662 keV\n10 MHz'])
            if col==0:ax.set_ylabel('Percent')
    handles,_=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,ncol=3,loc='upper center',bbox_to_anchor=(.5,.96),fontsize=9)
    fig.suptitle('New independent native test: alpha and gamma both matter; two windows per condition',y=1.01)
    fig.tight_layout(rect=[0,0,1,.94])
    for ext in ('png','pdf'):fig.savefig(artifacts/f'heldout-three-methods.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)
    cal=json.loads((folder/'calibration.json').read_text());bank=CalibratedBank(cal);t=np.arange(.25,12000.,2.)
    old=json.loads((previous/'controlled.json').read_text())['rows']
    new=json.loads((artifacts/'controlled-prototype.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(11,6.5))
    for i,(gap,pair) in enumerate(((10.,['alpha','gamma']),(200.,['gamma','alpha']))):
        a=next(r for r in old if r['gap_ns']==gap and r['pair']==pair and r['fit_stride']==2)
        b=next(r for r in new if r['gap_ns']==gap and r['pair']==pair)
        types=[0 if p=='alpha' else 1 for p in pair]
        y=10.+bank.matrix(t,[5000.,5000.+gap],types,[3000.,1000.])@np.array([3000.,1000.])
        y+=np.random.default_rng(a['seed']).normal(0,.1,len(t))
        mask=(t>=4975)&(t<=5400)
        axes[i,0].plot(t[mask]-5000,y[mask],color='0.7',lw=1,label='Control waveform')
        for r,color,label in ((a,colors[1],'Previous sparse fit'),(b,colors[2],'Joint fit')):
            c=r['fit']['all_components'];k=[0 if s=='alpha' else 1 for s in c['particle_name']]
            model=10.+bank.matrix(t,c['times_ns'],k,c['shape_charge_pe'])@np.array(c['charge_pe'])
            axes[i,0].plot(t[mask]-5000,model[mask],color=color,label=label)
            axes[i,1].plot(t[mask]-5000,(y-model)[mask],color=color,label=label)
        axes[i,0].set_title(f'{pair[0]} then {pair[1]}, gap={gap:g} ns')
        axes[i,0].set_ylabel('Voltage (mV)');axes[i,1].set_ylabel('Residual (mV)')
        axes[i,1].axhline(0,color='0.7',lw=.5)
    axes[0,0].legend(fontsize=8)
    for ax in axes[-1]:ax.set_xlabel('Time relative to first pulse (ns)')
    fig.suptitle('Numerical controls: same template + 0.1 mV white noise; not physical resolution')
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(artifacts/f'controlled-failure-repair.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('folder','artifacts','previous'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();summarize(a.folder,a.artifacts);figures(a.folder,a.artifacts,a.previous)
