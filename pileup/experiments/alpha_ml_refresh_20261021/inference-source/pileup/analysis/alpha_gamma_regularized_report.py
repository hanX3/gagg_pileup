"""Freeze validation choices and summarize an untouched mixed holdout."""
from __future__ import annotations
from dataclasses import asdict
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import numpy as np
from .alpha_gamma_regularized_fit import RegularizedConfig
from .alpha_rate_campaign import digest,save

PROJECT=Path(__file__).resolve().parents[2]


def metrics(rows,particle,inputs=None):
    data=[r['score']['by_particle'][particle] for r in rows]
    keys=('target_count','matched_count','correct_type_count','recovered_count','missed_count','wrong_type_count',
          'central_injected','central_nonzero_charge')
    counts={k:sum(r[k] for r in data) for k in keys}
    predictions=sum(sum(p==particle and r['gate_ns'][0]<=t<r['gate_ns'][1]
                        for p,t in zip(r['fit']['particle_name'],r['fit']['times_ns'])) for r in rows)
    errors=np.array([x for d in data for x in d['relative_charge_errors']])
    timing=np.array([x for d in data for x in d['time_errors_ns']])
    n=counts['target_count'];correct=counts['recovered_count']
    pure=counts['correct_type_count'];precise=correct
    # Selection keeps the predeclared target-count F1 in the frozen record.
    # Final prediction purity uses the prediction gate, including matched guard
    # truths and excluding target matches whose fitted time is beyond the gate.
    if inputs is not None:
        pure=precise=0
        for r in rows:
            truth=dict(np.load(Path(inputs)/'truth'/f"{Path(r['row']['root']).stem}.npz"))
            for i,j in r['score']['all_matched_pairs']:
                if (r['gate_ns'][0]<=r['fit']['times_ns'][j]<r['gate_ns'][1] and
                    r['fit']['particle_name'][j]==particle and truth['particle_name'][i]==particle):
                    pure+=1
                    precise+=int(abs(r['fit']['charge_pe'][j]/truth['charge_pe'][i]-1)<=.2)
    recall=correct/n if n else None;precision=precise/predictions if predictions else None
    f1=(2*recall*precision/(recall+precision) if recall+precision else 0.) if recall is not None and precision is not None else None
    return dict(**counts,prediction_count=predictions,pure_prediction_count=pure,precise_prediction_count=precise,
        full_recovery=recall,type_purity=pure/predictions if predictions else None,
        full_precision=precision,
        full_f1=f1 if inputs is not None else (2*correct/(n+predictions) if n+predictions else None),
        charge_bias=float(errors.mean()) if len(errors) else None,
        charge_std=float(errors.std(ddof=1)) if len(errors)>1 else None,
        charge_p16_p50_p84=np.quantile(errors,[.16,.5,.84]).tolist() if len(errors) else None,
        time_rmse_ns=float(np.sqrt(np.mean(timing**2))) if len(timing) else None,
        unmatched_predictions=sum(r['score']['unmatched_predictions_by_class'][particle] for r in rows))


def freeze(folder,artifacts):
    folder,artifacts=Path(folder),Path(artifacts)
    frozen=folder/'frozen-fit-settings.json'
    if frozen.exists():raise FileExistsError(frozen)
    rule=json.loads((folder/'selection-rule.json').read_text());choices=[]
    for penalty in rule['penalties']:
        files=sorted((artifacts/f'validation_mixed/regularized_{penalty:g}').glob('*.json'))
        assert len(files)==4
        rows=[json.loads(f.read_text()) for f in files]
        assert len({r['fitter_sha256'] for r in rows})==1
        cells=[dict(gamma_model=r['row']['gamma_model'],rate_hz=r['row']['rate_hz'],particle=p,
                    **metrics([r],p)) for r in rows for p in ('alpha','gamma')]
        values=[c['full_f1'] for c in cells if c['full_f1'] is not None]
        choices.append(dict(penalty=penalty,macro_full_f1=float(np.mean(values)),cells=cells,
                            result_sha256={str(p):digest(p) for p in files}))
    selected=max(choices,key=lambda x:(x['macro_full_f1'],x['penalty']))
    source=PROJECT/'pileup/analysis/alpha_gamma_regularized_fit.py'
    cal=artifacts/'calibration.json'
    result=dict(created=datetime.now().astimezone().isoformat(),config=asdict(RegularizedConfig(
        selection_penalty=selected['penalty'])),selection_rule=rule,validation_choices=choices,
        fitter_sha256=digest(source),calibration_sha256=digest(cal),test_opened=False,
        environment=dict(python=platform.python_version(),platform=platform.platform(),
            **{k:importlib.metadata.version(k) for k in ('numpy','scipy','uproot','awkward','matplotlib')}))
    # Preserve the initially run launcher before correcting its frozen import
    # context. Physical C++ and all generated macros remain byte-for-byte.
    old=folder/'source/pileup/analysis/alpha_gamma_regularized_campaign.py'
    shutil.copy2(old,folder/'initial-campaign-source.py.txt')
    for package in ('pileup/analysis','analysis'):
        for p in (PROJECT/package).glob('*.py'):
            dst=folder/'source'/package/p.name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
    for package in ('pileup','pileup/analysis','analysis'):
        target=folder/'source'/package/'__init__.py'
        if not target.exists():target.write_text('')
    cfg=folder/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json'
    cfg.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(folder/'source/brems/g4/configs'/cfg.name,cfg)
    shutil.copy2(cal,folder/'calibration.json')
    shutil.copy2(PROJECT/'jupyter/run_records/alpha_gamma_fit_20261015_093500/calibration-settings.json',
                 folder/'baseline-settings.json')
    result['source_sha256']={str(p.relative_to(folder)):digest(p) for p in (folder/'source').rglob('*') if p.is_file()}
    save(frozen,result)
    (folder/'.gitignore').write_text('bin/\n__pycache__/\n')
    (folder/'run_analysis.py').write_text('''from pathlib import Path
import runpy,sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
sys.path.insert(1,str(folder.parents[2]))
runpy.run_module('pileup.analysis.alpha_gamma_regularized_study',run_name='__main__')
''')
    (folder/'run_campaign.py').write_text('''from pathlib import Path
import importlib.util,sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
sys.path.insert(1,str(folder.parents[2]))
import pileup.analysis.alpha_gamma_regularized_campaign as module
module.PROJECT=folder.parents[2]
module.main()
''')
    print('FROZEN',selected['penalty'],[(c['penalty'],c['macro_full_f1']) for c in choices],flush=True)
    return result


def summarize(folder,artifacts):
    folder,artifacts=Path(folder),Path(artifacts)
    frozen=json.loads((folder/'frozen-fit-settings.json').read_text());penalty=frozen['config']['selection_penalty']
    methods=['baseline',f'regularized_{penalty:g}'];groups=[];all_rows=[]
    for method in methods:
        files=sorted((artifacts/'test_mixed'/method).glob('*.json'));assert len(files)==8
        rows=[json.loads(f.read_text()) for f in files];all_rows+=rows
        for model in ('brems','mono662'):
            for rate in (1000000,10000000):
                subset=[r for r in rows if r['row']['gamma_model']==model and r['row']['rate_hz']==rate]
                assert len(subset)==2
                for particle in ('alpha','gamma'):
                    groups.append(dict(method=method,gamma_model=model,alpha_injection_MHz=rate/1e6,
                        gamma_injection_MHz=10.,particle=particle,windows=2,**metrics(subset,particle,artifacts/'inputs')))
    result=dict(groups=groups,methods=methods,selected_penalty=penalty,heldout_windows=8,
        definition='Q truth >=500 pe, class-blind one-to-one time match <=10 ns, correct type, charge error <=20%',
        error_conditioning='Time-matched targets, including misclassified events; charge is not deposited energy',
        purity_gate='Predicted time within central gate; correct-type matches to all above-threshold truths, including guards',
        recovery_gate='Truth emission within central gate, even if fitted time is just beyond the gate',
        fit_seconds={m:sum(r['fit_seconds'] for r in all_rows if r['method']==m) for m in methods},
        diagnostics={m:{k:sum(r['score'][k] for r in all_rows if r['method']==m)
                     for k in ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags')}
                     for m in methods},
        sparse_solver=dict(success_flags=sum(r['fit']['selection']['success'] for r in all_rows if r['method']!='baseline'),
            max_projected_gradient=max(r['fit']['selection']['projected_gradient_max'] for r in all_rows if r['method']!='baseline'),
            interpretation='Relative objective stopping is not strict KKT convergence; selected support is approximate'),
        source_sha256=frozen['fitter_sha256'],calibration_sha256=frozen['calibration_sha256'],
        limitations=['Two windows per condition; no precise efficiency or maximum-rate claim',
            'Central isolated calibration, only 13 nonzero gamma calibration examples; no occupancy/position calibration',
            'AR(1) covariance and charge-dependent mean template are approximations',
            'Only sparse support is selected; discrete/time refinements are local, not guaranteed global',
            'L1 solver success flags do not imply strict KKT tolerance; saved projected-gradient diagnostics must be inspected',
            'No calibrated classification posterior or rejection policy; no deposited-energy calibration'])
    save(artifacts/'summary.json',result);save(artifacts/'test-window-results.json',all_rows)
    print(json.dumps(groups,ensure_ascii=False,indent=2))
    return result


def make_figures(folder, artifacts):
    """Plot saved predictions, without changing any fit or selection setting."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator, ScalarFormatter
    from .alpha_gamma_regularized_fit import CalibratedBank
    folder, artifacts = Path(folder), Path(artifacts)
    summary = json.loads((artifacts/'summary.json').read_text())
    cal = json.loads((artifacts/'calibration.json').read_text())
    old = json.loads((folder/'baseline-settings.json').read_text())
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    figures = []
    def save_figure(fig, name):
        for extension in ('png', 'pdf'):
            path = artifacts/f'{name}.{extension}'
            fig.savefig(path, dpi=170, bbox_inches='tight')
            figures.append(str(path))
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, particle, color in zip(axes, ('alpha','gamma'), ('#1565a0','#c15f12')):
        rows = [r for r in cal['isolated_rows'] if r['particle']==particle and r['used']]
        ax.scatter([r['charge_pe'] for r in rows], [r['slow_fraction'] for r in rows],
                   color=color, label=f'Independent singles (n={len(rows)})')
        par = cal['classes'][particle]
        q = np.geomspace(*par['charge_range_pe'], 100)
        ax.plot(q, par['slow_intercept']+par['slow_log_slope']*np.log(q/1000), color=color, label='Calibrated mean')
        ax.axhline(old['templates'][particle]['slow_fraction'], linestyle='--', color='0.45', label='Old example template')
        ax.set(xscale='log', xlabel='Realized charge Q (p.e.)', ylabel='Slow-component fraction', title=particle)
        ax.set_xticks([2000,5000,10000] if particle=='alpha' else [300,600,1000,2000])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_locator(NullLocator())
        ax.legend(fontsize=8)
    fig.suptitle('Isolated central-source calibration; not an occupancy calibration', fontsize=11)
    fig.tight_layout(); save_figure(fig, 'calibration-shapes')

    conditions = [(g,r) for g in ('brems','mono662') for r in (1.,10.)]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True)
    for col, particle in enumerate(('alpha','gamma')):
        for row, metric in enumerate(('full_recovery','type_purity')):
            ax = axes[row,col]
            for i, (method, color, label) in enumerate(zip(summary['methods'], ('#1565a0','#c15f12'), ('Baseline','Calibrated sparse fit'))):
                vals = [next(g for g in summary['groups'] if g['method']==method and g['particle']==particle and
                             (g['gamma_model'],g['alpha_injection_MHz'])==c) for c in conditions]
                x = np.arange(4)+(i-.5)*.36
                heights = [100*v[metric] if v[metric] is not None else 0 for v in vals]
                ax.bar(x, heights, width=.34, color=color, label=label)
                for xp, height, v in zip(x, heights, vals):
                    numerator = v['recovered_count'] if metric=='full_recovery' else v['pure_prediction_count']
                    denominator = v['target_count'] if metric=='full_recovery' else v['prediction_count']
                    ax.text(xp,height+2+(5*i if metric=='type_purity' else 0),f'{numerator}/{denominator}',ha='center',fontsize=8,rotation=0)
            ax.set_ylim(0,118);ax.set_yticks([0,25,50,75,100]);ax.grid(axis='y',alpha=.2)
            ax.set_title(f'{particle}: '+('full recovery' if metric=='full_recovery' else 'type purity'))
            ax.set_xticks(np.arange(4),['Brems\n1 MHz','Brems\n10 MHz','662 keV\n1 MHz','662 keV\n10 MHz'])
            if col==0:ax.set_ylabel('Percent')
            if row==1:ax.set_xlabel('Gamma spectrum / alpha source injection rate')
    axes[0,0].legend(fontsize=8,loc='lower left')
    fig.suptitle('New held-out native waveforms: two windows per condition; gamma source = 10 MHz',fontsize=11)
    fig.tight_layout();save_figure(fig,'heldout-comparison')

    rows = json.loads((artifacts/'test-window-results.json').read_text())
    fig, axes = plt.subplots(4, 2, figsize=(12, 10), sharex=True)
    examples = []
    for i,(model, rate) in enumerate(conditions):
        pairs = {method:sorted([r for r in rows if r['method']==method and r['row']['gamma_model']==model and
                                r['row']['rate_hz']==rate*1e6],key=lambda r:r['row']['label'])[0]
                 for method in summary['methods']}
        root = pairs['baseline']['row']['root'];stem=Path(root).stem
        wave = dict(np.load(artifacts/'inputs/waveforms'/f'{stem}.npz'));t=wave['time_ns']/1000
        mask=(t>=5.5)&(t<=7.5)
        axes[i,0].plot(t[mask],wave['analog_mV'][mask],color='0.65',lw=.8,label='Native noisy voltage')
        for method,color,label in zip(summary['methods'],('#1565a0','#c15f12'),('Baseline','Calibrated sparse fit')):
            fitted=np.load(artifacts/'test_mixed'/method/f'{stem}.npz')['model_mV']
            axes[i,0].plot(t[mask],fitted[mask],color=color,lw=1,label=label)
            axes[i,1].plot(t[mask],(wave['analog_mV']-fitted)[mask],color=color,lw=.7,label=label)
        axes[i,0].set_ylabel(f'{model}, alpha {rate:g} MHz\nVoltage (mV)')
        axes[i,1].set_ylabel('Residual (mV)');axes[i,1].axhline(0,color='0.6',lw=.5)
        examples.append(dict(gamma_model=model,alpha_MHz=rate,root=root,range_us=[5.5,7.5]))
    for ax in axes[-1]:ax.set_xlabel('Time (us)')
    axes[0,0].legend(fontsize=8);axes[0,1].legend(fontsize=8)
    fig.suptitle('Fixed examples: first test window of each condition, same 5.5–7.5 us interval',fontsize=11)
    fig.tight_layout();save_figure(fig,'heldout-waveforms')
    save(artifacts/'figure-record.json',dict(files=figures,examples=examples,
         selection='First planned test label of each condition; fixed 5.5–7.5 us, no result-based selection'))
    return figures


def comparison_tables(artifacts):
    """Create Chinese tables from the unchanged scored JSON, including failures."""
    artifacts=Path(artifacts)
    summary=json.loads((artifacts/'summary.json').read_text())
    windows=json.loads((artifacts/'test-window-results.json').read_text())
    pct=lambda x:'—' if x is None else f'{100*x:.1f}%'
    def lookup(model,rate,method,particle):
        return next(g for g in summary['groups'] if (g['gamma_model'],g['alpha_injection_MHz'],g['method'],g['particle'])==
                    (model,rate,method,particle))
    main=['| γ 来源 | α 注入率 / MHz | α 完整恢复：旧 → 新 | γ 完整恢复：旧 → 新 |',
          '| --- | ---: | --- | --- |']
    details=['| γ / α MHz | 方法 | 类别 | 时间匹配 / 目标 | 分型正确 / 预测（纯度） | 漏检 | 错类型 | 电荷偏差 / 条件标准差 | 电荷误差 P16 / P50 / P84 | 时间 RMSE / ns |',
             '| --- | --- | --- | ---: | --- | ---: | ---: | --- | --- | ---: |']
    paired=[]
    for model in ('brems','mono662'):
        for rate in (1.,10.):
            name='轫致辐射' if model=='brems' else '661.657 keV'
            counts=[]
            for particle in ('alpha','gamma'):
                vals=[lookup(model,rate,m,particle) for m in summary['methods']]
                counts.append(' → '.join(f"{v['recovered_count']}/{v['target_count']}（{pct(v['full_recovery'])}）" for v in vals))
            main.append(f'| {name} | {rate:g} | '+ ' | '.join(counts)+' |')
            for method in summary['methods']:
                for particle in ('alpha','gamma'):
                    v=lookup(model,rate,method,particle);qs=v['charge_p16_p50_p84']
                    details.append(f"| {model} / {rate:g} | {'旧' if method=='baseline' else '新'} | {particle} | {v['matched_count']}/{v['target_count']} | "
                        f"{v['pure_prediction_count']}/{v['prediction_count']}（{pct(v['type_purity'])}） | {v['missed_count']} | {v['wrong_type_count']} | "
                        f"{pct(v['charge_bias'])} / {pct(v['charge_std'])} | "+(' / '.join(map(pct,qs)) if qs else '—')+' | '+
                        ('—' if v['time_rmse_ns'] is None else f"{v['time_rmse_ns']:.2f}")+' |')
            subset=[r for r in windows if r['method']=='baseline' and r['row']['gamma_model']==model and r['row']['rate_hz']==rate*1e6]
            for particle in ('alpha','gamma'):
                counts=dict(both=0,baseline_only=0,regularized_only=0,neither=0)
                for old in subset:
                    new=next(r for r in windows if r['method']!='baseline' and r['row']['root']==old['row']['root'])
                    truth=dict(np.load(artifacts/'inputs/truth'/f"{Path(old['row']['root']).stem}.npz"))
                    def recovered(r):
                        return {i for i,j in r['score']['matched_pairs'] if r['fit']['particle_name'][j]==particle and
                                truth['particle_name'][i]==particle and abs(r['fit']['charge_pe'][j]/truth['charge_pe'][i]-1)<=.2}
                    a,b=recovered(old),recovered(new)
                    targets={i for i in old['score']['target_indices'] if truth['particle_name'][i]==particle}
                    for key,s in (('both',a&b),('baseline_only',a-b),('regularized_only',b-a),('neither',targets-(a|b))):
                        counts[key]+=len(s)
                paired.append(dict(gamma_model=model,alpha_MHz=rate,particle=particle,**counts))
    diagnostic=['| 方法 | 目标 | 时间匹配 | 漏检 | 额外预测 | 误拆分标志 | 误合并标志 |',
                '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for method in summary['methods']:
        d=summary['diagnostics'][method]
        diagnostic.append('| '+('旧' if method=='baseline' else '新')+' | '+' | '.join(str(d[k]) for k in
            ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags'))+' |')
    perwindow=['| 测试窗 | 方法 | α 完整恢复 | γ 完整恢复 | 额外预测 | 单窗拟合 / s |',
               '| --- | --- | ---: | ---: | ---: | ---: |']
    for r in sorted(windows,key=lambda r:(r['row']['label'],r['method'])):
        scores=r['score']['by_particle']
        perwindow.append('| '+r['row']['label']+' | '+('旧' if r['method']=='baseline' else '新')+' | '+
            ' | '.join(f"{scores[p]['recovered_count']}/{scores[p]['target_count']}" for p in ('alpha','gamma'))+
            f" | {r['score']['unmatched_predictions']} | {r['fit_seconds']:.2f} |")
    result=dict(main='\n'.join(main),details='\n'.join(details),diagnostics='\n'.join(diagnostic),
                per_window='\n'.join(perwindow),paired=paired)
    save(artifacts/'comparison-tables.json',result)
    return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('freeze','summarize','report'));p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--artifacts',type=Path,required=True);a=p.parse_args()
    (freeze if a.stage=='freeze' else summarize)(a.folder,a.artifacts)
    if a.stage=='report':
        make_figures(a.folder,a.artifacts)
        comparison_tables(a.artifacts)
