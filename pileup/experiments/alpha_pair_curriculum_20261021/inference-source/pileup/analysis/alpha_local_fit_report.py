"""Postprocessing and notebook delivery for the frozen local-alpha experiment."""
from __future__ import annotations

import ast
import argparse
import base64
from datetime import datetime
import json
from pathlib import Path
import re
import shutil
import sys

project=next(p for p in (Path.cwd(),*Path.cwd().parents) if (p/'pileup/analysis/alpha_local_fit.py').exists())
sys.path.insert(0,str(project))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
import numpy as np

from pileup.analysis.alpha_local_fit_study import PROJECT,OUTPUT,EXPERIMENT,RECORD,settings
from pileup.analysis.alpha_full_signal_study import make_template
from pileup.analysis.alpha_rate_campaign import digest,save
from pileup.analysis.alpha_rate_data import open_root
from pileup.analysis.alpha_pileup_fit import score


def load():
    return json.loads((OUTPUT/'test-summary.json').read_text())


def conditional_metrics(summary):
    result=[]
    for rate in sorted({r['rate_hz'] for r in summary['windows']}):
        rows=[r for r in summary['windows'] if r['rate_hz']==rate]
        item=dict(rate_hz=rate,gained=0,lost=0,both_correct=0,both_failed=0,
            gained_previously_unmatched=0,gained_charge_corrected=0,
            lost_unmatched=0,lost_charge_failed=0)
        for method in ('baseline','local'):
            q=np.array([e for r in rows for e in r['methods'][method]['score']['relative_charge_errors']])
            tt=np.array([e for r in rows for e in r['methods'][method]['score']['time_errors_ns']])
            item[method]=dict(matched_charge_mean_percent=float(q.mean()*100) if len(q) else None,
                matched_charge_std_percent=float(q.std(ddof=1)*100) if len(q)>1 else None,
                matched_time_mean_ns=float(tt.mean()) if len(tt) else None,
                matched_time_std_ns=float(tt.std(ddof=1)) if len(tt)>1 else None)
        for row in rows:
            good={};matched={}
            for method in ('baseline','local'):
                s=row['methods'][method]['score']
                good[method]={int(pair[0]) for pair,e in zip(s['matched_pairs'],s['relative_charge_errors']) if abs(e)<=.2}
                matched[method]={int(pair[0]) for pair in s['matched_pairs']}
            old,new=good['baseline'],good['local'];n=row['methods']['baseline']['score']['target_count']
            item['gained']+=len(new-old);item['lost']+=len(old-new)
            item['both_correct']+=len(old&new);item['both_failed']+=n-len(old|new)
            item['gained_previously_unmatched']+=len((new-old)-matched['baseline'])
            item['gained_charge_corrected']+=len((new-old)&matched['baseline'])
            item['lost_unmatched']+=len((old-new)-matched['local'])
            item['lost_charge_failed']+=len((old-new)&matched['local'])
        result.append(item)
    return result


def shifted_controls(summary):
    """Post-fit time mismatches; diagnostic reference, not corrected efficiency."""
    path=OUTPUT/'time-shift-controls.json'
    if path.exists():return json.loads(path.read_text())
    shifts=(-1739.,-887.,887.,1739.)
    rows=[]
    for row in summary['windows']:
        stem=Path(row['root']).stem
        with np.load(OUTPUT/'inputs/truth'/f'{stem}.npz') as z:truth={k:z[k] for k in ('first_active_ns','emission_ns','charge_pe')}
        cases={}
        for method in ('baseline','local'):
            fit=row['methods'][method]['fit'];cases[method]=[]
            for shift in shifts:
                shifted={**fit,'times_ns':(np.asarray(fit['times_ns'])+shift).tolist()}
                scored=score(shifted,truth,[5000.,15000.])
                cases[method].append(dict(shift_ns=shift,target_count=scored['target_count'],
                    recovered_count=scored['recovered_count'],reconstructed_count=scored['reconstructed_count']))
        rows.append(dict(root=row['root'],rate_hz=row['rate_hz'],methods=cases))
    rates=[]
    for rate in sorted({r['rate_hz'] for r in rows}):
        group=[r for r in rows if r['rate_hz']==rate];item=dict(rate_hz=rate)
        n=sum(r['methods']['baseline'][0]['target_count'] for r in group)
        for method in ('baseline','local'):
            values=[sum(r['methods'][method][i]['recovered_count'] for r in group) for i in range(4)]
            item[method]=dict(target_count=n,recovered_counts=values,
                mean_fraction=float(np.mean(values)/n) if n else None,
                min_fraction=min(values)/n if n else None,max_fraction=max(values)/n if n else None)
        rates.append(item)
    result=dict(shifts_ns=list(shifts),rows=rows,rates=rates,
        limitation='Shift all reconstructed onsets, keep charges; five-us guards cover these shifts. Stationary-source diagnostic only; not an unbiased correction, confidence interval, or strict physical limit. Never used to choose frozen inference parameters.')
    save(path,result);return result


def plots(summary=None):
    plt.rcParams.update({'font.size':10,'axes.spines.right':False,'axes.spines.top':False})
    controls=json.loads((OUTPUT/'controls/controlled.json').read_text())['rows']
    fig,axes=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
    for ax,method,color in zip(axes,['baseline','local'],['#4b6a87','#b3552a']):
        for ratio in (.125,.25,.5,1.,2.,4.):
            values=[];gaps=[]
            for gap in (4.,8.,12.,20.,40.,80.,200.):
                cases=[r for r in controls if r['kind']=='double' and r['gap_ns']==gap and r['ratio']==ratio]
                values.append(100*np.mean([r['scores'][method]['recovered_count']==2 and r['scores'][method]['unmatched_predictions']==0 for r in cases]));gaps.append(gap)
            ax.plot(gaps,values,'o-',lw=1,label=f'Q2/Q1={ratio:g}')
        ax.set_xscale('log');ax.set_xticks(gaps,labels=[f'{g:g}' for g in gaps]);ax.set_ylim(-5,105)
        ax.set(xlabel='True separation (ns)',ylabel='Both recovered, no extra (%)',title=('(a) ' if method=='baseline' else '(b) ')+method)
        ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.2)
    fig.suptitle('Figure 17-1: Matched-template double controls; 3 electronic-noise trials per point')
    fig.savefig(OUTPUT/'local_controls.png',dpi=160);plt.close(fig)
    if summary is None:return

    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    rates=summary['rates']['baseline'];x=np.arange(len(rates));labels=[f"{r['rate_hz']/1e6:g}" for r in rates]
    for method,color in [('baseline','#4b6a87'),('local','#b3552a')]:
        rr=summary['rates'][method]
        axes[0,0].plot(x,[100*r['recovery_fraction'] if r['recovery_fraction'] is not None else np.nan for r in rr],'o-',color=color,label=method)
        axes[0,1].plot(x,[r['unmatched_predictions'] for r in rr],'o-',color=color,label=method)
        for i,r in enumerate(rr):
            a=r['charge_relative_p16_p50_p84']
            if a:
                lo,mid,hi=np.array(a)*100
                axes[1,0].errorbar(i+(-.08 if method=='baseline' else .08),mid,[[mid-lo],[hi-mid]],fmt='o',capsize=3,color=color)
        axes[1,1].plot(x,[r['missed_count'] for r in rr],'o-',color=color,label=method)
    axes[0,0].set(ylabel='Complete event recovery (%)',ylim=(-3,106))
    axes[0,1].set(ylabel='Unmatched predictions (count)')
    axes[1,0].set(ylabel='Matched charge error: median, P16-P84 (%)');axes[1,0].axhline(0,color='.6',ls=':')
    axes[1,1].set(ylabel='Targets without a time match (count)')
    for ax,title in zip(axes.flat,['(a) Complete recovery','(b) Unmatched outputs','(c) Conditional charge errors','(d) Missed targets']):ax.set_title(title)
    for ax in axes.flat:ax.set_xticks(x,labels);ax.set_xlabel('Source alpha injection rate (MHz)');ax.grid(alpha=.2)
    axes[0,0].legend();fig.suptitle('Figure 17-2: 23 new native windows; time <=10 ns and charge error <=20%')
    fig.savefig(OUTPUT/'local_native_rates.png',dpi=160);plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for index,(axrow,rate) in enumerate(zip(axes,(1e7,4e7))):
        row=next(r for r in summary['windows'] if r['rate_hz']==rate);stem=Path(row['root']).stem
        with np.load(OUTPUT/'inputs/waveforms'/f'{stem}.npz') as z:t,y=z['time_ns'],z['analog_mV']
        axrow[0].plot(t/1000,y,color='.6',lw=1,label='Observed voltage')
        for method,color in [('baseline','#4b6a87'),('local','#b3552a')]:
            with np.load(OUTPUT/'test/selected'/f'{stem}_{method}.npz') as z:axrow[0].plot(t/1000,z['model_mV'],color=color,lw=.9,label=method)
            fit=row['methods'][method]['fit'];xx=np.array(fit['times_ns']);q=np.array(fit['charge_pe']);gate=(xx>=5000)&(xx<7000)
            axrow[1].scatter(xx[gate]/1000,q[gate],marker='x' if method=='baseline' else '+',color=color,label=method)
        with np.load(OUTPUT/'inputs/truth'/f'{stem}.npz') as z:
            tt,qq=z['first_active_ns'],z['charge_pe'];mask=(tt>=5000)&(tt<7000)&(qq>=500)
        axrow[1].scatter(tt[mask]/1000,qq[mask],facecolors='none',edgecolors='.2',label='Truth for display only')
        axrow[0].set(xlim=(5,7),ylabel='Voltage (mV)',title=f'({chr(97+index*2)}) {rate/1e6:g} MHz; first predeclared test window')
        axrow[1].set(xlim=(5,7),ylabel='Realized avalanche charge (p.e.)',title=f'({chr(98+index*2)}) Individual event charges')
        for ax in axrow:ax.grid(alpha=.2);ax.set_xlabel('Record time (us)')
    axes[0,0].legend(fontsize=8);axes[0,1].legend(fontsize=8)
    fig.suptitle('Figure 17-3: Fixed 5-7 us display interval; total waveform and individual-event charge')
    fig.savefig(OUTPUT/'local_native_examples.png',dpi=160);plt.close(fig)


def validate(summary):
    frozen=json.loads((OUTPUT/'freeze.json').read_text());manifest=json.loads((EXPERIMENT/'manifest.json').read_text())
    for name,expected in frozen['frozen_files'].items():assert digest(PROJECT/name)==expected,name
    checks={name:json.loads((RECORD/name).read_text()) for name in
            ('numerical-tests.json','blind-replay-validation.json','read-cache-validation.json')}
    assert all(value['passed'] for value in checks.values())
    assert checks['numerical-tests.json']['source_sha256']==digest(PROJECT/'pileup/analysis/alpha_local_fit.py')
    assert checks['numerical-tests.json']['test_source_sha256']==digest(PROJECT/'pileup/analysis/test_alpha_local_fit.py')
    assert checks['blind-replay-validation.json']['scoring_truth_arrays_opened'] is False
    recovery=json.loads((RECORD/'redundant-time-validation.json').read_text())
    assert recovery['passed'] and recovery['reference_passed']
    assert len(recovery['all_arrays_exactly_equal'])==16
    assert recovery['recovered_checks']['one_to_one_complete_coverage']
    assert digest(PROJECT/'pileup/analysis/root_primary_time_recovery.py')==recovery['adapter_source_sha256']
    recovered_stem=Path(recovery['recovered_root']).stem
    with np.load(RECORD/'before-redundant-time-recovery/waveforms'/f'{recovered_stem}.npz') as old,np.load(OUTPUT/'inputs/waveforms'/f'{recovered_stem}.npz') as now:
        assert set(old.files)==set(now.files)
        for key in old.files:assert np.array_equal(old[key],now[key],equal_nan=old[key].dtype.kind in 'fc'),key
    capacity=json.loads((RECORD/'capacity-equivalence.json').read_text())
    assert capacity['passed'] and len(capacity['all_arrays_exactly_equal'])==16
    capacity_path=PROJECT/manifest['capacity_extension']['experiment']
    capacity_manifest=json.loads((capacity_path/'manifest.json').read_text())
    assert capacity_manifest['equivalence_passed']
    for name,expected in capacity_manifest['source_sha256'].items():assert digest(capacity_path/name)==expected,name
    original_config=json.loads((EXPERIMENT/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json').read_text())
    enlarged=json.loads((PROJECT/manifest['capacity_extension']['configuration']).read_text())
    assert enlarged['sensor']['max_candidates']==30_000_000
    enlarged['sensor']['max_candidates']=10_000_000
    assert enlarged==original_config
    for row in manifest['failed_attempts']:
        assert (PROJECT/row['root']).stat().st_size==row['root_bytes']
        assert (PROJECT/row['macro']).exists() and (PROJECT/row['log']).exists()
    assert len(manifest['runs'])==23 and len(summary['windows'])==23
    assert len({r['root'] for r in manifest['runs']})==23
    assert len({r['seed'] for r in manifest['runs']})==23
    cfg=settings();template=make_template(cfg);checked=[]
    records=json.loads((OUTPUT/'input-records.json').read_text())
    cached_checks=0
    for category in ('calibration_single','development_single'):
        for row in records[category]:
            for key,expected in row.get('cached_sha256',{}).items():
                assert digest(PROJECT/row[key])==expected,(category,row['root'],key)
                cached_checks+=1
    for row in records['development_pileup']:
        for name,expected in row['input_sha256'].items():assert digest(PROJECT/name)==expected,name
    frozen_at=datetime.fromisoformat(frozen['frozen_at'])
    prior_path=RECORD/'prechecked-roots.json'
    prior=json.loads(prior_path.read_text())['checked'] if prior_path.exists() else []
    prior={r['root']:r for r in prior}
    for native in manifest['runs']:
        assert native['status']=='completed' and datetime.fromisoformat(native['started'])>frozen_at
        assert digest(PROJECT/native['macro'])==native['macro_sha256']
        path=PROJECT/native['root'];stat=path.stat();previous=prior.get(native['root'])
        if previous and (previous['size'],previous['mtime_ns'])==(stat.st_size,stat.st_mtime_ns):
            assert previous['sha256']==native['root_sha256']
        else:
            assert digest(path)==native['root_sha256']
        with open_root(PROJECT/native['root']) as f:
            info=f['RunInfo'].arrays(library='np')
            assert int(info['random_seed'][0])==native['seed']
            assert int(info['sipm_avalanche_base_seed'][0])==native['sipm_seed']
            assert int(info['sipm_readout_base_seed'][0])==native['readout_seed']
            assert float(info['alpha_p11b_rate_hz'][0])==native['rate_hz']
            assert float(info['t_length_ps'][0])==20_000_000
            assert not info['brems_enabled'][0] and not info['c12_capture_enabled'][0]
            actual_config=json.loads(info['sipm_config_json'][0])
            config_path=PROJECT/native['config'] if 'config' in native else EXPERIMENT/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json'
            assert actual_config==json.loads(config_path.read_text())
            assert f['WaveformEvent'].num_entries==1 and f['SiPMEvent'].num_entries==1
        row=next(r for r in summary['windows'] if r['root']==native['root']);stem=Path(row['root']).stem
        metadata=json.loads((OUTPUT/'inputs/metadata'/f'{stem}.json').read_text())
        assert Path(metadata['root']).resolve()==(PROJECT/native['root']).resolve()
        assert metadata['root_sha256']==native['root_sha256']
        assert metadata['macro_contents']==(PROJECT/native['macro']).read_text()
        with np.load(OUTPUT/'inputs/truth'/f'{stem}.npz') as z:truth={k:z[k] for k in ('first_active_ns','emission_ns','charge_pe')}
        with np.load(OUTPUT/'inputs/waveforms'/f'{stem}.npz') as z:
            t=z['time_ns'];assert len(t)==10000 and np.isfinite(z['analog_mV']).all()
        for method in ('baseline','local'):
            fit=row['methods'][method]['fit'];rescored=score(fit,truth,[5000.,15000.])
            for name in ('target_count','reconstructed_count','matched_count','recovered_count','missed_count','unmatched_predictions','merge_flags','split_flags'):
                assert rescored[name]==row['methods'][method]['score'][name],(stem,method,name)
            if fit['status']=='fitted':
                comp=fit['all_components'] if method=='local' else fit
                times=np.asarray(comp['times_ns']);aa=np.asarray(comp['peak_mV'])
                model=10+template(t[:,None]-times[None,:])@aa
                with np.load(OUTPUT/'test/selected'/f'{stem}_{method}.npz') as z:np.testing.assert_allclose(model,z['model_mV'],atol=1e-8,rtol=1e-10)
            else:
                with np.load(OUTPUT/'test/selected'/f'{stem}_{method}.npz') as z:np.testing.assert_allclose(z['model_mV'],10.,atol=1e-12)
        checked.append(dict(root=native['root'],seed=native['seed'],passed=True,
            derived_sha256={kind:digest(OUTPUT/'inputs'/kind/f'{stem}.npz') for kind in ('waveforms','truth')}))
    save(RECORD/'data-validation.json',dict(passed=True,checked=checked,frozen_source_files=len(frozen['frozen_files']),
        waveform_models_checked=46,postfreeze_independent_runs=23,scoring_comparisons=46,
        calibration_and_development_cache_hashes_checked=cached_checks,
        numerical_tests_passed=checks['numerical-tests.json']['tests_run'],
        blind_replay_without_truth_passed=True,read_cache_16_arrays_exactly_equal=True,
        root_hash_policy='Full canonical-file SHA256; earlier closed-file checks reused only when size/mtime_ns remain identical',
        capacity_extension_equivalence_passed=True,capacity_extension_source_files=len(capacity_manifest['source_sha256']),
        aborted_attempts_preserved=len(manifest['failed_attempts']),
        scoring_redundant_field_recovery=recovery,blind_voltage_before_recovery_exactly_preserved=True,
        read_cache_transition=manifest.get('read_cache_transition')))


def delivery(summary):
    p=PROJECT/'jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb'
    before=nbformat.read(RECORD/'codex_0003_Alpha_Pileup_Rate_Fitting.ipynb.before',as_version=4)
    notebook=nbformat.read(p,as_version=4)
    notebook.cells=[c for c in notebook.cells if c.id not in ('local-alpha-progress-summary','local-alpha-provisional-20261020')]
    assert notebook==before,'Notebook changed during run; review before appending.'
    dev=json.loads((OUTPUT/'development-summary.json').read_text());frozen=json.loads((OUTPUT/'freeze.json').read_text())
    pen=json.loads((OUTPUT/'calibration/penalty.json').read_text());ctrl=json.loads((OUTPUT/'controls/controlled.json').read_text())
    shifts=json.loads((OUTPUT/'time-shift-controls.json').read_text())
    high_shift=next(r for r in shifts['rates'] if r['rate_hz']==1e8)
    usable_cal=[r for r in pen['rows'] if not r.get('no_candidate',False)]
    control_totals={}
    for kind in ('empty','single','double'):
        rows=[r for r in ctrl['rows'] if r['kind']==kind]
        control_totals[kind]={method:dict(targets=sum(r['scores'][method]['target_count'] for r in rows),
            recovered=sum(r['scores'][method]['recovered_count'] for r in rows),
            extra=sum(r['scores'][method]['unmatched_predictions'] for r in rows),
            exact_cases=sum(r['scores'][method]['recovered_count']==r['scores'][method]['target_count'] and r['scores'][method]['unmatched_predictions']==0 for r in rows)) for method in ('baseline','local')}
    lines=['| 源注入率/MHz | 新窗口数 | 原基线完整恢复 | 局部搜索完整恢复 | 原→新额外预测 | 原→新漏检 |',
           '| ---: | ---: | --- | --- | --- | --- |']
    for old,new in zip(summary['rates']['baseline'],summary['rates']['local']):
        def fraction(r):return f"{r['recovered_count']}/{r['target_count']}（{100*r['recovery_fraction']:.1f}%）" if r['target_count'] else '0/0'
        lines.append(f"| {old['rate_hz']/1e6:g} | {old['windows']} | {fraction(old)} | {fraction(new)} | {old['unmatched_predictions']}→{new['unmatched_predictions']} | {old['missed_count']}→{new['missed_count']} |")
    table='\n'.join(lines)
    metrics=conditional_metrics(summary)
    ten=next(r for r in metrics if r['rate_hz']==1e7)
    forty=next(r for r in metrics if r['rate_hz']==4e7)
    error_table=['| 源率/MHz | 原电荷误差均值±σ/% | 新电荷误差均值±σ/% | 新增正确事件 | 原正确但新失败 |',
                 '| ---: | --- | --- | ---: | ---: |']
    for row in metrics:
        def fmt(method):
            x=row[method];a=x['matched_charge_mean_percent'];b=x['matched_charge_std_percent']
            return f'{a:.2f} ± {b:.2f}' if a is not None and b is not None else '样本不足'
        error_table.append(f"| {row['rate_hz']/1e6:g} | {fmt('baseline')} | {fmt('local')} | {row['gained']} | {row['lost']} |")
    error_table='\n'.join(error_table)
    totals={}
    for method in ('baseline','local'):
        rr=summary['rates'][method];totals[method]={k:sum(r[k] for r in rr) for k in ('target_count','recovered_count','unmatched_predictions','missed_count','merge_flags','split_flags')}
    cells=[]
    def md(text,id):cells.append(nbformat.v4.new_markdown_cell(text.strip()+'\n',id=id))
    md(f'''
<a id="pure-alpha-local-model-selection-20261020"></a>
## 17. 主动拆峰、局部脉冲数搜索和多初值联合重拟合：实际实现

**研究问题**：落实第16节建议，在固定双衰减模板和完整动态范围下，比较原盲拟合与主动局部模型搜索能否恢复更多纯α事件。

**完成范围**：独立单α标定48窗、开发单α24窗；三个已公开的1/10/40 MHz窗仅作开发。冻结参数后生成23个新种子原生窗口，八档源注入率0.1/1/10/20/40/60/80/100 MHz，各8/3/2/2/2/2/2/2窗。另有174条固定均值模板软件控制。单窗20 μs、中央[5,15) μs评价；p–¹¹B有效单α源谱，无γ背景，原生共享微胞状态、2 ns采样和0.1 mV电子噪声保留。输入是完整含噪analog_mV；ADC位数、量化和削顶不参与拒绝。

**盲输入边界**：两种方法均未知真实K、到达时刻和电荷，已知α类型和独立固定模板95/351 ns及慢权重。新方法从原基线的观测波形拟合作初值；祖先/事件真值只在推断完成后打开评分。原始ROOT、旧结果和前16节不改。电荷为实际堆积微胞状态下的满恢复单胞等效量，尚未换算为沉积能量。

### 17.1 实现与标定

保留原基线，对所有重数都执行近邻局部联合时间优化；固定时刻时以非负最小二乘消去幅度，用多个时间初值搜索。主动试拆一个候选为两个，另有残差增峰、删除和近邻合并。起点最小间隔3 ns、局部调时±24 ns、两轮扫描和每窗最多160次增量是预算设置，不是器件分辨极限；默认局部优化使用1400 ns尾部，接受比较在完整记录电压域重新计算，保留外部脉冲尾部。

增加一个分量需满足局部全记录SSE改善大于 `factor × (0.5 mV² + c × Qref)`。Qref来自本轮局部拟合电荷并在该模型比较中固定，不读真值。独立单α标定的最大归一化拆峰收益给出 c={pen['coefficient_mV2_per_pe']:.8f} mV²/p.e.；开发比较factor=0.5/1/2，先最小化24个开发单α窗额外预测，再最大化三个开发率的完整恢复F1宏平均，平手选较大代价。本轮选择factor={dev['selected']['factor']:g}。这是经验选择规则，不是校准的χ²检验、p值或全局最优保证。

48个标定窗中，{len(usable_cal)}个有原基线候选并参与归一化拆峰收益计算，其拟合电荷范围{min(r['fitted_charge_pe'] for r in usable_cal):.1f}–{max(r['fitted_charge_pe'] for r in usable_cal):.1f} p.e.；其余无候选窗仍保留，不能把全部48窗称为48个阈上可用脉冲。该孤立标定不覆盖所有高占用下的形状/噪声条件。

内部允许100 p.e.以上干扰分量，最终仍只报告Q≥500 p.e.；这样可拟合阈下尾部，但也可能让弱分量分走目标电荷。快慢寿命和形状权重均固定，本轮不同时引入条件模板/GLS，便于识别搜索策略的效果。

开发首版因每次局部更新重建全模板矩阵而中断；已保存原源码及26个完成结果到development_v1记录，改为只更新局部模板列后重跑开发。该改动发生在冻结前；五项数值检查覆盖梯度、隐藏双脉冲、空/单误拆分、冷/热启动一致性及超过64分量的局部更新。

参数冻结时间：{frozen['frozen_at']}。实际Python/C++/配置快照、二进制身份、宏、三套种子和运行清单保存在新实验目录；源注入率不等于探测信号率或p–¹¹B反应率。

共享盘读取大ROOT时，mmap导出耗时过长；在同一60 MHz文件上验证冻结导出器的16个数组逐项一致后，末段只读导出改为读取普通本机字节副本。生产ROOT始终直接写入并保留在共享盘g4/data；复制时检查SHA256和源文件稳定性，元数据仍指向原始ROOT。该I/O调整不改变导出算法、物理、拟合源码或冻结参数，记录见read-cache-validation.json及manifest的read_cache_transition。

100 MHz、20 μs的首次生产触发了`max_candidates=10000000`计算保护，未保存完整响应，另一进程随后由控制器停止；两份不完整ROOT和日志均保留在failed_attempts，未计作恢复失败。独立alpha_local_capacity_20261020源码/配置快照只将处理候选数保护上限提高至30000000，C++/二进制和其他物理/读出参数不变；相同1 MHz源/SiPM/读出种子重放的16个数组逐项一致后，重试原计划的两份100 MHz种子。容量验证重放不加入23窗独立测试分母，拟合参数始终冻结。该上限是计算预算，不是ADC动态范围、物理分辨率或可接受计数率。
''','local-alpha-method-and-data')
    md('''
### 17.2 固定均值模板控制

174条：24空、24单脉冲，另7种间隔×6种幅比×3次电子噪声，共126条双脉冲；第一电荷4000 p.e.，第二/第一电荷比0.125–4。成功要求两个目标均满足时间±10 ns、电荷±20%，且无额外输出。生成器有真值，拟合不接收K/t/Q。这里未加入光子/微胞涨落，与原生性能严格分开；3次重复不足以确定分辨阈值。

幅比0.125时第二脉冲恰为500 p.e.报告阈值，噪声/拟合涨落可使输出跨阈，因此即使长间隔仍可能漏报；该曲线不能单独解释为双脉冲分辨能力。
''','local-alpha-controls-description')
    md(f'''原→新：空记录额外预测{control_totals['empty']['baseline']['extra']}→{control_totals['empty']['local']['extra']}；单脉冲完整恢复{control_totals['single']['baseline']['recovered']}/24→{control_totals['single']['local']['recovered']}/24，额外预测均为0。双脉冲合计252个目标，完整恢复事件数{control_totals['double']['baseline']['recovered']}→{control_totals['double']['local']['recovered']}；两事件都恢复且无额外输出的案例{control_totals['double']['baseline']['exact_cases']}/126→{control_totals['double']['local']['exact_cases']}/126。''','local-alpha-controls-counts')
    def image_cell(file,id):
        code=f'from pathlib import Path\nfrom IPython.display import display, Image\nproject = next(p for p in (Path.cwd(), *Path.cwd().parents) if (p / "pileup/analysis/alpha_local_fit.py").exists())\ndisplay(Image(filename=str(project / "jupyter/artifacts/alpha_local_model_selection_20261020/{file}")))'
        encoded=base64.b64encode((OUTPUT/file).read_bytes()).decode()
        c=nbformat.v4.new_code_cell(code,id=id,execution_count=1,
            outputs=[nbformat.v4.new_output('display_data',data={'image/png':encoded,'text/plain':'<IPython.core.display.Image object>'},metadata={})])
        # Execute the actual display expression with a collector; verify that
        # the resulting Image reads exactly the embedded artifact bytes.
        seen=[]
        def capture(obj):seen.append(obj)
        import IPython.display
        old_display=IPython.display.display;IPython.display.display=capture
        try:exec(compile(code,f'<{id}>','exec'),{})
        finally:IPython.display.display=old_display
        assert len(seen)==1 and seen[0].data==(OUTPUT/file).read_bytes()
        cells.append(c)
    image_cell('local_controls.png','local-alpha-figure-17-1')
    md('<a id="codex-figure-17-1"></a>\n**图 17-1：固定均值模板双脉冲的完整恢复控制**\n\n(a) 左：原基线；(b) 右：局部模型搜索。横轴为真间隔，曲线为电荷比；仅电子白噪声控制。','local-alpha-caption-17-1')
    md('### 17.3 新原生窗口的同窗比较\n\n'+table+'\n\n完整恢复沿用Q≥500 p.e.、一对一时间匹配±10 ns及电荷误差≤20%。漏检指未时间匹配；已匹配但电荷不合格另计，不能与完整恢复率互补。额外预测、误拆分/合并标志及电荷/时间误差均保留在逐窗JSON和summary。','local-alpha-native-results')
    md(error_table+'\n\n电荷均值和σ只对各方法时间匹配的样本计算，匹配集合可能不同；它们不是全体事件的能量偏差或能量分辨率。新增正确/原正确但新失败按同一真值事件逐一比较，避免只看净恢复数掩盖交换的事件。时间误差均值/σ另存report.json，同样受±10 ns匹配门截断。'
       +f"\n\n40 MHz新增{forty['gained']}个完整正确事件中，{forty['gained_previously_unmatched']}个原来无时间匹配，{forty['gained_charge_corrected']}个原来已匹配但电荷不合格。10 MHz丢失的{ten['lost']}个原正确事件中，{ten['lost_unmatched']}个变成无时间匹配，{ten['lost_charge_failed']}个仍匹配但电荷超出±20%；完整恢复的变化同时包含检出和电荷估计，不全部等同于新增发现脉冲。以上是测试后的描述性分解，没有用于选参。",'local-alpha-conditional-errors')
    image_cell('local_native_rates.png','local-alpha-figure-17-2')
    md('<a id="codex-figure-17-2"></a>\n**图 17-2：八档源注入率的独立原生恢复与失败比较**\n\n(a) 左上：完整恢复率；(b) 右上：未匹配输出数；(c) 左下：时间匹配后的电荷误差中位数及16–84百分位；(d) 右下：无时间匹配的漏检数。计数型指标随曝光和事件数变化，不直接跨率当比例比较。','local-alpha-caption-17-2')
    image_cell('local_native_examples.png','local-alpha-figure-17-3')
    md('<a id="codex-figure-17-3"></a>\n**图 17-3：预定10/40 MHz首个测试窗的总波形和逐事件电荷**\n\n(a)(b) 第一行10 MHz，(c)(d) 第二行40 MHz；左列为观测/拟合总电压，右列为真值及拟合事件电荷。固定显示记录5–7 μs，无按成功与否选图；真值只用于推断后显示。','local-alpha-caption-17-3')
    md(f'''
### 17.4 限制、验证和下一步

本轮按率报告，不把23个短窗中所有事件当独立重复来给窄置信区间。低率目标少，高率每档仅2窗；不能据此宣布最高可用率、通用最小间隔或实物性能。高密度下±10 ns的一对一时间匹配也有偶然配对成分。恢复的是电荷而非已标定沉积能量。

推断后另将全部重建时间错位−1739/−887/+887/+1739 ns，保持电荷不变，以原评分估计偶然配对的参考水平；5 μs保护区覆盖这些错位。100 MHz下四种错位的平均完整匹配比例为原基线{100*high_shift['baseline']['mean_fraction']:.2f}%、新方法{100*high_shift['local']['mean_fraction']:.2f}%，须与17.3的未错位成功率同时阅读。该对照不是无偏效率校正、置信区间或严格物理极限，也未用于调参。所有率的结果见time-shift-controls.json。

全部目标合计{totals['baseline']['target_count']}个，仅作这组非均衡曝光的描述：原→新完整恢复{totals['baseline']['recovered_count']}→{totals['local']['recovered_count']}，额外预测{totals['baseline']['unmatched_predictions']}→{totals['local']['unmatched_predictions']}，漏检{totals['baseline']['missed_count']}→{totals['local']['missed_count']}；误合并标志{totals['baseline']['merge_flags']}→{totals['local']['merge_flags']}、误拆分标志{totals['baseline']['split_flags']}→{totals['local']['split_flags']}。这些标志是操作定义，不是完备因果归类。

尚未达到更高恢复率的可能来源包括：最初候选遗漏、局部优化有限预算、近邻幅度/时间不可辨识、固定平均模板失配、相关的光子/雪崩涨落，以及弱干扰分量分走电荷。即使残差更小，逐事件电荷也可能变差；保留所有退步条件。下一步可在新的开发/测试划分下消融相关残差加权和独立条件形状标定，不用本轮新测试再选参。

五项数值检查、无真值输入的冷重放一致性、只读缓存及容量扩展的16个数组等价性，以及23个ROOT的RunInfo/实际宏/种子/哈希与冻结时序、46组重评分及模型重构见验证记录。图像及表仅加载已保存结果，未重跑旧计算。新图号17-1至17-3固定，原75个单元、原图及输出逐项保留。

其中一份100 MHz ROOT的`PrimaryPhoton.photon_arrival_t_ps`无法解压，独立重读仍失败；其余所需主光子分支及WaveformEvent/SiPMEvent可读。评分导出以全局WaveformEvent中同一光子的冗余时间，通过存储的x/y/能量浮点位精确连接恢复：要求键唯一、每条主光子记录精确且一对一对应、覆盖全部全局光子，否则拒绝恢复。先在另一份完整100 MHz ROOT上证明16个导出数组与原流程逐项一致，再应用于该文件；总波形的4个数组另与恢复前逐项核对。该适配仅恢复评分字段，原始ROOT、错误记录和原冻结导出源码保留，拟合仍只接收总电压；详见[冗余字段恢复验证](run_records/alpha_local_model_selection_20261020/redundant-time-validation.json)。该文件原分支仍不可直接读取，重放评分导出需使用保存的适配器或已验证的缓存。

[新测试完整结果](artifacts/alpha_local_model_selection_20261020/test-summary.json) · [冻结参数及源码哈希](artifacts/alpha_local_model_selection_20261020/freeze.json) · [原生运行清单](../pileup/experiments/alpha_local_model_selection_20261020/manifest.json) · [验证记录](run_records/alpha_local_model_selection_20261020/validation.json)

复用逻辑见[局部拟合实现](../pileup/analysis/alpha_local_fit.py)，下列命令从项目根目录调用保存的冻结源码重拟合，使用全新的输出目录保留本轮结果。软件控制将`test`替换为`controls`；只重画当前三张图可运行`python -m pileup.analysis.alpha_local_fit_report --figures-only`，同样通过项目Python入口调用。

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
  bash tools/run_gagg_python.sh \\
  pileup/experiments/alpha_local_model_selection_20261020/replay_fit.py \\
  test --workers 2 --output jupyter/artifacts/alpha_local_replay_new
```
''','local-alpha-limits-and-validation')
    for c in cells:
        if c.cell_type=='code':ast.parse(c.source)
        else:
            for link in re.findall(r'\]\(([^)]+)\)',c.source):assert not re.search(r'\.(pdf|ipynb|md)(?:#.*)?$',link,re.I)
    notice=nbformat.v4.new_markdown_cell('''**2026-10-20最新方法比较**：[第17节](#pure-alpha-local-model-selection-20261020)实现主动拆峰、局部脉冲数比较和多初值联合拟合，并在冻结后生成的23个新原生窗口上与原基线同窗比较。原第14节的36窗属于历史基线；评价新版请使用第17.3节的共同测试表，保留不同率的改善和退步。

新增固定图号：[图17-1：软件双脉冲控制](#codex-figure-17-1)、[图17-2：新原生逐率比较](#codex-figure-17-2)、[图17-3：预定波形与事件电荷](#codex-figure-17-3)。原有单元、图号及输出保持。\n''',id='local-alpha-current-summary')
    notebook.cells.insert(1,notice);notebook.cells.extend(cells);nbformat.validate(notebook);nbformat.write(notebook,p)
    actual=nbformat.read(p,as_version=4)
    old_ids={c.id for c in before.cells};assert [c for c in actual.cells if c.id in old_ids]==before.cells
    doc=PROJECT/'docs/literature-use.md'
    assert doc.read_bytes()==(RECORD/'literature-use.md.before').read_bytes()
    addition=f'''

<a id="pure-alpha-local-model-selection-20261020"></a>
### L40：2026-10-20——纯α局部脉冲数搜索的实际实现与独立测试

落实[L39](#pure-alpha-fit-methods-20261020)的主动拆峰、局部K比较与多初值联合时间/非负电荷重拟合；保持95/351 ns和原慢权重，完整含噪电压输入，无K/t/Q或源率真值输入。以独立48个单α窗最大拆峰SSE收益/Q标定代价系数{pen['coefficient_mV2_per_pe']:.8f}，24个开发单α及三个已公开1/10/40 MHz窗选择factor={dev['selected']['factor']:g}，于{frozen['frozen_at']}冻结后生成23个新原生测试窗。四篇文献条件/审核边界沿用L39；局部搜索/代价不是原论文逐式复现或校准的统计检验。

八档源率0.1/1/10/20/40/60/80/100 MHz分别8/3/2/2/2/2/2/2窗，20 μs/窗、中央10 μs。原生共享微胞、2 ns采样、0.1 mV噪声不变；Q≥500 p.e.、±10 ns和±20%完整恢复。逐率结果、额外预测、漏检及条件电荷误差见[Notebook0003第17节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#pure-alpha-local-model-selection-20261020)及[完整结果](../jupyter/artifacts/alpha_local_model_selection_20261020/test-summary.json)。非均衡曝光合计{totals['baseline']['target_count']}个目标，原→新完整恢复{totals['baseline']['recovered_count']}→{totals['local']['recovered_count']}、额外预测{totals['baseline']['unmatched_predictions']}→{totals['local']['unmatched_predictions']}、漏检{totals['baseline']['missed_count']}→{totals['local']['missed_count']}；不能替代逐率结果或窄区间推断。

另有174条匹配均值软件控制，仅电子白噪声，不当实物时间分辨率。保留开发首版源码及26个完成结果；局部列更新在冻结前通过五项数值检查。新[实验快照](../pileup/experiments/alpha_local_model_selection_20261020/README.md)含实际C++/配置/二进制身份、Python源码、宏和种子；[验证记录](../jupyter/run_records/alpha_local_model_selection_20261020/validation.json)核对冻结时序、ROOT/宏/RunInfo、46组重评分/模型及原Notebook输出。未改变物理、做GLS/条件形状消融或完成沉积能标；每个高率仅2窗，局部搜索不保证全局最优，低残差不等于正确事件数。

100 MHz首次因处理候选数保护上限10000000而中止，保留两份不完整ROOT/日志；在独立alpha_local_capacity_20261020快照仅将计算上限提高至30000000，以相同1 MHz种子的16数组逐项一致验证后重试原计划种子。此工程容量调整不改变物理/读出概率或冻结拟合参数，验证重放不计入测试分母；已知失败不混入恢复率。

一份100 MHz ROOT的主光子时间分支无法解压；以WaveformEvent的冗余时间及唯一的x/y/能量浮点位精确一对一连接恢复评分字段，要求覆盖全部光子。另一完整100 MHz文件的16数组等价验证及恢复前后4个总波形数组相同见[恢复验证](../jupyter/run_records/alpha_local_model_selection_20261020/redundant-time-validation.json)。原ROOT保留且该原分支仍不可直接读取，评分重导出使用保存的适配器；盲拟合输入、源码及参数不变。
'''
    doc.write_text(doc.read_text()+addition)
    validation=json.loads((RECORD/'data-validation.json').read_text())
    validation.update(original_notebook_cells=len(before.cells),new_notebook_cells=len(cells)+1,
        all_old_cells_and_outputs_exactly_preserved=True,figure_numbers=['17-1','17-2','17-3'],
        embedded_image_cells_executed_and_verified=True,nbformat_passed=True,
        document_old_prefix_preserved=doc.read_bytes().startswith((RECORD/'literature-use.md.before').read_bytes()),
        result_totals=totals,source_note='New tests never used to select parameters')
    save(RECORD/'validation.json',validation)
    save(OUTPUT/'report.json',dict(totals=totals,table_markdown=table,conditional_metrics=metrics,
        error_table_markdown=error_table,control_totals=control_totals,time_shift_rates=shifts['rates'],
        figures=['local_controls.png','local_native_rates.png','local_native_examples.png']))
    print(table,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--figures-only',action='store_true',help='Regenerate saved figures without editing notebooks or documents')
    p.add_argument('--validate-only',action='store_true')
    args=p.parse_args();summary=load()
    if args.figures_only:plots(summary);return
    validate(summary)
    if args.validate_only:return
    shifted_controls(summary);plots(summary);delivery(summary)


if __name__=='__main__':main()
