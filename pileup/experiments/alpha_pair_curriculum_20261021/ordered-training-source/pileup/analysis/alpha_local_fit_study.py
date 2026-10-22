"""Calibrate/freeze and evaluate local pure-alpha fits, from the project root."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, replace
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import time

import numpy as np

from .alpha_local_fit import LocalConfig, fit_local, optimize_local
from .alpha_pileup_fit import FitConfig, fit_waveform, score
from .alpha_full_signal_study import make_template, summarize
from .alpha_rate_campaign import digest, save

PROJECT = Path(__file__).resolve().parents[2]
EXPERIMENT = PROJECT/'pileup/experiments/alpha_local_model_selection_20261020'
OUTPUT = PROJECT/'jupyter/artifacts/alpha_local_model_selection_20261020'
RECORD = PROJECT/'jupyter/run_records/alpha_local_model_selection_20261020'


def settings():
    return json.loads((OUTPUT/'baseline-settings.json').read_text())


def prepare():
    assert not OUTPUT.exists()
    for p in (OUTPUT, RECORD, EXPERIMENT): p.mkdir(parents=True, exist_ok=True)
    for name in ('calibration', 'development', 'controls', 'test', 'inputs'):
        (OUTPUT/name).mkdir()
    for name in ('jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb',
                 'jupyter/Run_Catalog.ipynb', 'docs/literature-use.md', 'README.md', 'pileup/README.md'):
        shutil.copyfile(PROJECT/name, RECORD/(Path(name).name+'.before'))
    old = PROJECT/'jupyter/artifacts/alpha_rate_fit_20261013_211541/analysis_config.json'
    shutil.copyfile(old, OUTPUT/'baseline-settings.json')
    groups = json.loads((PROJECT/'jupyter/artifacts/alpha_gamma_shared_tau_20261020/input-records.json').read_text())
    cal = sorted([r for r in groups['calibration'] if r['particle']=='alpha'], key=lambda r:r['seed'])
    dev = sorted([r for r in groups['development'] if r['particle']=='alpha'], key=lambda r:r['seed'])
    rng = np.random.default_rng(262023001)
    cal = [cal[i] for i in sorted(rng.choice(len(cal), min(48,len(cal)), replace=False))]
    dev = [dev[i] for i in sorted(rng.choice(len(dev), min(24,len(dev)), replace=False))]
    old_summary = json.loads((PROJECT/'jupyter/artifacts/alpha_full_signal_20261014_143000/summary.json').read_text())
    mixed = [next(r for r in old_summary['windows'] if r['rate_hz']==rate) for rate in (1e6,1e7,4e7)]
    save(OUTPUT/'input-records.json', dict(calibration_single=cal,development_single=dev,
        development_pileup=mixed, exposed_regression_windows=36,
        selection='Deterministic seed 262023001; pure-alpha development rates 1/10/40 MHz first predeclared old window'))
    save(OUTPUT/'protocol.json', dict(created=datetime.now().astimezone().isoformat(),
        purpose='Pure alpha local K search; preserve fixed double-decay template and full voltage input',
        candidate_penalty_multipliers=[.5,1.,2.],
        selection='Lowest extra predictions on native development singles first; then maximal macro complete-recovery F1 over three old development rate windows; higher penalty wins ties',
        scoring=dict(min_charge_pe=500,time_tolerance_ns=10,relative_charge_tolerance=.2),
        calibration='Max observed positive split SSE gain / fitted charge; independently labelled single-alpha calibration only',
        test_plan={'100000':8,'1000000':3,'10000000':2,'20000000':2,'40000000':2,'60000000':2,'80000000':2,'100000000':2},
        software_control='Frozen mean template + 0.1 mV white noise; separate from native detector performance',
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        environment=dict(python=platform.python_version(),platform=platform.platform(),
            packages={k:importlib.metadata.version(k) for k in ('numpy','scipy','uproot','awkward','matplotlib')})))
    (RECORD/'git-status-before.txt').write_text(subprocess.check_output(
        ['git','-c','core.fileMode=false','status','--short','--untracked-files=no'],cwd=PROJECT,text=True))
    # Retain the actual existing edits as context without committing them.
    (RECORD/'working-changes.patch').write_bytes(subprocess.check_output(
        ['git','-c','core.fileMode=false','diff','--binary'],cwd=PROJECT))
    print('PREPARED',len(cal),'calibration singles',len(dev),'development singles',flush=True)


def load_voltage(row, *, crop_single=False):
    with np.load(PROJECT/row['waveform']) as a:
        t,y=a['time_ns'].copy(),a['analog_mV'].copy()
    if crop_single:
        keep=(t>=3000)&(t<9000);t,y=t[keep]-3000,y[keep]
    return t,y


def calibrate():
    d=settings();template=make_template(d);base=FitConfig(**d['fit_config']);cfg=LocalConfig()
    rows=json.loads((OUTPUT/'input-records.json').read_text())['calibration_single']
    result=[]
    for row in rows:
        t,y=load_voltage(row,crop_single=True)
        initial=fit_waveform(t,y,template,base)
        if not initial['times_ns']:
            result.append(dict(root=row['root'],no_candidate=True));continue
        center=initial['times_ns'][int(np.argmax(initial['charge_pe']))]
        lo,hi=max(t[0],center-24),min(t[-1],center+24)
        mask=(t>=lo-4)&(t<=hi+cfg.tail_ns)
        one,_=optimize_local(t[mask],y[mask]-10,[[center]],template,lo,hi,cfg)
        seeds=[]
        for gap in cfg.split_gaps_ns:
            seeds.extend(([center-gap/2,center+gap/2],[center,center+gap],[center-gap,center]))
        two,_=optimize_local(t[mask],y[mask]-10,seeds,template,lo,hi,cfg)
        gain=0.
        if two is not None and np.min(two['amplitude']/(base.area_per_pe_mV_ns*template.unit_area_peak))>=cfg.component_charge_pe:
            r1=y-10-template(t[:,None]-one['times'])@one['amplitude']
            r2=y-10-template(t[:,None]-two['times'])@two['amplitude']
            gain=max(0.,float(r1@r1-r2@r2))
        q=max(500.,float(np.sum(one['amplitude'])/(base.area_per_pe_mV_ns*template.unit_area_peak)))
        result.append(dict(root=row['root'],seed=row['seed'],fitted_charge_pe=q,
            best_split_gain_mV2=gain,normalized_gain_mV2_per_pe=gain/q))
        print('CALIBRATION',len(result),round(q,1),round(gain/q,6),flush=True)
    coefficient=max([r.get('normalized_gain_mV2_per_pe',0) for r in result]+[1e-6])
    save(OUTPUT/'calibration/penalty.json',dict(rows=result,coefficient_mV2_per_pe=coefficient,
        rule='Maximum independent-single split gain divided by fitted charge, floor 1e-6; no pileup test truth used',
        baseline_settings_sha256=digest(OUTPUT/'baseline-settings.json')))


def common_fit(row, config, *, single=False, prefix='development', tag='selected'):
    folder=OUTPUT/prefix/tag;folder.mkdir(parents=True,exist_ok=True)
    dest=folder/(Path(row['root']).stem+'.json')
    if dest.exists():return json.loads(dest.read_text())
    d=settings();template=make_template(d);base=FitConfig(**d['fit_config'])
    if single:
        t,y=load_voltage(row,crop_single=True);gate=[0.,6000.];shift=3000.
    else:
        source=Path(row['source']) if 'source' in row else OUTPUT/'inputs'
        if not source.is_absolute():source=PROJECT/source
        stem=Path(row['root']).stem
        with np.load(source/'waveforms'/f'{stem}.npz') as a:t,y=a['time_ns'].copy(),a['analog_mV'].copy()
        row={**row,'truth':str((source/'truth'/f'{stem}.npz').relative_to(PROJECT))}
        gate=[5000.,15000.];shift=0.
    cache=OUTPUT/'baseline-cache';cache.mkdir(exist_ok=True)
    stem=Path(row['root']).stem+('_single' if single else '')
    cached=cache/(stem+'.json')
    if cached.exists():
        previous=json.loads(cached.read_text());initial=previous['fit'];baseline_seconds=previous['seconds']
        with np.load(cache/(stem+'.npz')) as z:initial['model_mV']=z['model_mV'].copy()
    else:
        tick=time.monotonic();initial=fit_waveform(t,y,template,base);baseline_seconds=time.monotonic()-tick
        frozen_initial={k:v for k,v in initial.items() if k!='model_mV'}
        np.savez_compressed(cache/(stem+'.npz'),model_mV=initial['model_mV'])
        save(cached,dict(fit=frozen_initial,seconds=baseline_seconds))
    tick=time.monotonic();new=fit_local(t,y,template,base,config,initial=initial);seconds=time.monotonic()-tick
    # Deliberately open scoring truth only after both blind fits are finished.
    with np.load(PROJECT/row['truth']) as z:
        truth={k:z[k].copy() for k in ('first_active_ns','emission_ns','charge_pe')}
    truth['first_active_ns']-=shift;truth['emission_ns']-=shift
    methods={}
    for key,fitted,elapsed in [('baseline',initial,baseline_seconds),('local',new,baseline_seconds+seconds)]:
        model=fitted.pop('model_mV')
        np.savez_compressed(folder/(Path(row['root']).stem+f'_{key}.npz'),model_mV=model)
        methods[key]=dict(fit=fitted,score=score(fitted,truth,gate),seconds=elapsed)
    result=dict(root=row['root'],seed=row.get('seed'),rate_hz=row.get('rate_hz',0),
        waveform=row.get('waveform'),truth=row['truth'],single=single,methods=methods,
        refinement_seconds=seconds,input_sha256=dict(voltage=digest(PROJECT/row['waveform'])) if single else {},
        inference_inputs='time_ns, total analog_mV, independent fixed template/config; baseline observed-waveform warm start')
    save(dest,result)
    print('FIT',prefix,tag,row['root'],[(k,v['score']['recovered_count'],v['score']['target_count'],v['score']['unmatched_predictions']) for k,v in methods.items()],round(seconds,2),flush=True)
    return result


def develop(workers=2):
    data=json.loads((OUTPUT/'input-records.json').read_text());p=json.loads((OUTPUT/'calibration/penalty.json').read_text())
    candidates=[]
    for factor in (.5,1.,2.):
        cfg=replace(LocalConfig(),addition_cost_mV2_per_pe=p['coefficient_mV2_per_pe'],penalty_multiplier=factor)
        tag=f'penalty_{factor:g}'
        tasks=[(r,True) for r in data['development_single']]+[(r,False) for r in data['development_pileup']]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures=[pool.submit(common_fit,r,cfg,single=single,tag=tag) for r,single in tasks]
            rows=[f.result() for f in as_completed(futures)]
        single_rows=[r for r in rows if r['single']];pileup=[r for r in rows if not r['single']]
        extras=sum(r['methods']['local']['score']['unmatched_predictions'] for r in single_rows)
        f1=[]
        for r in pileup:
            z=r['methods']['local']['score'];n=z['target_count'];rec=z['recovered_count'];pred=z['reconstructed_count']
            f1.append(2*rec/(n+pred) if n+pred else 1.)
        candidates.append(dict(factor=factor,tag=tag,config=asdict(cfg),single_extra_predictions=extras,
            macro_complete_recovery_f1=float(np.mean(f1)),rows=len(rows)))
    selected=max(candidates,key=lambda x:(-x['single_extra_predictions'],x['macro_complete_recovery_f1'],x['factor']))
    save(OUTPUT/'development-summary.json',dict(candidates=candidates,selected=selected,rule='Predeclared protocol; no new test inspected'))


def freeze():
    chosen=json.loads((OUTPUT/'development-summary.json').read_text())['selected']
    folder=EXPERIMENT/'inference-source'
    assert not folder.exists()
    originals=list((PROJECT/'pileup/analysis').glob('*.py'))+list((PROJECT/'analysis').glob('*.py'))
    for path in originals:
        dest=folder/path.relative_to(PROJECT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
    for package in ('analysis','pileup','pileup/analysis'):
        (folder/package/'__init__.py').touch(exist_ok=True)
    files=[OUTPUT/'baseline-settings.json',OUTPUT/'calibration/penalty.json',OUTPUT/'protocol.json',OUTPUT/'input-records.json']
    files+=list(folder.rglob('*.py'))+originals
    files+=list((EXPERIMENT/'source').rglob('*'))+[EXPERIMENT/'bin/gagg']
    files=[p for p in files if p.is_file()]
    record=dict(frozen_at=datetime.now().astimezone().isoformat(),selected=chosen,
        frozen_files={str(p.relative_to(PROJECT)):digest(p) for p in files},
        source='Actual implementation and dependencies copied; no checkout needed for replay')
    save(OUTPUT/'freeze.json',record);save(EXPERIMENT/'frozen-fit-settings.json',record)
    print('FROZEN',record['frozen_at'],chosen['factor'],flush=True)


def controls():
    f=json.loads((OUTPUT/'freeze.json').read_text());cfg=LocalConfig(**f['selected']['config'])
    d=settings();template=make_template(d);base=FitConfig(**d['fit_config'])
    rng=np.random.default_rng(262023002);t=np.arange(0.,6000.,2.)+.137;rows=[]
    plans=[('empty',[],[],0.,0.)]*24
    plans += [('single',[1000.3],[q],0.,0.) for q in (500.,1000.,4000.,16000.) for _ in range(6)]
    plans += [('double',[1000.3,1000.3+gap],[4000.,4000.*ratio],gap,ratio)
              for gap in (4.,8.,12.,20.,40.,80.,200.) for ratio in (.125,.25,.5,1.,2.,4.) for _ in range(3)]
    for i,(kind,onsets,charges,gap,ratio) in enumerate(plans):
        a=np.asarray(charges)*base.area_per_pe_mV_ns*template.unit_area_peak
        y=10+template(t[:,None]-np.asarray(onsets)[None,:])@a+rng.normal(0,.1,len(t))
        initial=fit_waveform(t,y,template,base);new=fit_local(t,y,template,base,cfg,initial=initial)
        truth=dict(first_active_ns=np.array(onsets),emission_ns=np.array(onsets),charge_pe=np.array(charges))
        rows.append(dict(kind=kind,gap_ns=gap,ratio=ratio,scores={name:score(fit,truth,[500.,4000.]) for name,fit in [('baseline',initial),('local',new)]}))
    save(OUTPUT/'controls/controlled.json',dict(seed=262023002,rows=rows,
        limitation='Matched mean-template software controls only, white electronic noise; not native shared-cell double-pulse resolution'))
    print('CONTROLS',len(rows),flush=True)


def test(workers=2):
    frozen=json.loads((OUTPUT/'freeze.json').read_text())
    for name,expected in frozen['frozen_files'].items():assert digest(PROJECT/name)==expected,name
    cfg=LocalConfig(**frozen['selected']['config'])
    manifest=json.loads((EXPERIMENT/'manifest.json').read_text())
    assert all(r['status']=='completed' for r in manifest['runs'])
    assert len(manifest['runs'])==23
    with ThreadPoolExecutor(max_workers=workers) as pool:
        fs=[pool.submit(common_fit,r,cfg,prefix='test') for r in manifest['runs']]
        rows=[f.result() for f in as_completed(fs)]
    rates={}
    for method in ('baseline','local'):
        packaged=[]
        for row in rows:
            m=row['methods'][method]
            packaged.append(dict(rate_hz=row['rate_hz'],fit=m['fit'],score=m['score']))
        rates[method]=summarize(packaged)
    save(OUTPUT/'test-summary.json',dict(rates=rates,windows=sorted(rows,key=lambda r:(r['rate_hz'],r['seed'])),
        frozen_at=frozen['frozen_at'],test_windows=23,cohort='New native seeds, generated after algorithm freeze'))


def main():
    global OUTPUT
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','calibrate','develop','freeze','controls','test']);p.add_argument('--workers',type=int,default=2)
    p.add_argument('--output',type=Path,help='New project-relative output for frozen test/control replay; preserves original results')
    args=p.parse_args();assert 1<=args.workers<=4
    if args.output is not None:
        assert args.stage in ('test','controls'),'Only frozen test/control replay accepts --output'
        target=(PROJECT/args.output).resolve();target.relative_to(PROJECT)
        assert not target.exists(),'Replay output must be new'
        target.mkdir(parents=True)
        for name in ('freeze.json','baseline-settings.json'):
            shutil.copyfile(OUTPUT/name,target/name)
        for name in ('inputs','calibration'):
            shutil.copytree(OUTPUT/name,target/name)
        (target/'test').mkdir();(target/'controls').mkdir();OUTPUT=target
    if args.stage in ('develop','test'):globals()[args.stage](args.workers)
    else:globals()[args.stage]()


if __name__=='__main__':main()
