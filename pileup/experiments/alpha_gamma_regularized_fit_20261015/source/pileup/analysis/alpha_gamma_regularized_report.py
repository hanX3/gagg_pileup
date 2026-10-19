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


def metrics(rows,particle):
    data=[r['score']['by_particle'][particle] for r in rows]
    keys=('target_count','matched_count','correct_type_count','recovered_count','missed_count','wrong_type_count',
          'central_injected','central_nonzero_charge')
    counts={k:sum(r[k] for r in data) for k in keys}
    predictions=sum(sum(p==particle and r['gate_ns'][0]<=t<r['gate_ns'][1]
                        for p,t in zip(r['fit']['particle_name'],r['fit']['times_ns'])) for r in rows)
    errors=np.array([x for d in data for x in d['relative_charge_errors']])
    timing=np.array([x for d in data for x in d['time_errors_ns']])
    n=counts['target_count'];correct=counts['recovered_count']
    return dict(**counts,prediction_count=predictions,full_recovery=correct/n if n else None,
        type_purity=counts['correct_type_count']/predictions if predictions else None,
        full_precision=correct/predictions if predictions else None,
        full_f1=2*correct/(n+predictions) if n+predictions else None,
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
                        gamma_injection_MHz=10.,particle=particle,windows=2,**metrics(subset,particle)))
    result=dict(groups=groups,methods=methods,selected_penalty=penalty,heldout_windows=8,
        definition='Q truth >=500 pe, class-blind one-to-one time match <=10 ns, correct type, charge error <=20%',
        error_conditioning='Time-matched targets, including misclassified events; charge is not deposited energy',
        fit_seconds={m:sum(r['fit_seconds'] for r in all_rows if r['method']==m) for m in methods},
        diagnostics={m:{k:sum(r['score'][k] for r in all_rows if r['method']==m)
                     for k in ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags')}
                     for m in methods},
        source_sha256=frozen['fitter_sha256'],calibration_sha256=frozen['calibration_sha256'],
        limitations=['Two windows per condition; no precise efficiency or maximum-rate claim',
            'Central isolated calibration, only 13 nonzero gamma calibration examples; no occupancy/position calibration',
            'AR(1) covariance and charge-dependent mean template are approximations',
            'Only sparse support is selected; discrete/time refinements are local, not guaranteed global',
            'No calibrated classification posterior or rejection policy; no deposited-energy calibration'])
    save(artifacts/'summary.json',result);save(artifacts/'test-window-results.json',all_rows)
    print(json.dumps(groups,ensure_ascii=False,indent=2))
    return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('freeze','summarize'));p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--artifacts',type=Path,required=True);a=p.parse_args()
    (freeze if a.stage=='freeze' else summarize)(a.folder,a.artifacts)
