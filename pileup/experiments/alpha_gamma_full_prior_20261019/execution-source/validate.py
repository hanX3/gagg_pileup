from collections import defaultdict
from datetime import datetime
from pathlib import Path
import ast
import hashlib
import json
import re
import subprocess
import sys
import numpy as np
import nbformat

project=Path(__file__).resolve().parents[3];record=Path(__file__).resolve().parent
folder=project/'pileup/experiments/alpha_gamma_full_prior_20261019'
work=Path('/tmp/gagg-full-prior-20261019')
sys.path.insert(0,str(work/'source'));sys.path.insert(0,str(work/'dependency-source'))
from pileup.analysis.alpha_gamma_regularized_fit import CalibratedBank
from pileup.analysis.alpha_gamma_pileup_fit import score_mixed
from pileup.analysis.alpha_gamma_prior_study import measure
from pileup.analysis.alpha_rate_data import open_root

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text())
verified=parsed=0
for name in ('source-manifest.json','dependency-source-manifest.json','test-source-manifest.json',
             'native-source-manifest.json','postprocessing-source-manifest.json','postprocessing-source-manifest-v2.json',
             'postprocessing-source-manifest-v3.json','execution-source-manifest.json'):
    for relative,expected in load(folder/name)['files'].items():
        path=folder/relative;assert sha(path)==expected,relative;verified+=1
        if path.suffix=='.py':ast.parse(path.read_text());parsed+=1
protocol=load(folder/'protocol.json');settings=load(folder/'frozen-fit-settings.json')
assert sha(folder/'source-manifest.json')==protocol['source_manifest_sha256']
assert sha(folder/'protocol.json')==settings['development_protocol_sha256']
assert sha(work/'outputs/summary.json')==settings['development_summary_sha256']
assert sha(work/'models/unet/model.pt')==settings['neural_model_sha256']
assert sha(folder/'calibration.json')==settings['calibration_sha256']
assert all(value==0 for value in load(record/'inference-preflight.json')['dense_max_abs_differences'].values())
assert all(value==0 for value in load(record/'inference-preflight.json')['full_end_to_end_fit_max_abs_differences'].values())
assert load(record/'objective-state-check.json')['explicit_charge_penalty']==0.

manifest=load(folder/'manifest.json');runs=manifest['runs']
assert len(runs)==16 and all(r['status']=='completed' for r in runs)
assert {r['seed'] for r in runs}==set(range(261091200,261091216))
starts=[datetime.fromisoformat(r['started']) for r in runs]
assert all((b-a).total_seconds()>=1. for a,b in zip(starts,starts[1:]))
assert min(starts)>datetime.fromisoformat(settings['created'])
for row in runs:
    root=project/row['root'];stem=root.stem
    assert re.fullmatch(r'\d{8}_\d{2}h\d{2}m\d{2}s\.root',root.name)
    assert sha(root)==row['root_sha256']
    assert sha(project/row['macro'])==row['macro_sha256']
    with open_root(root) as data:
        assert all(data[name].num_entries==1 for name in ('RunInfo','WaveformEvent','SiPMEvent'))
    metadata=load(work/'test-inputs/metadata'/f'{stem}.json')
    assert metadata['seed']==row['seed'] and metadata['rate_hz']==row['rate_hz']
    slot='brems' if row['gamma_model']=='brems' else 'c12_capture'
    assert metadata['source_rates_hz'][slot]==row['gamma_rate_hz']
    assert metadata['macro_contents']==(project/row['macro']).read_text()
    assert metadata['gate_ns']==[5000.,15000.]
    wave=np.load(work/'test-inputs/waveforms'/f'{stem}.npz')
    assert len(wave['time_ns'])==10000 and 0.<=wave['time_ns'][0]<2.
    assert np.allclose(np.diff(wave['time_ns']),2.,atol=1e-8,rtol=0)
    assert np.isfinite(wave['analog_mV']).all()
    case=load(work/'test-cases'/f'{stem}.json')
    assert set(case)=={'stem','waveform','waveform_sha256'}
    assert sha(work/case['waveform'])==case['waveform_sha256']
    audit=load(work/'test-logs'/f'{stem}.file-access.json')
    assert audit['root_truth_guard_enabled'] and audit['fit_completed']
    assert not any(p.endswith('.root') or '/truth/' in p for p in audit['opened_files'])

cal=load(folder/'calibration.json');bank=CalibratedBank(cal)
model_errors=[];monotone=0;subproblem_status=defaultdict(int)
paths=list((work/'outputs/fits/class_only').glob('*.json'))
paths+=list((work/'outputs/fits/replay_none').glob('*.json'))+list((work/'outputs/fits/replay_both').glob('*.json'))
paths+=[p for name in ('hybrid_none','hybrid_both','class_only') for p in (work/'test-outputs/fits'/name).glob('*.json')]
assert len(paths)==70
for path in paths:
    item=load(path);fit=item['fit'];assert fit['status']=='fitted' and fit['fit_truth_inputs']==[]
    history=fit['weighted_cost_history'];assert np.isfinite(history).all() and np.all(np.diff(history)<=1e-6)
    monotone+=1
    for status,count in fit['optimizer_status_counts'].items():subproblem_status[status]+=count
    component=fit['all_components']
    names=component['particle_name'];kinds=[0 if name=='alpha' else 1 for name in names]
    assert set(names)<={'alpha','gamma'} and all(q>0 for q in component['charge_pe'])
    report=np.asarray(component['charge_pe'])>=500.
    assert fit['times_ns']==np.asarray(component['times_ns'])[report].tolist()
    assert fit['charge_pe']==np.asarray(component['charge_pe'])[report].tolist()
    saved=np.load(path.with_suffix('.npz'));t=saved['time_ns']
    model=bank.matrix(t,component['times_ns'],kinds,component['shape_charge_pe'])@np.asarray(component['charge_pe'])+cal['baseline_mV']
    error=float(np.max(abs(model-saved['model_mV'])))
    assert error<1e-8,error;model_errors.append(error)

rescored=0
for rows,output,native in ((protocol['scoring_rows'],work/'outputs',False),(runs,work/'test-outputs',True)):
    summary=load(output/'summary.json');by_method=defaultdict(list)
    for row in rows:
        stem=Path(row['root']).stem
        truthpath=work/'test-inputs/truth'/f'{stem}.npz' if native else project/row['truth']
        metapath=work/'test-inputs/metadata'/f'{stem}.json' if native else project/row['metadata']
        with np.load(truthpath) as data:
            truth={key:data[key] for key in ('first_active_ns','emission_ns','particle_name','charge_pe')}
        gate=load(metapath)['gate_ns']
        for method in summary['methods']:
            if native or method=='class_only':path=output/'fits'/method/f'{stem}.json'
            else:
                old={'baseline':'baseline','hybrid_none':'hybrid_unet_0','hybrid_both':'hybrid_unet_10'}[method]
                path=project/protocol['parent_artifacts']/'results/validation'/old/f'{stem}.json'
            fitted=load(path)['fit'];score=score_mixed(fitted,truth,gate)
            assert score==load(output/'scores'/method/f'{stem}.json')['score']
            by_method[method].append(dict(row=row,fit=fitted,truth=truth,score=score,gate_ns=gate));rescored+=1
    for stored in summary['pooled']:
        measured=measure(by_method[stored['method']],stored['particle'])
        assert measured=={key:stored[key] for key in measured}
        grouped=[r for r in summary['groups'] if (r['method'],r['particle'])==(stored['method'],stored['particle'])]
        assert sum(r['target_count'] for r in grouped)==stored['target_count']
        assert sum(r['recovered_count'] for r in grouped)==stored['recovered_count']
    events=load(output/'events.json')
    for pooled in summary['pooled']:
        count=sum(e['particle']==pooled['particle'] and e['outcomes'][pooled['method']]['outcome']=='recovered' for e in events)
        assert count==pooled['recovered_count']
assert rescored==160

artifact=project/'jupyter/artifacts/alpha_gamma_full_prior_20261019'
archive=load(record/'final-archive.json')
for relative,expected in archive['files'].items():assert sha(artifact/relative)==expected
notebook=project/'jupyter/codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb'
book=nbformat.read(notebook,as_version=4);before=nbformat.read(record/(notebook.name+'.before'),as_version=4)
nbformat.validate(book);assert book.cells[:49]==before.cells and len(before.cells)==49
new=book.cells[49:]
assert all(c.get('execution_count') is not None and not any(o.output_type=='error' for o in c.get('outputs',[])) for c in new if c.cell_type=='code')
links=0
for cell in new:
    if cell.cell_type!='markdown':continue
    for link in re.findall(r'\[[^\]]*\]\(([^)]+)\)',cell.source):
        target=link.split('#')[0]
        if not target or '://' in target:continue
        assert (notebook.parent/target).resolve().is_file(),target;links+=1
catalog=nbformat.read(project/'jupyter/Run_Catalog.ipynb',as_version=4)
oldcatalog=nbformat.read(record/'Run_Catalog.ipynb.before',as_version=4)
nbformat.validate(catalog);assert catalog.cells[:-2]==oldcatalog.cells
subprocess.run(['git','-c','core.fileMode=false','diff','--check','--','README.md','docs/literature-use.md'],cwd=project,check=True)
result=dict(status='passed',created=datetime.now().astimezone().isoformat(),source_entries_verified=verified,python_files_parsed=parsed,
    native_roots_verified=16,all_test_fits_guarded=16,blind_method_results=80,score_records_recomputed=rescored,
    monotone_hybrid_fits=monotone,physical_models_reconstructed=len(model_errors),max_model_abs_difference=max(model_errors),
    nonlinear_subproblem_status_counts=dict(subproblem_status),budget_status_zero_is_not_global_convergence=True,
    unit_tests=dict(passed=3,scope='Previously run on the identical snapshotted full-prior code'),
    original_notebook_cells_preserved=49,new_notebook_cells=len(new),new_local_links_checked=links,
    original_catalog_cells_preserved=len(oldcatalog.cells),notebook_sha256=sha(notebook),
    development_summary_sha256=sha(artifact/'development/summary.json'),test_summary_sha256=sha(artifact/'test/summary.json'),
    settings_sha256=sha(folder/'frozen-fit-settings.json'),
    limits=['Finite search budgets and 4 native windows per condition','No real detector, deposited-energy or domain-shift validation',
            'Whole convenience replay wrappers not rerun; saved actual runners, endpoint controls and all rescores verified'])
(record/'validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result),flush=True)
