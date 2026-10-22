from pathlib import Path
from datetime import datetime
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys

project=Path(__file__).resolve().parents[3]
record=Path(__file__).resolve().parent
folder=project/'pileup/experiments/alpha_gamma_full_prior_20261019'
work=Path('/tmp/gagg-full-prior-20261019')
prior=project/'pileup/experiments/alpha_gamma_prior_ablation_20261019'
parent=json.loads((prior/'protocol.json').read_text())
folder.mkdir(parents=True,exist_ok=False);work.mkdir(exist_ok=False)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

for name in ('README.md','docs/literature-use.md','jupyter/codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb'):
    shutil.copy2(project/name,record/(Path(name).name+'.before'))
status=subprocess.check_output(['git','-c','core.fileMode=false','status','--short'],cwd=project,text=True)
save(record/'git-status-before.json',dict(status=status))
snapshot={}
for path in sorted((project/'pileup/analysis').glob('*.py')):
    target=folder/'source'/path.relative_to(project)
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
    snapshot[str(target.relative_to(folder))]=sha(target)
cal=project/'pileup/experiments/alpha_gamma_ml_20261016/calibration.json'
target=folder/'source'/cal.relative_to(project);target.parent.mkdir(parents=True,exist_ok=True)
shutil.copy2(cal,target);snapshot[str(target.relative_to(folder))]=sha(target)
shutil.copy2(cal,folder/'calibration.json')
shutil.copytree(prior/'dependency-source',folder/'dependency-source')
shutil.copy2(prior/'dependency-source-manifest.json',folder/'dependency-source-manifest.json')
save(folder/'source-manifest.json',dict(created=datetime.now().astimezone().isoformat(),files=snapshot))
shutil.copytree(prior/'cases',folder/'cases')
reference={}
for method,old in {'baseline':'baseline','hybrid_none':'hybrid_unet_0','hybrid_both':'hybrid_unet_10'}.items():
    reference[method]={Path(row['root']).stem:sha(project/parent['parent_artifacts']/'results/validation'/old/(Path(row['root']).stem+'.json')) for row in parent['scoring_rows']}
sys.path.insert(0,str(folder/'source'));sys.path.insert(0,str(folder/'dependency-source'))
from dataclasses import asdict
from pileup.analysis.alpha_gamma_full_prior_fit import FullPriorConfig
protocol=dict(created=datetime.now().astimezone().isoformat(),state='fixed_before_full_pipeline_development_results',
    git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=project,text=True).strip(),
    git_status_before=status, question='Does removing only the explicit Q penalty help both classes in the full floating-time and residual-addition hybrid?',
    parent_experiment=parent['parent_experiment'],parent_artifacts=parent['parent_artifacts'],
    ablation_parent=str(prior.relative_to(project)),
    data_role='20 disclosed development windows; old final tests not reused for selection',
    scoring_rows=parent['scoring_rows'],calibration_sha256=sha(folder/'calibration.json'),
    neural_model_sha256='40ecf7d65b33df04f85381c808b407e6829e2ea6c4e9664d848c726d7ad5b7a9',
    fit_config=asdict(FullPriorConfig()),baseline_mV=json.loads(cal.read_text())['baseline_mV'],
    reference_sha256=reference,source_manifest_sha256=sha(folder/'source-manifest.json'),
    bootstrap_seed=261091101,workers=4,threads_per_worker=1,
    development_gate=dict(reference='baseline',pooled_recovery='Each particle recovered count >= baseline',
        precision_margin_absolute=.02,extra_multiplier=1.25,extra_offset=5,
        cell_recovery_margin_absolute=.05,cell_minimum_targets=20,
        interpretation='Practical predeclared screening tolerances, not evidence of statistical equivalence',
        action_if_failed='Record failure and diagnose; do not spend fresh final-test data on this unchanged candidate',
        action_if_passed='Freeze candidate, then generate new independent native test windows before scoring'),
    changes=dict(charge_strength=[10,0],class_strength=[10,10],
        retained=['Neural+blind candidate union','Neural Q initialization and Q-dependent shape',
                  'Frozen baseline whitening','Two time/type/Q passes','Up to 8 residual additions',
                  'Complete analog voltage and adequate dynamic range']),
    blind_case_policy='Packets contain waveform, cached blind fit/network outputs and hashes; truth opened only after fits',
    new_g4=False,training=False)
save(folder/'protocol.json',protocol)
wrapper='''from pathlib import Path
import runpy
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
sys.path.insert(0,str(folder/'dependency-source'))
runpy.run_module('pileup.analysis.alpha_gamma_full_prior_study',run_name='__main__')
'''
(folder/'run_analysis.py').write_text(wrapper)
(folder/'.gitignore').write_text('__pycache__/\n*.pyc\n')
for name in ('source','dependency-source','cases'):
    shutil.copytree(folder/name,work/name)
for name in ('calibration.json','run_analysis.py'):
    shutil.copy2(folder/name,work/name)
local_prior=Path('/tmp/gagg-prior-ablation-20261019')
for path in sorted((folder/'cases').glob('*.json')):
    case=json.loads(path.read_text());stem=path.stem
    for item in case['inputs'].values():
        source=local_prior/item['path']
        if not source.is_file():source=project/'jupyter/artifacts/alpha_gamma_prior_ablation_20261019/blind-inputs'/stem/Path(item['path']).name
        assert sha(source)==item['sha256']
        target=work/item['path'];target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
save(record/'preparation.json',dict(created=datetime.now().astimezone().isoformat(),work=str(work),
    experiment=str(folder.relative_to(project)),protocol_sha256=sha(folder/'protocol.json'),
    runtime=dict(executable=sys.executable,python=platform.python_version(),platform=platform.platform(),
        versions={name:importlib.metadata.version(name) for name in ('numpy','scipy','matplotlib','uproot','nbformat')}),
    frozen_source_entries=len(snapshot),old_notebook_sha256=sha(record/'codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb.before')))
print(json.dumps(dict(work=str(work),protocol_sha256=sha(folder/'protocol.json'),source_entries=len(snapshot))),flush=True)
