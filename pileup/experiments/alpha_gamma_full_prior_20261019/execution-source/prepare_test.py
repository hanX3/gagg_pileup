from dataclasses import asdict
from datetime import datetime
from pathlib import Path
import hashlib
import json
import shutil
import sys

project=Path(__file__).resolve().parents[3];record=Path(__file__).resolve().parent
folder=project/'pileup/experiments/alpha_gamma_full_prior_20261019'
parent=project/'pileup/experiments/alpha_gamma_ml_20261016'
work=Path('/tmp/gagg-full-prior-20261019')
summary=json.loads((work/'outputs/summary.json').read_text());assert summary['development_gate']['passed']
protocol=json.loads((folder/'protocol.json').read_text())
old=json.loads((parent/'manifest.json').read_text())
sys.path.insert(0,str(work/'source'));sys.path.insert(0,str(work/'dependency-source'))
from pileup.analysis.alpha_gamma_full_prior_fit import FullPriorConfig
from pileup.analysis.alpha_gamma_ml_hybrid import MixedHybridConfig

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')

test_source=folder/'test-source';test_source.mkdir(exist_ok=False)
module=project/'pileup/analysis/alpha_gamma_full_prior_test.py'
shutil.copy2(module,test_source/module.name)
save(folder/'test-source-manifest.json',dict(files={'test-source/'+module.name:sha(test_source/module.name)}))
native_files={}
(folder/'bin').mkdir(exist_ok=False)
for model in ('brems','mono662'):
    source=parent/f'source/{model}/g4';destination=folder/f'source/{model}/g4'
    for path in source.rglob('*'):
        if path.is_file():assert sha(path)==old['source_sha256'][str(path.relative_to(parent))]
    shutil.copytree(source,destination)
    for path in destination.rglob('*'):
        if path.is_file():native_files[str(path.relative_to(folder))]=sha(path)
    assert sha(parent/f'bin/gagg-{model}')==old['parents'][model]['binary_sha256']
    shutil.copy2(parent/f'bin/gagg-{model}',folder/f'bin/gagg-{model}')
save(folder/'native-source-manifest.json',dict(parent=str(parent.relative_to(project)),files=native_files))
(folder/'.gitignore').write_text('__pycache__/\n*.pyc\n/bin/\n')
known=set()
for path in (project/'pileup/experiments').glob('*/manifest.json'):
    obj=json.loads(path.read_text())
    for item in obj.get('runs',[])+obj.get('planned_runs',[]):
        for name in ('seed','sipm_seed','readout_seed'):
            if name in item:known.add(item[name])
plans=[]
for model in ('brems','mono662'):
    for rate in (1000000,10000000):
        for i in range(4):
            seed=261091200+len(plans)
            assert not ({seed,seed+100000,seed+200000}&known)
            plans.append(dict(label=f'test_full_prior_{model}_{rate}_{i:02}',role='test_mixed',
                gamma_model=model,rate_hz=rate,gamma_rate_hz=10000000,seed=seed,preflight=i==0))
settings=dict(created=datetime.now().astimezone().isoformat(),state='frozen_before_new_test_generation',
    development_protocol_sha256=sha(folder/'protocol.json'),development_summary_sha256=sha(work/'outputs/summary.json'),
    development_gate=protocol['development_gate'],baseline_mV=protocol['baseline_mV'],
    calibration_sha256=protocol['calibration_sha256'],neural_model_sha256=protocol['neural_model_sha256'],
    neural_thresholds=[.4,.55],source_manifest_sha256=sha(folder/'source-manifest.json'),
    test_source_manifest_sha256=sha(folder/'test-source-manifest.json'),
    native_source_manifest_sha256=sha(folder/'native-source-manifest.json'),
    hybrid_configs=dict(hybrid_none=asdict(FullPriorConfig(charge_strength=0.,class_strength=0.)),
                        hybrid_both=asdict(FullPriorConfig(charge_strength=10.,class_strength=10.)),
                        class_only=asdict(FullPriorConfig())),
    physical_baseline='Frozen expanded fit_model_selection with addition_cost=1',
    common_methods=['baseline','unet','hybrid_none','hybrid_both','class_only'],
    test_bootstrap_seed=261091102,test_plan=plans,test_labels_seen=False,test_roots_generated=False,
    candidate_parameters_changed_after_development=False,
    test_rule='Apply the same practical screen and report all strata/paired-window uncertainty without retuning')
save(folder/'frozen-fit-settings.json',settings)
manifest=dict(created=datetime.now().astimezone().isoformat(),git_head=protocol['git_head'],
    purpose='New native independent test of separated explicit neural penalties',window_us=20.,gate_us=[5.,15.],
    max_workers=4,parents=old['parents'],physical_parent=str(parent.relative_to(project)),physics_changed=False,
    planned_runs=plans,runs=[],preflight_validated=False,input_policy='Full analog voltage; adequate dynamic range',
    frozen_fit_settings_sha256=sha(folder/'frozen-fit-settings.json'),source_sha256=native_files,
    git_status_at_development_freeze=protocol['git_status_before'])
save(folder/'manifest.json',manifest)
native='''from pathlib import Path
import argparse
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'));sys.path.insert(0,str(folder/'dependency-source'))
from pileup.analysis import alpha_gamma_regularized_campaign as campaign
p=argparse.ArgumentParser()
p.add_argument('--stage',choices=['preflight','test_mixed'],required=True)
p.add_argument('--inputs',type=Path,required=True)
p.add_argument('--workers',type=int,default=4)
a=p.parse_args();campaign.PROJECT=folder.parents[2]
campaign.run(folder,a.stage,a.inputs,a.workers)
'''
(folder/'run_native.py').write_text(native)
runner='''from pathlib import Path
import runpy
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'));sys.path.insert(0,str(folder/'dependency-source'))
runpy.run_path(str(folder/'test-source/alpha_gamma_full_prior_test.py'),run_name='__main__')
'''
(folder/'run_test.py').write_text(runner)
shutil.copytree(test_source,work/'test-source')
shutil.copy2(folder/'run_test.py',work/'run_test.py')
shutil.copy2(folder/'frozen-fit-settings.json',work/'test-settings.json')
destination=work/'models/unet';destination.mkdir(parents=True,exist_ok=False)
for name in ('model.pt','training-config.json','frozen-model.json'):
    shutil.copy2(project/'jupyter/artifacts/alpha_gamma_ml_20261016/unet'/name,destination/name)
assert sha(destination/'model.pt')==settings['neural_model_sha256']
save(record/'test-preparation.json',dict(created=datetime.now().astimezone().isoformat(),
    settings_sha256=sha(folder/'frozen-fit-settings.json'),native_source_files=len(native_files),
    used_previous_seeds_checked=len(known),planned_windows=16,preflight_windows=4,
    test_model_sha256=sha(destination/'model.pt'),test_source_sha256=sha(test_source/module.name)))
print(json.dumps(dict(frozen=True,settings_sha256=sha(folder/'frozen-fit-settings.json'),windows=16)),flush=True)
