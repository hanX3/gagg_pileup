"""Independent native calibration, validation and held-out mixed windows.

Run from the repository root. Each file is one serial G4 event. Frozen source,
configuration, executable identity and macros are retained in the experiment.
"""
from __future__ import annotations
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from .alpha_rate_campaign import digest, save
from .alpha_rate_data import export_window, open_root

PROJECT = Path(__file__).resolve().parents[2]


def prepare(folder):
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=False)
    for part in ('macros', 'logs', 'bin', 'source/pileup/analysis'):
        (folder/part).mkdir(parents=True, exist_ok=True)
    parents = {}
    for model in ('brems', 'mono662'):
        old = PROJECT/f'pileup/experiments/alpha_gamma_{model}_fit_20261015'
        record = json.loads((old/'manifest.json').read_text())
        assert digest(old/'build/gagg') == record['binary_sha256']
        shutil.copytree(old/'source/g4', folder/f'source/{model}/g4')
        shutil.copy2(old/'build/gagg', folder/f'bin/gagg-{model}')
        parents[model] = dict(experiment=str(old.relative_to(PROJECT)),
                             binary_sha256=record['binary_sha256'], environment=record['environment'])
    shutil.copy2(__file__, folder/'source/pileup/analysis'/Path(__file__).name)
    shutil.copy2(Path(__file__).with_name('alpha_rate_campaign.py'),
                 folder/'source/pileup/analysis/alpha_rate_campaign.py')
    plans = []
    for role, repetitions in (('calibration', {'alpha':4,'gamma':8}),
                              ('validation_single', {'alpha':2,'gamma':2})):
        for particle, energies in (('alpha',(2.,3.5,5.486,8.)),
                                    ('gamma',(.05,.1,.2,.661657,1.5))):
            for energy in energies:
                for i in range(repetitions[particle]):
                    plans.append(dict(label=f'{role}_{particle}_{energy:g}_{i:02}',
                        role=role, particle=particle, energy_MeV=energy, gamma_model='brems',
                        rate_hz=0, gamma_rate_hz=0, seed=261030000+len(plans)))
    # Four cheap coverage checks must finish before the rest of the campaign.
    selected = [next(p for p in plans if p['role']=='calibration' and
                     p['particle']==particle and p['energy_MeV']==energy)
                for particle,energy in (('alpha',2.),('alpha',8.),('gamma',.1),('gamma',.661657))]
    for p in selected:p['preflight']=True
    for role, repetitions, base in (('validation_mixed',1,261033000),('test_mixed',2,261034000)):
        for model in ('brems','mono662'):
            for rate in (1000000,10000000):
                for i in range(repetitions):
                    plans.append(dict(label=f'{role}_{model}_{rate}_{i:02}',role=role,
                        gamma_model=model,rate_hz=rate,gamma_rate_hz=10000000,
                        seed=base+len([p for p in plans if p['role']==role])))
    save(folder/'manifest.json',dict(created=datetime.now().astimezone().isoformat(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        purpose='Independent native calibration/validation and frozen-parameter held-out mixed test',
        window_us=20.,gate_us=[5.,15.],max_workers=4,parents=parents,planned_runs=plans,runs=[],
        preflight_validated=False,input_policy='Full analog voltage; adequate dynamic range',
        source_sha256={str(p.relative_to(folder)):digest(p) for p in (folder/'source').rglob('*') if p.is_file()}))
    (folder/'git-status-before.txt').write_text(subprocess.check_output(
        ['git','-c','core.fileMode=false','status','--short'],cwd=PROJECT,text=True))
    (folder/'working-changes.patch').write_bytes(subprocess.check_output(
        ['git','-c','core.fileMode=false','diff','--binary'],cwd=PROJECT))
    (folder/'run_campaign.py').write_text('''from pathlib import Path
import importlib.util, sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder.parents[2]))
spec=importlib.util.spec_from_file_location('pileup.analysis.frozen_regularized_campaign',folder/'source/pileup/analysis/alpha_gamma_regularized_campaign.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.PROJECT=folder.parents[2]
module.main()
''')
    return folder


def macro(folder, row, output):
    cfg=folder/f"source/{row['gamma_model']}/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json"
    lines=f'''/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/gagg/run/windowLength 20 us
/gagg/run/preWindow 5 us
/gagg/run/postWindow 5 us
/gagg/run/seed {row['seed']}
/gagg/run/outputFile {output}
/gagg/sipm/configFile {cfg}
/gagg/sipm/seed {row['seed']+100000}
/gagg/sipm/readoutSeed {row['seed']+200000}
/gagg/sipm/enable true
/gagg/step/setGaggMaxStep 1 um
/gagg/step/setMylarFrontMaxStep 0.5 um
/gagg/step/setMylarSideMaxStep 1 um
/run/initialize
'''
    if 'particle' in row:
        lines+=f'''/gagg/pileup/enable false
/gagg/source/useDiskConeSource false
/gun/particle {row['particle']}
/gun/energy {row['energy_MeV']} MeV
/gun/position 0 0 0 mm
/gun/direction 0 0 1
/gun/time 5 us
'''
    else:
        lines+=f'''/gagg/source/useDiskConeSource true
/gagg/source/aimAtGagg true
/gagg/source/diskRadius 10 mm
/gagg/pileup/enable true
/gagg/source/brems/enable {'true' if row['gamma_model']=='brems' else 'false'}
/gagg/source/brems/rateHz {row['gamma_rate_hz'] if row['gamma_model']=='brems' else 0}
/gagg/source/c12Capture/enable {'true' if row['gamma_model']=='mono662' else 'false'}
/gagg/source/c12Capture/rateHz {row['gamma_rate_hz'] if row['gamma_model']=='mono662' else 0}
/gagg/source/alphaP11B/enable true
/gagg/source/alphaP11B/rateHz {row['rate_hz']}
'''
    return lines+'/run/beamOn 1\n'


def run(folder, stage, inputs, workers=4):
    folder,inputs=Path(folder).resolve(),Path(inputs).resolve()
    record=json.loads((folder/'manifest.json').read_text())
    assert 1<=workers<=record['max_workers']
    if stage!='preflight':assert record['preflight_validated']
    if stage=='test_mixed':assert (folder/'frozen-fit-settings.json').exists(), 'Freeze fit before test generation'
    completed={r['label'] for r in record['runs'] if r['status']=='completed'}
    assert all(r['status']=='completed' for r in record['runs']), 'Inspect prior failed/interrupted runs'
    plans=[p for p in record['planned_runs'] if p['label'] not in completed and
           (p.get('preflight',False) if stage=='preflight' else p['role']==stage)]
    env=dict(os.environ);env.update(record['parents']['brems']['environment']['g4_dataset_environment'])
    active=[];last=-float('inf')
    def persist():save(folder/'manifest.json',record)
    try:
        while plans or active:
            if plans and len(active)<workers and time.monotonic()-last>=1.05:
                output=PROJECT/'g4/data'/datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss.root')
                if output.exists():time.sleep(.2);continue
                p=plans.pop(0);row=dict(p,root=str(output.relative_to(PROJECT)),status='running',
                    sipm_seed=p['seed']+100000,readout_seed=p['seed']+200000,
                    started=datetime.now().astimezone().isoformat())
                path=folder/'macros'/(p['label']+'.mac');log=folder/'logs'/(p['label']+'.log')
                assert not path.exists() and not log.exists()
                path.write_text(macro(folder,row,output));stream=log.open('w')
                binary=folder/f"bin/gagg-{p['gamma_model']}"
                assert digest(binary)==record['parents'][p['gamma_model']]['binary_sha256']
                start=time.monotonic();proc=subprocess.Popen([str(binary),str(path)],cwd=PROJECT/'g4/build',
                    env=env,stdout=stream,stderr=subprocess.STDOUT)
                row.update(macro=str(path.relative_to(PROJECT)),macro_sha256=digest(path),
                    log=str(log.relative_to(PROJECT)),pid=proc.pid,command=[str(binary),str(path)],cwd='g4/build')
                record['runs'].append(row);active.append((proc,stream,row,start));last=start;persist()
                print('START',row['label'],row['root'],flush=True)
            for task in list(active):
                proc,stream,row,start=task
                if proc.poll() is None:continue
                stream.close();active.remove(task);log=(PROJECT/row['log']).read_text()
                root=PROJECT/row['root'];ok=proc.returncode==0 and root.exists() and 'ROOT file is saved' in log
                ok=ok and not any(s in log for s in ('COMMAND NOT FOUND','Error in <TFile','Error in <TTree','Error in <TBranch'))
                row.update(returncode=proc.returncode,elapsed_s=time.monotonic()-start,
                    status='completed' if ok else 'failed',finished=datetime.now().astimezone().isoformat())
                if root.exists():row.update(root_sha256=digest(root),root_bytes=root.stat().st_size)
                persist()
                if not ok:raise RuntimeError(row['log'])
                meta=export_window(root,inputs,allowed_particles=('alpha','gamma'))
                assert meta['seed']==row['seed'] and meta['macro_contents']==(PROJECT/row['macro']).read_text()
                with open_root(root) as f:
                    import awkward as ak
                    info=ak.to_list(f['RunInfo'].arrays())[0]
                    assert info['sipm_avalanche_base_seed']==row['sipm_seed']
                    assert info['sipm_readout_base_seed']==row['readout_seed']
                row['validation']=dict(macro_verified=True,seeds_verified=True,charge_conserved=True,
                    primaries=meta['primaries'],charge_pe=meta['charge_pe'])
                persist();print('DONE',row['label'],round(row['elapsed_s'],1),'s',flush=True)
            if plans or active:time.sleep(.3)
        if stage=='preflight':record['preflight_validated']=True;persist()
    except BaseException:
        for proc,stream,row,start in active:
            proc.terminate();proc.wait();stream.close();row.update(status='interrupted',returncode=proc.returncode)
        persist();raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder',type=Path,required=True);p.add_argument('--inputs',type=Path)
    p.add_argument('--prepare',action='store_true')
    p.add_argument('--stage',choices=('preflight','calibration','validation_single','validation_mixed','test_mixed'))
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args()
    if a.prepare:prepare(a.folder)
    elif a.stage and a.inputs:run(a.folder,a.stage,a.inputs,a.workers)
    else:p.error('Select prepare or stage with inputs')


if __name__=='__main__':main()
