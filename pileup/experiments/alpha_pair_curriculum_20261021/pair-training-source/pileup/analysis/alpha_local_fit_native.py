"""New native pure-alpha windows, after inference freeze; serial G4 per ROOT."""
from __future__ import annotations

import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import awkward as ak
import numpy as np

from .alpha_local_fit_study import PROJECT, EXPERIMENT, OUTPUT
from .alpha_rate_campaign import digest, save
from .alpha_rate_data import export_window, open_root


def prepare():
    assert not (EXPERIMENT/'manifest.json').exists()
    parent=PROJECT/'pileup/experiments/alpha_gamma_shared_tau_20261020'
    previous=json.loads((parent/'manifest.json').read_text())
    shutil.copytree(parent/'source/brems/g4',EXPERIMENT/'source/g4')
    (EXPERIMENT/'bin').mkdir();(EXPERIMENT/'macros').mkdir();(EXPERIMENT/'logs').mkdir()
    shutil.copyfile(parent/'bin/gagg-brems',EXPERIMENT/'bin/gagg')
    assert digest(EXPERIMENT/'bin/gagg')==previous['binary_sha256']
    plans=[]
    for rate,n in [(100000,8),(1000000,3),(10000000,2),(20000000,2),(40000000,2),(60000000,2),(80000000,2),(100000000,2)]:
        for i in range(n):
            plans.append(dict(label=f'alpha_{rate:09d}_{i:02}',rate_hz=rate,
                seed=262024000+len(plans),role='test_pileup',preflight=(rate in (1000000,10000000) and i==0)))
    seeds={r['seed'] for r in plans}
    for path in (PROJECT/'pileup/experiments').glob('*/manifest.json'):
        data=json.loads(path.read_text())
        for r in data.get('runs',[])+data.get('planned_runs',[]):
            assert r.get('seed') not in seeds,path
    record=dict(created=datetime.now().astimezone().isoformat(),planned_runs=plans,runs=[],max_workers=4,
        preflight_validated=False,parent=str(parent.relative_to(PROJECT)),
        binary_sha256=previous['binary_sha256'],environment=previous['parents']['brems']['environment'],
        source_sha256={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'source').rglob('*') if p.is_file()},
        source_model='Current p-11B effective single-alpha spectrum; no gamma background',
        rate_definition='Source-injected alpha/s; distinct from detector signal rate and fusion reaction rate',
        window_us=20.,gate_us=[5.,15.],input_policy='Full analog voltage; adequate dynamic range; no rail rejection')
    save(EXPERIMENT/'manifest.json',record)
    (EXPERIMENT/'replay_fit.py').write_text('''from pathlib import Path
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'inference-source'))
from pileup.analysis import alpha_local_fit_study as module
module.PROJECT=folder.parents[2]
module.EXPERIMENT=folder
module.OUTPUT=module.PROJECT/'jupyter/artifacts/alpha_local_model_selection_20261020'
module.RECORD=module.PROJECT/'jupyter/run_records/alpha_local_model_selection_20261020'
if __name__=='__main__':module.main()
''')
    print('NATIVE_PREPARED',len(plans),flush=True)


def macro(row,root):
    return f'''/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/gagg/run/windowLength 20 us
/gagg/run/preWindow 5 us
/gagg/run/postWindow 5 us
/gagg/run/seed {row['seed']}
/gagg/run/outputFile {root}
/gagg/sipm/configFile {EXPERIMENT/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json'}
/gagg/sipm/seed {row['seed']+100000}
/gagg/sipm/readoutSeed {row['seed']+200000}
/gagg/sipm/enable true
/gagg/step/setGaggMaxStep 1 um
/gagg/step/setMylarFrontMaxStep 0.5 um
/gagg/step/setMylarSideMaxStep 1 um
/run/initialize
/gagg/source/useDiskConeSource true
/gagg/source/aimAtGagg true
/gagg/source/diskRadius 10 mm
/gagg/pileup/enable true
/gagg/source/brems/enable false
/gagg/source/brems/rateHz 0
/gagg/source/c12Capture/enable false
/gagg/source/c12Capture/rateHz 0
/gagg/source/alphaP11B/enable true
/gagg/source/alphaP11B/rateHz {row['rate_hz']}
/run/beamOn 1
'''


def run(workers=4):
    assert 1<=workers<=4
    manifest=json.loads((EXPERIMENT/'manifest.json').read_text());frozen=json.loads((OUTPUT/'freeze.json').read_text())
    for name,expected in frozen['frozen_files'].items():assert digest(PROJECT/name)==expected,name
    assert all(r['status']=='completed' for r in manifest['runs'])
    completed={r['label'] for r in manifest['runs']}
    plans=[r for r in manifest['planned_runs'] if r['label'] not in completed]
    plans.sort(key=lambda r:(not r['preflight'],r['seed']))
    cache=Path('/tmp/gagg-alpha-local-native-20261020')/manifest['binary_sha256'];cache.mkdir(parents=True,exist_ok=True)
    binary=cache/'gagg'
    if not binary.exists():shutil.copyfile(EXPERIMENT/'bin/gagg',binary);binary.chmod(0o755)
    assert digest(binary)==manifest['binary_sha256']
    env=dict(os.environ);env.update(manifest['environment']['g4_dataset_environment'])
    active=[];last=-float('inf')
    def persist():save(EXPERIMENT/'manifest.json',manifest)
    try:
        while plans or active:
            if plans and len(active)<workers and (manifest['preflight_validated'] or plans[0]['preflight']) and time.monotonic()-last>=1.05:
                root=PROJECT/'g4/data'/datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss.root')
                if root.exists():time.sleep(.15);continue
                p=plans.pop(0);row={**p,'status':'running','root':str(root.relative_to(PROJECT)),
                    'sipm_seed':p['seed']+100000,'readout_seed':p['seed']+200000,
                    'started':datetime.now().astimezone().isoformat()}
                mac=EXPERIMENT/'macros'/(p['label']+'.mac');log=EXPERIMENT/'logs'/(p['label']+'.log');local_log=cache/(root.stem+'.log')
                assert not mac.exists() and not log.exists() and not local_log.exists()
                mac.write_text(macro(row,root));stream=local_log.open('x')
                proc=subprocess.Popen([str(binary),str(mac)],cwd=PROJECT/'g4/build',env=env,stdout=stream,stderr=subprocess.STDOUT)
                row.update(macro=str(mac.relative_to(PROJECT)),macro_sha256=digest(mac),log=str(log.relative_to(PROJECT)),
                    command=[str(binary),str(mac)],cwd='g4/build',actual_stdout_file=str(local_log),pid=proc.pid,
                    binary_sha256=manifest['binary_sha256'],inference_frozen_at=frozen['frozen_at'])
                active.append((proc,stream,row,time.monotonic()));manifest['runs'].append(row);last=time.monotonic();persist()
                print('NATIVE_START',row['label'],row['root'],flush=True)
            for entry in active.copy():
                proc,stream,row,tick=entry
                if proc.poll() is None:continue
                stream.close();shutil.copyfile(row['actual_stdout_file'],PROJECT/row['log'])
                row.update(elapsed_s=time.monotonic()-tick,returncode=proc.returncode,status='completed' if proc.returncode==0 else 'failed')
                root=PROJECT/row['root']
                if proc.returncode or not root.exists() or 'ROOT file is saved' not in (PROJECT/row['log']).read_text():
                    row['status']='failed';persist();raise RuntimeError(row)
                meta=export_window(root,OUTPUT/'inputs',expected=row)
                with open_root(root) as f:
                    assert all(f[name].num_entries==1 for name in ('RunInfo','WaveformEvent','SiPMEvent'))
                    info=ak.to_list(f['RunInfo'].arrays())[0]
                    assert info['alpha_p11b_rate_hz']==row['rate_hz'] and not info['brems_enabled'] and not info['c12_capture_enabled']
                with np.load(OUTPUT/'inputs/waveforms'/f'{root.stem}.npz') as a:
                    assert len(a['time_ns'])==10000 and np.isfinite(a['analog_mV']).all()
                row.update(root_sha256=meta['root_sha256'],root_bytes=root.stat().st_size,
                    export_metadata=str((OUTPUT/'inputs/metadata'/f'{root.stem}.json').relative_to(PROJECT)))
                active.remove(entry);persist();print('NATIVE_DONE',row['label'],round(row['elapsed_s'],2),flush=True)
            if not manifest['preflight_validated'] and sum(r['status']=='completed' and r['preflight'] for r in manifest['runs'])==2:
                manifest['preflight_validated']=True;manifest['preflight_passed_at']=datetime.now().astimezone().isoformat();persist();print('PREFLIGHT_PASSED',flush=True)
            time.sleep(.15)
    finally:
        for proc,stream,row,_ in active:
            if proc.poll() is None:proc.terminate();proc.wait()
            stream.close();row['status']='interrupted'
        persist()


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run']);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    prepare() if a.stage=='prepare' else run(a.workers)


if __name__=='__main__':main()
