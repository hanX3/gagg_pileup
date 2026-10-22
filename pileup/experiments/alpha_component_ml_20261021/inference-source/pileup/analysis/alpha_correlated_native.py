"""Independent native test windows after correlated/GRU inference freeze."""
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

from .alpha_correlated_study import PROJECT,EXPERIMENT,OUTPUT,RECORD
from .alpha_rate_campaign import digest,save
from .alpha_rate_data import export_window,open_root


def prepare():
    assert not (EXPERIMENT/'manifest.json').exists()
    parent=PROJECT/'pileup/experiments/alpha_local_capacity_20261020'
    capacity=json.loads((parent/'manifest.json').read_text());assert capacity['equivalence_passed']
    original=json.loads((PROJECT/'pileup/experiments/alpha_local_model_selection_20261020/manifest.json').read_text())
    shutil.copytree(parent/'source/g4',EXPERIMENT/'source/g4')
    for sub in ('bin','macros','logs'):(EXPERIMENT/sub).mkdir()
    shutil.copyfile(parent/'bin/gagg',EXPERIMENT/'bin/gagg')
    assert digest(EXPERIMENT/'bin/gagg')==capacity['binary_sha256']
    plans=[]
    for rate,n in [(100000,8),(1000000,3),(10000000,2),(20000000,2),(40000000,2),(60000000,2),(80000000,2),(100000000,2)]:
        for i in range(n):plans.append(dict(label=f'alpha_{rate:09d}_{i:02}',rate_hz=rate,seed=262026000+len(plans),role='new_test',preflight=rate in (1000000,10000000) and i==0))
    seeds={p['seed'] for p in plans}
    for path in (PROJECT/'pileup/experiments').glob('*/manifest.json'):
        data=json.loads(path.read_text())
        for row in data.get('runs',[])+data.get('planned_runs',[]):assert row.get('seed') not in seeds,path
    save(EXPERIMENT/'manifest.json',dict(created=datetime.now().astimezone().isoformat(),planned_runs=plans,runs=[],
        preflight_validated=False,max_workers=4,binary_sha256=capacity['binary_sha256'],environment=original['environment'],
        parent=str(parent.relative_to(PROJECT)),source_sha256={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'source').rglob('*') if p.is_file()},
        source_model='Unchanged p-11B effective single-alpha spectrum, no gamma',
        capacity_guard='Previously validated 30M candidate implementation; physical/readout parameters unchanged',
        window_us=20.,gate_us=[5.,15.],adequate_dynamic_range=True,root_rule='One serial G4 event per ROOT; no overwrite',
        source_rate_definition='Injected alpha/s, not detector interaction or reaction rate'))
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


def export(row,cache):
    root=PROJECT/row['root'];local=cache/root.name
    assert not local.exists()
    shutil.copyfile(root,local)
    try:
        original_hash=digest(root);assert digest(local)==original_hash
        try:metadata=export_window(local,OUTPUT/'inputs',expected=row)
        except Exception as error:
            if 'decompress' not in str(error).lower():raise
            # Exact redundant-field adapter was independently validated in the
            # preceding experiment. Never guess ownership/time or alter voltage.
            from .root_primary_time_recovery import scoring_time_recovery
            checks=dict(original_error=repr(error))
            with scoring_time_recovery(export_window,checks):metadata=export_window(local,OUTPUT/'inputs',expected=row)
            save(RECORD/(root.stem+'-redundant-time.json'),checks)
            row['scoring_time_recovery']=str((RECORD/(root.stem+'-redundant-time.json')).relative_to(PROJECT))
        assert metadata['root_sha256']==original_hash
        metadata['root']=str(root.relative_to(PROJECT));metadata['io_note']='Verified regular byte copy to local temporary file; no shared mmap or link'
        save(OUTPUT/'inputs/metadata'/f'{root.stem}.json',metadata)
        with open_root(local) as f:
            info=ak.to_list(f['RunInfo'].arrays())[0]
            assert all(f[k].num_entries==1 for k in ('RunInfo','WaveformEvent','SiPMEvent'))
            assert info['random_seed']==row['seed'] and info['alpha_p11b_rate_hz']==row['rate_hz']
            assert not info['brems_enabled'] and not info['c12_capture_enabled']
        with np.load(OUTPUT/'inputs/waveforms'/f'{root.stem}.npz') as z:
            assert len(z['time_ns'])==10000 and np.isfinite(z['analog_mV']).all()
        row.update(root_sha256=original_hash,root_bytes=root.stat().st_size,
            export_metadata=str((OUTPUT/'inputs/metadata'/f'{root.stem}.json').relative_to(PROJECT)))
    finally:
        # Only this function's verified temporary copy; original ROOT retained.
        local.unlink(missing_ok=True)


def run(workers=4):
    assert 1<=workers<=4
    manifest=json.loads((EXPERIMENT/'manifest.json').read_text());frozen=json.loads((OUTPUT/'freeze.json').read_text())
    for name,expected in frozen['frozen_files'].items():assert digest(PROJECT/name)==expected,name
    assert all(r['status']=='completed' for r in manifest['runs'])
    completed={r['label'] for r in manifest['runs']};plans=[r for r in manifest['planned_runs'] if r['label'] not in completed]
    plans.sort(key=lambda r:(not r['preflight'],r['seed']))
    cache=Path('/tmp/gagg-correlated-native-20261020')/manifest['binary_sha256'];cache.mkdir(parents=True,exist_ok=True)
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
                    'started':datetime.now().astimezone().isoformat(),'sipm_seed':p['seed']+100000,'readout_seed':p['seed']+200000}
                mac=EXPERIMENT/'macros'/(p['label']+'.mac');log=EXPERIMENT/'logs'/(p['label']+'.log');local_log=cache/(root.stem+'.log')
                assert not mac.exists() and not log.exists() and not local_log.exists()
                mac.write_text(macro(row,root));stream=local_log.open('x')
                proc=subprocess.Popen([str(binary),str(mac)],cwd=PROJECT/'g4/build',env=env,stdout=stream,stderr=subprocess.STDOUT)
                row.update(macro=str(mac.relative_to(PROJECT)),macro_sha256=digest(mac),log=str(log.relative_to(PROJECT)),
                    actual_stdout_file=str(local_log),command=[str(binary),str(mac)],cwd='g4/build',pid=proc.pid,
                    binary_sha256=manifest['binary_sha256'],inference_frozen_at=frozen['frozen_at'])
                manifest['runs'].append(row);active.append((proc,stream,row,time.monotonic()));last=time.monotonic();persist()
                print('NATIVE_START',row['label'],row['root'],flush=True)
            for item in active.copy():
                proc,stream,row,tick=item
                if proc.poll() is None:continue
                stream.close();shutil.copyfile(row['actual_stdout_file'],PROJECT/row['log'])
                row.update(returncode=proc.returncode,elapsed_s=time.monotonic()-tick,status='completed' if proc.returncode==0 else 'failed')
                if proc.returncode or 'ROOT file is saved' not in (PROJECT/row['log']).read_text():
                    row['status']='failed';persist();raise RuntimeError(row)
                export(row,cache);active.remove(item);persist();print('NATIVE_DONE',row['label'],round(row['elapsed_s'],1),flush=True)
            if not manifest['preflight_validated'] and sum(r['status']=='completed' and r['preflight'] for r in manifest['runs'])==2:
                manifest['preflight_validated']=True;manifest['preflight_passed_at']=datetime.now().astimezone().isoformat();persist();print('PREFLIGHT_PASSED',flush=True)
            time.sleep(.15)
    finally:
        for proc,stream,row,_ in active:
            if proc.poll() is None:proc.terminate();proc.wait()
            stream.close();row['status']='interrupted'
        persist()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','run']);parser.add_argument('--workers',type=int,default=4);a=parser.parse_args()
    prepare() if a.stage=='prepare' else run(a.workers)
