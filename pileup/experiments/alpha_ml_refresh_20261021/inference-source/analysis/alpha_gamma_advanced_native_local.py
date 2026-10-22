"""Frozen-test native runner with ordinary local executable/log caches.

Model/config/seeds and ROOT destinations are unchanged. Local caches avoid
thousands of shared-folder writes from Geant4 initialization diagnostics.
"""
from datetime import datetime
import os
from pathlib import Path
import shutil
import subprocess
import time
import json

from pileup.analysis import alpha_gamma_regularized_campaign as native
from .alpha_gamma_short_window_study import digest, save

PROJECT = Path(__file__).resolve().parents[1]


def run(folder, workers=4):
    folder = Path(folder).resolve()
    manifest = json.loads((folder/'manifest.json').read_text())
    output = PROJECT/manifest['analysis_output']
    frozen = json.loads((output/'freeze.json').read_text())
    assert 1 <= workers <= min(8, manifest['max_workers']) and manifest['preflight_validated']
    for relative, expected in frozen['frozen_files'].items():
        assert digest(PROJECT/relative) == expected, relative
    assert all(r['status']=='completed' for r in manifest['runs'])
    complete = {r['label'] for r in manifest['runs']}
    plans = [r for r in manifest['planned_runs'] if r['role']=='test_single' and r['label'] not in complete]
    expected = manifest['parents']['brems']['binary_sha256']
    cache = Path('/tmp/gagg-advanced-native-cache-20261020')/expected
    cache.mkdir(parents=True, exist_ok=True)
    binary = cache/'gagg-brems'
    if not binary.exists():
        shutil.copy2(folder/'bin/gagg-brems', binary)
    assert digest(binary) == digest(folder/'bin/gagg-brems') == expected
    manifest['test_runner'] = {'source_sha256': digest(__file__),
        'source': 'inference-source/analysis/alpha_gamma_advanced_native_local.py',
        'local_binary_cache': str(binary), 'binary_sha256': expected,
        'stdout_policy': 'Regular local log during run, byte copy to retained experiment log after exit',
        'no_symbolic_links': True, 'workers': workers}
    save(folder/'frozen-fit-settings.json', frozen)
    env = dict(os.environ)
    env.update(manifest['parents']['brems']['environment']['g4_dataset_environment'])
    native.PROJECT = PROJECT
    active = []; last = -float('inf')
    def persist():
        save(folder/'manifest.json', manifest)
    try:
        while plans or active:
            if plans and len(active)<workers and time.monotonic()-last >= 1.05:
                root = PROJECT/'g4/data'/datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss.root')
                if root.exists():
                    time.sleep(.2); continue
                planned = plans.pop(0)
                row = {**planned, 'root': str(root.relative_to(PROJECT)), 'status':'running',
                       'sipm_seed': planned['seed']+100000, 'readout_seed': planned['seed']+200000,
                       'started': datetime.now().astimezone().isoformat()}
                macro = folder/'macros'/(row['label']+'.mac')
                log = folder/'logs'/(row['label']+'.log')
                local_log = cache/(root.stem+'.log')
                assert not macro.exists() and not log.exists() and not local_log.exists()
                text = native.macro(folder,row,root).replace('/gun/position 0 0 0 mm',
                       f"/gun/position {row['x_mm']} {row['y_mm']} 0 mm")
                macro.write_text(text)
                stream = local_log.open('x')
                proc = subprocess.Popen([str(binary),str(macro)],cwd=PROJECT/'g4/build',env=env,
                                        stdout=stream,stderr=subprocess.STDOUT)
                row.update(macro=str(macro.relative_to(PROJECT)), macro_sha256=digest(macro),
                           log=str(log.relative_to(PROJECT)), actual_stdout_file=str(local_log),
                           pid=proc.pid, command=[str(binary),str(macro)], cwd='g4/build',
                           binary_sha256=expected)
                tick = time.monotonic(); last = tick
                manifest['runs'].append(row); active.append((proc,stream,row,tick)); persist()
                print('START_LOCAL_LOG',row['label'],row['root'],flush=True)
            for task in list(active):
                proc,stream,row,tick = task
                if proc.poll() is None:
                    continue
                stream.close(); active.remove(task)
                local_log = Path(row['actual_stdout_file'])
                shutil.copy2(local_log,PROJECT/row['log'])
                text = local_log.read_text(); root = PROJECT/row['root']
                ok = proc.returncode==0 and root.exists() and 'ROOT file is saved' in text
                ok = ok and not any(s in text for s in ('COMMAND NOT FOUND','Error in <TFile','Error in <TTree','Error in <TBranch'))
                row.update(returncode=proc.returncode,elapsed_s=time.monotonic()-tick,
                           status='completed' if ok else 'failed',finished=datetime.now().astimezone().isoformat())
                if root.exists():
                    row.update(root_sha256=digest(root),root_bytes=root.stat().st_size)
                persist(); assert ok,row['log']
                meta = native.export_window(root,output/'native-inputs',allowed_particles=('alpha','gamma'))
                assert meta['seed']==row['seed'] and meta['macro_contents']==(PROJECT/row['macro']).read_text()
                row['validation'] = {'macro_verified':True,'seeds_verified':True,'charge_conserved':True,
                                     'primaries':meta['primaries'],'charge_pe':meta['charge_pe']}
                persist(); print('DONE_LOCAL_LOG',row['label'],round(row['elapsed_s'],2),'s',flush=True)
            if plans or active:
                time.sleep(.2)
    except BaseException:
        for proc,stream,row,tick in active:
            proc.terminate();proc.wait();stream.close()
            shutil.copy2(row['actual_stdout_file'],PROJECT/row['log'])
            row.update(status='interrupted',returncode=proc.returncode)
        persist();raise
    from . import alpha_gamma_advanced_native as campaign
    campaign.PROJECT = PROJECT
    campaign.run(folder,'test_single',workers=workers)
