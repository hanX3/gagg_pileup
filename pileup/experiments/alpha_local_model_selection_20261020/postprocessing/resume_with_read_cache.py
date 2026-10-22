"""Resume completed physics exports and two remaining runs using read caches.

Production ROOT files are still written directly to shared g4/data. Only
read-only analysis uses regular local byte copies; frozen inference is intact.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import awkward as ak
import numpy as np

PROJECT=Path.cwd();EXPERIMENT=PROJECT/'pileup/experiments/alpha_local_model_selection_20261020'
OUTPUT=PROJECT/'jupyter/artifacts/alpha_local_model_selection_20261020';RECORD=PROJECT/'jupyter/run_records/alpha_local_model_selection_20261020'
sys.path.insert(0,str(EXPERIMENT/'inference-source'))
from pileup.analysis import alpha_local_fit_study as study
study.PROJECT=PROJECT;study.EXPERIMENT=EXPERIMENT;study.OUTPUT=OUTPUT;study.RECORD=RECORD
from pileup.analysis import alpha_local_fit_native as native
from pileup.analysis.alpha_rate_data import export_window,open_root
from pileup.analysis.alpha_rate_campaign import digest,save
native.PROJECT=PROJECT;native.EXPERIMENT=EXPERIMENT;native.OUTPUT=OUTPUT


def local_copy(source):
    cache=Path('/tmp/gagg-alpha-local-read-cache');cache.mkdir(parents=True,exist_ok=True)
    target=cache/source.name;sidecar=cache/(source.stem+'-read-cache.json')
    s=source.stat();fingerprint=dict(size=s.st_size,mtime_ns=s.st_mtime_ns)
    if target.exists() and sidecar.exists():
        record=json.loads(sidecar.read_text());assert record['source_fingerprint']==fingerprint
        assert digest(target)==record['sha256'];return target,record
    if target.exists():
        preserved=cache/(target.name+'.unfinished-'+datetime.now().strftime('%H%M%S'))
        target.rename(preserved)
    tick=time.monotonic();checksum=hashlib.sha256()
    with source.open('rb') as src,target.open('xb') as dst:
        while block:=src.read(8*1024*1024):
            checksum.update(block);dst.write(block)
    after=source.stat();assert after.st_size==s.st_size and after.st_mtime_ns==s.st_mtime_ns
    value=checksum.hexdigest();assert digest(target)==value
    record=dict(source=str(source),source_fingerprint=fingerprint,local_cache=str(target),
        sha256=value,copy_seconds=time.monotonic()-tick,
        method='Canonical-file SHA256 computed during sequential byte copy, then independently check local copy')
    save(sidecar,record);return target,record


def cached_export(row):
    source=PROJECT/row['root'];local,copy=local_copy(source);stem=source.stem
    if row.get('root_sha256'):assert copy['sha256']==row['root_sha256']
    meta_path=OUTPUT/'inputs/metadata'/f'{stem}.json'
    if not meta_path.exists():
        # Preserve incomplete cache outputs made by the interrupted read-only
        # exporter. No production ROOT or useful completed result is removed.
        for sub in ('waveforms','truth'):
            old=OUTPUT/'inputs'/sub/f'{stem}.npz'
            if old.exists():
                backup=RECORD/'interrupted-export'/sub/old.name
                backup.parent.mkdir(parents=True,exist_ok=True)
                assert not backup.exists();old.rename(backup)
    tick=time.monotonic()
    meta=export_window(local,OUTPUT/'inputs',expected={**row,'root_sha256':copy['sha256']})
    export_seconds=time.monotonic()-tick
    assert meta['root_sha256']==copy['sha256']
    meta.update(root=str(source),read_cache=str(local),read_cache_source_sha256=copy['sha256'],
                cache_copy_seconds=copy['copy_seconds'],cache_export_seconds=export_seconds)
    save(meta_path,meta)
    with open_root(local) as f:
        info=ak.to_list(f['RunInfo'].arrays())[0]
        assert int(info['random_seed'])==row['seed']
        assert int(info['sipm_avalanche_base_seed'])==row['sipm_seed']
        assert int(info['sipm_readout_base_seed'])==row['readout_seed']
        assert info['alpha_p11b_rate_hz']==row['rate_hz'] and info['t_length_ps']==20_000_000
        assert not info['brems_enabled'] and not info['c12_capture_enabled']
        assert f['WaveformEvent'].num_entries==f['SiPMEvent'].num_entries==1
    with np.load(OUTPUT/'inputs/waveforms'/f'{stem}.npz') as z:
        assert len(z['time_ns'])==10000 and np.isfinite(z['analog_mV']).all()
    return dict(root_sha256=copy['sha256'],root_bytes=source.stat().st_size,read_cache=str(local),
        read_cache_record=str(local.parent/(source.stem+'-read-cache.json')),
        export_metadata=str(meta_path.relative_to(PROJECT)),cache_copy_seconds=copy['copy_seconds'],
        cache_export_seconds=export_seconds,status='completed')


def main():
    frozen=json.loads((OUTPUT/'freeze.json').read_text())
    for name,expected in frozen['frozen_files'].items():assert digest(PROJECT/name)==expected,name
    proof=json.loads((RECORD/'read-cache-validation.json').read_text());assert proof['passed']
    manifest=json.loads((EXPERIMENT/'manifest.json').read_text())
    assert len(manifest['runs'])==21
    assert all(r['status'] in ('completed','interrupted','running') for r in manifest['runs'])
    # An exited analysis process can leave its manifest at "running" after
    # physics has closed the ROOT file. Recover only the explicitly checked
    # completed-physics rows; never restart or overwrite their production.
    for row in manifest['runs']:
        if row['status']=='completed':continue
        log=Path(row['actual_stdout_file']).read_text(errors='replace')
        assert 'ROOT file is saved' in log and 'Graphics systems deleted.' in log
        assert (PROJECT/row['root']).is_file()
    shutil.copyfile(EXPERIMENT/'manifest.json',RECORD/'manifest-after-export-interrupt.json')
    existing={r['label'] for r in manifest['runs']}
    plans=[r for r in manifest['planned_runs'] if r['label'] not in existing]
    assert len(plans)==2 and all(r['rate_hz']==100_000_000 for r in plans)
    binary=Path('/tmp/gagg-alpha-local-native-20261020')/manifest['binary_sha256']/'gagg'
    assert digest(binary)==manifest['binary_sha256']
    env=dict(os.environ);env.update(manifest['environment']['g4_dataset_environment'])
    manifest['read_cache_transition']=dict(at=datetime.now().astimezone().isoformat(),
        reason='Shared-folder mmap export slow; identical frozen exporter on regular local byte copies',
        equivalence_record=str((RECORD/'read-cache-validation.json').relative_to(PROJECT)),
        original_frozen_inference_unchanged=True,production_root_policy='Still write directly to shared g4/data',
        controller=str(Path(__file__).relative_to(PROJECT)),controller_sha256=digest(Path(__file__)),
        g4_workers_after_transition=2,export_workers=1)
    active=[];exports={};last=-float('inf')
    def persist():save(EXPERIMENT/'manifest.json',manifest)
    with ThreadPoolExecutor(max_workers=1) as pool:
        for row in manifest['runs']:
            if row['status']=='completed':continue
            log=Path(row['actual_stdout_file']).read_text(errors='replace')
            assert 'ROOT file is saved' in log and 'Graphics systems deleted.' in log
            target=PROJECT/row['log']
            if not target.exists():shutil.copyfile(row['actual_stdout_file'],target)
            row['status_before_cache_transition']=row['status'];row['physics_finished_before_export_interrupt']=True
            row['status']='exporting';exports[pool.submit(cached_export,dict(row))]=row
        persist()
        try:
            while plans or active or exports:
                if plans and len(active)<2 and time.monotonic()-last>=1.05:
                    root=PROJECT/'g4/data'/datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss.root')
                    if root.exists():time.sleep(.2);continue
                    planned=plans.pop(0);row={**planned,'status':'running','root':str(root.relative_to(PROJECT)),
                        'sipm_seed':planned['seed']+100000,'readout_seed':planned['seed']+200000,
                        'started':datetime.now().astimezone().isoformat(),'inference_frozen_at':frozen['frozen_at']}
                    mac=EXPERIMENT/'macros'/(row['label']+'.mac');log=EXPERIMENT/'logs'/(row['label']+'.log')
                    local_log=binary.parent/(root.stem+'.log')
                    assert not mac.exists() and not log.exists() and not local_log.exists()
                    mac.write_text(native.macro(row,root));stream=local_log.open('x')
                    command=[str(binary),str(mac)]
                    proc=subprocess.Popen(command,cwd=PROJECT/'g4/build',env=env,stdout=stream,stderr=subprocess.STDOUT)
                    row.update(macro=str(mac.relative_to(PROJECT)),macro_sha256=digest(mac),log=str(log.relative_to(PROJECT)),
                        actual_stdout_file=str(local_log),command=command,cwd='g4/build',pid=proc.pid,binary_sha256=manifest['binary_sha256'])
                    manifest['runs'].append(row);active.append((proc,stream,row,time.monotonic()));last=time.monotonic();persist()
                    print('NATIVE_START',row['label'],row['root'],flush=True)
                for entry in active.copy():
                    proc,stream,row,tick=entry
                    if proc.poll() is None:continue
                    stream.close();shutil.copyfile(row['actual_stdout_file'],PROJECT/row['log'])
                    row.update(returncode=proc.returncode,elapsed_s=time.monotonic()-tick)
                    assert proc.returncode==0 and 'ROOT file is saved' in (PROJECT/row['log']).read_text()
                    row['status']='exporting';exports[pool.submit(cached_export,dict(row))]=row
                    active.remove(entry);persist();print('PHYSICS_DONE',row['label'],flush=True)
                for future in list(exports):
                    if not future.done():continue
                    row=exports.pop(future);row.update(future.result());persist()
                    print('CACHE_EXPORT_DONE',row['label'],round(row['cache_copy_seconds'],2),round(row['cache_export_seconds'],2),flush=True)
                time.sleep(.2)
        finally:
            for proc,stream,row,_ in active:
                if proc.poll() is None:proc.terminate();proc.wait()
                stream.close();row['status']='interrupted'
            persist()
    assert len(manifest['runs'])==23 and all(r['status']=='completed' for r in manifest['runs'])
    print('NATIVE_ALL_COMPLETE',len(manifest['runs']),flush=True)


if __name__=='__main__':main()
