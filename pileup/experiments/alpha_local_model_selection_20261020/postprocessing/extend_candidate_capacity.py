"""Preserve aborted attempts; validate and retry with a larger guard only."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

P=Path.cwd();R=P/'jupyter/run_records/alpha_local_model_selection_20261020'
E=P/'pileup/experiments/alpha_local_model_selection_20261020'
C=P/'pileup/experiments/alpha_local_capacity_20261020'
O=P/'jupyter/artifacts/alpha_local_model_selection_20261020'
spec=importlib.util.spec_from_file_location('cached_native',R/'resume_with_read_cache.py')
io=importlib.util.module_from_spec(spec);spec.loader.exec_module(io)
digest,save=io.digest,io.save


def main():
    m=json.loads((E/'manifest.json').read_text());frozen=json.loads((O/'freeze.json').read_text())
    for name,h in frozen['frozen_files'].items():assert digest(P/name)==h,name
    assert not C.exists()
    failed=[r for r in m['runs'] if r['status']!='completed']
    assert len(failed)==2 and all(r['rate_hz']==100_000_000 for r in failed)
    assert 'max_candidates exceeded' in Path(failed[0]['actual_stdout_file']).read_text()
    save(R/'manifest-after-capacity-abort.json',m)
    C.mkdir();shutil.copytree(E/'source',C/'source');(C/'bin').mkdir();(C/'macros').mkdir();(C/'logs').mkdir()
    shutil.copyfile(E/'bin/gagg',C/'bin/gagg');assert digest(C/'bin/gagg')==m['binary_sha256']
    old_config=E/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json'
    config=C/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json'
    original=json.loads(old_config.read_text());changed=json.loads(config.read_text())
    assert original==changed and original['sensor']['max_candidates']==10_000_000
    changed['sensor']['max_candidates']=30_000_000
    config.write_text(json.dumps(changed,ensure_ascii=False,indent=2)+'\n')
    comparison=json.loads(config.read_text());comparison['sensor']['max_candidates']=10_000_000
    assert comparison==original
    extension=dict(created=datetime.now().astimezone().isoformat(),
        parent=str(E.relative_to(P)),reason='100 MHz / 20 us aborts at the 10-million processed-candidate guard',
        change=dict(field='sensor.max_candidates',old=10_000_000,new=30_000_000),
        physical_and_readout_parameters_identical=True,binary_sha256=m['binary_sha256'],
        code_evidence='source/g4/src/SiPMResponse.cc:205; guard checks count, does not change RNG/probabilities below the limit',
        source_sha256={str(p.relative_to(C)):digest(p) for p in (C/'source').rglob('*') if p.is_file()},
        runs=[],validation_runs=[])
    (C/'README.md').write_text('''# SiPM候选处理容量扩展（2026-10-20）

父实验的100 MHz、20 μs窗口触发了sensor.max_candidates=10000000保护并未保存完整响应。原失败ROOT和日志保留在父manifest.failed_attempts。此独立源码/配置快照只将该计算保护上限提高至30000000；C++和二进制身份、物理概率、微胞/读出及拟合参数均不变。先以相同1 MHz种子重放验证导出的16个数组逐项一致，再重试原计划100 MHz的两套种子。正式ROOT仍直接写入共享g4/data，普通本机副本仅用于只读导出；实际宏、日志和运行记录见manifest.json。此配置扩展不构成新物理模型，也不提高物理动态范围。
''')
    for row in failed:
        row['status']='failed_capacity_guard' if row.get('returncode')==-6 else 'interrupted_after_peer_failure'
        row['root_bytes']=(P/row['root']).stat().st_size
        row['failure_recorded_at']=datetime.now().astimezone().isoformat()
        shutil.copyfile(row['actual_stdout_file'],P/row['log'])
    m['failed_attempts']=m.get('failed_attempts',[])+failed
    m['runs']=[r for r in m['runs'] if r['status']=='completed']
    m['capacity_extension']=dict(experiment=str(C.relative_to(P)),configuration=str(config.relative_to(P)),
        configuration_sha256=digest(config),change=extension['change'],fit_parameters_unchanged=True,
        validation_record=str((R/'capacity-equivalence.json').relative_to(P)))
    save(E/'manifest.json',m);save(C/'manifest.json',extension)
    environment=dict(os.environ);environment.update(m['environment']['g4_dataset_environment'])
    binary=Path('/tmp/gagg-alpha-local-native-20261020')/m['binary_sha256']/'gagg'
    assert digest(binary)==m['binary_sha256']
    def persist():save(E/'manifest.json',m);save(C/'manifest.json',extension)
    def start(plan,role):
        root=P/'g4/data'/datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss.root')
        while root.exists():
            time.sleep(1.05);root=P/'g4/data'/datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss.root')
        row={**plan,'role':role,'status':'running','root':str(root.relative_to(P)),
            'sipm_seed':plan['seed']+100000,'readout_seed':plan['seed']+200000,
            'started':datetime.now().astimezone().isoformat(),'inference_frozen_at':frozen['frozen_at'],
            'capacity_extension':str(C.relative_to(P)),'config':str(config.relative_to(P)),
            'config_sha256':digest(config),'binary_sha256':m['binary_sha256']}
        mac=C/'macros'/(row['label']+'.mac');log=C/'logs'/(row['label']+'.log');local_log=binary.parent/(root.stem+'.log')
        assert not mac.exists() and not log.exists() and not local_log.exists()
        mac.write_text(io.native.macro(row,root).replace(str(old_config),str(config)))
        stream=local_log.open('x');command=[str(binary),str(mac)]
        process=subprocess.Popen(command,cwd=P/'g4/build',env=environment,stdout=stream,stderr=subprocess.STDOUT)
        row.update(macro=str(mac.relative_to(P)),macro_sha256=digest(mac),log=str(log.relative_to(P)),
            actual_stdout_file=str(local_log),command=command,cwd='g4/build',pid=process.pid)
        print('CAPACITY_START',row['label'],row['root'],flush=True)
        return process,stream,row,time.monotonic()
    # Reuse a completed low-rate seed solely as a paired capacity validation;
    # it is not added to the independent-test denominator.
    reference=next(r for r in m['runs'] if r['rate_hz']==1_000_000)
    process,stream,proof,tick=start(dict(label='capacity_equivalence_001MHz',rate_hz=1_000_000,seed=reference['seed']), 'capacity_validation')
    extension['validation_runs'].append(proof);persist()
    process.wait();stream.close();shutil.copyfile(proof['actual_stdout_file'],P/proof['log'])
    proof.update(returncode=process.returncode,elapsed_s=time.monotonic()-tick,status='physics_completed')
    assert process.returncode==0 and 'ROOT file is saved' in (P/proof['log']).read_text()
    local,copy=io.local_copy(P/proof['root'])
    validation_folder=O/'capacity-equivalence'
    meta=io.export_window(local,validation_folder,expected=proof)
    stem=Path(proof['root']).stem;refstem=Path(reference['root']).stem;keys=[]
    for kind in ('waveforms','truth'):
        with np.load(O/'inputs'/kind/f'{refstem}.npz') as a,np.load(validation_folder/kind/f'{stem}.npz') as b:
            assert set(a.files)==set(b.files)
            for key in a.files:
                assert a[key].dtype==b[key].dtype
                assert np.array_equal(a[key],b[key],equal_nan=a[key].dtype.kind in 'fc'),(kind,key)
                keys.append(kind+'/'+key)
    meta.update(root=str(P/proof['root']),read_cache=str(local));save(validation_folder/'metadata'/f'{stem}.json',meta)
    proof.update(status='completed',root_sha256=meta['root_sha256'],root_bytes=(P/proof['root']).stat().st_size)
    save(R/'capacity-equivalence.json',dict(passed=True,reference_root=reference['root'],
        replay_root=proof['root'],seed=reference['seed'],change=extension['change'],
        all_arrays_exactly_equal=keys,physical_configuration_identical_except_guard=True,
        binary_unchanged=True,test_cohort_entry=False))
    extension['equivalence_passed']=True;persist();print('CAPACITY_EQUIVALENCE_PASSED',len(keys),flush=True)
    plans=[next(p for p in m['planned_runs'] if p['label']==row['label']) for row in failed]
    active=[];exports={};last=-float('inf')
    with ThreadPoolExecutor(max_workers=1) as pool:
        try:
            while plans or active or exports:
                if plans and len(active)<2 and time.monotonic()-last>=1.05:
                    entry=start(plans.pop(0),'test_pileup');active.append(entry)
                    m['runs'].append(entry[2]);extension['runs'].append(entry[2]);last=time.monotonic();persist()
                for entry in active.copy():
                    proc,stream,row,tick=entry
                    if proc.poll() is None:continue
                    stream.close();shutil.copyfile(row['actual_stdout_file'],P/row['log'])
                    row.update(returncode=proc.returncode,elapsed_s=time.monotonic()-tick,
                        physics_finished=datetime.now().astimezone().isoformat())
                    assert proc.returncode==0 and 'ROOT file is saved' in (P/row['log']).read_text()
                    row['status']='exporting';exports[pool.submit(io.cached_export,dict(row))]=row
                    active.remove(entry);persist();print('CAPACITY_PHYSICS_DONE',row['label'],flush=True)
                for future in list(exports):
                    if not future.done():continue
                    row=exports.pop(future);row.update(future.result());persist()
                    print('CAPACITY_EXPORT_DONE',row['label'],flush=True)
                time.sleep(.2)
        finally:
            for proc,stream,row,_ in active:
                if proc.poll() is None:proc.terminate();proc.wait()
                stream.close();row['status']='interrupted';row['returncode']=proc.returncode
            persist()
    assert len(m['runs'])==23 and all(r['status']=='completed' for r in m['runs'])
    print('CAPACITY_RETRY_ALL_COMPLETE',flush=True)


if __name__=='__main__':main()
