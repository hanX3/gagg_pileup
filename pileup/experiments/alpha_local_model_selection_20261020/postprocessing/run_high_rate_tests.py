"""Run the six frozen high-rate fits in independent Python processes.

One fit runs alongside native production; after G4 completion at most four fit
processes run. This changes scheduling only, not the frozen inference code.
"""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import time

PROJECT=Path.cwd();EXPERIMENT=PROJECT/'pileup/experiments/alpha_local_model_selection_20261020'
OUTPUT=PROJECT/'jupyter/artifacts/alpha_local_model_selection_20261020';RECORD=PROJECT/'jupyter/run_records/alpha_local_model_selection_20261020'
sys.path.insert(0,str(EXPERIMENT/'inference-source'))
from pileup.analysis import alpha_local_fit_study as module
module.PROJECT=PROJECT;module.EXPERIMENT=EXPERIMENT;module.OUTPUT=OUTPUT;module.RECORD=RECORD


def main():
    p=argparse.ArgumentParser();p.add_argument('--one');args=p.parse_args()
    frozen=json.loads((OUTPUT/'freeze.json').read_text())
    for name,expected in frozen['frozen_files'].items():assert module.digest(PROJECT/name)==expected,name
    config=module.LocalConfig(**frozen['selected']['config'])
    if args.one:
        data=json.loads((EXPERIMENT/'manifest.json').read_text())
        row=next(r for r in data['runs'] if r['root']==args.one)
        assert row['status']=='completed' and row['rate_hz']>=60_000_000
        module.common_fit(row,config,prefix='test');return
    environment=dict(os.environ);environment.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    active=[];seen=set();history=[]
    try:
        while len(seen)<6:
            try:data=json.loads((EXPERIMENT/'manifest.json').read_text())
            except json.JSONDecodeError:time.sleep(.5);continue
            limit=4 if len(data['runs'])==23 and all(r['status']=='completed' for r in data['runs']) else 1
            for proc,stream,item in active.copy():
                if proc.poll() is None:continue
                stream.close();item.update(returncode=proc.returncode,finished=datetime.now().astimezone().isoformat())
                assert proc.returncode==0,item
                seen.add(item['root']);active.remove((proc,stream,item))
                print('HIGH_FIT_DONE',item['root'],flush=True)
            ready=[r for r in data['runs'] if r['status']=='completed' and r['rate_hz']>=60_000_000 and r['root'] not in seen and r['root'] not in {a[2]['root'] for a in active}]
            for row in ready:
                if len(active)>=limit:break
                destination=OUTPUT/'test/selected'/(Path(row['root']).stem+'.json')
                if destination.exists():seen.add(row['root']);continue
                log=RECORD/(Path(row['root']).stem+'-high-fit.log')
                assert not log.exists(),log
                stream=log.open('x');command=[sys.executable,str(Path(__file__).resolve()),'--one',row['root']]
                proc=subprocess.Popen(command,cwd=PROJECT,env=environment,stdout=stream,stderr=subprocess.STDOUT)
                item=dict(root=row['root'],command=command,log=str(log.relative_to(PROJECT)),pid=proc.pid,
                          started=datetime.now().astimezone().isoformat(),workers_limit=limit)
                active.append((proc,stream,item));history.append(item)
                print('HIGH_FIT_START',row['root'],flush=True)
            module.save(RECORD/'high-rate-execution.json',dict(frozen_at=frozen['frozen_at'],jobs=history,
                script_sha256=module.digest(Path(__file__)),policy='1 fit while native G4 runs; max 4 independent fit processes after native completion'))
            time.sleep(.5)
    finally:
        for proc,stream,item in active:
            if proc.poll() is None:proc.terminate();proc.wait()
            stream.close()
    print('HIGH_RATE_TEST_COMPLETE',len(seen),flush=True)


if __name__=='__main__':main()
