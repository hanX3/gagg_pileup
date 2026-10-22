"""Finish independent frozen inference stages after physical assets are ready."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
import json
import subprocess
import sys
import time

EXPERIMENT=Path(__file__).resolve().parent
PROJECT=EXPERIMENT.parents[2]
RECORD=PROJECT/'jupyter/run_records/alpha_pair_resolution_20261022'
RUN=EXPERIMENT/'run.py'


def stage(name,arguments):
    command=[sys.executable,str(RUN),*arguments]
    info=dict(name=name,command=command,started=datetime.now().astimezone().isoformat())
    print('STAGE_START',name,flush=True)
    with (EXPERIMENT/'logs'/(name+'.log')).open('w') as out:
        result=subprocess.run(command,cwd=PROJECT,stdout=out,stderr=subprocess.STDOUT)
    info.update(returncode=result.returncode,finished=datetime.now().astimezone().isoformat())
    (RECORD/('stage-'+name+'.json')).write_text(json.dumps(info,indent=2)+'\n')
    print('STAGE_END',name,result.returncode,flush=True)
    if result.returncode:raise RuntimeError('Stage failed: '+name)


def shared():
    stage('shared-predict',['fit','predict','--dataset','shared','--workers','4'])
    stage('shared-diagnostics',['fit','diagnostics','--dataset','shared'])


def controls():
    stage('linear-stage',['linear'])
    for name in ('linear','shared-noise'):
        stage(name+'-predict',['fit','predict','--dataset',name,'--workers','2'])
        stage(name+'-diagnostics',['fit','diagnostics','--dataset',name])


def main():
    while True:
        marker=RECORD/'physical-io-adapter.json'
        try:ready=marker.exists() and json.loads(marker.read_text()).get('passed')
        except json.JSONDecodeError:ready=False
        if ready:break
        location=RECORD/'physical-local-cache.json'
        if location.exists():
            cache=Path(json.loads(location.read_text())['directory'])
            print('WAIT_PHYSICAL',len(list((cache/'cases').glob('*.npz'))),'/ 3480 checkpoints',flush=True)
        time.sleep(45)
    failures=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(shared),pool.submit(controls)]
        for future in as_completed(futures):
            try:future.result()
            except Exception as error:failures.append(repr(error))
    if failures:raise RuntimeError(failures)
    print('ALL_FROZEN_PREDICTIONS_AND_DIAGNOSTICS_COMPLETE',flush=True)


if __name__=='__main__':main()
