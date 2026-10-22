"""Complete predeclared data/evaluation/report stages, stop on any error."""
import os
from pathlib import Path
import subprocess
import sys
import time

project=Path(__file__).resolve().parents[3]
record=project/'jupyter/run_records/alpha_correlated_global_20261020'
environment=dict(os.environ)
environment.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
    PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR='/tmp/gagg-correlated-mpl')
while not (record/'remaining-production.json').exists() or not (record/'recovered-export.json').exists():
    time.sleep(2.)

stages=[
    ('finish_native_exports.py',False),
    ('evaluate_after_native.py',False),
    ('-m pileup.analysis.alpha_correlated_evaluate neural',True),
    ('-m pileup.analysis.alpha_correlated_evaluate shifts',False),
    ('-m pileup.analysis.alpha_correlated_report',False),
    ('validate_results.py',False),
    ('write_notebooks.py',False),
    ('validate_notebooks.py',False),
    ('write_docs.py',False),
]
for stage,torch_required in stages:
    arguments=stage.split() if stage.startswith('-m ') else [str(record/stage)]
    env=dict(environment)
    if torch_required:env['PYTHONPATH']='/tmp/gagg-correlated-ml-deps-20261020'
    log=record/('pipeline-'+arguments[-1].replace('/','_')+'.log')
    print('PIPELINE_STAGE_START',stage,flush=True)
    with log.open('x') as stream:
        completed=subprocess.run([sys.executable,*arguments],cwd=project,env=env,stdout=stream,stderr=subprocess.STDOUT)
    print('PIPELINE_STAGE_DONE',stage,completed.returncode,flush=True)
    assert completed.returncode==0,f'{stage} failed; see {log}'
print('COMPLETE_PREDECLARED_STUDY',flush=True)
