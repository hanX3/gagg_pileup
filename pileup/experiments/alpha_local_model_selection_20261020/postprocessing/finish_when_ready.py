"""Finish the frozen analysis once all planned native fits have completed."""
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


P=Path.cwd();R=P/'jupyter/run_records/alpha_local_model_selection_20261020'
E=P/'pileup/experiments/alpha_local_model_selection_20261020'
O=P/'jupyter/artifacts/alpha_local_model_selection_20261020'


def main():
    previous=None
    while True:
        try:
            manifest=json.loads((E/'manifest.json').read_text())
            files=[O/'test/selected'/(Path(row['root']).stem+'.json') for row in manifest['runs']]
            # Count directory entries while production/export is unfinished;
            # avoid rereading every saved fit JSON on the shared drive.
            present={path.name for path in (O/'test/selected').glob('*.json')}
            count=sum(path.name in present for path in files)
        except json.JSONDecodeError:
            time.sleep(15);continue
        assert all(row['status'] in ('completed','exporting','running') for row in manifest['runs'])
        state=(len(manifest['runs']),sum(row['status']=='completed' for row in manifest['runs']),count)
        if state!=previous:print('DELIVERY_WAIT total/native/scored',state,flush=True);previous=state
        if state==(23,23,23):
            try:parsed=[json.loads(path.read_text()) for path in files]
            except json.JSONDecodeError:time.sleep(15);continue
            break
        time.sleep(15)
    # Save the actual postprocessing/operational sources independently of the
    # immutable inference snapshot. No files protected by freeze are edited.
    for dest in (P/'pileup/analysis/alpha_local_fit_report.py',E/'postprocessing/alpha_local_fit_report.py'):
        shutil.copyfile(R/'build_report.py',dest)
    for name in ('resume_with_read_cache.py','profile_local_root_read.py','replay_blind_check.py',
                 'verify_native_roots.py','update_catalog.py','audit_delivery.py','finalize_docs.py',
                 'extend_candidate_capacity.py','finish_when_ready.py','run_high_rate_tests.py',
                 'record_progress.py','stream_low_rate_test.py','stream_mid_rate_test.py'):
        shutil.copyfile(R/name,E/'postprocessing'/name)
    paths=[P/'pileup/analysis/alpha_local_fit_report.py',*sorted((E/'postprocessing').glob('*.py'))]
    identity={str(path.relative_to(P)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    (R/'postprocessing-source.json').write_text(json.dumps(identity,indent=2)+'\n')
    env=dict(os.environ);env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
        PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR='/tmp/gagg-alpha-local-report')
    stages=[('aggregate',[sys.executable,str(E/'replay_fit.py'),'test','--workers','2']),
        ('report',[sys.executable,'-m','pileup.analysis.alpha_local_fit_report']),
        ('catalog',[sys.executable,str(R/'update_catalog.py')]),
        ('readmes',[sys.executable,str(R/'finalize_docs.py')]),
        ('audit',[sys.executable,str(R/'audit_delivery.py')]),
        ('diff-check',['git','-c','core.fileMode=false','diff','--check','--','README.md','docs/literature-use.md','pileup/README.md'])]
    jobs=[]
    for name,command in stages:
        log=R/(name+'-delivery.log');assert not log.exists(),log
        row=dict(stage=name,command=command,started=datetime.now().astimezone().isoformat(),log=str(log.relative_to(P)))
        jobs.append(row);print('DELIVERY_STAGE_START',name,flush=True)
        with log.open('x') as stream:
            proc=subprocess.run(command,cwd=P,env=env,stdout=stream,stderr=subprocess.STDOUT)
        row.update(returncode=proc.returncode,finished=datetime.now().astimezone().isoformat())
        (R/'finish-execution.json').write_text(json.dumps(jobs,indent=2)+'\n')
        assert proc.returncode==0,(name,log.read_text()[-5000:])
        print('DELIVERY_STAGE_PASS',name,flush=True)
    result=dict(completed_at=datetime.now().astimezone().isoformat(),passed=True,
        test_windows=23,auxiliary_windows=3,source_snapshot='freeze.json plus independent candidate-capacity extension',
        notebook='jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb',
        notebook_section=17,figures=['17-1','17-2','17-3'],stage_records=str((R/'finish-execution.json').relative_to(P)))
    (R/'completion.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('DELIVERY_ALL_CHECKS_COMPLETE',flush=True)


if __name__=='__main__':main()
