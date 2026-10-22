from pathlib import Path
from datetime import datetime
import hashlib
import json
import shutil

project=Path(__file__).resolve().parents[3];record=Path(__file__).resolve().parent
folder=project/'pileup/experiments/alpha_gamma_full_prior_20261019'
work=Path('/tmp/gagg-full-prior-20261019')
artifact=project/'jupyter/artifacts/alpha_gamma_full_prior_20261019'
assert len(json.loads((record/'independent-fit-runs.json').read_text()))==16
assert json.loads((work/'test-outputs/summary.json').read_text())['windows']==16
for source,destination in [('test-outputs','test'),('test-inputs','test-inputs'),('test-cases','test-cases'),('models','models')]:
    shutil.copytree(work/source,artifact/destination)
shutil.copytree(work/'test-logs',record/'test-logs')
actual=folder/'execution-source';actual.mkdir(exist_ok=False)
names=['prepare.py','prepare_test.py','run_fits.py','run_independent.py','recover_independent.py','check_inference.py',
       'audit_fit.py','trace_objective.py','append_catalog.py','record_analysis_interruption.py',
       'archive_test.py','finish.py','write_docs.py','write_notebook.py','validate.py']
files={}
for name in names:
    shutil.copy2(record/name,actual/name)
    files[str((actual/name).relative_to(folder))]=hashlib.sha256((actual/name).read_bytes()).hexdigest()
shutil.copy2(work/'blind_test_runner.py',actual/'blind_test_runner.py')
files['execution-source/blind_test_runner.py']=hashlib.sha256((actual/'blind_test_runner.py').read_bytes()).hexdigest()
(folder/'execution-source-manifest.json').write_text(json.dumps(dict(files=files,
    scope='Exact scripts used in this run; use portable replay wrappers for a new output location'),indent=2)+'\n')
archive={}
for path in sorted(artifact.rglob('*')):
    if path.is_file():archive[str(path.relative_to(artifact))]=hashlib.sha256(path.read_bytes()).hexdigest()
for source,destination in [('test-outputs','test'),('test-inputs','test-inputs'),('test-cases','test-cases'),('models','models')]:
    for path in (work/source).rglob('*'):
        if path.is_file():
            target=artifact/destination/path.relative_to(work/source)
            assert hashlib.sha256(path.read_bytes()).hexdigest()==hashlib.sha256(target.read_bytes()).hexdigest()
(record/'final-archive.json').write_text(json.dumps(dict(created=datetime.now().astimezone().isoformat(),
    artifact=str(artifact.relative_to(project)),files=archive,source=str(work)),indent=2)+'\n')
print('final_artifacts_archived',len(archive),'execution_sources',len(files),flush=True)
