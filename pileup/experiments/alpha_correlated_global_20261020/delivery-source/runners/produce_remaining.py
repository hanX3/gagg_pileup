"""Generate the four still-unstarted, predeclared seeds; defer scoring export.

Frozen G4 binary/source/config/macro and inference remain unchanged. Two G4
processes limit memory contention on the four-core VM. Original failed scoring
export is retained; no existing seed is rerun or ROOT overwritten.
"""
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

project=Path(__file__).resolve().parents[3];sys.path.insert(0,str(project))
from pileup.analysis.alpha_correlated_evaluate import frozen
from pileup.analysis.alpha_correlated_native import macro
from pileup.analysis.alpha_correlated_study import EXPERIMENT,OUTPUT,RECORD
from pileup.analysis.alpha_rate_campaign import digest,save

configuration=frozen();manifest=json.loads((EXPERIMENT/'manifest.json').read_text())
for row in manifest['runs']:
    assert row['returncode']==0 and 'ROOT file is saved' in (project/row['log']).read_text(),row
remaining=[r for r in manifest['planned_runs'] if r['label'] not in {x['label'] for x in manifest['runs']}]
assert len(remaining)==4
cache=Path('/tmp/gagg-correlated-native-20261020')/manifest['binary_sha256']
binary=cache/'gagg';assert digest(binary)==manifest['binary_sha256']
env=dict(os.environ);env.update(manifest['environment']['g4_dataset_environment'])
active=[];last=-float('inf')
manifest['production_resume']=dict(reason='Scoring-only ROOT branch decoding failure after 19 successful G4 runs; preserve all original ROOTs and seeds',
    inference_changed=False,physics_changed=False,workers=2,exports_deferred_until_g4_complete=True,
    runner=str(Path(__file__).relative_to(project)),runner_sha256=digest(Path(__file__)))
save(EXPERIMENT/'manifest.json',manifest)
try:
    while remaining or active:
        if remaining and len(active)<2 and time.monotonic()-last>=1.05:
            root=project/'g4/data'/datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss.root')
            if root.exists():time.sleep(.15);continue
            plan=remaining.pop(0);row={**plan,'status':'running','root':str(root.relative_to(project)),
                'started':datetime.now().astimezone().isoformat(),'sipm_seed':plan['seed']+100000,'readout_seed':plan['seed']+200000}
            mac=EXPERIMENT/'macros'/(plan['label']+'.mac');log=EXPERIMENT/'logs'/(plan['label']+'.log');local_log=cache/(root.stem+'.log')
            assert not mac.exists() and not log.exists() and not local_log.exists()
            mac.write_text(macro(row,root));stream=local_log.open('x')
            command=[str(binary),str(mac)]
            process=subprocess.Popen(command,cwd=project/'g4/build',env=env,stdout=stream,stderr=subprocess.STDOUT)
            row.update(macro=str(mac.relative_to(project)),macro_sha256=digest(mac),log=str(log.relative_to(project)),
                actual_stdout_file=str(local_log),command=command,cwd='g4/build',pid=process.pid,
                binary_sha256=manifest['binary_sha256'],inference_frozen_at=configuration['frozen_at'])
            manifest['runs'].append(row);active.append((process,stream,row,time.monotonic()));last=time.monotonic()
            save(EXPERIMENT/'manifest.json',manifest);print('REMAINING_NATIVE_START',row['label'],row['root'],flush=True)
        for item in active.copy():
            process,stream,row,tick=item
            if process.poll() is None:continue
            stream.close();shutil.copyfile(row['actual_stdout_file'],project/row['log'])
            row.update(returncode=process.returncode,elapsed_s=time.monotonic()-tick,
                g4_finished=datetime.now().astimezone().isoformat(),status='generated' if process.returncode==0 else 'failed')
            assert process.returncode==0 and 'ROOT file is saved' in (project/row['log']).read_text(),row
            row.update(root_sha256=digest(project/row['root']),root_bytes=(project/row['root']).stat().st_size)
            active.remove(item);save(EXPERIMENT/'manifest.json',manifest)
            print('REMAINING_NATIVE_GENERATED',row['label'],round(row['elapsed_s'],1),flush=True)
        time.sleep(.15)
finally:
    for process,stream,row,_ in active:
        if process.poll() is None:process.terminate();process.wait()
        stream.close();row['status']='interrupted'
    save(EXPERIMENT/'manifest.json',manifest)
save(RECORD/'remaining-production.json',dict(runs=[r for r in manifest['runs'] if r['status']=='generated'],
    frozen_at=configuration['frozen_at'],inference_changed=False,physics_changed=False))
print('REMAINING_NATIVE_GENERATION_COMPLETE',flush=True)
