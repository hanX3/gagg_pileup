from pathlib import Path
import json
import shutil
import sys

record=Path(__file__).resolve().parent
source=Path('/tmp/gagg-full-prior-20261019')
work=Path('/tmp/gagg-full-prior-audit-20261019')
work.mkdir(exist_ok=False)
stem=sorted((source/'cases').glob('*.json'))[0].stem
shutil.copytree(source/'inputs'/stem,work/'inputs'/stem)
(work/'cases').mkdir();shutil.copy2(source/'cases'/f'{stem}.json',work/'cases'/f'{stem}.json')
shutil.copy2(source/'calibration.json',work/'calibration.json')
access=[]

def audit(event,args):
    if event=='open' and isinstance(args[0],(str,bytes)):
        path=str(args[0])
        if path.endswith('.root') or '/truth/' in path:
            raise RuntimeError('Truth/ROOT access forbidden during blind reconstruction: '+path)
        access.append(path)

sys.addaudithook(audit)
sys.path.insert(0,str(source/'source'));sys.path.insert(0,str(source/'dependency-source'))
from pileup.analysis.alpha_gamma_full_prior_study import fit_case
fit_case(work,stem)
value=dict(stem=stem,work=str(work),root_or_truth_access_blocked=True,
           fit_completed=True,opened_files=sorted(set(access)),number_of_file_opens=len(access))
(record/'fit-file-access.json').write_text(json.dumps(value,indent=2)+'\n')
print(json.dumps(dict(stem=stem,opens=len(access),blind_fit_completed=True)),flush=True)
