"""Score the 13 predeclared low-rate tests while higher-rate G4 runs finish."""
import json
from pathlib import Path
import sys
import time

PROJECT=Path.cwd()
EXPERIMENT=PROJECT/'pileup/experiments/alpha_local_model_selection_20261020'
OUTPUT=PROJECT/'jupyter/artifacts/alpha_local_model_selection_20261020'
sys.path.insert(0,str(EXPERIMENT/'inference-source'))
from pileup.analysis import alpha_local_fit_study as module
module.PROJECT=PROJECT;module.EXPERIMENT=EXPERIMENT;module.OUTPUT=OUTPUT
module.RECORD=PROJECT/'jupyter/run_records/alpha_local_model_selection_20261020'
frozen=json.loads((OUTPUT/'freeze.json').read_text())
for name,expected in frozen['frozen_files'].items():assert module.digest(PROJECT/name)==expected,name
config=module.LocalConfig(**frozen['selected']['config'])
seen=set()
while len(seen)<13:
    try:data=json.loads((EXPERIMENT/'manifest.json').read_text())
    except json.JSONDecodeError:time.sleep(.5);continue
    rows=[r for r in data['runs'] if r['status']=='completed' and r['rate_hz']<=10_000_000 and r['root'] not in seen]
    if not rows:time.sleep(.5);continue
    for row in rows:
        module.common_fit(row,config,prefix='test')
        seen.add(row['root'])
print('LOW_RATE_TEST_COMPLETE',len(seen),flush=True)
