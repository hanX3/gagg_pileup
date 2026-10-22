"""Replay one saved native macro to a fresh ROOT without changing its physics."""
import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

folder = Path(__file__).resolve().parent
project = folder.parents[2]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--label', required=True)
p.add_argument('--prepare-only', action='store_true')
a = p.parse_args()
manifest = json.loads((folder/'manifest.json').read_text())
row = next(r for r in manifest['runs'] if r['label'] == a.label and r['status'] == 'completed')
digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
binary = folder/'bin/gagg-brems'
assert digest(binary) == manifest['parents']['brems']['binary_sha256']
original = project/row['macro']
assert digest(original) == row['macro_sha256']
text = original.read_text()
# Rebase the original absolute project paths when the retained project moves.
config = re.search(r'^/gagg/sipm/configFile (.+)$', text, re.M).group(1)
old_project = config.split('/pileup/experiments/')[0]
assert old_project != config
text = text.replace(old_project+'/', str(project)+'/')
assert Path(re.search(r'^/gagg/sipm/configFile (.+)$', text, re.M).group(1)).is_file()
while True:
    stamp = datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')
    root = project/'g4/data'/(stamp+'.root')
    record = project/'jupyter/run_records'/('advanced_native_replay_'+stamp)
    if not root.exists() and not record.exists():
        record.mkdir()
        break
    time.sleep(1.05)
text, replaced = re.subn(r'^/gagg/run/outputFile .+$', '/gagg/run/outputFile '+str(root), text, flags=re.M)
assert replaced == 1 and not root.exists()
macro = record/'actual.mac'; macro.write_text(text)
command = [str(binary), str(macro)]
replay = {'original_label': a.label, 'original_root': row['root'], 'original_macro_sha256': row['macro_sha256'],
          'output_root': str(root.relative_to(project)), 'actual_macro': str(macro.relative_to(project)),
          'actual_macro_sha256': digest(macro), 'seed': row['seed'], 'command': command,
          'binary_sha256': digest(binary), 'same_seed_replay_not_independent_test': True,
          'status': 'prepared_only' if a.prepare_only else 'running'}
destination = record/'replay.json'
destination.write_text(json.dumps(replay, indent=2)+'\n')
if not a.prepare_only:
    env = dict(os.environ)
    env.update(manifest['parents']['brems']['environment']['g4_dataset_environment'])
    with (record/'native.log').open('w') as log:
        result = subprocess.run(command, cwd=project/'g4/build', env=env, stdout=log, stderr=subprocess.STDOUT)
    replay.update(returncode=result.returncode, status='completed' if result.returncode == 0 and root.exists() else 'failed')
    if root.exists():
        replay['root_sha256'] = digest(root)
    destination.write_text(json.dumps(replay, indent=2)+'\n')
    assert replay['status'] == 'completed'
print(json.dumps(replay, indent=2))
