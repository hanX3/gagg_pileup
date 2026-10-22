"""Verify closed production ROOT bytes while the final two runs compute."""
from datetime import datetime
import hashlib
import json
from pathlib import Path


project=Path.cwd()
exp=project/'pileup/experiments/alpha_local_model_selection_20261020'
record=project/'jupyter/run_records/alpha_local_model_selection_20261020'
manifest=json.loads((exp/'manifest.json').read_text())
rows=[r for r in manifest['runs'] if r['status']=='completed']
checks=[]
for row in rows:
    path=project/row['root'];before=path.stat();h=hashlib.sha256()
    with path.open('rb') as f:
        while block:=f.read(8*1024*1024):h.update(block)
    after=path.stat()
    assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
    assert h.hexdigest()==row['root_sha256'],row['root']
    checks.append(dict(root=row['root'],sha256=h.hexdigest(),size=after.st_size,
        mtime_ns=after.st_mtime_ns,checked_at=datetime.now().astimezone().isoformat()))
    (record/'prechecked-roots.json').write_text(json.dumps(dict(passed=True,checked=checks),indent=2)+'\n')
    print('ROOT_SHA_OK',len(checks),len(rows),row['root'],flush=True)
