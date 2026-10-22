"""Verify the two closed 100 MHz files without rereading prior checks."""
from datetime import datetime
import hashlib
import json
from pathlib import Path

P=Path.cwd();R=P/'jupyter/run_records/alpha_local_model_selection_20261020'
m=json.loads((P/'pileup/experiments/alpha_local_model_selection_20261020/manifest.json').read_text())
prior=R/'prechecked-roots.json';record=json.loads(prior.read_text())
backup=R/'prechecked-roots.before-additional.json';assert not backup.exists();backup.write_bytes(prior.read_bytes())
for row in m['runs']:
    if row['rate_hz']!=100_000_000:continue
    path=P/row['root'];before=path.stat();h=hashlib.sha256()
    cache=Path('/tmp/gagg-alpha-local-read-cache')/(path.stem+'-read-cache.json')
    expected=json.loads(cache.read_text())['sha256']
    with path.open('rb') as stream:
        while block:=stream.read(8*1024*1024):h.update(block)
    after=path.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
    assert h.hexdigest()==expected
    record['checked'].append(dict(root=row['root'],sha256=expected,size=after.st_size,
        mtime_ns=after.st_mtime_ns,checked_at=datetime.now().astimezone().isoformat()))
    temp=prior.with_suffix('.pending.json');temp.write_text(json.dumps(record,indent=2)+'\n');temp.replace(prior)
    print('ADDITIONAL_ROOT_SHA_OK',row['root'],flush=True)
