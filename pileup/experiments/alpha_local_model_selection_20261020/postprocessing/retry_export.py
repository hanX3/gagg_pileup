"""Retry the same closed ROOT bytes with the unchanged frozen exporter."""
import importlib.util
import json
from pathlib import Path
import shutil
import traceback

P=Path.cwd();R=P/'jupyter/run_records/alpha_local_model_selection_20261020'
E=P/'pileup/experiments/alpha_local_model_selection_20261020'
O=P/'jupyter/artifacts/alpha_local_model_selection_20261020'
spec=importlib.util.spec_from_file_location('cached_native',R/'resume_with_read_cache.py')
io=importlib.util.module_from_spec(spec);spec.loader.exec_module(io)
m=json.loads((E/'manifest.json').read_text())
row=next(r for r in m['runs'] if r['label']=='alpha_100000000_01')
stem=Path(row['root']).stem;local=Path('/tmp/gagg-alpha-local-read-cache')/(stem+'.root')
sidecar=json.loads(local.with_name(stem+'-read-cache.json').read_text())
assert io.digest(local)==sidecar['sha256']
folder=O/'export-retry-diagnostic';assert not folder.exists()
try:
    meta=io.export_window(local,folder,expected={**row,'root_sha256':sidecar['sha256']})
except Exception as error:
    io.save(R/'export-retry-result.json',dict(passed=False,root=row['root'],
        local_copy_sha256=sidecar['sha256'],exception=repr(error),traceback=traceback.format_exc(),
        source='Unchanged frozen exporter; independent retry in a fresh Python process'))
    raise
else:
    for kind in ('waveforms','truth'):
        dest=O/'inputs'/kind/(stem+'.npz')
        if dest.exists():
            old=R/'export-recovery-preexisting'/kind/dest.name;old.parent.mkdir(parents=True,exist_ok=True)
            assert not old.exists();shutil.copyfile(dest,old)
        shutil.copyfile(folder/kind/(stem+'.npz'),dest)
    meta.update(root=str(P/row['root']),read_cache=str(local),read_cache_source_sha256=sidecar['sha256'])
    io.save(O/'inputs/metadata'/(stem+'.json'),meta)
    row.update(status='completed',root_sha256=meta['root_sha256'],root_bytes=(P/row['root']).stat().st_size,
        export_metadata=str((O/'inputs/metadata'/(stem+'.json')).relative_to(P)))
    io.save(E/'manifest.json',m)
    capacity=P/m['capacity_extension']['experiment']/'manifest.json';cm=json.loads(capacity.read_text())
    cm['runs']=[r for r in m['runs'] if r['rate_hz']==100_000_000];io.save(capacity,cm)
    io.save(R/'export-retry-result.json',dict(passed=True,root=row['root'],
        local_copy_sha256=sidecar['sha256'],source='Unchanged frozen exporter; independent retry in a fresh Python process'))
    print('ROOT_EXPORT_RETRY_PASSED',row['root'],flush=True)
