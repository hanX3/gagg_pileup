"""Replay low-rate diagnostics with their frozen source and scorer."""
from pathlib import Path
import hashlib,json,runpy,sys
folder=Path(__file__).resolve().parent
project=folder.parents[2]
manifest=json.loads((folder/'manifest.json').read_text())
for relative,expected in manifest['source_sha256'].items():
    with (folder/relative).open('rb') as f:
        assert hashlib.file_digest(f,'sha256').hexdigest()==expected,relative
old=project/manifest['reference_experiment']
with (old/'analysis-source.json').open('rb') as f:
    assert hashlib.file_digest(f,'sha256').hexdigest()==manifest['reference_analysis_source_sha256']
source=json.loads((old/'analysis-source.json').read_text())
for relative,expected in source['files'].items():
    with (old/relative).open('rb') as f:
        assert hashlib.file_digest(f,'sha256').hexdigest()==expected,relative
sys.path.insert(0,str(old/'source'))
sys.path.insert(0,str(folder/'source'))
sys.path.append(str(project))
for option,key in (('--reference','reference_results'),('--experiment','reference_experiment'),('--widths','width_reference')):
    if option not in sys.argv:sys.argv.extend([option,manifest[key]])
if __name__=='__main__':runpy.run_module('pileup.analysis.alpha_ml_low_rate',run_name='__main__')
