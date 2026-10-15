"""Replay the frozen new neural/physical comparison from the project root."""
from pathlib import Path
import hashlib
import json
import sys

FOLDER=Path(__file__).resolve().parent
PROJECT=FOLDER.parents[2]
record=json.loads((FOLDER/'analysis-source.json').read_text())
for relative,expected in record['files'].items():
    with (FOLDER/relative).open('rb') as stream:
        assert hashlib.file_digest(stream,'sha256').hexdigest()==expected,relative
sys.path.insert(0,str(FOLDER/'source'));sys.path.append(str(PROJECT))
from pileup.analysis.alpha_ml_multimethod import main
if __name__=='__main__':
    if '--experiment' not in sys.argv:sys.argv.extend(['--experiment',str(FOLDER.relative_to(PROJECT))])
    main()
