"""Replay the frozen full-voltage analysis from the project root."""
from pathlib import Path
import hashlib
import json
import sys

FOLDER = Path(__file__).resolve().parent
PROJECT = FOLDER.parents[2]
manifest = json.loads((FOLDER / 'manifest.json').read_text())
for relative, expected in manifest['source_sha256'].items():
    with (FOLDER / relative).open('rb') as stream:
        assert hashlib.file_digest(stream, 'sha256').hexdigest() == expected, relative
sys.path.insert(0, str(FOLDER / 'source'))
sys.path.append(str(PROJECT))
from pileup.analysis.alpha_full_signal_study import main
if __name__ == '__main__':
    if '--inputs' not in sys.argv:
        sys.argv.extend(['--inputs', *manifest['input_folders']])
    main()
