"""Replay one prefix using ONLY an experiment's saved inference sources."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--experiment',type=Path,required=True);p.add_argument('--root',type=Path,required=True)
p.add_argument('--pre',type=float,required=True);p.add_argument('--post',type=float,required=True)
p.add_argument('--output',type=Path,required=True);a=p.parse_args()
folder=a.experiment.resolve();project=folder.parents[2]
sys.path.insert(0,str(folder/'source'))
from analysis import alpha_gamma_short_window as engine
import numpy as np
assert Path(engine.__file__).resolve().is_relative_to(folder/'source')
manifest=json.loads((folder/'manifest.json').read_text())
output=project/manifest['analysis_output']
settings=json.loads((folder/'analysis-freeze.json').read_text())
digest=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
assert digest(Path(engine.__file__))==settings['source_sha256']['analysis/alpha_gamma_short_window.py']
assert digest(output/'calibration/calibration-arrays.npz')==settings['calibration_arrays_sha256']
assert digest(output/'calibration/calibration.json')==settings['calibration_sha256']
cal=engine.PrefixCalibration.load(output/'calibration')
with np.load(output/'native-inputs/waveforms'/f'{a.root.stem}.npz') as wave:
    result=engine.fit_prefix(wave['time_ns'],wave['analog_mV'],cal,pre_ns=a.pre,post_ns=a.post)
if a.output.exists():raise FileExistsError(a.output)
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print('FROZEN_PREFIX_REPLAY',a.root,a.pre,a.post)
