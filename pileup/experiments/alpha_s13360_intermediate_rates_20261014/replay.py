"""Replay one saved run with frozen code/config and a new ROOT filename.

Example (from project root):
  bash tools/run_gagg_python.sh pileup/experiments/alpha_s13360_intermediate_rates_20261014/replay.py alpha_020000000_00

The same seed reproduces the same physical event and is not a new independent
sample. ROOT files are never overwritten. Rebuild the saved source if necessary.
"""
import argparse
from datetime import datetime
import importlib.util
import json
import os
from pathlib import Path

folder=Path(__file__).resolve().parent
project=folder.parents[2]
parser=argparse.ArgumentParser();parser.add_argument('label');args=parser.parse_args()
record=json.loads((folder/'manifest.json').read_text())
if any(r['status']=='running' for r in record['runs']):
    raise RuntimeError('Wait for production to finish before changing its manifest')
for key,value in record['environment']['g4_dataset_environment'].items():
    if key not in os.environ and Path(value).exists():os.environ[key]=value
row=next(r for r in record['runs'] if r['label']==args.label)
spec=importlib.util.spec_from_file_location('saved_runner',folder/'source/pileup/analysis/alpha_rate_campaign.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
runner.PROJECT=project
label='replay_'+args.label+'_'+datetime.now().strftime('%Y%m%d_%H%M%S')
runner.run_one(folder,label,row['rate_hz'],row['seed'],duration=row['window_us'],guard=row['gate_us'][0])
record=json.loads((folder/'manifest.json').read_text())
replayed=record['runs'].pop()
replayed['replay_of']=args.label
# Keep deterministic replays out of the independent comparison sample.
record.setdefault('replays',[]).append(replayed)
runner.save(folder/'manifest.json',record)
