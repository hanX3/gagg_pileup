"""Replay frozen inference into a NEW project-relative output, from repo root.

Working analysis modules need not match an old Git commit. Actual frozen source
and calibration/weights are used; original ROOTs and verified exports are needed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import sys

folder=Path(__file__).resolve().parent
project=folder.parents[2]
parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['one','all','controls','optimizer_controls'])
parser.add_argument('--output',type=Path,required=True);parser.add_argument('--root',help='ROOT stem for one-window physical replay')
a=parser.parse_args()
target=(project/a.output).resolve();target.relative_to(project)
assert not target.exists(),'Replay output must be new'
target.mkdir(parents=True)
assets=folder/'frozen-analysis'
for item in assets.iterdir():
    if item.is_dir():shutil.copytree(item,target/item.name)
    else:shutil.copyfile(item,target/item.name)
origin=project/'jupyter/artifacts/alpha_correlated_global_20261020'
frozen=json.loads((folder/'frozen-fit-settings.json').read_text())
mutable_prefixes=('pileup/analysis/','analysis/',str(origin.relative_to(project))+'/')
frozen['frozen_files']={k:v for k,v in frozen['frozen_files'].items() if not k.startswith(mutable_prefixes)}
frozen['replay_origin_sha256']=hashlib.sha256((folder/'frozen-fit-settings.json').read_bytes()).hexdigest()
(target/'freeze.json').write_text(json.dumps(frozen,ensure_ascii=False,indent=2)+'\n')
sys.path.insert(0,str(folder/'inference-source'))
from pileup.analysis import alpha_correlated_study as study
study.PROJECT=project;study.EXPERIMENT=folder;study.OUTPUT=target;study.OLD=target/'local-reference'
from pileup.analysis import alpha_correlated_evaluate as evaluate
evaluate.frozen()
if a.stage=='one':
    rows=json.loads((folder/'manifest.json').read_text())['runs']
    row=next(r for r in rows if Path(r['root']).stem==a.root) if a.root else rows[0]
    row={**row,'source':str((origin/'inputs').relative_to(project))}
    evaluate.physical_window(row,frozen['selected'])
elif a.stage in ('controls','optimizer_controls'):
    getattr(evaluate,a.stage)()
else:
    # Copy only small waveform/truth/metadata exports, not the original ROOTs.
    shutil.copytree(origin/'inputs',target/'inputs')
    evaluate.physical(2);evaluate.neural();evaluate.controls();evaluate.optimizer_controls();evaluate.shifts()
