"""Retrain a frozen ordered model into a new project-local output directory."""
from pathlib import Path
import argparse
import os
import shutil
import sys

folder=Path(__file__).resolve().parent
project=folder.parents[2]
os.environ['GAGG_PAIR_PROJECT']=str(project)
sys.path.insert(0,str(folder/'inference-source'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--variant',choices=('pair_ordered_direct','pair_ordered_physical'),required=True)
    p.add_argument('--seed',type=int,choices=(262050001,262050002),required=True)
    p.add_argument('--destination',type=Path,required=True)
    a=p.parse_args()
    from pileup.analysis import alpha_pair_ordered as model
    original=model.OUTPUT
    destination=a.destination.resolve()
    if project not in destination.parents:
        raise ValueError('Choose a new output directory inside the project for relative provenance paths')
    destination.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(original/'single-selection.json',destination/'single-selection.json')
    shutil.copyfile(original/'train-shared-pair.npz',destination/'train-shared-pair.npz')
    chosen=model.load(original/'single-selection.json')['variant']
    pretrained=destination/'models'/f'{chosen}_{a.seed}'
    pretrained.mkdir(parents=True)
    shutil.copyfile(original/'models'/f'{chosen}_{a.seed}'/'model.pt',pretrained/'model.pt')
    # Only output/pretrained-copy paths move. Native train/validation pools,
    # calibration and validation controls keep their original immutable paths.
    model.OUTPUT=destination
    model.train(a.variant,a.seed)


if __name__=='__main__':main()
