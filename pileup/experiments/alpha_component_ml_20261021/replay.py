"""Replay inference or training from this experiment's frozen Python snapshot.

Dependencies and original local waveform/training-cache/model files are needed;
no Git checkout or modification of the original results is required.
"""
from pathlib import Path
import argparse
import json
import os
import sys

project=Path(__file__).resolve().parents[3]
os.environ['GAGG_COMPONENT_PROJECT']=str(project)
sys.path.insert(0,str(Path(__file__).resolve().parent/'inference-source'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('infer','train'))
    p.add_argument('--waveform',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--variant')
    p.add_argument('--seed',type=int)
    p.add_argument('--threads',type=int,default=2)
    args=p.parse_args()
    from pileup.analysis import alpha_component_ml as ml
    from pileup.analysis import alpha_component_evaluate as evaluation
    from pileup.analysis.alpha_rate_campaign import digest,save
    if args.output.exists():raise FileExistsError('Choose a new output directory')
    evaluation.check_freeze()
    if args.stage=='train':
        ml.OUTPUT=args.output.resolve()
        ml.train(args.variant,args.seed,steps=2400,threads=args.threads)
        return
    if args.waveform is None:p.error('--waveform is required for infer')
    import numpy as np
    nets=evaluation.networks()
    with np.load(args.waveform) as f:
        t,y=f['time_ns'].copy(),f['analog_mV'].copy()
    fits,seconds,dense=evaluation.infer(t,y,nets)
    args.output.mkdir(parents=True)
    for name,fit in fits.items():
        if 'model_mV' in fit:
            np.savez_compressed(args.output/(name+'.npz'),model_mV=fit.pop('model_mV'))
    for name,pred in dense.items():np.savez_compressed(args.output/(name+'_dense.npz'),**pred)
    save(args.output/'predictions.json',dict(waveform=str(args.waveform),waveform_sha256=digest(args.waveform),
        fits=fits,seconds=seconds,truth_inputs=[],rate_inputs=[],source='Frozen experiment snapshot'))
    print(args.output/'predictions.json')


if __name__=='__main__':main()
