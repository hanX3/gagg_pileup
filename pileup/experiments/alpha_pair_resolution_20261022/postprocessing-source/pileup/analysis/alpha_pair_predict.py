"""Inference-only command: a complete 1024-sample local voltage frame."""
from pathlib import Path
import argparse

import numpy as np
import torch

from .alpha_pair_registry import load_model,predict
from .alpha_rate_campaign import save,digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',required=True)
    p.add_argument('--waveforms',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    torch.set_num_threads(1)
    with np.load(a.waveforms,allow_pickle=False) as f:
        # Load only the voltage field, never iterate over possible label arrays.
        voltage=f['voltage_mV'] if 'voltage_mV' in f.files else f['analog_mV']
    if voltage.shape[-1]!=1024 or not np.isfinite(voltage).all():
        raise ValueError('This model requires1024 samples at2ns; it is not a20us high-multiplicity decoder')
    model,cfg,selected=load_model(a.model)
    fits,arrays=predict(model,voltage,selected.get('threshold',.5))
    a.output.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(a.output/'predictions.npz',**arrays)
    save(a.output/'events.json',dict(fits=fits,truth_inputs=[],model=a.model,
         waveform_sha256=digest(a.waveforms),known_single=cfg['single_known_K'],
         source_kind='Only observed voltage and frozen calibration/weights'))


if __name__=='__main__':main()
