"""Serialization repair for the frozen controls; predictor/source stay unchanged.

The original controls completed numerical work but tried to JSON-encode the
physical model arrays. Move only those saved arrays to NPZ, then save JSON.
"""
from pathlib import Path
import os
import sys

project=Path(__file__).resolve().parents[3]
os.environ['GAGG_COMPONENT_PROJECT']=str(project)
sys.path.insert(0,str(Path(__file__).resolve().parent/'inference-source'))
from pileup.analysis import alpha_component_evaluate as evaluation
from pileup.analysis.alpha_rate_campaign import save,digest
import numpy as np


def array_safe_save(path,data):
    if Path(path).name=='controls.json':
        folder=evaluation.OUTPUT/'control-models';folder.mkdir(exist_ok=True)
        for row in data['cases']:
            for name,result in row['methods'].items():
                if 'model_mV' not in result['fit']:continue
                target=folder/f'{row["noise"]}_{row["label"]}_{name}.npz'
                np.savez_compressed(target,model_mV=result['fit'].pop('model_mV'))
                result['saved_model']=str(target.relative_to(project))
    save(path,data)


if __name__=='__main__':
    evaluation.check_freeze()
    assert not (evaluation.OUTPUT/'controls.json').exists()
    evaluation.save=array_safe_save
    evaluation.controls()
    record=project/'jupyter/run_records/alpha_component_ml_20261021'
    save(record/'control-serialization-repair.json',dict(
        cause='Original frozen controls tried to JSON-encode model_mV ndarrays after all62 cases.',
        correction='Same62 cases and original seed regenerated; only move model_mV to NPZ before JSON serialization.',
        inference_source_sha256=digest(evaluation.__file__),adapter_sha256=digest(__file__),
        inference_parameters_changed=False,case_or_seed_replacements=0))
