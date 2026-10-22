"""A small GRU event head on the same ROOT split/features as the frozen CNN.

This is a project adaptation, not a reproduction of Regadio2024. Keep the
previous supervised objective so architecture is the principal new variable.
It is offline sliding-window inference, with no claim of streaming latency.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime
import json
from pathlib import Path
import time

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

from .alpha_correlated_study import PROJECT,EXPERIMENT,OUTPUT
from .alpha_ml_models import windowed_prediction
from .alpha_ml_multimethod import choice
from .alpha_pileup_ml import load_pilot,loss_function
from .alpha_rate_campaign import digest,save


class GRUEvent(nn.Module):
    """8 ns recurrent tokens plus a 2 ns local branch; three event heads."""
    def __init__(self):
        super().__init__()
        self.embed=nn.Conv1d(5,12,9,stride=4,padding=4)
        self.gru=nn.GRU(12,24,batch_first=True)
        self.local=nn.Sequential(nn.Conv1d(5,8,5,padding=2),nn.SiLU())
        self.head=nn.Sequential(nn.Conv1d(32,24,5,padding=2),nn.SiLU(),nn.Conv1d(24,3,1))

    def forward(self,x):
        tokens=F.silu(self.embed(x)).transpose(1,2)
        recurrent,_=self.gru(tokens)
        context=F.interpolate(recurrent.transpose(1,2),size=x.shape[-1],mode='linear',align_corners=False)
        return self.head(torch.cat((context,self.local(x)),dim=1))


def train(steps=2400,threads=2):
    folder=OUTPUT/'gru';folder.mkdir(exist_ok=False)
    torch.set_num_threads(threads);torch.manual_seed(262025009);rng=np.random.default_rng(262025009)
    parent=PROJECT/'pileup/experiments/alpha_ml_comparison_20261014'
    settings=json.loads((parent/'calibration-settings.json').read_text())
    split=json.loads((parent/'data-split.json').read_text())
    items=[load_pilot(r,settings) for r in split['windows']]
    training=[r for r in items if r['record']['split']=='train'];validation=[r for r in items if r['record']['split']=='validation']
    assert len(training)==26 and len(validation)==10
    assert not {r['record']['root_sha256'] for r in training}&{r['record']['root_sha256'] for r in validation}
    # No new test waveform or label is opened by this module.
    save(folder/'data-split.json',split)
    config=dict(architecture='Small unidirectional GRU on 8 ns tokens with native-grid convolution head',
        steps=steps,batch_size=16,seed=262025009,optimizer='AdamW',learning_rate=.001,weight_decay=.0001,
        crop_samples=1024,training_core_samples=500,guard_samples=262,threads=threads,
        loss='Unchanged soft BCE + 2 log-charge Huber + .5 time-offset Huber; no target-delay relaxation',
        inference='Sliding 1024-point windows, each guarded 500-point core emitted once',
        input_features='Same five frozen voltage-derived features as prior CNN',
        source_sha256=digest(Path(__file__)),calibration=settings,
        calibration_sha256=digest(parent/'calibration-settings.json'),split_sha256=digest(parent/'data-split.json'),
        parameters=sum(p.numel() for p in GRUEvent().parameters()),torch=torch.__version__,
        model_selection='Validation-only macro complete-recovery F1 over thresholds .15 to .85',
        paper_reproduction=False,truth_input_at_inference=[],rate_input_at_inference=False,
        comparison_limit='Same training ROOTs, features, objective, steps and batch size; different seed and architecture; not an equal wall-time benchmark')
    save(folder/'config.json',config)
    groups=defaultdict(list)
    for item in training:groups[item['record']['rate_hz']].append(item)
    rates=sorted(groups);model=GRUEvent();optimizer=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001)
    history=[];best=None;tick=time.monotonic()
    print('GRU_TRAIN',config['parameters'],flush=True)
    for step in range(1,steps+1):
        x=[];labels=[]
        for _ in range(16):
            group=groups[rates[int(rng.integers(len(rates)))]];item=group[int(rng.integers(len(group)))];start=int(rng.integers(2250,6739))
            x.append(item['features'][:,start:start+1024]);labels.append(item['targets'][:,start:start+1024])
        model.train();optimizer.zero_grad(set_to_none=True)
        loss,parts=loss_function(model(torch.from_numpy(np.stack(x))),torch.from_numpy(np.stack(labels)))
        if not torch.isfinite(loss):raise FloatingPointError('Nonfinite GRU loss')
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
        if step==1 or step%200==0 or step==steps:
            model.eval();selected=choice(model,validation)
            row=dict(step=step,loss=float(loss.detach()),loss_terms=[float(p.detach()) for p in parts],
                validation_objective=selected['objective'],threshold=selected['threshold'],elapsed_s=time.monotonic()-tick)
            history.append(row);save(folder/'history.json',history)
            if best is None or selected['objective']>best['objective']:
                best=dict(selected,step=step);torch.save(model.state_dict(),folder/'model.pt');save(folder/'validation-selection.json',best)
            print('GRU_STEP',step,round(row['loss'],4),round(row['validation_objective'],4),round(row['elapsed_s'],1),flush=True)
    frozen=dict(frozen_at=datetime.now().astimezone().isoformat(),weights_sha256=digest(folder/'model.pt'),
        config_sha256=digest(folder/'config.json'),split_sha256=digest(folder/'data-split.json'),
        threshold=best['threshold'],step=best['step'],validation_objective=best['objective'],
        training_seconds=time.monotonic()-tick,parameters=config['parameters'],source_sha256=digest(Path(__file__)))
    save(folder/'frozen-model.json',frozen);print('GRU_FROZEN',frozen['step'],frozen['threshold'],flush=True)


def load():
    folder=OUTPUT/'gru';frozen=json.loads((folder/'frozen-model.json').read_text())
    assert digest(folder/'model.pt')==frozen['weights_sha256']
    assert digest(folder/'config.json')==frozen['config_sha256']
    assert digest(Path(__file__))==frozen['source_sha256']
    model=GRUEvent();model.load_state_dict(torch.load(folder/'model.pt',map_location='cpu',weights_only=True));model.eval()
    return model,frozen,json.loads((folder/'config.json').read_text())['calibration']


if __name__=='__main__':train()
