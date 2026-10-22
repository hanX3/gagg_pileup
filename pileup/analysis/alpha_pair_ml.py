"""Single-pulse calibration and two-component supervised reconstruction.

This is a new project experiment, not a paper reproduction. Two independent
event slots replace nearest-bin labels/8 ns peak suppression. The direct head
partitions positive observed voltage into two components and a nuisance part;
the physical head partitions its integrated charge and uses fixed dual-tau
responses. Counts are inferred through presence probabilities. Truth is used
only by training losses and by separate scoring functions.
"""
from __future__ import annotations

import itertools
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

from .alpha_component_ml import PhysicalDecoder

SAMPLES = 1024
DT_NS = 2.


class ComponentNet(nn.Module):
    """Observed voltage -> presence, first-active times, Q and component curves.

    Charge integration is an explicit inductive bias, in the same avalanche-pe
    units as the fixed electronics calibration. Finite-window tail correction
    uses the mean response and learned time, never a true arrival or true Q.
    Direct curves are voltage masks, not predetermined exponential curves.
    """
    def __init__(self, settings, *, slots=2, head='direct', samples=SAMPLES):
        super().__init__()
        if slots not in (1, 2) or head not in ('direct', 'physical'):
            raise ValueError((slots, head))
        self.slots, self.head, self.samples = slots, head, samples
        self.baseline = float(settings['fit_config']['baseline_mV'])
        self.area_per_pe = float(settings['fit_config']['area_per_pe_mV_ns'])
        self.decoder = PhysicalDecoder(settings)
        self.register_buffer('grid', torch.arange(samples)*DT_NS + .25)
        self.encoder = nn.Sequential(
            nn.Conv1d(4, 16, 9, stride=2, padding=4), nn.SiLU(),
            nn.Conv1d(16, 24, 7, stride=2, padding=3), nn.SiLU(),
            nn.Conv1d(24, 32, 5, stride=2, padding=2), nn.SiLU(),
            nn.Conv1d(32, 32, 5, padding=2), nn.SiLU())
        self.fine = nn.Sequential(nn.Conv1d(36, 24, 5, padding=2), nn.SiLU(),
                                  nn.Conv1d(24, 24, 5, padding=2), nn.SiLU())
        self.time_head = nn.Conv1d(24, slots, 1)
        self.mask_head = nn.Conv1d(24, slots+1, 1)
        self.global_head = nn.Sequential(nn.Linear(32*16+1, 64), nn.SiLU(),
                                         nn.Linear(64, 2*slots+2))

    def forward(self, voltage):
        if voltage.ndim != 2 or voltage.shape[1] != self.samples:
            raise ValueError('Expected a batch of complete local voltage frames')
        y = voltage - self.baseline
        positive = y.clamp_min(0.)
        derivative = F.pad((y[:, 2:]-y[:, :-2])/4., (1, 1))
        smooth = F.avg_pool1d(y[:, None], 33, stride=1, padding=16)[:, 0]
        coord = (self.grid/self.grid[-1]*2-1).expand_as(y)
        x = torch.stack((torch.asinh(y/500.), torch.asinh(derivative/10.),
                         torch.asinh(smooth/500.), coord), dim=1)
        encoded = self.encoder(x)
        fine = self.fine(torch.cat((F.interpolate(encoded, size=self.samples,
                            mode='linear', align_corners=False), x), dim=1))
        logits = self.time_head(fine)
        density = F.softmax(logits, dim=-1)
        times = (density*self.grid).sum(-1)
        raw_area_q = positive.sum(-1)*DT_NS/self.area_per_pe
        pooled = F.adaptive_avg_pool1d(encoded, 16).flatten(1)
        global_z = self.global_head(torch.cat((pooled, torch.log1p(raw_area_q/3000.)[:, None]), 1))
        presence_logits = global_z[:, :self.slots]
        unit = self.decoder.response(self.grid[None, None, :] - times[:, :, None])
        fraction = (unit.sum(-1)*DT_NS/self.area_per_pe).clamp_min(.02)
        if self.slots == 1:
            # Known-single diagnostic. The bounded residual corrects native
            # nuisance charge and finite-window/mean-shape differences.
            q = raw_area_q[:, None]/fraction * torch.exp(.25*torch.tanh(global_z[:, -1:]))
            curves = q[:, :, None]*unit
        elif self.head == 'direct':
            masks = F.softmax(self.mask_head(fine), dim=1)
            curves = masks[:, :self.slots]*positive[:, None]
            q = curves.sum(-1)*DT_NS/self.area_per_pe/fraction
        else:
            weights = F.softmax(global_z[:, self.slots:self.slots+self.slots+1], dim=1)
            q = raw_area_q[:, None]*weights[:, :self.slots]/fraction
            curves = q[:, :, None]*unit
        return dict(presence_logits=presence_logits, time_logits=logits,
                    times_ns=times, charge_pe=q.clamp_min(.01), components_mV=curves)


def set_loss(pred, target, settings, *, normalization='fixed', component_weight=1., single=False):
    """Exactly enumerate two assignments; no label/count is an inference input.

    Fixed normalization is a training-pool global constant. Relative uses each
    true component norm; both are independently compared, not assumed unbiased.
    The time cross entropy is a training target, never test-time realignment.
    """
    slots = pred['charge_pe'].shape[1]
    permutations = list(itertools.permutations(range(slots)))
    costs, details = [], []
    grid = torch.arange(pred['time_logits'].shape[-1], device=pred['time_logits'].device)*2.+.25
    for permutation in permutations:
        order = list(permutation)
        active = target['active'][:, order]
        true_q = target['charge_pe'][:, order]
        true_t = target['times_ns'][:, order]
        actual = target['components_mV'][:, order]
        pl = F.binary_cross_entropy_with_logits(pred['presence_logits'], active, reduction='none')
        ql = F.smooth_l1_loss(torch.log(pred['charge_pe']/true_q.clamp_min(1.)),
                             torch.zeros_like(true_q), reduction='none')*active
        time_target = torch.exp(-.5*((grid[None, None]-true_t[:, :, None])/2.)**2)
        time_target = time_target/time_target.sum(-1, keepdim=True).clamp_min(1e-20)
        tl = -(time_target*F.log_softmax(pred['time_logits'], -1)).sum(-1)*active
        error = (pred['components_mV']-actual).square().sum(-1)
        fixed = float(settings['component_fixed_norm'])
        if normalization == 'relative':
            norm = torch.where(active.bool(), actual.square().sum(-1).clamp_min(1.),
                               torch.full_like(active, fixed))
        elif normalization == 'fixed':
            norm = torch.full_like(active, fixed)
        else:
            raise ValueError(normalization)
        cl = error/norm
        # Known-single stage reports K=1 explicitly and does not train a count.
        cost = (0. if single else 1.)*pl + 6.*ql + .5*tl + component_weight*cl
        costs.append(cost.sum(-1))
        details.append(torch.stack((pl.sum(-1), ql.sum(-1), tl.sum(-1), cl.sum(-1)), -1))
    matrix = torch.stack(costs, -1)
    assignment = matrix.detach().argmin(-1)
    index = torch.arange(len(assignment), device=assignment.device)
    parts = torch.stack(details, 1)[index, assignment].mean(0)
    return matrix[index, assignment].mean(), dict(zip(('presence', 'logq', 'time_ce', 'component'), parts))


def blind_predict(model, voltage_mV, *, threshold=.5):
    """Only the local observed voltage, model and frozen threshold are accepted."""
    y = np.asarray(voltage_mV, np.float32)
    one = y.ndim == 1
    if one:
        y = y[None]
    with torch.inference_mode():
        output = model(torch.from_numpy(y))
    arrays = {key:value.cpu().numpy() for key,value in output.items()}
    presence = 1/(1+np.exp(-np.clip(arrays['presence_logits'], -40, 40)))
    result = []
    for i in range(len(y)):
        keep = np.ones(model.slots, bool) if model.slots == 1 else presence[i] >= threshold
        if model.slots != 1:
            keep &= arrays['charge_pe'][i] >= 500.
        ids = np.flatnonzero(keep)
        ids = ids[np.argsort(arrays['times_ns'][i, ids])]
        result.append(dict(status='predicted', times_ns=arrays['times_ns'][i, ids].astype(float).tolist(),
                           charge_pe=arrays['charge_pe'][i, ids].astype(float).tolist(),
                           slot_ids=ids.tolist(), presence=presence[i].astype(float).tolist(),
                           known_single=model.slots == 1, truth_inputs=[]))
    return (result[0] if one else result), arrays


def tensors(batch):
    return {key:torch.from_numpy(np.asarray(value, np.float32)) for key,value in batch.items()}
