"""Component-supervised neural training with a differentiable physical decoder.

Run from the repository root. This is a GAGG adaptation inspired by SAEFit,
not a reproduction. Training can use native ancestry. Blind inference accepts
only sampled voltage and frozen calibration/weights. Unknown multiplicity and
the old 2 ns/8 ns dense event decoder are retained for a controlled comparison.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import replace
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

from .alpha_correlated_gru import GRUEvent
from .alpha_full_signal_study import make_template
from .alpha_ml_models import UNetEvent, windowed_prediction
from .alpha_pileup_ml import waveform_features, supervised_targets, loss_function, decode_prediction
from .alpha_pileup_fit import score
from .alpha_rate_campaign import digest, save
from .sipm_readout import ReadoutConfig, synthesize_readout
from .sipm_response import AvalancheResult, ELECTRON_CHARGE_C

PROJECT = Path(os.environ.get('GAGG_COMPONENT_PROJECT', Path(__file__).resolve().parents[2])).resolve()
EXPERIMENT = PROJECT / 'pileup/experiments/alpha_component_ml_20261021'
OUTPUT = PROJECT / 'jupyter/artifacts/alpha_component_ml_20261021'
OFFSETS = np.arange(-16., 1200.1, 2.)
THRESHOLDS = (.15, .25, .35, .45, .55, .65, .75, .85, .9, .95)
VARIANTS = ('gru_base', 'gru_component', 'unetgru_base', 'unetgru_component')
SEEDS = (262041001, 262041002)


class UNetGRUEvent(UNetEvent):
    """Old U-Net with a unidirectional GRU at its 32 ns bottleneck."""
    def __init__(self):
        super().__init__()
        self.recurrent = nn.GRU(64, 64, batch_first=True)

    def forward(self, x):
        skips = []
        for i, encoder in enumerate(self.encoders):
            if i:
                x = F.avg_pool1d(x, 2)
            x = encoder(x)
            skips.append(x)
        x, _ = self.recurrent(x.transpose(1, 2))
        x = x.transpose(1, 2)
        for decoder, skip in zip(self.decoders, reversed(skips[:-1])):
            x = F.interpolate(x, size=skip.shape[-1], mode='linear', align_corners=False)
            x = decoder(torch.cat((x, skip), dim=1))
        return self.head(x)


def model_for(variant):
    if variant not in VARIANTS:
        raise ValueError(variant)
    return GRUEvent() if variant.startswith('gru_') else UNetGRUEvent()


class PhysicalDecoder(nn.Module):
    """Fixed dual-decay, per-charge voltage response with differentiable shifts.

    Linear interpolation is piecewise differentiable in event time. No truth,
    source rate or fitted count is stored in this decoder.
    """
    def __init__(self, settings):
        super().__init__()
        template = make_template(settings)
        grid = np.arange(-24., 6008.001, .25)
        values = template(grid) * settings['fit_config']['area_per_pe_mV_ns'] * template.unit_area_peak
        self.register_buffer('kernel', torch.tensor(values, dtype=torch.float32))
        self.register_buffer('offsets', torch.tensor(OFFSETS, dtype=torch.float32))
        self.grid_start, self.grid_step = -24., .25

    def response(self, relative_ns):
        p = (relative_ns - self.grid_start) / self.grid_step
        lo = torch.floor(p).long()
        valid = (lo >= 0) & (lo < self.kernel.numel() - 1)
        lo = lo.clamp(0, self.kernel.numel() - 2)
        w = p - lo.to(p.dtype)
        return ((1-w) * self.kernel[lo] + w * self.kernel[lo+1]) * valid

    def components(self, charge, offset):
        return charge[:, None] * self.response(self.offsets[None, :] - offset[:, None])


def native_component_cache(record, settings, destination):
    """Cache actual native components, without redigitizing independent pulses.

    Original avalanche times/charges retain shared-cell history. Only the linear
    electronics are replayed, with exact grid-aligned local origins. Full-Q is
    checked before restricting the display/loss support. Loss covers 1216 ns;
    that support is neither an acquisition-window recommendation nor Q integral.
    """
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    with np.load(PROJECT / record['waveforms']) as f:
        t, y = f['time_ns'].copy(), f['analog_mV'].copy()
    metadata = json.loads((PROJECT / record['metadata']).read_text())
    assert metadata['root_sha256'] == record['root_sha256']
    with np.load(PROJECT / record['truth']) as f:
        truth = {k: f[k] for k in ('first_active_ns', 'emission_ns', 'charge_pe',
                                  'avalanche_time_ns', 'avalanche_charge_C', 'avalanche_owner')}
    if not np.allclose(t, np.arange(len(t))*2. + .25, atol=1e-8, rtol=0):
        raise ValueError('Native aperture/grid differs from the declared 2 ns configuration')
    q0 = metadata['config']['sensor']['gain_electrons'] * ELECTRON_CHARGE_C
    readout = replace(ReadoutConfig(**metadata['config']['readout']), baseline_mV=0., noise_rms_mV=0.)
    owners = truth['avalanche_owner']
    order = np.argsort(owners, kind='stable')
    records = np.empty(len(order), dtype=[('time_ns', 'f8'), ('charge_C', 'f8')])
    records['time_ns'] = truth['avalanche_time_ns'][order]
    records['charge_C'] = truth['avalanche_charge_C'][order]
    bounds = np.r_[0, np.cumsum(np.bincount(owners+1, minlength=len(truth['charge_pe'])+1))]
    times, charges = truth['first_active_ns'], truth['charge_pe']
    ids = np.flatnonzero(np.isfinite(times) & (charges >= 500.) & (times >= 4900.) & (times < 15100.))
    anchors = np.rint((times[ids] - t[0])/2.).astype(int)
    unique, counts = np.unique(anchors, return_counts=True)
    collided = set(unique[counts > 1].tolist())
    valid = np.array([a not in collided for a in anchors], dtype=bool)
    collision_count = int((~valid).sum())
    ids, anchors = ids[valid], anchors[valid]
    curves = np.empty((len(ids), len(OFFSETS)), dtype=np.float32)
    max_charge_error, max_local_error = 0., 0.
    for j, (i, anchor) in enumerate(zip(ids, anchors)):
        block = records[bounds[i+1]:bounds[i+2]]
        max_charge_error = max(max_charge_error, abs(float(block['charge_C'].sum()/q0-charges[i])))
        origin = t[anchor] - 16. - readout.aperture_ns/2
        if len(block) and block['time_ns'].min() < origin:
            raise ValueError('A local component would omit its prehistory')
        local = block[block['time_ns'] < origin + 1232.].copy()
        local['time_ns'] -= origin
        signal = synthesize_readout(AvalancheResult(local, 1232., q0, np.array([]), {}), readout, 0)
        curves[j] = signal.arrays['signal_mV'][np.arange(len(OFFSETS))]
        # One full-record reference per ROOT verifies the local-origin shortcut.
        if j == 0:
            full = synthesize_readout(AvalancheResult(block, float(metadata['duration_ns']), q0,
                                                      np.array([]), {}), readout, 0)
            reference = full.arrays['signal_mV'][anchor-8+np.arange(len(OFFSETS))]
            max_local_error = float(np.max(abs(curves[j].astype(float)-reference)))
    # Float32 target storage tolerance is separate from exact float64 readout.
    if max_charge_error > 1e-7 or max_local_error > .003:
        raise ValueError((max_charge_error, max_local_error))
    features = waveform_features(t, y, settings)
    targets = supervised_targets(t, truth)
    decoder = PhysicalDecoder(settings)
    with torch.no_grad():
        mean = decoder.components(torch.tensor(charges[ids], dtype=torch.float32),
                                  torch.tensor(times[ids]-t[anchors], dtype=torch.float32)).numpy()
        norm = np.sum((charges[ids, None] * decoder.response(decoder.offsets).numpy())**2, axis=1)
    floor = np.sum((mean-curves)**2, axis=1)/np.maximum(norm, 1e-8)
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, t=t, features=features, targets=targets,
        times=times, charges=charges, emission_ns=truth['emission_ns'],
        event_ids=ids, anchors=anchors, curves=curves, component_norm=norm.astype(np.float32),
        mean_decoder_floor=floor.astype(np.float32))
    return dict(events=len(ids), same_bin_excluded_from_component_loss=collision_count,
                evaluation_exclusions=0, charge_max_error_pe=max_charge_error,
                local_full_max_error_mV=max_local_error,
                mean_template_floor_median=float(np.median(floor)) if len(floor) else None,
                mean_template_floor_p90=float(np.quantile(floor,.9)) if len(floor) else None)


def batch_components(items, starts):
    batches, anchors, curves, norms = [], [], [], []
    # Unique closest-bin events only. Other base labels and scoring remain intact.
    for b, (item, start) in enumerate(zip(items, starts)):
        keep = (item['anchors'] >= start+262) & (item['anchors'] < start+762)
        n = int(keep.sum())
        batches.extend([b]*n)
        anchors.extend((item['anchors'][keep]-start).tolist())
        curves.extend(item['curves'][keep])
        norms.extend(item['component_norm'][keep])
    n = len(batches)
    return dict(batch=torch.tensor(batches, dtype=torch.long),
                anchor=torch.tensor(anchors, dtype=torch.long),
                curve=torch.tensor(np.asarray(curves, np.float32).reshape(n, len(OFFSETS))),
                norm=torch.tensor(norms, dtype=torch.float32))


def component_loss(z, targets, decoder):
    if targets['batch'].numel() == 0:
        zero = z.sum()*0
        return zero, zero
    selected = z[targets['batch'], :, targets['anchor']]
    q = 3000. * torch.exp(selected[:, 1].clamp(-5., 5.))
    offset = 4. * torch.tanh(selected[:, 2])
    predicted = decoder.components(q, offset)
    # Per-component normalization prevents the strongest event dominating.
    individual = ((predicted-targets['curve']).square().sum(1)/targets['norm'].clamp_min(1e-8)).mean()
    # Conditional sum on the common native 2 ns grid. These centers are training
    # labels, never supplied at inference; this is not an unsupervised full-y loss.
    width = z.shape[-1] + len(OFFSETS)
    indices = (targets['batch'][:,None]*width + targets['anchor'][:,None]
               + torch.arange(len(OFFSETS))[None,:]).flatten()
    out = torch.zeros(z.shape[0]*width, dtype=z.dtype).index_add(0, indices, predicted.flatten())
    truth = torch.zeros_like(out).index_add(0, indices, targets['curve'].flatten())
    out, truth = out.reshape(z.shape[0],width), truth.reshape(z.shape[0],width)
    nonempty = truth.square().sum(1) > 0
    total = ((out-truth).square().sum(1)/truth.square().sum(1).clamp_min(1.))[nonempty].mean()
    return individual, total


def validation_grid(model, items):
    dense = [windowed_prediction(model, item['features']) for item in items]
    choices = []
    for threshold in THRESHOLDS:
        groups = defaultdict(lambda: np.zeros(4))
        for item, pred in zip(items, dense):
            fit = decode_prediction(item['t'], pred, threshold)
            s = score(fit, item['truth'], [5000.,15000.])
            groups[item['record']['rate_hz']] += [s['target_count'], s['reconstructed_count'],
                                                 s['recovered_count'], s['unmatched_predictions']]
        objective = float(np.mean([2*v[2]/max(v[0]+v[1],1) for v in groups.values()]))
        choices.append(dict(threshold=threshold, objective=objective,
                            extra_count=int(sum(v[3] for v in groups.values())),
                            recovered_count=int(sum(v[2] for v in groups.values())),
                            by_rate={str(k):v.tolist() for k,v in groups.items()}))
    return max(choices, key=lambda c:(c['objective'],c['threshold'])), choices


def load_items():
    split = json.loads((EXPERIMENT/'data-split.json').read_text())
    items = []
    for record in split['windows']:
        with np.load(PROJECT/record['component_cache']) as f:
            item = {k:f[k] for k in f.files}
        item.update(record=record, truth=dict(first_active_ns=item['times'],charge_pe=item['charges'],
                                            emission_ns=item['emission_ns']))
        items.append(item)
    return items


def train(variant, seed, steps=2400, threads=2):
    if variant not in VARIANTS or seed not in SEEDS or not 1 <= threads <= 4:
        raise ValueError('Undeclared training run')
    destination = OUTPUT/'models'/f'{variant}_{seed}'
    destination.mkdir(parents=True, exist_ok=False)
    torch.set_num_threads(threads)
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    settings = json.loads((EXPERIMENT/'calibration-settings.json').read_text())
    items = load_items()
    training = [i for i in items if i['record']['split']=='train']
    validation = [i for i in items if i['record']['split']=='validation']
    groups = defaultdict(list)
    for item in training:
        groups[item['record']['rate_hz']].append(item)
    rates = sorted(groups)
    model, decoder = model_for(variant), PhysicalDecoder(settings)
    optimizer = torch.optim.AdamW(model.parameters(), lr=.001, weight_decay=.0001)
    config = dict(variant=variant,seed=seed,steps=steps,batch_size=16,threads=threads,
        parameters=sum(p.numel() for p in model.parameters()),torch=torch.__version__,
        source_sha256=digest(__file__),split_sha256=digest(EXPERIMENT/'data-split.json'),
        calibration_sha256=digest(EXPERIMENT/'calibration-settings.json'),
        component_weight=1. if variant.endswith('_component') else 0., conditional_sum_weight=.1,
        component_loss='Native ancestry voltage vs fixed physical decoder at supervised nearest-bin anchors',
        supervised_center_use='Training loss only; unknown event count at blind inference',
        inference_truth_inputs=[],rate_input_at_inference=False, candidate_separation_ns=8.,
        scope='Fixed 95/351 ns; Q is realized avalanche-equivalent pe, not deposited energy')
    save(destination/'config.json', config)
    best, history, started = None, [], time.monotonic()
    print('TRAIN_BEGIN',variant,seed,config['parameters'],flush=True)
    for step in range(1, steps+1):
        selected_items, starts = [], []
        for _ in range(16):
            group = groups[rates[int(rng.integers(len(rates)))]]
            selected_items.append(group[int(rng.integers(len(group)))])
            starts.append(int(rng.integers(2250,6739)))
        x = torch.tensor(np.stack([i['features'][:,s:s+1024] for i,s in zip(selected_items,starts)]))
        labels = torch.tensor(np.stack([i['targets'][:,s:s+1024] for i,s in zip(selected_items,starts)]))
        model.train()
        optimizer.zero_grad(set_to_none=True)
        z = model(x)
        base, parts = loss_function(z, labels)
        comp, total = z.sum()*0, z.sum()*0
        if config['component_weight']:
            comp, total = component_loss(z, batch_components(selected_items,starts), decoder)
        loss = base + config['component_weight']*(comp + .1*total)
        if not torch.isfinite(loss):
            raise FloatingPointError('Nonfinite component training loss')
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(),5.)
        optimizer.step()
        if step == 1 or step % 400 == 0 or step == steps:
            model.eval()
            choice, grid = validation_grid(model,validation)
            row = dict(step=step,loss=float(loss.detach()),base_loss=float(base.detach()),
                       component_loss=float(comp.detach()),conditional_sum_loss=float(total.detach()),
                       validation_objective=choice['objective'],threshold=choice['threshold'],
                       elapsed_s=time.monotonic()-started)
            history.append(row)
            save(destination/'history.json',history)
            if best is None or choice['objective'] > best['objective']:
                best = dict(choice,step=step)
                torch.save(model.state_dict(),destination/'model.pt')
                save(destination/'selection.json',dict(best=best,grid=grid))
            print('TRAIN_STEP',variant,seed,step,round(row['loss'],4),round(choice['objective'],4),
                  round(row['elapsed_s'],1),flush=True)
    frozen = dict(frozen_at=datetime.now().astimezone().isoformat(),weights_sha256=digest(destination/'model.pt'),
        config_sha256=digest(destination/'config.json'),source_sha256=digest(__file__),
        threshold=best['threshold'],step=best['step'],validation_objective=best['objective'],
        training_seconds=time.monotonic()-started,variant=variant,seed=seed)
    save(destination/'frozen.json',frozen)
    print('TRAIN_FROZEN',variant,seed,best['step'],best['threshold'],flush=True)


def load_model(name):
    folder=OUTPUT/'models'/name
    frozen=json.loads((folder/'frozen.json').read_text())
    assert frozen['weights_sha256']==digest(folder/'model.pt')
    assert frozen['config_sha256']==digest(folder/'config.json')
    assert frozen['source_sha256']==digest(__file__)
    model=model_for(frozen['variant'])
    model.load_state_dict(torch.load(folder/'model.pt',map_location='cpu',weights_only=True))
    return model.eval(),frozen


def blind_predict(model, t_ns, voltage_mV, settings, threshold):
    """Inference boundary: no record metadata, rate, ancestry or truth accepted."""
    dense = windowed_prediction(model, waveform_features(t_ns, voltage_mV, settings))
    return decode_prediction(np.asarray(t_ns), dense, threshold), dense
