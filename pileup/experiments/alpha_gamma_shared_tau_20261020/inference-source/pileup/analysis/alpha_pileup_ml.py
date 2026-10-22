"""Exploratory dense neural event reconstruction on native shared-cell SiPM data.

This is a small project baseline inspired by time-sequence unfolding, not a
reproduction of a published network. At inference only complete sampled noisy
voltage, its grid, frozen calibration and learned weights are accepted. Entire
ROOTs are assigned to train/validation before random crops. Truth is used for
supervised training and separate evaluation, never for blind prediction.

Optional dependency: CPU PyTorch (see requirements-ml.txt). Run from project root.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import replace
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import time

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import torch
from torch import nn
from torch.nn import functional as F

from .alpha_rate_campaign import digest, save
from .alpha_pileup_fit import FitConfig, candidate_peaks, fit_waveform, score
from .alpha_full_signal_study import make_template, summarize


class EventNet(nn.Module):
    """Five calibrated waveform features -> probability, charge and time offset.

    Seven dilated residual blocks, 24 channels, no fixed pulse multiplicity.
    Receptive field is 517 samples = 1034 ns at this study's 2 ns sampling.
    The inverse-filter feature is computed on the full voltage record; this is
    offline reconstruction, not a causal streaming/latency benchmark.
    """
    def __init__(self, channels=24):
        super().__init__()
        self.stem = nn.Conv1d(5, channels, 9, padding=4)
        self.blocks = nn.ModuleList([
            nn.Sequential(nn.Conv1d(channels, channels, 5, padding=2*d, dilation=d),
                          nn.SiLU(), nn.Conv1d(channels, channels, 1))
            for d in (1, 2, 4, 8, 16, 32, 64)])
        self.head = nn.Conv1d(channels, 3, 1)

    def forward(self, x):
        x = F.silu(self.stem(x))
        for block in self.blocks:
            x = F.silu(x + .3*block(x))
        return self.head(x)


def waveform_features(t_ns, voltage_mV, settings):
    """All channels depend only on measured voltage and frozen calibration."""
    t, y = np.asarray(t_ns, float), np.asarray(voltage_mV, float)
    if t.ndim != 1 or y.shape != t.shape or not np.isfinite(y).all():
        raise ValueError('Invalid voltage samples')
    if not np.allclose(np.diff(t), 2., atol=1e-8, rtol=0):
        raise ValueError('This checkpoint is calibrated for 2 ns samples')
    cfg, template = FitConfig(**settings['fit_config']), make_template(settings)
    values = y-cfg.baseline_mV
    _, inv8 = candidate_peaks(t, values, template, cfg)
    _, inv2 = candidate_peaks(t, values, template, replace(cfg, deconvolution_sigma_ns=2.))
    features = np.array([values/500., np.gradient(values, 2.)/10., inv8/200.,
                         inv2/200., gaussian_filter1d(values, 50.)/500.], dtype=np.float32)
    # Fixed invertible compression; no test-set normalization or rate oracle.
    return np.arcsinh(features)


def supervised_targets(t, truth):
    """Soft event peaks and nearest-event charge/time labels, q >= 500 pe.

    Events within one 2 ns bin cannot have separate output coordinates; their
    ambiguity remains in scoring. Low-charge physical signals stay in input.
    """
    times, charges = np.asarray(truth['first_active_ns']), np.asarray(truth['charge_pe'])
    keep = np.isfinite(times) & (charges >= 500.)
    times, charges = times[keep], charges[keep]
    order = np.argsort(times)
    times, charges = times[order], charges[order]
    n = len(t)
    labels = np.zeros((4, n), dtype=np.float32)
    if not len(times):
        return labels
    # Nearest physical arrival; charge is not used to select nearest labels.
    j = np.searchsorted(times, t)
    a, b = np.maximum(j-1, 0), np.minimum(j, len(times)-1)
    ix = np.where(abs(t-times[a]) <= abs(t-times[b]), a, b)
    distance = t-times[ix]
    labels[0] = np.exp(-.5*(distance/2.)**2)
    labels[1] = np.log(charges[ix]/3000.)
    labels[2] = (times[ix]-t)/4.
    labels[3] = (abs(distance) <= 4.).astype(np.float32)
    return labels


def dense_prediction(model, features):
    with torch.inference_mode():
        z = model(torch.from_numpy(features[None]))[0].cpu().numpy()
    return dict(probability=1/(1+np.exp(-np.clip(z[0], -30., 30.))),
                charge_pe=3000.*np.exp(np.clip(z[1], -5., 5.)),
                time_offset_ns=4.*np.tanh(z[2]))


def decode_prediction(t, dense, threshold):
    peaks, _ = find_peaks(dense['probability'], height=threshold, distance=4)
    keep = dense['charge_pe'][peaks] >= 500.
    peaks = peaks[keep]
    times = t[peaks]+dense['time_offset_ns'][peaks]
    order = np.argsort(times)
    return dict(status='predicted', times_ns=times[order].astype(float).tolist(),
                charge_pe=dense['charge_pe'][peaks][order].astype(float).tolist(),
                candidates=len(peaks), input_policy='complete_sampled_voltage',
                threshold=float(threshold))


def load_pilot(record, settings):
    with np.load(record['waveforms']) as f:
        t, voltage = f['time_ns'], f['analog_mV']
    with np.load(record['truth']) as f:
        truth = {k: f[k] for k in ('first_active_ns', 'emission_ns', 'charge_pe')}
    metadata = json.loads(Path(record['metadata']).read_text())
    return dict(record=record, t=t, features=waveform_features(t, voltage, settings),
                targets=supervised_targets(t, truth), truth=truth, gate_ns=metadata['gate_ns'])


def loss_function(prediction, labels):
    # Train only middle 500 samples: 262-sample guards exceed 258-sample radius.
    z, labels = prediction[:, :, 262:-262], labels[:, :, 262:-262]
    event, charge, offset, mask = labels.unbind(1)
    classification = F.binary_cross_entropy_with_logits(z[:, 0], event,
                         pos_weight=torch.tensor(12., device=z.device))
    denominator = mask.sum().clamp(min=1.)
    qloss = (F.smooth_l1_loss(z[:, 1], charge, reduction='none')*mask).sum()/denominator
    tloss = (F.smooth_l1_loss(torch.tanh(z[:, 2]), offset, reduction='none')*mask).sum()/denominator
    return classification + 2.*qloss + .5*tloss, (classification, qloss, tloss)


def validation_choice(model, validation):
    dense = [dense_prediction(model, item['features']) for item in validation]
    choices = []
    for threshold in (.15, .25, .35, .45, .55, .65, .75, .85):
        groups = defaultdict(lambda: np.zeros(3))
        scores = []
        for item, prediction in zip(validation, dense):
            fitted = decode_prediction(item['t'], prediction, threshold)
            s = score(fitted, item['truth'], item['gate_ns'])
            groups[item['record']['rate_hz']] += [s['target_count'], s['reconstructed_count'], s['recovered_count']]
            scores.append(dict(label=item['record']['label'], score=s))
        # Macro charge-qualified F1: 2 successful / (targets + predictions).
        objective = float(np.mean([2*v[2]/max(v[0]+v[1], 1.) for v in groups.values()]))
        choices.append(dict(threshold=threshold, objective=objective, windows=scores))
    return max(choices, key=lambda r: (r['objective'], r['threshold']))


def train(experiment, output, steps=2400, threads=2):
    if not 1 <= threads <= 4 or steps < 1:
        raise ValueError('Invalid CPU/step budget')
    output.mkdir(parents=True, exist_ok=False)
    torch.set_num_threads(threads)
    torch.manual_seed(2610141620)
    rng = np.random.default_rng(2610141620)
    settings = json.loads((experiment/'calibration-settings.json').read_text())
    split = json.loads((experiment/'data-split.json').read_text())
    assert not {r['root_sha256'] for r in split['windows'] if r['split']=='train'} & {
        r['root_sha256'] for r in split['windows'] if r['split']=='validation'}
    for r in split['windows']:
        r['file_sha256'] = {k:digest(r[k]) for k in ('waveforms','truth','metadata')}
    save(output/'data-split.json', split)
    save(output/'analysis_config.json', dict(settings, model='Small dilated residual 1D CNN, calibrated waveform features',
        supervised_truth_fields=['first_active_ns','charge_pe'], truth_input_at_inference=[],
        rate_input_at_inference=False, candidate_count_fixed=False,
        normalization='Fixed asinh compression of voltage, derivative and inverse-filter channels',
        label_sigma_ns=2.,charge_scale_pe=3000.,min_separation_ns=8.,
        loss='BCE soft event peaks, positive weight 12; 2x masked log-charge Huber; 0.5x masked time-offset Huber',
        model_selection='Max validation macro charge-qualified F1 over declared threshold grid; no new test data accessed',
        waveform_crop_samples=1024, scored_training_core_samples=500,
        training_seed=2610141620,steps=steps,batch_size=16,cpu_threads=threads,
        environment=dict(python=platform.python_version(), torch=torch.__version__,
                         numpy=np.__version__, gpu=torch.cuda.is_available())))
    items = [load_pilot(r, settings) for r in split['windows']]
    training = [r for r in items if r['record']['split']=='train']
    validation = [r for r in items if r['record']['split']=='validation']
    by_rate = defaultdict(list)
    for item in training:
        by_rate[item['record']['rate_hz']].append(item)
    rates = sorted(by_rate)
    model = EventNet()
    optimizer = torch.optim.AdamW(model.parameters(), lr=.001, weight_decay=.0001)
    history, best, started = [], None, time.monotonic()
    print('TRAIN', len(training), 'VALIDATION', len(validation),
          'parameters',sum(p.numel() for p in model.parameters()),flush=True)
    for step in range(1, steps+1):
        x, labels = [], []
        for _ in range(16):
            group = by_rate[rates[int(rng.integers(len(rates)))]]
            item = group[int(rng.integers(len(group)))]
            start = int(rng.integers(2250, 6739))  # Core remains in the central evaluation interval.
            x.append(item['features'][:,start:start+1024])
            labels.append(item['targets'][:,start:start+1024])
        model.train()
        optimizer.zero_grad(set_to_none=True)
        loss, terms = loss_function(model(torch.from_numpy(np.stack(x))), torch.from_numpy(np.stack(labels)))
        if not torch.isfinite(loss):
            raise FloatingPointError('Nonfinite training loss')
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.)
        optimizer.step()
        if step == 1 or step % 200 == 0 or step == steps:
            model.eval()
            selected = validation_choice(model, validation)
            row = dict(step=step,loss=float(loss.detach()),
                loss_terms=[float(v.detach()) for v in terms], validation_objective=selected['objective'],
                threshold=selected['threshold'],elapsed_s=time.monotonic()-started)
            history.append(row)
            if best is None or selected['objective'] > best['objective']:
                best = dict(selected, step=step)
                torch.save(model.state_dict(), output/'model.pt')
                save(output/'validation-selection.json', best)
            save(output/'training-history.json', history)
            print('STEP', step, 'LOSS',round(row['loss'],4),'VAL',round(row['validation_objective'],4),
                  'THRESHOLD',row['threshold'],'seconds',round(row['elapsed_s'],1),flush=True)
    checkpoint = dict(created=datetime.now().astimezone().isoformat(),state='frozen_before_test',
        weights_sha256=digest(output/'model.pt'),config_sha256=digest(output/'analysis_config.json'),
        data_split_sha256=digest(output/'data-split.json'),threshold=best['threshold'],step=best['step'],
        validation_objective=best['objective'],training_seconds=time.monotonic()-started,
        parameters=sum(p.numel() for p in model.parameters()),
        model_source_sha256=digest(Path(__file__)))
    save(output/'frozen-model.json', checkpoint)
    print('MODEL FROZEN', checkpoint,flush=True)


def evaluate(experiment, output, data):
    """Run the frozen neural model and unchanged blind fit on the same new tests."""
    frozen = json.loads((output/'frozen-model.json').read_text())
    assert frozen['weights_sha256']==digest(output/'model.pt')
    assert frozen['config_sha256']==digest(output/'analysis_config.json')
    assert frozen['model_source_sha256']==digest(Path(__file__))
    settings = json.loads((output/'analysis_config.json').read_text())
    torch.set_num_threads(1)
    model = EventNet()
    model.load_state_dict(torch.load(output/'model.pt', map_location='cpu',weights_only=True))
    model.eval()
    for sub in ('neural','fit'):
        (output/sub).mkdir(exist_ok=True)
    manifest = json.loads((experiment/'manifest.json').read_text())
    rows = []
    training_seeds = {r['seed'] for r in json.loads((output/'data-split.json').read_text())['windows']}
    for run in manifest['runs']:
        if run['status']!='completed' or run['seed'] in training_seeds:
            raise ValueError('Incomplete or overlapping test ROOT')
        stem = Path(run['root']).stem
        result_path=output/'fit'/f'{stem}.json'
        if result_path.exists():
            rows.append(json.loads(result_path.read_text()))
            continue
        with np.load(data/'waveforms'/f'{stem}.npz') as f:
            t, voltage = f['time_ns'], f['analog_mV']
        start = time.monotonic()
        dense = dense_prediction(model, waveform_features(t, voltage, settings))
        neural = decode_prediction(t, dense, frozen['threshold'])
        nn_seconds = time.monotonic()-start
        start = time.monotonic()
        fitted = fit_waveform(t, voltage, make_template(settings), FitConfig(**settings['fit_config']))
        fit_seconds = time.monotonic()-start
        np.savez_compressed(output/'fit'/f'{stem}.npz',model_mV=fitted.pop('model_mV'))
        np.savez_compressed(output/'neural'/f'{stem}.npz', **dense)
        # Open test truth only after both algorithms have completed blind inference.
        with np.load(data/'truth'/f'{stem}.npz') as f:
            truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe')}
        meta=json.loads((data/'metadata'/f'{stem}.json').read_text())
        row=dict(label=run['label'],root=run['root'],root_sha256=run['root_sha256'],
            seed=run['seed'],rate_hz=run['rate_hz'],gate_ns=meta['gate_ns'],
            neural_fit=neural,neural_score=score(neural,truth,meta['gate_ns']),
            fit=fitted,score=score(fitted,truth,meta['gate_ns']),
            neural_seconds=nn_seconds,fit_seconds=fit_seconds,
            waveform_path=str(data/'waveforms'/f'{stem}.npz'),
            truth_path=str(data/'truth'/f'{stem}.npz'),
            waveform_sha256=digest(data/'waveforms'/f'{stem}.npz'),
            truth_sha256=digest(data/'truth'/f'{stem}.npz'))
        save(result_path,row)
        rows.append(row)
        save(output/'test-progress.json',dict(completed=len(rows),planned=len(manifest['planned_runs'])))
        print('TEST',run['label'],'NN',row['neural_score']['recovered_count'],
              'FIT',row['score']['recovered_count'],'/',row['score']['target_count'],flush=True)
    assert len(rows)==len(manifest['planned_runs'])
    save(output/'summary.json',dict(state='completed',windows=rows,
        fit_rates=summarize(rows),neural_rates=summarize([dict(r,fit=r['neural_fit'],score=r['neural_score']) for r in rows]),
        frozen=frozen,input_policy='complete_sampled_voltage',
        caveats=['Small exploratory pilot; new ROOT seeds independent of training/validation',
                 'Supervised simulation training, no measured detector domain test',
                 'Actual shared-cell charge q>=500 pe, not reconstructed deposited-energy spectrum',
                 'No test threshold/weight tuning; 100 MHz has no pilot validation ROOT',
                 '8 ns candidate separation and 2 ns output grid limit very close pairs',
                 'Fixed pure-alpha family, fixed sampling/noise/shape; no background or rate oracle',
                 'Adequate acquisition dynamic range; no clipping or quantization considered']))
    print('COMPLETE',output/'summary.json',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('train','evaluate'))
    p.add_argument('--experiment',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--data',type=Path)
    p.add_argument('--steps',type=int,default=2400)
    p.add_argument('--threads',type=int,default=2)
    args=p.parse_args()
    if args.stage=='train':
        train(args.experiment,args.output,args.steps,args.threads)
    else:
        if args.data is None:
            p.error('--data is required for evaluate')
        evaluate(args.experiment,args.output,args.data)


if __name__=='__main__':
    main()
