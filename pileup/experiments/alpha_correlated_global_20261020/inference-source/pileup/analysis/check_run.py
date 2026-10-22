#!/usr/bin/env python3
"""Inspect pure-alpha simulation output; this is not a reconstruction algorithm."""
import argparse
import json
from pathlib import Path
import numpy as np
import uproot
from uproot.source.file import MemmapSource


def check(path, min_photons=20):
    with uproot.open(path, handler=MemmapSource) as root:
        info = {name: root['RunInfo'][name].array(library='np').tolist()
                for name in root['RunInfo'].keys()}
        if root['RunInfo'].num_entries != 1:
            raise ValueError('Expected one RunInfo row.')
        scalar = lambda key: info[key][0]
        window = scalar('t_length_ps')
        start, end = scalar('analysis_start_ps'), scalar('analysis_end_ps')
        if not 0 <= start < end <= window:
            raise ValueError('Invalid analysis gate.')
        if not scalar('pileup_enabled') or not scalar('alpha_p11b_enabled'):
            raise ValueError('Expected alpha pile-up mode.')
        if scalar('brems_enabled') or scalar('c12_capture_enabled'):
            raise ValueError('Photon backgrounds are enabled.')
        events = root['WaveformEvent'].arrays(
            ['event_id', 'n_primary', 'n_optical_photons_arrived_total', 'edep_total_MeV'], library='np')
        primaries = root['PrimaryPhoton'].arrays(
            ['event_id', 'primary_id', 'primary_particle_name', 'primary_response_name',
             'primary_energy_MeV', 'primary_time_ps', 'edep_total_MeV', 'n_optical_photons_arrived'], library='np')
        n_events = len(events['event_id'])
        if n_events == 0:
            raise ValueError('No waveform windows were recorded.')
        if n_events != scalar('n_events') or not np.array_equal(events['event_id'], np.arange(n_events)):
            raise ValueError('Incomplete or unexpected event IDs.')
        names = set(primaries['primary_particle_name'])
        if not names <= {'alpha'}:
            raise ValueError(f'Non-alpha primaries: {names}')
        p_event = primaries['event_id'].astype(int)
        if len(p_event) and (p_event.min() < 0 or p_event.max() >= n_events):
            raise ValueError('Primary event IDs are outside the recorded events.')
        counts = np.bincount(p_event, minlength=n_events)
        arrivals = np.bincount(p_event, weights=primaries['n_optical_photons_arrived'], minlength=n_events)
        edep = np.bincount(p_event, weights=primaries['edep_total_MeV'], minlength=n_events)
        if not np.array_equal(counts, events['n_primary']):
            raise ValueError('Primary multiplicities disagree between trees.')
        if not np.array_equal(arrivals, events['n_optical_photons_arrived_total']):
            raise ValueError('Photon counts disagree between trees.')
        if not np.allclose(edep, events['edep_total_MeV'], rtol=3e-4, atol=2e-5):
            raise ValueError('Deposited energies disagree beyond float-accumulation tolerance.')
        times = primaries['primary_time_ps']
        if np.any(times >= window):
            raise ValueError('Emission time outside source window.')
        for event_id in np.unique(p_event):
            mask = p_event == event_id
            if np.any(np.diff(times[mask].astype(np.int64)) < 0):
                raise ValueError('Primary times are not sorted.')
        central = (times >= start) & (times < end)
        signal = central & (primaries['edep_total_MeV'] > 0) & (primaries['n_optical_photons_arrived'] >= min_photons)
        duration = n_events * window * 1e-12
        central_duration = n_events * (end - start) * 1e-12
        report = {'file': str(path), 'status': 'validated', 'n_windows': n_events,
                  'primary_particles': int(counts.sum()), 'primary_types': sorted(names),
                  'empty_windows': int((counts == 0).sum()), 'window_us': window*1e-6,
                  'analysis_gate_us': [start*1e-6, end*1e-6],
                  'configured_injection_rate_hz': scalar('alpha_p11b_rate_hz'),
                  'observed_injection_rate_hz': float(counts.sum()/duration),
                  'central_source_count': int(central.sum()),
                  'central_signal_count': int(signal.sum()),
                  'central_signal_rate_hz': float(signal.sum()/central_duration),
                  'signal_definition': f'Edep > 0 and arrived photons >= {min_photons}; central source-emission time',
                  'mean_primaries_per_window': float(counts.mean()),
                  'expected_mean_primaries_per_window': scalar('alpha_p11b_rate_hz')*window*1e-12,
                  'random_seed': scalar('random_seed'),
                  'note': 'Data-integrity and source-rate check only; no pile-up reconstruction has been performed.'}
        manifest_path = path.parent/'manifest.json'
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text())
            if n_events != manifest['events'] or scalar('random_seed') != manifest['seed']:
                raise ValueError('ROOT event count or seed differs from the run manifest.')
            if scalar('alpha_p11b_rate_hz') != manifest['point']['rate_hz']:
                raise ValueError('ROOT rate differs from the run manifest.')
            config = manifest['config']
            if not np.isclose(window, config['window_us']*1e6, rtol=0, atol=0.5):
                raise ValueError('ROOT full window differs from the run manifest.')
            expected_gate = [config['pre_window_us']*1e6,
                             (config['window_us']-config['post_window_us'])*1e6]
            if not np.allclose([start,end], expected_gate, rtol=0, atol=0.5):
                raise ValueError('ROOT analysis gate differs from the run manifest.')
        return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root_file', type=Path)
    parser.add_argument('--min-photons', type=int, default=20)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.min_photons < 1:
        parser.error('--min-photons must be positive.')
    report = check(args.root_file.resolve(), args.min_photons)
    text = json.dumps(report, ensure_ascii=False, indent=2) + '\n'
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(text)
    print(text, end='')


if __name__ == '__main__':
    main()
