"""Validate frozen results; omit only the measured wall time in exact replay.

The archived validator compares a whole WLS result dictionary, which includes
``refinement_seconds``. The numerical reconstruction must replay exactly, but a
new execution's wall time cannot. The adapter exempts only that profiling field from exact replay, evaluates
the same physical validation sum in 128-row blocks, and records progress with
garbage collection between windows. It preserves the archived file and all
numerical assertions, and also checks original ROOT sample branches.
"""
from __future__ import annotations

import difflib
import gc
import json
import os
from pathlib import Path
import sys
import resource

PROJECT = Path(__file__).resolve().parents[3]
EXPERIMENT = Path(__file__).resolve().parent
RECORD = PROJECT / 'jupyter/run_records/alpha_component_ml_20261021'
OUTPUT = PROJECT / 'jupyter/artifacts/alpha_component_ml_20261021'
os.environ['GAGG_COMPONENT_PROJECT'] = str(PROJECT)
sys.path.insert(0, str(EXPERIMENT / 'inference-source'))

import awkward as ak
import numpy as np
from pileup.analysis import alpha_component_validate as archived
from pileup.analysis.alpha_rate_campaign import digest, save


def numerical_fit(fit):
    """Only this explicit profiling field is exempt from exact replay."""
    return {key: value for key, value in fit.items() if key != 'refinement_seconds'}


def progress(stage, identity):
    row=dict(stage=stage,identity=str(identity),max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    save(RECORD/'validation-progress.json',row)
    print('VALIDATION_PROGRESS',stage,identity,'RSS_KiB',row['max_rss_kib'],flush=True)
    gc.collect()


def reconstructed_model(time, times, amplitudes, baseline, template):
    """Same per-row physical sum, with bounded temporary interpolation arrays."""
    output=np.empty_like(time,dtype=float)
    for start in range(0,len(time),128):
        output[start:start+128]=baseline+template(time[start:start+128,None]-times[None,:])@amplitudes
    return output


def main():
    source_path = Path(archived.__file__)
    source = source_path.read_text()
    before = "assert compatible(fit)==first['methods'][name]['fit'],name"
    after = "assert numerical_fit(compatible(fit))==numerical_fit(first['methods'][name]['fit']),name"
    assert source.count(before) == 1
    adapted = source.replace(before, after)
    dense_model="actual=settings['fit_config']['baseline_mV']+template(t[:,None]-tt[None,:])@aa"
    block_model="actual=reconstructed_model(t,tt,aa,settings['fit_config']['baseline_mV'],template)"
    assert adapted.count(dense_model)==1
    adapted=adapted.replace(dense_model,block_model)
    for before,after in [
        ("root_checks.append(dict(root=r['root'],seed=r['seed'],actual_runinfo_verified=True))",
         "root_checks.append(dict(root=r['root'],seed=r['seed'],actual_runinfo_verified=True))\n        progress('RunInfo',r['root'])"),
        ("stem=Path(r['root']).stem", "stem=Path(r['root']).stem\n        progress('scores-and-models',stem)")]:
        assert adapted.count(before)==1
        adapted=adapted.replace(before,after)
    delta = ''.join(difflib.unified_diff(source.splitlines(keepends=True),
                                       adapted.splitlines(keepends=True),
                                       fromfile='frozen-validator', tofile='timer-exempt-validator'))
    (RECORD / 'validation-timer-adapter.diff').write_text(delta)
    namespace = dict(vars(archived), numerical_fit=numerical_fit,
                     reconstructed_model=reconstructed_model,progress=progress,
                     __name__='component_validation_adapter')
    exec(compile(adapted, str(source_path) + ':timer-exemption', 'exec'), namespace)
    namespace['validate']()

    summary = json.loads((OUTPUT / 'test-summary.json').read_text())
    declared_gallery = json.loads((RECORD / 'gallery-choice-before-test-inference.json').read_text())
    actual_gallery = json.loads((OUTPUT / 'gallery-selection.json').read_text())
    # The record's field may contain only the selected model and rationale;
    # compare the selection itself, rather than incidental metadata wording.
    assert declared_gallery['model'] == actual_gallery['model']
    replay = json.loads((RECORD / 'cold-replay/predictions.json').read_text())
    cold = next(row for row in summary['windows'] if row['waveform'] == replay['waveform'])
    assert cold['waveform_sha256'] == replay['waveform_sha256']
    timer_differences = []
    for name, fit in replay['fits'].items():
        expected = cold['methods'][name]['fit']
        assert numerical_fit(fit) == numerical_fit(expected), name
        differing = [key for key in fit if fit[key] != expected.get(key)]
        assert set(differing) <= {'refinement_seconds'}, (name, differing)
        timer_differences.extend((name, key) for key in differing)

    rows = []
    for row in summary['windows']:
        with archived.open_root(PROJECT / row['root']) as f:
            tree = f['SiPMEvent']
            sample_time = ak.to_numpy(tree['sample_time_ns'].array()[0])
            analog = ak.to_numpy(tree['analog_mV'].array()[0])
        with np.load(PROJECT / row['waveform']) as f:
            assert np.array_equal(sample_time, f['time_ns'])
            assert np.array_equal(analog, f['analog_mV'])
        rows.append(dict(root=row['root'],samples=len(sample_time),
                         native_voltage_and_time_exact=True))
    preserved_galleries = []
    for previous in (RECORD / 'partial-gallery-pass1').glob('*/components.npz'):
        current = OUTPUT / 'gallery' / previous.parent.name / previous.name
        with np.load(previous) as a, np.load(current) as b:
            assert a.files == b.files
            assert all(np.array_equal(a[key], b[key]) for key in a.files)
        preserved_galleries.append(previous.parent.name)
    result = dict(status='passed', actual_cold_cli_methods=len(replay['fits']),
                  cold_replay_non_numeric_differences=timer_differences,
                  replay_exempt_profiling_fields=['refinement_seconds'],
                  frozen_validator_sha256=digest(source_path),
                  validator_adapter_sha256=digest(__file__),
                  initial_validation_exit_code=143,
                  validation_reconstruction_rows_per_block=128,
                  inference_source_unchanged=True,
                  rerendered_gallery_arrays_exact_to_completed_first_pass=preserved_galleries,
                  original_root_sample_checks=rows)
    save(RECORD / 'additional-validation.json', result)
    print('ADDITIONAL_VALIDATION', len(rows), 'original ROOT waveforms,',
          len(replay['fits']), 'cold CLI numerical fits exact', flush=True)


if __name__ == '__main__':
    main()
