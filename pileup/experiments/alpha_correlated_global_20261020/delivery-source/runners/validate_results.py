"""Independent consistency audit of fixed inference, exports and saved scores.

Not a new evaluation protocol or a source of inference parameters.
"""
from datetime import datetime
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile
from unittest.mock import patch

import awkward as ak
import numpy as np
import uproot

project = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project))
from pileup.analysis.alpha_correlated_evaluate import frozen
from pileup.analysis.alpha_correlated_fit import CorrelatedConfig, InnovationsMetric, ResidualCalibration, fit_correlated, projected
from pileup.analysis.alpha_correlated_study import OUTPUT, EXPERIMENT, RECORD, OLD
from pileup.analysis.alpha_full_signal_study import make_template
from pileup.analysis.alpha_pileup_fit import FitConfig, fit_waveform, score
from pileup.analysis.alpha_rate_campaign import digest, save
from pileup.analysis.alpha_rate_data import open_root


def main():
    f = frozen(); old = json.loads((OLD/'freeze.json').read_text())
    io = json.loads((EXPERIMENT/'io-recovery.json').read_text())
    for name, sha in io['source_sha256'].items():
        assert digest(project/name) == sha, name
    for name, sha in io['snapshot_sha256'].items():
        assert digest(project/name) == sha, name
    budget = json.loads((OUTPUT/'optimizer-budget-controls.json').read_text())
    assert len(budget['rows']) == 108
    assert digest(project/'pileup/analysis/alpha_optimizer_budget.py') == budget['protocol']['source_sha256']
    assert digest(OUTPUT/'calibration/covariance.json') == budget['protocol']['covariance_sha256']
    assert digest(OUTPUT/'baseline-settings.json') == budget['protocol']['baseline_settings_sha256']
    for row in budget['rows']:
        assert all(v['profile_evaluations'] == 500 for v in row['results'].values())
        for fitted in row['results'].values():
            correct = np.all(abs(np.asarray(fitted['time_errors_ns'])) <= 10) and np.all(abs(np.asarray(fitted['relative_charge_errors'])) <= .2)
            assert bool(correct) == fitted['all_components_complete']
    for name, sha in old['frozen_files'].items():
        assert digest(project/name) == sha, name
    print('FROZEN_HASHES', len(f['frozen_files']), len(old['frozen_files']), flush=True)
    manifest = json.loads((EXPERIMENT/'manifest.json').read_text())
    test = json.loads((OUTPUT/'test-summary.json').read_text())
    assert len(manifest['runs']) == len(test['windows']) == 23
    assert manifest['preflight_validated']
    before = json.loads((OUTPUT/'input-records.json').read_text())
    development = before['calibration_single']+before['development_single']+before['development_pileup']
    split = json.loads((OUTPUT/'gru/data-split.json').read_text())['windows']
    exposed = {r['root'] for r in development+split}
    old_tests = json.loads((OLD/'test-summary.json').read_text())['windows']
    exposed |= {r['root'] for r in old_tests}
    # This checks identities against all listed old neural windows, including
    # previously exposed tests; no fitting data are silently reclassified.
    assert not exposed & {r['root'] for r in manifest['runs']}
    old_seeds = {r['seed'] for r in development+split+old_tests if 'seed' in r}
    assert not old_seeds & {r['seed'] for r in manifest['runs']}
    starts = sorted(datetime.fromisoformat(r['started']) for r in manifest['runs'])
    assert min(starts) > datetime.fromisoformat(f['frozen_at'])
    assert min((b-a).total_seconds() for a, b in zip(starts[:-1], starts[1:])) >= 1.
    for seed_key in ('seed', 'sipm_seed', 'readout_seed'):
        assert len({r[seed_key] for r in manifest['runs']}) == 23
    settings = json.loads((OUTPUT/'baseline-settings.json').read_text())
    base = FitConfig(**settings['fit_config']); template = make_template(settings)
    cal = ResidualCalibration(**json.loads((OUTPUT/'calibration/covariance.json').read_text())['parameters'])
    for case in (0,17,35):
        with np.load(OUTPUT/'optimizer-budget'/f'case-{case:02d}.npz') as z:
            bt,by,true_t,true_q=(z[k].copy() for k in ('time_ns','analog_mV','true_onsets_ns','true_charge_pe'))
        initial=fit_waveform(bt,by,template,base)
        metric=InnovationsMetric(initial['model_mV']-10,cal)
        for row in (r for r in budget['rows'] if r['case']==case):
            for fitted in row['results'].values():
                cost,_,amplitude=projected(bt,by-10,fitted['times_ns'],template,metric)
                np.testing.assert_allclose(cost,fitted['cost'],rtol=1e-11,atol=1e-8)
                q=amplitude/(base.area_per_pe_mV_ns*template.unit_area_peak)
                np.testing.assert_allclose(q,fitted['charge_pe'],rtol=1e-11,atol=1e-8)
                np.testing.assert_allclose(np.asarray(fitted['times_ns'])-true_t,fitted['time_errors_ns'],atol=1e-11)
                np.testing.assert_allclose(q/true_q-1,fitted['relative_charge_errors'],atol=1e-11)
    rows = {r['root']: r for r in test['windows']}
    charge_errors = []; model_errors = []; audit_rows = []; failed_initial_models = 0
    for run in manifest['runs']:
        assert run['status'] == 'completed' and run['returncode'] == 0
        assert re.fullmatch(r'g4/data/\d{8}_\d{2}h\d{2}m\d{2}s\.root', run['root'])
        assert (project/run['root']).stat().st_size == run['root_bytes']
        assert digest(project/run['macro']) == run['macro_sha256']
        meta = json.loads((project/run['export_metadata']).read_text())
        assert meta['root'] == run['root'] and meta['root_sha256'] == run['root_sha256']
        if 'scoring_field_recovery' in run:
            recovery = run['scoring_field_recovery']
            assert recovery['one_to_one_complete_coverage'] and recovery['full_stored_common_fields_equal']
            assert recovery['global_photons'] == recovery['primary_photons'] == meta['photons']
            assert not recovery['inference_or_model_inputs_used']
        assert meta['seed'] == run['seed'] and meta['rate_hz'] == run['rate_hz']
        assert meta['source_rates_hz']['brems'] == meta['source_rates_hz']['c12_capture'] == 0
        assert meta['macro_contents'] == (project/run['macro']).read_text()
        # Small RunInfo branches use ordinary buffered reads on the shared
        # original, avoiding a full ROOT copy for this per-file seed audit.
        with uproot.open(project/run['root'], handler=uproot.source.file.MultithreadedFileSource,
                         num_workers=1, decompression_executor=uproot.TrivialExecutor(),
                         interpretation_executor=uproot.TrivialExecutor()) as actual_root:
            actual_seeds = actual_root['RunInfo'].arrays(
                ['random_seed','sipm_avalanche_base_seed','sipm_readout_base_seed'],library='np')
            recorded_seeds = {k:int(v[0]) for k,v in actual_seeds.items()}
        assert recorded_seeds == dict(random_seed=run['seed'], sipm_avalanche_base_seed=run['sipm_seed'],
                                      sipm_readout_base_seed=run['readout_seed'])
        assert meta['gate_ns'] == [5000., 15000.] and meta['duration_ns'] == 20000.
        assert meta['config']['readout']['sample_period_ns'] == 2.
        assert meta['config']['sensor']['dark_rate_hz'] == 2000000.
        assert meta['config']['sensor']['max_candidates'] == 30000000
        row = rows[run['root']]
        assert row['inference_truth_inputs'] == []
        assert digest(project/row['waveform']) == row['waveform_sha256']
        assert digest(project/row['truth']) == row['truth_sha256']
        with np.load(project/row['waveform']) as z:
            t, y = z['time_ns'].copy(), z['analog_mV'].copy()
        with np.load(project/row['truth']) as z:
            truth = {k: z[k].copy() for k in ('first_active_ns', 'emission_ns', 'charge_pe')}
            charge = z['avalanche_charge_C']; owners = z['avalanche_owner']
            q0 = meta['config']['sensor']['gain_electrons']*1.602176634e-19
            light = owners >= 0
            # Independent aggregation from saved ancestry, retaining dark chains.
            reconstructed = np.bincount(owners[light], weights=charge[light]/q0, minlength=len(truth['charge_pe']))
            np.testing.assert_allclose(reconstructed, truth['charge_pe'], rtol=2e-11, atol=1e-7)
            dark = float(charge[~light].sum()/q0)
            np.testing.assert_allclose(dark, meta['dark_charge_pe'], rtol=1e-12, atol=1e-7)
            np.testing.assert_allclose(reconstructed.sum()+dark, meta['charge_pe'], rtol=2e-11, atol=1e-6)
            charge_errors.append(float(np.max(abs(reconstructed-truth['charge_pe']), initial=0)))
        for name, record in row['methods'].items():
            # JSON stores matched (truth, prediction) tuples as two-element
            # lists. Normalize representation before checking every value.
            rescored = json.loads(json.dumps(score(record['fit'], truth, row['gate_ns'])))
            assert rescored == record['score'], (run['root'], name)
        assert len(row['methods']) == 7
        stem = Path(run['root']).stem
        with np.load(OUTPUT/'pilots'/(stem+'.npz')) as z: pilot_mean = z['model_mV'].copy()-base.baseline_mV
        for name in ('diagonal', 'correlated', 'global'):
            fit = row['methods'][name]['fit']
            if fit['status'] != 'fitted':
                assert fit['correlated_status'] == 'initial_fit_failed'
                assert fit['times_ns'] == fit['charge_pe'] == []
                with np.load(OUTPUT/'test'/name/(stem+'.npz')) as z:
                    np.testing.assert_array_equal(z['model_mV'], pilot_mean+base.baseline_mV)
                failed_initial_models += 1
                continue
            all_fit = fit['all_components']
            model = base.baseline_mV+template(t[:, None]-np.asarray(all_fit['times_ns']))@np.asarray(all_fit['peak_mV'])
            with np.load(OUTPUT/'test'/name/(stem+'.npz')) as z:
                difference = float(np.max(abs(model-z['model_mV'])))
                np.testing.assert_allclose(model, z['model_mV'], rtol=1e-12, atol=2e-8)
            model_errors.append(difference)
            metric = InnovationsMetric(pilot_mean, cal, diagonal=name == 'diagonal')
            r = metric.apply(y-model)
            np.testing.assert_allclose(r@r, fit['whitened_cost'], rtol=1e-10, atol=1e-5)
            np.testing.assert_allclose(metric.logdet, fit['logdet_covariance'], rtol=1e-13)
            assert fit['fit_truth_inputs'] == []
            assert fit['global_calls_used'] <= (12 if name == 'global' else 0)
        audit_rows.append(dict(root=run['root'], rate_hz=run['rate_hz'], seed=run['seed'],
                               run_info_base_seeds=recorded_seeds,
                               rescored_methods=7, physical_models_recomputed=sum(row['methods'][m]['fit']['status']=='fitted' for m in ('diagonal','correlated','global')), charge_balanced=True))
    print('SCORES_MODELS_ANCESTRY', len(audit_rows), flush=True)
    # Read actual ROOTs again, through ordinary local byte copies. The nominal
    # config stored in these RunInfo records, rather than today's macro, is checked.
    root_checks = []
    for rate in (100000, 10000000, 100000000):
        run = min((r for r in manifest['runs'] if r['rate_hz'] == rate), key=lambda r: r['seed'])
        row = rows[run['root']]
        with tempfile.TemporaryDirectory(prefix='gagg-correlated-root-audit-', dir='/tmp') as temporary:
            local = Path(temporary)/Path(run['root']).name
            shutil.copyfile(project/run['root'], local)
            assert digest(local) == digest(project/run['root']) == run['root_sha256']
            with open_root(local) as root:
                assert all(root[name].num_entries == 1 for name in ('RunInfo', 'WaveformEvent', 'SiPMEvent'))
                info = ak.to_list(root['RunInfo'].arrays())[0]
                assert info['random_seed'] == run['seed'] and info['alpha_p11b_rate_hz'] == rate
                assert info['macro_contents'] == (project/run['macro']).read_text()
                sipm = root['SiPMEvent'].arrays(['sample_time_ns', 'analog_mV', 'signal_mV'])
                with np.load(project/row['waveform']) as z:
                    assert np.array_equal(ak.to_numpy(sipm['sample_time_ns'][0]), z['time_ns'])
                    assert np.array_equal(ak.to_numpy(sipm['analog_mV'][0]), z['analog_mV'])
                with np.load(project/row['truth']) as z:
                    assert np.array_equal(ak.to_numpy(sipm['signal_mV'][0]), z['signal_mV'])
        root_checks.append(dict(root=run['root'], hash_and_runinfo=True, voltage_and_signal_bitwise_equal=True))
        print('ROOT_RECHECK', run['root'], flush=True)
    # A cold fit from measured arrays must work with all numpy file loading
    # denied. This is stronger than the public function having no truth argument.
    chosen = min((r for r in manifest['runs'] if r['rate_hz'] == 1000000), key=lambda r: r['seed'])
    row = rows[chosen['root']]
    with np.load(project/row['waveform']) as z: t, y = z['time_ns'].copy(), z['analog_mV'].copy()
    with patch('numpy.load', side_effect=AssertionError('Inference attempted file input')):
        pilot = fit_waveform(t, y, template, base)
        replays = {}
        for name, chosen_config in f['selected'].items():
            fit = fit_correlated(t, y, template, base, cal, CorrelatedConfig(**chosen_config['config']), initial=pilot)
            reference = row['methods'][name]['fit']
            assert np.array_equal(fit['times_ns'], reference['times_ns'])
            assert np.array_equal(fit['charge_pe'], reference['charge_pe'])
            replays[name] = dict(times_and_charge_bitwise_equal=True, numpy_file_loading_denied=True)
    result = dict(checked_at=datetime.now().astimezone().isoformat(),
        current_frozen_files=len(f['frozen_files']), old_frozen_files_unchanged=len(old['frozen_files']),
        io_source_and_snapshot_sha256_checked=True,
        supplemental_search_runs=108, supplemental_equal_profile_budget=500,
        supplemental_search_predictions_independently_recomputed=18,
        run_info_three_base_seeds_directly_checked_roots=23,
        scoring_field_recovery_roots=[r['root'] for r in manifest['runs'] if 'scoring_field_recovery' in r],
        test_roots=23, saved_scores_recomputed=161, saved_weighted_models_recomputed=len(model_errors),
        initial_failure_models_checked=failed_initial_models,
        starts_after_freeze=True, minimum_start_interval_s=min((b-a).total_seconds() for a,b in zip(starts[:-1], starts[1:])),
        known_training_calibration_development_and_exposed_roots_disjoint=True,
        max_ancestry_charge_error_pe=max(charge_errors), max_model_difference_mV=max(model_errors),
        rows=audit_rows, direct_root_rechecks=root_checks, blind_cold_replay=dict(root=chosen['root'], methods=replays))
    save(RECORD/'scientific-validation.json', result)
    print('SCIENTIFIC_VALIDATION_PASSED', flush=True)


if __name__ == '__main__':
    main()
