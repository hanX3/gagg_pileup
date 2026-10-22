"""Blind full-pipeline prior comparison, followed by development-only scoring."""
from collections import defaultdict
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
import argparse
import json
import time
import numpy as np
from .alpha_gamma_full_prior_fit import FullPriorConfig, fit_full_prior
from .alpha_gamma_prior_ablation import decode_cached
from .alpha_gamma_prior_study import digest, save, measure
from .alpha_gamma_pileup_fit import score_mixed


CONFIGS = {
    'class_only': FullPriorConfig(),
    'replay_none': FullPriorConfig(charge_strength=0., class_strength=0.),
    'replay_both': FullPriorConfig(charge_strength=10., class_strength=10.),
}


def fit_case(work, stem, controls=False):
    work = Path(work)
    case = json.loads((work/'cases'/f'{stem}.json').read_text())
    for name, item in case['inputs'].items():
        assert digest(work/item['path']) == item['sha256'], name
    assert digest(work/'calibration.json') == case['calibration_sha256']
    wave = dict(np.load(work/case['inputs']['waveform']['path']))
    dense = dict(np.load(work/case['inputs']['dense']['path']))
    baseline = json.loads((work/case['inputs']['baseline']['path']).read_text())
    neural = json.loads((work/case['inputs']['neural']['path']).read_text())
    for cached in (baseline, neural):
        assert cached['waveform_sha256'] == case['inputs']['waveform']['sha256']
        assert cached['calibration_sha256'] == case['calibration_sha256']
    baseline['fit']['model_mV'] = np.load(work/case['inputs']['baseline_model']['path'])['model_mV']
    t, voltage = wave['time_ns'], wave['analog_mV']
    assert decode_cached(t, dense, neural['fit']['thresholds']) == neural['fit']['all_components']
    calibration = json.loads((work/'calibration.json').read_text())
    results = []
    for method in CONFIGS if controls else ['class_only']:
        start = time.perf_counter()
        fitted = fit_full_prior(t, voltage, calibration, neural['fit'], dense, baseline['fit'], CONFIGS[method])
        seconds = time.perf_counter()-start
        model = fitted.pop('model_mV')
        destination = work/'outputs/fits'/method/f'{stem}.json'
        if destination.exists():
            raise FileExistsError(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(destination.with_suffix('.npz'), time_ns=t, model_mV=model)
        save(destination, dict(stem=stem, method=method, fit=fitted, seconds=seconds,
            calibration_sha256=case['calibration_sha256'], input_sha256={k:v['sha256'] for k,v in case['inputs'].items()},
            algorithm_sha256={n:digest(Path(__file__).with_name(n)) for n in
                              ('alpha_gamma_full_prior_fit.py','alpha_gamma_full_prior_study.py')},
            fit_truth_inputs=[], truth_file_access=False))
        results.append(dict(method=method, seconds=seconds, reported=len(fitted['times_ns']),
                            additions=len(fitted['residual_additions']), pair_trials=fitted['pair_trials']))
    print(json.dumps(dict(stem=stem, results=results)), flush=True)


def analyze(project, folder, output):
    project, folder, output = map(Path, (project, folder, output))
    protocol = json.loads((folder/'protocol.json').read_text())
    rows = protocol['scoring_rows']
    assert len(rows) == 20
    assert all((output/'fits/class_only'/(Path(r['root']).stem+'.json')).exists() for r in rows)
    old_methods = {'baseline':'baseline', 'hybrid_none':'hybrid_unet_0', 'hybrid_both':'hybrid_unet_10'}
    methods = list(old_methods)+['class_only']
    by_method = defaultdict(list); strata = defaultdict(list); events = []; numerical = []
    for row in rows:
        stem = Path(row['root']).stem
        strata[(row['gamma_model'],row['rate_hz'])].append(stem)
        assert digest(project/row['truth']) == row['input_sha256']['truth']
        assert digest(project/row['metadata']) == row['input_sha256']['metadata']
        with np.load(project/row['truth']) as data:
            truth = {key:data[key] for key in ('first_active_ns','emission_ns','particle_name','charge_pe')}
        gate = json.loads((project/row['metadata']).read_text())['gate_ns']
        wave = dict(np.load(project/row['waveform']))
        for method in methods:
            if method == 'class_only':
                record = json.loads((output/'fits'/method/f'{stem}.json').read_text())
                assert record['fit']['config'] == protocol['fit_config']
                history = record['fit']['weighted_cost_history']
                numerical.append(dict(stem=stem, cost_nonincreasing=bool(np.all(np.diff(history) <= 1e-6)),
                    objective_history=history, pair_trials=record['fit']['pair_trials'],
                    optimizer_status_counts=record['fit']['optimizer_status_counts'],
                    additions=len(record['fit']['residual_additions'])))
            else:
                record = json.loads((project/protocol['parent_artifacts']/'results/validation'/old_methods[method]/f'{stem}.json').read_text())
                expected = protocol['reference_sha256'][method][stem]
                assert digest(project/protocol['parent_artifacts']/'results/validation'/old_methods[method]/f'{stem}.json') == expected
            fitted = record['fit']; score = score_mixed(fitted, truth, gate)
            by_method[method].append(dict(row=row, fit=fitted, truth=truth, score=score, gate_ns=gate))
            save(output/'scores'/method/f'{stem}.json', dict(row=row, method=method, score=score, gate_ns=gate))
        physical = np.flatnonzero(np.isfinite(truth['first_active_ns']) & (truth['charge_pe'] >= 50.))
        for i in by_method['baseline'][-1]['score']['target_indices']:
            onset, charge = float(truth['first_active_ns'][i]), float(truth['charge_pe'][i])
            others = physical[physical != i]
            gap = float(np.min(abs(truth['first_active_ns'][others]-onset))) if len(others) else None
            mask = (wave['time_ns'] >= onset-24.) & (wave['time_ns'] < onset-4.)
            history = max(0.,float(wave['analog_mV'][mask].mean()-protocol['baseline_mV'])) if mask.any() else 0.
            event = dict(stem=stem, truth_index=i, particle=str(truth['particle_name'][i]), true_charge_pe=charge,
                preceding_voltage_mV=history, nearest_other_q50_ns=gap, strata=dict(
                    charge='500-1000' if charge < 1000 else '1000-2000' if charge < 2000 else '>=2000',
                    history='<50' if history < 50 else '50-200' if history < 200 else '>=200',
                    spacing='<50' if gap is not None and gap < 50 else '50-150' if gap is not None and gap < 150 else
                            '150-500' if gap is not None and gap < 500 else '>=500'), outcomes={})
            for method in methods:
                result = by_method[method][-1]
                pair = dict(result['score']['matched_pairs']); prediction = None
                outcome = 'missed_at_report_threshold'
                if i in pair:
                    j = pair[i]; fit = result['fit']
                    kind, q = fit['particle_name'][j], fit['charge_pe'][j]
                    outcome = 'wrong_type' if kind != event['particle'] else 'bad_charge' if abs(q/charge-1) > .2 else 'recovered'
                    prediction = dict(index=j, time_ns=fit['times_ns'][j], particle=kind, charge_pe=q)
                event['outcomes'][method] = dict(outcome=outcome, prediction=prediction)
            events.append(event)
    pooled = [dict(method=m, particle=p, **measure(by_method[m],p)) for m in methods for p in ('alpha','gamma')]
    groups = []
    for method in methods:
        for model, rate in sorted(strata):
            subset = [r for r in by_method[method] if (r['row']['gamma_model'],r['row']['rate_hz']) == (model,rate)]
            groups.extend(dict(method=method, gamma_model=model, rate_hz=rate, particle=p,
                               windows=len(subset), **measure(subset,p)) for p in ('alpha','gamma'))
    indexed = {m:{Path(r['row']['root']).stem:r for r in by_method[m]} for m in methods}
    rng = np.random.default_rng(protocol['bootstrap_seed'])
    draws = [[str(s) for cell in strata.values() for s in rng.choice(cell,len(cell),replace=True)] for _ in range(2000)]
    paired = []
    for reference in old_methods:
        for particle in ('alpha','gamma'):
            counts = {s:indexed['class_only'][s]['score']['by_particle'][particle]['recovered_count']-
                        indexed[reference][s]['score']['by_particle'][particle]['recovered_count'] for s in indexed[reference]}
            n = {s:indexed[reference][s]['score']['by_particle'][particle]['target_count'] for s in indexed[reference]}
            deltas = [sum(counts[s] for s in draw)/sum(n[s] for s in draw) for draw in draws if sum(n[s] for s in draw)]
            paired.append(dict(a='class_only', b=reference, particle=particle, delta_count=sum(counts.values()),
                delta_recovery=sum(counts.values())/sum(n.values()),
                paired_window_p025_p50_p975=np.quantile(deltas,[.025,.5,.975]).tolist()))
    diagnostics = {m:{k:sum(r['score'][k] for r in by_method[m]) for k in
        ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags')} for m in methods}
    checks = []
    margin = protocol['development_gate']['precision_margin_absolute']
    for particle in ('alpha','gamma'):
        new = next(r for r in pooled if r['method']=='class_only' and r['particle']==particle)
        old = next(r for r in pooled if r['method']=='baseline' and r['particle']==particle)
        checks.extend([dict(check=particle+'_pooled_recovery', passed=new['recovered_count'] >= old['recovered_count'],
                            candidate=new['recovered_count'], reference=old['recovered_count']),
                       dict(check=particle+'_precision_margin', passed=new['full_precision'] >= old['full_precision']-margin,
                            candidate=new['full_precision'], reference=old['full_precision'], margin=margin)])
    limit = diagnostics['baseline']['unmatched_predictions']*protocol['development_gate']['extra_multiplier']+protocol['development_gate']['extra_offset']
    checks.append(dict(check='extra_predictions_budget', passed=diagnostics['class_only']['unmatched_predictions'] <= limit,
                       candidate=diagnostics['class_only']['unmatched_predictions'], limit=limit))
    for new in (r for r in groups if r['method']=='class_only' and r['target_count'] >= 20):
        old = next(r for r in groups if r['method']=='baseline' and
                   (r['particle'],r['gamma_model'],r['rate_hz']) == (new['particle'],new['gamma_model'],new['rate_hz']))
        margin = protocol['development_gate']['cell_recovery_margin_absolute']
        checks.append(dict(check=f"cell_{new['gamma_model']}_{new['rate_hz']}_{new['particle']}",
            passed=new['full_recovery'] >= old['full_recovery']-margin,
            candidate=new['full_recovery'], reference=old['full_recovery'], margin=margin))
    bins = defaultdict(lambda: defaultdict(int))
    for event in events:
        for method,item in event['outcomes'].items():
            for feature,interval in event['strata'].items():
                bins[(method,event['particle'],feature,interval)]['n'] += 1
                bins[(method,event['particle'],feature,interval)][item['outcome']] += 1
    summary = dict(created=datetime.now().astimezone().isoformat(), windows=20, methods=methods, pooled=pooled,
        groups=groups, paired_contrasts=paired, diagnostics=diagnostics, numerical_diagnostics=numerical,
        event_strata=[dict(method=m,particle=p,feature=f,bin=b,**v) for (m,p,f,b),v in sorted(bins.items())],
        development_gate=dict(passed=all(r['passed'] for r in checks), checks=checks,
                              interpretation='Predeclared project screening rule, not a statistical noninferiority test'),
        scope='Twenty disclosed development windows; full time/class/Q search and residual additions; no new final test',
        target='Q>=500 pe; time match<=10 ns; correct class; |dQ/Q|<=20%; source emission [5,15) us',
        explicit_change='Charge penalty 10->0; class penalty remains 10; same candidates, initial shape Q and search budgets',
        limits=['Neural Q still initializes template shape and influences candidates',
                'Class scores remain uncalibrated','Local search and limited gamma counts',
                'Q is realized charge, not deposited energy','No fresh final test, retraining or G4 generation'])
    save(output/'summary.json',summary); save(output/'events.json',events)
    print(json.dumps(dict(development_gate=summary['development_gate'], pooled=pooled)),flush=True)
    return summary


def figures(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = Path(output); summary = json.loads((output/'summary.json').read_text())
    labels = ['Physical\nfit','Hybrid\nno penalty','Hybrid\nQ + class','Hybrid\nclass only']
    fig, axes = plt.subplots(1,2,figsize=(10.5,4.2))
    for ax,particle in zip(axes,('alpha','gamma')):
        rows = [next(r for r in summary['pooled'] if r['method']==m and r['particle']==particle) for m in summary['methods']]
        x = np.arange(len(rows))
        ax.bar(x-.18,[100*r['full_recovery'] for r in rows],.36,label='Full recovery')
        ax.bar(x+.18,[100*r['full_precision'] for r in rows],.36,label='Full precision')
        ax.set(xticks=x,xticklabels=labels,title=particle,ylabel='Percent',ylim=(0,100))
        ax.grid(axis='y',alpha=.2)
    axes[0].legend(); fig.suptitle('Full optimizer: 20 disclosed development windows'); fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(output/f'full-prior-comparison.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)
    fig, axes = plt.subplots(1,2,figsize=(10.5,4.2))
    outcomes = ['recovered','bad_charge','wrong_type','missed_at_report_threshold']
    colors = ['#3a9b62','#e6b441','#d25c57','#999999']
    events = json.loads((output/'events.json').read_text())
    for ax,particle in zip(axes,('alpha','gamma')):
        subset = [e for e in events if e['particle']==particle]; bottom=np.zeros(4)
        for outcome,color in zip(outcomes,colors):
            counts = np.array([sum(e['outcomes'][m]['outcome']==outcome for e in subset) for m in summary['methods']])
            ax.bar(np.arange(4),counts,bottom=bottom,color=color,label=outcome.replace('_',' '));bottom+=counts
        ax.set(xticks=np.arange(4),xticklabels=labels,title=particle,ylabel='Target events')
    axes[1].legend(loc='upper left',bbox_to_anchor=(1,1));fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(output/f'full-prior-outcomes.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action',required=True)
    fit = sub.add_parser('fit'); fit.add_argument('--work',type=Path,required=True); fit.add_argument('--stem',required=True)
    fit.add_argument('--controls',action='store_true')
    report = sub.add_parser('report')
    for name in ('project','folder','output'):report.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args()
    if args.action == 'fit':fit_case(args.work,args.stem,args.controls)
    else:analyze(args.project,args.folder,args.output);figures(args.output)


if __name__ == '__main__':main()
