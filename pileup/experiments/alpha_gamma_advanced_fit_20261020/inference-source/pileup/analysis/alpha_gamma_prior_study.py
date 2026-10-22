"""Replay cached blind outputs, then score the fixed-candidate prior ablation."""
from pathlib import Path
from collections import defaultdict
from datetime import datetime
import argparse
import hashlib
import json
import time
import numpy as np
from .alpha_gamma_prior_ablation import prepare, fit_fixed, decode_cached, ARMS
from .alpha_gamma_pileup_fit import score_mixed


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def fit_case(work, stem):
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
    decoded = decode_cached(t, dense, neural['fit']['thresholds'])
    assert decoded == neural['fit']['all_components'], 'Dense cache and saved neural events disagree'
    start = time.perf_counter()
    common = prepare(t, voltage, json.loads((work/'calibration.json').read_text()),
                     neural['fit'], dense, baseline['fit'])
    prepare_seconds = time.perf_counter()-start
    save(work/'outputs/common'/f'{stem}.json', dict(common['metadata'], candidate_sha256=common['candidate_sha256'],
        waveform_sha256=case['inputs']['waveform']['sha256'], prepare_seconds=prepare_seconds))
    results = []
    # This stage's case files contain no truth, source rate, seed, or known count.
    for arm, config in ARMS.items():
        start = time.perf_counter(); fitted = fit_fixed(common, config); seconds = time.perf_counter()-start
        model = fitted.pop('model_mV')
        destination = work/'outputs/fits'/arm/f'{stem}.json'
        if destination.exists():
            raise FileExistsError(destination)
        npz = destination.with_suffix('.npz'); npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(npz, time_ns=t, model_mV=model)
        saved = dict(stem=stem, arm=arm, fit=fitted, seconds=seconds,
                     calibration_sha256=case['calibration_sha256'], input_sha256={k:v['sha256'] for k,v in case['inputs'].items()},
                     candidate_sha256=common['candidate_sha256'], fit_truth_inputs=[], truth_file_access=False,
                     algorithm_sha256={n:digest(Path(__file__).with_name(n)) for n in
                                       ('alpha_gamma_prior_ablation.py','alpha_gamma_prior_study.py')})
        save(destination, saved); results.append(dict(arm=arm, seconds=seconds, reported=len(fitted['times_ns']),
                                                      passes=fitted['passes'], objective=fitted['objective_history'][-1]))
    print(json.dumps(dict(stem=stem, candidates=len(common['times']), results=results)), flush=True)


def measure(rows, particle):
    keys = ('target_count','matched_count','correct_type_count','recovered_count','missed_count','wrong_type_count')
    counts = {key:sum(r['score']['by_particle'][particle][key] for r in rows) for key in keys}
    predictions = pure = precise = 0
    errors = []; correct_errors = []
    for r in rows:
        fit, truth, score = r['fit'], r['truth'], r['score']
        gate = r['gate_ns']
        mask = [p == particle and gate[0] <= t < gate[1] for p,t in zip(fit['particle_name'],fit['times_ns'])]
        predictions += sum(mask)
        for i,j in score['all_matched_pairs']:
            if mask[j] and truth['particle_name'][i] == particle:
                pure += 1; precise += int(abs(fit['charge_pe'][j]/truth['charge_pe'][i]-1) <= .2)
        for i,j in score['matched_pairs']:
            if truth['particle_name'][i] == particle:
                error = float(fit['charge_pe'][j]/truth['charge_pe'][i]-1)
                errors.append(error)
                if fit['particle_name'][j] == particle:
                    correct_errors.append(error)
    recall = counts['recovered_count']/counts['target_count'] if counts['target_count'] else None
    precision = precise/predictions if predictions else 0.
    def moments(values):
        return dict(n=len(values), bias=float(np.mean(values)) if values else None,
                    std=float(np.std(values,ddof=1)) if len(values)>1 else None)
    return dict(**counts, prediction_count=predictions, pure_prediction_count=pure, precise_prediction_count=precise,
                full_recovery=recall, full_precision=precision, type_purity=pure/predictions if predictions else None,
                full_f1=2*recall*precision/(recall+precision) if recall and precision else 0.,
                conditional_charge= moments(errors), conditional_correct_type_charge=moments(correct_errors))


def analyze(project, folder, output):
    project, folder, output = map(Path, (project, folder, output))
    protocol = json.loads((folder/'protocol.json').read_text())
    rows = protocol['scoring_rows']
    assert len(rows) == 20
    assert all((output/'fits'/arm/(Path(r['root']).stem+'.json')).exists() for r in rows for arm in ARMS)
    methods = list(ARMS)+['original_baseline','original_hybrid_0','original_hybrid_10']
    by_method = defaultdict(list); fixed_events = []; reported_events = []
    strata = defaultdict(list)
    for row in rows:
        stem = Path(row['root']).stem
        strata[(row['gamma_model'],row['rate_hz'])].append(stem)
        assert digest(project/row['truth']) == row['input_sha256']['truth']
        with np.load(project/row['truth']) as data:
            truth = {key:data[key] for key in ('first_active_ns','emission_ns','particle_name','charge_pe')}
        meta = json.loads((project/row['metadata']).read_text()); gate = meta['gate_ns']
        wave = dict(np.load(project/row['waveform']))
        common = json.loads((output/'common'/f'{stem}.json').read_text())
        dummy = dict(times_ns=common['times_ns'], particle_name=['alpha']*len(common['times_ns']),
                     charge_pe=[500.]*len(common['times_ns']))
        association = score_mixed(dummy, truth, gate); fixed_pair = dict(association['matched_pairs'])
        physical = np.flatnonzero(np.isfinite(truth['first_active_ns']) & (truth['charge_pe'] >= 50.))
        fits = {}
        for method in methods:
            if method in ARMS:
                record = json.loads((output/'fits'/method/f'{stem}.json').read_text())
                assert record['candidate_sha256'] == common['candidate_sha256']
                fitted = record['fit']
                assert fitted['all_components']['times_ns'] == common['times_ns']
                assert fitted['all_components']['shape_charge_pe'] == common['shape_charge_pe']
            else:
                old = {'original_baseline':'baseline','original_hybrid_0':'hybrid_unet_0','original_hybrid_10':'hybrid_unet_10'}[method]
                record = json.loads((project/protocol['parent_artifacts']/'results/validation'/old/f'{stem}.json').read_text())
                fitted = record['fit']
            score = score_mixed(fitted, truth, gate)
            by_method[method].append(dict(row=row, fit=fitted, truth=truth, score=score, gate_ns=gate))
            fits[method] = fitted
            save(output/'scores'/method/f'{stem}.json', dict(row=row, method=method, score=score, gate_ns=gate))
        for i in association['target_indices']:
            particle = str(truth['particle_name'][i]); onset = float(truth['first_active_ns'][i]); charge = float(truth['charge_pe'][i])
            others = physical[physical != i]
            gap = float(np.min(abs(truth['first_active_ns'][others]-onset))) if len(others) else None
            mask = (wave['time_ns'] >= onset-24.) & (wave['time_ns'] < onset-4.)
            history = max(0., float(wave['analog_mV'][mask].mean()-10.)) if mask.any() else 0.
            descriptors = dict(charge='500-1000' if charge < 1000 else '1000-2000' if charge < 2000 else '>=2000',
                history='<50' if history < 50 else '50-200' if history < 200 else '>=200',
                spacing='<50' if gap is not None and gap < 50 else '50-150' if gap is not None and gap < 150 else
                        '150-500' if gap is not None and gap < 500 else '>=500')
            event = dict(stem=stem, particle=particle, truth_index=i, true_charge_pe=charge,
                         preceding_voltage_mV=history, nearest_other_q50_ns=gap, strata=descriptors, outcomes={})
            for arm in ARMS:
                outcome = 'no_candidate'; prediction = None
                if i in fixed_pair:
                    j = fixed_pair[i]; fitted = fits[arm]['all_components']
                    q = fitted['charge_pe'][j]; kind = fitted['particle_name'][j]
                    outcome = ('below_report' if q < 500 else 'wrong_type' if kind != particle else
                               'bad_charge' if abs(q/charge-1) > .2 else 'recovered')
                    prediction = dict(candidate_index=j, particle=kind, charge_pe=q,
                                      neural_prior_charge_pe=common['neural_prior_charge_pe'][j],
                                      raw_class_cost_per_strength=common['class_cost_per_strength'][j])
                event['outcomes'][arm] = dict(outcome=outcome, prediction=prediction)
            fixed_events.append(event)
    pooled = [dict(method=m, particle=p, **measure(by_method[m],p)) for m in methods for p in ('alpha','gamma')]
    grouped = []
    for method in methods:
        for model, rate in sorted(strata):
            subset = [r for r in by_method[method] if (r['row']['gamma_model'],r['row']['rate_hz']) == (model,rate)]
            grouped.extend(dict(method=method, gamma_model=model, rate_hz=rate, particle=p,
                                windows=len(subset), **measure(subset,p)) for p in ('alpha','gamma'))
    bins = defaultdict(lambda:defaultdict(int))
    for e in fixed_events:
        for arm,item in e['outcomes'].items():
            for feature,interval in e['strata'].items():
                bins[(arm,e['particle'],feature,interval)]['n'] += 1
                bins[(arm,e['particle'],feature,interval)][item['outcome']] += 1
    contrasts = [('charge','none'),('class','none'),('both','charge'),('both','class'),('both','none')]
    paired = []
    indexed = {m:{Path(r['row']['root']).stem:r for r in by_method[m]} for m in ARMS}
    rng = np.random.default_rng(261091001)
    draws = [[str(s) for cell in strata.values() for s in rng.choice(cell,len(cell),replace=True)] for _ in range(2000)]
    for a,b in contrasts:
        for p in ('alpha','gamma'):
            deltas = []
            for sample in draws:
                n = sum(indexed[a][s]['score']['by_particle'][p]['target_count'] for s in sample)
                if n:
                    deltas.append(sum(indexed[a][s]['score']['by_particle'][p]['recovered_count']-
                                      indexed[b][s]['score']['by_particle'][p]['recovered_count'] for s in sample)/n)
            numerator = sum(indexed[a][s]['score']['by_particle'][p]['recovered_count']-
                            indexed[b][s]['score']['by_particle'][p]['recovered_count'] for s in indexed[a])
            n = sum(indexed[a][s]['score']['by_particle'][p]['target_count'] for s in indexed[a])
            paired.append(dict(a=a,b=b,particle=p,delta_count=numerator,delta_recovery=numerator/n,
                               paired_window_p025_p50_p975=np.quantile(deltas,[.025,.5,.975]).tolist()))
    summary = dict(created=datetime.now().astimezone().isoformat(), windows=20, methods=methods,
        pooled=pooled, groups=grouped, paired_contrasts=paired,
        diagnostics={m:{k:sum(r['score'][k] for r in by_method[m]) for k in
                     ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags')} for m in methods},
        fixed_association_bins=[dict(arm=a,particle=p,feature=f,bin=b,**counts) for (a,p,f,b),counts in sorted(bins.items())],
        scope='Twenty disclosed development windows; fixed blind candidate onsets/counts, templates and whitening; no time search/additions',
        priors='Uncalibrated scores and charge penalties at the previous strength 10; no tuning on these results',
        target='Q>=500 pe; time match<=10 ns; correct class; |dQ/Q|<=20%; source emission [5,15) us',
        fixed_association='Separate common time-only mapping before report-threshold selection; truth used only at scoring',
        limits=['Development diagnostic, not new held-out evidence','Limited gamma counts and local discrete minima',
                'None arm retains common neural candidates, initial labels and shape charges',
                'Fixed onsets cannot reproduce production hybrid arrival-time or residual-addition effects',
                'Q is realized charge, not calibrated deposited energy','No new training, probability calibration or G4 generation'])
    save(output/'summary.json',summary); save(output/'fixed-association-events.json',fixed_events)
    print(json.dumps(dict(windows=20,pooled=[{k:r[k] for k in ('method','particle','target_count','recovered_count','wrong_type_count','full_precision')} for r in pooled])),flush=True)
    return summary


def figures(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output=Path(output); summary=json.loads((output/'summary.json').read_text())
    arms=list(ARMS); labels=['No explicit\nprior','Charge\nprior','Class\nprior','Both\npriors']
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,p in zip(axes,('alpha','gamma')):
        rows=[next(r for r in summary['pooled'] if r['method']==a and r['particle']==p) for a in arms]
        x=np.arange(4)
        ax.bar(x-.18,[100*r['full_recovery'] for r in rows],.36,label='Full recovery')
        ax.bar(x+.18,[100*r['full_precision'] for r in rows],.36,label='Full precision')
        ax.set(xticks=x,xticklabels=labels,title=p,ylim=(0,100),ylabel='Percent')
        ax.grid(axis='y',alpha=.2)
    axes[0].legend();fig.suptitle('Fixed blind candidates: 20 disclosed development windows')
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(output/f'prior-comparison.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)
    events=json.loads((output/'fixed-association-events.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    outcomes=['recovered','bad_charge','wrong_type','below_report','no_candidate']
    colors=['#3a9b62','#e6b441','#d25c57','#8c78ad','#999999']
    for ax,p in zip(axes,('alpha','gamma')):
        subset=[e for e in events if e['particle']==p];bottom=np.zeros(4)
        for outcome,color in zip(outcomes,colors):
            counts=np.array([sum(e['outcomes'][a]['outcome']==outcome for e in subset) for a in arms])
            ax.bar(np.arange(4),counts,bottom=bottom,color=color,label=outcome.replace('_',' '));bottom+=counts
        ax.set(xticks=np.arange(4),xticklabels=labels,title=p+' common time association',ylabel='Target events')
    axes[1].legend(loc='upper left',bbox_to_anchor=(1,1));fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(output/f'fixed-association-outcomes.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=('fit','report'))
    p.add_argument('--work',type=Path);p.add_argument('--stem');p.add_argument('--project',type=Path)
    p.add_argument('--folder',type=Path);p.add_argument('--output',type=Path);args=p.parse_args()
    if args.stage=='fit':fit_case(args.work,args.stem)
    else:analyze(args.project,args.folder,args.output);figures(args.output)
