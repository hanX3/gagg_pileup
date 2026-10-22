"""Refine the original two-event free-lifetime window scan.

This measures same-event parameter stability against a saved 2500 ns fit.
It does not estimate lifetime accuracy, coverage, or a universal rate limit.
No truth label, source time, charge, or lifetime is passed to the voltage fit.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess

import numpy as np

from .alpha_gamma_waveform_fit import fit_voltage, voltage_model

PROJECT = Path(__file__).resolve().parents[1]
REFERENCE = Path('jupyter/artifacts/alpha_gamma_short_windows_20261013_205215/summary.json')
WAVEFORMS = Path('jupyter/artifacts/alpha_gamma_fit_20261013_165844/fit_arrays.npz')
ENDS_NS = tuple(float(t) for t in range(500, 751, 25))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n')


def stable_from(scans, reference, fields, tolerance):
    """Earliest tested endpoint for which every later tested endpoint passes."""
    tests = [not fit['active_bounds'] and all(
        abs(fit['params'][key]/reference[key]-1.) <= tolerance for key in fields)
        for _, fit in scans]
    return next((time for j, (time, _) in enumerate(scans) if all(tests[j:])), None)


def run(output, records):
    output, records = Path(output).resolve(), Path(records).resolve()
    output.mkdir(parents=True, exist_ok=True)
    records.mkdir(parents=True, exist_ok=True)
    if (output/'summary.json').exists():
        raise FileExistsError('Completed analysis already exists; preserve it.')
    original = json.loads((PROJECT/REFERENCE).read_text())
    for relative, expected in original['source_sha256'].items():
        if sha(PROJECT/relative) != expected:
            raise ValueError(f'Original fitting dependency changed: {relative}')
    source_paths = [Path(__file__).resolve(), *(PROJECT/path for path in original['source_sha256'])]
    protocol = dict(created=datetime.now().astimezone().isoformat(),
        state='frozen_before_dense_scan', endpoints_ns=ENDS_NS, reference_end_ns=2500.,
        tolerances=[.05,.10,.20], lifetime_fields=['tau_fast_ns','tau_slow_ns'],
        lifetime_rule='Both relative errors <= tolerance, no active bound, and every later tested endpoint also passes',
        area_rule='Relative area difference <= tolerance, no active bound, every later tested endpoint passes',
        interpretation='Exploratory stability of the same two known records; long-fit reference is not truth',
        inference_truth_inputs=[], localization='Original full-record observed peak and disjoint pre-pulse baseline',
        pre_crossing_fit_ns=200., model='Original six-parameter free-lifetime convolution fit, original bounds and four starts',
        sources={p.relative_to(PROJECT).as_posix():sha(p) for p in source_paths},
        inputs={REFERENCE.as_posix():sha(PROJECT/REFERENCE),WAVEFORMS.as_posix():sha(PROJECT/WAVEFORMS)},
        environment=dict(python=platform.python_version(), **{p:importlib.metadata.version(p) for p in ('numpy','scipy','matplotlib')}),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        git_status=subprocess.check_output(['git','-c','core.fileMode=false','status','--short'],cwd=PROJECT,text=True))
    write_json(records/'protocol.json', protocol)
    for p in source_paths:
        destination=records/'source'/p.relative_to(PROJECT)
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(p,destination)
    for p in (PROJECT/'jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb', PROJECT/'docs/literature-use.md'):
        destination=records/('notebook-before.ipynb' if p.suffix=='.ipynb' else 'literature-use-before.md')
        if not destination.exists():shutil.copyfile(p,destination)

    combined={}; endpoint_checks=[]
    with np.load(PROJECT/WAVEFORMS) as arrays:
        for label, group in original['results'].items():
            t,y = arrays[f'{label}_time_ns'],arrays[f'{label}_adc_mV']
            readout=group['recorded_config']['readout']
            fits={float(s['end_ns']):s['fits']['free'] for s in group['scans']}
            for end in ENDS_NS:
                fit=fit_voltage(t,y,readout,stop_ns=end,taus=None)
                if not fit['success']:raise RuntimeError((label,end,fit['optimizer_message']))
                ref=fit['reference']; selected=(t>=ref['anchor_ns'])&(t<=ref['anchor_ns']+2000.)
                prediction=voltage_model(t[selected],fit['params'],readout)
                fit['evaluation_0_2000_nrmse_percent_peak']=float(100*np.sqrt(np.mean((y[selected]-prediction)**2))/ref['peak_mV_above_baseline'])
                if end in fits:
                    expected=fits[end]
                    relative=max(abs(fit['params'][key]/expected['params'][key]-1.) for key in fit['params'])
                    sse_relative=abs(fit['sse_mV2']/expected['sse_mV2']-1.)
                    if relative>1e-6 or sse_relative>1e-8:
                        raise AssertionError((label,end,relative,sse_relative))
                    endpoint_checks.append(dict(particle=label,end_ns=end,max_relative_parameter_difference=relative,sse_relative_difference=sse_relative))
                fits[end]=fit
                write_json(output/f'{label}_{int(end):04d}ns.json',fit)
                print(label,int(end),'taus',round(fit['params']['tau_fast_ns'],3),round(fit['params']['tau_slow_ns'],3),'bounds',fit['active_bounds'],flush=True)
            combined[label]=dict(metadata=group['metadata'],scans=[dict(end_ns=end,fit=fit) for end,fit in sorted(fits.items())])

    minima=[]; rows=[]
    for label, group in combined.items():
        scans=[(s['end_ns'],s['fit']) for s in group['scans']]
        reference=next(f['params'] for end,f in scans if end==2500.)
        for tolerance in (.05,.1,.2):
            minima.append(dict(particle=label,tolerance=tolerance,
                both_lifetimes_stable_from_ns=stable_from(scans,reference,('tau_fast_ns','tau_slow_ns'),tolerance),
                area_stable_from_ns=stable_from(scans,reference,('area_mV_ns',),tolerance)))
        for end,fit in scans:
            p=fit['params']
            rows.append(dict(particle=label,end_ns=end,**p,
                tau_fast_relative_to_2500=p['tau_fast_ns']/reference['tau_fast_ns']-1.,
                tau_slow_relative_to_2500=p['tau_slow_ns']/reference['tau_slow_ns']-1.,
                area_relative_to_2500=p['area_mV_ns']/reference['area_mV_ns']-1.,
                active_bounds=fit['active_bounds']))
    report=dict(created=datetime.now().astimezone().isoformat(),protocol='protocol.json',
        new_fits=2*len(ENDS_NS),endpoint_regression=endpoint_checks,results=combined,
        parameter_rows=rows,minimum_stability=minima,
        reference_is_truth=False,independent_particle_events=2,new_ROOT=0,
        limits=['One alpha and one selected partial-deposition gamma only',
                'Same-record stability, not parameter accuracy or confidence coverage',
                'Full-record localization and pre-pulse baseline are inherited; not online latency',
                'Window-dependent original area bounds retained; all bounds reported',
                'Dense 500–750 ns grid, not a continuous minimum or all-charge population threshold'])
    write_json(output/'summary.json',report)
    make_figure(rows,output)
    print(json.dumps(minima,ensure_ascii=False,indent=2),flush=True)
    return report


def make_figure(rows, output):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/gagg-minimum-window-mpl-20261019')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.3),layout='constrained')
    for ax,label in zip(axes,('alpha','gamma')):
        selected=[r for r in rows if r['particle']==label and 500<=r['end_ns']<=1000]
        ax.axhspan(-10,10,color='#d2e7db',alpha=.6,label='±10% stability band')
        ax.axhline(0,color='black',lw=.7)
        for level in (-5,5):ax.axhline(level,color='gray',ls=':',lw=.8)
        for field,color,marker,title in (
            ('tau_fast_relative_to_2500','#2466a8','o','Fast lifetime'),
            ('tau_slow_relative_to_2500','#c67d23','s','Slow lifetime')):
            ax.plot([r['end_ns'] for r in selected],[100*r[field] for r in selected],marker+'-',color=color,label=title)
        ax.set(title=label,xlabel='Fit end T after observed crossing / ns',
               ylabel='Difference from same-event 2500 ns fit / %',xlim=(490,1010))
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Free lifetimes: two original isolated records; stability does not imply accuracy')
    fig.savefig(output/'minimum_window_stability.png',dpi=170)
    fig.savefig(output/'minimum_window_stability.pdf')
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--records',type=Path,required=True)
    args=parser.parse_args()
    run(args.output,args.records)


if __name__=='__main__':main()
