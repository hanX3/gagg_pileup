"""Short-window sensitivity of the existing isolated GAGG/SiPM fits.

Run from the project root. No new Geant4 events are generated. The full saved
record locates the isolated peak and pre-trigger baseline; this is not a causal
online trigger or a latency measurement. Each window includes 200 ns before the
observed threshold crossing and T ns after it.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess

import numpy as np

from .alpha_gamma_waveform_fit import (
    DEFAULT_INPUTS, PROJECT, digest, fit_voltage, load_event, voltage_model, write_json,
)

WINDOW_ENDS_NS = (100., 200., 300., 500., 750., 1000., 1500., 2000., 2500.)


def analyze_window_scan(events, ends=WINDOW_ENDS_NS):
    """Fit both particle hypotheses to every observed ADC window."""
    ends = tuple(float(end) for end in ends)
    if not ends or any(end <= 0 for end in ends) or len(set(ends)) != len(ends):
        raise ValueError("Window ends must be distinct positive times")
    results = {}
    for label, event in events.items():
        t, y, readout = event['time_ns'], event['adc_mV'], event['config']['readout']
        scans = []
        for end in sorted(ends):
            fits = {
                'fixed': fit_voltage(t, y, readout, stop_ns=end),
                'free': fit_voltage(t, y, readout, stop_ns=end, taus=None),
                'hyp_alpha': fit_voltage(t, y, readout, stop_ns=end, slow_fraction=.70),
                'hyp_gamma': fit_voltage(t, y, readout, stop_ns=end, slow_fraction=.21),
            }
            ref = fits['fixed']['reference']
            full = (t >= ref['anchor_ns']) & (t <= ref['anchor_ns']+2000.)
            for fit in fits.values():
                prediction = voltage_model(t[full], fit['params'], readout)
                fit['evaluation_0_2000_nrmse_percent_peak'] = float(
                    100*np.sqrt(np.mean((y[full]-prediction)**2))/ref['peak_mV_above_baseline'])
            alpha_sse, gamma_sse = (fits['hyp_'+particle]['sse_mV2'] for particle in ('alpha','gamma'))
            delta = gamma_sse-alpha_sse
            scans.append({
                'end_ns': end, 'fit_start_ns_relative': -200.,
                'interval_duration_ns': end+200.,
                'post_crossing_sample_count': int(np.sum((t>=ref['anchor_ns']) & (t<=ref['anchor_ns']+end))),
                'fits': fits, 'delta_sse_mV2': delta,
                'relative_sse_contrast': delta/(gamma_sse+alpha_sse),
                'preferred_hypothesis': 'alpha' if delta>0 else 'gamma',
            })
        results[label] = {'metadata': event['metadata'], 'recorded_config': event['config'], 'scans': scans}
    return results


def make_figures(events, results, output):
    import os
    os.environ.setdefault('MPLCONFIGDIR', str(PROJECT/'jupyter/artifacts/.matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    colors = {'alpha':'#c64a39','gamma':'#2466a8'}
    with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
        fig, axes = plt.subplots(1,3,figsize=(15,4.7),layout='constrained')
        for label, result in results.items():
            scans = result['scans']
            x = [s['end_ns'] for s in scans]
            axes[0].plot(x,[s['fits']['fixed']['params']['slow_fraction'] for s in scans],
                         'o-', color=colors[label],label=label)
            axes[1].plot(x,[s['fits']['free']['params']['tau_slow_ns'] for s in scans],
                         'o-', color=colors[label],label=label)
            axes[2].plot(x,[s['relative_sse_contrast'] for s in scans],
                         'o-',color=colors[label],label='observed '+label)
            # Mark any constrained parameter, including area or baseline.
            for axis, key, field in ((axes[0],'fixed','slow_fraction'),(axes[1],'free','tau_slow_ns')):
                limited = [s for s in scans if s['fits'][key]['active_bounds']]
                axis.plot([s['end_ns'] for s in limited],
                          [s['fits'][key]['params'][field] for s in limited],
                          'x',color='black',ms=9,mew=1.5,ls='none')
        axes[1].axhline(351.,color='0.4',ls='--',label='simulation input: 351 ns')
        axes[2].axhline(0.,color='0.4',ls='--')
        axes[0].set(title='Fixed lifetimes: 95 / 351 ns',ylabel='Fitted slow area fraction',ylim=(0,1))
        axes[1].set(title='Both lifetimes free',ylabel='Effective slow lifetime (ns)',yscale='log')
        axes[2].set(title='Two fixed-shape hypotheses',ylabel='Relative SSE contrast (not probability)',ylim=(-1.05,1.05))
        ticks = [100,200,500,1000,2500]
        for axis in axes:
            axis.set_xscale('log')
            axis.set_xticks(ticks)
            axis.xaxis.set_major_formatter(ScalarFormatter())
            axis.minorticks_off()
            axis.set_xlabel('Fit end T after observed crossing (ns)')
            axis.grid(alpha=.2)
            axis.legend(fontsize=8)
        fig.suptitle('Same two events, window [-200 ns, T]; black crosses: at least one parameter hits a bound')
        fig.savefig(output/'short_window_scan.png',dpi=170)
        plt.close(fig)

        fig,axes=plt.subplots(2,2,figsize=(12,7),sharex=True,layout='constrained')
        for row,(label,event) in enumerate(events.items()):
            for col,end in enumerate((100.,500.)):
                scan=next(s for s in results[label]['scans'] if s['end_ns']==end)
                ref=scan['fits']['fixed']['reference']
                x=event['time_ns']-ref['anchor_ns']
                mask=(x>=-30.) & (x<=1000.)
                scale=ref['peak_mV_above_baseline']
                ax=axes[row,col]
                ax.plot(x[mask],(event['adc_mV'][mask]-ref['baseline_mV'])/scale,
                        color='0.65',lw=1,label='Observed ADC')
                for key,color,title in (('fixed',colors[label],'Fixed lifetimes'),('free','#b1872d','Free lifetimes')):
                    prediction=voltage_model(event['time_ns'][mask],scan['fits'][key]['params'],event['config']['readout'])
                    fit_region=x[mask]<=end
                    for region,style in ((fit_region,'-'),(~fit_region,'--')):
                        ax.plot(x[mask][region],((prediction-ref['baseline_mV'])/scale)[region],
                                style,color=color,lw=1.4,label=title if style=='-' else None)
                ax.axvline(end,color='0.3',ls=':',lw=1)
                ax.axvspan(end,1000,color='0.8',alpha=.18)
                ax.set(title=f'{label}: fit ends at {end:g} ns',ylabel='Voltage / observed peak',ylim=(-.025,1.13))
                if row==1: ax.set_xlabel('Time from observed threshold crossing (ns)')
                ax.legend(fontsize=8)
        fig.suptitle('Solid: fitted interval; dashed/shaded: extrapolation within the same event, not independent validation')
        fig.savefig(output/'short_window_extrapolation.png',dpi=170)
        plt.close(fig)
    return ['short_window_scan.png','short_window_extrapolation.png']


def run_window_scan(output, inputs=None):
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=True)
    inputs=inputs or DEFAULT_INPUTS
    events={label:load_event(PROJECT/path) for label,path in inputs.items()}
    results=analyze_window_scan(events)
    figures=make_figures(events,results,output)
    dependencies=[Path(__file__).resolve(),PROJECT/'analysis/alpha_gamma_waveform_fit.py',
                  PROJECT/'analysis/single_particle_sipm.py']
    report={
        'created_utc':datetime.now(timezone.utc).isoformat(),
        'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        'source_sha256':{str(p.relative_to(PROJECT)):digest(p) for p in dependencies},
        'environment':{'python':platform.python_version(),
            **{p:importlib.metadata.version(p) for p in ('numpy','scipy','matplotlib','uproot','awkward')}},
        'window_ends_ns':list(WINDOW_ENDS_NS),'results':results,'figures':figures,
        'score_definition':'(SSE_gamma-SSE_alpha)/(SSE_gamma+SSE_alpha), descriptive only; not probability',
        'fit_truth_inputs':[],
        'timing_scope':'Full-record observed peak and 1500 ns pre-pulse baseline; not an online latency estimate',
        'limits':'One independent alpha and one selected gamma, no PSD efficiency or confidence intervals',
    }
    write_json(output/'summary.json',report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=run_window_scan(args.output)
    for label,result in report['results'].items():
        for scan in result['scans']:
            fixed,free=scan['fits']['fixed'],scan['fits']['free']
            print(label,scan['end_ns'],'ws',round(fixed['params']['slow_fraction'],4),
                  'free taus',round(free['params']['tau_fast_ns'],2),round(free['params']['tau_slow_ns'],2),
                  'bounds',fixed['active_bounds'],free['active_bounds'],
                  'preferred',scan['preferred_hypothesis'],'contrast',round(scan['relative_sse_contrast'],4),flush=True)


if __name__=='__main__': main()
