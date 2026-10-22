"""Matched-template double-pulse controls, explicitly separate from native tests.

These intentionally challenge an at-most-one-pulse prefix hypothesis. They do
not simulate shared-cell recovery and cannot establish physical pileup recovery.
"""
import argparse
from pathlib import Path
import numpy as np

from .alpha_gamma_short_window import PrefixCalibration, fit_prefix
from .alpha_gamma_short_window_study import save


def run(output):
    output=Path(output);cal=PrefixCalibration.load(output/'calibration')
    t=np.arange(.25,9500.,2.);rows=[]
    for first,second in [('alpha','gamma'),('gamma','alpha')]:
        for gap in (50.,100.,200.,400.):
            for ratio in (.25,1.,2.):
                seed=261099000+len(rows)
                q1,q2=4000.,4000.*ratio
                y=(cal.meta['baseline_mV']+q1*cal.template(t-5000.,first,q1)
                   +q2*cal.template(t-5000.-gap,second,q2)
                   +np.random.default_rng(seed).normal(0.,.1,len(t)))
                for post in (100.,200.,300.):
                    fit=fit_prefix(t,y,cal,pre_ns=50.,post_ns=post,methods=('physics_ols','template_gls'))
                    for result in fit['methods'].values():
                        result['relative_first_charge_error']=result['charge_pe']/q1-1
                    rows.append({'first_particle':first,'second_particle':second,'gap_ns':gap,'charge_ratio':ratio,
                        'post_ns':post,'seed':seed,'fit':fit,
                        'second_signal_can_enter_prefix':bool(fit['detected'] and fit['last_observed_ns']>=5000.+gap-2.)})
    summary={'rows':rows,'scope':'Matched independent mean-template sum plus electronics noise; NOT native shared-cell pileup',
        'inference_truth_inputs':[],'estimator_hypothesis':'At most one pulse; two pulses are a deliberate model mismatch',
        'known_control_truth':'Pair, times and Q specify generation/evaluation only; never passed to fit_prefix'}
    save(output/'controlled-doubles.json',summary)
    print('CONTROLLED_PREFIXES',len(rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output)
