"""Information lost by truncating an ideal single-exponential photon process.

This is an analytic design diagnostic, not an ADC fitter or a GAGG performance
prediction. Arrival times are observed exactly in [0, T]; onset is known; there
is no background, readout, SiPM nonlinearity or second decay component. N is the
expected count over [0, infinity), and is unknown to the hypothetical estimator.
The example value N=10,000 is illustrative, not inferred from a project event.

For lambda(t)=N exp(-t/tau)/tau, the Fisher matrix for (log N, log tau)
is N F [[1, g], [g, g*g + v]], with F=1-exp(-x), x=T/tau,
g=-x/expm1(x), v=Var(t/tau | t<=T)=1-[x/(2 sinh(x/2))]**2.
Profiling N gives information N F v about log tau. At x<<1 it is
N x**3/12. Quoted errors are local Cramer-Rao bounds for locally unbiased
estimators under this ideal model, not eventwise error guarantees.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad


def exponential_bounds(window_ns: float, tau_ns: float, expected_total: float) -> dict:
    """Return fractional local uncertainty bounds; all inputs must be positive."""
    values = (window_ns, tau_ns, expected_total)
    if any(not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError('Window, lifetime and expected photon count must be positive and finite')
    x = window_ns / tau_ns
    fraction = -math.expm1(-x)
    g = -x / math.expm1(x) if x < 700 else 0.
    variance = (x*x/12 - x**4/240 + x**6/6048) if x < .01 else (
        1 - (x/(2*math.sinh(x/2)))**2 if x < 1400 else 1.)
    count = expected_total * fraction
    covariance = np.array([[variance + g*g, -g], [-g, 1.]]) / (count*variance)
    fisher = count * np.array([[1., g], [g, variance + g*g]])
    return {
        'window_ns': float(window_ns), 'tau_ns': float(tau_ns), 'expected_total': float(expected_total),
        'expected_observed': count, 'observed_fraction': fraction,
        'relative_sigma_N_tau_fixed': 1/math.sqrt(count),
        'relative_sigma_N_tau_free': math.sqrt(float(covariance[0, 0])),
        'relative_sigma_tau_N_free': math.sqrt(float(covariance[1, 1])),
        'log_parameter_fisher': fisher.tolist(), 'log_parameter_covariance_bound': covariance.tolist(),
    }


def validate_formula() -> dict:
    """Cross-check against independent quadrature of Poisson score products."""
    maximum = 0.
    for x in (.001, .05, .1, .5, 1., 3., 10.):
        result = exponential_bounds(x, 1., 10_000.)
        def entry(i, j):
            return 10_000.*quad(lambda z: math.exp(-z)*(z-1)**(i+j), 0., x,
                                epsabs=1e-13, epsrel=1e-12)[0]
        numeric = np.array([[entry(i,j) for j in range(2)] for i in range(2)])
        analytic = np.array(result['log_parameter_fisher'])
        maximum = max(maximum, float(np.max(np.abs(numeric-analytic))))
        np.testing.assert_allclose(analytic, numeric, rtol=1e-11, atol=1e-10)
        np.testing.assert_allclose(np.linalg.inv(numeric), result['log_parameter_covariance_bound'],
                                   rtol=3e-8, atol=1e-10)
    for x in (.001, .002):
        result = exponential_bounds(x, 1., 10_000.)
        information = 1/result['relative_sigma_tau_N_free']**2
        assert abs(information/(10_000.*x**3/12)-1) < .002
    short, long = [exponential_bounds(t, 351., 10_000.) for t in (100., 500.)]
    assert short['relative_sigma_tau_N_free'] > long['relative_sigma_tau_N_free']
    assert short['relative_sigma_N_tau_free'] > short['relative_sigma_N_tau_fixed']
    return {'status': 'passed', 'independent_fisher_integrals': 28,
            'inverted_matrix_checks': 7, 'small_window_limit_checks': 2,
            'max_absolute_fisher_difference': maximum}


def run(output: Path) -> dict:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output.mkdir(parents=True, exist_ok=True)
    validation = validate_formula()
    times = (50., 100., 150., 200., 300., 500., 750., 1000.)
    rows = [exponential_bounds(t, 351., 10_000.) for t in times]
    x = np.linspace(25., 1000., 250)
    grid = [exponential_bounds(t, 351., 10_000.) for t in x]
    with plt.rc_context({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False}):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
        axes[0].plot(x, [100*r['relative_sigma_N_tau_fixed'] for r in grid], label='Lifetime fixed')
        axes[0].plot(x, [100*r['relative_sigma_N_tau_free'] for r in grid], label='Lifetime also estimated')
        axes[0].set(ylabel='Relative sigma bound on total N (%)', title='Total photon count estimation')
        for n in (1000., 10_000.):
            axes[1].plot(x, [100*exponential_bounds(t, 351., n)['relative_sigma_tau_N_free'] for t in x],
                         label=f'Total N={n:,.0f}; N unknown to estimator')
        axes[1].set(ylabel='Relative sigma bound on lifetime (%)', title='Lifetime estimation')
        for axis in axes:
            axis.set(xlabel='Observation after known onset (ns)', yscale='log')
            axis.grid(alpha=.2)
            axis.legend(fontsize=8)
        fig.suptitle('Ideal one-exponential Poisson example: lifetime 351 ns; no ADC or baseline noise')
        fig.savefig(output/'ideal_window_information.png', dpi=170)
        fig.savefig(output/'ideal_window_information.pdf')
        plt.close(fig)
    report = {'created_utc': datetime.now(timezone.utc).isoformat(),
              'scope': __doc__, 'rows': rows, 'validation': validation,
              'figures': ['ideal_window_information.png', 'ideal_window_information.pdf'],
              'new_G4_events': 0, 'new_ADC_fits': 0, 'literature_performance_claim': False}
    (output/'information.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = run(args.output)
    print(json.dumps(report['validation']))
    for row in report['rows']:
        print(row['window_ns'], 'sigma N fixed/free and sigma tau (%):',
              [round(100*row[k], 3) for k in ('relative_sigma_N_tau_fixed', 'relative_sigma_N_tau_free', 'relative_sigma_tau_N_free')])


if __name__ == '__main__':
    main()
