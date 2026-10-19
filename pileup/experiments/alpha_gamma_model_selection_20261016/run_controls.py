"""Replay the disclosed ideal-template regression controls with frozen code."""
from pathlib import Path
import argparse
import json
import sys

FOLDER = Path(__file__).resolve().parent
sys.path.insert(0, str(FOLDER/'source'))
sys.path.insert(1, str(FOLDER.parents[2]))
import numpy as np
from pileup.analysis.alpha_gamma_regularized_fit import CalibratedBank
from pileup.analysis.alpha_gamma_model_selection_fit import fit_model_selection, ModelSelectionConfig
from pileup.analysis.alpha_gamma_pileup_fit import score_mixed
from pileup.analysis.alpha_rate_campaign import digest, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    frozen = json.loads((FOLDER/'frozen-fit-settings.json').read_text())
    assert digest(FOLDER/'source/pileup/analysis/alpha_gamma_model_selection_fit.py') == frozen['fitter_sha256']
    assert digest(FOLDER/'calibration.json') == frozen['calibration_sha256']
    cal = json.loads((FOLDER/'calibration.json').read_text()); bank = CalibratedBank(cal)
    previous = json.loads(args.previous.read_text())
    times = np.arange(.25, 12000., 2.); rows = []
    for r in previous['rows']:
        if r['fit_stride'] != 2:
            continue
        starts = np.array([5000., 5000.+r['gap_ns']]); charges = np.array(r['truth_charge_pe'])
        kinds = [0 if p == 'alpha' else 1 for p in r['pair']]
        voltage = cal['baseline_mV'] + bank.matrix(times, starts, kinds, charges) @ charges
        voltage += np.random.default_rng(r['seed']).normal(0., .1, len(times))
        fit = fit_model_selection(times, voltage, cal, ModelSelectionConfig(**frozen['config']))
        fit.pop('model_mV')
        truth = dict(first_active_ns=starts, emission_ns=starts, charge_pe=charges, particle_name=bank.names[kinds])
        rows.append(dict(gap_ns=r['gap_ns'], pair=r['pair'], seed=r['seed'], fit=fit,
                         score=score_mixed(fit, truth, [4000., 8000.])))
    save(args.output, dict(rows=rows, fitter_sha256=frozen['fitter_sha256'],
        calibration_sha256=frozen['calibration_sha256'], launcher_sha256=digest(Path(__file__)),
        previous_sha256=digest(args.previous), purpose='Previously exposed ideal controls; not independent native performance'))
    print('CONTROLS', len(rows), 'recovered', sum(r['score']['charge_and_type_recovered_count'] for r in rows), flush=True)


if __name__ == '__main__':
    main()
