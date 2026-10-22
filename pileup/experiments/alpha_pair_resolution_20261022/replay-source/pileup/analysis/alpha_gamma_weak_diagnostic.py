"""Truth-assisted weak-gamma subtraction; never a blind recovery benchmark."""
from pathlib import Path
import argparse
import json
from dataclasses import replace
import numpy as np
import matplotlib.pyplot as plt
from .alpha_gamma_study import calibration
from .alpha_gamma_pileup_fit import fit_waveform, score_mixed
from .alpha_rate_campaign import save, digest
from .sipm_response import AvalancheResult, AVALANCHE_DTYPE
from .sipm_readout import ReadoutConfig, synthesize_readout


def diagnose(settings, inputs, original_result, output):
    inputs, output = Path(inputs), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    original = json.loads(Path(original_result).read_text())
    stem = Path(original['row']['root']).stem
    truth = dict(np.load(inputs/'truth'/f'{stem}.npz'))
    wave = dict(np.load(inputs/'waveforms'/f'{stem}.npz'))
    meta = json.loads((inputs/'metadata'/f'{stem}.json').read_text())
    _, bank, config = calibration(settings)
    readout = replace(ReadoutConfig(**meta['config']['readout']), noise_rms_mV=0.)
    q0 = meta['config']['sensor']['gain_electrons']*1.602176634e-19
    records = np.zeros(len(truth['avalanche_time_ns']), dtype=AVALANCHE_DTYPE)
    records['time_ns'] = truth['avalanche_time_ns']
    records['charge_C'] = truth['avalanche_charge_C']
    def electronics(selected):
        a = AvalancheResult(records[selected], meta['duration_ns'], q0, np.array([]), {})
        return synthesize_readout(a, readout, 0).arrays['signal_mV']
    replay = electronics(np.ones(len(records), dtype=bool))
    replay_error = float(np.max(np.abs(replay-truth['signal_mV'])))
    assert replay_error < 1e-6, replay_error
    owner = truth['avalanche_owner']
    weak = (truth['particle_name']=='gamma') & (truth['charge_pe']<500.)
    remove = (owner>=0) & weak[np.maximum(owner, 0)]
    weak_signal = electronics(remove)
    corrected = wave['analog_mV']-weak_signal
    fit = fit_waveform(wave['time_ns'], corrected, bank, config)
    scored = score_mixed(fit, truth, meta['gate_ns'])
    model = fit.pop('model_mV')
    result = dict(root=original['row']['root'], role='Truth-assisted preflight diagnostic, excluded from formal test',
        procedure='Subtract only linear readout of realized subthreshold-gamma avalanches; preserve all actual shared-cell charges, dark chains and native voltage noise',
        limitation='Uses truth; not an achievable blind method, and not a simulation of lower gamma rate because microcell history is unchanged',
        settings_sha256=digest(settings), all_avalanche_readout_max_error_mV=replay_error,
        removed_charge_pe=float(records['charge_C'][remove].sum()/q0),
        original_score=original['score'], diagnostic_score=scored, diagnostic_fit=fit)
    save(output/'weak-gamma-diagnostic.json',result)
    np.savez_compressed(output/'weak-gamma-diagnostic.npz',time_ns=wave['time_ns'],
        original_mV=wave['analog_mV'],subtracted_mV=corrected,weak_signal_mV=weak_signal,model_mV=model)
    fig,axes=plt.subplots(2,1,figsize=(11,6),sharex=True)
    t=wave['time_ns']/1000
    axes[0].plot(t,wave['analog_mV'],label='Original mixture',lw=.8)
    axes[0].plot(t,corrected,label='Truth-assisted subtraction',lw=.8)
    axes[0].plot(t,model,label='Fit after subtraction',lw=.8)
    axes[1].plot(t,weak_signal,lw=.8,color='tab:orange')
    axes[0].legend();axes[0].set_ylabel('Voltage / mV')
    axes[1].set_ylabel('Removed weak gamma / mV');axes[1].set_xlabel('Time / us')
    axes[1].set_xlim(5,15);fig.suptitle('Preflight only: same shared-cell history, truth-assisted subtraction')
    fig.tight_layout()
    for extension in ('png','pdf'):fig.savefig(output/f'weak_gamma_diagnostic.{extension}',dpi=160)
    plt.close(fig)
    print('DIAGNOSTIC',replay_error, {p:(original['score']['by_particle'][p]['recovered_count'],scored['by_particle'][p]['recovered_count']) for p in ('alpha','gamma')},flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('settings','inputs','original-result','output'):parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args();diagnose(args.settings,args.inputs,args.original_result,args.output)
