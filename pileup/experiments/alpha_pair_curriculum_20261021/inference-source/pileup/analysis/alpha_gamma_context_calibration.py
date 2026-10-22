"""Independent native calibration, including position and shared-cell history.

Ancestry is used ONLY here to measure the component shapes of labelled
calibration events. Blind reconstruction never imports calibration truth.
"""
from dataclasses import replace
from datetime import datetime
from pathlib import Path
import argparse
import copy
import json
import numpy as np
from scipy.optimize import least_squares, nnls
from .alpha_pileup_fit import AlphaTemplate
from .alpha_rate_campaign import digest, save
from .sipm_response import AVALANCHE_DTYPE, AvalancheResult, ELECTRON_CHARGE_C
from .sipm_readout import ReadoutConfig, synthesize_readout


def preceding_voltage(time_ns, voltage_mV, onset_ns, baseline_mV=10.):
    """An observable, causal pre-onset level in mV, not a microcell occupancy."""
    t = np.asarray(time_ns)
    mask = (t >= onset_ns-24.) & (t < onset_ns-4.)
    return max(0., float(np.mean(np.asarray(voltage_mV)[mask])-baseline_mV)) if mask.any() else 0.


def component_signal(truth, meta, owner=None):
    mask = np.ones(len(truth['avalanche_time_ns']), dtype=bool) if owner is None else truth['avalanche_owner'] == owner
    records = np.zeros(int(mask.sum()), dtype=AVALANCHE_DTYPE)
    records['time_ns'] = truth['avalanche_time_ns'][mask]
    records['charge_C'] = truth['avalanche_charge_C'][mask]
    q0 = meta['config']['sensor']['gain_electrons'] * ELECTRON_CHARGE_C
    avalanches = AvalancheResult(records, meta['duration_ns'], q0, np.array([]), {})
    cfg = replace(ReadoutConfig(**meta['config']['readout']), noise_rms_mV=0.)
    return synthesize_readout(avalanches, cfg, 0).arrays['signal_mV']


def extract(folder, inputs, old_calibration, output, roles):
    folder, inputs, old_calibration, output = map(Path, (folder, inputs, old_calibration, output))
    if output.exists():
        raise FileExistsError(output)
    rows = [r for r in json.loads((folder/'manifest.json').read_text())['runs'] if r['role'] in roles]
    assert rows and all(r['status'] == 'completed' for r in rows)
    cal = json.loads(old_calibration.read_text())
    area = cal['area_per_pe_mV_ns']
    basis = [AlphaTemplate(cal['readout'], f) for f in (0., 1.)]
    observations, checks, hashes = [], [], {}
    for row in rows:
        stem = Path(row['root']).stem
        paths = {part: inputs/part/(stem+suffix) for part, suffix in
                 [('waveforms','.npz'), ('truth','.npz'), ('metadata','.json')]}
        wave = dict(np.load(paths['waveforms'])); truth = dict(np.load(paths['truth']))
        meta = json.loads(paths['metadata'].read_text())
        t, y = wave['time_ns'], wave['analog_mV']
        # Independently verify that the ancestry reconstruction preserves the
        # native readout, including finite aperture and continuous arrival times.
        replay = component_signal(truth, meta)
        maximum = float(np.max(abs(replay-truth['signal_mV'])))
        assert maximum < 2e-7, (stem, maximum)
        checks.append(dict(root=row['root'],max_readout_difference_mV=maximum))
        for j, q in enumerate(truth['charge_pe']):
            onset = float(truth['first_active_ns'][j])
            item = dict(root=row['root'],label=row['label'],role=row['role'],seed=row['seed'],truth_index=j,
                particle=str(truth['particle_name'][j]),charge_pe=float(q),
                source_energy_MeV=float(truth['source_energy_MeV'][j]),edep_MeV=float(truth['edep_MeV'][j]),
                first_active_ns=onset if np.isfinite(onset) else None,
                position_mm=(row['schedule'][j]['position_mm'] if 'schedule' in row else row.get('position_mm',[0.,0.,0.])),
                is_probe=bool('schedule' in row and j>0),gap_ns=row.get('gap_ns'),
                used=bool(q>=50 and np.isfinite(onset)),
                preceding_voltage_mV=preceding_voltage(t,y,onset,cal['baseline_mV']) if np.isfinite(onset) else None)
            if item['used']:
                values = component_signal(truth, meta, j)
                keep = (t>=onset-20.) & (t<=onset+3000.)
                tx, target = t[keep]-onset, values[keep]
                def model(z):
                    f, delay = z
                    return q*area*((1-f)*basis[0](tx-delay)*basis[0].unit_area_peak+
                                   f*basis[1](tx-delay)*basis[1].unit_area_peak)
                opt = least_squares(lambda z:model(z)-target,[.7 if item['particle']=='alpha' else .21,0.],
                    bounds=([.02,-10.],[.98,10.]),x_scale=[.2,2.],max_nfev=100)
                assert opt.success, (row['label'],j,opt.message)
                item.update(slow_fraction=float(opt.x[0]),delay_ns=float(opt.x[1]),
                    solver_success=bool(opt.success),fit_rms_mV=float(np.sqrt(np.mean(opt.fun**2))),
                    component_peak_mV=float(values.max()),shape_source='Native ancestral avalanches; shared-cell losses retained')
            observations.append(item)
        hashes.update({str(p):digest(p) for p in paths.values()})
    result = dict(rows=observations,readout_checks=checks,source_sha256=digest(Path(__file__)),
        manifest_sha256=digest(folder/'manifest.json'),input_sha256=hashes,
        truth_policy='Labelled calibration/diagnostic only; no truth passed to the blind fitter')
    save(output,result)
    for p in ('alpha','gamma'):
        used = [r for r in observations if r['particle']==p and r['used']]
        print('OBSERVATIONS',p,len(used),'Q',np.quantile([r['charge_pe'] for r in used],[0,.5,1]).tolist() if used else [],flush=True)
    return result


def regress(rows, context=False):
    q = np.array([r['charge_pe'] for r in rows])
    y = np.array([r['slow_fraction'] for r in rows])
    X = np.c_[np.ones(len(q)),np.log(q/1000.)]
    if context:
        X = np.c_[X,np.array([r['preceding_voltage_mV'] for r in rows])/100.]
    w = np.sqrt(np.minimum(q,4000.))
    coeff = np.linalg.lstsq(X*w[:,None],y*w,rcond=None)[0]
    residual = y-X@coeff
    return coeff, residual, q, w


def calibrate(observation_file, old_calibration, output):
    observation_file, old_calibration, output = map(Path,(observation_file,old_calibration,output))
    if output.exists():raise FileExistsError(output)
    data = json.loads(observation_file.read_text())
    # No development or validation rows may be used here, even accidentally.
    rows = [r for r in data['rows'] if r['used'] and r['role'] in ('calibration_single','calibration_history')]
    assert rows
    cal = copy.deepcopy(json.loads(old_calibration.read_text()))
    for key in ('isolated_rows','manifest_sha256','waveform_sha256','truth_sha256'):
        cal.pop(key,None)
    cal['created']=datetime.now().astimezone().isoformat()
    cal['classes']={};cal['context_classes']={}
    for p in ('alpha','gamma'):
        isolated = [r for r in rows if r['particle']==p and r['role']=='calibration_single']
        group = isolated+[r for r in rows if r['particle']==p and r['role']=='calibration_history' and r['is_probe']]
        assert len(isolated)>=12 and len(group)>len(isolated)
        coeff,residual,q,w = regress(isolated)
        cal['classes'][p]=dict(slow_intercept=float(coeff[0]),slow_log_slope=float(coeff[1]),
            delay_ns=float(np.average([r['delay_ns'] for r in isolated],weights=w*w)),
            charge_range_pe=[float(q.min()),float(q.max())],used_samples=len(isolated),
            rms_fraction_residual=float(np.sqrt(np.mean(residual**2))))
        cc,rr,cq,cw=regress(group,True)
        # Heteroscedastic residual variance: nonnegative constant + 1/Q term.
        vc=nnls(np.c_[np.ones(len(cq)),1000./cq],rr**2)[0]
        cal['context_classes'][p]=dict(slow_intercept=float(cc[0]),slow_log_slope=float(cc[1]),
            slow_per_100mV=float(cc[2]),charge_range_pe=[float(cq.min()),float(cq.max())],
            history_range_mV=[0.,max(r['preceding_voltage_mV'] for r in group)],
            delay_ns=cal['classes'][p]['delay_ns'],used_samples=len(group),
            rms_fraction_residual=float(np.sqrt(np.mean(rr**2))),
            fraction_variance_coeff=vc.tolist(),
            fraction_residual_p16_p50_p84=np.quantile(rr,[.16,.5,.84]).tolist())
    cal.update(model='Expanded-position/low-charge mean; optional voltage-history linear slow-fraction regression',
        observation_file=str(observation_file),observation_sha256=digest(observation_file),
        parent_calibration_sha256=digest(old_calibration),physics_changed=False,
        source_roles=['calibration_single','calibration_history'],
        noise_policy='Keep previous independent AR(1) noise parameters fixed to isolate template changes',
        context_feature='Causal mean analog voltage in [onset-24,onset-4) ns minus baseline; clamp at zero',
        limitations=['Position is unobserved at inference and marginalized by equal planned position sampling',
                    'Voltage history is an occupancy proxy; it cannot determine local microcell history',
                    'Source energy, deposited energy and realized charge are distinct quantities',
                    'Class/position-dependent detection acceptance affects usable gamma calibration',
                    'Finite calibration residual RMS includes photon statistics, transport and shape-fit variation'])
    save(output,cal)
    print('CALIBRATION',cal['classes'],cal['context_classes'],flush=True)
    return cal


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['extract','calibrate'])
    for arg in ('folder','inputs','old-calibration','output','observations'):p.add_argument('--'+arg,type=Path)
    p.add_argument('--roles',nargs='+')
    a=p.parse_args()
    if a.stage=='extract':extract(a.folder,a.inputs,a.old_calibration,a.output,a.roles)
    else:calibrate(a.observations,a.old_calibration,a.output)
