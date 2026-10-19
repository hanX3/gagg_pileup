"""Fixed empty/single/four-type pair tests: matched mean templates, white noise.

These tests deliberately omit photon and cell fluctuations. They diagnose
algorithmic failures and equal-charge class cues, not detector resolution.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import torch
from .alpha_rate_campaign import save,digest
from .alpha_gamma_ml_study import load_model
from .alpha_gamma_ml import voltage_features,predict_features,decode
from .alpha_gamma_ml_hybrid import fit_hybrid,MixedHybridConfig
from .alpha_gamma_model_selection_fit import fit_model_selection,ModelSelectionConfig
from .alpha_gamma_regularized_fit import CalibratedBank
from .alpha_gamma_pileup_fit import score_mixed


def run(folder,output):
    folder,output=Path(folder),Path(output);torch.set_num_threads(1)
    frozen=json.loads((folder/'frozen-fit-settings.json').read_text())
    cal=json.loads((folder/'calibration.json').read_text());bank=CalibratedBank(cal)
    networks={kind:load_model(output/kind) for kind in ('cnn','unet')}
    cases=[dict(label='empty',types=[],times=[],charges=[])]
    for k in (0,1):
        for q in (500.,1000.,3000.,10000.):
            cases.append(dict(label=f'single_{bank.names[k]}_{q:g}',types=[k],times=[8000.],charges=[q]))
    for gap in (10.,20.,40.,80.,200.):
        for pair in ((0,0),(0,1),(1,0),(1,1)):
            for charge in ([1000.,1000.],[3000.,600.]):
                label=f'pair_{bank.names[pair[0]]}_{bank.names[pair[1]]}_{gap:g}_{charge[0]:g}_{charge[1]:g}'
                cases.append(dict(label=label,types=list(pair),times=[8000.,8000.+gap],charges=charge))
    destination=output/'synthetic-controls';destination.mkdir(exist_ok=True)
    rng=np.random.default_rng(261084500);rows=[];t=np.arange(10000)*2.;gate=[5000,15000]
    for case in cases:
        q=np.array(case['charges']);types=np.array(case['types'],int);times=np.array(case['times'])
        voltage=cal['baseline_mV']+bank.matrix(t,times,types,q)@q+rng.normal(0,.1,len(t))
        truth=dict(first_active_ns=times,emission_ns=times,particle_name=bank.names[types],charge_pe=q)
        path=destination/(case['label']+'.json')
        if path.exists():rows.append(json.loads(path.read_text()));continue
        fits={};seconds={};dense={};features=voltage_features(t,voltage,cal)
        for kind,(network,cfg,freeze) in networks.items():
            start=time.perf_counter();dense[kind]=predict_features(network,features,cfg)
            fits[kind]=decode(t,dense[kind],freeze['selected']['thresholds'],cfg);seconds[kind]=time.perf_counter()-start
        start=time.perf_counter();fits['baseline']=fit_model_selection(t,voltage,cal,ModelSelectionConfig(addition_cost=1.))
        seconds['baseline']=time.perf_counter()-start
        kind=frozen['hybrid_kind'];method=frozen['hybrid_method'];start=time.perf_counter()
        fits[method]=fit_hybrid(t,voltage,cal,fits[kind],dense[kind],fits['baseline'],MixedHybridConfig(**frozen['hybrid_config']))
        seconds[method]=time.perf_counter()-start
        records={}
        for name,fit in fits.items():
            score=score_mixed(fit,truth,gate);fit=dict(fit);fit.pop('model_mV',None)
            records[name]=dict(fit=fit,score=score,seconds=seconds[name])
        row=dict(case=case,methods=records,calibration_sha256=digest(folder/'calibration.json'),
            model='Matched mean calibration templates + independent 0.1 mV white electronics noise; no photon/microcell fluctuations',
            oracle_inputs_to_fits=[],source_sha256=digest(Path(__file__)))
        save(path,row);rows.append(row)
        print('CONTROL',case['label'],{m:r['score']['charge_and_type_recovered_count'] for m,r in records.items()},flush=True)
    save(output/'synthetic-control-summary.json',dict(rows=rows,cases=len(rows),seed=261084500,
        scope='Unknown-count/class blind inference on synthetic controls; not physical resolution or native gamma-gamma validation'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.folder,a.output)
