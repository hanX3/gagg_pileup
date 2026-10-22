"""Controlled 0/1/2 photon injections, jointly digitized with shared SiPM cells.

This is offline redigitization of grouped Geant4 optical transports, with new
digitizer/readout seeds. It is distinct both from linear waveform sums and from
new complete C++ native G4 runs. No event truth is passed to an inference model.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from dataclasses import replace
from pathlib import Path
import shutil
import tempfile

import numpy as np

from . import alpha_pair_data as data
from .alpha_pair_data import PROJECT,EXPERIMENT,OUTPUT,load
from .alpha_rate_campaign import digest,save
from .sipm_response import PhotonHits,SensorConfig,read_photon_hits,simulate_avalanches,map_cells
from .sipm_readout import ReadoutConfig,synthesize_readout


def donor_cache(role,records):
    folder=OUTPUT/'photons'/role
    folder.mkdir(parents=True,exist_ok=True)
    result=[]
    for i,r in enumerate(records):
        path=folder/f'{i:04}.npz'
        if not path.exists():
            with tempfile.TemporaryDirectory(prefix='gagg-pair-donor-') as td:
                local=Path(td)/Path(r['root']).name
                shutil.copyfile(PROJECT/r['root'],local)
                assert digest(local)==load(PROJECT/r['metadata'])['root_sha256']
                hits=read_photon_hits(local)
            config=load(PROJECT/r['metadata'])['config']
            sensor=SensorConfig(**config['sensor'])
            cells=map_cells(hits,sensor)
            first=float(hits.time_ns[cells>=0].min()) if np.any(cells>=0) else 5000.
            np.savez_compressed(path,time_ns=hits.time_ns-first,x_mm=hits.x_mm,y_mm=hits.y_mm,energy_eV=hits.energy_eV)
        with np.load(path) as f:
            result.append({k:f[k].copy() for k in f.files})
        if (i+1)%25==0:
            print('PHOTON_DONORS',role,i+1,len(records),flush=True)
    save(folder/'manifest.json',dict(records=records,files={p.name:digest(p) for p in folder.glob('*.npz')},
        scope='Original transport arrays, first-active-photon-relative; no photon-amplitude scaling'))
    return result


def one_case(donors,config,seed):
    rng=np.random.default_rng(seed)
    count=int(rng.choice((0,1,2),p=(.2,.3,.5)))
    ids=rng.choice(len(donors),count,replace=count>len(donors)).tolist()
    first=float(rng.uniform(96.,288.));gap=float(rng.choice(data.GAPS))
    window=4096.
    columns={k:[] for k in ('time_ns','x_mm','y_mm','energy_eV')};owners=[]
    for j,i in enumerate(ids):
        times=donors[i]['time_ns']+first+j*gap
        keep=(times>=0.)&(times<window)
        for k in columns:
            columns[k].append((times if k=='time_ns' else donors[i][k])[keep])
        owners.extend([j]*int(keep.sum()))
    arrays={k:np.concatenate(v) if v else np.array([]) for k,v in columns.items()}
    owners=np.asarray(owners,int);order=np.argsort(arrays['time_ns'],kind='stable')
    arrays={k:v[order] for k,v in arrays.items()};owners=owners[order]
    hits=PhotonHits(**arrays,window_ns=window)
    sensor=SensorConfig(**config['sensor']);readout=ReadoutConfig(**config['readout'])
    avalanches=simulate_avalanches(hits,sensor,seed+100000)
    output=synthesize_readout(avalanches,readout,seed+200000)
    rec=avalanches.records
    avalanche_owner=np.full(len(rec),-1,int)
    direct=rec['kind']==0
    avalanche_owner[direct]=owners[rec['photon_index'][direct]]
    for i in np.flatnonzero(rec['parent_id']>=0):
        assert rec['parent_id'][i]<i
        avalanche_owner[i]=avalanche_owner[rec['parent_id'][i]]
    q=np.zeros(2);times=np.full(2,first);curves=np.zeros((2,data.SAMPLES))
    replay=replace(readout,baseline_mV=0.,noise_rms_mV=0.)
    signal_sum=np.zeros_like(output.arrays['signal_mV'])
    cells=map_cells(hits,sensor)
    from .sipm_response import AvalancheResult
    charge_error=0.
    for owner in range(-1,count):
        block=rec[avalanche_owner==owner]
        contribution=synthesize_readout(AvalancheResult(block,window,avalanches.full_charge_C,np.array([]),{}),replay,0).arrays['signal_mV']
        signal_sum+=contribution
        if owner>=0:
            q[owner]=block['charge_C'].sum()/avalanches.full_charge_C
            mask=(owners==owner)&(cells>=0)
            times[owner]=float(hits.time_ns[mask].min()) if mask.any() else first+owner*gap
            curves[owner]=contribution[:data.SAMPLES]
            charge_error+=block['charge_C'].sum()
    charge_error=abs(charge_error+rec['charge_C'][avalanche_owner<0].sum()-rec['charge_C'].sum())
    replay_error=float(np.max(abs(signal_sum-output.arrays['signal_mV'])))
    assert charge_error<1e-20 and replay_error<1e-8
    active=(q>=500.).astype(np.float32)
    curves*=active[:,None]
    batch=dict(voltage_mV=output.arrays['analog_mV'][:data.SAMPLES],active=active,
               charge_pe=np.maximum(q,1.),times_ns=times,components_mV=curves)
    info=dict(seed=seed,sipm_seed=seed+100000,readout_seed=seed+200000,parents=ids,
              injected_count=count,qualified_count=int(active.sum()),gap_ns=gap if count==2 else None,
              window_ns=window,observed_ns=data.SAMPLES*2.,shared_cells=True,
              charge_balance_error_C=charge_error,component_sum_max_error_mV=replay_error)
    return batch,info


def create(role,*,count,seed,workers=2,records=None):
    if not 1<=workers<=4:
        raise ValueError('Use at most4 processes/threads')
    destination=OUTPUT/f'{role}-shared-pair.npz'
    if destination.exists():
        raise FileExistsError(destination)
    if records is None:
        records=load(EXPERIMENT/'data-split.json')[role]
    donors=donor_cache(role,records)
    config=load(PROJECT/records[0]['metadata'])['config']
    output=[None]*count;info=[None]*count
    # Sensor sampling is Python-heavy; threads also avoid duplicating optical
    # arrays. Maximum CPU concurrency is recorded alongside every generation.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(one_case,donors,config,seed+i):i for i in range(count)}
        for future in as_completed(futures):
            i=futures[future];output[i],info[i]=future.result()
            if sum(x is not None for x in output)%32==0:
                print('SHARED_PAIR',role,sum(x is not None for x in output),count,flush=True)
    batch={k:np.asarray([r[k] for r in output],np.float32) for k in output[0]}
    data.save_batch(destination,batch,info)
    save(destination.with_suffix('.provenance.json'),dict(role=role,parents=len(records),count=count,seed=seed,
        workers=workers,scope='Joint offline SiPM redigitization; shared pixels, PDE, XT/dark, full analog noise; same model laws, new RNG',
        transport_inputs='Grouped independent native single-alpha ROOTs; artificial controlled shifts of photon times',
        source_sha256={str(p.relative_to(PROJECT)):digest(p) for p in (Path(__file__),PROJECT/'pileup/analysis/sipm_response.py',PROJECT/'pileup/analysis/sipm_readout.py')}))
    return batch


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('role',choices=('train','validation'))
    p.add_argument('--workers',type=int,default=2)
    a=p.parse_args()
    create(a.role,count=512 if a.role=='train' else 128,seed=262052000 if a.role=='train' else 262053000,workers=a.workers)


if __name__=='__main__':main()
