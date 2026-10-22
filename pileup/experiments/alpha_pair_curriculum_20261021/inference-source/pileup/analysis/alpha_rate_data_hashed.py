"""Export native SiPM samples and, separately, photon-parent truth for scoring.

Read-only I/O optimization: full records use checked hash lookup rather than
repeated structured binary search; rare hash collisions fall back to the original.
Truth attribution follows accepted photon indices and secondary-parent links.
It does not redigitize individual pulses: shared-cell recovery remains that of
the actual pileup. Photon labels never enter the blind waveform input file.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import awkward as ak
import uproot
from uproot.source.file import MemmapSource

from .alpha_rate_campaign import digest,save
from .root_redundant_photon_fields import record_hash


def open_root(path):
    return uproot.open(path,handler=MemmapSource,num_workers=1,
        decompression_executor=uproot.TrivialExecutor(),interpretation_executor=uproot.TrivialExecutor())


def photon_key(t,x):
    # Equality of integer picosecond time and the exact stored float32 x bits.
    # Uniqueness is checked before use; ambiguity is an error, never guessed.
    return (np.asarray(t,dtype=np.uint64)<<np.uint64(32))|np.asarray(x,dtype=np.float32).view(np.uint32).astype(np.uint64)


def photon_full_key(t,x,y,energy):
    """Exact saved photon record, including coordinates and energy float bits.

    Rounded picosecond times and boundary coordinates can coincide for distinct
    photons. Never break a remaining full-record tie using row order or truth.
    """
    key=np.empty(len(t),dtype=[('t_ps','<i8'),('x_bits','<u4'),('y_bits','<u4'),('energy_bits','<u4')])
    key['t_ps']=np.asarray(t,dtype=np.int64)
    for name,values in (('x_bits',x),('y_bits',y),('energy_bits',energy)):
        key[name]=np.asarray(values,dtype=np.float32).view(np.uint32)
    return key


def export_window(root_path,output, *, expected=None, allowed_particles=('alpha',)):
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    root_path=Path(root_path); stem=root_path.stem
    for name in ('waveforms','truth','metadata'):(output/name).mkdir(exist_ok=True)
    if (output/'metadata'/f'{stem}.json').exists():
        cached=json.loads((output/'metadata'/f'{stem}.json').read_text())
        with np.load(output/'truth'/f'{stem}.npz') as truth:
            if 'particle_name' in truth:
                unexpected=set(np.unique(truth['particle_name']))-set(allowed_particles)
                if unexpected:raise ValueError(f'Unexpected cached primary particles: {unexpected}')
            elif set(allowed_particles)!={'alpha'}:
                raise ValueError('Legacy export has no particle labels; use a new mixed-export directory')
        if expected:
            assert cached['seed']==expected['seed'] and cached['rate_hz']==expected['rate_hz']
            if 'root_sha256' in expected:assert cached['root_sha256']==expected['root_sha256']
            if 'gamma_rate_hz' in expected:
                slot='brems' if expected['gamma_model']=='brems' else 'c12_capture'
                assert cached['source_rates_hz'][slot]==expected['gamma_rate_hz']
        return cached
    with open_root(root_path) as f:
        assert all(f[k].num_entries==1 for k in ('RunInfo','WaveformEvent','SiPMEvent'))
        info=ak.to_list(f['RunInfo'].arrays())[0]
        assert info['sipm_enabled']
        if expected:
            assert info['random_seed']==expected['seed']
            assert info['alpha_p11b_rate_hz']==expected['rate_hz']
            if set(allowed_particles)=={'alpha'}:
                assert not info['brems_enabled'] and not info['c12_capture_enabled']
            if 'gamma_rate_hz' in expected:
                slot='brems' if expected['gamma_model']=='brems' else 'c12_capture'
                assert info[f'{slot}_rate_hz']==expected['gamma_rate_hz']
        config=json.loads(info['sipm_config_json'])
        s=f['SiPMEvent'].arrays(library='ak',how=dict)
        av=lambda k:ak.to_numpy(s[k][0])
        scalar=lambda k:ak.to_list(s[k])[0]
        t=av('sample_time_ns'); analog=av('analog_mV');codes=av('adc_counts')
        r=config['readout']; lsb=(r['adc_max_mV']-r['adc_min_mV'])/2**r['adc_bits']
        # The fit sees rail codes, not simulation's analog clipping truth flag.
        rails=(codes==0)|(codes==2**r['adc_bits']-1)
        np.savez_compressed(output/'waveforms'/f'{stem}.npz',time_ns=t,
            adc_mV=r['adc_min_mV']+(codes+.5)*lsb,analog_mV=analog,
            adc_rail_mask=rails)

        # This section is exclusively a scoring export, kept in a different file.
        total=f['WaveformEvent']
        pt=ak.to_numpy(total['all_photon_arrival_t_ps'].array()[0])
        px=ak.to_numpy(total['all_photon_arrival_x_mm'].array()[0])
        keys=photon_key(pt,px);order=np.argsort(keys);sorted_keys=keys[order]
        collisions=int(np.count_nonzero(sorted_keys[1:]==sorted_keys[:-1]))
        full_key=bool(collisions or np.any(pt<0) or np.any(pt>=2**32))
        hashed_join=False
        if full_key:
            del keys,order,sorted_keys
            py=ak.to_numpy(total['all_photon_arrival_y_mm'].array()[0])
            pe=ak.to_numpy(total['photon_energy_eV'].array()[0])
            keys=photon_full_key(pt,px,py,pe);del py,pe
            hashes=record_hash(keys);order=np.argsort(hashes);sorted_keys=hashes[order]
            hashed_join=not np.any(sorted_keys[1:]==sorted_keys[:-1])
            if hashed_join:full_keys=keys
            else:
                # Rare hash collisions use the original exact structured join.
                order=np.argsort(keys);sorted_keys=keys[order]
                if len(keys)>1 and np.any(sorted_keys[1:]==sorted_keys[:-1]):
                    raise ValueError('Ambiguous full photon record; ancestry cannot be assigned uniquely')
            del hashes
        del keys,pt,px
        owner=np.full(len(order),-1,dtype=np.int32)
        rows=[]
        branches=['primary_id','primary_particle_name','primary_energy_MeV','primary_time_ps',
                  'edep_total_MeV','photon_arrival_t_ps','photon_arrival_x_mm','photon_arrival_y_mm']
        if full_key:branches.append('photon_energy_eV')
        for chunk in f['PrimaryPhoton'].iterate(branches,step_size='32 MB',library='ak'):
            for p in chunk:
                particle=str(p['primary_particle_name'])
                if particle not in allowed_particles:
                    raise ValueError(f'Unexpected primary particle {particle}')
                tt=ak.to_numpy(p['photon_arrival_t_ps']);xx=ak.to_numpy(p['photon_arrival_x_mm']);yy=ak.to_numpy(p['photon_arrival_y_mm'])
                k=photon_full_key(tt,xx,yy,ak.to_numpy(p['photon_energy_eV'])) if full_key else photon_key(tt,xx)
                lookup=record_hash(k) if hashed_join else k
                loc=np.searchsorted(sorted_keys,lookup)
                assert np.all(loc<len(sorted_keys)) and np.array_equal(sorted_keys[loc],lookup)
                if hashed_join:assert np.array_equal(full_keys[order[loc]],k)
                assert np.all(owner[order[loc]]==-1)
                owner[order[loc]]=len(rows)
                inside=(xx>=config['sensor']['x_min_mm'])&(xx<config['sensor']['x_max_mm'])&(yy>=config['sensor']['y_min_mm'])&(yy<config['sensor']['y_max_mm'])
                rows.append(dict(primary_id=int(p['primary_id']),particle_name=particle,
                    emission_ns=float(p['primary_time_ps'])/1000,
                    first_active_ns=float(tt[inside].min())/1000 if inside.any() else np.nan,
                    source_energy_MeV=float(p['primary_energy_MeV']),edep_MeV=float(p['edep_total_MeV']),
                    active_photons=int(inside.sum())))
        assert np.all(owner>=0)
        del order,sorted_keys
        if hashed_join:del full_keys
        kind=av('avalanche_kind'); pi=av('avalanche_photon_index');parent=av('avalanche_parent')
        charge=av('avalanche_charge_C');at=av('avalanche_time_ns')
        parents=np.full(len(kind),-1,dtype=np.int32)
        direct=kind==0;parents[direct]=owner[pi[direct]]
        del owner
        secondary=np.flatnonzero(parent>=0)
        assert np.all(parent[secondary]<secondary)
        # Each secondary inherits its root photon's primary; dark chains stay -1.
        for j in secondary:parents[j]=parents[parent[j]]
        q0=config['sensor']['gain_electrons']*1.602176634e-19
        light=parents>=0
        q=np.bincount(parents[light],weights=charge[light],minlength=len(rows))/q0
        dark_q=float(charge[~light].sum()/q0)
        assert np.isclose(q.sum()+dark_q,charge.sum()/q0,rtol=1e-11,atol=1e-6)
        data={k:np.asarray([x[k] for x in rows]) for k in rows[0]} if rows else {
            k:np.array([]) for k in ('primary_id','emission_ns','first_active_ns','source_energy_MeV','edep_MeV','active_photons')}
        if not rows:data['particle_name']=np.array([],dtype='<U5')
        data['charge_pe']=q
        # Actual waveform and avalanche ancestry are evaluation-only quantities.
        np.savez_compressed(output/'truth'/f'{stem}.npz',**data,signal_mV=av('signal_mV'),
            avalanche_time_ns=at,avalanche_charge_C=charge,avalanche_owner=parents)
        meta=dict(root=str(root_path),root_sha256=digest(root_path),seed=int(info['random_seed']),
            rate_hz=float(info['alpha_p11b_rate_hz']),
            source_rates_hz={k:float(info[f'{k}_rate_hz']) for k in ('alpha_p11b','brems','c12_capture')},
            gate_ns=[info['analysis_start_ps']/1000,info['analysis_end_ps']/1000],
            duration_ns=info['t_length_ps']/1000,primaries=len(rows),
            photons=int(scalar('n_input_photons')),avalanche_count=len(kind),
            charge_pe=float(charge.sum()/q0),dark_charge_pe=dark_q,
            clipped_samples=int(scalar('n_clipped_samples')),rail_samples=int(rails.sum()),
            charge_balance_error_C=float(scalar('charge_balance_error_C')),
            config=config,macro_contents=info['macro_contents'],
            photon_key_mode='time_x_y_energy_exact_bits' if full_key else 'time_x_exact_bits',
            photon_time_x_collisions=collisions,
            photon_join='unique_hash_lookup_plus_exact_full_record_check' if hashed_join else 'original_exact_record_join',
            truth_definition='Realized charge in shared-cell pileup, with photon-descendant XT; excludes dark chains. Time is first active-plane photon, not source emission.')
        save(output/'metadata'/f'{stem}.json',meta)
        return meta
