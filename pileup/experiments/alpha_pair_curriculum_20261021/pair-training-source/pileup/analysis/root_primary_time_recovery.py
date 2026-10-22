"""Recover a duplicated primary-photon time field by an exact record join.

Used only for scoring exports when one ROOT branch cannot be decompressed.
All global photon records must have unique stored (x, y, energy) float bits;
each primary record must map exactly once and cover the global record set.
No pulse time, charge, particle count or fit result is inferred by this adapter.
"""
from contextlib import contextmanager

import awkward as ak
import numpy as np


def coordinate_energy_key(x,y,energy):
    key=np.empty(len(x),dtype=[('x','<u4'),('y','<u4'),('energy','<u4')])
    for name,values in zip(key.dtype.names,(x,y,energy)):
        key[name]=np.asarray(values,dtype=np.float32).view(np.uint32)
    return key


def lookup_hash(key):
    """Fast lookup only: hash uniqueness and exact full records are checked."""
    value=key['x'].astype(np.uint64)|(key['y'].astype(np.uint64)<<np.uint64(32))
    value^=key['energy'].astype(np.uint64)*np.uint64(0x9E3779B97F4A7C15)
    value=(value^(value>>np.uint64(30)))*np.uint64(0xBF58476D1CE4E5B9)
    value=(value^(value>>np.uint64(27)))*np.uint64(0x94D049BB133111EB)
    return value^(value>>np.uint64(31))


class PrimaryTimeView:
    def __init__(self,tree,total,checks):
        self.tree,self.total,self.checks=tree,total,checks

    def __getattr__(self,name):return getattr(self.tree,name)

    def iterate(self,branches,**options):
        read=[name for name in branches if name!='photon_arrival_t_ps']
        if 'photon_energy_eV' not in read:read.append('photon_energy_eV')
        get=lambda name:ak.to_numpy(self.total[name].array()[0])
        times=get('all_photon_arrival_t_ps')
        keys=coordinate_energy_key(get('all_photon_arrival_x_mm'),
            get('all_photon_arrival_y_mm'),get('photon_energy_eV'))
        hashes=lookup_hash(keys);order=np.argsort(hashes);sorted_hashes=hashes[order];del hashes
        assert not np.any(sorted_hashes[1:]==sorted_hashes[:-1]),'Ambiguous/colliding lookup keys; recovery refused'
        used=np.zeros(len(order),dtype=bool);count=0
        for chunk in self.tree.iterate(read,**options):
            lengths=ak.to_numpy(ak.num(chunk['photon_arrival_x_mm']))
            flat=lambda name:ak.to_numpy(ak.flatten(chunk[name]))
            incoming=coordinate_energy_key(flat('photon_arrival_x_mm'),flat('photon_arrival_y_mm'),flat('photon_energy_eV'))
            where=np.searchsorted(sorted_hashes,lookup_hash(incoming))
            assert np.all(where<len(sorted_hashes)),'Missing photon record; recovery refused'
            indices=order[where]
            assert np.array_equal(keys[indices],incoming),'Inexact photon record; recovery refused'
            assert len(np.unique(indices))==len(indices) and not used[indices].any(),'Repeated photon ownership; recovery refused'
            used[indices]=True;count+=len(indices)
            recovered=ak.unflatten(times[indices],lengths)
            yield ak.with_field(chunk,recovered,'photon_arrival_t_ps')
        assert used.all() and count==len(times),'Incomplete photon coverage; recovery refused'
        self.checks.update(global_photons=len(times),mapped_primary_photons=count,
            exact_unique_float_bits=True,one_to_one_complete_coverage=True)


class RootView:
    def __init__(self,root,checks):self.root,self.checks=root,checks
    def __getattr__(self,name):return getattr(self.root,name)
    def __getitem__(self,name):
        if name=='PrimaryPhoton':return PrimaryTimeView(self.root[name],self.root['WaveformEvent'],self.checks)
        return self.root[name]


@contextmanager
def scoring_time_recovery(export_function,checks):
    """Replace the file view only; leave frozen export/inference source intact."""
    globals_=export_function.__globals__;original=globals_['open_root']
    @contextmanager
    def open_view(path):
        with original(path) as root:yield RootView(root,checks)
    globals_['open_root']=open_view
    try:yield
    finally:globals_['open_root']=original
