"""Exact-field scoring recovery must be bijective or refuse the input."""
import unittest
from unittest.mock import patch
import zlib

import awkward as ak
import numpy as np

from .root_redundant_photon_fields import FIELDS, recover, _Root


class Branch:
    def __init__(self, values, bad=False): self.values, self.bad = values, bad
    def array(self):
        if self.bad: raise zlib.error('Error -5 while decompressing data: incomplete or truncated stream')
        return self.values


class Tree:
    def __init__(self, data, bad): self.data, self.bad = data, bad
    @property
    def num_entries(self): return len(next(iter(self.data.values())))
    def __getitem__(self, name): return Branch(self.data[name], name in self.bad)
    def iterate(self, names, **kwargs):
        yield ak.Array({name:self.data[name] for name in names})


def fixture(global_bad=(), primary_bad=(), ambiguous=False, mismatch=False):
    values = dict(time=np.array([10,20,30,40,50,60], dtype=np.int64),
        x=np.array([.1,.2,.3,.4,.5,.6],dtype=np.float32),
        y=np.array([.7,.8,.9,1.,1.1,1.2],dtype=np.float32),
        energy=np.array([2.,2.1,2.2,2.3,2.4,2.5],dtype=np.float32))
    if ambiguous:
        values['y'][1]=values['y'][0];values['energy'][1]=values['energy'][0]
    order=np.array([5,1,3,0,2,4]);lengths=np.array([3,0,3])
    global_data={g:ak.unflatten(values[name],np.array([6])) for name,g,_,_ in FIELDS}
    primary_data={p:ak.unflatten(values[name][order].copy(),lengths) for name,_,p,_ in FIELDS}
    if mismatch:
        a=ak.to_numpy(ak.flatten(primary_data['photon_arrival_y_mm'])).copy();a[0]+=1
        primary_data['photon_arrival_y_mm']=ak.unflatten(a,lengths)
    root={'WaveformEvent':Tree(global_data,set(global_bad)), 'PrimaryPhoton':Tree(primary_data,set(primary_bad))}
    return root,values,order


class RecoveryTests(unittest.TestCase):
    def test_missing_global_coordinate(self):
        root,values,_=fixture(global_bad=('all_photon_arrival_x_mm',));checks={}
        view=_Root(root,recover(root,checks))
        np.testing.assert_array_equal(ak.to_numpy(view['WaveformEvent']['all_photon_arrival_x_mm'].array()[0]),values['x'])
        self.assertTrue(checks['one_to_one_complete_coverage'])

    def test_missing_primary_time(self):
        root,values,order=fixture(primary_bad=('photon_arrival_t_ps',));checks={}
        view=_Root(root,recover(root,checks));chunk=next(view['PrimaryPhoton'].iterate(['photon_arrival_t_ps','photon_arrival_x_mm']))
        np.testing.assert_array_equal(ak.to_numpy(ak.flatten(chunk['photon_arrival_t_ps'])),values['time'][order])
        np.testing.assert_array_equal(ak.to_numpy(ak.num(chunk['photon_arrival_t_ps'])),[3,0,3])

    def test_two_missing_fields_unique_common_records(self):
        root,values,order=fixture(global_bad=('all_photon_arrival_x_mm',),primary_bad=('photon_arrival_t_ps',));checks={}
        recovered=recover(root,checks)
        self.assertEqual(checks['common_exact_fields'],['y','energy'])
        np.testing.assert_array_equal(ak.to_numpy(recovered['WaveformEvent']['all_photon_arrival_x_mm'][0]),values['x'])
        np.testing.assert_array_equal(ak.to_numpy(ak.flatten(recovered['PrimaryPhoton']['photon_arrival_t_ps'])),values['time'][order])

    def test_ambiguous_common_fields_refused(self):
        root,_,_=fixture(global_bad=('all_photon_arrival_x_mm',),primary_bad=('photon_arrival_t_ps',),ambiguous=True)
        with self.assertRaisesRegex(AssertionError,'Non-unique'):recover(root,{})

    def test_mismatching_photon_refused(self):
        root,_,_=fixture(global_bad=('all_photon_arrival_x_mm',),mismatch=True)
        with self.assertRaises(AssertionError):recover(root,{})

    def test_both_copies_unreadable_or_hash_collision_refused(self):
        root,_,_=fixture(global_bad=('all_photon_arrival_t_ps',),primary_bad=('photon_arrival_t_ps',))
        with self.assertRaisesRegex(AssertionError,'Both copies'):recover(root,{})
        root,_,_=fixture(global_bad=('all_photon_arrival_x_mm',))
        with patch('pileup.analysis.root_redundant_photon_fields.record_hash',side_effect=lambda x:np.zeros(len(x),dtype=np.uint64)):
            with self.assertRaisesRegex(AssertionError,'Non-unique'):recover(root,{})


if __name__=='__main__':unittest.main()
