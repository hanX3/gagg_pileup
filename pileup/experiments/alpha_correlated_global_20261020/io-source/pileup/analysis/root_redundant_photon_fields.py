"""Exact recovery of unreadable duplicated photon fields, for scoring only.

Global WaveformEvent and PrimaryPhoton store the same photon records. Join on
every readable common field, require unique global keys and complete bijection,
and copy missing values from the other representation. No time/position/owner
is guessed. The original ROOT and measured SiPM voltage are never changed.
"""
from contextlib import contextmanager

import awkward as ak
import numpy as np


FIELDS = (
    ('time', 'all_photon_arrival_t_ps', 'photon_arrival_t_ps', '<i8'),
    ('x', 'all_photon_arrival_x_mm', 'photon_arrival_x_mm', '<u4'),
    ('y', 'all_photon_arrival_y_mm', 'photon_arrival_y_mm', '<u4'),
    ('energy', 'photon_energy_eV', 'photon_energy_eV', '<u4'),
)


def record_hash(records):
    """Lookup accelerator; every selected record is checked for exact equality."""
    h = np.zeros(len(records), dtype=np.uint64)
    for i, name in enumerate(records.dtype.names):
        v = records[name].astype(np.uint64)
        v += np.uint64(0x9E3779B97F4A7C15) ^ np.uint64(i)
        v = (v ^ (v >> np.uint64(30)))*np.uint64(0xBF58476D1CE4E5B9)
        v = (v ^ (v >> np.uint64(27)))*np.uint64(0x94D049BB133111EB)
        h ^= v ^ (v >> np.uint64(31))
        h = (h << np.uint64(17)) | (h >> np.uint64(47))
    return h


def key(values, common):
    result = np.empty(len(values[common[0]]), dtype=[(name, dtype) for name, _, _, dtype in FIELDS if name in common])
    for name in result.dtype.names:
        result[name] = values[name] if name == 'time' else np.asarray(values[name], np.float32).view(np.uint32)
    return result


def recover(root, checks):
    global_values = {}; primary_values = {}; lengths = None
    failures = {}; n = None
    for name, global_branch, primary_branch, _ in FIELDS:
        try:
            values = ak.to_numpy(root['WaveformEvent'][global_branch].array()[0])
            global_values[name] = values
            if n is None: n = len(values)
            assert len(values) == n
        except Exception as error:
            if 'decompress' not in str(error).lower() and 'truncated stream' not in str(error).lower(): raise
            failures['WaveformEvent/'+global_branch] = repr(error)
        try:
            jagged = root['PrimaryPhoton'][primary_branch].array()
            counts = ak.to_numpy(ak.num(jagged))
            if lengths is None: lengths = counts
            else: np.testing.assert_array_equal(lengths, counts)
            primary_values[name] = ak.to_numpy(ak.flatten(jagged))
        except Exception as error:
            if 'decompress' not in str(error).lower() and 'truncated stream' not in str(error).lower(): raise
            failures['PrimaryPhoton/'+primary_branch] = repr(error)
    assert failures, 'Recovery is only permitted for an observed decoding failure'
    assert n is not None and lengths is not None and int(lengths.sum()) == n
    assert all(len(a) == n for a in primary_values.values())
    common = [name for name, *_ in FIELDS if name in global_values and name in primary_values]
    assert common, 'No readable common photon fields; exact recovery refused'
    for name, *_ in FIELDS:
        assert name in global_values or name in primary_values, f'Both copies unreadable: {name}; recovery refused'
    gkey = key(global_values, common); pkey = key(primary_values, common)
    hashes = record_hash(gkey); order = np.argsort(hashes); ordered = hashes[order]
    # Hash collisions are rejected here, not silently resolved. Equality of all
    # common stored fields is checked after lookup, in addition to uniqueness.
    assert not np.any(ordered[1:] == ordered[:-1]), 'Non-unique/colliding common keys; recovery refused'
    locations = np.searchsorted(ordered, record_hash(pkey))
    assert np.all(locations < n), 'Missing photon key; recovery refused'
    indices = order[locations]
    assert np.array_equal(gkey[indices], pkey), 'Inexact photon join; recovery refused'
    used = np.zeros(n, dtype=bool); used[indices] = True
    assert used.all(), 'Duplicate/incomplete primary coverage; recovery refused'
    # Equal source/destination counts plus complete coverage implies a bijection.
    global_recovered = {}; primary_recovered = {}
    for name, global_branch, primary_branch, _ in FIELDS:
        if name not in global_values:
            values = np.empty(n, dtype=primary_values[name].dtype)
            values[indices] = primary_values[name]
            global_recovered[global_branch] = ak.unflatten(values, np.array([n], dtype=np.int64))
        if name not in primary_values:
            primary_recovered[primary_branch] = ak.unflatten(global_values[name][indices], lengths)
    checks.update(failures=failures, common_exact_fields=common, global_photons=n,
                  primary_photons=int(lengths.sum()), unique_lookup_keys=True,
                  full_stored_common_fields_equal=True, one_to_one_complete_coverage=True,
                  recovered_global_fields=list(global_recovered), recovered_primary_fields=list(primary_recovered),
                  inference_or_model_inputs_used=False)
    return {'WaveformEvent':global_recovered, 'PrimaryPhoton':primary_recovered}


class _Branch:
    def __init__(self, values): self.values = values

    def array(self, *args, **kwargs):
        assert not args and kwargs.get('library', 'ak') == 'ak'
        return self.values[kwargs.get('entry_start', 0):kwargs.get('entry_stop', None)]


class _Tree:
    def __init__(self, original, replacements): self.original, self.replacements = original, replacements
    def __getattr__(self, name): return getattr(self.original, name)
    def __getitem__(self, name):
        return _Branch(self.replacements[name]) if name in self.replacements else self.original[name]

    def iterate(self, branches, **options):
        read = [name for name in branches if name not in self.replacements]
        offset = 0
        for chunk in self.original.iterate(read, **options):
            count = len(chunk)
            for name in branches:
                if name in self.replacements:
                    chunk = ak.with_field(chunk, self.replacements[name][offset:offset+count], name)
            yield chunk
            offset += count
        assert offset == self.original.num_entries


class _Root:
    def __init__(self, original, replacements): self.original, self.replacements = original, replacements
    def __getattr__(self, name): return getattr(self.original, name)
    def __getitem__(self, name):
        return _Tree(self.original[name], self.replacements[name]) if name in self.replacements else self.original[name]


@contextmanager
def exact_redundant_fields(export_function, checks):
    """Adapt only ROOT field reads; keep the export and measured voltage intact."""
    namespace = export_function.__globals__; original = namespace['open_root']
    @contextmanager
    def opened(path):
        with original(path) as root:
            replacements = recover(root, checks)
            yield _Root(root, replacements)
    namespace['open_root'] = opened
    try: yield
    finally: namespace['open_root'] = original
