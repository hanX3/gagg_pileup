"""Post-inference diagnostic of teacher anchors and the reported dense peaks.

This is an oracle diagnostic of saved outputs, never a reconstruction method.
No network, fit, threshold, scoring pair or recorded prediction is changed.
Run from the project root, or supply --project explicitly.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(project):
    project=Path(project).resolve()
    out=project/'jupyter/artifacts/alpha_component_ml_20261021'
    summary=json.loads((out/'test-summary.json').read_text())
    assert len(summary['windows'])==60
    methods=[f'{v}_{s}' for v in ('gru_base','gru_component','unetgru_base','unetgru_component')
             for s in (262041001,262041002)]
    groups={name:dict(target_count=0,anchor_q_qualified=0,matched_count=0,
                     final_q_qualified=0,matched_anchor_q_qualified=0,
                     final_q_failed_but_anchor_q_qualified=0,
                     final_q_qualified_but_anchor_q_failed=0,
                     reported_peak_at_nearest_true_bin=0,
                     anchor_relative_q_errors=[],final_relative_q_errors=[],
                     matched_anchor_relative_q_errors=[],peak_bin_delta=[])
            for name in methods}
    windows=[]
    for row in summary['windows']:
        with np.load(project/row['truth']) as f:
            truth_t=f['first_active_ns'].copy();truth_q=f['charge_pe'].copy()
        with np.load(project/row['waveform']) as f:
            time=f['time_ns'].copy()
        stem=Path(row['root']).stem
        for name in methods:
            result=row['methods'][name];fit=result['fit'];score=result['score']
            with np.load(out/'test/dense'/f'{stem}_{name}.npz') as f:
                probability=f['probability'].copy();dense_q=f['charge_pe'].copy()
                offset=f['time_offset_ns'].copy()
            peaks,_=find_peaks(probability,height=fit['threshold'],distance=4)
            peaks=peaks[dense_q[peaks]>=500.]
            order=np.argsort(time[peaks]+offset[peaks]);peaks=peaks[order]
            assert np.array_equal((time[peaks]+offset[peaks]),fit['times_ns'])
            assert np.array_equal(dense_q[peaks].astype(float),fit['charge_pe'])
            ids=np.asarray(score['target_indices'],dtype=int)
            anchors=np.rint((truth_t[ids]-time[0])/2.).astype(int)
            assert np.all((anchors>=0)&(anchors<len(time)))
            anchor_error=dense_q[anchors]/truth_q[ids]-1.
            pairs=np.asarray(score['matched_pairs'],dtype=int).reshape(-1,2)
            matched_ids,js=pairs.T
            matched_anchor=np.rint((truth_t[matched_ids]-time[0])/2.).astype(int)
            matched_anchor_error=dense_q[matched_anchor]/truth_q[matched_ids]-1.
            final_error=np.asarray(fit['charge_pe'])[js]/truth_q[matched_ids]-1.
            at_anchor=abs(matched_anchor_error)<=.2;final=abs(final_error)<=.2
            assert int(final.sum())==score['recovered_count']
            delta=peaks[js]-matched_anchor
            counts=dict(target_count=len(ids),anchor_q_qualified=int((abs(anchor_error)<=.2).sum()),
                        matched_count=len(js),final_q_qualified=int(final.sum()),
                        matched_anchor_q_qualified=int(at_anchor.sum()),
                        final_q_failed_but_anchor_q_qualified=int((~final&at_anchor).sum()),
                        final_q_qualified_but_anchor_q_failed=int((final&~at_anchor).sum()),
                        reported_peak_at_nearest_true_bin=int((delta==0).sum()))
            g=groups[name]
            for key,value in counts.items():g[key]+=value
            g['anchor_relative_q_errors'].extend(anchor_error.tolist())
            g['final_relative_q_errors'].extend(final_error.tolist())
            g['matched_anchor_relative_q_errors'].extend(matched_anchor_error.tolist())
            g['peak_bin_delta'].extend(delta.tolist())
            windows.append(dict(seed=row['seed'],rate_hz=row['rate_hz'],method=name,**counts))
    rows=[]
    for name,g in groups.items():
        for key in ('anchor_relative_q_errors','matched_anchor_relative_q_errors',
                    'final_relative_q_errors','peak_bin_delta'):
            values=g.pop(key)
            g[key+'_p16_p50_p84']=np.quantile(values,[.16,.5,.84]).tolist() if values else None
        rows.append(dict(method=name,**g))
    record=dict(status='completed',scope='Exploratory oracle diagnostic only; no inference or scoring changes',
                diagnostic_source_sha256=sha(__file__),test_summary_sha256=sha(out/'test-summary.json'),
                truth_use='Known nearest true bin after saved blind prediction; not achievable blind recovery',
                all_test_targets_retained=True,windows=windows,totals=rows)
    path=out/'anchor-head-diagnostic.json'
    assert not path.exists()
    path.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    print('ANCHOR_HEAD_DIAGNOSTIC',len(windows),'saved cases; inference unchanged',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',default='.')
    main(parser.parse_args().project)
