"""Export new test voltage and ancestry into separate files, without model tuning."""
from pathlib import Path
import argparse
import json
import sys
import time

FOLDER=Path(__file__).resolve().parent
PROJECT=FOLDER.parents[2]
sys.path.insert(0,str(FOLDER/'source'))
sys.path.append(str(PROJECT))
from pileup.analysis.alpha_rate_campaign import digest,save
from pileup.analysis.alpha_rate_data import export_window


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--wait',action='store_true')
    args=p.parse_args()
    while True:
        manifest=json.loads((FOLDER/'manifest.json').read_text())
        if any(r['status'] in ('failed','interrupted') for r in manifest['runs']):
            raise RuntimeError('Acquisition failed: inspect manifest')
        if len(manifest['runs'])==len(manifest['planned_runs']) and all(r['status']=='completed' for r in manifest['runs']):
            break
        if not args.wait:
            raise RuntimeError('Finish acquisition first, or use --wait')
        time.sleep(5)
    for name,expected in manifest['source_sha256'].items():
        assert digest(FOLDER/name)==expected,name
    records=[]
    for run in manifest['runs']:
        root=PROJECT/run['root']
        assert digest(root)==run['root_sha256']
        assert digest(PROJECT/run['macro'])==run['macro_sha256']
        meta=export_window(root,args.output,expected=run)
        assert meta['duration_ns']==20000. and meta['gate_ns']==[5000.,15000.]
        assert meta['macro_contents']==(PROJECT/run['macro']).read_text()
        assert abs(meta['charge_balance_error_C'])<1e-15
        records.append(dict(label=run['label'],root=run['root'],root_sha256=run['root_sha256'],
                            seed=run['seed'],rate_hz=run['rate_hz'],metadata=meta))
        save(args.output/'export-progress.json',dict(completed=len(records),planned=len(manifest['planned_runs'])))
        print('EXPORTED',run['label'],len(records),flush=True)
    save(args.output/'export-manifest.json',dict(state='completed',records=records,
        truth_policy='Supervised labels/evaluation stored separately; no clipping rejection or test model tuning'))
    print('EXPORT COMPLETE',len(records),flush=True)


if __name__=='__main__':
    main()
