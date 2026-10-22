"""Harvest completed runs after an I/O exception; retry the same failed seeds."""
from pathlib import Path
from datetime import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

PROJECT=next(p for p in Path(__file__).resolve().parents if (p/'g4/data').is_dir()
             and (p/'pileup/experiments/alpha_pair_resolution_20261022').is_dir())
EXPERIMENT=PROJECT/'pileup/experiments/alpha_pair_resolution_20261022'
for key in ('GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(PROJECT)
sys.path.insert(0,str(EXPERIMENT/'replay-source'))


def main():
    import awkward as ak
    from pileup.analysis.alpha_pair_resolution_data import OUTPUT,RECORD,load,save,digest,check_freeze
    from pileup.analysis.alpha_rate_data import export_window,open_root
    from pileup.analysis.alpha_gamma_regularized_campaign import macro
    from pileup.analysis.sipm_response import read_photon_hits
    check_freeze()
    manifest=load(EXPERIMENT/'manifest.json')
    backup=RECORD/'manifest.before_native_recovery.json'
    assert not backup.exists(), 'Recovery already attempted; inspect saved state before resuming.'
    shutil.copyfile(EXPERIMENT/'manifest.json',backup)
    old={r['label']:r for r in manifest['runs']}
    assert len(list((EXPERIMENT/'macros').glob('alpha_*.mac')))==240
    cache=Path(tempfile.mkdtemp(prefix='gagg-resolution-recovery-'))
    binary=cache/'gagg-brems';shutil.copyfile(EXPERIMENT/'bin/gagg-brems',binary);binary.chmod(0o755)
    env=os.environ.copy();env.update(manifest['environment']['g4_dataset_environment'])
    attempts=[];recovered=[];last_launch=0.

    def validate(plan,root,macro_path,log_path,started,export_folder):
        body=macro_path.read_text();log=log_path.read_text()
        assert 'ROOT file is saved.' in log
        assert 'Error in <TFile' not in log, 'ROOT reported a write error.'
        local=cache/root.name;shutil.copyfile(root,local)
        sha=digest(root);assert sha==digest(local)
        with open_root(local) as f:
            info=ak.to_list(f['RunInfo'].arrays(library='ak'))[0]
            assert f['WaveformEvent'].num_entries==f['SiPMEvent'].num_entries==1
            assert info['macro_contents']==body and info['random_seed']==plan['seed']
            assert info['sipm_avalanche_base_seed']==plan['seed']+100000
            assert info['sipm_readout_base_seed']==plan['seed']+200000
        # Read all donor photon arrays now, before admitting this ROOT to tests.
        hits=read_photon_hits(local);assert len(hits.time_ns)>0
        folder=OUTPUT/export_folder
        metadata=folder/'metadata'/(root.stem+'.json')
        if metadata.exists():
            meta=load(metadata)
            assert meta['root_sha256']==sha
            for part in ('waveforms','truth'):
                assert (folder/part/(root.stem+'.npz')).exists()
        else:
            meta=export_window(local,folder,allowed_particles=('alpha',))
            meta.update(root=str(root.relative_to(PROJECT)),io_note='Hash-verified regular local recovery copy; no links.')
            save(metadata,meta)
        assert meta['primaries']==1 and meta['root_sha256']==sha
        # Save a separate provenance record; do not overwrite an existing export.
        save(OUTPUT/'recovered-provenance'/(root.stem+'.json'),dict(root=str(root.relative_to(PROJECT)),
             root_sha256=sha,macro=str(macro_path.relative_to(PROJECT)),macro_contents_verified=True,
             all_photon_arrays_read=True,io_note='Hash-verified regular local copy, no links.'))
        return dict(**plan,root=str(root.relative_to(PROJECT)),started=started,
            started_basis='Original macro creation time / retry Popen time',
            sipm_seed=plan['seed']+100000,readout_seed=plan['seed']+200000,
            macro=str(macro_path.relative_to(PROJECT)),macro_sha256=digest(macro_path),
            log=str(log_path.relative_to(PROJECT)),status='completed',returncode=0,
            root_sha256=sha,root_bytes=root.stat().st_size,
            finished=datetime.now().astimezone().isoformat(),export_folder=export_folder,
            validation=dict(runinfo_macro_verified=True,seed_verified=True,all_photon_arrays_read=True,
                            charge_balance_error_C=meta['charge_balance_error_C']))

    for plan in manifest['planned_runs']:
        if plan['label'] in old:
            recovered.append(old[plan['label']]);continue
        macro_path=EXPERIMENT/'macros'/(plan['label']+'.mac')
        log_path=EXPERIMENT/'logs'/(plan['label']+'.log')
        root=Path(re.search(r'^/gagg/run/outputFile (.+)$',macro_path.read_text(),re.M)[1])
        started=datetime.fromtimestamp(macro_path.stat().st_mtime).astimezone().isoformat()
        try:
            row=validate(plan,root,macro_path,log_path,started,'native-inputs')
        except Exception as exc:
            attempts.append(dict(label=plan['label'],root=str(root.relative_to(PROJECT)),
                                 macro=str(macro_path.relative_to(PROJECT)),log=str(log_path.relative_to(PROJECT)),
                                 retained=True,reason=repr(exc),seed=plan['seed']))
            save(RECORD/'native-recovery-attempts.json',attempts)
            # Same seed and configuration, serial retry, new second-resolution ROOT.
            while True:
                remaining=1.05-(time.monotonic()-last_launch)
                if remaining>0:time.sleep(remaining)
                root=PROJECT/'g4/data'/(datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root')
                if not root.exists():break
                time.sleep(1.05)
            macro_path=EXPERIMENT/'macros'/('recovery_'+plan['label']+'_01.mac')
            log_path=EXPERIMENT/'logs'/('recovery_'+plan['label']+'_01.log')
            assert not macro_path.exists() and not log_path.exists()
            macro_path.write_text(macro(EXPERIMENT,dict(plan,root=str(root.relative_to(PROJECT))),root))
            started=datetime.now().astimezone().isoformat();last_launch=time.monotonic()
            with log_path.open('w') as out:
                subprocess.run([str(binary),str(macro_path)],cwd=PROJECT/'g4/build',
                               env=env,stdout=out,stderr=subprocess.STDOUT,check=True)
            row=validate(plan,root,macro_path,log_path,started,'recovery-inputs')
            row['same_seed_retry']=True;row['original_attempt']=attempts[-1]
        recovered.append(row)
        manifest['runs']=recovered+list(old[k] for k in old if k not in {r['label'] for r in recovered})
        manifest['retained_failed_attempts']=attempts
        save(EXPERIMENT/'manifest.json',manifest)
        if len(recovered)%10==0:print('RECOVERED',len(recovered),'/ 240',flush=True)
    assert len(recovered)==240 and len({r['seed'] for r in recovered})==240
    manifest['runs']=recovered;manifest['retained_failed_attempts']=attempts
    save(EXPERIMENT/'manifest.json',manifest)
    records=[]
    for row in sorted(recovered,key=lambda r:r['seed']):
        stem=Path(row['root']).stem;folder=OUTPUT/row.get('export_folder','native-inputs')
        records.append(dict(root=row['root'],seed=row['seed'],energy_MeV=row['energy_MeV'],
            replicate=row['replicate'],particle='alpha',root_original_sha256=row['root_sha256'],
            source_manifest=str((EXPERIMENT/'manifest.json').relative_to(PROJECT)),
            waveform=str((folder/'waveforms'/(stem+'.npz')).relative_to(PROJECT)),
            truth=str((folder/'truth'/(stem+'.npz')).relative_to(PROJECT)),
            metadata=str((folder/'metadata'/(stem+'.json')).relative_to(PROJECT))))
    save(OUTPUT/'native-records.json',dict(test_single=records))
    save(RECORD/'native-recovery.json',dict(passed=True,roots=240,harvested=len(recovered)-len(old),
        original_recorded=len(old),retries=len(attempts),all_original_attempts_retained=True,
        seeds_or_physics_changed=False,reason='Original future exception exited collection; other submitted G4 jobs completed.'))
    print('RECOVERY_COMPLETE',len(records),'same-seed retries',len(attempts),flush=True)


if __name__=='__main__':main()
