"""Frozen native alpha/gamma experiments, one complete SiPM window per ROOT.

Run from the project root. The monoenergetic variant reuses the C12 source slot
in its own source snapshot; it is explicitly NOT a C12 capture simulation.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import time

from pileup.analysis.alpha_rate_campaign import digest, save

PROJECT = Path(__file__).resolve().parents[2]


def prepare(folder, gamma_model):
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=False)
    for sub in ('source/g4/include', 'source/g4/src', 'source/g4/tests', 'source/g4/configs', 'macros', 'logs'):
        (folder/sub).mkdir(parents=True, exist_ok=True)
    paths = [Path('g4/CMakeLists.txt'), Path('g4/main.cc'),
             Path('g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json'),
             Path('pileup/analysis/alpha_gamma_campaign.py'), Path('pileup/analysis/alpha_rate_campaign.py')]
    paths += list(Path('g4/include').glob('*')) + list(Path('g4/src').glob('*')) + list(Path('g4/tests').glob('*.cc'))
    for path in paths:
        if path.is_file():
            target = folder/'source'/path
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(PROJECT/path, target)
    modifications = []
    if gamma_model == 'mono662':
        header = folder/'source/g4/include/Constants.hh'
        original = header.read_text()
        (folder/'Constants-before-mono662.hh.txt').write_text(original)
        for name, old in [('GAMMA0', '16.1'), ('GAMMA1', '11.7'), ('GAMMA_STAR', '4.439')]:
            before = f'C12_CAPTURE_{name}_ENERGY = {old}*MeV;'
            after = f'C12_CAPTURE_{name}_ENERGY = 0.661657*MeV;'
            assert original.count(before) == 1
            original = original.replace(before, after)
            modifications.append(dict(parameter=f'C12_CAPTURE_{name}_ENERGY', before_MeV=float(old),
                                      after_MeV=.661657, reason='User-requested monoenergetic gamma comparison'))
        header.write_text('// EXPERIMENT ONLY: C12 source slot repurposed as monoenergetic 661.657 keV gamma.\n'+original)
        generator = folder/'source/g4/src/PrimaryGeneratorAction.cc'
        content = generator.read_text()
        for name in ('gamma0', 'gamma1', 'gammaStar'):
            content = content.replace(f'"gamma_c12_capture_{name}"', '"gamma_mono_661657_keV"')
        generator.write_text(content)
    base = 261020000 if gamma_model == 'brems' else 261021000
    planned = [dict(label=f'{gamma_model}_alpha_{rate:08d}_{i:02d}', rate_hz=rate,
                    gamma_rate_hz=10000000, gamma_model=gamma_model, seed=base+100+group*20+i)
               for group, rate in enumerate((1000000, 10000000)) for i in range(8)]
    reference = json.loads((PROJECT/'pileup/experiments/alpha_ml_multimethod_20261014/manifest.json').read_text())
    versions = {p: importlib.metadata.version(p) for p in ('numpy','scipy','uproot','awkward','matplotlib')}
    versions.update(python=platform.python_version(), platform=platform.platform(),
                    geant4=subprocess.check_output(['geant4-config','--version'],text=True).strip(),
                    root=subprocess.check_output(['root-config','--version'],text=True).strip(),
                    g4_dataset_environment=reference['environment']['g4_dataset_environment'])
    record = dict(created=datetime.now().astimezone().isoformat(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        purpose='Blind fitting of unknown-count alpha/gamma total analog waveforms; labels only for scoring',
        gamma_model=gamma_model, gamma_slot='brems' if gamma_model=='brems' else 'c12Capture',
        gamma_description='Existing thermal brems spectrum, kTe=100 keV, 10 keV–5 MeV' if gamma_model=='brems' else
                          '661.657 keV independent photons; C12 slot repurposed in this frozen snapshot only',
        physical_modifications=modifications, alpha_model='Existing effective p–11B single-alpha spectrum',
        rate_definition='Independent source-particle injection rates; neither detected event nor reaction rates',
        window_us=20., gate_us=[5.,15.], max_parallel_processes=2, checks_validated=False,
        planned_checks=[dict(label=f'{gamma_model}_preflight', rate_hz=1000000,
                             gamma_rate_hz=10000000, gamma_model=gamma_model, seed=base)],
        planned_runs=planned, checks=[], runs=[], environment=versions,
        source_sha256={str(p.relative_to(folder)):digest(p) for p in (folder/'source').rglob('*') if p.is_file()},
        input_policy='Complete analog_mV including sample aperture/electronic noise; adequate ADC range assumed')
    (folder/'git-status-before.txt').write_text(subprocess.check_output(
        ['git','-c','core.fileMode=false','status','--short'],cwd=PROJECT,text=True))
    (folder/'source-changes.patch').write_bytes(subprocess.check_output(['git','-c','core.fileMode=false','diff','--binary'],cwd=PROJECT))
    save(folder/'manifest.json',record)
    wrapper = '''from pathlib import Path
import importlib.util
import sys
folder = Path(__file__).resolve().parent
project = folder.parents[2]
sys.path.insert(0, str(project))
sys.path.insert(0, str(folder/"source"))
spec = importlib.util.spec_from_file_location("frozen_mixed_campaign", folder/"source/pileup/analysis/alpha_gamma_campaign.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.PROJECT = project
module.main(default_folder=folder)
'''
    (folder/'run_campaign.py').write_text(wrapper)
    return record


def macro(folder, plan, output):
    seed=plan['seed']; gamma=plan['gamma_rate_hz']; alpha=plan['rate_hz']
    brems=gamma if plan['gamma_model']=='brems' else 0
    mono=gamma if plan['gamma_model']=='mono662' else 0
    return f'''# Independent alpha + {plan['gamma_model']} gamma sources; C12 slot is mono662 only in that snapshot.
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/gagg/run/windowLength 20 us
/gagg/run/preWindow 5 us
/gagg/run/postWindow 5 us
/gagg/run/seed {seed}
/gagg/run/outputFile {output}
/gagg/sipm/configFile {folder}/source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json
/gagg/sipm/seed {seed+100000}
/gagg/sipm/readoutSeed {seed+200000}
/gagg/sipm/enable true
/gagg/step/setGaggMaxStep 1 um
/gagg/step/setMylarFrontMaxStep 0.5 um
/gagg/step/setMylarSideMaxStep 1 um
/run/initialize
/gagg/source/useDiskConeSource true
/gagg/source/aimAtGagg true
/gagg/source/diskRadius 10 mm
/gagg/pileup/enable true
/gagg/source/brems/enable {'true' if brems else 'false'}
/gagg/source/brems/rateHz {brems}
/gagg/source/c12Capture/enable {'true' if mono else 'false'}
/gagg/source/c12Capture/rateHz {mono}
/gagg/source/alphaP11B/enable true
/gagg/source/alphaP11B/rateHz {alpha}
/run/beamOn 1
'''


def run_stage(folder, stage, workers):
    folder=Path(folder).resolve(); record=json.loads((folder/'manifest.json').read_text())
    if not 1 <= workers <= record['max_parallel_processes']: raise ValueError('Invalid workers')
    if stage=='production' and not record['checks_validated']: raise RuntimeError('Validate checks first')
    key='checks' if stage=='checks' else 'runs'
    existing={r['label']:r for r in record[key]}
    if any(r['status']!='completed' for r in existing.values()): raise RuntimeError('Inspect incomplete runs')
    pending=[p for p in record['planned_checks' if key=='checks' else 'planned_runs'] if p['label'] not in existing]
    binary=folder/'build/gagg'
    if digest(binary)!=record['binary_sha256']: raise RuntimeError('Binary changed')
    environment=dict(os.environ);environment.update(record['environment']['g4_dataset_environment'])
    active=[];last=-float('inf')
    def persist():
        save(folder/'manifest.json.tmp',record); (folder/'manifest.json.tmp').replace(folder/'manifest.json')
    try:
        while pending or active:
            if pending and len(active)<workers and time.monotonic()-last>=1.05:
                output=PROJECT/'g4/data'/(datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root')
                if output.exists(): time.sleep(1.05); continue
                plan=pending.pop(0); path=folder/'macros'/(plan['label']+'.mac');log=folder/'logs'/(plan['label']+'.log')
                if path.exists() or log.exists(): raise RuntimeError('Existing run artifacts')
                path.write_text(macro(folder,plan,output));stream=log.open('w')
                command=[str(binary),str(path)]
                start=time.monotonic();proc=subprocess.Popen(command,cwd=folder/'build',env=environment,
                                    stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
                row=dict(plan,sipm_seed=plan['seed']+100000,readout_seed=plan['seed']+200000,
                    window_us=20.,gate_us=[5.,15.],root=str(output.relative_to(PROJECT)),
                    macro=str(path.relative_to(PROJECT)),macro_sha256=digest(path),log=str(log.relative_to(PROJECT)),
                    command=command,cwd=str((folder/'build').relative_to(PROJECT)),pid=proc.pid,
                    started=datetime.now().astimezone().isoformat(),status='running')
                record[key].append(row);persist();active.append((proc,stream,row,start,output));last=start
                print('START',row['label'],row['root'],flush=True)
            for task in list(active):
                proc,stream,row,start,output=task; code=proc.poll()
                if code is None:continue
                stream.close();text=(PROJECT/row['log']).read_text()
                write_errors=[s for s in ('Error in <TFile::WriteBuffer>', 'Error in <TBranchElement::WriteBasketImpl>',
                                          'Error in <TTree::Fill>') if s in text]
                ok=(code==0 and output.exists() and 'ROOT file is saved' in text and
                    'COMMAND NOT FOUND' not in text and not write_errors)
                row.update(returncode=code,elapsed_s=time.monotonic()-start,
                    finished=datetime.now().astimezone().isoformat(),status='completed' if ok else 'failed')
                if write_errors:row['root_write_errors']=write_errors
                if output.exists():row.update(root_sha256=digest(output),root_bytes=output.stat().st_size)
                active.remove(task);persist();print('DONE',row['label'],row['status'],round(row['elapsed_s'],1),flush=True)
                if not ok:raise RuntimeError(f"Failed run: {row['log']}")
            if pending or active:time.sleep(.5)
    except BaseException:
        for proc,stream,row,start,output in active:
            if proc.poll() is None:os.killpg(proc.pid,signal.SIGTERM);proc.wait()
            stream.close();row.update(status='interrupted',returncode=proc.returncode,elapsed_s=time.monotonic()-start)
        persist();raise
    print('STAGE COMPLETE',stage,len(record[key]),flush=True)


def main(default_folder=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',type=Path,default=default_folder,required=default_folder is None)
    parser.add_argument('--prepare',choices=('brems','mono662'))
    parser.add_argument('--stage',choices=('checks','production'))
    parser.add_argument('--workers',type=int,default=1)
    args=parser.parse_args()
    if args.prepare:prepare(args.folder,args.prepare)
    elif args.stage:run_stage(args.folder,args.stage,args.workers)
    else:parser.error('Select --prepare or --stage')


if __name__=='__main__':main()
