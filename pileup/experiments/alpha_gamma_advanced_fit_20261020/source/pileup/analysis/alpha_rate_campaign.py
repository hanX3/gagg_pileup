"""Small native G4/S13360 pure-alpha campaign; one ROOT window per process.

Run from the project root. Frozen sources/config/macros and the runner are saved
with each new experiment; no old campaign is modified or used as a test set.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import time

PROJECT = Path(__file__).resolve().parents[2]


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def prepare(folder):
    folder.mkdir(parents=True, exist_ok=False)
    for sub in ('macros', 'logs', 'source/g4/include', 'source/g4/src', 'source/g4/tests'):
        (folder/sub).mkdir(parents=True, exist_ok=True)
    paths = [Path('g4/CMakeLists.txt'), Path('g4/main.cc')]
    paths += list(Path('g4/include').glob('*')) + list(Path('g4/src').glob('*'))
    paths += list(Path('g4/tests').glob('*.cc'))
    paths += [Path('pileup/analysis/alpha_rate_campaign.py')]
    for path in paths:
        if path.is_file():
            dest = folder/'source'/path
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dest)
    config = json.loads((PROJECT/'g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json').read_text())
    # A computational guard, not a change to sensor physics or ADC dynamic range.
    config['sensor']['max_candidates'] = 80000000
    config['configuration_notes']['rate_study'] = 'Only candidate guard raised 10M -> 80M; device/readout unchanged.'
    (folder/'source/g4/configs').mkdir()
    save(folder/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json', config)
    (folder/'git-status-before.txt').write_text(subprocess.check_output(['git','status','--short'], text=True))
    (folder/'source-changes.patch').write_bytes(subprocess.check_output(['git','diff','--binary']))
    record = {
        'created': datetime.now().astimezone().isoformat(),
        'git_head': subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
        'source_sha256': {str(p.relative_to(folder)):digest(p) for p in (folder/'source').rglob('*') if p.is_file()},
        'environment': {'python':platform.python_version(), 'platform':platform.platform(),
                        **{k:importlib.metadata.version(k) for k in ('numpy','scipy','uproot','awkward','matplotlib')}},
        'purpose':'Pure alpha pilot for blind pulse number/time/charge reconstruction; no gamma background',
        'source_model':'Current p-11B effective single-alpha spectrum; directed disk source radius 10 mm',
        'rate_definition':'Source-injected alpha particles per second, not detected events or fusion reactions',
        'max_parallel_processes':1, 'runs':[]}
    save(folder/'manifest.json', record)
    return record


def macro(folder, seed, rate, output, duration=20., guard=5.):
    config = folder/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json'
    return f'''# Pure alpha rate study; independent source, SiPM and readout seeds.
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/gagg/run/windowLength {duration:g} us
/gagg/run/preWindow {guard:g} us
/gagg/run/postWindow {guard:g} us
/gagg/run/seed {seed}
/gagg/run/outputFile {output}
/gagg/sipm/configFile {config}
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
/gagg/source/brems/enable false
/gagg/source/brems/rateHz 0
/gagg/source/c12Capture/enable false
/gagg/source/c12Capture/rateHz 0
/gagg/source/alphaP11B/enable true
/gagg/source/alphaP11B/rateHz {rate:g}
/run/beamOn 1
'''


def run_one(folder, label, rate, seed, *, duration=20., guard=5.):
    record = json.loads((folder/'manifest.json').read_text())
    if any(x['label']==label for x in record['runs']):
        raise ValueError('Run label already used; preserve its output')
    while True:
        output = PROJECT/'g4/data'/(datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root')
        if not output.exists(): break
        time.sleep(1.05)
    path = folder/'macros'/(label+'.mac')
    path.write_text(macro(folder, seed, rate, output, duration, guard))
    row = dict(label=label, rate_hz=rate, seed=seed, sipm_seed=seed+100000,
               readout_seed=seed+200000, window_us=duration, gate_us=[guard,duration-guard],
               macro=str(path.relative_to(PROJECT)), macro_sha256=digest(path),
               root=str(output.relative_to(PROJECT)), started=datetime.now().astimezone().isoformat(),
               status='running')
    record['runs'].append(row); save(folder/'manifest.json', record)
    start = time.monotonic()
    command = [str(folder/'build/gagg'),str(path)]
    log = folder/'logs'/(label+'.log')
    with log.open('w') as stream:
        result = subprocess.run(command, cwd=folder/'build', stdout=stream, stderr=subprocess.STDOUT)
    row.update(elapsed_s=time.monotonic()-start, returncode=result.returncode,
               command=command, cwd=str((folder/'build').relative_to(PROJECT)),
               log=str(log.relative_to(PROJECT)), status='completed' if result.returncode==0 else 'failed')
    if output.exists(): row.update(root_sha256=digest(output), root_bytes=output.stat().st_size)
    save(folder/'manifest.json', record)
    if result.returncode or not output.exists() or 'ROOT file is saved' not in log.read_text():
        raise RuntimeError(f'Run failed; inspect {log}')
    print(label, row['root'], f"{row['elapsed_s']:.1f} s", flush=True)
    time.sleep(max(0.,1.1-(time.monotonic()-start)))
    return row


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('folder',type=Path)
    parser.add_argument('--prepare',action='store_true')
    parser.add_argument('--label')
    parser.add_argument('--rate',type=float)
    parser.add_argument('--seed',type=int)
    parser.add_argument('--duration',type=float,default=20.)
    parser.add_argument('--guard',type=float,default=5.)
    args=parser.parse_args()
    folder=args.folder.resolve()
    if args.prepare: prepare(folder)
    else: run_one(folder,args.label,args.rate,args.seed,duration=args.duration,guard=args.guard)
