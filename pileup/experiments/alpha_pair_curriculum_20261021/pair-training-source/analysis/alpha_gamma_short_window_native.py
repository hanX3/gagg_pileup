"""Frozen native held-out single-particle campaign for prefix fitting."""
from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import shutil
import subprocess

from pileup.analysis import alpha_gamma_regularized_campaign as native
from .alpha_gamma_short_window_study import digest, save

PROJECT=Path(__file__).resolve().parents[1]


def prepare(folder, analysis_output, *, only_gamma=False, gamma_repetitions=4, seed_start=261095000):
    folder=Path(folder).resolve();analysis_output=Path(analysis_output).resolve()
    if folder.exists():
        assert not (folder/'manifest.json').exists(), 'Completed preparation is immutable'
        assert not list((folder/'macros').glob('*.mac')), 'Do not resume after native generation'
    else:
        folder.mkdir(parents=True)
    parent=PROJECT/'pileup/experiments/alpha_gamma_regularized_fit_20261015'
    old=json.loads((parent/'manifest.json').read_text())
    for sub in ('macros','logs','bin','source/analysis','source/pileup/analysis'):
        (folder/sub).mkdir(parents=True,exist_ok=True)
    assert digest(parent/'bin/gagg-brems')==old['parents']['brems']['binary_sha256']
    physics_target=folder/'source/brems/g4'
    if physics_target.exists():
        for p in (parent/'source/brems/g4').rglob('*'):
            if p.is_file():assert digest(p)==digest(physics_target/p.relative_to(parent/'source/brems/g4'))
    else:
        shutil.copytree(parent/'source/brems/g4',physics_target)
    shutil.copy2(parent/'bin/gagg-brems',folder/'bin/gagg-brems')
    source_files=list((PROJECT/'analysis').glob('*.py'))
    source_files += [PROJECT/'pileup/analysis'/name for name in
        ('__init__.py','alpha_gamma_regularized_campaign.py','alpha_rate_campaign.py','alpha_rate_data.py',
         'sipm_response.py','sipm_readout.py')]
    source_files += [PROJECT/'requirements-analysis.lock.txt']
    for source in source_files:
        target=folder/'source'/source.relative_to(PROJECT)
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists() and digest(target)!=digest(source):
            draft=folder/'preparation-drafts'/source.relative_to(PROJECT)
            draft.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(target,draft)
        shutil.copy2(source,target)
    (folder/'source/pileup/__init__.py').write_text('')
    plans=[]
    specifications=([('alpha',(1.5,1.8,2.5,5.486,8.),2)] if not only_gamma else [])
    specifications += [('gamma',(.1,.2,.661657,1.5),gamma_repetitions)]
    for particle,energies,repetitions in specifications:
        for energy in energies:
            for x,y in ((0.,0.),(4.,0.),(4.,4.)):
                for repetition in range(repetitions):
                    plans.append({'label':f'test_single_{particle}_{energy:g}_{x:g}_{y:g}_{repetition:02}',
                        'role':'test_single','particle':particle,'energy_MeV':energy,'x_mm':x,'y_mm':y,
                        'gamma_model':'brems','rate_hz':0,'gamma_rate_hz':0,'seed':seed_start+len(plans)})
    for particle,energy in [('alpha',1.5),('alpha',8.),('gamma',.1),('gamma',.661657)]:
        candidates=[r for r in plans if r['particle']==particle and r['energy_MeV']==energy and r['x_mm']==r['y_mm']==0.]
        if candidates:candidates[0]['preflight']=True
    save(folder/'manifest.json',{'created':datetime.now().astimezone().isoformat(),
        'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        'purpose':'New held-out isolated native windows after prefix algorithm freeze; same source/readout physics as calibration',
        'parents':{'brems':old['parents']['brems']},'planned_runs':plans,'runs':[],
        'window_us':20.,'gate_us':[5.,15.],'max_workers':4,'preflight_validated':False,
        'analysis_output':str(analysis_output.relative_to(PROJECT)),
        'input_policy':'Full analog voltage, adequate range; all gamma zero-deposition runs retained',
        'rate_note':'One fixed primary per ROOT, not a Poisson rate experiment; rate fields unused by single mode',
        'source_sha256':{str(p.relative_to(folder)):digest(p) for p in (folder/'source').rglob('*') if p.is_file()}})
    (folder/'run_frozen.py').write_text('''from pathlib import Path
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
from analysis import alpha_gamma_short_window_native as module
module.PROJECT=folder.parents[2]
module.main()
''')
    (folder/'analyze_frozen.py').write_text('''from pathlib import Path
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
from analysis import alpha_gamma_short_window_study as module
module.PROJECT=folder.parents[2]
module.main()
''')
    print('NATIVE_PREPARED',len(plans),flush=True)


def run(folder,stage,inputs):
    folder=Path(folder).resolve();manifest=json.loads((folder/'manifest.json').read_text())
    output=PROJECT/manifest['analysis_output'];settings=json.loads((output/'frozen-settings.json').read_text())
    for path,expected in settings['source_sha256'].items():
        snapshot=folder/'source'/path
        assert digest(snapshot)==expected, str(snapshot)
    assert (output/'freeze-complete.json').exists(), 'Finish development and snapshot verification before new tests'
    if not (folder/'frozen-fit-settings.json').exists():
        save(folder/'frozen-fit-settings.json',settings)
    native.PROJECT=PROJECT
    base_macro=native.macro
    native.macro=lambda f,r,o:base_macro(f,r,o).replace('/gun/position 0 0 0 mm',
        f"/gun/position {r['x_mm']} {r['y_mm']} 0 mm")
    native.run(folder,stage,inputs,workers=4)
    manifest=json.loads((folder/'manifest.json').read_text())
    rows=[]
    for row in manifest['runs']:
        if row['status']!='completed':continue
        stem=Path(row['root']).stem
        rows.append({'root':row['root'],'seed':row['seed'],'particle':row['particle'],
            **{name:str((Path(inputs)/sub/(stem+suffix)).relative_to(PROJECT))
               for name,sub,suffix in [('waveform','waveforms','.npz'),('truth','truth','.npz'),('metadata','metadata','.json')]}})
    save(output/'test-input-records.json',rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=('prepare','preflight','test_single'))
    p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path)
    p.add_argument('--inputs',type=Path);p.add_argument('--only-gamma',action='store_true')
    p.add_argument('--gamma-repetitions',type=int,default=4);p.add_argument('--seed-start',type=int,default=261095000)
    a=p.parse_args()
    if a.stage=='prepare':prepare(a.folder,a.output,only_gamma=a.only_gamma,
        gamma_repetitions=a.gamma_repetitions,seed_start=a.seed_start)
    else:run(a.folder,a.stage,a.inputs.resolve())


if __name__=='__main__':main()
