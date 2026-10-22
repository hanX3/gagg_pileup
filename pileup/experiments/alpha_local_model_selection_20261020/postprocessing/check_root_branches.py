"""Read each required primary branch independently to localize a read error."""
import json
from pathlib import Path
import traceback
import uproot

R=Path('jupyter/run_records/alpha_local_model_selection_20261020')
root=Path('/tmp/gagg-alpha-local-read-cache/20261020_19h15m59s.root')
rows=[]
with uproot.open(root,handler=uproot.source.file.MemmapSource,num_workers=1,
                 decompression_executor=uproot.TrivialExecutor(),interpretation_executor=uproot.TrivialExecutor()) as f:
    tree=f['PrimaryPhoton']
    for name in ['primary_id','primary_particle_name','primary_energy_MeV','primary_time_ps',
                 'edep_total_MeV','photon_arrival_t_ps','photon_arrival_x_mm','photon_arrival_y_mm']:
        try:
            array=tree[name].array(library='ak')
            row=dict(branch=name,passed=True,entries=len(array))
            del array
        except Exception as error:
            row=dict(branch=name,passed=False,error=repr(error),traceback=traceback.format_exc())
        rows.append(row);print('BRANCH',name,row['passed'],flush=True)
        (R/'root-branch-diagnostic.json').write_text(json.dumps(dict(root=str(root),branches=rows),indent=2)+'\n')
