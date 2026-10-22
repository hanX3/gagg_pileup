"""Evaluation-only outcome transitions; never imported by reconstruction."""
from collections import Counter
import argparse
from pathlib import Path
import json
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
output=parser.parse_args().output
events=json.loads((output/'events.json').read_text())
diagnostics=[]
for reference in ('baseline','hybrid_none','hybrid_both'):
    for particle in ('alpha','gamma'):
        subset=[e for e in events if e['particle']==particle]
        lost=[e for e in subset if e['outcomes'][reference]['outcome']=='recovered' and e['outcomes']['class_only']['outcome']!='recovered']
        gained=[e for e in subset if e['outcomes'][reference]['outcome']!='recovered' and e['outcomes']['class_only']['outcome']=='recovered']
        diagnostics.append(dict(reference=reference,particle=particle,lost=len(lost),gained=len(gained),
            lost_to=dict(Counter(e['outcomes']['class_only']['outcome'] for e in lost)),
            gained_from=dict(Counter(e['outcomes'][reference]['outcome'] for e in gained)),
            lost_events=[dict(stem=e['stem'],truth_index=e['truth_index']) for e in lost],
            gained_events=[dict(stem=e['stem'],truth_index=e['truth_index']) for e in gained]))

weak=[e for e in events if e['particle']=='gamma' and e['true_charge_pe'] < 1000.]
strata=[]
for group,subset in [('weak_gamma_500_1000',weak),
                     ('gamma_history_ge200',[e for e in events if e['particle']=='gamma' and e['preceding_voltage_mV']>=200.]),
                     ('gamma_nearest_q50_lt50',[e for e in events if e['particle']=='gamma' and e['nearest_other_q50_ns'] is not None and e['nearest_other_q50_ns']<50.])]:
    strata.append(dict(group=group,n=len(subset),outcomes={m:dict(Counter(e['outcomes'][m]['outcome'] for e in subset))
        for m in ('baseline','hybrid_none','hybrid_both','class_only')}))

# Final fitted Q and frozen shape-Q belong to the same all-component row.
# Their difference is observable and does not need true charge. It documents
# remaining template dependence, not a causal estimate of template bias.
shape_gap=[]
for path in sorted((output/'fits/class_only').glob('*.json')):
    fit=json.loads(path.read_text())['fit']['all_components']
    for q,shape,kind,onset in zip(fit['charge_pe'],fit['shape_charge_pe'],fit['particle_name'],fit['times_ns']):
        if q>=500.:
            shape_gap.append(dict(stem=path.stem,time_ns=onset,particle=kind,charge_pe=q,
                                  fixed_shape_charge_pe=shape,relative_gap=shape/q-1.))
shape_summary={}
for kind in ('alpha','gamma'):
    gaps=[r['relative_gap'] for r in shape_gap if r['particle']==kind]
    shape_summary[kind]=dict(n=len(gaps),absolute_gap_over_20_percent=sum(abs(x)>.2 for x in gaps),
                             gap_p16_p50_p84=np.quantile(gaps,[.16,.5,.84]).tolist() if gaps else None)
result=dict(scope='Post-fit evaluation-only transitions on disclosed development windows',
    transitions=diagnostics,strata=strata,frozen_shape_vs_final_fitted_charge=shape_summary,
    caution='Shape/fitted-Q disagreement is descriptive, not proof that changing shape improves recovery')
(output/'transition-diagnostics.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result),flush=True)
