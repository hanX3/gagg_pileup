"""Post-fit charge-prior diagnostics; truth is accessed only for evaluation."""
from pathlib import Path
from collections import Counter, defaultdict
import argparse
import json
import numpy as np


def run(project, output):
    project, output = Path(project), Path(output)
    events = json.loads((output/'fixed-association-events.json').read_text())
    common = {p.stem:json.loads(p.read_text()) for p in (output/'common').glob('*.json')}
    neural = {stem:json.loads((output/'blind-inputs'/stem/'neural.json').read_text())['fit']['all_components'] for stem in common}
    transitions = []
    raw = []
    for particle in ('alpha','gamma'):
        subset = [e for e in events if e['particle'] == particle]
        for arm in ('none','charge','class','both'):
            lost = [e for e in subset if e['outcomes']['none']['outcome']=='recovered' and e['outcomes'][arm]['outcome']!='recovered']
            gained = [e for e in subset if e['outcomes']['none']['outcome']!='recovered' and e['outcomes'][arm]['outcome']=='recovered']
            transitions.append(dict(particle=particle,arm=arm,lost=len(lost),gained=len(gained),
                                    loss_outcomes=dict(Counter(e['outcomes'][arm]['outcome'] for e in lost))))
            raw.append(dict(particle=particle,arm=arm,targets=len(subset),
                class_wrong_at_any_charge=sum(e['outcomes'][arm]['prediction'] is not None and
                    e['outcomes'][arm]['prediction']['particle'] != particle for e in subset),
                no_candidate=sum(e['outcomes'][arm]['prediction'] is None for e in subset)))
    charge_errors=defaultdict(list)
    for e in events:
        pred=e['outcomes']['none']['prediction']
        if pred is None or pred['neural_prior_charge_pe'] <= 0:continue
        t=common[e['stem']]['times_ns'][pred['candidate_index']]
        nn=neural[e['stem']];m=int(np.argmin(abs(np.asarray(nn['times_ns'])-t)))
        kind=nn['particle_name'][m]
        charge_errors[(e['particle'],kind,e['strata']['charge'])].append(pred['neural_prior_charge_pe']/e['true_charge_pe']-1.)
    charge_groups=[dict(true_particle=p,decoded_neural_particle=k,charge_bin=b,n=len(v),
        median_relative_q_error=float(np.median(v)),p16_p84=np.quantile(v,[.16,.84]).tolist())
        for (p,k,b),v in sorted(charge_errors.items())]
    losses=[e for e in events if e['particle']=='gamma' and e['outcomes']['none']['outcome']=='recovered'
            and e['outcomes']['charge']['outcome']=='wrong_type']
    example=sorted(losses,key=lambda e:(e['stem'],e['truth_index']))[0] if losses else None
    result=dict(transitions=transitions,raw_classification=raw,decoded_charge_groups=charge_groups,
                failure_example=example,selection='Post-hoc first gamma event recovered by none but misclassified by charge, sorted by stem/truth index',
                interpretation='Descriptive development diagnostics; no probability calibration or independent-test claim')
    (output/'prior-diagnostics.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    if example is not None:plot_example(project,output,example)
    return result


def plot_example(project,output,event):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    stem=event['stem'];r=json.loads((output/'scores/none'/f'{stem}.json').read_text())
    wave=dict(np.load(project/r['row']['waveform']))
    with np.load(project/r['row']['truth']) as truth:onset=float(truth['first_active_ns'][event['truth_index']])
    t=wave['time_ns']-onset;keep=(t>=-100)&(t<=700)
    fig,axes=plt.subplots(1,2,figsize=(12,4.2))
    axes[0].plot(t[keep],wave['analog_mV'][keep],color='black',lw=.9,label='Observed voltage')
    q=[event['true_charge_pe'],event['outcomes']['none']['prediction']['neural_prior_charge_pe']]
    names=['True\ngamma','Decoded\nNN']
    for arm,color in zip(('none','charge','class','both'),('#777777','#df9a25','#2878b5','#a05293')):
        model=np.load(output/'fits'/arm/f'{stem}.npz')['model_mV']
        axes[0].plot(t[keep],model[keep],color=color,lw=1.,label=arm)
        pred=event['outcomes'][arm]['prediction'];q.append(pred['charge_pe']);names.append(arm+'\n'+pred['particle'])
    axes[0].axvline(0,color='gray',ls=':');axes[0].set(xlabel='Time relative to true gamma onset / ns',ylabel='Voltage / mV')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    axes[1].bar(np.arange(len(q)),q,color=['#3a9b62','#d25c57','#777777','#df9a25','#2878b5','#a05293'])
    axes[1].set(xticks=np.arange(len(q)),xticklabels=names,ylabel='Realized / reconstructed charge / p.e.')
    for i,value in enumerate(q):axes[1].text(i,value+8,f'{value:.0f}',ha='center',fontsize=9)
    axes[1].set_ylim(0,max(q)*1.2);axes[1].grid(axis='y',alpha=.2)
    fig.suptitle('Post-hoc diagnostic: a single decoded charge prior is reused across class hypotheses')
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(output/f'charge-induced-class-flip.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.project,a.output)
