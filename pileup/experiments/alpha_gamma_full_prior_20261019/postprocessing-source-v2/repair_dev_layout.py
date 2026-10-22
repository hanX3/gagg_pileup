"""Wrap development-figure tick labels; retain the original rendered files."""
from pathlib import Path
import argparse
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
output=p.parse_args().output
old=output/'figure-layout-v1';old.mkdir(exist_ok=False)
for extension in ('png','pdf'):
    shutil.copy2(output/f'full-prior-outcomes.{extension}',old/f'full-prior-outcomes.{extension}')
summary=json.loads((output/'summary.json').read_text());events=json.loads((output/'events.json').read_text())
methods=summary['methods']
labels=['Physical\nfit','Hybrid\nno\npenalty','Hybrid\nQ +\nclass','Hybrid\nclass\nonly']
fig,axes=plt.subplots(1,2,figsize=(11.5,4.8))
for ax,particle in zip(axes,('alpha','gamma')):
    subset=[e for e in events if e['particle']==particle];bottom=np.zeros(4)
    for outcome,color in zip(['recovered','bad_charge','wrong_type','missed_at_report_threshold'],
                             ['#3a9b62','#e6b441','#d25c57','#999999']):
        counts=np.array([sum(e['outcomes'][m]['outcome']==outcome for e in subset) for m in methods])
        ax.bar(np.arange(4),counts,bottom=bottom,color=color,label=outcome.replace('_',' '));bottom+=counts
    ax.set(xticks=np.arange(4),xticklabels=labels,title=particle,ylabel='Target events')
axes[1].legend(loc='upper left',bbox_to_anchor=(1,1))
fig.tight_layout()
for extension in ('png','pdf'):fig.savefig(output/f'full-prior-outcomes.{extension}',dpi=175,bbox_inches='tight')
plt.close(fig)
print('development_plot_labels_wrapped_original_retained',flush=True)
