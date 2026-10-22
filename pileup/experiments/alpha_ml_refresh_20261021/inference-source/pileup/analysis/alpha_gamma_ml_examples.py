"""Predeclared first-window voltage/event illustrations; post-fit labels only."""
from pathlib import Path
import argparse,json
import numpy as np

def run(folder,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder,output=Path(folder),Path(output)
    project=next(p for p in [Path.cwd(),*Path.cwd().parents] if (p/'g4/data').is_dir())
    frozen=json.loads((folder/'frozen-fit-settings.json').read_text())
    manifest=json.loads((folder/'manifest.json').read_text())
    cases=[]
    for spectrum in ('brems','mono662'):
        for rate in (1000000,10000000):
            row=next(r for r in manifest['runs'] if r['role']=='test_mixed' and
                     r['gamma_model']==spectrum and r['rate_hz']==rate)
            stem=Path(row['root']).stem;wave=dict(np.load(output/'inputs/waveforms'/f'{stem}.npz'))
            truth=dict(np.load(output/'inputs/truth'/f'{stem}.npz'))
            # Independent of every method's success or failure: earliest gamma
            # above the common threshold in the central truth gate.
            central=(truth['emission_ns']>=5000)&(truth['emission_ns']<15000)
            gamma=central&(truth['particle_name']=='gamma')&(truth['charge_pe']>=500)
            onset=float(np.min(truth['first_active_ns'][gamma])) if gamma.any() else 10000.
            lo=max(5000.,min(onset-500.,13000.));hi=lo+2000.
            fig,axes=plt.subplots(2,1,figsize=(12,6),sharex=True,gridspec_kw={'height_ratios':[2,1.5]})
            mask=(wave['time_ns']>=lo)&(wave['time_ns']<=hi)
            axes[0].plot(wave['time_ns'][mask]/1000,wave['analog_mV'][mask],color='#333333',lw=.8,label='Observed native voltage')
            for method,color in [('baseline','#8888cc'),(frozen['hybrid_method'],'#299b79')]:
                model=np.load(output/'results/test'/method/f'{stem}.npz')['model_mV']
                axes[0].plot(wave['time_ns'][mask]/1000,model[mask],color=color,lw=1.,label=method)
            axes[0].set_ylabel('Voltage / mV');axes[0].legend(fontsize=8)
            series=[('Truth',dict(times_ns=truth['first_active_ns'],particle_name=truth['particle_name'],charge_pe=truth['charge_pe']))]
            for method in frozen['common_methods']:
                fit=json.loads((output/'results/test'/method/f'{stem}.json').read_text())['fit'];series.append((method,fit))
            for level,(name,events) in enumerate(series):
                for particle,color,marker in [('alpha','#de8737','^'),('gamma','#2878b5','o')]:
                    times=np.asarray(events['times_ns']);charge=np.asarray(events['charge_pe']);names=np.asarray(events['particle_name'])
                    use=(times>=lo)&(times<=hi)&(charge>=500)&(names==particle)
                    axes[1].scatter(times[use]/1000,np.full(use.sum(),level),s=14+20*np.sqrt(charge[use]/1000),
                        c=color,marker=marker,label=particle if level==0 else None,alpha=.8)
            axes[1].set_yticks(range(len(series)),[x[0] for x in series]);axes[1].invert_yaxis()
            axes[1].set_xlabel('Time / us');axes[1].legend(loc='upper right',fontsize=8)
            axes[1].grid(axis='x',alpha=.2)
            fig.suptitle(f'{spectrum}; alpha source {rate/1e6:g} MHz; gamma source 10 MHz\nFirst declared test window; marker area grows with Q, Q >= 500 p.e.',fontsize=11)
            fig.tight_layout();name=f'example-{spectrum}-{rate}'
            for ext in ('png','pdf'):fig.savefig(output/f'{name}.{ext}',dpi=170,bbox_inches='tight')
            plt.close(fig)
            cases.append(dict(root=row['root'],gamma_model=spectrum,alpha_rate_hz=rate,range_ns=[lo,hi],figure=name))
    (output/'example-selection.json').write_text(json.dumps(dict(cases=cases,
        rule='First declared test per condition; first central gamma Q>=500; independent of fit performance',
        labels_used_after_fitting=True),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.folder,a.output)
