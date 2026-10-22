"""Fit isolated GAGG/S13360 waveforms; ROOT truth is kept outside voltage fits.

Run from the project root:
    python -m analysis.alpha_gamma_waveform_fit --output jupyter/artifacts/<new-run>

Voltage objectives are unweighted least squares, not calibrated chi-square:
avalanche fluctuations create correlated residuals. Two saved events illustrate
model behavior; repeated seeds do not constitute a PSD efficiency sample.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import subprocess

import numpy as np
from scipy.optimize import least_squares
from scipy.special import xlogy

from .single_particle_sipm import PROJECT, array, digest, load_root, write_json

DEFAULT_INPUTS = {
    "alpha": "g4/data/20261013_09h58m53s.root",
    "gamma": "g4/data/20261013_09h59m09s.root",
}
FROZEN_CONSTANTS = "jupyter/run_records/sipm_native_20261013_095833/source/g4/include/Constants.hh"


def exponential_cdf(t, tau):
    u = np.maximum(np.asarray(t, dtype=float), 0.)
    return -np.expm1(-u/tau)


def chain_cdf(t, tau, readout):
    """CDF of E(tau)*E(rise)*E(decay)*E(lowpass), each with unit area.

    The present saved configuration has lowpass == rise. Its repeated pole is
    evaluated analytically, avoiding near-coincident numerical partial fractions.
    Distinct positive poles and zero lowpass are also supported. Scintillation
    poles must not coincide with a readout pole in this parameterization.
    """
    r, d, low = (float(readout[k]) for k in ("rise_ns", "decay_ns", "lowpass_ns"))
    if min(tau, r, d) <= 0 or low < 0 or r == d:
        raise ValueError("Require positive distinct rise/decay and nonnegative lowpass")
    if min(abs(tau-r), abs(tau-d)) < 1e-5:
        raise ValueError("Scintillation and readout poles are too close")
    u = np.maximum(np.asarray(t, dtype=float), 0.)
    if np.isclose(low, r, atol=1e-12, rtol=0):
        a = tau**3 / ((tau-d)*(tau-r)**2)
        b = d**3 / ((d-tau)*(d-r)**2)
        e = r**2 / ((r-tau)*(r-d))
        c = 1.-a-b-e
        x = u/r
        result = (a*exponential_cdf(u, tau) + b*exponential_cdf(u, d)
                  + c*exponential_cdf(u, r) + e*(-np.expm1(-x)-x*np.exp(-x)))
    else:
        poles = [tau, r, d] + ([low] if low else [])
        if any(abs(x-y) < 1e-5 for i,x in enumerate(poles) for y in poles[i+1:]):
            raise ValueError("Unsupported repeated pole")
        result = np.zeros_like(u)
        for i, pole in enumerate(poles):
            coefficient = pole**(len(poles)-1) / np.prod([pole-other for j,other in enumerate(poles) if j!=i])
            result += coefficient*exponential_cdf(u, pole)
    return np.clip(result, 0., 1.)


def sampled_component(t, t0, tau, readout):
    """Average a unit-area continuous response over the recorded ADC aperture."""
    width = float(readout["aperture_ns"])
    if width <= 0:
        raise ValueError("A positive sample aperture is required")
    u = np.asarray(t)-t0
    return np.maximum((chain_cdf(u+width/2, tau, readout)
                       - chain_cdf(u-width/2, tau, readout))/width, 0.)


def voltage_model(t, params, readout):
    fast = sampled_component(t, params["t0_ns"], params["tau_fast_ns"], readout)
    slow = sampled_component(t, params["t0_ns"], params["tau_slow_ns"], readout)
    w = params["slow_fraction"]
    return params["baseline_mV"] + params["area_mV_ns"]*((1-w)*fast+w*slow)


def observed_reference(t, y):
    """Locate the isolated pulse and baseline from voltage samples only."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    if len(t) != len(y) or len(t) < 20 or not np.all(np.diff(t) > 0):
        raise ValueError("Invalid waveform grid")
    if not np.isfinite(y).all():
        raise ValueError("Nonfinite waveform")
    peak = int(np.argmax(y))
    pre = (t >= t[peak]-2000.) & (t < t[peak]-500.)
    if pre.sum() < 30:
        raise ValueError("Need at least 30 clean pre-pulse samples")
    baseline = float(np.median(y[pre]))
    sigma = float(1.4826*np.median(np.abs(y[pre]-baseline)))
    amplitude = float(y[peak]-baseline)
    if amplitude < max(20*sigma, 1.):
        raise ValueError("No strong isolated pulse for this shape-fit study")
    threshold = baseline + max(.05*amplitude, 8*sigma)
    search = np.flatnonzero((t >= t[peak]-250.) & (t <= t[peak]) & (y > threshold))
    candidates = [i for i in search if i+2<len(y) and np.all(y[i:i+3]>threshold)]
    if not candidates:
        raise ValueError("No sustained leading-edge threshold crossing")
    return {"anchor_ns": float(t[candidates[0]]), "baseline_mV": baseline,
            "baseline_mad_sigma_mV": sigma, "peak_mV_above_baseline": amplitude,
            "baseline_sample_count": int(pre.sum()), "peak_time_ns": float(t[peak])}


def fit_voltage(t, y, readout, *, taus=(95.,351.), slow_fraction=None,
                single=False, stop_ns=2000., starts=None):
    """Fit the observed voltage, without accepting truth labels or truth times.

    taus=None fits both scintillation lifetimes with positive ordered separation;
    single=True fits one exponential followed by the same recorded readout.
    The baseline is bounded using a disjoint pre-pulse interval. All hypotheses
    use the same data-derived window and baseline constraint.
    """
    ref = observed_reference(t, y)
    time = np.asarray(t, float)-ref["anchor_ns"]
    values = np.asarray(y, float)
    mask = (time >= -200.) & (time <= stop_ns)
    x, observed = time[mask], values[mask]
    baseline = ref["baseline_mV"]
    area0 = max(float(np.trapezoid(np.maximum(observed-baseline,0), x)), 1.)
    bspan = max(.5, 5*ref["baseline_mad_sigma_mV"])
    lower = [math.log(area0/5), baseline-bspan, -30.]
    upper = [math.log(area0*5), baseline+bspan, 15.]
    names = ["log_area", "baseline_mV", "t0_relative_ns"]
    if single:
        lower += [math.log(30.)]; upper += [math.log(1500.)]; names += ["log_tau"]
    else:
        if slow_fraction is None:
            lower += [0.]; upper += [1.]; names += ["slow_fraction"]
        if taus is None:
            lower += [math.log(30.),math.log(30.)]
            upper += [math.log(220.),math.log(2000.)]
            names += ["log_tau_fast","log_tau_separation"]

    def unpack(z):
        p = {"area_mV_ns": math.exp(z[0]), "baseline_mV": float(z[1]), "t0_ns": float(z[2])}
        if single:
            p.update(tau_fast_ns=math.exp(z[3]), tau_slow_ns=math.exp(z[3]), slow_fraction=0.)
        else:
            index = 3
            w = slow_fraction
            if w is None:
                w = float(z[index]); index += 1
            tf, ts = taus if taus is not None else (math.exp(z[index]), math.exp(z[index])+math.exp(z[index+1]))
            p.update(tau_fast_ns=float(tf), tau_slow_ns=float(ts), slow_fraction=float(w))
        return p

    def residual(z):
        return (voltage_model(x, unpack(z), readout)-observed)/ref["peak_mV_above_baseline"]

    solutions = []
    guesses = starts or ([(.5,95.,351.)] if slow_fraction is not None else
                        [(w,tf,ts) for w in (.2,.7) for tf,ts in ((80.,300.),(140.,600.))])
    if single:
        guesses = [(0.,tf,tf) for tf in (80.,200.,400.)]
    for w, tf, ts in guesses:
        initial = [math.log(area0),baseline,-3.]
        if single:
            initial += [math.log(tf)]
        else:
            if slow_fraction is None: initial += [w]
            if taus is None: initial += [math.log(tf),math.log(ts-tf)]
        fit = least_squares(residual, initial, bounds=(lower,upper), x_scale="jac",
                            ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=1500)
        solutions.append(fit)
    best = min(solutions,key=lambda f: float(f.fun@f.fun))
    params = unpack(best.x)
    predicted = voltage_model(x,params,readout)
    residual_mV = observed-predicted
    roi = x >= 0
    rms = float(np.sqrt(np.mean(residual_mV[roi]**2)))
    sse = float(residual_mV@residual_mV)
    raw = residual_mV[roi]-np.mean(residual_mV[roi])
    acf1 = float(raw[:-1]@raw[1:] / (raw@raw)) if raw@raw else 0.
    relative = params.pop("t0_ns")
    params["t0_ns"] = relative+ref["anchor_ns"]
    return {"params":params,"reference":ref,"stop_ns":stop_ns,
            "baseline_bounds_mV":[baseline-bspan,baseline+bspan],
            "sse_mV2":sse,"rms_mV":rms,"nrmse_percent_peak":100*rms/ref["peak_mV_above_baseline"],
            "residual_lag1_correlation":acf1,"n_samples":len(x),"n_parameters":len(best.x),
            "success":bool(best.success),"optimizer_message":str(best.message),
            "active_bounds":[name for name,active in zip(names,best.active_mask) if active],
            "multistart_sse_mV2":[float(f.fun@f.fun)*ref["peak_mV_above_baseline"]**2 for f in solutions],
            "objective":"Unweighted voltage SSE; correlated avalanche residuals, no chi-square p value"}


def fit_photon_tail(arrival_ns, *, taus=(95.,351.), start_ns=20., stop_ns=2000., bin_ns=10.):
    """Bin-integrated Poisson fit to one optical event, with known optical onset.

    This diagnostic uses the first active-area photon. It is not a blind ADC
    estimator, and the start is not passed to any voltage fit.
    """
    origin = float(np.min(arrival_ns))
    edges = np.arange(start_ns,stop_ns+bin_ns/2,bin_ns)
    observed,_ = np.histogram(np.asarray(arrival_ns)-origin,edges)
    n = len(arrival_ns)
    lower=[math.log(n/5),0.]; upper=[math.log(n*5),1.]
    if taus is None:
        lower += [math.log(30.),math.log(30.)]; upper += [math.log(220.),math.log(2000.)]

    def unpack(z):
        tf,ts=taus if taus is not None else (math.exp(z[2]),math.exp(z[2])+math.exp(z[3]))
        return math.exp(z[0]),float(z[1]),float(tf),float(ts)

    def expected(z):
        total,w,tf,ts=unpack(z)
        probabilities=((1-w)*(np.exp(-edges[:-1]/tf)-np.exp(-edges[1:]/tf))
                       +w*(np.exp(-edges[:-1]/ts)-np.exp(-edges[1:]/ts)))
        return np.maximum(total*probabilities,1e-12)

    def residual(z):
        mu=expected(z)
        terms=np.maximum(2*(mu-observed+xlogy(observed,observed/mu)),0.)
        return np.sign(mu-observed)*np.sqrt(terms)

    solutions=[]
    for w,tf,ts in [(w,tf,ts) for w in (.2,.7) for tf,ts in ((80.,300.),(140.,600.))]:
        z=[math.log(n),w]+([math.log(tf),math.log(ts-tf)] if taus is None else [])
        solutions.append(least_squares(residual,z,bounds=(lower,upper),x_scale="jac",
                                      ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=1500))
    best=min(solutions,key=lambda f:float(f.fun@f.fun))
    total,w,tf,ts=unpack(best.x)
    return {"total_extrapolated_photons":total,"slow_fraction":w,"tau_fast_ns":tf,
            "tau_slow_ns":ts,"origin_first_active_photon_ns":origin,"start_ns":start_ns,
            "stop_ns":stop_ns,"bin_ns":bin_ns,"fitted_photons":int(observed.sum()),
            "poisson_deviance":float(best.fun@best.fun),"n_bins":len(observed),
            "n_parameters":len(best.x),"success":bool(best.success),
            "active_bounds":np.flatnonzero(best.active_mask).tolist(),
            "model_note":"One event; optical transport retained; onset truth-assisted; no PSD efficiency"}


def optical_prediction(edges, fit):
    tf,ts,w=fit["tau_fast_ns"],fit["tau_slow_ns"],fit["slow_fraction"]
    return fit["total_extrapolated_photons"]*((1-w)*(np.exp(-edges[:-1]/tf)-np.exp(-edges[1:]/tf))
                                            +w*(np.exp(-edges[:-1]/ts)-np.exp(-edges[1:]/ts)))


def load_event(path):
    data=load_root(path)
    info, opt, primary, digi=(data[k] for k in ("RunInfo","WaveformEvent","PrimaryPhoton","SiPMEvent"))
    if int(opt["n_primary"]) != 1:
        raise ValueError("This analysis requires one known isolated primary")
    cfg=json.loads(str(info["sipm_config_json"]))
    if cfg["sensor"]["name"]!="Hamamatsu_S13360-6050PE":
        raise ValueError("Only the current named SiPM is in scope")
    if int(digi["n_clipped_samples"]):
        raise ValueError("Clipped samples need a censored-data model; this study uses unclipped runs")
    x,y=(array(opt,f"all_photon_arrival_{axis}_mm") for axis in ("x","y"))
    s=cfg["sensor"]
    active=(x>=s["x_min_mm"])&(x<=s["x_max_mm"])&(y>=s["y_min_mm"])&(y<=s["y_max_mm"])
    arrival=array(opt,"all_photon_arrival_t_ps")[active].astype(float)/1000.
    r=cfg["readout"]
    lsb=(r["adc_max_mV"]-r["adc_min_mV"])/2**r["adc_bits"]
    out={"time_ns":array(digi,"sample_time_ns"),"analog_mV":array(digi,"analog_mV"),
         "adc_mV":r["adc_min_mV"]+(array(digi,"adc_counts")+.5)*lsb,
         "signal_mV":array(digi,"signal_mV"),"arrival_ns":arrival,"config":cfg,
         "metadata":{"root":str(Path(path).resolve().relative_to(PROJECT)),"root_sha256":digest(path),
                     "particle_truth":str(primary["primary_particle_name"]),
                     "source_energy_MeV":float(primary["primary_energy_MeV"]),
                     "deposited_energy_MeV":float(opt["edep_total_MeV"]),
                     "source_time_ns_truth":float(primary["primary_time_ps"])/1000.,
                     "n_active_photons":len(arrival),"n_clipped_samples":int(digi["n_clipped_samples"]),
                     "g4_seed":int(info["random_seed"]),"sipm_base_seed":int(info["sipm_avalanche_base_seed"]),
                     "readout_base_seed":int(info["sipm_readout_base_seed"]),
                     "macro_file":str(info["macro_file"]),"macro_contents":str(info["macro_contents"]),
                     "sipm_model_version":str(info["sipm_model_version"]),"adc_lsb_mV":lsb}}
    return out


def analyze_events(events):
    """Run fixed/free shape fits and explicit simulation-informed hypotheses."""
    results={}
    for label,event in events.items():
        r=event["config"]["readout"]
        result={"metadata":event["metadata"],"recorded_config":event["config"],
                "optical_fixed":fit_photon_tail(event["arrival_ns"]),
                "optical_free":fit_photon_tail(event["arrival_ns"],taus=None),"voltage":{}}
        for stage in ("analog","adc"):
            y=event[stage+"_mV"]
            fits={"two_fixed":fit_voltage(event["time_ns"],y,r),
                  "two_free":fit_voltage(event["time_ns"],y,r,taus=None),
                  "single":fit_voltage(event["time_ns"],y,r,single=True),
                  "hyp_alpha":fit_voltage(event["time_ns"],y,r,slow_fraction=.70),
                  "hyp_gamma":fit_voltage(event["time_ns"],y,r,slow_fraction=.21)}
            fits["shape_score_delta_sse_mV2"]=fits["hyp_gamma"]["sse_mV2"]-fits["hyp_alpha"]["sse_mV2"]
            fits["preferred_simulation_hypothesis"]="alpha" if fits["shape_score_delta_sse_mV2"]>0 else "gamma"
            result["voltage"][stage]=fits
        result["adc_window_scan"]=[fit_voltage(event["time_ns"],event["adc_mV"],r,stop_ns=stop)
                                   for stop in (1000.,1500.,2500.)]
        result["adc_free_window_scan"]=[fit_voltage(event["time_ns"],event["adc_mV"],r,taus=None,stop_ns=stop)
                                        for stop in (1000.,1500.,2500.)]
        result["optical_start_scan"]=[fit_photon_tail(event["arrival_ns"],start_ns=start)
                                      for start in (10.,50.,100.)]
        results[label]=result
    return results


def make_figures(events,results,output):
    import os
    os.environ.setdefault("MPLCONFIGDIR",str(PROJECT/"jupyter/artifacts/.matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors={"alpha":"#c64a39","gamma":"#2466a8"}
    figures=[]
    with plt.rc_context({"font.size":10,"axes.spines.top":False,"axes.spines.right":False}):
        fig,axs=plt.subplots(2,2,figsize=(12,7),gridspec_kw={"height_ratios":[3,1]},sharex="col")
        for j,(label,e) in enumerate(events.items()):
            f=results[label]["optical_fixed"]
            edges=np.arange(f["start_ns"],f["stop_ns"]+f["bin_ns"]/2,f["bin_ns"])
            n,_=np.histogram(e["arrival_ns"]-f["origin_first_active_photon_ns"],edges)
            centres=(edges[:-1]+edges[1:])/2
            axs[0,j].stairs(n,edges,color=colors[label],alpha=.7,label="Active-area photons")
            for name,style in (("optical_fixed","-"),("optical_free","--")):
                mu=optical_prediction(edges,results[label][name])
                axs[0,j].plot(centres,mu,style,color="#222222" if name.endswith("fixed") else "#d3922a",label=name.replace("optical_","")+" lifetimes")
            mu=optical_prediction(edges,f)
            axs[1,j].plot(centres,(n-mu)/np.sqrt(mu),".",ms=3,color=colors[label])
            axs[1,j].axhline(0,color="0.4",lw=.8)
            axs[0,j].set(title=f"{label}: before SiPM",ylabel="Photons / 10 ns",yscale="log",ylim=(.2,max(n)*1.4))
            axs[1,j].set(xlabel="Time from first active photon (ns)",ylabel="Pearson residual",xlim=(0,2000))
            axs[0,j].legend(fontsize=9)
        fig.suptitle("One optical event per particle; photon onset used only in this diagnostic")
        fig.tight_layout();p=output/"optical_fits.png";fig.savefig(p,dpi=170);plt.close(fig);figures.append(p.name)

        fig,axs=plt.subplots(2,2,figsize=(12,7),gridspec_kw={"height_ratios":[3,1]},sharex="col")
        style={"two_fixed":("#202020","-","Two decays, fixed lifetimes"),
               "two_free":("#d3922a","--","Two decays, free lifetimes"),
               "single":("#378767",":","Single decay")}
        for j,(label,e) in enumerate(events.items()):
            fits=results[label]["voltage"]["adc"];ref=fits["two_fixed"]["reference"]
            t=e["time_ns"];x=t-ref["anchor_ns"];mask=(x>=-100)&(x<=2000)
            axs[0,j].plot(x[mask],e["adc_mV"][mask]-ref["baseline_mV"],color=colors[label],lw=.8,alpha=.7,label="ADC voltage")
            for name,(color,ls,title) in style.items():
                model=voltage_model(t[mask],fits[name]["params"],e["config"]["readout"])
                axs[0,j].plot(x[mask],model-ref["baseline_mV"],ls,color=color,lw=1.4,label=title)
                axs[1,j].plot(x[mask],e["adc_mV"][mask]-model,ls,color=color,lw=.8,alpha=.85)
            axs[1,j].axhline(0,color="0.4",lw=.8)
            axs[0,j].set(title=f"{label}: after SiPM and ADC",ylabel="Voltage above baseline (mV)")
            axs[1,j].set(xlabel="Time from observed threshold crossing (ns)",ylabel="Residual (mV)",xlim=(-100,2000))
            axs[0,j].legend(fontsize=8)
        fig.suptitle("Full recorded readout response and 0.5 ns aperture; voltage-only timing")
        fig.tight_layout();p=output/"adc_fits_residuals.png";fig.savefig(p,dpi=170);plt.close(fig);figures.append(p.name)

        fig,axs=plt.subplots(1,2,figsize=(12,4.4))
        for j,(label,e) in enumerate(events.items()):
            fits=results[label]["voltage"]["adc"];ref=fits["two_fixed"]["reference"]
            t=e["time_ns"];x=t-ref["anchor_ns"];mask=(x>=-20)&(x<=1000)
            scale=ref["peak_mV_above_baseline"]
            axs[j].plot(x[mask],(e["adc_mV"][mask]-ref["baseline_mV"])/scale,color="0.65",lw=1,label="Observed ADC")
            for particle in ("alpha","gamma"):
                f=fits["hyp_"+particle]
                model=voltage_model(t[mask],f["params"],e["config"]["readout"])
                axs[j].plot(x[mask],(model-ref["baseline_mV"])/scale,color=colors[particle],label=particle+" hypothesis")
            axs[j].set(title=f"Observed {label}: preferred {fits['preferred_simulation_hypothesis']}",
                       xlabel="Time from observed threshold crossing (ns)",ylabel="Voltage / observed peak")
            axs[j].legend()
        fig.suptitle("Simulation-informed fixed shapes; two examples, not a PSD efficiency measurement")
        fig.tight_layout();p=output/"template_hypotheses.png";fig.savefig(p,dpi=170);plt.close(fig);figures.append(p.name)

        fig,axs=plt.subplots(1,2,figsize=(12,4.4))
        for label in events:
            r=results[label]
            windows=sorted(r["adc_window_scan"]+[r["voltage"]["adc"]["two_fixed"]],key=lambda f:f["stop_ns"])
            axs[0].plot([f["stop_ns"] for f in windows],[f["params"]["slow_fraction"] for f in windows],"o-",color=colors[label],label=label+" ADC fit")
            pf=r["optical_fixed"]["slow_fraction"]
            axs[0].axhline(pf,color=colors[label],ls="--",alpha=.7,label=label+" photon fit")
            nrmse=[r["voltage"]["adc"][name]["nrmse_percent_peak"] for name in ("single","two_fixed","two_free")]
            axs[1].plot([0,1,2],nrmse,"o-",color=colors[label],label=label)
        axs[0].set(xlabel="ADC fit end from observed crossing (ns)",ylabel="Fitted slow area fraction",ylim=(0,1))
        axs[1].set(xticks=[0,1,2],xticklabels=["Single decay","Two: fixed lifetimes","Two: free lifetimes"],ylabel="Residual RMS / peak (%)")
        for ax in axs:ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.suptitle("Window sensitivity and descriptive fit error; no confidence interval inferred")
        fig.tight_layout();p=output/"fit_sensitivity.png";fig.savefig(p,dpi=170);plt.close(fig);figures.append(p.name)
    return figures


def run_analysis(output,inputs=None):
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=True)
    inputs=inputs or DEFAULT_INPUTS
    events={label:load_event(PROJECT/path) for label,path in inputs.items()}
    if len({json.dumps(e["config"],sort_keys=True) for e in events.values()})!=1:
        raise ValueError("Compare only identical recorded response configurations")
    results=analyze_events(events)
    figures=make_figures(events,results,output)
    source_hashes={str(Path(__file__).resolve().relative_to(PROJECT)):digest(__file__),
                   FROZEN_CONSTANTS:digest(PROJECT/FROZEN_CONSTANTS)}
    report={"schema_version":1,"purpose":"Two independent optical events; isolated waveform fitting, not efficiency",
            "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=PROJECT,text=True).strip(),
            "source_sha256":source_hashes,"environment":{"python":platform.python_version(),
                **{p:importlib.metadata.version(p) for p in ("numpy","scipy","matplotlib","uproot","awkward")}},
            "simulation_shape_reference":{"source":FROZEN_CONSTANTS,"tau_fast_ns":95.,"tau_slow_ns":351.,
                                          "slow_fractions":{"alpha":.70,"gamma":.21}},
            "voltage_fit_truth_inputs":[],"known_condition":"One isolated pulse per input, selected signal gamma",
            "results":results,"figures":figures}
    write_json(output/"summary.json",report)
    arrays={}
    for label,e in events.items():
        arrays[label+"_time_ns"]=e["time_ns"]
        arrays[label+"_adc_mV"]=e["adc_mV"]
        for name in ("single","two_fixed","two_free","hyp_alpha","hyp_gamma"):
            arrays[f"{label}_{name}_prediction_mV"]=voltage_model(e["time_ns"],results[label]["voltage"]["adc"][name]["params"],e["config"]["readout"])
    np.savez_compressed(output/"fit_arrays.npz",**arrays)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    report=run_analysis(args.output)
    for label,r in report["results"].items():
        print(label,"optical fixed/free:",r["optical_fixed"]["slow_fraction"],r["optical_free"])
        for name in ("single","two_fixed","two_free","hyp_alpha","hyp_gamma"):
            f=r["voltage"]["adc"][name]
            print(name,f["params"],"NRMSE%",f["nrmse_percent_peak"],"bounds",f["active_bounds"])
    print("Saved",args.output)


if __name__=="__main__":
    main()
