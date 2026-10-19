"""Voltage-only history conditioning of independently calibrated pulse shapes."""
from dataclasses import asdict
import numpy as np
from .alpha_gamma_context_calibration import preceding_voltage
from .alpha_gamma_model_selection_fit import fit_model_selection, ModelSelectionConfig
from .alpha_gamma_shape_profile import fit_shape_profile, ShapeConfig


def shape_parameters(time_ns, voltage_mV, times, charges, calibration):
    history=np.array([preceding_voltage(time_ns,voltage_mV,t,calibration['baseline_mV']) for t in times])
    centers=np.empty((len(times),2));sigmas=np.empty_like(centers)
    for k,particle in enumerate(('alpha','gamma')):
        p=calibration['context_classes'][particle]
        q=np.clip(np.maximum(charges,1.),*p['charge_range_pe'])
        h=np.clip(history,*p['history_range_mV'])
        centers[:,k]=np.clip(p['slow_intercept']+p['slow_log_slope']*np.log(q/1000.)+p['slow_per_100mV']*h/100.,.02,.98)
        if 'fraction_variance_coeff' in p:
            c=p['fraction_variance_coeff']
            sigmas[:,k]=np.sqrt(np.maximum(c[0]+c[1]*1000./q,.003**2))
        else:sigmas[:,k]=p['rms_fraction_residual']
    return centers,sigmas,history


def fit_context(time_ns,voltage_mV,calibration,config=ShapeConfig(sigma_scale=.5),*,initial=None):
    if initial is None:
        initial=fit_model_selection(time_ns,voltage_mV,calibration,ModelSelectionConfig(addition_cost=1.))
    if initial['status']!='fitted' or not initial['all_components']['times_ns']:
        return initial
    c=initial['all_components']
    centers,sigmas,history=shape_parameters(time_ns,voltage_mV,np.array(c['times_ns']),
        np.array(c['shape_charge_pe']),calibration)
    result=fit_shape_profile(time_ns,voltage_mV,calibration,config,initial=initial,
                            center_grid=centers,sigma_grid=sigmas)
    result['context_calibration']=dict(initial_times_ns=c['times_ns'],preceding_voltage_mV=history.tolist(),
        slow_fraction_centers=centers.tolist(),slow_fraction_rms=sigmas.tolist(),
        class_order=['alpha','gamma'],position_input=False,
        feature='Causal observed voltage mean [t-24,t-4) ns; nonnegative baseline-subtracted mV')
    result['objective']='Fixed voltage-only pilot times/weights; independent charge/history-dependent class shapes and RMS penalty'
    result['fit_truth_inputs']=[]
    return result
