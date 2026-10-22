"""Equivalent offline digitizer: bulk independent photon RNG and cached PDE interpolation.

No parameters or proposal ordering change. Adoption requires the recorded exact
array comparison against the frozen scalar implementation.
"""
from bisect import bisect_right
import heapq
import math
import numpy as np
from .sipm_response import (map_cells, trigger_fraction, _neighbours, HC_EV_NM,
                            KINDS, AVALANCHE_DTYPE, AvalancheResult)

def simulate_avalanches(hits, config, seed):
    """Merge photons, primary dark proposals and recursive secondary proposals.

    At full recovery, dark_rate_hz is the primary dark avalanche rate. Each
    avalanche proposes at most one XT and one AP, with probabilities p*q/q0.
    Their target cell accepts with trigger_fraction(recovery). All mechanisms
    share cell state. Table PDE is sampled directly once, without pre-thinning.
    """
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("Seed must be a nonnegative integer")
    if hits.energy_eV is None and config.pde_mode != "constant":
        raise ValueError("Energy-resolved PDE requires photon_energy_eV")
    cells = map_cells(hits, config)
    rng_photon, rng_transit, rng_dark, rng_secondary, rng_trigger = [
        np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(5)]
    times = hits.time_ns.copy()
    if config.transit_delay_sigma_ns:
        shape = (config.transit_delay_mean_ns / config.transit_delay_sigma_ns) ** 2
        scale = config.transit_delay_sigma_ns ** 2 / config.transit_delay_mean_ns
        times += rng_transit.gamma(shape, scale, len(times))
    else:
        times += config.transit_delay_mean_ns
    if not np.isfinite(times).all():
        raise ValueError("Nonfinite transit delays")
    # Original photon indices remain usable after transit delays reorder photons.
    order = np.argsort(times, kind="stable")
    order = order[(cells[order] >= 0) & (times[order] < hits.window_ns)]
    photon_uniform = rng_photon.random(len(order))
    if config.pde_mode == "table":
        wavelengths = HC_EV_NM / hits.energy_eV
        w = np.asarray(config.pde_wavelength_nm)
        volts = np.asarray(config.pde_overvoltage_V)
        table = np.asarray(config.pde_probabilities)
        if len(order) and (np.any(wavelengths[order] < w[0]) or np.any(wavelengths[order] > w[-1])):
            raise ValueError("PDE wavelength table does not cover all photons on the sensor")
        wavelength_index = np.clip(np.searchsorted(w, wavelengths, side="right") - 1, 0, len(w) - 2)
        wavelength_fraction = (wavelengths - w[wavelength_index]) / (w[wavelength_index + 1] - w[wavelength_index])
        pde_by_voltage = table[:, wavelength_index] * (1 - wavelength_fraction)[None, :] + table[:, wavelength_index + 1] * wavelength_fraction[None, :]
    last = np.full(config.n_cells, -np.inf)
    queue, records = [], []
    sequence = 0
    counters = {"input_photons": len(times), "outside_footprint": int(np.sum(cells < 0)),
                "photon_candidates_after_window": int(np.sum((cells >= 0) & (times >= hits.window_ns))),
                "processed_candidates": 0, "dark_proposals": 0,
                "crosstalk_proposals": 0, "afterpulse_proposals": 0,
                "crosstalk_after_window": 0, "afterpulse_after_window": 0}

    def enqueue(t, cell, kind, parent):
        nonlocal sequence
        if not math.isfinite(t):
            raise ValueError("Nonfinite candidate time")
        if t >= hits.window_ns:
            if kind in (2, 3):
                counters[KINDS[kind] + "_after_window"] += 1
            return
        sequence += 1
        heapq.heappush(queue, (t, sequence, cell, kind, parent))

    def next_dark(t):
        if config.dark_rate_hz:
            t += rng_dark.exponential(1e9 / config.dark_rate_hz)
            if t < hits.window_ns:
                enqueue(t, int(rng_dark.integers(config.n_cells)), 1, -1)

    next_dark(0.0)
    pointer = 0
    while pointer < len(order) or queue:
        next_time = times[order[pointer]] if pointer < len(order) else math.inf
        # Queued prompt secondaries precede remaining photons with equal times.
        if queue and queue[0][0] <= next_time:
            t, _, cell, kind, parent = heapq.heappop(queue)
            photon_index = -1
            if kind == 1:
                counters["dark_proposals"] += 1
                next_dark(t)
        else:
            photon_index = int(order[pointer])
            pointer += 1
            t, cell, kind, parent = float(times[photon_index]), int(cells[photon_index]), 0, -1
        counters["processed_candidates"] += 1
        if counters["processed_candidates"] > config.max_candidates:
            raise RuntimeError("Candidate limit reached; increase max_candidates explicitly after checking the model")
        r = (-math.expm1(-(t - last[cell]) / config.recovery_ns)
             if config.recovery_enabled else 1.0)
        if kind == 0:
            if config.pde_mode == "constant":
                probability = config.constant_pde * trigger_fraction(r, config)
            else:
                voltage = r * config.overvoltage_V
                vi = min(bisect_right(volts, voltage) - 1, len(volts) - 2)
                vi = max(vi, 0)
                vf = (voltage - volts[vi]) / (volts[vi + 1] - volts[vi])
                wi, wf = wavelength_index[photon_index], wavelength_fraction[photon_index]
                low = pde_by_voltage[vi, photon_index]
                high = pde_by_voltage[vi + 1, photon_index]
                probability = low * (1 - vf) + high * vf
            accepted = photon_uniform[pointer - 1] < probability
        else:
            accepted = rng_trigger.random() < trigger_fraction(r, config)
        if not accepted:
            continue
        aid = len(records)
        records.append((t, cell, config.full_charge_C * r, kind, parent, photon_index, r))
        last[cell] = t
        if config.crosstalk_probability and rng_secondary.random() < config.crosstalk_probability * r:
            neighbours = _neighbours(cell, config)
            if neighbours:
                target = neighbours[int(rng_secondary.integers(len(neighbours)))]
                delay = rng_secondary.exponential(config.crosstalk_delay_ns) if config.crosstalk_delay_ns else 0.0
                counters["crosstalk_proposals"] += 1
                enqueue(t + delay, target, 2, aid)
        if config.afterpulse_probability and rng_secondary.random() < config.afterpulse_probability * r:
            counters["afterpulse_proposals"] += 1
            enqueue(t + rng_secondary.exponential(config.afterpulse_delay_ns), cell, 3, aid)
    data = np.array(records, dtype=AVALANCHE_DTYPE)
    gate = (data["time_ns"] >= hits.gate_start_ns) & (data["time_ns"] < hits.gate_end_ns)
    counters.update({f"{name}_avalanches": int(np.sum(data["kind"] == k)) for k, name in KINDS.items()})
    return AvalancheResult(data, hits.window_ns, config.full_charge_C, last,
                           {"seed": int(seed), "counters": counters,
                            "initial_state": "fully recovered, no pre-window secondary history",
                            "gate_ns": [hits.gate_start_ns, hits.gate_end_ns],
                            "gate_avalanche_count": int(gate.sum()),
                            "gate_avalanche_charge_C": float(data["charge_C"][gate].sum()),
                            "avalanche_charge_C": float(data["charge_C"].sum()),
                            "kind_codes": KINDS})
