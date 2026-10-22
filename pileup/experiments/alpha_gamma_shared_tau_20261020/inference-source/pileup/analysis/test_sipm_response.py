"""Physical limits, event ordering, units and reproducibility of the SiPM model."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import tempfile
import unittest

import awkward as ak
import numpy as np
import uproot

from pileup.analysis.sipm_response import (
    AVALANCHE_DTYPE, HC_EV_NM, AvalancheResult, PhotonHits, SensorConfig,
    digitize_root, map_cells, read_photon_hits, sha256, simulate_avalanches,
)
from pileup.analysis.sipm_readout import ReadoutConfig, synthesize_readout


def sensor(**changes):
    base = SensorConfig(
        name="test only", parameter_status="hypothetical", parameter_source="synthetic test",
        temperature_C=20, cells_x=10, cells_y=10, x_min_mm=0, x_max_mm=1,
        y_min_mm=0, y_max_mm=1, overvoltage_V=3, gain_electrons=1e6,
        recovery_ns=10, trigger_scale_V=1, constant_pde=1,
    )
    return replace(base, **changes)


def plateau_sensor(**changes):
    # Deterministic acceptance for any appreciably recovered cell; synthetic only.
    return sensor(pde_mode="table", constant_pde=None,
                  pde_wavelength_nm=(400, 500, 600, 650), pde_overvoltage_V=(0, 1e-8, 3),
                  pde_probabilities=((0, 0, 0, 0), (1, 1, 0, 0), (1, 1, 0, 0)), **changes)


def photons(times=(), *, x=None, y=None, wavelengths=None, window=1000, gate=0):
    n = len(times)
    return PhotonHits(times, np.full(n, .25) if x is None else x,
                      np.full(n, .25) if y is None else y,
                      HC_EV_NM / np.full(n, 500.) if wavelengths is None else HC_EV_NM / np.asarray(wavelengths),
                      window, gate, window)


def readout(**changes):
    return replace(ReadoutConfig(model="double_exponential", internal_dt_ns=.25,
                                 sample_period_ns=2, aperture_ns=.5, rise_ns=2, decay_ns=20,
                                 transimpedance_ohm=50, lowpass_ns=0, baseline_mV=0,
                                 noise_rms_mV=0, adc_bits=12, adc_min_mV=-20, adc_max_mV=20), **changes)


def avalanches(times=(), *, fractions=None, window=1000):
    q0 = sensor().full_charge_C
    if fractions is None:
        fractions = np.ones(len(times))
    data = np.array([(t, 0, q0 * q, 0, -1, i, q) for i, (t, q) in enumerate(zip(times, fractions))],
                    dtype=AVALANCHE_DTYPE)
    return AvalancheResult(data, window, q0, np.array([]), {})


class AvalanchePhysicsChecks(unittest.TestCase):
    def test_dark_disabled_zero_light(self):
        result = simulate_avalanches(photons(), sensor(), 1)
        self.assertEqual(len(result.records), 0)
        self.assertEqual(result.metadata["avalanche_charge_C"], 0)

    def test_simultaneous_same_cell_saturates_but_other_cells_do_not(self):
        same = simulate_avalanches(photons([0, 0]), sensor(cells_x=2, cells_y=1), 1)
        different = simulate_avalanches(photons([0, 0], x=[.25, .75]), sensor(cells_x=2, cells_y=1), 1)
        self.assertEqual(len(same.records), 1)
        self.assertEqual(len(different.records), 2)
        np.testing.assert_allclose(different.records["charge_C"], sensor().full_charge_C)

    def test_recovery_charge_and_gate_do_not_reset_state(self):
        h = photons([0, 10, 30], gate=10)
        a = simulate_avalanches(h, plateau_sensor(), 1)
        np.testing.assert_allclose(a.records["recovery_fraction"], [1, 1-np.exp(-1), 1-np.exp(-2)])
        self.assertEqual(a.metadata["gate_avalanche_count"], 2)

    def test_rejected_photon_does_not_restart_recovery(self):
        a = simulate_avalanches(photons([0, 10, 20], wavelengths=[500, 600, 500]), plateau_sensor(), 1)
        np.testing.assert_array_equal(a.records["photon_index"], [0, 2])
        self.assertAlmostEqual(a.records["recovery_fraction"][-1], 1-np.exp(-2))

    def test_linear_baseline_has_binomial_mean_and_variance(self):
        n, p, trials = 80, .23, 1500
        h, s = photons(np.zeros(n)), sensor(constant_pde=p, recovery_enabled=False)
        counts = np.array([len(simulate_avalanches(h, s, i).records) for i in range(trials)])
        variance = n*p*(1-p)
        self.assertLess(abs(counts.mean()-n*p), 6*np.sqrt(variance/trials))
        self.assertLess(abs(counts.var(ddof=1)-variance), .15*variance)

    def test_two_dimensional_pde_interpolation(self):
        n = 12000
        s = sensor(recovery_enabled=False, overvoltage_V=1.5, pde_mode="table", constant_pde=None,
                   pde_wavelength_nm=(400, 600), pde_overvoltage_V=(0, 3),
                   pde_probabilities=((0, 0), (.4, .8)))
        h = photons(np.zeros(n))  # 500 nm, 1.5 V gives p=0.3.
        count = len(simulate_avalanches(h, s, 9).records)
        self.assertLess(abs(count - n*.3), 6*np.sqrt(n*.3*.7))

    def test_short_flash_occupancy_fixed_photon_number(self):
        rng = np.random.default_rng(52)
        n, p, nc, trials = 600, .5, 100, 250
        s = sensor(constant_pde=p)
        counts = []
        for i in range(trials):
            h = photons(np.zeros(n), x=rng.random(n), y=rng.random(n))
            counts.append(len(simulate_avalanches(h, s, i).records))
        expected = nc*(1-(1-p/nc)**n)
        self.assertLess(abs(np.mean(counts)-expected), .5)
        self.assertLessEqual(max(counts), nc)

    def test_short_flash_occupancy_poisson_illumination(self):
        rng = np.random.default_rng(61)
        mean, p, nc = 500, .3, 100
        s = sensor(constant_pde=p)
        counts = []
        for i in range(300):
            n = rng.poisson(mean)
            h = photons(np.zeros(n), x=rng.random(n), y=rng.random(n))
            counts.append(len(simulate_avalanches(h, s, i).records))
        self.assertLess(abs(np.mean(counts)-nc*(1-np.exp(-mean*p/nc))), .7)

    def test_dark_rate_is_primary_rate_at_full_recovery(self):
        s = sensor(recovery_enabled=False, dark_rate_hz=1e6)
        results = [simulate_avalanches(photons(window=10000), s, i).records for i in range(1000)]
        counts = np.array([len(a) for a in results])
        self.assertLess(abs(counts.mean()-10), .6)
        self.assertLess(abs(counts.var(ddof=1)-10), 1.5)
        times = np.concatenate([a["time_ns"] for a in results])
        self.assertLess(abs(times.mean()-5000), 100)

    def test_prompt_crosstalk_respects_neighbours_and_recovery(self):
        s = sensor(cells_x=2, cells_y=1, crosstalk_probability=.8)
        counts = []
        for i in range(200):
            a = simulate_avalanches(photons([0]), s, i).records
            counts.append(len(a))
            self.assertLessEqual(len(a), 2)
            if len(a) == 2:
                self.assertEqual(a[1]["kind"], 2)
                self.assertNotEqual(a[0]["cell_id"], a[1]["cell_id"])
                self.assertEqual(a[1]["parent_id"], 0)
        self.assertGreater(np.mean(counts), 1.65)

    def test_crosstalk_chain_multiplicity_in_linear_limit(self):
        s = sensor(cells_x=2, cells_y=1, recovery_enabled=False, crosstalk_probability=.2)
        counts = [len(simulate_avalanches(photons([0]), s, i).records) for i in range(1500)]
        self.assertLess(abs(np.mean(counts)-1/(1-.2)), .08)

    def test_afterpulse_delay_and_multiplicity_in_linear_limit(self):
        s = sensor(recovery_enabled=False, afterpulse_probability=.7, afterpulse_delay_ns=40)
        delays, counts = [], []
        for i in range(1200):
            a = simulate_avalanches(photons([0], window=10000), s, i).records
            counts.append(len(a))
            ap = a[a["kind"] == 3]
            delays.extend(ap["time_ns"] - a[ap["parent_id"]]["time_ns"])
        self.assertLess(abs(np.mean(delays)-40), 3)
        self.assertLess(abs(np.mean(counts)-1/(1-.7)), .25)

    def test_all_noise_types_share_chronological_microcell_state(self):
        s = sensor(cells_x=2, cells_y=1, dark_rate_hz=1e8,
                   crosstalk_probability=.2, crosstalk_delay_ns=.5,
                   afterpulse_probability=.5, afterpulse_delay_ns=5)
        a = simulate_avalanches(photons(np.arange(0, 500, 5), window=1000), s, 81).records
        self.assertTrue(np.all(np.diff(a["time_ns"]) >= 0))
        self.assertEqual(set(a["kind"]), {0, 1, 2, 3})
        last = {}
        for i, event in enumerate(a):
            cell, t = int(event["cell_id"]), event["time_ns"]
            expected = 1 if cell not in last else 1-np.exp(-(t-last[cell])/s.recovery_ns)
            self.assertAlmostEqual(event["charge_C"]/s.full_charge_C, expected, places=12)
            if event["parent_id"] >= 0:
                self.assertLess(event["parent_id"], i)
            last[cell] = t
        self.assertTrue(np.any(a["recovery_fraction"][a["kind"] == 3] < .9))

    def test_late_secondary_candidates_are_reported(self):
        s = sensor(afterpulse_probability=.8, afterpulse_delay_ns=1000)
        late = 0
        for i in range(100):
            a = simulate_avalanches(photons([99.9], window=100), s, i)
            self.assertTrue(np.all(a.records["time_ns"] < 100))
            late += a.metadata["counters"]["afterpulse_after_window"]
        self.assertGreater(late, 60)

    def test_causal_transit_spread_is_sorted_with_original_indices(self):
        n = 10000
        h = photons(np.arange(n)*.01, window=1000)
        s = sensor(recovery_enabled=False, transit_delay_mean_ns=5, transit_delay_sigma_ns=2)
        a = simulate_avalanches(h, s, 17).records
        delays = a["time_ns"] - h.time_ns[a["photon_index"]]
        self.assertTrue(np.all(delays >= 0))
        self.assertTrue(np.all(np.diff(a["time_ns"]) >= 0))
        self.assertFalse(np.array_equal(a["photon_index"], np.arange(n)))
        self.assertLess(abs(delays.mean()-5), .08)
        self.assertLess(abs(delays.std()-2), .08)

    def test_seed_reproducibility(self):
        h = photons(np.arange(500), window=1000)
        s = sensor(constant_pde=.4, dark_rate_hz=1e7, afterpulse_probability=.1)
        a, b = simulate_avalanches(h, s, 123), simulate_avalanches(h, s, 123)
        np.testing.assert_array_equal(a.records, b.records)
        self.assertEqual(a.metadata, b.metadata)

    def test_footprint_edges_and_outside_photons(self):
        h = photons([0]*6, x=[0, .49, .5, 1, -.001, 1.001])
        s = sensor(cells_x=2, cells_y=1)
        np.testing.assert_array_equal(map_cells(h, s), [0, 0, 1, 1, -1, -1])
        self.assertEqual(simulate_avalanches(h, s, 1).metadata["counters"]["outside_footprint"], 2)

    def test_invalid_input_and_uncovered_spectrum_are_rejected(self):
        for times in ([1, 0], [-1], [1000], [np.nan]):
            with self.assertRaises(ValueError):
                photons(times)
        with self.assertRaises(ValueError):
            PhotonHits([0], [], [0], [2], 1000)
        with self.assertRaises(ValueError):
            simulate_avalanches(photons([0], wavelengths=[700]), plateau_sensor(), 1)
        with self.assertRaises(ValueError):
            sensor(crosstalk_probability=.6, afterpulse_probability=.5)
        with self.assertRaises(RuntimeError):
            simulate_avalanches(photons([0, 1, 2]), sensor(max_candidates=2), 1)


class ReadoutPhysicsChecks(unittest.TestCase):
    def test_off_grid_avalanches_conserve_integrated_charge(self):
        a = avalanches([.137, 10.123], fractions=[1, .5], window=1000)
        w = synthesize_readout(a, readout(), 1)
        recovered = w.arrays["current_A"].sum()*.25e-9
        self.assertAlmostEqual(recovered/a.records["charge_C"].sum(), 1, places=12)
        self.assertAlmostEqual(w.arrays["signal_mV"][0],
                               w.arrays["current_A"][:2].mean()*50*1000, places=12)

    def test_window_tail_matches_analytic_integral(self):
        a = avalanches([19.8], window=20)
        w = synthesize_readout(a, readout(), 1)
        tail_fraction = (20*np.exp(-.2/20)-2*np.exp(-.2/2))/18
        observed = w.metadata["in_window_current_charge_C"]/a.full_charge_C
        self.assertAlmostEqual(observed, 1-tail_fraction, places=12)
        self.assertLess(observed, .001)

    def test_linear_superposition_with_readout_filter(self):
        r = readout(lowpass_ns=3)
        first = synthesize_readout(avalanches([10.12]), r, 1)
        second = synthesize_readout(avalanches([17.65], fractions=[.4]), r, 1)
        together = synthesize_readout(avalanches([10.12, 17.65], fractions=[1, .4]), r, 1)
        np.testing.assert_allclose(together.arrays["signal_mV"],
                                   first.arrays["signal_mV"]+second.arrays["signal_mV"], rtol=1e-12, atol=1e-15)

    def test_readout_time_step_converges(self):
        a = avalanches([10.137], window=200)
        signals = [synthesize_readout(a, readout(internal_dt_ns=dt, aperture_ns=1, lowpass_ns=3), 1)
                   .arrays["signal_mV"] for dt in (.5, .25, .125, .0625)]
        errors = [np.linalg.norm(s-signals[-1]) for s in signals[:-1]]
        self.assertTrue(errors[0] > errors[1] > errors[2])
        self.assertLess(errors[1]/np.linalg.norm(signals[-1]), .005)

    def test_measured_template_scaling_area_and_causality(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"spe.csv"
            path.write_text("time_ns,voltage_mV\n0,0\n1,1\n2,0\n")
            r = readout(model="measured_voltage", template_csv=str(path),
                        sample_period_ns=.25, aperture_ns=.25)
            a = avalanches([10.13], fractions=[.4], window=100)
            w = synthesize_readout(a, r, 1)
            v, t = w.arrays["signal_mV"], w.arrays["time_ns"]
            self.assertAlmostEqual(v.sum()*.25, .4, places=12)
            self.assertLess(np.max(np.abs(v[t < 10.13])), 1e-14)
            self.assertAlmostEqual(np.sum(t*v)/v.sum(), 11.25, places=10)
            self.assertEqual(w.arrays["current_A"].size, 0)
            with self.assertRaises(ValueError):
                replace(r, lowpass_ns=1)

    def test_adc_quantization_and_both_clipping_limits(self):
        a = avalanches(window=100)
        r = readout(adc_bits=2, adc_min_mV=0, adc_max_mV=20, baseline_mV=12)
        w = synthesize_readout(a, r, 1)
        np.testing.assert_array_equal(w.arrays["adc_counts"], 2)
        self.assertEqual(w.metadata["saturated_samples"], 0)
        for baseline, code in ((-1, 0), (20, 3)):
            w = synthesize_readout(a, replace(r, baseline_mV=baseline), 1)
            np.testing.assert_array_equal(w.arrays["adc_counts"], code)
            self.assertEqual(w.metadata["saturated_samples"], 50)

    def test_readout_noise_and_seed_are_separate_from_avalanches(self):
        a = avalanches([10], window=20000)
        before = a.records.copy()
        r = readout(sample_period_ns=1, noise_rms_mV=.2)
        first, again, other = [synthesize_readout(a, r, seed) for seed in (11, 11, 12)]
        residual = first.arrays["analog_mV"]-first.arrays["signal_mV"]
        self.assertLess(abs(residual.mean()), .01)
        self.assertLess(abs(residual.std()-.2), .01)
        np.testing.assert_array_equal(first.arrays["analog_mV"], again.arrays["analog_mV"])
        np.testing.assert_array_equal(first.arrays["signal_mV"], other.arrays["signal_mV"])
        self.assertFalse(np.array_equal(first.arrays["analog_mV"], other.arrays["analog_mV"]))
        np.testing.assert_array_equal(before, a.records)

    def test_invalid_sampling_is_explicit(self):
        with self.assertRaises(ValueError):
            readout(sample_period_ns=.3)
        with self.assertRaises(ValueError):
            synthesize_readout(avalanches(window=101), readout(), 1)


class RootInterfaceChecks(unittest.TestCase):
    def write_root(self, path, energy=True):
        with uproot.recreate(path) as root:
            root["RunInfo"] = {"t_length_ps": np.array([100000]),
                               "analysis_start_ps": np.array([20000]), "analysis_end_ps": np.array([80000])}
            branches = {"event_id": np.array([0]), "n_optical_photons_arrived_total": np.array([3]),
                        "all_photon_arrival_t_ps": ak.Array([[0, 20000, 50000]]),
                        "all_photon_arrival_x_mm": ak.Array([[.1, .2, .3]]),
                        "all_photon_arrival_y_mm": ak.Array([[.2, .2, .2]])}
            if energy:
                branches["photon_energy_eV"] = ak.Array([[2.4, 2.5, 2.6]])
            root["WaveformEvent"] = branches

    def test_energy_alignment_and_explicit_legacy_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"input.root"
            self.write_root(path)
            h = read_photon_hits(path)
            np.testing.assert_array_equal(h.time_ns, [0, 20, 50])
            np.testing.assert_allclose(h.energy_eV, [2.4, 2.5, 2.6])
            self.assertEqual((h.gate_start_ns, h.gate_end_ns), (20, 80))
            legacy = Path(directory)/"legacy.root"
            self.write_root(legacy, energy=False)
            with self.assertRaises(ValueError):
                read_photon_hits(legacy)
            h = read_photon_hits(legacy, allow_missing_energy=True)
            self.assertIsNone(h.energy_eV)
            with self.assertRaises(ValueError):
                simulate_avalanches(h, plateau_sensor(), 1)

    def test_output_provenance_input_preservation_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            path, config, output = d/"input.root", d/"config.json", d/"result"
            self.write_root(path)
            config.write_text(json.dumps({"schema_version": 1, "sensor": asdict(sensor()),
                                         "readout": asdict(readout())}))
            original = sha256(path)
            manifest = digitize_root(path, config, output, 1, 2)
            self.assertEqual(manifest["status"], "completed")
            self.assertEqual(sha256(path), original)
            self.assertEqual(manifest["input_sha256"], original)
            for name, expected in manifest["output_sha256"].items():
                self.assertEqual(sha256(output/name), expected)
            with np.load(output/"avalanches.npz", allow_pickle=False) as data:
                self.assertEqual(len(data["records"]), 3)
            with self.assertRaises(FileExistsError):
                digitize_root(path, config, output, 1, 2)

    def test_failed_run_is_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            path, config, output = d/"legacy.root", d/"config.json", d/"result"
            self.write_root(path, energy=False)
            config.write_text(json.dumps({"schema_version": 1, "sensor": asdict(sensor()),
                                         "readout": asdict(readout())}))
            with self.assertRaises(ValueError):
                digitize_root(path, config, output, 1, 2)
            manifest = json.loads((output/"manifest.json").read_text())
            self.assertEqual(manifest["status"], "failed")
            self.assertIn("photon_energy_eV", manifest["error"])

    def test_legacy_constant_pde_run_is_marked_as_approximation(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            path, config = d/"legacy.root", d/"config.json"
            self.write_root(path, energy=False)
            config.write_text(json.dumps({"schema_version": 1, "sensor": asdict(sensor()),
                                         "readout": asdict(readout())}))
            result = digitize_root(path, config, d/"result", 1, 2, allow_missing_energy=True)
            self.assertTrue(result["input_metadata"]["missing_energy_approximation"])
            self.assertEqual(result["status"], "completed")

    def test_relative_measured_template_is_copied_and_replayable(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            path, config, template = d/"input.root", d/"config.json", d/"spe.csv"
            self.write_root(path)
            template.write_text("time_ns,voltage_mV\n0,0\n1,1\n3,0\n")
            config.write_text(json.dumps({"schema_version": 1, "sensor": asdict(sensor()),
                                         "readout": asdict(readout(model="measured_voltage", template_csv="spe.csv"))}))
            first = digitize_root(path, config, d/"first", 1, 2)
            second = digitize_root(path, d/"first/config.json", d/"second", 1, 2)
            self.assertEqual(first["template_sha256"], sha256(template))
            self.assertEqual(first["template_sha256"], second["template_sha256"])
            with np.load(d/"first/waveform.npz") as a, np.load(d/"second/waveform.npz") as b:
                for key in a.files:
                    np.testing.assert_array_equal(a[key], b[key])


if __name__ == "__main__":
    unittest.main()
