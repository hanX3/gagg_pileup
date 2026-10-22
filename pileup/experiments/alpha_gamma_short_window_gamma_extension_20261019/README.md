# Independent causal-prefix isolated-pulse benchmark

2026-10-19. This directory preserves 96 native single-event runs. Combine it with the other `alpha_gamma_short_window*_20261019` directory: the complete test is 30 alpha + 144 gamma windows. Algorithms and calibration were frozen before new generation. These are single-particle controls, not a source-rate or blind multi-pulse experiment.

Each ROOT is one 20 us native Geant4/SiPM window, one primary at 5 us. Native jobs are serial internally, at most four independent jobs, starts at least one second apart. Files remain under project `g4/data/`, named to seconds. All zero-signal gamma windows are retained. Full analog voltage is used with adequate acquisition range; historical ADC clipping does not reject a window.

`manifest.json` records actual macro/ROOT paths, hashes, seeds, positions, source energies, status, and run validation. Source energy is not deposited or visible energy. Actual physical C++ sources/configs, binary, Python sources, macros and logs are retained here; the original physical model is unchanged. `analysis-freeze.json` records frozen source and calibration hashes. External Geant4/ROOT/Python dependencies must still be available.

From project root, the existing native runner resumes missing planned runs without overwriting completed ROOT files:

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_gamma_short_window_gamma_extension_20261019/run_frozen.py test_single --folder pileup/experiments/alpha_gamma_short_window_gamma_extension_20261019 --inputs jupyter/artifacts/alpha_gamma_short_window_implementation_20261019/native-inputs
```

To replay a prefix with saved inference sources, after obtaining the local calibration and exported waveform arrays:

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_gamma_short_window_gamma_extension_20261019/replay_prefix.py --experiment pileup/experiments/alpha_gamma_short_window_gamma_extension_20261019 --root g4/data/20261019_16h19m17s.root --pre 20 --post 100 --output /tmp/gagg-prefix-alpha_gamma_short_window_gamma_extension_20261019.json
```

`replay_prefix.py` checks the saved module and calibration hashes, uses `source/` and avoids dependence on the current analysis code. The destination must be new. The generated `analyze_frozen.py` preparation helper is not the standalone replay entry: process workers can derive a different project root. Use the checked serial replay above for snapshot reproduction. Frozen sources remain untouched; later reporting/validation sources are under `postprocessing-source/`.

Results and limits are in the local [research notebook](../../../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-implementation-20261019), [complete report](../../../jupyter/artifacts/alpha_gamma_short_window_implementation_20261019/report.json), and [ROOT validation](../../../jupyter/run_records/alpha_gamma_short_window_implementation_20261019/native-data-validation.json). The combined list is `test-input-records-combined.json`, not the per-campaign `test-input-records.json`.

Four voltage-only fits assume at most one pulse. Primary results use Q >= 500 pe and charge error <=5%; missed signals remain failures. Q is realized photon-descendant avalanche charge, not energy. Finite test fractions do not certify a universal shortest window. This study fixes response shape information and does not measure the minimum window for free decay constants or a maximum count rate.

Git HEAD at preparation: `6ae005296b26a7b52f20fef18da9e5e3fe8996ec` with uncommitted local analysis changes; full snapshots and hashes preserve those changes. Python/dependency versions are in local `freeze-complete.json`. ROOT, notebooks, artifacts, and run records are local and not included in a push. An optional `/tmp` binary cache attempt failed at symbolic-link creation; original binaries remained in place and all production ran normally.
