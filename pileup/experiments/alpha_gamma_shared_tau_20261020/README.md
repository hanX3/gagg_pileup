# Independent shared-lifetime short-window study — 2026-10-20

Compare fixed95/351 ns, independently shared, calibration-constrained and event-free double lifetimes. All physical flows profile nonnegative fast/slow amplitudes and a constrained baseline, and fit the onset using the same prefix and pilot-frozen covariance. The readout and sample aperture are convolved with both components. A calibrated full-template GLS is the external control.

## Calibration and frozen design

317 existing calibration windows, 129 existing development windows; no previously public tests used as new development or independent scoring.113 usable calibration signals jointly estimate effective taus91.743114/364.822935 ns. Each calibration event has its own amplitudes/baseline/onset. Shared-lifetime uncertainty uses2000 linearized ROOT bootstrap replicates plus0.02 log-SD floor. Development selects prior inflation3 from1/3/10. This is a Gaussian penalty, not a full hierarchical marginal likelihood.

All decisions, calibration, inference source, unchanged G4 source/configuration and executable are frozen before test production. New168 native single-primary windows:24 alpha and144 monoenergetic gamma, three positions,20 us duration,2 ns full analog samples. Source energies are alpha1.5/2.5/5.486/8 MeV and gamma0.1/0.2/0.661657/1.5 MeV. Seeds262021000–262021167, avalanche/readout seeds+100000/+200000. Actual macros and ROOT RunInfo are authoritative. No ADC rail rejection. Native rate fields are unused single-mode settings, not a rate scan.

## Results

Targets have realized photon-descendant avalanche charge Q>=500 pe. This is not deposited-energy resolution. Joint success means charge within5% and correct class; misses remain failures.

| 触发后T/ns | 方法 | α Q±5% | γ Q±5% | α Q＋类型 | γ Q＋类型 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 100 | fixed_tau | 3/17 | 3/12 | 3/17 | 3/12 |
| 100 | shared_tau | 5/17 | 2/12 | 5/17 | 2/12 |
| 100 | constrained_tau | 5/17 | 1/12 | 5/17 | 1/12 |
| 100 | free_tau | 3/17 | 2/12 | 3/17 | 2/12 |
| 100 | template_gls | 11/17 | 8/12 | 11/17 | 8/12 |
| 150 | fixed_tau | 5/17 | 6/12 | 5/17 | 6/12 |
| 150 | shared_tau | 7/17 | 5/12 | 7/17 | 5/12 |
| 150 | constrained_tau | 9/17 | 3/12 | 9/17 | 3/12 |
| 150 | free_tau | 5/17 | 1/12 | 5/17 | 1/12 |
| 150 | template_gls | 15/17 | 8/12 | 15/17 | 8/12 |
| 200 | fixed_tau | 8/17 | 7/12 | 8/17 | 7/12 |
| 200 | shared_tau | 9/17 | 4/12 | 9/17 | 4/12 |
| 200 | constrained_tau | 13/17 | 4/12 | 13/17 | 4/12 |
| 200 | free_tau | 1/17 | 0/12 | 1/17 | 0/12 |
| 200 | template_gls | 15/17 | 10/12 | 15/17 | 10/12 |
| 300 | fixed_tau | 16/17 | 9/12 | 16/17 | 9/12 |
| 300 | shared_tau | 15/17 | 6/12 | 15/17 | 6/12 |
| 300 | constrained_tau | 15/17 | 5/12 | 15/17 | 5/12 |
| 300 | free_tau | 2/17 | 3/12 | 2/17 | 3/12 |
| 300 | template_gls | 16/17 | 11/12 | 16/17 | 11/12 |
| 500 | fixed_tau | 17/17 | 11/12 | 17/17 | 11/12 |
| 500 | shared_tau | 17/17 | 11/12 | 17/17 | 11/12 |
| 500 | constrained_tau | 17/17 | 11/12 | 17/17 | 11/12 |
| 500 | free_tau | 9/17 | 8/12 | 9/17 | 5/12 |
| 500 | template_gls | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | fixed_tau | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | shared_tau | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | constrained_tau | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | free_tau | 16/17 | 10/12 | 16/17 | 9/12 |
| 1000 | template_gls | 17/17 | 12/12 | 17/17 | 12/12 |

The complete summary retains charge bands, bias/SD/RMSE, missed/wrong-class events, optimizer/boundary/interval failures, and zero-signal gamma windows. Exploratory earliest passing windows are descriptive grid values, not population95% certification. Native fitted taus are effective response parameters; agreement with microscopic95/351 ns is not guaranteed by this approximate voltage model. Fixed/shared taus are external information, not eventwise lifetime measurements.

## Controls and validation

24 independent matched-model calibration signals plus48 tests use exact double-exponential/readout means and0.1mV white Gaussian electronics, without native photon or microcell fluctuations. Their Q is a mean-area parameter, not realized avalanche charge. The observation origin is preset; this is not a trigger-latency test. Separate charge5% and both-tau10% results are retained.24 additive mean double pulses,48 prefixes, expose the single-pulse assumption; no second-pulse recovery or native shared-cell performance is claimed.

9 independent numerical tests,168 native ROOT/RunInfo/three-seed/macro/config/ancestry checks,360 independently rescored summary groups,4 future-NaN/time-shift checks and2 saved-source replays passed. Old code, ROOTs and notebook outputs are preserved.

## Run records and reproduction

Run from the repository root with the project Python environment. Existing stage outputs are validated and reused; raw production never overwrites a ROOT. Native tasks are serial internally, up to4 independent processes, starts at least1 second apart. Ordinary local executable/log caches avoid shared-folder initialization I/O; ROOT remains under`g4/data/`. No symlinks are created.

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_gamma_shared_tau_20261020/analyze_frozen.py test --output jupyter/artifacts/alpha_gamma_shared_tau_20261020 --workers 4
bash tools/run_gagg_python.sh jupyter/run_records/alpha_gamma_shared_tau_20261020/frozen_replay.py --root REPLACE_WITH_RECORDED_ROOT --output /tmp/shared-tau-replay.json
```

`source/brems/g4/` and`bin/` preserve actual native source/config and executable bytes;`inference-source/` preserves the actual Python analysis and dependency lock.`macros/`,`logs/`,`manifest.json` preserve all ROOT mappings, seeds and commands.`git-status-before.txt` and`working-changes.patch` capture pre-existing dirty work. Replaying a native macro requires a fresh second-resolution ROOT output path; preserve the recorded macro and create a new one rather than overwriting it. External G4/ROOT datasets/libraries remain dependencies.

Detailed Chinese process and fixed figures12-1–12-4 are local in`jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb`. Artifacts are in`jupyter/artifacts/alpha_gamma_shared_tau_20261020/`; validation/backup/workflow logs in`jupyter/run_records/alpha_gamma_shared_tau_20261020/`. The shared literature library was read-only; project assumptions and limits are recorded in`docs/literature-use.md` L38.
