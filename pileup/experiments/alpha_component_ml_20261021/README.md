# Component-supervised physical training and U-Net–GRU comparison

Run from the repository root. The local research record is Notebook 0004,
sections 23–28. Training, frozen native evaluation and numerical checks have completed. This is a GAGG adaptation inspired by SAEFit and time-sequence
unfolding, not a reproduction of their detector models or published efficiency.

The four workflows are GRU/original loss, GRU/component loss,
U-Net–GRU/original loss, and U-Net–GRU/component loss. Each uses two declared
training seeds and 2400 steps on the same 49 training ROOTs and 26 validation
ROOTs. Previously published 23 tests are now training data; the last 16
confirmations are validation data. They are no longer held-out tests.

The additional loss compares native ancestry components with a fixed physical
decoder, plus a conditional sum of those supervised components. True centers
are used only in training. Inference retains the original dense three heads,
2 ns grid and 8 ns peak selection; multiplicity is unknown. The decoder fixes
95/351 ns scintillation decays and includes the calibrated electronics. Its
1216 ns loss support is not a minimum acquisition-window result. Residuals at
true Q/time are reference model mismatch, not a mathematical lower bound.

The frozen 60-window plan has 40/8/2/2/2/2/2/2 windows at source injection
0.1/1/10/20/40/60/80/100 MHz. The low-rate counts were expanded before the
all-model freeze; all 23 initially chosen seeds were retained. The old plan and
revision are recorded locally. Weights and thresholds were selected on
development/validation data and frozen before the first new ROOT. The final
test plan is immutable.

The C++ source and executable are copied from the validated prior experiment.
Each ROOT contains one serial G4 WaveformEvent, a 20 us record and a central
10 us scoring gate. The SiPM retains shared-cell history, dark chains and
electronic noise. Analysis uses complete analog voltage with adequate dynamic
range. Q is realized avalanche-equivalent charge, not deposited energy.

`source/g4`, `bin`, `macros`, `assets`, `inference-source`, and `manifest.json`
retain actual source/configuration/run inputs. Local models, inputs and results
remain under `jupyter/artifacts/alpha_component_ml_20261021`; runtime records
are under `jupyter/run_records/alpha_component_ml_20261021`. Reading the saved
Notebook needs the analysis environment; retraining or inference also needs
the optional PyTorch dependency in `requirements-ml.txt` and the retained local
input files. No Git checkout is needed.

Use `replay.py infer --waveform <saved-waveform.npz> --output <new-directory>`
to replay frozen inference; `replay.py train --variant <name> --seed <seed>
--output <new-directory>` retrains from saved development inputs. Both refuse
an existing output directory. Original production uses `run_frozen.py` stages.

The frozen software-control routine completed numerical work but initially
failed to JSON-encode physical model arrays. `save_controls.py` regenerates the
same 62 cases with the original seed and moves those arrays to NPZ before JSON
serialization. Frozen inference source, model weights and fit parameters are
unchanged; the repair has its own source hash and local record.

`stream_evaluate.py` schedules unchanged frozen inference on each completed
predeclared record while remaining simulations run. `validate_saved.py` retains
the archived validator, exempting only the `refinement_seconds` wall-clock field
from exact numerical replay. It also checks original ROOT time/voltage samples
against saved arrays and the independently executed cold CLI replay. Neither
adapter changes inference parameters, waveforms, truth or scoring rules.

An initial gallery export was interrupted with exit 143 and no Python
traceback. Its partial images/data are retained in the local run records.
`export_gallery.py` replays archived gallery calculations with only the
Matplotlib legend position fixed to upper right, avoiding expensive automatic
placement. It keeps all components, samples, selected cases and figure IDs.

Per-rate recovery, Q errors, misses, split/merge flags, extras, time-shift controls
and complete-record component plots are reported together. Two training
seeds on the same ROOTs are not independent physical repeats; the high-rate
sample is still exploratory. Contrastive real-data pretraining is outside this
round because no measured waveform dataset is being supplied here.

## Completed results

All 60 predeclared native windows are included, with 4253 threshold-qualified targets. Complete recovery requires time agreement within 10 ns and Q error within 20%. The two seeds share the same test ROOTs.

| Workflow | Complete recovery, seed 1 / seed 2 | Extra reports, seed 1 / seed 2 |
| --- | ---: | ---: |
| gru_base | 1041/1100 | 525/498 |
| gru_component | 856/832 | 512/521 |
| unetgru_base | 1124/1164 | 462/455 |
| unetgru_component | 922/830 | 475/489 |

Current WLS: 812 complete, 300 extras. Read per-rate results, charge failures and paired event gains/losses before judging a tradeoff. The remaining density, response-mismatch and event-decoding limitations apply.

Component supervision reduced complete recovery for both architectures and both seeds. Keep the original-loss U-Net–GRU as a candidate, not as a demonstrated breakthrough over the previous frozen U-Net. The post-inference, true-bin charge-head diagnostic also retains the negative charge bias; it is an oracle diagnostic, not a new blind reconstruction. Neither decoder-position errors nor mean-response mismatch alone have been isolated as the cause. The reporting-only source is retained in `reporting-source`.
