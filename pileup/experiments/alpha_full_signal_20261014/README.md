# 完整电压下的纯 α 堆积恢复（2026-10-14）

用户明确本轮假设采集动态范围足够，研究计数率、脉冲长度与盲拟合恢复的关系。全部预定窗口统一使用保存的总 `analog_mV`（含电子学噪声、采样和原有 SiPM 响应），不输入 ADC 码值或触顶标记，也不增加新的量化模型。没有重跑 G4。

该目录保存独立分析快照；旧 G4 实验、ROOT 与 ADC 结果保留。输入来自两个已有结果目录，分别含 24 和 12 个独立窗口；开发窗及同种子环境对照不纳入。本次不调整拟合参数、阈值或模板，不把真实脉冲数/时间/电荷交给拟合器。

从项目根目录运行（输出目录须不存在）：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/gagg-full-signal-mpl \
  bash tools/run_gagg_python.sh pileup/experiments/alpha_full_signal_20261014/replay.py \
  --output jupyter/artifacts/alpha_full_signal_replay --workers 2
```

运行所用源码、输入结果哈希、环境与输出路径在 `manifest.json` 中记录。原宏、源率、三组种子和每个 ROOT 的哈希由输入结果及其 G4 实验 manifest 追溯。该结果是动态范围充足且忽略量化的模拟基准，不代表已验证某款 ADC；恢复的是实际输出电荷，尚未连接沉积能量标定。

研究问题、表图、限制和下一步统一记录于本地 [纯 α 笔记第14节](../../../jupyter/codex_Alpha_Pileup_Rate_Fitting.ipynb#full-signal-rates-20261014)。
