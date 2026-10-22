# 轫致辐射 γ 与 α 混合

2026-10-15 的 α/γ 原生混合波形拟合先导。α 沿用有效 p–¹¹B 单粒子谱，源注入率 1、10 MHz；γ 源注入率固定 10 MHz。每条件 8 个独立 20 μs 窗，中央 [5,15) μs 评价；另有 1 个预检查窗，不计入正式分母。源注入率不是探测事件率或物理反应率。

采用 S13360-6050PE v2，所有粒子的光子在完整窗共享微胞历史。拟合只读取完整含噪 `analog_mV`，假定动态范围充足，不施加 ADC 削顶或量化限制。固定模板来自旧的两条孤立示例；脉冲数、时刻、电荷、类型均未知。真值仅用于后续评价。电荷单位是完全恢复单胞等效 p.e.，尚未换算为沉积能量。

实际 Geant4 11.4.2、ROOT 6.40.04、GCC 16.2.1；G4 内部串行，本实验最多两个独立进程，两种 γ 合计最多四个。实际配置、三个种子、宏、源码哈希和 ROOT 映射见 [manifest.json](manifest.json)。`source/` 保留生产时的源码和运行器，`analysis_source/` 保留实际拟合源码，重跑不依赖切换 Git 提交。原始 ROOT、分析产物和 Notebook 按项目约定只在本地保留。

从项目根目录、加载本机 Geant4 / ROOT 环境后运行（将 EXP 替换为本目录路径）：

```bash
cmake -S EXP/source/g4 -B EXP/build -DCMAKE_BUILD_TYPE=Release
cmake --build EXP/build -j 2
ctest --test-dir EXP/build --output-on-failure
bash tools/run_gagg_python.sh EXP/run_campaign.py --stage checks --workers 1
# 检查实际 ROOT、原生电荷守恒和小样本拟合，记录 checks_validated 后再生产。
bash tools/run_gagg_python.sh EXP/run_campaign.py --stage production --workers 2
bash tools/run_gagg_python.sh EXP/run_analysis.py --stage production \
  --settings EXP/calibration-settings.json \
  --inputs jupyter/artifacts/alpha_gamma_fit_inputs_20261015_093500 \
  --output jupyter/artifacts/alpha_gamma_fit_20261015_093500 \
  --experiments pileup/experiments/alpha_gamma_brems_fit_20261015 pileup/experiments/alpha_gamma_mono662_fit_20261015
```

已完成的计划行会复用；不要修改旧宏、配置或覆盖 ROOT。新物理配置应创建新的实验。重拟合使用新输出路径即可保留原结果；当前结果缓存检查模板设置与拟合源码哈希。

完整过程见本地 [α/γ 混合拟合 Notebook](../../../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb)，实际文件目录见 [Run_Catalog](../../../jupyter/Run_Catalog.ipynb#alpha-gamma-mixed-20261015)，依据与限制见 [文献采用记录](../../../docs/literature-use.md#alpha-gamma-mixed-20261015)。完整恢复要求 Q≥500 p.e.、一对一时间误差≤10 ns、类型正确及电荷误差≤20%；同时报告漏检、额外候选、误拆分/误合并、条件电荷/时间误差和按独立窗重采样的不确定度。

实际写出异常与重放：γ 使用现有 exp(−E/kT)/E 源，kT=100 keV，10 keV–5 MeV，物理源码没有更改。原 20261015_09h47m13s.root 的 PrimaryPhoton/photon_energy_eV 分支写入失败；同种子重放为 20261015_10h29m25s.root。所有可读数组一致，且拟合输入与评分数组逐项相同，保留对应拟合结果并记录复用证据。 原文件保留并标记 failed_readback，替代文件不增加独立样本数。ROOT 曾在写入错误后返回退出码 0；最终采用完整分支读取核验。生产时运行器快照保留原样，可复用模块已增加 ROOT 写入错误日志识别。

`run_analysis.py` 固定 OpenBLAS / OpenMP 单线程，与本轮实际命令一致；并行应在独立窗口之间安排，不用数值库的隐式多线程改变本轮配置。

数值加速试验未用于正式统计：Cholesky缩减保持相同线性目标，但在高率完整非线性流程中，舍入误差影响后续局部搜索；虽恢复计数、分型、配对相同，最大时刻差0.103 ns、电荷差24.81 p.e.未达到预设0.05 ns/1 p.e.严格容差。原结果与试验均保留，32个正式窗口统一使用原直接NNLS。对照记录位于本地 `jupyter/run_records/alpha_gamma_fit_20261015_093500/solver-equivalence-full-windows.json`。
