# Pure-alpha ML comparison with current fitting

本轮重做本地 Notebook 0004 的共同比较。保留 CNN、U-Net、Transformer、GRU 的已有权重和验证阈值，加入当前对角 WLS，以及新的 CNN 候选＋残差模型搜索 WLS。旧 CNN＋物理拟合继续作为对照。

新的混合方法把原观测波形 pilot 与 CNN 候选合并；CNN 候选保留训练所需的262点边缘上下文，250 pe仅作候选门。最终电荷由完整波形的非负 WLS 重新求解，没有神经电荷先验；经验方差固定于原观测 pilot，局部搜索允许增峰、拆分、合并和删除。阈值、K代价和搜索预算继承已有验证选择，没有用本轮测试再选参数。这是多个改动的整体比较，不单独归因于某一改动。

- 输入为20 μs完整含噪模拟电压、2 ns采样；纯p–¹¹B有效单α谱、无γ背景、固定95/351 ns、实际共享微胞。假设采集动态范围充分。
- 中央[5,15) μs评价，保护区保留；Q≥500 p.e.目标，时间±10 ns一对一匹配及电荷±20%共同通过。Q为实际雪崩等效电荷，尚无沉积能量标定。
- 原10个验证ROOT仅作冻结配置诊断。已有23个公开窗口作为共同复评，不当第二批独立重复。
- 确认测试在源码、权重、配置冻结后生成16个新种子ROOT；0.1/1/10/20/40/60/80/100 MHz为8/2/1/1/1/1/1/1窗。每档高率只有一个新窗口，不能据此确定稳定排名或最高可用率。
- 未重新训练网络。100 MHz仍没有同率神经验证集；模板、噪声及模型失配的边界保留。

`inference-source/`保存实际分析源码，`assets/`保存冻结物理配置，`models/`保存全部实际权重、阈值及训练溯源；`source/g4/`、`bin/gagg`、`macros/`和`manifest.json`保留本轮原生生产条件。复跑不需要切换Git提交。ROOT统一在项目根目录`g4/data/`，输出和Notebook只在本地Jupyter目录维护。

从项目根目录运行已准备实验的阶段（已完成结果会校验后续用）：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_refresh_20261021/replay.py validation
bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_refresh_20261021/replay.py common
bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_refresh_20261021/replay.py confirmation
```

真正重新推断、拟合到一个不存在的新输出目录：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_refresh_20261021/replay_analysis.py \
  --output jupyter/artifacts/alpha_ml_refresh_replay_new --phase confirmation
```

可加`--one <ROOT文件名主干>`只复跑一窗。上述分析需要本机CPU PyTorch及原本地输入；不启动G4、不训练、不改原ROOT。已完成评分会与原结果逐项核对。

实际记录包括`jupyter/artifacts/alpha_ml_refresh_20261021/`中的协议、冻结哈希、逐窗预测/模型、评分、控制与图像，以及`jupyter/run_records/alpha_ml_refresh_20261021/`中的预检查、原工作区差异、日志和交付验证。文献依据沿用项目L21、L40–L42及其人工审核边界；不把其他装置中的性能数字移植为GAGG能力。
