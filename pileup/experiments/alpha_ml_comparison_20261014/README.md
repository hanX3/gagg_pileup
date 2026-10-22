# 纯 α pile-up：时序卷积网络与盲拟合先导比较

本实验按完整 ROOT 划分原有36个窗口为26训练、10验证，再生成17个新种子测试窗口。输入为完整含噪 `analog_mV`，假定采集动态范围足够；沿用实际共享微胞 SiPM 状态，不考虑ADC削顶或量化。网络和拟合都从电压估计未知事件数、时间和电荷，真值仅用于监督训练和分别评分。

完整过程与图表在本地 [独立 Notebook](../../../jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb)；[采用记录](../../../docs/literature-use.md#ml-pileup-comparison-20261014) 说明文献方向和本项目的实际结果。它是自建小型1D膨胀残差CNN，含固定逆滤波等电压特征，不是论文U-Net/SAEFit完整实现的复现。

## 数据与配置

- [data-split.json](data-split.json)：按完整ROOT、物理/SiPM/电子学种子划分，先于切片；原先200个保留最终测试ROOT没有使用。100 MHz只有一个旧窗，全部训练，缺少同率验证。
- [manifest.json](manifest.json)：原生G4源码、器件配置、二进制哈希、实际宏/日志、种子、ROOT和资源记录。新测试为0.1 MHz 8窗，1/10 MHz各2窗，其余20/40/60/80/100 MHz各1窗。每窗20 μs，中央[5,15) μs评价。
- [calibration-settings.json](calibration-settings.json)：与原完整电压拟合一致的独立模板、固定基线与评分配置。
- [source](source)：实际C++、Python和配置快照；[analysis-source.json](analysis-source.json) 保存校验哈希。沿用已通过同种子检查的G4二进制，完整源码可独立重编译，不依赖切换Git提交。外部Geant4、ROOT及项目材料/能谱输入仍需提供。
- 模型输出是实际电荷pe与首有效面光子时间，尚未标定为沉积能量；q≥500 pe、时间误差≤10 ns、相对电荷误差≤20%才计成功。密集波形的同一2 ns网格及8 ns候选最小间隔仍可能无法区分事件。

## 重跑分析

从项目根目录执行。基础环境沿用 `tools/` 的各系统本机环境，ML是可选依赖：

```bash
bash tools/run_gagg_python.sh -m pip install -r requirements-ml.txt --index-url https://download.pytorch.org/whl/cpu
```

本轮实际PyTorch 2.14.1+cpu安装于本机临时目录 `/tmp/gagg-ml-deps-20261014`，使用 `PYTHONPATH=/tmp/gagg-ml-deps-20261014`；重跑可在本机项目环境安装以上可选依赖，不要共享跨系统虚拟环境。版本记录位于本地 `jupyter/run_records/alpha_ml_comparison_20261014_162000/environment.json`。

训练读取已冻结的划分和输入路径。新输出目录须不存在，禁止覆盖旧模型/结果：

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bash tools/run_gagg_python.sh \
  pileup/experiments/alpha_ml_comparison_20261014/replay.py train \
  --output jupyter/artifacts/alpha_ml_comparison_REPLAY --steps 2400 --threads 2

bash tools/run_gagg_python.sh \
  pileup/experiments/alpha_ml_comparison_20261014/export_tests.py \
  --output jupyter/artifacts/alpha_ml_test_inputs_REPLAY

env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bash tools/run_gagg_python.sh \
  pileup/experiments/alpha_ml_comparison_20261014/replay.py evaluate \
  --output jupyter/artifacts/alpha_ml_comparison_REPLAY \
  --data jupyter/artifacts/alpha_ml_test_inputs_REPLAY

bash tools/run_gagg_python.sh -m pileup.analysis.alpha_ml_report \
  jupyter/artifacts/alpha_ml_comparison_REPLAY
```

训练每200步在验证集选择全率通用检出阈值，保存 `frozen-model.json`、权重与配置哈希。独立测试前冻结，测试数据和标签不调参。读取Notebook中保存的表/图只需基础环境，不会自动重训或启动G4。

## 控制与限制

本地记录保留初始化失败：首个过大的整数种子被G4拒绝，未生成ROOT；更正为261015000系列后重新启动。缺sympy的依赖探测也未训练权重；补齐依赖后才正式训练。新1 MHz同种子检查与原物理分支完全一致，该窗不作为独立测试。

另外保存空、单、低阈、超训练电荷幅度、10/20/50/100/200 ns双脉冲控制。它们是平均模板加噪的软件检查，不等同真实共享微胞双脉冲分辨率；网络失败也保留。测试后用中心预测事件的1–9 μs循环时间移位检验偶然时间匹配，保持预测数量与电荷，不用于调参。

高率每档仅1个独立新窗、训练总量有限，结果是模拟先导，不是实物性能或通用最大工作率；未知背景、模板/噪声失配及电荷到能量映射待验证。时间打乱的多次抽样不是新增独立ROOT，不给窗口间置信区间。推断计时含特征预处理，不含ROOT读取、训练和初始库加载；双侧上下文的离线运行不能解释为在线吞吐保证。
