# 五方法纯α堆积模拟比较（2026-10-14）

在动态范围足够、完整含噪电压可记录的前提下，实际比较原盲模板拟合、冻结CNN、小型U-Net、局部Transformer和CNN＋物理拟合。26训练/10验证ROOT，新生成22个独立测试ROOT，旧17个公开窗仅作回归。未知K，真值只监督/评分，不输入推断。原200个保留最终测试未用。

[实际Notebook第二轮](../../../jupyter/codex_Alpha_Pileup_ML_Comparison_20261014.ipynb#ml-multimethod-20261014) · [L21采用记录](../../../docs/literature-use.md#ml-multimethod-20261014) · [逐窗结果](../../../jupyter/artifacts/alpha_ml_multimethod_20261014_180000/new_test/summary.json) · [验证](../../../jupyter/run_records/alpha_ml_multimethod_20261014_180000/validation-results.json)

## 已完成的方法与条件

CNN保留第一轮25,707参数和原权重，U-Net134,355参数、Transformer30,731参数；新网络各2400步、相同特征/损失/切片序列，验证选择800/0.75及2400/0.85。U-Net未复刻文献LSTM，Transformer为32维、4 heads、两层和局部窗口，不是SDD文献完整实现。新网络1024点窗输出中心500点，CNN全窗，容量/上下文/计算时间不同。

混合：CNN未知数量候选（250 pe门）、冻结95/351 ns模板联合非负幅度及±6 ns时间调整；500 pe输出门、最终合并<8 ns。验证只在λ=0/0.1/1中选择电荷先验，选中0.1。小矩阵NNLS，大矩阵有界Gram优化，保存数值停止/KKT；不向算法传入真值时间或数量。不能补出没有神经候选的事件。

每ROOT恰好一个20 μs `WaveformEvent`、2 ns采样，中央[5,15) μs评分；含光子、共享微胞恢复/串扰、读出孔径和0.1 mV电子噪声，无外部物理背景。q为实际归属α的满恢复雪崩当量；尚无沉积能量标定。所有方法只用`analog_mV`，不根据ADC设置拒绝/量化。时间匹配±10 ns，电荷误差≤20%。

| 源注入率 / MHz | 窗数 | 阈上目标 | 拟合 / % | CNN / % | U-Net / % | Transformer / % | CNN＋物理 / % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 8 | 8 | 100.00 | 100.00 | 87.50 | 87.50 | 100.00 |
| 1 | 2 | 10 | 100.00 | 90.00 | 80.00 | 90.00 | 100.00 |
| 10 | 2 | 148 | 57.43 | 60.81 | 57.43 | 62.84 | 62.16 |
| 20 | 2 | 283 | 38.87 | 44.52 | 43.46 | 45.58 | 45.94 |
| 40 | 2 | 594 | 17.00 | 31.82 | 35.35 | 28.96 | 32.15 |
| 60 | 2 | 761 | 10.12 | 23.39 | 25.10 | 20.37 | 23.65 |
| 80 | 2 | 1104 | 4.89 | 16.12 | 19.38 | 13.68 | 16.30 |
| 100 | 2 | 1344 | 2.46 | 17.86 | 18.38 | 15.03 | 13.17 |

每个高率两个窗，100 MHz无同率验证，仍为小样本先导；误检/误合并/误拆分、条件误差、时间移位和幅度外推失败在Notebook中完整保留。不能据此给装置最大工作率、能谱FWHM或实时吞吐。

## 源码、数据与运行溯源

`source/`保存55份训练/推断/物理及报告源码，`analysis-source.json`逐文件校验；`replay.py`优先加载冻结源码。模型核心6份源码及权重哈希在本地结果`frozen-methods.json`保存，方法冻结时间为2026-10-14T19:25:05.063410+08:00，测试无调参。`post-test-diagnostics-source.json`另保存事后固定区间/残差诊断源码，其区域为6.0–6.6 μs、20/80 MHz第一窗，未用于方法选择。

`data-split.json`记录26/10划分及22个测试种子；`manifest.json`记录全部实际宏、ROOT哈希、三个种子、版本及未提交改动。新增1个同种子环境对照＋22测试ROOT；G4串行、最多2独立任务、启动间隔≥1秒。Geant411.4.2、ROOT6.40.04、GCC16.2.1；二进制从已验证同物理版本复制，新1 MHz同种子检查物理字段逐项一致。所有ROOT继续只在项目`g4/data/`，分析产物在本地Jupyter。

共享目录原子替换清单失败后：60 MHz 00已正常保存，经原始RunInfo/实际宏/波形复核恢复完成状态；80 MHz 00原宏未启动且没有ROOT，保留该宏，使用`_retry`新标签和原种子。无有效数据删除或覆盖。故障前清单和核验在本地[恢复记录](../../../jupyter/run_records/alpha_ml_multimethod_20261014_180000/acquisition_manifest_permission_failure/recovery.json)。

## 重跑分析

从项目根目录运行，保留本地ROOT、标签导出和原CNN权重；使用新输出目录，训练拒绝覆盖已有模型。基础依赖见项目`requirements-analysis.txt`，可选CPU PyTorch见`dependency-requirements.txt`：

```bash
bash tools/run_gagg_python.sh -m pip install -r requirements-ml.txt --index-url https://download.pytorch.org/whl/cpu

mkdir -p jupyter/artifacts/alpha_ml_multimethod_REPLAY/models
cp -a jupyter/artifacts/alpha_ml_multimethod_20261014_180000/models/cnn jupyter/artifacts/alpha_ml_multimethod_REPLAY/models/

bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_multimethod_20261014/replay.py train --model unet --steps 2400 --threads 2 --output jupyter/artifacts/alpha_ml_multimethod_REPLAY
bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_multimethod_20261014/replay.py train --model transformer --steps 2400 --threads 2 --output jupyter/artifacts/alpha_ml_multimethod_REPLAY
bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_multimethod_20261014/replay.py select-hybrid --output jupyter/artifacts/alpha_ml_multimethod_REPLAY

bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_multimethod_20261014/export_tests.py --output jupyter/artifacts/alpha_ml_multimethod_inputs_REPLAY
bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_multimethod_20261014/replay.py evaluate --benchmark new_test --output jupyter/artifacts/alpha_ml_multimethod_REPLAY --data jupyter/artifacts/alpha_ml_multimethod_inputs_REPLAY
```

本轮实际PyTorch2.14.1+cpu安装于本机临时目录`/tmp/gagg-ml-deps-20261014`，运行设置`PYTHONPATH`，没有放入共享虚拟环境。换系统需在本机安装可选依赖；以上输入和权重只在本地维护，GitHub克隆不含Jupyter。仅重评估时复制本地完整模型目录、`frozen-methods.json`和`hybrid-config.json`到新结果目录后执行`evaluate`，不要重新选择参数。

报告可用冻结副本：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_multimethod_20261014/replay_report.py report --output jupyter/artifacts/alpha_ml_multimethod_REPLAY
bash tools/run_gagg_python.sh pileup/experiments/alpha_ml_multimethod_20261014/replay_report.py diagnostics --output jupyter/artifacts/alpha_ml_multimethod_REPLAY
```

运行器不覆盖已完成ROOT。另做物理模型或固定配置变化时新建实验目录，不改本快照；换机器重建冻结C++源码仍需外部Geant4/ROOT及对应物理数据集。
