# 纯 α / S13360 中间计数率补充扫描

研究入口：[纯 α 堆积拟合笔记](../../../jupyter/codex_Alpha_Pileup_Rate_Fitting.ipynb)。
源注入率为 2、4、6、8 × 10⁷ Hz，每档预定 3 个独立窗口。
每个 ROOT 保存一个 20 μs 窗，中央 [5,15) μs 评价，关闭外部 γ 背景。
注入率不是探测计数率或物理反应产生率。

`source/` 保存实际 C++、配置和分析源码；物理与读出逐字节沿用
`alpha_s13360_rate_fit_20261013`，包括 1 V ADC 上限和 8000 万候选计算保护。
该保护值不改变器件概率。本轮仅增加中间率，不修改旧实验。
`macros/`、`logs/` 和 `manifest.json` 保存宏、种子、ROOT 映射、版本、资源消耗和状态。

上次数据来自 Debian / ROOT 6.34.10，本轮使用 Fedora / ROOT 6.40.04、
Geant4 11.4.2、GCC 16.2.1。先重放旧 1 MHz 开发窗和第一个 10 MHz 比较窗，
检查实际 RunInfo、输出结构和同种子物理/读出结果；这两窗保存在 `checks`，
不计为新的独立比较样本。正式运行最多两个独立进程，G4 内部串行，
启动间隔至少一秒；并发数根据检查窗内存决定。

从项目根目录构建和检查（依赖路径按本机调整）：

```bash
cmake -S pileup/experiments/alpha_s13360_intermediate_rates_20261014/source/g4 \
  -B pileup/experiments/alpha_s13360_intermediate_rates_20261014/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DGeant4_DIR=/opt/geant4/11.4.2/lib64/cmake/Geant4 \
  -DROOT_DIR=/opt/root/6.40.04/cmake
cmake --build pileup/experiments/alpha_s13360_intermediate_rates_20261014/build -j 2
ctest --test-dir pileup/experiments/alpha_s13360_intermediate_rates_20261014/build --output-on-failure
bash tools/run_gagg_python.sh pileup/experiments/alpha_s13360_intermediate_rates_20261014/run_campaign.py --stage checks
```

检查通过并记录 `checks_validated` 后，运行同一脚本的 `--stage production --workers 2`。
已完成标签会跳过，失败或中断记录要求先检查，输出不会被静默覆盖。
重放已有比较窗可使用同目录 `replay.py <label>`；重放结果写入新 ROOT 并单独列入
`replays`，同种子重复不增加独立样本量。

分析沿用原固定模板、阈值和时间/电荷成功标准，同时报告 ADC 削顶、拒绝窗比例、
恢复率、漏检和误拆分/合并。未削顶模拟电压仅作诊断，不能替代真实 ADC 的恢复结果。
运行进度与验证以 manifest 和研究 Notebook 为准。

分析和出图（从项目根目录运行）：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bash tools/run_gagg_python.sh pileup/experiments/alpha_s13360_intermediate_rates_20261014/analyze.py --wait
bash tools/run_gagg_python.sh pileup/experiments/alpha_s13360_intermediate_rates_20261014/plot_results.py
```

`--wait`只等待已经排定的12份输出，不启动仿真。同种子环境检查不计入结果；
部分结果标为进行中，全部完成后才产生最终合并结果。

本轮第二个20 MHz窗口发现两个光子的整数ps时间和横坐标相同，原归属检查
保护性中止。`analysis_source/`另存修复后的评分导出器：必要时用完整时间、
x/y和能量的精确保存位值做唯一匹配，仍有歧义则报错，不按行顺序猜归属。
原导出器保留于`source/`，新版本哈希见`analysis-source.json`；盲拟合函数、
固定配置和G4物理未变。波形输入与评分真值仍分开保存。

结果在`jupyter/artifacts/alpha_intermediate_rates_20261014_114000/`，
验证、修改前Notebook副本和交付记录在同名`jupyter/run_records/`目录。
结果使用实现的共享微胞输出电荷，不直接解释为沉积能量分辨率。

GitHub 上传遵守项目的本地笔记约定：`jupyter/` 及此快照内的
`source/jupyter/` 校准输出副本只在本地保留。现有本地快照及其哈希不变；
从 GitHub 克隆后执行 `analyze.py`，还需提供原有本地 ROOT 和校准输出，
将 `jupyter/artifacts/alpha_gamma_fit_20261013_165844/summary.json` 的已核验副本
放回快照中对应的 `source/jupyter/artifacts/` 路径，其 SHA256 须与
`manifest.json` 的 `source_sha256` 记录一致。
