# 纯 α / S13360 计数率拟合先导试验

新实验，不修改 2026-09-23 的冻结扫描。研究入口为
[codex_Alpha_Pileup_Rate_Fitting.ipynb](../../../jupyter/codex_Alpha_Pileup_Rate_Fitting.ipynb)。

`source/` 保存本次实际 G4 C++、配置和运行器副本；`macros/` 为实际运行宏，
`logs/` 为完整运行输出；`manifest.json` 保存 ROOT 映射、种子、哈希、版本、环境和状态。
物理模型沿用当前主项目，只有 S13360 配置的计算保护 `max_candidates` 从 1000 万提高到 8000 万；
器件、PDE、恢复时间、电子学和 ADC 量程没有改变。增大这一保护值不改变模拟概率。

每个 ROOT 是一个 20 μs 的独立窗口，评价门为 [5,15) μs；G4 串行、进程也串行。
源为当前 p–¹¹B 有效单 α 谱、10 mm 圆盘朝向晶体，未加入 γ 或其他源。
计数率指注入 α/s，不等于探测信号率或聚变反应率。

开发窗 `pilot_1e6` 不计入比较结果。比较样本预定为：10⁵ Hz 16 窗、10⁶ Hz 4 窗、
10⁷ Hz 3 窗、10⁸ Hz 1 窗。均为小样本探索，不是最终测试集。
原来 400 窗的数据划分和保留测试集没有用于本次拟合。

## 重建与重放

从项目根目录构建保存的源码，不需要切换 Git 提交：

```bash
cmake -S pileup/experiments/alpha_s13360_rate_fit_20261013/source/g4 \
  -B pileup/experiments/alpha_s13360_rate_fit_20261013/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DGeant4_DIR=/opt/geant4/11.4.2/lib/cmake/Geant4 \
  -DROOT_DIR=/opt/root/6.34.10/cmake
cmake --build pileup/experiments/alpha_s13360_rate_fit_20261013/build -j 2
ctest --test-dir pileup/experiments/alpha_s13360_rate_fit_20261013/build --output-on-failure
.venv/bin/python pileup/experiments/alpha_s13360_rate_fit_20261013/replay.py alpha_001000000_00
```

依赖路径按当前系统记录；移机时只调整依赖位置。重放写入新的秒级 ROOT 文件名，
保留原输出并记录 `replay_of`。同种子重放不能增加独立样本量。
原宏含本次运行的绝对输出路径，不要直接再次执行原宏覆盖输出。

分析默认只读 ROOT。可复用实现位于 `pileup/analysis/alpha_pileup_fit.py` 等模块；
本次分析源码副本、初版计算上限对照、验证和 Notebook 生成脚本位于
`jupyter/run_records/alpha_rate_fit_20261013_211541/`。
