# alpha_p11b_v1

本目录是可独立编译和重跑的实验快照。`g4/`、宏、配置、运行器和检查脚本均有实际副本，不依赖项目根目录的后续修改，也不要求切换 Git 提交。复用旧结果时保留整个实验目录；Geant4、ROOT 与编译环境仍须可用，版本和依赖路径见环境及构建记录。

在本目录执行：

```bash
python3 run.py list
python3 run.py verify
python3 run.py build --jobs 2
python3 run.py run --point alpha_1e5Hz --events 2
# 确认成本后，运行完整一档：
python3 run.py run --point alpha_1e5Hz
# 或运行全部条件：
python3 run.py run --all
# 在新的输出目录重跑某次历史运行：
python3 run.py replay runs/<条件>/<运行目录>/manifest.json
python3 analysis/check_run.py runs/<条件>/<运行目录>/waveform.root
```

每次运行生成新的 `runs/<条件>/<UTC时间_随机后缀>/`，保存 ROOT、实际 `run.mac`、标准输出日志及 `manifest.json`。`--events` 与 `--seed` 的覆盖值会记录在宏和清单中。重放使用历史运行的实际事件数和种子；ROOT 文件本身有时间戳／UUID，比较重跑结果应比较物理数据分支，而不是要求文件字节一致。

默认每档 100 个独立窗，完整窗 20 μs，中央区 [5, 15) μs。只启用当前 p–¹¹B 的有效单 α 源；源注入率不等于融合反应率或形成可测信号的事件率。前后保护区不是已验证收敛的最终选择。

本阶段仅提供仿真数据入口和一致性检查；未知脉冲数的 α 堆积重建及精度评价将在阅读用户提供的文献结果后实现。生产扫描由用户提交。

冻结文件改变后 `verify` 和运行器会拒绝继续复用该快照。需要改模型／宏参数时，从主项目重新生成另一个实验目录。
