# 纯 α 随机堆积实验

当前生产采用 **G4 串行、一窗一 ROOT、最多同时运行 8 个独立任务**。每档 100 个窗对应 100 个 ROOT 文件。新入口为 [submit_scan.py](submit_scan.py)，实际数据统一写入项目的 `g4/data/`；文件名为 `YYYYMMDD_HHhMMmSSs.root`，启动至少间隔 1 秒，同名则等待。

第一步准备纯 α 输入数据，研究对象是给定计数率下随机出现的多脉冲。沿用主工程的 p–¹¹B 有效单 α 能谱和光学参数，暂不加入 brems / ¹²C γ 背景。未知脉冲数的重建将在精读用户提供的文献结果后实现；这里的检查脚本只读取仿真真值做数据一致性检查。

## 首轮条件

当前提交器使用 `alpha_p11b_v2`；原 `alpha_p11b_v1` 快照仍可独立重跑。v2 仅修改 RootIO 的文件命名和停写重复的 `data.log`，物理源码及配置保持一致，默认宏改为每次一个窗。每档计划 100 个独立随机窗，每窗完整生成和记录 20 μs，中央评价区为 `[5,15)` μs。前 5 μs 提供此前粒子的响应，后 5 μs 保留中央区粒子的尾部。保护区内照常生成粒子并记录光子，后续算法应接收完整波形；边界收敛尚未验证，不能先将保护区裁掉再拟合。

| 条件 | α 注入率 / Hz | 每完整窗期望 α 数 | 每中央区期望 α 数 | 随机种子 |
| --- | ---: | ---: | ---: | ---: |
| alpha_1e5Hz | 100000 | 2 | 1 | 260923101 |
| alpha_1e6Hz | 1000000 | 20 | 10 | 260923102 |
| alpha_1e7Hz | 10000000 | 200 | 100 | 260923103 |
| alpha_1e8Hz | 100000000 | 2000 | 1000 | 260923104 |

源在完整窗内抽样 `N ~ Poisson(R × T)`，再抽样并排序均匀发射时间；100 MHz 表示平均间隔 10 ns，不是每隔 10 ns 固定发射。当前每个窗使用不同种子的独立串行进程，种子规则见下文；多窗历史实验则在一次运行中依次使用同一随机流。

注入率是生成器发射的 α 粒子数率，不是融合反应率，也不是实际形成信号的事件率。当前源是原有有效单 α 能谱，不生成一次 p–¹¹B 反应的三个关联 α。源能量、晶体入口能量、沉积能量和由光输出标定的能量需要分别处理。

100 窗是每档的初始总规模，可用新提交器的 `--files-per-rate` 调整；每次 G4 运行固定 `/run/beamOn 1`。最高速率默认期望生成 20 万个 α，完整光学输运成本尚未测定；先运行小样本再决定正式统计量。

## 文件与结果

```text
pileup/
  configs/alpha_p11b_v1.json    # 准备新实验的输入
  prepare.py                  # 只创建新实验，不运行 Geant4
  tools/run_experiment.py      # 复制到实验目录成为 run.py
  analysis/check_run.py        # 数据一致性与源速率检查
  experiments/alpha_p11b_v1/
    g4/                       # 实际源码副本，含光产额表和原有宏
    macros/alpha_*.mac         # 四档自包含的扫描宏
    config.json
    bundle.json               # 源版本、工作区状态、全部冻结文件 SHA256
    environment_prepare.json
    source_changes.patch
    run.py
    analysis/check_run.py
    build/                    # 本实验的构建目录，Git 忽略
    runs/<条件>/<运行编号>/    # 早先测试的宏、日志与清单；ROOT 已迁至 g4/data
  submit_scan.py              # 一窗一文件的独立任务提交器
  results/submissions/        # 各批次的脚本副本、实际宏、日志、清单和状态，Git 忽略
```

每次运行创建新的输出目录；旧 ROOT 不被覆盖。运行清单保存种子、实际事件数、完整配置、宏与二进制哈希、依赖版本、耗时、返回码和 ROOT 哈希。`RunInfo` 另存实际时间窗、中央评价区、源启用状态/速率及随机引擎状态。

## 运行与重跑

需要 CMake、C++17 编译器、Geant4 和 ROOT；本地验证环境为 Geant4 11.3.2、ROOT 6.34.10。加载依赖环境后，从项目根目录执行：

```bash
python3 pileup/experiments/alpha_p11b_v2/run.py verify
python3 pileup/experiments/alpha_p11b_v2/run.py build --jobs 2
python3 pileup/submit_scan.py --files-per-rate 100 --workers 8
```

一次只启动一个生产批次，以保持总并发上限。当前已经提交的批次在本机 tmux 后台运行；不需要再运行上面的提交命令。每个任务调用冻结实验的串行 `gagg`，独立设置种子、秒级输出路径和 `/run/beamOn 1`。冻结目录内早先的 100 窗宏保留用于历史记录，当前实际运行宏由新提交器保存。

提交器打印 `pileup/results/submissions/<批次>/`。其中 `submission.json` 保存任务安排，`status.json` 持续更新，`runs/<条件>/<文件序号>/` 保存实际 `run.mac`、`manifest.json`、运行日志和检查结果；ROOT 本体统一位于 `g4/data/`。每份完成后自动检查恰好一个 WaveformEvent、时间窗、种子、速率与真值树一致性。

默认种子为各速率的原始种子加 `200 + 1000 × 文件序号`。新增独立批次应换用不重叠的 `--seed-offset`；相同种子的重跑是复现检查，不增加独立统计量。若要保留在后台运行，可先加 `--prepare-only`，再在 tmux 中执行所打印目录下的 `run_scan.py --submission <该目录>`。在批次目录建立 `STOP_AFTER_CURRENT` 文件可让正在运行的单窗完成后停止后续任务。

继续一个已经停止的单窗批次，只提交尚未完成的任务，保留原文件序号和种子：

```bash
python3 pileup/submit_scan.py --resume pileup/results/submissions/<已停止批次> --workers 8
```

续跑建立新提交目录，并通过 `campaign_id` 将前后结果关联；已完成窗口不重复运行。跨源码快照续跑仅允许已核对的 RootIO 改动，物理源码和配置必须一致。此次文件整理前已完成 99 个窗口，余下 301 个续跑，共同组成四档各 100 个窗口。相同种子的 v1/v2 小样本已核对物理分支完全一致。

重跑一个已完成的单窗文件，指定其实际清单；输出仍用新的秒级文件名写入 `g4/data/`：

```bash
python3 pileup/submit_scan.py --replay pileup/results/submissions/<批次>/runs/<条件>/<文件序号>/manifest.json --workers 1
```

不同快照需同时指定 `--experiment`。ROOT 包含时间戳、UUID 和路径，复现应比较物理分支，不能只比较文件哈希。

本地 `jupyter/Run_Catalog.ipynb` 统一列出 ROOT、实际宏、注入率、种子、窗口、状态和说明，CSV／JSON 导出位于 `jupyter/run_records/`。执行索引单元即可刷新；历史文件缺失原始宏时明确标注，不以当前同名宏替代。经用户授权，已删除被替代的多窗测试、中断 ROOT 及空 ROOT；删除清单、旧名到新名的映射和 `data.log` 原文均在 Jupyter 中归档，实际宏及运行清单保留。单粒子图与模板位于 `jupyter/artifacts/`，`g4/data/` 只存 ROOT。

## 保留旧实验，继续修改新实验

`g4/` 主工程继续用于开发。需要改变模型、时间窗、保护区或扫描配置时，修改主工程／配置文件，再生成一个新名称：

```bash
python3 pileup/prepare.py --config pileup/configs/alpha_p11b_v2.json --name alpha_p11b_v3
```

准备器拒绝覆盖已有实验。每个实验内的源码、宏、运行器、检查脚本都有实际副本并校验哈希；其运行不读取主工程，因而不要求 Git 回退。以后另存整个实验目录及其 `runs/` 即可保留配置与结果；换机器时重新构建，Geant4、ROOT 及所需物理数据集仍须可用，环境记录不是这些外部依赖的副本。

## 数据含义与当前限制

- 完整波形是理想 SiPM 面的到达光子序列，没有 PDE、电子学卷积、采样噪声或饱和模型。
- 中央区的真值计数目前按**源发射时间**划分；将来评价算法时需明确输运时间偏移和匹配容差。
- 检查报告暂将 `Edep > 0 且到达光子数 ≥ 20` 作为“有信号”的诊断定义，可通过 `--min-photons` 调整；这不是测量阈值、算法检出率或可恢复性的判断。
- 真值树用于检查与后续评分。重建算法的输入应只有总波形及预先建立的响应模型，不提供真实粒子数、时间或能量。
- 5 μs 保护区、统计量、时间分箱和光输出到 α 能量的标定仍需验证。尚未得到重建效率、漏检／误拆分／误合并、能量分辨率或可恢复计数率的结论。

研究过程和实际验证记录保存在本地 `jupyter/Note.ipynb` 与 `jupyter/Research_Plan_20260921.ipynb`。数据索引见本地 `jupyter/Run_Catalog.ipynb`。源码可按完成阶段做本地 Git 提交；上传须由用户另行明确发起。
