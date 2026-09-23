# 纯 α 随机堆积实验

第一步准备纯 α 输入数据，研究对象是给定计数率下随机出现的多脉冲。沿用主工程的 p–¹¹B 有效单 α 能谱和光学参数，暂不加入 brems / ¹²C γ 背景。未知脉冲数的重建将在精读用户提供的文献结果后实现；这里的检查脚本只读取仿真真值做数据一致性检查。

## 首轮条件

实验名为 `alpha_p11b_v1`。每档默认 100 个独立随机窗，每窗完整生成和记录 20 μs，中央评价区为 `[5,15)` μs。前 5 μs 提供此前粒子的响应，后 5 μs 保留中央区粒子的尾部。保护区内照常生成粒子并记录光子，后续算法应接收完整波形；边界收敛尚未验证，不能先将保护区裁掉再拟合。

| 条件 | α 注入率 / Hz | 每完整窗期望 α 数 | 每中央区期望 α 数 | 随机种子 |
| --- | ---: | ---: | ---: | ---: |
| alpha_1e5Hz | 100000 | 2 | 1 | 260923101 |
| alpha_1e6Hz | 1000000 | 20 | 10 | 260923102 |
| alpha_1e7Hz | 10000000 | 200 | 100 | 260923103 |
| alpha_1e8Hz | 100000000 | 2000 | 1000 | 260923104 |

源在完整窗内抽样 `N ~ Poisson(R × T)`，再抽样并排序均匀发射时间；100 MHz 表示平均间隔 10 ns，不是每隔 10 ns 固定发射。各窗依次使用同一随机流，分别独立抽样，不在每个窗开始时重置种子。

注入率是生成器发射的 α 粒子数率，不是融合反应率，也不是实际形成信号的事件率。当前源是原有有效单 α 能谱，不生成一次 p–¹¹B 反应的三个关联 α。源能量、晶体入口能量、沉积能量和由光输出标定的能量需要分别处理。

100 窗是初始运行规模，可用 `--events` 覆盖并自动记录。最高速率默认期望生成 20 万个 α，完整光学输运成本尚未测定；先运行小样本再决定正式统计量。

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
    runs/<条件>/<运行编号>/    # ROOT、实际宏、日志、运行清单，Git 忽略
```

每次运行创建新的输出目录；旧 ROOT 不被覆盖。运行清单保存种子、实际事件数、完整配置、宏与二进制哈希、依赖版本、耗时、返回码和 ROOT 哈希。`RunInfo` 另存实际时间窗、中央评价区、源启用状态/速率及随机引擎状态。

## 运行与重跑

需要 CMake、C++17 编译器、Geant4 和 ROOT；本地验证环境为 Geant4 11.3.2、ROOT 6.34.10。加载依赖环境后，从项目根目录执行：

```bash
python3 pileup/experiments/alpha_p11b_v1/run.py list
python3 pileup/experiments/alpha_p11b_v1/run.py verify
python3 pileup/experiments/alpha_p11b_v1/run.py build --jobs 2
python3 pileup/experiments/alpha_p11b_v1/run.py run --point alpha_1e5Hz --events 2
```

确认日志和小样本后，由用户提交完整一档或四档：

```bash
python3 pileup/experiments/alpha_p11b_v1/run.py run --point alpha_1e5Hz
python3 pileup/experiments/alpha_p11b_v1/run.py run --all
```

上述两个命令分别是单档与全部条件的运行选项。`--all` 按顺序运行各档。正式运行使用冻结宏中的种子；增加独立重复样本时，为单档指定新的 `--seed`。相同种子和事件数的重复运行是复现检查，不能当作新增独立统计量。

记录运行器打印的结果目录。重跑时，把下面路径中的 `<条件>/<运行编号>` 换成那次运行的实际路径：

```bash
python3 pileup/experiments/alpha_p11b_v1/run.py replay pileup/experiments/alpha_p11b_v1/runs/<条件>/<运行编号>/manifest.json
python3 pileup/experiments/alpha_p11b_v1/analysis/check_run.py pileup/experiments/alpha_p11b_v1/runs/<条件>/<运行编号>/waveform.root
```

`replay` 使用历史运行的实际种子和事件数，输出到新目录。ROOT 包含时间戳、UUID 和运行路径，复现判断应比较物理数据分支，不能只比较文件哈希。检查脚本需要 numpy、uproot；可用 `--output 检查结果.json` 保存报告，已有同名报告不会被覆盖。

## 保留旧实验，继续修改新实验

`g4/` 主工程继续用于开发。需要改变模型、时间窗、保护区或扫描配置时，修改主工程／配置文件，再生成一个新名称：

```bash
python3 pileup/prepare.py --config pileup/configs/alpha_p11b_v1.json --name alpha_p11b_v2
```

准备器拒绝覆盖已有实验。每个实验内的源码、宏、运行器、检查脚本都有实际副本并校验哈希；其运行不读取主工程，因而不要求 Git 回退。以后另存整个实验目录及其 `runs/` 即可保留配置与结果；换机器时重新构建，Geant4、ROOT 及所需物理数据集仍须可用，环境记录不是这些外部依赖的副本。

## 数据含义与当前限制

- 完整波形是理想 SiPM 面的到达光子序列，没有 PDE、电子学卷积、采样噪声或饱和模型。
- 中央区的真值计数目前按**源发射时间**划分；将来评价算法时需明确输运时间偏移和匹配容差。
- 检查报告暂将 `Edep > 0 且到达光子数 ≥ 20` 作为“有信号”的诊断定义，可通过 `--min-photons` 调整；这不是测量阈值、算法检出率或可恢复性的判断。
- 真值树用于检查与后续评分。重建算法的输入应只有总波形及预先建立的响应模型，不提供真实粒子数、时间或能量。
- 5 μs 保护区、统计量、时间分箱和光输出到 α 能量的标定仍需验证。尚未得到重建效率、漏检／误拆分／误合并、能量分辨率或可恢复计数率的结论。

研究过程和实际验证记录保存在本地 `jupyter/Note.ipynb` 与 `jupyter/Research_Plan_20260921.ipynb`。源码可按完成阶段做本地 Git 提交；上传须由用户另行明确发起。
