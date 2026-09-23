# GAGG 探测器响应与脉冲堆积模拟

本项目使用 Geant4 模拟 GAGG(Ce) 中的能量沉积、闪烁发光及光子输运，记录到达虚拟 SiPM 的光子时间序列，再用 Python 分析波形。研究对象是 p–¹¹B 相关的 α 粒子、热轫致辐射背景和 ¹²C 俘获/退激 γ，目标是研究计数率升高时的随机堆积，以及 α 数量和能谱的恢复能力。

当前波形的纵轴是每个时间 bin 的到达光子数；实际 SiPM 的 PDE、雪崩、串扰、后脉冲及读出电子学尚未建模。已有分粒子波形来自仿真真值，自动盲重建算法仍待实现与评估。

## 目录与版本

| 目录 | 内容 |
| --- | --- |
| [g4](g4/) | 当前 Geant4 工程，基于原 g4_260623；包含源码、宏、分析脚本和历史小型结果图。 |
| [11BH_alpha_spectra](11BH_alpha_spectra/) | 简化 p–¹¹B 单 α 能谱、生成脚本、CSV 和检查图。 |
| [bremsstrahlung](bremsstrahlung/) | 参数化热轫致辐射光子谱、生成脚本、CSV 和检查图。 |
| [stopping_power](stopping_power/) | 阻止本领表、Birks 光产额曲线和生成程序。 |

整个项目使用一个 Git 仓库。旧工程统一以 `g4/` 为路径，使用以下标签回溯：

| 标签 | 对应内容 |
| --- | --- |
| `g4_260621` | 早期光学响应、波长谱检查版本。 |
| `g4_260622` | 粒子相关光产额、α/γ 单粒子波形与 PSD 分析版本。 |
| `g4_260623` | 随机多源堆积版本。 |

标签保存的是 **2026-09-23 迁移时各旧目录的工作区快照**，包括当时尚未提交的源码和分析结果；标签名称沿用旧目录名，不表示已完整恢复那个日期的运行环境。原有提交也作为历史父提交保留，路径统一到 `g4/` 并排除本地笔记，因此提交号发生变化；提交说明中的 `Original-Commit` 记录原编号。

`main` 在最新仿真版本基础上汇集早期可复用的光学、单粒子分析脚本及结果图。C++ 源码、物理参数和当前宏保持原 g4_260623 的内容。历史图的来源见 [分析目录说明](g4/analysis/README.md)。

查看差异或在独立目录运行旧版本：

```bash
git log --first-parent --oneline --decorate
git diff g4_260622 g4_260623 -- g4/
git worktree add --detach .worktrees/g4_260622 g4_260622
```

旧版本的数据文件需另行取得。Git 标签不包含 ROOT 原始数据、构建产物和本地笔记。`AGENTS.md`、`jupyter/`、所有 `.ipynb`、`papers/` 与本地归档通过 `.gitignore` 保留在本机。

## 当前模型

模拟流程：入射粒子 → GAGG 内能量沉积 → 闪烁光子产生与传播 → SiPM 到达时间 → ROOT → 时间直方图。

- GAGG 尺寸为 10 × 10 × 1 mm，入口位于 z = 250 mm；前 Mylar 膜厚 5 μm，侧膜厚 20 μm。
- 背面虚拟 SiPM 尺寸为 10 × 10 × 0.01 mm，记录进入的光子并终止其轨迹。
- 堆积模式下，一个 Geant4 event 是一个时间窗。各源独立抽样 N ~ Poisson(R × T)，再在窗内均匀抽样发射时间。
- 当前源码 T = 10 μs。部分历史宏、帮助文字及 ROOT 树标题仍写着 1 ms；解释历史数据时使用其 RunInfo.t_length_ps。
- 常用宏启用圆盘源并将粒子瞄准晶体前表面；速率表示有效注入粒子/光子数率，不直接等于等离子体总产额或探测事件率。
- α 使用有效单粒子能谱；¹²C γ 使用单光子谱线混合，当前没有生成完整三 α 关联事件或 γ 级联符合对。

几何、时间窗、光学与固定能谱参数见 [Constants.hh](g4/include/Constants.hh)；源抽样见 [PrimaryGeneratorAction.cc](g4/src/PrimaryGeneratorAction.cc)；光产额与边界设置见 [DetectorConstruction.cc](g4/src/DetectorConstruction.cc)。

## 编译与运行

需要 CMake、C++ 编译器、带 UI/可视化组件的 Geant4，以及 ROOT 的 Core/RIO/Tree。迁移前工作区使用 Geant4 11.3.2、ROOT 6.34.10；这不代表所有历史数据的生成环境已被完整追溯。

从项目根目录运行，先按本机安装方式加载 Geant4/ROOT 环境：

```bash
cmake -S g4 -B g4/build
cmake --build g4/build -j
cd g4/build
./gagg ../macros/run_all_pileup.mac
```

程序使用相对输出路径 `../data`，因此应在 `g4/build/` 中启动；新结果位于 `g4/data/`。无参数运行 `./gagg` 会进入交互可视化模式。不要复用旧目录中含绝对路径的 CMake 缓存。

当前 [run_all_pileup.mac](g4/macros/run_all_pileup.mac) 配置为轫致辐射 10⁶ Hz、¹²C γ 10 Hz、α 10⁶ Hz，生成 10 个时间窗。每窗期望初级粒子数约 20。另有单源宏与 [低速率三源宏](g4/macros/run_all_p11b_sources_pileup.mac)。大规模运行前应先用独立输出目录和少量事件检查配置。

## 分析已有数据

Python 绘图依赖 numpy、matplotlib、uproot、awkward；源谱辅助脚本还使用 pandas。已核对的本地读取环境为 Python 3.10.12、uproot 5.7.4、awkward 2.9.1、numpy 1.26.4、matplotlib 3.5.1；尚未建立依赖锁定文件。

取得对应 ROOT 数据后，从项目根目录运行：

```bash
cd g4/analysis
python3 plot_pileup_waveform.py ../data/gagg_waveform_20260624_16h28m34s.root --bin-width-ns 1.0
```

脚本目前只画第一个 event，PNG 保存到启动命令所在目录。它会读取完整事件树，大文件分析需要预留内存。本机 uproot 默认文件读取曾出现等待，显式使用 `uproot.source.file.MemmapSource` 可读取；此迁移未改变脚本的读取算法。

单粒子比较示例和各历史图的来源见 [分析目录说明](g4/analysis/README.md)。

## ROOT 数据与运行记录

| 树 | 每行含义 | 主要字段 |
| --- | --- | --- |
| RunInfo | 一次运行的信息 | 时间窗 ps、事件数、随机种子、步长限制、宏路径、源配置标签。 |
| WaveformEvent | 一个完整时间窗 | 初级粒子数、总沉积能量 MeV、全部到达光子的时间 ps 与位置 mm。 |
| PrimaryPhoton | 一个初级粒子及其后代贡献，包含 α | 粒子/源标签、入射能量和时间、沉积能量、闪烁光子数、到达光子的时间与位置。 |

ps 除以 1000 得到 ns。`edep_total_MeV` 是沉积能量真值，不是淬灭及光收集后的可见能量。模板形状、绝对光输出和能量标定须分别处理。

本地现有多源样本均只有一个 10 μs 窗；下表来自各自 ROOT 的 RunInfo 和 PrimaryPhoton，不由当前宏推断。完整文件名为 `gagg_waveform_<时间戳>.root`，这些 ROOT 文件不随 Git 仓库分发。

| 时间戳 | brems / ¹²C γ / α 注入率（Hz） | 初级粒子数 | 实际 brems / ¹²C γ / α 数 |
| --- | --- | ---: | --- |
| 20260624_16h28m34s | 10⁶ / 10 / 10⁶ | 19 | 10 / 0 / 9 |
| 20260624_16h28m48s | 10⁷ / 10² / 10⁷ | 206 | 102 / 0 / 104 |
| 20260624_16h29m17s | 10⁸ / 10³ / 10⁸ | 1993 | 977 / 0 / 1016 |

单窗示例不能给出可靠的重建效率或误差统计。后续主实验应使用多段独立随机流，处理窗前历史和窗末截断，评价 α 效率、漏检、误拆分/误合并、假 α 率、能量偏差/分辨率和整体能谱偏差。

每次运行记录代码提交号及未提交改动、完整宏、随机种子、事件数、源谱和注入率、时间窗、关键参数、依赖环境和输出路径。`data/data.log` 与 RunInfo 不能代替完整源码和参数记录；现有历史数据的几何、光学参数及环境仍有追溯缺口。
