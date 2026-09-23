# GAGG 探测器响应与脉冲堆积模拟

本项目使用 Geant4 模拟 GAGG(Ce) 中的能量沉积、闪烁发光及光子输运，记录到达虚拟 SiPM 的光子时间序列，再用 Python 分析波形。研究对象是 p–¹¹B 相关的 α 粒子、热轫致辐射背景和 ¹²C 俘获/退激 γ，目标是研究计数率升高时的随机堆积，以及 α 数量和能谱的恢复能力。

当前波形的纵轴是每个时间 bin 的到达光子数；实际 SiPM 的 PDE、雪崩、串扰、后脉冲及读出电子学尚未建模。已有分粒子波形来自仿真真值，自动盲重建算法仍待实现与评估。

## 项目目录

| 目录 | 内容 |
| --- | --- |
| [g4](g4/) | Geant4 工程，包含源码、宏、分析脚本和结果图。 |
| [pileup](pileup/) | 纯 α 计数率实验：配置、可独立重跑的源码快照、运行器和数据检查。 |
| [11BH_alpha_spectra](11BH_alpha_spectra/) | 简化 p–¹¹B 单 α 能谱、生成脚本、CSV 和检查图。 |
| [bremsstrahlung](bremsstrahlung/) | 参数化热轫致辐射光子谱、生成脚本、CSV 和检查图。 |
| [stopping_power](stopping_power/) | 阻止本领表、Birks 光产额曲线和生成程序。 |

ROOT 原始数据、构建产物和本地研究笔记不随仓库分发，排除规则见 [.gitignore](.gitignore)。

## 当前模型

模拟流程：入射粒子 → GAGG 内能量沉积 → 闪烁光子产生与传播 → SiPM 到达时间 → ROOT → 时间直方图。

- GAGG 尺寸为 10 × 10 × 1 mm，入口位于 z = 250 mm；前 Mylar 膜厚 5 μm，侧膜厚 20 μm。
- 背面虚拟 SiPM 尺寸为 10 × 10 × 0.01 mm，记录进入的光子并终止其轨迹。
- 堆积模式下，一个 Geant4 event 是一个时间窗。各源独立抽样 N ~ Poisson(R × T)，再在窗内均匀抽样发射时间。
- 默认 T = 10 μs，可用 `/gagg/run/windowLength` 配置。纯 α 首轮实验使用完整 20 μs 窗和中央 10 μs 评价区。部分历史 ROOT 树标题写着 1 ms；解释历史数据时使用其 RunInfo.t_length_ps。
- 常用宏启用圆盘源并将粒子瞄准晶体前表面；速率表示有效注入粒子/光子数率，不直接等于等离子体总产额或探测事件率。
- α 使用有效单粒子能谱；¹²C γ 使用单光子谱线混合，当前没有生成完整三 α 关联事件或 γ 级联符合对。

几何、时间窗、光学与固定能谱参数见 [Constants.hh](g4/include/Constants.hh)；源抽样见 [PrimaryGeneratorAction.cc](g4/src/PrimaryGeneratorAction.cc)；光产额与边界设置见 [DetectorConstruction.cc](g4/src/DetectorConstruction.cc)。

## 编译与运行

需要 CMake、C++ 编译器、带 UI/可视化组件的 Geant4，以及 ROOT 的 Core/RIO/Tree。已验证的本地环境为 Geant4 11.3.2、ROOT 6.34.10。

从项目根目录运行，先按本机安装方式加载 Geant4/ROOT 环境：

```bash
cmake -S g4 -B g4/build
cmake --build g4/build -j
cd g4/build
./gagg ../macros/run_all_pileup.mac
```

未指定输出文件时，程序使用相对目录 `../data`，因此上述运行方式将新结果写入 `g4/data/`。宏可用 `/gagg/run/outputFile` 指定相对于启动目录的 ROOT 路径；已有同名文件会被拒绝覆盖。无参数运行 `./gagg` 会进入交互可视化模式。

当前 [run_all_pileup.mac](g4/macros/run_all_pileup.mac) 配置为轫致辐射 10⁶ Hz、¹²C γ 10 Hz、α 10⁶ Hz，生成 10 个时间窗。每窗期望初级粒子数约 20。另有单源宏与 [低速率三源宏](g4/macros/run_all_p11b_sources_pileup.mac)。大规模运行前应先用独立输出目录和少量事件检查配置。

## 纯 α 堆积实验

首轮沿用当前 p–¹¹B 有效单 α 能谱，扫描注入率 10⁵、10⁶、10⁷、10⁸ Hz，关闭光子背景。准备与运行方法见 [pileup/README.md](pileup/README.md)。每个实验目录保存自己的源码、宏、运行器和配置；以后修改主工程时，旧实验仍能自行编译重跑，无需切换 Git 提交。Geant4、ROOT 等外部依赖须保持可用。

运行接口支持以下命令，均应在 `/run/beamOn` 前设置；`preWindow` 和 `postWindow` 定义中央评价区，不会删除保护区内的粒子或光子：

```text
/gagg/run/seed 260923101
/gagg/run/windowLength 20 us
/gagg/run/preWindow 5 us
/gagg/run/postWindow 5 us
/gagg/run/outputFile waveform.root
```

未知脉冲数的 pileup 重建尚未实现；当前阶段先准备可追溯的输入数据，待精读用户提供的文献结果后再实现算法。

## 分析已有数据

Python 绘图依赖 numpy、matplotlib、uproot、awkward；源谱辅助脚本还使用 pandas。已核对的本地读取环境为 Python 3.10.12、uproot 5.7.4、awkward 2.9.1、numpy 1.26.4、matplotlib 3.5.1；尚未建立依赖锁定文件。

取得对应 ROOT 数据后，从项目根目录运行：

```bash
cd g4/analysis
python3 plot_pileup_waveform.py ../data/gagg_waveform_20260624_16h28m34s.root --bin-width-ns 1.0
```

脚本目前只画第一个 event，PNG 保存到启动命令所在目录。它会读取完整事件树，大文件分析需要预留内存。本机 uproot 默认文件读取曾出现等待，显式使用 `uproot.source.file.MemmapSource` 可读取。

单粒子比较示例与分析工具说明见 [分析目录说明](g4/analysis/README.md)。

## ROOT 数据与运行记录

| 树 | 每行含义 | 主要字段 |
| --- | --- | --- |
| RunInfo | 一次运行的信息 | 时间窗和中央评价区 ps、事件数、实际种子与随机引擎状态、步长、宏路径与内容、各源启用状态及注入率。 |
| WaveformEvent | 一个完整时间窗 | 初级粒子数、总沉积能量 MeV、全部到达光子的时间 ps 与位置 mm。 |
| PrimaryPhoton | 一个初级粒子及其后代贡献，包含 α | 粒子/源标签、源能量和发射时间、沉积能量、闪烁光子数、到达光子的时间与位置。 |

ps 除以 1000 得到 ns。`edep_total_MeV` 是沉积能量真值，不是淬灭及光收集后的可见能量。模板形状、绝对光输出和能量标定须分别处理。

新运行将发射与光子到达时间向下量化到 1 ps，记录范围为 `[0,T)`；旧版使用最近整数舍入。差异小于 1 ps，物理输运和发光模型不变。`primary_time_ps` 是源发射时间，并非进入晶体或首次沉积的时刻。旧文件可能没有新增的 RunInfo 字段。

本地现有多源样本均只有一个 10 μs 窗；下表来自各自 ROOT 的 RunInfo 和 PrimaryPhoton，不由当前宏推断。完整文件名为 `gagg_waveform_<时间戳>.root`，这些 ROOT 文件不随 Git 仓库分发。

| 时间戳 | brems / ¹²C γ / α 注入率（Hz） | 初级粒子数 | 实际 brems / ¹²C γ / α 数 |
| --- | --- | ---: | --- |
| 20260624_16h28m34s | 10⁶ / 10 / 10⁶ | 19 | 10 / 0 / 9 |
| 20260624_16h28m48s | 10⁷ / 10² / 10⁷ | 206 | 102 / 0 / 104 |
| 20260624_16h29m17s | 10⁸ / 10³ / 10⁸ | 1993 | 977 / 0 / 1016 |

单窗示例不能给出可靠的重建效率或误差统计。后续主实验应使用多段独立随机流，处理窗前历史和窗末截断，评价 α 效率、漏检、误拆分/误合并、假 α 率、能量偏差/分辨率和整体能谱偏差。

每次运行记录代码提交号及未提交改动、完整宏、随机种子、事件数、源谱和注入率、时间窗、关键参数、依赖环境和输出路径。`data/data.log` 与 RunInfo 不能代替完整源码和参数记录；现有历史数据的几何、光学参数及环境仍有追溯缺口。
