# GAGG 探测器响应与脉冲堆积模拟

本项目使用 Geant4 模拟 GAGG(Ce) 中的能量沉积、闪烁发光及光子输运，记录到达虚拟 SiPM 的光子时间序列，再用 Python 分析波形。研究对象是 p–¹¹B 相关的 α 粒子、热轫致辐射背景和 ¹²C 俘获/退激 γ，目标是研究计数率升高时的随机堆积，以及 α 数量和能谱的恢复能力。

Geant4 保留到达光子，并可通过 `G4VDigitizerModule` / `G4DigiManager` 在窗末生成 SiPM 雪崩、电流、电压和 ADC，保存到同一 ROOT 的 `SiPMEvent`。当前单粒子默认器件为滨松 **S13360-6050PE**（25°C、3 V 过压），采用数据表典型参数和近似读图 PDE；恢复时间使用同系列代理值，光耦合与读出仍需标定。[SiPM 响应说明](docs/sipm-response.md) 包含统一的原生与离线入口、逐项器件依据和验证结果。原有光子计数波形与数字化电压分别保留；已有分粒子波形来自仿真真值，自动盲重建算法仍待实现与评估。

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

模拟流程：入射粒子 → GAGG 内能量沉积 → 闪烁光子产生与传播 → 完整时间窗的光子排序 → 可选 SiPM 数字化 → ROOT。后续可分析电信号，也可用保留的光子输入重新离线扫描器件参数。

- GAGG 尺寸为 10 × 10 × 1 mm，入口位于 z = 250 mm；前 Mylar 膜厚 5 μm，侧膜厚 20 μm。
- 背面虚拟光子记录层为 10 × 10 × 0.01 mm，记录进入的光子并终止其轨迹；S13360-6050PE 数字化只接收中央 6 × 6 mm 有效面内的光子，区域外光子仍保留在 ROOT。当前记录层采用理想折射率匹配，尚未加入真实封装及光学脂。
- 堆积模式下，一个 Geant4 event 是一个时间窗。各源独立抽样 N ~ Poisson(R × T)，再在窗内均匀抽样发射时间。
- 默认 T = 10 μs，可用 `/gagg/run/windowLength` 配置。纯 α 首轮实验使用完整 20 μs 窗和中央 10 μs 评价区。部分历史 ROOT 树标题写着 1 ms；解释历史数据时使用其 RunInfo.t_length_ps。
- 常用宏启用圆盘源并将粒子瞄准晶体前表面；速率表示有效注入粒子/光子数率，不直接等于等离子体总产额或探测事件率。
- α 使用有效单粒子能谱；¹²C γ 使用单光子谱线混合，当前没有生成完整三 α 关联事件或 γ 级联符合对。

几何、时间窗、光学与固定能谱参数见 [Constants.hh](g4/include/Constants.hh)；源抽样见 [PrimaryGeneratorAction.cc](g4/src/PrimaryGeneratorAction.cc)；光产额与边界设置见 [DetectorConstruction.cc](g4/src/DetectorConstruction.cc)。

## 编译与运行

需要 CMake、C++17 编译器、带 UI/可视化组件的 Geant4、ROOT 的 Core/RIO/Tree，以及 nlohmann_json ≥ 3.2（读取 SiPM 配置；优先使用独立 CMake 包，也支持 ROOT 自带头文件）。既有生产数据使用 Geant4 11.3.2、ROOT 6.34.10；2026-10-09 本机使用 Geant4 11.4.2、ROOT 6.40.04，主工程、两个冻结实验和阻止本领工具已重新构建。环境检查记录保存在本地 `jupyter/Run_Catalog.ipynb`，历史结果仍按各自的运行记录解释。

从项目根目录运行，先按本机安装方式加载 Geant4/ROOT 环境：

```bash
cmake -S g4 -B g4/build
cmake --build g4/build -j 2
cd g4/build
./gagg ../macros/run_all_pileup.mac
```

未指定输出文件时，程序使用相对目录 `../data` 和秒级文件名 `YYYYMMDD_HHhMMmSSs.root`，因此上述运行方式将新结果写入 `g4/data/`；同名时等待下一秒。宏可用 `/gagg/run/outputFile` 指定相对于启动目录的 ROOT 路径；显式指定的已有同名文件会被拒绝覆盖。无参数运行 `./gagg` 会进入交互可视化模式。

当前 [run_all_pileup.mac](g4/macros/run_all_pileup.mac) 配置为轫致辐射 10⁶ Hz、¹²C γ 10 Hz、α 10⁶ Hz，生成 10 个时间窗。每窗期望初级粒子数约 20。另有单源宏与 [低速率三源宏](g4/macros/run_all_p11b_sources_pileup.mac)。大规模运行前应先用独立输出目录和少量事件检查配置。

单 α / γ 电信号的入口如下。两个宏分别发射一个 5.486 MeV α 和一个 661.657 keV γ，20 μs 窗内在 5 μs 发射，中央正入射；每次只产生一个 ROOT。单γ宏选用固定32例样本中有信号的一例种子，便于查看电信号；该宏不是探测效率样本。固定样本运行器仍从预定种子生成并保留所有γ，包括仅含暗计数和电子学噪声的无沉积事件。

```bash
ctest --test-dir g4/build --output-on-failure
cd g4/build
./gagg ../macros/run_single_alpha_sipm.mac
./gagg ../macros/run_single_gamma_sipm.mac
```

从项目根目录运行 `.venv/bin/python -m analysis.single_particle_sipm run --alpha 4 --gamma 32` 可保存固定样本及三组随机数对照、实际源码与宏、ROOT 索引、读出交叉检查和图。默认数字化关闭，原生 SiPM 运行启用后要求 `/run/beamOn 1`。具体参数、削顶标记与模型限制见 [原生数字化](docs/sipm-response.md#native-sipm)。 本次完整改动与实际结果另见本地独立笔记 [codex_SiPM_Digitization.ipynb](jupyter/codex_SiPM_Digitization.ipynb)。

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

当前生产采用 G4 串行、每个 ROOT 一个波形窗、最多 8 个独立任务。ROOT 统一存入 `g4/data/`，使用原有秒级文件名，任务启动至少间隔 1 秒。具体入口与运行记录见 [pileup 说明](pileup/README.md)。

截至 2026-09-30，首轮四档各 100 个独立时间窗已全部完成并通过逐窗数据检查，生产 ROOT 合计约 54.27 GB。未知脉冲数的 pileup 重建尚未实现；下一步核验已有 α 模板及光输出—沉积能量标定，再建立纯 α 盲重建和计数率性能曲线。方法依据与适用条件见 [文献采用记录](docs/literature-use.md)；运行明细与研究计划在本地 Jupyter 笔记中维护。

## 分析已有数据

Python 分析及 Notebook 依赖见 [requirements-analysis.txt](requirements-analysis.txt)，包括 numpy、matplotlib、uproot、awkward、pandas、scipy 和 Notebook 执行组件。本机使用项目根目录的 `.venv/`，从项目根目录激活后运行分析：

```bash
source .venv/bin/activate
python -m unittest pileup.analysis.test_alpha_response -v
```

需要补充 Python 依赖时执行 `python -m pip install -r requirements-analysis.txt`。新机器先用其 Python 的 venv/pip 组件创建虚拟环境。2026-10-09 的本机 Python 为 3.14.4，实际包版本及单窗 G4、历史 ROOT 读取和 Notebook 执行检查见本地 [环境记录](jupyter/run_records/environment_setup_20261009_170541/environment.json)。完整 Python 版本清单在同目录的 `requirements.freeze.txt`；已在 `.venv` 中注册 `GAGG pileup (.venv)` 内核。此次仅核验运行与数据接口，未验证新旧 Geant4 的物理响应一致性。此前读取环境为 Python 3.10.12、uproot 5.7.4、awkward 2.9.1、numpy 1.26.4、matplotlib 3.5.1，保留为历史结果的环境信息。

本机继续使用经典 Jupyter Notebook 界面，五本研究笔记统一选择 `GAGG pileup (.venv)` 内核。在项目根目录运行 `bash jupyter/start_gagg_jupyter.sh`，即可打开 SiPM 笔记的经典页面。启动器使用项目 `.venv` 中的 NbClassic、Jupyter Server 和计算内核；界面与分析依赖统一由 `requirements-analysis.txt` 安装。启动配置保存在本地 `jupyter/classic_notebook_config.json`。

虚拟环境不能直接跨系统 Python 版本搬用。2026-10-13 当前系统为 Python 3.13.5，原 3.14.4 环境因缺少对应标准库而无法启动，已备份后重建。旧环境、原笔记与图像备份及本次验证记录位于本地 `jupyter/run_records/notebook_repair_20261013_112957/`。2026-10-09/10 的版本与验证记录描述当时的环境，继续保留用于追溯；旧 `.cache/jupyter-classic/` 组件不再是当前启动器的依赖。

在新机器恢复环境时，先将失效的 `.venv` 移至独立备份目录，再从项目根目录执行下列命令。共享目录不支持符号链接时，预先创建普通 `lib64` 目录，并用 `--copies` 复制解释器：

```bash
mkdir -p .venv/lib64
python3 -m venv --copies .venv
.venv/bin/python -m pip install -r requirements-analysis.txt
.venv/bin/python -m ipykernel install --prefix .venv --name gagg-pileup --display-name 'GAGG pileup (.venv)'
bash jupyter/start_gagg_jupyter.sh
```

取得对应 ROOT 数据后，从项目根目录运行：

```bash
cd g4/analysis
python3 plot_pileup_waveform.py ../data/20260624_16h28m34s.root --bin-width-ns 1.0
```

脚本目前只画第一个 event，PNG 保存到启动命令所在目录。它会读取完整事件树，大文件分析需要预留内存。本机 uproot 默认文件读取曾出现等待，显式使用 `uproot.source.file.MemmapSource` 可读取。

单粒子比较示例与分析工具说明见 [分析目录说明](g4/analysis/README.md)。

纯 α 响应核验与条件能标模块见 [pileup/analysis/alpha_response.py](pileup/analysis/alpha_response.py)，使用方法见 [pileup 说明](pileup/README.md#α-响应核验与能量标定)。首批 100 个校准/验证窗的检查、长尾模板及能标比较已完成，实际输出记录在本地 `jupyter/Alpha_Response_Calibration.ipynb`；未知脉冲数的重建仍待实现。

离线 SiPM 数字化入口为 [sipm_response.py](pileup/analysis/sipm_response.py)，使用同一 [S13360-6050PE 配置](g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json)，当前示例、验证入口和笔记均采用该器件。用完整时间窗推进所有粒子共享的 pixel 状态，另存雪崩表、波形、种子、配置及源码记录；不修改输入 ROOT。命令、概率定义、实测单雪崩模板接口和限制见 [使用说明](docs/sipm-response.md)。历史无能量字段的 ROOT 仅在明确启用常数有效 PDE 近似时接受。

## ROOT 数据与运行记录

| 树 | 每行含义 | 主要字段 |
| --- | --- | --- |
| RunInfo | 一次运行的信息 | 时间窗和中央评价区 ps、事件数、实际种子与随机引擎状态、步长、宏路径与内容、各源启用状态及注入率。 |
| WaveformEvent | 一个完整时间窗 | 初级粒子数、总沉积能量 MeV、全部到达光子的时间 ps、位置 mm 与光子能量 eV。 |
| PrimaryPhoton | 一个初级粒子及其后代贡献，包含 α | 粒子/源标签、源能量和发射时间、沉积能量、闪烁光子数、到达光子的时间、位置与光子能量 eV。 |
| SiPMEvent（启用时） | 整个时间窗对应的一个 Geant4 digi | 雪崩时间 ns、微胞、类型、父编号、原光子索引、电荷 C、恢复系数；电流 A、电压 mV、ADC码与削顶标记。 |

`RunInfo.sipm_config_json` 保存实际完整参数；另记录模块版本、随机数引擎说明和独立基础种子，`SiPMEvent` 记录该事件派生后的两个实际种子。新增电信号的粒子标签不参与数字化计算。

ps 除以 1000 得到 ns。`edep_total_MeV` 是沉积能量真值，不是淬灭及光收集后的可见能量。模板形状、绝对光输出和能量标定须分别处理。

2026-10-10 起，主工程在两棵树中新增 `photon_energy_eV`（`vector<float>`）分支，记录进入虚拟 SiPM 时 pre-step 的光子能量。`WaveformEvent` 中它与 `all_photon_arrival_t_ps/x_mm/y_mm` 按索引对应，`PrimaryPhoton` 中与 `photon_arrival_t_ps/x_mm/y_mm` 对应；能量随同一光子一起按到达时间排序。真空波长可用 `wavelength_nm = 1239.842 / photon_energy_eV` 换算。历史 ROOT 和冻结实验没有该字段，分析时先检查分支是否存在；`primary_energy_MeV` 仍表示初级粒子能量。

新运行将发射与光子到达时间向下量化到 1 ps，记录范围为 `[0,T)`；旧版使用最近整数舍入。相对连续时间的量化误差小于 1 ps，与旧版整数记录最多相差 1 ps；物理输运和发光模型不变。`primary_time_ps` 是源发射时间，并非进入晶体或首次沉积的时刻。旧文件可能没有新增的 RunInfo 字段。

本地现有多源样本均只有一个 10 μs 窗；下表来自各自 ROOT 的 RunInfo 和 PrimaryPhoton，不由当前宏推断。完整文件名为 `<时间戳>.root`，这些 ROOT 文件不随 Git 仓库分发。

| 时间戳 | brems / ¹²C γ / α 注入率（Hz） | 初级粒子数 | 实际 brems / ¹²C γ / α 数 |
| --- | --- | ---: | --- |
| 20260624_16h28m34s | 10⁶ / 10 / 10⁶ | 19 | 10 / 0 / 9 |
| 20260624_16h28m48s | 10⁷ / 10² / 10⁷ | 206 | 102 / 0 / 104 |
| 20260624_16h29m17s | 10⁸ / 10³ / 10⁸ | 1993 | 977 / 0 / 1016 |

单窗示例不能给出可靠的重建效率或误差统计。后续主实验应使用多段独立随机流，处理窗前历史和窗末截断，评价 α 效率、漏检、误拆分/误合并、假 α 率、能量偏差/分辨率和整体能谱偏差。

每次运行记录代码提交号及未提交改动、完整宏、随机种子、事件数、源谱和注入率、时间窗、关键参数、依赖环境和输出路径。实际宏和 manifest 保存在对应提交目录，ROOT 与宏的索引记录在本地 `jupyter/Run_Catalog.ipynb`。原 `data.log` 已归档到 `jupyter/run_records/`，新代码不再生成这个重复日志；RunInfo 继续保存运行条件和实际宏内容。现有历史数据的几何、光学参数及环境仍有追溯缺口。
