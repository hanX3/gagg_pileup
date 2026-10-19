# GAGG 探测器响应与脉冲堆积模拟

本项目使用 Geant4 模拟 GAGG(Ce) 中的能量沉积、闪烁发光及光子输运，记录到达虚拟 SiPM 的光子时间序列，再用 Python 分析波形。研究对象是 p–¹¹B 相关的 α 粒子、热轫致辐射背景和 ¹²C 俘获/退激 γ，目标是研究计数率升高时的随机堆积，以及 α 数量和能谱的恢复能力。

Geant4 保留到达光子，并可通过 `G4VDigitizerModule` / `G4DigiManager` 在窗末生成 SiPM 雪崩、电流、电压和 ADC，保存到同一 ROOT 的 `SiPMEvent`。当前单粒子默认器件为滨松 **S13360-6050PE**（25°C、3 V 过压），采用数据表典型参数和近似读图 PDE；恢复时间使用同系列代理值，光耦合与读出仍需标定。[SiPM 响应说明](docs/sipm-response.md) 包含统一的原生与离线入口、逐项器件依据和验证结果。原有光子计数波形与数字化电压分别保留；已有分粒子波形来自仿真真值。当前已加入纯 α 盲重建先导基线，拟合仅输入观测波形，真值单独用于评价；能量标定与独立样本评价仍待完成。

## 项目目录

| 目录 | 内容 |
| --- | --- |
| [g4](g4/) | Geant4 工程，包含源码、宏、分析脚本和结果图。 |
| [pileup](pileup/) | 纯 α 与 α/γ 混合计数率实验：配置、可独立重跑的源码快照、运行器和数据检查。 |
| [analysis](analysis/) | 单粒子 SiPM 分析、孤立 α/γ 波形拟合和拟合窗口扫描。 |
| [tools](tools/) | 各系统独立 Python 环境、Notebook 启动及本机内核注册工具。 |
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

从项目根目录运行 `bash tools/run_gagg_python.sh -m analysis.single_particle_sipm run --alpha 4 --gamma 32` 可保存固定样本及三组随机数对照、实际源码与宏、ROOT 索引、读出交叉检查和图。默认数字化关闭，原生 SiPM 运行启用后要求 `/run/beamOn 1`。具体参数、削顶标记与模型限制见 [原生数字化](docs/sipm-response.md#native-sipm)。 本次完整改动与实际结果另见本地独立笔记 [codex_SiPM_Digitization.ipynb](jupyter/codex_SiPM_Digitization.ipynb)。

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

截至 2026-09-30，首轮四档各 100 个独立时间窗已全部完成并通过逐窗数据检查，生产 ROOT 合计约 54.27 GB。2026-10-13 另建 [S13360 纯 α 计数率先导实验](pileup/experiments/alpha_s13360_rate_fit_20261013/README.md)，保留源码、配置、宏、种子和运行清单；四档比较样本分别为 16、4、3、1 个独立窗，未使用首轮保留测试集。2026-10-14 补充 20、40、60、80 MHz 各 3 窗，并按用户指定的“采集动态范围足够”前提，对全部 36 窗进行 [完整电压统一重分析](pileup/experiments/alpha_full_signal_20261014/README.md)。当前 [盲重建基线](pileup/analysis/alpha_pileup_fit.py) 从完整含噪电压估计脉冲数、到达时间和电荷；本轮不考虑 ADC 位数、量化和削顶，旧 ADC 整窗拒绝结果只作为历史记录。同日完成 [时序 CNN 与拟合的比较](pileup/experiments/alpha_ml_comparison_20261014/README.md)：26 窗训练、10 窗验证、17 个新种子测试窗，20–100 MHz 的事件恢复率改善，误检增加；另保存近邻/幅度外推失败及高率偶然时间匹配对照。完整过程记录在本地 [独立 ML Notebook](jupyter/codex_Alpha_Pileup_ML_Comparison_20261014.ipynb)。这些小样本结果不构成最终能谱或实物性能验证。下一步仍需条件能量标定、模型失配与扩大独立样本评价。方法依据与适用条件见 [文献采用记录](docs/literature-use.md)；运行明细与研究计划在本地 Jupyter 笔记中维护。

<a id="ml-multimethod-20261014"></a>
同日扩展为[五种流程的共同测试](pileup/experiments/alpha_ml_multimethod_20261014/README.md)：新增U-Net、局部Transformer及CNN＋物理拟合，保留原CNN与盲拟合；两网络各训练2400步，22个新独立窗口覆盖八档源率，原17窗作回归。算法参数仅由验证集选择，完整电压输入和足够动态范围前提保持。不同率的改善并不一致，受控近邻、幅度外推失败及高密度偶然匹配均原样记录；详见[独立Notebook第二轮](jupyter/codex_Alpha_Pileup_ML_Comparison_20261014.ipynb#ml-multimethod-20261014)和[采用记录L21](docs/literature-use.md#ml-multimethod-20261014)。

<a id="alpha-gamma-mixed-20261015"></a>
2026-10-15 新增 [α/轫致辐射 γ](pileup/experiments/alpha_gamma_brems_fit_20261015/README.md) 与 [α/661.657 keV γ](pileup/experiments/alpha_gamma_mono662_fit_20261015/README.md) 的盲拟合：α 注入 1/10 MHz、γ 均为 10 MHz，各条件 8 个独立原生共享微胞窗口。拟合未知脉冲数、时刻、电荷和类型，保持完整含噪电压输入；另有147次空/单/双脉冲控制。以下恢复要求真值Q≥500 p.e.、时间误差≤10 ns、类型正确及电荷误差≤20%，不是能量分辨率。

| γ 来源 | α 注入率 / MHz | α 完整恢复 | γ 完整恢复 |
| --- | ---: | ---: | ---: |
| 轫致辐射谱 | 1 | 44/64（68.8%） | 2/30（6.7%） |
| 轫致辐射谱 | 10 | 270/571（47.3%） | 3/42（7.1%） |
| 661.657 keV | 1 | 49/58（84.5%） | 29/32（90.6%） |
| 661.657 keV | 10 | 282/605（46.6%） | 12/31（38.7%） |

详细的检出/分型、α样本纯度、漏检及误拆分/误合并、条件电荷误差和偶然匹配对照见本地 [新 Notebook](jupyter/codex_Alpha_Gamma_Pileup_Fitting_20261015.ipynb)。固定模板仅来自两条旧孤立示例，结果受阈下信号和模型失配影响，不能推出通用最大计数率。采用依据与待验证事项见 [L22](docs/literature-use.md#alpha-gamma-mixed-20261015)。

2026-10-15–16 实施 [独立标定与稀疏加权拟合](pileup/experiments/alpha_gamma_regularized_fit_20261015/README.md)：另生成56个标定、18个单粒子开发检查、4个混合验证及8个冻结后测试窗，加入电荷相关模板、固定相关噪声权重、稀疏选择与无L1幅度重拟合，并保留弱分量。相同新测试窗中，α源10 MHz时，α完整恢复在轫致辐射背景下由81/153增至97/153，在单能γ背景下由69/155增至91/155；α源1 MHz时两种背景均退步，γ未稳定改善。额外预测53→19、漏检60→61、误合并标志24→31，不能称为全面提高。实际结果、失败控制、纯度和条件电荷误差见 [Notebook第13节](jupyter/codex_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-regularized-20261015)及 [L24](docs/literature-use.md#alpha-gamma-regularized-20261015)；每条件仅两窗，尚无占用依赖标定、沉积能谱或实物验证。

<a id="alpha-gamma-model-selection-20261016"></a>
2026-10-16 继续[联合时刻/类型搜索与残差模型比较](pileup/experiments/alpha_gamma_model_selection_20261016/README.md)，按α、γ同等重要选择方案并另生成8个测试窗。同模板控制24/40→40/40；新原生低α率轫致辐射条件相对上一版为alpha 12/16→14/16、gamma 7/12→11/12。高率为brems：α 81/142→78/142，γ 0/6→1/6；mono662：α 85/143→86/143，γ 0/6→1/6，仍不足以可靠恢复两类。真值辅助诊断和标定覆盖检查说明候选搜索之外的问题仍需处理；完整同窗三方法比较见[Notebook第14节](jupyter/codex_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-model-selection-20261016)和[L25](docs/literature-use.md#alpha-gamma-model-selection-20261016)。

<a id="alpha-gamma-ml-20261016"></a>
2026-10-16—17 完成[α＋γ多任务CNN、U-Net与网络＋物理拟合](pileup/experiments/alpha_gamma_ml_20261016/README.md)，另生成16个冻结后测试窗。混合流程相对扩展标定拟合的α完整恢复374/641→416/641，γ35/75→28/75，额外预测41→64；未达到两类同时改善，保留现有拟合默认流程。弱γ错分仍突出，匹配模板控制也显示网络尚未掌握稳定分离。完整结果、失败分层及独立性核验见[Notebook第8节](jupyter/codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-ml-implementation-20261016)和[L29](docs/literature-use.md#alpha-gamma-ml-20261016)。

## 分析已有数据

Python 分析及 Notebook 依赖见 [requirements-analysis.txt](requirements-analysis.txt)，包括 numpy、matplotlib、uproot、awkward、pandas、scipy 和 Notebook 执行组件。[requirements-analysis.lock.txt](requirements-analysis.lock.txt) 固定已验证的依赖版本；各系统分别安装对应的二进制包，不共享虚拟环境。环境工具保存在随仓库分发的 [tools](tools/) 中；`jupyter/` 下的研究笔记、图像和运行记录仍只在本地维护，GitHub 克隆不包含这些文件。

项目位于共享文件夹时，**Debian 和 Fedora 各用本机用户目录中的环境**。目录按项目路径、发行版、架构和 Python 次版本区分，默认位于 `~/.local/share/gagg-pileup/envs/`，也支持本机 `XDG_DATA_HOME`。七本当前研究笔记统一选择 **GAGG pileup (本机环境)** 内核（`gagg-pileup`）。切换系统时，Notebook 与 ROOT 路径保持相对引用；Python、依赖、内核注册和临时连接文件留在本机。

已接入项目内核的本机，日常继续使用原有 Notebook 启动入口，无需切换到项目目录：

```bash
bash ~/start_jupyter.sh
```

本机入口先检查 GAGG 环境，自动准备依赖，并把 `gagg-pileup` 内核注册到原有 Notebook 环境，再按原设置打开 `/media/sf_works`。直接打开 `GAGG_pileup/jupyter/` 下的笔记，即可使用已保存的项目内核设置；原有其他内核仍可使用。环境可用时直接复用，首次使用、依赖清单更新或系统 Python 次版本变化时自动准备，无需提前执行 setup。首次准备需要下载依赖；系统需提供 Python 3 和 venv/ensurepip 组件。安装记录保存在 `jupyter/run_records/`，中断后可重试同一个本机启动入口。

新机器只需在本机启动脚本执行 Jupyter 之前接入一次以下调用。`notebook_dir` 指向共享目录，`venv_dir` 指向该启动脚本原来使用的本机 Jupyter 环境；按该机器实际目录调整这两个变量。此后环境检查与注册由原脚本自动完成：

```bash
gagg_environment="$notebook_dir/GAGG_pileup/tools/gagg_environment.py"
if [[ -f "$gagg_environment" ]]; then
    "${GAGG_BASE_PYTHON:-/usr/bin/python3}" "$gagg_environment" register --prefix "$venv_dir"
fi
```

2026-10-14 已接入当前 Fedora 上的 `~/start_jupyter.sh` 和 `~/software/jupyter/start_jupyter.sh`；原文件备份和实际启动验证见本地 [入口记录](jupyter/run_records/local_jupyter_entry_20261014_093600/local_server_checks.json)。旧的 `jupyter/gagg_environment.py` 入口继续转交给 `tools/`，保留现有启动方式。Debian 的本机脚本尚需切换系统后接入并实测。

从项目根目录运行 `bash tools/start_gagg_jupyter.sh` 可独立启动经典 NbClassic 界面，也适用于没有原启动脚本的新机器。存在本地 SiPM 笔记时默认打开该笔记，否则显示项目文件列表。两种入口均可传入 `--no-browser`、`--ServerApp.port=8890` 等 Jupyter 参数。查看本机环境位置或执行分析脚本（分析入口同样会自动准备环境）：

```bash
python3 tools/gagg_environment.py info
bash tools/run_gagg_python.sh -m unittest pileup.analysis.test_alpha_response -v
```

项目内备用启动器按脚本位置定位项目，允许从其他目录启动；分析脚本仍按各自说明选择工作目录。自定义解释器可设置本机 `GAGG_BASE_PYTHON`，自定义环境可设置本机 `GAGG_ENV_DIR`；自定义环境应放在本机文件系统中。本机入口保留原来的 Notebook 配置，备用 NbClassic 入口使用 [classic_notebook_config.json](tools/classic_notebook_config.json)，如存在本地 `jupyter/classic_notebook_config.json` 则优先使用它。

共享 `.venv` 曾分别基于 Python 3.14.4 和 3.13.5 创建，切换系统后出现缺少标准库的内核启动错误；它不再是 Jupyter 默认入口。2026-10-09/10 与 2026-10-13 的环境记录仍按当时条件保留，用于追溯历史结果。2026-10-13 的旧环境和笔记备份位于 `jupyter/run_records/notebook_repair_20261013_112957/`；本次跨系统环境调整、原笔记和输出备份位于 `jupyter/run_records/notebook_repair_20261014_084256/`。本次只调整运行环境与笔记执行入口，不据此推断新旧 Geant4 的物理响应一致性。 本次已在 Fedora 44 / Python 3.14.3 实际执行七本 Notebook，35项分析测试及经典页面/内核检查通过；Debian目录选择已检查，实际执行待切换系统确认。[验证结果](jupyter/run_records/notebook_repair_20261014_084256/final_validation.json) 保留本次环境和逐本执行状态。

取得对应 ROOT 数据后，从项目根目录运行：

```bash
cd g4/analysis
bash ../../tools/run_gagg_python.sh plot_pileup_waveform.py ../data/20260624_16h28m34s.root --bin-width-ns 1.0
```

脚本只画第一个 event，PNG 保存到启动命令所在目录。它使用 `uproot.source.file.MemmapSource`，只加载首窗绘图需要的光子时间和粒子标签；不加载无关的位置、能量或后续窗光子数组。高计数率下绘制全部真值分量仍需预留内存。

单粒子比较示例与分析工具说明见 [分析目录说明](g4/analysis/README.md)。

纯 α 响应核验与条件能标模块见 [pileup/analysis/alpha_response.py](pileup/analysis/alpha_response.py)，使用方法见 [pileup 说明](pileup/README.md#α-响应核验与能量标定)。首批 100 个校准/验证窗的检查、长尾模板及能标比较已完成，实际输出记录在本地 `jupyter/Alpha_Response_Calibration.ipynb`；当前电荷重建基线尚未完成与此条件能标的连接。

孤立 α/γ 波形拟合与短窗口扫描分别见 [alpha_gamma_waveform_fit.py](analysis/alpha_gamma_waveform_fit.py) 和 [alpha_gamma_fit_window_scan.py](analysis/alpha_gamma_fit_window_scan.py)。默认示例需要另行提供对应 ROOT 和冻结运行记录；相关本地路径在模块中列明。两条孤立模拟波形只能展示模型行为，不能据此推断分类效率。

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
