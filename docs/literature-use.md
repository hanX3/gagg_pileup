# 本项目文献采用记录

本文件记录 GAGG(Ce) α 能谱与堆积重建项目实际采用的文献结论、证据位置、适用条件和待验证事项。推导、计算与图像继续保存在 [研究计划](../jupyter/Research_Plan_20260921.ipynb) 和 [仿真与分析笔记](../jupyter/Note.ipynb)，在具体采用条目中链接相关章节。

## 知识库入口与使用顺序

共享知识库根目录为 `/media/sf_works/papers`，从本项目只读访问。以下链接按本项目与 `papers/` 为同级目录组织。

1. 从 [topic-index.md](../../papers/topic-index.md) 定位主题，并从 [review-index.md](../../papers/review-index.md) 查询已有审核意见。
2. 优先读取与当前问题相关的主题页：
   - [数字脉冲处理与堆积恢复](../../papers/methods/数字脉冲处理与堆积恢复.md)
   - [脉冲形状甄别与粒子鉴别](../../papers/methods/脉冲形状甄别与粒子鉴别.md)
   - [闪烁体的光产额与能量分辨率](../../papers/concepts/闪烁体的光产额与能量分辨率.md)
3. 沿主题页链接读取相关 `sources/` 精读笔记；已有记录足以回答的问题直接复用。
4. 关键公式、参数或引用证据有缺口时，按具体问题定点回查原文，保留 DOI、页码、公式、图表等定位信息。

## 采用记录

### 2026-09-30：纯 α 随机堆积第一版方法与评价口径

**研究问题**：沿用当前 p–¹¹B 有效单 α 谱，如何从总光子到达波形中估计未知的脉冲数、时间、光输出及沉积能谱，并比较四档注入率？已完成的 400 窗只证明数据准备和一致性检查完成，尚无重建性能结果。

**本次阅读范围**：从主题与审核索引进入三个优先主题页，再复用下列来源笔记的相关方法、结果、限制和用户意见；核对了 Morad 第 553 项实际人工意见。本轮没有重新阅读全文 PDF、运行论文代码或验证论文性能。旧 PLAN 第 9 节的“只读公开摘录”属于 2026-09-22 的历史阅读范围，现在补入的是共享库已有的全文精读记录，不冒充本轮原文复核。

**本项目决策**：先核验已有 α 响应和光输出标定，再做纯 α 的候选发现、非负幅度/时间联合拟合与残差增删候选。先在 10⁶ Hz 独立测试窗闭环，10⁵ Hz 作低率对照，再比较 10⁷、10⁸ Hz。算法结构是本项目待验证方案；不直接移植其他探测器的性能数字。执行顺序见 [PLAN 第 12 节](../jupyter/Research_Plan_20260921.ipynb#alpha-next-20260930)。

#### L01：Zhang 2024——模板误差与隐藏脉冲

- **来源位置**：[数字脉冲主题 B23-052](../../papers/methods/数字脉冲处理与堆积恢复.md)；[来源笔记第 3、4、6 节](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)；原文 PDF2–5，式 2–5、图 2–7；[DOI](https://doi.org/10.1016/j.nima.2024.169273)。
- **文献条件**：2 英寸 × 2 英寸 LaBr₃(Ce)、R6233 PMT、500 MHz/14 bit、32 μs 离线记录；模拟采用与拟合同一模板、417/1042/1667 keV 三离散能量、高斯噪声 σ=7 ADC 道（约 1.5 keV）、指数间隔、10⁶–10⁸ s⁻¹。实验快中子区域不超过 4×10⁶ s⁻¹。
- **结论与采用**：固定模板失配可诱发假脉冲；作者按残差尝试增峰并重拟合。本项目先查 α 形状随能量、对齐及尾部截断的变化，再决定是否需要分能量模板；增峰阈值和停止规则用独立验证集校准。幅度依赖、非负约束及光子计数的噪声模型均需按本项目定义，不称逐式复现。
- **验证限制**：论文约 10 ns 间隔及建议不超过 6×10⁷ s⁻¹来自特定模拟；实验“恢复超过 98%”借助模拟效率推算。均不能作为 GAGG 的分辨间隔、最大计数率或已测效率。式 5 修正系数坐标及增峰标准仍有缺口。
- **审核/项目状态**：来源为“AI 精读草稿，待人工审核”；本项目仅采用方法方向，尚未实现或验证。

#### L02：Adrich 2021——峰数不等于脉冲数

- **来源位置**：[数字脉冲主题 B23-044](../../papers/methods/数字脉冲处理与堆积恢复.md)；[来源笔记第 3、4、6 节](../../papers/sources/2021_Adrich_NIMA_Least-squares-fitting-algorithm-for-peak-pile-up-correction-in-gamma-ray-spectroscopy.md)；原文 PDF2–5，式 1–2、表 1、图 3–6；[DOI](https://doi.org/10.1016/j.nima.2020.164962)。
- **文献条件**：1 英寸 × 1 英寸 LaBr₃:Ce、R6231-100 PMT、有源底座、CAEN DT5730，¹³⁷Cs γ；约 10⁵ 个低率单脉冲建模板，高率实测约 80 kcps、400 万事例。当前复用笔记未列采样率；平滑宽度、噪声/协方差和 p 值阈值未完整给出。另将实测单脉冲两两叠加，保存能量与偏移作真值。
- **结论与采用**：一个观测峰也可能包含多个脉冲；单模板拟合失败后允许增加分量。项目拟合的 K 为未知数，不固定为 2，也不等于 R×T；用独立单脉冲检验假拆分，用少量双脉冲辅助定位失败。
- **验证限制**：表 1 的 99.68% 是已识别双候选的接受结果；92.3(10)% 是求和峰变化估计；8.5 ns 约 50%、17 ns 超过 90% 是两能量均通过 ±3σ 的合成测试。分母、成功定义和样本不同，不能互换。高拟合优度不等于正确恢复真实脉冲数。
- **审核/项目状态**：“AI 精读草稿，待人工审核”；采用候选数与真值数分离的原则，项目阈值、匹配和能量检验尚待实现。

#### L03：Lopatin 2012——发现候选与估计能量分开

- **来源位置**：[数字脉冲主题 B23-024](../../papers/methods/数字脉冲处理与堆积恢复.md)；[来源笔记第 3 节](../../papers/sources/2012_Lopatin_IEEE_Pileup-attenuation-for-spectroscopic-signals-using-a-sparse-reconstruction.md)；原文 PDF2–4，式 1–6、算法 1。
- **文献条件**：HPGe、¹³⁷Cs/²⁴¹Am 混合源；模型使用单位积分字典、Poisson 到达和独立高斯噪声，实际采样率、源强、计数率及噪声幅度未给出完整信息。证据主要为实测能谱的定性改善，仍有未消除的求和结构。
- **结论与采用**：论文使用较强正则化寻找起点、较弱正则化估计能量。若项目后续使用稀疏候选，拟合能量时另做无 L1 收缩的非负幅度重拟合，并检验偏差；这是项目调整，不是论文原算法。
- **验证限制**：原式 4 未约束非负，式 5 对同一项列出相互冲突的阈值条件。保留疑点，不自行修成“作者算法”。目前仅列为候选发现的后备方法，不与第一版模板法同时铺开。
- **审核/项目状态**：“AI 精读草稿，待人工审核”；未复现，未选择为当前必需模块。

#### L04：GAGG 响应——形状模板与能量标定分别核验

- **来源位置**：[PSD 主题 B23-021/B23-031](../../papers/methods/脉冲形状甄别与粒子鉴别.md)、[光产额主题 B23-031](../../papers/concepts/闪烁体的光产额与能量分辨率.md)；[Kobayashi 2012 第 3–4 节](../../papers/sources/2012_Kobayashi_NIMA_Significantly-different-pulse-shapes-for-gamma-and-alpha-rays-in-Gd3Al2Ga3O12-Ce3%2B-scintillating-crystals.md)，PDF2–3、表 2，[DOI](https://doi.org/10.1016/j.nima.2012.07.055)；[Tamagawa 2015 第 3 节](../../papers/sources/2015_Tamagawa_NIMA_Alpha-gamma-pulse-shape-discrimination-in-Gd3Al2Ga3O12-GAGG-Ce3%2B-crystal-scintillator-using-shape-indicator.md)，PDF2–3、图 2–3，[DOI](https://doi.org/10.1016/j.nima.2015.05.052)。
- **文献条件**：Kobayashi 的 5 mm GAGG 立方体/H2256 PMT、5 GHz 采样、512 波形平均、0–1000 ns 衰减拟合，比较 0.662 MeV γ 与约 5.47 MeV α；另用 1 cm 立方体/双 PMT、净 1 μs 积分测相对光量。Tamagawa 为 1×1×3 cm³ GAGG/R6233K PMT、8 bit/500 MHz、2048 ns，外部约 5.48 MeV α、γ 标能范围 0.66–1.84 MeV。未提供可移植到本项目的噪声模型或高率 pileup 性能。
- **结论与采用**：平均形状差异不是逐事件堆积效率；有限积分、表面收光及 α 淬灭影响能标，γ 积分线性也不证明 α 动能线性。旧 α 模板的首光子对齐、1 μs 截断及逐事件归一化须核验；Qα(E_dep) 需保留绝对光输出与响应展宽，不能仅由单位面积模板得到。
- **项目依据**：[Note 第 4.2 节](../jupyter/Note.ipynb#single-alpha-gamma-waveforms) 及 [模板复算记录](../jupyter/artifacts/note_single_particle_4c4ltqti/run_record.json)。现模型 α 为共享 95/351 ns、光量权重 0.30/0.70；Kobayashi 表 2 α 为 145/764 ns、0.83/0.17。这是已知模型差异，本轮不替换冻结实验参数；先验证现模型下的算法，再另建物理失配对照。
- **审核/项目状态**：两篇来源均为“AI 精读草稿，待人工审核”，状态未改变。2026-09-30 已完成第一批长尾、能量分组与条件能标检查，见下面的项目验证；真实 GAGG 标定、未知脉冲数重建与物理失配仍待验证。

#### L05：Morad 2024——把漏检和条件性能分开报告

- **来源位置**：[数字脉冲主题中的 Morad 2024 条目](../../papers/methods/数字脉冲处理与堆积恢复.md)；[来源笔记第 3–5 节](../../papers/sources/2024_Morad_TNS_Model-Based-Deep-Learning-Algorithm-for-Detection-and-Classification-at-High-Event-Rates.md)；原文 PDF2–3、7、10 表 II；[DOI](https://doi.org/10.1109/TNS.2024.3371573)。
- **文献条件**：25 mm × 25 mm CLYC、PMT 800 V、1 GS/s、17 mV 触发；实测单脉冲库和模拟库，Poisson 到达及线性叠加构造高率 n/γ 信号，模拟加 AWGN。文中“0 dB”采用其本底噪声幅度定义，不是通用信号功率定义。6 Mcps 是构造信号的测试条件，不是持续在线 DAQ 证明。
- **结论与采用**：第 553 项确认的 6 Mcps 结果是条件分类 89.1%、含漏检的总准确率 72.9%、另列虚警 14.1%。项目只借鉴分母分离：匹配成功子集的能量误差必须与全体目标 α 的恢复效率一起报告，假脉冲数/有效时间另列；不把这些数值作为 α 重建性能，也不据此启动深度学习。
- **审核状态**：**“第553项的限定采用意见已获人工确认，整篇其余内容仍待审核”**；已读取 [实际人工意见](../../papers/global-review.md#human-review-553-morad2024-confirmed)，确认范围没有扩大。论文虚警分母、训练/测试独立性、公式及其他问题仍保留。
- **项目状态**：评价原则纳入 PLAN；α 的匹配规则、信号定义、误拆分/误合并统计与置信区间尚待实现。

<a id="sipm-use-20261010"></a>
### 2026-10-10：SiPM 光子响应与数字化接口（通用方案）

**问题与决定**：用户尚未选型。采用“光子到达 → 概率触发/微胞恢复 → 单雪崩电荷与波形 → 电子学/ADC”的分层。完整公式、字段和边界探测注意事项见 [PLAN 第 13 节](../jupyter/Research_Plan_20260921.ipynb#sipm-response-20261010)。2026-10-10 用户授权实现后已建立离线模型与假设参数验证，实际器件尚未标定，冻结实验未改动。

**实现与验证（2026-10-10）**：建议的演变见 [PLAN 第13.9节](../jupyter/Research_Plan_20260921.ipynb#sipm-implementation-plan-20261010)，完成记录见 [第13.10节](../jupyter/Research_Plan_20260921.ipynb#sipm-implementation-results-20261010)。[离线实现](../pileup/analysis/sipm_response.py) 使用共享微胞队列，保留线性基线，支持二维 PDE 表、恢复及相关噪声；读出支持参数化电流或实测电压模板。软件/合成极限检查和已有单窗后处理不等同于实物验证。单候选串扰/AP分支、其 q/q0 缩放、Gamma 触发延迟和 ADC 孔径是明确标注的项目子模型，未声称直接复现某篇论文的完整代码。配置中所有数值均标为假设，概率定义和适用限制见 [说明](sipm-response.md)。

**证据范围**：只读检索共享库主题与审核索引、主题页和三篇来源笔记；定点回查下列 PDF 页。三篇均保留“AI 精读草稿，待人工审核”，未找到针对这三篇的后续人工确认。没有重新全文精读或递归阅读其引用；另核对官方 Geant4 11.4 文档与本机 11.4.2 源码。

#### L06：Acerbi 与 Gundacker 2019——PDE、单雪崩电荷与响应叠加

- **来源**：*Understanding and simulating SiPMs*，NIMA 926, 16–35；[DOI](https://doi.org/10.1016/j.nima.2018.11.118)；[笔记](../../papers/sources/2019_Acerbi_NIMA_Understanding-and-simulating-SiPMs.md)；[原文](../../papers/raw/2019_Acerbi_NIMA_Understanding-and-simulating-SiPMs.pdf)。本轮核 PDF4 式3、PDF5 式5、PDF7 恢复/电路、PDF13 式18、PDF15 式19–24；PDF15另看图像。
- **采用**：PDE 含触发概率和微胞填充因子，雪崩电荷约为过压乘电容；低占用可叠加单次雪崩响应。本项目拟用面积归一双指数核记账电荷；原式19按峰高归一，该转换是项目建模选择。
- **条件与缺口**：综述的 TOF-PET 前沿计时示例在特定低暗率下省略部分长尾等效应，不能据此忽略高率 GAGG 的恢复、基线及相关噪声。核和噪声参数须由选定器件/读出链标定；当前95/351 ns是已包含在到达时间中的闪烁参数，不作为 SiPM 常数。

- **流程补核（2026-10-10）**：本次实际查看 PDF16 图26（刊页31），核实其包含 PDE/LTE 筛选、SPTR、串扰/暗事件、单胞信号叠加、低通与阈值定时；适用范围仍是指定 TOF-PET 示例。此前的 SiPMDigitizer 是项目拟建模块名，尚无已验证的同名实现；后续新增文献已定位G4SiPM专门方法论文，见L09；流程依据与项目接口分界见 [PLAN 第13.8节](../jupyter/Research_Plan_20260921.ipynb#sipm-literature-flow-20261010)。

#### L07：Jaime Rosado 2019——恢复同时影响触发概率和电荷

- **来源**：*Modeling the Nonlinear Response of Silicon Photomultipliers*，IEEE Sensors Journal 19(24), 12031–12039；[DOI](https://doi.org/10.1109/JSEN.2019.2938018)；[笔记](../../papers/sources/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.md)；[原文](../../papers/raw/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.pdf)。本轮核 PDF1 范围/式1–3、PDF2 式4–6与种子定义；PDF2另看图像。
- **采用**：逐微胞按上次触发后的实际过压计算 PDE 和电荷；所有粒子的光子共享同一器件历史。短脉冲占用式只作极限校验，不能对长 GAGG 脉冲逐 bin 独立施加硬饱和。
- **条件与缺口**：据既有笔记，实验为 S13360-1325CS/1350CS、3×3×20 mm³ LYSO、400 ns积分、约300–2100 keV及连续 LED；采样率/全部测试温度未完整列出。平均电荷式含两种子等近似，不替代强非线性逐事件波形/方差；恢复时间、PDE参数和相关噪声不移植成项目默认值。

#### L08：Ogawara 与 Ishikawa 2016——光学输运与实测读出响应分层

- **来源**：*Signal pulse emulation for scintillation detectors using Geant4 Monte Carlo with light tracking simulation*，RSI 87, 075114；[DOI](https://doi.org/10.1063/1.4959186)；[笔记](../../papers/sources/2016_Ogawara_RSI_Signal-pulse-emulation-for-scintillation-detectors-using-Geant4-Monte-Carlo-with-light-tracking-simulation.md)；[原文](../../papers/raw/2016_Ogawara_RSI_Signal-pulse-emulation-for-scintillation-detectors-using-Geant4-Monte-Carlo-with-light-tracking-simulation.pdf)。本轮核 PDF4 图1/方法/仪器、PDF6 §III.B。
- **采用及条件**：借鉴 Geant4 到达时间与读出响应卷积分层。原装置为 R6231-100 PMT、5 GS/s、600 MHz；BC408 提供包含有限闪烁和读出链的有效响应，并非 SiPM SPE。指定晶体平均662 keV脉冲验证不能作为当前高率 pileup 性能。

**工程核查**：调研时的历史 ROOT 只有光子时间和位置，支持假设常数有效 PDE 的首版后处理。2026-10-10 随后已在主工程两棵数据树增加 photon_energy_eV，记录光子进入 SiPM 时的能量，与时间/位置同步排序；新输出可据此查光谱 PDE，历史数据仍缺此字段。构建及同种子小样本回归通过，见 [PLAN 第13.6节](../jupyter/Research_Plan_20260921.ipynb#photon-energy-20261010)。真实封装反射/返回路径也不能靠离线缩放补出。Geant4 EFFICIENCY 是边界吸收后的探测概率；Digitize 需要显式调用。依据：[官方边界文档](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/TrackingAndPhysics/physicsProcess.html#boundary-process)、[数字化文档](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/Detector/digitization.html)及本机 G4OpBoundaryProcess 的 DoAbsorption/InvokeSD、G4DigiManager 接口。加入器件响应后重新核验模板和光输出能标；旧 η=0.78348 不作 PDE。

<a id="sipm-p10-use-20261010"></a>
### 2026-10-10：新增 P10 SiPM 文献与已实现模型的对照

**范围与状态**：新增18篇题录/导读已筛查，重点回读以下4篇完整来源笔记并定点核PDF；其余仅作导航。SiPM审核检索及作者/题名补查未定位这4篇后续人工意见，仍保留AI草稿待人工审核。完整公式、证据表、代码对照与后续优先级见 [PLAN第13.11节](../jupyter/Research_Plan_20260921.ipynb#sipm-p10-20261010)。库内文件只读；本次更新文献记录，不修改器件模型、假设参数或既有仿真结果。

#### L09：Niggemann等2015——G4SiPM的逐胞时间队列

- **来源**：*G4SiPM: A novel silicon photomultiplier simulation package for Geant4*；[DOI 10.1016/j.nima.2015.01.067](https://doi.org/10.1016/j.nima.2015.01.067)；[笔记](../../papers/sources/2015_Niggemann_NIMA_G4SiPM-A-novel-silicon-photomultiplier-simulation-package-for-Geant4.md)；[原文](../../papers/raw/2015_Niggemann_NIMA_G4SiPM-A-novel-silicon-photomultiplier-simulation-package-for-Geant4.pdf)
- **核验/采用**：回查PDF1–3，另看PDF2图1/式3–9原页。直接支持光子输运后建立时间—微胞队列、加入暗事件、递归相关触发和恢复。四邻居Binomial串扰、双分量Poisson后脉冲是候选扩展；当前单子候选模型并非该论文算法的完整复现。
- **条件/限制**：作者的5%符合限于S10362-33-050、3600胞、485 nm、12 ns光脉冲和250 ns积分；过压/温度与全部参数未列全。式4的PDE参考面、式5–6的概率语义及触发恢复假设不能无条件移植。未获取/运行作者包，当前兼容性未核。

#### L10：Garutti、Klanner、Rolph、Schwandt 2021——低占用Geiger Array与单雪崩电流

- **来源**：*Simulation of the response of SiPMs; Part I: Without saturation effects*；[DOI 10.1016/j.nima.2021.165853](https://doi.org/10.1016/j.nima.2021.165853)；[笔记](../../papers/sources/2021_Garutti_NIMA_Simulation-of-the-response-of-SiPMs%3B-Part-I-Without-saturation-effects.md)；[原文](../../papers/raw/2021_Garutti_NIMA_Simulation-of-the-response-of-SiPMs%3B-Part-I-Without-saturation-effects.pdf)
- **核验/采用**：PDF2–4，PDF4原页核式1–3/AP段。采用“放电数组→单位面积单胞核→积分/电路卷积”的分层及门前历史要求。其快慢衰减核是指数之和；现有双指数差核是另一现象学选择，未替换代码。
- **条件/限制**：低占用，不自行计算PDE、绝对增益或电路推导响应；输入主放电不能与到面光子混用。AP的0/1候选也是作者明确的简化。KETEK谱对照条件与证据层级见PLAN，未把其定性比较当作当前GAGG高率验证。

#### L11：Klanner 2024——饱和、概率/电荷恢复与公共负载

- **来源**：*Simulation of the response of SiPMs Part II: With saturation effects*；[DOI 10.1016/j.nima.2023.169018](https://doi.org/10.1016/j.nima.2023.169018)；[笔记](../../papers/sources/2024_Klanner_NIMA_Simulation-of-the-response-of-SiPMs-Part-II-With-saturation-effects.md)；[原文](../../papers/raw/2024_Klanner_NIMA_Simulation-of-the-response-of-SiPMs-Part-II-With-saturation-effects.pdf)
- **核验/采用**：PDF3–6，另看PDF4式3–6与接受步骤。采用高占用必须保留恢复历史和检查负载压降的物理依据；继续保留当前x/y逐胞模拟，未改成论文的均匀单胞平均场近似。
- **条件/限制**：平均响应验证不能自动支持方差、时间分辨或高DCR性能；据笔记PDF7–11仍有亚ns/12 ns区失配及高DCR验证缺口，本轮未重查这些结果页。式5恢复基线和AP条件化存在需源码澄清之处，不照抄有歧义的式子。当前跨阻输出与ADC削顶不包含公共偏置反馈。

#### L12：Corsi等2007——单胞电荷与前端电路

- **来源**：*Modelling a silicon photomultiplier (SiPM) as a signal source for optimum front-end design*；[DOI 10.1016/j.nima.2006.10.219](https://doi.org/10.1016/j.nima.2006.10.219)；[笔记](../../papers/sources/2007_Corsi_NIMA_Modelling-a-silicon-photomultiplier-%28SiPM%29-as-a-signal-source-for-optimum-front-end-design.md)；[原文](../../papers/raw/2007_Corsi_NIMA_Modelling-a-silicon-photomultiplier-%28SiPM%29-as-a-signal-source-for-optimum-front-end-design.pdf)
- **核验/采用**：PDF2–3文本核电荷关系与前放描述，图形细节沿用完整笔记。单胞电荷约为过压乘有效电容，端口脉形取决于器件电容、淬灭和前端阻抗/带宽。优先以同读出链SPE或可验证电路约束波形。
- **条件/限制**：单胞暗脉冲电路拟合，不包含完整探测概率或高占用噪声；未移植旧器件参数、未使用有疑点的反演式。

**本次检索时的项目判断**：保留Geant4→ROOT→离线逐胞雪崩→读出的接口。优先器件/SPE/双脉冲/暗态标定，再比较独立prompt/delayed串扰、多个后脉冲时间分量与增益涨落，检查前史长度及高率负载。各项仍是待验证的模型选择，未因新增文献改写已有35项软件检查的意义。记录核查结果见 [清单](../jupyter/run_records/sipm_literature_p10_20261010.json)。

<a id="sipm-native-use-20261010"></a>
### 2026-10-10：原生 Geant4 数字化的采用记录

用户随后选择从官方 `G4VDigitizerModule` / `G4DigiManager` 实现，先检查单α和单γ电信号。现已增加 [SiPMDigitizer](../g4/src/SiPMDigitizer.cc) 和C++数值核心，在完整窗光子排序后显式调用，保存当前G4Event的digi集合并写入 `SiPMEvent`；Python离线入口保留用于参数扫描与独立读出核对。此实现采用L06/L07的PDE—恢复—电荷关系、L09的统一时序队列方向和L10的雪崩数组/电路分层。L09四邻居多候选、双分量后脉冲以及L11公共偏置反馈尚未实现；参数仍为本项目假设，文献AI待审核状态不变。

官方机制依据：[Geant4 digitization](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/Detector/digitization.html)及本机11.4.2头文件。未引入第三方G4SiPM，亦不将项目 `SiPMDigitizer` 称为Geant4提供的完整SiPM物理模型。单粒子条件、配置和验证入口见 [响应说明](sipm-response.md#native-sipm)，完整改动与实际结果见独立笔记 [codex_SiPM_Digitization.ipynb](../jupyter/codex_SiPM_Digitization.ipynb)，研究路线见 [PLAN第13.12节](../jupyter/Research_Plan_20260921.ipynb#sipm-native-20261010)。这些软件/合成核验不替代真实器件或高率实验验证。

<a id="hamamatsu-s13360-20261013"></a>
## L13：滨松 S13360-6050PE 的数据表配置（2026-10-13）

用户要求将假设 SiPM 换成一个具体滨松型号；本次选择单通道 **S13360-6050PE，25°C、3 V 过压**。最新入口为 [v2 配置](../g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json)，完整依据、对照结果和限制写在 [codex_SiPM_Digitization.ipynb 器件参数节](../jupyter/codex_SiPM_Digitization.ipynb#hamamatsu-s13360-20261013)。`datasheet_based` 表示典型参数与注明的近似，不能称为实物标定。

- **厂家原始资料**：Hamamatsu Photonics, *MPPC S13360 series*，文件名 `s13360_series_kapd1052e.pdf`，Cat.No.KAPD1052E07，Jul.2025；无DOI。[官方产品页](https://www.hamamatsu.com.cn/cn/zh-cn/product/optical-sensors/mppc/mppc_mppc-array/S13360-6050PE.html)与[官方PDF](https://www.hamamatsu.com.cn/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s13360_series_kapd1052e.pdf)。项目下载副本、SHA256和页面图保存在本地 `jupyter/run_records/codex_hamamatsu_s13360_sources/`，未写入共享文献库。核对PDF文本及2/11页型号/几何；实际查看3–5页原页图像，核参数表、PE光谱曲线KAPDB0322EB及过压曲线KAPDB0324EA。
- **采用的典型输入**：有效面6×6 mm、50 μm pitch、14400 pixel、增益1.7×10⁶、DCR 2 Mcps（0.5 p.e.阈值）、串扰3%；峰值PDE 40%指定450 nm，不直接作为GAGG全光谱效率。PDE不含串扰/AP，已含74%填充因子，不重复乘该因子。端电容1280 pF、击穿电压53±5 V、工作电压温度系数54 mV/°C和环氧折射率1.55仅作为资料记录；当前没有完整封装、电路或温度反馈模型。
- **近似读图与插值**：只配置400–650 nm、0–3 V，覆盖当前GAGG发光表。3 V光谱点按原图近似读取，说明性读图精度约±1个百分点；2 V按450 nm处0.315/0.400缩放，2 V以下线性连零。波长/过压可分离与低过压行为属项目近似，不是厂家完整二维实测表。当前α有效面光子的满恢复加权PDE为0.34943，属于项目计算。
- **恢复时间代理**：Jaime Rosado, *Performance of SiPMs in the nonlinear region*, NIMA912(2018)39–42，DOI [10.1016/j.nima.2017.10.031](https://doi.org/10.1016/j.nima.2017.10.031)。文件名 `2018_Rosado_NIMA_Performance-of-SiPMs-in-the-nonlinear-region.pdf`；[共享笔记](../../papers/sources/2018_Rosado_NIMA_Performance-of-SiPMs-in-the-nonlinear-region.md)。本轮实际查看原PDF第3页表1，仅采用同系列50 μm的1350CS/3050CS行所列29 ns作代理。该表转引[1]，本轮未重新核查[1]；它不测6050PE，亦未提供可移植温度/误差。AI精读草稿仍待人工审核，未扩大已有用户采用范围；表中25 μm行的pixel数疑点不用于当前配置。
- **剩余未知项**：AP和SPTR暂不附加，零参数表示未建模；低恢复过压PDE、噪声触发曲线、29 ns代理、SPE/前端传递函数及真实光耦合需标定。50 Ω、2/20 ns核、采样/ADC均为项目演示读出参数。拟合恢复代理取20/29/40 ns时，同一α光子历史的离线电荷相对29 ns变化+2.47%/0/−2.76%，是单例敏感性，不是实物误差范围。

本轮运行保留11个新ROOT：500 mV器件诊断6窗、1000 mV量程复跑5窗；每个ROOT仍仅1个20 μs窗。v1的α峰值508.018 mV、4点削顶，v2无削顶；同种子量程对照的光子、雪崩和模拟波形一致。源能/沉积能量仍分别记录，γ展示沿用已知有沉积种子，不能用于效率估计。实现内部一致性与真实器件性能验证分开记录。

### 2026-10-13：当前入口统一为 S13360-6050PE

按用户要求，当前配置目录移除40000 pixel假设器件，笔记主体、原生与离线示例、C++测试配置和器件分析统一使用S13360-6050PE。器件数值沿用已保存的v2配置；不新增未经标定的AP/SPTR或封装参数。原配置、旧笔记与图像归档，历史ROOT、RunInfo和冻结源码不改写。

本次通过 `paper-search` 只读复用共享库索引、Rosado2018完整精读笔记及本节已有选型记录；29 ns继续标为1350CS/3050CS代理，AI草稿待人工审核状态保留，未重新精读PDF或扩展引用证据。新图与分析只读取S13360-6050PE的两批记录；旧器件对照退出当前入口。变更与验证见本地 [统一记录](../jupyter/run_records/sipm_device_only_20261013_143939/validation.json)。

<a id="psd-fit-use-20261013"></a>
### 2026-10-13：α / γ 波形拟合函数的文献核对

用户希望用波形拟合作GAGG + S13360-6050PE粒子鉴别。本轮按 `paper-search` 只读检索本地库，复用来源笔记并定点核原PDF；完整函数、参数表、书目信息/DOI、原文链接、条件和推导已集中迁入 [独立拟合笔记第2节](../jupyter/codex_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#psd-fit-functions-20261013)，核验清单及PDF哈希在 [文献核验记录](../jupyter/run_records/psd_fit_literature_20261013_163822/sources.json)。文献阶段仅核对函数；随后按用户要求完成现有两条独立模拟波形的拟合与模板判别示例，见下方项目计算更新。未重新运行Geant4或更改物理输入。

#### L14：GAGG双指数——分别保留α和γ的衰减常数

- **来源与证据**：延续L04的 [Kobayashi2012笔记](../../papers/sources/2012_Kobayashi_NIMA_Significantly-different-pulse-shapes-for-gamma-and-alpha-rays-in-Gd3Al2Ga3O12-Ce3%2B-scintillating-crystals.md)（DOI 10.1016/j.nima.2012.07.055）与 [Tamagawa2015笔记](../../papers/sources/2015_Tamagawa_NIMA_Alpha-gamma-pulse-shape-discrimination-in-Gd3Al2Ga3O12-GAGG-Ce3%2B-crystal-scintillator-using-shape-indicator.md)（DOI 10.1016/j.nima.2015.05.052）；本轮核两篇PDF2方法及PDF3原图，分别为图1/表2和图3(b)，DOI核各PDF1。
- **可用结论**：衰减段可用 `B + Af exp(-u/tauf) + As exp(-u/taus)`。Kobayashi γ为95.4/351 ns、α145/764 ns；Tamagawa γ91.2/367 ns、α138/599 ns。幅度比例不等于面积比例，后者须乘各自寿命；基线不计入光量。纯衰减段不能在两幅度均自由时独立辨识起点。
- **适用条件与缺口**：前者5 mm立方/H2256 PMT、512平均、1000 ns拟合；后者1×1×3 cm³/R6233K PMT、500 MHz/8 bit、2048 ns；约5.47–5.48 MeV外α与662 keV γ。两者均非S13360读出，无可移植的参数误差/噪声协方差。Tamagawa逐事件鉴别实际用SI，不是逐事件自由双指数拟合。
- **项目采用状态**：作为候选尾部函数和参数初值/失配对照。当前α/γ共享95/351 ns、α快光量0.30的模型仍保留；这一输入与上述实验参数组有已知差别。两篇AI草稿待人工审核，本轮没有新增人工认可或实物验证。

#### L15：SiPM后电压——发光与器件/前端响应卷积

- **来源与证据**：[Wang2019完整笔记](../../papers/sources/2019_Wang_JINST_Improvement-on-pulse-shape-discrimination-performance-of-SiPM-array-coupled-Cs2LiYCl6Ce3%2B-detector.md)，DOI 10.1088/1748-0221/14/11/P11005，本轮实看PDF6式3.2；以及L10的 [Garutti2021完整笔记](../../papers/sources/2021_Garutti_NIMA_Simulation-of-the-response-of-SiPMs%3B-Part-I-Without-saturation-effects.md)，DOI 10.1016/j.nima.2021.165853，本轮实看PDF4式1–3及§2.11。
- **条件与用途**：Wang为CLYC/SensL 64元阵列/OPA657、500 MS/s、12 bit、20 μs，热中子约3.2 MeVee及宇宙线/γ代理；Garutti为低占据SiPM模型，按单雪崩核叠加并与电子响应卷积。它们提供响应分层依据，不提供本项目GAGG α/γ实测常数。
- **项目推导**：双指数发光经单RC核后的每项为 `[exp(-u/tau)-exp(-u/taue)]/(tau-taue)`，按面积归一；完整候选保留光学、SPE、前端及采样孔径。`taue`不等于29 ns恢复代理。当前2/20 ns电流核及2 ns低通应显式卷积或用校准电压模板，不能直接冒充单RC实物模型。
- **边界与状态**：Wang原印R0C0约10.50 μs与“far less than 1 μs”措辞冲突，保留原文疑点且不迁移该数值。有限微胞/恢复非线性可能改变形状，ADC无削顶不代表SiPM线性。两篇AI草稿待人工审核；函数归一、卷积代数和现有模拟样例拟合已检查，未完成器件标定或完整非线性响应反演。

#### L16：两类模板分别拟合——以拟合优度差作PSD候选

- **来源与证据**：S. Marrone等/n_TOF Collaboration，*Pulse shape analysis of liquid scintillators for neutron studies*，NIMA490(2002)299–307；[完整笔记](../../papers/sources/2002_Marrone_NIMA_Pulse-shape-analysis-of-liquid-scintillators-for-neutron-studies.md) / [PDF](../../papers/raw/2002_Marrone_NIMA_Pulse-shape-analysis-of-liquid-scintillators-for-neutron-studies.pdf)。本轮核PDF3/4的式1/2和参数固定流程、PDF6χ²差及PDF7图5；实看3/6/7。DOI原PDF未印，本轮未核实，不由PII猜补。
- **原方法与条件**：NE213/XP2041和C6D6/EMI9823QKA，1 GS/s、8 bit、250 MHz；60Co/AmBe、约200 keVee硬阈值。先标定形状并约束参数，再对同一事件用两类固定模板分别拟合幅度、比较归一χ²。属于液闪n/γ先例，不转用时间常数、阈值或性能数字。
- **项目候选**：在同一采样窗、相同自由参数/噪声规则下，α/γ各拟合面积、观测起点及基线，定义 `D = chi2_gamma - chi2_alpha`；正值倾向α。阈值由独立验证样本确定；协方差随假设变化时需采用完整似然。当前展示的一对波形、首光子真值对齐和不同光量不足以证明分类性能。
- **原文疑点与状态**：式1慢项符号与式2正慢尾不一致，PDF7图5轴差值方向与PDF6正文相反；本项目使用明确的新定义，不宣称逐式复现。来源AI草稿待人工审核，无本函数新增人工采用意见。后续先标定孤立脉冲，再按重建光量分箱报告γ接受率、α误留率及ROC。

<a id="psd-fit-results-20261013"></a>
### 项目计算更新：2026-10-13 孤立 α / γ 波形拟合

以下是当前模拟的计算结果，不是上述文献或真实器件的性能结论。过程、拟合式、静态参数表、残差图和执行输出在 [codex_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb](../jupyter/codex_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#fit-inputs)；可复用实现为 [alpha_gamma_waveform_fit.py](../analysis/alpha_gamma_waveform_fit.py)，结果与来源信息见 [summary.json](../jupyter/artifacts/alpha_gamma_fit_20261013_165844/summary.json)。

- 采用无ADC削顶的S13360-6050PE ROOT `20261013_09h58m53s`（α）和 `09h59m09s`（γ），各一个已知孤立脉冲。源能5.486/0.661657 MeV，晶体沉积能4.939493/0.202202 MeV；γ是已选定的部分沉积有信号事件。参数从RunInfo与原运行冻结源码读取，未用当前宏推断历史条件。重复种子和仅改变电子噪声的记录不计为独立样本。
- SiPM前对有效面光子到达直方图作bin积分Poisson尾部拟合；SiPM后显式卷积2/20 ns单雪崩核、2 ns低通和0.5 ns采样孔径，拟合analog与ADC。光子诊断使用首光子对齐；电压拟合仅由观测波形估计起点/基线，不传入事件粒子标签、真值时刻或沉积能。
- 固定模拟的95/351 ns寿命，光子侧慢面积占比α/γ为0.70697/0.19464，ADC侧为0.74017/0.20474；ADC残差RMS/峰高分别0.889%/0.833%，单指数对照为2.566%/1.144%。1–2.5 μs拟合窗下固定寿命慢占比变化较小。自由寿命拟合可降低残差，但γ的ADC有效慢常数约217 ns，不能当作351 ns内禀输入的准确恢复。
- 分别用模拟输入定义的α/γ形状拟合同一ADC记录，两例的较小SSE假设都与评价标签一致。当前使用 `D_SSE = SSE_gamma - SSE_alpha`，并非已校准χ²；相邻样点残差相关系数约0.963/0.941，未建立完整协方差。两条不同光量的样例不支持鉴别率、FOM或ROC，也不能据此证明实物适用性。
- 独立状态方程与解析卷积、无噪声参数回收、32次电压优化的多初值收敛、量化与哈希/重复事件检查已通过，见 [数值验证记录](../jupyter/run_records/alpha_gamma_fit_20261013_165844/validation_analysis.json)。后续需要独立校准/验证/测试事件、共同光量区间和实物SPE/前端标定；文献审核状态保持不变。

- **短窗口追加检查**：将ADC拟合终点扩展到阈值后100、200、300、500、750 ns，并保留1–2.5 μs对照，起点仍为阈值前200 ns。固定95/351 ns时，500 ns的慢占比α/γ为0.7365/0.2012，接近2 μs的0.7402/0.2047；100–300 ns自由寿命出现大幅变化、参数触界或分量权重趋零。两类固定模板在这两个样例的全部窗口仍偏向对应标签，不构成最短窗口或鉴别效率验证。100 ns短窗拟合的α自由模型在0–2 μs共同评价区间RMS/峰高约8.48%，显示窗内残差小不保证长尾外推准确。9档窗口共72次拟合已核对样点、残差和2 μs回归，见 [拟合笔记第5.2节](../jupyter/codex_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-scan-20261013)及 [追加记录](../jupyter/run_records/alpha_gamma_short_windows_20261013_205215/validation_analysis.json)。仍需独立样本决定窗口；已有长记录定位和预触发基线需求没有消失。

<a id="alpha-rate-fitting-20261013"></a>
### 项目计算更新：2026-10-13 纯 α / S13360 计数率与单脉冲恢复

研究与全部图表见新笔记 [codex_Alpha_Pileup_Rate_Fitting.ipynb](../jupyter/codex_Alpha_Pileup_Rate_Fitting.ipynb)。本轮只读复用主题/审核索引、数字脉冲处理、PSD、光产额主题及 L01/L02/L03 来源笔记，并核对 Morad 第553项的实际意见；没有新增论文下载或原文逐页复核。实际采用 L01 的残差发现遗漏成分、L02 的脉冲数与峰数分离、L05 的条件误差与总体效率分母分离。当前逆滤波候选、非负幅度/浮点时间重拟合和停止规则是本项目方案，不称复现作者未完整给出的阈值/χ²检验；L03 的 LASSO 没有作为本次算法实现。

- **数据与模型**：新冻结实验 [alpha_s13360_rate_fit_20261013](../pileup/experiments/alpha_s13360_rate_fit_20261013/README.md)，当前 p–¹¹B 有效单 α 谱、无 γ/其他物理源、S13360-6050PE 原生波长 PDE/共享微胞恢复/串扰/暗计数及 ADC。10⁵/10⁶/10⁷/10⁸ Hz 分别16/4/3/1个独立20 μs窗，中心门[5,15) μs；另1个开发窗排除于结果。原有400窗及保留测试集未用于本次性能分析。SiPM只将计算候选保护1000万→8000万，物理/读出参数未改；实际宏、RunInfo、源码、种子、环境和哈希全部保存。
- **真值隔离与拟合函数**：拟合器只接收总电压、采样时刻、可观察ADC端点码和固定α模板，不接收注入率、真实K/时刻/能量。模板由前次独立孤立α的ADC结果固定95/351 ns、有效慢面积占比0.740174，并卷积现有读出。当前G4发光权重仍0.30/0.70，L04记录的真实GAGG模型差异未消除。实际雪崩按光子祖先归属仅用于评分，保留堆积中的微胞状态，不把独立重数字化的电荷当作本次真实贡献。
- **实际结果**：目标为中心源发射且电荷≥500满恢复雪崩等效数；一对一时间容差±10 ns，恢复还要求电荷误差≤20%。四档分别恢复10/10、35/36、142/204、0/639。10⁷ Hz 时间匹配177/204、漏检27、匹配后电荷未通过35、额外候选13，近邻误合并标志8；零误拆分标志仅限定义的±10 ns邻域。10⁸ Hz这一个中心窗全部ADC到顶，当前整窗拒绝策略使恢复数为0；未削顶模拟analog另作诊断也只通过23/639，不替代ADC性能。模拟analog峰值约2281 mV，当前ADC上限1000 mV。
- **条件误差与限制**：前三档已匹配事件的相对电荷误差中位数约+0.79%、+0.59%、+3.27%，时间误差中位数约−1.73、−0.98、−0.78 ns（相对首个有效面光子，不是源时间）。未校准沉积能量，不能称能量分辨率；电荷/光子数随高率下降含恢复及串扰，也不是厂家PDE测量。小样本探索、已知α模板、固定基线/形状、阈值选择及ADC拒绝政策均限制外推，未给出装置最大计数率或最终测试集性能。
- **验证与开发记录**：构建/SiPM CTest、25份实际RunInfo/宏和ROOT哈希、光子祖先唯一映射、电荷守恒及264个简化控制检查通过，见 [验证记录](../jupyter/run_records/alpha_rate_fit_20261013_211541/validation_analysis.json)。模板+电子噪声双脉冲在本次4–8 ns合并、12–160 ns的每点8个试验通过，不含光子/微胞涨落，不作为实物12 ns能力。最初10⁷ Hz的261–275个候选触发拟合器250保护值，已保存原结果并提高到500；这项数值保护失败不解释为物理不可分。

下一步为独立单脉冲/分能量模板与基线/相关噪声标定、更多独立窗口、增大量程或降低跨阻的新配置对照，并研究削顶区外的局部恢复和非线性电荷反演。论文来源的“AI草稿待审核”及 Morad 的限定确认范围均未改变。

## 新文献整理需求待办

下列是后续定点核验或文献整理需求；共享库的修改在 `papers` 项目中处理。本轮未使用存在歧义的公式实现算法，因此不以重新全文精读作为启动项目响应检查的前提。

| 发现日期 | 当前问题与已有来源 | 缺少的公式、参数或证据 | 需回查的原文位置／新增整理内容 | 状态 |
| --- | --- | --- | --- | --- |
| 2026-09-30 | L01，残差增峰与幅度相关模板 | 式 5 局部坐标；噪声权重、增峰显著性、复杂度及终止规则的完整实现 | Zhang PDF2–4、式 4–5及方法段；若原文未给出，后续按需查看作者代码并记录版本，不把自选参数当文献参数 | 待实现前定点核验；本轮采用方法方向 |
| 2026-09-30 | L02，模型选择与假拆分 | p 值阈值、平滑后相关噪声及拟合区间边界；合成与真实噪声差别 | Adrich PDF2–5，方法、图 6；未给出项由项目独立验证集校准 | 待核；不影响先作单脉冲失配检查 |
| 2026-09-30 | L03，后备稀疏候选方法 | 式 5 冲突条件、分组定义、实际采样和计数率缺口 | Lopatin PDF3 式 5/算法 1；需原文、勘误或作者实现证据，不自行补下标 | 后备路线启用时处理 |
| 2026-09-30 | L04，当前能区的真实 α 响应 | 连续 α 能区的光输出、长尾、位置依赖及与现模型不同的完整参数组 | 先核现有单粒子 ROOT 和训练数据；需要新的低能 α 标定文献时，在 papers 项目单列整理 | 待项目响应检查后确定缺口 |
| 2026-10-10 | L06–L08，通用 SiPM 数字化方案 | 选定型号/过压/温度后的到达光谱 PDE、增益、pixel布局与恢复、SPTR/DCR/串扰/后脉冲、同读出链 SPE 与 ADC 参数 | 已核S13360-6050PE厂家资料；还需实物恢复/SPTR/AP、低过压光谱PDE和SPE/ADC读出链标定；新精读在 papers 项目处理 | 2026-10-13已建立数据表驱动配置，见L13；实物标定待完成 |
| 2026-10-10 | L09–L12，新增P10与现有离线实现 | 串扰候选/观测概率转换、AP分量与条件化、真实SPE和恢复、公共负载参数、前史收敛 | 优先实物标定；必要时核G4SiPM式4–8、Part I表2与Part II式4–6的原实现，不把AI辨析当作者勘误 | 本轮已核方法；参数/代码复现与器件验证未完成 |
| 2026-10-13 | L14–L16，GAGG + S13360波形拟合PSD | 同晶体/温度/偏压/前端的α/γ形状与SPE标定、光量依赖、独立标签和噪声协方差 | 先用独立单脉冲校准/验证数据检查；新增实物拟合文献若需入库，在 papers 项目另行处理 | 函数文献与两条现有模拟波形拟合已完成；独立分类性能、参数不确定度及实物标定待完成 |
| 2026-10-13 | L01/L02/L05，纯α计数率恢复先导试验 | 有SiPM恢复/相关涨落时的候选阈值、局部削顶恢复和电荷→沉积能量标定；最终独立测试的效率不确定度 | 优先本项目独立校准/验证；若需采用论文完整噪声/模型选择算法，再定点核原实现，不把本次自选规则当文献参数 | 新24窗盲拟合与264个简化控制已完成；真实器件与正式性能验证待做 |


### 项目验证更新：2026-09-30 α 响应与数据接口

以下是本项目计算，不是上述文献的性能结论。过程、图表和执行输出保存在 [Alpha_Response_Calibration.ipynb](../jupyter/Alpha_Response_Calibration.ipynb)，结果及限制见 [计算摘要](../jupyter/artifacts/alpha_response_20260930_v2/summary.json)。

- 四档各按完整 ROOT 划分为 30 校准、20 验证、50 最终测试；本轮只重新读取 10⁵、10⁶ Hz 的 100 个校准/验证窗。逐文件生产 SHA256、元数据、光子及能量一致性核验通过；总波形和逐 primary 光子贡献在 1/2/5 ns 分箱相等。200 个最终测试 ROOT 未用于本轮分析。
- 以中央发射、剩余记录至少 5 μs、E_dep>0、Q_5us≥20 的条件取得 283 个校准响应与 163 个验证响应。逐事件首光子对齐后，1–5 μs 光输出占比均值 4.0327%；旧模板的源起始 1 μs 记录还包含飞行时间，不能把该均值直接当作每个旧事件的精确补偿系数。
- 线性/二次能标在低沉积能区偏差明显。增加实际冻结 α 累积光产额表乘有效比例的候选，仅在校准集拟合比例 η=0.78348。它与仿真共享物理输入，属于模型内标定；没有更改淬灭参数，也不是 SiPM PDE 测量或低能光产额的独立实验验证。
- 在三种能标共同支持的同一 157 个验证响应上，线性、二次、冻结表候选的能量 RMSE 分别为 80.57、72.31、45.28 keV。它们使用真值分离后的 Q，包含响应涨落，仅是条件能标误差，不能称为 pileup 能量分辨率。全部 163 个验证响应及各方法超范围样本另行保留；不将外推值夹到端点。
- 两个由真值筛出的独立 K=1 验证窗已用总波形做积分/模板检查。未知 K 的识别、假拆分/合并、计数率性能与保护区历史收敛尚未完成。1/2/5 ns 光子守恒也不代表重建性能对分箱已收敛。

当前代码见 [alpha_response.py](../pileup/analysis/alpha_response.py)。与本次结果配套的源码副本、划分清单、输入哈希和环境记录保存在同一结果目录。文献审核状态与项目计算状态分别维护。
