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

官方机制依据：[Geant4 digitization](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/Detector/digitization.html)及本机11.4.2头文件。未引入第三方G4SiPM，亦不将项目 `SiPMDigitizer` 称为Geant4提供的完整SiPM物理模型。单粒子条件、配置和验证入口见 [响应说明](sipm-response.md#native-sipm)，完整改动与实际结果见独立笔记 [codex_0001_SiPM_Digitization.ipynb](../jupyter/codex_0001_SiPM_Digitization.ipynb)，研究路线见 [PLAN第13.12节](../jupyter/Research_Plan_20260921.ipynb#sipm-native-20261010)。这些软件/合成核验不替代真实器件或高率实验验证。

<a id="hamamatsu-s13360-20261013"></a>
## L13：滨松 S13360-6050PE 的数据表配置（2026-10-13）

用户要求将假设 SiPM 换成一个具体滨松型号；本次选择单通道 **S13360-6050PE，25°C、3 V 过压**。最新入口为 [v2 配置](../g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json)，完整依据、对照结果和限制写在 [codex_0001_SiPM_Digitization.ipynb 器件参数节](../jupyter/codex_0001_SiPM_Digitization.ipynb#hamamatsu-s13360-20261013)。`datasheet_based` 表示典型参数与注明的近似，不能称为实物标定。

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

用户希望用波形拟合作GAGG + S13360-6050PE粒子鉴别。本轮按 `paper-search` 只读检索本地库，复用来源笔记并定点核原PDF；完整函数、参数表、书目信息/DOI、原文链接、条件和推导已集中迁入 [独立拟合笔记第2节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#psd-fit-functions-20261013)，核验清单及PDF哈希在 [文献核验记录](../jupyter/run_records/psd_fit_literature_20261013_163822/sources.json)。文献阶段仅核对函数；随后按用户要求完成现有两条独立模拟波形的拟合与模板判别示例，见下方项目计算更新。未重新运行Geant4或更改物理输入。

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

以下是当前模拟的计算结果，不是上述文献或真实器件的性能结论。过程、拟合式、静态参数表、残差图和执行输出在 [codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#fit-inputs)；可复用实现为 [alpha_gamma_waveform_fit.py](../analysis/alpha_gamma_waveform_fit.py)，结果与来源信息见 [summary.json](../jupyter/artifacts/alpha_gamma_fit_20261013_165844/summary.json)。

- 采用无ADC削顶的S13360-6050PE ROOT `20261013_09h58m53s`（α）和 `09h59m09s`（γ），各一个已知孤立脉冲。源能5.486/0.661657 MeV，晶体沉积能4.939493/0.202202 MeV；γ是已选定的部分沉积有信号事件。参数从RunInfo与原运行冻结源码读取，未用当前宏推断历史条件。重复种子和仅改变电子噪声的记录不计为独立样本。
- SiPM前对有效面光子到达直方图作bin积分Poisson尾部拟合；SiPM后显式卷积2/20 ns单雪崩核、2 ns低通和0.5 ns采样孔径，拟合analog与ADC。光子诊断使用首光子对齐；电压拟合仅由观测波形估计起点/基线，不传入事件粒子标签、真值时刻或沉积能。
- 固定模拟的95/351 ns寿命，光子侧慢面积占比α/γ为0.70697/0.19464，ADC侧为0.74017/0.20474；ADC残差RMS/峰高分别0.889%/0.833%，单指数对照为2.566%/1.144%。1–2.5 μs拟合窗下固定寿命慢占比变化较小。自由寿命拟合可降低残差，但γ的ADC有效慢常数约217 ns，不能当作351 ns内禀输入的准确恢复。
- 分别用模拟输入定义的α/γ形状拟合同一ADC记录，两例的较小SSE假设都与评价标签一致。当前使用 `D_SSE = SSE_gamma - SSE_alpha`，并非已校准χ²；相邻样点残差相关系数约0.963/0.941，未建立完整协方差。两条不同光量的样例不支持鉴别率、FOM或ROC，也不能据此证明实物适用性。
- 独立状态方程与解析卷积、无噪声参数回收、32次电压优化的多初值收敛、量化与哈希/重复事件检查已通过，见 [数值验证记录](../jupyter/run_records/alpha_gamma_fit_20261013_165844/validation_analysis.json)。后续需要独立校准/验证/测试事件、共同光量区间和实物SPE/前端标定；文献审核状态保持不变。

- **短窗口追加检查**：将ADC拟合终点扩展到阈值后100、200、300、500、750 ns，并保留1–2.5 μs对照，起点仍为阈值前200 ns。固定95/351 ns时，500 ns的慢占比α/γ为0.7365/0.2012，接近2 μs的0.7402/0.2047；100–300 ns自由寿命出现大幅变化、参数触界或分量权重趋零。两类固定模板在这两个样例的全部窗口仍偏向对应标签，不构成最短窗口或鉴别效率验证。100 ns短窗拟合的α自由模型在0–2 μs共同评价区间RMS/峰高约8.48%，显示窗内残差小不保证长尾外推准确。9档窗口共72次拟合已核对样点、残差和2 μs回归，见 [拟合笔记第5.2节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-scan-20261013)及 [追加记录](../jupyter/run_records/alpha_gamma_short_windows_20261013_205215/validation_analysis.json)。仍需独立样本决定窗口；已有长记录定位和预触发基线需求没有消失。

<a id="alpha-rate-fitting-20261013"></a>
### 项目计算更新：2026-10-13 纯 α / S13360 计数率与单脉冲恢复

研究与全部图表见新笔记 [codex_0003_Alpha_Pileup_Rate_Fitting.ipynb](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb)。本轮只读复用主题/审核索引、数字脉冲处理、PSD、光产额主题及 L01/L02/L03 来源笔记，并核对 Morad 第553项的实际意见；没有新增论文下载或原文逐页复核。实际采用 L01 的残差发现遗漏成分、L02 的脉冲数与峰数分离、L05 的条件误差与总体效率分母分离。当前逆滤波候选、非负幅度/浮点时间重拟合和停止规则是本项目方案，不称复现作者未完整给出的阈值/χ²检验；L03 的 LASSO 没有作为本次算法实现。

- **数据与模型**：新冻结实验 [alpha_s13360_rate_fit_20261013](../pileup/experiments/alpha_s13360_rate_fit_20261013/README.md)，当前 p–¹¹B 有效单 α 谱、无 γ/其他物理源、S13360-6050PE 原生波长 PDE/共享微胞恢复/串扰/暗计数及 ADC。10⁵/10⁶/10⁷/10⁸ Hz 分别16/4/3/1个独立20 μs窗，中心门[5,15) μs；另1个开发窗排除于结果。原有400窗及保留测试集未用于本次性能分析。SiPM只将计算候选保护1000万→8000万，物理/读出参数未改；实际宏、RunInfo、源码、种子、环境和哈希全部保存。
- **真值隔离与拟合函数**：拟合器只接收总电压、采样时刻、可观察ADC端点码和固定α模板，不接收注入率、真实K/时刻/能量。模板由前次独立孤立α的ADC结果固定95/351 ns、有效慢面积占比0.740174，并卷积现有读出。当前G4发光权重仍0.30/0.70，L04记录的真实GAGG模型差异未消除。实际雪崩按光子祖先归属仅用于评分，保留堆积中的微胞状态，不把独立重数字化的电荷当作本次真实贡献。
- **实际结果**：目标为中心源发射且电荷≥500满恢复雪崩等效数；一对一时间容差±10 ns，恢复还要求电荷误差≤20%。四档分别恢复10/10、35/36、142/204、0/639。10⁷ Hz 时间匹配177/204、漏检27、匹配后电荷未通过35、额外候选13，近邻误合并标志8；零误拆分标志仅限定义的±10 ns邻域。10⁸ Hz这一个中心窗全部ADC到顶，当前整窗拒绝策略使恢复数为0；未削顶模拟analog另作诊断也只通过23/639，不替代ADC性能。模拟analog峰值约2281 mV，当前ADC上限1000 mV。
- **条件误差与限制**：前三档已匹配事件的相对电荷误差中位数约+0.79%、+0.59%、+3.27%，时间误差中位数约−1.73、−0.98、−0.78 ns（相对首个有效面光子，不是源时间）。未校准沉积能量，不能称能量分辨率；电荷/光子数随高率下降含恢复及串扰，也不是厂家PDE测量。小样本探索、已知α模板、固定基线/形状、阈值选择及ADC拒绝政策均限制外推，未给出装置最大计数率或最终测试集性能。
- **验证与开发记录**：构建/SiPM CTest、25份实际RunInfo/宏和ROOT哈希、光子祖先唯一映射、电荷守恒及264个简化控制检查通过，见 [验证记录](../jupyter/run_records/alpha_rate_fit_20261013_211541/validation_analysis.json)。模板+电子噪声双脉冲在本次4–8 ns合并、12–160 ns的每点8个试验通过，不含光子/微胞涨落，不作为实物12 ns能力。最初10⁷ Hz的261–275个候选触发拟合器250保护值，已保存原结果并提高到500；这项数值保护失败不解释为物理不可分。

下一步为独立单脉冲/分能量模板与基线/相关噪声标定、更多独立窗口、增大量程或降低跨阻的新配置对照，并研究削顶区外的局部恢复和非线性电荷反演。论文来源的“AI草稿待审核”及 Morad 的限定确认范围均未改变。

<a id="literature-rate-limits-20261014"></a>
### L17：2026-10-14——拟合解堆积的计数率与恢复质量

**检索与核验范围**：本轮按 `paper-search` 只读检索本地共享库，核对既有精读笔记与审核记录，并定点查原PDF关键文字/图表。完整题录、证据表、装置/噪声条件和待核事项集中在 [纯α堆积笔记第11节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#literature-rate-limits-20261014)，不重复维护完整报告。以下来源仍为AI草稿待人工审核；`global-review.md` B23归并未新增人工确认。本轮不联网、不修改共享库，不运行论文代码或新增仿真。

- [Zhang2024](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)，[DOI](https://doi.org/10.1016/j.nima.2024.169273)：原PDF4–6、图6–7/9/11。同模板、低噪声LaBr₃模拟建议≤6×10⁷ s⁻¹，但本轮粗读图7(a)，10⁷/6×10⁷/10⁸ s⁻¹全能峰效率仅约90%/30%/8%；不是本项目逐事件电荷/时间成功率。实测快中子区≤4×10⁶ s⁻¹的>98%来自模拟推算；γ-flash约6.67×10⁷ s⁻¹由恢复间隔拟合获得，缺同率独立效率真值。
- [Khilkevitch2020](../../papers/sources/2020_Khilkevitch_NIMA_Advanced-algorithms-for-signal-processing-scintillation-gamma-ray-detectors-at-high-counting-rates.md)，[DOI](https://doi.org/10.1016/j.nima.2020.164309)：原PDF8/10/11，§4.5、表2、图17。LaBr₃+PMT实测5.1×10⁶ s⁻¹时fitting的4.44 MeV峰FWHM 0.098 MeV；模拟3×10⁷ s⁻¹时拟合补权后的总率仍偏差约33%。统计补权不恢复每个事件的时间和能量，离线处理不等于在线吞吐。
- [Löher2012](../../papers/sources/2012_Loher_NIMA_High-count-rate-gamma-ray-spectroscopy-with-LaBr3-Ce-scintillation-detectors.md)，[DOI](https://doi.org/10.1016/j.nima.2012.05.051)：原PDF5–6、图7–8/结论；模板矩阵解叠在约0.5–10 Mcps明显改善部分峰计数，最高14 Mcps仍有峰宽恶化，不能称未知脉冲数的全波形拟合或全事件恢复。
- [Shevelev2016](../../papers/sources/2016_Shevelev_NIMA_High-performance-gamma-ray-spectrometer-for-runaway-electron-studies-on-the-FT-2-tokamak.md)，[DOI](https://doi.org/10.1016/j.nima.2016.05.075)：原PDF7、表1/结论；固定形状拟合的FT-2谱仪结论至10⁷ s⁻¹，约2×10⁷ s⁻¹已观察到谱失真。

**本项目采用**：把“能看到恢复后的能谱峰”“给定阈值的逐事件恢复率”“能量偏差/分辨率”和“在线吞吐”分开定义；不移植上述LaBr₃+PMT数值作为GAGG+SiPM α上限。当前95/351 ns模板、微胞恢复和ADC量程，以及源注入率/探测事件率的差别均保留。研究主线、物理参数及拟合器不变；新增记录只说明证据与评价口径。

<a id="pulse-duration-rate-20261014"></a>
### L18：2026-10-14——用信号长度解释计数率与拟合难度

**项目计算与采用口径**：在 [纯α堆积笔记第12节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#pulse-duration-rate-20261014) 中，沿用独立单α拟合的95/351 ns模板和冻结读出，计算电压FWHM 181.5 ns、10%–90%上升时间25.81 ns、90%/99%面积时长0.728/1.536 μs；这是当前平均模拟模板的结果，不是材料或实物脉宽测量。明确采用泊松占用量 `μ=R_sig W` 和下一邻居重叠概率 `1−exp(−μ)`，另以 `1−exp(−2μ)` 表示两侧任一邻居；`R_sig`为信号事件率，与源注入率及算法恢复率分开。重叠概率不作为拟合失败率，宽度不作为死时间。

**方法条件与实际结果**：已知恰好两个脉冲及其时间、固定基线/正确模板、白噪声时，幅度标准差相对孤立脉冲放大 `1/sqrt(1−ρ²)`；20/50/100 ns间隔对应3.67/1.89/1.35倍，不是盲分辨效率。本次只对已保存盲拟合作真值评价：10⁷ Hz注入组204个阈上目标全有±1.536 μs内有电荷α邻居，仍恢复142个；最近邻<12 ns为3/23、≥50 ns为83/99。3个独立窗的事后分组受幅度、多个邻居和共享微胞状态混杂，未给出通用分辨时间或最大计数率。原漏检/假脉冲、条件电荷误差及10⁸ Hz削顶限制仍保留，未新增能量标定或G4生产。

**文献来源与迁移限制**：复用 [Zhang2024](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)，定点查原PDF2正文和PDF3图2：LaBr₃+R6233 PMT的FWHM主带约20个2 ns样点，即约40 ns，且随幅度变化。当前模板/40 ns的FWHM比约4.54，仅用于相同信号率下的占用量比较；不将文献60 Mcps按此比例换算为GAGG工作上限。来源仍为AI草稿待人工审核，本项目解析推导和数值验证不构成对文献的人工批准。

**补核原文讨论**：见 [第12.7节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#duration-literature-direct-20261014)。[Löher2012](../../papers/sources/2012_Loher_NIMA_High-count-rate-gamma-ray-spectroscopy-with-LaBr3-Ce-scintillation-detectors.md) PDF4 §3.3/图5直接扫描间隔与恢复误差，PDF5式(5)/图8使用`nτ`及60/200 ns时间尺度联系计数率；[Khilkevitch2020](../../papers/sources/2020_Khilkevitch_NIMA_Advanced-algorithms-for-signal-processing-scintillation-gamma-ray-detectors-at-high-counting-rates.md) PDF2明确PMT、ADC带宽和电缆使实测脉冲长于材料衰减；Zhang2024 PDF2/4区分峰/尾堆积并用式(6)指数间隔和图6解释近邻漏失；Shevelev2016 PDF2以旧BGO约320 ns的高率堆积说明更换快LaBr₃的动机。这些支持时间尺度与率联合评价，不提供GAGG通用上限；当前模板数值、近邻分组和ρ误差公式仍明确为本项目结果/推导。保留Löher式(5)对P(0)的文字解释与公式不一致这一限制，审核状态不变。

**验证与待办**：宽度网格/尾部收敛、50万组泊松间隔及每种条件4000次孤立/已知时间双脉冲白噪声最小二乘检查已通过；新Notebook单元已执行，原31单元及输出保留。见 [结果](../jupyter/artifacts/pulse_duration_rate_20261014_105829/summary.json) 和 [验证记录](../jupyter/run_records/pulse_duration_rate_20261014_105829/validation.json)。后续需在独立样本上分别改变慢尾/前沿，明确固定电荷或固定峰高、采样和噪声，并检验幅度比、相关噪声、微胞恢复和ADC量程；当前单一宽度数据不足以测出 `R_max(W)`。

<a id="intermediate-rates-20261014"></a>
### 项目计算更新：2026-10-14 中间源注入率补充

沿用L17/L18的率、效率与时间尺度区分，新增20/40/60/80 MHz各3个独立20 μs窗口，中心门每档共30 μs；物理与0–1000 mV ADC配置不变。实际ADC逐事件恢复为20 MHz：51/443、40 MHz：0/858、60 MHz：0/1223、80 MHz：0/1640；中心ADC端点样点比例依次0.81%、24.82%、65.80%、96.71%。任何端点触及都会按原政策整窗拒绝，故其效率下降不能单独解释为拟合的物理极限。未削顶模拟电压仅对拒绝窗作诊断，不当作实际ADC结果。三窗/档及单一脉宽仍不足以定义最大工作率或`R_max(W)`。

新环境的两份同种子检查物理分支完全一致；14份新ROOT（12独立+2对照）的实际宏/种子/哈希及电荷守恒通过。高密度光子时间+x碰撞改用精确y/能量完成评分归属，原冻结导出器保留；盲拟合不接收真值、源率或已知脉冲数。完整结果、失败分解和输入限制见[纯α笔记第13节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#intermediate-rates-20261014)及[验证记录](../jupyter/run_records/alpha_intermediate_rates_20261014_114000/validation_results.json)。本轮没有新增文献结论或修改共享库审核状态，后续优先区分量程、局部削顶恢复和固定宽度下的拟合失配。

<a id="full-signal-rates-20261014"></a>
### 项目计算更新：2026-10-14 动态范围足够的完整电压基准

**用户指定的当前假设**：采集动态范围足够，完整信号可被记录；本轮不以ADC位数、量化或削顶限制堆积恢复。统一用已有36个独立窗口的含噪 `analog_mV` 重新盲拟合，保留原SiPM状态、采样/电子学响应、固定模板和评分标准；不传入ADC端点掩码、真值或源率。旧ADC结果和“触顶即整窗拒绝”仅作历史记录，先前提出的量程扫描/局部削顶恢复不作为本轮前置工作。

**实际计算**：八档36窗均完成拟合；新增四档为20 MHz：152/443（34.31%）、40 MHz：140/858（16.32%）、60 MHz：119/1223（9.73%）、80 MHz：68/1640（4.15%）。本次20 MHz包含全部3窗，消除旧analog仅分析拒绝窗的条件选择。其余率、失败分解、条件电荷误差及固定脉宽占用量集中于[纯α笔记第14节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#full-signal-rates-20261014)。没有新G4运行；原ROOT/宏/冻结实验保留。

**验证与边界**：4项完整电压控制、36窗输入与评分检查通过，见[验证记录](../jupyter/run_records/alpha_full_signal_20261014_143000/validation_results.json)。这是假设动态范围足够且忽略量化的模拟研究，未验证某款ADC；输出电荷尚未转为沉积能谱，小样本与单一宽度不足以定义最大工作率。下一步聚焦近邻间隔、幅度比、模板失配及受控脉宽扫描。本条是项目假设与计算更新，没有新增文献结论，共享库审核状态不变。

<a id="ml-pileup-literature-20261014"></a>
### L19：2026-10-14——机器学习／人工智能的pileup候选方法

**检索与核验**：按`paper-search`只读共享库，优先主题／审核记录并全文补查，选15篇重点来源。各篇题名、作者、DOI、原PDF／笔记链接、页码、训练／测试条件和限制集中于[纯α笔记第15节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#ml-pileup-literature-20261014)及[证据表](../jupyter/run_records/ml_pileup_literature_20261014_153426/evidence.json)。复用已有完整精读，15个DOI与原PDF首页匹配；定点原页文字／图像核对范围逐篇记录，没有再次全文精读所有论文。共享库保持只读，未联网或训练模型。

**当前采用范围**：将SAEFit分量监督参数回归、1D卷积自编码器和U-Net事件序列展开作为后续候选；Morad的未知到达数／时间检出可参考，Jeon的候选峰高回归需保留检出依赖。Huang2025 Sensors的统计谱／真实率先验与TIM的时序能量分类分开；Lin的识别／计数、Fu／Liu的尾部PSD收益不等同逐事件时间＋电荷恢复。Ma2019已知4／9重数、Kafaee／Wang2022已知双脉冲只用于受控基准。没有据此替换当前算法、移植器件参数或认定GAGG高率性能。

**审核与研究条件**：第545、547、549、553、555项局部意见按原对象与范围保留；其余AI笔记仍待审核，不因本项目引用获得新认可。沿用动态范围足够的完整电压口径，后续同条件比较未知K、时间／电荷成功率、漏检、误合并／误拆分和条件误差，能量标定后另报能谱指标。网络分类99%或谱形接近不直接作为当前恢复率；局部窗容量、长尾历史和GAGG共享微胞状态都要进入训练／测试条件。

<a id="ml-pileup-comparison-20261014"></a>
### L20：2026-10-14——ML与拟合的实际模拟先导比较

**采用方向与实现**：依据L19的Regadío时序事件展开方向，自建含固定电压/导数/逆滤波特征的小型1D膨胀残差CNN，25,707参数、未知K事件概率/时间/电荷输出；不是论文U-Net/LSTM或SAEFit的完整复现。共享库只读，原人工审核范围及待审核状态不变。

**训练与独立评价**：原36个原生共享微胞SiPM窗口先按完整ROOT/种子分为26训练、10验证，随后抽切片。仿真真值仅监督/评分，不输入测试推断。新生成17个独立测试ROOT（0.1 MHz 8窗、1/10 MHz各2窗，其余五档各1窗），每窗20 μs、中央10 μs门；首轮200个保留最终测试未用。2400步CPU训练后由验证集选择第1200步/阈值0.65，保存哈希再测试，无测试调参。沿用完整含噪电压/足够动态范围前提，不使用ADC拒绝或量化限制；拟合保持原独立标定和算法参数。

**本项目计算**：同一目标按q≥500 pe、±10 ns一对一时间匹配且相对电荷误差≤20%计成功；电荷是实际共享微胞雪崩当量，不是沉积能量。

| 源注入率 / MHz | 新测试窗 | 阈上目标数 | 拟合成功 / % | CNN成功 / % |
| --- | ---: | ---: | ---: | ---: |
| 0.1 | 8 | 4 | 4 / 100.00 | 4 / 100.00 |
| 1 | 2 | 14 | 14 / 100.00 | 14 / 100.00 |
| 10 | 2 | 140 | 93 / 66.43 | 91 / 65.00 |
| 20 | 1 | 143 | 58 / 40.56 | 76 / 53.15 |
| 40 | 1 | 281 | 45 / 16.01 | 86 / 30.60 |
| 60 | 1 | 447 | 36 / 8.05 | 91 / 20.36 |
| 80 | 1 | 531 | 30 / 5.65 | 109 / 20.53 |
| 100 | 1 | 654 | 16 / 2.45 | 100 / 15.29 |

**失败与边界**：空/单/阈下和20–200 ns平均模板双脉冲控制正常，但10 ns双脉冲CNN恢复1/2（拟合2/2），40000 pe超训练范围单脉冲CNN0/1且有2个额外预测（拟合1/1）。训练电荷最大约7989 pe，未验证大幅度外推。CNN网格峰距8 ns后允许时间偏移，而拟合还合并最终<8 ns分量，约束不完全相同。高率时间移位对照保留预测数量/电荷，检验偶然匹配；其结果与漏检/误拆分/误合并、条件电荷/时间误差集中于[独立Notebook](../jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb#ml-pileup-comparison-20261014)。不能将高密度时间匹配和电荷先验贡献全部称为真实事件身份恢复。

**验证与下一步**：4项模型检查、同种子物理一致性及17窗输入/哈希/评分/仅电压网络重放通过，见[验证](../jupyter/run_records/alpha_ml_comparison_20261014_162000/validation-results.json)。训练约8.60分钟，离线单窗CPU中位耗时CNN 0.1147 s、拟合 4.723 s；不等同实时吞吐。高率一窗/档、100 MHz无同率验证、已知纯α家族与固定噪声/响应，仅为模拟先导。后续扩大独立曝光、比较神经候选+物理refit、宽度/形状/背景失配及能量标定，不能据此定义最大工作率。

<a id="ml-multimethod-20261014"></a>
### L21：2026-10-14——U-Net、Transformer与CNN＋物理拟合共同测试

**采用方向**：在L19/L20基础上新增小型1D U-Net（四级跳接，无论文LSTM瓶颈）及局部Transformer（32维、4 heads、两层，8 ns token并保留2 ns分支）；不是完整论文复现，没有实施SAEFit。复用[Regadío2021](../../papers/sources/2021_Regadio_NIMA_Unfolding-using-deep-learning-and-its-application-on-pulse-height-analysis-and-pile-up-management.md)、[DOI](https://doi.org/10.1016/j.nima.2021.165403)与[Wang2025](../../papers/sources/2025_Wang_NIMA_Transformer-algorithm-for-pile-up-correction-in-synchrotron-radiation-spectroscopic-detection-experiments.md)、[DOI](https://doi.org/10.1016/j.nima.2025.170215)，沿用已有精读的装置/采样/计数率条件与审核状态，不把其结论移植为GAGG高率能力。共享库只读。

**实现与选参**：保持原CNN和原盲拟合。新网络仍用26训练/10验证完整ROOT、相同特征/损失/随机切片及2400步，由验证集选中U-Net第800步/0.75、Transformer第2400步/0.85。两者输出1024点窗的500点中心，原CNN维持全窗推断，容量和计算时间不完全相同。混合方法使用CNN≥250 pe候选，在实际总电压上联合非负幅度及±6 ns时间调整，最终要求≥500 pe并合并<8 ns分量；`||H A−y||²+λΣ||h_j||²(A_j−A_CNN,j)²`中的λ只由验证集在0/0.1/1中选为0.1。冻结后不再调参，推断不输入真值、源率或已知K。

**本项目计算**：另生成22个新种子共同测试ROOT，低率8窗，其余各2窗，20 μs/窗，中心[5,15) μs。原17个已公开测试仅作回归、未用于选择新模型；原200个保留最终测试仍未使用。保留原生共享微胞状态和电子噪声，用完整电压、足够动态范围前提，不讨论ADC量化/削顶。统一q≥500 pe、±10 ns一对一时间匹配及电荷误差≤20%：

| 源注入率 / MHz | 窗数 | 阈上目标 | 拟合 / % | CNN / % | U-Net / % | Transformer / % | CNN＋物理 / % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 8 | 8 | 100.00 | 100.00 | 87.50 | 87.50 | 100.00 |
| 1 | 2 | 10 | 100.00 | 90.00 | 80.00 | 90.00 | 100.00 |
| 10 | 2 | 148 | 57.43 | 60.81 | 57.43 | 62.84 | 62.16 |
| 20 | 2 | 283 | 38.87 | 44.52 | 43.46 | 45.58 | 45.94 |
| 40 | 2 | 594 | 17.00 | 31.82 | 35.35 | 28.96 | 32.15 |
| 60 | 2 | 761 | 10.12 | 23.39 | 25.10 | 20.37 | 23.65 |
| 80 | 2 | 1104 | 4.89 | 16.12 | 19.38 | 13.68 | 16.30 |
| 100 | 2 | 1344 | 2.46 | 17.86 | 18.38 | 15.03 | 13.17 |

**失败、控制与限制**：九组平均模板控制中10 ns双脉冲U-Net/Transformer/拟合2/2，CNN1/2、混合0/2；40000 pe超训练范围单脉冲三网络及混合均0/1、拟合1/1。均不作为实物分辨能力。每档高率仅两窗；100 MHz无同率验证；时间移位显示高密度偶然匹配贡献。条件电荷/时间P16/50/84、漏检、额外预测、误拆分/误合并标志、固定区间示例、残差与离线耗时集中于[独立Notebook第二轮](../jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb#ml-multimethod-20261014)。较低波形RMS不保证逐事件电荷正确，混合不能补出完全没有神经候选的事件。尚无沉积能量标定、背景/形状失配和实物泛化，不给最大工作率或能谱分辨率。

**验证与待办**：5项新算法检查、同种子物理字段一致性、39窗哈希/独立ROOT与种子/仅电压精确重放/统一评分/物理电压检查通过，见[验证记录](../jupyter/run_records/alpha_ml_multimethod_20261014_180000/validation-results.json)。混合优化117次，终止成功117次、达到预算或其它非成功停止0次；终止成功不是事件恢复成功。实际源码、宏、种子、环境及共享目录清单故障恢复记录保留在[新冻结实验](../pileup/experiments/alpha_ml_multimethod_20261014/README.md)。后续扩大独立曝光并验证宽度/形状/噪声/背景失配；候选增峰或幅度范围扩展另建实验，不用本测试调参。

<a id="ml-low-rates-20261015"></a>
### 项目计算更新：2026-10-15——五方法低率0.1/1/10 MHz复查

用户指定只分析已有三档；复用L21的12个新低率ROOT，另12个公开回归窗分别列出，没有G4新运行或重新训练/调参。保持完整电压、共享微胞状态、中央[5,15) μs、q≥500 pe、±10 ns及电荷误差≤20%。真值仅用于推断后的重评分/分组/失败解释。

| 源注入率/MHz | 目标数 | 盲拟合 | CNN | U-Net | Transformer | CNN＋物理拟合 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 8 | 8/8（100.00%） | 8/8（100.00%） | 7/8（87.50%） | 7/8（87.50%） | 8/8（100.00%） |
| 1 | 10 | 10/10（100.00%） | 9/10（90.00%） | 8/10（80.00%） | 9/10（90.00%） | 10/10（100.00%） |
| 10 | 148 | 85/148（57.43%） | 90/148（60.81%） | 85/148（57.43%） | 93/148（62.84%） | 92/148（62.16%） |

**失败性质**：新0.1 MHz五方法时间匹配均8/8，U-Net/Transformer各1次电荷失败；1 MHz只有Transformer漏检1次，其余方法时间匹配10/10，CNN有1、U-Net有2次电荷失败。此前笼统的低率“漏失更多”已纠正。0.1 MHz孤立655.54 pe脉冲被两网络高估42.14%/23.43%；1 MHz618.19 pe近邻81.16 ns时CNN/U-Net高估24.65%/34.63%，Transformer峰得分0.84670未达冻结0.85；另有7622.61 pe脉冲被U-Net低估25.13%，仍在训练范围内。阈值没有按测试调整。

**与长度关联**：复用L18的当前模板FWHM181.505 ns、W90 728.024 ns。三档实测有电荷α率0.1/0.65/8.7 MHz，`R_charge W90`为0.0728/0.4732/6.3338；全部有电荷α（含保护区、阈下）查两侧近邻，W90内目标为0/8、7/10、148/148。Poisson参考概率不当失败率，分组受幅度/共享微胞混杂，不给因果分辨时间或最高工作率。0.1/1 MHz的8/8与10/10只有少量目标，逐ROOT删除范围不是置信区间，回归/新测试不合并。

**证据和边界**：[同一本Notebook第11–14节](../jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb#ml-low-rates-20261015)保存条件电荷/时间P16/50/84、全部漏检/额外预测/合并拆分、同目标配对、原生失败图及逐窗敏感性；[检查](../jupyter/artifacts/alpha_ml_low_rate_20261015_082901/validation.json)覆盖24低率窗五方法原评分、哈希和计数分解。快照及重跑见[后续实验](../pileup/experiments/alpha_ml_low_rate_20261015/README.md)。没有新增论文结论、改变Regadío待审核或Wang第549项的限定范围；共享库只读，已有文献器件/采样/噪声条件不能直接移植为GAGG低率性能。尚未沉积能量标定或实物泛化；本轮未模拟低于0.1 MHz。

<a id="alpha-gamma-mixed-20261015"></a>
### L22：2026-10-15——α/γ 原生混合波形的拟合先导

**已有依据与审核状态**：只读复用 L14–L16 的 Kobayashi 2012（GAGG/PMT 平均 α/γ 形状）、Marrone 2002（液闪/PMT 形状拟合、实测 γ 波形软件叠加）及相关 SiPM 卷积笔记，完整出处、页码与适用条件见 [既有采用记录](#psd-fit-use-20261013)。共享库的 AI 草稿待人工审核状态保持。没有文献证据可把这些装置下的性能直接视为当前 GAGG/S13360 混合逐事件准确度。

**本项目假设与实现**：两条旧孤立模拟电压固定 α/γ 慢分量权重 0.7414998/0.2047319，快/慢常数95/351 ns，加同一读出核；未知K、时刻、电荷、类型，双模板反卷积给候选，非负最小二乘及近邻离散类型/连续时刻联合搜索。算法在理想空/单/双脉冲控制中修正，正式标签只作评价。有限局部搜索不保证全局最优。另存的矩阵缩减加速试验未达到完整流程的严格数值一致性容差，没有混入正式统计。

**实际输入**：两个独立冻结实验，各自16个新正式20 μs窗，α有效p–¹¹B源注入1/10 MHz、γ注入10 MHz，中央[5,15) μs，每条件8窗。γ分别为原有轫致辐射exp(−E/100 keV)/E（10 keV–5 MeV）与661.657 keV单能源；后者仅在独立快照改写C12源槽三条能量和响应标签，不是C12俘获谱。S13360-6050PE v2，完整窗共享微胞状态，2 ns采样、0.5 ns孔径、0.1 mV电子噪声。使用完整analog电压，遵循动态范围充足前提。

**本次计算结果（不是文献性能）**：完整恢复定义为真值Q≥500 p.e.、一对一类型无关时间匹配≤10 ns、类型正确、电荷误差≤20%。γ注入率均为10 MHz。

| γ 来源 | α 注入率 / MHz | α 完整恢复 | γ 完整恢复 |
| --- | ---: | ---: | ---: |
| 轫致辐射谱 | 1 | 44/64（68.8%） | 2/30（6.7%） |
| 轫致辐射谱 | 10 | 270/571（47.3%） | 3/42（7.1%） |
| 661.657 keV | 1 | 49/58（84.5%） | 29/32（90.6%） |
| 661.657 keV | 10 | 282/605（46.6%） | 12/31（38.7%） |

同时保存检出率、条件分型率、选择纯度、条件电荷/时间误差、未匹配候选、邻接图误拆分/误合并及按独立完整窗bootstrap区间。147次理想控制的孤立18/18成功；双脉冲中的逐脉冲恢复在10 ns为32/36，20–1000 ns各36/36（每档18个双脉冲波形），严格同时受算法8 ns最小间隔限制；不含光子和微胞涨落，不能视为真实分辨时间。五个时间平移对照用于检查高密度偶然匹配，不作为p值或直接效率修正。

**失配诊断与限制**：不计入正式统计的轫致辐射预检查窗中，用真值扣除阈下γ的线性电压，保留实际共享恢复电荷、暗链与噪声，α完整恢复6/10→9/10、γ仍1/4。该单窗辅助诊断支持弱γ干扰α，但不能解释全部分型错误，也不是盲处理或降低γ率的重模拟。电荷是实际单胞等效量，未建立两类沉积能量条件标定，不能称为能量分辨率或能谱已恢复。样本量与模板标定均有限。

**可追溯记录**：[新Notebook](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb)、[轫致辐射实验](../pileup/experiments/alpha_gamma_brems_fit_20261015/README.md)、[单能γ实验](../pileup/experiments/alpha_gamma_mono662_fit_20261015/README.md)。保留36个新ROOT，两个写入失败原文件同种子重放后排除；正式分母仍32窗。全部有效文件进行各树/分支读取、宏/种子/哈希核验。未修改共享文献库。

**待验证**：扩大独立窗；建立同器件、不同沉积能量和位置/微胞占用下的独立模板与条件电荷能量标定；扫描α/γ率比；比较连续慢分量与弱背景模型。若需新文献精读或入库，仍在papers项目处理。

<a id="alpha-gamma-failure-diagnosis-20261015"></a>
**2026-10-15 失败机制补查与改进方向**：正式32窗中，α源1 MHz时，轫致辐射/单能γ的非零电荷率分别8.600/0.650 MHz，而不是相同的10 MHz探测率；轫致辐射30个过阈值γ中25个已时间匹配，其中22个被判α。其阈下γ占该类总电荷78.1%（α源10 MHz时70.0%），最终只保留大分量的模型存在明显解释缺口。新真值辅助诊断只用两条预检查波形：给定全部非零粒子的时刻与类型、保留小分量并从原电压拟合幅度，轫致辐射窗α/γ完整恢复由盲拟合6/10、1/4变为10/10、3/4。该结果使用真值，不是盲恢复成绩或严格物理上限。放松类型为同一时刻的双模板连续系数后，残差2.821→2.645 mV，α却变为5/10，支持“残差更小不保证分型正确”的本次计算结论。原来的弱γ真值扣除实验仍单独保留。

改进顺序为：分开内部建模阈值与报告阈值，联合描述弱γ并约束复杂度；使用独立孤立样本标定模板分布和信号相关噪声；再比较混合样本CNN初值＋物理联合拟合。保持500 p.e./10 ns/20%评价，同时检查纯度、漏检、假候选、误拆分/合并及电荷偏差，阈值拒判须同时报告效率。此次建议为项目待验证方案，不能从辅助真值实验承诺盲算法增益。

文献依据据既有精读：[Tamagawa等2015完整笔记](../../papers/sources/2015_Tamagawa_NIMA_Alpha-gamma-pulse-shape-discrimination-in-Gd3Al2Ga3O12-GAGG-Ce3%2B-crystal-scintillator-using-shape-indicator.md)，PDF2–4，GAGG/R6233K PMT、500 MHz/2048 ns、选定能窗的平均模板SI及低能误判边界，[出版社/DOI](https://doi.org/10.1016/j.nima.2015.05.052)；[Marrone等2002完整笔记](../../papers/sources/2002_Marrone_NIMA_Pulse-shape-analysis-of-liquid-scintillators-for-neutron-studies.md)，PDF3–4、6–7，NE213/C6D6液闪/PMT标定形状拟合及堆积，[出版社/DOI](https://doi.org/10.1016/S0168-9002(02)01063-X)。本轮核对出版社题录与DOI，未重新全文或PDF精读，AI待审核状态及已有局部采用意见保持。以上文献没有证明当前GAGG/S13360多MHz原生混合流的性能。

详细过程：[现有Notebook第10节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-failure-diagnosis-20261015)；[诊断结果](../jupyter/artifacts/alpha_gamma_fit_20261015_093500/oracle_diagnostic/preflight-onset-types.json)；[可复用诊断模块](../pileup/analysis/alpha_gamma_oracle_diagnostic.py)。原有25个Notebook单元及输出保留，新段落记录实际计算与方案边界。

<a id="alpha-gamma-weakfit-pilot-20261015"></a>
**随后实际尝试的盲改进及负结果**：新增[内部门槛与报告门槛分离实现](../pileup/analysis/alpha_gamma_weakfit.py)，弱分量保留在完整电压模型中，输出/评价仍为500 p.e./10 ns/20%。拟合只输入波形和固定标定。两条预检查窗的三次探索中，轫致辐射原始α6/10、γ1/4；内部200 p.e.变为9/10、0/4，内部100 p.e.变为5/10、1/4。单能γ原始α5/7、γ4/6；内部100 p.e.变为1/7、5/6。其γ选择纯度同时从5/5变为5/10，说明只看该类恢复增加会掩盖α泄漏。残差均降低，却未取得同时稳定改善，故不替换原基线、不扩大到32窗宣传新准确度。这些是小样本探索，不是独立性能检验。

当前证据支持下一步同时做弱分量复杂度约束、独立模板/噪声标定与分型置信度，而非单独继续降低门槛。尚未实现正则背景或混合CNN。新控制确认小分量进入模型却不进入过阈值报告，5项既有混合拟合测试与3份保存模型重构检查通过。[Notebook第11节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-weakfit-pilot-20261015)保留波形、残差、效率、纯度、假候选及运行时间；[结果索引](../jupyter/artifacts/alpha_gamma_fit_20261015_093500/weakfit_pilot/summary.json)。没有修改文献结论或共享库。


<a id="pileup-fit-optimization-literature-20261015"></a>
### L23：2026-10-15——面向混合α/γ的拟合优化文献

**问题与实际证据**：L22的三次低内部门槛盲试验降低了残差，却未稳定同时改善α/γ恢复与纯度。本次按`paper-search`只读共享库，复用完整精读笔记并定点核5篇原文；没有修改拟合器、启动新仿真或宣称准确度提高。完整题录、笔记/PDF链接、原页位置、仪器/噪声/计数率条件和优先级见[现有Notebook第12节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#pileup-fit-optimization-literature-20261015)及[证据记录](../jupyter/run_records/pileup_fit_literature_20261015/evidence.json)。

**采用的方法原则与限制**：

- [Zhang等2024](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)，[DOI](https://doi.org/10.1016/j.nima.2024.169273)，本轮PDF3–4、实际看PDF3图2–3/§3.2：幅度相关模板及残差补峰，明确模板失配可制造假脉冲。LaBr3/PMT、500 MHz及同模板合成条件保留；当前已有残差补峰，新增重点是独立模板分布与接受准则，不直接移植10 ns/60 Mcps性能。
- [Lopatin等2012](../../papers/sources/2012_Lopatin_IEEE_Pileup-attenuation-for-spectroscopic-signals-using-a-sparse-reconstruction.md)，5页会议论文，DOI未独立核；本轮PDF3原图/式(4)–(6)：较强LASSO选起点、较弱LASSO估能量，原式没有非负约束且式(5)条件冲突。借鉴两阶段原则；非负/分组候选、离散类型比较和NNLS减偏是本项目拟议改写，不是原文完整复现。HPGe实测谱缺完整采样/率参数，不据此给GAGG上限。
- [Fixsen等2004](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)，[DOI](https://doi.org/10.1016/j.nima.2003.11.313)，本轮PDF2原图/式(1)：独立标定均值与完整协方差，逆协方差加权。仅TES非线性/非平稳噪声模型验证，不是GAGG混合恢复证据；协方差随参数变化时，完整高斯似然还需log det项，协方差正则化及共享微胞相关性须本项目验证。
- [Luo等2018](../../papers/sources/2018_Luo_NIMA_Pulse-pile-up-identification-and-reconstruction-for-liquid-scintillator-based-neutron-detectors.md)，[DOI](https://doi.org/10.1016/j.nima.2018.03.078)，本轮PDF4文字：液闪n/γ四组合重建仅作联合类型模型先例；当前已实现近邻四组合。原判据为有符号差求和后取绝对值，可能抵消，不直接作为本项目似然；受控双脉冲结果不推广至随机高重数。
- [Novak与Vencelj2009](../../papers/sources/2009_Novak_TNS_Gauss-Seidel-Iterative-Method-as-a-Real-Time-Pile-Up-Solver-of-Scintillation-Pulses.md)，[DOI](https://doi.org/10.1109/TNS.2009.2030725)，本轮PDF2–3文字：已知可靠时刻、固定线性形状下用局部小系统解幅度。可借鉴缓存、局部块及暖启动；不解决漏候选/错类型，求解核吞吐不当完整恢复性能。

**项目下一步（待实现）**：独立原生孤立/空窗的模板与噪声标定 → 加权稀疏/复杂度候选选择 → 保留受约束弱分量并减小幅度收缩 → 类型/数量模型比较与置信度 → 局部求解加速。这个组合是项目方案，不是五篇论文共同验证的算法；保留500 p.e./10 ns/20%及完整电压口径，报告效率、纯度、漏失、误拆分/合并和电荷误差，能量分辨率另需标定。当前32窗已公开分析，新方法须使用独立标定/验证并冻结后另留测试；真值不进入盲推断。

**审核及核验范围**：五篇仍为AI草稿待人工审核；global-review B23-024/B23-052的原式与同模板验证限制保留，其他论文的限定认可不转移。有限网页题录核对没有补足Lopatin的可读出版社DOI证据；未下载新论文。原页文字/图像、共享来源SHA256与验证记录保存在本地jupyter。Trigano2015和Kafaee/Goodarzi2020仅查到引文/队列线索，列为待补原件，未把书目关系当已核算法依据。共享库未写入。


<a id="alpha-gamma-regularized-20261015"></a>
### L24：2026-10-15–16——独立标定、稀疏加权与条件减偏的实际实现

**采用范围**：落实L23中Zhang2024的幅度相关模板方向、Lopatin2012的选择/估计两阶段原则、Fixsen2004的独立均值/噪声标定；近邻四类型模型比较沿用已有实现。所有来源、原页与装置条件见[L23](#pileup-fit-optimization-literature-20261015)及其Notebook链接，五篇的AI草稿待人工审核状态不变。本轮没有新增文献结论或修改共享库，也没有使用Lopatin式(5)的歧义条件。本实现是项目改写，不是原论文完整复现，性能均为项目模拟计算。

**数据与模型**：另建[冻结实验](../pileup/experiments/alpha_gamma_regularized_fit_20261015/README.md)，56个原生孤立标定窗、18个单粒子开发检查窗、4个混合选参窗、8个冻结后新测试窗。保留两个父实验的G4物理、S13360-6050PE v2、2 ns采样/0.5 ns孔径/0.1 mV噪声与完整共享微胞历史；没有按测试改变物理参数。盲拟合只输入完整含噪电压和冻结参数，真值只在标定/评价中使用。保持用户指定的足够动态范围。

**由旧值到新值**：原固定α/γ慢分量0.741500/0.204732，改为在95/351 ns基函数上的独立 `f(Q)=a+b log(Q/1000 p.e.)`；α(a,b)=(0.680950,0.032393)、γ=(0.243212,0.004669)，时延分别−1.9104/−0.3973 ns。56个入射标定事件中仅16个α与13个γ达到Q≥50 p.e.，Q覆盖1545.6–16860.2和256.9–1996.7 p.e.，范围外形状夹住而电荷不夹住。独立残差给出基线方差0.0100869 mV²、信号方差系数0.614075 mV、2 ns相关系数0.899770。中心单粒子不足以标定位置/占用依赖，模板失配也进入残差，不能把全部残差方差称作本征统计噪声。

**算法与验证集选择**：标准化α/γ字典列的非负L1选候选，条件无L1 NNLS减小幅度收缩，保留弱分量，局部时间/类型与近邻四组合重拟合。候选80 p.e.，报告500 p.e.；原2 ns发现、4 ns拟合、2 ns保存模型。候选权重使用80 ns平滑电压，最终协方差由纯电压初步模型产生后固定，是可行加权最小二乘而非完整参数相关似然。λ=3/8只在4个混合验证窗选择，预声明的宏平均完整恢复F1为0.59490/0.47698，选λ=3后于2026-10-15 17:14:52 +08:00冻结，晚于此才生成测试。没有完整后验、拒判或全局最优保证；形状在每轮内固定，条件减偏不保证最终电荷严格无偏。

**同一新测试样本的结果**：α源1/10 MHz，γ源均10 MHz；中央[5,15) μs，真值Q≥500 p.e.、一对一时间匹配≤10 ns、类型正确、电荷误差≤20%。电荷为实际完全恢复单胞等效量，不是沉积能量。

| γ 来源 | α源 / MHz | α完整恢复：旧 → 新 | γ完整恢复：旧 → 新 |
| --- | ---: | --- | --- |
| 轫致辐射 | 1 | 14/16 → 11/16 | 4/14 → 3/14 |
| 轫致辐射 | 10 | 81/153 → 97/153 | 1/11 → 1/11 |
| 661.657 keV | 1 | 15/16 → 13/16 | 1/2 → 2/2 |
| 661.657 keV | 10 | 69/155 → 91/155 | 1/2 → 0/2 |

额外预测53→19、漏检60→61、误拆分标志0→0、误合并标志24→31。高率α改善伴随低率退步、γ不稳定，某些条件的电荷偏差/标准差变差；每条件仅2窗，单能γ每档仅2个目标，不能定义稳定效率或最大率。旧基线只替换带KKT检查的等价线性NNLS缩减，非线性方法不变；单窗与原直接NNLS回查的数量、类别、配对和恢复数一致，但不承诺所有局部搜索逐位一致。原32窗统计与本轮分开。

**保留的失败**：首版权重给噪声下行涨落过高权，单粒子开发检查λ=3只有α2/8、γ1/3恢复且额外13；改后8/8、3/3、额外0，与旧基线相同。开发集已用于修正，不当独立测试。40次同模板双脉冲检查中，2/4 ns网格恢复计数相同；每网格10/20/40/80/200 ns分别为0/8、2/8、8/8、8/8、6/8，含10/20 ns漏候选与200 ns γα局部错解。该控制没有光子/微胞涨落，不能宣称40 ns实物分辨能力。

**追溯与待验证**：已有Notebook[第13节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-regularized-20261015)记录全过程、纯度、完整条件误差、固定区间图和逐目标配对。86个ROOT全树/分支、宏/种子/哈希与冻结时序核验，16次测试拟合共同重评分、8个新方法模型重构、一窗只波形重放和13项检查见[验证记录](../jupyter/run_records/alpha_gamma_regularized_20261015/validation.json)。8次L1成功终止只说明满足优化器停止规则，最大投影梯度0.15813（λ=3）仍高于严格1e−5梯度容差；候选选择不是已证实的精确最优解。尚未消融各改动、完成γ位置/占用标定、分型概率校准、沉积能标或实物验证。下一步在新实验验证求解精度、受约束残差增峰、单/双模型竞争与多初值，处理控制中已见的漏候选/局部解；不得以本次已读测试再选参。

<a id="alpha-gamma-model-selection-20261016"></a>
### L25：2026-10-16——α/γ等权目标下的联合模型搜索与独立复验

**采用依据**：沿用[L23](#pileup-fit-optimization-literature-20261015)的Zhang2024残差候选/模板失配警示、Lopatin2012两阶段原则，以及Luo2018联合类别模型先例；原页、装置和采样/噪声限制不变，AI草稿待人工审核状态不变。本次只读复用已有资料，没有新增论文性能结论或修改共享库。以下均为项目自己的实现和计算；复杂度代价、KKT求解和搜索预算不是论文给定参数。

**实现及划分**：用户要求α/γ逐事件恢复同样重要。新增[独立实验](../pileup/experiments/alpha_gamma_model_selection_20261016/README.md)：保持非负L1目标，用活跃二次子问题检查1e−5投影KKT；近邻500 ns内联合两次到达时间、电荷与四类型组合；残差补峰通过固定加权SSE收益与复杂度代价比较。复用18个单粒子和4个混合开发窗，类别/条件等权完整恢复F1为上一版0.59490、增峰代价1为0.62757、代价9为0.62557，冻结代价1、λ=3。09:14:48冻结，09:15:09后生成8个新种子测试窗，每条件两窗，三方法同窗比较。未用新测试择参。

**计算结果**：已公开的同模板/仅电子白噪声控制，10/20/40/80/200 ns × 四类型组合，24/40→40/40脉冲完整恢复；冻结源码重放一致，不是物理10 ns分辨率。新原生低α率轫致辐射条件，相比上一版为alpha 12/16→14/16、gamma 7/12→11/12。高率为brems：α 81/142→78/142，γ 0/6→1/6；mono662：α 85/143→86/143，γ 0/6→1/6。高率未同时可靠恢复两类，不能称为问题已解决。

| γ来源 | α源/MHz | 类别 | 原始拟合 | 上一版稀疏 | 本轮联合 |
| --- | ---: | --- | --- | --- | --- |
| brems | 1 | alpha | 9/16（56.2%） | 12/16（75.0%） | 14/16（87.5%） |
| brems | 1 | gamma | 3/12（25.0%） | 7/12（58.3%） | 11/12（91.7%） |
| brems | 10 | alpha | 61/142（43.0%） | 81/142（57.0%） | 78/142（54.9%） |
| brems | 10 | gamma | 0/6（0.0%） | 0/6（0.0%） | 1/6（16.7%） |
| mono662 | 1 | alpha | 17/19（89.5%） | 17/19（89.5%） | 17/19（89.5%） |
| mono662 | 1 | gamma | 2/2（100.0%） | 2/2（100.0%） | 2/2（100.0%） |
| mono662 | 10 | alpha | 84/143（58.7%） | 85/143（59.4%） | 86/143（60.1%） |
| mono662 | 10 | gamma | 3/6（50.0%） | 0/6（0.0%） | 1/6（16.7%） |

仍用Q≥500 p.e.、类型无关匹配≤10 ns、类型正确且电荷误差≤20%；Q是实际单胞等效电荷，不是沉积能量。源注入率α1/10 MHz、γ10 MHz不等于阈上探测率。G4物理、共享微胞历史和读出未改变，输入完整含噪电压、动态范围足够。

**失败诊断和下一步**：本轮高α率12个阈上γ中，2个完整恢复、1个漏检、8个类型错、1个仅电荷不合格；高率α被误判为γ也从上一版9个增至26个，必须结合纯度判断。在旧开发窗刻意提供真实数量/时刻、用于形状和协方差的真实电荷及协方差用真实类型，另比较固定真实类型/双类型NNLS。高率brems即使固定真实类型也只恢复α53/73、γ1/4；该诊断不是盲算法，也不是严格物理上限。这提示当前平均模板及噪声模型可能受形状失配、涨落或重叠病态性影响，尚未做完整因果分离。开发目标78/161个α低于α标定下限1545.6 p.e.、2/12个γ高于γ标定上限1996.7 p.e.；应优先补齐低Qα、γ及位置/前史占用标定，再验证受约束形状分布模型，不能只靠进一步补峰或降低阈值。

**追溯**：现有Notebook[第14节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-model-selection-20261016)保存逐类纯度、电荷条件误差、漏检/误拆分/误合并及图像；[验证记录](../jupyter/run_records/alpha_gamma_model_selection_20261016/validation.json)核对8个ROOT全分支、132份冻结源码、24个拟合的共同评分/模型、一窗盲重放及17项算法检查。本轮线性选择达到KKT容差，但非线性搜索有有限预算，不保证全局最优；每条件两窗，未验证稳定γ效率、沉积能谱或实物泛化。

<a id="alpha-gamma-profile-fit-20261016"></a>
### L26：2026-10-16——全窗电荷重估、相关噪声与受约束类内形状

**文献复用与条件**：只读复用[Fixsen2004精读](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)的独立平均模板/协方差标定原则（原PDF2 §2–3；TES模型、非平稳噪声，不能移植其性能为GAGG结论），以及[Novak2009精读](../../papers/sources/2009_Novak_TNS_Gauss-Seidel-Iterative-Method-as-a-Real-Time-Pile-Up-Solver-of-Scintillation-Pulses.md)关于固定模板、已知时刻线性系统的条件及有限邻域/模型误差（原PDF1–3，NaI(Tl)离线例与理想模型/FPGA吞吐分开）。两份笔记仍是AI草稿待人工审核。本轮没有新增原文精读或修改共享库；来源哈希见[literature-reuse.json](../jupyter/run_records/alpha_gamma_profile_fit_20261016/literature-reuse.json)。以下算法与数值是项目实现，不是论文复现。

**实现与选择**：新增全窗非负电荷重估与时刻搜索；原AR(1)/对角/均匀及SPE卷积带状协方差比较；电压初始前沿约束；受标定RMS约束的类内慢分量变化。共12组设置，完整失败结果保留。带状协方差的经验乘数2.07796由29条独立孤立标定得到，包含模型失配，不是器件本征参数；未精确描述共享微胞相关性。类内形状均值、时刻与权重固定，两边界模板非负组合，正则化是项目经验规则，不是校准的Bayesian后验。

开发采用原4混合窗，α/γ×四条件等权完整恢复F1：原V3 0.62756657，shape_0.5为0.63009504，其余11候选更低。开发α92/161→93/161，γ5/12不变，额外预测12→11；18个孤立开发窗的阈上α8/8、γ3/3保持，无额外预测。冻结0.5倍原标定RMS、三倍尺度边界，未根据后续测试改参。原G4/SiPM和完整analog_mV保持，Q为实际单胞等效电荷；源注入率与探测率分别理解。

自由时刻开发优化48次中17次正常停止、26次达到60次迭代上限、5次ABNORMAL。负结果只限于当前初始化和预算，不代表该类方法充分收敛后的极限；选中shape_0.5后处理固定时刻。见[停止状态](../jupyter/artifacts/alpha_gamma_profile_fit_20261016/optimizer-status.json)。

**独立复验**：10:29:36冻结，10:30:40后生成8个新原生窗；两γ谱×α源1/10 MHz×每条件2窗，γ源10 MHz，Q≥500 p.e.、时间匹配≤10 ns、类型正确且电荷误差≤20%。两方法同窗、真值在拟合后评分。

本轮同一批8个新原生窗比较，受约束可变形状拟合使α完整恢复从170/290（58.62%）升到173/290（59.66%），γ仍为18/38（47.37%）。额外预测均17个，漏检49→48，误拆分标志均0，误合并标志18→17。高α源10 MHz条件合计：α148/262→151/262，γ5/16不变。改进幅度很小，且每条件只有两窗，不能据此声称稳定提高了实际效率。

保留可变形状作为实验候选，暂不替换上一版比较基线。γ正确分型数21→22，但同时满足类型和电荷要求的预测比例由18/49降为18/51；γ的完整恢复F1从0.41379降到0.40449。因此，本轮没有达到α和γ都可靠逐事件恢复的目标。12组试验也表明，在当前标定与数据条件下，更多拟合自由度、更小残差或更复杂的噪声权重不能自动转化为更好的两类恢复。

**诊断与局限**：使用祖先真值的形状描述量在开发α/γ中位数分别0.7283/0.2321，不能把所有失败归因于两类本来相同；这也不保证叠加电压可稳定分解。显式提供真数量/时刻/Q以及协方差用真实类型后，带状协方差且固定真实类型的高率brems仍仅恢复α52/73、γ2/4；是条件诊断而非盲性能或严格物理上限。下一步优先补齐α约500–1500 p.e.、γ标定上端及不同位置/前史条件的独立响应标定，建立条件化平均模板和形状散布；用含真实光子统计与共享微胞恢复的受控双脉冲扫描间隔和幅度比，检验模板失配、随机涨落与分解病态性各自的作用。当前证据尚不足以把失败归因于单一原因，也不能宣称已达到物理分辨极限。该补充标定与因果扫描本轮尚未实施。

**追溯与验证**：[冻结实验](../pileup/experiments/alpha_gamma_profile_fit_20261016/README.md)、[Notebook第15节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-profile-fit-20261016)、[验证](../jupyter/run_records/alpha_gamma_profile_fit_20261016/validation.json)。24项测试、8个新ROOT全分支、144份冻结源码、86个原生拟合模型/评分、40个控制评分与一窗只输入电压的完整重放通过。两方法在20组同模板仅电子白噪声控制均40/40，不是实测10 ns分辨率。旧52个Notebook单元及输出保留。

## 新文献整理需求待办

下列是后续定点核验或文献整理需求；共享库的修改在 `papers` 项目中处理。本轮未使用存在歧义的公式实现算法，因此不以重新全文精读作为启动项目响应检查的前提。

| 发现日期 | 当前问题与已有来源 | 缺少的公式、参数或证据 | 需回查的原文位置／新增整理内容 | 状态 |
| --- | --- | --- | --- | --- |
| 2026-10-15 | L23，两阶段/迭代稀疏重建的期刊方法 | Trigano/Gildin/Sepulcre2015（DOI线索10.1109/LSP.2015.2406911）及Kafaee/Goodarzi2020（TNS67,858–862；DOI未核）的独立原文/完整笔记 | 当前只有他文引文；后者已有共享队列键Kafaee-Goodarzi-2020-TNS67-858。需要时在papers项目核同身份原件并安排精读，不由本项目修改队列 | 书目线索，未读原文；不作为已验证性能证据 |
| 2026-10-15 | L23，模板、噪声与模型选择 | 同GAGG/SiPM的独立模板散布/协方差，弱分量复杂度、非规则增峰检验及类型不确定性 | Zhang PDF3–4、Lopatin PDF3、Fixsen PDF2已定点核；具体实现缺口与L01/L03待办合并追踪。收缩/带状协方差和减偏NNLS由项目独立验证，不冒充原文算法 | L24/L25/L26已记录独立新测试及失败；L26已试受约束类内形状，位置/占用标定与概率校准仍待完成 |
| 2026-10-14 | L19，ML／AI pileup候选 | 未知K事件输出、分量损失／尺度／排序、GAGG训练标签与独立ROOT划分；代码／权重／数据许可与真实泛化 | 优先定点核SAEFit、Regadío U-Net及Morad实现，按需要在papers项目安排新精读；本项目另建同条件空／单／双及多事件基准 | 文献检索完成；L20/L21已实施CNN、U-Net、Transformer、CNN＋物理拟合共同模拟先导，论文完整架构及实物泛化未验证 |
| 2026-09-30 | L01，残差增峰与幅度相关模板 | 式 5 局部坐标；噪声权重、增峰显著性、复杂度及终止规则的完整实现 | Zhang PDF2–4、式 4–5及方法段；若原文未给出，后续按需查看作者代码并记录版本，不把自选参数当文献参数 | 待实现前定点核验；本轮采用方法方向 |
| 2026-09-30 | L02，模型选择与假拆分 | p 值阈值、平滑后相关噪声及拟合区间边界；合成与真实噪声差别 | Adrich PDF2–5，方法、图 6；未给出项由项目独立验证集校准 | 待核；不影响先作单脉冲失配检查 |
| 2026-09-30 | L03，后备稀疏候选方法 | 式 5 冲突条件、分组定义、实际采样和计数率缺口 | Lopatin PDF3 式 5/算法 1；需原文、勘误或作者实现证据，不自行补下标 | 后备路线启用时处理 |
| 2026-09-30 | L04，当前能区的真实 α 响应 | 连续 α 能区的光输出、长尾、位置依赖及与现模型不同的完整参数组 | 先核现有单粒子 ROOT 和训练数据；需要新的低能 α 标定文献时，在 papers 项目单列整理 | 待项目响应检查后确定缺口 |
| 2026-10-10 | L06–L08，通用 SiPM 数字化方案 | 选定型号/过压/温度后的到达光谱 PDE、增益、pixel布局与恢复、SPTR/DCR/串扰/后脉冲、同读出链 SPE 与 ADC 参数 | 已核S13360-6050PE厂家资料；还需实物恢复/SPTR/AP、低过压光谱PDE和SPE/ADC读出链标定；新精读在 papers 项目处理 | 2026-10-13已建立数据表驱动配置，见L13；实物标定待完成 |
| 2026-10-10 | L09–L12，新增P10与现有离线实现 | 串扰候选/观测概率转换、AP分量与条件化、真实SPE和恢复、公共负载参数、前史收敛 | 优先实物标定；必要时核G4SiPM式4–8、Part I表2与Part II式4–6的原实现，不把AI辨析当作者勘误 | 本轮已核方法；参数/代码复现与器件验证未完成 |
| 2026-10-13 | L14–L16，GAGG + S13360波形拟合PSD | 同晶体/温度/偏压/前端的α/γ形状与SPE标定、光量依赖、独立标签和噪声协方差 | 先用独立单脉冲校准/验证数据检查；新增实物拟合文献若需入库，在 papers 项目另行处理 | 函数文献与两条现有模拟波形拟合已完成；独立分类性能、参数不确定度及实物标定待完成 |
| 2026-10-13 | L01/L02/L05，纯α计数率恢复先导试验 | 有SiPM恢复/相关涨落时的候选阈值、局部削顶恢复和电荷→沉积能量标定；最终独立测试的效率不确定度 | 优先本项目独立校准/验证；若需采用论文完整噪声/模型选择算法，再定点核原实现，不把本次自选规则当文献参数 | 新24窗盲拟合与264个简化控制已完成；真实器件与正式性能验证待做 |
| 2026-10-14 | L17，文献计数率与本项目逐事件恢复率的比较 | 统一峰效率/逐事件成功分母，同探测器条件的独立模板、噪声与饱和验证；高率独立真值效率 | 已核四篇关键原页；若采用具体算法，仍需其阈值/噪声/模型选择实现；新精读在papers项目处理 | 本轮完成条件核对；未推定GAGG+SiPM最大计数率 |
| 2026-10-14 | L18，计数率与信号长度 | 不同宽度定义、前沿/慢尾独立变化、固定电荷/峰高及相关噪声下的同指标恢复曲线 | Zhang PDF2–3宽度口径已核；优先项目独立单/双脉冲标定，若需新器件模板文献则在papers项目整理 | 当前模板、泊松占用和原样本近邻分组已分析；宽度因果扫描与最大工作率未验证 |


### 项目验证更新：2026-09-30 α 响应与数据接口

以下是本项目计算，不是上述文献的性能结论。过程、图表和执行输出保存在 [Alpha_Response_Calibration.ipynb](../jupyter/Alpha_Response_Calibration.ipynb)，结果及限制见 [计算摘要](../jupyter/artifacts/alpha_response_20260930_v2/summary.json)。

- 四档各按完整 ROOT 划分为 30 校准、20 验证、50 最终测试；本轮只重新读取 10⁵、10⁶ Hz 的 100 个校准/验证窗。逐文件生产 SHA256、元数据、光子及能量一致性核验通过；总波形和逐 primary 光子贡献在 1/2/5 ns 分箱相等。200 个最终测试 ROOT 未用于本轮分析。
- 以中央发射、剩余记录至少 5 μs、E_dep>0、Q_5us≥20 的条件取得 283 个校准响应与 163 个验证响应。逐事件首光子对齐后，1–5 μs 光输出占比均值 4.0327%；旧模板的源起始 1 μs 记录还包含飞行时间，不能把该均值直接当作每个旧事件的精确补偿系数。
- 线性/二次能标在低沉积能区偏差明显。增加实际冻结 α 累积光产额表乘有效比例的候选，仅在校准集拟合比例 η=0.78348。它与仿真共享物理输入，属于模型内标定；没有更改淬灭参数，也不是 SiPM PDE 测量或低能光产额的独立实验验证。
- 在三种能标共同支持的同一 157 个验证响应上，线性、二次、冻结表候选的能量 RMSE 分别为 80.57、72.31、45.28 keV。它们使用真值分离后的 Q，包含响应涨落，仅是条件能标误差，不能称为 pileup 能量分辨率。全部 163 个验证响应及各方法超范围样本另行保留；不将外推值夹到端点。
- 两个由真值筛出的独立 K=1 验证窗已用总波形做积分/模板检查。未知 K 的识别、假拆分/合并、计数率性能与保护区历史收敛尚未完成。1/2/5 ns 光子守恒也不代表重建性能对分箱已收敛。

当前代码见 [alpha_response.py](../pileup/analysis/alpha_response.py)。与本次结果配套的源码副本、划分清单、输入哈希和环境记录保存在同一结果目录。文献审核状态与项目计算状态分别维护。


<a id="alpha-gamma-context-calibration-20261016"></a>
### 2026-10-16：L27——低电荷、位置及前史标定的实际验证

- **复用依据与状态**：[Fixsen2004精读](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)，PDF2 §2–3、PDF3 §4（TES模拟，条件平均形状/协方差标定，[DOI](https://doi.org/10.1016/j.nima.2003.11.313)）；[Rosado2019精读](../../papers/sources/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.md)，PDF2–4恢复和非线性响应、PDF5–7两种Hamamatsu器件与LYSO/LED验证，[DOI](https://doi.org/10.1109/JSEN.2019.2938018)。两篇均为“AI精读草稿，待人工审核”，无新增人工确认。仅采用检查条件响应和恢复效应的原则，不移植TES性能或SiPM平均响应为逐事件pileup指标。
- **本项目实现**：独立新增141个标定窗，覆盖低电荷α、三点位置及四种原生双脉冲顺序；新孤立γ只有8个可用例，拟合前完整合并旧56个独立标定窗，所有零信号保留。时间表接口只改变受控源的指定方式，物理和共享微胞/读出参数不改；逐事件祖先标签仅用于标定形状及评价。旧样本集中在中央，合并后的平均不能解释成连续均匀位置的精确平均。
- **方法假设**：新增均值标定；再试电荷/可观测前史依赖的慢分量及常数+1/Q散布约束，0.5/1两尺度。回归和惩罚形式为项目方案，非论文给出的GAGG算法。H是估计起点前电压的代理量，位置/真实微胞占用不输入盲拟合；AR(1)噪声沿用旧独立标定。系数原值/新值和适用范围在[Notebook第16节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-context-calibration-20261016)及[冻结calibration](../pileup/experiments/alpha_gamma_context_calibration_20261016/calibration.json)对照。
- **实际计算**：冻结后的16个新混合窗口中，expanded相对原V3的α完整恢复为374/656（57.01%）→378/656（57.62%），γ为30/64（46.88%）→31/64（48.44%）；额外预测55→44，漏检120→117，误合并标志58→56。α完整预测比例373/563→377/555，γ30/91→31/91，γ预测仍有很大误分类/误恢复负担。各条件并不一致：高α率轫致辐射背景的γ反而4/16→3/16，单能γ背景为2/17→3/17；两种高α率条件合计γ仍仅6/33（18.18%）。按条件分层、成对重采样窗口得到恢复率差的95%百分位区间：α为−0.45至+1.65个百分点，γ为−3.13至+5.80个百分点，均跨零。该区间只包含本次窗口波动，不含标定和物理系统误差，不能据此确认稳定提升。 保留扩展标定及全部失败候选，暂不替换V3默认流程。新增标定补齐了原低电荷α范围，但均值修正尚不足以解决密集混合波形：测试后分组中，500–1000 p.e.的α只恢复38/157；与最近另一Q≥50 p.e.脉冲相距小于50 ns的γ仍为9/36，其中18例被错分型；前史电压≥200 mV的γ仅3/21，其中13例错分型。这些是真值辅助的相关性描述，不能当成因果证明或物理时间分辨极限。两个前史形状候选在开发集的α恢复均从92/161降至83/161，γ仍5/12，因此未进入本轮独立测试。每条件仅4个混合窗，原生双脉冲验证中γ目标也只有2例。 开发按两类和四条件等权F1选择，另外27孤立/16原生双脉冲验证；16新混合窗在冻结后生成。完整漏检/误检/分型/条件Q偏差与散布及窗口bootstrap见[实验记录](../pileup/experiments/alpha_gamma_context_calibration_20261016/README.md)。Q为实现p.e.，不称能量分辨率。
- **验证/待办**：26项Python测试及C++ SiPM检查通过，209个新ROOT、220份拟合复核、一个完整电压盲重跑。下一步优先补齐强前史下的γ条件标定：本次仅4个有效γ前史标定，前史电压上限约172 mV，未覆盖上述≥200 mV失败组。可从独立标定库选取多份有效γ光子输入，与前脉冲在同一微胞阵列按新数字化种子重放，系统覆盖前史、间隔和位置；这属于条件响应检验，不能拿筛选后的样本估计γ探测效率。在独立重复样本上检验条件均值和协方差是否足以描述模板涨落，再在开发集比较类型推断与电荷分配，并用新的冻结后测试判断是否改善。这些是后续方案，本轮尚未实施；本轮16个新测试窗没有用于调参。 仍需处理位置/前史覆盖不足、标定不确定性、当前共享95/351 ns基底与真实器件差异；没有把本轮结果当成物理不可辨识极限。


<a id="alpha-gamma-ml-literature-20261016"></a>
### L28：2026-10-16——α＋γ混合波形的机器学习解堆积候选

**范围与追溯**：使用本地paper-search，只读复用8篇来源笔记，6篇定点核原页；没有联网、新训练或新仿真。
完整书目、DOI、装置/采样/噪声/率条件、原页链接、项目采用范围见[新Notebook](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb)和[evidence.json](../jupyter/run_records/alpha_gamma_ml_literature_20261016/evidence.json)。
本次检索没有核实同GAGG＋共享SiPM条件的α＋γ逐事件完整恢复基准；n/γ文献作为算法类比，不能转移其数值性能。

**采用方向与边界**：Kim2023（DOI 10.1016/j.apradiso.2023.110880，PDF3–6）支持共享CNN的脉高/类型多任务，但以已找到峰为前提；12 ns双脉冲MT分类84.80%，9重降至60.22%，式3/5口径有歧义，不采作能量准确率/召回证据。
Morad2024（10.1109/TNS.2024.3371573，PDF3–7、10）支持数量/到时检测与上下文分类分开；6 Mcps的89.1%是检出后分类，含漏检总准确率72.9%，FA14.1%另列，无逐事件电荷要求。
仅该范围已有[第553项人工限定确认](../../papers/global-review.md#human-review-553-morad2024-confirmed)，其余问题仍待审核；其他7篇保持AI草稿状态。
Deltoro2025（10.1007/s41365-024-01606-y，PDF4–7）支持双解码器分离对照，限250 MHz EJ301的合成nγ/γn、12–160 ns；0.988波形相关不是98.8%事件成功，类别保持约83–85%。
Regadío2021（10.1016/j.nima.2021.165403，PDF2–3、6）支持全波形到事件序列的U-Net方向，没有α/γ分型。
Yang2017（10.1109/TNS.2017.2691729，PDF2–3、5–6）直接做BaF₂ α/γ神经网络甄别但未解堆积，0.21%是符合γ误判α比例；仅用于孤立PSD基准。
Fu2018、Peng2022本轮仅复用笔记：分别保留尾部预测在<60 ns退化、组合分类排除<100 ns的条件，不能作为当前近邻弱γ已可恢复的依据。
Tamagawa2015（10.1016/j.nima.2015.05.052，PDF3）GAGG＋PMT的α138/599 ns与当前共用95/351 ns基底不同，提示域偏移，未更改物理参数。

**项目方案，尚未实施**：优先比较带α/γ事件、时刻、电荷输出的共享1D CNN/U-Net，以及网络/残差候选并集＋允许改类的联合物理拟合；双解码器先验证四种双脉冲与空/单脉冲。固定同一候选的PSD分类作诊断对照。训练需补弱γ、近邻及强前史，保留幅度并用Q匹配对照检查能谱捷径；使用原生共享微胞数字化，按ROOT/种子/父光子身份在增强前分组。α/γ等权开发，同时报告两类完整恢复率、precision/F1、条件电荷偏差/宽度、漏检、错分及误拆分/合并，不以更小波形误差替代事件成功。模板约束保留消融，不能假定混合拟合必优于纯网络。

**现有事实与待办**：L27扩展拟合的高α源10 MHz条件γ仅6/33完整恢复，其中19错分；<50 ns组9/36恢复，≥200 mV前史组3/21，后两组相互重叠、仅作关联。训练覆盖不足与模板失配仍需检验，不能宣称物理不可分。该16窗今后作为已公开回归数据，新的ML最终测试需在网络和阈值冻结后独立评估；源注入率不替代实际探测率，Q不替代沉积能量。本轮Notebook仅汇总旧JSON并生成图表，未给出混合ML性能。共享文献源/原PDF哈希和旧Notebook保持验证见[validation.json](../jupyter/run_records/alpha_gamma_ml_literature_20261016/validation.json)。

**新增文献整理需求（在papers项目处理）**：若精确复现，定点补Kim式3/5与峰检分母，Morad FA分母/NoA先验/基础脉冲划分，Deltoro表2与正文差异/基础事件独立性；同GAGG/SiPM原生α＋γ逐事件标签及模型失配证据仍缺。本项目只记录缺口，不修改库内审核或队列。


<a id="alpha-gamma-ml-20261016"></a>
### L29：2026-10-16—17——α＋γ多任务网络与物理拟合的实际比较

**采用与边界**：落实L28的A/B路线：扩张CNN、U-Net的事件/类型/Q/时间输出，以及神经候选＋可改类的联合物理拟合。文献装置条件、原页与审核状态沿用[L28](#alpha-gamma-ml-literature-20261016)，没有新增人工确认或论文性能外推；这是项目实现，未实现双解码器/幅度单独PSD对照。

**输入和独立性**：原生共享SiPM完整20 μs电压，2 ns采样；52混合＋141孤立/前史训练窗、20混合验证窗，按ROOT/种子先划分。新加12训练窗γ源50 MHz；全部方法冻结后另生成16测试窗（两种γ×α源1/10 MHz×4窗，γ源10 MHz）。两类和四条件等权F1选参数；16旧混合只作回归，43原生控制是开发资料。物理和动态范围假设不变。

**实际结果**：完整恢复定义仍为真值Q≥500 p.e.、时间≤10 ns、类型正确及Q误差≤20%，Q不是沉积能量。

| 方法 | α完整恢复 | γ完整恢复 | α完整precision | γ完整precision | 额外预测 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 扩展标定拟合 | 374/641（58.3%） | 35/75（46.7%） | 70.9% | 37.6% | 41 |
| 多任务CNN | 342/641（53.4%） | 25/75（33.3%） | 67.3% | 41.7% | 20 |
| 多任务U-Net | 375/641（58.5%） | 21/75（28.0%） | 69.3% | 51.2% | 27 |
| U-Net＋物理拟合 | 416/641（64.9%） | 28/75（37.3%） | 69.0% | 38.9% | 64 |

**采用决定**：未达到两类恢复同时改善，不替换已有拟合默认流程。混合相对拟合α+6.55、γ−9.33个百分点；按条件分层的成对窗口bootstrap95%区间分别+2.62至+10.49、−16.88至−1.23个百分点，只含有限窗口波动。弱γ500–1000 p.e.恢复19/49→10/49、错分15/49→26/49；前史≥200 mV的γ2/13→0/13，各分层重叠。高α率γ仍4/31，未解决密集γ恢复。匹配模板双脉冲80事件拟合80/80、CNN46/80、U-Net45/80、混合71/80，提示当前学习和组合仍有缺陷，不能解释成物理极限。

**追溯、验证与待办**：[冻结实验](../pileup/experiments/alpha_gamma_ml_20261016/README.md)、[Notebook第8节](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-ml-implementation-20261016)、[评分/重放复核](../jupyter/run_records/alpha_gamma_ml_20261016/analysis-validation.json)。9项算法检查、28新有效ROOT、420份重评分、四方法电压盲重放；原生控制γ仅3个。保留序列化失败及VirtualBox受阻记录，按原种子在本机缓存完成训练后归档。下一步Q匹配及前史条件原生训练、固定候选PSD/幅度消融与先验校准，开发约束逐类恢复和误报，另用新测试；真实器件、物理域偏移和能量标定仍未验证，不修改共享文献库。

<a id="alpha-gamma-prior-ablation-20261019"></a>
### L30：2026-10-19——固定候选下拆分网络电荷与类别惩罚

**依据与范围**：沿用L28/L29中Kim2023、Morad2024的多任务及检出/分类分开评价方向和原适用条件；本轮只复用既有精读笔记与项目代码，没有重新核原PDF或新增论文性能结论。Morad第553项限定确认、其余AI草稿状态保持。将两个显式惩罚拆分是本项目的受控诊断，文献没有验证当前归一化网络分数为校准后验，也没有验证当前单个Q先验可跨类别共用。

**设计与实现**：20个已公开开发窗（两γ背景×α源1/10 MHz×5窗，γ源10 MHz），791个α/86个γ阈上目标。固定盲候选、初始类别/Q、到达时刻、模板Q及噪声权重，无到时搜索或增峰；仅比较(λQ,λC)=(0,0)/(10,0)/(0,10)/(10,10)。强度沿用旧开发选择，没有新调参。“无显式先验”仍有共同网络候选与初始形状，不能当作纯拟合。盲输入仅完整电压与保存的盲输出，真值在全部拟合后评分。

**实际计算**：无/仅电荷/仅类别/两者的α完整恢复为483/791、517/791、503/791、520/791；γ为42/86、37/86、42/86、38/86；额外预测53/83/59/83。仅电荷项使6个原成功γ失败（5个错分α、1个Q超差），补回1个；γ差−5.81个百分点，成对整窗分层重采样95%范围−10.34至−2.15个百分点。仅类别项α增加20个、γ净变化0（4失4得），γ完整precision37.2%→54.5%，但α完整precision略降、额外预测增加。区间仅描述已公开开发窗的波动，不作为新独立性能或实物结论。

**解释与采用决定**：当前把已解码事件的一个Q中心同时用于α、γ假设，Q惩罚能够影响分型。事后首个错分示例真值γ为554.14 p.e.、网络Q967.02 p.e.；无先验得584.67 p.e.的γ，仅电荷后得948.95 p.e.的α。该例用于解释全窗优化的变化，不把全部失败或物理可辨识性归因于一个因素。暂保留默认流程；下一轮优先检验完整流程中的候选/类别提示与物理Q重估，再做独立的类型/电荷/可观测前史条件Q校准与Q匹配训练，冻结后另留新测试。

**记录与验证**：[实验和复算](../pileup/experiments/alpha_gamma_prior_ablation_20261019/README.md)、[现有Notebook第9节](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-prior-ablation-20261019)、[本地结果](../jupyter/artifacts/alpha_gamma_prior_ablation_20261019/summary.json)。4项数值测试、20份密集输出解码、80组共同候选/目标单调性及预定首窗四组重放纳入验证。全部满足至多6轮的停止容差，但不是全局最优保证。没有新G4、训练、概率校准或能量标定；ROOT和旧结果保留，共享文献库未写入。

<a id="alpha-gamma-full-prior-20261019"></a>
### L31：2026-10-19——仅类别惩罚的完整混合流程及新独立测试

**依据与边界**：沿用L28–L30的既有精读与原审核状态，没有新增原文复核或论文性能结论。拆分惩罚及筛选容差是项目方案；归一化网络事件分数仍不当作已校准后验。当前Q相关模板仍由初始网络/物理电荷决定，去掉显式Q项不等于去掉全部网络Q影响。

**实现与冻结**：仅λQ 10→0、λC保留10，候选/形状初值、固定噪声权重、时刻/类别/电荷搜索及残差补峰预算沿用旧完整流程。20公开开发窗相对物理拟合α471/791→505/791、γ38/86→42/86，额外预测45→57，达到事先设置的逐类恢复、precision与额外预测筛选。工程容差不作为文献阈值或统计等效性结论。其后冻结旧U-Net权重、阈值[0.4,0.55]、标定和拟合配置，另生成16新种子原生窗，两γ背景×α源1/10 MHz×各4窗，γ源10 MHz。物理与完整电压/足够动态范围前提不变。

**独立计算结果**：同窗物理拟合→仅类别混合，α345/643→382/643，γ30/79→29/79；额外预测40→54。完整恢复仍要求Q≥500 p.e.、匹配≤10 ns、类型正确、Q误差≤20%，不是沉积能量分辨率。本批独立数据未通过与开发相同的筛选，不替换默认拟合流程。未通过项：gamma_pooled_recovery, cell_brems_1000000_alpha, cell_brems_10000000_gamma。全部五方法、分条件、错分/漏检与条件Q误差保存在本地summary。- 相对扩展标定物理拟合的alpha恢复数变化+37，恢复率差+5.75个百分点；成对整窗95%范围[+4.10, +7.72]个百分点。 - 相对扩展标定物理拟合的gamma恢复数变化-1，恢复率差-1.27个百分点；成对整窗95%范围[-6.06, +3.70]个百分点。

**解释、验证与待办**：γ损失8个原成功事件、补回7个；弱γ为扩展标定物理拟合 16/52、旧混合：Q＋类别 13/52、新混合：仅类别 13/52。窗口重采样只含有限统计，各分层重叠；若干非线性子问题达到评价预算，不声称全局最优。3项数值检查、旧端点和冻结网络密集图/五方法数值复现、16窗重建文件访问检查、ROOT/RunInfo/原宏与重评分纳入验证。首次缺输出目录在生成ROOT前修复，日志保留；分析外层143中断后仅补跑缺结果的2窗，14个完整结果保持原哈希，不新增样本或调参。下一步在独立标定/开发资料中检验类型/Q/可观测前史条件的模板与不确定度，再设计新实验；此批测试不再用于调参。未验证真实器件、沉积能量与物理域偏移，不修改共享文献库。

**追溯**：[冻结实验](../pileup/experiments/alpha_gamma_full_prior_20261019/README.md)、[Notebook第10节](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-full-prior-20261019)、[本地验证](../jupyter/run_records/alpha_gamma_full_prior_20261019/validation.json)。

**补充分型对照**：同批完整混合由无显式惩罚到仅类别提示，α367/643→382/643、γ32/79→29/79，额外预测49→54；丢失5个原成功γ且均错分α，补回2个。相对物理拟合损失的8个原成功γ也全部错分α。类别提示仍有取舍，下一步优先在训练/独立标定资料中做Q匹配PSD、类别分数校准和可观测前史/Q条件模板检查，不用本批测试选参。γ相对物理拟合的窗口区间跨零，不声称总体显著下降。

<a id="alpha-gamma-recovery-strategy-20261019"></a>
### L32：2026-10-19——针对弱γ错分的文献复核与后续方案

**检索与证据层级**：按`paper-search`只读共享库，复用主题页、完整精读笔记和global-review原始答复。7篇重点来源的题名、作者、DOI、真实笔记/PDF链接、装置/采样/噪声/率条件及逐篇核验范围集中于[Notebook第11节](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-recovery-strategy-20261019)和[本地证据表](../jupyter/run_records/alpha_gamma_recovery_strategy_20261019/evidence.json)。DOI核PDF首页；实际原页图为Tamagawa3–4、Kim6、Morad10、Fixsen2、Klanner4，其余仅定点文字/既有笔记，不声明本轮重新全文精读。Morad第553项仍仅限检出/分类指标和构造条件，其余AI草稿状态不变。共享库未修改，未联网或追读引用。

**采用范围**：Tamagawa2015（10.1016/j.nima.2015.05.052，PDF3–4）用于总面积归一形状判别对照；99.95% γ保留/0.22% α误留限GAGG/PMT孤立源事件0.5–0.75 MeVee，不是本项目密集堆积成绩。Luo2018（10.1016/j.nima.2018.03.078，PDF4–6）用于多类型假设比较的先例，限BC-501A双脉冲；不采用其会残差抵消的绝对残差总和。Kim2023（10.1016/j.apradiso.2023.110880，PDF6表2/3/5）保留多任务取舍：12 ns分类MT84.80%低于ST87.91%，9重MT总体分类60.22%；原标Acc.H/Recall疑点不消除。Morad2024（10.1109/TNS.2024.3371573，PDF3–7、10）用于检出、上下文分类分开；6 Mcps条件分类89.1%与含漏检72.9%分列，不证明热图分数已校准。

Rosado2019（10.1109/JSEN.2019.2938018，PDF2–5、7）及Klanner2024（10.1016/j.nima.2023.169018，PDF3–4、9）支持SiPM恢复/光时间结构影响响应的机制，验证以指定器件平均电荷为主，不是GAGG逐事件解堆积；不直接移植参数或有疑点的公式。Fixsen2004（10.1016/j.nima.2003.11.313，PDF2–3）只借鉴条件均值/完整协方差标定，限TES模型的局部最小二乘；条件高斯log det及形状分布归一化是本项目拟议统计模型要求，不冒充该文已实现的波形pileup方法。

**项目核查与轻量计算**：重新汇总L31已公开79个γ，物理拟合为成功30/错类型35/正确类型但Q超差4/漏失10；仅类别混合为29/40/2/8。相对物理拟合丢失的8个成功γ全部变错类型；统计区间仍跨零，不宣称总体显著下降。当前混合只用Q条件模板且形状Q初值固定，不使用保存的context_classes；旧物理基线已有Q形状更新，旧context拟合也已试过，不能称本轮才首次想到。旧γ前史参数来自25个有效样本、最大172.20 mV，而公开失败分层有25个≥200 mV目标；这是标定覆盖缺口，不是器件上限。按当前目标，λC=10、假设scoreα/scoreγ=4，会给γ增加27.73白化SSE代价；这是代数量级，不是实测失败分数分布。

**拟议改进，尚未实施**：先固定盲候选，用观测减盲估邻居的残差做独立SI/白化模板/PSD对照；不对已经生成的固定类型模板再分型。用原生Q/可观测前史匹配资料检验幅度捷径，检出和分型输出分开，在独立按ROOT分组的标定集校准融合；同一电压的物理/网络分数不作为独立证据相乘。补弱γ强前史条件响应与不确定度，再检验Q/时刻/形状自洽的双假设比较。自由快慢比例会吸收类间差；若采用条件概率须处理形状分布归一化、条件协方差和形状涨落重复计入，而不是仅增加自由度或给旧经验惩罚机械补项。以上为本项目推导/方案，不是文献性能保证。

**验证与边界**：3个轻量Notebook单元实际执行，旧65单元/图像保留；格式、链接、计数恒等式、来源/输入哈希检查见[本地验证](../jupyter/run_records/alpha_gamma_recovery_strategy_20261019/validation.json)。没有改重建器、训练模型或生成ROOT，默认物理拟合不变。旧公开测试只作失败描述/回归，不择参；新方案冻结后再做独立共同测试，继续要求α/γ完整恢复和误报同时受控。

**新增文献/标定待办**：同GAGG/SiPM原生混合逐事件实测证据、弱γ强前史的条件形状与噪声分解、融合分数校准/重复计入控制仍缺。Kim式3/5、Morad FA分母/NoA先验、Klanner实现异文沿原记录保留；需要新增精读或代码核验时在papers项目另行处理，不从本项目修改共享库。

<a id="free-lifetime-minimum-window-20261019"></a>
### L33：2026-10-19——自由双衰减时间的最短观察窗口诊断

**依据与范围**：沿用L14–L16的完整响应卷积和原S13360两个孤立样例，不新增文献性能结论或物理参数。自由τ是当前电压模型的有效参数；GAGG/PMT平均波形的文献寿命不等于本器件逐事件可辨识精度。原文来源/审核条件保持，原始文献没有规定100/500 ns通用门槛。

**实际计算**：复用原ADC缓存、原约束和四初值，500–750 ns每25 ns扫描，新增22次拟合；两个粒子×500/750 ns共4个历史端点数值回归。T指观测阈值后长度，另用前200 ns且定位/基线依赖更长记录。按两τ均与同事件2.5 μs自由解相差≤5%/10%/20%、无触边且全部更长已测T均满足，两例共同最早已测试T分别为550/550/525 ns。门槛是本轮探索性定义，非置信区间、实物精度或连续意义下最短值；独立粒子仍只有1α+1γ。

**解释与边界**：500 ns自由拟合面积已与长窗相差约3%以内，但γ慢τ仍高约23%；固定τ的100 ns γ慢占比0.3335，相比2 μs的0.2047并不稳定。长窗自由解α80.36/348.05、γ81.41/221.63 ns仍不同于仿真输入95/351 ns，因此本节只回答相对长窗解稳定，未确定真实寿命准确恢复所需T。Poisson前向无新到达概率exp(−RdetT)仅作计数率联系，Rdet非源率、不代表恢复效率，不能由1/T给最大率；后续需独立分Q样本和受控间隔/高率联合拟合。

**追溯与验证**：[原Notebook第8节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#free-lifetime-minimum-window-20261019)、[新汇总](../jupyter/artifacts/alpha_gamma_window_minimum_20261019/summary.json)、[协议与冻结源码](../jupyter/run_records/alpha_gamma_window_minimum_20261019/protocol.json)、[验证](../jupyter/run_records/alpha_gamma_window_minimum_20261019/validation.json)。原29单元/图像保留，未生成ROOT或替换旧结果；有效寿命与本征发光寿命、观察长度与最大率分开，未将同一事件多个截窗当独立样本。

**固定τ补充对照**：同日仅复算原9档固定95/351 ns结果，仍自由拟合A/B/t0/ws，没有新增拟合。相对各自2500 ns固定τ解，要求本点及所有更长已测点的数值均合格，两类共同最早已测档位为面积≤10%：200 ns、面积≤5%：300 ns、面积和ws各≤5%/10%：500 ns。100 ns的γ面积+11.49%、ws+62.88%，不能说固定τ后100 ns普遍足够。α的1000 ns拟合有基线上界触边；若严格加入上一节的无触边条件，共同起点为1500 ns。必须将这个基线/模型质量问题与数值稳定分开，不能把前述300/500 ns称为已通过全部质量条件的最短窗，也不能把1500 ns解释为物理必需长度。两事件、粗网格及长基线依赖不变，无新增文献或总体精度结论。[Notebook固定τ对照](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#fixed-tau-window-comparison-20261019)、[复算与两种标准](../jupyter/artifacts/alpha_gamma_fixed_tau_window_20261019/summary.json)、[验证](../jupyter/run_records/alpha_gamma_fixed_tau_window_20261019/validation.json)。

<a id="short-window-literature-strategy-20261019"></a>
### L34：2026-10-19——缩短拟合观察窗的方法与信息量诊断

**证据与审核范围**：按`paper-search`只读共享库，复用Belli2008、Radeka1968、Fixsen2004、Nakhostin2019 JINST、Yang2022、Saxena2020六篇完整已有笔记，定点核原页；实际看图Belli3–4、Radeka6–7、Fixsen2、Nakhostin4、Yang4、Saxena2，其他位置分清本轮文字/旧笔记层级。方法与适用条件见[Notebook第9节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-literature-strategy-20261019)；完整题名、作者、DOI/未核状态、来源位置和审核记录见[证据表](../jupyter/run_records/alpha_gamma_short_window_strategy_20261019/evidence.json)。两次导航脚本45 s超时后用rg/索引，不据无返回说无文献。六篇仍为AI草稿；第554项认可的是另一篇Nakhostin TNS背景抑制，第658项不转移Fixsen；共享库和审核状态未改。

**采用与限制**：Belli（10.1016/j.nima.2008.06.045，PDF3–4图4–6）支持独立单脉冲截断补尾的误差验证设计；NE213固定7.5/45 ns、200 MS/s条件下7点35 ns/加缓冲50 ns峰距及均值/SD约±8%，不能移为GAGG最短窗。Radeka（DOI本轮未核，PDF6–7/13）用于有限窗口权重及噪声代价原则，串/并联白噪声/冲量信号前提保持；式16–17原印冲突不直接实现。Fixsen（10.1016/j.nima.2003.11.313，PDF2式1–2）用于独立均值/完整协方差标定原则，不作当前SiPM短窗性能证明。Nakhostin（10.1088/1748-0221/14/04/P04012，PDF4/9–10）的NaI离线795 kHz缩尾到20 ns后峰宽16%，100 ns窗含20 ns基线，噪声代价及式2.2/2.3不自洽保留，不照抄为已核代码。Yang（10.1016/j.nima.2022.166376，PDF2–5）的双平均窗尾部预测仍拒绝0.75 μs近邻；Saxena（10.1016/j.nima.2018.09.123，PDF2–5）的HPGe实时自适应成形为1 Mcps输入/370 kcps输出，0.7 μs平顶/2.4 μs复位等待不能照搬。两者只借鉴残尾估计/窗口自适应思想，不据此承诺100 ns逐事件准确恢复。

**项目推导与实际计算**：拟议固定τ后以Af/As线性幅度与少数非线性参数profile搜索、独立基线及条件协方差；相同模型/目标下参数改写不增加信息，不能单独保证精度改善。短窗先截协方差再求逆，不能截长窗逆矩阵；原定位使用整段峰和长前窗，不是严格因果前缀。另实现理想单指数Poisson过程的N/τ Fisher诊断：给定已知起点、任意示例N=10000/τ=351 ns，100 ns总量在τ固定时相对1σ界约2.01%，τ自由时总量约21.25%、τ约24.47%。这是本项目解析说明，非双指数GAGG/ADC下界或实测精度；28项独立积分、7次逆矩阵和2个短窗极限检查通过。

**下一步及未验证项**：先实现完整响应两幅度拟合与只依赖前缀的触发/基线，再检验条件均值/协方差和α/γ形状不确定度，缩尾作有正则的对照；逐事件双τ自由与跨事件共享τ分开。预定扫Tpre=20–200 ns、Tpost=50–1000 ns，按ROOT分标定/开发/新测试，真值只评价；电荷≤5%且事件合格率≥95%/两τ各≤10%只是拟议工程目标。报告两类检测/分类、偏差、分辨率、漏失/误拆/误并及选择效率，校准自适应停止后的覆盖率，不用两条波形的稳定性作“正确”证明。本轮没有新增GAGG拟合、训练或ROOT。[协议（尚未执行）](../jupyter/run_records/alpha_gamma_short_window_strategy_20261019/proposed-protocol.json)、[独立模块](../analysis/short_window_information.py)、[计算结果](../jupyter/artifacts/alpha_gamma_short_window_strategy_20261019/information.json)、[本轮验证](../jupyter/run_records/alpha_gamma_short_window_strategy_20261019/validation.json)。原41 Notebook单元及旧拟合器保留。待补本器件短窗独立原生精度、条件噪声/形状标定、因果总延迟，以及缩尾原式代码/勘误证据；新增精读在papers项目另行处理。


**2026-10-20笔记展示整理**：按用户要求，Notebook 0002聚焦双指数模型，移除单衰减拟合对照及理想信息量诊断的正文／图9-1，合并重复叙述并按研究问题排序；保留其余图号和双指数数值。历史推导、原始结果与验证仍保留在上述独立模块、产物及运行记录中。第9节的第一阶段方案已在下方L35执行；此次整理不新增拟合、文献采用或物理结论。

<a id="short-window-implementation-20261019"></a>
### L35：2026-10-19——因果短前缀拟合的实现与独立原生检验

**文献依据与采用边界**：落实[L34](#short-window-literature-strategy-20261019)的有限窗补尾、线性幅度profile和独立条件均值/协方差原则；Belli NE213固定7.5/45 ns、Radeka理想噪声、Fixsen TES的装置/采样/模型条件与AI草稿审核状态均沿L34，不引入新的文献性能或人工确认。Radeka/Nakhostin原式疑点未消除，未照抄缩尾公式；本轮没有修改共享papers库。α/γ完整模板、条件pilot、50%收缩和停止规则为本项目实现，不冒充论文给出的GAGG性能。

**实现与独立性**：149个旧孤立标定窗、45个旧开发窗，70个Q≥250 pe（52α/18γ）用于三档类别均值/726维协方差；基线均值和电子噪声独立标定，σB=0.2 mV及50%对角收缩为工程设定。四方法是固定95/351 ns两幅度OLS/GLS及独立完整α/γ形状OLS/GLS，输入同一完整原生analog_mV。触发用过去滚动中位数的连续三点2 mV门槛，窗口只读[触发−Tpre,触发+Tpost]；先截协方差再分解，未使用未来峰值/尾部定位。固定类别完整形状比固定双τ且自由比例约束更强，参数改写本身不提高信息量。

配置、标定、4个核心源码哈希冻结后，新增174个独立原生20 μs单事件ROOT（30α/144γ），预定包含γ扩展并保留零信号；单能γ0.1/0.2/0.661657/1.5 MeV、α1.5/1.8/2.5/5.486/8 MeV，三点位置。真值只在该事件44个窗口全部拟合完成后打开，不作为推断输入。Q≥500 pe目标数为α19、γ12；这是实现光子后代雪崩电荷，不是源/沉积/可见能量。

**实际结果**：主方法template_gls。开发选择的前20＋后300 ns候选在新测试通过本批点比例门槛（α19/19（100.0%），γ12/12（100.0%））；目标为两类各≥95%事件Q误差≤5%。测试网格仅作探索描述：前20 ns：后300 ns；前50 ns：后300 ns；前100 ns：后300 ns；前200 ns：后300 ns。小样本点比例和逐条件Wilson区间不能认证总体95%门槛；同事件多个截窗不当作独立样本，未据测试重选算法/停止阈值。前20 ns的关键比较如下。

| 触发后T / ns | 方法 | α电荷误差≤5% | γ电荷误差≤5% |
| ---: | --- | ---: | ---: |
| 100 | physics_ols | 3/19（15.8%） | 0/12（0.0%） |
| 100 | physics_gls | 4/19（21.1%） | 1/12（8.3%） |
| 100 | template_ols | 14/19（73.7%） | 6/12（50.0%） |
| 100 | template_gls | 14/19（73.7%） | 7/12（58.3%） |
| 200 | physics_ols | 10/19（52.6%） | 4/12（33.3%） |
| 200 | physics_gls | 10/19（52.6%） | 7/12（58.3%） |
| 200 | template_ols | 15/19（78.9%） | 3/12（25.0%） |
| 200 | template_gls | 18/19（94.7%） | 12/12（100.0%） |
| 300 | physics_ols | 11/19（57.9%） | 6/12（50.0%） |
| 300 | physics_gls | 18/19（94.7%） | 10/12（83.3%） |
| 300 | template_ols | 17/19（89.5%） | 7/12（58.3%） |
| 300 | template_gls | 19/19（100.0%） | 12/12（100.0%） |
| 500 | physics_ols | 11/19（57.9%） | 3/12（25.0%） |
| 500 | physics_gls | 19/19（100.0%） | 11/12（91.7%） |
| 500 | template_ols | 14/19（73.7%） | 8/12（66.7%） |
| 500 | template_gls | 19/19（100.0%） | 12/12（100.0%） |
| 1000 | physics_ols | 15/19（78.9%） | 6/12（50.0%） |
| 1000 | physics_gls | 19/19（100.0%） | 12/12（100.0%） |
| 1000 | template_ols | 18/19（94.7%） | 10/12（83.3%） |
| 1000 | template_gls | 19/19（100.0%） | 12/12（100.0%） |

**停止与负对照**：开发最坏窗口误差得到k=54.45149，残差阈值2.10599，固定前50 ns；新测试α接受0/19、γ接受0/12，拒绝保留为无条件失败。局部σ和经验放大不能视作已校准顺序置信区间。24个标定均值相加的双脉冲（72次前缀检查）只检查单脉冲假设遇邻居的失效，不是新原生共享微胞双脉冲恢复验证。

**验证与限制**：4项有独立数值参照/未来污染的算法检查通过；174个ROOT的实际RunInfo、宏、种子、配置、哈希与祖先电荷核验，1408个汇总组独立重评分，4次原生未来NaN/平移检查及2次快照源码重放通过。[Notebook第10节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-implementation-20261019)保留原50单元和输出；[实验](../pileup/experiments/alpha_gamma_short_window_20261019/README.md)、[测试汇总](../jupyter/artifacts/alpha_gamma_short_window_implementation_20261019/test-summary.json)、[数据验证](../jupyter/run_records/alpha_gamma_short_window_implementation_20261019/native-data-validation.json)。本轮假定至多一个孤立脉冲、不估计数量，尚不能报告盲多脉冲误拆/误并，也不回答free双τ最短窗。旧两条ADC的长窗数值稳定与本轮独立模拟真值精度口径不同，不直接比较缩短倍数。待补有效γ/条件前史/标定不确定度、原生双脉冲、跨事件共享τ与freeτ、实物域偏移及高率联合拟合；不以1/T替代计数率能力。新文献整理仍在papers项目处理。

<a id="advanced-fit-literature-20261020"></a>
### L36：2026-10-20——短窗口GLS之后的条件响应、形状模态与标定统计

**问题与核验层级**：用户询问是否有更好的拟合方法。按`paper-search`只读本地库，复用导航、完整来源笔记及global-review原始答复；两次knowledge_index导航各45秒无输出超时，改用rg，不以无返回判断无文献。五篇正式期刊来源与逐篇实际文本／原页范围、文件哈希、作者、DOI或未核状态、设备／能区／采样／噪声／率条件汇总于[证据表](../jupyter/run_records/alpha_gamma_advanced_fit_literature_20261020/evidence.json)。这是定点证据核对，不是五篇新全文精读或入库；共享库未写入，未联网。过程与候选追加到[笔记0002第9.1–9.2节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#advanced-fit-literature-20261020)，原代码输出和图号保持。

| 来源 | 实际依据与采用范围 | 限制／审核 |
| --- | --- | --- |
| D. J. Fixsen, S. H. Moseley等，2004，*Pulse estimation in nonlinear detectors with nonstationary noise*；[精读](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)／[原文](../../papers/raw/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.pdf)；[DOI](https://doi.org/10.1016/j.nima.2003.11.313) | PDF2式1–2（原页实看）、PDF3文字；条件均值和逆协方差随标定能量插值，用于提出自洽条件拟合。 | TES模型0.5–5 eV、每点5000次，采样／率未给；原式无log det。AI待审核，无本篇人工确认。 |
| Chris H. Whitford，2005，*Implementation of optimal filters for pulse height measurement from non-linear detectors with non-stationary noise*；[精读](../../papers/sources/2005_Whitford_NIMA_Implementation-of-optimal-filters-for-pulse-height-measurement-from-non-linear-detectors-with-non-stationary-noise.md)／[原文](../../papers/raw/2005_Whitford_NIMA_Implementation-of-optimal-filters-for-pulse-height-measurement-from-non-linear-detectors-with-non-stationary-noise.pdf)；[DOI](https://doi.org/10.1016/j.nima.2005.09.001) | PDF2式4–5、PDF3图3（两页实看）、PDF4式9–10文字；局部dS/dE与端点自洽插值，并保留时刻／基线代价。 | TES 2–20 keV；20 keV改善19%相对未知相位峰取值，仍差下界14%；40 eV是指定2 keV间距的能标偏差。无堆积或硬件验证，AI待审。 |
| S. E. Busch, J. S. Adams等，2016，*Progress Towards Improved Analysis of TES X-ray Data Using Principal Component Analysis*；[精读](../../papers/sources/2016_Busch_JLTP_Progress-towards-improved-analysis-of-TES-X-ray-data-using-principal-component-analysis.md)／[原文](../../papers/raw/2016_Busch_JLTP_Progress-towards-improved-analysis-of-TES-X-ray-data-using-principal-component-analysis.pdf)；[DOI](https://doi.org/10.1007/s10909-015-1357-z) | PDF3与6原页实看，PDF5／7文字；少量波形变化分量及能量方向标定，作为结构化形状／协方差候选。 | Mo/Au TES、1 MS/s；Mn Gaussian FWHM3.19±0.14→1.56±0.13 eV，Al无改善；仅第658项限定采用，独立划分／权重与4317/4984计数差仍待核。 |
| Roger Barlow, Christine Beeston，1993，*Fitting using finite Monte Carlo samples*；[精读](../../papers/sources/1993_Barlow_CPC_Fitting-using-finite-Monte-Carlo-samples.md)／[原文](../../papers/raw/1993_Barlow_CPC_Fitting-using-finite-Monte-Carlo-samples.pdf)；DOI未核 | PDF2式8–9原页实看，PDF3条件最大化文字；有限MC模板统计进入辅助似然，借鉴到有限标定误差传播。 | 原文为Poisson分箱计数，不能直接用于相关电压；模板均值profile不等于边缘积分。DOI未核、PII不改写为DOI；原式7等疑点保留，AI待审。 |
| Jaime Rosado，2019，*Modeling the Nonlinear Response of Silicon Photomultipliers*；[精读](../../papers/sources/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.md)／[原文](../../papers/raw/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.pdf)；[DOI](https://doi.org/10.1109/JSEN.2019.2938018) | PDF4式19原页实看；PDF2–3／5／7–8文字；SiPM恢复、触发概率与电荷改变，为条件响应提供物理动机。 | S13360-1325CS/1350CS，LYSO3×3×20 mm／LED，300–2100 keV、400 ns电荷积分；仅平均电荷／电流模型，非6050PE波形逆算法，AI待审。 |


**与当前实现的差别**：已核`analysis/alpha_gamma_short_window.py`，现有流程已有线性幅度／基线profile、连续时刻搜索、三Q档类别模板，以及`C=Qpilot*K(tpilot)+sigma_e²I`和50%对角收缩；不是完全固定且与Q无关的协方差。具体缺口是模板和C按单次pilot冻结，拟合Q／t0不再自洽更新。有效标定70条、γ每档6条与726维协方差提示应优先验证统计稳定性及标定范围。L26已有形状试验α170/290→173/290、γ18/38不变；不把增加自由度写成已经改善。

**采用的研究建议／尚未验证**：优先比较自洽条件均值／噪声，随后比较开发集选择的低秩结构与有限标定不确定度。若采用参数相关条件高斯似然，`rᵀ C(theta)⁻¹r + log det C(theta)`是本项目统计建模，不冒充Fixsen原式；旧同前缀共用固定C的log det为常数。PCA旋转不自动白化；形状模态若已在残差C内，不重复添加；有限标定误差会被多事件共享。固定形状Q恢复、带标定先验／共享τ和逐事件完全free双τ分开验证，不能由300 ns小样本100%推出自由寿命最短窗或高率能力。

预定先用开发数据比较100/150/200/300 ns，在公开旧测试上只回归；方法冻结后另设新测试，逐类报告Q误差、分型、漏检、接受度与区间覆盖，保留失败。未实现新拟合器、未重新拟合、未生成ROOT，尚不能排序本项目真实性能。Scoullar2008会议PPR仅复用既有笔记／原修订记录，公开似然及事件数选择细节不足，未作为可直接复现的首选。

**文献整理待办**：如后续采用，先另行确认相关电压的有限模板联合统计模型、类别／Q／时刻条件协方差、PCA能量方向与分量选择、自由双τ可辨识／覆盖。Barlow DOI和原印疑点、Whitford详细TES依赖等仍待核；补充文献整理在papers项目按独立授权处理，不在本项目修改共享库或自动追读。


<a id="advanced-fit-implementation-20261020"></a>
### L37：2026-10-20——自洽条件似然、低秩协方差与标定均值误差的实现

**文献采用与边界**：落实[L36](#advanced-fit-literature-20261020)中Fixsen/Whitford的条件响应思想、Busch的低维变化方向及Barlow–Beeston的有限模板统计原则，原探测器/采样/噪声/能区/率条件和审核状态全部沿L36，没有新文献性能或人工确认。本轮不重读或修改共享库。参数相关高斯目标`rᵀC⁻¹r+log det C`、残差对角的采样解释和32次分层ROOT bootstrap高斯均值误差是本项目方案；不是Fixsen原式或Poisson Barlow–Beeston的原样电压实现。PCA替代相关协方差，不新增无限制形状幅度，不重复计入同一涨落。均值标定误差跨事件共享，当前未完整重采样分档、协方差、器件系统误差。

**实现和分组**：七条流程含原标定GLS、扩充标定GLS、自洽密集条件似然、PCA8/PCA16和两种均值误差组合。只拟合至多一个孤立脉冲的未知类型/Q/连续时刻/基线，完整双分量发光物理保持。基线解析profile及0.2 mV共同软约束；模板随当前log Q插值，协方差随当前Q/时刻更新，先截前缀再分解。新增420个原生20 μs单粒子ROOT，4个预检查计入标定；168标定＋84开发＋168冻结后测试，保留所有零信号γ，不拒绝ADC触顶。α源1.5/2.5/5.486/8 MeV、单能γ源0.1/0.2/0.661657/1.5 MeV、三点位置，保持原物理/SiPM与2 ns analog_mV输入。

加入旧149标定和45开发后，为317/129窗；有效标定Q≥250 pe从70增加到113，γ各档从[6, 6, 6]增加到[14, 14, 13]。旧174个已公开测试不用于择法或本轮评分。开发按两类等权、100/150/200/300 ns等权的Q±5%且类型正确选中`expanded_gls`，前窗20 ns；开发95%点比例规则选择后300 ns。主短窗目标150 ns预先固定，方法/标定/源码冻结后才生成168个新测试；按ROOT/三个种子隔离，真值在全部方法/窗口推断后打开。

**独立结果**：目标为实现光子后代雪崩Q≥500 pe，而非源/沉积/可见能量。本批测试α16、γ15个目标；保留零电荷γ108个。150 ns联合成功由原标定GLS的alpha16/16、gamma11/15，经扩充标定GLS的alpha16/16、gamma9/15，到所选方法的alpha16/16、gamma9/15。完整关键窗口如下，同事件截窗不当作独立样本。

| 触发后T / ns | 同批方法 | α Q±5%且类型正确 | γ Q±5%且类型正确 |
| ---: | --- | ---: | ---: |
| 100 | legacy_gls | 14/16 | 10/15 |
| 100 | expanded_gls | 12/16 | 10/15 |
| 150 | legacy_gls | 16/16 | 11/15 |
| 150 | expanded_gls | 16/16 | 9/15 |
| 200 | legacy_gls | 15/16 | 15/15 |
| 200 | expanded_gls | 15/16 | 15/15 |
| 300 | legacy_gls | 16/16 | 15/15 |
| 300 | expanded_gls | 16/16 | 15/15 |
| 500 | legacy_gls | 16/16 | 15/15 |
| 500 | expanded_gls | 16/16 | 15/15 |

150 ns配对比较：alpha对legacy_gls差+0.00个百分点、配对ROOT重采样95%范围[+0.00,+0.00]；alpha对expanded_gls差+0.00个百分点、配对ROOT重采样95%范围[+0.00,+0.00]；gamma对legacy_gls差-13.33个百分点、配对ROOT重采样95%范围[-40.00,+13.33]；gamma对expanded_gls差+0.00个百分点、配对ROOT重采样95%范围[+0.00,+0.00]。范围仅表示指定模拟设计下有限ROOT统计，不是含标定/器件系统误差的总体保证。所选方法在新测试网格中两类≥95%且所有更长档通过的探索性最早点为后300 ns；未据此重新选择方法或窗口，小样本100%不能认证总体≥95%。电荷偏差/SD为匹配条件统计，联合成功和覆盖保留漏检为失败。

**区间与失败控制**：所选方法局部名义95%区间：alpha后150 ns覆盖14/16、有效区间16/16；gamma后150 ns覆盖10/15、有效区间15/15；alpha后300 ns覆盖7/16、有效区间16/16；gamma后300 ns覆盖12/15、有效区间15/15。区间条件于当前类别和模型，不声称已覆盖类别切换/模型失配或顺序停止。所选方法在150 ns、108个零信号窗上有0个触发、0个Q≥500 pe误报。优化器状态、起点/Q触边、无效区间和分Q结果全部保留。另有24条标定均值相加＋白噪声双脉冲、48个前缀负对照，只评分第一个脉冲，未模拟共享微胞非线性或恢复第二个脉冲，不是原生盲解堆积。

**验证和结论适用性**：9项独立数值参照检查、4个旧开发端点回归通过；420个ROOT的RunInfo/实际宏/三个种子/配置/原文件哈希/祖先电荷核验，280组独立重评分，4组原生未来NaN/时移/缓存重放和2组冻结源码重放通过。原短窗源码及全部历史结果哈希保持，Notebook原51单元的代码/输出/图像保留，新增图11-1/11-2/11-3。只验证同G4/SiPM模型、有限能区/位置的孤立固定形状Q；未给free双τ最短窗、实物精度、原生高率误拆/误并或最高源率，不用1/T替代计数率能力。下一步仍需条件噪声/标定系统误差、跨事件共享τ、原生双脉冲及混合高率的独立开发与测试。

**复现与记录**：[实验及冻结源码](../pileup/experiments/alpha_gamma_advanced_fit_20261020/README.md)、[Notebook第11节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#advanced-fit-implementation-20261020)、[完整报告](../jupyter/artifacts/alpha_gamma_advanced_fit_20261020/report.json)、[验证](../jupyter/run_records/alpha_gamma_advanced_fit_20261020/native-data-validation.json)。文献追加整理仍在papers项目按独立授权处理。


<a id="shared-tau-implementation-20261020"></a>
### L38：2026-10-20——跨事件共享双τ与独立标定约束的实现

**文献与项目方案边界**：延续L36–L37中Fixsen2004/Whitford2005条件响应与统计标定的有限借鉴范围；原TES设备、模拟/采样/率条件、缺口及AI待审核状态保持，不新增人工确认或文献性能。本轮跨事件共享τ、log τ高斯约束、逐ROOT线性化bootstrap与慢分量比例分类器是本项目统计方案，不声称原论文已经给出此GAGG方法或普适短窗长度。共享papers库未修改，也未新增原文核验。

**方法改变**：保持双分量发光和原生SiPM物理不变；同样的单位面积双指数与读出/孔径卷积，快慢幅度非负，基线/时刻逐事件拟合。比较固定95/351 ns、独立标定联合共享τ、逐事件带log τ约束、逐事件free双τ，完整类别模板GLS另作外部对照。四种物理方法共用同一前缀的电压pilot冻结协方差，先截矩阵后求解；variable projection消去两幅度/基线。free仍有fast25–250 ns、slow250.5–2000 ns和触发相对时刻−20至＋2 ns范围，不叫无界拟合。原有模型/结果不变，新图12-1至12-4。

**独立标定与选择**：沿用317标定、129开发，不使用之前174或168个公开测试作本轮择参/评分。113条Q≥250 pe标定在触发后1000 ns联合估计共享τ＝91.743114/364.822935 ns，每事件各自幅度/基线/时刻；Q与类型真值只用于独立标定样本和后置分类器，事件推断不读真值。该τ属于有效电压响应模型，不直接等于材料发光寿命。2000次逐ROOT profile-score/information线性化bootstrap近似标定log τ协方差，并加0.02 log标准差下限；不是完整bootstrap refit或全系统误差。约束宽度1/3/10倍仅由开发按两类等权、100/150/200/300 ns联合Q±5%且类型正确选择，冻结为3倍。分类器按各流程独立标定的三Q档慢比例均值/SD，先拟合后分型，未把真类型输入事件拟合。

**新测试与实际结果**：冻结源码/标定/物理快照后生成168个原生20 μs单粒子ROOT（24α/144γ），源α1.5/2.5/5.486/8 MeV、γ0.1/0.2/0.661657/1.5 MeV，三点位置、5 μs发射、2 ns analog电压；保持所有零信号γ和完整动态范围，不拒绝削顶。仍假设至多一个孤立脉冲，不估计数量。下表以实现光子后代雪崩Q≥500 pe为分母，不是能量分辨率；同事件各截窗不作独立样本。

| 触发后T/ns | 方法 | α Q±5% | γ Q±5% | α Q＋类型 | γ Q＋类型 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 100 | fixed_tau | 3/17 | 3/12 | 3/17 | 3/12 |
| 100 | shared_tau | 5/17 | 2/12 | 5/17 | 2/12 |
| 100 | constrained_tau | 5/17 | 1/12 | 5/17 | 1/12 |
| 100 | free_tau | 3/17 | 2/12 | 3/17 | 2/12 |
| 100 | template_gls | 11/17 | 8/12 | 11/17 | 8/12 |
| 150 | fixed_tau | 5/17 | 6/12 | 5/17 | 6/12 |
| 150 | shared_tau | 7/17 | 5/12 | 7/17 | 5/12 |
| 150 | constrained_tau | 9/17 | 3/12 | 9/17 | 3/12 |
| 150 | free_tau | 5/17 | 1/12 | 5/17 | 1/12 |
| 150 | template_gls | 15/17 | 8/12 | 15/17 | 8/12 |
| 200 | fixed_tau | 8/17 | 7/12 | 8/17 | 7/12 |
| 200 | shared_tau | 9/17 | 4/12 | 9/17 | 4/12 |
| 200 | constrained_tau | 13/17 | 4/12 | 13/17 | 4/12 |
| 200 | free_tau | 1/17 | 0/12 | 1/17 | 0/12 |
| 200 | template_gls | 15/17 | 10/12 | 15/17 | 10/12 |
| 300 | fixed_tau | 16/17 | 9/12 | 16/17 | 9/12 |
| 300 | shared_tau | 15/17 | 6/12 | 15/17 | 6/12 |
| 300 | constrained_tau | 15/17 | 5/12 | 15/17 | 5/12 |
| 300 | free_tau | 2/17 | 3/12 | 2/17 | 3/12 |
| 300 | template_gls | 16/17 | 11/12 | 16/17 | 11/12 |
| 500 | fixed_tau | 17/17 | 11/12 | 17/17 | 11/12 |
| 500 | shared_tau | 17/17 | 11/12 | 17/17 | 11/12 |
| 500 | constrained_tau | 17/17 | 11/12 | 17/17 | 11/12 |
| 500 | free_tau | 9/17 | 8/12 | 9/17 | 5/12 |
| 500 | template_gls | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | fixed_tau | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | shared_tau | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | constrained_tau | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | free_tau | 16/17 | 10/12 | 16/17 | 9/12 |
| 1000 | template_gls | 17/17 | 12/12 | 17/17 | 12/12 |

**寿命与最短窗的解释**：固定/shared的τ来自设定/独立标定，不记为单事件寿命恢复。原生free τ与G495/351 ns的偏差混合有限统计、不可辨识、优化触边及前向响应失配；不能仅据长窗收敛认证材料寿命正确。每类各≥95%并且所有更长档合格的探索性最早网格、分Q结果、偏差/SD/RMSE、类型错误、漏失、触边及局部区间覆盖全部保存；不据测试重新选择窗口/方法，小样本点比例不认证总体95%。

**控制与验证**：另有24个匹配均值模型标定及48个独立测试，单位面积双指数卷积＋0.1 mV白高斯电子噪声，无光子/微胞涨落；Q为均值面积参数，观察参考时刻预设，不混作原生实现Q或真实触发延迟。这一对照分开评价Q±5%和两个逐事件τ±10%。另24条均值相加双脉冲、48个前缀只检查单脉冲模型失效，不恢复第二个事件，也无共享微胞非线性。9项独立数值检查、168个ROOT的RunInfo/三个种子/实际宏/配置/祖先电荷、360组独立重评分、4组原生未来NaN/时间平移和2次冻结源码重放通过。

**采用范围和待办**：本次如实保存成功与失败，不把共享τ当作每事件free测量或认为约束必然改善。后续需拆开未观察尾部的随机涨落、模型失配与类别退化，校验真实SiPM前向响应/占用依赖及标定系统误差，再研究原生混合高率；逐事件能量标定、实物、盲多脉冲误拆/误并和最高源率未验证，不用1/T替代计数率能力。新文献精读/入库仍在papers项目按独立授权处理。

**复现**：[Notebook第12节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#shared-tau-implementation-20261020)、[实验及冻结源码](../pileup/experiments/alpha_gamma_shared_tau_20261020/README.md)、[完整结果](../jupyter/artifacts/alpha_gamma_shared_tau_20261020/report.json)、[独立验证](../jupyter/run_records/alpha_gamma_shared_tau_20261020/native-data-validation.json)。


**2026-10-20笔记整体复核与整理**：Notebook0002将当前结论、三批独立样本和评价定义前置，区分两事件长窗稳定性、原生实现电荷精度、有效电压τ及匹配均值对照；把第9节旧待实施计划更新为第10–12节已完成／尚未完成状态，精简重复解释，修正条件均值和约束逆矩阵公式排版。全部32个代码单元、保存输出与17个固定图号保持。此次仅整理展示，不新增拟合、仿真、文献采用或科学性能结论。
