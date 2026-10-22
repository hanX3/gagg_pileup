# Pure-alpha correlated fitting, local population search and light GRU

本实验实现笔记0003第18节的三个建议：把无权残差改成完整相关协方差加权；加入有预算的局部差分进化（DE）与SLSQP比较；增加轻量GRU与旧冻结CNN的对照。过程和图19-1至19-5记录在本地 `jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb`。本目录编号沿用2026-10-20准备日期，共同冻结与新原生测试实际发生于2026-10-21。

## 模型和评价

- p–¹¹B有效单α源谱，无γ源。固定发光双衰减95/351 ns、慢分量0.7401741260169167；物理/光学/读出沿用父实验，不研究free-τ或最短窗口。
- 实际共享微胞原生响应，20 μs记录、2 ns采样，中央[5,15) μs评分。采集动态范围充分，盲重建用完整 `analog_mV`，不按ADC触顶拒绝。
- 源率0.1/1/10/20/40/60/80/100 MHz，不是探测信号率或反应率。每率新ROOT数为8/3/2/2/2/2/2/2，共23。
- Q是归属初级粒子的实际雪崩等效电荷，不是沉积能量。目标Q≥500 pe，完整恢复要求一对一时间差≤10 ns及相对Q误差≤20%。同时保留漏检、额外输出、误合并/误拆分操作标志、条件Q/时间误差和五种循环错位对照。
- 已知单α标定48窗（其中33窗有观测候选用于协方差）；开发24单α及1/10/40 MHz各一窗。协方差和K代价只在这些数据中估计/选择。旧23窗测试已公开，不再视为未见测试。
- GRU使用原CNN的26训练/10验证ROOT和相同五特征、损失、2400步、batch16；架构/种子不同，非等墙钟训练。7435参数，8 ns循环token及2 ns局部卷积头；阈值/权重只由旧验证集选择。

原神经划分的100 MHz只有1个训练ROOT、没有同率验证，最高验证率80 MHz；100 MHz按外推测试解释，不用新测试重选阈值。

协方差是本项目的经验残差模型：

```text
C[i,j] = s[i] s[j] rho^abs(i-j) + sigma_e^2 delta[i,j]
s[i]^2 = v_dark + eta max(mu_observed_pilot[i], 0)
```

整窗C由不读真值的原始观测pilot固定，所有K/时间/幅度方案比较使用同一个C，logdet为常数。innovation白化与非负幅度变量投影改变了拟合目标；局部方案使用边缘协方差、最终接受使用完整记录残差。增删K的经验代价重新标定，非χ²显著性检验。`v_dark=0`是非负回归边界，原生2 MHz暗计数仍启用。

该AR(1)/线性方差近似包含模板失配，尚非真实SiPM微观协方差。白化不足也作为结果保留。盲原生搜索实际提出局部一/二分量；最多12次DE按采集时间分布，不能保证全局最优。另做的已知K=2/3诊断有更强条件，不能算作未知K三脉冲恢复。

## 实际快照和本地数据

- `source/g4/`、`bin/gagg`：沿用已验证30M候选容量版，哈希在manifest记录；该容量保护已由父实验验证，不改变物理参数。
- `macros/`、`logs/`、`manifest.json`：23个实际串行G4任务、三个种子、命令、环境和状态，单ROOT单WaveformEvent；启动间隔至少1秒。先通过1/10 MHz两个预检查，再启动其余，预检查包含在23窗中。
- `inference-source/`：推断当时实际Python源码。`frozen-analysis/`保存配置、标定、开发选择、GRU权重和原局部K配置；`reference-cnn/`保存旧CNN权重/阈值。无须切换Git提交来复现推断。
- `frozen-fit-settings.json`：共同冻结时间及源代码/配置/权重SHA256，测试启动晚于冻结。
- 本地 `jupyter/artifacts/alpha_correlated_global_20261020/`：电压与评分真值分开导出、预测/模型/评分、软件控制、图像和完整表。
- 本地 `jupyter/run_records/alpha_correlated_global_20261020/`：原工作区补丁、环境记录、源码加速前后等价检查、验证脚本及Notebook原文备份。
- 新ROOT原件仍保存在项目 `g4/data/`。读取使用经过SHA256核对的普通本机字节副本，没有共享盘临时链接。真值的冗余时间恢复仅在原分支解压确实失败且完整精确记录可唯一一对一匹配时允许，发生情况单独记录。

运行使用项目本机Python环境。训练/神经网络推断另需记录版本的CPU PyTorch 2.14.1；本次从已有本机wheel缓存恢复到普通临时目录，未改动共享虚拟环境。基础拟合、读保存表和绘图不依赖PyTorch。ROOT、Notebook、运行记录和权重等本地数据并不随Git仓库分发，复跑需保有上述数据及运行时依赖。

## I/O故障与恢复

19个G4任务正常结束后，60 MHz第二个ROOT的全局光子x分支解压失败。按共同时间/y/能量的完整存储值，将16,534,845条记录与PrimaryPhoton做唯一一对一映射，恢复另一份已有的x值；原ROOT不改写，电压4数组逐位不变。只读适配在公共字段不唯一、双份均坏或覆盖不完整时拒绝恢复，6项恢复/拒绝测试已通过。

`io-source/`和`io-recovery.json`单独保留故障后的评分I/O适配及其哈希，不修改原冻结推断源码/配置/权重。checked hash定位之后仍检查完整记录，相对原导出器在正常0.1/10/40 MHz共48个数组上逐位相同；没有依据模型填补评分真值。其余解码异常的实际恢复字段逐ROOT记录在manifest/metadata。

原19个种子不重跑、不更换；剩余四个预定80/100 MHz种子沿用冻结物理/二进制/宏，由2个G4进程生成后顺序导出，降低4核虚拟机的内存争用。所有23个预定窗口均进入比较，历史ROOT和失败导出记录保留。

## 安全复跑到新目录

从项目根目录执行，`--output`必须是不存在的项目相对目录。

```bash
# 使用冻结源码/配置，对首个1 MHz测试窗重跑五种物理方法。
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_correlated_global_20261020/replay_analysis.py one --output jupyter/artifacts/alpha_correlated_replay_new

# 指定原始ROOT文件名主干；不复制或覆盖ROOT。
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_correlated_global_20261020/replay_analysis.py one --root 20261021_08h25m09s --output jupyter/artifacts/alpha_correlated_replay_named_new

# 对全部原生导出、神经推断及软件诊断重跑，需要PyTorch。
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_correlated_global_20261020/replay_analysis.py all --output jupyter/artifacts/alpha_correlated_replay_all_new

# 仅重跑软件诊断，可分别使用新的输出目录。
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_correlated_global_20261020/replay_analysis.py controls --output jupyter/artifacts/alpha_correlated_controls_new
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_correlated_global_20261020/replay_analysis.py optimizer_controls --output jupyter/artifacts/alpha_correlated_optimizer_new

# 只读取本轮已完成结果、重画固定编号图，不训练、不新增ROOT。
bash jupyter/run_gagg_python.sh -m pileup.analysis.alpha_correlated_report
```

开发、冻结、原生运行器分别为 `alpha_correlated_study.py`、`alpha_correlated_native.py`，保留原始协议；不要对既有输出重新执行prepare/freeze或改动既有实验快照。新调参/物理条件须新建实验和冻结后新测试。

## 软件诊断的边界

216窗匹配模板控制分白噪声与独立生成的相关噪声，每种12空窗、12单脉冲、84双脉冲。间隔4/8/12/20/40/80/200 ns，Q1=4000 pe，Q2/Q1=0.125/0.25/0.5/1/2/4，各2次。这些控制不含原生光子/共享微胞统计，不能宣布探测器分离极限。

36个已知K=2/3条件诊断比较相同触发区域和三组初值下的SLSQP与DE＋SLSQP；真时刻/幅度不作输入，DE花费额外评估，不是等计算预算比较。七种盲原生方法仍为原SSE、原局部K-SSE、对角WLS、相关GLS、GLS＋DE、旧CNN和轻量GRU。

之后单独冻结补充协议并生成36个新噪声案例，每例3个优化随机种子，在严格相同的500次实际白化NNLS变量投影评估下，对照随机多初值SLSQP与DE/局部混合。非法时间向量开销另记，不声称严格同墙钟；54次/K共享18个噪声实现，不当54个独立噪声样本。这是已知K条件精度的补充，没有改变原生盲流程或阈值，实际源码保存于`supplemental-search-source/`。

高率每档仅2窗，恢复总数受不均衡曝光影响；条件误差含选择效应，高密度有偶然时间配对。未做沉积能量标定、free-τ窗口极限或实物验证。文献依据与审核范围见 `docs/literature-use.md` L41/L42；Cleland的已知模板协方差、Huang的条件种群搜索和Regadío的GRU只作方法启发，不转移原文装置中的性能数字。

## 本轮实际结果

| 方法 | 完整恢复/目标 | 漏检（无时间匹配） | 额外输出 | 合并/拆分标志 |
|---|---:|---:|---:|---:|
| 原始SSE | 486/4333 | 2480 | 253 | 871/0 |
| 局部K-SSE | 726/4333 | 2118 | 282 | 978/7 |
| 对角WLS | 776/4333 | 2043 | 269 | 1010/4 |
| 相关GLS | 636/4333 | 2243 | 114 | 962/1 |
| GLS＋局部DE | 635/4333 | 2244 | 116 | 962/1 |
| 旧CNN | 1038/4333 | 1448 | 604 | 960/113 |
| 轻量GRU | 982/4333 | 1837 | 199 | 1006/38 |


上述非均衡曝光合计只描述本次样本，逐率及条件误差见本地 `report-tables.json`、`test-summary.json`与固定图。100 MHz实际/五种错位均值为：局部K-SSE 4.73% / 2.79%；相关GLS 3.49% / 2.31%；GLS＋局部DE 3.35% / 2.33%；旧CNN 17.60% / 10.09%；轻量GRU 15.93% / 7.04%。不能以低残差、条件Q标准差或偶然时间匹配代替完整事件效率。

已知K软件诊断：K=2全部分量完整恢复1/18→13/18，评估次数944→4163；K=3全部分量完整恢复1/18→7/18，评估次数1515→6174。与盲原生结果分开，DE额外评估。白化后相邻残差相关0.403，经验C不足，不根据这批测试继续调参。

追加新的36个噪声案例、每例3个优化种子，在相同500次实际NNLS profile评估预算下比较随机多初值SLSQP与DE混合：K=2：随机多初值SLSQP 36/54，DE混合44/54；K=3：随机多初值SLSQP 20/54，DE混合26/54。每K的54次运行共享18个噪声实现，不视为54个独立样本；该补充明确已知K，不重新选择原生流程参数。

验证重新计算161份保存评分、69个加权模型/目标、23窗祖先电荷守恒，复读三档ROOT、核对新旧冻结源码和冻结前后隔离；禁止NumPy文件加载的冷盲拟合及本目录冻结源码新路径复跑均通过。Notebook只执行本次追加单元，历史94单元与输出不变。ROOT—实际宏索引与RunInfo核验随本地Run_Catalog保存。
