# 低电荷、位置和前史标定后的 α/γ 拟合

旧α标定最低约1546 p.e.，4个已公开开发窗中78/161个α目标在其下方。本实验扩展标定覆盖，比较更新均值模板和按电荷/可观测前史约束形状的拟合。α与γ同等重要，选择指标为两类×四条件的完整恢复F1等权平均。

冻结后的16个新混合窗口中，expanded相对原V3的α完整恢复为374/656（57.01%）→378/656（57.62%），γ为30/64（46.88%）→31/64（48.44%）；额外预测55→44，漏检120→117，误合并标志58→56。α完整预测比例373/563→377/555，γ30/91→31/91，γ预测仍有很大误分类/误恢复负担。各条件并不一致：高α率轫致辐射背景的γ反而4/16→3/16，单能γ背景为2/17→3/17；两种高α率条件合计γ仍仅6/33（18.18%）。按条件分层、成对重采样窗口得到恢复率差的95%百分位区间：α为−0.45至+1.65个百分点，γ为−3.13至+5.80个百分点，均跨零。该区间只包含本次窗口波动，不含标定和物理系统误差，不能据此确认稳定提升。

保留扩展标定及全部失败候选，暂不替换V3默认流程。新增标定补齐了原低电荷α范围，但均值修正尚不足以解决密集混合波形：测试后分组中，500–1000 p.e.的α只恢复38/157；与最近另一Q≥50 p.e.脉冲相距小于50 ns的γ仍为9/36，其中18例被错分型；前史电压≥200 mV的γ仅3/21，其中13例错分型。这些是真值辅助的相关性描述，不能当成因果证明或物理时间分辨极限。两个前史形状候选在开发集的α恢复均从92/161降至83/161，γ仍5/12，因此未进入本轮独立测试。每条件仅4个混合窗，原生双脉冲验证中γ目标也只有2例。

## 数据与物理条件

新209个ROOT：9个预检查、93个孤立标定、48个原生前史标定、27个孤立验证、16个双脉冲验证和16个冻结后混合测试。另完整复用旧独立标定56窗，训练合计197窗；没有将旧验证或测试转入训练。新孤立γ48次注入只有8个Q≥50 p.e.，触发最少12例保护检查后，按[方案调整](calibration-protocol-amendment.json)完整加入旧标定。零沉积事件全部保留。

孤立α源能量1.5/1.8/2.5/5.486/8 MeV，γ为0.1/0.2/0.661657/1.5 MeV；源位置(0,0,0)、(4,0,0)、(4,4,0) mm，沿+z入射。前史覆盖AA/AG/GA/GG，发射间隔30/100/300/1000 ns，前者中央，后者中央或边缘；独立验证为60/200 ns。实际首光子间隔受飞行/输运影响。旧标定均为中央，池化回归不是连续位置均匀分布的精确平均。

G4物理、光产额、95/351 ns基底、α/γ发光慢分量0.70/0.21、S13360-6050PE v2、恢复参数和读出全部沿用。只在本目录冻结的brems生成器中增加`/gagg/source/scheduleFile`接口；时间表格式示例见[schedules](schedules)。时间表与每次实际宏、源码和哈希一起保存，所有粒子在同一原生共享微胞窗口数字化。单脉冲新旧入口同种子导出数组逐值一致。每个ROOT一个20 μs窗，最多8个独立串行G4，文件名秒级且启动间隔≥1秒；ROOT仅写根目录g4/data。

盲重建仅输入时间、完整analog_mV和独立冻结标定；没有ADC触顶拒绝规则。Q是实现的满恢复微胞等效电荷p.e.，不等于源能量、沉积能量或已标定可见能量。标定阶段使用祖先真值分离形状；这些标签不进入未知混合波形的推断。

## 方法及结果

`expanded`沿用V3搜索和旧AR(1)噪声，只更新池化标定的电荷相关均值及平均延迟。`context_0.5`和`context_1`增加前史电压H（估计起点前[-24,-4) ns的平均电压减基线，最低0），采用`f = b0 + bQ log(Q/1000) + bH H/100mV`和非负的常数+1/Q残差方差。可变形状约束为中心±3倍缩放散布及[0.02,0.98]；Q/H形状外推取边界，估计电荷本身不截断。两候选尺度0.5/1，固定盲初解时刻/权重，用非负端点模板组合重估类型和电荷。位置不作为推断输入，前史电压也不能唯一确定局部微胞占用。

| 开发方案 | 等权完整恢复F1 | α完整恢复 | γ完整恢复 | 额外预测 |
| --- | ---: | ---: | ---: | ---: |
| baseline | 0.62757 | 92/161 | 5/12 | 12 |
| expanded | 0.64358 | 93/161 | 5/12 | 11 |
| context_0.5 | 0.58422 | 83/161 | 5/12 | 7 |
| context_1 | 0.58319 | 83/161 | 5/12 | 8 |

[冻结记录](frozen-fit-settings.json)选择`expanded`，早于16个新测试生成。两种γ（轫致辐射、661.657 keV）×α源1/10 MHz×4窗；γ源均10 MHz。源注入率不是探测率。评价发射[5,15) μs、真值Q≥500 p.e.；一对一时间匹配≤10 ns、类型正确且Q误差≤20%。

| γ来源 | α源/MHz | 类别 | 原V3 | 新标定候选 |
| --- | ---: | --- | --- | --- |
| brems | 1 | alpha | 33/39 (84.6%) | 33/39 (84.6%) |
| brems | 1 | gamma | 12/18 (66.7%) | 13/18 (72.2%) |
| brems | 10 | alpha | 158/296 (53.4%) | 161/296 (54.4%) |
| brems | 10 | gamma | 4/16 (25.0%) | 3/16 (18.8%) |
| mono662 | 1 | alpha | 29/35 (82.9%) | 29/35 (82.9%) |
| mono662 | 1 | gamma | 12/13 (92.3%) | 12/13 (92.3%) |
| mono662 | 10 | alpha | 154/286 (53.8%) | 155/286 (54.2%) |
| mono662 | 10 | gamma | 2/17 (11.8%) | 3/17 (17.6%) |

| 方法 | 类别 | 完整恢复 | 类型纯度 | 完整预测比例 | Q偏差/标准差 | 漏检 | 错分型 |
| --- | --- | --- | --- | --- | --- | ---: | ---: |
| baseline | alpha | 374/656 | 496/563 | 373/563 | 6.7% / 51.5% | 113 | 46 |
| baseline | gamma | 30/64 | 33/91 | 30/91 | 27.2% / 41.6% | 7 | 24 |
| expanded | alpha | 378/656 | 498/555 | 377/555 | 5.4% / 43.2% | 111 | 46 |
| expanded | gamma | 31/64 | 35/91 | 31/91 | 26.2% / 42.6% | 6 | 23 |

| 方法 | 目标 | 时间匹配 | 漏检 | 额外预测 | 误拆分标志 | 误合并标志 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 720 | 600 | 120 | 55 | 0 | 58 |
| expanded | 720 | 603 | 117 | 44 | 0 | 56 |

Q偏差/标准差是时间匹配子集的条件统计，不是能量分辨率。按条件分层、成对重采样完整窗口4000次的差异区间保存在summary；只反映小样本窗口波动，不含标定不确定性及物理系统误差。context计算时间包括expanded初始化。若context选中，其expanded中间结果作为预先声明的附加消融保存，不按测试结果重新选型。

## 复现与记录

从项目根目录使用本机GAGG Python环境。环境身份见冻结JSON，代码和配置在[source](source)，开发阶段版本在development-source-v1/v2。原始ROOT、导出NPZ、运行记录及Notebook按项目约定只在本地保留。完整生成/分析驱动脚本位于本地`jupyter/run_records/alpha_gamma_context_calibration_20261016/`；实际运行行、宏、随机数和输出哈希在[manifest](manifest.json)。

```bash
# 验证本机Python及库可用；下列命令均从项目根目录执行。
bash tools/run_gagg_python.sh pileup/experiments/alpha_gamma_context_calibration_20261016/run_analysis.py pileup.analysis.alpha_gamma_context_report --folder pileup/experiments/alpha_gamma_context_calibration_20261016 --artifacts jupyter/artifacts/alpha_gamma_context_calibration_20261016 --development-inputs jupyter/artifacts/alpha_gamma_regularized_20261015/inputs

# 使用保存的ROOT/NPZ和冻结模型重跑盲拟合；新输出目录避免覆盖历史结果。
bash tools/run_gagg_python.sh pileup/experiments/alpha_gamma_context_calibration_20261016/run_analysis.py pileup.analysis.alpha_gamma_context_study --folder pileup/experiments/alpha_gamma_context_calibration_20261016 --inputs jupyter/artifacts/alpha_gamma_context_calibration_20261016/inputs --calibration pileup/experiments/alpha_gamma_context_calibration_20261016/calibration.json --old-calibration pileup/experiments/alpha_gamma_context_calibration_20261016/old-calibration.json --output jupyter/artifacts/context_replay --role test_mixed --methods baseline expanded
```

报告重生成仅改本轮图表，原始ROOT与拟合结果保持。`context_replay`若已存在，会检查源码/输入身份后复用；需要从头复算时选另一个新目录。G4可由本目录source/brems/g4和source/mono662/g4独立CMake构建，不依赖切换Git提交。各实际宏包含原始绝对路径以对应RunInfo；异机重跑时按新项目根重定位配置和时间表，将outputFile改成新的秒级文件名，；同名则等待下一秒重新命名，不覆盖原文件。每个宏仍只执行beamOn 1。不要直接清空原manifest以重跑，应保存独立运行记录。

26项Python测试、1项C++ SiPM测试通过；209个ROOT和220份原生拟合结果复核，另对一个新测试窗完成仅电压重跑。详见本地[验证记录](../../../jupyter/run_records/alpha_gamma_context_calibration_20261016/validation.json)、[Notebook第16节](../../../jupyter/codex_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-context-calibration-20261016)、[Run Catalog](../../../jupyter/Run_Catalog.ipynb#alpha-gamma-context-calibration-20261016)和[文献采用L27](../../../docs/literature-use.md#alpha-gamma-context-calibration-20261016)。

下一步优先补齐强前史下的γ条件标定：本次仅4个有效γ前史标定，前史电压上限约172 mV，未覆盖上述≥200 mV失败组。可从独立标定库选取多份有效γ光子输入，与前脉冲在同一微胞阵列按新数字化种子重放，系统覆盖前史、间隔和位置；这属于条件响应检验，不能拿筛选后的样本估计γ探测效率。在独立重复样本上检验条件均值和协方差是否足以描述模板涨落，再在开发集比较类型推断与电荷分配，并用新的冻结后测试判断是否改善。这些是后续方案，本轮尚未实施；本轮16个新测试窗没有用于调参。

研究限制：三点位置和有限前史网格、弱γ标定/测试统计量、未直接测量的SiPM恢复参数、理想光耦合、局部拟合最优性及真实GAGG模板失配。当前α/γ共用95/351 ns是既有冻结模型设定；文献α衰减不同的问题见采用记录L04。本轮未声称物理极限或实物准确率，未远程上传。
