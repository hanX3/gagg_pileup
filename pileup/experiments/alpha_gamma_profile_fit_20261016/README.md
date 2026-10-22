# α/γ全窗拟合与受约束形状试验，2026-10-16

本实验继续按α、γ逐事件恢复同等重要的目标优化。比较12组开发设置：全窗电荷重估和时刻搜索、噪声协方差、前沿正则化、类内形状变化。完整失败结果保留。开发中只有shape_0.5略优于原V3：等权完整恢复F1 0.62756657→0.63009504；α92/161→93/161、γ5/12不变、额外预测12→11。

冻结时间2026-10-16 10:29:36 +08:00；首个新原生测试10:30:40启动。保持原G4物理、S13360-6050PE v2、2 ns采样、0.5 ns孔径、0.1 mV电子噪声；拟合使用完整analog_mV并假设动态范围足够。Q是实际完全恢复单胞等效电荷，不是沉积能量。报告Q≥500 p.e.，10 ns内类型无关一对一匹配、类型正确且电荷误差≤20%。

## 本轮结果

本轮同一批8个新原生窗比较，受约束可变形状拟合使α完整恢复从170/290（58.62%）升到173/290（59.66%），γ仍为18/38（47.37%）。额外预测均17个，漏检49→48，误拆分标志均0，误合并标志18→17。高α源10 MHz条件合计：α148/262→151/262，γ5/16不变。改进幅度很小，且每条件只有两窗，不能据此声称稳定提高了实际效率。

保留可变形状作为实验候选，暂不替换上一版比较基线。γ正确分型数21→22，但同时满足类型和电荷要求的预测比例由18/49降为18/51；γ的完整恢复F1从0.41379降到0.40449。因此，本轮没有达到α和γ都可靠逐事件恢复的目标。12组试验也表明，在当前标定与数据条件下，更多拟合自由度、更小残差或更复杂的噪声权重不能自动转化为更好的两类恢复。

完整逐类恢复率、纯度、条件电荷误差、漏检/误拆分/误合并及图表见本地[Notebook第15节](../../../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-profile-fit-20261016)和[summary](../../../jupyter/artifacts/alpha_gamma_profile_fit_20261016/summary.json)。文献原则及适用条件见[L26](../../../docs/literature-use.md#alpha-gamma-profile-fit-20261016)。新方法需要先运行V3再做形状后处理，不能把后处理时间当成完整流程时间。

## 数据与源码边界

- 原56个独立孤立标定窗中29条用于标定；原4个混合和18个单粒子开发窗复用。见父实验[manifest](../alpha_gamma_regularized_fit_20261015/manifest.json)。本轮不修改原标定或旧结果。
- 新8个20 μs原生窗：brems/661.657 keV × α源1/10 MHz × 每条件2窗；γ源10 MHz，中央评价10 μs。种子261052000–261052007；源注入率不等于探测率，单能γ借用C12源槽，不是C12俘获谱。
- [manifest.json](manifest.json)保存实际ROOT、宏、种子、二进制及环境。源代码位于[source](source/)，144份科学源码哈希在[frozen-fit-settings.json](frozen-fit-settings.json)；本轮没有改动C++，相同二进制的既有预检查复用。每ROOT一个WaveformEvent，最多4个串行G4进程，启动间隔≥1秒。
- 四个开发阶段的实际新增源码分别位于development-source-v1/v2/v3/v4；v1三种基础权重，v2固定时刻与带状协方差，v3前沿约束，v4类内形状。run_development.py在冻结依赖上覆盖对应阶段源码，避免用新源文件冒充旧运行版本。
- shape_0.5采用原标定慢分量RMS的0.5倍：α0.010084、γ0.024353；边界为均值±3倍所选尺度。固定初始时刻、形状均值和噪声权重；弱分量保留。正则化依赖同一波形的初始解，不能解释为独立先验或校准后的分类概率。
- physical-calibration.json仅用于开发噪声试验；经验方差乘数2.07796含模型失配。正式选中shape_0.5仍采用原AR(1)权重。
- 形状祖先诊断和physical-oracle-development明确使用真值，只是开发诊断，不是盲性能或物理上限。新测试拟合完成后才打开真值，两类分别评分，未用新测试调参。

## 复算

从项目根目录运行，Python环境由tools/run_gagg_python.sh选择。本地ROOT、缓存NPZ、Notebook和图表不随Git仓库分发，需另行提供本实验和父实验的数据。以下复算写入新的replay_test目录，不覆盖原结果：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bash tools/run_gagg_python.sh pileup/experiments/alpha_gamma_profile_fit_20261016/run_analysis.py pileup.analysis.alpha_gamma_model_selection_study --folder pileup/experiments/alpha_gamma_profile_fit_20261016 --inputs jupyter/artifacts/alpha_gamma_profile_fit_20261016/inputs --calibration pileup/experiments/alpha_gamma_profile_fit_20261016/calibration.json --output jupyter/artifacts/alpha_gamma_profile_fit_20261016/replay_test --role test_mixed --penalties 1
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bash tools/run_gagg_python.sh pileup/experiments/alpha_gamma_profile_fit_20261016/run_analysis.py pileup.analysis.alpha_gamma_shape_study --initial jupyter/artifacts/alpha_gamma_profile_fit_20261016/replay_test/model_selection_1 --inputs jupyter/artifacts/alpha_gamma_profile_fit_20261016/inputs --calibration pileup/experiments/alpha_gamma_profile_fit_20261016/calibration.json --output jupyter/artifacts/alpha_gamma_profile_fit_20261016/replay_test --scale 0.5
```

开发阶段复算使用`run_development.py 阶段号 模块名`，其余CLI参数与对应阶段study脚本相同，并指定新输出目录。例如阶段1为`1 pileup.analysis.alpha_gamma_profile_study`，阶段4为`4 pileup.analysis.alpha_gamma_shape_study`；开发初始解来自`jupyter/artifacts/alpha_gamma_model_selection_20261016/validation_mixed/model_selection_1`，输入来自`jupyter/artifacts/alpha_gamma_regularized_20261015/inputs`。shot变体使用physical-calibration.json，其余使用calibration.json。

实际新ROOT生产入口为run_campaign.py，实际双方法分析入口为[run_test_pipeline.py](run_test_pipeline.py)。后者从完成且已核验导出的窗口取波形，在两方法拟合结束后评分；入口哈希见test-driver.json。既有manifest仅补充未完成计划；追加新物理/配置/种子实验应另建目录。

## 验证与限制

24项数值检查通过；20组同模板电子噪声控制两方法均40/40，无额外预测，不能宣称10 ns物理分辨率。8个新ROOT全分支读取，宏/种子/哈希和冻结时序核对；86个原生拟合重新评分并重构模型，40个控制重评分；新窗一例只输入电压完整重放。详见[validation](../../../jupyter/run_records/alpha_gamma_profile_fit_20261016/validation.json)及[tests.log](../../../jupyter/run_records/alpha_gamma_profile_fit_20261016/tests.log)。

每条件两窗，γ阈上数量小；当前标定主要是中心孤立响应，低Q α、位置及前史占用仍覆盖不足。有限离散/连续搜索不保证全局最优。数值收敛、较小残差和理想控制通过都不能替代原生逐事件效率。

自由时刻开发试验的48次局部优化中，17次正常停止，26次达到60次迭代上限，5次返回ABNORMAL。负结果只适用于所用初始化和预算，不能冒充充分收敛后的方法极限；正常停止也不保证全局最优。选中的可变形状后处理固定时刻，不运行这项优化。逐组状态见[optimizer-status.json](../../../jupyter/artifacts/alpha_gamma_profile_fit_20261016/optimizer-status.json)。
