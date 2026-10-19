# α/γ等权恢复：联合模型搜索的独立测试

2026-10-16。针对上一版10/20 ns控制漏候选、200 ns γα局部错解和L1提前停止，实现更严格的线性选择、联合两脉冲时间/类型/电荷搜索及带复杂度代价的残差补峰。α、γ逐事件恢复同样重要，不用总准确率隐藏类别差异。来源及限制见[L25](../../../docs/literature-use.md#alpha-gamma-model-selection-20261016)；研究过程见本地[Notebook第14节](../../../jupyter/codex_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-model-selection-20261016)。

## 数据、划分与冻结

复用父实验56个单粒子标定的冻结模型、18个单粒子开发窗及4个混合开发窗；本轮只新增8个原生测试ROOT。每窗20 μs、每文件一个WaveformEvent、中央评价[5,15) μs。两种γ谱（轫致辐射、661.657 keV）× α源1/10 MHz × 2窗，γ源均10 MHz；这是源注入率。源α为有效p–¹¹B单α谱，源能量不等于入射/沉积能量。单能γ快照借用C12槽，不是C12俘获谱。

[selection-rule.json](selection-rule.json)预先定义类别/条件等权完整恢复F1和与上一版比较的开发门槛。代价1/9结果0.62757/0.62557，上一版0.59490；选择1，固定L1强度3。单粒子开发检查两候选均为α8/8、γ3/3、额外0；7个其他γ入射无阈上目标。[frozen-fit-settings.json](frozen-fit-settings.json)于09:14:48冻结，首个新测试09:15:09启动；新种子261042000–261042007与父实验不重叠。新测试未用于调参。

[manifest.json](manifest.json)是本轮8个新ROOT、实际宏、G4/SiPM/读出种子、状态和哈希的权威清单；[analysis-manifest.json](analysis-manifest.json)是测试前的开发输入快照，当时new_native_runs为空。Git基准f4c54d4f21a707206df0b44992c4f1efc62a81ac，未提交内容由git-status-before.txt、working-changes.patch及实际source快照记录。冻结132份源码包含原生C++和Python分析；不需切换Git提交即可重拟合。

物理保持两个父实验相同C++与二进制；Geant4 11.4.2、ROOT 6.40.04，S13360-6050PE v2、2 ns采样/0.5 ns孔径/0.1 mV电子噪声，共享微胞历史。复用同二进制的4个已核预检查及86文件全读证据。G4内部串行、最多4个独立任务、启动间隔≥1秒。完整含噪analog_mV进入拟合，保持足够动态范围；真值在盲拟合后才打开评分。

## 方法与原生同窗比较

原非负L1目标不变，在L-BFGS-B后用活跃二次子问题校正，投影KKT阈值1e−5。近邻500 ns内四类型比较，联合两次到达时间与电荷；每次局部最小二乘最多35次函数求值、最小间隔8 ns，两轮扫描。固定平均模板形状/权重条件下比较单/双模型，SSE收益须>1，两个分量均≥50 p.e.；每轮至多16个候选、总新增≤12。保留弱分量，报告Q≥500 p.e.。没有完整参数相关似然、全局搜索、已校准p值或后验置信度。

沿用父实验标定：95/351 ns基函数、Q相关慢分量及固定AR(1)近似。仅中心孤立标定，α Q范围1545.6–16860.2 p.e.，γ256.9–1996.7 p.e.，形状Q超界夹住，估计电荷不夹住。候选用2 ns采样、拟合4 ns，最终模型回原2 ns网格；保存全部分量与shape_charge_pe以精确重构。

完整恢复要求真值Q≥500 p.e.、类型无关一对一匹配≤10 ns、类别正确且电荷误差≤20%。Q是实际单胞等效电荷，尚未换算为沉积能量。三种方法处理同一组8个新波形：

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

| 方法 | 目标 | 时间匹配 | 漏检 | 额外预测 | 误拆分标志 | 误合并标志 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 原始拟合 | 346 | 290 | 56 | 54 | 1 | 24 |
| 上一版稀疏 | 346 | 284 | 62 | 25 | 0 | 24 |
| 本轮联合 | 346 | 295 | 51 | 17 | 0 | 21 |

相比上一版，低α率brems为alpha 12/16→14/16、gamma 7/12→11/12。高率为brems：α 81/142→78/142，γ 0/6→1/6；mono662：α 85/143→86/143，γ 0/6→1/6。逐类纯度、条件电荷偏差/标准差、漏检、逐窗结果见本地Notebook与jupyter/artifacts/alpha_gamma_model_selection_20261016。每条件仅两窗，不支持稳定效率或最大计数率结论。

## 控制、诊断和验证

- 已公开的20组理想双脉冲（四类型×10/20/40/80/200 ns，3000/1000 p.e.，仅0.1 mV白噪声、同模板生成/反演）完整恢复24/40→40/40，无额外预测。它们用于开发，不是独立性能；冻结重放与开发记录逐字段相同。不能将10 ns当真实物理分辨率。
- 故意使用真值的开发诊断：提供全部有电荷事件数量/时刻、真实电荷用于形状/协方差、真实类型用于协方差，再比较固定真实类型或双类型NNLS。高率brems固定类别仍只有α53/73、γ1/4完整恢复。该诊断不是盲性能或严格物理上限；提示平均模板失配、涨落或重叠病态性，未做完全分离。
- 开发阈上目标78/161个α低于α标定Q下限、2/12个γ高于γ上限。位置/高占用标定仍缺。下一轮优先补齐两类形状分布及电荷范围，再独立验证受约束形状模型；不只继续降低候选门槛。
- 8个新ROOT全部树/分支、宏/种子/哈希、冻结时序及独立性检查；132份科学源码未变；24个拟合共同重评分/模型重构，一窗只电压重放、20组控制重放及17项算法/数值检查。新方法线性选择达到KKT容差，但有限预算的非线性搜索不保证全局最优。记录见[validation.json](../../../jupyter/run_records/alpha_gamma_model_selection_20261016/validation.json)。

下列是最多4个独立进程下任务墙钟耗时之和，含并发资源竞争，不是总时间、受控加速比或实时吞吐：

| 方法 | 8窗任务耗时之和/s |
| --- | ---: |
| baseline | 2648.47 |
| regularized_3 | 107.40 |
| model_selection_1 | 517.45 |

## 复算

从项目根目录运行。Notebook、ROOT与artifacts按项目约定仅本地维护，新克隆需要取得相应本地数据。每系统使用自己的Python环境；重拟合不依赖G4二进制。

```bash
gagg_exp=pileup/experiments/alpha_gamma_model_selection_20261016
gagg_art=jupyter/artifacts/alpha_gamma_model_selection_20261016
# 三种方法重拟合；另存输出，不覆盖正式结果。
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bash tools/run_gagg_python.sh \
  "$gagg_exp/run_test.py" --inputs "$gagg_art/inputs" \
  --output "$gagg_art/replay_test" --workers 4
# 用冻结方法重放已公开的理想控制。
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bash tools/run_gagg_python.sh \
  "$gagg_exp/run_controls.py" \
  --previous jupyter/artifacts/alpha_gamma_regularized_20261015/controlled.json \
  --output "$gagg_art/replay_controls_new.json"
# 从正式结果复算表格和PNG/PDF。
env MPLCONFIGDIR=/tmp/gagg-model-report OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  bash tools/run_gagg_python.sh "$gagg_exp/run_report.py" --artifacts "$gagg_art" \
  --previous jupyter/artifacts/alpha_gamma_regularized_20261015
```

run_test.py/run_controls.py优先导入source冻结源码；run_validation.py可用旧开发输入重放选参，参数见--help。run_campaign.py按实际保存宏与固定计划运行并跳过已完成文件；改变物理/固定配置或使用新编译且哈希不同的G4程序须另建复现实验，不能改写本次身份。postprocessing单独保存测试后的统计/绘图和哈希，不改source；bin是忽略的本机二进制。共享文献库未写入，本轮未上传远端。
