# α/γ 混合堆积：独立标定与稀疏加权拟合

2026-10-15 建立与冻结，2026-10-16 完成同窗基线和记录。针对旧的两示例固定模板、阈下贡献及相关噪声问题，实现电荷相关平均模板、固定协方差加权、非负稀疏候选选择、弱分量保留和无L1条件幅度重拟合。方法借鉴已有文献原则，属于项目改写，未完整复现任何单篇论文。依据和限制见 [L23/L24](../../../docs/literature-use.md#alpha-gamma-regularized-20261015)，过程见本地 [Notebook第13节](../../../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-regularized-20261015)。

## 独立数据与适用范围

全部86个新原生ROOT各含一个20 μs波形：56个单粒子标定（包括4个预检查）、18个单粒子开发检查、4个混合验证、8个冻结后测试。原32个已公开混合窗没有作为本次未接触测试。实际ROOT、宏、三个种子、源码身份和状态见 [manifest.json](manifest.json)，基准Git为 `f4c54d4f21a707206df0b44992c4f1efc62a81ac`，当时未提交改动在 `working-changes.patch` 与 `git-status-before.txt`。独立G4最多4个串行任务，启动间隔至少1秒。

G4物理没有修改，两个父实验的实际C++源码分别保存在 `source/brems/g4` 和 `source/mono662/g4`；使用原冻结二进制，SHA256见manifest。Geant4 11.4.2、ROOT 6.40.04，S13360-6050PE v2；所有粒子共享微胞恢复历史，2 ns采样、0.5 ns孔径、0.1 mV电子噪声。单能γ在独立快照中借用C12槽的661.657 keV光子，不代表C12俘获谱。

孤立标定源为中央单粒子：α为2/3.5/5.486/8 MeV，γ为0.05/0.1/0.2/0.661657/1.5 MeV。实际用于模板的Q≥50 p.e.信号只有16个α与13个γ；覆盖1545.6–16860.2和256.9–1996.7 p.e.，没有位置或微胞占用标定。保存的是源能量，拟合电荷尚未换算为沉积能量。

混合采用有效p–¹¹B单α源1/10 MHz，γ源10 MHz，两种γ谱均比较，每条件2个新测试窗，中央[5,15) μs。源注入率不等于探测信号率。推断只读取完整含噪 `analog_mV`、采样网格和冻结标定，真值在拟合之后用于评分。遵循动态范围足够前提，不以ADC削顶拒绝窗口。

## 冻结方法与结果

模板固定95/351 ns基函数与原读出，标定 `f(Q)=a+b log(Q/1000 p.e.)` 和平均时延。形状使用的Q在标定范围内截取，但拟合电荷本身不截取。标定残差估计异方差AR(1)；候选使用80 ns平滑电压定权，最终权重由纯电压初步模型产生后冻结。因此是可行加权最小二乘，未声称完整参数相关高斯似然。两个类别列统一加权范数后做非负L1，随后无L1 NNLS与局部时刻/类型搜索。保留弱分量，最终只报告Q≥500 p.e.。

4个混合验证窗按 [selection-rule.json](selection-rule.json) 选择λ=3（候选3/8，宏平均完整恢复F1为0.59490/0.47698）。[frozen-fit-settings.json](frozen-fit-settings.json) 于2026-10-15 17:14:52 +08:00创建，早于首个测试窗17:16:24；保存完整配置、环境和源码哈希。选择/重拟合网格4 ns，发现候选和保存模型使用原2 ns；两轮refit、时刻±12 ns、最小间隔8 ns、近邻四类型组合距离100 ns。形状在每次坐标循环内固定，保存 `shape_charge_pe` 以复现；没有全局最优保证、概率校准或自动拒判。

完整恢复要求真值Q≥500 p.e.、类型无关一对一时间匹配≤10 ns、类型正确且相对电荷误差≤20%。Q为完全恢复单胞等效电荷，不是能量分辨率。下表两方法使用相同8窗：

| γ 来源 | α源 / MHz | α：旧 → 新 | γ：旧 → 新 |
| --- | ---: | --- | --- |
| 轫致辐射 | 1 | 14/16（87.5%）→11/16（68.8%） | 4/14（28.6%）→3/14（21.4%） |
| 轫致辐射 | 10 | 81/153（52.9%）→97/153（63.4%） | 1/11（9.1%）→1/11（9.1%） |
| 661.657 keV | 1 | 15/16（93.8%）→13/16（81.2%） | 1/2（50.0%）→2/2（100.0%） |
| 661.657 keV | 10 | 69/155（44.5%）→91/155（58.7%） | 1/2（50.0%）→0/2（0.0%） |

全部369个阈上目标中，时间匹配309→308、漏检60→61；额外预测53→19，误拆分标志0→0，误合并标志24→31。高α率下α改善，低率α退步，γ不稳定。部分条件电荷偏差/标准差也变坏；完整条件误差、纯度、逐目标配对转移和固定波形图保存在Notebook及本地 `jupyter/artifacts/alpha_gamma_regularized_20261015/summary.json`。每条件仅2窗，尤其单能γ总共4个阈上目标，不能据此给稳定效率或最大工作率。

旧基线使用原模板、发现/阈值与非线性搜索，只将NNLS改为带原空间KKT检查的等价Gram/Cholesky或QR缩减。局部非线性搜索不承诺逐位一致；新低率单窗与原直接NNLS回查，数量、类型、配对和恢复数相同，最大时刻差1.30×10⁻⁶ ns、电荷差6.37×10⁻⁶ p.e.。此前32窗原结果保持。新方法总拟合23.93 s；旧方法8窗任务耗时相加2149.60 s。两者运行日期和并行负载不同，不能将比值称为受控加速比或实时能力。

## 失败与实际验证

- 第一版权重对噪声下行涨落敏感，λ=3孤立开发目标只恢复α2/8、γ1/3，额外预测13；修改权重后为8/8、3/3、额外0，与旧基线相同。原失败结果在 `validation_single/`，改后在 `validation_single_v2/`。该检查集已用于开发，不能当独立测试。
- 40次同模板双脉冲控制（不含光子/微胞涨落），2/4 ns拟合网格计数相同。每网格下10/20/40/80/200 ns分别完整恢复0/8、2/8、8/8、8/8、6/8。10/20 ns漏候选，200 ns γα出现局部错解；不能宣称40 ns实物分辨时间。
- 8个测试窗的L1求解器均以成功状态停止，但最大投影梯度0.15813（λ=3）未达到严格1e−5梯度容差。相对目标变化停止不等于KKT完全收敛，支撑选择与后续局部拟合仍是近似解；诊断原样保留，没有按测试重调精度后择优替换。
- 全86个ROOT所有树/分支读回、宏/种子/哈希、独立划分、冻结时序和源码核验；16个测试拟合共同重评分，8个新模型含弱分量的精确重构，一窗只电压重放及13项算法/线性代数检查。证据见本地 [validation.json](../../../jupyter/run_records/alpha_gamma_regularized_20261015/validation.json)。
- 没有单独消融各项改进，没有位置/占用校准、完整分型后验、沉积能标或实物泛化。下一轮优先验证残差增峰/单双模型竞争和多初值，再扩大γ标定与独立测试；不依据已读测试重调本实验。

## 重跑

以下从项目根目录运行。`jupyter/`和ROOT按约定只保存在本地；在新克隆中需取得对应本地数据。各系统独立安装Python/本机G4/ROOT，不共享虚拟环境。

```bash
gagg_exp=pileup/experiments/alpha_gamma_regularized_fit_20261015
gagg_art=jupyter/artifacts/alpha_gamma_regularized_20261015
# 用冻结模型重建8个测试窗；新输出路径保留旧结果。
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bash tools/run_gagg_python.sh \
  "$gagg_exp/run_analysis.py" --stage test_mixed --folder "$gagg_exp" \
  --inputs "$gagg_art/inputs" --calibration "$gagg_exp/calibration.json" \
  --output "$gagg_art/replay_test" --no-baseline
# 原方法在相同新输出目录处理这8窗，每窗独立进程。
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bash tools/run_gagg_python.sh \
  "$gagg_exp/run_baseline_parallel.py" --inputs "$gagg_art/inputs" \
  --output "$gagg_art/replay_test" --workers 4
# 从原保存结果重新生成统计与PNG/PDF，保持拟合和冻结源码不变。
env MPLCONFIGDIR=/tmp/gagg-report-mpl OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  bash tools/run_gagg_python.sh "$gagg_exp/run_report.py" --artifacts "$gagg_art"
```

`source/`保存测试前冻结科学代码。`postprocessing/`另存测试完成后的绘图/表格代码，未覆盖冻结版本。`run_analysis.py`、`run_campaign.py`优先导入本实验快照；初次运行器修正前版本另保存在 `initial-campaign-source.py.txt`。单粒子预检查/标定/验证和混合阶段运行器会跳过已完成计划，旧ROOT不覆盖；改物理或固定配置须另建目录。

若需要重新编译原生程序，可分别对 `source/brems/g4`、`source/mono662/g4` 执行CMake配置与构建（本机须安装manifest记录的依赖/数据集）。现有运行器会核对原二进制SHA256；新编译的不同二进制应建立新的复现实验记录，不改写当前manifest的身份来伪装原运行。`bin/`是本地二进制，已忽略；重拟合已有波形不依赖G4二进制。
