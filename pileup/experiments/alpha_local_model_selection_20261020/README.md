# 纯α：主动拆峰与局部脉冲数搜索

本实验落实纯α笔记0003第16节提出的拟合升级：在原反卷积/非负幅度/循环调时基线上，增加主动单→双拆分、残差候选、删除/合并及多初值局部联合调时。保持原95/351 ns双衰减和慢分量权重、完整含噪电压输入；真实K、时间、电荷及源注入率不进入拟合。方法是本项目实现，文献依据及其器件/采样/噪声/审核范围沿用 docs/literature-use.md 的L39。

## 数据划分与冻结

独立单α标定48窗、开发单α24窗，复用其他任务中已经独立划分的原生孤立事件；三个已公开的1/10/40 MHz窗口仅用于开发。标定单脉冲最大拆峰SSE收益除以观测拟合电荷，建立局部复杂度代价的电荷系数；开发比较0.5/1/2三个倍率，优先减少孤立窗额外预测，再比较三档开发完整恢复F1宏平均。开发首版的完成结果及实际源码另存本地run_records，不覆盖有效旧数据。

参数选择见本地 artifacts/alpha_local_model_selection_20261020/development-summary.json；正式freeze.json和frozen-fit-settings.json产生后才允许生成新测试。inference-source/保存实际Python依赖，source/g4保存实际C++及配置，bin/gagg与父实验二进制SHA一致。父实验使用Geant4 11.4.2、ROOT 6.40.04；本轮不修改物理模型，不需要切换Git提交来重拟合。

新原生测试23窗：源注入率0.1/1/10/20/40/60/80/100 MHz，各8/3/2/2/2/2/2/2窗；源是当前p–¹¹B有效单α谱，没有γ背景。每个ROOT只含一个20 μs WaveformEvent，中央[5,15) μs评价，保留完整共享微胞历史。G4内部串行，独立进程不超过4个；启动至少间隔1秒，ROOT在g4/data/使用秒级唯一文件名。先验证1/10 MHz两个新预检查窗，再放行其余生产；预检查计入23窗。

宏、源/SiPM/读出种子、RunInfo、ROOT及二进制哈希、启动/结束、命令、日志和状态由manifest.json及实际文件追溯。可执行文件和stdout使用普通本地临时文件，随后复制保留日志，不在共享盘创建临时链接。

## 方法与评价

局部可变K，非负幅度使用variable projection，SLSQP同时优化局部到达时间；局部优化后在完整记录上比较SSE收益。3 ns最小间隔、±24 ns局部调时、1400 ns局部尾部、两轮及160次增量上限是工程预算，不是物理分辨极限。内部干扰分量最低100 p.e.，报告阈值500 p.e.。所有重数均允许局部联合更新；原基线仅≤64候选时有最后一轮全局联合非线性优化。

完整恢复沿用真实电荷≥500 p.e.、一对一匹配±10 ns及电荷误差≤20%。同时保存漏检、额外预测、误拆分/合并标志及条件电荷/时间误差。Q是实际堆积微胞状态下的等效雪崩电荷，尚不是已标定沉积能量。经验SSE代价不是p值、后验概率或全局最优保证。新方法时间包含初始化的原基线时间；并发任务墙钟之和不当作实时吞吐或严格加速比。

174条固定平均模板控制（24空、24单、126双）仅含0.1 mV电子白噪声，单列于本地controls/controlled.json。每间隔/幅比仅3次噪声重复，不作为原生/实物双脉冲分辨能力。新原生同窗比较以test-summary.json为准；小样本及高密度偶然匹配限制须保留，不能只凭总波形残差较小认定逐事件恢复正确。

## 重拟合

从项目根目录运行；ROOT、Notebook和artifacts按项目约定仅在本地维护。每系统使用自己的Python环境。下列命令调用冻结源码，输出须为新的项目内目录，保留原始结果：

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  bash tools/run_gagg_python.sh \
  pileup/experiments/alpha_local_model_selection_20261020/replay_fit.py \
  test --workers 2 --output jupyter/artifacts/alpha_local_replay_test
```

同样可用controls阶段重放软件控制。原始生产宏包含当时输出路径，不能直接重复执行覆盖旧ROOT；需在新实验中使用保存的源码/配置重建、分配新的输出文件并记录重放身份。正式结果、图17-1至17-3、限制及验证统一记录于本地codex_0003_Alpha_Pileup_Rate_Fitting.ipynb第17节。未经本次新的明确指示，不上传Notebook、数据或其他修改。


## 本轮完成结果

参数于2026-10-20T17:34:07.113436+08:00冻结；选择倍率1，电荷项系数为0.07574128 mV²/p.e.。23份原生ROOT及两方法拟合全部完成；20/40/100 MHz完整恢复分别为121/291→150/291、90/526→140/526、33/1348→56/1348，10 MHz为90/136→88/136。40/100 MHz额外预测为50→66、43→48；不能只按净恢复数评价。匹配后的电荷误差、同一目标新增/丢失、漏检与拆分/合并标志、错位时间对照保存在本地test-summary.json和report.json。

只读导出的I/O调整见manifest.read_cache_transition：生产ROOT始终直接写在共享g4/data；经60 MHz文件16个数组逐项等价检查后，末段使用校验过SHA256的普通/tmp字节副本读取。原始ROOT、未完成导出的中间文件和日志均保留，拟合/物理/冻结参数没有改变。该调度与I/O改动的实际代码保存在postprocessing/，推断仍使用inference-source/。

五项数值检查、无真值输入冷重放、23个原生ROOT/宏/RunInfo/种子/冻结时序、46组模型重构和重评分、旧Notebook输出保留由本地run_records/alpha_local_model_selection_20261020的验证JSON追溯。绘图/交付脚本为postprocessing/alpha_local_fit_report.py及对应可复用模块；从项目根目录仅重画三张图使用：

```bash
env MPLCONFIGDIR=/tmp/gagg-alpha-local-report \
  bash tools/run_gagg_python.sh -m pileup.analysis.alpha_local_fit_report --figures-only
```

完整的文档交付入口对旧笔记/文档有保护检查，用于本轮一次性交付；已交付后使用--figures-only或--validate-only。重拟合使用前述新输出目录命令。高率每档两窗、局部搜索预算、固定平均形状失配及偶然匹配限制均保留；不据此给出普适时间分辨率、最大计数率、实物或实时吞吐结论。

100 MHz首次触发max_candidates=10000000计算保护，原失败ROOT及日志保留在manifest.failed_attempts；独立alpha_local_capacity_20261020配置快照仅将该保护上限提高至30000000。相同1 MHz种子的16数组逐项一致验证通过后，重试原计划两套100 MHz种子；拟合、物理/读出概率和二进制保持。该配对重放不计入23窗测试分母，三份附加ROOT/宏映射也登记在本地Run_Catalog。

一份100 MHz ROOT的主光子到达时间分支无法解压，原文件保留且该分支仍不可直接读取。评分导出用全局WaveformEvent中的冗余时间，按唯一的x/y/能量浮点位精确一对一连接并核对完整覆盖；另一完整100 MHz文件16个导出数组逐项一致，恢复前后总波形4数组也相同。恢复代码为postprocessing/root_primary_time_recovery.py及recover_primary_times.py，验证在本地redundant-time-validation.json；冻结拟合源码/参数和总电压输入不变。重放拟合直接使用已验证的缓存；重新导出该文件评分字段需启用此适配器。
