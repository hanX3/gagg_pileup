# 纯α双事件条件分辨曲线

本实验将“多近不能解堆积”改为指定Q、噪声、误报预算和恢复容差下的条件曲线。固定95/351 ns双发光衰减，2 ns采样，2.048 μs局部模拟电压；充分动态范围，不按ADC触顶剔除。Q为实际雪崩等效电荷，不是沉积能量。受控单/双事件不对应Poisson源率。

`protocol.json`在新测试前冻结：17个间隔（含严格同刻控制）、7组均值条件、5组共同微胞条件，每个条件/间隔40对。240个新单α光学母ROOT，源能量2.5/5.486/8 MeV各80个；同一母对跨间隔复用，不能算作独立曝光。五组共同微胞对由新光子输运离线联合数字化，不是新增原生双事件ROOT。另比较独立原生分量的线性叠加和同一共同微胞实现的1 mV噪声控制。

盲方法为原变投影K=0/1/2拟合、两套已训练有序物理头、两套ML＋局部拟合。本轮不重训、不重选阈值。观察文件仅含`voltage_mV`；真实K/t只在单独的诊断阶段打开。已知K拟合是特定求解器诊断，不是理论上界。Fisher/GLS信息只在准确均值与已知电子学噪声假设下成立。

严格双事件成功：各Q相对误差≤20%，各时间误差≤min(10 ns, Δ/4)，两个真事件一对一匹配，无未匹配输出。Δ=0不评个体分离率。Δ90为从该点到所有之后已测正间隔点均达到90%的最小网格点；若不存在不报告阈值。空和相应单事件控制各自要求≤5%帧有未匹配输出；逐点Wilson95区间另外报告，非同时置信带。网络训练间隔4–200 ns，2/400/800 ns为外推。

## 冻结与文件

- `source/`、`bin/`：实际G4源码/运行二进制。G4内部串行，4个独立进程；ROOT启动至少相隔1秒，每ROOT单个20 μs WaveformEvent。
- `replay-source/`：实际数据、推断和继承算法源码；`models/`为两套权重/配置/阈值；`freeze.json`记录哈希。
- `manifest.json`、`macros/`、`logs/`：实际ROOT路径、三类种子、宏和RunInfo校验。
- `run.py`：实际可用入口。原`replay.py`的spawn命令分派问题由该入口修复；原冻结文件保留。第一次均值推断失败时尚无预测，失败日志保留在`logs/mean-predict.log`。
- `recovery-tools/recover_native.py`：本轮共享盘ROOT短写后的恢复程序。原240个任务均已运行；其中1个文件压缩流截断，保留原文件，以相同种子/物理配置串行写入新ROOT。另重新核对并登记了原异常之后未被收集器写入清单的有效输出，最终240个有效母样本，未重选测试种子。
- `recovery-tools/physical_cached.py`：实际离线数字化适配器，本地普通文件缓存、逐例检查点及4个工作进程；原缓存准备后长期无进度的阶段以130中断，原因未确定，尚无正式波形数据集落盘。适配器经历2进程执行检查后恢复4进程，所有完成的检查点保留。
- `optimization-source/`：独立光子RNG批量生成、PDE预插值和等价二分查找；不改变模型参数或候选顺序。运行记录中的六组比较逐项验证雪崩/微胞末态/波形完全一致，CPU时间合计约加快1.50倍；该小样本计时不是在线吞吐测量。原冻结文件未覆盖，实际适配器显式记录所用优化哈希。
- `postprocessing-source/`、`report.py`：完成后保存的实际报告源码及入口，独立于冻结推断；不会覆盖原冻结文件。
- 分析数组/预测/图：项目本地`jupyter/artifacts/alpha_pair_resolution_20261022/`。
- 验证、环境、启动器旧版本和Notebook备份：项目本地`jupyter/run_records/alpha_pair_resolution_20261022/`。
- 研究记录：本地Notebook0004第36–40节；ROOT完整映射在本地Run_Catalog。项目采用说明在`docs/literature-use.md` L49。

## 重放

以下在项目根目录运行，需要本批记录的G4数据路径和Python依赖（NumPy、SciPy、uproot、awkward、PyTorch等）。实际环境及版本保存于运行记录；模型和输入均须保留。脚本不自动下载依赖。`PYTHONPATH=/tmp/gagg-correlated-ml-deps-20261020`仅为本机本次PyTorch依赖位置，跨系统应使用自己的等效安装。

```bash
python pileup/experiments/alpha_pair_resolution_20261022/run.py mean
python pileup/experiments/alpha_pair_resolution_20261022/run.py native --workers 4
python pileup/experiments/alpha_pair_resolution_20261022/run.py physical --workers 4
python pileup/experiments/alpha_pair_resolution_20261022/run.py linear

python pileup/experiments/alpha_pair_resolution_20261022/run.py fit predict --dataset mean --workers 4
python pileup/experiments/alpha_pair_resolution_20261022/run.py fit predict --dataset shared --workers 4
python pileup/experiments/alpha_pair_resolution_20261022/run.py fit predict --dataset linear --workers 4
python pileup/experiments/alpha_pair_resolution_20261022/run.py fit predict --dataset shared-noise --workers 4

python pileup/experiments/alpha_pair_resolution_20261022/run.py fit diagnostics --dataset mean
python pileup/experiments/alpha_pair_resolution_20261022/run.py fit diagnostics --dataset shared
python pileup/experiments/alpha_pair_resolution_20261022/run.py fit diagnostics --dataset linear
python pileup/experiments/alpha_pair_resolution_20261022/run.py fit diagnostics --dataset shared-noise
python pileup/experiments/alpha_pair_resolution_20261022/report.py all
```

本次实际数字化使用`python pileup/experiments/alpha_pair_resolution_20261022/recovery-tools/physical_cached.py`替代上面的参考`physical`命令，依赖保存的优化一致性检查和本地缓存记录。`post_physical.py`在该阶段成功后接续线性控制、三域盲推断及单独诊断，同时最多6个拟合工作进程。运行记录保存每个后续阶段的实际命令、开始/结束时间和退出码。

数据生成命令拒绝覆盖既有资产；已有预测可按原worker数续跑，诊断拒绝覆盖。仅重画图用`report.py figures`；重算评分和图用`report.py all`。不同实验须复制源码并建立新输出路径，不删除或静默覆盖这批有效数据。所有后续重放直接调用实验快照，不需要切换Git提交。

均值信息图的数值差分敏感性检查在报告阶段完成：实际插值模板导数与冻结网格梯度近似在极短间隔会有可见差异，两者都保留；不改变任何盲预测或评分。线性叠加与共同微胞使用相同光学母样本，但微胞实现不同，不能把差值全部解释为占用；加噪声版本才是同一微胞实现的电子学噪声消融。

验证结果以`science-validation.json`和`notebook-validation.json`为准。新报告单元只加载保存结果，不启动训练或G4。原Notebook计算及图像均保留；本次新图固定编号0004-37-1、37-2、38-1、38-2、38-3、38-4、39-1。

首次Notebook写入进程以143结束，主笔记出现零字节文件；已从运行前备份逐字节恢复，Run_Catalog未受影响，信号发送原因未确定。更新器随后改为临时文件完整写入、fsync、字节与nbformat核验后原子替换；恢复记录见`notebook-write-recovery.json`，最终历史单元/图像保留情况由`notebook-validation.json`验证。此中断不涉及原始数据或拟合结果。

## 实际结果

中能共同微胞条件，盲拟合与两套ML＋拟合在40 ns均29/40严格成功，60–800 ns各已测点均40/40，且相关空/单控制通过，合格Δ90为60 ns。均值参考Q4000/4000 pe为盲拟合16 ns、混合方法6 ns，两域电荷/响应不同，不能据此单独估计某项物理效应。

原始网络在中能8 ns分别17/40、35/40，接固定模板拟合均降至0/40；物理重拟合并非必然改善。盲拟合高Q单事件额外输出46/80，含高Q条件未通过误报预算；网络/混合各空/单控制0/80。新图38-4直接比较原始ML与混合输出，未修改评分或阈值。完整条件阈值、独立性和置信限制见本地报告；本实验不支持一个普适的正数不可分间隔。
