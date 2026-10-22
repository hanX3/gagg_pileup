# 单脉冲校准与直接双分量机器学习

从项目根目录运行，使用项目分析Python环境及本地CPU PyTorch。实际源、训练与推断快照分别在single-training-source、pair-training-source、ordered-training-source、inference-source；原G4源码、二进制、宏及种子另存本目录。保留数据、权重和外部依赖，无需切换Git提交。

当前任务为2ns采样、1024点局部0/1/2事件；不是任意高率20µs解码器。三种单脉冲损失各两个种子，无序/有序×直接/物理四种双分量方案各两个种子，共14次训练。有序方法同时加入计数分类器和全局条件解码，不单独归因于分量顺序。固定95/351ns、模拟电压充分动态范围；真值用于训练和评分，盲推断只读观测。Q为雪崩等效pe。

结果和限制见本地Notebook0004第30–34节，结果表与验证位于jupyter/artifacts/alpha_pair_curriculum_20261021。运行过程、备份、原始失败/修复和阶段日志保留在jupyter/run_records/alpha_pair_curriculum_20261021及本目录logs。

单脉冲/无序双分量重训使用replay.py train --variant <已记录方案> --seed <已记录种子> --destination <新模型目录>；有序方案使用replay_ordered.py --variant <pair_ordered_direct或pair_ordered_physical> --seed <已记录种子> --destination <项目内新运行目录>，复制对应冻结单脉冲权重，复用原训练/验证输入。禁止覆盖既有结果。

冻结推断使用replay_predict.py --model <模型名> --waveforms <只含观测电压的NPZ> --output <新目录>，入口自动定位项目与inference-source。NPZ键为voltage_mV、形状N×1024，校准为2ns采样/10mV基线；不接受真实K/t/Q或源率。原始replay.py保留冻结时版本，后补的有序重训/推断便捷入口及哈希记录在本地运行记录extra-entrypoints.json。

24个原生ROOT生成在权重/选择/源码冻结之后；native-records和manifest保存实际宏/三种子/RunInfo核对。native旧视图为真时刻居中诊断，jitter-native补充为按实际源触发时刻随机裁帧，两者共享24母ROOT、未重训或改阈值，不相加成双倍曝光。补充设计在最初推断完成后、审阅其准确度前记录，仍属有外部触发的局部任务。线性控制与联合光子重数字化分开标注；重复案例不等于独立ROOT。known_K_fit_diagnostic仅为已知脉冲数的特定求解器诊断，不是理论上限。普通及严格时间评价均保留。

实际校验器修复了JSON将tuple转list造成的比较错误；原冻结副本保留，实际修正版及检查说明在运行记录validate_actual.py/implementation-notes.json，未修改评分/推断函数。报告绘图为结果显示层，实际版本留在report_actual.py。Notebook只读取保存结果，不会自动触发训练或G4。
