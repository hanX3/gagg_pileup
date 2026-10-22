"""Record the completed experiment without changing historical conclusions."""
import json
from pathlib import Path

project = Path(__file__).resolve().parents[3]
output = project/'jupyter/artifacts/alpha_correlated_global_20261020'
record = project/'jupyter/run_records/alpha_correlated_global_20261020'
experiment = project/'pileup/experiments/alpha_correlated_global_20261020'


def append(path, marker, text):
    before = path.read_text()
    assert marker not in before, f'Already recorded: {marker}'
    path.write_text(before.rstrip()+'\n\n'+text.strip()+'\n')


def main():
    assert (record/'scientific-validation.json').exists() and (record/'notebook-update.json').exists()
    report = json.loads((output/'report-tables.json').read_text())
    diagnostic = json.loads((output/'calibration/residual-diagnostics.json').read_text())
    audit = json.loads((record/'scientific-validation.json').read_text())
    frozen = json.loads((output/'freeze.json').read_text())
    gru = json.loads((output/'gru/frozen-model.json').read_text())
    totals = {r['method']:r for r in report['aggregate']}
    rates = {(r['rate_hz'],r['method']):r for r in report['native']}
    labels = {'baseline':'原始SSE', 'local':'局部K-SSE', 'diagonal':'对角WLS',
              'correlated':'相关GLS', 'global':'GLS＋局部DE', 'cnn':'旧CNN', 'gru':'轻量GRU'}
    table = '| 方法 | 完整恢复/目标 | 漏检（无时间匹配） | 额外输出 | 合并/拆分标志 |\n|---|---:|---:|---:|---:|\n'
    for method in ('baseline','local','diagonal','correlated','global','cnn','gru'):
        row=totals[method]
        table += f"| {labels[method]} | {row['recovered_count']}/{row['target_count']} | {row['missed_count']} | {row['unmatched_predictions']} | {row['merge_flags']}/{row['split_flags']} |\n"
    sentence = f"同一批23个新原生窗中，旧局部K、对角WLS、相关GLS、GLS＋局部DE、旧CNN和轻量GRU的完整恢复数依次为{totals['local']['recovered_count']}、{totals['diagonal']['recovered_count']}、{totals['correlated']['recovered_count']}、{totals['global']['recovered_count']}、{totals['cnn']['recovered_count']}、{totals['gru']['recovered_count']}（共同目标{totals['local']['target_count']}个）。这是非均衡曝光总计，逐率结果及额外输出须同时阅读。"
    chance='；'.join(f"{labels[m]} {100*rates[(100000000,m)]['recovery_fraction']:.2f}% / {100*rates[(100000000,m)]['shifted_fraction_mean']:.2f}%" for m in ('local','correlated','global','cnn','gru'))
    known = []
    for k in (2,3):
        entries = {r['method']:r for r in report['known_count'] if r['known_count']==k}
        a,b=entries['slsqp'],entries['population_plus_slsqp']
        known.append(f"K={k}全部分量完整恢复{a['all_components_complete']}/{a['cases']}→{b['all_components_complete']}/{b['cases']}，评估次数{a['total_evaluations']}→{b['total_evaluations']}")
    append(project/'docs/literature-use.md','id="alpha-correlated-global-gru-20261021"',rf'''
<a id="alpha-correlated-global-gru-20261021"></a>
## L42. 相关残差加权、局部种群搜索与轻量GRU的实际对照（2026-10-21）

**研究问题与采用依据**：实施L41提出的三个方向，沿用其定点原文证据及审核边界，不重复整理精读报告。Cleland/Stern（1994，PDF7–9式21–44、图4；DOI `10.1016/0168-9002(94)91332-3`仅本地题录）启发全协方差权重，其液氩/已知模板结论不等于未知K恢复；Huang等（2017，PDF3、8及原笔记表1–4，DOI `10.1007/s41365-016-0161-z`）启发有预算种群搜索，但原文没有全局最优证书或完整随机流未知K验证；Regadío等（2024，PDF3/7、原笔记PDF8–10，DOI `10.1016/j.apradiso.2023.111137`）启发GRU对照，本轮未复现CR–RC模型或允许标签延迟的MMSE。三篇的完整sources路径和适用条件沿用L41；均保持AI草稿/待人工审核，既有其他论文的审核认可不转移。

**本项目假设和实现**：固定95/351 ns双发光衰减、慢分量0.740174；纯p–¹¹B有效单α谱，无γ。源注入率0.1–100 MHz，不是探测/反应率。原生共享微胞20 μs、2 ns完整模拟电压，中央[5,15) μs评分，动态范围充分，无ADC整窗拒绝。

采用经验$C_{{ij}}=s_i s_j\rho^{{|i-j|}}+\sigma_e^2\delta_{{ij}}$、$s_i^2=v_d+\eta\max(\mu_{{pilot,i}},0)$。C由原始观测pilot整窗固定，不读真值/源率，不随候选K/时间/Q改变，logdet因此为常数。innovation白化后非负幅度变量投影，局部边缘协方差提出方案、全记录固定C接受。K增删代价在单粒子标定重新定标，并非χ²显著性；与旧SSE差异不能只归于C一项。对角WLS为非对角相关项的消融。

48个已知单α标定中33窗有观测候选，得σe=0.1 mV、vd=0 mV²、η=0.7976401192 mV、ρ=0.9819124013。vd落在回归非负边界，原生2 MHz暗计数仍开；残差包含模板失配，非真实微观SiPM协方差。白化后相邻相关{diagnostic['whitened_acf_lag1']:.3f}、二阶矩{diagnostic['whitened_second_moment']:.3f}，保留模型不充分的证据。24个单α及1/10/40 MHz三个开发窗按预定误报/F1规则在0.5/1/2中均选0.5；没有在新测试上重选C或K代价。

局部DE＋SLSQP每窗最多12次DE，预算分布于采集时间；盲流水线实际提出一/二分量，已知K=2/3诊断单列，未声称未知K三分量全局恢复。GRU {gru['parameters']}参数，原CNN的26训练/10验证ROOT、五特征与损失、2400步/batch16；第{gru['step']}步、阈值{gru['threshold']}由旧验证选择，训练{gru['training_seconds']:.1f}秒。架构/种子不同，非等训练墙钟；离线滑动窗口，不宣称实时硬件延迟。

**新计算结果**：共同冻结`{frozen['frozen_at']}`后生成23个新ROOT，与标定/开发、神经训练/验证及已公开测试身份/种子隔离。0.1/1/10/20/40/60/80/100 MHz各8/3/2/2/2/2/2/2窗。目标Q≥500 pe，时间≤10 ns且Q相对误差≤20%才完整恢复。Q为实际祖先归属雪崩等效电荷，不是沉积能量。

{table}

{sentence} 条件Q偏差/标准差、时间误差、逐真值ID新增/损失事件、墙钟在同一本笔记中；并发任务影响墙钟，不作严格硬件吞吐比较。100 MHz实际/五种循环错位平均完整匹配为：{chance}。错位保留输出数量/电荷，不是精确FDR、无偏效率校正或真实性能上限。

216个匹配模板盲软件窗含白/独立相关噪声的空/单/双脉冲，双间隔4/8/12/20/40/80/200 ns、比0.125–4，各2次；没有原生光子/共享微胞涨落，不用于宣称探测器分离极限。36个已知K条件诊断中DE在{report['known_count_objective_improved']}/36降低目标；{'；'.join(known)}。DE额外评估，且K已知，不进入盲效率表。

**记录与验证**：本地 `jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb` 第19节、固定图19-1至19-5，原94单元/图输出不变；原生索引同步本地Run_Catalog。实现位于 `pileup/analysis/alpha_correlated_*.py`，实际源/宏/二进制/标定/权重快照和无需Git切换的重跑入口在 `pileup/experiments/alpha_correlated_global_20261020/`。结果与独立验证分别位于本地 `jupyter/artifacts/alpha_correlated_global_20261020/` 和 `jupyter/run_records/alpha_correlated_global_20261020/`。

**I/O故障单独记录**：19个G4任务全部正常结束后，60 MHz第二个ROOT的全局光子x分支解压失败。由另一份PrimaryPhoton记录按共同时间/y/能量完整值唯一一对一映射16,534,845条记录恢复重复x字段，原ROOT不改写，4个电压数组逐位不变。新增适配只用于评分导出，缺唯一键/双份字段均坏/覆盖不足即拒绝；6项恢复/拒绝测试通过。checked hash加速导出在正常0.1/10/40 MHz的48个数组上与原导出器逐位相同。I/O源码/验证另存io-source及io-recovery.json，不追写原推断冻结；原种子、物理和拟合均不改变，剩余四个预定80/100 MHz由2进程继续、生成后顺序导出。23个预定窗口全部保留，不因算法表现换种子或删测试。

通过6项针对性数值检查、白化快慢核及完整拟合逐位等价；重新计算161份保存评分、{audit['saved_weighted_models_recomputed']}个模型/加权代价和祖先电荷守恒，核验新旧源码SHA、冻结时序/间隔和数据隔离，复读0.1/10/100 MHz实际ROOT，禁用NumPy文件加载的冷盲拟合及冻结源码新路径复跑。原ROOT保留，读取使用SHA相同普通本机字节副本，无临时链接。

**限制与待办**：高率每档仅2短窗、低率目标少；相关模型不足、固定均值模板/局部K预算、电荷选择效应和高密度偶然匹配仍存在。原pilot初始候选>500则退出、局部新增最多160，属于计算预算，不是物理上限；本批初始候选保护触发{report['candidate_limits']['initial_candidate_limit_windows']}/23窗，状态逐窗保留。不宣布最高可用率/通用最小间隔/沉积能谱分辨率；free-τ最短窗及α/γ混合不由本节回答。后续需在新开发协议研究更符合原生涨落的C及未知K候选/联合搜索，再冻结并另作独立测试；这轮测试已公开，不能继续用来选参。NIM、Kim accuracy、Regadío率轴等L41原文缺口仍待文献库独立维护，本轮未改papers。
''')
    append(project/'README.md','alpha_correlated_global_20261020/README.md',f'''
2026-10-21 完成[相关GLS、局部种群搜索与轻量GRU的独立对照](pileup/experiments/alpha_correlated_global_20261020/README.md)：固定双τ和充分动态范围，三种新方法与旧SSE/局部K/CNN在23个冻结后新原生窗同窗比较。{sentence} 协方差白化仍不充分，已知K种群搜索诊断与未知K盲效率分开报告；不构成通用最高率或能量分辨率结论。[L42采用记录与验证边界](docs/literature-use.md#alpha-correlated-global-gru-20261021)对应本地Notebook0003第19节及固定图19-1至19-5，原输出保留。
''')
    append(project/'pileup/README.md','## 相关GLS、局部种群搜索与GRU',f'''
## 相关GLS、局部种群搜索与GRU（2026-10-21）

[加权拟合实现](analysis/alpha_correlated_fit.py)将固定观测pilot的全协方差innovation白化接入局部K/非负幅度变量投影，增加对角消融及每窗12次预算的局部DE＋SLSQP；[轻量GRU](analysis/alpha_correlated_gru.py)使用原CNN的训练/验证ROOT、特征和损失，验证选择权重/阈值。48标定、24单α＋三个率窗开发后共同冻结，另生成0.1–100 MHz八档23个新原生ROOT，完整电压、固定95/351 ns和原生共享微胞条件不变。

{sentence} 完整表同时报告漏检/额外输出/拆并标志、条件Q与时间误差、事件得失和五种循环错位参考。216个盲软件控制、36个已知K=2/3优化诊断单列；条件优化成功不能直接解释为未知K恢复，且DE使用额外评估。

[实际快照及新目录复跑](experiments/alpha_correlated_global_20261020/README.md)保存实际源码、配置、宏、二进制、权重与冻结文件身份，不必切换Git提交；本地Notebook0003追加第19节和图19-1至19-5。[L42](../docs/literature-use.md#alpha-correlated-global-gru-20261021)保留文献审核边界、经验协方差失配和实际验证。Q仍为实现雪崩等效电荷，不是沉积能量；高率仅两窗，不能认证最高可用率或实物分辨率。
''')
    append(experiment/'README.md','## 本轮实际结果',f'''
## 本轮实际结果

{table}

上述非均衡曝光合计只描述本次样本，逐率及条件误差见本地 `report-tables.json`、`test-summary.json`与固定图。100 MHz实际/五种错位均值为：{chance}。不能以低残差、条件Q标准差或偶然时间匹配代替完整事件效率。

已知K软件诊断：{'；'.join(known)}。与盲原生结果分开，DE额外评估。白化后相邻残差相关{diagnostic['whitened_acf_lag1']:.3f}，经验C不足，不根据这批测试继续调参。

验证重新计算161份保存评分、{audit['saved_weighted_models_recomputed']}个加权模型/目标、23窗祖先电荷守恒，复读三档ROOT、核对新旧冻结源码和冻结前后隔离；禁止NumPy文件加载的冷盲拟合及本目录冻结源码新路径复跑均通过。Notebook只执行本次追加单元，历史94单元与输出不变。ROOT—实际宏索引与RunInfo核验随本地Run_Catalog保存。
''')
    print('DOCS_UPDATED_L42')


if __name__=='__main__':
    main()
