"""Append the completed section 19 and ROOT catalog, preserving prior cells.

Execute only the appended cells in a fresh Python namespace and embed their
actual display outputs. This does not rerun the historical Notebook cells.
"""
from datetime import datetime
import json
import os
from pathlib import Path
import sys

import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output

project = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project))
from pileup.analysis.alpha_rate_campaign import save

record = project/'jupyter/run_records/alpha_correlated_global_20261020'
output = project/'jupyter/artifacts/alpha_correlated_global_20261020'
experiment = project/'pileup/experiments/alpha_correlated_global_20261020'


def read(name):
    return json.loads((output/name).read_text())


def markdown(text):
    return nbformat.v4.new_markdown_cell(text.strip()+'\n')


def code(text):
    return nbformat.v4.new_code_cell(text.strip()+'\n')


def fraction(row):
    return f"{row['recovered_count']}/{row['target_count']} ({100*row['recovery_fraction']:.1f}%)"


def execute_appended(notebook, start, label):
    os.environ.setdefault('IPYTHONDIR', '/tmp/gagg-correlated-ipython')
    Path(os.environ['IPYTHONDIR']).mkdir(parents=True, exist_ok=True)
    InteractiveShell.instance()
    namespace = {'__name__': '__main__'}
    counts = [c.get('execution_count') or 0 for c in notebook.cells[:start] if c.cell_type == 'code']
    count = max(counts, default=0)
    executed = []
    for index, cell in enumerate(notebook.cells[start:], start=start):
        if cell.cell_type != 'code':
            continue
        with capture_output(stdout=True, stderr=True, display=True) as captured:
            exec(compile(cell.source, f'{label}:cell-{index}', 'exec'), namespace)
        outputs = []
        if captured.stdout:
            outputs.append(nbformat.v4.new_output('stream', name='stdout', text=captured.stdout))
        if captured.stderr:
            outputs.append(nbformat.v4.new_output('stream', name='stderr', text=captured.stderr))
        for item in captured.outputs:
            outputs.append(nbformat.v4.new_output('display_data', data=item.data, metadata=item.metadata))
        count += 1
        cell.execution_count = count
        cell.outputs = outputs
        executed.append(dict(index=index, outputs=len(outputs), no_error=True))
    return executed


def main():
    notebook_path = project/'jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb'
    before_path = record/'before__jupyter__codex_0003_Alpha_Pileup_Rate_Fitting.ipynb.before'
    notebook = nbformat.read(notebook_path, as_version=4)
    original = nbformat.read(before_path, as_version=4)
    assert notebook.cells == original.cells, 'Preserve existing work; do not overwrite an intervening edit'
    assert len(notebook.cells) == 94
    report = read('report-tables.json'); covariance = read('calibration/covariance.json')
    diagnostic = read('calibration/residual-diagnostics.json')
    validation = json.loads((record/'scientific-validation.json').read_text())
    budget_diagnostic = read('optimizer-budget-controls.json')
    f = read('freeze.json'); gru = read('gru/frozen-model.json')
    method_names = ('baseline', 'local', 'diagonal', 'correlated', 'global', 'cnn', 'gru')
    labels = dict(zip(method_names, ('原始SSE', '局部K-SSE', '对角WLS', '相关GLS', 'GLS＋局部DE', '旧CNN', '轻量GRU')))
    totals = {r['method']: r for r in report['aggregate']}
    native = {(r['rate_hz'], r['method']): r for r in report['native']}
    methods_table = '| 方法 | 残差/搜索 | 模型选择或监督 |\n|---|---|---|\n'
    methods_table += '| 原始SSE | 原始候选检测与无权最小二乘 | 原冻结配置 |\n'
    methods_table += '| 局部K-SSE | 第17节局部拆分/合并、增删K；无权残差 | 原标定/开发集选定配置 |\n'
    methods_table += '| 对角WLS | 同一观测pilot给方差；不含非对角项 | 本轮标定K代价，开发集选择倍率 |\n'
    methods_table += '| 相关GLS | 全协方差白化＋非负幅度变量投影 | 同上；协方差逐窗固定 |\n'
    methods_table += '| GLS＋局部DE | GLS基础上增加有预算的差分进化，再SLSQP | 每窗最多12次DE；其余配置冻结 |\n'
    methods_table += '| 旧CNN | 保留原冻结权重/阈值 | 原26训练＋10验证ROOT，未重新训练 |\n'
    methods_table += '| 轻量GRU | 8 ns recurrent token＋2 ns局部卷积头 | 同26/10 ROOT、特征、损失、2400步；独立种子 |\n'
    result_table = '| 源注入率 (MHz) | 目标数 | '+' | '.join(labels[m] for m in method_names)+' |\n'
    result_table += '|---:|---:|'+'---:|'*7+'\n'
    for rate in sorted({r['rate_hz'] for r in report['native']}):
        result_table += f"| {rate/1e6:g} | {native[(rate,'local')]['target_count']} | "+' | '.join(fraction(native[(rate,m)]) for m in method_names)+' |\n'
    aggregate_table = '| 方法 | 完整恢复/目标 | 无时间匹配漏检 | 额外输出 | 合并标志 | 拆分标志 | 相对局部K新增/损失的完整事件 |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for m in method_names:
        r = totals[m]
        aggregate_table += f"| {labels[m]} | {fraction(r)} | {r['missed_count']} | {r['unmatched_predictions']} | {r['merge_flags']} | {r['split_flags']} | +{r['gained_vs_local']} / −{r['lost_vs_local']} |\n"
    chance = '| 方法 | 100 MHz实际完整匹配 | 五种错位均值 |\n|---|---:|---:|\n'
    for m in method_names:
        r = native[(100000000,m)]
        chance += f"| {labels[m]} | {100*r['recovery_fraction']:.2f}% | {100*r['shifted_fraction_mean']:.2f}% |\n"
    conditional = '| 已知K | 局部搜索：全部分量恢复 | DE＋局部：全部分量恢复 | 函数评估次数：局部 / DE＋局部 |\n|---:|---:|---:|---:|\n'
    for k in (2,3):
        a = next(r for r in report['known_count'] if r['known_count'] == k and r['method'] == 'slsqp')
        b = next(r for r in report['known_count'] if r['known_count'] == k and r['method'] == 'population_plus_slsqp')
        conditional += f"| {k} | {a['all_components_complete']}/{a['cases']} | {b['all_components_complete']}/{b['cases']} | {a['total_evaluations']} / {b['total_evaluations']} |\n"
    cells = [markdown(r'''
<a id="correlated-global-gru-20261021"></a>
## 19. 相关噪声加权、局部种群搜索与GRU：实现及独立测试（2026-10-21）

**研究问题**：第18节文献建议能否增加纯α盲pile-up的逐事件恢复？本节实际实现并冻结三条路线，再生成新原生波形测试，不用此前已看过的测试窗重新宣称独立验证。

仍使用p–¹¹B有效单α源谱、GAGG–SiPM共享微胞响应，**没有γ源**。双发光衰减固定95/351 ns，慢分量0.740174；不是free-τ或短窗研究。使用完整20 μs模拟电压、2 ns采样，中央[5,15) μs评分；动态范围充分，不按ADC触顶拒绝。

Q表示实际共享微胞历史下归属该初级粒子的**雪崩等效电荷（pe）**，不是已标定沉积能量。源注入率也不等于探测事件率或物理反应率。真值仅用于监督训练、构造软件对照和推断后的评价。
'''), markdown('### 19.1 七种方法及真正改变的部分\n\n'+methods_table+r'''

相关GLS（广义最小二乘）改了目标函数。定义观测pilot均值$\mu_0$，固定
$$ C_{ij}=s_i s_j\rho^{|i-j|}+\sigma_e^2\delta_{ij},\qquad s_i^2=v_d+\eta\max(\mu_{0,i},0). $$
用状态空间innovation递推计算$W$，满足$WCW^T=I$；固定候选时刻时解$\min_{A\ge0}\|W(H A-(y-b))\|^2$，再优化到达时刻。每个记录的C在初始观测pilot之后固定，比较不同K/时刻/幅度时不随候选重算；$\log\det C$因此为常数。拟合不读取真值、源率或真实脉冲数。

局部提出方案使用最长1400 ns尾部的边缘协方差，**最终接受用完整记录的同一个C**。增加/删除分量按完整加权残差增益与经验代价$f(1+aQ_{ref})\Delta K$比较；该代价由孤立标定控制，不是χ²显著性检验，也不声称求得全局MAP。候选非负幅度可低至100 pe作干扰分量，最终仅报告≥500 pe。

本轮同时重新标定了K的增删代价。与旧SSE的差异不能全部归因于协方差单项；对角WLS与相关GLS是更直接的非对角协方差对照。两者仍为有预算的局部K启发式搜索。

差分进化（DE）搜索局部归一化时间盒，在SLSQP前额外探索。12次预算按采集时间分布，不能保证最优；盲原生流水线实际提出的是局部一/二分量模型。三分量搜索只在19.4明确已知K的诊断中实施，不能写成已实现未知K三脉冲全局恢复。

原观测pilot仍保留“初始候选>500则退出”的计算保护，局部方法每窗最多接受160次新增。这些是本实现的预算限制，不是探测器物理上限；若pilot失败，依赖它的物理改进也返回初始化失败，不从真值补候选。
'''), markdown(f'''
### 19.2 标定、开发和冻结边界

沿用独立48个已知单α标定窗，其中33窗有观测候选，可用于残差矩估计；另24个单α开发窗及1/10/40 MHz三个开发窗选择K代价倍率。三种新物理方法都从0.5/1/2中选得0.5，先最少孤立额外输出，再比较开发集宏平均完整恢复F1。旧局部K配置保持第17节冻结值。所有新测试种子与这些数据及神经网络训练/验证ROOT分开。

经验协方差为σe=0.1 mV、vd={covariance['parameters']['dark_variance_mV2']:.6g} mV²、η={covariance['parameters']['shot_coefficient_mV']:.6g} mV、ρ={covariance['parameters']['rho']:.6g}。**vd=0是回归落在非负边界，不是将原生2 MHz暗计数关闭**。该残差包含光子/雪崩涨落及模板失配；AR(1)和线性方差是本项目假设，不是文献证明的GAGG真实协方差。

白化后标定残差仍有相邻相关{diagnostic['whitened_acf_lag1']:.3f}，二阶矩{diagnostic['whitened_second_moment']:.3f}，应与理想0和1比较。图19-1把经验相关与含电子白噪声的模型预测同时显示；本轮保留失配，不根据新测试再改C。

GRU有{gru['parameters']}个参数，CPU训练{gru['training_seconds']:.1f}秒，验证集选择第{gru['step']}步与阈值{gru['threshold']:.2f}。训练使用相同特征和损失、2400步、batch16，但架构/随机种子不同，不是等训练墙钟比较，也未复现文献的目标延迟MMSE。滑动1024点，带保护的500点核心逐段输出，属于离线推断。

原数据划分在100 MHz仅有1个训练ROOT，没有同率验证ROOT，80 MHz是最高验证率。因此100 MHz的神经网络结果属于该划分下的外推测试，不把新测试拿来重新选择阈值。

共同冻结时间：`{f['frozen_at']}`。随后生成23个新ROOT，0.1/1/10/20/40/60/80/100 MHz分别8/3/2/2/2/2/2/2窗。先通过1/10 MHz两个原生预检查再运行其余；每ROOT单个WaveformEvent、G4内部串行，最多4个G4进程，分析与生产合计不超过8核。
'''), markdown('''
<a id="codex-figure-19-1"></a>
**图 19-1（0003-19-1）：单α标定残差与经验协方差的适用程度。** (a) 残差二阶矩与均值；(b) 方差标准化后的经验相关及模型预测；(c) 完整白化后剩余相关。这里只看标定数据，不把残差协方差当作已验证的微观噪声模型。
'''), code('''
from pathlib import Path
import sys
from IPython.display import display, Image
correlated_project = next(p for p in (Path.cwd(), *Path.cwd().parents) if (p/'pileup/analysis/alpha_correlated_fit.py').exists())
if str(correlated_project) not in sys.path:
    sys.path.insert(0, str(correlated_project))
from pileup.analysis.alpha_correlated_report import calibration_figure
correlated_artifacts = correlated_project/'jupyter/artifacts/alpha_correlated_global_20261020'
calibration_figure()
display(Image(filename=str(correlated_artifacts/'covariance_calibration.png')))
'''), markdown('''
### 19.3 新原生数据：盲重建结果

“完整恢复”同时要求一对一时间匹配≤10 ns、相对电荷误差≤20%；目标为中央注入且实际Q≥500 pe的α。时间匹配不按电荷选配；额外输出是中央未匹配预测，漏检是没有时间匹配的目标。已匹配但电荷失败另列。合并/拆分只是接近关系的操作标志，高密度下不作完整因果归因。

'''+result_table+'\n以下总计仅描述这组非均衡曝光，不把所有事件视为独立重复来产生窄置信区间。\n\n'+aggregate_table), markdown('''
<a id="codex-figure-19-2"></a>
**图 19-2（0003-19-2）：七种冻结方法在新原生波形上的恢复与失败。** (a) 完整恢复率；(b) 额外输出占中央输出的比例；(c) 无时间匹配的漏检比例；(d) 匹配后电荷失败比例。原始计数/分母保存在上表和report-tables.json。
'''), code('''
import json
from pileup.analysis.alpha_correlated_report import native_figures
correlated_report = json.loads((correlated_artifacts/'report-tables.json').read_text())
native_figures(correlated_report)
display(Image(filename=str(correlated_artifacts/'native_recovery.png')))
'''), markdown('''
<a id="codex-figure-19-3"></a>
**图 19-3（0003-19-3）：条件误差、事件得失、计算成本与错位对照。** (a) 时间匹配后的电荷误差中位数；(b) 电荷误差样本标准差；(c) 时间误差样本标准差；(d) 每20 μs记录的推断墙钟；(e) 相对旧局部K按同一真值ID计算新增/损失的完整事件；(f) 100 MHz实际与五种循环错位的平均匹配率。条件误差有选择效应，不是全部入射α的能量分辨率。
'''), code('''
display(Image(filename=str(correlated_artifacts/'native_diagnostics.png')))
'''), markdown('''
同一20 μs完整记录的物理方法计时包含观测pilot，几种改进共享其结果；CNN/GRU从电压特征开始计时。运行期间存在其他生产/分析任务，墙钟只表示本次离线实现的量级，不是严格隔离的硬件吞吐或实时延迟基准。

循环错位预先取1379/2933/4871/7039/8911 ns，保留中央输出数量和电荷、保护区预测不动，仅作为偶然匹配参考；不用于调参，不直接相减作真实效率。100 MHz下：

'''+chance), markdown('''
### 19.4 受控双脉冲与已知K的优化诊断

另做216个盲软件对照：白噪声/独立生成的相关噪声各108窗，每种含12空窗、12孤立、84双脉冲；双脉冲间隔4/8/12/20/40/80/200 ns，Q1=4000 pe，Q2/Q1=0.125/0.25/0.5/1/2/4，每组2次。均值模板匹配，噪声独立采样；没有原生光子数、共享微胞饱和或真实能谱，不能给出探测器时间分辨极限。

此外做36个**已知K=2/3**的条件诊断。两方法使用相同固定940–1200 ns触发区域和三组起点，不输入真实时刻/幅度；DE额外花费种群评估，因而不视为等计算预算对照。目标与残差噪声同上。

'''+conditional+f"\nDE在{report['known_count_objective_improved']}/36案例找到严格更低目标值，但是否恢复每个分量仍由时间和电荷共同评价。已知K的成功不计入19.3的盲重建效率。\n"), markdown('''
<a id="codex-figure-19-4"></a>
**图 19-4（0003-19-4）：软件诊断区分脉冲发现与局部优化。** (a)(b) 未提供K/时刻的双脉冲恢复；(c) 已知K=2/3时两个优化器的加权残差；(d) 已知K时全部分量同时恢复的案例数。下排条件更强，不能和上排或原生结果混为同一任务。
'''), code('''
from pileup.analysis.alpha_correlated_report import controls_figure
correlated_controls = json.loads((correlated_artifacts/'controls-summary.json').read_text())
correlated_optimizer = json.loads((correlated_artifacts/'optimizer-controls.json').read_text())
controls_figure(correlated_controls, correlated_optimizer)
display(Image(filename=str(correlated_artifacts/'software_controls.png')))
'''), markdown('''
<a id="codex-figure-19-5"></a>
**图 19-5（0003-19-5）：预定10/40 MHz首个测试窗的固定5.0–5.8 μs片段。** 左列为观测及四种物理拟合的总电压，右列为七种方法逐事件电荷；×为拟合结束后才显示的真值。按预先种子顺序选窗，不按重建好坏选图。总波形贴合并不保证其分解正确。
'''), code('''
from pileup.analysis.alpha_correlated_report import examples_figure
correlated_test = json.loads((correlated_artifacts/'test-summary.json').read_text())
examples_figure(correlated_test)
display(Image(filename=str(correlated_artifacts/'native_examples.png')))
'''), markdown('''
### 19.5 复现、验证与限制

可复用实现：`pileup/analysis/alpha_correlated_fit.py`、`alpha_correlated_gru.py`、`alpha_correlated_evaluate.py`及`alpha_correlated_report.py`。协议/标定/开发选择/权重/源码的实际副本随`pileup/experiments/alpha_correlated_global_20261020`保存；冻结副本独立复跑，无须切换Git提交。新绘图单元只读保存的结果，不依赖PyTorch；训练/神经推断才需要记录版本的PyTorch。

从项目根目录运行：
```bash
bash jupyter/run_gagg_python.sh -m pileup.analysis.alpha_correlated_report
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_correlated_global_20261020/replay_analysis.py one --output jupyter/artifacts/alpha_correlated_replay_new
```
输出目录必须新建，避免覆盖原结果。ROOT、实际宏、三个种子、RunInfo、配置和输出映射另保存在本次root-macro-index.json。共享盘读取使用经过SHA核对的普通本机字节副本，没有临时链接；原ROOT保留在g4/data。

数值验证覆盖innovation与稠密Cholesky、对角极限、独立噪声样本、变量投影梯度、时间平移和空窗。加速白化核与原核在矩阵和完整拟合上逐位相同。161份保存评分重新从预测/真值计算，69份加权模型重新构造并核对残差目标；检查祖先电荷守恒、旧冻结源码不变、新测试冻结时序/身份隔离，并在禁用NumPy文件读取时冷启动盲拟合。直接复读预定0.1/10/100 MHz ROOT核验电压和RunInfo。实际验证详情随运行记录保存。

高率每档仅2个短窗；协方差近似、固定均值模板和局部K预算均可能限制恢复。低率事件数也少，不能宣布最高可用率、通用最小间隔或实物能量分辨率。最短free-τ窗口及α/γ混合属于其他研究问题，不从本节结果推断。

文献采用边界延续第18节：Cleland/Stern启发相关噪声权重，但没有未知K多事件恢复；Huang启发种群搜索，但没有全局最优保证；Regadío启发GRU架构对照，但本实现不是其CR–RC数据和延迟MMSE的复现。新增文献仍保留AI草稿/待人工审核状态。
''')]
    findings = '### 19.6 本轮实际判断\n\n'
    findings += f"同一批新原生数据中，旧局部K完整恢复{fraction(totals['local'])}，对角WLS为{fraction(totals['diagonal'])}，相关GLS为{fraction(totals['correlated'])}，GLS＋局部DE为{fraction(totals['global'])}。这些总数只描述本次不均衡曝光，逐率判断以上表为准。\n\n"
    findings += f"加入预算DE相对不含DE的GLS，完整事件数变化为{totals['global']['recovered_count']-totals['correlated']['recovered_count']:+d}，额外输出变化为{totals['global']['unmatched_predictions']-totals['correlated']['unmatched_predictions']:+d}；它的已知K受控优势不能替代未知K时候选发现、拆并选择和模型失配的处理。\n\n"
    findings += f"旧CNN为{fraction(totals['cnn'])}，轻量GRU为{fraction(totals['gru'])}，额外输出分别{totals['cnn']['unmatched_predictions']}和{totals['gru']['unmatched_predictions']}。本轮只考察这一种轻量GRU和既定训练预算；不能外推为GRU架构的一般性能排序。\n\n"
    findings += '这轮把算法改动实际落地，并给出了失败位置：协方差白化不充分；已知K时改善优化目标仍不保证未知K的完整恢复；高率还有电荷分解与偶然时间配对。后续应针对这些失配和未知K候选问题设计新的开发/测试协议，不能继续用本轮测试集挑超参数。\n'
    cells.append(markdown(findings))
    budget_table='| 已知K | 每方法运行数 | 多随机初值SLSQP：全部分量合格 | DE混合：全部分量合格 | 每方法总profile评估 |\n|---:|---:|---:|---:|---:|\n'
    for known in (2,3):
        a=next(r for r in budget_diagnostic['summary'] if r['known_count']==known and r['strategy']=='random_multistart_slsqp')
        b=next(r for r in budget_diagnostic['summary'] if r['known_count']==known and r['strategy']=='de_then_slsqp_and_restart_fill')
        assert a['total_profile_evaluations']==b['total_profile_evaluations']==500*a['runs']
        budget_table+=f"| {known} | {a['runs']}（18噪声×3种子） | {a['all_components_complete']} | {b['all_components_complete']} | {a['total_profile_evaluations']} |\n"
    budget_text='''### 19.4.1 补充：相同profile评估预算、三个优化随机种子

前面的36个案例使用三组固定粗初值，DE额外计算；该诊断不能证明同算力下优于充分多初值的局部优化。为此在看过该初步结果后追加本检查，先另存协议/源码SHA，再生成新的合成噪声；不是原生盲测试的重新选参。

新36个已知K=2/3案例采用相同模板/噪声、940–1200 ns触发区域、间隔和电荷比，每个噪声案例配3个优化种子。两策略均包含同三组粗初值，并各消耗**恰好500次实际白化NNLS变量投影评估**：一方用随机多初值SLSQP；另一方将公共初值之后剩余预算的70%用于DE，再局部抛光并以重启填满预算。所有已评估的合法候选取最小值；不输入真时刻/电荷。三种子共享同一噪声，108次运行不是108个独立噪声样本。

'''+budget_table+'''
预算包括公共初值与抛光中的实际profile计算，另记无效时间向量开销；不是严格相等墙钟。表中“合格”是已知K时的时间±10 ns及Q±20%条件精度，不是未知K检出/报告效率。原生七方法、阈值和全部23个种子保持原冻结；补充记录在optimizer-budget-controls.json，图19-4仍展示前一组粗初值/追加DE案例。
'''
    position=next(i for i,c in enumerate(cells) if c.cell_type=='markdown' and 'id="codex-figure-19-5"' in c.source)
    cells.insert(position,markdown(budget_text))
    candidate_failures = report['candidate_limits']['initial_candidate_limit_windows']
    cells.append(markdown(f'本批原生测试中，500个初始候选的计算保护实际触发{candidate_failures}/23窗。各率/方法的拟合状态及初始候选数保存在report-tables.json；解释高率结果时须同时考虑这个状态和新增预算。'))
    for cell in cells:
        if cell.cell_type == 'markdown':
            cell.source = cell.source.replace('69份加权模型', f"{validation['saved_weighted_models_recomputed']}份加权模型")
    cells.append(markdown('''
#### 本轮数据I/O故障及精确恢复边界

19个G4任务正常结束后，60 MHz第二个ROOT的`WaveformEvent/all_photon_arrival_x_mm`解压失败，电压字段仍可读取。使用另一份PrimaryPhoton记录，按时间、y坐标和能量的完整存储值建立16,534,845条记录的唯一一对一映射，取回重复保存的x值；不猜测坐标、祖先、脉冲时刻或电荷，原ROOT不改写。恢复前后4个电压输入数组逐位相同。

新增只读I/O适配`root_redundant_photon_fields.py`在缺少唯一公共键、任一字段双份均坏、全记录不等或覆盖不完整时拒绝恢复。6项数值/拒绝测试通过；正常0.1/10/40 MHz ROOT的48个导出数组与原导出器逐位相同。加速路径只以hash定位，再检查完整字段相等；不是用hash碰巧相等当证据。

这项I/O适配在解码故障后补充，单独记录源码SHA与验证，不改变此前冻结的拟合/网络、物理、原定种子或评价条件。剩余四个预定80/100 MHz种子按原宏/二进制继续，以2个G4进程生成完再顺序导出，减少内存争用。所有原件/失败记录保留；全部23个预定窗口进入评价，没有因重建好坏筛掉测试。实际故障与恢复字段随manifest/metadata保存。
'''))
    start = len(notebook.cells); notebook.cells.extend(cells)
    executions = execute_appended(notebook, start, 'Notebook0003-section19')
    assert notebook.cells[:start] == original.cells
    nbformat.validate(notebook); nbformat.write(notebook, notebook_path)

    manifest = json.loads((experiment/'manifest.json').read_text())
    catalog_rows = []
    for row in sorted(manifest['runs'], key=lambda r: (r['rate_hz'], r['seed'])):
        stem = Path(row['root']).stem
        meta = json.loads((project/row['export_metadata']).read_text())
        catalog_rows.append({**row, 'window_ns':meta['duration_ns'], 'gate_ns':meta['gate_ns'],
            'run_info_base_seeds':next(v['run_info_base_seeds'] for v in validation['rows'] if v['root']==row['root']),
            'source_rates_hz':meta['source_rates_hz'], 'primaries':meta['primaries'],
            'waveform':str((output/'inputs/waveforms'/(stem+'.npz')).relative_to(project)),
            'truth':str((output/'inputs/truth'/(stem+'.npz')).relative_to(project)),
            'physical_model':'Unchanged native p-11B effective single-alpha spectrum and shared-cell SiPM',
            'fit_input':'analog_mV; adequate dynamic range; no rail rejection'})
    save(output/'root-macro-index.json', dict(rows=catalog_rows, inference_frozen_at=f['frozen_at'],
        source_snapshot=str((experiment/'source').relative_to(project)), binary_sha256=manifest['binary_sha256']))
    catalog_path = project/'jupyter/Run_Catalog.ipynb'
    catalog = nbformat.read(catalog_path, as_version=4)
    catalog_before = nbformat.read(record/'before__jupyter__Run_Catalog.ipynb.before', as_version=4)
    assert catalog.cells == catalog_before.cells
    catalog_start = len(catalog.cells)
    catalog.cells.extend([markdown(f'''
<a id="alpha-correlated-global-20261021"></a>
## 2026-10-21：相关GLS／局部种群搜索／GRU的23个独立纯α测试窗

共同冻结`{f['frozen_at']}`后生产；源注入率0.1/1/10/20/40/60/80/100 MHz，分别8/3/2/2/2/2/2/2窗。源为p–¹¹B有效单α谱，无γ背景；20 μs、2 ns采样、中央[5,15) μs评价。双发光衰减95/351 ns及原生共享微胞读出保持；离线使用完整模拟电压，动态范围充分，不按ADC触顶剔除。

每ROOT单个WaveformEvent，G4串行、最多4个独立进程，启动至少间隔1秒；1/10 MHz预检查包含在23窗内。源率不是探测事件率或反应率。新ROOT仍保存于g4/data，实际源码/配置/二进制、宏、日志、三个种子及冻结源码/权重副本保存在实验目录。共享盘读取通过SHA相同的普通本机字节副本，没有临时链接。

[完整ROOT—实际宏索引](artifacts/alpha_correlated_global_20261020/root-macro-index.json) · [运行清单](../pileup/experiments/alpha_correlated_global_20261020/manifest.json) · [冻结清单](artifacts/alpha_correlated_global_20261020/freeze.json) · [161份重评分与原生一致性验证](run_records/alpha_correlated_global_20261020/scientific-validation.json)

216个盲软件空/孤立/双脉冲及36个已知K初步诊断不是G4 ROOT；后续追加36个新噪声案例、每例3个优化种子的相同500次profile预算对照，也不计入23窗。原生与软件结果记录在笔记0003第19节。
'''), code('''
from pathlib import Path
import json
import pandas as pd
from IPython.display import display
correlated_catalog_project = next(p for p in (Path.cwd(), *Path.cwd().parents) if (p/'g4/data').is_dir())
correlated_catalog = json.loads((correlated_catalog_project/'jupyter/artifacts/alpha_correlated_global_20261020/root-macro-index.json').read_text())
correlated_catalog_df = pd.DataFrame(correlated_catalog['rows'])
assert len(correlated_catalog_df) == 23 and (correlated_catalog_df['status'] == 'completed').all()
display(correlated_catalog_df[['rate_hz','root','macro','seed','sipm_seed','readout_seed','status']])
''')])
    catalog_execution = execute_appended(catalog, catalog_start, 'RunCatalog-correlated')
    assert catalog.cells[:catalog_start] == catalog_before.cells
    nbformat.validate(catalog); nbformat.write(catalog, catalog_path)
    save(record/'notebook-update.json', dict(at=datetime.now().astimezone().isoformat(),
        old_cells_preserved=start, appended_cells=len(cells), appended_code_executed=executions,
        previous_catalog_cells_preserved=catalog_start, catalog_executed=catalog_execution,
        execution_scope='Only appended cells executed in fresh Python namespaces; historical cells not rerun',
        figure_ids=[f'0003-19-{i}' for i in range(1,6)]))
    print('NOTEBOOKS_UPDATED', len(notebook.cells), len(catalog.cells))


if __name__ == '__main__':
    main()
