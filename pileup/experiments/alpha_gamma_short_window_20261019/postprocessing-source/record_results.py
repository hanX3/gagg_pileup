"""Append executed local research records while preserving all prior cells."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys

import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output

project=Path(__file__).resolve().parents[3]
record=Path(__file__).resolve().parent
output=project/'jupyter/artifacts/alpha_gamma_short_window_implementation_20261019'
sys.path.insert(0,str(project))
from analysis.alpha_gamma_short_window_report import selected_group

summary=json.loads((output/'test-summary.json').read_text())
report=json.loads((output/'report.json').read_text())
settings=json.loads((output/'frozen-settings.json').read_text())
validation=json.loads((record/'native-data-validation.json').read_text())
assert summary['independent_events']==174 and validation['status']=='passed'

def pct(k,n):
    return f'{k}/{n}（{100*k/n:.1f}%）' if n else '无目标'

def table(pre=20.,methods=('physics_ols','physics_gls','template_ols','template_gls')):
    rows=['| 触发后T / ns | 方法 | α电荷误差≤5% | γ电荷误差≤5% |',
          '| ---: | --- | ---: | ---: |']
    for post in (100.,200.,300.,500.,1000.):
        for method in methods:
            a,g=[selected_group(summary,method,p,pre,post) for p in ('alpha','gamma')]
            rows.append(f"| {post:g} | {method} | {pct(a['charge_5pct_success'],a['targets'])} | {pct(g['charge_5pct_success'],g['targets'])} |")
    return '\n'.join(rows)

choice=report['development_choice']
chosen_pre=float(choice['chosen']['pre_ns']);chosen_post=float(choice['chosen']['post_ns'])
chosen=[selected_group(summary,'template_gls',p,chosen_pre,chosen_post) for p in ('alpha','gamma')]
choice_pass=all(g['targets'] and g['charge_5pct_fraction']>=.95 for g in chosen)
choice_statement=('通过本批点比例门槛' if choice_pass else '未通过本批点比例门槛')
primary_grid=[r for r in report['earliest_grid_points'] if r['method']=='template_gls']
grid_statement='；'.join(f"前{r['pre_ns']:g} ns："+(f"后{r['earliest_post_ns']:g} ns" if r['earliest_post_ns'] is not None else '所有已测档位均未满足') for r in primary_grid)

anchor='short-window-implementation-20261019'
nb_path=project/'jupyter/codex_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb'
notebook=nbformat.read(nb_path,as_version=4)
before=nbformat.read(record/'codex_GAGG_SiPM_Alpha_Gamma_Fitting-before.ipynb',as_version=4)
assert len(before.cells)==50 and notebook.cells==before.cells, 'Do not append twice or replace existing work'
shell=InteractiveShell.instance()
shell.run_cell("from pathlib import Path\nimport os\nos.chdir("+repr(str(project))+')',store_history=False,silent=True)
execution=max([int(c.get('execution_count') or 0) for c in notebook.cells]+[0])
new_cells=[]

def md(text):
    new_cells.append(nbformat.v4.new_markdown_cell(text))

def code(source):
    global execution
    execution+=1
    with capture_output() as captured:
        result=shell.run_cell(source,store_history=False,silent=False)
    if result.error_before_exec or result.error_in_exec:
        raise RuntimeError(result.error_before_exec or result.error_in_exec)
    outputs=[]
    if captured.stdout:outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
    if captured.stderr:outputs.append(nbformat.v4.new_output('stream',name='stderr',text=captured.stderr))
    for item in captured.outputs:
        outputs.append(nbformat.v4.new_output('display_data',data=item.data,metadata=item.metadata))
    new_cells.append(nbformat.v4.new_code_cell(source,execution_count=execution,outputs=outputs))

md(f'''<a id="{anchor}"></a>
## 10. 实现与独立测试：用短前缀估计孤立α/γ的完整电荷

2026-10-19，落实第9节的第一阶段。研究问题是：**在完整形状或衰减基底由独立标定约束时，只看有限的波形前缀，能否把完整电荷Q估到±5%？** 主指标为两类各自≥95%的目标成功率，阈值与扫描网格在新测试前冻结。这里的Q是光子后代雪崩的实现电荷，以完全恢复单微胞电荷为1 pe；**不是源能量、沉积能量或已标定的可见能量，也不是逐事件自由双τ拟合**。

旧第5/8节仅有1α＋1γ、ADC输入、与同事件长窗解比较，且使用长记录定位/基线。本轮另用完整原生含噪模拟电压、因果触发、独立实现Q真值；所有四方法使用相同输入。因此不能把新旧表直接称作“300 ns缩短到某值”的成对提升。

[可复用实现](../analysis/alpha_gamma_short_window.py)、[数据划分和运行入口](../analysis/alpha_gamma_short_window_study.py)、[冻结设置](artifacts/alpha_gamma_short_window_implementation_20261019/frozen-settings.json)、[实验快照](../pileup/experiments/alpha_gamma_short_window_20261019/README.md)、[采用记录L35](../docs/literature-use.md#{anchor})。''')

md('''### 10.1 数据、标定和盲输入

按ROOT/随机种子分组，复用149个独立孤立标定窗与45个开发窗；同一事件的全部截窗始终在同一组。149个标定窗中Q≥250 pe且有起点的70个（α52、γ18）用于三电荷档的类别均值和残差协方差，所有149个用于独立静态基线/电子噪声估计。每个γ形状档仅6例；这是明显的标定限制。基线均值来自标定，基线软约束σB=0.2 mV、协方差50%对角收缩为事先设定的工程参数，不是已验证的实物漂移分布。726维时间协方差来自70条样本，不能视作精确噪声模型。

新测试在代码、标定和自适应阈值冻结后生成，共174个原生20 μs窗口：30α、144γ。α源能量1.5/1.8/2.5/5.486/8 MeV，γ为0.1/0.2/0.661657/1.5 MeV；位置(x,y)=(0,0)/(4,0)/(4,4) mm，中央正入射，单初级粒子在5 μs发射。γ是四个单能控制，**不是两种连续背景谱的高率实验**。能量和位置配额是预定测试设计，不能据此计算真实源谱的效率。

推断只接收采样时刻、电压和冻结标定，未知类型/时刻/Q；**显式假定至多一个孤立脉冲，不估计脉冲数**。标定时允许使用独立的类型、Q和首个有效光子时刻；测试真值在该事件44个窗口的全部拟合结束后才打开。Q≥500 pe仅用于评价分母，不用于选波形或输入拟合。所有零信号γ保留。原生共享微胞、读出和物理配置保持冻结版本，使用analog_mV，不作ADC触顶拒绝。''')

code("""import json
from pathlib import Path
import pandas as pd
from IPython.display import display, Image
prefix_project=next(p for p in (Path.cwd(),*Path.cwd().parents) if (p/'analysis').is_dir() and (p/'g4').is_dir())
prefix_output=prefix_project/'jupyter/artifacts/alpha_gamma_short_window_implementation_20261019'
prefix_summary=json.loads((prefix_output/'test-summary.json').read_text())
prefix_report=json.loads((prefix_output/'report.json').read_text())
prefix_validation=json.loads((prefix_project/'jupyter/run_records/alpha_gamma_short_window_implementation_20261019/native-data-validation.json').read_text())
display(pd.DataFrame(prefix_summary['truth_counts']).T.rename(columns={'injected':'注入窗口','above_500':'评价Q≥500 pe','zero_charge':'光子后代Q=0'}))
print('独立新窗口:',prefix_summary['independent_events'],'；ROOT/实际宏/种子/祖先电荷检查:',prefix_validation['new_native_roots'])
print('冻结前开发选择:',prefix_report['development_choice'])
""")

md('''### 10.2 四种方法与窗口定义

- `physics_ols`：原95/351 ns固定衰减基底，与SPE、前端和采样响应卷积；两个非负幅度Af/As、基线B，profile搜索起点。快慢比例仍自由，无独立类别输出。
- `physics_gls`：相同模型，采用独立原生残差协方差加权。
- `template_ols`：分别拟合独立标定的α、γ完整平均形状，每类一个Q幅度，比较残差选择类型。
- `template_gls`（预定主方法）：类别模板＋协方差加权。两类使用同一个由电压pilot估计的Q/t0条件协方差；本轮不是参数变化下归一化的高斯似然。模板的Q插值也使用pilot，不输入真值Q。

两幅度求解枚举非负活跃集合、消去基线，起点在有限网格加局部搜索。协方差**先截当前窗口，再分解/求逆**。OLS和GLS的基线先验使用一致的电压单位。参数改写本身不增加信息；模板方法额外约束了完整形状，不能把其改善全归因于求解器。

触发是相对过去滚动中位数连续3个采样超过2 mV的第三个采样，2 ns采样的确认延迟4 ns包含在Tpre中；不读取全记录最大值。拟合只读`[观测触发−Tpre, 观测触发+Tpost]`。Tpre=20/50/100/200 ns，Tpost=50/75/100/150/200/250/300/400/500/750/1000 ns；75 ns在2 ns网格实际止于74 ns。**Tpost从观测触发算，既不是全局0 ns，也不是已知真实事件起点。** Tpre+Tpost是拟合数据跨度，不是CPU运行时间或完整采集记录长度；没有把数字化前端的延迟抹掉。未来电压替换为NaN以及全局时刻平移的原生重放均通过。''')

code("""prefix_rows=[]
for _post in (100.,200.,300.,500.,1000.):
    for _method in ('physics_ols','physics_gls','template_ols','template_gls'):
        _row={'Tpre/ns':20,'Tpost/ns':int(_post),'method':_method}
        for _particle in ('alpha','gamma'):
            _g=next(g for g in prefix_summary['groups'] if g['pre_ns']==20 and g['post_ns']==_post and g['method']==_method and g['particle']==_particle and g['band']=='all_above_500')
            _row[_particle+' Q±5%']=f"{_g['charge_5pct_success']}/{_g['targets']} ({100*_g['charge_5pct_fraction']:.1f}%)"
        prefix_rows.append(_row)
display(pd.DataFrame(prefix_rows))
""")
code("display(Image(filename=str(prefix_output/'prefix_success.png'),width=1100))")
code("display(Image(filename=str(prefix_output/'prefix_errors.png'),width=1100))")
code("display(Image(filename=str(prefix_output/'prefix_window_map.png'),width=1100))")

md(f'''### 10.3 独立结果的解释

开发集Q≥500 pe仅18α＋4γ。以两类各自至少95%的Q±5%成功且本点及全部更长已测窗口均合格为规则，开发主方法四个Tpre均从Tpost=300 ns起合格；开发确定较短的**前{chosen_pre:g}＋后{chosen_post:g}={chosen_pre+chosen_post:g} ns**作为待验证候选。该候选选择文件在新测试评分前写入；部分原生预检查ROOT在文件写入前已生成，但未用来选窗。

在全部新独立测试中，该候选**{choice_statement}**：α {pct(chosen[0]['charge_5pct_success'],chosen[0]['targets'])}，γ {pct(chosen[1]['charge_5pct_success'],chosen[1]['targets'])}。这是固定候选的测试结论，未据测试重调形状、协方差、阈值或算法。

前20＋后100 ns时，主方法类型正确为α14/19、γ9/12，电荷合格为14/19、7/12；短窗失败同时包含类型混淆与正确类型下的Q超差。200 ns的α18/19（94.74%）低于预定95%点比例门槛，250 ns又降为17/19，因此不能只挑200 ns的较好结果宣布足够。300 ns起，两类的Q和类型同时全合格。这里的错误分解只描述本批测试，未用于再次选参。

整张测试网格的探索性描述为：{grid_statement}。这些是本批样本的离散扫描结果，同事件截窗高度相关；即使点比例≥95%，也不代表真实总体成功率已用95%置信度证明≥95%，更不是连续意义下精确最短窗。图中Wilson区间只描述逐条件有限事件统计，不含标定、物理模型和多窗口选择误差。

电荷偏差/散布/90%绝对误差为匹配成功事件的条件统计；Q±5%成功率以全部Q≥500 pe目标为分母，漏检保留为失败。下一表同时列出分型及“Q和类型都正确”，避免把电荷准确率当成粒子逐事件完整恢复。两幅度物理基线没有分类器，不给其凭空记分类失败。''')

code("""prefix_detail=[]
for _post in (100.,200.,300.,500.,1000.):
    for _particle in ('alpha','gamma'):
        _g=next(g for g in prefix_summary['groups'] if g['pre_ns']==20 and g['post_ns']==_post and g['method']=='template_gls' and g['particle']==_particle and g['band']=='all_above_500')
        prefix_detail.append({'Tpost/ns':int(_post),'particle':_particle,'targets':_g['targets'],'匹配':_g['detected'],'类型正确':_g['correct_type'],'Q和类型均正确':_g['joint_charge_type_success'],'条件偏差/%':100*_g['charge_bias'],'条件相对误差SD/%':100*_g['charge_error_sd'],'Q±5% Wilson95':_g['charge_5pct_wilson95']})
display(pd.DataFrame(prefix_detail))
display(pd.DataFrame([g for g in prefix_summary['groups'] if g['method']=='template_gls' and g['pre_ns']==20 and g['post_ns']==300 and g['band']!='all_above_500'])[['particle','band','targets','detected','charge_5pct_success','charge_bias','charge_error_sd']])
print('Q=0事件上的触发及误报（保留全部零信号）：')
display(pd.DataFrame(prefix_report['false_triggers_on_no_signal']).T)
""")

md('''### 10.4 自适应停止和双脉冲负对照

自适应规则在开发集冻结：Tpre=50 ns，首次满足`k·局部σQ/Q≤5%`、残差门槛及无起点触边时停止。k取每个开发事件跨窗口最坏标准化误差的类内95%分位数，再取两类最大，为54.4515；残差上限2.106。局部曲率不包含类型切换/标定不确定度，开发错误导致k很大，因此规则可能几乎不接受事件。**不因测试中接受率低而降低k**；拒绝/未停止仍计入无条件失败。这只是已有前缀的离线重放，没有验证顺序停止后的置信覆盖率或真实实时执行能力。

另有24条合成受控双脉冲（两顺序、50/100/200/400 ns间隔、3种幅度比），各看100/200/300 ns，共72次前缀检查。它们由标定均值相加＋白噪声生成，**不是新的G4共享微胞双脉冲**，不用于调参或推断实物分辨极限；只用来检查“单脉冲短窗算法”遇到邻居后的失效，不能将它当作已经完成的盲解pile-up。''')

code("""display(pd.DataFrame(prefix_report['adaptive']).T)
prefix_controls=json.loads((prefix_output/'controlled-doubles.json').read_text())
print('合成双脉冲前缀检查:',len(prefix_controls['rows']))
_double_rows=[]
for _r in prefix_controls['rows']:
    for _method,_f in _r['fit']['methods'].items():
        _double_rows.append({'gap/ns':_r['gap_ns'],'Tpost/ns':_r['post_ns'],'method':_method,'邻居进入窗口':_r['second_signal_can_enter_prefix'],'第一Q误差/%':100*_f['relative_first_charge_error'],'第一Q±5%':abs(_f['relative_first_charge_error'])<=.05})
display(pd.DataFrame(_double_rows).groupby(['gap/ns','Tpost/ns','method','邻居进入窗口'],as_index=False).agg(样本数=('第一Q±5%','size'),第一Q合格数=('第一Q±5%','sum'),第一Q平均误差百分比=('第一Q误差/%','mean')))
""")

md('''### 10.5 当前结论、限制与后续

上表直接回答固定完整响应下的孤立电荷精度；**不能回答逐事件free双τ的最短窗**。固定95/351 ns仍自由快慢幅度与固定α/γ全形状也是不同的先验强度；若同时要求ws或两τ准确，必须另外制定并验证目标。当前没有改动GAGG光产额、衰减时间、Birks或SiPM硬件参数，没有将Q精度称为能量分辨率。

本轮只覆盖少数源能量/三点位置，γ有效样本与标定均有限；模板区间外的Q会夹在端点形状。协方差使用混合样本、pilot固定和工程收缩，未做独立噪声模型校验。原生共有物理模型和“同域”标定仍会低估真实器件/温漂/位置外推风险。无时间线性漂移项、邻居/暗脉冲数估计，也未验证高率下共享微胞前史。漏检、假触发、错类型均可报告；由于至多一脉冲假设，不能从本轮给出盲多脉冲的误拆分/误合并效率。

下一阶段先在独立开发数据上补有效γ、条件Q/前史和标定不确定度，检验固定候选失败原因；如要改变方法，需另设冻结后测试。再做原生受控双脉冲与混合高率联合拟合，并分别比较自由τ与跨事件共享τ。实际探测事件率与源注入率分开，不能用1/(Tpre+Tpost)作为最大计数率。

[完整逐窗结果](artifacts/alpha_gamma_short_window_implementation_20261019/test-summary.json)、[报告](artifacts/alpha_gamma_short_window_implementation_20261019/report.json)、[174个新ROOT索引](Run_Catalog.ipynb#short-window-native-20261019)、[数据与重放验证](run_records/alpha_gamma_short_window_implementation_20261019/native-data-validation.json)。原50单元与全部原有计算/图像保持；新增计算和图表为本轮实际执行。''')

notebook.cells.extend(new_cells)
nbformat.validate(notebook)
nbformat.write(notebook,nb_path)

catalog_path=project/'jupyter/Run_Catalog.ipynb'
catalog=nbformat.read(catalog_path,as_version=4)
old_catalog=nbformat.read(record/'Run_Catalog-before.ipynb',as_version=4)
assert catalog.cells==old_catalog.cells
catalog.cells.append(nbformat.v4.new_markdown_cell('''<a id="short-window-native-20261019"></a>
## 2026-10-19：因果短窗口拟合的174个独立孤立测试

30α＋144γ，代码/标定/自适应阈值冻结后新种子生成；包含6个预检查，不重复计数。每ROOT一个20 μs事件、5 μs发射一个初级粒子，三点位置，2 ns模拟电压采样。单能γ四档0.1/0.2/0.661657/1.5 MeV，α五档1.5/1.8/2.5/5.486/8 MeV；这里MeV为源能量、不是沉积/可见能量。每次G4内部串行、最多4个独立任务、实际启动至少间隔1秒；不是Poisson源率扫描，清单rate_hz=0只是未使用的占位字段。全部零沉积γ保留。

与独立旧149标定/45开发按ROOT和种子隔离。物理源码/二进制、宏、配置与日志分别冻结在[主实验](../pileup/experiments/alpha_gamma_short_window_20261019/README.md)和[预定γ扩展](../pileup/experiments/alpha_gamma_short_window_gamma_extension_20261019/README.md)，两个目录必须合并成174窗使用。所有RunInfo/实际宏/随机种子/配置、ROOT哈希、祖先实现电荷已逐文件核验，见[数据验证](run_records/alpha_gamma_short_window_implementation_20261019/native-data-validation.json)。拟合输入完整analog_mV、动态范围足够，不拒绝ADC触顶。

下表来自各次实际manifest及已核RunInfo；truth Q仅作评价。完整算法/结果见[主题Notebook第10节](codex_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-implementation-20261019)。'''))
catalog_source="""from pathlib import Path
import json
import pandas as pd
from IPython.display import display
short_catalog_project=next(p for p in (Path.cwd(),*Path.cwd().parents) if (p/'g4/data').is_dir())
short_catalog_rows=[]
for _name in ('alpha_gamma_short_window_20261019','alpha_gamma_short_window_gamma_extension_20261019'):
    _runs=json.loads((short_catalog_project/'pileup/experiments'/_name/'manifest.json').read_text())['runs']
    for _r in _runs:
        _stem=Path(_r['root']).stem
        _q=json.loads((short_catalog_project/'jupyter/artifacts/alpha_gamma_short_window_implementation_20261019/test'/(_stem+'.json')).read_text())['truth_for_evaluation_only']['charge_pe']
        short_catalog_rows.append({k:_r[k] for k in ('root','macro','particle','energy_MeV','x_mm','y_mm','seed','sipm_seed','readout_seed','status')}|{'Q_truth_pe':_q})
assert len(short_catalog_rows)==174
with pd.option_context('display.max_rows',180,'display.max_colwidth',105):
    display(pd.DataFrame(short_catalog_rows))
"""
with capture_output() as cap:
    result=shell.run_cell(catalog_source,store_history=False)
if result.error_in_exec or result.error_before_exec:raise RuntimeError(result.error_in_exec or result.error_before_exec)
catalog_outputs=[]
if cap.stdout:catalog_outputs.append(nbformat.v4.new_output('stream',name='stdout',text=cap.stdout))
for item in cap.outputs:catalog_outputs.append(nbformat.v4.new_output('display_data',data=item.data,metadata=item.metadata))
catalog.cells.append(nbformat.v4.new_code_cell(catalog_source,execution_count=max(int(c.get('execution_count') or 0) for c in catalog.cells)+1,outputs=catalog_outputs))
nbformat.validate(catalog);nbformat.write(catalog,catalog_path)

literature=project/'docs/literature-use.md'
original=(record/'literature-use-before.md').read_text()
assert literature.read_text()==original
literature.write_text(original+f'''

<a id="{anchor}"></a>
### L35：2026-10-19——因果短前缀拟合的实现与独立原生检验

**文献依据与采用边界**：落实[L34](#short-window-literature-strategy-20261019)的有限窗补尾、线性幅度profile和独立条件均值/协方差原则；Belli NE213固定7.5/45 ns、Radeka理想噪声、Fixsen TES的装置/采样/模型条件与AI草稿审核状态均沿L34，不引入新的文献性能或人工确认。Radeka/Nakhostin原式疑点未消除，未照抄缩尾公式；本轮没有修改共享papers库。α/γ完整模板、条件pilot、50%收缩和停止规则为本项目实现，不冒充论文给出的GAGG性能。

**实现与独立性**：149个旧孤立标定窗、45个旧开发窗，70个Q≥250 pe（52α/18γ）用于三档类别均值/726维协方差；基线均值和电子噪声独立标定，σB=0.2 mV及50%对角收缩为工程设定。四方法是固定95/351 ns两幅度OLS/GLS及独立完整α/γ形状OLS/GLS，输入同一完整原生analog_mV。触发用过去滚动中位数的连续三点2 mV门槛，窗口只读[触发−Tpre,触发+Tpost]；先截协方差再分解，未使用未来峰值/尾部定位。固定类别完整形状比固定双τ且自由比例约束更强，参数改写本身不提高信息量。

配置、标定、4个核心源码哈希冻结后，新增174个独立原生20 μs单事件ROOT（30α/144γ），预定包含γ扩展并保留零信号；单能γ0.1/0.2/0.661657/1.5 MeV、α1.5/1.8/2.5/5.486/8 MeV，三点位置。真值只在该事件44个窗口全部拟合完成后打开，不作为推断输入。Q≥500 pe目标数为α{summary['truth_counts']['alpha']['above_500']}、γ{summary['truth_counts']['gamma']['above_500']}；这是实现光子后代雪崩电荷，不是源/沉积/可见能量。

**实际结果**：主方法template_gls。开发选择的前{chosen_pre:g}＋后{chosen_post:g} ns候选在新测试{choice_statement}（α{pct(chosen[0]['charge_5pct_success'],chosen[0]['targets'])}，γ{pct(chosen[1]['charge_5pct_success'],chosen[1]['targets'])}）；目标为两类各≥95%事件Q误差≤5%。测试网格仅作探索描述：{grid_statement}。小样本点比例和逐条件Wilson区间不能认证总体95%门槛；同事件多个截窗不当作独立样本，未据测试重选算法/停止阈值。前20 ns的关键比较如下。

{table()}

**停止与负对照**：开发最坏窗口误差得到k={settings['adaptive_multiplier']:.5f}，残差阈值{settings['adaptive_residual_limit']:.5f}，固定前50 ns；新测试α接受{report['adaptive']['alpha']['accepted']}/{report['adaptive']['alpha']['targets']}、γ接受{report['adaptive']['gamma']['accepted']}/{report['adaptive']['gamma']['targets']}，拒绝保留为无条件失败。局部σ和经验放大不能视作已校准顺序置信区间。24个标定均值相加的双脉冲（72次前缀检查）只检查单脉冲假设遇邻居的失效，不是新原生共享微胞双脉冲恢复验证。

**验证与限制**：4项有独立数值参照/未来污染的算法检查通过；174个ROOT的实际RunInfo、宏、种子、配置、哈希与祖先电荷核验，1408个汇总组独立重评分，4次原生未来NaN/平移检查及2次快照源码重放通过。[Notebook第10节](../jupyter/codex_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#{anchor})保留原50单元和输出；[实验](../pileup/experiments/alpha_gamma_short_window_20261019/README.md)、[测试汇总](../jupyter/artifacts/alpha_gamma_short_window_implementation_20261019/test-summary.json)、[数据验证](../jupyter/run_records/alpha_gamma_short_window_implementation_20261019/native-data-validation.json)。本轮假定至多一个孤立脉冲、不估计数量，尚不能报告盲多脉冲误拆/误并，也不回答free双τ最短窗。旧两条ADC的长窗数值稳定与本轮独立模拟真值精度口径不同，不直接比较缩短倍数。待补有效γ/条件前史/标定不确定度、原生双脉冲、跨事件共享τ与freeτ、实物域偏移及高率联合拟合；不以1/T替代计数率能力。新文献整理仍在papers项目处理。
''')

names=('alpha_gamma_short_window_20261019','alpha_gamma_short_window_gamma_extension_20261019')
for name in names:
    folder=project/'pileup/experiments'/name
    shutil.copy2(record/'frozen_replay.py',folder/'replay_prefix.py')
    post=folder/'postprocessing-source';post.mkdir(exist_ok=True)
    for source in (project/'analysis/alpha_gamma_short_window_report.py',project/'analysis/alpha_gamma_short_window_controls.py',record/'validate_results.py',Path(__file__),record/'frozen_replay.py'):
        shutil.copy2(source,post/source.name)
    manifest=json.loads((folder/'manifest.json').read_text())
    example=manifest['runs'][0]['root']
    (folder/'README.md').write_text(f'''# Independent causal-prefix isolated-pulse benchmark

2026-10-19. This directory preserves {len(manifest['runs'])} native single-event runs. Combine it with the other `alpha_gamma_short_window*_20261019` directory: the complete test is 30 alpha + 144 gamma windows. Algorithms and calibration were frozen before new generation. These are single-particle controls, not a source-rate or blind multi-pulse experiment.

Each ROOT is one 20 us native Geant4/SiPM window, one primary at 5 us. Native jobs are serial internally, at most four independent jobs, starts at least one second apart. Files remain under project `g4/data/`, named to seconds. All zero-signal gamma windows are retained. Full analog voltage is used with adequate acquisition range; historical ADC clipping does not reject a window.

`manifest.json` records actual macro/ROOT paths, hashes, seeds, positions, source energies, status, and run validation. Source energy is not deposited or visible energy. Actual physical C++ sources/configs, binary, Python sources, macros and logs are retained here; the original physical model is unchanged. `analysis-freeze.json` records frozen source and calibration hashes. External Geant4/ROOT/Python dependencies must still be available.

From project root, the existing native runner resumes missing planned runs without overwriting completed ROOT files:

```bash
bash tools/run_gagg_python.sh pileup/experiments/{name}/run_frozen.py test_single --folder pileup/experiments/{name} --inputs jupyter/artifacts/alpha_gamma_short_window_implementation_20261019/native-inputs
```

To replay a prefix with saved inference sources, after obtaining the local calibration and exported waveform arrays:

```bash
bash tools/run_gagg_python.sh pileup/experiments/{name}/replay_prefix.py --experiment pileup/experiments/{name} --root {example} --pre 20 --post 100 --output /tmp/gagg-prefix-{name}.json
```

`replay_prefix.py` checks the saved module and calibration hashes, uses `source/` and avoids dependence on the current analysis code. The destination must be new. The generated `analyze_frozen.py` preparation helper is not the standalone replay entry: process workers can derive a different project root. Use the checked serial replay above for snapshot reproduction. Frozen sources remain untouched; later reporting/validation sources are under `postprocessing-source/`.

Results and limits are in the local [research notebook](../../../jupyter/codex_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-implementation-20261019), [complete report](../../../jupyter/artifacts/alpha_gamma_short_window_implementation_20261019/report.json), and [ROOT validation](../../../jupyter/run_records/alpha_gamma_short_window_implementation_20261019/native-data-validation.json). The combined list is `test-input-records-combined.json`, not the per-campaign `test-input-records.json`.

Four voltage-only fits assume at most one pulse. Primary results use Q >= 500 pe and charge error <=5%; missed signals remain failures. Q is realized photon-descendant avalanche charge, not energy. Finite test fractions do not certify a universal shortest window. This study fixes response shape information and does not measure the minimum window for free decay constants or a maximum count rate.

Git HEAD at preparation: `{manifest['git_head']}` with uncommitted local analysis changes; full snapshots and hashes preserve those changes. Python/dependency versions are in local `freeze-complete.json`. ROOT, notebooks, artifacts, and run records are local and not included in a push. An optional `/tmp` binary cache attempt failed at symbolic-link creation; original binaries remained in place and all production ran normally.
''')
    hashes={str(p.relative_to(folder)):hashlib.sha256(p.read_bytes()).hexdigest() for p in post.glob('*.py')}
    hashes['replay_prefix.py']=hashlib.sha256((folder/'replay_prefix.py').read_bytes()).hexdigest()
    (folder/'postprocessing-sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')

save_record={'created_utc':datetime.now(timezone.utc).isoformat(),'old_notebook_cells':50,
             'new_notebook_cells':len(new_cells),'old_catalog_cells':len(old_catalog.cells),
             'new_catalog_cells':2,'development_choice_passed_test_point_target':choice_pass}
(record/'recording-complete.json').write_text(json.dumps(save_record,ensure_ascii=False,indent=2)+'\n')
print('RECORDED_PREFIX_RESEARCH',save_record,flush=True)
