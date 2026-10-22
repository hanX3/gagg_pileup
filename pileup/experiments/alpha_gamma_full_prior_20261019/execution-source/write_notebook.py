from pathlib import Path
import ast
import contextlib
import hashlib
import io
import json
import platform
import nbformat

ROOT=Path(__file__).resolve().parents[3];RECORD=Path(__file__).resolve().parent
PATH=ROOT/'jupyter/codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb'
original=nbformat.read(RECORD/(PATH.name+'.before'),as_version=4)
book=nbformat.read(PATH,as_version=4)
assert len(original.cells)==49 and book.cells==original.cells
artifact=ROOT/'jupyter/artifacts/alpha_gamma_full_prior_20261019'
dev=json.loads((artifact/'development/summary.json').read_text())
test=json.loads((artifact/'test/summary.json').read_text())
diagnostic=json.loads((artifact/'test/transition-diagnostics.json').read_text())
cells=[]
def md(value):cells.append(nbformat.v4.new_markdown_cell(value))
def code(value):cells.append(nbformat.v4.new_code_cell(value))

md('''<a id="alpha-gamma-full-prior-20261019"></a>
## 10. 2026-10-19：去掉显式电荷惩罚的完整混合流程与新独立测试

第9节固定候选消融提示跨类别共用单个网络Q的惩罚可能损害γ。本节把仅类别惩罚接回完整流程，恢复到达时刻搜索和残差补峰，先用20个已公开开发窗检查，再冻结并生成16个新种子原生窗口。α与γ逐事件恢复同等重要；开发与独立测试分别报告。

[实验说明与复算](../pileup/experiments/alpha_gamma_full_prior_20261019/README.md)、[开发协议](../pileup/experiments/alpha_gamma_full_prior_20261019/protocol.json)、[测试冻结设置](../pileup/experiments/alpha_gamma_full_prior_20261019/frozen-fit-settings.json)、[采用记录L31](../docs/literature-use.md#alpha-gamma-full-prior-20261019)。第1–9节及其计算输出保留。''')
md('''### 10.1 唯一的算法改动与事先固定的开发筛选

相对旧完整混合流程仅将显式Q惩罚λQ从10改为0，类别强度λC保持10；共同的神经/物理候选并集、初始电荷与Q相关模板形状、固定噪声权重及搜索预算均沿用。每次最多两轮时刻/类型/电荷搜索，单脉冲时刻搜索范围±8 ns、100 ns以内近邻沿用旧版联合时刻约束，单次非线性子问题最多20次函数评价；残差最多补8个候选。类型可以改变，数量和时间均不是已知真值。

“无显式Q惩罚”仍保留网络Q对候选和模板形状初始化的影响，不能称为完全摆脱网络Q。类别分数仍未概率校准；有限局部优化也不保证全局最优。

开发结果出来前固定筛选：两类总体恢复数均不低于原物理拟合；两类完整precision最多下降2个百分点；额外预测≤1.25×物理拟合额外预测+5；有≥20个目标的源条件/粒子组，恢复率最多下降5个百分点。这些容差是本轮工程筛选选择，不是文献给出的阈值或统计等效性界限。''')
code('''from pathlib import Path
import json
import pandas as pd
_full_project = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p/'pileup/analysis').is_dir())
_full_artifact = _full_project/'jupyter/artifacts/alpha_gamma_full_prior_20261019'
_full_dev = json.loads((_full_artifact/'development/summary.json').read_text())
_full_test = json.loads((_full_artifact/'test/summary.json').read_text())
_full_names = {'baseline':'物理拟合','unet':'U-Net','hybrid_none':'混合：无显式惩罚',
               'hybrid_both':'旧混合：Q＋类别','class_only':'新混合：仅类别'}
def _full_table(summary):
    rows=[]
    for method in summary['methods']:
        items={r['particle']:r for r in summary['pooled'] if r['method']==method}
        rows.append({'方法':_full_names[method],
            'α完整恢复':f"{items['alpha']['recovered_count']}/{items['alpha']['target_count']} ({100*items['alpha']['full_recovery']:.1f}%)",
            'γ完整恢复':f"{items['gamma']['recovered_count']}/{items['gamma']['target_count']} ({100*items['gamma']['full_recovery']:.1f}%)",
            'α完整precision/%':round(100*items['alpha']['full_precision'],1),
            'γ完整precision/%':round(100*items['gamma']['full_precision'],1),
            '额外预测':summary['diagnostics'][method]['unmatched_predictions']})
    return pd.DataFrame(rows)
_full_table(_full_dev)''')
md('''开发集新流程α505/791、γ42/86，相对物理拟合471/791、38/86；额外预测45→57，α完整precision70.3%→68.9%，γ35.2%→53.8%。通过预定筛选，但γ仅增加4个，成对整窗重采样的恢复率差95%范围为−2.44至+12.25个百分点，包含零。相对旧混合两项惩罚，α520→505，γ38→42，额外预测81→57；这仍是两类间的取舍，不能只看合并指标。

完整恢复沿用真值Q≥500 p.e.、类型无关一对一时间匹配≤10 ns、类型正确且|ΔQ/Q|≤20%，中央发射区[5,15) μs。完整precision分母是重建时刻落入中央区的该类报告预测，分子要求与正电荷真值匹配、类型及Q均正确；保留低于目标阈值与保护区真值的配对，避免将边界/阈下真实脉冲直接记为误检。因此precision分子不必与目标完整恢复数完全相同。Q为实现电荷p.e.，不是沉积能量或能量分辨率。''')
md('''### 10.2 冻结后的独立原生测试

16窗：两种γ背景（轫致辐射谱、661.657 keV）×α源1/10 MHz×每条件4窗；γ源10 MHz。源注入率与阈上探测率区分。每窗完整20 μs、2 ns采样，一个ROOT一个串行G4事件，4个独立进程；4个预检窗计入16窗，不额外计样本。与前轮完全相同的物理源码、配置和原生共享SiPM读出，足够动态范围，使用完整模拟电压，不依据ADC削顶拒绝。

U-Net权重、阈值[0.4,0.55]、扩展标定及全部拟合设置在新波形生成前固定，无重新训练。新种子261091200–261091215，与已有生产种子核对无重用。五种方法在共同新波形上计算：物理拟合、U-Net、无显式惩罚混合、两项惩罚旧混合、仅类别新混合。所有盲拟合结束后才统一读真值评分；重建子进程的文件访问检查禁止打开ROOT或truth路径。

最初启动漏建macros/logs目录，在启动G4之前失败，未产生ROOT；补齐目录后重试，原失败日志保留，未改变物理或拟合参数。另一次外层运行以143中断时已有14窗全部盲结果与访问记录，后2窗尚无结果；保持冻结代码/配置和原ROOT，仅补跑这2窗，不增加独立样本数。首次14窗的外层返回码未落盘，以完整结果和访问记录恢复；补跑日志单独保留。原始ROOT/宏/RunInfo映射见[Run Catalog](Run_Catalog.ipynb#alpha-gamma-full-prior-native-20261019)与[实际manifest](../pileup/experiments/alpha_gamma_full_prior_20261019/manifest.json)。''')
code('''_full_table(_full_test)''')
code('''pd.DataFrame([{'方法':_full_names[m],
    '目标数':_full_test['diagnostics'][m]['target_count'],
    '时间漏检':_full_test['diagnostics'][m]['missed_count'],
    '额外预测':_full_test['diagnostics'][m]['unmatched_predictions'],
    '误拆分标志':_full_test['diagnostics'][m]['split_flags'],
    '误合并标志':_full_test['diagnostics'][m]['merge_flags']} for m in _full_test['methods']])''')
code('''_full_charge=[]
for r in _full_test['pooled']:
    q=r['conditional_correct_type_charge']
    _full_charge.append({'方法':_full_names[r['method']],'粒子':r['particle'],
        '阈上目标':r['target_count'],'时间漏检':r['missed_count'],'错类型':r['wrong_type_count'],
        '条件Q统计样本':q['n'],
        '正确类型条件Q偏差/%':round(100*q['bias'],2) if q['bias'] is not None else None,
        '正确类型条件Q标准差/%':round(100*q['std'],2) if q['std'] is not None else None})
pd.DataFrame(_full_charge)''')
code('''_full_timing=json.loads((_full_artifact/'test/timing-summary.json').read_text())
pd.DataFrame([{'方法':_full_names[r['method']], '执行组':r['execution_group'], '窗口数':r['windows'],
    '本方法步骤中位耗时/s':round(r['fit_median_s'],3),
    '含前置步骤中位耗时/s':round(r['pipeline_median_s'],3),
    '含前置步骤P95耗时/s':round(r['pipeline_p95_s'],3)} for r in _full_timing['groups']])''')
md('''耗时为本机处理完整20 μs输入的观测墙钟时间，包含并行竞争；首14完成窗为4窗并行，因外层中断补跑的后2窗为2窗并行，每进程BLAS/Torch各1线程，分组展示。混合流水线含物理基线与网络步骤，不含进程/模型启动和文件导出。它是本次离线处理成本记录，不作为硬件吞吐率基准。''')
code('''from IPython.display import Image
Image(filename=str(_full_artifact/'test/independent-test-comparison.png'))''')
code('''_full_delta=[]
for r in _full_test['paired_contrasts']:
    lo,mid,hi=r['paired_window_p025_p50_p975']
    _full_delta.append({'比较':'仅类别 − '+_full_names[r['b']], '粒子':r['particle'],
        '恢复数变化':r['delta_count'],'恢复率差/百分点':round(100*r['delta_recovery'],2),
        '成对窗口95%范围/百分点':f'[{100*lo:.2f}, {100*hi:.2f}]'})
pd.DataFrame(_full_delta)''')
code('''Image(filename=str(_full_artifact/'test/independent-test-conditions.png'))''')

failed=[r['check'] for r in test['frozen_screen']['checks'] if not r['passed']]
decision=('本批新数据通过同一套筛选规则。' if not failed else '本批新数据未通过同一套筛选规则；未通过项：'+', '.join(failed)+'。')
line=[]
for particle,label in [('alpha','α'),('gamma','γ')]:
    new=next(r for r in test['pooled'] if r['method']=='class_only' and r['particle']==particle)
    old=next(r for r in test['pooled'] if r['method']=='baseline' and r['particle']==particle)
    transition=next(r for r in diagnostic['transitions'] if r['reference']=='baseline' and r['particle']==particle)
    line.append(f"{label}完整恢复{old['recovered_count']}/{old['target_count']}→{new['recovered_count']}/{new['target_count']}；逐事件损失{transition['lost']}个原成功事件，补回{transition['gained']}个。")
md('### 10.3 独立结果与取舍\n\n'+decision+' '.join(line)+
   ' 不能把净增加解释为每个事件都改善，也不能从小样本推断普遍有效。以下按原条件分层的2000次成对整窗bootstrap只描述窗口波动，不含模板、器件与物理模型系统误差。测试结果未用于改阈值或重新选模型。')
code('''_full_diag=json.loads((_full_artifact/'test/transition-diagnostics.json').read_text())
_full_strata=[]
for group in _full_diag['strata']:
    row={'γ分层':group['group'],'目标数':group['n']}
    for method in ('baseline','hybrid_both','class_only'):
        counts=group['outcomes'][method]
        row[_full_names[method]+'完整恢复']=counts.get('recovered',0)
        row[_full_names[method]+'错类型']=counts.get('wrong_type',0)
    _full_strata.append(row)
pd.DataFrame(_full_strata)''')
md('''### 10.4 追溯、验证、限制与下一步

新优化器的3项针对性数值检查包括与旧强度0/10端点的复现、浮动时刻与残差补峰，以及非法惩罚/采样拒绝。首个真实开发窗的旧端点、重新运行冻结U-Net的密集图和端到端五方法均数值一致；新测试保存逐窗输入哈希、访问记录、目标历史及优化停止状态。部分近邻子问题达到20次函数评价预算；目标下降不等同于全局收敛。实际验证范围以[validation.json](run_records/alpha_gamma_full_prior_20261019/validation.json)为准。

每源条件仅4窗，γ数量有限；同一物理模型下的新种子并不验证真实器件域偏移。未完成可见电荷到沉积能量的条件标定，不以这些Q误差报告能谱分辨率。源码、原宏、配置、模型哈希与独立输出保留；旧ROOT与旧结果未覆盖。Notebook只执行本节新增的汇总代码，旧49个单元及输出保留；查看本节不需要Torch，重新推理需要固定CPU Torch依赖。

下一步仍应针对剩余的弱γ错分和强前史/近邻失败，先在独立标定中检验类型、电荷与可观测前史下的模板/Q偏差及不确定度，再决定Q条件模板迭代或类型条件软约束。不能直接把两个未充分监督的Q图当成已标定先验，也不继续用本批测试调参。文献条件与审核状态沿用L28–L30；本节是项目计算结果，未新增论文性能结论或修改共享文献库。''')

# The common execution block is appended from the previously validated writer.
namespace={'__name__':'__full_prior_notebook__'}
count=max((c.get('execution_count') or 0 for c in original.cells if c.cell_type=='code'),default=0)
executed=0
for cell in cells:
    if cell.cell_type!='code':continue
    count+=1;executed+=1
    tree=ast.parse(cell.source);last=tree.body.pop() if tree.body and isinstance(tree.body[-1],ast.Expr) else None
    stdout=io.StringIO();stderr=io.StringIO()
    with contextlib.redirect_stdout(stdout),contextlib.redirect_stderr(stderr):
        exec(compile(tree,'<full-prior-cell>','exec'),namespace)
        value=eval(compile(ast.Expression(last.value),'<full-prior-expression>','eval'),namespace) if last else None
    outputs=[]
    if stdout.getvalue():outputs.append(nbformat.v4.new_output('stream',name='stdout',text=stdout.getvalue()))
    if stderr.getvalue():outputs.append(nbformat.v4.new_output('stream',name='stderr',text=stderr.getvalue()))
    if value is not None:
        data={'text/plain':repr(value)};metadata={}
        if hasattr(value,'_repr_mimebundle_'):
            bundle=value._repr_mimebundle_()
            rich,metadata=bundle if isinstance(bundle,tuple) else (bundle,{})
            data.update(rich)
        elif hasattr(value,'_repr_html_'):
            html=value._repr_html_()
            if html:data['text/html']=html
        outputs.append(nbformat.v4.new_output('execute_result',data=data,metadata=metadata,execution_count=count))
    cell.execution_count=count;cell.outputs=outputs
book.cells.extend(cells)
book.metadata['full_prior_execution']={'runtime':platform.python_version(),'mode':'New cell Python code and final expressions evaluated in project runtime; old cells preserved without reexecution'}
nbformat.validate(book)
nbformat.write(book,PATH)
after=nbformat.read(PATH,as_version=4)
assert after.cells[:len(original.cells)]==original.cells
result={'original_cells':len(original.cells),'appended_cells':len(cells),'new_code_cells_executed':executed,
        'embedded_images':sum('image/png' in o.get('data',{}) for c in cells for o in c.get('outputs',[])),
        'all_original_cells_outputs_preserved':True,'notebook_sha256':hashlib.sha256(PATH.read_bytes()).hexdigest()}
(RECORD/'notebook-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
